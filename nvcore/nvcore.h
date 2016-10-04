// This code is in the public domain -- Ignacio Castaño <castano@gmail.com>
#pragma once

#include <vector>
#include <assert.h>
#include <stdarg.h> // va_list
#include <stdint.h>
#include <stdio.h>

#ifdef _MSC_VER
// Ignore gcc attributes.
#define __attribute__(X)
#define restrict
#define NV_FORCEINLINE __forceinline
#else
#define restrict __restrict__
#define NV_FORCEINLINE  inline __attribute__((always_inline))
#endif

#define nvCheck(exp)     if (!(exp)) { nvDebugPrint("%s %s %s\n", #exp, __FILE__, __LINE__); }
#define nvDebugCheck(exp) assert(exp)
#define nvDebug(...)    nvDebugPrint(__VA_ARGS__)
void nvDebugPrint( const char *msg, ... ) __attribute__((format (printf, 1, 2)));

// Just in case. Grrr.
#undef min
#undef max

#define NV_UINT32_MAX   0xffffffff
#define NV_FLOAT_MAX    3.402823466e+38F

namespace nv {

/// Swap two values.
template <typename T>
inline void swap(T &a, T &b)
{
	T temp(a);
	a = b;
	b = temp;
}

/// Return the maximum of the two arguments. For floating point values, it returns the second value if the first is NaN.
template <typename T>
//inline const T & max(const T & a, const T & b)
inline T max(const T &a, const T &b)
{
	return (b < a) ? a : b;
}

/// Return the maximum of the three arguments.
template <typename T>
//inline const T & max3(const T & a, const T & b, const T & c)
inline T max3(const T &a, const T &b, const T &c)
{
	return max(a, max(b, c));
}

/// Return the minimum of two values.
template <typename T>
//inline const T & min(const T & a, const T & b)
inline T min(const T &a, const T &b)
{
	return (a < b) ? a : b;
}

/// Return the maximum of the three arguments.
template <typename T>
//inline const T & min3(const T & a, const T & b, const T & c)
inline T min3(const T &a, const T &b, const T &c)
{
	return min(a, min(b, c));
}

/// Clamp between two values.
template <typename T>
//inline const T & clamp(const T & x, const T & a, const T & b)
inline T clamp(const T &x, const T &a, const T &b)
{
	return min(max(x, a), b);
}

inline float saturate(float f)
{
	return clamp(f, 0.0f, 1.0f);
}

/** Return the next power of two.
* @see http://graphics.stanford.edu/~seander/bithacks.html
* @warning Behaviour for 0 is undefined.
* @note isPowerOfTwo(x) == true -> nextPowerOfTwo(x) == x
* @note nextPowerOfTwo(x) = 2 << log2(x-1)
*/
inline uint32_t nextPowerOfTwo(uint32_t x)
{
	nvDebugCheck( x != 0 );
#if 1	// On modern CPUs this is supposed to be as fast as using the bsr instruction.
	x--;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	return x + 1;
#else
	uint32_t p = 1;
	while ( x > p ) {
		p += p;
	}
	return p;
#endif
}

inline uint64_t nextPowerOfTwo(uint64_t x)
{
	nvDebugCheck(x != 0);
	uint32_t p = 1;
	while (x > p) {
		p += p;
	}
	return p;
}

// Simple bit array.
class BitArray
{
public:

	BitArray() {}
	BitArray(uint32_t sz)
	{
		resize(sz);
	}

	uint32_t size() const
	{
		return m_size;
	}
	void clear()
	{
		resize(0);
	}

	void resize(uint32_t new_size)
	{
		m_size = new_size;
		m_wordArray.resize( (m_size + 31) >> 5 );
	}

	/// Get bit.
	bool bitAt(uint32_t b) const
	{
		nvDebugCheck( b < m_size );
		return (m_wordArray[b >> 5] & (1 << (b & 31))) != 0;
	}

	// Set a bit.
	void setBitAt(uint32_t idx)
	{
		nvDebugCheck(idx < m_size);
		m_wordArray[idx >> 5] |=  (1 << (idx & 31));
	}

	// Toggle a bit.
	void toggleBitAt(uint32_t idx)
	{
		nvDebugCheck(idx < m_size);
		m_wordArray[idx >> 5] ^= (1 << (idx & 31));
	}

	// Set a bit to the given value. @@ Rename modifyBitAt?
	void setBitAt(uint32_t idx, bool b)
	{
		nvDebugCheck(idx < m_size);
		m_wordArray[idx >> 5] = setBits(m_wordArray[idx >> 5], 1 << (idx & 31), b);
		nvDebugCheck(bitAt(idx) == b);
	}

	// Clear all the bits.
	void clearAll()
	{
		memset(m_wordArray.data(), 0, m_wordArray.size() * sizeof(uint32_t ));
	}

	// Set all the bits.
	void setAll()
	{
		memset(m_wordArray.data(), 0xFF, m_wordArray.size() * sizeof(uint32_t ));
	}

private:
	// See "Conditionally set or clear bits without branching" at http://graphics.stanford.edu/~seander/bithacks.html
	inline uint32_t setBits(uint32_t w, uint32_t m, bool b)
	{
		return (w & ~m) | (-int(b) & m);
	}

	// Number of bits stored.
	uint32_t m_size;

	// Array of bits.
	std::vector<uint32_t> m_wordArray;
};

inline uint32_t sdbmHash(const void *data_in, uint32_t size, uint32_t h = 5381)
{
	const uint8_t *data = (const uint8_t *) data_in;
	uint32_t i = 0;
	while (i < size) {
		h = (h << 16) + (h << 6) - h + (uint32_t ) data[i++];
	}
	return h;
}

// Note that this hash does not handle NaN properly.
inline uint32_t sdbmFloatHash(const float *f, uint32_t count, uint32_t h = 5381)
{
	for (uint32_t i = 0; i < count; i++) {
		union {
			float f;
			uint32_t i;
		} x = { f[i] };
		if (x.i == 0x80000000) x.i = 0;
		h = sdbmHash(&x, 4, h);
	}
	return h;
}

template <typename T>
inline uint32_t hash(const T &t, uint32_t h = 5381)
{
	return sdbmHash(&t, sizeof(T), h);
}

template <>
inline uint32_t hash(const float &f, uint32_t h)
{
	return sdbmFloatHash(&f, 1, h);
}

// Functors for hash table:
template <typename Key> struct Hash
{
	uint32_t operator()(const Key &k) const
	{
		return hash(k);
	}
};

template <typename Key> struct Equal
{
	bool operator()(const Key &k0, const Key &k1) const
	{
		return k0 == k1;
	}
};

// Based on Pierre Terdiman's and Michael Herf's source code.
// http://www.codercorner.com/RadixSortRevisited.htm
// http://www.stereopsis.com/radix.html
class RadixSort
{
public:
	// Constructor/Destructor
	RadixSort();
	~RadixSort();

	RadixSort &sort(const float *input, uint32_t count);
	inline RadixSort &sort(const std::vector<float> &input)
	{
		return sort(input.data(), input.size());
	}

	// Access to results. m_ranks is a list of indices in sorted order, i.e. in the order you may further process your data
	inline const uint32_t *ranks() const
	{
		nvDebugCheck(m_validRanks);
		return m_ranks;
	}
	inline uint32_t *ranks()
	{
		nvDebugCheck(m_validRanks);
		return m_ranks;
	}

private:
	uint32_t m_size;
	uint32_t *m_ranks;
	uint32_t *m_ranks2;
	bool m_validRanks;

	// Internal methods
	template <typename T> void insertionSort(const T *input, uint32_t count);
	template <typename T> void radixSort(const T *input, uint32_t count);
};

} // namespace nv
