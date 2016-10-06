// This code is in the public domain -- castanyo@yahoo.es

#pragma once
#ifndef ATLAS_H
#define ATLAS_H

#include <algorithm>
#include <unordered_map>
#include <vector>
#include <assert.h>
#include <stdarg.h> // va_list
#include <stdint.h>
#include <stdio.h>
#include <cmath>
#include <float.h>
#include <math.h>
#include <time.h>

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

#ifndef PI
#define PI                  float(3.1415926535897932384626433833)
#endif

#define NV_EPSILON          (0.0001f)
#define NV_NORMAL_EPSILON   (0.001f)

namespace nv
{
inline int align(int x, int a)
{
	return (x + a - 1) & ~(a - 1);
}

inline bool isAligned(int x, int a)
{
	return (x & (a - 1)) == 0;
}

/// Return the maximum of the three arguments.
template <typename T>
//inline const T & max3(const T & a, const T & b, const T & c)
inline T max3(const T &a, const T &b, const T &c)
{
	return std::max(a, std::max(b, c));
}

/// Return the maximum of the three arguments.
template <typename T>
//inline const T & min3(const T & a, const T & b, const T & c)
inline T min3(const T &a, const T &b, const T &c)
{
	return std::min(a, std::min(b, c));
}

/// Clamp between two values.
template <typename T>
//inline const T & clamp(const T & x, const T & a, const T & b)
inline T clamp(const T &x, const T &a, const T &b)
{
	return std::min(std::max(x, a), b);
}

inline float saturate(float f)
{
	return clamp(f, 0.0f, 1.0f);
}

// Robust floating point comparisons:
// http://realtimecollisiondetection.net/blog/?p=89
inline bool equal(const float f0, const float f1, const float epsilon = NV_EPSILON)
{
	//return fabs(f0-f1) <= epsilon;
	return fabs(f0 - f1) <= epsilon * max3(1.0f, fabsf(f0), fabsf(f1));
}

NV_FORCEINLINE int ftoi_floor(float val)
{
	return (int)val;
}

NV_FORCEINLINE int ftoi_ceil(float val)
{
	return (int)ceilf(val);
}

NV_FORCEINLINE int ftoi_round(float f)
{
	return int(floorf(f + 0.5f));
}

inline bool isZero(const float f, const float epsilon = NV_EPSILON)
{
	return fabs(f) <= epsilon;
}

inline float lerp(float f0, float f1, float t)
{
	const float s = 1.0f - t;
	return f0 * s + f1 * t;
}

inline float square(float f)
{
	return f * f;
}

inline int square(int i)
{
	return i * i;
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

class Vector2
{
public:
	typedef Vector2 const &Arg;

	Vector2();
	explicit Vector2(float f);
	Vector2(float x, float y);
	Vector2(Vector2::Arg v);

	//template <typename T> explicit Vector2(const T & v) : x(v.x), y(v.y) {}
	//template <typename T> operator T() const { return T(x, y); }

	const Vector2 &operator=(Vector2::Arg v);

	const float *ptr() const;

	void set(float x, float y);

	Vector2 operator-() const;
	void operator+=(Vector2::Arg v);
	void operator-=(Vector2::Arg v);
	void operator*=(float s);
	void operator*=(Vector2::Arg v);

	friend bool operator==(Vector2::Arg a, Vector2::Arg b);
	friend bool operator!=(Vector2::Arg a, Vector2::Arg b);

	union {
		struct {
			float x, y;
		};
		float component[2];
	};
};

class Vector3
{
public:
	typedef Vector3 const &Arg;

	Vector3();
	explicit Vector3(float x);
	//explicit Vector3(int x) : x(float(x)), y(float(x)), z(float(x)) {}
	Vector3(float x, float y, float z);
	Vector3(Vector2::Arg v, float z);
	Vector3(Vector3::Arg v);

	//template <typename T> explicit Vector3(const T & v) : x(v.x), y(v.y), z(v.z) {}
	//template <typename T> operator T() const { return T(x, y, z); }

	const Vector3 &operator=(Vector3::Arg v);

	Vector2 xy() const;

	const float *ptr() const;

	void set(float x, float y, float z);

	Vector3 operator-() const;
	void operator+=(Vector3::Arg v);
	void operator-=(Vector3::Arg v);
	void operator*=(float s);
	void operator/=(float s);
	void operator*=(Vector3::Arg v);
	void operator/=(Vector3::Arg v);

	friend bool operator==(Vector3::Arg a, Vector3::Arg b);
	friend bool operator!=(Vector3::Arg a, Vector3::Arg b);

	union {
		struct {
			float x, y, z;
		};
		float component[3];
	};
};

// Vector2
inline Vector2::Vector2() {}
inline Vector2::Vector2(float f) : x(f), y(f) {}
inline Vector2::Vector2(float x, float y) : x(x), y(y) {}
inline Vector2::Vector2(Vector2::Arg v) : x(v.x), y(v.y) {}

inline const Vector2 &Vector2::operator=(Vector2::Arg v)
{
	x = v.x;
	y = v.y;
	return *this;
}

inline const float *Vector2::ptr() const
{
	return &x;
}

inline void Vector2::set(float x, float y)
{
	this->x = x;
	this->y = y;
}

inline Vector2 Vector2::operator-() const
{
	return Vector2(-x, -y);
}

inline void Vector2::operator+=(Vector2::Arg v)
{
	x += v.x;
	y += v.y;
}

inline void Vector2::operator-=(Vector2::Arg v)
{
	x -= v.x;
	y -= v.y;
}

inline void Vector2::operator*=(float s)
{
	x *= s;
	y *= s;
}

inline void Vector2::operator*=(Vector2::Arg v)
{
	x *= v.x;
	y *= v.y;
}

inline bool operator==(Vector2::Arg a, Vector2::Arg b)
{
	return a.x == b.x && a.y == b.y;
}
inline bool operator!=(Vector2::Arg a, Vector2::Arg b)
{
	return a.x != b.x || a.y != b.y;
}


// Vector3
inline Vector3::Vector3() {}
inline Vector3::Vector3(float f) : x(f), y(f), z(f) {}
inline Vector3::Vector3(float x, float y, float z) : x(x), y(y), z(z) {}
inline Vector3::Vector3(Vector2::Arg v, float z) : x(v.x), y(v.y), z(z) {}
inline Vector3::Vector3(Vector3::Arg v) : x(v.x), y(v.y), z(v.z) {}

inline const Vector3 &Vector3::operator=(Vector3::Arg v)
{
	x = v.x;
	y = v.y;
	z = v.z;
	return *this;
}


inline Vector2 Vector3::xy() const
{
	return Vector2(x, y);
}

inline const float *Vector3::ptr() const
{
	return &x;
}

inline void Vector3::set(float x, float y, float z)
{
	this->x = x;
	this->y = y;
	this->z = z;
}

inline Vector3 Vector3::operator-() const
{
	return Vector3(-x, -y, -z);
}

inline void Vector3::operator+=(Vector3::Arg v)
{
	x += v.x;
	y += v.y;
	z += v.z;
}

inline void Vector3::operator-=(Vector3::Arg v)
{
	x -= v.x;
	y -= v.y;
	z -= v.z;
}

inline void Vector3::operator*=(float s)
{
	x *= s;
	y *= s;
	z *= s;
}

inline void Vector3::operator/=(float s)
{
	float is = 1.0f / s;
	x *= is;
	y *= is;
	z *= is;
}

inline void Vector3::operator*=(Vector3::Arg v)
{
	x *= v.x;
	y *= v.y;
	z *= v.z;
}

inline void Vector3::operator/=(Vector3::Arg v)
{
	x /= v.x;
	y /= v.y;
	z /= v.z;
}

inline bool operator==(Vector3::Arg a, Vector3::Arg b)
{
	return a.x == b.x && a.y == b.y && a.z == b.z;
}
inline bool operator!=(Vector3::Arg a, Vector3::Arg b)
{
	return a.x != b.x || a.y != b.y || a.z != b.z;
}

// Functions


// Vector2

inline Vector2 operator+(Vector2::Arg a, Vector2::Arg b)
{
	return Vector2(a.x + b.x, a.y + b.y);
}

inline Vector2 operator-(Vector2::Arg a, Vector2::Arg b)
{
	return Vector2(a.x - b.x, a.y - b.y);
}

inline Vector2 operator*(Vector2::Arg v, float s)
{
	return Vector2(v.x * s, v.y * s);
}

inline Vector2 operator*(Vector2::Arg v1, Vector2::Arg v2)
{
	return Vector2(v1.x * v2.x, v1.y * v2.y);
}

inline Vector2 operator/(Vector2::Arg v, float s)
{
	return Vector2(v.x / s, v.y / s);
}

inline Vector2 lerp(Vector2::Arg v1, Vector2::Arg v2, float t)
{
	const float s = 1.0f - t;
	return Vector2(v1.x * s + t * v2.x, v1.y * s + t * v2.y);
}

inline float dot(Vector2::Arg a, Vector2::Arg b)
{
	return a.x * b.x + a.y * b.y;
}

inline float lengthSquared(Vector2::Arg v)
{
	return v.x * v.x + v.y * v.y;
}

inline float length(Vector2::Arg v)
{
	return sqrtf(lengthSquared(v));
}

inline float distance(Vector2::Arg a, Vector2::Arg b)
{
	return length(a - b);
}

inline bool isNormalized(Vector2::Arg v, float epsilon = NV_NORMAL_EPSILON)
{
	return equal(length(v), 1, epsilon);
}

inline Vector2 normalize(Vector2::Arg v, float epsilon = NV_EPSILON)
{
	float l = length(v);
	nvDebugCheck(!isZero(l, epsilon));
	Vector2 n = v * (1.0f / l);
	nvDebugCheck(isNormalized(n));
	return n;
}

inline Vector2 normalizeSafe(Vector2::Arg v, Vector2::Arg fallback, float epsilon = NV_EPSILON)
{
	float l = length(v);
	if (isZero(l, epsilon)) {
		return fallback;
	}
	return v * (1.0f / l);
}

inline bool equal(Vector2::Arg v1, Vector2::Arg v2, float epsilon = NV_EPSILON)
{
	return equal(v1.x, v2.x, epsilon) && equal(v1.y, v2.y, epsilon);
}

inline Vector2 max(Vector2::Arg a, Vector2::Arg b)
{
	return Vector2(std::max(a.x, b.x), std::max(a.y, b.y));
}

inline bool isFinite(Vector2::Arg v)
{
	return std::isfinite(v.x) && std::isfinite(v.y);
}

// Note, this is the area scaled by 2!
inline float triangleArea(Vector2::Arg v0, Vector2::Arg v1)
{
	return (v0.x * v1.y - v0.y * v1.x); // * 0.5f;
}
inline float triangleArea(Vector2::Arg a, Vector2::Arg b, Vector2::Arg c)
{
	// IC: While it may be appealing to use the following expression:
	//return (c.x * a.y + a.x * b.y + b.x * c.y - b.x * a.y - c.x * b.y - a.x * c.y); // * 0.5f;
	// That's actually a terrible idea. Small triangles far from the origin can end up producing fairly large floating point
	// numbers and the results becomes very unstable and dependent on the order of the factors.
	// Instead, it's preferable to subtract the vertices first, and multiply the resulting small values together. The result
	// in this case is always much more accurate (as long as the triangle is small) and less dependent of the location of
	// the triangle.
	//return ((a.x - c.x) * (b.y - c.y) - (a.y - c.y) * (b.x - c.x)); // * 0.5f;
	return triangleArea(a - c, b - c);
}

inline float triangleArea2(Vector2::Arg v1, Vector2::Arg v2, Vector2::Arg v3)
{
	return 0.5f * (v3.x * v1.y + v1.x * v2.y + v2.x * v3.y - v2.x * v1.y - v3.x * v2.y - v1.x * v3.y);
}

template <>
inline uint32_t hash(const Vector2 &v, uint32_t h)
{
	return sdbmFloatHash(v.component, 2, h);
}



// Vector3

inline Vector3 add(Vector3::Arg a, Vector3::Arg b)
{
	return Vector3(a.x + b.x, a.y + b.y, a.z + b.z);
}
inline Vector3 add(Vector3::Arg a, float b)
{
	return Vector3(a.x + b, a.y + b, a.z + b);
}
inline Vector3 operator+(Vector3::Arg a, Vector3::Arg b)
{
	return add(a, b);
}
inline Vector3 operator+(Vector3::Arg a, float b)
{
	return add(a, b);
}

inline Vector3 sub(Vector3::Arg a, Vector3::Arg b)
{
	return Vector3(a.x - b.x, a.y - b.y, a.z - b.z);
}
inline Vector3 sub(Vector3::Arg a, float b)
{
	return Vector3(a.x - b, a.y - b, a.z - b);
}
inline Vector3 operator-(Vector3::Arg a, Vector3::Arg b)
{
	return sub(a, b);
}
inline Vector3 operator-(Vector3::Arg a, float b)
{
	return sub(a, b);
}

inline Vector3 cross(Vector3::Arg a, Vector3::Arg b)
{
	return Vector3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

inline Vector3 operator*(Vector3::Arg v, float s)
{
	return Vector3(v.x * s, v.y * s, v.z * s);
}

inline Vector3 operator*(float s, Vector3::Arg v)
{
	return Vector3(v.x * s, v.y * s, v.z * s);
}

inline Vector3 operator*(Vector3::Arg v, Vector3::Arg s)
{
	return Vector3(v.x * s.x, v.y * s.y, v.z * s.z);
}

inline Vector3 operator/(Vector3::Arg v, float s)
{
	return v * (1.0f / s);
}

inline Vector3 lerp(Vector3::Arg v1, Vector3::Arg v2, float t)
{
	const float s = 1.0f - t;
	return Vector3(v1.x * s + t * v2.x, v1.y * s + t * v2.y, v1.z * s + t * v2.z);
}

inline float dot(Vector3::Arg a, Vector3::Arg b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline float lengthSquared(Vector3::Arg v)
{
	return v.x * v.x + v.y * v.y + v.z * v.z;
}

inline float length(Vector3::Arg v)
{
	return sqrtf(lengthSquared(v));
}

inline float distance(Vector3::Arg a, Vector3::Arg b)
{
	return length(a - b);
}

inline float distanceSquared(Vector3::Arg a, Vector3::Arg b)
{
	return lengthSquared(a - b);
}

inline bool isNormalized(Vector3::Arg v, float epsilon = NV_NORMAL_EPSILON)
{
	return equal(length(v), 1, epsilon);
}

inline Vector3 normalize(Vector3::Arg v, float epsilon = NV_EPSILON)
{
	float l = length(v);
	nvDebugCheck(!isZero(l, epsilon));
	Vector3 n = v * (1.0f / l);
	nvDebugCheck(isNormalized(n));
	return n;
}

inline Vector3 normalizeSafe(Vector3::Arg v, Vector3::Arg fallback, float epsilon = NV_EPSILON)
{
	float l = length(v);
	if (isZero(l, epsilon)) {
		return fallback;
	}
	return v * (1.0f / l);
}

inline bool equal(Vector3::Arg v1, Vector3::Arg v2, float epsilon = NV_EPSILON)
{
	return equal(v1.x, v2.x, epsilon) && equal(v1.y, v2.y, epsilon) && equal(v1.z, v2.z, epsilon);
}

inline Vector3 min(Vector3::Arg a, Vector3::Arg b)
{
	return Vector3(std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z));
}

inline Vector3 max(Vector3::Arg a, Vector3::Arg b)
{
	return Vector3(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z));
}

inline Vector3 clamp(Vector3::Arg v, float min, float max)
{
	return Vector3(clamp(v.x, min, max), clamp(v.y, min, max), clamp(v.z, min, max));
}

inline Vector3 saturate(Vector3::Arg v)
{
	return Vector3(saturate(v.x), saturate(v.y), saturate(v.z));
}

inline Vector3 floor(Vector3::Arg v)
{
	return Vector3(floorf(v.x), floorf(v.y), floorf(v.z));
}

inline Vector3 ceil(Vector3::Arg v)
{
	return Vector3(ceilf(v.x), ceilf(v.y), ceilf(v.z));
}

inline bool isFinite(Vector3::Arg v)
{
	return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

template <>
inline uint32_t hash(const Vector3 &v, uint32_t h)
{
	return sdbmFloatHash(v.component, 3, h);
}

/// Basis class to compute tangent space basis, ortogonalizations and to
/// transform vectors from one space to another.
class Basis
{
public:

	/// Create a null basis.
	Basis() : tangent(0, 0, 0), bitangent(0, 0, 0), normal(0, 0, 0) {}

	void buildFrameForDirection(Vector3::Arg d, float angle = 0)
	{
		nvCheck(isNormalized(d));
		normal = d;
		// Choose minimum axis.
		if (fabsf(normal.x) < fabsf(normal.y) && fabsf(normal.x) < fabsf(normal.z)) {
			tangent = Vector3(1, 0, 0);
		} else if (fabsf(normal.y) < fabsf(normal.z)) {
			tangent = Vector3(0, 1, 0);
		} else {
			tangent = Vector3(0, 0, 1);
		}
		// Ortogonalize
		tangent -= normal * dot(normal, tangent);
		tangent = normalize(tangent);
		bitangent = cross(normal, tangent);
		// Rotate frame around normal according to angle.
		if (angle != 0.0f) {
			float c = cosf(angle);
			float s = sinf(angle);
			Vector3 tmp = c * tangent - s * bitangent;
			bitangent = s * tangent + c * bitangent;
			tangent = tmp;
		}
	}


	Vector3 tangent;
	Vector3 bitangent;
	Vector3 normal;
};

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

/// Bit map. This should probably be called BitImage.
class BitMap
{
public:
	BitMap() : m_width(0), m_height(0) {}
	BitMap(uint32_t w, uint32_t h) : m_width(w), m_height(h), m_bitArray(w * h) {}

	uint32_t width() const
	{
		return m_width;
	}
	uint32_t height() const
	{
		return m_height;
	}

	void resize(uint32_t w, uint32_t h, bool initValue)
	{
		BitArray tmp(w * h);
		if (initValue) tmp.setAll();
		else tmp.clearAll();
		// @@ Copying one bit at a time. This could be much faster.
		for (uint32_t y = 0; y < m_height; y++) {
			for (uint32_t x = 0; x < m_width; x++) {
				//tmp.setBitAt(y*w + x, bitAt(x, y));
				if (bitAt(x, y) != initValue) tmp.toggleBitAt(y * w + x);
			}
		}
		std::swap(m_bitArray, tmp);
		m_width = w;
		m_height = h;
	}


	bool bitAt(uint32_t x, uint32_t y) const
	{
		nvDebugCheck(x < m_width && y < m_height);
		return m_bitArray.bitAt(y * m_width + x);
	}

	void setBitAt(uint32_t x, uint32_t y)
	{
		nvDebugCheck(x < m_width && y < m_height);
		m_bitArray.setBitAt(y * m_width + x);
	}

	void clearAll()
	{
		m_bitArray.clearAll();
	}

private:
	uint32_t m_width;
	uint32_t m_height;
	BitArray m_bitArray;
};

// Axis Aligned Bounding Box.
class Box
{
public:

	inline Box() {}
	inline Box(const Box &b) : minCorner(b.minCorner), maxCorner(b.maxCorner) {}
	inline Box(const Vector3 &mins, const Vector3 &maxs) : minCorner(mins), maxCorner(maxs) {}

	operator const float *() const
	{
		return reinterpret_cast<const float *>(this);
	}

	// Clear the bounds.
	void clearBounds()
	{
		minCorner.set(FLT_MAX, FLT_MAX, FLT_MAX);
		maxCorner.set(-FLT_MAX, -FLT_MAX, -FLT_MAX);
	}

	// Return extents of the box.
	Vector3 Box::extents() const
	{
		return (maxCorner - minCorner) * 0.5f;
	}

	// Add a point to this box.
	void addPointToBounds(const Vector3 &p)
	{
		minCorner = min(minCorner, p);
		maxCorner = max(maxCorner, p);
	}

	// Get the volume of the box.
	float volume() const
	{
		Vector3 d = extents();
		return 8.0f * (d.x * d.y * d.z);
	}

	const Vector3 &corner(int i) const
	{
		return (&minCorner)[i];
	}

	Vector3 minCorner;
	Vector3 maxCorner;
};

namespace fit {
Vector3 computeCentroid(int n, const Vector3 *points);

Vector3 computeCovariance(int n, const Vector3 *points, float *covariance);

bool isPlanar(int n, const Vector3 *points, float epsilon = NV_EPSILON);

bool eigenSolveSymmetric3(const float matrix[6], float eigenValues[3], Vector3 eigenVectors[3]);
} // namespace fit

/// Fixed size vector class.
class FullVector
{
public:
	FullVector(uint32_t dim) { m_array.resize(dim); }
	FullVector(const FullVector &v) : m_array(v.m_array) {}

	const FullVector &operator=(const FullVector &v)
	{
		nvCheck(dimension() == v.dimension());
		m_array = v.m_array;
		return *this;
	}

	uint32_t dimension() const { return m_array.size(); }
	const float &operator[]( uint32_t index ) const { return m_array[index]; }
	float &operator[] ( uint32_t index ) { return m_array[index]; }

	void fill(float f)
	{
		const uint32_t dim = dimension();
		for (uint32_t i = 0; i < dim; i++) {
			m_array[i] = f;
		}
	}

	void operator+=(const FullVector &v)
	{
		nvDebugCheck(dimension() == v.dimension());
		const uint32_t dim = dimension();
		for (uint32_t i = 0; i < dim; i++) {
			m_array[i] += v.m_array[i];
		}
	}

	void operator-=(const FullVector &v)
	{
		nvDebugCheck(dimension() == v.dimension());
		const uint32_t dim = dimension();
		for (uint32_t i = 0; i < dim; i++) {
			m_array[i] -= v.m_array[i];
		}
	}

	void operator*=(const FullVector &v)
	{
		nvDebugCheck(dimension() == v.dimension());
		const uint32_t dim = dimension();
		for (uint32_t i = 0; i < dim; i++) {
			m_array[i] *= v.m_array[i];
		}
	}

	void operator+=(float f)
	{
		const uint32_t dim = dimension();
		for (uint32_t i = 0; i < dim; i++) {
			m_array[i] += f;
		}
	}

	void operator-=(float f)
	{
		const uint32_t dim = dimension();
		for (uint32_t i = 0; i < dim; i++) {
			m_array[i] -= f;
		}
	}

	void operator*=(float f)
	{
		const uint32_t dim = dimension();
		for (uint32_t i = 0; i < dim; i++) {
			m_array[i] *= f;
		}
	}

private:
	std::vector<float> m_array;
};

namespace HalfEdge {
class Face;
class Vertex;

class Edge
{
public:
	uint32_t id;
	Edge *next;
	Edge *prev;	// This is not strictly half-edge, but makes algorithms easier and faster.
	Edge *pair;
	Vertex *vertex;
	Face *face;

	// Default constructor.
	Edge(uint32_t id) : id(id), next(NULL), prev(NULL), pair(NULL), vertex(NULL), face(NULL)
	{
	}

	// Vertex queries.
	const Vertex *from() const
	{
		return vertex;
	}

	Vertex *from()
	{
		return vertex;
	}

	const Vertex *to() const
	{
		return pair->vertex;    // This used to be 'next->vertex', but that changed often when the connectivity of the mesh changes.
	}

	Vertex *to()
	{
		return pair->vertex;
	}

	// Edge queries.
	void setNext(Edge *e)
	{
		next = e;
		if (e != NULL) e->prev = this;
	}
	void setPrev(Edge *e)
	{
		prev = e;
		if (e != NULL) e->next = this;
	}

	// @@ It would be more simple to only check m_pair == NULL
	// Face queries.
	bool isBoundary() const
	{
		return !(face && pair->face);
	}

	// @@ This is not exactly accurate, we should compare the texture coordinates...
	bool isSeam() const
	{
		return vertex != pair->next->vertex || next->vertex != pair->vertex;
	}

	bool isValid() const
	{
		// null face is OK.
		if (next == NULL || prev == NULL || pair == NULL || vertex == NULL) return false;
		if (next->prev != this) return false;
		if (prev->next != this) return false;
		if (pair->pair != this) return false;
		return true;
	}

	// Geometric queries.
	Vector3 midPoint() const;

	float length() const;

	// Return angle between this edge and the previous one.
	float angle() const;
};

class Vertex
{
public:
	uint32_t id;
	Edge *edge;
	Vertex *next;
	Vertex *prev;
	Vector3 pos;
	Vector3 nor;
	Vector2 tex;

	Vertex(uint32_t id) : id(id), edge(NULL), pos(0.0f), nor(0.0f), tex(0.0f)
	{
		next = this;
		prev = this;
	}

	// Set first edge of all colocals.
	void setEdge(Edge *e)
	{
		for (VertexIterator it(colocals()); !it.isDone(); it.advance()) {
			it.current()->edge = e;
		}
	}

	// Update position of all colocals.
	void setPos(const Vector3 &p)
	{
		for (VertexIterator it(colocals()); !it.isDone(); it.advance()) {
			it.current()->pos = p;
		}
	}

	uint32_t colocalCount() const
	{
		uint32_t count = 0;
		for (ConstVertexIterator it(colocals()); !it.isDone(); it.advance()) {
			++count;
		}
		return count;
	}

	uint32_t valence() const
	{
		uint32_t count = 0;
		for (ConstEdgeIterator it(edges()); !it.isDone(); it.advance()) {
			++count;
		}
		return count;
	}

	bool isFirstColocal() const
	{
		return firstColocal() == this;
	}

	const Vertex *firstColocal() const
	{
		uint32_t firstId = id;
		const Vertex *vertex = this;
		for (ConstVertexIterator it(colocals()); !it.isDone(); it.advance()) {
			if (it.current()->id < firstId) {
				firstId = vertex->id;
				vertex = it.current();
			}
		}
		return vertex;
	}

	Vertex *firstColocal()
	{
		Vertex *vertex = this;
		uint32_t firstId = id;
		for (VertexIterator it(colocals()); !it.isDone(); it.advance()) {
			if (it.current()->id < firstId) {
				firstId = vertex->id;
				vertex = it.current();
			}
		}
		return vertex;
	}

	bool isColocal(const Vertex *v) const
	{
		if (this == v) return true;
		if (pos != v->pos) return false;
		for (ConstVertexIterator it(colocals()); !it.isDone(); it.advance()) {
			if (v == it.current()) {
				return true;
			}
		}
		return false;
	}

	void linkColocal(Vertex *v)
	{
		next->prev = v;
		v->next = next;
		next = v;
		v->prev = this;
	}
	void unlinkColocal()
	{
		next->prev = prev;
		prev->next = next;
		next = this;
		prev = this;
	}

	// @@ Note: This only works if linkBoundary has been called.
	bool isBoundary() const
	{
		return (edge && !edge->face);
	}

	// Iterator that visits the edges around this vertex in counterclockwise order.
	class EdgeIterator //: public Iterator<Edge *>
	{
	public:
		EdgeIterator(Edge *e) : m_end(NULL), m_current(e) { }

		virtual void advance()
		{
			if (m_end == NULL) m_end = m_current;
			m_current = m_current->pair->next;
			//m_current = m_current->prev->pair;
		}

		virtual bool isDone() const
		{
			return m_end == m_current;
		}
		virtual Edge *current() const
		{
			return m_current;
		}
		Vertex *vertex() const
		{
			return m_current->vertex;
		}

	private:
		Edge *m_end;
		Edge *m_current;
	};

	EdgeIterator edges()
	{
		return EdgeIterator(edge);
	}
	EdgeIterator edges(Edge *e)
	{
		return EdgeIterator(e);
	}

	// Iterator that visits the edges around this vertex in counterclockwise order.
	class ConstEdgeIterator //: public Iterator<Edge *>
	{
	public:
		ConstEdgeIterator(const Edge *e) : m_end(NULL), m_current(e) { }
		ConstEdgeIterator(EdgeIterator it) : m_end(NULL), m_current(it.current()) { }

		virtual void advance()
		{
			if (m_end == NULL) m_end = m_current;
			m_current = m_current->pair->next;
			//m_current = m_current->prev->pair;
		}

		virtual bool isDone() const
		{
			return m_end == m_current;
		}
		virtual const Edge *current() const
		{
			return m_current;
		}
		const Vertex *vertex() const
		{
			return m_current->to();
		}

	private:
		const Edge *m_end;
		const Edge *m_current;
	};

	ConstEdgeIterator edges() const
	{
		return ConstEdgeIterator(edge);
	}
	ConstEdgeIterator edges(const Edge *e) const
	{
		return ConstEdgeIterator(e);
	}

	// Iterator that visits all the colocal vertices.
	class VertexIterator //: public Iterator<Edge *>
	{
	public:
		VertexIterator(Vertex *v) : m_end(NULL), m_current(v) { }

		virtual void advance()
		{
			if (m_end == NULL) m_end = m_current;
			m_current = m_current->next;
		}

		virtual bool isDone() const
		{
			return m_end == m_current;
		}
		virtual Vertex *current() const
		{
			return m_current;
		}

	private:
		Vertex *m_end;
		Vertex *m_current;
	};

	VertexIterator colocals()
	{
		return VertexIterator(this);
	}

	// Iterator that visits all the colocal vertices.
	class ConstVertexIterator //: public Iterator<Edge *>
	{
	public:
		ConstVertexIterator(const Vertex *v) : m_end(NULL), m_current(v) { }

		virtual void advance()
		{
			if (m_end == NULL) m_end = m_current;
			m_current = m_current->next;
		}

		virtual bool isDone() const
		{
			return m_end == m_current;
		}
		virtual const Vertex *current() const
		{
			return m_current;
		}

	private:
		const Vertex *m_end;
		const Vertex *m_current;
	};

	ConstVertexIterator colocals() const
	{
		return ConstVertexIterator(this);
	}

};

class Face
{
public:
	uint32_t id;
	uint16_t group;
	uint16_t material;
	Edge *edge;

	Face(uint32_t id) : id(id), group(~0), material(~0), edge(NULL) {}

	float area() const
	{
		float area = 0;
		const Vector3 &v0 = edge->from()->pos;
		for (ConstEdgeIterator it(edges(edge->next)); it.current() != edge->prev; it.advance()) {
			const Edge *e = it.current();
			const Vector3 &v1 = e->vertex->pos;
			const Vector3 &v2 = e->next->vertex->pos;
			area += length(cross(v1 - v0, v2 - v0));
		}
		return area * 0.5f;
	}

	float parametricArea() const
	{
		float area = 0;
		const Vector2 &v0 = edge->from()->tex;
		for (ConstEdgeIterator it(edges(edge->next)); it.current() != edge->prev; it.advance()) {
			const Edge *e = it.current();
			const Vector2 &v1 = e->vertex->tex;
			const Vector2 &v2 = e->next->vertex->tex;
			area += triangleArea(v0, v1, v2);
		}
		return area * 0.5f;
	}

	Vector3 normal() const
	{
		Vector3 n(0);
		const Vertex *vertex0 = NULL;
		for (ConstEdgeIterator it(edges()); !it.isDone(); it.advance()) {
			const Edge *edge = it.current();
			nvCheck(edge != NULL);
			if (vertex0 == NULL) {
				vertex0 = edge->vertex;
			} else if (edge->next->vertex != vertex0) {
				const HalfEdge::Vertex *vertex1 = edge->from();
				const HalfEdge::Vertex *vertex2 = edge->to();
				const Vector3 &p0 = vertex0->pos;
				const Vector3 &p1 = vertex1->pos;
				const Vector3 &p2 = vertex2->pos;
				Vector3 v10 = p1 - p0;
				Vector3 v20 = p2 - p0;
				n += cross(v10, v20);
			}
		}
		return normalizeSafe(n, Vector3(0, 0, 1), 0.0f);
	}

	Vector3 centroid() const
	{
		Vector3 sum(0.0f);
		uint32_t count = 0;
		for (ConstEdgeIterator it(edges()); !it.isDone(); it.advance()) {
			const Edge *edge = it.current();
			sum += edge->from()->pos;
			count++;
		}
		return sum / float(count);
	}

	bool isValid() const
	{
		uint32_t count = 0;
		for (ConstEdgeIterator it(edges()); !it.isDone(); it.advance()) {
			const Edge *edge = it.current();
			if (edge->face != this) return false;
			if (!edge->isValid()) return false;
			if (!edge->pair->isValid()) return false;
			count++;
		}
		if (count < 3) return false;
		return true;
	}

	bool contains(const Edge *e) const
	{
		for (ConstEdgeIterator it(edges()); !it.isDone(); it.advance()) {
			if (it.current() == e) return true;
		}
		return false;
	}

	uint32_t edgeCount() const
	{
		uint32_t count = 0;
		for (ConstEdgeIterator it(edges()); !it.isDone(); it.advance()) {
			++count;
		}
		return count;
	}

	// The iterator that visits the edges of this face in clockwise order.
	class EdgeIterator //: public Iterator<Edge *>
	{
	public:
		EdgeIterator(Edge *e) : m_end(NULL), m_current(e) { }

		virtual void advance()
		{
			if (m_end == NULL) m_end = m_current;
			m_current = m_current->next;
		}

		virtual bool isDone() const
		{
			return m_end == m_current;
		}
		virtual Edge *current() const
		{
			return m_current;
		}
		Vertex *vertex() const
		{
			return m_current->vertex;
		}

	private:
		Edge *m_end;
		Edge *m_current;
	};

	EdgeIterator edges()
	{
		return EdgeIterator(edge);
	}
	EdgeIterator edges(Edge *e)
	{
		nvDebugCheck(contains(e));
		return EdgeIterator(e);
	}

	// The iterator that visits the edges of this face in clockwise order.
	class ConstEdgeIterator //: public Iterator<const Edge *>
	{
	public:
		ConstEdgeIterator(const Edge *e) : m_end(NULL), m_current(e) { }
		ConstEdgeIterator(const EdgeIterator &it) : m_end(NULL), m_current(it.current()) { }

		virtual void advance()
		{
			if (m_end == NULL) m_end = m_current;
			m_current = m_current->next;
		}

		virtual bool isDone() const
		{
			return m_end == m_current;
		}
		virtual const Edge *current() const
		{
			return m_current;
		}
		const Vertex *vertex() const
		{
			return m_current->vertex;
		}

	private:
		const Edge *m_end;
		const Edge *m_current;
	};

	ConstEdgeIterator edges() const
	{
		return ConstEdgeIterator(edge);
	}
	ConstEdgeIterator edges(const Edge *e) const
	{
		nvDebugCheck(contains(e));
		return ConstEdgeIterator(e);
	}
};

/// Simple half edge mesh designed for dynamic mesh manipulation.
class Mesh
{
public:

	Mesh();
	Mesh(const Mesh *mesh);
	~Mesh();

	void clear();

	Vertex *addVertex(const Vector3 &pos);
	//Vertex * addVertex(uint32_t id, const Vector3 & pos);
	//void addVertices(const Mesh * mesh);

	void linkColocals();
	void linkColocalsWithCanonicalMap(const std::vector<uint32_t> &canonicalMap);

	Face *addFace();
	Face *addFace(uint32_t v0, uint32_t v1, uint32_t v2);
	Face *addFace(uint32_t v0, uint32_t v1, uint32_t v2, uint32_t v3);
	Face *addFace(const std::vector<uint32_t> &indexArray);
	Face *addFace(const std::vector<uint32_t> &indexArray, uint32_t first, uint32_t num);
	//void addFaces(const Mesh * mesh);

	// These functions disconnect the given element from the mesh and delete it.
	void disconnect(Edge *edge);

	void remove(Edge *edge);
	void remove(Vertex *vertex);
	void remove(Face *face);

	// Remove holes from arrays and reassign indices.
	void compactEdges();
	void compactVertices();
	void compactFaces();

	void triangulate();

	void linkBoundary();

	bool splitBoundaryEdges(); // Returns true if any split was made.

	// Sew the boundary that starts at the given edge, returns one edge that still belongs to boundary, or NULL if boundary closed.
	HalfEdge::Edge *sewBoundary(Edge *startEdge);


	// Vertices
	uint32_t vertexCount() const
	{
		return m_vertexArray.size();
	}
	const Vertex *vertexAt(int i) const
	{
		return m_vertexArray[i];
	}
	Vertex *vertexAt(int i)
	{
		return m_vertexArray[i];
	}

	uint32_t colocalVertexCount() const
	{
		return m_colocalVertexCount;
	}

	// Faces
	uint32_t faceCount() const
	{
		return m_faceArray.size();
	}
	const Face *faceAt(int i) const
	{
		return m_faceArray[i];
	}
	Face *faceAt(int i)
	{
		return m_faceArray[i];
	}

	// Edges
	uint32_t edgeCount() const
	{
		return m_edgeArray.size();
	}
	const Edge *edgeAt(int i) const
	{
		return m_edgeArray[i];
	}
	Edge *edgeAt(int i)
	{
		return m_edgeArray[i];
	}

	class ConstVertexIterator;

	class VertexIterator
	{
		friend class ConstVertexIterator;
	public:
		VertexIterator(Mesh *mesh) : m_mesh(mesh), m_current(0) { }

		virtual void advance()
		{
			m_current++;
		}
		virtual bool isDone() const
		{
			return m_current == m_mesh->vertexCount();
		}
		virtual Vertex *current() const
		{
			return m_mesh->vertexAt(m_current);
		}

	private:
		HalfEdge::Mesh *m_mesh;
		uint32_t m_current;
	};
	VertexIterator vertices()
	{
		return VertexIterator(this);
	}

	class ConstVertexIterator
	{
	public:
		ConstVertexIterator(const Mesh *mesh) : m_mesh(mesh), m_current(0) { }
		ConstVertexIterator(class VertexIterator &it) : m_mesh(it.m_mesh), m_current(it.m_current) { }

		virtual void advance()
		{
			m_current++;
		}
		virtual bool isDone() const
		{
			return m_current == m_mesh->vertexCount();
		}
		virtual const Vertex *current() const
		{
			return m_mesh->vertexAt(m_current);
		}

	private:
		const HalfEdge::Mesh *m_mesh;
		uint32_t m_current;
	};
	ConstVertexIterator vertices() const
	{
		return ConstVertexIterator(this);
	}

	class ConstFaceIterator;

	class FaceIterator
	{
		friend class ConstFaceIterator;
	public:
		FaceIterator(Mesh *mesh) : m_mesh(mesh), m_current(0) { }

		virtual void advance()
		{
			m_current++;
		}
		virtual bool isDone() const
		{
			return m_current == m_mesh->faceCount();
		}
		virtual Face *current() const
		{
			return m_mesh->faceAt(m_current);
		}

	private:
		HalfEdge::Mesh *m_mesh;
		uint32_t m_current;
	};
	FaceIterator faces()
	{
		return FaceIterator(this);
	}

	class ConstFaceIterator
	{
	public:
		ConstFaceIterator(const Mesh *mesh) : m_mesh(mesh), m_current(0) { }
		ConstFaceIterator(const FaceIterator &it) : m_mesh(it.m_mesh), m_current(it.m_current) { }

		virtual void advance()
		{
			m_current++;
		}
		virtual bool isDone() const
		{
			return m_current == m_mesh->faceCount();
		}
		virtual const Face *current() const
		{
			return m_mesh->faceAt(m_current);
		}

	private:
		const HalfEdge::Mesh *m_mesh;
		uint32_t m_current;
	};
	ConstFaceIterator faces() const
	{
		return ConstFaceIterator(this);
	}

	class ConstEdgeIterator;

	class EdgeIterator
	{
		friend class ConstEdgeIterator;
	public:
		EdgeIterator(Mesh *mesh) : m_mesh(mesh), m_current(0) { }

		virtual void advance()
		{
			m_current++;
		}
		virtual bool isDone() const
		{
			return m_current == m_mesh->edgeCount();
		}
		virtual Edge *current() const
		{
			return m_mesh->edgeAt(m_current);
		}

	private:
		HalfEdge::Mesh *m_mesh;
		uint32_t m_current;
	};
	EdgeIterator edges()
	{
		return EdgeIterator(this);
	}

	class ConstEdgeIterator
	{
	public:
		ConstEdgeIterator(const Mesh *mesh) : m_mesh(mesh), m_current(0) { }
		ConstEdgeIterator(const EdgeIterator &it) : m_mesh(it.m_mesh), m_current(it.m_current) { }

		virtual void advance()
		{
			m_current++;
		}
		virtual bool isDone() const
		{
			return m_current == m_mesh->edgeCount();
		}
		virtual const Edge *current() const
		{
			return m_mesh->edgeAt(m_current);
		}

	private:
		const HalfEdge::Mesh *m_mesh;
		uint32_t m_current;
	};
	ConstEdgeIterator edges() const
	{
		return ConstEdgeIterator(this);
	}

	// @@ Add half-edge iterator.

	bool isValid() const;

public:

	// Error status:
	mutable uint32_t errorCount;
	mutable uint32_t errorIndex0;
	mutable uint32_t errorIndex1;

private:

	bool canAddFace(const std::vector<uint32_t> &indexArray, uint32_t first, uint32_t num) const;
	bool canAddEdge(uint32_t i, uint32_t j) const;
	Edge *addEdge(uint32_t i, uint32_t j);

	Edge *findEdge(uint32_t i, uint32_t j) const;

	void linkBoundaryEdge(Edge *edge);
	Vertex *splitBoundaryEdge(Edge *edge, float t, const Vector3 &pos);
	void splitBoundaryEdge(Edge *edge, Vertex *vertex);

private:

	std::vector<Vertex *> m_vertexArray;
	std::vector<Edge *> m_edgeArray;
	std::vector<Face *> m_faceArray;

	struct Key {
		Key() {}
		Key(const Key &k) : p0(k.p0), p1(k.p1) {}
		Key(uint32_t v0, uint32_t v1) : p0(v0), p1(v1) {}
		void operator=(const Key &k)
		{
			p0 = k.p0;
			p1 = k.p1;
		}
		bool operator==(const Key &k) const
		{
			return p0 == k.p0 && p1 == k.p1;
		}

		uint32_t p0;
		uint32_t p1;
	};
	friend struct Hash<Mesh::Key>;

	std::unordered_map<Key, Edge *, Hash<Key>, Equal<Key> > m_edgeMap;

	uint32_t m_colocalVertexCount;

};

class MeshTopology
{
public:
	MeshTopology(const Mesh *mesh)
	{
		buildTopologyInfo(mesh);
	}

	/// Determine if the mesh is connected.
	bool isConnected() const
	{
		return m_connectedCount == 1;
	}

	/// Determine if the mesh is closed. (Each edge is shared by two faces)
	bool isClosed() const
	{
		return m_boundaryCount == 0;
	}

	/// Return true if the mesh has the topology of a disk.
	bool isDisk() const
	{
		return isConnected() && m_boundaryCount == 1/* && m_eulerNumber == 1*/;
	}

private:
	void buildTopologyInfo(const Mesh *mesh)
	{
		const uint32_t vertexCount = mesh->colocalVertexCount();
		const uint32_t faceCount = mesh->faceCount();
		const uint32_t edgeCount = mesh->edgeCount();
		nvDebug( "--- Building mesh topology:\n" );
		std::vector<uint32_t> stack(faceCount);
		BitArray bitFlags(faceCount);
		bitFlags.clearAll();
		// Compute connectivity.
		nvDebug( "---   Computing connectivity.\n" );
		m_connectedCount = 0;
		for (uint32_t f = 0; f < faceCount; f++ ) {
			if ( bitFlags.bitAt(f) == false ) {
				m_connectedCount++;
				stack.push_back( f );
				while ( !stack.empty() ) {
					const uint32_t top = stack.back();
					nvCheck(top != uint32_t(~0));
					stack.pop_back();
					if ( bitFlags.bitAt(top) == false ) {
						bitFlags.setBitAt(top);
						const Face *face = mesh->faceAt(top);
						const Edge *firstEdge = face->edge;
						const Edge *edge = firstEdge;
						do {
							const Face *neighborFace = edge->pair->face;
							if (neighborFace != NULL) {
								stack.push_back(neighborFace->id);
							}
							edge = edge->next;
						} while (edge != firstEdge);
					}
				}
			}
		}
		nvCheck(stack.empty());
		nvDebug( "---   %d connected components.\n", m_connectedCount );
		// Count boundary loops.
		nvDebug( "---   Counting boundary loops.\n" );
		m_boundaryCount = 0;
		bitFlags.resize(edgeCount);
		bitFlags.clearAll();
		// Don't forget to link the boundary otherwise this won't work.
		for (uint32_t e = 0; e < edgeCount; e++) {
			const Edge *startEdge = mesh->edgeAt(e);
			if (startEdge != NULL && startEdge->isBoundary() && bitFlags.bitAt(e) == false) {
				nvDebugCheck(startEdge->face != NULL);
				nvDebugCheck(startEdge->pair->face == NULL);
				startEdge = startEdge->pair;
				m_boundaryCount++;
				const Edge *edge = startEdge;
				do {
					bitFlags.setBitAt(edge->id / 2);
					edge = edge->next;
				} while (startEdge != edge);
			}
		}
		nvDebug("---   %d boundary loops found.\n", m_boundaryCount );
		// Compute euler number.
		m_eulerNumber = vertexCount - edgeCount + faceCount;
		nvDebug("---   Euler number: %d.\n", m_eulerNumber);
		// Compute genus. (only valid on closed connected surfaces)
		m_genus = -1;
		if ( isClosed() && isConnected() ) {
			m_genus = (2 - m_eulerNumber) / 2;
			nvDebug("---   Genus: %d.\n", m_genus);
		}
	}

private:

	///< Number of boundary loops.
	int m_boundaryCount;

	///< Number of connected components.
	int m_connectedCount;

	///< Euler number.
	int m_eulerNumber;

	/// Mesh genus.
	int m_genus;
};

float computeSurfaceArea(const Mesh *mesh);
float computeParametricArea(const Mesh *mesh);
uint32_t countMeshTriangles(const Mesh *mesh);
Mesh *unifyVertices(const Mesh *inputMesh);
Mesh *triangulate(const Mesh *inputMesh);

} //  namespace HalfEdge

/// Mersenne twister random number generator.
class MTRand
{
public:

	enum time_e { Time };
	enum { N = 624 };       // length of state vector
	enum { M = 397 };

	/// Constructor that uses the current time as the seed.
	MTRand( time_e )
	{
		seed((uint32_t )time(NULL));
	}

	/// Constructor that uses the given seed.
	MTRand( uint32_t s = 0 )
	{
		seed(s);
	}

	/// Provide a new seed.
	void seed( uint32_t s )
	{
		initialize(s);
		reload();
	}

	/// Get a random number between 0 - 65536.
	uint32_t get()
	{
		// Pull a 32-bit integer from the generator state
		// Every other access function simply transforms the numbers extracted here
		if ( left == 0 ) {
			reload();
		}
		left--;
		uint32_t s1;
		s1 = *next++;
		s1 ^= (s1 >> 11);
		s1 ^= (s1 <<  7) & 0x9d2c5680U;
		s1 ^= (s1 << 15) & 0xefc60000U;
		return ( s1 ^ (s1 >> 18) );
	};

	/// Get a random number on [0, max] interval.
	uint32_t getRange( uint32_t max )
	{
		if (max == 0) return 0;
		if (max == NV_UINT32_MAX) return get();
		const uint32_t np2 = nextPowerOfTwo( max + 1 ); // @@ This fails if max == NV_UINT32_MAX
		const uint32_t mask = np2 - 1;
		uint32_t n;
		do {
			n = get() & mask;
		} while ( n > max );
		return n;
	}


private:

	void initialize( uint32_t seed )
	{
		// Initialize generator state with seed
		// See Knuth TAOCP Vol 2, 3rd Ed, p.106 for multiplier.
		// In previous versions, most significant bits (MSBs) of the seed affect
		// only MSBs of the state array.  Modified 9 Jan 2002 by Makoto Matsumoto.
		uint32_t *s = state;
		uint32_t *r = state;
		int i = 1;
		*s++ = seed & 0xffffffffUL;
		for ( ; i < N; ++i ) {
			*s++ = ( 1812433253UL * ( *r ^ (*r >> 30) ) + i ) & 0xffffffffUL;
			r++;
		}
	}

	void reload()
	{
		// Generate N new values in state
		// Made clearer and faster by Matthew Bellew (matthew.bellew@home.com)
		uint32_t *p = state;
		int i;
		for ( i = N - M; i--; ++p )
			*p = twist( p[M], p[0], p[1] );
		for ( i = M; --i; ++p )
			*p = twist( p[M - N], p[0], p[1] );
		*p = twist( p[M - N], p[0], state[0] );
		left = N, next = state;
	}

	uint32_t hiBit( uint32_t u ) const
	{
		return u & 0x80000000U;
	}
	uint32_t loBit( uint32_t u ) const
	{
		return u & 0x00000001U;
	}
	uint32_t loBits( uint32_t u ) const
	{
		return u & 0x7fffffffU;
	}
	uint32_t mixBits( uint32_t u, uint32_t v ) const
	{
		return hiBit(u) | loBits(v);
	}
	uint32_t twist( uint32_t m, uint32_t s0, uint32_t s1 ) const
	{
		return m ^ (mixBits(s0, s1) >> 1) ^ ((~loBit(s1) + 1) & 0x9908b0dfU);
	}

private:

	uint32_t state[N];	// internal state
	uint32_t *next;	// next value to get from state
	int left;		// number of values left before reload needed

};

namespace morton {
// Code from ryg:
// http://fgiesen.wordpress.com/2009/12/13/decoding-morton-codes/

// Inverse of part1By1 - "delete" all odd-indexed bits
inline uint32_t compact1By1(uint32_t x)
{
	x &= 0x55555555;                  // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
	x = (x ^ (x >>  1)) & 0x33333333; // x = --fe --dc --ba --98 --76 --54 --32 --10
	x = (x ^ (x >>  2)) & 0x0f0f0f0f; // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
	x = (x ^ (x >>  4)) & 0x00ff00ff; // x = ---- ---- fedc ba98 ---- ---- 7654 3210
	x = (x ^ (x >>  8)) & 0x0000ffff; // x = ---- ---- ---- ---- fedc ba98 7654 3210
	return x;
}

// Inverse of part1By2 - "delete" all bits not at positions divisible by 3
inline uint32_t compact1By2(uint32_t x)
{
	x &= 0x09249249;                  // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
	x = (x ^ (x >>  2)) & 0x030c30c3; // x = ---- --98 ---- 76-- --54 ---- 32-- --10
	x = (x ^ (x >>  4)) & 0x0300f00f; // x = ---- --98 ---- ---- 7654 ---- ---- 3210
	x = (x ^ (x >>  8)) & 0xff0000ff; // x = ---- --98 ---- ---- ---- ---- 7654 3210
	x = (x ^ (x >> 16)) & 0x000003ff; // x = ---- ---- ---- ---- ---- --98 7654 3210
	return x;
}

inline uint32_t decodeMorton2X(uint32_t code)
{
	return compact1By1(code >> 0);
}

inline uint32_t decodeMorton2Y(uint32_t code)
{
	return compact1By1(code >> 1);
}

inline uint32_t decodeMorton3X(uint32_t code)
{
	return compact1By2(code >> 0);
}

inline uint32_t decodeMorton3Y(uint32_t code)
{
	return compact1By2(code >> 1);
}

inline uint32_t decodeMorton3Z(uint32_t code)
{
	return compact1By2(code >> 2);
}
} // namespace morton

// A simple, dynamic proximity grid based on Jon's code.
// Instead of storing pointers here I store indices.
struct ProximityGrid 
{
	void init(const Box &box, uint32_t count)
	{
		cellArray.clear();
		// Determine grid size.
		float cellWidth;
		Vector3 diagonal = box.extents() * 2.f;
		float volume = box.volume();
		if (equal(volume, 0)) {
			// Degenerate box, treat like a quad.
			Vector2 quad;
			if (diagonal.x < diagonal.y && diagonal.x < diagonal.z) {
				quad.x = diagonal.y;
				quad.y = diagonal.z;
			} else if (diagonal.y < diagonal.x && diagonal.y < diagonal.z) {
				quad.x = diagonal.x;
				quad.y = diagonal.z;
			} else {
				quad.x = diagonal.x;
				quad.y = diagonal.y;
			}
			float cellArea = quad.x * quad.y / count;
			cellWidth = sqrtf(cellArea); // pow(cellArea, 1.0f / 2.0f);
		} else {
			// Ideally we want one cell per point.
			float cellVolume = volume / count;
			cellWidth = powf(cellVolume, 1.0f / 3.0f);
		}
		nvDebugCheck(cellWidth != 0);
		sx = std::max(1, ftoi_ceil(diagonal.x / cellWidth));
		sy = std::max(1, ftoi_ceil(diagonal.y / cellWidth));
		sz = std::max(1, ftoi_ceil(diagonal.z / cellWidth));
		invCellSize.x = float(sx) / diagonal.x;
		invCellSize.y = float(sy) / diagonal.y;
		invCellSize.z = float(sz) / diagonal.z;
		cellArray.resize(sx * sy * sz);
		corner = box.minCorner; // @@ Align grid better?
	}

	inline int index_x(float x) const
	{
		return clamp(ftoi_floor((x - corner.x) * invCellSize.x),  0, sx - 1);
	}

	inline int index_y(float y) const
	{
		return clamp(ftoi_floor((y - corner.y) * invCellSize.y),  0, sy - 1);
	}

	inline int index_z(float z) const
	{
		return clamp(ftoi_floor((z - corner.z) * invCellSize.z),  0, sz - 1);
	}

	inline int index(int x, int y, int z) const
	{
		nvDebugCheck(x >= 0 && x < sx);
		nvDebugCheck(y >= 0 && y < sy);
		nvDebugCheck(z >= 0 && z < sz);
		int idx = (z * sy + y) * sx + x;
		nvDebugCheck(idx >= 0 && uint32_t(idx) < cellArray.size());
		return idx;
	}

	uint32_t mortonCount() const
	{
		uint64_t s = uint64_t(max3(sx, sy, sz));
		s = nextPowerOfTwo(s);
		if (s > 1024) {
			return uint32_t(s * s * min3(sx, sy, sz));
		}
		return uint32_t(s * s * s);
	}

	int mortonIndex(uint32_t code) const
	{
		uint32_t x, y, z;
		uint32_t s = uint32_t(max3(sx, sy, sz));
		if (s > 1024) {
			// Use layered two-dimensional morton order.
			s = nextPowerOfTwo(s);
			uint32_t layer = code / (s * s);
			code = code % (s * s);
			uint32_t layer_count = uint32_t(min3(sx, sy, sz));
			if (sx == layer_count) {
				x = layer;
				y = morton::decodeMorton2X(code);
				z = morton::decodeMorton2Y(code);
			} else if (sy == layer_count) {
				x = morton::decodeMorton2Y(code);
				y = layer;
				z = morton::decodeMorton2X(code);
			} else { /*if (sz == layer_count)*/
				x = morton::decodeMorton2X(code);
				y = morton::decodeMorton2Y(code);
				z = layer;
			}
		} else {
			x = morton::decodeMorton3X(code);
			y = morton::decodeMorton3Y(code);
			z = morton::decodeMorton3Z(code);
		}
		if (x >= uint32_t(sx) || y >= uint32_t(sy) || z >= uint32_t(sz)) {
			return -1;
		}
		return index(x, y, z);
	}

	inline void add(const Vector3 &pos, uint32_t key)
	{
		int x = index_x(pos.x);
		int y = index_y(pos.y);
		int z = index_z(pos.z);
		uint32_t idx = index(x, y, z);
		cellArray[idx].indexArray.push_back(key);
	}

	// Gather all points inside the given sphere.
	// Radius is assumed to be small, so we don't bother culling the cells.
	void gather(const Vector3 &position, float radius, std::vector<uint32_t> &indexArray)
	{
		int x0 = index_x(position.x - radius);
		int x1 = index_x(position.x + radius);
		int y0 = index_y(position.y - radius);
		int y1 = index_y(position.y + radius);
		int z0 = index_z(position.z - radius);
		int z1 = index_z(position.z + radius);
		for (int z = z0; z <= z1; z++) {
			for (int y = y0; y <= y1; y++) {
				for (int x = x0; x <= x1; x++) {
					int idx = index(x, y, z);
					indexArray.insert(indexArray.begin(), cellArray[idx].indexArray.begin(), cellArray[idx].indexArray.end());
				}
			}
		}
	}

	struct Cell {
		std::vector<uint32_t> indexArray;
	};

	std::vector<Cell> cellArray;

	Vector3 corner;
	Vector3 invCellSize;
	int sx, sy, sz;
};

// Based on Pierre Terdiman's and Michael Herf's source code.
// http://www.codercorner.com/RadixSortRevisited.htm
// http://www.stereopsis.com/radix.html
class RadixSort
{
public:
	RadixSort() : m_size(0), m_ranks(NULL), m_ranks2(NULL), m_validRanks(false) {}
	~RadixSort()
	{
		// Release everything
		free(m_ranks2);
		free(m_ranks);
	}

	RadixSort &sort(const float *input, uint32_t count)
	{
		if (input == NULL || count == 0) return *this;
		// Resize lists if needed
		if (count != m_size) {
			if (count > m_size) {
				m_ranks2 = (uint32_t *)realloc(m_ranks2, sizeof(uint32_t ) * count);
				m_ranks = (uint32_t *)realloc(m_ranks, sizeof(uint32_t ) * count);
			}
			m_size = count;
			m_validRanks = false;
		}
		if (count < 32) {
			insertionSort(input, count);
		} else {
			// @@ Avoid touching the input multiple times.
			for (uint32_t i = 0; i < count; i++) {
				FloatFlip((uint32_t &)input[i]);
			}
			radixSort<uint32_t>((const uint32_t *)input, count);
			for (uint32_t i = 0; i < count; i++) {
				IFloatFlip((uint32_t &)input[i]);
			}
		}
		return *this;
	}

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

	inline void FloatFlip(uint32_t &f)
	{
		int32_t mask = (int32_t(f) >> 31) | 0x80000000; // Warren Hunt, Manchor Ko.
		f ^= mask;
	}

	inline void IFloatFlip(uint32_t &f)
	{
		uint32_t mask = ((f >> 31) - 1) | 0x80000000; // Michael Herf.
		f ^= mask;
	}

	template<typename T>
	void createHistograms(const T *buffer, uint32_t count, uint32_t *histogram)
	{
		const uint32_t bucketCount = sizeof(T); // (8 * sizeof(T)) / log2(radix)
		// Init bucket pointers.
		uint32_t *h[bucketCount];
		for (uint32_t i = 0; i < bucketCount; i++) {
			h[i] = histogram + 256 * i;
		}
		// Clear histograms.
		memset(histogram, 0, 256 * bucketCount * sizeof(uint32_t ));
		// @@ Add support for signed integers.
		// Build histograms.
		const uint8_t *p = (const uint8_t *)buffer;  // @@ Does this break aliasing rules?
		const uint8_t *pe = p + count * sizeof(T);
		while (p != pe) {
			h[0][*p++]++, h[1][*p++]++, h[2][*p++]++, h[3][*p++]++;
			if (bucketCount == 8) h[4][*p++]++, h[5][*p++]++, h[6][*p++]++, h[7][*p++]++;
		}
	}

	template <typename T> void insertionSort(const T *input, uint32_t count)
	{
		if (!m_validRanks) {
			m_ranks[0] = 0;
			for (uint32_t i = 1; i != count; ++i) {
				int rank = m_ranks[i] = i;
				uint32_t j = i;
				while (j != 0 && input[rank] < input[m_ranks[j - 1]]) {
					m_ranks[j] = m_ranks[j - 1];
					--j;
				}
				if (i != j) {
					m_ranks[j] = rank;
				}
			}
			m_validRanks = true;
		} else {
			for (uint32_t i = 1; i != count; ++i) {
				int rank = m_ranks[i];
				uint32_t j = i;
				while (j != 0 && input[rank] < input[m_ranks[j - 1]]) {
					m_ranks[j] = m_ranks[j - 1];
					--j;
				}
				if (i != j) {
					m_ranks[j] = rank;
				}
			}
		}
	}

	template <typename T> void radixSort(const T *input, uint32_t count)
	{
		const uint32_t P = sizeof(T); // pass count
		// Allocate histograms & offsets on the stack
		uint32_t histogram[256 * P];
		uint32_t *link[256];
		createHistograms(input, count, histogram);
		// Radix sort, j is the pass number (0=LSB, P=MSB)
		for (uint32_t j = 0; j < P; j++) {
			// Pointer to this bucket.
			const uint32_t *h = &histogram[j * 256];
			const uint8_t *inputBytes = (const uint8_t *)input; // @@ Is this aliasing legal?
			inputBytes += j;
			if (h[inputBytes[0]] == count) {
				// Skip this pass, all values are the same.
				continue;
			}
			// Create offsets
			link[0] = m_ranks2;
			for (uint32_t i = 1; i < 256; i++) link[i] = link[i - 1] + h[i - 1];
			// Perform Radix Sort
			if (!m_validRanks) {
				for (uint32_t i = 0; i < count; i++) {
					*link[inputBytes[i * P]]++ = i;
				}
				m_validRanks = true;
			} else {
				for (uint32_t i = 0; i < count; i++) {
					const uint32_t idx = m_ranks[i];
					*link[inputBytes[idx * P]]++ = idx;
				}
			}
			// Swap pointers for next pass. Valid indices - the most recent ones - are in m_ranks after the swap.
			std::swap(m_ranks, m_ranks2);
		}
		// All values were equal, generate linear ranks.
		if (!m_validRanks) {
			for (uint32_t i = 0; i < count; i++) {
				m_ranks[i] = i;
			}
			m_validRanks = true;
		}
	}
};

namespace raster {
class ClippedTriangle
{
public:
	ClippedTriangle(Vector2::Arg a, Vector2::Arg b, Vector2::Arg c)
	{
		m_numVertices = 3;
		m_activeVertexBuffer = 0;
		m_verticesA[0] = a;
		m_verticesA[1] = b;
		m_verticesA[2] = c;
		m_vertexBuffers[0] = m_verticesA;
		m_vertexBuffers[1] = m_verticesB;
	}

	uint32_t vertexCount()
	{
		return m_numVertices;
	}

	const Vector2 *vertices()
	{
		return m_vertexBuffers[m_activeVertexBuffer];
	}

	inline void clipHorizontalPlane(float offset, float clipdirection)
	{
		Vector2 *v  = m_vertexBuffers[m_activeVertexBuffer];
		m_activeVertexBuffer ^= 1;
		Vector2 *v2 = m_vertexBuffers[m_activeVertexBuffer];
		v[m_numVertices] = v[0];
		float dy2,   dy1 = offset - v[0].y;
		int   dy2in, dy1in = clipdirection * dy1 >= 0;
		uint32_t  p = 0;
		for (uint32_t k = 0; k < m_numVertices; k++) {
			dy2   = offset - v[k + 1].y;
			dy2in = clipdirection * dy2 >= 0;
			if (dy1in) v2[p++] = v[k];
			if ( dy1in + dy2in == 1 ) { // not both in/out
				float dx = v[k + 1].x - v[k].x;
				float dy = v[k + 1].y - v[k].y;
				v2[p++] = Vector2(v[k].x + dy1 * (dx / dy), offset);
			}
			dy1 = dy2;
			dy1in = dy2in;
		}
		m_numVertices = p;
		//for (uint32_t k=0; k<m_numVertices; k++) printf("(%f, %f)\n", v2[k].x, v2[k].y); printf("\n");
	}

	inline void clipVerticalPlane(float offset, float clipdirection )
	{
		Vector2 *v  = m_vertexBuffers[m_activeVertexBuffer];
		m_activeVertexBuffer ^= 1;
		Vector2 *v2 = m_vertexBuffers[m_activeVertexBuffer];
		v[m_numVertices] = v[0];
		float dx2,   dx1   = offset - v[0].x;
		int   dx2in, dx1in = clipdirection * dx1 >= 0;
		uint32_t  p = 0;
		for (uint32_t k = 0; k < m_numVertices; k++) {
			dx2 = offset - v[k + 1].x;
			dx2in = clipdirection * dx2 >= 0;
			if (dx1in) v2[p++] = v[k];
			if ( dx1in + dx2in == 1 ) { // not both in/out
				float dx = v[k + 1].x - v[k].x;
				float dy = v[k + 1].y - v[k].y;
				v2[p++] = Vector2(offset, v[k].y + dx1 * (dy / dx));
			}
			dx1 = dx2;
			dx1in = dx2in;
		}
		m_numVertices = p;
	}

	void computeAreaCentroid()
	{
		Vector2 *v  = m_vertexBuffers[m_activeVertexBuffer];
		v[m_numVertices] = v[0];
		m_area = 0;
		float centroidx = 0, centroidy = 0;
		for (uint32_t k = 0; k < m_numVertices; k++) {
			// http://local.wasp.uwa.edu.au/~pbourke/geometry/polyarea/
			float f = v[k].x * v[k + 1].y - v[k + 1].x * v[k].y;
			m_area += f;
			centroidx += f * (v[k].x + v[k + 1].x);
			centroidy += f * (v[k].y + v[k + 1].y);
		}
		m_area = 0.5f * fabsf(m_area);
		if (m_area == 0) {
			m_centroid = Vector2(0.0f);
		} else {
			m_centroid = Vector2(centroidx / (6 * m_area), centroidy / (6 * m_area));
		}
	}

	void clipAABox(float x0, float y0, float x1, float y1)
	{
		clipVerticalPlane  ( x0, -1);
		clipHorizontalPlane( y0, -1);
		clipVerticalPlane  ( x1,  1);
		clipHorizontalPlane( y1,  1);
		computeAreaCentroid();
	}

	Vector2 centroid()
	{
		return m_centroid;
	}

	float area()
	{
		return m_area;
	}

private:
	Vector2 m_verticesA[7 + 1];
	Vector2 m_verticesB[7 + 1];
	Vector2 *m_vertexBuffers[2];
	uint32_t    m_numVertices;
	uint32_t    m_activeVertexBuffer;
	float   m_area;
	Vector2 m_centroid;
};

/// A triangle vertex.
struct Vertex
{
	Vector2 pos;	// Position.
	Vector3 tex;	// Texcoord. (Barycentric coordinate)
};

/// A callback to sample the environment. Return false to terminate rasterization.
typedef bool (* SamplingCallback)(void *param, int x, int y, Vector3::Arg bar, Vector3::Arg dx, Vector3::Arg dy, float coverage);

/// A triangle for rasterization.
struct Triangle
{
	Triangle(Vector2::Arg v0, Vector2::Arg v1, Vector2::Arg v2, Vector3::Arg t0, Vector3::Arg t1, Vector3::Arg t2)
	{
		// Init vertices.
		this->v1 = v0;
		this->v2 = v2;
		this->v3 = v1;
		// Set barycentric coordinates.
		this->t1 = t0;
		this->t2 = t2;
		this->t3 = t1;
		// make sure every triangle is front facing.
		flipBackface();
		// Compute deltas.
		valid = computeDeltas();
		computeUnitInwardNormals();
	}

	/// Compute texture space deltas.
	/// This method takes two edge vectors that form a basis, determines the
	/// coordinates of the canonic vectors in that basis, and computes the
	/// texture gradient that corresponds to those vectors.
	bool computeDeltas()
	{
		Vector2 e0 = v3 - v1;
		Vector2 e1 = v2 - v1;
		Vector3 de0 = t3 - t1;
		Vector3 de1 = t2 - t1;
		float denom = 1.0f / (e0.y * e1.x - e1.y * e0.x);
		if (!std::isfinite(denom)) {
			return false;
		}
		float lambda1 = - e1.y * denom;
		float lambda2 = e0.y * denom;
		float lambda3 = e1.x * denom;
		float lambda4 = - e0.x * denom;
		dx = de0 * lambda1 + de1 * lambda2;
		dy = de0 * lambda3 + de1 * lambda4;
		return true;
	}

	bool draw(const Vector2 &extents, bool enableScissors, SamplingCallback cb, void *param)
	{
		// 28.4 fixed-point coordinates
		const int Y1 = ftoi_round(16.0f * v1.y);
		const int Y2 = ftoi_round(16.0f * v2.y);
		const int Y3 = ftoi_round(16.0f * v3.y);
		const int X1 = ftoi_round(16.0f * v1.x);
		const int X2 = ftoi_round(16.0f * v2.x);
		const int X3 = ftoi_round(16.0f * v3.x);
		// Deltas
		const int DX12 = X1 - X2;
		const int DX23 = X2 - X3;
		const int DX31 = X3 - X1;
		const int DY12 = Y1 - Y2;
		const int DY23 = Y2 - Y3;
		const int DY31 = Y3 - Y1;
		// Fixed-point deltas
		const int FDX12 = DX12 << 4;
		const int FDX23 = DX23 << 4;
		const int FDX31 = DX31 << 4;
		const int FDY12 = DY12 << 4;
		const int FDY23 = DY23 << 4;
		const int FDY31 = DY31 << 4;
		int minx, miny, maxx, maxy;
		if (enableScissors) {
			int frustumX0 =  0 << 4;
			int frustumY0 =  0 << 4;
			int frustumX1 =  (int)extents.x << 4;
			int frustumY1 =  (int)extents.y << 4;
			// Bounding rectangle
			minx = (std::max(min3(X1, X2, X3), frustumX0) + 0xF) >> 4;
			miny = (std::max(min3(Y1, Y2, Y3), frustumY0) + 0xF) >> 4;
			maxx = (std::min(max3(X1, X2, X3), frustumX1) + 0xF) >> 4;
			maxy = (std::min(max3(Y1, Y2, Y3), frustumY1) + 0xF) >> 4;
		} else {
			// Bounding rectangle
			minx = (min3(X1, X2, X3) + 0xF) >> 4;
			miny = (min3(Y1, Y2, Y3) + 0xF) >> 4;
			maxx = (max3(X1, X2, X3) + 0xF) >> 4;
			maxy = (max3(Y1, Y2, Y3) + 0xF) >> 4;
		}
		// Block size, standard 8x8 (must be power of two)
		const int q = 8;
		// @@ This won't work when minx,miny are negative. This code path is not used. Leaving as is for now.
		nvCheck(minx >= 0);
		nvCheck(miny >= 0);
		// Start in corner of 8x8 block
		minx &= ~(q - 1);
		miny &= ~(q - 1);
		// Half-edge constants
		int C1 = DY12 * X1 - DX12 * Y1;
		int C2 = DY23 * X2 - DX23 * Y2;
		int C3 = DY31 * X3 - DX31 * Y3;
		// Correct for fill convention
		if (DY12 < 0 || (DY12 == 0 && DX12 > 0)) C1++;
		if (DY23 < 0 || (DY23 == 0 && DX23 > 0)) C2++;
		if (DY31 < 0 || (DY31 == 0 && DX31 > 0)) C3++;
		// Loop through blocks
		for (int y = miny; y < maxy; y += q) {
			for (int x = minx; x < maxx; x += q) {
				// Corners of block
				int x0 = x << 4;
				int x1 = (x + q - 1) << 4;
				int y0 = y << 4;
				int y1 = (y + q - 1) << 4;
				// Evaluate half-space functions
				bool a00 = C1 + DX12 * y0 - DY12 * x0 > 0;
				bool a10 = C1 + DX12 * y0 - DY12 * x1 > 0;
				bool a01 = C1 + DX12 * y1 - DY12 * x0 > 0;
				bool a11 = C1 + DX12 * y1 - DY12 * x1 > 0;
				int a = (a00 << 0) | (a10 << 1) | (a01 << 2) | (a11 << 3);
				bool b00 = C2 + DX23 * y0 - DY23 * x0 > 0;
				bool b10 = C2 + DX23 * y0 - DY23 * x1 > 0;
				bool b01 = C2 + DX23 * y1 - DY23 * x0 > 0;
				bool b11 = C2 + DX23 * y1 - DY23 * x1 > 0;
				int b = (b00 << 0) | (b10 << 1) | (b01 << 2) | (b11 << 3);
				bool c00 = C3 + DX31 * y0 - DY31 * x0 > 0;
				bool c10 = C3 + DX31 * y0 - DY31 * x1 > 0;
				bool c01 = C3 + DX31 * y1 - DY31 * x0 > 0;
				bool c11 = C3 + DX31 * y1 - DY31 * x1 > 0;
				int c = (c00 << 0) | (c10 << 1) | (c01 << 2) | (c11 << 3);
				// Skip block when outside an edge
				if (a == 0x0 || b == 0x0 || c == 0x0) continue;
				// Accept whole block when totally covered
				if (a == 0xF && b == 0xF && c == 0xF) {
					Vector3 texRow = t1 + dy * (y0 - v1.y) + dx * (x0 - v1.x);
					for (int iy = y; iy < y + q; iy++) {
						Vector3 tex = texRow;
						for (int ix = x; ix < x + q; ix++) {
							//Vector3 tex = t1 + dx * (ix - v1.x) + dy * (iy - v1.y);
							if (!cb(param, ix, iy, tex, dx, dy, 1.0)) {
								// early out.
								return false;
							}
							tex += dx;
						}
						texRow += dy;
					}
				} else { // Partially covered block
					int CY1 = C1 + DX12 * y0 - DY12 * x0;
					int CY2 = C2 + DX23 * y0 - DY23 * x0;
					int CY3 = C3 + DX31 * y0 - DY31 * x0;
					Vector3 texRow = t1 + dy * (y0 - v1.y) + dx * (x0 - v1.x);
					for (int iy = y; iy < y + q; iy++) {
						int CX1 = CY1;
						int CX2 = CY2;
						int CX3 = CY3;
						Vector3 tex = texRow;
						for (int ix = x; ix < x + q; ix++) {
							if (CX1 > 0 && CX2 > 0 && CX3 > 0) {
								if (!cb(param, ix, iy, tex, dx, dy, 1.0)) {
									// early out.
									return false;
								}
							}
							CX1 -= FDY12;
							CX2 -= FDY23;
							CX3 -= FDY31;
							tex += dx;
						}
						CY1 += FDX12;
						CY2 += FDX23;
						CY3 += FDX31;
						texRow += dy;
					}
				}
			}
		}
		return true;
	}

	// extents has to be multiple of BK_SIZE!!
	bool drawAA(const Vector2 &extents, bool enableScissors, SamplingCallback cb, void *param)
	{
		const float PX_INSIDE = 1.0f/sqrt(2.0f);
		const float PX_OUTSIDE = -1.0f/sqrt(2.0f);
		const float BK_SIZE = 8;
		const float BK_INSIDE = sqrt(BK_SIZE*BK_SIZE/2.0f);
		const float BK_OUTSIDE = -sqrt(BK_SIZE*BK_SIZE/2.0f);

		float minx, miny, maxx, maxy;
		if (enableScissors) {
			// Bounding rectangle
			minx = floorf(std::max(min3(v1.x, v2.x, v3.x), 0.0f));
			miny = floorf(std::max(min3(v1.y, v2.y, v3.y), 0.0f));
			maxx = ceilf( std::min(max3(v1.x, v2.x, v3.x), extents.x - 1.0f));
			maxy = ceilf( std::min(max3(v1.y, v2.y, v3.y), extents.y - 1.0f));
		} else {
			// Bounding rectangle
			minx = floorf(min3(v1.x, v2.x, v3.x));
			miny = floorf(min3(v1.y, v2.y, v3.y));
			maxx = ceilf( max3(v1.x, v2.x, v3.x));
			maxy = ceilf( max3(v1.y, v2.y, v3.y));
		}
		// There's no reason to align the blocks to the viewport, instead we align them to the origin of the triangle bounds.
		minx = floorf(minx);
		miny = floorf(miny);
		//minx = (float)(((int)minx) & (~((int)BK_SIZE - 1))); // align to blocksize (we don't need to worry about blocks partially out of viewport)
		//miny = (float)(((int)miny) & (~((int)BK_SIZE - 1)));
		minx += 0.5;
		miny += 0.5; // sampling at texel centers!
		maxx += 0.5;
		maxy += 0.5;
		// Half-edge constants
		float C1 = n1.x * (-v1.x) + n1.y * (-v1.y);
		float C2 = n2.x * (-v2.x) + n2.y * (-v2.y);
		float C3 = n3.x * (-v3.x) + n3.y * (-v3.y);
		// Loop through blocks
		for (float y0 = miny; y0 <= maxy; y0 += BK_SIZE) {
			for (float x0 = minx; x0 <= maxx; x0 += BK_SIZE) {
				// Corners of block
				float xc = (x0 + (BK_SIZE - 1) / 2.0f);
				float yc = (y0 + (BK_SIZE - 1) / 2.0f);
				// Evaluate half-space functions
				float aC = C1 + n1.x * xc + n1.y * yc;
				float bC = C2 + n2.x * xc + n2.y * yc;
				float cC = C3 + n3.x * xc + n3.y * yc;
				// Skip block when outside an edge
				if ( (aC <= BK_OUTSIDE) || (bC <= BK_OUTSIDE) || (cC <= BK_OUTSIDE) ) continue;
				// Accept whole block when totally covered
				if ( (aC >= BK_INSIDE) && (bC >= BK_INSIDE) && (cC >= BK_INSIDE) ) {
					Vector3 texRow = t1 + dy * (y0 - v1.y) + dx * (x0 - v1.x);
					for (float y = y0; y < y0 + BK_SIZE; y++) {
						Vector3 tex = texRow;
						for (float x = x0; x < x0 + BK_SIZE; x++) {
							if (!cb(param, (int)x, (int)y, tex, dx, dy, 1.0f)) {
								return false;
							}
							tex += dx;
						}
						texRow += dy;
					}
				} else { // Partially covered block
					float CY1 = C1 + n1.x * x0 + n1.y * y0;
					float CY2 = C2 + n2.x * x0 + n2.y * y0;
					float CY3 = C3 + n3.x * x0 + n3.y * y0;
					Vector3 texRow = t1 + dy * (y0 - v1.y) + dx * (x0 - v1.x);
					for (float y = y0; y < y0 + BK_SIZE; y++) { // @@ This is not clipping to scissor rectangle correctly.
						float CX1 = CY1;
						float CX2 = CY2;
						float CX3 = CY3;
						Vector3 tex = texRow;
						for (float x = x0; x < x0 + BK_SIZE; x++) { // @@ This is not clipping to scissor rectangle correctly.
							if (CX1 >= PX_INSIDE && CX2 >= PX_INSIDE && CX3 >= PX_INSIDE) {
								// pixel completely covered
								Vector3 tex = t1 + dx * (x - v1.x) + dy * (y - v1.y);
								if (!cb(param, (int)x, (int)y, tex, dx, dy, 1.0f)) {
									return false;
								}
							} else if ((CX1 >= PX_OUTSIDE) && (CX2 >= PX_OUTSIDE) && (CX3 >= PX_OUTSIDE)) {
								// triangle partially covers pixel. do clipping.
								ClippedTriangle ct(v1 - Vector2(x, y), v2 - Vector2(x, y), v3 - Vector2(x, y));
								ct.clipAABox(-0.5, -0.5, 0.5, 0.5);
								Vector2 centroid = ct.centroid();
								float area = ct.area();
								if (area > 0.0f) {
									Vector3 texCent = tex - dx * centroid.x - dy * centroid.y;
									//nvCheck(texCent.x >= -0.1f && texCent.x <= 1.1f); // @@ Centroid is not very exact...
									//nvCheck(texCent.y >= -0.1f && texCent.y <= 1.1f);
									//nvCheck(texCent.z >= -0.1f && texCent.z <= 1.1f);
									//Vector3 texCent2 = t1 + dx * (x - v1.x) + dy * (y - v1.y);
									if (!cb(param, (int)x, (int)y, texCent, dx, dy, area)) {
										return false;
									}
								}
							}
							CX1 += n1.x;
							CX2 += n2.x;
							CX3 += n3.x;
							tex += dx;
						}
						CY1 += n1.y;
						CY2 += n2.y;
						CY3 += n3.y;
						texRow += dy;
					}
				}
			}
		}
		return true;
	}

	void flipBackface()
	{
		// check if triangle is backfacing, if so, swap two vertices
		if ( ((v3.x - v1.x) * (v2.y - v1.y) - (v3.y - v1.y) * (v2.x - v1.x)) < 0 ) {
			Vector2 hv = v1;
			v1 = v2;
			v2 = hv; // swap pos
			Vector3 ht = t1;
			t1 = t2;
			t2 = ht; // swap tex
		}
	}

	// compute unit inward normals for each edge.
	void computeUnitInwardNormals()
	{
		n1 = v1 - v2;
		n1 = Vector2(-n1.y, n1.x);
		n1 = n1 * (1.0f / sqrtf(n1.x * n1.x + n1.y * n1.y));
		n2 = v2 - v3;
		n2 = Vector2(-n2.y, n2.x);
		n2 = n2 * (1.0f / sqrtf(n2.x * n2.x + n2.y * n2.y));
		n3 = v3 - v1;
		n3 = Vector2(-n3.y, n3.x);
		n3 = n3 * (1.0f / sqrtf(n3.x * n3.x + n3.y * n3.y));
	}

	// Vertices.
	Vector2 v1, v2, v3;
	Vector2 n1, n2, n3; // unit inward normals
	Vector3 t1, t2, t3;

	// Deltas.
	Vector3 dx, dy;

	float sign;
	bool valid;
};

enum Mode
{
	Mode_Nearest,
	Mode_Antialiased
};

// Process the given triangle. Returns false if rasterization was interrupted by the callback.
bool drawTriangle(Mode mode, Vector2::Arg extents, bool enableScissors, const Vector2 v[3], SamplingCallback cb, void *param);

// Process the given quad. Returns false if rasterization was interrupted by the callback.
bool drawQuad(Mode mode, Vector2::Arg extents, bool enableScissors, const Vector2 v[4], SamplingCallback cb, void *param);
} // namespace raster

namespace sparse {
// Full and sparse vector and matrix classes. BLAS subset.
// Pseudo-BLAS interface.
void saxpy(float a, const FullVector &x, FullVector &y);   // y = a * x + y
void copy(const FullVector &x, FullVector &y);
void scal(float a, FullVector &x);
float dot(const FullVector &x, const FullVector &y);


enum Transpose {
	NoTransposed = 0,
	Transposed = 1
};

/**
* Sparse matrix class. The matrix is assumed to be sparse and to have
* very few non-zero elements, for this reason it's stored in indexed
* format. To multiply column vectors efficiently, the matrix stores
* the elements in indexed-column order, there is a list of indexed
* elements for each row of the matrix. As with the FullVector the
* dimension of the matrix is constant.
**/
class Matrix
{
public:
	// An element of the sparse array.
	struct Coefficient
	{
		uint32_t x;  // column
		float v; // value
	};

	Matrix(uint32_t d) : m_width(d) { m_array.resize(d); }
	Matrix(uint32_t w, uint32_t h) : m_width(w) { m_array.resize(h); }
	Matrix(const Matrix &m) : m_width(m.m_width) { m_array = m.m_array; }

	const Matrix &operator=(const Matrix &m)
	{
		nvCheck(width() == m.width());
		nvCheck(height() == m.height());
		m_array = m.m_array;
		return *this;
	}

	uint32_t width() const { return m_width; }
	uint32_t height() const { return m_array.size(); }
	bool isSquare() const { return width() == height(); }

	// x is column, y is row
	float getCoefficient(uint32_t x, uint32_t y) const
	{
		nvDebugCheck( x < width() );
		nvDebugCheck( y < height() );
		const uint32_t count = m_array[y].size();
		for (uint32_t i = 0; i < count; i++) {
			if (m_array[y][i].x == x) return m_array[y][i].v;
		}
		return 0.0f;
	}

	void setCoefficient(uint32_t x, uint32_t y, float f)
	{
		nvDebugCheck( x < width() );
		nvDebugCheck( y < height() );
		const uint32_t count = m_array[y].size();
		for (uint32_t i = 0; i < count; i++) {
			if (m_array[y][i].x == x) {
				m_array[y][i].v = f;
				return;
			}
		}
		if (f != 0.0f) {
			Coefficient c = { x, f };
			m_array[y].push_back( c );
		}
	}

	void addCoefficient(uint32_t x, uint32_t y, float f)
	{
		nvDebugCheck( x < width() );
		nvDebugCheck( y < height() );
		if (f != 0.0f) {
			const uint32_t count = m_array[y].size();
			for (uint32_t i = 0; i < count; i++) {
				if (m_array[y][i].x == x) {
					m_array[y][i].v += f;
					return;
				}
			}
			Coefficient c = { x, f };
			m_array[y].push_back( c );
		}
	}

	void mulCoefficient(uint32_t x, uint32_t y, float f)
	{
		nvDebugCheck( x < width() );
		nvDebugCheck( y < height() );
		const uint32_t count = m_array[y].size();
		for (uint32_t i = 0; i < count; i++) {
			if (m_array[y][i].x == x) {
				m_array[y][i].v *= f;
				return;
			}
		}
		if (f != 0.0f) {
			Coefficient c = { x, f };
			m_array[y].push_back( c );
		}
	}


	float sumRow(uint32_t y) const
	{
		nvDebugCheck( y < height() );
		const uint32_t count = m_array[y].size();
		float sum = 0;
		for (uint32_t i = 0; i < count; i++) {
			sum += m_array[y][i].v;
		}
		return sum;
	}

	float dotRow(uint32_t y, const FullVector &v) const
	{
		nvDebugCheck( y < height() );
		const uint32_t count = m_array[y].size();
		float sum = 0;
		for (uint32_t i = 0; i < count; i++) {
			sum += m_array[y][i].v * v[m_array[y][i].x];
		}
		return sum;
	}

	void madRow(uint32_t y, float alpha, FullVector &v) const
	{
		nvDebugCheck(y < height());
		const uint32_t count = m_array[y].size();
		for (uint32_t i = 0; i < count; i++) {
			v[m_array[y][i].x] += alpha * m_array[y][i].v;
		}
	}


	void clearRow(uint32_t y)
	{
		nvDebugCheck( y < height() );
		m_array[y].clear();
	}

	void scaleRow(uint32_t y, float f)
	{
		nvDebugCheck( y < height() );
		const uint32_t count = m_array[y].size();
		for (uint32_t i = 0; i < count; i++) {
			m_array[y][i].v *= f;
		}
	}

	void normalizeRow(uint32_t y)
	{
		nvDebugCheck( y < height() );
		float norm = 0.0f;
		const uint32_t count = m_array[y].size();
		for (uint32_t i = 0; i < count; i++) {
			float f = m_array[y][i].v;
			norm += f * f;
		}
		scaleRow(y, 1.0f / sqrtf(norm));
	}


	void clearColumn(uint32_t x)
	{
		nvDebugCheck(x < width());
		for (uint32_t y = 0; y < height(); y++) {
			const uint32_t count = m_array[y].size();
			for (uint32_t e = 0; e < count; e++) {
				if (m_array[y][e].x == x) {
					m_array[y][e].v = 0.0f;
					break;
				}
			}
		}
	}

	void scaleColumn(uint32_t x, float f)
	{
		nvDebugCheck(x < width());
		for (uint32_t y = 0; y < height(); y++) {
			const uint32_t count = m_array[y].size();
			for (uint32_t e = 0; e < count; e++) {
				if (m_array[y][e].x == x) {
					m_array[y][e].v *= f;
					break;
				}
			}
		}
	}

	const std::vector<Coefficient> &getRow(uint32_t y) const { return m_array[y]; }

	bool isSymmetric() const
	{
		for (uint32_t y = 0; y < height(); y++) {
			const uint32_t count = m_array[y].size();
			for (uint32_t e = 0; e < count; e++) {
				const uint32_t x = m_array[y][e].x;
				if (x > y) {
					float v = m_array[y][e].v;
					if (!equal(getCoefficient(y, x), v)) {  // @@ epsilon
						return false;
					}
				}
			}
		}
		return true;
	}

private:

	/// Number of columns.
	const uint32_t m_width;

	/// Array of matrix elements.
	std::vector< std::vector<Coefficient> > m_array;

};

void transpose(const Matrix &A, Matrix &B);

void mult(const Matrix &M, const FullVector &x, FullVector &y);
void mult(Transpose TM, const Matrix &M, const FullVector &x, FullVector &y);

// y = alpha*A*x + beta*y
void sgemv(float alpha, const Matrix &A, const FullVector &x, float beta, FullVector &y);
void sgemv(float alpha, Transpose TA, const Matrix &A, const FullVector &x, float beta, FullVector &y);

void mult(const Matrix &A, const Matrix &B, Matrix &C);
void mult(Transpose TA, const Matrix &A, Transpose TB, const Matrix &B, Matrix &C);

// C = alpha*A*B + beta*C
void sgemm(float alpha, const Matrix &A, const Matrix &B, float beta, Matrix &C);
void sgemm(float alpha, Transpose TA, const Matrix &A, Transpose TB, const Matrix &B, float beta, Matrix &C);

} // namespace sparse

namespace solver {
// Linear solvers.
bool LeastSquaresSolver(const sparse::Matrix &A, const FullVector &b, FullVector &x, float epsilon = 1e-5f);
bool LeastSquaresSolver(const sparse::Matrix &A, const FullVector &b, FullVector &x, const uint32_t *lockedParameters, uint32_t lockedCount, float epsilon = 1e-5f);
} // namespace solver

class Atlas;
class Chart;

namespace param {
bool computeLeastSquaresConformalMap(HalfEdge::Mesh *mesh);
bool computeOrthogonalProjectionMap(HalfEdge::Mesh *mesh);
void computeSingleFaceMap(HalfEdge::Mesh *mesh);

// Dummy implementation of a priority queue using sort at insertion.
// - Insertion is o(n)
// - Smallest element goes at the end, so that popping it is o(1).
// - Resorting is n*log(n)
// @@ Number of elements in the queue is usually small, and we'd have to rebalance often. I'm not sure it's worth implementing a heap.
// @@ Searcing at removal would remove the need for sorting when priorities change.
struct PriorityQueue
{
	PriorityQueue(uint32_t size = UINT_MAX) : maxSize(size) {}

	void push(float priority, uint32_t face)
	{
		uint32_t i = 0;
		const uint32_t count = pairs.size();
		for (; i < count; i++) {
			if (pairs[i].priority > priority) break;
		}
		Pair p = { priority, face };
		pairs.insert(pairs.begin() + i, p);
		if (pairs.size() > maxSize) {
			pairs.erase(pairs.begin());
		}
	}

	// push face out of order, to be sorted later.
	void push(uint32_t face)
	{
		Pair p = { 0.0f, face };
		pairs.push_back(p);
	}

	uint32_t pop()
	{
		uint32_t f = pairs.back().face;
		pairs.pop_back();
		return f;
	}

	void sort()
	{
		//sort(pairs); // @@ My intro sort appears to be much slower than it should!
		std::sort(pairs.begin(), pairs.end());
	}

	void clear()
	{
		pairs.clear();
	}

	uint32_t count() const
	{
		return pairs.size();
	}

	float firstPriority() const
	{
		return pairs.back().priority;
	}

	const uint32_t maxSize;

	struct Pair
	{
		bool operator <(const Pair &p) const
		{
			return priority > p.priority;    // !! Sort in inverse priority order!
		}

		float priority;
		uint32_t face;
	};

	std::vector<Pair> pairs;
};

struct ChartBuildData
{
	ChartBuildData(int id) : id(id)
	{
		planeNormal = Vector3(0);
		centroid = Vector3(0);
		coneAxis = Vector3(0);
		coneAngle = 0;
		area = 0;
		boundaryLength = 0;
		normalSum = Vector3(0);
		centroidSum = Vector3(0);
	}

	int id;

	// Proxy info:
	Vector3 planeNormal;
	Vector3 centroid;
	Vector3 coneAxis;
	float coneAngle;

	float area;
	float boundaryLength;
	Vector3 normalSum;
	Vector3 centroidSum;

	std::vector<uint32_t> seeds;  // @@ These could be a pointers to the HalfEdge faces directly.
	std::vector<uint32_t> faces;
	PriorityQueue candidates;
};

struct SegmentationSettings
{
	SegmentationSettings()
	{
		// Charts have no area or boundary limits right now.
		maxChartArea = NV_FLOAT_MAX;
		maxBoundaryLength = NV_FLOAT_MAX;
		proxyFitMetricWeight = 1.0f;
		roundnessMetricWeight = 0.1f;
		straightnessMetricWeight = 0.25f;
		normalSeamMetricWeight = 1.0f;
		textureSeamMetricWeight = 0.1f;
	}

	float maxChartArea;
	float maxBoundaryLength;

	float proxyFitMetricWeight;
	float roundnessMetricWeight;
	float straightnessMetricWeight;
	float normalSeamMetricWeight;
	float textureSeamMetricWeight;
};

struct AtlasBuilder
{
	AtlasBuilder(const HalfEdge::Mesh *m);
	~AtlasBuilder();

	void markUnchartedFaces(const std::vector<uint32_t> &unchartedFaces);

	void computeShortestPaths();

	void placeSeeds(float threshold, uint32_t maxSeedCount);
	void createRandomChart(float threshold);

	void addFaceToChart(ChartBuildData *chart, uint32_t f, bool recomputeProxy = false);

	bool growCharts(float threshold, uint32_t faceCount);
	bool growChart(ChartBuildData *chart, float threshold, uint32_t faceCount);

	void resetCharts();

	void updateCandidates(ChartBuildData *chart, uint32_t face);

	void updateProxies();
	void updateProxy(ChartBuildData *chart);

	bool relocateSeeds();
	bool relocateSeed(ChartBuildData *chart);

	void updatePriorities(ChartBuildData *chart);

	float evaluatePriority(ChartBuildData *chart, uint32_t face);
	float evaluateProxyFitMetric(ChartBuildData *chart, uint32_t face);
	float evaluateDistanceToBoundary(ChartBuildData *chart, uint32_t face);
	float evaluateDistanceToSeed(ChartBuildData *chart, uint32_t face);
	float evaluateRoundnessMetric(ChartBuildData *chart, uint32_t face, float newBoundaryLength, float newChartArea);
	float evaluateStraightnessMetric(ChartBuildData *chart, uint32_t face);

	float evaluateNormalSeamMetric(ChartBuildData *chart, uint32_t f);
	float evaluateTextureSeamMetric(ChartBuildData *chart, uint32_t f);
	float evaluateSeamMetric(ChartBuildData *chart, uint32_t f);

	float evaluateChartArea(ChartBuildData *chart, uint32_t f);
	float evaluateBoundaryLength(ChartBuildData *chart, uint32_t f);
	Vector3 evaluateChartNormalSum(ChartBuildData *chart, uint32_t f);
	Vector3 evaluateChartCentroidSum(ChartBuildData *chart, uint32_t f);

	Vector3 computeChartCentroid(const ChartBuildData *chart);


	void fillHoles(float threshold);
	void mergeCharts();

	// @@ Cleanup.
	struct Candidate {
		uint32_t face;
		ChartBuildData *chart;
		float metric;
	};

	const Candidate &getBestCandidate() const;
	void removeCandidate(uint32_t f);
	void updateCandidate(ChartBuildData *chart, uint32_t f, float metric);

	void mergeChart(ChartBuildData *owner, ChartBuildData *chart, float sharedBoundaryLength);


	uint32_t chartCount() const
	{
		return chartArray.size();
	}
	const std::vector<uint32_t> &chartFaces(uint32_t i) const;

	const HalfEdge::Mesh *mesh;
	uint32_t facesLeft;
	std::vector<int> faceChartArray;
	std::vector<ChartBuildData *> chartArray;
	std::vector<float> shortestPaths;

	std::vector<float> edgeLengths;
	std::vector<float> faceAreas;

	std::vector<Candidate> candidateArray; //
	std::vector<uint32_t> faceCandidateArray; // Map face index to candidate index.

	MTRand rand;

	SegmentationSettings settings;
};

struct AtlasPacker
{
	AtlasPacker(Atlas *atlas) : m_atlas(atlas), m_bitmap(256, 256), m_width(0), m_height(0) {}
	void packCharts(int quality, float texelArea, bool blockAligned, bool conservative);

	float computeAtlasUtilization() const
	{
		const uint32_t w = m_width;
		const uint32_t h = m_height;
		nvDebugCheck(w <= m_bitmap.width());
		nvDebugCheck(h <= m_bitmap.height());
		uint32_t count = 0;
		for (uint32_t y = 0; y < h; y++) {
			for (uint32_t x = 0; x < w; x++) {
				count += m_bitmap.bitAt(x, y);
			}
		}
		return float(count) / (w * h);
	}

private:
	void findChartLocation(int quality, const BitMap *bitmap, Vector2::Arg extents, int w, int h, int *best_x, int *best_y, int *best_w, int *best_h, int *best_r);
	void findChartLocation_bruteForce(const BitMap *bitmap, Vector2::Arg extents, int w, int h, int *best_x, int *best_y, int *best_w, int *best_h, int *best_r);
	void findChartLocation_random(const BitMap *bitmap, Vector2::Arg extents, int w, int h, int *best_x, int *best_y, int *best_w, int *best_h, int *best_r, int minTrialCount);

	void drawChartBitmapDilate(const Chart *chart, BitMap *bitmap, int padding);
	void drawChartBitmap(const Chart *chart, BitMap *bitmap, const Vector2 &scale, const Vector2 &offset);

	bool canAddChart(const BitMap *bitmap, int atlas_w, int atlas_h, int offset_x, int offset_y, int r)
	{
		nvDebugCheck(r == 0 || r == 1);
		// Check whether the two bitmaps overlap.
		const int w = bitmap->width();
		const int h = bitmap->height();
		if (r == 0) {
			for (int y = 0; y < h; y++) {
				int yy = y + offset_y;
				if (yy >= 0) {
					for (int x = 0; x < w; x++) {
						int xx = x + offset_x;
						if (xx >= 0) {
							if (bitmap->bitAt(x, y)) {
								if (xx < atlas_w && yy < atlas_h) {
									if (m_bitmap.bitAt(xx, yy)) return false;
								}
							}
						}
					}
				}
			}
		} else if (r == 1) {
			for (int y = 0; y < h; y++) {
				int xx = y + offset_x;
				if (xx >= 0) {
					for (int x = 0; x < w; x++) {
						int yy = x + offset_y;
						if (yy >= 0) {
							if (bitmap->bitAt(x, y)) {
								if (xx < atlas_w && yy < atlas_h) {
									if (m_bitmap.bitAt(xx, yy)) return false;
								}
							}
						}
					}
				}
			}
		}
		return true;
	}

	void addChart(const BitMap *bitmap, int atlas_w, int atlas_h, int offset_x, int offset_y, int r)
	{
		nvDebugCheck(r == 0 || r == 1);
		// Check whether the two bitmaps overlap.
		const int w = bitmap->width();
		const int h = bitmap->height();
		if (r == 0) {
			for (int y = 0; y < h; y++) {
				int yy = y + offset_y;
				if (yy >= 0) {
					for (int x = 0; x < w; x++) {
						int xx = x + offset_x;
						if (xx >= 0) {
							if (bitmap->bitAt(x, y)) {
								if (xx < atlas_w && yy < atlas_h) {
									nvDebugCheck(m_bitmap.bitAt(xx, yy) == false);
									m_bitmap.setBitAt(xx, yy);
								}
							}
						}
					}
				}
			}
		} else if (r == 1) {
			for (int y = 0; y < h; y++) {
				int xx = y + offset_x;
				if (xx >= 0) {
					for (int x = 0; x < w; x++) {
						int yy = x + offset_y;
						if (yy >= 0) {
							if (bitmap->bitAt(x, y)) {
								if (xx < atlas_w && yy < atlas_h) {
									nvDebugCheck(m_bitmap.bitAt(xx, yy) == false);
									m_bitmap.setBitAt(xx, yy);
								}
							}
						}
					}
				}
			}
		}
	}

	static bool checkBitsCallback(void *param, int x, int y, Vector3::Arg, Vector3::Arg, Vector3::Arg, float)
	{
		BitMap *bitmap = (BitMap * )param;
		nvDebugCheck(bitmap->bitAt(x, y) == false);
		return true;
	}

	static bool setBitsCallback(void *param, int x, int y, Vector3::Arg, Vector3::Arg, Vector3::Arg, float area)
	{
		BitMap *bitmap = (BitMap * )param;
		if (area > 0.0) {
			bitmap->setBitAt(x, y);
		}
		return true;
	}

	Atlas *m_atlas;
	BitMap m_bitmap;
	RadixSort m_radix;
	uint32_t m_width;
	uint32_t m_height;
	MTRand m_rand;
};

// Estimate quality of existing parameterization.
class ParameterizationQuality
{
public:
	ParameterizationQuality()
	{
		m_totalTriangleCount = 0;
		m_flippedTriangleCount = 0;
		m_zeroAreaTriangleCount = 0;
		m_parametricArea = 0.0f;
		m_geometricArea = 0.0f;
		m_stretchMetric = 0.0f;
		m_maxStretchMetric = 0.0f;
		m_conformalMetric = 0.0f;
		m_authalicMetric = 0.0f;
	}

	ParameterizationQuality(const HalfEdge::Mesh *mesh)
	{
		nvDebugCheck(mesh != NULL);
		m_totalTriangleCount = 0;
		m_flippedTriangleCount = 0;
		m_zeroAreaTriangleCount = 0;
		m_parametricArea = 0.0f;
		m_geometricArea = 0.0f;
		m_stretchMetric = 0.0f;
		m_maxStretchMetric = 0.0f;
		m_conformalMetric = 0.0f;
		m_authalicMetric = 0.0f;
		const uint32_t faceCount = mesh->faceCount();
		for (uint32_t f = 0; f < faceCount; f++) {
			const HalfEdge::Face *face = mesh->faceAt(f);
			const HalfEdge::Vertex *vertex0 = NULL;
			Vector3 p[3];
			Vector2 t[3];
			for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance()) {
				const HalfEdge::Edge *edge = it.current();
				if (vertex0 == NULL) {
					vertex0 = edge->vertex;
					p[0] = vertex0->pos;
					t[0] = vertex0->tex;
				} else if (edge->to() != vertex0) {
					p[1] = edge->from()->pos;
					p[2] = edge->to()->pos;
					t[1] = edge->from()->tex;
					t[2] = edge->to()->tex;
					processTriangle(p, t);
				}
			}
		}
		if (m_flippedTriangleCount + m_zeroAreaTriangleCount == faceCount) {
			// If all triangles are flipped, then none is.
			m_flippedTriangleCount = 0;
		}
		nvDebugCheck(std::isfinite(m_parametricArea) && m_parametricArea >= 0);
		nvDebugCheck(std::isfinite(m_geometricArea) && m_geometricArea >= 0);
		nvDebugCheck(std::isfinite(m_stretchMetric));
		nvDebugCheck(std::isfinite(m_maxStretchMetric));
		nvDebugCheck(std::isfinite(m_conformalMetric));
		nvDebugCheck(std::isfinite(m_authalicMetric));
	}

	bool isValid() const
	{
		return m_flippedTriangleCount == 0; // @@ Does not test for self-overlaps.
	}

	float rmsStretchMetric() const
	{
		if (m_geometricArea == 0) return 0.0f;
		float normFactor = sqrtf(m_parametricArea / m_geometricArea);
		return sqrtf(m_stretchMetric / m_geometricArea) * normFactor;
	}

	float maxStretchMetric() const
	{
		if (m_geometricArea == 0) return 0.0f;
		float normFactor = sqrtf(m_parametricArea / m_geometricArea);
		return m_maxStretchMetric * normFactor;
	}

	float rmsConformalMetric() const
	{
		if (m_geometricArea == 0) return 0.0f;
		return sqrtf(m_conformalMetric / m_geometricArea);
	}

	float maxAuthalicMetric() const
	{
		if (m_geometricArea == 0) return 0.0f;
		return sqrtf(m_authalicMetric / m_geometricArea);
	}

	void operator+=(const ParameterizationQuality &pq)
	{
		m_totalTriangleCount += pq.m_totalTriangleCount;
		m_flippedTriangleCount += pq.m_flippedTriangleCount;
		m_zeroAreaTriangleCount += pq.m_zeroAreaTriangleCount;
		m_parametricArea += pq.m_parametricArea;
		m_geometricArea += pq.m_geometricArea;
		m_stretchMetric += pq.m_stretchMetric;
		m_maxStretchMetric = std::max(m_maxStretchMetric, pq.m_maxStretchMetric);
		m_conformalMetric += pq.m_conformalMetric;
		m_authalicMetric += pq.m_authalicMetric;
	}

private:
	void processTriangle(Vector3 q[3], Vector2 p[3])
	{
		m_totalTriangleCount++;
		// Evaluate texture stretch metric. See:
		// - "Texture Mapping Progressive Meshes", Sander, Snyder, Gortler & Hoppe
		// - "Mesh Parameterization: Theory and Practice", Siggraph'07 Course Notes, Hormann, Levy & Sheffer.
		float t1 = p[0].x;
		float s1 = p[0].y;
		float t2 = p[1].x;
		float s2 = p[1].y;
		float t3 = p[2].x;
		float s3 = p[2].y;
		float geometricArea = length(cross(q[1] - q[0], q[2] - q[0])) / 2;
		float parametricArea = ((s2 - s1) * (t3 - t1) - (s3 - s1) * (t2 - t1)) / 2;
		if (isZero(parametricArea)) {
			m_zeroAreaTriangleCount++;
			return;
		}
		Vector3 Ss = (q[0] * (t2 - t3) + q[1] * (t3 - t1) + q[2] * (t1 - t2)) / (2 * parametricArea);
		Vector3 St = (q[0] * (s3 - s2) + q[1] * (s1 - s3) + q[2] * (s2 - s1)) / (2 * parametricArea);
		float a = dot(Ss, Ss); // E
		float b = dot(Ss, St); // F
		float c = dot(St, St); // G
		// Compute eigen-values of the first fundamental form:
		float sigma1 = sqrtf(0.5f * std::max(0.0f, a + c - sqrtf(square(a - c) + 4 * square(b)))); // gamma uppercase, min eigenvalue.
		float sigma2 = sqrtf(0.5f * std::max(0.0f, a + c + sqrtf(square(a - c) + 4 * square(b)))); // gamma lowercase, max eigenvalue.
		nvCheck(sigma2 >= sigma1);
		// isometric: sigma1 = sigma2 = 1
		// conformal: sigma1 / sigma2 = 1
		// authalic: sigma1 * sigma2 = 1
		float rmsStretch = sqrtf((a + c) * 0.5f);
		float rmsStretch2 = sqrtf((square(sigma1) + square(sigma2)) * 0.5f);
		nvDebugCheck(equal(rmsStretch, rmsStretch2, 0.01f));
		if (parametricArea < 0.0f) {
			// Count flipped triangles.
			m_flippedTriangleCount++;
			parametricArea = fabsf(parametricArea);
		}
		m_stretchMetric += square(rmsStretch) * geometricArea;
		m_maxStretchMetric = std::max(m_maxStretchMetric, sigma2);
		if (!isZero(sigma1, 0.000001f)) {
			// sigma1 is zero when geometricArea is zero.
			m_conformalMetric += (sigma2 / sigma1) * geometricArea;
		}
		m_authalicMetric += (sigma1 * sigma2) * geometricArea;
		// Accumulate total areas.
		m_geometricArea += geometricArea;
		m_parametricArea += parametricArea;
		//triangleConformalEnergy(q, p);
	}

	uint32_t m_totalTriangleCount;
	uint32_t m_flippedTriangleCount;
	uint32_t m_zeroAreaTriangleCount;
	float m_parametricArea;
	float m_geometricArea;
	float m_stretchMetric;
	float m_maxStretchMetric;
	float m_conformalMetric;
	float m_authalicMetric;
};
} // namespace param

} // namespace nv

#endif // ATLAS_H
