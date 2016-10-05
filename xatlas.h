// This code is in the public domain -- castanyo@yahoo.es

#pragma once
#ifndef ATLAS_H
#define ATLAS_H

#include <algorithm>
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

} // namespace nv

#endif // ATLAS_H
