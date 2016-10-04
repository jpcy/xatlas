// This code is in the public domain -- castanyo@yahoo.es

#pragma once
#ifndef NV_MATH_H
#define NV_MATH_H

#include "nvcore/nvcore.h"

#include <cmath>
#include <float.h>
#include <math.h>
#include <time.h>

#ifndef PI
#define PI                  float(3.1415926535897932384626433833)
#endif

#define NV_EPSILON          (0.0001f)
#define NV_NORMAL_EPSILON   (0.001f)

namespace nv
{
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

void convexHull(const std::vector<Vector2> &input, std::vector<Vector2> &output, float epsilon = 0);

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
struct ProximityGrid {
	ProximityGrid();

	void init(const Box &box, uint32_t count);

	int index_x(float x) const;
	int index_y(float y) const;
	int index_z(float z) const;
	int index(int x, int y, int z) const;

	uint32_t mortonCount() const;
	int mortonIndex(uint32_t code) const;

	void add(const Vector3 &pos, uint32_t key);

	void gather(const Vector3 &pos, float radius, std::vector<uint32_t> &indices);

	struct Cell {
		std::vector<uint32_t> indexArray;
	};

	std::vector<Cell> cellArray;

	Vector3 corner;
	Vector3 invCellSize;
	int sx, sy, sz;
};

inline int ProximityGrid::index_x(float x) const
{
	return clamp(ftoi_floor((x - corner.x) * invCellSize.x),  0, sx - 1);
}

inline int ProximityGrid::index_y(float y) const
{
	return clamp(ftoi_floor((y - corner.y) * invCellSize.y),  0, sy - 1);
}

inline int ProximityGrid::index_z(float z) const
{
	return clamp(ftoi_floor((z - corner.z) * invCellSize.z),  0, sz - 1);
}

inline int ProximityGrid::index(int x, int y, int z) const
{
	nvDebugCheck(x >= 0 && x < sx);
	nvDebugCheck(y >= 0 && y < sy);
	nvDebugCheck(z >= 0 && z < sz);
	int idx = (z * sy + y) * sx + x;
	nvDebugCheck(idx >= 0 && uint32_t(idx) < cellArray.size());
	return idx;
}

inline void ProximityGrid::add(const Vector3 &pos, uint32_t key)
{
	int x = index_x(pos.x);
	int y = index_y(pos.y);
	int z = index_z(pos.z);
	uint32_t idx = index(x, y, z);
	cellArray[idx].indexArray.push_back(key);
}

class SparseMatrix;
class FullVector;

namespace solver {
// Linear solvers.
bool LeastSquaresSolver(const SparseMatrix &A, const FullVector &b, FullVector &x, float epsilon = 1e-5f);
bool LeastSquaresSolver(const SparseMatrix &A, const FullVector &b, FullVector &x, const uint32_t *lockedParameters, uint32_t lockedCount, float epsilon = 1e-5f);
} // namespace solver

} // namespace nv

#endif // NV_MATH_H
