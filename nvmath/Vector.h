// This code is in the public domain -- castanyo@yahoo.es

#pragma once
#ifndef NV_MATH_VECTOR_H
#define NV_MATH_VECTOR_H

#include <cmath>
#include "nvmath.h"
#include "nvcore/nvcore.h" // min, max

namespace nv
{
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

inline Vector2 add(Vector2::Arg a, Vector2::Arg b)
{
	return Vector2(a.x + b.x, a.y + b.y);
}
inline Vector2 operator+(Vector2::Arg a, Vector2::Arg b)
{
	return add(a, b);
}

inline Vector2 sub(Vector2::Arg a, Vector2::Arg b)
{
	return Vector2(a.x - b.x, a.y - b.y);
}
inline Vector2 operator-(Vector2::Arg a, Vector2::Arg b)
{
	return sub(a, b);
}

inline Vector2 scale(Vector2::Arg v, float s)
{
	return Vector2(v.x * s, v.y * s);
}

inline Vector2 scale(Vector2::Arg v, Vector2::Arg s)
{
	return Vector2(v.x * s.x, v.y * s.y);
}

inline Vector2 operator*(Vector2::Arg v, float s)
{
	return scale(v, s);
}

inline Vector2 operator*(Vector2::Arg v1, Vector2::Arg v2)
{
	return Vector2(v1.x * v2.x, v1.y * v2.y);
}

inline Vector2 operator*(float s, Vector2::Arg v)
{
	return scale(v, s);
}

inline Vector2 operator/(Vector2::Arg v, float s)
{
	return scale(v, 1.0f / s);
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

inline float inverseLength(Vector2::Arg v)
{
	return 1.0f / sqrtf(lengthSquared(v));
}

inline bool isNormalized(Vector2::Arg v, float epsilon = NV_NORMAL_EPSILON)
{
	return equal(length(v), 1, epsilon);
}

inline Vector2 normalize(Vector2::Arg v, float epsilon = NV_EPSILON)
{
	float l = length(v);
	nvDebugCheck(!isZero(l, epsilon));
	Vector2 n = scale(v, 1.0f / l);
	nvDebugCheck(isNormalized(n));
	return n;
}

inline Vector2 normalizeSafe(Vector2::Arg v, Vector2::Arg fallback, float epsilon = NV_EPSILON)
{
	float l = length(v);
	if (isZero(l, epsilon)) {
		return fallback;
	}
	return scale(v, 1.0f / l);
}

// Safe, branchless normalization from Andy Firth. All error checking ommitted.
// http://altdevblogaday.com/2011/08/21/practical-flt-point-tricks/
inline Vector2 normalizeFast(Vector2::Arg v)
{
	const float very_small_float = 1.0e-037f;
	float l = very_small_float + length(v);
	return scale(v, 1.0f / l);
}

inline bool equal(Vector2::Arg v1, Vector2::Arg v2, float epsilon = NV_EPSILON)
{
	return equal(v1.x, v2.x, epsilon) && equal(v1.y, v2.y, epsilon);
}

inline Vector2 min(Vector2::Arg a, Vector2::Arg b)
{
	return Vector2(min(a.x, b.x), min(a.y, b.y));
}

inline Vector2 max(Vector2::Arg a, Vector2::Arg b)
{
	return Vector2(max(a.x, b.x), max(a.y, b.y));
}

inline Vector2 clamp(Vector2::Arg v, float min, float max)
{
	return Vector2(clamp(v.x, min, max), clamp(v.y, min, max));
}

inline Vector2 saturate(Vector2::Arg v)
{
	return Vector2(saturate(v.x), saturate(v.y));
}

inline bool isFinite(Vector2::Arg v)
{
	return std::isfinite(v.x) && std::isfinite(v.y);
}

inline Vector2 validate(Vector2::Arg v, Vector2::Arg fallback = Vector2(0.0f))
{
	if (!isFinite(v)) return fallback;
	Vector2 vf = v;
	nv::floatCleanup(vf.component, 2);
	return vf;
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

inline Vector3 scale(Vector3::Arg v, float s)
{
	return Vector3(v.x * s, v.y * s, v.z * s);
}

inline Vector3 scale(Vector3::Arg v, Vector3::Arg s)
{
	return Vector3(v.x * s.x, v.y * s.y, v.z * s.z);
}

inline Vector3 operator*(Vector3::Arg v, float s)
{
	return scale(v, s);
}

inline Vector3 operator*(float s, Vector3::Arg v)
{
	return scale(v, s);
}

inline Vector3 operator*(Vector3::Arg v, Vector3::Arg s)
{
	return scale(v, s);
}

inline Vector3 operator/(Vector3::Arg v, float s)
{
	return scale(v, 1.0f / s);
}

/*inline Vector3 add_scaled(Vector3::Arg a, Vector3::Arg b, float s)
{
    return Vector3(a.x + b.x * s, a.y + b.y * s, a.z + b.z * s);
}*/

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

inline float inverseLength(Vector3::Arg v)
{
	return 1.0f / sqrtf(lengthSquared(v));
}

inline bool isNormalized(Vector3::Arg v, float epsilon = NV_NORMAL_EPSILON)
{
	return equal(length(v), 1, epsilon);
}

inline Vector3 normalize(Vector3::Arg v, float epsilon = NV_EPSILON)
{
	float l = length(v);
	nvDebugCheck(!isZero(l, epsilon));
	Vector3 n = scale(v, 1.0f / l);
	nvDebugCheck(isNormalized(n));
	return n;
}

inline Vector3 normalizeSafe(Vector3::Arg v, Vector3::Arg fallback, float epsilon = NV_EPSILON)
{
	float l = length(v);
	if (isZero(l, epsilon)) {
		return fallback;
	}
	return scale(v, 1.0f / l);
}

// Safe, branchless normalization from Andy Firth. All error checking ommitted.
// http://altdevblogaday.com/2011/08/21/practical-flt-point-tricks/
inline Vector3 normalizeFast(Vector3::Arg v)
{
	const float very_small_float = 1.0e-037f;
	float l = very_small_float + length(v);
	return scale(v, 1.0f / l);
}

inline bool equal(Vector3::Arg v1, Vector3::Arg v2, float epsilon = NV_EPSILON)
{
	return equal(v1.x, v2.x, epsilon) && equal(v1.y, v2.y, epsilon) && equal(v1.z, v2.z, epsilon);
}

inline Vector3 min(Vector3::Arg a, Vector3::Arg b)
{
	return Vector3(min(a.x, b.x), min(a.y, b.y), min(a.z, b.z));
}

inline Vector3 max(Vector3::Arg a, Vector3::Arg b)
{
	return Vector3(max(a.x, b.x), max(a.y, b.y), max(a.z, b.z));
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

inline Vector3 validate(Vector3::Arg v, Vector3::Arg fallback = Vector3(0.0f))
{
	if (!isFinite(v)) return fallback;
	Vector3 vf = v;
	nv::floatCleanup(vf.component, 3);
	return vf;
}

inline Vector3 reflect(Vector3::Arg v, Vector3::Arg n)
{
	return v - (2 * dot(v, n)) * n;
}

template <>
inline uint32_t hash(const Vector3 &v, uint32_t h)
{
	return sdbmFloatHash(v.component, 3, h);
}

} // nv namespace

#endif // NV_MATH_VECTOR_H
