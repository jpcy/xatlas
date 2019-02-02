/*
Copyright (c) 2018 Jonathan Young
Copyright (c) 2013 Thekla, Inc
Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include <memory>
#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#undef min
#undef max
#include "xatlas.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define XA_STR(x) #x
#define XA_XSTR(x) XA_STR(x)

#ifndef XA_ASSERT
#define XA_ASSERT(exp) if (!(exp)) { XA_PRINT(0, "\rASSERT: %s %s %d\n", XA_XSTR(exp), __FILE__, __LINE__); }
#endif

#ifndef XA_DEBUG_ASSERT
#define XA_DEBUG_ASSERT(exp) assert(exp)
#endif

#define XA_ALLOC(type) (type *)internal::Realloc(NULL, sizeof(type), __FILE__, __LINE__)
#define XA_ALLOC_ARRAY(type, num) (type *)internal::Realloc(NULL, sizeof(type) * num, __FILE__, __LINE__)
#define XA_REALLOC(ptr, type, num) (type *)internal::Realloc(ptr, sizeof(type) * num, __FILE__, __LINE__)
#define XA_FREE(ptr) internal::Realloc(ptr, 0, __FILE__, __LINE__)
#define XA_NEW(type, ...) new (XA_ALLOC(type)) type(__VA_ARGS__)

#ifndef XA_PRINT
#define XA_PRINT(flags, ...) \
	if (xatlas::internal::s_print && (flags == 0 || (xatlas::internal::s_printFlags & flags) != 0)) \
		xatlas::internal::s_print(__VA_ARGS__);
#endif

#define XA_EPSILON          (0.0001f)
#define XA_NORMAL_EPSILON   (0.001f)

namespace xatlas {
namespace internal {

static ReallocFunc s_realloc = realloc;
static int s_printFlags = 0;
static PrintFunc s_print = printf;

#define XA_DEBUG_HEAP 0
#define XA_DEBUG_EXPORT_OBJ 0
#define XA_DEBUG_EXPORT_OBJ_INDIVIDUAL_CHARTS 0

#if XA_DEBUG_HEAP
struct AllocHeader
{
	size_t size;
	const char *file;
	int line;
	AllocHeader *prev, *next;
};

static AllocHeader *s_allocRoot = NULL;
static size_t s_allocTotalSize = 0;
static size_t s_allocPeakSize = 0;

static void *Realloc(void *ptr, size_t size, const char *file, int line)
{
	if (!size && !ptr)
		return NULL;
	uint8_t *realPtr = NULL;
	AllocHeader *header = NULL;
	if (ptr) {
		realPtr = ((uint8_t *)ptr) - sizeof(AllocHeader);
		header = (AllocHeader *)realPtr;
	}
	if (!size || realPtr) {
		// free or realloc, either way, remove.
		s_allocTotalSize -= header->size;
		if (header->prev)
			header->prev->next = header->next;
		else
			s_allocRoot = header->next;
		if (header->next)
			header->next->prev = header->prev;
	}
	if (!size)
		return s_realloc(realPtr, 0); // free
	size += sizeof(AllocHeader);
	uint8_t *newPtr = (uint8_t *)s_realloc(realPtr, size);
	if (!newPtr)
		return NULL;
	header = (AllocHeader *)newPtr;
	header->size = size;
	header->file = file;
	header->line = line;
	if (!s_allocRoot) {
		s_allocRoot = header;
		header->prev = header->next = 0;
	} else {
		header->prev = NULL;
		header->next = s_allocRoot;
		s_allocRoot = header;
		header->next->prev = header;
	}
	s_allocTotalSize += size;
	if (s_allocTotalSize > s_allocPeakSize)
		s_allocPeakSize = s_allocTotalSize;
	return newPtr + sizeof(AllocHeader);
}

static void ReportAllocs()
{
	AllocHeader *header = s_allocRoot;
	while (header) {
		printf("Leak: %d bytes %s %d\n", header->size, header->file, header->line);
		header = header->next;
	}
	printf("%0.2fMB peak memory usage\n", s_allocPeakSize / 1024.0f / 1024.0f);
}
#else
static void *Realloc(void *ptr, size_t size, const char * /*file*/, int /*line*/)
{
	void *mem = s_realloc(ptr, size);
	if (size > 0) {
		XA_DEBUG_ASSERT(mem);
	}
	return mem;
}
#endif

static int align(int x, int a)
{
	return (x + a - 1) & ~(a - 1);
}

/// Return the maximum of the three arguments.
template <typename T>
static T max3(const T &a, const T &b, const T &c)
{
	return std::max(a, std::max(b, c));
}

/// Return the maximum of the three arguments.
template <typename T>
static T min3(const T &a, const T &b, const T &c)
{
	return std::min(a, std::min(b, c));
}

/// Clamp between two values.
template <typename T>
static T clamp(const T &x, const T &a, const T &b)
{
	return std::min(std::max(x, a), b);
}

// Robust floating point comparisons:
// http://realtimecollisiondetection.net/blog/?p=89
static bool equal(const float f0, const float f1, const float epsilon = XA_EPSILON)
{
	//return fabs(f0-f1) <= epsilon;
	return fabs(f0 - f1) <= epsilon * max3(1.0f, fabsf(f0), fabsf(f1));
}

static int ftoi_ceil(float val)
{
	return (int)ceilf(val);
}

static int ftoi_round(float f)
{
	return int(floorf(f + 0.5f));
}

static bool isZero(const float f, const float epsilon = XA_EPSILON)
{
	return fabs(f) <= epsilon;
}

static float square(float f)
{
	return f * f;
}

/** Return the next power of two.
* @see http://graphics.stanford.edu/~seander/bithacks.html
* @warning Behaviour for 0 is undefined.
* @note isPowerOfTwo(x) == true -> nextPowerOfTwo(x) == x
* @note nextPowerOfTwo(x) = 2 << log2(x-1)
*/
static uint32_t nextPowerOfTwo(uint32_t x)
{
	XA_DEBUG_ASSERT( x != 0 );
	// On modern CPUs this is supposed to be as fast as using the bsr instruction.
	x--;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	return x + 1;
}

static uint32_t sdbmHash(const void *data_in, uint32_t size, uint32_t h = 5381)
{
	const uint8_t *data = (const uint8_t *) data_in;
	uint32_t i = 0;
	while (i < size) {
		h = (h << 16) + (h << 6) - h + (uint32_t ) data[i++];
	}
	return h;
}

template <typename T>
static uint32_t hash(const T &t, uint32_t h = 5381)
{
	return sdbmHash(&t, sizeof(T), h);
}

// Functors for hash table:
template <typename Key> struct Hash
{
	uint32_t operator()(const Key &k) const { return hash(k); }
};

template <typename Key> struct Equal
{
	bool operator()(const Key &k0, const Key &k1) const { return k0 == k1; }
};

class Vector2
{
public:
	typedef Vector2 const &Arg;

	Vector2() {}
	explicit Vector2(float f) : x(f), y(f) {}
	Vector2(float x, float y): x(x), y(y) {}
	Vector2(Vector2::Arg v) : x(v.x), y(v.y) {}

	const Vector2 &operator=(Vector2::Arg v)
	{
		x = v.x;
		y = v.y;
		return
			*this;
	}

	Vector2 operator-() const
	{
		return Vector2(-x, -y);
	}

	void operator+=(Vector2::Arg v)
	{
		x += v.x;
		y += v.y;
	}

	void operator-=(Vector2::Arg v)
	{
		x -= v.x;
		y -= v.y;
	}

	void operator*=(float s)
	{
		x *= s;
		y *= s;
	}

	void operator*=(Vector2::Arg v)
	{
		x *= v.x;
		y *= v.y;
	}

	friend bool operator==(Vector2::Arg a, Vector2::Arg b)
	{
		return a.x == b.x && a.y == b.y;
	}

	friend bool operator!=(Vector2::Arg a, Vector2::Arg b)
	{
		return a.x != b.x || a.y != b.y;
	}

	float x, y;
};

static Vector2 operator+(Vector2::Arg a, Vector2::Arg b)
{
	return Vector2(a.x + b.x, a.y + b.y);
}

static Vector2 operator-(Vector2::Arg a, Vector2::Arg b)
{
	return Vector2(a.x - b.x, a.y - b.y);
}

static Vector2 operator*(Vector2::Arg v, float s)
{
	return Vector2(v.x * s, v.y * s);
}

static Vector2 operator*(Vector2::Arg v1, Vector2::Arg v2)
{
	return Vector2(v1.x * v2.x, v1.y * v2.y);
}

static Vector2 lerp(Vector2::Arg v1, Vector2::Arg v2, float t)
{
	const float s = 1.0f - t;
	return Vector2(v1.x * s + t * v2.x, v1.y * s + t * v2.y);
}

static float dot(Vector2::Arg a, Vector2::Arg b)
{
	return a.x * b.x + a.y * b.y;
}

static float lengthSquared(Vector2::Arg v)
{
	return v.x * v.x + v.y * v.y;
}

static float length(Vector2::Arg v)
{
	return sqrtf(lengthSquared(v));
}

#ifdef _DEBUG
static bool isNormalized(Vector2::Arg v, float epsilon = XA_NORMAL_EPSILON)
{
	return equal(length(v), 1, epsilon);
}
#endif

static Vector2 normalize(Vector2::Arg v, float epsilon = XA_EPSILON)
{
	float l = length(v);
	XA_DEBUG_ASSERT(!isZero(l, epsilon));
#ifdef NDEBUG
	epsilon = epsilon; // silence unused parameter warning
#endif
	Vector2 n = v * (1.0f / l);
	XA_DEBUG_ASSERT(isNormalized(n));
	return n;
}

static bool equal(Vector2::Arg v1, Vector2::Arg v2, float epsilon = XA_EPSILON)
{
	return equal(v1.x, v2.x, epsilon) && equal(v1.y, v2.y, epsilon);
}

static Vector2 min(Vector2::Arg a, Vector2::Arg b)
{
	return Vector2(std::min(a.x, b.x), std::min(a.y, b.y));
}

static Vector2 max(Vector2::Arg a, Vector2::Arg b)
{
	return Vector2(std::max(a.x, b.x), std::max(a.y, b.y));
}

static bool isFinite(Vector2::Arg v)
{
	return std::isfinite(v.x) && std::isfinite(v.y);
}

// Note, this is the area scaled by 2!
static float triangleArea(Vector2::Arg v0, Vector2::Arg v1)
{
	return (v0.x * v1.y - v0.y * v1.x); // * 0.5f;
}

static float triangleArea(Vector2::Arg a, Vector2::Arg b, Vector2::Arg c)
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

static bool pointInTriangle(const Vector2 &p, const Vector2 &a, const Vector2 &b, const Vector2 &c)
{
	return triangleArea(a, b, p) >= 0.00001f && triangleArea(b, c, p) >= 0.00001f && triangleArea(c, a, p) >= 0.00001f;
}

class Vector3
{
public:
	typedef Vector3 const &Arg;

	Vector3() {}
	explicit Vector3(float f) : x(f), y(f), z(f) {}
	Vector3(float x, float y, float z) : x(x), y(y), z(z) {}
	Vector3(Vector2::Arg v, float z) : x(v.x), y(v.y), z(z) {}
	Vector3(Vector3::Arg v) : x(v.x), y(v.y), z(v.z) {}

	const Vector3 &operator=(Vector3::Arg v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
		return *this;
	}

	Vector2 xy() const
	{
		return Vector2(x, y);
	}

	Vector3 operator-() const
	{
		return Vector3(-x, -y, -z);
	}

	void operator+=(Vector3::Arg v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
	}

	void operator-=(Vector3::Arg v)
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;
	}

	void operator*=(float s)
	{
		x *= s;
		y *= s;
		z *= s;
	}

	void operator/=(float s)
	{
		float is = 1.0f / s;
		x *= is;
		y *= is;
		z *= is;
	}

	void operator*=(Vector3::Arg v)
	{
		x *= v.x;
		y *= v.y;
		z *= v.z;
	}

	void operator/=(Vector3::Arg v)
	{
		x /= v.x;
		y /= v.y;
		z /= v.z;
	}

	friend bool operator==(Vector3::Arg a, Vector3::Arg b)
	{
		return a.x == b.x && a.y == b.y && a.z == b.z;
	}

	friend bool operator!=(Vector3::Arg a, Vector3::Arg b)
	{
		return a.x != b.x || a.y != b.y || a.z != b.z;
	}

	float x, y, z;
};

static Vector3 add(Vector3::Arg a, Vector3::Arg b)
{
	return Vector3(a.x + b.x, a.y + b.y, a.z + b.z);
}

static Vector3 operator+(Vector3::Arg a, Vector3::Arg b)
{
	return add(a, b);
}

static Vector3 sub(Vector3::Arg a, Vector3::Arg b)
{
	return Vector3(a.x - b.x, a.y - b.y, a.z - b.z);
}

static Vector3 operator-(Vector3::Arg a, Vector3::Arg b)
{
	return sub(a, b);
}

static Vector3 cross(Vector3::Arg a, Vector3::Arg b)
{
	return Vector3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

static Vector3 operator*(Vector3::Arg v, float s)
{
	return Vector3(v.x * s, v.y * s, v.z * s);
}

static Vector3 operator*(float s, Vector3::Arg v)
{
	return Vector3(v.x * s, v.y * s, v.z * s);
}

static Vector3 operator/(Vector3::Arg v, float s)
{
	return v * (1.0f / s);
}

static Vector3 lerp(Vector3::Arg v1, Vector3::Arg v2, float t)
{
	const float s = 1.0f - t;
	return Vector3(v1.x * s + t * v2.x, v1.y * s + t * v2.y, v1.z * s + t * v2.z);
}

static float dot(Vector3::Arg a, Vector3::Arg b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

static float lengthSquared(Vector3::Arg v)
{
	return v.x * v.x + v.y * v.y + v.z * v.z;
}

static float length(Vector3::Arg v)
{
	return sqrtf(lengthSquared(v));
}

static bool isNormalized(Vector3::Arg v, float epsilon = XA_NORMAL_EPSILON)
{
	return equal(length(v), 1, epsilon);
}

static Vector3 normalize(Vector3::Arg v, float epsilon = XA_EPSILON)
{
	float l = length(v);
	XA_DEBUG_ASSERT(!isZero(l, epsilon));
#ifdef NDEBUG
	epsilon = epsilon; // silence unused parameter warning
#endif
	Vector3 n = v * (1.0f / l);
	XA_DEBUG_ASSERT(isNormalized(n));
	return n;
}

static Vector3 normalizeSafe(Vector3::Arg v, Vector3::Arg fallback, float epsilon = XA_EPSILON)
{
	float l = length(v);
	if (isZero(l, epsilon)) {
		return fallback;
	}
	return v * (1.0f / l);
}

static bool equal(const Vector3 &v0, const Vector3 &v1, float epsilon = XA_EPSILON)
{
	return fabs(v0.x - v1.x) <= epsilon && fabs(v0.y - v1.y) <= epsilon && fabs(v0.z - v1.z) <= epsilon;
}

#ifdef _DEBUG
bool isFinite(Vector3::Arg v)
{
	return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}
#endif

struct Vector3Hash
{
	uint32_t operator()(const Vector3 &v) const
	{
		int32_t data[3];
		data[0] = (int32_t)(v.x * 100.0f);
		data[1] = (int32_t)(v.y * 100.0f);
		data[2] = (int32_t)(v.z * 100.0f);
		return sdbmHash(data, sizeof(data));
	}
};

struct Vector3Equal
{
	bool operator()(const Vector3 &v0, const Vector3 &v1) const
	{
		return equal(v0, v1);
	}
};

template <typename T>
static void construct_range(T * ptr, uint32_t new_size, uint32_t old_size) {
	for (uint32_t i = old_size; i < new_size; i++) {
		new(ptr+i) T; // placement new
	}
}

template <typename T>
static void construct_range(T * ptr, uint32_t new_size, uint32_t old_size, const T & elem) {
	for (uint32_t i = old_size; i < new_size; i++) {
		new(ptr+i) T(elem); // placement new
	}
}

template <typename T>
static void construct_range(T * ptr, uint32_t new_size, uint32_t old_size, const T * src) {
	for (uint32_t i = old_size; i < new_size; i++) {
		new(ptr+i) T(src[i]); // placement new
	}
}

template <typename T>
static void destroy_range(T * ptr, uint32_t new_size, uint32_t old_size) {
	for (uint32_t i = new_size; i < old_size; i++) {
		(ptr+i)->~T(); // Explicit call to the destructor
	}
}

/**
* Replacement for std::vector that is easier to debug and provides
* some nice foreach enumerators. 
*/
template<typename T>
class Array {
public:
	typedef uint32_t size_type;

	Array() : m_buffer(NULL), m_capacity(0), m_size(0) {}

	Array(const Array & a) : m_buffer(NULL), m_capacity(0), m_size(0) {
		copy(a.m_buffer, a.m_size);
	}

	Array(const T * ptr, uint32_t num) : m_buffer(NULL), m_capacity(0), m_size(0) {
		copy(ptr, num);
	}

	explicit Array(uint32_t capacity) : m_buffer(NULL), m_capacity(0), m_size(0) {
		setArrayCapacity(capacity);
	}

	~Array() {
		clear();
		XA_FREE(m_buffer);
	}

	const T & operator[]( uint32_t index ) const
	{
		XA_DEBUG_ASSERT(index < m_size);
		return m_buffer[index];
	}
	
	T & operator[] ( uint32_t index )
	{
		XA_DEBUG_ASSERT(index < m_size);
		return m_buffer[index];
	}

	uint32_t size() const { return m_size; }
	const T * data() const { return m_buffer; }
	T * data() { return m_buffer; }
	T * begin() { return m_buffer; }
	T * end() { return m_buffer + m_size; }
	const T * begin() const { return m_buffer; }
	const T * end() const { return m_buffer + m_size; }
	bool isEmpty() const { return m_size == 0; }

	void push_back( const T & val )
	{
		XA_DEBUG_ASSERT(&val < m_buffer || &val >= m_buffer+m_size);
		uint32_t old_size = m_size;
		uint32_t new_size = m_size + 1;
		setArraySize(new_size);
		construct_range(m_buffer, new_size, old_size, val);
	}

	void pop_back()
	{
		XA_DEBUG_ASSERT( m_size > 0 );
		resize( m_size - 1 );
	}

	const T & back() const
	{
		XA_DEBUG_ASSERT( m_size > 0 );
		return m_buffer[m_size-1];
	}

	T & back()
	{
		XA_DEBUG_ASSERT( m_size > 0 );
		return m_buffer[m_size-1];
	}

	const T & front() const
	{
		XA_DEBUG_ASSERT( m_size > 0 );
		return m_buffer[0];
	}

	T & front()
	{
		XA_DEBUG_ASSERT( m_size > 0 );
		return m_buffer[0];
	}

	// Remove the element at the given index. This is an expensive operation!
	void removeAt(uint32_t index)
	{
		XA_DEBUG_ASSERT(index >= 0 && index < m_size);
		if (m_size == 1) {
			clear();
		}
		else {
			m_buffer[index].~T();
			memmove(m_buffer+index, m_buffer+index+1, sizeof(T) * (m_size - 1 - index));
			m_size--;
		}
	}

	// Insert the given element at the given index shifting all the elements up.
	void insertAt(uint32_t index, const T & val = T())
	{
		XA_DEBUG_ASSERT( index >= 0 && index <= m_size );
		setArraySize(m_size + 1);
		if (index < m_size - 1) {
			memmove(m_buffer+index+1, m_buffer+index, sizeof(T) * (m_size - 1 - index));
		}
		// Copy-construct into the newly opened slot.
		new(m_buffer+index) T(val);
	}

	void append(const Array<T> & other)
	{
		append(other.m_buffer, other.m_size);
	}

	void resize(uint32_t new_size)
	{
		uint32_t old_size = m_size;
		// Destruct old elements (if we're shrinking).
		destroy_range(m_buffer, new_size, old_size);
		setArraySize(new_size);
		// Call default constructors
		construct_range(m_buffer, new_size, old_size);
	}

	void resize(uint32_t new_size, const T & elem)
	{
		XA_DEBUG_ASSERT(&elem < m_buffer || &elem > m_buffer+m_size);
		uint32_t old_size = m_size;
		// Destruct old elements (if we're shrinking).
		destroy_range(m_buffer, new_size, old_size);
		setArraySize(new_size);
		// Call copy constructors
		construct_range(m_buffer, new_size, old_size, elem);
	}

	void clear()
	{
		// Destruct old elements
		destroy_range(m_buffer, 0, m_size);
		m_size = 0;
	}

	void reserve(uint32_t desired_size)
	{
		if (desired_size > m_capacity) {
			setArrayCapacity(desired_size);
		}
	}

	void copy(const T * data, uint32_t count)
	{
		destroy_range(m_buffer, 0, m_size);
		setArraySize(count);
		construct_range(m_buffer, count, 0, data);
	}

	Array<T> & operator=( const Array<T> & a )
	{
		copy(a.m_buffer, a.m_size);
		return *this;
	}

	friend void swap(Array<T> & a, Array<T> & b)
	{
		std::swap(a.m_buffer, b.m_buffer);
		std::swap(a.m_capacity, b.m_capacity);
		std::swap(a.m_size, b.m_size);
	}

protected:
	void setArraySize(uint32_t new_size)
	{
		m_size = new_size;
		if (new_size > m_capacity) {
			uint32_t new_buffer_size;
			if (m_capacity == 0) {
				// first allocation is exact
				new_buffer_size = new_size;
			}
			else {
				// following allocations grow array by 25%
				new_buffer_size = new_size + (new_size >> 2);
			}
			setArrayCapacity( new_buffer_size );
		}
	}
	void setArrayCapacity(uint32_t new_capacity)
	{
		XA_DEBUG_ASSERT(new_capacity >= m_size);
		if (new_capacity == 0) {
			// free the buffer.
			if (m_buffer != NULL) {
				XA_FREE(m_buffer);
				m_buffer = NULL;
			}
		}
		else {
			// realloc the buffer
			m_buffer = XA_REALLOC(m_buffer, T, new_capacity);
		}
		m_capacity = new_capacity;
	}

	T * m_buffer;
	uint32_t m_capacity;
	uint32_t m_size;
};

/// Basis class to compute tangent space basis, ortogonalizations and to
/// transform vectors from one space to another.
class Basis
{
public:
	/// Create a null basis.
	Basis() : tangent(0, 0, 0), bitangent(0, 0, 0), normal(0, 0, 0) {}

	void buildFrameForDirection(Vector3::Arg d, float angle = 0)
	{
		XA_ASSERT(isNormalized(d));
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
	BitArray() : m_size(0) {}
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
		XA_DEBUG_ASSERT( b < m_size );
		return (m_wordArray[b >> 5] & (1 << (b & 31))) != 0;
	}

	// Set a bit.
	void setBitAt(uint32_t idx)
	{
		XA_DEBUG_ASSERT(idx < m_size);
		m_wordArray[idx >> 5] |=  (1 << (idx & 31));
	}

	// Toggle a bit.
	void toggleBitAt(uint32_t idx)
	{
		XA_DEBUG_ASSERT(idx < m_size);
		m_wordArray[idx >> 5] ^= (1 << (idx & 31));
	}

	// Set a bit to the given value. @@ Rename modifyBitAt?
	void setBitAt(uint32_t idx, bool b)
	{
		XA_DEBUG_ASSERT(idx < m_size);
		m_wordArray[idx >> 5] = setBits(m_wordArray[idx >> 5], 1 << (idx & 31), b);
		XA_DEBUG_ASSERT(bitAt(idx) == b);
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
	uint32_t setBits(uint32_t w, uint32_t m, bool b)
	{
		return (w & ~m) | (-int(b) & m);
	}

	// Number of bits stored.
	uint32_t m_size;

	// Array of bits.
	Array<uint32_t> m_wordArray;
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
		XA_DEBUG_ASSERT(x < m_width && y < m_height);
		return m_bitArray.bitAt(y * m_width + x);
	}

	void setBitAt(uint32_t x, uint32_t y)
	{
		XA_DEBUG_ASSERT(x < m_width && y < m_height);
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

class Fit
{
public:
	static Vector3 computeCentroid(int n, const Vector3 * points)
	{
		Vector3 centroid(0.0f);
		for (int i = 0; i < n; i++) {
			centroid += points[i];
		}
		centroid /= float(n);
		return centroid;
	}

	static Vector3 computeCovariance(int n, const Vector3 * points, float * covariance)
	{
		// compute the centroid
		Vector3 centroid = computeCentroid(n, points);
		// compute covariance matrix
		for (int i = 0; i < 6; i++) {
			covariance[i] = 0.0f;
		}
		for (int i = 0; i < n; i++) {
			Vector3 v = points[i] - centroid;
			covariance[0] += v.x * v.x;
			covariance[1] += v.x * v.y;
			covariance[2] += v.x * v.z;
			covariance[3] += v.y * v.y;
			covariance[4] += v.y * v.z;
			covariance[5] += v.z * v.z;
		}
		return centroid;
	}

	static bool isPlanar(int n, const Vector3 *points, float epsilon = XA_EPSILON)
	{
		// compute the centroid and covariance
		float matrix[6];
		computeCovariance(n, points, matrix);
		float eigenValues[3];
		Vector3 eigenVectors[3];
		if (!eigenSolveSymmetric3(matrix, eigenValues, eigenVectors)) {
			return false;
		}
		return eigenValues[2] < epsilon;
	}

	// Tridiagonal solver from Charles Bloom.
	// Householder transforms followed by QL decomposition.
	// Seems to be based on the code from Numerical Recipes in C.
	static bool eigenSolveSymmetric3(const float matrix[6], float eigenValues[3], Vector3 eigenVectors[3])
	{
		XA_DEBUG_ASSERT(matrix != NULL && eigenValues != NULL && eigenVectors != NULL);
		float subd[3];
		float diag[3];
		float work[3][3];
		work[0][0] = matrix[0];
		work[0][1] = work[1][0] = matrix[1];
		work[0][2] = work[2][0] = matrix[2];
		work[1][1] = matrix[3];
		work[1][2] = work[2][1] = matrix[4];
		work[2][2] = matrix[5];
		EigenSolver3_Tridiagonal(work, diag, subd);
		if (!EigenSolver3_QLAlgorithm(work, diag, subd)) {
			for (int i = 0; i < 3; i++) {
				eigenValues[i] = 0;
				eigenVectors[i] = Vector3(0);
			}
			return false;
		}
		for (int i = 0; i < 3; i++) {
			eigenValues[i] = (float)diag[i];
		}
		// eigenvectors are the columns; make them the rows :
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				(&eigenVectors[j].x)[i] = (float) work[i][j];
			}
		}
		// shuffle to sort by singular value :
		if (eigenValues[2] > eigenValues[0] && eigenValues[2] > eigenValues[1]) {
			std::swap(eigenValues[0], eigenValues[2]);
			std::swap(eigenVectors[0], eigenVectors[2]);
		}
		if (eigenValues[1] > eigenValues[0]) {
			std::swap(eigenValues[0], eigenValues[1]);
			std::swap(eigenVectors[0], eigenVectors[1]);
		}
		if (eigenValues[2] > eigenValues[1]) {
			std::swap(eigenValues[1], eigenValues[2]);
			std::swap(eigenVectors[1], eigenVectors[2]);
		}
		XA_DEBUG_ASSERT(eigenValues[0] >= eigenValues[1] && eigenValues[0] >= eigenValues[2]);
		XA_DEBUG_ASSERT(eigenValues[1] >= eigenValues[2]);
		return true;
	}

private:
	static void EigenSolver3_Tridiagonal(float mat[3][3], float *diag, float *subd)
	{
		// Householder reduction T = Q^t M Q
		//   Input:
		//     mat, symmetric 3x3 matrix M
		//   Output:
		//     mat, orthogonal matrix Q
		//     diag, diagonal entries of T
		//     subd, subdiagonal entries of T (T is symmetric)
		const float epsilon = 1e-08f;
		float a = mat[0][0];
		float b = mat[0][1];
		float c = mat[0][2];
		float d = mat[1][1];
		float e = mat[1][2];
		float f = mat[2][2];
		diag[0] = a;
		subd[2] = 0.f;
		if (fabsf(c) >= epsilon) {
			const float ell = sqrtf(b * b + c * c);
			b /= ell;
			c /= ell;
			const float q = 2 * b * e + c * (f - d);
			diag[1] = d + c * q;
			diag[2] = f - c * q;
			subd[0] = ell;
			subd[1] = e - b * q;
			mat[0][0] = 1;
			mat[0][1] = 0;
			mat[0][2] = 0;
			mat[1][0] = 0;
			mat[1][1] = b;
			mat[1][2] = c;
			mat[2][0] = 0;
			mat[2][1] = c;
			mat[2][2] = -b;
		} else {
			diag[1] = d;
			diag[2] = f;
			subd[0] = b;
			subd[1] = e;
			mat[0][0] = 1;
			mat[0][1] = 0;
			mat[0][2] = 0;
			mat[1][0] = 0;
			mat[1][1] = 1;
			mat[1][2] = 0;
			mat[2][0] = 0;
			mat[2][1] = 0;
			mat[2][2] = 1;
		}
	}

	static bool EigenSolver3_QLAlgorithm(float mat[3][3], float *diag, float *subd)
	{
		// QL iteration with implicit shifting to reduce matrix from tridiagonal
		// to diagonal
		const int maxiter = 32;
		for (int ell = 0; ell < 3; ell++) {
			int iter;
			for (iter = 0; iter < maxiter; iter++) {
				int m;
				for (m = ell; m <= 1; m++) {
					float dd = fabsf(diag[m]) + fabsf(diag[m + 1]);
					if ( fabsf(subd[m]) + dd == dd )
						break;
				}
				if ( m == ell )
					break;
				float g = (diag[ell + 1] - diag[ell]) / (2 * subd[ell]);
				float r = sqrtf(g * g + 1);
				if ( g < 0 )
					g = diag[m] - diag[ell] + subd[ell] / (g - r);
				else
					g = diag[m] - diag[ell] + subd[ell] / (g + r);
				float s = 1, c = 1, p = 0;
				for (int i = m - 1; i >= ell; i--) {
					float f = s * subd[i], b = c * subd[i];
					if ( fabsf(f) >= fabsf(g) ) {
						c = g / f;
						r = sqrtf(c * c + 1);
						subd[i + 1] = f * r;
						c *= (s = 1 / r);
					} else {
						s = f / g;
						r = sqrtf(s * s + 1);
						subd[i + 1] = g * r;
						s *= (c = 1 / r);
					}
					g = diag[i + 1] - p;
					r = (diag[i] - g) * s + 2 * b * c;
					p = s * r;
					diag[i + 1] = g + p;
					g = c * r - b;
					for (int k = 0; k < 3; k++) {
						f = mat[k][i + 1];
						mat[k][i + 1] = s * mat[k][i] + c * f;
						mat[k][i] = c * mat[k][i] - s * f;
					}
				}
				diag[ell] -= p;
				subd[ell] = g;
				subd[m] = 0;
			}
			if ( iter == maxiter )
				// should not get here under normal circumstances
				return false;
		}
		return true;
	}
};

/// Fixed size vector class.
class FullVector
{
public:
	FullVector(uint32_t dim) { m_array.resize(dim); }
	FullVector(const FullVector &v) : m_array(v.m_array) {}

	const FullVector &operator=(const FullVector &v)
	{
		XA_ASSERT(dimension() == v.dimension());
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

private:
	Array<float> m_array;
};

template<typename Key, typename Value, typename H = Hash<Key>, typename E = Equal<Key> >
class HashMap
{
public:
	HashMap()
	{
		alloc(4096);
	}

	explicit HashMap(uint32_t size)
	{
		alloc(size > 0 ? size : 4096);
	}

	~HashMap()
	{
		if (m_slots)
			XA_FREE(m_slots);
	}

	void add(const Key &key, const Value &value)
	{
		const uint32_t hash = computeHash(key);
		Element element;
		element.key = key;
		element.value = value;
		element.next = m_slots[hash];
		m_slots[hash] = m_elements.size();
		m_elements.push_back(element);
	}

	bool remove(const Key &key)
	{
		const uint32_t hash = computeHash(key);
		uint32_t i = m_slots[hash];
		Element *prevElement = NULL;
		E equal;
		while (i != UINT32_MAX) {
			Element *element = &m_elements[i];
			if (equal(element->key, key)) {
				if (prevElement)
					prevElement->next = element->next;
				else
					m_slots[hash] = element->next;
				// Don't remove from m_elements, that would mess up Element::next indices.
				return true;
			}
			prevElement = element;
			i = element->next;
		}
		return false;
	}

	struct Element
	{
		Key key;
		Value value;
		uint32_t next;
	};

	bool removeElement(const Element *element)
	{
		const uint32_t hash = computeHash(element->key);
		uint32_t i = m_slots[hash];
		Element *prevElement = NULL;
		while (i != UINT32_MAX) {
			Element *e = &m_elements[i];
			if (e == element) {
				if (prevElement)
					prevElement->next = e->next;
				else
					m_slots[hash] = e->next;
				// Don't remove from m_elements, that would mess up Element::next indices.
				return true;
			}
			prevElement = e;
			i = e->next;
		}
		return false;
	}

	const Element *get(const Key &key) const
	{
		const uint32_t hash = computeHash(key);
		uint32_t i = m_slots[hash];
		E equal;
		while (i != UINT32_MAX) {
			const Element *element = &m_elements[i];
			if (equal(element->key, key))
				return element;
			i = element->next;
		}
		return NULL;
	}

	const Element *getNext(const Element *current) const
	{
		uint32_t i = current->next;
		E equal;
		while (i != UINT32_MAX) {
			const Element *element = &m_elements[i];
			if (equal(element->key, current->key))
				return element;
			i = element->next;
		}
		return NULL;
	}

private:
	void alloc(uint32_t size)
	{
		m_numSlots = (uint32_t)(size * 1.3);
		m_slots = XA_ALLOC_ARRAY(uint32_t, m_numSlots);
		for (uint32_t i = 0; i < m_numSlots; i++)
			m_slots[i] = UINT32_MAX;
		m_elements.reserve(size);
	}

	uint32_t computeHash(const Key &key) const
	{
		H hash;
		return hash(key) % m_numSlots;
	}

	uint32_t m_numSlots;
	uint32_t *m_slots;
	internal::Array<Element> m_elements;
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
		if (max == UINT32_MAX) return get();
		const uint32_t np2 = nextPowerOfTwo( max + 1 ); // @@ This fails if max == UINT32_MAX
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

	uint32_t state[N];	// internal state
	uint32_t *next;	// next value to get from state
	int left;		// number of values left before reload needed
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
		XA_FREE(m_ranks2);
		XA_FREE(m_ranks);
	}

	RadixSort &sort(const float *input, uint32_t count)
	{
		if (input == NULL || count == 0) return *this;
		// Resize lists if needed
		if (count != m_size) {
			if (count > m_size) {
				m_ranks2 = XA_REALLOC(m_ranks2, uint32_t, count);
				m_ranks = XA_REALLOC(m_ranks, uint32_t, count);
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

	RadixSort &sort(const Array<float> &input)
	{
		return sort(input.data(), input.size());
	}

	// Access to results. m_ranks is a list of indices in sorted order, i.e. in the order you may further process your data
	const uint32_t *ranks() const
	{
		XA_DEBUG_ASSERT(m_validRanks);
		return m_ranks;
	}

	uint32_t *ranks()
	{
		XA_DEBUG_ASSERT(m_validRanks);
		return m_ranks;
	}

private:
	uint32_t m_size;
	uint32_t *m_ranks;
	uint32_t *m_ranks2;
	bool m_validRanks;

	void FloatFlip(uint32_t &f)
	{
		int32_t mask = (int32_t(f) >> 31) | 0x80000000; // Warren Hunt, Manchor Ko.
		f ^= mask;
	}

	void IFloatFlip(uint32_t &f)
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
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4127)
#endif
			if (bucketCount == 8) h[4][*p++]++, h[5][*p++]++, h[6][*p++]++, h[7][*p++]++;
#ifdef _MSC_VER
#pragma warning(pop)
#endif
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

struct Edge
{
	uint32_t relativeIndex; // absolute: face.firstIndex + relativeIndex
	uint32_t face;
	uint32_t index0;
	uint32_t index1;
};

struct Face
{
	uint32_t firstIndex; // Index into Mesh::m_indices.
	uint32_t nIndices;
};

struct FaceFlags
{
	enum
	{
		Ignore = 1<<0
	};
};

class Mesh
{
public:
	Mesh(uint32_t approxVertexCount = 0, uint32_t approxFaceCount = 0, uint32_t id = UINT32_MAX) : m_id(id), m_colocalVertexCount(0), m_edgeMap(approxFaceCount * 3), m_vertexToEdgeMap(approxVertexCount)
	{
		m_edges.reserve(approxFaceCount * 3);
		m_faces.reserve(approxFaceCount);
		m_faceFlags.reserve(approxFaceCount);
		m_faceGroups.reserve(approxFaceCount);
		m_indices.reserve(approxFaceCount * 3);
		m_positions.reserve(approxVertexCount);
		m_normals.reserve(approxVertexCount);
		m_texcoords.reserve(approxVertexCount);
	}

	uint32_t id() const { return m_id; }

	void addVertex(const Vector3 &pos, const Vector3 &normal = Vector3(0.0f), const Vector2 &texcoord = Vector2(0.0f))
	{
		XA_DEBUG_ASSERT(isFinite(pos));
		m_positions.push_back(pos);
		m_normals.push_back(normal);
		m_texcoords.push_back(texcoord);
	}

	void addFace(uint32_t v0, uint32_t v1, uint32_t v2, uint32_t flags = 0)
	{
		uint32_t indexArray[3];
		indexArray[0] = v0;
		indexArray[1] = v1;
		indexArray[2] = v2;
		addFace(indexArray, 3, flags);
	}

	void addFace(const Array<uint32_t> &indexArray, uint32_t flags = 0)
	{
		addFace(indexArray.data(), indexArray.size(), flags);
	}

	void addFace(const uint32_t *indexArray, uint32_t indexCount, uint32_t flags = 0)
	{
		Face face;
		face.firstIndex = m_indices.size();
		face.nIndices = indexCount;
		m_faces.push_back(face);
		m_faceFlags.push_back(flags);
		m_faceGroups.push_back(UINT32_MAX);
		for (uint32_t i = 0; i < indexCount; i++)
			m_indices.push_back(indexArray[i]);
		for (uint32_t i = 0; i < indexCount; i++) {
			Edge edge;
			edge.face = m_faces.size() - 1;
			edge.relativeIndex = i;
			edge.index0 = face.firstIndex + i;
			edge.index1 = face.firstIndex + (i + 1) % face.nIndices;
			m_edges.push_back(edge);
			const uint32_t edgeIndex = m_edges.size() - 1;
			const uint32_t vertex0 = m_indices[edge.index0];
			const uint32_t vertex1 = m_indices[edge.index1];
			m_vertexToEdgeMap.add(vertex0, edgeIndex);
			m_vertexToEdgeMap.add(vertex1, edgeIndex);
			EdgeKey key(vertex0, vertex1);
			//if (!m_edgeMap.get(key))
			m_edgeMap.add(key, edgeIndex);
		}
	}

	void createColocals()
	{
		XA_PRINT(PrintFlags::MeshCreation, "--- Linking colocals:\n");
		const uint32_t vertexCount = m_positions.size();
		typedef HashMap<Vector3, uint32_t, Vector3Hash, Vector3Equal> PositionHashMap;
		PositionHashMap positionHashMap(vertexCount);
		for (uint32_t i = 0; i < m_positions.size(); i++)
			positionHashMap.add(m_positions[i], i);
		Array<uint32_t> colocals;
		m_nextColocalVertex.resize(vertexCount, UINT32_MAX);
		for (uint32_t i = 0; i < vertexCount; i++) {
			if (m_nextColocalVertex[i] != UINT32_MAX)
				continue; // Already done.
			colocals.clear();
			const PositionHashMap::Element *ele = positionHashMap.get(m_positions[i]);
			while (ele) {
				colocals.push_back(ele->value);
				ele = positionHashMap.getNext(ele);
			}
			if (colocals.size() == 1) {
				// No colocals for this vertex.
				m_nextColocalVertex[i] = i;
				continue; 
			}
			std::sort(colocals.begin(), colocals.end());
			for (uint32_t j = 0; j < colocals.size(); j++)
				m_nextColocalVertex[colocals[j]] = colocals[(j + 1) % colocals.size()];
			XA_DEBUG_ASSERT(m_nextColocalVertex[i] != UINT32_MAX);
		}
	}

	// Check if the face duplicates any edges of any face already in the group.
	bool faceDuplicatesGroupEdge(uint32_t group, uint32_t face) const
	{
		for (FaceEdgeIterator edgeIt(this, face); !edgeIt.isDone(); edgeIt.advance()) {
			for (ColocalEdgeIterator colocalEdgeIt(this, edgeIt.vertex0(), edgeIt.vertex1()); !colocalEdgeIt.isDone(); colocalEdgeIt.advance()) {
				if (m_faceGroups[m_edges[colocalEdgeIt.edge()].face] == group)
					return true;
			}
		}
		return false;
	}

	// Check if the face mirrors any face already in the group.
	// i.e. don't want two-sided faces in the same group.
	// A face mirrors another face if all edges match with opposite winding.
	bool faceMirrorsGroupFace(uint32_t group, uint32_t face) const
	{
		FaceEdgeIterator edgeIt(this, face);
		for (ColocalEdgeIterator colocalEdgeIt(this, edgeIt.vertex1(), edgeIt.vertex0()); !colocalEdgeIt.isDone(); colocalEdgeIt.advance()) {
			const uint32_t candidateFace = m_edges[colocalEdgeIt.edge()].face;
			if (m_faceGroups[candidateFace] == group) {
				// Found a match for mirrored first edge, try the other edges.
				bool match = false;
				for (; !edgeIt.isDone(); edgeIt.advance()) {
					match = false;
					for (ColocalEdgeIterator colocalEdgeIt2(this, edgeIt.vertex1(), edgeIt.vertex0()); !colocalEdgeIt2.isDone(); colocalEdgeIt2.advance()) {
						if (m_edges[colocalEdgeIt2.edge()].face == candidateFace) {
							match = true;
							break;
						}
					}
					if (!match)
						break;
				}
				if (match)
					return true; // All edges are mirrored in this face.
				// Try the next face.
				edgeIt = FaceEdgeIterator(this, candidateFace);
			}
		}
		return false;
	}

	void createFaceGroups()
	{
		const uint32_t faceCount = m_faces.size();
		uint32_t group = 0;
		Array<uint32_t> growFaces;
		for (;;) {
			// Find an unassigned face.
			uint32_t face = UINT32_MAX;
			for (uint32_t f = 0; f < faceCount; f++) {
				if (m_faceGroups[f] == UINT32_MAX && !(m_faceFlags[f] & FaceFlags::Ignore)) {
					face = f;
					break;
				}
			}
			if (face == UINT32_MAX)
				break; // All faces assigned to a group (except ignored faces).
			m_faceGroups[face] = group;
			growFaces.clear();
			growFaces.push_back(face);
			// Find faces connected to the face and assign them to the same group as the face, unless they are already assigned to another group.
			for (;;) {
				if (growFaces.isEmpty())
					break;
				const uint32_t f = growFaces.back();
				growFaces.pop_back();
				for (FaceEdgeIterator edgeIt(this, f); !edgeIt.isDone(); edgeIt.advance()) {
					// Iterate opposite edges. There may be more than one - non-manifold geometry can have duplicate edges.
					// Prioritize the one with exact vertex match, not just colocal.
					// If *any* of the opposite edges are already assigned to this group, don't do anything.
					bool alreadyAssignedToThisGroup = false;
					uint32_t bestConnectedFace = UINT32_MAX;
					for (ColocalEdgeIterator oppositeEdgeIt(this, edgeIt.vertex1(), edgeIt.vertex0()); !oppositeEdgeIt.isDone(); oppositeEdgeIt.advance()) {
						const Edge &oppositeEdge = m_edges[oppositeEdgeIt.edge()];
						if (m_faceFlags[oppositeEdge.face] & FaceFlags::Ignore)
							continue; // Don't add ignored faces to group.
						if (m_faceGroups[oppositeEdge.face] == group) {
							alreadyAssignedToThisGroup = true;
							break;
						}
						if (m_faceGroups[oppositeEdge.face] != UINT32_MAX)
							continue; // Connected face is already assigned to another group.
						if (faceDuplicatesGroupEdge(group, oppositeEdge.face))
							continue; // Don't want duplicate edges in a group.
						if (faceMirrorsGroupFace(group, oppositeEdge.face))
							continue; // Don't want two-sided faces in a group.
						const uint32_t oppositeVertex0 = m_indices[oppositeEdge.index0];
						const uint32_t oppositeVertex1 = m_indices[oppositeEdge.index1];
						if (bestConnectedFace == UINT32_MAX || (oppositeVertex0 == edgeIt.vertex1() && oppositeVertex1 == edgeIt.vertex0()))
							bestConnectedFace = oppositeEdge.face;
					}
					if (!alreadyAssignedToThisGroup && bestConnectedFace != UINT32_MAX) {
						m_faceGroups[bestConnectedFace] = group;
						growFaces.push_back(bestConnectedFace);
					}
				}
			}
			group++;
		}
	}

	void createBoundaries()
	{
		XA_PRINT(PrintFlags::MeshProcessing, "--- Creating boundaries:\n");
		const uint32_t edgeCount = m_edges.size();
		const uint32_t faceCount = m_faces.size();
		const uint32_t vertexCount = m_positions.size();
		m_oppositeEdges.resize(edgeCount);
		m_boundaryVertices.resize(vertexCount);
		for (uint32_t i = 0; i < edgeCount; i++)
			m_oppositeEdges[i] = UINT32_MAX;
		for (uint32_t i = 0; i < vertexCount; i++)
			m_boundaryVertices[i] = false;
		uint32_t nBoundaryEdges = 0;
		for (uint32_t i = 0; i < faceCount; i++) {
			if (m_faceFlags[i] & FaceFlags::Ignore)
				continue;
			const Face &face = m_faces[i];
			for (uint32_t j = 0; j < face.nIndices; j++) {
				const uint32_t vertex0 = m_indices[face.firstIndex + j];
				const uint32_t vertex1 = m_indices[face.firstIndex + (j + 1) % face.nIndices];
				// If there is an edge with opposite winding to this one, the edge isn't on a boundary.
				const Edge *oppositeEdge = findEdge(m_faceGroups[i], vertex1, vertex0);
				if (oppositeEdge) {
					XA_DEBUG_ASSERT(m_faceGroups[oppositeEdge->face] == m_faceGroups[i]);
					XA_DEBUG_ASSERT(!(m_faceFlags[oppositeEdge->face] & FaceFlags::Ignore));
					m_oppositeEdges[face.firstIndex + j] = m_faces[oppositeEdge->face].firstIndex + oppositeEdge->relativeIndex;
				} else {
					m_boundaryVertices[vertex0] = m_boundaryVertices[vertex1] = true;
					nBoundaryEdges++;
				}
			}
		}
		XA_PRINT(PrintFlags::MeshProcessing, "---   %d boundary edges.\n", nBoundaryEdges);
	}

	void linkBoundaries()
	{
		XA_PRINT(PrintFlags::MeshProcessing, "---   Linking boundaries:\n");
		const uint32_t edgeCount = m_edges.size();
		m_nextBoundaryEdges.resize(edgeCount);
		for (uint32_t i = 0; i < edgeCount; i++)
			m_nextBoundaryEdges[i] = UINT32_MAX;
		uint32_t numBoundaryLoops = 0, numUnclosedBoundaries = 0;
		BitArray bitFlags(edgeCount);
		bitFlags.clearAll();
		for (;;) {
			uint32_t firstEdge = UINT32_MAX;
			for (uint32_t i = 0; i < edgeCount; i++) {
				if (m_oppositeEdges[i] == UINT32_MAX && !bitFlags.bitAt(i)) {
					firstEdge = i;
					break;
				}
			}
			if (firstEdge == UINT32_MAX)
				break;
			uint32_t currentEdge = firstEdge;
			for (;;) {
				const Edge &edge = m_edges[currentEdge];
				// Find the next boundary edge. The first vertex will be the same as (or colocal to) the current edge second vertex.
				const uint32_t startVertex = m_indices[edge.index1];
				uint32_t bestNextEdge = UINT32_MAX;
				for (ColocalVertexIterator it(this, startVertex); !it.isDone(); it.advance()) {
					const VertexToEdgeMap::Element *ele = m_vertexToEdgeMap.get(it.vertex());
					while (ele) {
						const Edge &otherEdge = m_edges[ele->value];
						if (m_oppositeEdges[ele->value] != UINT32_MAX)
							goto next; // Not a boundary edge.
						if (bitFlags.bitAt(ele->value))
							goto next; // Already linked.
						if (m_faceGroups[edge.face] != m_faceGroups[otherEdge.face])
							goto next; // Don't cross face groups.
						if (m_faceFlags[otherEdge.face] & FaceFlags::Ignore)
							goto next; // Face is ignored.
						if (m_indices[otherEdge.index0] != it.vertex())
							goto next; // Edge contains the vertex, but it's the wrong one.
						// First edge has the lowest priority, don't want to close the boundary loop prematurely.
						// Non-colocal vertex has the highest.
						if (bestNextEdge == UINT32_MAX || bestNextEdge == firstEdge || it.vertex() == startVertex)
							bestNextEdge = ele->value;
					next:
						ele = m_vertexToEdgeMap.getNext(ele);
					}
				}
				if (bestNextEdge == UINT32_MAX) {
					numUnclosedBoundaries++;
					break; // Can't find a next edge.
				}
				m_nextBoundaryEdges[currentEdge] = bestNextEdge;
				bitFlags.setBitAt(bestNextEdge);
				currentEdge = bestNextEdge;
				if (currentEdge == firstEdge) {
					numBoundaryLoops++;
					break; // Closed the boundary loop.
				}
			}
		}
		XA_PRINT(PrintFlags::MeshProcessing, "---   %d boundary loops.\n", numBoundaryLoops);
		XA_PRINT(PrintFlags::MeshProcessing, "---   %d unclosed boundaries.\n", numUnclosedBoundaries);
	}

	/// Find edge, test all colocals.
	const Edge *findEdge(uint32_t faceGroup, uint32_t vertex0, uint32_t vertex1) const
	{
		const Edge *result = NULL;
		if (m_nextColocalVertex.isEmpty()) {
			EdgeKey key(vertex0, vertex1);
			const EdgeMap::Element *ele = m_edgeMap.get(key);
			while (ele) {
				const Edge *edge = &m_edges[ele->value];
				// Don't find edges of ignored faces.
				if ((faceGroup == UINT32_MAX || m_faceGroups[edge->face] == faceGroup) && !(m_faceFlags[edge->face] & FaceFlags::Ignore)) {
					XA_DEBUG_ASSERT(!result); // duplicate edge
					result = edge;
#if NDEBUG
					return result;
#endif
				}
				ele = m_edgeMap.getNext(ele);
			}
		} else {
			for (ColocalVertexIterator it0(this, vertex0); !it0.isDone(); it0.advance()) {
				for (ColocalVertexIterator it1(this, vertex1); !it1.isDone(); it1.advance()) {
					EdgeKey key(it0.vertex(), it1.vertex());
					const EdgeMap::Element *ele = m_edgeMap.get(key);
					while (ele) {
						const Edge *edge = &m_edges[ele->value];
						// Don't find edges of ignored faces.
						if ((faceGroup == UINT32_MAX || m_faceGroups[edge->face] == faceGroup) && !(m_faceFlags[edge->face] & FaceFlags::Ignore)) {
							XA_DEBUG_ASSERT(!result); // duplicate edge
							result = edge;
#if NDEBUG
							return result;
#endif
						}
						ele = m_edgeMap.getNext(ele);
					}
				}
			}
		}
		return result;
	}

#if XA_DEBUG_EXPORT_OBJ
	void writeObjVertices(FILE *file) const
	{
		for (uint32_t i = 0; i < m_positions.size(); i++)
			fprintf(file, "v %g %g %g\n", m_positions[i].x, m_positions[i].y, m_positions[i].z);
		for (uint32_t i = 0; i < m_normals.size(); i++)
			fprintf(file, "vn %g %g %g\n", m_normals[i].x, m_normals[i].y, m_normals[i].z);
		for (uint32_t i = 0; i < m_texcoords.size(); i++)
			fprintf(file, "vt %g %g\n", m_texcoords[i].x, m_texcoords[i].y);
	}

	void writeObjFace(FILE *file, uint32_t face) const
	{
		const Face &f = m_faces[face];
		fprintf(file, "f ");
		for (uint32_t j = 0; j < f.nIndices; j++) {
			const uint32_t index = m_indices[f.firstIndex + j] + 1; // 1-indexed
			fprintf(file, "%d/%d/%d%c", index, index, index, j == f.nIndices - 1 ? '\n' : ' ');
		}
	}

	void writeObjBoundaryEges(FILE *file) const
	{
		fprintf(file, "o boundary_edges\n");
		for (uint32_t i = 0; i < m_edges.size(); i++) {
			if (m_oppositeEdges[i] != UINT32_MAX)
				continue;
			const Edge &edge = m_edges[i];
			fprintf(file, "l %d %d\n", m_indices[edge.index0] + 1, m_indices[edge.index1] + 1); // 1-indexed
		}
	}

	void writeObjLinkedBoundaries(FILE *file) const
	{
		if (m_nextBoundaryEdges.size() == 0)
			return; // Boundaries aren't linked.
		const uint32_t edgeCount = m_edges.size();
		BitArray bitFlags(edgeCount);
		bitFlags.clearAll();
		uint32_t boundary = 0;
		for (;;) {
			uint32_t firstEdge = UINT32_MAX;
			for (uint32_t i = 0; i < edgeCount; i++) {
				if (m_oppositeEdges[i] == UINT32_MAX && !bitFlags.bitAt(i)) {
					firstEdge = i;
					break;
				}
			}
			if (firstEdge == UINT32_MAX)
				break;
			uint32_t edge = firstEdge;
			fprintf(file, "o boundary_%0.4d\n", boundary);
			fprintf(file, "l");
			for (;;) {
				bitFlags.setBitAt(edge);
				const uint32_t vertex0 = m_indices[m_edges[edge].index0];
				const uint32_t vertex1 = m_indices[m_edges[edge].index1];
				fprintf(file, " %d", vertex0 + 1); // 1-indexed
				edge = m_nextBoundaryEdges[edge];
				if (edge == firstEdge || edge == UINT32_MAX) {
					fprintf(file, " %d\n", vertex1 + 1); // 1-indexed
					break;
				}

			}
			boundary++;
		}
	}

	void writeObj() const
	{
		char filename[256];
		sprintf(filename, "debug_mesh_%0.3u.obj", m_id);
		FILE *file = fopen(filename, "w");
		if (!file)
			return;
		writeObjVertices(file);
		// groups
		uint32_t numGroups = 0;
		for (uint32_t i = 0; i < m_faceGroups.size(); i++) {
			if (m_faceGroups[i] != UINT32_MAX)
				numGroups = std::max(numGroups, m_faceGroups[i] + 1);
		}
		for (uint32_t i = 0; i < numGroups; i++) {
			fprintf(file, "o group_%0.4d\n", i);
			fprintf(file, "s off\n");
			for (uint32_t f = 0; f < m_faceGroups.size(); f++) {
				if (m_faceGroups[f] == i)
					writeObjFace(file, f);
			}
		}
		fprintf(file, "o group_ignored\n");
		fprintf(file, "s off\n");
		for (uint32_t f = 0; f < m_faceGroups.size(); f++) {
			if (m_faceGroups[f] == UINT32_MAX)
				writeObjFace(file, f);
		}
		writeObjBoundaryEges(file);
		fclose(file);
	}

	void writeSimpleObj(const char *filename) const
	{
		FILE *file = fopen(filename, "w");
		if (!file)
			return;
		writeObjVertices(file);
		fprintf(file, "s off\n");
		fprintf(file, "o object\n");
		for (uint32_t i = 0; i < m_faces.size(); i++)
			writeObjFace(file, i);
		writeObjBoundaryEges(file);
		writeObjLinkedBoundaries(file);
		fclose(file);
	}
#endif

	float computeSurfaceArea() const
	{
		float area = 0;
		for (uint32_t f = 0; f < m_faces.size(); f++)
			area += faceArea(f);
		XA_DEBUG_ASSERT(area >= 0);
		return area;
	}

	float computeParametricArea() const
	{
		float area = 0;
		for (uint32_t f = 0; f < m_faces.size(); f++)
			area += faceParametricArea(f);
		XA_DEBUG_ASSERT(area >= 0);
		return fabsf(area);
	}

	uint32_t countTriangles() const
	{
		const uint32_t faceCount = m_faces.size();
		uint32_t triangleCount = 0;
		for (uint32_t f = 0; f < faceCount; f++) {
			const Face &face = m_faces[f];
			const uint32_t edgeCount = face.nIndices;
			XA_DEBUG_ASSERT(edgeCount > 2);
			triangleCount += edgeCount - 2;
		}
		return triangleCount;
	}

	float faceArea(uint32_t face) const
	{
		float area = 0;
		Vector3 firstPos(0.0f);
		for (FaceEdgeIterator  it(this, face); !it.isDone(); it.advance()) {
			if (it.relativeEdge() == 0)
				firstPos = it.position0();
			else
				area += length(cross(it.position0() - firstPos, it.position1() - firstPos));
		}
		return area * 0.5f;
	}

	Vector3 faceCentroid(uint32_t face) const
	{
		Vector3 sum(0.0f);
		uint32_t count = 0;
		for (FaceEdgeIterator  it(this, face); !it.isDone(); it.advance()) {
			sum += it.position0();
			count++;
		}
		return sum / float(count);
	}

	Vector3 faceNormal(uint32_t face) const
	{
		Vector3 n(0.0f);
		Vector3 p0(0.0f);
		for (FaceEdgeIterator  it(this, face); !it.isDone(); it.advance()) {
			if (it.relativeEdge() == 0) {
				p0 = it.position0();
			} else if (it.position1() != p0) {
				const Vector3 &p1 = it.position0();
				const Vector3 &p2 = it.position1();
				const Vector3 v10 = p1 - p0;
				const Vector3 v20 = p2 - p0;
				n += cross(v10, v20);
			}
		}
		return normalizeSafe(n, Vector3(0, 0, 1), 0.0f);
	}

	float faceParametricArea(uint32_t face) const
	{
		float area = 0;
		Vector2 firstTexcoord(0.0f);
		for (FaceEdgeIterator  it(this, face); !it.isDone(); it.advance()) {
			if (it.relativeEdge() == 0)
				firstTexcoord = it.texcoord0();
			else
				area += triangleArea(firstTexcoord, it.texcoord0(), it.texcoord1());
		}
		return area * 0.5f;
	}


	// Average of the edge midpoints weighted by the edge length.
	// I want a point inside the triangle, but closer to the cirumcenter.
	Vector3 triangleCenter(uint32_t face) const
	{
		const Vector3 &p0 = m_positions[m_indices[m_faces[face].firstIndex + 0]];
		const Vector3 &p1 = m_positions[m_indices[m_faces[face].firstIndex + 1]];
		const Vector3 &p2 = m_positions[m_indices[m_faces[face].firstIndex + 2]];
		const float l0 = length(p1 - p0);
		const float l1 = length(p2 - p1);
		const float l2 = length(p0 - p2);
		const Vector3 m0 = (p0 + p1) * l0 / (l0 + l1 + l2);
		const Vector3 m1 = (p1 + p2) * l1 / (l0 + l1 + l2);
		const Vector3 m2 = (p2 + p0) * l2 / (l0 + l1 + l2);
		return m0 + m1 + m2;
	}

	// Unnormalized face normal assuming it's a triangle.
	Vector3 triangleNormal(uint32_t face) const
	{
		return normalizeSafe(triangleNormalAreaScaled(face), Vector3(0), 0.0f);
	}

	Vector3 triangleNormalAreaScaled(uint32_t face) const
	{
		const Vector3 &p0 = m_positions[m_indices[m_faces[face].firstIndex + 0]];
		const Vector3 &p1 = m_positions[m_indices[m_faces[face].firstIndex + 1]];
		const Vector3 &p2 = m_positions[m_indices[m_faces[face].firstIndex + 2]];
		const Vector3 e0 = p2 - p0;
		const Vector3 e1 = p1 - p0;
		return cross(e0, e1);
	}

	// @@ This is not exactly accurate, we should compare the texture coordinates...
	bool isSeam(uint32_t edge) const
	{
		const uint32_t oppositeEdge = m_oppositeEdges[edge];
		if (oppositeEdge == UINT32_MAX)
			return true; // boundary edge
		const Edge &e = m_edges[edge];
		const Edge &oe = m_edges[oppositeEdge];
		return m_indices[e.index0] != m_indices[oe.index1] || m_indices[e.index1] != m_indices[oe.index0];
	}

	bool isNormalSeam(uint32_t edge) const
	{
		const uint32_t oppositeEdge = m_oppositeEdges[edge];
		if (oppositeEdge == UINT32_MAX)
			return true; // boundary edge
		const Edge &e = m_edges[edge];
		const Edge &oe = m_edges[oppositeEdge];
		return m_normals[m_indices[e.index0]] != m_normals[m_indices[oe.index1]] || m_normals[m_indices[e.index1]] != m_normals[m_indices[oe.index0]];
	}

	bool isTextureSeam(uint32_t edge) const
	{
		const uint32_t oppositeEdge = m_oppositeEdges[edge];
		if (oppositeEdge == UINT32_MAX)
			return true; // boundary edge
		const Edge &e = m_edges[edge];
		const Edge &oe = m_edges[oppositeEdge];
		return m_texcoords[m_indices[e.index0]] != m_texcoords[m_indices[oe.index1]] || m_texcoords[m_indices[e.index1]] != m_texcoords[m_indices[oe.index0]];
	}

	uint32_t firstColocal(uint32_t vertex) const
	{
		for (ColocalVertexIterator it(this, vertex); !it.isDone(); it.advance()) {
			if (it.vertex() < vertex)
				vertex = it.vertex();
		}
		return vertex;
	}

	bool areColocal(uint32_t vertex0, uint32_t vertex1) const
	{
		if (vertex0 == vertex1)
			return true;
		if (m_nextColocalVertex.isEmpty())
			return false;
		for (ColocalVertexIterator it(this, vertex0); !it.isDone(); it.advance()) {
			if (it.vertex() == vertex1)
				return true;
		}
		return false;
	}

	uint32_t edgeCount() const { return m_edges.size(); }
	const Edge *edgeAt(uint32_t edge) const { return &m_edges[edge]; }
	uint32_t oppositeEdge(uint32_t edge) const { return m_oppositeEdges[edge]; }
	bool isBoundaryEdge(uint32_t edge) const { return m_oppositeEdges[edge] == UINT32_MAX; }
	bool isBoundaryVertex(uint32_t vertex) const { return m_boundaryVertices[vertex]; }
	uint32_t colocalVertexCount() const { return m_colocalVertexCount; }
	uint32_t vertexCount() const { return m_positions.size(); }
	uint32_t vertexAt(uint32_t i) const { return m_indices[i]; }
	const Vector3 *positionAt(uint32_t i) const { return &m_positions[i]; }
	Vector3 *positionAt(uint32_t i) { return &m_positions[i]; }
	const Vector3 *normalAt(uint32_t i) const { return &m_normals[i]; }
	Vector3 *normalAt(uint32_t i) { return &m_normals[i]; }
	const Vector2 *texcoordAt(uint32_t i) const { return &m_texcoords[i]; }
	Vector2 *texcoordAt(uint32_t i) { return &m_texcoords[i]; }
	uint32_t faceCount() const { return m_faces.size(); }
	const Face *faceAt(uint32_t i) const { return &m_faces[i]; }
	Face *faceAt(uint32_t i) { return &m_faces[i]; }
	uint32_t faceFlagsAt(uint32_t i) const { return m_faceFlags[i]; }
	uint32_t faceGroupAt(uint32_t face) const { return m_faceGroups[face]; }

private:
	uint32_t m_id;
	Array<Edge> m_edges;
	Array<Face> m_faces;
	Array<uint32_t> m_faceFlags;
	Array<uint32_t> m_faceGroups;
	Array<uint32_t> m_indices;
	Array<Vector3> m_positions;
	Array<Vector3> m_normals;
	Array<Vector2> m_texcoords;

	// Populated by createColocals
	uint32_t m_colocalVertexCount;
	Array<uint32_t> m_nextColocalVertex; // In: vertex index. Out: the vertex index of the next colocal position.

	// Populated by createBoundaries
	Array<uint32_t> m_nextBoundaryEdges; // The index of the next boundary edge. UINT32_MAX if the edge is not a boundary edge.
	Array<bool> m_boundaryVertices;
	Array<uint32_t> m_oppositeEdges; // In: edge index. Out: the index of the opposite edge (i.e. wound the opposite direction). UINT32_MAX if the input edge is a boundary edge.

	struct EdgeKey
	{
		EdgeKey() {}
		EdgeKey(const EdgeKey &k) : v0(k.v0), v1(k.v1) {}
		EdgeKey(uint32_t v0, uint32_t v1) : v0(v0), v1(v1) {}

		void operator=(const EdgeKey &k)
		{
			v0 = k.v0;
			v1 = k.v1;
		}
		bool operator==(const EdgeKey &k) const
		{
			return v0 == k.v0 && v1 == k.v1;
		}

		uint32_t v0;
		uint32_t v1;
	};

	typedef HashMap<EdgeKey, uint32_t> EdgeMap;
	EdgeMap m_edgeMap;
	typedef HashMap<uint32_t, uint32_t> VertexToEdgeMap;
	VertexToEdgeMap m_vertexToEdgeMap;

public:
	class BoundaryEdgeIterator
	{
	public:
		BoundaryEdgeIterator(const Mesh *mesh, uint32_t edge) : m_mesh(mesh), m_first(UINT32_MAX), m_current(edge) {}

		void advance()
		{
			if (m_first == UINT32_MAX)
				m_first = m_current;
			m_current = m_mesh->m_nextBoundaryEdges[m_current];
		}

		bool isDone() const
		{
			return m_first == m_current || m_current == UINT32_MAX;
		}

		uint32_t edge() const
		{
			return m_current;
		}

		uint32_t nextEdge() const
		{
			return m_mesh->m_nextBoundaryEdges[m_current];
		}

	private:
		const Mesh *m_mesh;
		uint32_t m_first;
		uint32_t m_current;
	};

	class ColocalVertexIterator
	{
	public:
		ColocalVertexIterator(const Mesh *mesh, uint32_t v) : m_mesh(mesh), m_first(UINT32_MAX), m_current(v) {}

		void advance()
		{
			if (m_first == UINT32_MAX)
				m_first = m_current;
			if (!m_mesh->m_nextColocalVertex.isEmpty())
				m_current = m_mesh->m_nextColocalVertex[m_current];
		}

		bool isDone() const
		{
			return m_first == m_current;
		}

		uint32_t vertex() const
		{
			return m_current;
		}

		const Vector3 *pos() const
		{
			return &m_mesh->m_positions[m_current];
		}

	private:
		const Mesh *m_mesh;
		uint32_t m_first;
		uint32_t m_current;
	};

	class ColocalEdgeIterator
	{
	public:
		ColocalEdgeIterator(const Mesh *mesh, uint32_t vertex0, uint32_t vertex1) : m_mesh(mesh), m_vertex0It(mesh, vertex0), m_vertex1It(mesh, vertex1), m_vertex1(vertex1)
		{
			resetElement();
		}

		void advance()
		{
			advanceElement();
		}

		bool isDone() const
		{
			return m_vertex0It.isDone() && m_vertex1It.isDone() && !m_currentElement;
		}

		uint32_t edge() const
		{
			return m_currentElement->value;
		}

	private:
		void resetElement()
		{
			m_currentElement = m_mesh->m_edgeMap.get(Mesh::EdgeKey(m_vertex0It.vertex(), m_vertex1It.vertex()));
			for (;;) {
				if (!m_currentElement)
					break;
				if (!isIgnoredFace())
					break;
				m_currentElement = m_mesh->m_edgeMap.getNext(m_currentElement);
			}
			if (!m_currentElement)
				advanceVertex1();
		}

		void advanceElement()
		{
			for (;;) {
				m_currentElement = m_mesh->m_edgeMap.getNext(m_currentElement);
				if (!m_currentElement)
					break;
				if (!isIgnoredFace())
					break;
			}
			if (!m_currentElement)
				advanceVertex1();
		}

		void advanceVertex0()
		{
			m_vertex0It.advance();
			if (m_vertex0It.isDone())
				return;
			m_vertex1It = ColocalVertexIterator(m_mesh, m_vertex1);
			resetElement();
		}

		void advanceVertex1()
		{
			m_vertex1It.advance();
			if (m_vertex1It.isDone())
				advanceVertex0();
			else
				resetElement();
		}

		bool isIgnoredFace() const
		{
			const Edge *edge = &m_mesh->m_edges[m_currentElement->value];
			return (m_mesh->m_faceFlags[edge->face] & FaceFlags::Ignore) != 0;
		}

		const Mesh *m_mesh;
		ColocalVertexIterator m_vertex0It, m_vertex1It;
		const uint32_t m_vertex1;
		const Mesh::EdgeMap::Element *m_currentElement;
	};

	class FaceEdgeIterator 
	{
	public:
		FaceEdgeIterator (const Mesh *mesh, uint32_t face) : m_mesh(mesh), m_face(face), m_relativeEdge(0)
		{
			m_edge = m_mesh->m_faces[m_face].firstIndex;
		}

		void advance()
		{
			if (m_relativeEdge < m_mesh->m_faces[m_face].nIndices) {
				m_edge++;
				m_relativeEdge++;
			}
		}

		bool isDone() const
		{
			return m_relativeEdge == m_mesh->m_faces[m_face].nIndices;
		}

		bool isBoundary() const { return m_mesh->m_oppositeEdges[m_edge] == UINT32_MAX; }
		bool isSeam() const { return m_mesh->isSeam(m_edge); }
		bool isNormalSeam() const { return m_mesh->isNormalSeam(m_edge); }
		bool isTextureSeam() const { return m_mesh->isTextureSeam(m_edge); }
		uint32_t edge() const { return m_edge; }
		uint32_t relativeEdge() const { return m_relativeEdge; }
		uint32_t face() const { return m_face; }
		uint32_t oppositeEdge() const { return m_mesh->m_oppositeEdges[m_edge]; }
		
		uint32_t oppositeFace() const
		{
			const uint32_t oedge = m_mesh->m_oppositeEdges[m_edge];
			if (oedge == UINT32_MAX)
				return UINT32_MAX;
			return m_mesh->m_edges[oedge].face;
		}

		uint32_t vertex0() const
		{
			const Face &face = m_mesh->m_faces[m_face];
			return m_mesh->m_indices[face.firstIndex + m_relativeEdge];
		}

		uint32_t vertex1() const
		{
			const Face &face = m_mesh->m_faces[m_face];
			return m_mesh->m_indices[face.firstIndex + (m_relativeEdge + 1) % face.nIndices];
		}

		const Vector3 &position0() const { return m_mesh->m_positions[vertex0()]; }
		const Vector3 &position1() const { return m_mesh->m_positions[vertex1()]; }
		const Vector3 &normal0() const { return m_mesh->m_normals[vertex0()]; }
		const Vector3 &normal1() const { return m_mesh->m_normals[vertex1()]; }
		const Vector2 &texcoord0() const { return m_mesh->m_texcoords[vertex0()]; }
		const Vector2 &texcoord1() const { return m_mesh->m_texcoords[vertex1()]; }

	private:
		const Mesh *m_mesh;
		uint32_t m_face;
		uint32_t m_edge;
		uint32_t m_relativeEdge;
	};
};

/*
Fixing T-junctions.

- Find T-junctions. Find  vertices that are on an edge.
- This test is approximate.
- Insert edges on a spatial index to speedup queries.
- Consider only open edges, that is edges that have no pairs.
- Consider only vertices on boundaries.
- Close T-junction.
- Split edge.

*/
struct SplitEdge
{
	uint32_t vertex;
	uint32_t edge;
	float t;
};

static Mesh *meshSplitBoundaryEdges(const Mesh &inputMesh) // Returns NULL if no split was made.
{
	XA_PRINT(PrintFlags::MeshProcessing, "Fixing T-junctions:\n");
	Array<SplitEdge> splitEdges;
	const uint32_t vertexCount = inputMesh.vertexCount();
	const uint32_t edgeCount = inputMesh.edgeCount();
	for (uint32_t v = 0; v < vertexCount; v++) {
		if (!inputMesh.isBoundaryVertex(v))
			continue;
		// Find edges that this vertex overlaps with.
		const Vector3 x0 = *inputMesh.positionAt(v);
		for (uint32_t e = 0; e < edgeCount; e++) {
			if (!inputMesh.isBoundaryEdge(e))
				continue;
			const Edge *edge = inputMesh.edgeAt(e);
			const Vector3 x1 = *inputMesh.positionAt(inputMesh.vertexAt(edge->index0));
			const Vector3 x2 = *inputMesh.positionAt(inputMesh.vertexAt(edge->index1));
			if (x1 == x0 || x2 == x0)
				continue; // Vertex lies on either edge vertex.
			const Vector3 v01 = x0 - x1;
			const Vector3 v21 = x2 - x1;
			const float l = length(v21);
			const float d = length(cross(v01, v21)) / l;
			if (!isZero(d))
				continue;
			float t = dot(v01, v21) / (l * l);
			if (t < XA_EPSILON || t > 1.0f - XA_EPSILON)
				continue;
			//XA_DEBUG_ASSERT(lerp(x1, x2, t) == x0);
			SplitEdge splitEdge;
			splitEdge.vertex = v;
			splitEdge.edge = e;
			splitEdge.t = t;
			splitEdges.push_back(splitEdge);
		}
	}
	if (splitEdges.isEmpty())
		return NULL;
	const uint32_t faceCount = inputMesh.faceCount();
	Mesh *mesh = XA_NEW(Mesh, vertexCount + splitEdges.size(), faceCount);
	for (uint32_t v = 0; v < vertexCount; v++)
		mesh->addVertex(*inputMesh.positionAt(v), *inputMesh.normalAt(v), *inputMesh.texcoordAt(v));
	for (uint32_t se = 0; se < splitEdges.size(); se++) {
		const SplitEdge &splitEdge = splitEdges[se];
		const Edge *edge = inputMesh.edgeAt(splitEdge.edge);
		Vector3 normal = lerp(*inputMesh.normalAt(inputMesh.vertexAt(edge->index0)), *inputMesh.normalAt(inputMesh.vertexAt(edge->index1)), splitEdge.t);
		Vector2 texcoord = lerp(*inputMesh.texcoordAt(inputMesh.vertexAt(edge->index0)), *inputMesh.texcoordAt(inputMesh.vertexAt(edge->index1)), splitEdge.t);
		mesh->addVertex(*inputMesh.positionAt(splitEdge.vertex), normal, texcoord);
	}
	Array<uint32_t> indexArray;
	for (uint32_t f = 0; f < faceCount; f++) {
		indexArray.clear();
		for (Mesh::FaceEdgeIterator it(&inputMesh, f); !it.isDone(); it.advance()) {
			indexArray.push_back(it.vertex0());
			for (uint32_t se = 0; se < splitEdges.size(); se++) {
				const SplitEdge &splitEdge = splitEdges[se];
				if (splitEdge.edge == it.edge()) {
					indexArray.push_back(vertexCount + se);
					break;
				}
			}
		}
		mesh->addFace(indexArray, inputMesh.faceFlagsAt(f));
	}
	XA_PRINT(PrintFlags::MeshProcessing, " - %d edges split.\n", splitEdges.size());
	mesh->createColocals(); // Added new vertices, some may be colocal with existing vertices.
	return mesh;
}

// This is doing a simple ear-clipping algorithm that skips invalid triangles. Ideally, we should
// also sort the ears by angle, start with the ones that have the smallest angle and proceed in order.
static Mesh *meshTriangulate(const Mesh &inputMesh)
{
	if (inputMesh.faceCount() * 3 == inputMesh.edgeCount())
		return NULL;
	const uint32_t vertexCount = inputMesh.vertexCount();
	const uint32_t faceCount = inputMesh.faceCount();
	Mesh *mesh = XA_NEW(Mesh, vertexCount, faceCount);
	// Add all vertices.
	for (uint32_t v = 0; v < vertexCount; v++)
		mesh->addVertex(*inputMesh.positionAt(v), *inputMesh.normalAt(v), *inputMesh.texcoordAt(v));
	Array<uint32_t> polygonVertices;
	Array<float> polygonAngles;
	Array<Vector2> polygonPoints;
	for (uint32_t f = 0; f < faceCount; f++) {
		const Face *face = inputMesh.faceAt(f);
		const uint32_t edgeCount = face->nIndices;
		XA_DEBUG_ASSERT(edgeCount >= 3);
		polygonVertices.clear();
		polygonVertices.reserve(edgeCount);
		if (edgeCount == 3) {
			// Simple case for triangles.
			for (Mesh::FaceEdgeIterator it(&inputMesh, f); !it.isDone(); it.advance())
				polygonVertices.push_back(it.vertex0());
			mesh->addFace(polygonVertices[0], polygonVertices[1], polygonVertices[2]);
		} else {
			// Build 2D polygon projecting vertices onto normal plane.
			// Faces are not necesarily planar, this is for example the case, when the face comes from filling a hole. In such cases
			// it's much better to use the best fit plane.
			const Vector3 fn = inputMesh.faceNormal(f);
			Basis basis;
			basis.buildFrameForDirection(fn);
			polygonPoints.clear();
			polygonPoints.reserve(edgeCount);
			polygonAngles.clear();
			polygonAngles.reserve(edgeCount);
			for (Mesh::FaceEdgeIterator it(&inputMesh, f); !it.isDone(); it.advance()) {
				polygonVertices.push_back(it.vertex0());
				const Vector3 &pos = it.position0();
				polygonPoints.push_back(Vector2(dot(basis.tangent, pos), dot(basis.bitangent, pos)));
			}
			polygonAngles.resize(edgeCount);
			while (polygonVertices.size() > 2) {
				const uint32_t size = polygonVertices.size();
				// Update polygon angles. @@ Update only those that have changed.
				float minAngle = float(2.0f * M_PI);
				uint32_t bestEar = 0; // Use first one if none of them is valid.
				bool bestIsValid = false;
				for (uint32_t i = 0; i < size; i++) {
					uint32_t i0 = i;
					uint32_t i1 = (i + 1) % size; // Use Sean's polygon interation trick.
					uint32_t i2 = (i + 2) % size;
					Vector2 p0 = polygonPoints[i0];
					Vector2 p1 = polygonPoints[i1];
					Vector2 p2 = polygonPoints[i2];
					float d = clamp(dot(p0 - p1, p2 - p1) / (length(p0 - p1) * length(p2 - p1)), -1.0f, 1.0f);
					float angle = acosf(d);
					float area = triangleArea(p0, p1, p2);
					if (area < 0.0f)
						angle = float(2.0f * M_PI - angle);
					polygonAngles[i1] = angle;
					if (angle < minAngle || !bestIsValid) {
						// Make sure this is a valid ear, if not, skip this point.
						bool valid = true;
						for (uint32_t j = 0; j < size; j++) {
							if (j == i0 || j == i1 || j == i2)
								continue;
							Vector2 p = polygonPoints[j];
							if (pointInTriangle(p, p0, p1, p2)) {
								valid = false;
								break;
							}
						}
						if (valid || !bestIsValid) {
							minAngle = angle;
							bestEar = i1;
							bestIsValid = valid;
						}
					}
				}
				// Clip best ear:
				const uint32_t i0 = (bestEar + size - 1) % size;
				const uint32_t i1 = (bestEar + 0) % size;
				const uint32_t i2 = (bestEar + 1) % size;
				mesh->addFace(polygonVertices[i0], polygonVertices[i1], polygonVertices[i2]);
				polygonVertices.removeAt(i1);
				polygonPoints.removeAt(i1);
				polygonAngles.removeAt(i1);
			}
		}
	}
	mesh->createBoundaries();
	mesh->linkBoundaries();
	return mesh;
}

static Mesh *meshUnifyVertices(const Mesh &inputMesh)
{
	const uint32_t vertexCount = inputMesh.vertexCount();
	const uint32_t faceCount = inputMesh.faceCount();
	Mesh *mesh = XA_NEW(Mesh, vertexCount, faceCount);
	// Only add the first colocal.
	for (uint32_t v = 0; v < vertexCount; v++) {
		if (inputMesh.firstColocal(v) == v)
			mesh->addVertex(*inputMesh.positionAt(v), *inputMesh.normalAt(v), *inputMesh.texcoordAt(v));
	}
	Array<uint32_t> indexArray;
	// Add new faces pointing to first colocals.
	for (uint32_t f = 0; f < faceCount; f++) {
		indexArray.clear();
		for (Mesh::FaceEdgeIterator it(&inputMesh, f); !it.isDone(); it.advance())
			indexArray.push_back(inputMesh.firstColocal(it.vertex0()));
		mesh->addFace(indexArray, inputMesh.faceFlagsAt(f));
	}
	mesh->createBoundaries();
	mesh->linkBoundaries();
	return mesh;
}

// boundaryEdges are the first edges for each boundary loop.
static void meshGetBoundaryEdges(const Mesh &mesh, Array<uint32_t> &boundaryEdges)
{
	const uint32_t edgeCount = mesh.edgeCount();
	BitArray bitFlags(edgeCount);
	bitFlags.clearAll();
	boundaryEdges.clear();
	// Search for boundary edges. Mark all the edges that belong to the same boundary.
	for (uint32_t e = 0; e < edgeCount; e++) {
		if (bitFlags.bitAt(e) || !mesh.isBoundaryEdge(e))
			continue;
		for (Mesh::BoundaryEdgeIterator it(&mesh, e); !it.isDone(); it.advance())
			bitFlags.setBitAt(it.edge());
		boundaryEdges.push_back(e);
	}
}

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
		XA_PRINT(PrintFlags::ComputingCharts, "--- Building mesh topology:\n" );
		Array<uint32_t> stack(faceCount);
		BitArray bitFlags(faceCount);
		bitFlags.clearAll();
		// Compute connectivity.
		XA_PRINT(PrintFlags::ComputingCharts, "---   Computing connectivity.\n" );
		m_connectedCount = 0;
		for (uint32_t f = 0; f < faceCount; f++ ) {
			if (bitFlags.bitAt(f) == false) {
				m_connectedCount++;
				stack.push_back(f);
				while (!stack.isEmpty()) {
					const uint32_t top = stack.back();
					XA_ASSERT(top != uint32_t(~0));
					stack.pop_back();
					if (bitFlags.bitAt(top) == false) {
						bitFlags.setBitAt(top);
						for (Mesh::FaceEdgeIterator it(mesh, top); !it.isDone(); it.advance()) {
							const uint32_t oppositeFace = it.oppositeFace();
							if (oppositeFace != UINT32_MAX)
								stack.push_back(oppositeFace);
						}
					}
				}
			}
		}
		XA_ASSERT(stack.isEmpty());
		XA_PRINT(PrintFlags::ComputingCharts, "---   %d connected components.\n", m_connectedCount);
		// Count boundary loops.
		XA_PRINT(PrintFlags::ComputingCharts, "---   Counting boundary loops.\n" );
		m_boundaryCount = 0;
		bitFlags.resize(edgeCount);
		bitFlags.clearAll();
		// Don't forget to link the boundary otherwise this won't work.
		for (uint32_t e = 0; e < edgeCount; e++) {
			if (bitFlags.bitAt(e) || !mesh->isBoundaryEdge(e))
				continue;
			m_boundaryCount++;
			for (Mesh::BoundaryEdgeIterator it(mesh, e); !it.isDone(); it.advance())
				bitFlags.setBitAt(it.edge());
		}
		XA_PRINT(PrintFlags::ComputingCharts, "---   %d boundary loops found.\n", m_boundaryCount);
		// Compute euler number.
		m_eulerNumber = vertexCount - edgeCount + faceCount;
		XA_PRINT(PrintFlags::ComputingCharts, "---   Euler number: %d.\n", m_eulerNumber);
		// Compute genus. (only valid on closed connected surfaces)
		m_genus = -1;
		if (isClosed() && isConnected()) {
			m_genus = (2 - m_eulerNumber) / 2;
			XA_PRINT(PrintFlags::ComputingCharts, "---   Genus: %d.\n", m_genus);
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

	void clipHorizontalPlane(float offset, float clipdirection)
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
	}

	void clipVerticalPlane(float offset, float clipdirection )
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
		XA_ASSERT(minx >= 0);
		XA_ASSERT(miny >= 0);
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
								Vector3 tex2 = t1 + dx * (x - v1.x) + dy * (y - v1.y);
								if (!cb(param, (int)x, (int)y, tex2, dx, dy, 1.0f)) {
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
									//XA_ASSERT(texCent.x >= -0.1f && texCent.x <= 1.1f); // @@ Centroid is not very exact...
									//XA_ASSERT(texCent.y >= -0.1f && texCent.y <= 1.1f);
									//XA_ASSERT(texCent.z >= -0.1f && texCent.z <= 1.1f);
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
static bool drawTriangle(Mode mode, Vector2::Arg extents, bool enableScissors, const Vector2 v[3], SamplingCallback cb, void *param)
{
	Triangle tri(v[0], v[1], v[2], Vector3(1, 0, 0), Vector3(0, 1, 0), Vector3(0, 0, 1));
	// @@ It would be nice to have a conservative drawing mode that enlarges the triangle extents by one texel and is able to handle degenerate triangles.
	// @@ Maybe the simplest thing to do would be raster triangle edges.
	if (tri.valid) {
		if (mode == Mode_Antialiased) {
			return tri.drawAA(extents, enableScissors, cb, param);
		}
		if (mode == Mode_Nearest) {
			return tri.draw(extents, enableScissors, cb, param);
		}
	}
	return true;
}

} // namespace raster

// Full and sparse vector and matrix classes. BLAS subset.
// Pseudo-BLAS interface.
namespace sparse {
enum Transpose
{
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
		XA_ASSERT(width() == m.width());
		XA_ASSERT(height() == m.height());
		m_array = m.m_array;
		return *this;
	}

	uint32_t width() const { return m_width; }
	uint32_t height() const { return m_array.size(); }
	bool isSquare() const { return width() == height(); }

	// x is column, y is row
	float getCoefficient(uint32_t x, uint32_t y) const
	{
		XA_DEBUG_ASSERT( x < width() );
		XA_DEBUG_ASSERT( y < height() );
		const uint32_t count = m_array[y].size();
		for (uint32_t i = 0; i < count; i++) {
			if (m_array[y][i].x == x) return m_array[y][i].v;
		}
		return 0.0f;
	}

	void setCoefficient(uint32_t x, uint32_t y, float f)
	{
		XA_DEBUG_ASSERT( x < width() );
		XA_DEBUG_ASSERT( y < height() );
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

	float dotRow(uint32_t y, const FullVector &v) const
	{
		XA_DEBUG_ASSERT( y < height() );
		const uint32_t count = m_array[y].size();
		float sum = 0;
		for (uint32_t i = 0; i < count; i++) {
			sum += m_array[y][i].v * v[m_array[y][i].x];
		}
		return sum;
	}

	void madRow(uint32_t y, float alpha, FullVector &v) const
	{
		XA_DEBUG_ASSERT(y < height());
		const uint32_t count = m_array[y].size();
		for (uint32_t i = 0; i < count; i++) {
			v[m_array[y][i].x] += alpha * m_array[y][i].v;
		}
	}

	void clearRow(uint32_t y)
	{
		XA_DEBUG_ASSERT( y < height() );
		m_array[y].clear();
	}

	const Array<Coefficient> &getRow(uint32_t y) const { return m_array[y]; }

private:
	/// Number of columns.
	const uint32_t m_width;

	/// Array of matrix elements.
	Array< Array<Coefficient> > m_array;
};

// y = a * x + y
static void saxpy(float a, const FullVector &x, FullVector &y)
{
	XA_DEBUG_ASSERT(x.dimension() == y.dimension());
	const uint32_t dim = x.dimension();
	for (uint32_t i = 0; i < dim; i++) {
		y[i] += a * x[i];
	}
}

static void copy(const FullVector &x, FullVector &y)
{
	XA_DEBUG_ASSERT(x.dimension() == y.dimension());
	const uint32_t dim = x.dimension();
	for (uint32_t i = 0; i < dim; i++) {
		y[i] = x[i];
	}
}

static void scal(float a, FullVector &x)
{
	const uint32_t dim = x.dimension();
	for (uint32_t i = 0; i < dim; i++) {
		x[i] *= a;
	}
}

static float dot(const FullVector &x, const FullVector &y)
{
	XA_DEBUG_ASSERT(x.dimension() == y.dimension());
	const uint32_t dim = x.dimension();
	float sum = 0;
	for (uint32_t i = 0; i < dim; i++) {
		sum += x[i] * y[i];
	}
	return sum;
}

// y = M * x
static void mult(const Matrix &M, const FullVector &x, FullVector &y)
{
	uint32_t w = M.width();
	uint32_t h = M.height();
	XA_DEBUG_ASSERT( w == x.dimension() );
	XA_DEBUG_ASSERT( h == y.dimension() );
	for (uint32_t i = 0; i < h; i++)
		y[i] = M.dotRow(i, x);
#ifdef NDEBUG
	w = w; // silence unused parameter warning
#endif
}

static void sgemv(float alpha, Transpose TA, const Matrix &A, const FullVector &x, float beta, FullVector &y)
{
	uint32_t w = A.width();
	uint32_t h = A.height();
	if (TA == Transposed) {
		XA_DEBUG_ASSERT( h == x.dimension() );
		XA_DEBUG_ASSERT( w == y.dimension() );
		for (uint32_t i = 0; i < h; i++) {
			A.madRow(i, alpha * x[i], y);
		}
	} else {
		XA_DEBUG_ASSERT( w == x.dimension() );
		XA_DEBUG_ASSERT( h == y.dimension() );
		for (uint32_t i = 0; i < h; i++) {
			y[i] = alpha * A.dotRow(i, x) + beta * y[i];
		}
	}
#ifdef NDEBUG
	w = w; // silence unused parameter warning
	h = h;
#endif
}

// y = alpha*A*x + beta*y
static void sgemv(float alpha, const Matrix &A, const FullVector &x, float beta, FullVector &y)
{
	sgemv(alpha, NoTransposed, A, x, beta, y);
}

// dot y-row of A by x-column of B
static float dotRowColumn(int y, const Matrix &A, int x, const Matrix &B)
{
	const Array<Matrix::Coefficient> &row = A.getRow(y);
	const uint32_t count = row.size();
	float sum = 0.0f;
	for (uint32_t i = 0; i < count; i++) {
		const Matrix::Coefficient &c = row[i];
		sum += c.v * B.getCoefficient(x, c.x);
	}
	return sum;
}

static void transpose(const Matrix &A, Matrix &B)
{
	XA_DEBUG_ASSERT(A.width() == B.height());
	XA_DEBUG_ASSERT(B.width() == A.height());
	const uint32_t w = A.width();
	for (uint32_t x = 0; x < w; x++) {
		B.clearRow(x);
	}
	const uint32_t h = A.height();
	for (uint32_t y = 0; y < h; y++) {
		const Array<Matrix::Coefficient> &row = A.getRow(y);
		const uint32_t count = row.size();
		for (uint32_t i = 0; i < count; i++) {
			const Matrix::Coefficient &c = row[i];
			XA_DEBUG_ASSERT(c.x < w);
			B.setCoefficient(y, c.x, c.v);
		}
	}
}

static void sgemm(float alpha, const Matrix &A, const Matrix &B, float beta, Matrix &C)
{
	const uint32_t w = C.width();
	const uint32_t h = C.height();
#ifdef _DEBUG
	const uint32_t aw = A.width();
	const uint32_t ah = A.height();
	const uint32_t bw = B.width();
	const uint32_t bh = B.height();
	XA_DEBUG_ASSERT(aw == bh);
	XA_DEBUG_ASSERT(bw == ah);
	XA_DEBUG_ASSERT(w == bw);
	XA_DEBUG_ASSERT(h == ah);
#endif
	for (uint32_t y = 0; y < h; y++) {
		for (uint32_t x = 0; x < w; x++) {
			float c = beta * C.getCoefficient(x, y);
			// dot y-row of A by x-column of B.
			c += alpha * dotRowColumn(y, A, x, B);
			C.setCoefficient(x, y, c);
		}
	}
}

// C = A * B
static void mult(const Matrix &A, const Matrix &B, Matrix &C)
{
	sgemm(1.0f, A, B, 0.0f, C);
}

} // namespace sparse

class JacobiPreconditioner
{
public:
	JacobiPreconditioner(const sparse::Matrix &M, bool symmetric) : m_inverseDiagonal(M.width())
	{
		XA_ASSERT(M.isSquare());
		for (uint32_t x = 0; x < M.width(); x++) {
			float elem = M.getCoefficient(x, x);
			//XA_DEBUG_ASSERT( elem != 0.0f ); // This can be zero in the presence of zero area triangles.
			if (symmetric) {
				m_inverseDiagonal[x] = (elem != 0) ? 1.0f / sqrtf(fabsf(elem)) : 1.0f;
			} else {
				m_inverseDiagonal[x] = (elem != 0) ? 1.0f / elem : 1.0f;
			}
		}
	}

	void apply(const FullVector &x, FullVector &y) const
	{
		XA_DEBUG_ASSERT(x.dimension() == m_inverseDiagonal.dimension());
		XA_DEBUG_ASSERT(y.dimension() == m_inverseDiagonal.dimension());
		// @@ Wrap vector component-wise product into a separate function.
		const uint32_t D = x.dimension();
		for (uint32_t i = 0; i < D; i++) {
			y[i] = m_inverseDiagonal[i] * x[i];
		}
	}

private:
	FullVector m_inverseDiagonal;
};

// Linear solvers.
class Solver
{
public:
	// Solve the symmetric system: At·A·x = At·b
	static bool LeastSquaresSolver(const sparse::Matrix &A, const FullVector &b, FullVector &x, float epsilon = 1e-5f)
	{
		XA_DEBUG_ASSERT(A.width() == x.dimension());
		XA_DEBUG_ASSERT(A.height() == b.dimension());
		XA_DEBUG_ASSERT(A.height() >= A.width()); // @@ If height == width we could solve it directly...
		const uint32_t D = A.width();
		sparse::Matrix At(A.height(), A.width());
		sparse::transpose(A, At);
		FullVector Atb(D);
		sparse::mult(At, b, Atb);
		sparse::Matrix AtA(D);
		sparse::mult(At, A, AtA);
		return SymmetricSolver(AtA, Atb, x, epsilon);
	}

	// See section 10.4.3 in: Mesh Parameterization: Theory and Practice, Siggraph Course Notes, August 2007
	static bool LeastSquaresSolver(const sparse::Matrix &A, const FullVector &b, FullVector &x, const uint32_t *lockedParameters, uint32_t lockedCount, float epsilon = 1e-5f)
	{
		XA_DEBUG_ASSERT(A.width() == x.dimension());
		XA_DEBUG_ASSERT(A.height() == b.dimension());
		XA_DEBUG_ASSERT(A.height() >= A.width() - lockedCount);
		// @@ This is not the most efficient way of building a system with reduced degrees of freedom. It would be faster to do it on the fly.
		const uint32_t D = A.width() - lockedCount;
		XA_DEBUG_ASSERT(D > 0);
		// Compute: b - Al * xl
		FullVector b_Alxl(b);
		for (uint32_t y = 0; y < A.height(); y++) {
			const uint32_t count = A.getRow(y).size();
			for (uint32_t e = 0; e < count; e++) {
				uint32_t column = A.getRow(y)[e].x;
				bool isFree = true;
				for (uint32_t i = 0; i < lockedCount; i++) {
					isFree &= (lockedParameters[i] != column);
				}
				if (!isFree) {
					b_Alxl[y] -= x[column] * A.getRow(y)[e].v;
				}
			}
		}
		// Remove locked columns from A.
		sparse::Matrix Af(D, A.height());
		for (uint32_t y = 0; y < A.height(); y++) {
			const uint32_t count = A.getRow(y).size();
			for (uint32_t e = 0; e < count; e++) {
				uint32_t column = A.getRow(y)[e].x;
				uint32_t ix = column;
				bool isFree = true;
				for (uint32_t i = 0; i < lockedCount; i++) {
					isFree &= (lockedParameters[i] != column);
					if (column > lockedParameters[i]) ix--; // shift columns
				}
				if (isFree) {
					Af.setCoefficient(ix, y, A.getRow(y)[e].v);
				}
			}
		}
		// Remove elements from x
		FullVector xf(D);
		for (uint32_t i = 0, j = 0; i < A.width(); i++) {
			bool isFree = true;
			for (uint32_t l = 0; l < lockedCount; l++) {
				isFree &= (lockedParameters[l] != i);
			}
			if (isFree) {
				xf[j++] = x[i];
			}
		}
		// Solve reduced system.
		bool result = LeastSquaresSolver(Af, b_Alxl, xf, epsilon);
		// Copy results back to x.
		for (uint32_t i = 0, j = 0; i < A.width(); i++) {
			bool isFree = true;
			for (uint32_t l = 0; l < lockedCount; l++) {
				isFree &= (lockedParameters[l] != i);
			}
			if (isFree) {
				x[i] = xf[j++];
			}
		}
		return result;
	}

private:
	/**
	* Compute the solution of the sparse linear system Ab=x using the Conjugate
	* Gradient method.
	*
	* Solving sparse linear systems:
	* (1)		A·x = b
	*
	* The conjugate gradient algorithm solves (1) only in the case that A is
	* symmetric and positive definite. It is based on the idea of minimizing the
	* function
	*
	* (2)		f(x) = 1/2·x·A·x - b·x
	*
	* This function is minimized when its gradient
	*
	* (3)		df = A·x - b
	*
	* is zero, which is equivalent to (1). The minimization is carried out by
	* generating a succession of search directions p.k and improved minimizers x.k.
	* At each stage a quantity alfa.k is found that minimizes f(x.k + alfa.k·p.k),
	* and x.k+1 is set equal to the new point x.k + alfa.k·p.k. The p.k and x.k are
	* built up in such a way that x.k+1 is also the minimizer of f over the whole
	* vector space of directions already taken, {p.1, p.2, . . . , p.k}. After N
	* iterations you arrive at the minimizer over the entire vector space, i.e., the
	* solution to (1).
	*
	* For a really good explanation of the method see:
	*
	* "An Introduction to the Conjugate Gradient Method Without the Agonizing Pain",
	* Jonhathan Richard Shewchuk.
	*
	**/
	// Conjugate gradient with preconditioner.
	static bool ConjugateGradientSolver(const JacobiPreconditioner &preconditioner, const sparse::Matrix &A, const FullVector &b, FullVector &x, float epsilon)
	{
		XA_DEBUG_ASSERT( A.isSquare() );
		XA_DEBUG_ASSERT( A.width() == b.dimension() );
		XA_DEBUG_ASSERT( A.width() == x.dimension() );
		int i = 0;
		const int D = A.width();
		const int i_max = 4 * D;   // Convergence should be linear, but in some cases, it's not.
		FullVector r(D);    // residual
		FullVector p(D);    // search direction
		FullVector q(D);    //
		FullVector s(D);    // preconditioned
		float delta_0;
		float delta_old;
		float delta_new;
		float alpha;
		float beta;
		// r = b - A·x
		sparse::copy(b, r);
		sparse::sgemv(-1, A, x, 1, r);
		// p = M^-1 · r
		preconditioner.apply(r, p);
		delta_new = sparse::dot(r, p);
		delta_0 = delta_new;
		while (i < i_max && delta_new > epsilon * epsilon * delta_0) {
			i++;
			// q = A·p
			sparse::mult(A, p, q);
			// alpha = delta_new / p·q
			alpha = delta_new / sparse::dot(p, q);
			// x = alfa·p + x
			sparse::saxpy(alpha, p, x);
			if ((i & 31) == 0) { // recompute r after 32 steps
				// r = b - A·x
				sparse::copy(b, r);
				sparse::sgemv(-1, A, x, 1, r);
			} else {
				// r = r - alfa·q
				sparse::saxpy(-alpha, q, r);
			}
			// s = M^-1 · r
			preconditioner.apply(r, s);
			delta_old = delta_new;
			delta_new = sparse::dot( r, s );
			beta = delta_new / delta_old;
			// p = s + beta·p
			sparse::scal(beta, p);
			sparse::saxpy(1, s, p);
		}
		return delta_new <= epsilon * epsilon * delta_0;
	}

	static bool SymmetricSolver(const sparse::Matrix &A, const FullVector &b, FullVector &x, float epsilon = 1e-5f)
	{
		XA_DEBUG_ASSERT(A.height() == A.width());
		XA_DEBUG_ASSERT(A.height() == b.dimension());
		XA_DEBUG_ASSERT(b.dimension() == x.dimension());
		JacobiPreconditioner jacobi(A, true);
		return ConjugateGradientSolver(jacobi, A, b, x, epsilon);
	}
};

namespace param {
class Atlas;
class Chart;

// Fast sweep in 3 directions
static bool findApproximateDiameterVertices(Mesh *mesh, uint32_t *a, uint32_t *b)
{
	XA_DEBUG_ASSERT(a != NULL);
	XA_DEBUG_ASSERT(b != NULL);
	const uint32_t vertexCount = mesh->vertexCount();
	uint32_t minVertex[3];
	uint32_t maxVertex[3];
	minVertex[0] = minVertex[1] = minVertex[2] = UINT32_MAX;
	maxVertex[0] = maxVertex[1] = maxVertex[2] = UINT32_MAX;
	for (uint32_t v = 1; v < vertexCount; v++) {
		if (mesh->isBoundaryVertex(v)) {
			minVertex[0] = minVertex[1] = minVertex[2] = v;
			maxVertex[0] = maxVertex[1] = maxVertex[2] = v;
			break;
		}
	}
	if (minVertex[0] == UINT32_MAX) {
		// Input mesh has not boundaries.
		return false;
	}
	for (uint32_t v = 1; v < vertexCount; v++) {
		if (!mesh->isBoundaryVertex(v)) {
			// Skip interior vertices.
			continue;
		}
		const Vector3 *pos = mesh->positionAt(v);
		if (pos->x < mesh->positionAt(minVertex[0])->x)
			minVertex[0] = v;
		else if (pos->x > mesh->positionAt(maxVertex[0])->x)
			maxVertex[0] = v;
		if (pos->y < mesh->positionAt(minVertex[1])->y)
			minVertex[1] = v;
		else if (pos->y > mesh->positionAt(maxVertex[1])->y)
			maxVertex[1] = v;
		if (pos->z < mesh->positionAt(minVertex[2])->z)
			minVertex[2] = v;
		else if (pos->z > mesh->positionAt(maxVertex[2])->z)
			maxVertex[2] = v;
	}
	float lengths[3];
	for (int i = 0; i < 3; i++) {
		lengths[i] = length(*mesh->positionAt(minVertex[i]) - *mesh->positionAt(maxVertex[i]));
	}
	if (lengths[0] > lengths[1] && lengths[0] > lengths[2]) {
		*a = minVertex[0];
		*b = maxVertex[0];
	} else if (lengths[1] > lengths[2]) {
		*a = minVertex[1];
		*b = maxVertex[1];
	} else {
		*a = minVertex[2];
		*b = maxVertex[2];
	}
	return true;
}

// Conformal relations from Brecht Van Lommel (based on ABF):

static float vec_angle_cos(Vector3::Arg v1, Vector3::Arg v2, Vector3::Arg v3)
{
	Vector3 d1 = v1 - v2;
	Vector3 d2 = v3 - v2;
	return clamp(dot(d1, d2) / (length(d1) * length(d2)), -1.0f, 1.0f);
}

static float vec_angle(Vector3::Arg v1, Vector3::Arg v2, Vector3::Arg v3)
{
	float dot = vec_angle_cos(v1, v2, v3);
	return acosf(dot);
}

static void triangle_angles(Vector3::Arg v1, Vector3::Arg v2, Vector3::Arg v3, float *a1, float *a2, float *a3)
{
	*a1 = vec_angle(v3, v1, v2);
	*a2 = vec_angle(v1, v2, v3);
	*a3 = float(M_PI - *a2 - *a1);
}

static void setup_abf_relations(sparse::Matrix &A, int row, int id0, int id1, int id2, const Vector3 &p0, const Vector3 &p1, const Vector3 &p2)
{
	// @@ IC: Wouldn't it be more accurate to return cos and compute 1-cos^2?
	// It does indeed seem to be a little bit more robust.
	// @@ Need to revisit this more carefully!
	float a0, a1, a2;
	triangle_angles(p0, p1, p2, &a0, &a1, &a2);
	float s0 = sinf(a0);
	float s1 = sinf(a1);
	float s2 = sinf(a2);
	if (s1 > s0 && s1 > s2) {
		std::swap(s1, s2);
		std::swap(s0, s1);
		std::swap(a1, a2);
		std::swap(a0, a1);
		std::swap(id1, id2);
		std::swap(id0, id1);
	} else if (s0 > s1 && s0 > s2) {
		std::swap(s0, s2);
		std::swap(s0, s1);
		std::swap(a0, a2);
		std::swap(a0, a1);
		std::swap(id0, id2);
		std::swap(id0, id1);
	}
	float c0 = cosf(a0);
	float ratio = (s2 == 0.0f) ? 1.0f : s1 / s2;
	float cosine = c0 * ratio;
	float sine = s0 * ratio;
	// Note  : 2*id + 0 --> u
	//         2*id + 1 --> v
	int u0_id = 2 * id0 + 0;
	int v0_id = 2 * id0 + 1;
	int u1_id = 2 * id1 + 0;
	int v1_id = 2 * id1 + 1;
	int u2_id = 2 * id2 + 0;
	int v2_id = 2 * id2 + 1;
	// Real part
	A.setCoefficient(u0_id, 2 * row + 0, cosine - 1.0f);
	A.setCoefficient(v0_id, 2 * row + 0, -sine);
	A.setCoefficient(u1_id, 2 * row + 0, -cosine);
	A.setCoefficient(v1_id, 2 * row + 0, sine);
	A.setCoefficient(u2_id, 2 * row + 0, 1);
	// Imaginary part
	A.setCoefficient(u0_id, 2 * row + 1, sine);
	A.setCoefficient(v0_id, 2 * row + 1, cosine - 1.0f);
	A.setCoefficient(u1_id, 2 * row + 1, -sine);
	A.setCoefficient(v1_id, 2 * row + 1, -cosine);
	A.setCoefficient(v2_id, 2 * row + 1, 1);
}

static bool computeLeastSquaresConformalMap(Mesh *mesh)
{
	// For this to work properly, mesh should not have colocals that have the same
	// attributes, unless you want the vertices to actually have different texcoords.
	const uint32_t vertexCount = mesh->vertexCount();
	const uint32_t D = 2 * vertexCount;
	const uint32_t N = 2 * mesh->countTriangles();
	// N is the number of equations (one per triangle)
	// D is the number of variables (one per vertex; there are 2 pinned vertices).
	if (N < D - 4) {
		return false;
	}
	sparse::Matrix A(D, N);
	FullVector b(N);
	FullVector x(D);
	// Fill b:
	b.fill(0.0f);
	// Fill x:
	uint32_t v0, v1;
	if (!findApproximateDiameterVertices(mesh, &v0, &v1)) {
		// Mesh has no boundaries.
		return false;
	}
	if (*mesh->texcoordAt(v0) == *mesh->texcoordAt(v1)) {
		// LSCM expects an existing parameterization.
		return false;
	}
	for (uint32_t v = 0; v < vertexCount; v++) {
		const Vector2 *texcoord = mesh->texcoordAt(v);
		// Initial solution.
		x[2 * v + 0] = texcoord->x;
		x[2 * v + 1] = texcoord->y;
	}
	// Fill A:
	const uint32_t faceCount = mesh->faceCount();
	for (uint32_t f = 0, t = 0; f < faceCount; f++) {
		XA_DEBUG_ASSERT(mesh->faceAt(f)->nIndices == 3);
		uint32_t vertex0 = UINT32_MAX;
		for (Mesh::FaceEdgeIterator it(mesh, f); !it.isDone(); it.advance()) {
			if (vertex0 == UINT32_MAX) {
				vertex0 = it.vertex0();
			} else if (it.vertex1() != vertex0) {
				setup_abf_relations(A, t, vertex0, it.vertex0(), it.vertex1(), *mesh->positionAt(vertex0), it.position0(), it.position1());
				t++;
			}
		}
	}
	const uint32_t lockedParameters[] = {
		2 * v0 + 0,
		2 * v0 + 1,
		2 * v1 + 0,
		2 * v1 + 1
	};
	// Solve
	Solver::LeastSquaresSolver(A, b, x, lockedParameters, 4, 0.000001f);
	// Map x back to texcoords:
	for (uint32_t v = 0; v < vertexCount; v++) {
		*mesh->texcoordAt(v) = Vector2(x[2 * v + 0], x[2 * v + 1]);
	}
	return true;
}

static bool computeOrthogonalProjectionMap(Mesh *mesh)
{
	Vector3 axis[2];
	uint32_t vertexCount = mesh->vertexCount();
	Array<Vector3> points(vertexCount);
	points.resize(vertexCount);
	for (uint32_t i = 0; i < vertexCount; i++)
		points[i] = *mesh->positionAt(i);
	// Avoid redundant computations.
	float matrix[6];
	Fit::computeCovariance(vertexCount, points.data(), matrix);
	if (matrix[0] == 0 && matrix[3] == 0 && matrix[5] == 0) {
		return false;
	}
	float eigenValues[3];
	Vector3 eigenVectors[3];
	if (!Fit::eigenSolveSymmetric3(matrix, eigenValues, eigenVectors)) {
		return false;
	}
	axis[0] = normalize(eigenVectors[0]);
	axis[1] = normalize(eigenVectors[1]);
	// Project vertices to plane.
	for (uint32_t i = 0; i < vertexCount; i++)
		*mesh->texcoordAt(i) = Vector2(dot(axis[0], *mesh->positionAt(i)), dot(axis[1], *mesh->positionAt(i)));
	return true;
}

static void computeSingleFaceMap(Mesh *mesh)
{
	XA_DEBUG_ASSERT(mesh != NULL);
	XA_DEBUG_ASSERT(mesh->faceCount() == 1);
	Face *face = mesh->faceAt(0);
	XA_ASSERT(face != NULL);
	Vector3 p0 = *mesh->positionAt(mesh->vertexAt(face->firstIndex + 0));
	Vector3 p1 = *mesh->positionAt(mesh->vertexAt(face->firstIndex + 1));
	Vector3 X = normalizeSafe(p1 - p0, Vector3(0.0f), 0.0f);
	Vector3 Z = mesh->faceNormal(0);
	Vector3 Y = normalizeSafe(cross(Z, X), Vector3(0.0f), 0.0f);
	uint32_t i = 0;
	for (Mesh::FaceEdgeIterator it(mesh, 0); !it.isDone(); it.advance(), i++) {
		if (i == 0) {
			*mesh->texcoordAt(it.vertex0()) = Vector2(0);
		} else {
			Vector3 pn = it.position0();
			const float xn = dot((pn - p0), X);
			const float yn = dot((pn - p0), Y);
			*mesh->texcoordAt(it.vertex0()) = Vector2(xn, yn);
		}
	}
}

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
		pairs.insertAt(i, p);
		if (pairs.size() > maxSize)
			pairs.removeAt(0);
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

	Array<Pair> pairs;
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

	Array<uint32_t> seeds;
	Array<uint32_t> faces;
	PriorityQueue candidates;
};

struct AtlasBuilder
{
	AtlasBuilder(const Mesh *rm, const CharterOptions &options) : m_mesh(rm), m_facesLeft(rm->faceCount()), m_options(options)
	{
		const uint32_t faceCount = m_mesh->faceCount();
		m_faceChartArray.resize(faceCount, -1);
		m_faceCandidateArray.resize(faceCount, (uint32_t)-1);
		// @@ Floyd for the whole mesh is too slow. We could compute floyd progressively per patch as the patch grows. We need a better solution to compute most central faces.
		//computeShortestPaths();
		// Precompute edge lengths and face areas.
		const uint32_t edgeCount = m_mesh->edgeCount();
		m_edgeLengths.resize(edgeCount, 0.0f);
		m_faceAreas.resize(m_mesh->faceCount(), 0.0f);
		for (uint32_t f = 0; f < m_mesh->faceCount(); f++) {
			if ((m_mesh->faceFlagsAt(f) & FaceFlags::Ignore) != 0)
				continue;
			float &faceArea = m_faceAreas[f];
			Vector3 firstPos(0.0f);
			for (Mesh::FaceEdgeIterator it(m_mesh, f); !it.isDone(); it.advance()) {
				m_edgeLengths[it.edge()] = internal::length(it.position1() - it.position0());
				XA_DEBUG_ASSERT(m_edgeLengths[it.edge()] > 0.0f);
				if (it.relativeEdge() == 0)
					firstPos = it.position0();
				else
					faceArea += length(cross(it.position0() - firstPos, it.position1() - firstPos));
			}
			faceArea *= 0.5f;
			XA_DEBUG_ASSERT(faceArea > 0.0f);
		}
	}

	~AtlasBuilder()
	{
		const uint32_t chartCount = m_chartArray.size();
		for (uint32_t i = 0; i < chartCount; i++) {
			m_chartArray[i]->~ChartBuildData();
			XA_FREE(m_chartArray[i]);
		}
	}

	void placeSeeds(float threshold, uint32_t maxSeedCount)
	{
		// Instead of using a predefiened number of seeds:
		// - Add seeds one by one, growing chart until a certain treshold.
		// - Undo charts and restart growing process.
		// @@ How can we give preference to faces far from sharp features as in the LSCM paper?
		//   - those points can be found using a simple flood filling algorithm.
		//   - how do we weight the probabilities?
		for (uint32_t i = 0; i < maxSeedCount; i++) {
			if (m_facesLeft == 0) {
				// No faces left, stop creating seeds.
				break;
			}
			createRandomChart(threshold);
		}
	}

	void createRandomChart(float threshold)
	{
		const uint32_t randomFaceIdx = m_rand.getRange(m_facesLeft - 1);
		ChartBuildData *chart = XA_NEW(ChartBuildData, m_chartArray.size());
		m_chartArray.push_back(chart);
		// Pick random face that is not used by any chart yet.
		uint32_t i = 0;
		for (uint32_t f = 0; f != randomFaceIdx; f++, i++) {
			while (m_faceChartArray[i] != -1) i++;
		}
		while (m_faceChartArray[i] != -1) i++;
		chart->seeds.push_back(i);
		addFaceToChart(chart, i, true);
		// Grow the chart as much as possible within the given threshold.
		growChart(chart, threshold * 0.5f, m_facesLeft);
		//growCharts(threshold - threshold * 0.75f / chartCount(), facesLeft);
	}

	void addFaceToChart(ChartBuildData *chart, uint32_t f, bool recomputeProxy = false)
	{
		// Add face to chart.
		chart->faces.push_back(f);
		XA_DEBUG_ASSERT(m_faceChartArray[f] == -1);
		m_faceChartArray[f] = chart->id;
		m_facesLeft--;
		// Update area and boundary length.
		chart->area = evaluateChartArea(chart, f);
		chart->boundaryLength = evaluateBoundaryLength(chart, f);
		chart->normalSum = evaluateChartNormalSum(chart, f);
		chart->centroidSum = evaluateChartCentroidSum(chart, f);
		if (recomputeProxy) {
			// Update proxy and candidate's priorities.
			updateProxy(chart);
		}
		// Update candidates.
		removeCandidate(f);
		updateCandidates(chart, f);
		updatePriorities(chart);
	}

	// Returns true if any of the charts can grow more.
	bool growCharts(float threshold, uint32_t faceCount)
	{
		// Using one global list.
		faceCount = std::min(faceCount, m_facesLeft);
		for (uint32_t i = 0; i < faceCount; i++) {
			const Candidate &candidate = getBestCandidate();
			if (candidate.metric > threshold) {
				return false; // Can't grow more.
			}
			addFaceToChart(candidate.chart, candidate.face);
		}
		return m_facesLeft != 0; // Can continue growing.
	}

	bool growChart(ChartBuildData *chart, float threshold, uint32_t faceCount)
	{
		// Try to add faceCount faces within threshold to chart.
		for (uint32_t i = 0; i < faceCount; ) {
			if (chart->candidates.count() == 0 || chart->candidates.firstPriority() > threshold)
				return false;
			const uint32_t f = chart->candidates.pop();
			if (m_faceChartArray[f] == -1) {
				addFaceToChart(chart, f);
				i++;
			}
		}
		if (chart->candidates.count() == 0 || chart->candidates.firstPriority() > threshold)
			return false;
		return true;
	}

	void resetCharts()
	{
		const uint32_t faceCount = m_mesh->faceCount();
		for (uint32_t i = 0; i < faceCount; i++) {
			m_faceChartArray[i] = -1;
			m_faceCandidateArray[i] = (uint32_t)-1;
		}
		m_facesLeft = faceCount;
		m_candidateArray.clear();
		const uint32_t chartCount = m_chartArray.size();
		for (uint32_t i = 0; i < chartCount; i++) {
			ChartBuildData *chart = m_chartArray[i];
			const uint32_t seed = chart->seeds.back();
			chart->area = 0.0f;
			chart->boundaryLength = 0.0f;
			chart->normalSum = Vector3(0);
			chart->centroidSum = Vector3(0);
			chart->faces.clear();
			chart->candidates.clear();
			addFaceToChart(chart, seed);
		}
	}

	void updateCandidates(ChartBuildData *chart, uint32_t f)
	{
		// Traverse neighboring faces, add the ones that do not belong to any chart yet.
		for (Mesh::FaceEdgeIterator it(m_mesh, f); !it.isDone(); it.advance()) {
			if (it.oppositeEdge() == UINT32_MAX)
				continue;
			const Edge *oppositeEdge = m_mesh->edgeAt(it.oppositeEdge());
			if (m_faceChartArray[oppositeEdge->face] == -1) {
				XA_DEBUG_ASSERT(m_mesh->faceGroupAt(f) == m_mesh->faceGroupAt(oppositeEdge->face));
				chart->candidates.push(oppositeEdge->face);
			}
		}
	}

	void updateProxies()
	{
		const uint32_t chartCount = m_chartArray.size();
		for (uint32_t i = 0; i < chartCount; i++)
			updateProxy(m_chartArray[i]);
	}

	void updateProxy(ChartBuildData *chart) const
	{
		//#pragma message(NV_FILE_LINE "TODO: Use best fit plane instead of average normal.")
		chart->planeNormal = normalizeSafe(chart->normalSum, Vector3(0), 0.0f);
		chart->centroid = chart->centroidSum / float(chart->faces.size());
	}

	bool relocateSeeds()
	{
		bool anySeedChanged = false;
		const uint32_t chartCount = m_chartArray.size();
		for (uint32_t i = 0; i < chartCount; i++) {
			if (relocateSeed(m_chartArray[i])) {
				anySeedChanged = true;
			}
		}
		return anySeedChanged;
	}

	bool relocateSeed(ChartBuildData *chart) const
	{
		Vector3 centroid = computeChartCentroid(chart);
		const uint32_t N = 10;  // @@ Hardcoded to 10?
		PriorityQueue bestTriangles(N);
		// Find the first N triangles that fit the proxy best.
		const uint32_t faceCount = chart->faces.size();
		for (uint32_t i = 0; i < faceCount; i++) {
			float priority = evaluateProxyFitMetric(chart, chart->faces[i]);
			bestTriangles.push(priority, chart->faces[i]);
		}
		// Of those, choose the most central triangle.
		uint32_t mostCentral = 0;
		float maxDistance = -1;
		const uint32_t bestCount = bestTriangles.count();
		for (uint32_t i = 0; i < bestCount; i++) {
			Vector3 faceCentroid = m_mesh->triangleCenter(bestTriangles.pairs[i].face);
			float distance = length(centroid - faceCentroid);
			if (distance > maxDistance) {
				maxDistance = distance;
				mostCentral = bestTriangles.pairs[i].face;
			}
		}
		XA_DEBUG_ASSERT(maxDistance >= 0);
		// In order to prevent k-means cyles we record all the previously chosen seeds.
		for (uint32_t i = 0; i < chart->seeds.size(); i++) {
			if (chart->seeds[i] == mostCentral) {
				// Move new seed to the end of the seed array.
				uint32_t last = chart->seeds.size() - 1;
				std::swap(chart->seeds[i], chart->seeds[last]);
				return false;
			}
		}
		// Append new seed.
		chart->seeds.push_back(mostCentral);
		return true;
	}

	void updatePriorities(ChartBuildData *chart)
	{
		// Re-evaluate candidate priorities.
		uint32_t candidateCount = chart->candidates.count();
		for (uint32_t i = 0; i < candidateCount; i++) {
			PriorityQueue::Pair &pair = chart->candidates.pairs[i];
			pair.priority = evaluatePriority(chart, pair.face);
			if (m_faceChartArray[pair.face] == -1)
				updateCandidate(chart, pair.face, pair.priority);
		}
		// Sort candidates.
		chart->candidates.sort();
	}

	// Evaluate combined metric.
	float evaluatePriority(ChartBuildData *chart, uint32_t face) const
	{
		// Estimate boundary length and area:
		const float newBoundaryLength = evaluateBoundaryLength(chart, face);
		const float newChartArea = evaluateChartArea(chart, face);
		const float F = evaluateProxyFitMetric(chart, face);
		const float C = evaluateRoundnessMetric(chart, face, newBoundaryLength, newChartArea);
		const float P = evaluateStraightnessMetric(chart, face);
		// Penalize faces that cross seams, reward faces that close seams or reach boundaries.
		const float N = evaluateNormalSeamMetric(chart, face);
		const float T = evaluateTextureSeamMetric(chart, face);
		//float R = evaluateCompletenessMetric(chart, face);
		//float D = evaluateDihedralAngleMetric(chart, face);
		// @@ Add a metric based on local dihedral angle.
		// @@ Tweaking the normal and texture seam metrics.
		// - Cause more impedance. Never cross 90 degree edges.
		float cost =
			m_options.proxyFitMetricWeight * F +
			m_options.roundnessMetricWeight * C +
			m_options.straightnessMetricWeight * P +
			m_options.normalSeamMetricWeight * N +
			m_options.textureSeamMetricWeight * T;
		// Enforce limits strictly:
		if (newChartArea > m_options.maxChartArea)
			cost = FLT_MAX;
		if (newBoundaryLength > m_options.maxBoundaryLength)
			cost = FLT_MAX;
		// Make sure normal seams are fully respected:
		if (m_options.normalSeamMetricWeight >= 1000 && N != 0)
			cost = FLT_MAX;
		XA_DEBUG_ASSERT(std::isfinite(cost));
		return cost;
	}

	// Returns a value in [0-1].
	float evaluateProxyFitMetric(ChartBuildData *chart, uint32_t f) const
	{
		const Vector3 faceNormal = m_mesh->triangleNormal(f);
		// Use plane fitting metric for now:
		return 1 - dot(faceNormal, chart->planeNormal); // @@ normal deviations should be weighted by face area
	}

	float evaluateRoundnessMetric(ChartBuildData *chart, uint32_t /*face*/, float newBoundaryLength, float newChartArea) const
	{
		float roundness = square(chart->boundaryLength) / chart->area;
		float newRoundness = square(newBoundaryLength) / newChartArea;
		if (newRoundness > roundness) {
			return square(newBoundaryLength) / float(newChartArea * 4 * M_PI);
		} else {
			// Offer no impedance to faces that improve roundness.
			return 0;
		}
	}

	float evaluateStraightnessMetric(ChartBuildData *chart, uint32_t f) const
	{
		float l_out = 0.0f;
		float l_in = 0.0f;
		if (m_mesh->faceFlagsAt(f) & FaceFlags::Ignore)
			return 1.0f;
		for (Mesh::FaceEdgeIterator it(m_mesh, f); !it.isDone(); it.advance()) {
			float l = m_edgeLengths[it.edge()];
			if (it.isBoundary()) {
				l_out += l;
			} else {
				if (m_faceChartArray[it.oppositeFace()] != chart->id) {
					l_out += l;
				} else {
					l_in += l;
				}
			}
		}
		XA_DEBUG_ASSERT(l_in != 0.0f); // Candidate face must be adjacent to chart. @@ This is not true if the input mesh has zero-length edges.
		float ratio = (l_out - l_in) / (l_out + l_in);
		return std::min(ratio, 0.0f); // Only use the straightness metric to close gaps.
	}

	float evaluateNormalSeamMetric(ChartBuildData *chart, uint32_t f) const
	{
		float seamFactor = 0.0f;
		float totalLength = 0.0f;
		for (Mesh::FaceEdgeIterator it(m_mesh, f); !it.isDone(); it.advance()) {
			if (it.isBoundary())
				continue;
			if (m_faceChartArray[it.oppositeFace()] != chart->id)
				continue;
			float l = m_edgeLengths[it.edge()];
			totalLength += l;
			if (!it.isSeam())
				continue;
			// Make sure it's a normal seam.
			if (it.isNormalSeam()) {
				const Vector3 *n0 = m_mesh->normalAt(it.vertex0());
				const Vector3 *n1 = m_mesh->normalAt(it.vertex1());
				const Edge *oedge = m_mesh->edgeAt(it.oppositeEdge());
				const Vector3 *on0 = m_mesh->normalAt(m_mesh->vertexAt(oedge->index0));
				const Vector3 *on1 = m_mesh->normalAt(m_mesh->vertexAt(oedge->index1));
				const float d0 = clamp(dot(*n0, *on1), 0.0f, 1.0f);
				const float d1 = clamp(dot(*n1, *on0), 0.0f, 1.0f);
				l *= 1 - (d0 + d1) * 0.5f;
				seamFactor += l;
			}
		}
		if (seamFactor == 0)
			return 0.0f;
		return seamFactor / totalLength;
	}

	float evaluateTextureSeamMetric(ChartBuildData *chart, uint32_t f) const
	{
		float seamLength = 0.0f;
		float totalLength = 0.0f;
		for (Mesh::FaceEdgeIterator it(m_mesh, f); !it.isDone(); it.advance()) {
			if (it.isBoundary())
				continue;
			if (m_faceChartArray[it.oppositeFace()] != chart->id)
				continue;
			float l = m_edgeLengths[it.edge()];
			totalLength += l;
			if (!it.isSeam())
				continue;
			// Make sure it's a texture seam.
			if (it.isTextureSeam())
				seamLength += l;
		}
		if (seamLength == 0.0f)
			return 0.0f; // Avoid division by zero.
		return seamLength / totalLength;
	}

	float evaluateChartArea(ChartBuildData *chart, uint32_t f) const
	{
		return chart->area + m_faceAreas[f];
	}

	float evaluateBoundaryLength(ChartBuildData *chart, uint32_t f) const
	{
		float boundaryLength = chart->boundaryLength;
		// Add new edges, subtract edges shared with the chart.
		for (Mesh::FaceEdgeIterator it(m_mesh, f); !it.isDone(); it.advance()) {
			const float edgeLength = m_edgeLengths[it.edge()];
			if (it.isBoundary()) {
				boundaryLength += edgeLength;
			} else {
				if (m_faceChartArray[it.oppositeFace()] != chart->id)
					boundaryLength += edgeLength;
				else
					boundaryLength -= edgeLength;
			}
		}
		return std::max(0.0f, boundaryLength);  // @@ Hack!
	}

	Vector3 evaluateChartNormalSum(ChartBuildData *chart, uint32_t f) const
	{
		return chart->normalSum + m_mesh->triangleNormalAreaScaled(f);
	}

	Vector3 evaluateChartCentroidSum(ChartBuildData *chart, uint32_t f) const
	{
		return chart->centroidSum + m_mesh->faceCentroid(f);
	}

	Vector3 computeChartCentroid(const ChartBuildData *chart) const
	{
		Vector3 centroid(0);
		const uint32_t faceCount = chart->faces.size();
		for (uint32_t i = 0; i < faceCount; i++)
			centroid += m_mesh->triangleCenter(chart->faces[i]);
		return centroid / float(faceCount);
	}

	void fillHoles(float threshold)
	{
		while (m_facesLeft > 0)
			createRandomChart(threshold);
	}

	void mergeCharts()
	{
		Array<float> sharedBoundaryLengths;
		const uint32_t chartCount = m_chartArray.size();
		for (int c = chartCount - 1; c >= 0; c--) {
			sharedBoundaryLengths.clear();
			sharedBoundaryLengths.resize(chartCount, 0.0f);
			ChartBuildData *chart = m_chartArray[c];
			float externalBoundary = 0.0f;
			const uint32_t faceCount = chart->faces.size();
			for (uint32_t i = 0; i < faceCount; i++) {
				uint32_t f = chart->faces[i];
				for (Mesh::FaceEdgeIterator it(m_mesh, f); !it.isDone(); it.advance()) {
					const float l = m_edgeLengths[it.edge()];
					if (it.isBoundary()) {
						externalBoundary += l;
					} else {
						int neighborChart = m_faceChartArray[it.oppositeFace()];
						if (neighborChart != c) {
							if ((it.isSeam() && (it.isNormalSeam() || it.isTextureSeam())) || neighborChart == -2) {
								externalBoundary += l;
							} else {
								sharedBoundaryLengths[neighborChart] += l;
							}
						}
					}
				}
			}
			for (int cc = chartCount - 1; cc >= 0; cc--) {
				if (cc == c)
					continue;
				ChartBuildData *chart2 = m_chartArray[cc];
				if (chart2 == NULL)
					continue;
				if (sharedBoundaryLengths[cc] > 0.8 * std::max(0.0f, chart->boundaryLength - externalBoundary)) {
					// Try to avoid degenerate configurations.
					if (chart2->boundaryLength > sharedBoundaryLengths[cc]) {
						if (dot(chart2->planeNormal, chart->planeNormal) > -0.25) {
							mergeChart(chart2, chart, sharedBoundaryLengths[cc]);
							chart->~ChartBuildData();
							XA_FREE(chart);
							m_chartArray[c] = NULL;
							break;
						}
					}
				}
				if (sharedBoundaryLengths[cc] > 0.20 * std::max(0.0f, chart->boundaryLength - externalBoundary)) {
					// Compare proxies.
					if (dot(chart2->planeNormal, chart->planeNormal) > 0) {
						mergeChart(chart2, chart, sharedBoundaryLengths[cc]);
						chart->~ChartBuildData();
						XA_FREE(chart);
						m_chartArray[c] = NULL;
						break;
					}
				}
			}
		}
		// Remove deleted charts.
		for (int c = 0; c < int32_t(m_chartArray.size()); /*do not increment if removed*/) {
			if (m_chartArray[c] == NULL) {
				m_chartArray.removeAt(c);
				// Update m_faceChartArray.
				const uint32_t faceCount = m_faceChartArray.size();
				for (uint32_t i = 0; i < faceCount; i++) {
					XA_DEBUG_ASSERT(m_faceChartArray[i] != -1);
					XA_DEBUG_ASSERT(m_faceChartArray[i] != c);
					XA_DEBUG_ASSERT(m_faceChartArray[i] <= int32_t(m_chartArray.size()));
					if (m_faceChartArray[i] > c) {
						m_faceChartArray[i]--;
					}
				}
			} else {
				m_chartArray[c]->id = c;
				c++;
			}
		}
	}

	// @@ Cleanup.
	struct Candidate {
		uint32_t face;
		ChartBuildData *chart;
		float metric;
	};

	// @@ Get N best candidates in one pass.
	const Candidate &getBestCandidate() const
	{
		uint32_t best = 0;
		float bestCandidateMetric = FLT_MAX;
		const uint32_t candidateCount = m_candidateArray.size();
		XA_ASSERT(candidateCount > 0);
		for (uint32_t i = 0; i < candidateCount; i++) {
			const Candidate &candidate = m_candidateArray[i];
			if (candidate.metric < bestCandidateMetric) {
				bestCandidateMetric = candidate.metric;
				best = i;
			}
		}
		return m_candidateArray[best];
	}

	void removeCandidate(uint32_t f)
	{
		int c = m_faceCandidateArray[f];
		if (c != -1) {
			m_faceCandidateArray[f] = (uint32_t)-1;
			if (c == int(m_candidateArray.size() - 1)) {
				m_candidateArray.pop_back();
			} else {
				// Replace with last.
				m_candidateArray[c] = m_candidateArray[m_candidateArray.size() - 1];
				m_candidateArray.pop_back();
				m_faceCandidateArray[m_candidateArray[c].face] = c;
			}
		}
	}

	void updateCandidate(ChartBuildData *chart, uint32_t f, float metric)
	{
		if (m_faceCandidateArray[f] == (uint32_t)-1) {
			const uint32_t index = m_candidateArray.size();
			m_faceCandidateArray[f] = index;
			m_candidateArray.resize(index + 1);
			m_candidateArray[index].face = f;
			m_candidateArray[index].chart = chart;
			m_candidateArray[index].metric = metric;
		} else {
			const uint32_t c = m_faceCandidateArray[f];
			XA_DEBUG_ASSERT(c != (uint32_t)-1);
			Candidate &candidate = m_candidateArray[c];
			XA_DEBUG_ASSERT(candidate.face == f);
			if (metric < candidate.metric || chart == candidate.chart) {
				candidate.metric = metric;
				candidate.chart = chart;
			}
		}
	}

	void mergeChart(ChartBuildData *owner, ChartBuildData *chart, float sharedBoundaryLength)
	{
		const uint32_t faceCount = chart->faces.size();
		for (uint32_t i = 0; i < faceCount; i++) {
			uint32_t f = chart->faces[i];
			XA_DEBUG_ASSERT(m_faceChartArray[f] == chart->id);
			m_faceChartArray[f] = owner->id;
			owner->faces.push_back(f);
		}
		// Update adjacencies?
		owner->area += chart->area;
		owner->boundaryLength += chart->boundaryLength - sharedBoundaryLength;
		owner->normalSum += chart->normalSum;
		owner->centroidSum += chart->centroidSum;
		updateProxy(owner);
	}

	uint32_t facesLeft() const { return m_facesLeft;	}
	uint32_t chartCount() const { return m_chartArray.size(); }
	const Array<uint32_t> &chartFaces(uint32_t i) const { return m_chartArray[i]->faces; }

private:
	const Mesh *m_mesh;
	Array<float> m_edgeLengths;
	Array<float> m_faceAreas;
	uint32_t m_facesLeft;
	Array<int> m_faceChartArray;
	Array<ChartBuildData *> m_chartArray;
	Array<Candidate> m_candidateArray;
	Array<uint32_t> m_faceCandidateArray; // Map face index to candidate index.
	MTRand m_rand;
	CharterOptions m_options;
};

/// A chart is a connected set of faces with a certain topology (usually a disk).
class Chart
{
public:
	Chart(const Mesh *originalMesh, const Array<uint32_t> &faceArray) : atlasIndex(-1), m_blockAligned(true), m_chartMesh(NULL), m_unifiedMesh(NULL), m_isDisk(false)
	{
		// Copy face indices.
		m_faceArray = faceArray;
		const uint32_t meshVertexCount = originalMesh->vertexCount();
		m_chartMesh = XA_NEW(Mesh);
		m_unifiedMesh = XA_NEW(Mesh);
		Array<uint32_t> chartMeshIndices;
		chartMeshIndices.resize(meshVertexCount, (uint32_t)~0);
		Array<uint32_t> unifiedMeshIndices;
		unifiedMeshIndices.resize(meshVertexCount, (uint32_t)~0);
		// Add vertices.
		const uint32_t faceCount = faceArray.size();
		for (uint32_t f = 0; f < faceCount; f++) {
			for (Mesh::FaceEdgeIterator it(originalMesh, faceArray[f]); !it.isDone(); it.advance()) {
				const uint32_t vertex = it.vertex0();
				const uint32_t unifiedVertex = originalMesh->firstColocal(vertex);
				if (unifiedMeshIndices[unifiedVertex] == (uint32_t)~0) {
					unifiedMeshIndices[unifiedVertex] = m_unifiedMesh->vertexCount();
					XA_DEBUG_ASSERT(equal(it.position0(), *originalMesh->positionAt(unifiedVertex)));
					m_unifiedMesh->addVertex(it.position0());
				}
				if (chartMeshIndices[vertex] == (uint32_t)~0) {
					chartMeshIndices[vertex] = m_chartMesh->vertexCount();
					m_chartToOriginalMap.push_back(vertex);
					m_chartToUnifiedMap.push_back(unifiedMeshIndices[unifiedVertex]);
					m_chartMesh->addVertex(it.position0(), it.normal0(), it.texcoord0());
				}
			}
		}
		// This is ignoring the canonical map:
		// - Is it really necessary to link colocals?
		m_chartMesh->createColocals();
		Array<uint32_t> faceIndices;
		faceIndices.reserve(7);
		// Add faces.
		for (uint32_t f = 0; f < faceCount; f++) {
			const uint32_t faceFlags = originalMesh->faceFlagsAt(faceArray[f]);
			faceIndices.clear();
			for (Mesh::FaceEdgeIterator it(originalMesh, faceArray[f]); !it.isDone(); it.advance())
				faceIndices.push_back(chartMeshIndices[it.vertex0()]);
			m_chartMesh->addFace(faceIndices, faceFlags);
			faceIndices.clear();
			for (Mesh::FaceEdgeIterator it(originalMesh, faceArray[f]); !it.isDone(); it.advance()) {
				uint32_t unifiedVertex = originalMesh->firstColocal(it.vertex0());
				if (unifiedVertex == UINT32_MAX)
					unifiedVertex = it.vertex0();
				faceIndices.push_back(unifiedMeshIndices[unifiedVertex]);
			}
			m_unifiedMesh->addFace(faceIndices, faceFlags);
		}
		m_chartMesh->createBoundaries();
		m_chartMesh->linkBoundaries();
		m_unifiedMesh->createBoundaries();
		m_unifiedMesh->linkBoundaries();
		Mesh *splitUnifiedMesh = meshSplitBoundaryEdges(*m_unifiedMesh);
		if (splitUnifiedMesh) {
			m_unifiedMesh->~Mesh();
			XA_FREE(m_unifiedMesh);
			m_unifiedMesh = meshUnifyVertices(*splitUnifiedMesh);
			splitUnifiedMesh->~Mesh();
			XA_FREE(splitUnifiedMesh);
		}
		// Closing the holes is not always the best solution and does not fix all the problems.
		// We need to do some analysis of the holes and the genus to:
		// - Find cuts that reduce genus.
		// - Find cuts to connect holes.
		// - Use minimal spanning trees or seamster.
		bool closed = closeHoles();
#if XA_DEBUG_EXPORT_OBJ
		if (!closed)
			m_unifiedMesh->writeSimpleObj("debug_chart_not_closed.obj");
#endif
		XA_DEBUG_ASSERT(closed);
		closed = closed; // silence unused parameter warning;
		Mesh *triangulatedMesh = meshTriangulate(*m_unifiedMesh);
		if (triangulatedMesh) {
			m_unifiedMesh->~Mesh();
			XA_FREE(m_unifiedMesh);
			m_unifiedMesh = triangulatedMesh;
		}
		MeshTopology topology(m_unifiedMesh);
		m_isDisk = topology.isDisk();
#if XA_DEBUG_EXPORT_OBJ
		if (!m_isDisk)
			m_unifiedMesh->writeSimpleObj("debug_chart_not_disk.obj");
#endif
		XA_DEBUG_ASSERT(m_isDisk);
	}

	~Chart()
	{
		if (m_chartMesh) {
			m_chartMesh->~Mesh();
			XA_FREE(m_chartMesh);
		}
		if (m_unifiedMesh) {
			m_unifiedMesh->~Mesh();
			XA_FREE(m_unifiedMesh);
		}
	}

	bool isBlockAligned() const { return m_blockAligned; }
	bool isDisk() const { return m_isDisk; }
	uint32_t vertexCount() const { return m_chartMesh->vertexCount(); }
	uint32_t colocalVertexCount() const { return m_unifiedMesh->vertexCount(); }
	uint32_t faceCount() const { return m_faceArray.size(); }
	uint32_t mapFaceToSourceFace(uint32_t i) const { return m_faceArray[i]; }
	const Mesh *chartMesh() const { return m_chartMesh; }
	Mesh *chartMesh() { return m_chartMesh; }
	const Mesh *unifiedMesh() const { return m_unifiedMesh; }
	Mesh *unifiedMesh() { return m_unifiedMesh; }
	uint32_t mapChartVertexToOriginalVertex(uint32_t i) const { return m_chartToOriginalMap[i]; }

	// Transfer parameterization from unified mesh to chart mesh.
	void transferParameterization()
	{
		const uint32_t vertexCount = m_chartMesh->vertexCount();
		for (uint32_t v = 0; v < vertexCount; v++) {
			Vector2 *texcoord = m_chartMesh->texcoordAt(v);
			*texcoord = *m_unifiedMesh->texcoordAt(m_chartToUnifiedMap[v]);
		}
	}

	float computeSurfaceArea() const
	{
		return m_chartMesh->computeSurfaceArea();
	}

	float computeParametricArea() const
	{
		// This only makes sense in parameterized meshes.
		XA_DEBUG_ASSERT(m_isDisk);
		return m_chartMesh->computeParametricArea();
	}

	Vector2 computeParametricBounds() const
	{
		// This only makes sense in parameterized meshes.
		XA_DEBUG_ASSERT(m_isDisk);
		Vector2 minCorner(FLT_MAX, FLT_MAX);
		Vector2 maxCorner(-FLT_MAX, -FLT_MAX);
		const uint32_t vertexCount = m_chartMesh->vertexCount();
		for (uint32_t v = 0; v < vertexCount; v++) {
			const Vector2 *tex = m_chartMesh->texcoordAt(v);
			minCorner = min(minCorner, *tex);
			maxCorner = max(maxCorner, *tex);
		}
		return (maxCorner - minCorner) * 0.5f;
	}

	int32_t atlasIndex;

private:
	bool closeHoles()
	{
		Array<uint32_t> boundaryEdges;
		meshGetBoundaryEdges(*m_unifiedMesh, boundaryEdges);
		uint32_t boundaryCount = boundaryEdges.size();
		if (boundaryCount <= 1) {
			// Nothing to close.
			return true;
		}
		// Compute lengths.
		Array<float> boundaryLengths;
		for (uint32_t i = 0; i < boundaryCount; i++) {
			float boundaryLength = 0.0f;
			for (Mesh::BoundaryEdgeIterator it(m_unifiedMesh, boundaryEdges[i]); !it.isDone(); it.advance()) {
				const Edge *edge = m_unifiedMesh->edgeAt(it.edge());
				Vector3 t0 = *m_unifiedMesh->positionAt(m_unifiedMesh->vertexAt(edge->index0));
				Vector3 t1 = *m_unifiedMesh->positionAt(m_unifiedMesh->vertexAt(edge->index1));
				boundaryLength += length(t1 - t0);
			}
			boundaryLength = boundaryLength;
			boundaryLengths.push_back(boundaryLength);
		}
		// Find disk boundary.
		uint32_t diskBoundary = 0;
		float maxLength = boundaryLengths[0];
		for (uint32_t i = 1; i < boundaryCount; i++) {
			if (boundaryLengths[i] > maxLength) {
				maxLength = boundaryLengths[i];
				diskBoundary = i;
			}
		}
		// Close holes.
		for (uint32_t i = 0; i < boundaryCount; i++) {
			if (diskBoundary == i) {
				// Skip disk boundary.
				continue;
			}
			Array<uint32_t> vertexLoop;
			Array<const Edge *> edgeLoop;
			startOver:
			for (Mesh::BoundaryEdgeIterator it(m_unifiedMesh, boundaryEdges[i]); !it.isDone(); it.advance()) {
				const Edge *edge = m_unifiedMesh->edgeAt(it.edge());
				const uint32_t vertex = m_unifiedMesh->vertexAt(edge->index1);
				uint32_t j;
				for (j = 0; j < vertexLoop.size(); j++) {
					if (m_unifiedMesh->areColocal(vertex, vertexLoop[j]))
						break;
				}
				bool isCrossing = (j != vertexLoop.size());
				if (isCrossing) {
					// Close loop.
					edgeLoop.insertAt(0, edge);
					closeLoop(j + 1, edgeLoop);
					// Start over again.
					vertexLoop.clear();
					edgeLoop.clear();
					goto startOver; // HE mesh version is bugged, actually breaks at end of edge iteration instead.
				}
				vertexLoop.push_back(vertex);
				edgeLoop.insertAt(0, edge);
			}
			closeLoop(0, edgeLoop);
		}
		m_unifiedMesh->createBoundaries();
		m_unifiedMesh->linkBoundaries();
		meshGetBoundaryEdges(*m_unifiedMesh, boundaryEdges);
		boundaryCount = boundaryEdges.size();
		return boundaryCount == 1;
	}

	bool closeLoop(uint32_t startVertex, const Array<const Edge *> &loop)
	{
		const uint32_t vertexCount = loop.size() - startVertex;
		XA_DEBUG_ASSERT(vertexCount >= 3);
		if (vertexCount < 3)
			return false;
		// If the hole is planar, then we add a single face that will be properly triangulated later.
		// If the hole is not planar, we add a triangle fan with a vertex at the hole centroid.
		// This is still a bit of a hack. There surely are better hole filling algorithms out there.
		Array<Vector3> points;
		points.resize(vertexCount);
		for (uint32_t i = 0; i < vertexCount; i++)
			points[i] = *m_unifiedMesh->positionAt(m_unifiedMesh->vertexAt(loop[startVertex + i]->index0));
		const bool isPlanar = Fit::isPlanar(vertexCount, points.data());
		if (isPlanar) {
			Array<uint32_t> indices;
			indices.resize(vertexCount);
			for (uint32_t i = 0; i < vertexCount; i++)
				indices[i] = m_unifiedMesh->vertexAt(loop[startVertex + i]->index0);
			m_unifiedMesh->addFace(indices);
		} else {
			// If the polygon is not planar, we just cross our fingers, and hope this will work:
			// Compute boundary centroid:
			Vector3 centroidPos(0.0f);
			for (uint32_t i = 0; i < vertexCount; i++)
				centroidPos += points[i];
			centroidPos *= (1.0f / vertexCount);
			const uint32_t centroidVertex = m_unifiedMesh->vertexCount();
			m_unifiedMesh->addVertex(centroidPos);
			// Add one pair of edges for each boundary vertex.
			for (uint32_t j = vertexCount - 1, i = 0; i < vertexCount; j = i++) {
				const uint32_t vertex1 = m_unifiedMesh->vertexAt(loop[startVertex + j]->index0);
				const uint32_t vertex2 = m_unifiedMesh->vertexAt(loop[startVertex + i]->index0);
				m_unifiedMesh->addFace(centroidVertex, vertex1, vertex2);
			}
		}
		return true;
	}

	bool m_blockAligned;

	Mesh *m_chartMesh;
	Mesh *m_unifiedMesh;
	bool m_isDisk;

	// List of faces of the original mesh that belong to this chart.
	Array<uint32_t> m_faceArray;

	// Map vertices of the chart mesh to vertices of the original mesh.
	Array<uint32_t> m_chartToOriginalMap;

	Array<uint32_t> m_chartToUnifiedMap;
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

	ParameterizationQuality(const Mesh *mesh)
	{
		XA_DEBUG_ASSERT(mesh != NULL);
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
			uint32_t vertex0 = UINT32_MAX;
			Vector3 p[3];
			Vector2 t[3];
			for (Mesh::FaceEdgeIterator it(mesh, f); !it.isDone(); it.advance()) {
				if (vertex0 == UINT32_MAX) {
					vertex0 = it.vertex0();
					p[0] = it.position0();
					t[0] = it.texcoord0();
				} else if (it.vertex1() != vertex0) {
					p[1] = it.position0();
					p[2] = it.position1();
					t[1] = it.texcoord0();
					t[2] = it.texcoord1();
					processTriangle(p, t);
				}
			}
		}
		if (m_flippedTriangleCount + m_zeroAreaTriangleCount == faceCount) {
			// If all triangles are flipped, then none is.
			m_flippedTriangleCount = 0;
		}
		XA_DEBUG_ASSERT(std::isfinite(m_parametricArea) && m_parametricArea >= 0);
		XA_DEBUG_ASSERT(std::isfinite(m_geometricArea) && m_geometricArea >= 0);
		XA_DEBUG_ASSERT(std::isfinite(m_stretchMetric));
		XA_DEBUG_ASSERT(std::isfinite(m_maxStretchMetric));
		XA_DEBUG_ASSERT(std::isfinite(m_conformalMetric));
		XA_DEBUG_ASSERT(std::isfinite(m_authalicMetric));
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

	bool operator==(const ParameterizationQuality &pq)
	{
		if (m_totalTriangleCount != pq.m_totalTriangleCount || m_flippedTriangleCount != pq.m_flippedTriangleCount || m_zeroAreaTriangleCount != pq.m_zeroAreaTriangleCount)
			return false;
		if (!equal(m_parametricArea, pq.m_parametricArea))
			return false;
		if (!equal(m_geometricArea, pq.m_geometricArea))
			return false;
		if (!equal(m_stretchMetric, pq.m_stretchMetric))
			return false;
		if (!equal(m_maxStretchMetric, pq.m_maxStretchMetric))
			return false;
		if (!equal(m_conformalMetric, pq.m_conformalMetric))
			return false;
		if (!equal(m_authalicMetric, pq.m_authalicMetric))
			return false;
		return true;
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
		XA_ASSERT(sigma2 > sigma1 || equal(sigma1, sigma2));
		// isometric: sigma1 = sigma2 = 1
		// conformal: sigma1 / sigma2 = 1
		// authalic: sigma1 * sigma2 = 1
		float rmsStretch = sqrtf((a + c) * 0.5f);
		float rmsStretch2 = sqrtf((square(sigma1) + square(sigma2)) * 0.5f);
		XA_DEBUG_ASSERT(equal(rmsStretch, rmsStretch2, 0.01f));
#ifdef NDEBUG
		rmsStretch2 = rmsStretch2; // silence unused parameter warning
#endif
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

// Set of charts corresponding to mesh faces in the same face group.
class ChartGroup
{
public:
	ChartGroup(uint32_t id, const Mesh *sourceMesh, uint32_t faceGroup) : m_sourceId(sourceMesh->id()), m_id(id), m_isVertexMap(faceGroup == UINT32_MAX)
	{
		// Create new mesh from the source mesh, using faces that belong to this group.
		const uint32_t sourceFaceCount = sourceMesh->faceCount();
		for (uint32_t f = 0; f < sourceFaceCount; f++) {
			if (sourceMesh->faceGroupAt(f) == faceGroup)
				m_faceArray.push_back(f);
		}
		m_mesh = XA_NEW(Mesh);
		const uint32_t faceCount = m_faceArray.size();
		XA_DEBUG_ASSERT(faceCount > 0);
		Array<uint32_t> meshIndices;
		meshIndices.resize(sourceMesh->vertexCount(), (uint32_t)~0);
		for (uint32_t f = 0; f < faceCount; f++) {
			const Face *face = sourceMesh->faceAt(m_faceArray[f]);
			XA_DEBUG_ASSERT(face != NULL);
			for (uint32_t i = 0; i < face->nIndices; i++) {
				const uint32_t vertex = sourceMesh->vertexAt(face->firstIndex + i);
				if (meshIndices[vertex] == (uint32_t)~0) {
					meshIndices[vertex] = m_mesh->vertexCount();
					m_vertexToSourceVertexMap.push_back(vertex);
					m_mesh->addVertex(*sourceMesh->positionAt(vertex), *sourceMesh->normalAt(vertex), *sourceMesh->texcoordAt(vertex));
				}
			}
		}
		m_mesh->createColocals();
		Array<uint32_t> faceIndices;
		faceIndices.reserve(7);
		// Add faces.
		for (uint32_t f = 0; f < faceCount; f++) {
			const Face *face = sourceMesh->faceAt(m_faceArray[f]);
			faceIndices.clear();
			for (uint32_t i = 0; i < face->nIndices; i++) {
				const uint32_t vertex = sourceMesh->vertexAt(face->firstIndex + i);
				XA_DEBUG_ASSERT(meshIndices[vertex] != (uint32_t)~0);
				faceIndices.push_back(meshIndices[vertex]);
			}
			m_mesh->addFace(faceIndices);
		}
		m_mesh->createBoundaries();
		m_mesh->linkBoundaries();
#if XA_DEBUG_EXPORT_OBJ
		char filename[256];
		sprintf(filename, "debug_mesh_%0.3u_chartgroup_%0.3u.obj", m_sourceId, m_id);
		m_mesh->writeSimpleObj(filename);
#endif
	}

	~ChartGroup()
	{
		m_mesh->~Mesh();
		XA_FREE(m_mesh);
		for (uint32_t i = 0; i < m_chartArray.size(); i++) {
			m_chartArray[i]->~Chart();
			XA_FREE(m_chartArray[i]);
		}
	}

	uint32_t chartCount() const { return m_chartArray.size(); }
	Chart *chartAt(uint32_t i) const { return m_chartArray[i]; }
	bool isVertexMap() const { return m_isVertexMap; }
	uint32_t mapFaceToSourceFace(uint32_t face) const { return m_faceArray[face]; }
	uint32_t mapVertexToSourceVertex(uint32_t i) const { return m_vertexToSourceVertexMap[i]; }
	const Mesh *mesh() const { return m_mesh; }

	/*
	Compute charts using a simple segmentation algorithm.

	LSCM:
	- identify sharp features using local dihedral angles.
	- identify seed faces farthest from sharp features.
	- grow charts from these seeds.

	MCGIM:
	- phase 1: chart growth
	  - grow all charts simultaneously using dijkstra search on the dual graph of the mesh.
	  - graph edges are weighted based on planarity metric.
	  - metric uses distance to global chart normal.
	  - terminate when all faces have been assigned.
	- phase 2: seed computation:
	  - place new seed of the chart at the most interior face.
	  - most interior is evaluated using distance metric only.

	- method repeates the two phases, until the location of the seeds does not change.
	  - cycles are detected by recording all the previous seeds and chartification terminates.

	D-Charts:

	- Uniaxial conic metric:
	  - N_c = axis of the generalized cone that best fits the chart. (cone can a be cylinder or a plane).
	  - omega_c = angle between the face normals and the axis.
	  - Fitting error between chart C and tringle t: F(c,t) = (N_c*n_t - cos(omega_c))^2

	- Compactness metrics:
	  - Roundness:
		- C(c,t) = pi * D(S_c,t)^2 / A_c
		- S_c = chart seed.
		- D(S_c,t) = length of the shortest path inside the chart betwen S_c and t.
		- A_c = chart area.
	  - Straightness:
		- P(c,t) = l_out(c,t) / l_in(c,t)
		- l_out(c,t) = lenght of the edges not shared between C and t.
		- l_in(c,t) = lenght of the edges shared between C and t.

	- Combined metric:
	  - Cost(c,t) = F(c,t)^alpha + C(c,t)^beta + P(c,t)^gamma
	  - alpha = 1, beta = 0.7, gamma = 0.5

	Our basic approach:
	- Just one iteration of k-means?
	- Avoid dijkstra by greedily growing charts until a threshold is met. Increase threshold and repeat until no faces left.
	- If distortion metric is too high, split chart, add two seeds.
	- If chart size is low, try removing chart.

	Postprocess:
	- If topology is not disk:
	  - Fill holes, if new faces fit proxy.
	  - Find best cut, otherwise.
	- After parameterization:
	  - If boundary self-intersects:
		- cut chart along the closest two diametral boundary vertices, repeat parametrization.
		- what if the overlap is on an appendix? How do we find that out and cut appropiately?
		  - emphasize roundness metrics to prevent those cases.
	  - If interior self-overlaps: preserve boundary parameterization and use mean-value map.
	*/
	void computeCharts(const CharterOptions &options)
	{
		AtlasBuilder builder(m_mesh, options);
		if (builder.facesLeft() != 0) {
			// This seems a reasonable estimate.
			uint32_t maxSeedCount = std::max(6U, builder.facesLeft());
			// Create initial charts greedely.
			XA_PRINT(PrintFlags::ComputingCharts, "### Placing seeds\n");
			builder.placeSeeds(options.maxThreshold, maxSeedCount);
			XA_PRINT(PrintFlags::ComputingCharts, "###   Placed %d seeds (max = %d)\n", builder.chartCount(), maxSeedCount);
			builder.updateProxies();
			builder.mergeCharts();
#if 1
			XA_PRINT(PrintFlags::ComputingCharts, "### Relocating seeds\n");
			builder.relocateSeeds();
			XA_PRINT(PrintFlags::ComputingCharts, "### Reset charts\n");
			builder.resetCharts();
			XA_PRINT(PrintFlags::ComputingCharts, "### Growing charts\n");
			// Restart process growing charts in parallel.
			uint32_t iteration = 0;
			while (true) {
				if (!builder.growCharts(options.maxThreshold, options.growFaceCount)) {
					XA_PRINT(PrintFlags::ComputingCharts, "### Can't grow anymore\n");
					// If charts cannot grow more: fill holes, merge charts, relocate seeds and start new iteration.
					XA_PRINT(PrintFlags::ComputingCharts, "### Filling holes\n");
					builder.fillHoles(options.maxThreshold);
					XA_PRINT(PrintFlags::ComputingCharts, "###   Using %d charts now\n", builder.chartCount());
					builder.updateProxies();
					XA_PRINT(PrintFlags::ComputingCharts, "### Merging charts\n");
					builder.mergeCharts();
					XA_PRINT(PrintFlags::ComputingCharts, "###   Using %d charts now\n", builder.chartCount());
					XA_PRINT(PrintFlags::ComputingCharts, "### Reseeding\n");
					if (!builder.relocateSeeds()) {
						XA_PRINT(PrintFlags::ComputingCharts, "### Cannot relocate seeds anymore\n");
						// Done!
						break;
					}
					if (iteration == options.maxIterations) {
						XA_PRINT(PrintFlags::ComputingCharts, "### Reached iteration limit\n");
						break;
					}
					iteration++;
					XA_PRINT(PrintFlags::ComputingCharts, "### Reset charts\n");
					builder.resetCharts();
					XA_PRINT(PrintFlags::ComputingCharts, "### Growing charts\n");
				}
			}
#endif
			// Make sure no holes are left!
			XA_DEBUG_ASSERT(builder.facesLeft() == 0);
			const uint32_t chartCount = builder.chartCount();
			for (uint32_t i = 0; i < chartCount; i++) {
				Chart *chart = XA_NEW(Chart, m_mesh, builder.chartFaces(i));
				m_chartArray.push_back(chart);
#if XA_DEBUG_EXPORT_OBJ && XA_DEBUG_EXPORT_OBJ_INDIVIDUAL_CHARTS
				char filename[256];
				sprintf(filename, "debug_chart_%0.4d.obj", i);
				chart->chartMesh()->writeSimpleObj(filename);
				sprintf(filename, "debug_chart_%0.4d_unified.obj", i);
				chart->unifiedMesh()->writeSimpleObj(filename);
#endif
			}
#if XA_DEBUG_EXPORT_OBJ
			char filename[256];
			sprintf(filename, "debug_mesh_%0.3u_chartgroup_%0.3u_charts.obj", m_sourceId, m_id);
			FILE *file = fopen(filename, "w");
			if (file) {
				m_mesh->writeObjVertices(file);
				for (uint32_t i = 0; i < chartCount; i++) {
					fprintf(file, "o chart_%0.4d\n", i);
					fprintf(file, "s off\n");
					const Array<uint32_t> &faces = builder.chartFaces(i);
					for (uint32_t f = 0; f < faces.size(); f++)
						m_mesh->writeObjFace(file, faces[f]);
				}
				m_mesh->writeObjBoundaryEges(file);
				m_mesh->writeObjLinkedBoundaries(file);
				fclose(file);
			}
#endif
		}
	}

	void parameterizeCharts()
	{
		ParameterizationQuality globalParameterizationQuality;
		// Parameterize the charts.
		uint32_t diskCount = 0;
		const uint32_t chartCount = m_chartArray.size();
		for (uint32_t i = 0; i < chartCount; i++)
		{
			Chart *chart = m_chartArray[i];
			bool isValid = false;
			if (chart->isDisk())
			{
				diskCount++;
				ParameterizationQuality chartParameterizationQuality;
				if (chart->faceCount() == 1) {
					computeSingleFaceMap(chart->unifiedMesh());
					ParameterizationQuality quality = ParameterizationQuality(chart->unifiedMesh());
					chartParameterizationQuality = quality;
				} else {
					computeOrthogonalProjectionMap(chart->unifiedMesh());
					ParameterizationQuality orthogonalQuality(chart->unifiedMesh());
					computeLeastSquaresConformalMap(chart->unifiedMesh());
					ParameterizationQuality lscmQuality(chart->unifiedMesh());
					chartParameterizationQuality = lscmQuality;
				}
				isValid = chartParameterizationQuality.isValid();
				if (!isValid)
					XA_PRINT(PrintFlags::ParametizingCharts, "*** Invalid parameterization.\n");
				// @@ Check that parameterization quality is above a certain threshold.
				// @@ Detect boundary self-intersections.
				globalParameterizationQuality += chartParameterizationQuality;
			}
			// Transfer parameterization from unified mesh to chart mesh.
			chart->transferParameterization();
		}
		XA_PRINT(PrintFlags::ParametizingCharts, "  Parameterized %d/%d charts.\n", diskCount, chartCount);
		XA_PRINT(PrintFlags::ParametizingCharts, "  RMS stretch metric: %f\n", globalParameterizationQuality.rmsStretchMetric());
		XA_PRINT(PrintFlags::ParametizingCharts, "  MAX stretch metric: %f\n", globalParameterizationQuality.maxStretchMetric());
		XA_PRINT(PrintFlags::ParametizingCharts, "  RMS conformal metric: %f\n", globalParameterizationQuality.rmsConformalMetric());
		XA_PRINT(PrintFlags::ParametizingCharts, "  RMS authalic metric: %f\n", globalParameterizationQuality.maxAuthalicMetric());
	}

private:
	uint32_t m_sourceId, m_id;
	bool m_isVertexMap;
	Mesh *m_mesh;
	Array<uint32_t> m_faceArray; // List of faces of the source mesh that belong to this chart group.
	Array<uint32_t> m_vertexToSourceVertexMap; // Map vertices of the mesh to vertices of the source mesh.
	Array<Chart *> m_chartArray;
};

/// An atlas is a set of chart groups.
class Atlas
{
public:
	~Atlas()
	{
		for (uint32_t i = 0; i < m_chartGroups.size(); i++) {
			m_chartGroups[i]->~ChartGroup();
			XA_FREE(m_chartGroups[i]);
		}
	}

	uint32_t chartGroupCount(uint32_t mesh) const
	{
		uint32_t count = 0;
		for (uint32_t i = 0; i < m_chartGroups.size(); i++) {
			if (m_chartGroupSourceMeshes[i] == mesh)
				count++;
		}
		return count;
	}

	const ChartGroup *chartGroupAt(uint32_t mesh, uint32_t group) const
	{
		for (uint32_t c = 0; c < m_chartGroups.size(); c++) {
			if (m_chartGroupSourceMeshes[c] != mesh)
				continue;
			if (group == 0)
				return m_chartGroups[c];
			group--;
		}
		return NULL;
	}

	uint32_t chartCount() const
	{
		uint32_t count = 0;
		for (uint32_t i = 0; i < m_chartGroups.size(); i++)
			count += m_chartGroups[i]->chartCount();
		return count;
	}

	Chart *chartAt(uint32_t i)
	{
		for (uint32_t c = 0; c < m_chartGroups.size(); c++) {
			uint32_t count = m_chartGroups[c]->chartCount();
			if (i < count) {
				return m_chartGroups[c]->chartAt(i);
			}
			i -= count;
		}
		return NULL;
	}

	void computeCharts(const Array<Mesh *> &meshes, const CharterOptions &options, ProgressCallback progressCallback, void *progressCallbackUserData)
	{
		int progress = 0;
		if (progressCallback)
			progressCallback(ProgressCategory::ComputingCharts, 0, progressCallbackUserData);
		const uint32_t meshCount = meshes.size();
		for (uint32_t i = 0; i < meshCount; i++) {
			// Get list of face groups.
			const uint32_t faceCount = meshes[i]->faceCount();
			Array<uint32_t> faceGroups;
			for (uint32_t f = 0; f < faceCount; f++) {
				const uint32_t group = meshes[i]->faceGroupAt(f);
				bool exists = false;
				for (uint32_t g = 0; g < faceGroups.size(); g++) {
					if (faceGroups[g] == group) {
						exists = true;
						break;
					}
				}
				if (!exists)
					faceGroups.push_back(group);
			}
			// Create one chart group per face group.
			for (uint32_t g = 0; g < faceGroups.size(); g++) {
				ChartGroup *chartGroup = XA_NEW(ChartGroup, g, meshes[i], faceGroups[g]);
				if (!chartGroup->isVertexMap())
					chartGroup->computeCharts(options);
				m_chartGroups.push_back(chartGroup);
				m_chartGroupSourceMeshes.push_back(i);
				if (progressCallback) {
					const float groupProgess = (g + 1) / (float)faceGroups.size();
					const int newProgress = int(((i + groupProgess) / (float)meshCount) * 100.0f);
					if (newProgress != progress) {
						progress = newProgress;
						progressCallback(ProgressCategory::ComputingCharts, progress, progressCallbackUserData);
					}
				}
			}
		}
		if (progressCallback && progress != 100)
			progressCallback(ProgressCategory::ComputingCharts, 100, progressCallbackUserData);
	}

	void parameterizeCharts(ProgressCallback progressCallback, void *progressCallbackUserData)
	{
		int progress = 0;
		if (progressCallback)
			progressCallback(ProgressCategory::ParametizingCharts, 0, progressCallbackUserData);
		for (uint32_t i = 0; i < m_chartGroups.size(); i++) {
			if (!m_chartGroups[i]->isVertexMap())
				m_chartGroups[i]->parameterizeCharts();
			if (progressCallback) {
				const int newProgress = int((i + 1) / (float)m_chartGroups.size() * 100.0f);
				if (newProgress != progress) {
					progress = newProgress;
					progressCallback(ProgressCategory::ParametizingCharts, progress, progressCallbackUserData);
				}
			}
		}
		if (progressCallback && progress != 100)
			progressCallback(ProgressCategory::ParametizingCharts, 100, progressCallbackUserData);
	}

	void resetChartTexcoords()
	{
		const uint32_t nCharts = chartCount();
		if (m_originalChartTexcoords.isEmpty()) {
			// save
			m_originalChartTexcoords.resize(nCharts);
			for (uint32_t i = 0; i < nCharts; i++) {
				const Mesh *mesh = chartAt(i)->chartMesh();
				m_originalChartTexcoords[i].resize(mesh->vertexCount());
				for (uint32_t j = 0; j < mesh->vertexCount(); j++)
					m_originalChartTexcoords[i][j] = *mesh->texcoordAt(j);
			}
		} else {
			// restore
			for (uint32_t i = 0; i < nCharts; i++) {
				Mesh *mesh = chartAt(i)->chartMesh();
				for (uint32_t j = 0; j < mesh->vertexCount(); j++)
					*mesh->texcoordAt(j) = m_originalChartTexcoords[i][j];
			}
		}
	}

private:
	Array<ChartGroup *> m_chartGroups;
	Array<uint32_t> m_chartGroupSourceMeshes;
	Array<Array<Vector2> > m_originalChartTexcoords;
};

struct AtlasPacker
{
	AtlasPacker(Atlas *atlas) : m_atlas(atlas), m_width(0), m_height(0), m_texelsPerUnit(0)	{}

	~AtlasPacker()
	{
		for (uint32_t i = 0; i < m_bitmaps.size(); i++) {
			m_bitmaps[i]->~BitMap();
			XA_FREE(m_bitmaps[i]);
		}
	}

	uint32_t getWidth() const { return m_width; }
	uint32_t getHeight() const { return m_height; }
	uint32_t getNumAtlases() const { return m_bitmaps.size(); }
	float getTexelsPerUnit() const { return m_texelsPerUnit; }

	// Pack charts in the smallest possible rectangle.
	void packCharts(const PackerOptions &options, ProgressCallback progressCallback, void *progressCallbackUserData)
	{
		if (progressCallback)
			progressCallback(ProgressCategory::PackingCharts, 0, progressCallbackUserData);
		const uint32_t chartCount = m_atlas->chartCount();
		XA_PRINT(PrintFlags::PackingCharts, "   %u charts\n", chartCount);
		if (chartCount == 0) return;
		uint32_t resolution = options.resolution;
		m_texelsPerUnit = options.texelsPerUnit;
		if (resolution <= 0 || m_texelsPerUnit <= 0) {
			if (resolution <= 0 && m_texelsPerUnit <= 0)
				resolution = 1024;
			float meshArea = 0;
			for (uint32_t c = 0; c < chartCount; c++) {
				const Chart *chart = m_atlas->chartAt(c);
				if (!chart->isDisk())
					continue;
				meshArea += chart->computeSurfaceArea();
			}
			if (resolution <= 0) {
				// Estimate resolution based on the mesh surface area and given texel scale.
				const float texelCount = std::max(1.0f, meshArea * square(m_texelsPerUnit) / 0.75f); // Assume 75% utilization.
				resolution = nextPowerOfTwo(uint32_t(sqrtf(texelCount)));
				XA_PRINT(PrintFlags::PackingCharts, "      Estimating resolution as %d\n", resolution);
			}
			if (m_texelsPerUnit <= 0) {
				// Estimate a suitable texelsPerUnit to fit the given resolution.
				const float texelCount = std::max(1.0f, meshArea / 0.75f); // Assume 75% utilization.
				m_texelsPerUnit = sqrt((resolution * resolution) / texelCount);
				XA_PRINT(PrintFlags::PackingCharts, "      Estimating texelsPerUnit as %g\n", m_texelsPerUnit);
			}
		}
		m_rand = MTRand();
		Array<float> chartOrderArray;
		chartOrderArray.resize(chartCount);
		Array<Vector2> chartExtents;
		chartExtents.resize(chartCount);
		for (uint32_t c = 0; c < chartCount; c++) {
			Chart *chart = m_atlas->chartAt(c);
			if (!chart->isDisk()) {
				chartOrderArray[c] = 0;
				// Skip non-disks.
				continue;
			}
			Vector2 extents(0.0f);
			// Compute surface area to sort charts.
			float chartArea = chart->computeSurfaceArea();
			//chartOrderArray[c] = chartArea;
			// Compute chart scale
			float parametricArea = fabsf(chart->computeParametricArea());    // @@ There doesn't seem to be anything preventing parametric area to be negative.
			if (parametricArea < XA_EPSILON) {
				// When the parametric area is too small we use a rough approximation to prevent divisions by very small numbers.
				Vector2 bounds = chart->computeParametricBounds();
				parametricArea = bounds.x * bounds.y;
			}
			float scale = (chartArea / parametricArea) * m_texelsPerUnit;
			if (parametricArea == 0) { // < XA_EPSILON)
				scale = 0;
			}
			XA_ASSERT(std::isfinite(scale));
			// Compute bounding box of chart.
			Vector2 majorAxis, minorAxis, origin, end;
			computeBoundingBox(chart, &majorAxis, &minorAxis, &origin, &end);
			XA_ASSERT(isFinite(majorAxis) && isFinite(minorAxis) && isFinite(origin));
			// Sort charts by perimeter. @@ This is sometimes producing somewhat unexpected results. Is this right?
			//chartOrderArray[c] = ((end.x - origin.x) + (end.y - origin.y)) * scale;
			// Translate, rotate and scale vertices. Compute extents.
			Mesh *mesh = chart->chartMesh();
			const uint32_t vertexCount = mesh->vertexCount();
			for (uint32_t i = 0; i < vertexCount; i++) {
				Vector2 tmp;
				const Vector2 *texcoord = mesh->texcoordAt(i);
				tmp.x = dot(*texcoord, majorAxis);
				tmp.y = dot(*texcoord, minorAxis);
				tmp -= origin;
				tmp *= scale;
				if (tmp.x < 0 || tmp.y < 0) {
					XA_PRINT(PrintFlags::PackingCharts, "tmp: %f %f\n", tmp.x, tmp.y);
					XA_PRINT(PrintFlags::PackingCharts, "scale: %f\n", scale);
					XA_PRINT(PrintFlags::PackingCharts, "origin: %f %f\n", origin.x, origin.y);
					XA_PRINT(PrintFlags::PackingCharts, "majorAxis: %f %f\n", majorAxis.x, majorAxis.y);
					XA_PRINT(PrintFlags::PackingCharts, "minorAxis: %f %f\n", minorAxis.x, minorAxis.y);
					XA_DEBUG_ASSERT(false);
				}
				//XA_ASSERT(tmp.x >= 0 && tmp.y >= 0);
				XA_ASSERT(std::isfinite(tmp.x) && std::isfinite(tmp.y));
				*mesh->texcoordAt(i) = tmp;
				extents = max(extents, tmp);
			}
			XA_DEBUG_ASSERT(extents.x >= 0 && extents.y >= 0);
			// Limit chart size.
			if (extents.x > 1024 || extents.y > 1024) {
				float limit = std::max(extents.x, extents.y);
				scale = 1024 / (limit + 1);
				for (uint32_t i = 0; i < vertexCount; i++) {
					Vector2 *texcoord = mesh->texcoordAt(i);
					*texcoord *= scale;
				}
				extents *= scale;
				XA_DEBUG_ASSERT(extents.x <= 1024 && extents.y <= 1024);
			}
			// Scale the charts to use the entire texel area available. So, if the width is 0.1 we could scale it to 1 without increasing the lightmap usage and making a better
			// use of it. In many cases this also improves the look of the seams, since vertices on the chart boundaries have more chances of being aligned with the texel centers.
			float scale_x = 1.0f;
			float scale_y = 1.0f;
			float divide_x = 1.0f;
			float divide_y = 1.0f;
			if (extents.x > 0) {
				int cw = ftoi_ceil(extents.x);
				if (options.blockAlign && chart->isBlockAligned()) {
					// Align all chart extents to 4x4 blocks, but taking padding into account.
					if (options.conservative) {
						cw = align(cw + 2, 4) - 2;
					} else {
						cw = align(cw + 1, 4) - 1;
					}
				}
				scale_x = (float(cw) - XA_EPSILON);
				divide_x = extents.x;
				extents.x = float(cw);
			}
			if (extents.y > 0) {
				int ch = ftoi_ceil(extents.y);
				if (options.blockAlign && chart->isBlockAligned()) {
					// Align all chart extents to 4x4 blocks, but taking padding into account.
					if (options.conservative) {
						ch = align(ch + 2, 4) - 2;
					} else {
						ch = align(ch + 1, 4) - 1;
					}
				}
				scale_y = (float(ch) - XA_EPSILON);
				divide_y = extents.y;
				extents.y = float(ch);
			}
			for (uint32_t v = 0; v < vertexCount; v++) {
				Vector2 *texcoord = mesh->texcoordAt(v);
				texcoord->x /= divide_x;
				texcoord->y /= divide_y;
				texcoord->x *= scale_x;
				texcoord->y *= scale_y;
				XA_ASSERT(std::isfinite(texcoord->x) && std::isfinite(texcoord->y));
			}
			chartExtents[c] = extents;
			// Sort charts by perimeter.
			chartOrderArray[c] = extents.x + extents.y;
		}
		// @@ We can try to improve compression of small charts by sorting them by proximity like we do with vertex samples.
		// @@ How to do that? One idea: compute chart centroid, insert into grid, compute morton index of the cell, sort based on morton index.
		// @@ We would sort by morton index, first, then quantize the chart sizes, so that all small charts have the same size, and sort by size preserving the morton order.
		XA_PRINT(PrintFlags::PackingCharts, "      Sorting charts.\n");
		// Sort charts by area.
		m_radix = RadixSort();
		m_radix.sort(chartOrderArray);
		const uint32_t *ranks = m_radix.ranks();
		// Add sorted charts to bitmap.
		XA_PRINT(PrintFlags::PackingCharts, "      Rasterizing charts\n");
		int w = 0, h = 0;
		int progress = 0;
		for (uint32_t i = 0; i < chartCount; i++) {
			uint32_t c = ranks[chartCount - i - 1]; // largest chart first
			Chart *chart = m_atlas->chartAt(c);
			if (!chart->isDisk())
				continue;
			BitMap chart_bitmap;
			// @@ Add special cases for dot and line charts. @@ Lightmap rasterizer also needs to handle these special cases.
			// @@ We could also have a special case for chart quads. If the quad surface <= 4 texels, align vertices with texel centers and do not add padding. May be very useful for foliage.
			// @@ In general we could reduce the padding of all charts by one texel by using a rasterizer that takes into account the 2-texel footprint of the tent bilinear filter. For example,
			// if we have a chart that is less than 1 texel wide currently we add one texel to the left and one texel to the right creating a 3-texel-wide bitmap. However, if we know that the
			// chart is only 1 texel wide we could align it so that it only touches the footprint of two texels:
			//      |   |      <- Touches texels 0, 1 and 2.
			//    |   |        <- Only touches texels 0 and 1.
			// \   \ / \ /   /
			//  \   X   X   /
			//   \ / \ / \ /
			//    V   V   V
			//    0   1   2
			if (options.conservative) {
				// Init all bits to 0.
				chart_bitmap.resize(ftoi_ceil(chartExtents[c].x) + 1 + options.padding, ftoi_ceil(chartExtents[c].y) + 1 + options.padding, false);  // + 2 to add padding on both sides.
				// Rasterize chart and dilate.
				dchartBitmapDilate(chart, &chart_bitmap, options.padding);
			} else {
				// Init all bits to 0.
				chart_bitmap.resize(ftoi_ceil(chartExtents[c].x) + 1, ftoi_ceil(chartExtents[c].y) + 1, false);  // Add half a texels on each side.
				// Rasterize chart and dilate.
				dchartBitmap(chart, &chart_bitmap, Vector2(1), Vector2(0.5));
			}
			uint32_t currentBitmapIndex = 0;
			int best_x = 0, best_y = 0;
			int best_cw = 0, best_ch = 0;   // Includes padding now.
			int best_r = 0;
			for (;;)
			{
				bool firstChartInBitmap = false;
				if (currentBitmapIndex + 1 > m_bitmaps.size()) {
					// Chart doesn't fit in the current bitmap, create a new one.
					BitMap *bm = XA_NEW(BitMap);
					bm->clearAll();
					bm->resize(resolution, resolution, false);
					m_bitmaps.push_back(bm);
					firstChartInBitmap = true;
				}
				const bool foundLocation = findChartLocation(options.attempts, m_bitmaps[currentBitmapIndex], &chart_bitmap, chartExtents[c], w, h, &best_x, &best_y, &best_cw, &best_ch, &best_r, chart->isBlockAligned(), options.resolution <= 0);
				if (firstChartInBitmap && !foundLocation) {
					// Chart doesn't fit in an empty, newly allocated bitmap. texelsPerUnit must be too large for the resolution.
					XA_ASSERT(true && "chart doesn't fit");
					break;
				}
				if (options.resolution <= 0) {
					XA_DEBUG_ASSERT(foundLocation);
					break;
				}
				if (foundLocation)
					break;
				// Chart doesn't fit in the current bitmap, try the next one.
				currentBitmapIndex++;
			}
			/*if (w < best_x + best_cw || h < best_y + best_ch)
				XA_PRINT("Resize extents to (%d, %d).\n", best_x + best_cw, best_y + best_ch);
			*/
			// Update parametric extents.
			w = std::max(w, best_x + best_cw);
			h = std::max(h, best_y + best_ch);
			if (options.resolution <= 0) {
				// Resize bitmap if necessary.
				if (uint32_t(w) > m_bitmaps[0]->width() || uint32_t(h) > m_bitmaps[0]->height()) {
					m_bitmaps[0]->resize(nextPowerOfTwo(uint32_t(w)), nextPowerOfTwo(uint32_t(h)), false);
					XA_PRINT(PrintFlags::PackingCharts, "      Resize bitmap (%d, %d).\n", m_bitmaps[0]->width(), m_bitmaps[0]->height());
				}
			} else {
				w = std::min((int)options.resolution, w);
				h = std::min((int)options.resolution, h);
			}
			//XA_PRINT("Add chart at (%d, %d).\n", best_x, best_y);
			addChart(m_bitmaps[currentBitmapIndex], &chart_bitmap, w, h, best_x, best_y, best_r);
			chart->atlasIndex = (int32_t)currentBitmapIndex;
			//float best_angle = 2 * M_PI * best_r;
			// Translate and rotate chart texture coordinates.
			Mesh *mesh = chart->chartMesh();
			const uint32_t vertexCount = mesh->vertexCount();
			for (uint32_t v = 0; v < vertexCount; v++) {
				Vector2 *texcoord = mesh->texcoordAt(v);
				Vector2 t = *texcoord;
				if (best_r) std::swap(t.x, t.y);
				//vertex->tex.x = best_x + t.x * cosf(best_angle) - t.y * sinf(best_angle);
				//vertex->tex.y = best_y + t.x * sinf(best_angle) + t.y * cosf(best_angle);
				texcoord->x = best_x + t.x + 0.5f;
				texcoord->y = best_y + t.y + 0.5f;
				XA_ASSERT(texcoord->x >= 0 && texcoord->y >= 0);
				XA_ASSERT(std::isfinite(texcoord->x) && std::isfinite(texcoord->y));
			}
			if (progressCallback)
			{
				const int newProgress = int((i + 1) / (float)chartCount * 100.0f);
				if (newProgress != progress)
				{
					progress = newProgress;
					progressCallback(ProgressCategory::PackingCharts, progress, progressCallbackUserData);
				}
			}
		}
		//w -= padding - 1; // Leave one pixel border!
		//h -= padding - 1;
		m_width = std::max(0, w);
		m_height = std::max(0, h);
		if (options.resolution > 0)
			m_width = m_height = options.resolution;
		XA_PRINT(PrintFlags::PackingCharts, "      %dx%d resolution\n", m_width, m_height);
		if (progressCallback && progress != 100)
			progressCallback(ProgressCategory::PackingCharts, 0, progressCallbackUserData);
	}

	float computeAtlasUtilization(uint32_t atlasIndex) const
	{
		const uint32_t w = m_width;
		const uint32_t h = m_height;
		BitMap *bm = m_bitmaps[atlasIndex];
		XA_DEBUG_ASSERT(w <= bm->width());
		XA_DEBUG_ASSERT(h <= bm->height());
		uint32_t count = 0;
		for (uint32_t y = 0; y < h; y++) {
			for (uint32_t x = 0; x < w; x++) {
				count += bm->bitAt(x, y);
			}
		}
		return float(count) / (w * h);
	}

private:
	// IC: Brute force is slow, and random may take too much time to converge. We start inserting large charts in a small atlas. Using brute force is lame, because most of the space
	// is occupied at this point. At the end we have many small charts and a large atlas with sparse holes. Finding those holes randomly is slow. A better approach would be to
	// start stacking large charts as if they were tetris pieces. Once charts get small try to place them randomly. It may be interesting to try a intermediate strategy, first try
	// along one axis and then try exhaustively along that axis.
	bool findChartLocation(int attempts, const BitMap *atlasBitmap, const BitMap *chartBitmap, Vector2::Arg extents, int w, int h, int *best_x, int *best_y, int *best_w, int *best_h, int *best_r, bool blockAligned, bool resizableAtlas)
	{
		if (attempts <= 0 || attempts >= w * h)
			return findChartLocation_bruteForce(atlasBitmap, chartBitmap, extents, w, h, best_x, best_y, best_w, best_h, best_r, blockAligned, resizableAtlas);
		return findChartLocation_random(atlasBitmap, chartBitmap, extents, w, h, best_x, best_y, best_w, best_h, best_r, attempts, blockAligned, resizableAtlas);
	}

	bool findChartLocation_bruteForce(const BitMap *atlasBitmap, const BitMap *chartBitmap, Vector2::Arg /*extents*/, int w, int h, int *best_x, int *best_y, int *best_w, int *best_h, int *best_r, bool blockAligned, bool resizableAtlas)
	{
		bool result = false;
		const int BLOCK_SIZE = 4;
		int best_metric = INT_MAX;
		int step_size = blockAligned ? BLOCK_SIZE : 1;
		// Try two different orientations.
		for (int r = 0; r < 2; r++) {
			int cw = chartBitmap->width();
			int ch = chartBitmap->height();
			if (r & 1) std::swap(cw, ch);
			for (int y = 0; y <= h + 1; y += step_size) { // + 1 to extend atlas in case atlas full.
				for (int x = 0; x <= w + 1; x += step_size) { // + 1 not really necessary here.
					if (!resizableAtlas && (x > (int)atlasBitmap->width() - cw || y > (int)atlasBitmap->height() - ch))
						continue;
					// Early out.
					int area = std::max(w, x + cw) * std::max(h, y + ch);
					//int perimeter = max(w, x+cw) + max(h, y+ch);
					int extents = std::max(std::max(w, x + cw), std::max(h, y + ch));
					int metric = extents * extents + area;
					if (metric > best_metric) {
						continue;
					}
					if (metric == best_metric && std::max(x, y) >= std::max(*best_x, *best_y)) {
						// If metric is the same, pick the one closest to the origin.
						continue;
					}
					if (canAddChart(atlasBitmap, chartBitmap, w, h, x, y, r)) {
						result = true;
						best_metric = metric;
						*best_x = x;
						*best_y = y;
						*best_w = cw;
						*best_h = ch;
						*best_r = r;
						if (area == w * h) {
							// Chart is completely inside, do not look at any other location.
							goto done;
						}
					}
				}
			}
		}
	done:
		XA_DEBUG_ASSERT (best_metric != INT_MAX);
		return result;
	}

	bool findChartLocation_random(const BitMap *atlasBitmap, const BitMap *chartBitmap, Vector2::Arg /*extents*/, int w, int h, int *best_x, int *best_y, int *best_w, int *best_h, int *best_r, int minTrialCount, bool blockAligned, bool resizableAtlas)
	{
		bool result = false;
		const int BLOCK_SIZE = 4;
		int best_metric = INT_MAX;
		for (int i = 0; i < minTrialCount; i++) {
			int cw = chartBitmap->width();
			int ch = chartBitmap->height();
			int r = m_rand.getRange(1);
			if (r & 1)
				std::swap(cw, ch);
			// + 1 to extend atlas in case atlas full. We may want to use a higher number to increase probability of extending atlas.
			int xRange = w + 1;
			int yRange = h + 1;
			if (!resizableAtlas) {
				xRange = std::min(xRange, (int)atlasBitmap->width() - cw);
				yRange = std::min(yRange, (int)atlasBitmap->height() - ch);
			}
			int x = m_rand.getRange(xRange);
			int y = m_rand.getRange(yRange);
			if (blockAligned) {
				x = align(x, BLOCK_SIZE);
				y = align(y, BLOCK_SIZE);
				if (!resizableAtlas && (x > (int)atlasBitmap->width() - cw || y > (int)atlasBitmap->height() - ch))
					continue; // Block alignment pushed the chart outside the atlas.
			}
			// Early out.
			int area = std::max(w, x + cw) * std::max(h, y + ch);
			//int perimeter = max(w, x+cw) + max(h, y+ch);
			int extents = std::max(std::max(w, x + cw), std::max(h, y + ch));
			int metric = extents * extents + area;
			if (metric > best_metric) {
				continue;
			}
			if (metric == best_metric && std::min(x, y) > std::min(*best_x, *best_y)) {
				// If metric is the same, pick the one closest to the origin.
				continue;
			}
			if (canAddChart(atlasBitmap, chartBitmap, w, h, x, y, r)) {
				result = true;
				best_metric = metric;
				*best_x = x;
				*best_y = y;
				*best_w = cw;
				*best_h = ch;
				*best_r = r;
				if (area == w * h) {
					// Chart is completely inside, do not look at any other location.
					break;
				}
			}
		}
		return result;
	}

	void dchartBitmapDilate(const Chart *chart, BitMap *bitmap, int padding)
	{
		const int w = bitmap->width();
		const int h = bitmap->height();
		const Vector2 extents = Vector2(float(w), float(h));
		// Rasterize chart faces, check that all bits are not set.
		const uint32_t faceCount = chart->faceCount();
		for (uint32_t f = 0; f < faceCount; f++) {
			const Mesh *mesh = chart->chartMesh();
			Vector2 vertices[3];
			uint32_t edgeCount = 0;
			for (Mesh::FaceEdgeIterator it(mesh, f); !it.isDone(); it.advance()) {
				if (edgeCount < 3)
					vertices[edgeCount] = it.texcoord0() + Vector2(0.5f) + Vector2(float(padding), float(padding));
				edgeCount++;
			}
			XA_DEBUG_ASSERT(edgeCount == 3);
			raster::drawTriangle(raster::Mode_Antialiased, extents, true, vertices, AtlasPacker::setBitsCallback, bitmap);
		}
		// Expand chart by padding pixels. (dilation)
		BitMap tmp(w, h);
		for (int i = 0; i < padding; i++) {
			tmp.clearAll();
			for (int y = 0; y < h; y++) {
				for (int x = 0; x < w; x++) {
					bool b = bitmap->bitAt(x, y);
					if (!b) {
						if (x > 0) {
							b |= bitmap->bitAt(x - 1, y);
							if (y > 0) b |= bitmap->bitAt(x - 1, y - 1);
							if (y < h - 1) b |= bitmap->bitAt(x - 1, y + 1);
						}
						if (y > 0) b |= bitmap->bitAt(x, y - 1);
						if (y < h - 1) b |= bitmap->bitAt(x, y + 1);
						if (x < w - 1) {
							b |= bitmap->bitAt(x + 1, y);
							if (y > 0) b |= bitmap->bitAt(x + 1, y - 1);
							if (y < h - 1) b |= bitmap->bitAt(x + 1, y + 1);
						}
					}
					if (b) tmp.setBitAt(x, y);
				}
			}
			std::swap(tmp, *bitmap);
		}
	}

	void dchartBitmap(const Chart *chart, BitMap *bitmap, const Vector2 &scale, const Vector2 &offset)
	{
		const int w = bitmap->width();
		const int h = bitmap->height();
		const Vector2 extents = Vector2(float(w), float(h));
		static const Vector2 pad[4] = {
			Vector2(-0.5, -0.5),
			Vector2(0.5, -0.5),
			Vector2(-0.5, 0.5),
			Vector2(0.5, 0.5)
		};
		// Rasterize 4 times to add proper padding.
		for (int i = 0; i < 4; i++) {
			// Rasterize chart faces, check that all bits are not set.
			const uint32_t faceCount = chart->chartMesh()->faceCount();
			for (uint32_t f = 0; f < faceCount; f++) {
				const Mesh *mesh = chart->chartMesh();
				Vector2 vertices[3];
				uint32_t edgeCount = 0;
				for (Mesh::FaceEdgeIterator it(mesh, f); !it.isDone(); it.advance()) {
					if (edgeCount < 3) {
						vertices[edgeCount] = it.texcoord0() * scale + offset + pad[i];
						XA_ASSERT(ftoi_ceil(vertices[edgeCount].x) >= 0);
						XA_ASSERT(ftoi_ceil(vertices[edgeCount].y) >= 0);
						XA_ASSERT(ftoi_ceil(vertices[edgeCount].x) <= w);
						XA_ASSERT(ftoi_ceil(vertices[edgeCount].y) <= h);
					}
					edgeCount++;
				}
				XA_ASSERT(edgeCount == 3);
				raster::drawTriangle(raster::Mode_Antialiased, extents, /*enableScissors=*/true, vertices, AtlasPacker::setBitsCallback, bitmap);
			}
		}
		// Expand chart by padding pixels. (dilation)
		BitMap tmp(w, h);
		tmp.clearAll();
		for (int y = 0; y < h; y++) {
			for (int x = 0; x < w; x++) {
				bool b = bitmap->bitAt(x, y);
				if (!b) {
					if (x > 0) {
						b |= bitmap->bitAt(x - 1, y);
						if (y > 0) b |= bitmap->bitAt(x - 1, y - 1);
						if (y < h - 1) b |= bitmap->bitAt(x - 1, y + 1);
					}
					if (y > 0) b |= bitmap->bitAt(x, y - 1);
					if (y < h - 1) b |= bitmap->bitAt(x, y + 1);
					if (x < w - 1) {
						b |= bitmap->bitAt(x + 1, y);
						if (y > 0) b |= bitmap->bitAt(x + 1, y - 1);
						if (y < h - 1) b |= bitmap->bitAt(x + 1, y + 1);
					}
				}
				if (b) tmp.setBitAt(x, y);
			}
		}
		std::swap(tmp, *bitmap);
	}

	bool canAddChart(const BitMap *atlasBitmap, const BitMap *chartBitmap, int atlas_w, int atlas_h, int offset_x, int offset_y, int r)
	{
		XA_DEBUG_ASSERT(r == 0 || r == 1);
		// Check whether the two bitmaps overlap.
		const int w = chartBitmap->width();
		const int h = chartBitmap->height();
		if (r == 0) {
			for (int y = 0; y < h; y++) {
				int yy = y + offset_y;
				if (yy >= 0) {
					for (int x = 0; x < w; x++) {
						int xx = x + offset_x;
						if (xx >= 0) {
							if (chartBitmap->bitAt(x, y)) {
								if (xx < atlas_w && yy < atlas_h) {
									if (atlasBitmap->bitAt(xx, yy))
										return false;
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
							if (chartBitmap->bitAt(x, y)) {
								if (xx < atlas_w && yy < atlas_h) {
									if (atlasBitmap->bitAt(xx, yy))
										return false;
								}
							}
						}
					}
				}
			}
		}
		return true;
	}

	void addChart(BitMap *atlasBitmap, const BitMap *chartBitmap, int atlas_w, int atlas_h, int offset_x, int offset_y, int r)
	{
		XA_DEBUG_ASSERT(r == 0 || r == 1);
		// Check whether the two bitmaps overlap.
		const int w = chartBitmap->width();
		const int h = chartBitmap->height();
		if (r == 0) {
			for (int y = 0; y < h; y++) {
				int yy = y + offset_y;
				if (yy >= 0) {
					for (int x = 0; x < w; x++) {
						int xx = x + offset_x;
						if (xx >= 0) {
							if (chartBitmap->bitAt(x, y)) {
								if (xx < atlas_w && yy < atlas_h) {
									XA_DEBUG_ASSERT(atlasBitmap->bitAt(xx, yy) == false);
									atlasBitmap->setBitAt(xx, yy);
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
							if (chartBitmap->bitAt(x, y)) {
								if (xx < atlas_w && yy < atlas_h) {
									XA_DEBUG_ASSERT(atlasBitmap->bitAt(xx, yy) == false);
									atlasBitmap->setBitAt(xx, yy);
								}
							}
						}
					}
				}
			}
		}
	}

	static bool setBitsCallback(void *param, int x, int y, Vector3::Arg, Vector3::Arg, Vector3::Arg, float area)
	{
		BitMap *bitmap = (BitMap * )param;
		if (area > 0.0) {
			bitmap->setBitAt(x, y);
		}
		return true;
	}

	// Compute the convex hull using Graham Scan.
	static void convexHull(const Array<Vector2> &input, Array<Vector2> &output, float epsilon)
	{
		const uint32_t inputCount = input.size();
		Array<float> coords;
		coords.resize(inputCount);
		for (uint32_t i = 0; i < inputCount; i++) {
			coords[i] = input[i].x;
		}
		RadixSort radix;
		radix.sort(coords);
		const uint32_t *ranks = radix.ranks();
		Array<Vector2> top;
		top.reserve(inputCount);
		Array<Vector2> bottom;
		bottom.reserve(inputCount);
		Vector2 P = input[ranks[0]];
		Vector2 Q = input[ranks[inputCount - 1]];
		float topy = std::max(P.y, Q.y);
		float boty = std::min(P.y, Q.y);
		for (uint32_t i = 0; i < inputCount; i++) {
			Vector2 p = input[ranks[i]];
			if (p.y >= boty) top.push_back(p);
		}
		for (uint32_t i = 0; i < inputCount; i++) {
			Vector2 p = input[ranks[inputCount - 1 - i]];
			if (p.y <= topy) bottom.push_back(p);
		}
		// Filter top list.
		output.clear();
		output.push_back(top[0]);
		output.push_back(top[1]);
		for (uint32_t i = 2; i < top.size(); ) {
			Vector2 a = output[output.size() - 2];
			Vector2 b = output[output.size() - 1];
			Vector2 c = top[i];
			float area = triangleArea(a, b, c);
			if (area >= -epsilon) {
				output.pop_back();
			}
			if (area < -epsilon || output.size() == 1) {
				output.push_back(c);
				i++;
			}
		}
		uint32_t top_count = output.size();
		output.push_back(bottom[1]);
		// Filter bottom list.
		for (uint32_t i = 2; i < bottom.size(); ) {
			Vector2 a = output[output.size() - 2];
			Vector2 b = output[output.size() - 1];
			Vector2 c = bottom[i];
			float area = triangleArea(a, b, c);
			if (area >= -epsilon) {
				output.pop_back();
			}
			if (area < -epsilon || output.size() == top_count) {
				output.push_back(c);
				i++;
			}
		}
		// Remove duplicate element.
		XA_DEBUG_ASSERT(output.front() == output.back());
		output.pop_back();
	}

	// This should compute convex hull and use rotating calipers to find the best box. Currently it uses a brute force method.
	static void computeBoundingBox(const Chart *chart, Vector2 *majorAxis, Vector2 *minorAxis, Vector2 *minCorner, Vector2 *maxCorner)
	{
		// Compute list of boundary points.
		Array<Vector2> points;
		points.reserve(16);
		const Mesh *mesh = chart->chartMesh();
		const uint32_t vertexCount = mesh->vertexCount();
		uint32_t bp = 0;
		for (uint32_t v = 0; v < vertexCount; v++) {
			if (mesh->isBoundaryVertex(v)) {
				points.push_back(*mesh->texcoordAt(v));
				bp++;
			}
		}
		XA_DEBUG_ASSERT(points.size() > 0);
		Array<Vector2> hull;
		convexHull(points, hull, 0.00001f);
		// @@ Ideally I should use rotating calipers to find the best box. Using brute force for now.
		float best_area = FLT_MAX;
		Vector2 best_min(0);
		Vector2 best_max(0);
		Vector2 best_axis(0);
		const uint32_t hullCount = hull.size();
		for (uint32_t i = 0, j = hullCount - 1; i < hullCount; j = i, i++) {
			if (equal(hull[i], hull[j])) {
				continue;
			}
			Vector2 axis = normalize(hull[i] - hull[j], 0.0f);
			XA_DEBUG_ASSERT(isFinite(axis));
			// Compute bounding box.
			Vector2 box_min(FLT_MAX, FLT_MAX);
			Vector2 box_max(-FLT_MAX, -FLT_MAX);
			for (uint32_t v = 0; v < hullCount; v++) {
				Vector2 point = hull[v];
				float x = dot(axis, point);
				if (x < box_min.x) box_min.x = x;
				if (x > box_max.x) box_max.x = x;
				float y = dot(Vector2(-axis.y, axis.x), point);
				if (y < box_min.y) box_min.y = y;
				if (y > box_max.y) box_max.y = y;
			}
			// Compute box area.
			float area = (box_max.x - box_min.x) * (box_max.y - box_min.y);
			if (area < best_area) {
				best_area = area;
				best_min = box_min;
				best_max = box_max;
				best_axis = axis;
			}
		}
		// Consider all points, not only boundary points, in case the input chart is malformed.
		for (uint32_t i = 0; i < vertexCount; i++) {
			Vector2 point = *mesh->texcoordAt(i);
			float x = dot(best_axis, point);
			if (x < best_min.x) best_min.x = x;
			if (x > best_max.x) best_max.x = x;
			float y = dot(Vector2(-best_axis.y, best_axis.x), point);
			if (y < best_min.y) best_min.y = y;
			if (y > best_max.y) best_max.y = y;
		}
		*majorAxis = best_axis;
		*minorAxis = Vector2(-best_axis.y, best_axis.x);
		*minCorner = best_min;
		*maxCorner = best_max;
	}

	Atlas *m_atlas;
	Array<BitMap *> m_bitmaps;
	RadixSort m_radix;
	uint32_t m_width;
	uint32_t m_height;
	float m_texelsPerUnit;
	MTRand m_rand;
};

} // namespace param
} // namespace internal

struct Context
{
	Atlas atlas;
	internal::param::Atlas paramAtlas;
	internal::Array<internal::Mesh *> meshes;
};

Atlas *Create()
{
	Context *ctx = XA_NEW(Context);
	ctx->atlas.atlasCount = 0;
	ctx->atlas.chartCount = 0;
	ctx->atlas.height = 0;
	ctx->atlas.meshCount = 0;
	ctx->atlas.meshes = NULL;
	ctx->atlas.texelsPerUnit = 0;
	ctx->atlas.utilization = NULL;
	ctx->atlas.width = 0;
	return &ctx->atlas;
}

static void DestroyOutputMeshes(Context *ctx)
{
	if (!ctx->atlas.meshes)
		return;
	for (int i = 0; i < (int)ctx->atlas.meshCount; i++) {
		Mesh *mesh = ctx->atlas.meshes[i];
		for (uint32_t j = 0; j < mesh->chartCount; j++)
			XA_FREE(mesh->chartArray[j].indexArray);
		XA_FREE(mesh->chartArray);
		XA_FREE(mesh->vertexArray);
		XA_FREE(mesh->indexArray);
		XA_FREE(mesh);
	}
	XA_FREE(ctx->atlas.meshes);
	ctx->atlas.meshes = NULL;
}

void Destroy(Atlas *atlas)
{
	XA_DEBUG_ASSERT(atlas);
	Context *ctx = (Context *)atlas;
	if (atlas->utilization)
		XA_FREE(atlas->utilization);
	for (int i = 0; i < (int)ctx->meshes.size(); i++) {
		ctx->meshes[i]->~Mesh();
		XA_FREE(ctx->meshes[i]);
	}
	DestroyOutputMeshes(ctx);
	ctx->~Context();
	XA_FREE(ctx);
#if XA_DEBUG_HEAP
	internal::ReportAllocs();
#endif
}

static internal::Vector3 DecodePosition(const MeshDecl &meshDecl, uint32_t index)
{
	XA_DEBUG_ASSERT(meshDecl.vertexPositionData);
	XA_DEBUG_ASSERT(meshDecl.vertexPositionStride > 0);
	return *((const internal::Vector3 *)&((const uint8_t *)meshDecl.vertexPositionData)[meshDecl.vertexPositionStride * index]);
}

static internal::Vector3 DecodeNormal(const MeshDecl &meshDecl, uint32_t index)
{
	XA_DEBUG_ASSERT(meshDecl.vertexNormalData);
	XA_DEBUG_ASSERT(meshDecl.vertexNormalStride > 0);
	return *((const internal::Vector3 *)&((const uint8_t *)meshDecl.vertexNormalData)[meshDecl.vertexNormalStride * index]);
}

static internal::Vector2 DecodeUv(const MeshDecl &meshDecl, uint32_t index)
{
	XA_DEBUG_ASSERT(meshDecl.vertexUvData);
	XA_DEBUG_ASSERT(meshDecl.vertexUvStride > 0);
	return *((const internal::Vector2 *)&((const uint8_t *)meshDecl.vertexUvData)[meshDecl.vertexUvStride * index]);
}

static uint32_t DecodeIndex(IndexFormat::Enum format, const void *indexData, int32_t offset, uint32_t i)
{
	XA_DEBUG_ASSERT(indexData);
	if (format == IndexFormat::UInt16)
		return uint16_t((int32_t)((const uint16_t *)indexData)[i] + offset);
	return uint32_t((int32_t)((const uint32_t *)indexData)[i] + offset);
}

static float EdgeLength(internal::Vector3 pos1, internal::Vector3 pos2)
{
	return internal::length(pos2 - pos1);
}

AddMeshError::Enum AddMesh(Atlas *atlas, const MeshDecl &meshDecl)
{
	XA_DEBUG_ASSERT(atlas);
	XA_DEBUG_ASSERT(meshDecl.vertexCount > 0);
	XA_DEBUG_ASSERT(meshDecl.indexCount > 0);
	Context *ctx = (Context *)atlas;
	XA_PRINT(PrintFlags::MeshCreation, "Adding mesh %d: %u vertices, %u triangles\n", atlas->meshCount, meshDecl.vertexCount, meshDecl.indexCount / 3);
	// Expecting triangle faces.
	if ((meshDecl.indexCount % 3) != 0)
		return AddMeshError::InvalidIndexCount;
	// Check if any index is out of range.
	for (uint32_t i = 0; i < meshDecl.indexCount; i++) {
		const uint32_t index = DecodeIndex(meshDecl.indexFormat, meshDecl.indexData, meshDecl.indexOffset, i);
		if (index >= meshDecl.vertexCount)
			return AddMeshError::IndexOutOfRange;
	}
	internal::Mesh *mesh = XA_NEW(internal::Mesh, meshDecl.vertexCount, meshDecl.indexCount / 3, ctx->meshes.size());
	for (uint32_t i = 0; i < meshDecl.vertexCount; i++) {
		internal::Vector3 normal(0);
		internal::Vector2 texcoord(0);
		if (meshDecl.vertexNormalData)
			normal = DecodeNormal(meshDecl, i);
		if (meshDecl.vertexUvData)
			texcoord = DecodeUv(meshDecl, i);
		mesh->addVertex(DecodePosition(meshDecl, i), normal, texcoord);
	}
	mesh->createColocals();
	for (uint32_t i = 0; i < meshDecl.indexCount / 3; i++) {
		uint32_t tri[3];
		for (int j = 0; j < 3; j++)
			tri[j] = DecodeIndex(meshDecl.indexFormat, meshDecl.indexData, meshDecl.indexOffset, i * 3 + j);
		uint32_t faceFlags = 0;
		// Check for degenerate or zero length edges.
		for (int j = 0; j < 3; j++) {
			const uint32_t edges[6] = { 0, 1, 1, 2, 2, 0 };
			const uint32_t index1 = tri[edges[j * 2 + 0]];
			const uint32_t index2 = tri[edges[j * 2 + 1]];
			if (index1 == index2) {
				faceFlags |= internal::FaceFlags::Ignore;
				XA_PRINT(PrintFlags::MeshWarnings, "Mesh %d degenerate edge: index %d, index %d\n", (int)atlas->meshCount, index1, index2);
				break;
			}
			const internal::Vector3 pos1 = DecodePosition(meshDecl, index1);
			const internal::Vector3 pos2 = DecodePosition(meshDecl, index2);
			if (EdgeLength(pos1, pos2) <= 0.0f) {
				faceFlags |= internal::FaceFlags::Ignore;
				XA_PRINT(PrintFlags::MeshWarnings, "Mesh %d zero length edge: index %d position (%g %g %g), index %d position (%g %g %g)\n", (int)atlas->meshCount, index1, pos1.x, pos1.y, pos1.z, index2, pos2.x, pos2.y, pos2.z);
				break;
			}
		}
		// Check for zero area faces. Don't bother if a degenerate or zero length edge was already detected.
		if (!(faceFlags & internal::FaceFlags::Ignore))
		{
			const internal::Vector3 a = DecodePosition(meshDecl, tri[0]);
			const internal::Vector3 b = DecodePosition(meshDecl, tri[1]);
			const internal::Vector3 c = DecodePosition(meshDecl, tri[2]);
			const float area = internal::length(internal::cross(b - a, c - a)) * 0.5f;
			if (area <= 0.0f)
			{
				faceFlags |= internal::FaceFlags::Ignore;
				XA_PRINT(PrintFlags::MeshWarnings, "Mesh %d zero area face: %d, indices (%d %d %d)\n", (int)atlas->meshCount, i, tri[0], tri[1], tri[2]);
			}
		}
		if (meshDecl.faceIgnoreData && meshDecl.faceIgnoreData[i])
			faceFlags |= internal::FaceFlags::Ignore;
		mesh->addFace(tri[0], tri[1], tri[2], faceFlags);
	}
	mesh->createFaceGroups();
	mesh->createBoundaries();
#if XA_DEBUG_EXPORT_OBJ
	mesh->writeObj();
#endif
	ctx->meshes.push_back(mesh);
	atlas->meshCount++;
	return AddMeshError::Success;
}

void GenerateCharts(Atlas *atlas, CharterOptions charterOptions, ProgressCallback progressCallback, void *progressCallbackUserData)
{
	XA_DEBUG_ASSERT(atlas);
	Context *ctx = (Context *)atlas;
	if (ctx->meshes.isEmpty())
		return;
	atlas->atlasCount = 0;
	atlas->chartCount = 0;
	atlas->height = 0;
	atlas->texelsPerUnit = 0;
	atlas->width = 0;
	if (atlas->utilization) {
		XA_FREE(atlas->utilization);
		atlas->utilization = NULL;
	}
	DestroyOutputMeshes(ctx);
	// Chart meshes.
	XA_PRINT(PrintFlags::ComputingCharts, "Computing charts\n");
	ctx->paramAtlas.computeCharts(ctx->meshes, charterOptions, progressCallback, progressCallbackUserData);
	XA_PRINT(PrintFlags::ParametizingCharts, "Parameterizing charts\n");
	ctx->paramAtlas.parameterizeCharts(progressCallback, progressCallbackUserData);
	// Count charts.
	for (uint32_t i = 0; i < ctx->meshes.size(); i++) {
		for (uint32_t j = 0; j < ctx->paramAtlas.chartGroupCount(i); j++) {
			const internal::param::ChartGroup *chartGroup = ctx->paramAtlas.chartGroupAt(i, j);
			if (!chartGroup->isVertexMap())
				atlas->chartCount += chartGroup->chartCount();
		}
	}
}

void PackCharts(Atlas *atlas, PackerOptions packerOptions, ProgressCallback progressCallback, void *progressCallbackUserData)
{
	XA_DEBUG_ASSERT(atlas);
	Context *ctx = (Context *)atlas;
	atlas->atlasCount = 0;
	atlas->height = 0;
	atlas->texelsPerUnit = packerOptions.texelsPerUnit;
	atlas->width = 0;
	DestroyOutputMeshes(ctx);
	if (atlas->utilization) {
		XA_FREE(atlas->utilization);
		atlas->utilization = NULL;
	}
	if (atlas->chartCount <= 0)
		return;
	XA_PRINT(PrintFlags::PackingCharts, "Packing charts\n");
	ctx->paramAtlas.resetChartTexcoords();
	internal::param::AtlasPacker packer(&ctx->paramAtlas);
	packer.packCharts(packerOptions, progressCallback, progressCallbackUserData);
	atlas->atlasCount = packer.getNumAtlases();
	atlas->width = packer.getWidth();
	atlas->height = packer.getHeight();
	atlas->texelsPerUnit = packer.getTexelsPerUnit();
	atlas->utilization = XA_ALLOC_ARRAY(float, atlas->atlasCount);
	for (uint32_t i = 0; i < atlas->atlasCount; i++)
		atlas->utilization[i] = packer.computeAtlasUtilization(i);
	XA_PRINT(PrintFlags::BuildingOutputMeshes, "Building output meshes\n");
	int progress = 0;
	if (progressCallback)
		progressCallback(ProgressCategory::BuildingOutputMeshes, 0, progressCallbackUserData);
	atlas->meshes = XA_ALLOC_ARRAY(Mesh *, atlas->meshCount);
	for (uint32_t i = 0; i < atlas->meshCount; i++) {
		Mesh *outputMesh = atlas->meshes[i] = XA_ALLOC(Mesh);
		outputMesh->vertexCount = 0;
		outputMesh->indexCount = 0;
		for (uint32_t cg = 0; cg < ctx->paramAtlas.chartGroupCount(i); cg++) {
			const internal::param::ChartGroup *chartGroup = ctx->paramAtlas.chartGroupAt(i, cg);
			if (chartGroup->isVertexMap()) {
				outputMesh->vertexCount += chartGroup->mesh()->vertexCount();
				outputMesh->indexCount += chartGroup->mesh()->faceCount() * 3;
			} else {
				for (uint32_t c = 0; c < chartGroup->chartCount(); c++) {
					const internal::param::Chart *chart = chartGroup->chartAt(c);
					outputMesh->vertexCount += chart->vertexCount();
					outputMesh->indexCount += chart->faceCount() * 3;
				}
			}
		}
		// Vertices.
		outputMesh->vertexArray = XA_ALLOC_ARRAY(Vertex, outputMesh->vertexCount);
		uint32_t vertexOffset = 0;
		for (uint32_t cg = 0; cg < ctx->paramAtlas.chartGroupCount(i); cg++) {
			const internal::param::ChartGroup *chartGroup = ctx->paramAtlas.chartGroupAt(i, cg);
			if (chartGroup->isVertexMap()) {
				const internal::Mesh *mesh = chartGroup->mesh();
				for (uint32_t v = 0; v < mesh->vertexCount(); v++) {
					Vertex &vertex = outputMesh->vertexArray[vertexOffset++];
					vertex.atlasIndex = -1;
					const internal::Vector2 &uv = *mesh->texcoordAt(v);
					vertex.uv[0] = std::max(0.0f, uv.x);
					vertex.uv[1] = std::max(0.0f, uv.y);
					vertex.xref = chartGroup->mapVertexToSourceVertex(v);
				}
			} else {
				for (uint32_t c = 0; c < chartGroup->chartCount(); c++) {
					const internal::param::Chart *chart = chartGroup->chartAt(c);
					for (uint32_t v = 0; v < chart->vertexCount(); v++) {
						Vertex &vertex = outputMesh->vertexArray[vertexOffset++];
						XA_DEBUG_ASSERT(chart->atlasIndex >= 0);
						vertex.atlasIndex = chart->atlasIndex;
						const internal::Vector2 &uv = *chart->chartMesh()->texcoordAt(v);
						vertex.uv[0] = std::max(0.0f, uv.x);
						vertex.uv[1] = std::max(0.0f, uv.y);
						vertex.xref = chartGroup->mapVertexToSourceVertex(chart->mapChartVertexToOriginalVertex(v));
					}
				}
			}
		}
		// Indices.
		outputMesh->indexArray = XA_ALLOC_ARRAY(uint32_t, outputMesh->indexCount);
		vertexOffset = 0;
		for (uint32_t cg = 0; cg < ctx->paramAtlas.chartGroupCount(i); cg++) {
			const internal::param::ChartGroup *chartGroup = ctx->paramAtlas.chartGroupAt(i, cg);
			if (chartGroup->isVertexMap()) {
				const internal::Mesh *mesh = chartGroup->mesh();
				for (uint32_t f = 0; f < mesh->faceCount(); f++) {
					const internal::Face *face = mesh->faceAt(f);
					uint32_t indexOffset = chartGroup->mapFaceToSourceFace(f) * 3;
					for (uint32_t j = 0; j < 3; j++)
						outputMesh->indexArray[indexOffset++] = vertexOffset + mesh->vertexAt(face->firstIndex + j);
				}
			} else {
				for (uint32_t c = 0; c < chartGroup->chartCount(); c++) {
					const internal::param::Chart *chart = chartGroup->chartAt(c);
					const internal::Mesh *mesh = chart->chartMesh();
					for (uint32_t f = 0; f < chart->faceCount(); f++) {
						const internal::Face *face = mesh->faceAt(f);
						uint32_t indexOffset = chartGroup->mapFaceToSourceFace(chart->mapFaceToSourceFace(f)) * 3;
						for (uint32_t j = 0; j < 3; j++)
							outputMesh->indexArray[indexOffset++] = vertexOffset + mesh->vertexAt(face->firstIndex + j);
					}
					vertexOffset += chart->vertexCount();
				}
			}
		}
		// Charts.
		// Ignore vertex mapped charts.
		outputMesh->chartCount = 0;
		for (uint32_t j = 0; j < ctx->paramAtlas.chartGroupCount(i); j++) {
			const internal::param::ChartGroup *chartGroup = ctx->paramAtlas.chartGroupAt(i, j);
			if (!chartGroup->isVertexMap())
				outputMesh->chartCount += chartGroup->chartCount();
		}
		outputMesh->chartArray = XA_ALLOC_ARRAY(Chart, outputMesh->chartCount);
		vertexOffset = 0;
		uint32_t chartIndex = 0;
		for (uint32_t cg = 0; cg < ctx->paramAtlas.chartGroupCount(i); cg++) {
			const internal::param::ChartGroup *chartGroup = ctx->paramAtlas.chartGroupAt(i, cg);
			for (uint32_t c = 0; c < chartGroup->chartCount(); c++) {
				const internal::param::Chart *chart = chartGroup->chartAt(c);
				Chart *outputChart = &outputMesh->chartArray[chartIndex];
				XA_DEBUG_ASSERT(chart->atlasIndex >= 0);
				outputChart->atlasIndex = (uint32_t)chart->atlasIndex;
				const internal::Mesh *mesh = chart->chartMesh();
				outputChart->indexCount = mesh->faceCount() * 3;
				outputChart->indexArray = XA_ALLOC_ARRAY(uint32_t, outputChart->indexCount);
				for (uint32_t k = 0; k < mesh->faceCount(); k++) {
					const internal::Face *face = mesh->faceAt(k);
					outputChart->indexArray[3 * k + 0] = vertexOffset + mesh->vertexAt(face->firstIndex + 0);
					outputChart->indexArray[3 * k + 1] = vertexOffset + mesh->vertexAt(face->firstIndex + 1);
					outputChart->indexArray[3 * k + 2] = vertexOffset + mesh->vertexAt(face->firstIndex + 2);
				}
				vertexOffset += chart->vertexCount();
				chartIndex++;
			}
		}
		XA_PRINT(PrintFlags::BuildingOutputMeshes, "   mesh %u: %u vertices, %u triangles, %u charts\n", i, outputMesh->vertexCount, outputMesh->indexCount / 3, outputMesh->chartCount);
		if (progressCallback)
		{
			const int newProgress = int((i + 1) / (float)atlas->meshCount * 100.0f);
			if (newProgress != progress)
			{
				progress = newProgress;
				progressCallback(ProgressCategory::BuildingOutputMeshes, progress, progressCallbackUserData);
			}
		}
	}
	if (progressCallback && progress != 100)
		progressCallback(ProgressCategory::BuildingOutputMeshes, 0, progressCallbackUserData);
}

void SetRealloc(ReallocFunc reallocFunc)
{
	internal::s_realloc = reallocFunc;
}

void SetPrint(int flags, PrintFunc print)
{
	internal::s_printFlags = flags;
	internal::s_print = print ? print : printf;
}

const char *StringForEnum(AddMeshError::Enum error)
{
	if (error == AddMeshError::IndexOutOfRange)
		return "Index out of range";
	if (error == AddMeshError::InvalidIndexCount)
		return "Invalid index count";
	return "Success";
}

const char *StringForEnum(ProgressCategory::Enum category)
{
	if (category == ProgressCategory::ComputingCharts)
		return "Computing charts";
	if (category == ProgressCategory::ParametizingCharts)
		return "Parametizing charts";
	if (category == ProgressCategory::PackingCharts)
		return "Packing charts";
	if (category == ProgressCategory::BuildingOutputMeshes)
		return "Building output meshes";
	return "";
}

} // namespace xatlas
