// This code is in the public domain -- castanyo@yahoo.es
#pragma once
#ifndef XATLAS_RASTER_H
#define XATLAS_RASTER_H

namespace xatlas {
namespace internal {

class Vector2
{
public:
	typedef Vector2 const &Arg;

	Vector2();
	explicit Vector2(float f);
	Vector2(float x, float y);
	Vector2(Vector2::Arg v);
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

	union
	{
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4201)
#endif
		struct
		{
			float x, y;
		};
#ifdef _MSC_VER
#pragma warning(pop)
#endif

		float component[2];
	};
};

class Vector3
{
public:
	typedef Vector3 const &Arg;

	Vector3();
	explicit Vector3(float x);
	Vector3(float x, float y, float z);
	Vector3(Vector2::Arg v, float z);
	Vector3(Vector3::Arg v);
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

	union
	{
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4201)
#endif
		struct
		{
			float x, y, z;
		};
#ifdef _MSC_VER
#pragma warning(pop)
#endif

		float component[3];
	};
};

namespace raster {

/// A callback to sample the environment. Return false to terminate rasterization.
typedef bool (* SamplingCallback)(void *param, int x, int y, Vector3::Arg bar, Vector3::Arg dx, Vector3::Arg dy, float coverage);

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
} // namespace internal
} // namespace xatlas

#endif // XATLAS_RASTER_H
