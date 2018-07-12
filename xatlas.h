// This code is in the public domain -- castanyo@yahoo.es
#pragma once
#ifndef XATLAS_H
#define XATLAS_H

namespace xatlas {

struct PackMethod
{
	enum Enum
	{
		TexelArea, // texel_area determines resolution
		ApproximateResolution, // guess texel_area to approximately match desired resolution
		ExactResolution // run the packer multiple times to exactly match the desired resolution (slow)
	};
};

struct Options
{
	struct
	{
        float proxy_fit_metric_weight;
        float roundness_metric_weight;
        float straightness_metric_weight;
        float normal_seam_metric_weight;
        float texture_seam_metric_weight;
        float max_chart_area;
        float max_boundary_length;
    }
	charter;

	struct
	{
		PackMethod::Enum method;
        int packing_quality;
        float texel_area;       // This is not really texel area, but 1 / texel width?
		uint32_t resolution;
        bool block_align;       // Align charts to 4x4 blocks. 
        bool conservative;      // Pack charts with extra padding.
		int padding;
    }
	packer;

	void (*Print)(const char *, ...);
};

struct Input_Mesh
{
	uint32_t vertexCount;
    const void *vertexPositionData;
	uint32_t vertexPositionStride;
	const void *vertexNormalData;
	uint32_t vertexNormalStride;
	uint32_t indexCount;
    const uint32_t *indexData;
	const uint16_t *faceMaterialData; // optional. indexCount / 3 in length.
};

struct Output_Vertex
{
    float uv[2];
    uint32_t xref;   // Index of input vertex from which this output vertex originated.
};

struct Output_Mesh
{
	uint32_t vertexCount;
	uint32_t indexCount;
    Output_Vertex *vertexArray;
	uint32_t *indexArray;
};

enum Error
{
    Error_Success,
    Error_Invalid_Args,
    Error_Invalid_Options,
    Error_Invalid_Mesh,
    Error_Invalid_Mesh_Non_Manifold,
    Error_Not_Implemented,
};

struct Atlas
{
	Error error;
	int errorMeshIndex;
	int width;
	int height;
	int nCharts;
	int nMeshes;
	Output_Mesh **meshes;
};

void set_default_options(Options * options);
void add_mesh(const Input_Mesh *mesh);
Atlas atlas_generate(const Options *options);
void atlas_free(Atlas atlas);

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

#endif // XATLAS_H
