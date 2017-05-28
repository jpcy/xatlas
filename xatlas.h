// This code is in the public domain -- castanyo@yahoo.es
#pragma once
#ifndef XATLAS_H
#define XATLAS_H

#ifndef xaAssert
#define xaAssert(exp) if (!(exp)) { xaPrint("%s %s %s\n", #exp, __FILE__, __LINE__); }
#endif
#ifndef xaDebugAssert
#define xaDebugAssert(exp) assert(exp)
#endif
#ifndef xaPrint
#define xaInternalPrint
#define xaPrint(...) xatlas::internal::Print(__VA_ARGS__)
#endif

namespace xatlas {
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
        int packing_quality;
        float texel_area;       // This is not really texel area, but 1 / texel width?
        bool block_align;       // Align charts to 4x4 blocks. 
        bool conservative;      // Pack charts with extra padding.
		int padding;
    }
	packer;
};

struct Input_Vertex
{
    float position[3];
    float normal[3];
    float uv[2];
    int first_colocal;
};

struct Input_Face
{
    int vertex_index[3];
    int material_index;
};

struct Input_Mesh
{
    int vertex_count;
    int face_count;
    Input_Vertex * vertex_array;
    Input_Face * face_array;
};

struct Output_Vertex
{
    float uv[2];
    int xref;   // Index of input vertex from which this output vertex originated.
};

struct Output_Mesh
{
    int vertex_count;
    int index_count;
    Output_Vertex * vertex_array;
    int * index_array;
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
	int width;
	int height;
	int nMeshes;
	Output_Mesh **meshes;
};

void set_default_options(Options * options);
void add_mesh(const Input_Mesh *mesh);
Atlas atlas_generate(const Options *options);
void atlas_free(Atlas atlas);

#ifdef _MSC_VER
// Ignore gcc attributes.
#define __attribute__(X)
#endif

namespace internal {
#ifdef xaInternalPrint
void Print( const char *msg, ... ) __attribute__((format (printf, 1, 2)));
#endif

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
