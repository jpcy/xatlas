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
	const void *vertexNormalData; // optional
	uint32_t vertexNormalStride; // optional
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

} // namespace xatlas

#endif // XATLAS_H
