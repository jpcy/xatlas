/*
MIT License

Copyright (c) 2018-2019 Jonathan Young

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
/*
thekla_atlas
MIT License
https://github.com/Thekla/thekla_atlas
Copyright (c) 2013 Thekla, Inc
Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>
*/
#pragma once
#ifndef XATLAS_H
#define XATLAS_H
#include <stdint.h>

namespace xatlas {

struct Chart
{
	uint32_t atlasIndex;
	uint32_t *indexArray;
	uint32_t indexCount;
};

struct Vertex
{
	int32_t atlasIndex; // -1 if the vertex doesn't exist in any atlas.
	float uv[2]; // Not normalized, values are in Atlas width and height range.
	uint32_t xref; // Index of input vertex from which this output vertex originated.
};

struct Mesh
{
	Chart *chartArray;
	uint32_t chartCount;
	uint32_t *indexArray;
	uint32_t indexCount;
	Vertex *vertexArray;
	uint32_t vertexCount;
};

struct Atlas
{
	uint32_t width;
	uint32_t height;
	uint32_t atlasCount;
	uint32_t chartCount;
	uint32_t meshCount;
	Mesh *meshes;

	// Normalized atlas texel utilization. atlasCount in length.
	float *utilization;

	// Equal to PackOptions texelsPerUnit if texelsPerUnit > 0, otherwise an estimated value to try and match PackOptions resolution.
	float texelsPerUnit;
};

Atlas *Create();
void Destroy(Atlas *atlas);

struct AddMeshError
{
	enum Enum
	{
		Success,
		IndexOutOfRange,
		InvalidIndexCount // Not evenly divisible by 3 - expecting triangles.
	};
};

struct IndexFormat
{
	enum Enum
	{
		UInt16,
		UInt32
	};
};

struct MeshDecl
{
	uint32_t vertexCount = 0;
	const void *vertexPositionData = nullptr;
	uint32_t vertexPositionStride = 0;
	const void *vertexNormalData = nullptr; // optional
	uint32_t vertexNormalStride = 0; // optional
	const void *vertexUvData = nullptr; // optional. The input UVs are provided as a hint to the chart generator.
	uint32_t vertexUvStride = 0; // optional
	uint32_t indexCount = 0;
	const void *indexData = nullptr;
	int32_t indexOffset = 0; // optional. Add this offset to all indices.
	IndexFormat::Enum indexFormat = IndexFormat::UInt16;
	
	// optional. indexCount / 3 in length.
	// Don't atlas faces set to true. Faces will still exist in the output meshes, Vertex uv will be (0, 0) and Vertex atlasIndex will be -1.
	const bool *faceIgnoreData = nullptr;
};

AddMeshError::Enum AddMesh(Atlas *atlas, const MeshDecl &meshDecl);

struct ProgressCategory
{
	enum Enum
	{
		ComputeCharts,
		ParameterizeCharts,
		PackCharts,
		BuildOutputMeshes
	};
};

typedef void (*ProgressFunc)(ProgressCategory::Enum category, int progress, void *userData);

struct ChartOptions
{
	float proxyFitMetricWeight = 2.0f;
	float roundnessMetricWeight = 0.01f;
	float straightnessMetricWeight = 6.0f;
	float normalSeamMetricWeight = 4.0f; // If > 1000, normal seams are fully respected.
	float textureSeamMetricWeight = 0.5f;
	float maxChartArea = 0.0f; // Don't grow charts to be larger than this. 0 means no limit.
	float maxBoundaryLength = 0.0f; // Don't grow charts to have a longer boundary than this. 0 means no limit.
	float maxThreshold = 2.0f;
	uint32_t growFaceCount = 32;
	uint32_t maxIterations = 1;
};

typedef void (*ParameterizeFunc)(const float *positions, float *texcoords, uint32_t vertexCount, const uint32_t *indices, uint32_t indexCount, bool isPlanar);

struct PackOptions
{
	// The number of attempts to find a suitable random chart location.
	// 0 is brute force - very slow, but best results. Faster if blockAlign is true;
	int attempts = 4096;

	// Unit to texel scale. e.g. a 1x1 quad with texelsPerUnit of 32 will take up approximately 32x32 texels in the atlas.
	// If 0, an estimated value will be calculated to try and match the given resolution.
	// If resolution is also 0, the estimated value will try to match a 1024x1024 atlas.
	float texelsPerUnit = 0.0f;

	// If 0, generate a single atlas with texelsPerUnit determining the final resolution.
	// If not 0, generate 1 or more atlases with that exact resolution.
	uint32_t resolution = 0;

	// Charts larger than this will be scaled down.
	uint32_t maxChartSize = 1024;

	// Align charts to 4x4 blocks. 
	bool blockAlign = false;

	// Pack charts with extra padding.
	bool conservative = false;

	// Number of pixels to pad. conservative must be true.
	uint32_t padding = 0;
};

// Equivalent to calling ComputeCharts, ParameterizeCharts and PackCharts in sequence. Can be called multiple times to regenerate with different options.
void Generate(Atlas *atlas, ChartOptions chartOptions = ChartOptions(), ParameterizeFunc paramFunc = nullptr, PackOptions packOptions = PackOptions(), ProgressFunc progressFunc = nullptr, void *progressUserData = nullptr);

// Call after AddMesh. Can be called multiple times to recompute charts with different options.
void ComputeCharts(Atlas *atlas, ChartOptions chartOptions = ChartOptions(), ProgressFunc progressFunc = nullptr, void *progressUserData = nullptr);

// Call after ComputeCharts. Can be called multiple times to re-parameterize charts with a different ParameterizeFunc.
void ParameterizeCharts(Atlas *atlas, ParameterizeFunc func = nullptr, ProgressFunc progressFunc = nullptr, void *progressUserData = nullptr);

// Call after ParameterizeCharts. Can be called multiple times to re-pack charts with different options.
void PackCharts(Atlas *atlas, PackOptions packOptions = PackOptions(), ProgressFunc progressFunc = nullptr, void *progressUserData = nullptr);

typedef void *(*ReallocFunc)(void *, size_t);
void SetRealloc(ReallocFunc reallocFunc);

typedef int (*PrintFunc)(const char *, ...);
void SetPrint(PrintFunc print, bool verbose);

const char *StringForEnum(AddMeshError::Enum error);
const char *StringForEnum(ProgressCategory::Enum category);

} // namespace xatlas

#endif // XATLAS_H
