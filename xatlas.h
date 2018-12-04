// This code is in the public domain -- castanyo@yahoo.es
#pragma once
#ifndef XATLAS_H
#define XATLAS_H
#include <stdint.h>
#include <float.h> // FLT_MAX

namespace xatlas {

struct Atlas;

struct PrintFlags
{
	enum Enum
	{
		BuildingOutputMeshes = 1<<0,
		ComputingCharts = 1<<1,
		MeshCreation = 1<<2,
		MeshProcessing = 1<<3,
		MeshWarnings = 1<<4,
		PackingCharts = 1<<5,
		ParametizingCharts = 1<<6,
		All = ~0
	};
};

typedef int (*PrintFunc)(const char *, ...);

void SetPrint(int flags, PrintFunc print = NULL);
Atlas *Create();
void Destroy(Atlas *atlas);

struct AddMeshError
{
	enum Enum
	{
		Success,
		IndexOutOfRange, // index0 is the index
		InvalidIndexCount // not evenly divisible by 3 - expecting triangles
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

struct InputMesh
{
	uint32_t vertexCount;
	const void *vertexPositionData;
	uint32_t vertexPositionStride;
	const void *vertexNormalData; // optional
	uint32_t vertexNormalStride; // optional
	const void *vertexUvData; // optional. The input UVs are provided as a hint to the chart generator.
	uint32_t vertexUvStride;
	uint32_t indexCount;
	const void *indexData;
	IndexFormat::Enum indexFormat;
	
	// optional. indexCount / 3 in length.
	// Don't atlas faces set to true. Faces will still exist in the output meshes, uv will be (0, 0).
	const bool *faceIgnoreData;

	InputMesh()
	{
		vertexCount = 0;
		vertexPositionData = NULL;
		vertexPositionStride = 0;
		vertexNormalData = NULL;
		vertexNormalStride = 0;
		vertexUvData = NULL;
		vertexUvStride = 0;
		indexCount = 0;
		indexData = NULL;
		indexFormat = IndexFormat::UInt16;
		faceIgnoreData = NULL;
	}
};

// useColocalVertices - generates fewer charts (good), but is more sensitive to bad geometry.
AddMeshError::Enum AddMesh(Atlas *atlas, const InputMesh &mesh, bool useColocalVertices = true);

struct ProgressCategory
{
	enum Enum
	{
		ComputingCharts,
		ParametizingCharts,
		PackingCharts,
		BuildingOutputMeshes
	};
};

typedef void (*ProgressCallback)(ProgressCategory::Enum category, int progress, void *userData);

struct CharterOptions
{
	float proxyFitMetricWeight;
	float roundnessMetricWeight;
	float straightnessMetricWeight;
	float normalSeamMetricWeight;
	float textureSeamMetricWeight;
	float maxChartArea;
	float maxBoundaryLength;
	float maxThreshold;
	uint32_t growFaceCount;
	uint32_t maxIterations;

	CharterOptions()
	{
		// These are the default values we use on The Witness.
		proxyFitMetricWeight = 2.0f;
		roundnessMetricWeight = 0.01f;
		straightnessMetricWeight = 6.0f;
		normalSeamMetricWeight = 4.0f;
		textureSeamMetricWeight = 0.5f;
		/*
		proxyFitMetricWeight = 1.0f;
		roundnessMetricWeight = 0.1f;
		straightnessMetricWeight = 0.25f;
		normalSeamMetricWeight = 1.0f;
		textureSeamMetricWeight = 0.1f;
		*/
		maxChartArea = FLT_MAX;
		maxBoundaryLength = FLT_MAX;
		maxThreshold = 2;
		growFaceCount = 32;
		maxIterations = 4;
	}
};

void GenerateCharts(Atlas *atlas, CharterOptions charterOptions = CharterOptions(), ProgressCallback progressCallback = NULL, void *progressCallbackUserData = NULL);

struct PackerOptions
{
	// The number of attempts to find a suitable random chart location.
	// 0 is brute force - slowest, but best results.
	int attempts;

	// This is not really texel area, but 1 / texel width?
	// If 0, an estimated value will be calculated to try and match the given resolution.
	// If resolution is also 0, the estimated value will try to match a 1024x1024 atlas.
	float texelArea;

	// If 0, generate a single atlas with texelArea determining the final resolution.
	// If not 0, generate 1 or more atlases with that exact resolution.
	uint32_t resolution;

	// Align charts to 4x4 blocks. 
	bool blockAlign;

	// Pack charts with extra padding.
	bool conservative;

	// Number of pixels to pad. conservative must be true.
	int padding;

	PackerOptions()
	{
		attempts = 4096;
		texelArea = 0;
		resolution = 0;
		blockAlign = false;
		conservative = false;
		padding = 0;
	}
};

void PackCharts(Atlas *atlas, PackerOptions packerOptions = PackerOptions(), ProgressCallback progressCallback = NULL, void *progressCallbackUserData = NULL);

struct OutputChart
{
	uint32_t atlasIndex;
	uint32_t *indexArray;
	uint32_t indexCount;
};

struct OutputVertex
{
	int32_t atlasIndex; // -1 if the vertex doesn't exist in any atlas.
	float uv[2]; // not normalized
	uint32_t xref; // Index of input vertex from which this output vertex originated.
};

struct OutputMesh
{
	OutputChart *chartArray;
	uint32_t chartCount;
	uint32_t *indexArray;
	uint32_t indexCount;
	OutputVertex *vertexArray;
	uint32_t vertexCount;
};

uint32_t GetWidth(const Atlas *atlas);
uint32_t GetHeight(const Atlas *atlas);
uint32_t GetNumAtlases(const Atlas *atlas);
uint32_t GetNumCharts(const Atlas *atlas);
const OutputMesh * const *GetOutputMeshes(const Atlas *atlas);
const char *StringForEnum(AddMeshError::Enum error);
const char *StringForEnum(ProgressCategory::Enum category);

} // namespace xatlas

#endif // XATLAS_H
