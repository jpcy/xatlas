/*
MIT License

Copyright (c) 2018-2020 Jonathan Young

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
#ifndef XATLAS_C_H
#define XATLAS_C_H
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef XATLAS_EXPORT_API
#define XATLAS_EXPORT_API 0
#endif

#ifndef XATLAS_IMPORT_API
#define XATLAS_IMPORT_API 0
#endif

#ifndef XATLAS_API
	#if XATLAS_EXPORT_API
		#ifdef _MSC_VER
			#define XATLAS_API __declspec(dllexport)
		#else
			#define XATLAS_API __attribute__((visibility("default")))
		#endif
	#elif XATLAS_IMPORT_API
		#ifdef _MSC_VER
			#define XATLAS_API __declspec(dllimport)
		#else
			#define XATLAS_API
		#endif
	#else
		#define XATLAS_API
	#endif
#endif

typedef enum
{
	XATLAS_CHART_TYPE_PLANAR,
	XATLAS_CHART_TYPE_ORTHO,
	XATLAS_CHART_TYPE_LSCM,
	XATLAS_CHART_TYPE_PIECEWISE,
	XATLAS_CHART_TYPE_INVALID
}
xatlasChartType;

typedef struct
{
	uint32_t *faceArray;
	uint32_t atlasIndex;
	uint32_t faceCount;
	xatlasChartType type;
	uint32_t material;
}
xatlasChart;

typedef struct
{
	int32_t atlasIndex;
	int32_t chartIndex;
	float uv[2];
	uint32_t xref;
}
xatlasVertex;

typedef struct
{
	xatlasChart *chartArray;
	uint32_t *indexArray;
	xatlasVertex *vertexArray;
	uint32_t chartCount;
	uint32_t indexCount;
	uint32_t vertexCount;
}
xatlasMesh;

static const uint32_t xatlasImageChartIndexMask = 0x1FFFFFFF;
static const uint32_t xatlasImageHasChartIndexBit = 0x80000000;
static const uint32_t xatlasImageIsBilinearBit = 0x40000000;
static const uint32_t xatlasImageIsPaddingBit = 0x20000000;

typedef struct
{
	uint32_t *image;
	xatlasMesh *meshes;
	float *utilization;
	uint32_t width;
	uint32_t height;
	uint32_t atlasCount;
	uint32_t chartCount;
	uint32_t meshCount;
	float texelsPerUnit;
}
xatlasAtlas;

typedef enum
{
	XATLAS_INDEX_FORMAT_UINT16,
	XATLAS_INDEX_FORMAT_UINT32
}
xatlasIndexFormat;

typedef struct
{
	const void *vertexPositionData;
	const void *vertexNormalData;
	const void *vertexUvData;
	const void *indexData;
	const bool *faceIgnoreData;
	const uint32_t *faceMaterialData;
	const uint8_t *faceVertexCount;
	uint32_t vertexCount;
	uint32_t vertexPositionStride;
	uint32_t vertexNormalStride;
	uint32_t vertexUvStride;
	uint32_t indexCount;
	int32_t indexOffset;
	uint32_t faceCount;
	xatlasIndexFormat indexFormat;
	float epsilon;
}
xatlasMeshDecl;

typedef enum
{
	XATLAS_ADD_MESH_ERROR_SUCCESS,
	XATLAS_ADD_MESH_ERROR_ERROR,
	XATLAS_ADD_MESH_ERROR_INDEXOUTOFRANGE,
	XATLAS_ADD_MESH_ERROR_INVALIDFACEVERTEXCOUNT,
	XATLAS_ADD_MESH_ERROR_INVALIDINDEXCOUNT
}
xatlasAddMeshError;

typedef struct
{
	const void *vertexUvData;
	const void *indexData;
	const uint32_t *faceMaterialData;
	uint32_t vertexCount;
	uint32_t vertexStride;
	uint32_t indexCount;
	int32_t indexOffset;
	xatlasIndexFormat indexFormat;
}
xatlasUvMeshDecl;

typedef void (*xatlasParameterizeFunc)(const float *positions, float *texcoords, uint32_t vertexCount, const uint32_t *indices, uint32_t indexCount);

typedef struct
{
	xatlasParameterizeFunc paramFunc;
	float maxChartArea;
	float maxBoundaryLength;
	float normalDeviationWeight;
	float roundnessWeight;
	float straightnessWeight;
	float normalSeamWeight;
	float textureSeamWeight;
	float maxCost;
	uint32_t maxIterations;
	bool useInputMeshUvs;
	bool fixWinding;
}
xatlasChartOptions;

typedef struct
{
	uint32_t maxChartSize;
	uint32_t padding;
	float texelsPerUnit;
	uint32_t resolution;
	bool bilinear;
	bool blockAlign;
	bool bruteForce;
	bool createImage;
	bool rotateChartsToAxis;
	bool rotateCharts;
}
xatlasPackOptions;

typedef enum
{
	XATLAS_PROGRESS_CATEGORY_ADDMESH,
	XATLAS_PROGRESS_CATEGORY_COMPUTECHARTS,
	XATLAS_PROGRESS_CATEGORY_PACKCHARTS,
	XATLAS_PROGRESS_CATEGORY_BUILDOUTPUTMESHES
}
xatlasProgressCategory;

typedef bool (*xatlasProgressFunc)(xatlasProgressCategory category, int progress, void *userData);
typedef void *(*xatlasReallocFunc)(void *, size_t);
typedef void (*xatlasFreeFunc)(void *);
typedef int (*xatlasPrintFunc)(const char *, ...);

XATLAS_API xatlasAtlas *xatlasCreate();
XATLAS_API void xatlasDestroy(xatlasAtlas *atlas);
XATLAS_API xatlasAddMeshError xatlasAddMesh(xatlasAtlas *atlas, const xatlasMeshDecl *meshDecl, uint32_t meshCountHint);
XATLAS_API void xatlasAddMeshJoin(xatlasAtlas *atlas);
XATLAS_API xatlasAddMeshError xatlasAddUvMesh(xatlasAtlas *atlas, const xatlasUvMeshDecl *decl);
XATLAS_API void xatlasComputeCharts(xatlasAtlas *atlas, const xatlasChartOptions *chartOptions);
XATLAS_API void xatlasPackCharts(xatlasAtlas *atlas, const xatlasPackOptions *packOptions);
XATLAS_API void xatlasGenerate(xatlasAtlas *atlas, const xatlasChartOptions *chartOptions, const xatlasPackOptions *packOptions);
XATLAS_API void xatlasSetProgressCallback(xatlasAtlas *atlas, xatlasProgressFunc progressFunc, void *progressUserData);
XATLAS_API void xatlasSetAlloc(xatlasReallocFunc reallocFunc, xatlasFreeFunc freeFunc);
XATLAS_API void xatlasSetPrint(xatlasPrintFunc print, bool verbose);
XATLAS_API const char *xatlasAddMeshErrorString(xatlasAddMeshError error);
XATLAS_API const char *xatlasProgressCategoryString(xatlasProgressCategory category);
XATLAS_API void xatlasMeshDeclInit(xatlasMeshDecl *meshDecl);
XATLAS_API void xatlasUvMeshDeclInit(xatlasUvMeshDecl *uvMeshDecl);
XATLAS_API void xatlasChartOptionsInit(xatlasChartOptions *chartOptions);
XATLAS_API void xatlasPackOptionsInit(xatlasPackOptions *packOptions);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // XATLAS_C_H
