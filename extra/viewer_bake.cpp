/*
xatlas
https://github.com/jpcy/xatlas
Copyright (c) 2018 Jonathan Young

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#include <atomic>
#include <cstddef>
#include <mutex>
#include <thread>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <bx/os.h>
#include <bx/rng.h>
#include <embree3/rtcore.h>
#include <imgui/imgui.h>
#include <OpenImageDenoise/oidn.h>
#include <TaskScheduler_c.h>
#include "shaders/shared.h"
#include "../xatlas.h"
#include "viewer.h"

static bx::Vec3 operator+(bx::Vec3 a, bx::Vec3 b) { return bx::add(a, b); }
static bx::Vec3 operator-(bx::Vec3 a, bx::Vec3 b) { return bx::sub(a, b); }
static bx::Vec3 operator*(bx::Vec3 a, bx::Vec3 b) { return bx::mul(a, b); }
static bx::Vec3 operator*(bx::Vec3 a, float s) { return bx::mul(a, s); }

/***********************************************************
* A single header file OpenGL lightmapping library         *
* https://github.com/ands/lightmapper                      *
* no warranty implied | use at your own risk               *
* author: Andreas Mantler (ands) | last change: 12.06.2016 *
*                                                          *
* License:                                                 *
* This software is in the public domain.                   *
* Where that dedication is not recognized,                 *
* you are granted a perpetual, irrevocable license to copy *
* and modify this file however you want.                   *
***********************************************************/
typedef struct lm_vec2 { float x, y; } lm_vec2;
static lm_vec2  lm_v2i       (int     x, int     y) { lm_vec2 v = { (float)x, (float)y }; return v; }
static lm_vec2  lm_v2        (float   x, float   y) { lm_vec2 v = { x, y }; return v; }
static lm_vec2  lm_add2      (lm_vec2 a, lm_vec2 b) { return lm_v2(a.x + b.x, a.y + b.y); }
static lm_vec2  lm_sub2      (lm_vec2 a, lm_vec2 b) { return lm_v2(a.x - b.x, a.y - b.y); }
static lm_vec2  lm_scale2    (lm_vec2 a, float   b) { return lm_v2(a.x * b, a.y * b); }
static lm_vec2  lm_div2      (lm_vec2 a, float   b) { return lm_scale2(a, 1.0f / b); }
static lm_vec2  lm_min2      (lm_vec2 a, lm_vec2 b) { return lm_v2(bx::min(a.x, b.x), bx::min(a.y, b.y)); }
static lm_vec2  lm_max2      (lm_vec2 a, lm_vec2 b) { return lm_v2(bx::max(a.x, b.x), bx::max(a.y, b.y)); }
static lm_vec2  lm_floor2    (lm_vec2 a           ) { return lm_v2(floorf(a.x), floorf(a.y)); }
static lm_vec2  lm_ceil2     (lm_vec2 a           ) { return lm_v2(ceilf (a.x), ceilf (a.y)); }
static float    lm_dot2      (lm_vec2 a, lm_vec2 b) { return a.x * b.x + a.y * b.y; }
static float    lm_cross2    (lm_vec2 a, lm_vec2 b) { return a.x * b.y - a.y * b.x; } // pseudo cross product

static lm_vec2 lm_toBarycentric(lm_vec2 p1, lm_vec2 p2, lm_vec2 p3, lm_vec2 p)
{
	// http://www.blackpawn.com/texts/pointinpoly/
	// Compute vectors
	lm_vec2 v0 = lm_sub2(p3, p1);
	lm_vec2 v1 = lm_sub2(p2, p1);
	lm_vec2 v2 = lm_sub2(p, p1);
	// Compute dot products
	float dot00 = lm_dot2(v0, v0);
	float dot01 = lm_dot2(v0, v1);
	float dot02 = lm_dot2(v0, v2);
	float dot11 = lm_dot2(v1, v1);
	float dot12 = lm_dot2(v1, v2);
	// Compute barycentric coordinates
	float invDenom = 1.0f / (dot00 * dot11 - dot01 * dot01);
	float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	float v = (dot00 * dot12 - dot01 * dot02) * invDenom;
	return lm_v2(u, v);
}

static inline int lm_leftOf(lm_vec2 a, lm_vec2 b, lm_vec2 c)
{
	float x = lm_cross2(lm_sub2(b, a), lm_sub2(c, b));
	return x < 0 ? -1 : x > 0;
}

static bool lm_lineIntersection(lm_vec2 x0, lm_vec2 x1, lm_vec2 y0, lm_vec2 y1, lm_vec2* res)
{
	lm_vec2 dx = lm_sub2(x1, x0);
	lm_vec2 dy = lm_sub2(y1, y0);
	lm_vec2 d = lm_sub2(x0, y0);
	float dyx = lm_cross2(dy, dx);
	if (dyx == 0.0f)
		return false;
	dyx = lm_cross2(d, dx) / dyx;
	if (dyx <= 0 || dyx >= 1)
		return false;
	res->x = y0.x + dyx * dy.x;
	res->y = y0.y + dyx * dy.y;
	return true;
}

// this modifies the poly array! the poly array must be big enough to hold the result!
// res must be big enough to hold the result!
static int lm_convexClip(lm_vec2 *poly, int nPoly, const lm_vec2 *clip, int nClip, lm_vec2 *res)
{
	int nRes = nPoly;
	int dir = lm_leftOf(clip[0], clip[1], clip[2]);
	for (int i = 0, j = nClip - 1; i < nClip && nRes; j = i++)
	{
		if (i != 0)
			for (nPoly = 0; nPoly < nRes; nPoly++)
				poly[nPoly] = res[nPoly];
		nRes = 0;
		lm_vec2 v0 = poly[nPoly - 1];
		int side0 = lm_leftOf(clip[j], clip[i], v0);
		if (side0 != -dir)
			res[nRes++] = v0;
		for (int k = 0; k < nPoly; k++)
		{
			lm_vec2 v1 = poly[k], x;
			int side1 = lm_leftOf(clip[j], clip[i], v1);
			if (side0 + side1 == 0 && side0 && lm_lineIntersection(clip[j], clip[i], v0, v1, &x))
				res[nRes++] = x;
			if (k == nPoly - 1)
				break;
			if (side1 != -dir)
				res[nRes++] = v1;
			v0 = v1;
			side0 = side1;
		}
	}

	return nRes;
}

struct lm_context
{
	struct
	{
		bx::Vec3 p[3];
		lm_vec2 uv[3];
	} triangle;

	struct
	{
		int minx, miny;
		int maxx, maxy;
		int x, y;
	} rasterizer;

	struct
	{
		bx::Vec3 position;
		bx::Vec3 direction;
	} sample;
};

static bool lm_hasConservativeTriangleRasterizerFinished(lm_context *ctx)
{
	return ctx->rasterizer.y >= ctx->rasterizer.maxy;
}

static void lm_moveToNextPotentialConservativeTriangleRasterizerPosition(lm_context *ctx)
{
	if (++ctx->rasterizer.x >= ctx->rasterizer.maxx)
	{
		ctx->rasterizer.x = ctx->rasterizer.minx;
		++ctx->rasterizer.y;
	}
}

static bool lm_trySamplingConservativeTriangleRasterizerPosition(lm_context *ctx)
{
	if (lm_hasConservativeTriangleRasterizerFinished(ctx))
		return false;

	lm_vec2 pixel[16];
	pixel[0] = lm_v2i(ctx->rasterizer.x, ctx->rasterizer.y);
	pixel[1] = lm_v2i(ctx->rasterizer.x + 1, ctx->rasterizer.y);
	pixel[2] = lm_v2i(ctx->rasterizer.x + 1, ctx->rasterizer.y + 1);
	pixel[3] = lm_v2i(ctx->rasterizer.x, ctx->rasterizer.y + 1);

	lm_vec2 res[16];
	int nRes = lm_convexClip(pixel, 4, ctx->triangle.uv, 3, res);
	if (nRes > 0)
	{
		// do centroid sampling
		lm_vec2 centroid = res[0];
		float area = res[nRes - 1].x * res[0].y - res[nRes - 1].y * res[0].x;
		for (int i = 1; i < nRes; i++)
		{
			centroid = lm_add2(centroid, res[i]);
			area += res[i - 1].x * res[i].y - res[i - 1].y * res[i].x;
		}
		centroid = lm_div2(centroid, (float)nRes);
		area = bx::abs(area / 2.0f);

		if (area > 0.0f)
		{
			// calculate 3D sample position and orientation
			lm_vec2 uv = lm_toBarycentric(
				ctx->triangle.uv[0],
				ctx->triangle.uv[1],
				ctx->triangle.uv[2],
				centroid);

			// sample it only if its's not degenerate
			if (bx::isFinite(uv.x) && bx::isFinite(uv.y))
			{
				const bx::Vec3 &p0 = ctx->triangle.p[0];
				const bx::Vec3 &p1 = ctx->triangle.p[1];
				const bx::Vec3 &p2 = ctx->triangle.p[2];
				const bx::Vec3 v1 = p1 - p0;
				const bx::Vec3 v2 = p2 - p0;
				ctx->sample.position = p0 + v2 * uv.x + v1 * uv.y;
				ctx->sample.direction = bx::normalize(bx::cross(v1, v2));

				if (bx::isFinite(ctx->sample.position.x) && bx::isFinite(ctx->sample.position.y) && bx::isFinite(ctx->sample.position.z) &&
					bx::isFinite(ctx->sample.direction.x) && bx::isFinite(ctx->sample.direction.y) && bx::isFinite(ctx->sample.direction.z) &&
					bx::length(ctx->sample.direction) > 0.5f) // don't allow 0.0f. should always be ~1.0f
				{
					return true;
				}
			}
		}
	}
	return false;
}

// returns true if a sampling position was found and
// false if we finished rasterizing the current triangle
static bool lm_findFirstConservativeTriangleRasterizerPosition(lm_context *ctx)
{
	while (!lm_trySamplingConservativeTriangleRasterizerPosition(ctx))
	{
		lm_moveToNextPotentialConservativeTriangleRasterizerPosition(ctx);
		if (lm_hasConservativeTriangleRasterizerFinished(ctx))
			return false;
	}
	return true;
}

static bool lm_findNextConservativeTriangleRasterizerPosition(lm_context *ctx)
{
	lm_moveToNextPotentialConservativeTriangleRasterizerPosition(ctx);
	return lm_findFirstConservativeTriangleRasterizerPosition(ctx);
}

enum class BakeStatus
{
	Idle,
	InitEmbree,
	Rasterizing,
	Tracing,
	ThreadFinished,
	Finished,
	Cancelled,
	Error
};

enum class DenoiseStatus
{
	Idle,
	Working,
	ThreadFinished,
	Finished,
	Error
};

enum class UpdateStatus
{
	Idle,
	Pending,
	Updating
};

struct BakeOptions
{
	bool fitToWindow = true;
	float scale = 1.0f;
	bx::Vec3 skyColor = bx::Vec3(1.0f);
	int maxDepth = 8;
};

static const int kMaxDepth = 16;

struct SampleLocation
{
	bx::Vec3 pos;
	bx::Vec3 normal;
	uint32_t uv[2];
};

struct TexelData
{
	bx::Vec3 accumColor;
	uint32_t numColorSamples;
	int numPathsTraced; // Seed for randomDirHemisphere.
	float randomOffset[2];
};

struct
{
	std::mutex errorMessageMutex;
	char errorMessage[256];
	bool initialized = false;
	std::atomic<BakeStatus> status;
	BakeOptions options;
	enkiTaskScheduler *taskScheduler;
	void *embreeLibrary = nullptr;
	RTCDevice embreeDevice = nullptr;
	RTCScene embreeScene = nullptr;
	RTCGeometry embreeGeometry = nullptr;
	void *oidnLibrary = nullptr;
	std::vector<objzMaterial *> triMaterials;
	// lightmap
	bgfx::TextureHandle lightmap, denoisedLightmap;
	uint32_t lightmapWidth, lightmapHeight;
	// worker thread
	std::thread *workerThread = nullptr;
	std::atomic<bool> stopWorker;
	std::atomic<bool> cancelWorker;
	std::vector<SampleLocation> sampleLocations;
	std::vector<uint32_t> sampleLocationRanks;
	std::vector<TexelData> texels;
	std::vector<float> lightmapData;
	std::atomic<uint32_t> numTrianglesRasterized;
	std::atomic<uint32_t> samplesPerTexelCount;
	uint32_t numRaysTraced;
	bx::RngMwc rng;
	// denoise thread
	std::atomic<DenoiseStatus> denoiseStatus;
	std::thread *denoiseThread = nullptr;
	double denoiseProgress;
	std::vector<float> denoisedLightmapData;
	// lightmap update
	std::vector<float> updateData;
	std::atomic<UpdateStatus> updateStatus;
	uint32_t updateFinishedFrameNo;
}
s_bake;
	
static thread_local bx::RngMwc s_rng;

#define EMBREE_LIB "embree3.dll"

namespace embree
{
	typedef RTCDevice (*NewDeviceFunc)(const char* config);
	typedef void (*ReleaseDeviceFunc)(RTCDevice device);
	typedef void (*SetDeviceErrorFunctionFunc)(RTCDevice device, RTCErrorFunction error, void* userPtr);
	typedef ssize_t (*GetDevicePropertyFunc)(RTCDevice device, enum RTCDeviceProperty prop);
	typedef RTCScene (*NewSceneFunc)(RTCDevice device);
	typedef void (*ReleaseSceneFunc)(RTCScene scene);
	typedef unsigned int (*AttachGeometryFunc)(RTCScene scene, RTCGeometry geometry);
	typedef void (*CommitSceneFunc)(RTCScene scene);
	typedef RTCGeometry (*NewGeometryFunc)(RTCDevice device, enum RTCGeometryType type);
	typedef void (*ReleaseGeometryFunc)(RTCGeometry geometry);
	typedef void (*SetSharedGeometryBufferFunc)(RTCGeometry geometry, enum RTCBufferType type, unsigned int slot, enum RTCFormat format, const void* ptr, size_t byteOffset, size_t byteStride, size_t itemCount);
	typedef void (*CommitGeometryFunc)(RTCGeometry geometry);
	typedef void (*Intersect1Func)(RTCScene scene, struct RTCIntersectContext* context, struct RTCRayHit* rayhit);
	typedef void (*Intersect4Func)(const int* valid, RTCScene scene, struct RTCIntersectContext* context, struct RTCRayHit4* rayhit);
	typedef void (*Intersect8Func)(const int* valid, RTCScene scene, struct RTCIntersectContext* context, struct RTCRayHit8* rayhit);
	typedef void (*Intersect16Func)(const int* valid, RTCScene scene, struct RTCIntersectContext* context, struct RTCRayHit16* rayhit);
	typedef void (*Occluded1Func)(RTCScene scene, struct RTCIntersectContext* context, struct RTCRay* ray);
	typedef void (*Occluded4Func)(const int* valid, RTCScene scene, struct RTCIntersectContext* context, struct RTCRay4* ray);
	typedef void (*Occluded8Func)(const int* valid, RTCScene scene, struct RTCIntersectContext* context, struct RTCRay8* ray);
	typedef void (*Occluded16Func)(const int* valid, RTCScene scene, struct RTCIntersectContext* context, struct RTCRay16* ray);
	NewDeviceFunc NewDevice;
	ReleaseDeviceFunc ReleaseDevice;
	SetDeviceErrorFunctionFunc SetDeviceErrorFunction;
	GetDevicePropertyFunc GetDeviceProperty;
	NewSceneFunc NewScene;
	ReleaseSceneFunc ReleaseScene;
	AttachGeometryFunc AttachGeometry;
	CommitSceneFunc CommitScene;
	NewGeometryFunc NewGeometry;
	ReleaseGeometryFunc ReleaseGeometry;
	SetSharedGeometryBufferFunc SetSharedGeometryBuffer;
	CommitGeometryFunc CommitGeometry;
	Intersect1Func Intersect1;
	Intersect4Func Intersect4;
	Intersect8Func Intersect8;
	Intersect16Func Intersect16;
	Occluded1Func Occluded1;
	Occluded4Func Occluded4;
	Occluded8Func Occluded8;
	Occluded16Func Occluded16;
};

static void bakeEmbreeError(void* /*userPtr*/, enum RTCError /*code*/, const char* str)
{
	fprintf(stderr, "Embree error: %s\n", str);
	exit(EXIT_FAILURE);
}

static bool bakeInitEmbree()
{
	if (!s_bake.embreeLibrary) {
		s_bake.embreeLibrary = bx::dlopen(EMBREE_LIB);
		if (!s_bake.embreeLibrary) {
			std::lock_guard<std::mutex> lock(s_bake.errorMessageMutex);
			bx::snprintf(s_bake.errorMessage, sizeof(s_bake.errorMessage), "Embree not installed. Cannot open '%s'.", EMBREE_LIB);
			return false;
		}
		embree::NewDevice = (embree::NewDeviceFunc)bx::dlsym(s_bake.embreeLibrary, "rtcNewDevice");
		embree::ReleaseDevice = (embree::ReleaseDeviceFunc)bx::dlsym(s_bake.embreeLibrary, "rtcReleaseDevice");
		embree::SetDeviceErrorFunction = (embree::SetDeviceErrorFunctionFunc)bx::dlsym(s_bake.embreeLibrary, "rtcSetDeviceErrorFunction");
		embree::GetDeviceProperty = (embree::GetDevicePropertyFunc)bx::dlsym(s_bake.embreeLibrary, "rtcGetDeviceProperty");
		embree::NewScene = (embree::NewSceneFunc)bx::dlsym(s_bake.embreeLibrary, "rtcNewScene");
		embree::ReleaseScene = (embree::ReleaseSceneFunc)bx::dlsym(s_bake.embreeLibrary, "rtcReleaseScene");
		embree::AttachGeometry = (embree::AttachGeometryFunc)bx::dlsym(s_bake.embreeLibrary, "rtcAttachGeometry");
		embree::CommitScene = (embree::CommitSceneFunc)bx::dlsym(s_bake.embreeLibrary, "rtcCommitScene");
		embree::NewGeometry = (embree::NewGeometryFunc)bx::dlsym(s_bake.embreeLibrary, "rtcNewGeometry");
		embree::ReleaseGeometry = (embree::ReleaseGeometryFunc)bx::dlsym(s_bake.embreeLibrary, "rtcReleaseGeometry");
		embree::SetSharedGeometryBuffer = (embree::SetSharedGeometryBufferFunc)bx::dlsym(s_bake.embreeLibrary, "rtcSetSharedGeometryBuffer");
		embree::CommitGeometry = (embree::CommitGeometryFunc)bx::dlsym(s_bake.embreeLibrary, "rtcCommitGeometry");
		embree::Intersect1 = (embree::Intersect1Func)bx::dlsym(s_bake.embreeLibrary, "rtcIntersect1");
		embree::Intersect4 = (embree::Intersect4Func)bx::dlsym(s_bake.embreeLibrary, "rtcIntersect4");
		embree::Intersect8 = (embree::Intersect8Func)bx::dlsym(s_bake.embreeLibrary, "rtcIntersect8");
		embree::Intersect16 = (embree::Intersect16Func)bx::dlsym(s_bake.embreeLibrary, "rtcIntersect16");
		embree::Occluded1 = (embree::Occluded1Func)bx::dlsym(s_bake.embreeLibrary, "rtcOccluded1");
		embree::Occluded4 = (embree::Occluded4Func)bx::dlsym(s_bake.embreeLibrary, "rtcOccluded4");
		embree::Occluded8 = (embree::Occluded8Func)bx::dlsym(s_bake.embreeLibrary, "rtcOccluded8");
		embree::Occluded16 = (embree::Occluded16Func)bx::dlsym(s_bake.embreeLibrary, "rtcOccluded16");
	}
	if (!s_bake.embreeDevice) {
		s_bake.embreeDevice = embree::NewDevice(nullptr);
		if (!s_bake.embreeDevice) {
			std::lock_guard<std::mutex> lock(s_bake.errorMessageMutex);
			bx::snprintf(s_bake.errorMessage, sizeof(s_bake.errorMessage), "Error creating Embree device");
			return false;
		}
		embree::SetDeviceErrorFunction(s_bake.embreeDevice, bakeEmbreeError, nullptr);
		int version[3];
		for (int i = 0; i < 3; i++)
			version[i] = (int)embree::GetDeviceProperty(s_bake.embreeDevice, RTCDeviceProperty(RTC_DEVICE_PROPERTY_VERSION_MAJOR + i));
		printf("Embree version %d.%d.%d\n", version[0], version[1], version[2]);
		if (version[0] != RTC_VERSION_MAJOR || version[1] != RTC_VERSION_MINOR || version[2] != RTC_VERSION_PATCH)
			printf("   Expected version is %d.%d.%d\n", RTC_VERSION_MAJOR, RTC_VERSION_MINOR, RTC_VERSION_PATCH);
		if (embree::GetDeviceProperty(s_bake.embreeDevice, RTC_DEVICE_PROPERTY_NATIVE_RAY4_SUPPORTED))
			printf("   RTC_DEVICE_PROPERTY_NATIVE_RAY4_SUPPORTED\n");
		if (embree::GetDeviceProperty(s_bake.embreeDevice, RTC_DEVICE_PROPERTY_NATIVE_RAY8_SUPPORTED))
			printf("   RTC_DEVICE_PROPERTY_NATIVE_RAY8_SUPPORTED\n");
		if (embree::GetDeviceProperty(s_bake.embreeDevice, RTC_DEVICE_PROPERTY_NATIVE_RAY16_SUPPORTED))
			printf("   RTC_DEVICE_PROPERTY_NATIVE_RAY16_SUPPORTED\n");
		if (embree::GetDeviceProperty(s_bake.embreeDevice, RTC_DEVICE_PROPERTY_RAY_STREAM_SUPPORTED))
			printf("   RTC_DEVICE_PROPERTY_RAY_STREAM_SUPPORTED\n");
		s_bake.embreeGeometry = embree::NewGeometry(s_bake.embreeDevice, RTC_GEOMETRY_TYPE_TRIANGLE);
		const objzModel *model = modelGetData();
		embree::SetSharedGeometryBuffer(s_bake.embreeGeometry, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, model->vertices, offsetof(ModelVertex, pos), sizeof(ModelVertex), (size_t)model->numVertices);
		embree::SetSharedGeometryBuffer(s_bake.embreeGeometry, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, model->indices, 0, sizeof(uint32_t) * 3, (size_t)model->numIndices / 3);
		embree::CommitGeometry(s_bake.embreeGeometry);
		s_bake.embreeScene = embree::NewScene(s_bake.embreeDevice);
		embree::AttachGeometry(s_bake.embreeScene, s_bake.embreeGeometry);
		embree::CommitScene(s_bake.embreeScene);
	}
	return true;
}

static int compareSampleLocations(const void *a, const void *b)
{
	const uint32_t i0 = *(const uint32_t *)a;
	const uint32_t i1 = *(const uint32_t *)b;
	const uint32_t offset0 = s_bake.sampleLocations[i0].uv[0] + s_bake.sampleLocations[i0].uv[1] * s_bake.lightmapWidth;
	const uint32_t offset1 = s_bake.sampleLocations[i1].uv[0] + s_bake.sampleLocations[i1].uv[1] * s_bake.lightmapWidth;
	if (offset0 < offset1)
		return -1;
	else if (offset0 > offset1)
		return 1;
	return 0;
}

static bool bakeRasterize()
{
	// Only need to do this once, unless the atlas or models changes - bakeClear will be called.
	if (!s_bake.sampleLocations.empty())
		return true;
	s_bake.numTrianglesRasterized = 0;
	std::vector<ModelVertex> &modelVertices = *atlasGetVertices();
	std::vector<uint32_t> &modelIndices = *atlasGetIndices();
	// Don't allow duplicate samples at the same uv.
	std::vector<bool> sampleExists;
	sampleExists.resize(s_bake.lightmapWidth * s_bake.lightmapHeight);
	for (uint32_t i = 0; i < s_bake.lightmapWidth * s_bake.lightmapHeight; i++)
		sampleExists[i] = false;
	for (uint32_t tri = 0; tri < uint32_t(modelIndices.size() / 3); tri++) {
		// Skip emissive.
		const objzMaterial *mat = s_bake.triMaterials[tri];
		if (mat && (mat->emission[0] > 0.0f || mat->emission[1] > 0.0f || mat->emission[2] > 0.0f))
			continue;
		lm_context ctx;
		ctx.rasterizer.x = ctx.rasterizer.y = 0;
		lm_vec2 uvMin = lm_v2(FLT_MAX, FLT_MAX), uvMax = lm_v2(-FLT_MAX, -FLT_MAX);
		for (int i = 0; i < 3; i++) {
			const ModelVertex &vertex = modelVertices[modelIndices[tri * 3 + i]];
			ctx.triangle.p[i] = vertex.pos;
			ctx.triangle.uv[i].x = vertex.texcoord[2] * s_bake.lightmapWidth;
			ctx.triangle.uv[i].y = vertex.texcoord[3] * s_bake.lightmapHeight;
			// update bounds on lightmap
			uvMin = lm_min2(uvMin, ctx.triangle.uv[i]);
			uvMax = lm_max2(uvMax, ctx.triangle.uv[i]);
		}
		// Calculate area of interest (on lightmap) for conservative rasterization.
		lm_vec2 bbMin = lm_floor2(uvMin);
		lm_vec2 bbMax = lm_ceil2(uvMax);
		ctx.rasterizer.minx = ctx.rasterizer.x = bx::max((int)bbMin.x - 1, 0);
		ctx.rasterizer.miny = ctx.rasterizer.y = bx::max((int)bbMin.y - 1, 0);
		ctx.rasterizer.maxx = bx::min((int)bbMax.x + 1, (int)s_bake.lightmapWidth);
		ctx.rasterizer.maxy = bx::min((int)bbMax.y + 1, (int)s_bake.lightmapHeight);
		assert(ctx.rasterizer.minx <= ctx.rasterizer.maxx && ctx.rasterizer.miny <= ctx.rasterizer.maxy);
		if (lm_findFirstConservativeTriangleRasterizerPosition(&ctx)) {
			for (;;) {
				SampleLocation sample;
				sample.pos = ctx.sample.position;
				sample.normal = ctx.sample.direction;
				sample.uv[0] = ctx.rasterizer.x;
				sample.uv[1] = ctx.rasterizer.y;
				const uint32_t offset = sample.uv[0] + sample.uv[1] * s_bake.lightmapWidth;
				if (!sampleExists[offset]) {
					sampleExists[offset] = true;
					s_bake.sampleLocations.push_back(sample);
				}
				if (!lm_findNextConservativeTriangleRasterizerPosition(&ctx))
					break;
			}
		}
		s_bake.numTrianglesRasterized++;
		if (s_bake.cancelWorker)
			return false;
	}
	// Sort samples.
	s_bake.sampleLocationRanks.resize(s_bake.sampleLocations.size());
	for (uint32_t i = 0; i < s_bake.sampleLocationRanks.size(); i++)
		s_bake.sampleLocationRanks[i] = i;
	qsort(s_bake.sampleLocationRanks.data(), s_bake.sampleLocationRanks.size(), sizeof(uint32_t), compareSampleLocations);
	return true;
}

// https://en.wikipedia.org/wiki/Halton_sequence
static float halton(int index, int base)
{
	float result = 0;
	float f = 1;
	while (index > 0) {
		f /= base;
		result += f * (index % base);
		index = (int)bx::floor(index / (float)base);
	}
	return result;
}

// Precomputed Global Illumination in Frostbite (GDC 2018)
static bx::Vec3 randomDirHemisphere(int index, const float *offset, bx::Vec3 normal)
{
	// Generate 2 uniformly-distributed values in range 0 to 1
	float u = halton(index, 3);
	float v = halton(index, 5);
	// Apply per-texel randomization
	u = bx::fract(u + offset[0]);
	v = bx::fract(v + offset[1]);
	// Transform unit square sample to uniform hemisphere direction
	const float cosTheta = u * 2.0f - 1.0f;
	const float sinTheta = sqrtf(1.0f - cosTheta * cosTheta);
	const float vrad = v * bx::kPi2;
	const float sinPhi = bx::sin(vrad);
	const float cosPhi = bx::cos(vrad);
	bx::Vec3 dir(cosPhi * sinTheta, sinPhi * sinTheta, cosTheta);
	if (bx::dot(dir, normal) < 0.0f)
		return bx::neg(dir);
	return dir;
}

// https://github.com/aras-p/ToyPathTracer
static void bakeTraceRaysTask(uint32_t start, uint32_t end, uint32_t /*threadIndex*/, void * /*args*/)
{
	const float kNear = 0.01f * modelGetScale();
	// Initialize paths.
	for (uint32_t i = start; i < end; i++) {
		TexelData &texel = s_bake.texels[i];
		const SampleLocation &sample = s_bake.sampleLocations[s_bake.sampleLocationRanks[i]];
		bx::Vec3 color = bx::Vec3(0.0f);
		bx::Vec3 throughput = bx::Vec3(1.0f);
		bx::Vec3 rayOrigin = sample.pos;
		bx::Vec3 rayDir = randomDirHemisphere(texel.numPathsTraced, texel.randomOffset, sample.normal);
		for (int depth = 0; depth < s_bake.options.maxDepth; depth++) {
			RTCIntersectContext context;
			rtcInitIntersectContext(&context);
			alignas(16) RTCRayHit rh;
			rh.ray.org_x = rayOrigin.x;
			rh.ray.org_y = rayOrigin.y;
			rh.ray.org_z = rayOrigin.z;
			rh.ray.tnear = kNear;
			rh.ray.dir_x = rayDir.x;
			rh.ray.dir_y = rayDir.y;
			rh.ray.dir_z = rayDir.z;
			rh.ray.time = 0.0f;
			rh.ray.tfar = FLT_MAX;
			rh.ray.mask = UINT32_MAX;
			rh.ray.id = 0;
			rh.ray.flags = 0;
			rh.hit.primID = RTC_INVALID_GEOMETRY_ID;
			rh.hit.geomID = RTC_INVALID_GEOMETRY_ID;
			embree::Intersect1(s_bake.embreeScene, &context, &rh);
			s_bake.numRaysTraced++;
			texel.numPathsTraced++;
			if (rh.hit.geomID == RTC_INVALID_GEOMETRY_ID) {
				// Ray missed, use sky color.
				color = color + (throughput * s_bake.options.skyColor);
				break;
			}
			const uint32_t *indices = atlasGetIndices()->data();
			const ModelVertex *vertices = atlasGetVertices()->data();
			const ModelVertex &v0 = vertices[indices[rh.hit.primID * 3 + 0]];
			const ModelVertex &v1 = vertices[indices[rh.hit.primID * 3 + 1]];
			const ModelVertex &v2 = vertices[indices[rh.hit.primID * 3 + 2]];
			// we got a new ray bounced from the surface; recursively trace it
			bx::Vec3 diffuse = bx::Vec3(0.5f);
			bx::Vec3 emission(0.0f);
			const objzMaterial *mat = s_bake.triMaterials[rh.hit.primID];
			if (mat) {
				diffuse = bx::Vec3(mat->diffuse[0], mat->diffuse[1], mat->diffuse[2]);
				emission = bx::Vec3(mat->emission[0], mat->emission[1], mat->emission[2]);
				float uv[2];
				for (int j = 0; j < 2; j++)
					uv[j] = v0.texcoord[j] + (v1.texcoord[j] - v0.texcoord[j]) * rh.hit.u + (v2.texcoord[j] - v0.texcoord[j]) * rh.hit.v;
				bx::Vec3 texelColor;
				if (modelSampleMaterialDiffuse(mat, uv, &texelColor))
					diffuse = diffuse * texelColor;
				if (modelSampleMaterialEmission(mat, uv, &texelColor))
					emission = texelColor;
			}
			if (emission.x > 0.0f || emission.y > 0.0f || emission.z > 0.0f)
				color = color + (throughput * emission);
			else
				throughput = throughput * diffuse;
			// Russian Roulette
			// https://computergraphics.stackexchange.com/questions/2316/is-russian-roulette-really-the-answer
			const float p = bx::max(throughput.x, bx::max(throughput.y, throughput.z));
			if (bx::frnd(&s_rng) > p)
				break;
			throughput = throughput * (1.0f / p);
			if (depth + 1 < s_bake.options.maxDepth) {
				// Using barycentrics should be more precise than "origin + dir * rh.ray.tfar".
				rayOrigin = v0.pos + (v1.pos - v0.pos) * rh.hit.u + (v2.pos - v0.pos) * rh.hit.v;
				const bx::Vec3 normal = bx::normalize(bx::Vec3(rh.hit.Ng_x, rh.hit.Ng_y, rh.hit.Ng_z));
				rayDir = randomDirHemisphere(texel.numPathsTraced, texel.randomOffset, normal);
			}
		}
		texel.accumColor = texel.accumColor + color;
		texel.numColorSamples++;
	}
	// Copy texel data to lightmap.
	for (uint32_t i = start; i < end; i++) {
		TexelData &texel = s_bake.texels[i];
		const SampleLocation &sample = s_bake.sampleLocations[s_bake.sampleLocationRanks[i]];
		float *rgba = &s_bake.lightmapData[(sample.uv[0] + sample.uv[1] * s_bake.lightmapWidth) * 4];
		const float invn = 1.0f / (float)texel.numColorSamples;
		rgba[0] = texel.accumColor.x * invn;
		rgba[1] = texel.accumColor.y * invn;
		rgba[2] = texel.accumColor.z * invn;
		rgba[3] = 1.0f;
	}
}

static bool bakeTraceRays()
{
	s_bake.rng.reset();
	s_bake.numRaysTraced = 0;
	s_bake.lightmapData.resize(s_bake.lightmapWidth * s_bake.lightmapHeight * 4);
	memset(s_bake.lightmapData.data(), 0, s_bake.lightmapData.size() * sizeof(float));
	s_bake.updateData.resize(s_bake.lightmapWidth * s_bake.lightmapHeight * 4);
	s_bake.texels.resize(s_bake.sampleLocations.size());
	for (uint32_t i = 0; i < (uint32_t)s_bake.texels.size(); i++) {
		TexelData &texel = s_bake.texels[i];
		texel.accumColor = bx::Vec3(0.0f);
		texel.numColorSamples = 0;
		texel.numPathsTraced = 0;
		texel.randomOffset[0] = bx::frnd(&s_bake.rng);
		texel.randomOffset[1] = bx::frnd(&s_bake.rng);
	}
	enkiTaskSet *task = enkiCreateTaskSet(s_bake.taskScheduler, bakeTraceRaysTask);
	const clock_t start = clock();
	s_bake.samplesPerTexelCount = 0;
	for (;;) {
#if 0
		bakeTraceRaysTask(0, (int32_t)s_bake.sampleLocations.size(), 0, nullptr);
#else
		enkiAddTaskSetToPipe(s_bake.taskScheduler, task, 0, (int32_t)s_bake.sampleLocations.size());
		enkiWaitForTaskSet(s_bake.taskScheduler, task);
#endif
		if (s_bake.cancelWorker)
			break;
		s_bake.samplesPerTexelCount++;
		if (s_bake.stopWorker)
			break;
		// Handle lightmap updates.
		if (s_bake.updateStatus == UpdateStatus::Idle) {
			memcpy(s_bake.updateData.data(), s_bake.lightmapData.data(), s_bake.lightmapData.size() * sizeof(float));
			s_bake.updateStatus = UpdateStatus::Pending;
		}
	}
	enkiDeleteTaskSet(task);
	if (s_bake.cancelWorker)
		return false;
	const double elapsedSeconds = (clock() - start) / (double)CLOCKS_PER_SEC;
	printf("Finished tracing rays in %.2f seconds. %u samples per texel. %.2f Mrays/s.\n", elapsedSeconds, s_bake.samplesPerTexelCount.load(), s_bake.numRaysTraced / elapsedSeconds / 1000000.0);
	return true;
}

typedef bool (*FilterWritePredicate)(const float *data, uint32_t offset);
typedef bool (*FilterSamplePredicate)(const float *data, uint32_t offset, uint32_t chart);

static void bakeFilter(float *data, FilterWritePredicate writePredicate, FilterSamplePredicate samplePredicate)
{
	const int xOffsets[] = { -1, 0, 1, -1, 1, -1, 0, 1 };
	const int yOffsets[] = { -1, -1, -1, 0, 0, 1, 1, 1 };
	const uint32_t *atlasImage = atlasGetImage();
	for (uint32_t y = 0; y < s_bake.lightmapHeight; y++) {
		for (uint32_t x = 0; x < s_bake.lightmapWidth; x++) {
			const uint32_t offset = x + y * s_bake.lightmapWidth;
			if (!writePredicate(data, offset))
				continue;
			const uint32_t chart = atlasImage[offset] & xatlas::kImageChartIndexMask;
			float rgbSum[3] = { 0.0f, 0.0f, 0.0f }, n = 0;
			for (uint32_t si = 0; si < 8; si++) {
				const int sx = (int)x + xOffsets[si];
				const int sy = (int)y + yOffsets[si];
				if (sx < 0 || sy < 0 || sx >= (int)s_bake.lightmapWidth || sy >= (int)s_bake.lightmapHeight)
					continue; // Sample position is outside of atlas.
				const uint32_t sampleOffset = (uint32_t)sx + (uint32_t)sy * s_bake.lightmapWidth;
				if (atlasImage[offset] == 0)
					continue; // Empty space.
				if (!samplePredicate(data, sampleOffset, chart))
					continue;
				const float *rgba = &data[sampleOffset * 4];
				rgbSum[0] += (float)rgba[0];
				rgbSum[1] += (float)rgba[1];
				rgbSum[2] += (float)rgba[2];
				n++;
			}
			if (n != 0) {
				const float invn = 1.0f / (float)n;
				float *rgba = &data[offset * 4];
				rgba[0] = rgbSum[0] * invn;
				rgba[1] = rgbSum[1] * invn;
				rgba[2] = rgbSum[2] * invn;
				rgba[3] = 1.0f;
			}
		}
	}
}

static bool emptyFilterWritePredicate(const float *data, uint32_t offset)
{
	// Pixel in atlas image is a valid chart (but not bilinear or padding), but lightmap data was never sampled by the rasterizer (alpha is 0).
	const uint32_t imageData = atlasGetImage()[offset];
	const float *rgba = &data[offset * 4];
	return imageData & xatlas::kImageHasChartIndexBit && !(imageData & (xatlas::kImageIsBilinearBit | xatlas::kImageIsPaddingBit)) && rgba[3] <= 0.0f;
}

static bool emptyFilterSamplePredicate(const float * /*data*/, uint32_t offset, uint32_t chart)
{
	const uint32_t *atlasImage = atlasGetImage();
	if (atlasImage[offset] & (xatlas::kImageIsBilinearBit | xatlas::kImageIsPaddingBit))
		return false; // Don't sample from bilinear or padding.
	if ((atlasImage[offset] & xatlas::kImageChartIndexMask) != chart)
		return false; // Only sample from matching chart.
	return true;
}

static bool bilinearFilterWritePredicate(const float * /*data*/, uint32_t offset)
{
	// Only write to texels sampled by bilinear interpolation..
	return (atlasGetImage()[offset] & xatlas::kImageIsBilinearBit) != 0;
}

static bool bilinearFilterSamplePredicate(const float * /*data*/, uint32_t offset, uint32_t chart)
{
	const uint32_t *atlasImage = atlasGetImage();
	if (atlasImage[offset] & (xatlas::kImageIsBilinearBit | xatlas::kImageIsPaddingBit))
		return false; // Don't sample from bilinear or padding.
	if ((atlasImage[offset] & xatlas::kImageChartIndexMask) != chart)
		return false; // Only sample from matching chart.
	return true;
}

static void bakeFilter(float *data)
{
	bakeFilter(data, emptyFilterWritePredicate, emptyFilterSamplePredicate);
	bakeFilter(data, bilinearFilterWritePredicate, bilinearFilterSamplePredicate);
}

static void bakeWorkerThread()
{
	if (!bakeInitEmbree()) {
		s_bake.status = BakeStatus::Error;
		return;
	}
	s_bake.status = BakeStatus::Rasterizing;
	if (!bakeRasterize()) {
		s_bake.status = BakeStatus::Cancelled;
		return;
	}
	s_bake.status = BakeStatus::Tracing;
	if (!bakeTraceRays()) {
		s_bake.status = BakeStatus::Cancelled;
		return;
	}
	bakeFilter(s_bake.lightmapData.data());
	s_bake.status = BakeStatus::ThreadFinished;
}

static void shutdownWorkerThread()
{
	s_bake.cancelWorker = true;
	if (s_bake.workerThread) {
		if (s_bake.workerThread->joinable())
			s_bake.workerThread->join();
		delete s_bake.workerThread;
		s_bake.workerThread = nullptr;
	}
}

#if BX_ARCH_64BIT
#define OIDN_LIB "OpenImageDenoise.dll"

namespace oidn
{
	typedef OIDNDevice (*NewDeviceFunc)(OIDNDeviceType type);
	typedef void (*CommitDeviceFunc)(OIDNDevice device);
	typedef void (*ReleaseDeviceFunc)(OIDNDevice device);
	typedef void (*SetDevice1bFunc)(OIDNDevice device, const char* name, bool value);
	typedef void (*SetDeviceErrorFunctionFunc)(OIDNDevice device, OIDNErrorFunction func, void* userPtr);
	typedef OIDNFilter (*NewFilterFunc)(OIDNDevice device, const char* type);
	typedef void (*SetFilterProgressMonitorFunctionFunc)(OIDNFilter filter, OIDNProgressMonitorFunction func, void* userPtr);
	typedef void (*SetSharedFilterImageFunc)(OIDNFilter filter, const char* name, void* ptr, OIDNFormat format, size_t width, size_t height, size_t byteOffset, size_t bytePixelStride, size_t byteRowStride);
	typedef void (*SetFilter1bFunc)(OIDNFilter filter, const char* name, bool value);
	typedef void (*CommitFilterFunc)(OIDNFilter filter);
	typedef void (*ExecuteFilterFunc)(OIDNFilter filter);
	typedef void (*ReleaseFilterFunc)(OIDNFilter filter);
	NewDeviceFunc NewDevice;
	CommitDeviceFunc CommitDevice;
	ReleaseDeviceFunc ReleaseDevice;
	SetDevice1bFunc SetDevice1b;
	SetDeviceErrorFunctionFunc SetDeviceErrorFunction;
	NewFilterFunc NewFilter;
	SetFilterProgressMonitorFunctionFunc SetFilterProgressMonitorFunction;
	SetSharedFilterImageFunc SetSharedFilterImage;
	SetFilter1bFunc SetFilter1b;
	CommitFilterFunc CommitFilter;
	ExecuteFilterFunc ExecuteFilter;
	ReleaseFilterFunc ReleaseFilter;
};

static void bakeOidnError(void* /*userPtr*/, OIDNError code, const char* message)
{
	if (code == OIDN_ERROR_CANCELLED)
		return;
	fprintf(stderr, "OIDN error: %s\n", message);
	exit(EXIT_FAILURE);
}

static bool bakeOidnProgress(void* /*userPtr*/, double n)
{
	s_bake.denoiseProgress = n;
	return true;
}

static void bakeDenoiseThread()
{
	s_bake.denoiseProgress = 0;
	if (!s_bake.oidnLibrary) {
		s_bake.oidnLibrary = bx::dlopen(OIDN_LIB);
		if (!s_bake.oidnLibrary) {
			std::lock_guard<std::mutex> lock(s_bake.errorMessageMutex);
			bx::snprintf(s_bake.errorMessage, sizeof(s_bake.errorMessage), "OIDN not installed. Cannot open '%s'.", OIDN_LIB);
			fprintf(stderr, "OIDN not installed. Cannot open '%s'.\n", OIDN_LIB);
			s_bake.denoiseStatus = DenoiseStatus::Error;
			return;
		}
		oidn::NewDevice = (oidn::NewDeviceFunc)bx::dlsym(s_bake.oidnLibrary, "oidnNewDevice");
		oidn::CommitDevice = (oidn::CommitDeviceFunc)bx::dlsym(s_bake.oidnLibrary, "oidnCommitDevice");
		oidn::ReleaseDevice = (oidn::ReleaseDeviceFunc)bx::dlsym(s_bake.oidnLibrary, "oidnReleaseDevice");
		oidn::SetDevice1b = (oidn::SetDevice1bFunc)bx::dlsym(s_bake.oidnLibrary, "oidnSetDevice1b");
		oidn::SetDeviceErrorFunction = (oidn::SetDeviceErrorFunctionFunc)bx::dlsym(s_bake.oidnLibrary, "oidnSetDeviceErrorFunction");
		oidn::NewFilter = (oidn::NewFilterFunc)bx::dlsym(s_bake.oidnLibrary, "oidnNewFilter");
		oidn::SetFilterProgressMonitorFunction = (oidn::SetFilterProgressMonitorFunctionFunc)bx::dlsym(s_bake.oidnLibrary, "oidnSetFilterProgressMonitorFunction");
		oidn::SetSharedFilterImage = (oidn::SetSharedFilterImageFunc)bx::dlsym(s_bake.oidnLibrary, "oidnSetSharedFilterImage");
		oidn::SetFilter1b = (oidn::SetFilter1bFunc)bx::dlsym(s_bake.oidnLibrary, "oidnSetFilter1b");
		oidn::CommitFilter = (oidn::CommitFilterFunc)bx::dlsym(s_bake.oidnLibrary, "oidnCommitFilter");
		oidn::ExecuteFilter = (oidn::ExecuteFilterFunc)bx::dlsym(s_bake.oidnLibrary, "oidnExecuteFilter");
		oidn::ReleaseFilter = (oidn::ReleaseFilterFunc)bx::dlsym(s_bake.oidnLibrary, "oidnReleaseFilter");
	}
	s_bake.denoisedLightmapData.resize(s_bake.lightmapWidth * s_bake.lightmapHeight * 4);
	OIDNDevice device = oidn::NewDevice(OIDN_DEVICE_TYPE_DEFAULT);
	if (!device) {
		std::lock_guard<std::mutex> lock(s_bake.errorMessageMutex);
		bx::snprintf(s_bake.errorMessage, sizeof(s_bake.errorMessage), "Error creating OIDN device");
		fprintf(stderr, "Error creating OIDN device\n");
		s_bake.denoiseStatus = DenoiseStatus::Error;
		return;
	}
	oidn::SetDeviceErrorFunction(device, bakeOidnError, nullptr);
	oidn::SetDevice1b(device, "setAffinity", false);
	oidn::CommitDevice(device);
	OIDNFilter filter = oidn::NewFilter(device, "RTLightmap");
	oidn::SetFilterProgressMonitorFunction(filter, bakeOidnProgress, nullptr);
	oidn::SetSharedFilterImage(filter, "color", s_bake.lightmapData.data(), OIDN_FORMAT_FLOAT3, s_bake.lightmapWidth, s_bake.lightmapHeight, 0, sizeof(float) * 4, 0);
	oidn::SetSharedFilterImage(filter, "output", s_bake.denoisedLightmapData.data(), OIDN_FORMAT_FLOAT3, s_bake.lightmapWidth, s_bake.lightmapHeight, 0, sizeof(float) * 4, 0);
	oidn::CommitFilter(filter);
	oidn::ExecuteFilter(filter);
	oidn::ReleaseFilter(filter);
	oidn::ReleaseDevice(device);
	// Copy alpha channel.
	for (uint32_t i = 0; i < s_bake.lightmapWidth * s_bake.lightmapHeight; i++)
		s_bake.denoisedLightmapData[i * 4 + 3] = s_bake.lightmapData[i * 4 + 3];
	bakeFilter(s_bake.denoisedLightmapData.data());
	s_bake.denoiseStatus = DenoiseStatus::ThreadFinished;
}

static void bakeDenoise()
{
	if (s_bake.status != BakeStatus::Finished)
		return;
	s_bake.denoiseStatus = DenoiseStatus::Working;
	s_bake.denoiseThread = new std::thread(bakeDenoiseThread);
}

static void bakeShutdownDenoiseThread()
{
	if (s_bake.denoiseThread) {
		if (s_bake.denoiseThread->joinable())
			s_bake.denoiseThread->join();
		delete s_bake.denoiseThread;
		s_bake.denoiseThread = nullptr;
	}
}
#endif

void bakeInit()
{
	s_bake.status = BakeStatus::Idle;
	s_bake.denoiseStatus = DenoiseStatus::Idle;
	s_bake.updateStatus = UpdateStatus::Idle;
	s_bake.stopWorker = false;
	s_bake.taskScheduler = enkiNewTaskScheduler();
	enkiInitTaskSchedulerNumThreads(s_bake.taskScheduler, std::thread::hardware_concurrency());
}

void bakeShutdown()
{
#if BX_ARCH_64BIT
	bakeShutdownDenoiseThread();
#endif
	shutdownWorkerThread();
	enkiDeleteTaskScheduler(s_bake.taskScheduler);
	if (s_bake.embreeDevice) {
		embree::ReleaseGeometry(s_bake.embreeGeometry);
		embree::ReleaseScene(s_bake.embreeScene);
		embree::ReleaseDevice(s_bake.embreeDevice);
		s_bake.embreeDevice = nullptr;
	}
	if (s_bake.embreeLibrary)
		bx::dlclose(s_bake.embreeLibrary);
	if (s_bake.oidnLibrary)
		bx::dlclose(s_bake.oidnLibrary);
	if (s_bake.initialized) { 
		bgfx::destroy(s_bake.lightmap);
		bgfx::destroy(s_bake.denoisedLightmap);
	}
}

void bakeExecute()
{
	if (!(s_bake.status == BakeStatus::Idle || s_bake.status == BakeStatus::Finished || s_bake.status == BakeStatus::Cancelled || s_bake.status == BakeStatus::Error))
		return;
	{
		std::lock_guard<std::mutex> lock(s_bake.errorMessageMutex);
		s_bake.errorMessage[0] = 0;
	}
	// Re-create lightmap if atlas resolution has changed.
	const bool lightmapResolutionChanged = s_bake.lightmapWidth != atlasGetWidth() || s_bake.lightmapHeight != atlasGetHeight();
	if (!s_bake.initialized || lightmapResolutionChanged) {
		s_bake.lightmapWidth = atlasGetWidth();
		s_bake.lightmapHeight = atlasGetHeight();
		if (s_bake.initialized) {
			bgfx::destroy(s_bake.lightmap);
			bgfx::destroy(s_bake.denoisedLightmap);
		}
		s_bake.lightmap = bgfx::createTexture2D((uint16_t)s_bake.lightmapWidth, (uint16_t)s_bake.lightmapHeight, false, 1, bgfx::TextureFormat::RGBA32F, BGFX_SAMPLER_UVW_CLAMP);
		s_bake.denoisedLightmap = bgfx::createTexture2D((uint16_t)s_bake.lightmapWidth, (uint16_t)s_bake.lightmapHeight, false, 1, bgfx::TextureFormat::RGBA32F, BGFX_SAMPLER_UVW_CLAMP);
	}
	s_bake.initialized = true;
	g_options.shadeMode = ShadeMode::Lightmap;
	s_bake.status = BakeStatus::InitEmbree;
	s_bake.denoiseStatus = DenoiseStatus::Idle;
	s_bake.updateStatus = UpdateStatus::Idle;
	s_bake.stopWorker = false;
	s_bake.cancelWorker = false;
	s_bake.workerThread = new std::thread(bakeWorkerThread);
}

void bakeFrame(uint32_t frameNo)
{
	if (s_bake.status == BakeStatus::ThreadFinished) {
		shutdownWorkerThread();
		// Do a final update of the lightmap texture with the dilated result.
		bgfx::updateTexture2D(s_bake.lightmap, 0, 0, 0, 0, (uint16_t)s_bake.lightmapWidth, (uint16_t)s_bake.lightmapHeight, bgfx::makeRef(s_bake.lightmapData.data(), (uint32_t)s_bake.lightmapData.size() * sizeof(float)));
		s_bake.status = BakeStatus::Finished;
	} else if ((s_bake.status == BakeStatus::Cancelled || s_bake.status == BakeStatus::Error) && g_options.shadeMode == ShadeMode::Lightmap) {
		// Executing bake sets shade mode to lightmap. Set it back to charts if bake was cancelled or there was an error.
		g_options.shadeMode = ShadeMode::Charts;
	}
	if (s_bake.denoiseStatus == DenoiseStatus::ThreadFinished) {
#if BX_ARCH_64BIT
		bakeShutdownDenoiseThread();
#endif
		bgfx::updateTexture2D(s_bake.denoisedLightmap, 0, 0, 0, 0, (uint16_t)s_bake.lightmapWidth, (uint16_t)s_bake.lightmapHeight, bgfx::makeRef(s_bake.denoisedLightmapData.data(), (uint32_t)s_bake.denoisedLightmapData.size() * sizeof(float)));
		s_bake.denoiseStatus = DenoiseStatus::Finished;
	}
	if (s_bake.updateStatus == UpdateStatus::Pending) {
		bgfx::updateTexture2D(s_bake.lightmap, 0, 0, 0, 0, (uint16_t)s_bake.lightmapWidth, (uint16_t)s_bake.lightmapHeight, bgfx::makeRef(s_bake.updateData.data(), (uint32_t)s_bake.updateData.size() * sizeof(float)));
		s_bake.updateStatus = UpdateStatus::Updating;
		s_bake.updateFinishedFrameNo = frameNo + 2;
	} else if (s_bake.updateStatus == UpdateStatus::Updating) {
		if (frameNo == s_bake.updateFinishedFrameNo)
			s_bake.updateStatus = UpdateStatus::Idle;
	}
}

void bakeClear()
{
	s_bake.status = BakeStatus::Idle;
	s_bake.denoiseStatus = DenoiseStatus::Idle;
	g_options.shadeMode = atlasIsReady() ? ShadeMode::Charts : ShadeMode::Flat;
	if (s_bake.embreeDevice) {
		embree::ReleaseGeometry(s_bake.embreeGeometry);
		embree::ReleaseScene(s_bake.embreeScene);
		embree::ReleaseDevice(s_bake.embreeDevice);
		s_bake.embreeDevice = nullptr;
	}
	s_bake.sampleLocations.clear();
	// Need to lookup material from triangle index when processing ray hits.
	const objzModel *model = modelGetData();
	if (model) {
		s_bake.triMaterials.resize(model->numIndices / 3);
		uint32_t triIndex = 0;
		for (uint32_t i = 0; i < model->numMeshes; i++) {
			const objzMesh &mesh = model->meshes[i];
			for (uint32_t j = 0; j < mesh.numIndices / 3; j++)
				s_bake.triMaterials[triIndex++] = mesh.materialIndex == -1 ? nullptr : &model->materials[mesh.materialIndex];
		}
	}
}

void bakeShowGuiOptions()
{
	const ImVec2 buttonSize(ImVec2(ImGui::GetContentRegionAvailWidth() * 0.3f, 0.0f));
	const ImVec4 red(1.0f, 0.0f, 0.0f, 1.0f);
	ImGui::Text(ICON_FA_LIGHTBULB_O " Lightmap");
	ImGui::Spacing();
	if (atlasGetCount() > 1) {
		ImGui::PushStyleColor(ImGuiCol_Text, red);
		ImGui::Text("Baking doesn't support multiple atlases");
		ImGui::PopStyleColor();
		return;
	}
#if BX_ARCH_32BIT
	ImGui::PushStyleColor(ImGuiCol_Text, red);
	ImGui::Text("Baking not supported in 32-bit build");
	ImGui::PopStyleColor();
#else
	if ((s_bake.status == BakeStatus::Idle || s_bake.status == BakeStatus::Finished || s_bake.status == BakeStatus::Cancelled || s_bake.status == BakeStatus::Error) && (s_bake.denoiseStatus == DenoiseStatus::Idle || s_bake.denoiseStatus == DenoiseStatus::Finished || s_bake.denoiseStatus == DenoiseStatus::Error)) {
		if (ImGui::Button(ICON_FA_BOLT " Bake", buttonSize))
			bakeExecute();
		if (s_bake.status == BakeStatus::Finished && (s_bake.denoiseStatus == DenoiseStatus::Idle || s_bake.denoiseStatus == DenoiseStatus::Error)) {
			ImGui::SameLine();
			if (ImGui::Button(ICON_FA_VOLUME_OFF " Denoise", buttonSize))
				bakeDenoise();
		}
		ImGui::Columns(2, nullptr, false);
		guiColumnColorEdit("Sky color", "##skyColor", &s_bake.options.skyColor.x);
		guiColumnSliderInt("Max depth", "##maxDepth", &s_bake.options.maxDepth, 1, kMaxDepth);
		ImGui::Columns(1);
		std::lock_guard<std::mutex> lock(s_bake.errorMessageMutex);
		if (s_bake.errorMessage[0]) {
			ImGui::PushStyleColor(ImGuiCol_Text, red);
			ImGui::Text("%s", s_bake.errorMessage);
			ImGui::PopStyleColor();
		}
	} else if (s_bake.status == BakeStatus::InitEmbree) {
		ImGui::AlignTextToFramePadding();
		ImGui::Text("Initializing Embree");
		ImGui::SameLine();
		ImGui::Spinner("##rasterSpinner");
	} else if (s_bake.status == BakeStatus::Rasterizing) {
		ImGui::AlignTextToFramePadding();
		ImGui::Text("Rasterizing");
		ImGui::SameLine();
		ImGui::Spinner("##rasterSpinner");
		ImGui::ProgressBar(s_bake.numTrianglesRasterized / float(modelGetData()->numIndices / 3));
		if (ImGui::Button(ICON_FA_TIMES " Cancel"))
			shutdownWorkerThread();
	} else if (s_bake.status == BakeStatus::Tracing) {
		ImGui::AlignTextToFramePadding();
		ImGui::Text("Tracing rays [%u samples]", s_bake.samplesPerTexelCount.load());
		ImGui::SameLine();
		ImGui::Spinner("##traceSpinner");
		if (ImGui::Button(ICON_FA_STOP " Stop", buttonSize))
			s_bake.stopWorker = true;
		ImGui::SameLine();
		if (ImGui::Button(ICON_FA_TIMES " Cancel", buttonSize))
			shutdownWorkerThread();
	}
	if (s_bake.denoiseStatus == DenoiseStatus::Working) {
		ImGui::AlignTextToFramePadding();
		ImGui::Text("Denoising");
		ImGui::SameLine();
		ImGui::Spinner("##denoiseSpinner");
		ImGui::ProgressBar((float)s_bake.denoiseProgress);
	}
#endif
}

void bakeShowGuiWindow()
{
	if (s_bake.status == BakeStatus::Idle || s_bake.status == BakeStatus::InitEmbree || s_bake.status == BakeStatus::Error || !g_options.showLightmapWindow)
		return;
	if (ImGui::Begin("Lightmap", &g_options.showLightmapWindow, ImGuiWindowFlags_HorizontalScrollbar)) {
		ImGui::Checkbox("Fit to window", &s_bake.options.fitToWindow);
		if (!s_bake.options.fitToWindow) {
			ImGui::SameLine();
			ImGui::PushItemWidth(50.0f);
			ImGui::InputFloat("Scale", &s_bake.options.scale);
			ImGui::PopItemWidth();
		}
		GuiTexture texture;
		texture.bgfx.handle = bakeGetLightmap();
		texture.bgfx.flags = GuiTextureFlags::PointSampler;
		if (s_bake.options.fitToWindow) {
			// Fit to content while maintaining aspect ratio.
			ImVec2 size((float)s_bake.lightmapWidth, (float)s_bake.lightmapHeight);
			const float scale = bx::min(ImGui::GetContentRegionAvail().x / size.x, ImGui::GetContentRegionAvail().y / size.y);
			size.x *= scale;
			size.y *= scale;
			ImGui::Image(texture.imgui, size);
		}
		else
			ImGui::Image(texture.imgui, ImVec2((float)s_bake.lightmapWidth * s_bake.options.scale, (float)s_bake.lightmapHeight * s_bake.options.scale));
		ImGui::End();
	}
}

bgfx::TextureHandle bakeGetLightmap()
{
	return (s_bake.denoiseStatus == DenoiseStatus::Finished && g_options.useDenoisedLightmap) ? s_bake.denoisedLightmap : s_bake.lightmap;
}

uint32_t bakeGetLightmapSamplerFlags()
{
	uint32_t flags = BGFX_SAMPLER_UVW_CLAMP;
	if (g_options.lightmapPointSampling)
		flags |= BGFX_SAMPLER_POINT;
	return flags;
}

bool bakeIsLightmapReady()
{
	return !(s_bake.status == BakeStatus::Idle || s_bake.status == BakeStatus::InitEmbree || s_bake.status == BakeStatus::Error);
}

bool bakeIsDenoised()
{
	return s_bake.denoiseStatus == DenoiseStatus::Finished;
}
