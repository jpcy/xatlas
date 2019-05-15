/*
xatlas
https://github.com/jpcy/xatlas
Copyright (c) 2018 Jonathan Young

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#include <atomic>
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
#include "shaders/shared.h"
#include "viewer.h"

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
				bx::Vec3 p0 = ctx->triangle.p[0];
				bx::Vec3 p1 = ctx->triangle.p[1];
				bx::Vec3 p2 = ctx->triangle.p[2];
				bx::Vec3 v1 = bx::sub(p1, p0);
				bx::Vec3 v2 = bx::sub(p2, p0);
				ctx->sample.position = bx::add(p0, bx::add(bx::mul(v2, uv.x), bx::mul(v1, uv.y)));
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

static void lmImageDilate(const float *image, float *outImage, int w, int h, int c)
{
	assert(c > 0 && c <= 4);
	for (int y = 0; y < h; y++)
	{
		for (int x = 0; x < w; x++)
		{
			float color[4];
			bool valid = false;
			for (int i = 0; i < c; i++)
			{
				color[i] = image[(y * w + x) * c + i];
				valid |= color[i] > 0.0f;
			}
			if (!valid)
			{
				int n = 0;
				const int dx[] = { -1, 0, 1,  0 };
				const int dy[] = {  0, 1, 0, -1 };
				for (int d = 0; d < 4; d++)
				{
					int cx = x + dx[d];
					int cy = y + dy[d];
					if (cx >= 0 && cx < w && cy >= 0 && cy < h)
					{
						float dcolor[4];
						bool dvalid = false;
						for (int i = 0; i < c; i++)
						{
							dcolor[i] = image[(cy * w + cx) * c + i];
							dvalid |= dcolor[i] > 0.0f;
						}
						if (dvalid)
						{
							for (int i = 0; i < c; i++)
								color[i] += dcolor[i];
							n++;
						}
					}
				}
				if (n)
				{
					float in = 1.0f / n;
					for (int i = 0; i < c; i++)
						color[i] *= in;
				}
			}
			for (int i = 0; i < c; i++)
				outImage[(y * w + x) * c + i] = color[i];
		}
	}
}

struct BakeStatus
{
	enum Enum
	{
		Idle,
		InitEmbree,
		Rasterizing,
		Tracing,
		Denoising,
		ThreadFinished,
		Finished,
		Error
	};

	bool operator==(Enum value)
	{
		m_lock.lock();
		const bool result = m_value == value;
		m_lock.unlock();
		return result;
	}

	bool operator!=(Enum value)
	{
		m_lock.lock();
		const bool result = m_value != value;
		m_lock.unlock();
		return result;
	}

	BakeStatus &operator=(Enum value)
	{
		m_lock.lock();
		m_value = value;
		m_lock.unlock();
		return *this;
	}

private:
	std::mutex m_lock;
	Enum m_value = Idle;
};

struct UpdateStatus
{
	enum Enum
	{
		Idle,
		Pending,
		Updating
	};

	bool operator==(Enum value)
	{
		m_lock.lock();
		const bool result = m_value == value;
		m_lock.unlock();
		return result;
	}

	bool operator!=(Enum value)
	{
		m_lock.lock();
		const bool result = m_value != value;
		m_lock.unlock();
		return result;
	}

	UpdateStatus &operator=(Enum value)
	{
		m_lock.lock();
		m_value = value;
		m_lock.unlock();
		return *this;
	}

private:
	std::mutex m_lock;
	Enum m_value = Idle;
};

struct BakeOptions
{
	bool fitToWindow = true;
	bool denoise = true;
	bool sky = true;
	bx::Vec3 skyColor = bx::Vec3(0.5f, 0.5f, 0.5f);
	int numSamples = 16;
	int numBounces = 1;
};

struct SampleLocation
{
	bx::Vec3 pos;
	bx::Vec3 normal;
	uint32_t uv[2];
};

struct
{
	char errorMessage[256];
	bool initialized = false;
	BakeStatus status;
	BakeOptions options;
	void *embreeLibrary = nullptr;
	RTCDevice embreeDevice = nullptr;
	RTCScene embreeScene = nullptr;
	RTCGeometry embreeGeometry = nullptr;
	void *oidnLibrary = nullptr;
	// lightmap
	bgfx::TextureHandle lightmap;
	uint32_t lightmapWidth, lightmapHeight;
	// worker thread
	std::thread *workerThread = nullptr;
	bool cancelWorker = false;
	std::mutex cancelWorkerMutex;
	std::vector<SampleLocation> sampleLocations;
	std::vector<float> lightmapData;
	std::vector<float> denoisedLightmapData;
	std::vector<float> dilatedLightmapData;
	std::atomic<uint32_t> numTrianglesRasterized;
	std::atomic<uint32_t> numSampleLocationsProcessed;
	double denoiseProgress;
	bool denoiseSucceeded;
	bx::RngMwc rng;
	// lightmap update
	clock_t lastUpdateTime = 0;
	const double updateIntervalMs = 50;
	std::vector<float> updateData;
	UpdateStatus updateStatus;
	uint32_t updateFinishedFrameNo;
}
s_bake;

#define EMBREE_LIB "embree3.dll"

namespace embree
{
	typedef RTCDevice (*NewDeviceFunc)(const char* config);
	typedef void (*ReleaseDeviceFunc)(RTCDevice device);
	typedef void (*SetDeviceErrorFunctionFunc)(RTCDevice device, RTCErrorFunction error, void* userPtr);
	typedef RTCScene (*NewSceneFunc)(RTCDevice device);
	typedef void (*ReleaseSceneFunc)(RTCScene scene);
	typedef unsigned int (*AttachGeometryFunc)(RTCScene scene, RTCGeometry geometry);
	typedef void (*CommitSceneFunc)(RTCScene scene);
	typedef RTCGeometry (*NewGeometryFunc)(RTCDevice device, enum RTCGeometryType type);
	typedef void (*ReleaseGeometryFunc)(RTCGeometry geometry);
	typedef void (*SetSharedGeometryBufferFunc)(RTCGeometry geometry, enum RTCBufferType type, unsigned int slot, enum RTCFormat format, const void* ptr, size_t byteOffset, size_t byteStride, size_t itemCount);
	typedef void (*CommitGeometryFunc)(RTCGeometry geometry);
	typedef void (*Occluded1Func)(RTCScene scene, struct RTCIntersectContext* context, struct RTCRay* ray);
	typedef void (*Occluded4Func)(const int* valid, RTCScene scene, struct RTCIntersectContext* context, struct RTCRay4* ray);
	typedef void (*Occluded8Func)(const int* valid, RTCScene scene, struct RTCIntersectContext* context, struct RTCRay8* ray);
	typedef void (*Occluded16Func)(const int* valid, RTCScene scene, struct RTCIntersectContext* context, struct RTCRay16* ray);
	NewDeviceFunc NewDevice;
	ReleaseDeviceFunc ReleaseDevice;
	SetDeviceErrorFunctionFunc SetDeviceErrorFunction;
	NewSceneFunc NewScene;
	ReleaseSceneFunc ReleaseScene;
	AttachGeometryFunc AttachGeometry;
	CommitSceneFunc CommitScene;
	NewGeometryFunc NewGeometry;
	ReleaseGeometryFunc ReleaseGeometry;
	SetSharedGeometryBufferFunc SetSharedGeometryBuffer;
	CommitGeometryFunc CommitGeometry;
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
			bx::snprintf(s_bake.errorMessage, sizeof(s_bake.errorMessage), "Embree not installed. Cannot open '%s'.", EMBREE_LIB);
			return false;
		}
		embree::NewDevice = (embree::NewDeviceFunc)bx::dlsym(s_bake.embreeLibrary, "rtcNewDevice");
		embree::ReleaseDevice = (embree::ReleaseDeviceFunc)bx::dlsym(s_bake.embreeLibrary, "rtcReleaseDevice");
		embree::SetDeviceErrorFunction = (embree::SetDeviceErrorFunctionFunc)bx::dlsym(s_bake.embreeLibrary, "rtcSetDeviceErrorFunction");
		embree::NewScene = (embree::NewSceneFunc)bx::dlsym(s_bake.embreeLibrary, "rtcNewScene");
		embree::ReleaseScene = (embree::ReleaseSceneFunc)bx::dlsym(s_bake.embreeLibrary, "rtcReleaseScene");
		embree::AttachGeometry = (embree::AttachGeometryFunc)bx::dlsym(s_bake.embreeLibrary, "rtcAttachGeometry");
		embree::CommitScene = (embree::CommitSceneFunc)bx::dlsym(s_bake.embreeLibrary, "rtcCommitScene");
		embree::NewGeometry = (embree::NewGeometryFunc)bx::dlsym(s_bake.embreeLibrary, "rtcNewGeometry");
		embree::ReleaseGeometry = (embree::ReleaseGeometryFunc)bx::dlsym(s_bake.embreeLibrary, "rtcReleaseGeometry");
		embree::SetSharedGeometryBuffer = (embree::SetSharedGeometryBufferFunc)bx::dlsym(s_bake.embreeLibrary, "rtcSetSharedGeometryBuffer");
		embree::CommitGeometry = (embree::CommitGeometryFunc)bx::dlsym(s_bake.embreeLibrary, "rtcCommitGeometry");
		embree::Occluded1 = (embree::Occluded1Func)bx::dlsym(s_bake.embreeLibrary, "rtcOccluded1");
		embree::Occluded4 = (embree::Occluded4Func)bx::dlsym(s_bake.embreeLibrary, "rtcOccluded4");
		embree::Occluded8 = (embree::Occluded8Func)bx::dlsym(s_bake.embreeLibrary, "rtcOccluded8");
		embree::Occluded16 = (embree::Occluded16Func)bx::dlsym(s_bake.embreeLibrary, "rtcOccluded16");
	}
	if (!s_bake.embreeDevice) {
		s_bake.embreeDevice = embree::NewDevice(nullptr);
		if (!s_bake.embreeDevice) {
			bx::snprintf(s_bake.errorMessage, sizeof(s_bake.errorMessage), "Error creating Embree device");
			return false;
		}
		embree::SetDeviceErrorFunction(s_bake.embreeDevice, bakeEmbreeError, nullptr);
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

static bool bakeIsWorkerThreadCancelled()
{
	std::lock_guard<std::mutex> lock(s_bake.cancelWorkerMutex);
	return s_bake.cancelWorker;
}

static bool bakeRasterize()
{
	// Only need to do this once, unless the atlas or models changes - bakeClear will be called.
	if (!s_bake.sampleLocations.empty())
		return true;
	s_bake.numTrianglesRasterized = 0;
	std::vector<ModelVertex> &modelVertices = *atlasGetVertices();
	std::vector<uint32_t> &modelIndices = *atlasGetIndices();
	for (uint32_t tri = 0; tri < uint32_t(modelIndices.size() / 3); tri++) {
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
				s_bake.sampleLocations.push_back(sample);
				if (!lm_findNextConservativeTriangleRasterizerPosition(&ctx))
					break;
			}
		}
		s_bake.numTrianglesRasterized++;
		if (bakeIsWorkerThreadCancelled())
			return false;
	}
	return true;
}

struct Occluded4Functor
{
	void operator()(const int* valid, RTCScene scene, struct RTCIntersectContext* context, struct RTCRay4* ray)
	{
		embree::Occluded4(valid, scene, context, ray);
	}
};

struct Occluded8Functor
{
	void operator()(const int* valid, RTCScene scene, struct RTCIntersectContext* context, struct RTCRay8* ray)
	{
		embree::Occluded8(valid, scene, context, ray);
	}
};

struct Occluded16Functor
{
	void operator()(const int* valid, RTCScene scene, struct RTCIntersectContext* context, struct RTCRay16* ray)
	{
		embree::Occluded16(valid, scene, context, ray);
	}
};

template<typename Ray, typename Occluded, int N>
static int bakeCalculateRayPacketVisibility(RTCIntersectContext *context, const SampleLocation &sample, float near)
{
	Ray ray;
	int valid[N];
	for (int i = 0; i < N; i++) {
		ray.org_x[i] = sample.pos.x;
		ray.org_y[i] = sample.pos.y;
		ray.org_z[i] = sample.pos.z;
		const bx::Vec3 dir = bx::randUnitHemisphere(&s_bake.rng, sample.normal);
		ray.dir_x[i] = dir.x;
		ray.dir_y[i] = dir.y;
		ray.dir_z[i] = dir.z;
		ray.tnear[i] = near;
		ray.tfar[i] = FLT_MAX;
		ray.flags[i] = 0;
		valid[i] = -1;
	}
	Occluded occluded;
	occluded(valid, s_bake.embreeScene, context, &ray);
	int visible = 0;
	for (int i = 0; i < N; i++) {
		if (ray.tfar[i] > 0.0f)
			visible++;
	}
	return visible;
}

static bool bakeTraceRays()
{
	s_bake.rng.reset();
	s_bake.numSampleLocationsProcessed = 0;
	s_bake.lightmapData.resize(s_bake.lightmapWidth * s_bake.lightmapHeight * 4);
	memset(s_bake.lightmapData.data(), 0, s_bake.lightmapData.size() * sizeof(float));
	s_bake.updateData.resize(s_bake.lightmapWidth * s_bake.lightmapHeight * 4);
	RTCIntersectContext context;
	const float near = 0.01f * modelGetScale();
	for (uint32_t si = 0; si < (uint32_t)s_bake.sampleLocations.size(); si++) {
		rtcInitIntersectContext(&context);
		const SampleLocation &sample = s_bake.sampleLocations[si];
		float visibility = 0.0f;
		if (s_bake.options.numSamples == 4)
			visibility = bakeCalculateRayPacketVisibility<RTCRay4, Occluded4Functor, 4>(&context, sample, near) / 4.0f;
		else if (s_bake.options.numSamples == 8)
			visibility = bakeCalculateRayPacketVisibility<RTCRay8, Occluded8Functor, 8>(&context, sample, near) / 8.0f;
		else if (s_bake.options.numSamples == 16)
			visibility = bakeCalculateRayPacketVisibility<RTCRay16, Occluded16Functor, 16>(&context, sample, near) / 16.0f;
		else if (s_bake.options.numSamples == 256) {
			int visible = 0;
			for (int i = 0; i < 16; i++)
				visible += bakeCalculateRayPacketVisibility<RTCRay16, Occluded16Functor, 16>(&context, sample, near);
			visibility = visible / 256.0f;
		}
		else {
			RTCRay ray;
			ray.org_x = sample.pos.x;
			ray.org_y = sample.pos.y;
			ray.org_z = sample.pos.z;
			ray.dir_x = sample.normal.x;
			ray.dir_y = sample.normal.y;
			ray.dir_z = sample.normal.z;
			ray.tnear = near;
			ray.tfar = FLT_MAX;
			ray.flags = 0;
			embree::Occluded1(s_bake.embreeScene, &context, &ray);
			if (ray.tfar > 0.0f)
				visibility = 1.0f;
		}
		if (visibility > 0.0f) {
			float *rgba = &s_bake.lightmapData[(sample.uv[0] + sample.uv[1] * s_bake.lightmapWidth) * 4];
			rgba[0] = s_bake.options.skyColor.x * visibility;
			rgba[1] = s_bake.options.skyColor.y * visibility;
			rgba[2] = s_bake.options.skyColor.z * visibility;
			rgba[3] = 1.0f;
		}
		s_bake.numSampleLocationsProcessed++;
		if (bakeIsWorkerThreadCancelled())
			return false;
		// Handle lightmap updates.
		const double elapsedMs = (clock() - s_bake.lastUpdateTime) * 1000.0 / CLOCKS_PER_SEC;
		if (elapsedMs >= s_bake.updateIntervalMs) {
			if (s_bake.updateStatus == UpdateStatus::Idle) {
				memcpy(s_bake.updateData.data(), s_bake.lightmapData.data(), s_bake.lightmapData.size() * sizeof(float));
				s_bake.updateStatus = UpdateStatus::Pending;
			}
			s_bake.lastUpdateTime = clock();
		}
	}
	return true;
}

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
	if (bakeIsWorkerThreadCancelled())
		return false;
	s_bake.denoiseProgress = n;
	return true;
}

static void bakeDenoise()
{
	s_bake.denoiseSucceeded = false;
	s_bake.denoiseProgress = 0;
	if (!s_bake.options.denoise)
		return;
	if (!s_bake.oidnLibrary) {
		s_bake.oidnLibrary = bx::dlopen(OIDN_LIB);
		if (!s_bake.oidnLibrary) {
			fprintf(stderr, "OIDN not installed. Cannot open '%s'.", OIDN_LIB);
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
		fprintf(stderr, "Error creating OIDN device\n");
		return;
	}
	// Dilate into temp buffer.
	std::vector<float> temp;
	temp.resize(s_bake.lightmapWidth * s_bake.lightmapHeight * 4);
	lmImageDilate(s_bake.lightmapData.data(), temp.data(), (int)s_bake.lightmapWidth, (int)s_bake.lightmapHeight, 4);
	oidn::SetDeviceErrorFunction(device, bakeOidnError, nullptr);
	oidn::SetDevice1b(device, "setAffinity", false);
	oidn::CommitDevice(device);
	OIDNFilter filter = oidn::NewFilter(device, "RT");
	oidn::SetFilterProgressMonitorFunction(filter, bakeOidnProgress, nullptr);
	oidn::SetSharedFilterImage(filter, "color", temp.data(), OIDN_FORMAT_FLOAT3, s_bake.lightmapWidth, s_bake.lightmapHeight, 0, sizeof(float) * 4, 0);
	oidn::SetSharedFilterImage(filter, "output", s_bake.denoisedLightmapData.data(), OIDN_FORMAT_FLOAT3, s_bake.lightmapWidth, s_bake.lightmapHeight, 0, sizeof(float) * 4, 0);
	oidn::CommitFilter(filter);
	oidn::ExecuteFilter(filter);
	oidn::ReleaseFilter(filter);
	oidn::ReleaseDevice(device);
	if (bakeIsWorkerThreadCancelled())
		return;
	// Copy alpha channel.
	for (uint32_t i = 0; i < s_bake.lightmapWidth * s_bake.lightmapHeight; i++)
		s_bake.denoisedLightmapData[i * 4 + 3] = s_bake.lightmapData[i * 4 + 3];
	s_bake.denoiseSucceeded = true;
}

static void bakeWorkerThread()
{
	if (!bakeInitEmbree()) {
		s_bake.status = BakeStatus::Error;
		return;
	}
	s_bake.status = BakeStatus::Rasterizing;
	if (!bakeRasterize()) {
		s_bake.status = BakeStatus::Finished;
		return;
	}
	s_bake.status = BakeStatus::Tracing;
	if (!bakeTraceRays()) {
		s_bake.status = BakeStatus::Finished;
		return;
	}
	// Make sure there's a lightmap texture update before denoising.
	// Otherwise you're looking at an unfinished bake while the denoiser is running.
	if (s_bake.options.denoise) {
		while (s_bake.updateStatus != UpdateStatus::Idle) {}
		memcpy(s_bake.updateData.data(), s_bake.lightmapData.data(), s_bake.lightmapData.size() * sizeof(float));
		s_bake.updateStatus = UpdateStatus::Pending;
	}
	s_bake.status = BakeStatus::Denoising;
	bakeDenoise();
	if (bakeIsWorkerThreadCancelled()) {
		s_bake.status = BakeStatus::Finished;
		return;
	}
	// Dilate
	s_bake.dilatedLightmapData.resize(s_bake.lightmapWidth * s_bake.lightmapHeight * 4);
	lmImageDilate(s_bake.denoiseSucceeded ? s_bake.denoisedLightmapData.data() : s_bake.lightmapData.data(), s_bake.dilatedLightmapData.data(), (int)s_bake.lightmapWidth, (int)s_bake.lightmapHeight, 4);
	s_bake.status = BakeStatus::ThreadFinished;
}

static void bakeShutdownWorkerThread()
{
	s_bake.cancelWorkerMutex.lock();
	s_bake.cancelWorker = true;
	s_bake.cancelWorkerMutex.unlock();
	if (s_bake.workerThread) {
		s_bake.workerThread->join();
		delete s_bake.workerThread;
		s_bake.workerThread = nullptr;
	}
}

void bakeInit()
{
}

void bakeShutdown()
{
	bakeShutdownWorkerThread();
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
	if (s_bake.initialized)
		bgfx::destroy(s_bake.lightmap);
}

void bakeExecute()
{
	if (!(s_bake.status == BakeStatus::Idle || s_bake.status == BakeStatus::Finished || s_bake.status == BakeStatus::Error))
		return;
	// Re-create lightmap if atlas resolution has changed.
	const bool lightmapResolutionChanged = s_bake.lightmapWidth != atlasGetWidth() || s_bake.lightmapHeight != atlasGetHeight();
	if (!s_bake.initialized || lightmapResolutionChanged) {
		s_bake.lightmapWidth = atlasGetWidth();
		s_bake.lightmapHeight = atlasGetHeight();
		if (s_bake.initialized)
			bgfx::destroy(s_bake.lightmap);
		s_bake.lightmap = bgfx::createTexture2D((uint16_t)s_bake.lightmapWidth, (uint16_t)s_bake.lightmapHeight, false, 1, bgfx::TextureFormat::RGBA32F, BGFX_SAMPLER_UVW_CLAMP);
	}
	s_bake.initialized = true;
	g_options.shadeMode = ShadeMode::Lightmap;
	s_bake.status = BakeStatus::InitEmbree;
	s_bake.updateStatus = UpdateStatus::Idle;
	s_bake.cancelWorker = false;
	s_bake.workerThread = new std::thread(bakeWorkerThread);
}

void bakeFrame(uint32_t frameNo)
{
	if (s_bake.updateStatus == UpdateStatus::Pending) {
		bgfx::updateTexture2D(s_bake.lightmap, 0, 0, 0, 0, (uint16_t)s_bake.lightmapWidth, (uint16_t)s_bake.lightmapHeight, bgfx::makeRef(s_bake.updateData.data(), (uint32_t)s_bake.updateData.size() * sizeof(float)));
		s_bake.updateStatus = UpdateStatus::Updating;
		s_bake.updateFinishedFrameNo = frameNo + 2;
	} else if (s_bake.updateStatus == UpdateStatus::Updating) {
		if (frameNo == s_bake.updateFinishedFrameNo)
			s_bake.updateStatus = UpdateStatus::Idle;
	}
	if (s_bake.status == BakeStatus::ThreadFinished) {
		bakeShutdownWorkerThread();
		// Do a final update of the lightmap texture with the dilated result.
		bgfx::updateTexture2D(s_bake.lightmap, 0, 0, 0, 0, (uint16_t)s_bake.lightmapWidth, (uint16_t)s_bake.lightmapHeight, bgfx::makeRef(s_bake.dilatedLightmapData.data(), (uint32_t)s_bake.dilatedLightmapData.size() * sizeof(float)));
		s_bake.status = BakeStatus::Finished;
	} else if (s_bake.status == BakeStatus::Error && g_options.shadeMode == ShadeMode::Lightmap) {
		// Executing bake sets shade mode to lightmap. Set it back to charts if there was an error.
		g_options.shadeMode = ShadeMode::Charts;
	}
}

void bakeClear()
{
	s_bake.status = BakeStatus::Idle;
	g_options.shadeMode = atlasIsReady() ? ShadeMode::Charts : ShadeMode::Flat;
	if (s_bake.embreeDevice) {
		embree::ReleaseGeometry(s_bake.embreeGeometry);
		embree::ReleaseScene(s_bake.embreeScene);
		embree::ReleaseDevice(s_bake.embreeDevice);
		s_bake.embreeDevice = nullptr;
	}
	s_bake.sampleLocations.clear();
}

void bakeShowGuiOptions()
{
	ImGui::Text("Lightmap");
	if (atlasGetCount() > 1) {
		ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.0f, 0.0f, 1.0f));
		ImGui::Text("Baking doesn't support multiple atlases");
		ImGui::PopStyleColor();
		return;
	}
#if BX_ARCH_32BIT
	ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.0f, 0.0f, 1.0f));
	ImGui::Text("Baking not supported in 32-bit build");
	ImGui::PopStyleColor();
#else
	if (s_bake.status == BakeStatus::Idle || s_bake.status == BakeStatus::Finished || s_bake.status == BakeStatus::Error) {
		ImGui::Checkbox("Denoise", &s_bake.options.denoise);
		ImGui::Checkbox("Sky", &s_bake.options.sky);
		ImGui::SameLine();
		ImGui::ColorEdit3("Sky color", &s_bake.options.skyColor.x, ImGuiColorEditFlags_NoInputs);
		static const char * const sampleLabels[] = { "1", "4", "8", "16", "256" };
		static const int samples[] = { 1, 4, 8, 16, 256 };
		int sampleIndex = 0;
		for (int i = 0; i < (int)BX_COUNTOF(samples); i++) {
			if (samples[i] == s_bake.options.numSamples) {
				sampleIndex = i;
				break;
			}
		}
		ImGui::Combo("Samples", &sampleIndex, sampleLabels, (int)BX_COUNTOF(sampleLabels));
		s_bake.options.numSamples = samples[sampleIndex];
		//ImGui::SliderInt("Bounces", &s_bake.options.numBounces, 0, 4);
		const ImVec2 buttonSize(ImVec2(ImGui::GetContentRegionAvailWidth() * 0.3f, 0.0f));
		if (ImGui::Button("Bake", buttonSize))
			bakeExecute();
		if (s_bake.status == BakeStatus::Error) {
			ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.0f, 0.0f, 1.0f));
			ImGui::Text(s_bake.errorMessage);
			ImGui::PopStyleColor();
		}
	} else if (s_bake.status == BakeStatus::InitEmbree) {
		ImGui::Text("Initializing Embree...");
	} else if (s_bake.status == BakeStatus::Rasterizing) {
		ImGui::Text("Rasterizing...");
		ImGui::ProgressBar(s_bake.numTrianglesRasterized / float(modelGetData()->numIndices / 3));
		if (ImGui::Button("Cancel"))
			bakeShutdownWorkerThread();
	} else if (s_bake.status == BakeStatus::Tracing) {
		ImGui::Text("Tracing rays...");
		ImGui::ProgressBar(s_bake.numSampleLocationsProcessed / float(s_bake.sampleLocations.size()));
		if (ImGui::Button("Cancel"))
			bakeShutdownWorkerThread();
	} else if (s_bake.status == BakeStatus::Denoising) {
		ImGui::Text("Denoising...");
		ImGui::ProgressBar((float)s_bake.denoiseProgress);
		if (ImGui::Button("Cancel"))
			bakeShutdownWorkerThread();
	}
#endif
}

void bakeShowGuiWindow()
{
	if (s_bake.status == BakeStatus::Idle || s_bake.status == BakeStatus::InitEmbree || s_bake.status == BakeStatus::Error || !g_options.showLightmapWindow)
		return;
	ImGuiIO &io = ImGui::GetIO();
	const float size = 500;
	const float margin = 4.0f;
	ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x - size - margin, size + margin * 2.0f), ImGuiCond_FirstUseEver);
	ImGui::SetNextWindowSize(ImVec2(size, size), ImGuiCond_FirstUseEver);
	if (ImGui::Begin("Lightmap", &g_options.showLightmapWindow, ImGuiWindowFlags_HorizontalScrollbar)) {
		ImGui::Checkbox("Fit to window", &s_bake.options.fitToWindow);
		GuiTexture texture;
		texture.bgfx.handle = s_bake.lightmap;
		texture.bgfx.flags = GuiTextureFlags::PointSampler;
		if (s_bake.options.fitToWindow)
			ImGui::Image(texture.imgui, ImGui::GetContentRegionAvail());
		else
			ImGui::Image(texture.imgui, ImVec2((float)s_bake.lightmapWidth, (float)s_bake.lightmapHeight));
		ImGui::End();
	}
}

bgfx::TextureHandle bakeGetLightmap()
{
	return s_bake.lightmap;
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
