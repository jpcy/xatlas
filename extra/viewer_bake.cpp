/*
xatlas
https://github.com/jpcy/xatlas
Copyright (c) 2018 Jonathan Young

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#include <assert.h>
#include <time.h>
#include <mutex>
#include <thread>
#include <vector>
#include <bx/os.h>
#include <bx/rng.h>
#include <embree3/rtcore.h>
#include <imgui/imgui.h>
#include <OpenImageDenoise/oidn.h>
#include "shaders/shared.h"
#include "viewer.h"

#define DEBUG_RAY_BUNDLE 0

struct BakeStatus
{
	enum Enum
	{
		Idle,
		Rendering,
		ReadingLightmap,
		Denoising,
		WritingLightmap,
		Finished
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

struct BakeOptions
{
	bool fitToWindow = true;
	bool denoise = false;
	bool sky = true;
	bx::Vec3 skyColor = bx::Vec3(0.5f, 0.5f, 0.5f);
	int resolution = 512;
	int directionsPerFrame = 10;
	int numDirections = 300;
	int numBounces = 1;
};

struct Lightmap
{
	enum
	{
		CurrPass,
		PrevPass,
		Sum,
		Final,
		Num
	};
};

static constexpr int s_maxDirections = 5000;

struct
{
	const uint16_t rbDataTextureSize = 8192;
	float fsOrtho[16];
	bool enabled;
	bool initialized = false;
	BakeStatus status;
	BakeOptions options;
	int directionCount;
	int passCount;
	int resolution;
	uint32_t lightmapWidth, lightmapHeight;
	std::vector<float> lightmapData;
	std::vector<float> denoisedLightmapData;
	uint32_t lightmapDataReadyFrameNo;
	void *embreeLibrary = nullptr;
	RTCDevice embreeDevice = nullptr;
	RTCScene embreeScene = nullptr;
	RTCGeometry embreeGeometry = nullptr;
	std::thread *denoiseThread = nullptr;
	void *oidnLibrary = nullptr;
	bx::RngMwc rng;
	clock_t lastUpdateTime = 0;
	const double updateIntervalMs = 50;
	bx::Vec3 sampleDirections[s_maxDirections];
	// shaders
	bgfx::ShaderHandle fs_atomicCounterClear;
	bgfx::ShaderHandle fs_lightmapAverage;
	bgfx::ShaderHandle fs_lightmapOp;
	bgfx::ShaderHandle fs_rayBundleClear;
	bgfx::ShaderHandle fs_rayBundleIntegrate;
	bgfx::ShaderHandle fs_rayBundleLightmapClear;
	bgfx::ShaderHandle fs_rayBundleWrite;
	// programs
	bgfx::ProgramHandle atomicCounterClearProgram;
	bgfx::ProgramHandle lightmapAverageProgram;
	bgfx::ProgramHandle lightmapOpProgram;
	bgfx::ProgramHandle rayBundleClearProgram;
	bgfx::ProgramHandle rayBundleIntegrateProgram;
	bgfx::ProgramHandle rayBundleLightmapClearProgram;
	bgfx::ProgramHandle rayBundleWriteProgram;
	// uniforms
	bgfx::UniformHandle u_pass;
	bgfx::UniformHandle u_lightmapOp_bounceWeight;
	bgfx::UniformHandle u_lightmapSize_dataSize;
	bgfx::UniformHandle u_rayNormal;
	bgfx::UniformHandle u_skyColor_enabled;
	bgfx::UniformHandle s_atomicCounter;
	bgfx::UniformHandle s_rayBundleHeader;
	bgfx::UniformHandle s_rayBundleData;
	bgfx::UniformHandle s_rayBundleLightmap;
	bgfx::UniformHandle s_lightmap;
	bgfx::UniformHandle s_lightmapCurrPass;
	bgfx::UniformHandle s_lightmapPrevPass;
	bgfx::UniformHandle s_lightmapSum;
	// framebuffer targets
	bgfx::TextureHandle rayBundleTarget;
	bgfx::TextureHandle lightmapTarget;
	// atomic counter
	bgfx::FrameBufferHandle atomicCounterFb;
	bgfx::TextureHandle atomicCounterTexture;
	// ray bundle write
	bgfx::FrameBufferHandle rayBundleFb;
	bgfx::TextureHandle rayBundleHeader;
	bgfx::TextureHandle rayBundleData;
	// ray bundle integrate
	bgfx::FrameBufferHandle rayBundleIntegrateFb;
	bgfx::TextureHandle rayBundleLightmap; // Ray bundle is integrated into this. Cleared every ray bundle.
	bgfx::TextureHandle rayBundleIntegrateHack;
	// ray bundle lightmap clear
	bgfx::FrameBufferHandle rayBundleLightmapClearFb;
	// lightmap op
	bgfx::FrameBufferHandle lightmapOpFb;
	// lightmap average
	bgfx::FrameBufferHandle lightmapAverageFb;
	// lightmap textures
	bgfx::TextureHandle lightmaps[Lightmap::Num];
#if DEBUG_RAY_BUNDLE
	bgfx::TextureHandle rayBundleDebugWrite;
	bgfx::UniformHandle s_rayBundleDebugWrite;
	bgfx::TextureHandle rayBundleDebugIntegrate;
	bgfx::UniformHandle s_rayBundleDebugIntegrate;
#endif
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
};

namespace oidn
{
	typedef OIDNDevice (*NewDeviceFunc)(OIDNDeviceType type);
	typedef void (*CommitDeviceFunc)(OIDNDevice device);
	typedef void (*ReleaseDeviceFunc)(OIDNDevice device);
	typedef void (*SetDevice1bFunc)(OIDNDevice device, const char* name, bool value);
	typedef OIDNError (*GetDeviceErrorFunc)(OIDNDevice device, const char** outMessage);
	typedef OIDNFilter (*NewFilterFunc)(OIDNDevice device, const char* type);
	typedef void (*SetSharedFilterImageFunc)(OIDNFilter filter, const char* name, void* ptr, OIDNFormat format, size_t width, size_t height, size_t byteOffset, size_t bytePixelStride, size_t byteRowStride);
	typedef void (*SetFilter1bFunc)(OIDNFilter filter, const char* name, bool value);
	typedef void (*CommitFilterFunc)(OIDNFilter filter);
	typedef void (*ExecuteFilterFunc)(OIDNFilter filter);
	typedef void (*ReleaseFilterFunc)(OIDNFilter filter);
	NewDeviceFunc NewDevice;
	CommitDeviceFunc CommitDevice;
	ReleaseDeviceFunc ReleaseDevice;
	SetDevice1bFunc SetDevice1b;
	GetDeviceErrorFunc GetDeviceError;
	NewFilterFunc NewFilter;
	SetSharedFilterImageFunc SetSharedFilterImage;
	SetFilter1bFunc SetFilter1b;
	CommitFilterFunc CommitFilter;
	ExecuteFilterFunc ExecuteFilter;
	ReleaseFilterFunc ReleaseFilter;
};

static void bakeEmbreeError(void* /*userPtr*/, enum RTCError /*code*/, const char* str)
{
	fprintf(stderr, "Embree error: %s\n", str);
	exit(EXIT_FAILURE);
}

void bakeInit()
{
	if (!s_bake.embreeLibrary) {
		s_bake.embreeLibrary = bx::dlopen(EMBREE_LIB);
		if (!s_bake.embreeLibrary) {
			printf("Embree not installed.\n");
		} else {
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
		}
	}
	if (bgfx::getRendererType() != bgfx::RendererType::OpenGL) {
		s_bake.enabled = false;
		printf("Baking is disabled for all renderers except OpenGL. Use '--gl' command line argument.\n");
		return;
	}
	s_bake.enabled = (bgfx::getCaps()->supported & BGFX_CAPS_FRAMEBUFFER_RW) != 0;
	if (!s_bake.enabled) {
		printf("Read/Write frame buffer attachments are not supported. Baking is disabled.\n");
		return;
	}
	bx::mtxOrtho(s_bake.fsOrtho, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f, bgfx::getCaps()->homogeneousDepth);
}

void bakeShutdown()
{
	if (!s_bake.initialized)
		return;
	if (s_bake.embreeLibrary)
		bx::dlclose(s_bake.embreeLibrary);
	if (s_bake.oidnLibrary)
		bx::dlclose(s_bake.oidnLibrary);
	// shaders
	bgfx::destroy(s_bake.fs_atomicCounterClear);
	bgfx::destroy(s_bake.fs_lightmapOp);
	bgfx::destroy(s_bake.fs_lightmapAverage);
	bgfx::destroy(s_bake.fs_rayBundleClear);
	bgfx::destroy(s_bake.fs_rayBundleIntegrate);
	bgfx::destroy(s_bake.fs_rayBundleLightmapClear);
	bgfx::destroy(s_bake.fs_rayBundleWrite);
	// programs
	bgfx::destroy(s_bake.atomicCounterClearProgram);
	bgfx::destroy(s_bake.rayBundleClearProgram);
	bgfx::destroy(s_bake.rayBundleWriteProgram);
	bgfx::destroy(s_bake.rayBundleIntegrateProgram);
	bgfx::destroy(s_bake.rayBundleLightmapClearProgram);
	bgfx::destroy(s_bake.lightmapAverageProgram);
	bgfx::destroy(s_bake.lightmapOpProgram);
	// uniforms
	bgfx::destroy(s_bake.u_pass);
	bgfx::destroy(s_bake.u_lightmapOp_bounceWeight);
	bgfx::destroy(s_bake.u_lightmapSize_dataSize);
	bgfx::destroy(s_bake.u_rayNormal);
	bgfx::destroy(s_bake.u_skyColor_enabled);
	bgfx::destroy(s_bake.s_atomicCounter);
	bgfx::destroy(s_bake.s_rayBundleHeader);
	bgfx::destroy(s_bake.s_rayBundleData);
	bgfx::destroy(s_bake.s_rayBundleLightmap);
	bgfx::destroy(s_bake.s_lightmap);
	bgfx::destroy(s_bake.s_lightmapCurrPass);
	bgfx::destroy(s_bake.s_lightmapPrevPass);
	bgfx::destroy(s_bake.s_lightmapSum);
	// framebuffers
	bgfx::destroy(s_bake.rayBundleTarget);
	bgfx::destroy(s_bake.lightmapTarget);
	bgfx::destroy(s_bake.atomicCounterFb);
	bgfx::destroy(s_bake.rayBundleFb);
	bgfx::destroy(s_bake.rayBundleHeader);
	bgfx::destroy(s_bake.rayBundleData);
	bgfx::destroy(s_bake.rayBundleLightmap);
	bgfx::destroy(s_bake.rayBundleIntegrateFb);
	bgfx::destroy(s_bake.rayBundleIntegrateHack);
	bgfx::destroy(s_bake.rayBundleLightmapClearFb);
	bgfx::destroy(s_bake.lightmapOpFb);
	bgfx::destroy(s_bake.lightmapAverageFb);
	for (int i = 0; i < Lightmap::Num; i++)
		bgfx::destroy(s_bake.lightmaps[i]);
#if DEBUG_RAY_BUNDLE
	bgfx::destroy(s_bake.rayBundleDebugWrite);
	bgfx::destroy(s_bake.s_rayBundleDebugWrite);
	bgfx::destroy(s_bake.rayBundleDebugIntegrate);
	bgfx::destroy(s_bake.s_rayBundleDebugIntegrate);
#endif
}

static void setScreenSpaceQuadVertexBuffer()
{
	const uint32_t nVerts = 3;
	if (bgfx::getAvailTransientVertexBuffer(nVerts, PosVertex::decl) < nVerts)
		return;
	bgfx::TransientVertexBuffer vb;
	bgfx::allocTransientVertexBuffer(&vb, nVerts, PosVertex::decl);
	auto vertices = (PosVertex *)vb.data;
	vertices[0].pos[0] = -1.0f;
	vertices[0].pos[1] = 0.0f;
	vertices[1].pos[0] = 1.0f;
	vertices[1].pos[1] = 0.0f;
	vertices[2].pos[0] = 1.0f;
	vertices[2].pos[1] = 2.0f;
	bgfx::setVertexBuffer(0, &vb);
}

// From bx::generateSphereHammersley
static void generateHemisphereHammersley(void* _data, uint32_t _stride, uint32_t _num, float _scale = 1.0f)
{
	// Reference(s):
	// - Sampling with Hammersley and Halton Points
	//   https://web.archive.org/web/20190207230709/http://www.cse.cuhk.edu.hk/~ttwong/papers/udpoint/udpoints.html

	uint8_t* data = (uint8_t*)_data;

	for (uint32_t ii = 0; ii < _num; ii++)
	{
		float tt = 0.0f;
		float pp = 0.5;
		for (uint32_t jj = ii; jj; jj >>= 1)
		{
			tt += (jj & 1) ? pp : 0.0f;
			pp *= 0.5f;
		}

		tt = 2.0f * tt - 1.0f;

		const float phi = (ii + 0.5f) / _num;
		const float phirad = phi * bx::kPi;
		const float st = bx::sqrt(1.0f - tt * tt) * _scale;

		float* xyz = (float*)data;
		data += _stride;

		xyz[0] = st * bx::cos(phirad);
		xyz[1] = st * bx::sin(phirad);
		xyz[2] = tt * _scale;
	}
}

void bakeExecute()
{
	if (!(s_bake.status == BakeStatus::Idle || s_bake.status == BakeStatus::Finished))
		return;
	bakeClear();
	// Init embree.
	s_bake.embreeDevice = embree::NewDevice(nullptr);
	if (!s_bake.embreeDevice) {
		fprintf(stderr, "Error creating Embree device\n");
		return;
	}
	embree::SetDeviceErrorFunction(s_bake.embreeDevice, bakeEmbreeError, nullptr);
	s_bake.embreeGeometry = embree::NewGeometry(s_bake.embreeDevice, RTC_GEOMETRY_TYPE_TRIANGLE);
	const objzModel *model = modelGetData();
	embree::SetSharedGeometryBuffer(s_bake.embreeGeometry, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, model->vertices, offsetof(ModelVertex, pos), sizeof(ModelVertex), (size_t)model->numVertices);
	//embree::SetSharedGeometryBuffer(s_bake.embreeGeometry, RTC_BUFFER_TYPE_NORMAL, 0, RTC_FORMAT_FLOAT3, model->vertices, offsetof(ModelVertex, normal), sizeof(ModelVertex), (size_t)model->numVertices);
	embree::SetSharedGeometryBuffer(s_bake.embreeGeometry, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, model->indices, 0, sizeof(uint32_t), (size_t)model->numIndices / 3);
	embree::CommitGeometry(s_bake.embreeGeometry);
	s_bake.embreeScene = embree::NewScene(s_bake.embreeDevice);
	embree::AttachGeometry(s_bake.embreeScene, s_bake.embreeGeometry);
	embree::CommitScene(s_bake.embreeScene);
	generateHemisphereHammersley(s_bake.sampleDirections, sizeof(bx::Vec3), (uint32_t)s_bake.options.numDirections);
	if (!s_bake.initialized) {
		// shaders
		s_bake.u_pass = bgfx::createUniform("u_pass", bgfx::UniformType::Vec4);
		s_bake.u_lightmapOp_bounceWeight = bgfx::createUniform("u_lightmapOp_bounceWeight", bgfx::UniformType::Vec4);
		s_bake.u_lightmapSize_dataSize = bgfx::createUniform("u_lightmapSize_dataSize", bgfx::UniformType::Vec4);
		s_bake.u_rayNormal = bgfx::createUniform("u_rayNormal", bgfx::UniformType::Vec4);
		s_bake.u_skyColor_enabled = bgfx::createUniform("u_skyColor_enabled", bgfx::UniformType::Vec4);
		s_bake.s_atomicCounter = bgfx::createUniform("s_atomicCounter", bgfx::UniformType::Sampler);
		s_bake.s_rayBundleHeader = bgfx::createUniform("s_rayBundleHeader", bgfx::UniformType::Sampler);
		s_bake.s_rayBundleData = bgfx::createUniform("s_rayBundleData", bgfx::UniformType::Sampler);
		s_bake.s_rayBundleLightmap = bgfx::createUniform("s_rayBundleLightmap", bgfx::UniformType::Sampler);
		s_bake.s_lightmap = bgfx::createUniform("s_lightmap", bgfx::UniformType::Sampler);
		s_bake.s_lightmapCurrPass = bgfx::createUniform("s_lightmapCurrPass", bgfx::UniformType::Sampler);
		s_bake.s_lightmapPrevPass = bgfx::createUniform("s_lightmapPrevPass", bgfx::UniformType::Sampler);
		s_bake.s_lightmapSum = bgfx::createUniform("s_lightmapSum", bgfx::UniformType::Sampler);
#if DEBUG_RAY_BUNDLE
		s_bake.s_rayBundleDebugWrite = bgfx::createUniform("s_rayBundleDebugWrite", bgfx::UniformType::Sampler);
		s_bake.s_rayBundleDebugIntegrate = bgfx::createUniform("s_rayBundleDebugIntegrate", bgfx::UniformType::Sampler);
#endif
		s_bake.fs_atomicCounterClear = loadShader(ShaderId::fs_atomicCounterClear);
		s_bake.fs_lightmapAverage = loadShader(ShaderId::fs_lightmapAverage);
		s_bake.fs_lightmapOp = loadShader(ShaderId::fs_lightmapOp);
		s_bake.fs_rayBundleClear = loadShader(ShaderId::fs_rayBundleClear);
		s_bake.fs_rayBundleIntegrate = loadShader(ShaderId::fs_rayBundleIntegrate);
		s_bake.fs_rayBundleLightmapClear = loadShader(ShaderId::fs_rayBundleLightmapClear);
		s_bake.fs_rayBundleWrite = loadShader(ShaderId::fs_rayBundleWrite);
		s_bake.atomicCounterClearProgram = bgfx::createProgram(get_vs_position(), s_bake.fs_atomicCounterClear);
		s_bake.lightmapOpProgram = bgfx::createProgram(get_vs_position(), s_bake.fs_lightmapOp);
		s_bake.lightmapAverageProgram = bgfx::createProgram(get_vs_position(), s_bake.fs_lightmapAverage);
		s_bake.rayBundleClearProgram = bgfx::createProgram(get_vs_position(), s_bake.fs_rayBundleClear);
		s_bake.rayBundleIntegrateProgram = bgfx::createProgram(get_vs_position(), s_bake.fs_rayBundleIntegrate);
		s_bake.rayBundleLightmapClearProgram = bgfx::createProgram(get_vs_position(), s_bake.fs_rayBundleLightmapClear);
		s_bake.rayBundleWriteProgram = bgfx::createProgram(modelGet_vs_model(), s_bake.fs_rayBundleWrite);
		// Atomic counter
		bgfx::TextureHandle target = bgfx::createTexture2D(1, 1, false, 1, bgfx::TextureFormat::BGRA8, BGFX_TEXTURE_RT);
		s_bake.atomicCounterTexture = bgfx::createTexture2D(1, 1, false, 1, bgfx::TextureFormat::R32U, BGFX_TEXTURE_COMPUTE_WRITE | BGFX_SAMPLER_POINT | BGFX_SAMPLER_UVW_CLAMP);
		bgfx::Attachment attachments[2];
		attachments[0].init(target);
		attachments[1].init(s_bake.atomicCounterTexture, bgfx::Access::ReadWrite);
		s_bake.atomicCounterFb = bgfx::createFrameBuffer(BX_COUNTOF(attachments), attachments, true);
	}
	// Re-create ray bundle data if ray bundle resolution has changed.
	const bool rayBundleResolutionChanged = s_bake.resolution != s_bake.options.resolution;
	if (!s_bake.initialized || rayBundleResolutionChanged) {
		s_bake.resolution = s_bake.options.resolution;
		if (s_bake.initialized) {
			bgfx::destroy(s_bake.rayBundleTarget);
			bgfx::destroy(s_bake.rayBundleFb);
			bgfx::destroy(s_bake.rayBundleHeader);
			bgfx::destroy(s_bake.rayBundleData);
			bgfx::destroy(s_bake.rayBundleIntegrateHack);
#if DEBUG_RAY_BUNDLE
			bgfx::destroy(s_bake.rayBundleDebugWrite);
#endif
		}
		s_bake.rayBundleTarget = bgfx::createTexture2D((uint16_t)s_bake.resolution, (uint16_t)s_bake.resolution, false, 1, bgfx::TextureFormat::BGRA8, BGFX_TEXTURE_RT);
		s_bake.rayBundleHeader = bgfx::createTexture2D((uint16_t)s_bake.resolution, (uint16_t)s_bake.resolution, false, 1, bgfx::TextureFormat::R32U, BGFX_TEXTURE_COMPUTE_WRITE | BGFX_SAMPLER_POINT | BGFX_SAMPLER_UVW_CLAMP);
		s_bake.rayBundleData = bgfx::createTexture2D(s_bake.rbDataTextureSize, s_bake.rbDataTextureSize, false, 1, bgfx::TextureFormat::RGBA32U, BGFX_TEXTURE_COMPUTE_WRITE | BGFX_SAMPLER_POINT | BGFX_SAMPLER_UVW_CLAMP);
#if DEBUG_RAY_BUNDLE
		s_bake.rayBundleDebugWrite = bgfx::createTexture2D((uint16_t)s_bake.resolution, (uint16_t)s_bake.resolution, false, 1, bgfx::TextureFormat::RGBA8, BGFX_TEXTURE_COMPUTE_WRITE | BGFX_SAMPLER_POINT | BGFX_SAMPLER_UVW_CLAMP);
		bgfx::Attachment attachments[5];
		attachments[4].init(s_bake.rayBundleDebugWrite, bgfx::Access::ReadWrite);
#else
		bgfx::Attachment attachments[4];
#endif
		attachments[0].init(s_bake.rayBundleTarget);
		attachments[1].init(s_bake.atomicCounterTexture, bgfx::Access::ReadWrite);
		attachments[2].init(s_bake.rayBundleHeader, bgfx::Access::ReadWrite);
		attachments[3].init(s_bake.rayBundleData, bgfx::Access::ReadWrite);
		s_bake.rayBundleFb = bgfx::createFrameBuffer(BX_COUNTOF(attachments), attachments);
		s_bake.rayBundleIntegrateHack = bgfx::createTexture2D((uint16_t)s_bake.resolution, (uint16_t)s_bake.resolution, false, 1, bgfx::TextureFormat::BGRA8, BGFX_TEXTURE_COMPUTE_WRITE | BGFX_SAMPLER_POINT | BGFX_SAMPLER_UVW_CLAMP);
	}
	// Re-create lightmap if atlas resolution has changed.
	const bool lightmapResolutionChanged = s_bake.lightmapWidth != atlasGetWidth() || s_bake.lightmapHeight != atlasGetHeight();
	if (!s_bake.initialized || lightmapResolutionChanged) {
		if (s_bake.initialized) {
			bgfx::destroy(s_bake.lightmapTarget);
			bgfx::destroy(s_bake.lightmapOpFb);
			bgfx::destroy(s_bake.rayBundleLightmap);
			bgfx::destroy(s_bake.rayBundleLightmapClearFb);
			bgfx::destroy(s_bake.lightmapAverageFb);
			for (int i = 0; i < Lightmap::Num; i++)
				bgfx::destroy(s_bake.lightmaps[i]);
		}
		s_bake.lightmapWidth = atlasGetWidth();
		s_bake.lightmapHeight = atlasGetHeight();
		s_bake.lightmapTarget = bgfx::createTexture2D((uint16_t)s_bake.lightmapWidth, (uint16_t)s_bake.lightmapHeight, false, 1, bgfx::TextureFormat::BGRA8, BGFX_TEXTURE_RT);
		s_bake.rayBundleLightmap = bgfx::createTexture2D((uint16_t)s_bake.lightmapWidth * 4, (uint16_t)s_bake.lightmapHeight, false, 1, bgfx::TextureFormat::R32U, BGFX_TEXTURE_COMPUTE_WRITE | BGFX_SAMPLER_POINT | BGFX_SAMPLER_UVW_CLAMP);
		for (int i = 0; i < Lightmap::Num; i++)
			s_bake.lightmaps[i] = bgfx::createTexture2D((uint16_t)s_bake.lightmapWidth, (uint16_t)s_bake.lightmapHeight, false, 1, bgfx::TextureFormat::RGBA32F, BGFX_TEXTURE_COMPUTE_WRITE | BGFX_SAMPLER_UVW_CLAMP);
		{
			bgfx::Attachment attachments[2];
			attachments[0].init(s_bake.lightmapTarget);
			attachments[1].init(s_bake.rayBundleLightmap, bgfx::Access::ReadWrite);
			s_bake.rayBundleLightmapClearFb = bgfx::createFrameBuffer(BX_COUNTOF(attachments), attachments);
		}
		{
			bgfx::Attachment attachments[1 + Lightmap::Num];
			attachments[0].init(s_bake.lightmapTarget);
			for (int i = 0; i < Lightmap::Num; i++)
				attachments[1 + i].init(s_bake.lightmaps[i], bgfx::Access::ReadWrite);
			s_bake.lightmapOpFb = bgfx::createFrameBuffer(BX_COUNTOF(attachments), attachments);
		}
		{
			bgfx::Attachment attachments[3];
			attachments[0].init(s_bake.lightmapTarget);
			attachments[1].init(s_bake.rayBundleLightmap, bgfx::Access::Read);
			attachments[2].init(s_bake.lightmaps[Lightmap::CurrPass], bgfx::Access::ReadWrite);
			s_bake.lightmapAverageFb = bgfx::createFrameBuffer(BX_COUNTOF(attachments), attachments);
		}
	}
	// Re-create ray bundle integrate fb if ray bundle or atlas resolution has changed.
	if (!s_bake.initialized || rayBundleResolutionChanged || lightmapResolutionChanged) {
		if (s_bake.initialized) {
			bgfx::destroy(s_bake.rayBundleIntegrateFb);
#if DEBUG_RAY_BUNDLE
			bgfx::destroy(s_bake.rayBundleDebugIntegrate);
#endif
		}
#if DEBUG_RAY_BUNDLE
		s_bake.rayBundleDebugIntegrate = bgfx::createTexture2D((uint16_t)s_bake.lightmapWidth, (uint16_t)s_bake.lightmapHeight, false, 1, bgfx::TextureFormat::RGBA8, BGFX_TEXTURE_COMPUTE_WRITE | BGFX_SAMPLER_POINT | BGFX_SAMPLER_UVW_CLAMP);
		//s_bake.rayBundleDebugIntegrate = bgfx::createTexture2D((uint16_t)s_bake.resolution, (uint16_t)s_bake.resolution, false, 1, bgfx::TextureFormat::RGBA8, BGFX_TEXTURE_COMPUTE_WRITE | BGFX_SAMPLER_POINT | BGFX_SAMPLER_UVW_CLAMP);
		bgfx::Attachment attachments[5];
		attachments[4].init(s_bake.rayBundleDebugIntegrate, bgfx::Access::ReadWrite);
#else
		bgfx::Attachment attachments[5];
		attachments[4].init(s_bake.rayBundleIntegrateHack, bgfx::Access::ReadWrite);
#endif
		attachments[0].init(s_bake.rayBundleTarget);
		attachments[1].init(s_bake.rayBundleHeader, bgfx::Access::Read);
		attachments[2].init(s_bake.rayBundleData, bgfx::Access::Read);
		attachments[3].init(s_bake.rayBundleLightmap, bgfx::Access::ReadWrite);
		s_bake.rayBundleIntegrateFb = bgfx::createFrameBuffer(BX_COUNTOF(attachments), attachments);
	}
	s_bake.initialized = true;
	s_bake.status = BakeStatus::Rendering;
	s_bake.directionCount = 0;
	s_bake.passCount = 0;
	s_bake.rng.reset();
	g_options.shadeMode = ShadeMode::Lightmap;
}

static void bakeDenoise()
{
	if (!s_bake.oidnLibrary) {
		s_bake.oidnLibrary = bx::dlopen("OpenImageDenoise.dll");
		if (!s_bake.oidnLibrary) {
			s_bake.status = BakeStatus::Finished;
			return;
		}
		oidn::NewDevice = (oidn::NewDeviceFunc)bx::dlsym(s_bake.oidnLibrary, "oidnNewDevice");
		oidn::CommitDevice = (oidn::CommitDeviceFunc)bx::dlsym(s_bake.oidnLibrary, "oidnCommitDevice");
		oidn::ReleaseDevice = (oidn::ReleaseDeviceFunc)bx::dlsym(s_bake.oidnLibrary, "oidnReleaseDevice");
		oidn::SetDevice1b = (oidn::SetDevice1bFunc)bx::dlsym(s_bake.oidnLibrary, "oidnSetDevice1b");
		oidn::GetDeviceError = (oidn::GetDeviceErrorFunc)bx::dlsym(s_bake.oidnLibrary, "oidnGetDeviceError");
		oidn::NewFilter = (oidn::NewFilterFunc)bx::dlsym(s_bake.oidnLibrary, "oidnNewFilter");
		oidn::SetSharedFilterImage = (oidn::SetSharedFilterImageFunc)bx::dlsym(s_bake.oidnLibrary, "oidnSetSharedFilterImage");
		oidn::SetFilter1b = (oidn::SetFilter1bFunc)bx::dlsym(s_bake.oidnLibrary, "oidnSetFilter1b");
		oidn::CommitFilter = (oidn::CommitFilterFunc)bx::dlsym(s_bake.oidnLibrary, "oidnCommitFilter");
		oidn::ExecuteFilter = (oidn::ExecuteFilterFunc)bx::dlsym(s_bake.oidnLibrary, "oidnExecuteFilter");
		oidn::ReleaseFilter = (oidn::ReleaseFilterFunc)bx::dlsym(s_bake.oidnLibrary, "oidnReleaseFilter");
	}
	// OIDN_FORMAT_FLOAT4 not supported.
	std::vector<float> input, output;
	input.resize(s_bake.lightmapWidth * s_bake.lightmapHeight * 3);
	for (uint32_t i = 0; i < s_bake.lightmapWidth * s_bake.lightmapHeight; i++) {
		const float *rgbaIn = &s_bake.lightmapData[i * 4];
		float *rgbOut = &input[i * 3];
		rgbOut[0] = rgbaIn[0];
		rgbOut[1] = rgbaIn[1];
		rgbOut[2] = rgbaIn[2];
	}
	output.resize(s_bake.lightmapWidth * s_bake.lightmapHeight * 3);
	OIDNDevice device = oidn::NewDevice(OIDN_DEVICE_TYPE_DEFAULT);
	if (!device) {
		fprintf(stderr, "Error creating OIDN device\n");
		exit(EXIT_FAILURE);
	}
	oidn::SetDevice1b(device, "setAffinity", false);
	oidn::CommitDevice(device);
	s_bake.status = BakeStatus::WritingLightmap;
	OIDNFilter filter = oidn::NewFilter(device, "RT");
	oidn::SetSharedFilterImage(filter, "color", input.data(), OIDN_FORMAT_FLOAT3, s_bake.lightmapWidth, s_bake.lightmapHeight, 0, 0, 0);
	oidn::SetSharedFilterImage(filter, "output", output.data(), OIDN_FORMAT_FLOAT3, s_bake.lightmapWidth, s_bake.lightmapHeight, 0, 0, 0);
	oidn::CommitFilter(filter);
	oidn::ExecuteFilter(filter);
	const char *errorMessage;
	if (oidn::GetDeviceError(device, &errorMessage) != OIDN_ERROR_NONE)
		fprintf(stderr, "Denoiser error: %s\n", errorMessage);
	oidn::ReleaseFilter(filter);
	oidn::ReleaseDevice(device);
	s_bake.denoisedLightmapData.resize(s_bake.lightmapWidth * s_bake.lightmapHeight * 4);
	for (uint32_t i = 0; i < s_bake.lightmapWidth * s_bake.lightmapHeight; i++) {
		const float *rgbIn = &output[i * 3];
		float *rgbaOut = &s_bake.denoisedLightmapData[i * 4];
		rgbaOut[0] = rgbIn[0];
		rgbaOut[1] = rgbIn[1];
		rgbaOut[2] = rgbIn[2];
		rgbaOut[3] = 1.0f;
	}
}

static float bakeCalculateBounceWeight(int bounce, int numBounces)
{
	return bx::pow(2.0f, (numBounces - (bounce + 1.0f))) / (bx::pow(2.0f, (float)numBounces) - 1.0f);
}

static void bakeSubmitLightmapOp(bgfx::ViewId viewId, int op)
{
	bgfx::setViewFrameBuffer(viewId, s_bake.lightmapOpFb);
	bgfx::setViewRect(viewId, 0, 0, (uint16_t)s_bake.lightmapWidth, (uint16_t)s_bake.lightmapHeight);
	bgfx::setViewTransform(viewId, nullptr, s_bake.fsOrtho);
	bgfx::setTexture(1, s_bake.s_lightmapCurrPass, s_bake.lightmaps[Lightmap::CurrPass]);
	bgfx::setTexture(2, s_bake.s_lightmapPrevPass, s_bake.lightmaps[Lightmap::PrevPass]);
	bgfx::setTexture(3, s_bake.s_lightmapSum, s_bake.lightmaps[Lightmap::Sum]);
	bgfx::setTexture(4, s_bake.s_lightmap, s_bake.lightmaps[Lightmap::Final]);
	float weight = 1.0f;
	if (s_bake.passCount > 0)
		weight = bakeCalculateBounceWeight(s_bake.passCount - 1, s_bake.options.numBounces);
	const float opData[] = { (float)op, weight, 0.0f, 0.0f };
	bgfx::setUniform(s_bake.u_lightmapOp_bounceWeight, opData);
	setScreenSpaceQuadVertexBuffer();
	bgfx::setState(0);
	bgfx::submit(viewId, s_bake.lightmapOpProgram);
}

void bakeFrame(uint32_t bgfxFrame)
{
	bgfx::ViewId viewId = kFirstFreeView;
	if (s_bake.status == BakeStatus::Rendering) {
		bool finishedPass = false;
		for (uint32_t i = 0; i < (uint32_t)s_bake.options.directionsPerFrame; i++) {
			if (s_bake.directionCount == 0) {
				// Clear ray bundle lightmap.
				bgfx::setViewFrameBuffer(viewId, s_bake.rayBundleLightmapClearFb);
				bgfx::setViewRect(viewId, 0, 0, (uint16_t)s_bake.lightmapWidth, (uint16_t)s_bake.lightmapHeight);
				bgfx::setViewTransform(viewId, nullptr, s_bake.fsOrtho);
				bgfx::setTexture(1, s_bake.s_rayBundleLightmap, s_bake.rayBundleLightmap);
				setScreenSpaceQuadVertexBuffer();
				bgfx::setState(0);
				bgfx::submit(viewId, s_bake.rayBundleLightmapClearProgram);
				viewId++;
				if (s_bake.passCount == 0) {
					// Lightmap clear sum.
					bakeSubmitLightmapOp(viewId, LIGHTMAP_OP_CLEAR_SUM);
					viewId++;
				}
			}
			// Atomic counter clear.
			bgfx::setViewFrameBuffer(viewId, s_bake.atomicCounterFb);
			bgfx::setViewRect(viewId, 0, 0, 1, 1);
			bgfx::setViewTransform(viewId, nullptr, s_bake.fsOrtho);
			bgfx::setTexture(1, s_bake.s_atomicCounter, s_bake.atomicCounterTexture);
			setScreenSpaceQuadVertexBuffer();
			bgfx::setState(0);
			bgfx::submit(viewId, s_bake.atomicCounterClearProgram);
			viewId++;
			// Ray bundle clear.
			bgfx::setViewFrameBuffer(viewId, s_bake.rayBundleFb);
			bgfx::setViewRect(viewId, 0, 0, (uint16_t)s_bake.resolution, (uint16_t)s_bake.resolution);
			bgfx::setViewTransform(viewId, nullptr, s_bake.fsOrtho);
			bgfx::setTexture(2, s_bake.s_rayBundleHeader, s_bake.rayBundleHeader);
			setScreenSpaceQuadVertexBuffer();
			bgfx::setState(0);
			bgfx::submit(viewId, s_bake.rayBundleClearProgram);
			viewId++;
			// Ray bundle write.
			bx::Vec3 forward = s_bake.sampleDirections[s_bake.directionCount];
			bx::Vec3 right = bx::randUnitSphere(&s_bake.rng);
			bx::Vec3 up = bx::cross(forward, right);
#if DEBUG_RAY_BUNDLE
			forward = bx::Vec3(0, 0, -1);
			up = bx::Vec3(0, 1, 0);
#endif
			float view[16];
			bx::mtxLookAt(view, bx::Vec3(0.0f), forward, up, bx::Handness::Right);
			AABB aabb;
#if 1
			bx::Vec3 corners[8];
			modelGetAABB().getCorners(corners);
			for (uint32_t j = 0; j < 8; j++) {
				const float in[4] = { corners[j].x, corners[j].y, corners[j].z, 1.0f };
				float out[4];
				bx::vec4MulMtx(out, in, view);
				aabb.addPoint(bx::Vec3(out[0], out[1], out[2]));
			}
#else
			for (uint32_t j = 0; j < s_model.data->numVertices; j++) {
				const auto &v = ((ModelVertex *)s_model.data->vertices)[j];
				const float in[4] = { v.pos.x, v.pos.y, v.pos.z, 1.0f };
				float out[4];
				bx::vec4MulMtx(out, in, view);
				aabb.addPoint(bx::Vec3(out[0], out[1], out[2]));
			}
#endif
			float projection[16];
			bx::mtxOrtho(projection, aabb.min.x, aabb.max.x, aabb.min.y, aabb.max.y, -aabb.max.z, -aabb.min.z, 0.0f, bgfx::getCaps()->homogeneousDepth, bx::Handness::Right);
			bgfx::setViewFrameBuffer(viewId, s_bake.rayBundleFb);
			bgfx::setViewRect(viewId, 0, 0, (uint16_t)s_bake.resolution, (uint16_t)s_bake.resolution);
			bgfx::setViewTransform(viewId, view, projection);
			{
				const float passData[] = { (float)s_bake.passCount, 0.0f, 0.0f, 0.0f };
				const float sizes[] = { (float)s_bake.lightmapWidth, (float)s_bake.lightmapHeight, (float)s_bake.rbDataTextureSize, 0.0f };
				const objzModel *model = modelGetData();
				for (uint32_t j = 0; j < model->numMeshes; j++) {
					const objzMesh &mesh = model->meshes[j];
					bgfx::setIndexBuffer(atlasGetIb(), mesh.firstIndex, mesh.numIndices);
					bgfx::setVertexBuffer(0, atlasGetVb());
					bgfx::setState(0);
					bgfx::setTexture(0, s_bake.s_lightmapPrevPass, s_bake.lightmaps[Lightmap::PrevPass]);
					bgfx::setTexture(1, s_bake.s_atomicCounter, s_bake.atomicCounterTexture);
					bgfx::setTexture(2, s_bake.s_rayBundleHeader, s_bake.rayBundleHeader);
					bgfx::setTexture(3, s_bake.s_rayBundleData, s_bake.rayBundleData);
#if DEBUG_RAY_BUNDLE
					bgfx::setTexture(4, s_bake.s_rayBundleDebugWrite, s_bake.rayBundleDebugWrite);
#endif
					bgfx::setUniform(s_bake.u_pass, passData);
					bgfx::setUniform(s_bake.u_lightmapSize_dataSize, sizes);
					modelSetMaterialTexturesAndUniforms(mesh.materialIndex, 1);
					bgfx::submit(viewId, s_bake.rayBundleWriteProgram);
				}
			}
			viewId++;
			// Ray bundle integrate.
			bgfx::setViewFrameBuffer(viewId, s_bake.rayBundleIntegrateFb);
			bgfx::setViewRect(viewId, 0, 0, (uint16_t)s_bake.resolution, (uint16_t)s_bake.resolution);
			bgfx::setViewTransform(viewId, nullptr, s_bake.fsOrtho);
			bgfx::setTexture(1, s_bake.s_rayBundleHeader, s_bake.rayBundleHeader);
			bgfx::setTexture(2, s_bake.s_rayBundleData, s_bake.rayBundleData);
			bgfx::setTexture(3, s_bake.s_rayBundleLightmap, s_bake.rayBundleLightmap);
#if DEBUG_RAY_BUNDLE
			bgfx::setTexture(4, s_bake.s_rayBundleDebugIntegrate, s_bake.rayBundleDebugIntegrate);
#endif
			const float sizes[] = { (float)s_bake.lightmapWidth, (float)s_bake.lightmapHeight, (float)s_bake.rbDataTextureSize, 0.0f };
			bgfx::setUniform(s_bake.u_lightmapSize_dataSize, sizes);
			const float rayNormal[] = { -view[2], -view[6], -view[10], 0 };
			bgfx::setUniform(s_bake.u_rayNormal, rayNormal);
			const float sky[] = { s_bake.options.skyColor.x, s_bake.options.skyColor.y, s_bake.options.skyColor.z, s_bake.passCount == 0 && s_bake.options.sky ? 1.0f : 0.0f };
			bgfx::setUniform(s_bake.u_skyColor_enabled, sky);
			setScreenSpaceQuadVertexBuffer();
			bgfx::setState(0);
			bgfx::submit(viewId, s_bake.rayBundleIntegrateProgram);
			viewId++;
			// Finished with this direction.
			s_bake.directionCount++;
			finishedPass = s_bake.directionCount >= s_bake.options.numDirections;
			if (finishedPass)
				break;
		}
		const double elapsedMs = (clock() - s_bake.lastUpdateTime) * 1000.0 / CLOCKS_PER_SEC;
		if (finishedPass || elapsedMs >= s_bake.updateIntervalMs) {
			s_bake.lastUpdateTime = clock();
			// Lightmap clear average.
			bakeSubmitLightmapOp(viewId, LIGHTMAP_OP_CLEAR_CURR);
			viewId++;
			// Lightmap average.
			bgfx::setViewFrameBuffer(viewId, s_bake.lightmapAverageFb);
			bgfx::setViewRect(viewId, 0, 0, (uint16_t)s_bake.lightmapWidth, (uint16_t)s_bake.lightmapHeight);
			bgfx::setViewTransform(viewId, nullptr, s_bake.fsOrtho);
			bgfx::setTexture(1, s_bake.s_rayBundleLightmap, s_bake.rayBundleLightmap);
			bgfx::setTexture(2, s_bake.s_lightmap, s_bake.lightmaps[Lightmap::CurrPass]);
			setScreenSpaceQuadVertexBuffer();
			bgfx::setState(0);
			bgfx::submit(viewId, s_bake.lightmapAverageProgram);
			viewId++;
			// Lightmap write final.
			bakeSubmitLightmapOp(viewId, LIGHTMAP_OP_WRITE_FINAL);
			viewId++;
		}
		if (finishedPass) {
			// Lightmap finish pass.
			bakeSubmitLightmapOp(viewId, LIGHTMAP_OP_FINISH_PASS);
			viewId++;
			s_bake.passCount++;
			s_bake.directionCount = 0;
			// Finished all passes?
			if (s_bake.passCount == s_bake.options.numBounces + 1) {
				if (s_bake.options.denoise) {
					s_bake.status = BakeStatus::ReadingLightmap;
					s_bake.lightmapData.resize(s_bake.lightmapWidth * s_bake.lightmapHeight * 4 * sizeof(float));
					s_bake.lightmapDataReadyFrameNo = bgfx::readTexture(s_bake.lightmaps[Lightmap::Final], s_bake.lightmapData.data());
				}
				else {
					s_bake.status = BakeStatus::Finished;
				}
			}
		}
	} else if (s_bake.status == BakeStatus::ReadingLightmap) {
		if (bgfxFrame >= s_bake.lightmapDataReadyFrameNo) {
			s_bake.status = BakeStatus::Denoising;
			s_bake.denoiseThread = new std::thread(bakeDenoise);
		}
	} else if (s_bake.status == BakeStatus::WritingLightmap) {
		s_bake.denoiseThread->join();
		delete s_bake.denoiseThread;
		s_bake.denoiseThread = nullptr;
		bgfx::updateTexture2D(s_bake.lightmaps[Lightmap::Final], 0, 0, 0, 0, (uint16_t)s_bake.lightmapWidth, (uint16_t)s_bake.lightmapHeight, bgfx::makeRef(s_bake.denoisedLightmapData.data(), (uint32_t)s_bake.denoisedLightmapData.size() * sizeof(float)));
		s_bake.status = BakeStatus::Finished;
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
}

void bakeShowGuiOptions()
{
	if (!s_bake.enabled)
		return;
	ImGui::Text("Lightmap");
	if (atlasGetCount() > 1) {
		ImGui::Text("Baking doesn't support multiple atlases");
		return;
	}
	if (s_bake.status == BakeStatus::Idle || s_bake.status == BakeStatus::Finished) {
		ImGui::Checkbox("Denoise", &s_bake.options.denoise);
		ImGui::Checkbox("Sky", &s_bake.options.sky);
		ImGui::SameLine();
		ImGui::ColorEdit3("Sky color", &s_bake.options.skyColor.x, ImGuiColorEditFlags_NoInputs);
		static const char * const resolutionLabels[] = { "256x256", "512x512", "1024x1024" };
		static const int resolutions[] = { 256, 512, 1024 };
		int resolutionIndex = 0;
		for (int i = 0; i < (int)BX_COUNTOF(resolutions); i++) {
			if (resolutions[i] == s_bake.options.resolution) {
				resolutionIndex = i;
				break;
			}
		}
		ImGui::ListBox("Ray bundle resolution", &resolutionIndex, resolutionLabels, (int)BX_COUNTOF(resolutionLabels));
		s_bake.options.resolution = resolutions[resolutionIndex];
		ImGui::SliderInt("Ray bundle directions", &s_bake.options.numDirections, 1, s_maxDirections);
		ImGui::SliderInt("Directions per frame", &s_bake.options.directionsPerFrame, 1, 100);
		ImGui::SliderInt("Bounces", &s_bake.options.numBounces, 0, 4);
		const ImVec2 buttonSize(ImVec2(ImGui::GetContentRegionAvailWidth() * 0.3f, 0.0f));
		if (ImGui::Button("Bake", buttonSize))
			bakeExecute();
	}
	else {
		if (s_bake.directionCount < s_bake.options.numDirections && s_bake.passCount < s_bake.options.numBounces + 1) {
			ImGui::SliderInt("Directions per frame", &s_bake.options.directionsPerFrame, 1, 100);
			if (s_bake.passCount > 0)
				ImGui::Text("Baking bounce %d of %d...", s_bake.passCount, s_bake.options.numBounces);
			else
				ImGui::Text("Baking...");
			ImGui::ProgressBar(s_bake.directionCount / (float)s_bake.options.numDirections);
			if (ImGui::Button("Cancel"))
				s_bake.status = BakeStatus::Finished;
		}
		else {
			ImGui::AlignTextToFramePadding();
			ImGui::Text("Denoising...");
		}
	}
}

void bakeShowGuiWindow()
{
	if (s_bake.status == BakeStatus::Idle || !g_options.showLightmapWindow)
		return;
	ImGuiIO &io = ImGui::GetIO();
	const float size = 500;
	const float margin = 4.0f;
	ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x - size - margin, size + margin * 2.0f), ImGuiCond_FirstUseEver);
	ImGui::SetNextWindowSize(ImVec2(size, size), ImGuiCond_FirstUseEver);
	if (ImGui::Begin("Lightmap", &g_options.showLightmapWindow, ImGuiWindowFlags_HorizontalScrollbar)) {
#if DEBUG_RAY_BUNDLE
		ImVec2 imageSize(ImGui::GetContentRegionAvail().x * 0.45f, ImGui::GetContentRegionAvail().y * 0.45f);
		GuiTexture texture;
		texture.bgfx.flags = GuiTextureFlags::pointSampling;
		texture.bgfx.handle = s_bake.rayBundleDebugWrite;
		ImGui::Image(texture.imgui, imageSize);
		texture.bgfx.handle = s_bake.rayBundleDebugIntegrate;
		ImGui::SameLine();
		ImGui::Image(texture.imgui, imageSize);
		texture.bgfx.handle = s_bake.lightmaps[Lightmap::Final];
		ImGui::Image(texture.imgui, imageSize);
#else
		ImGui::Checkbox("Fit to window", &s_bake.options.fitToWindow);
		GuiTexture texture;
		texture.bgfx.handle = s_bake.lightmaps[Lightmap::Final];
		texture.bgfx.flags = GuiTextureFlags::PointSampler;
		if (s_bake.options.fitToWindow)
			ImGui::Image(texture.imgui, ImGui::GetContentRegionAvail());
		else
			ImGui::Image(texture.imgui, ImVec2((float)s_bake.lightmapWidth, (float)s_bake.lightmapHeight));
#endif
		ImGui::End();
	}
}

bgfx::TextureHandle bakeGetLightmap()
{
	return s_bake.lightmaps[Lightmap::Final];
}

uint32_t bakeGetLightmapSamplerFlags()
{
	uint32_t flags = BGFX_SAMPLER_UVW_CLAMP;
	if (g_options.lightmapPointSampling)
		flags |= BGFX_SAMPLER_POINT;
	return flags;
}

bool bakeIsIdle()
{
	return s_bake.status == BakeStatus::Idle;
}
