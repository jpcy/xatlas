/*
xatlas
https://github.com/jpcy/xatlas
Copyright (c) 2018 Jonathan Young

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#include <assert.h>
#include <mutex>
#include <thread>
#include <vector>
#include <bx/os.h>
#include <bx/rng.h>
#include <imgui/imgui.h>
#include <OpenImageDenoise/oidn.h>
#include "viewer.h"

struct BakeStatus
{
	enum Enum
	{
		Idle,
		Executing,
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
	bool showLightmap = true;
	bool denoise = false;
	bool sky = true;
	bx::Vec3 skyColor = bx::Vec3(1.0f, 1.0f, 1.0f);
	int resolution = 512;
	int directionsPerFrame = 10;
	int numDirections = 300;
};

struct
{
	const uint16_t rbDataTextureSize = 8192;
	float fsOrtho[16];
	bool enabled;
	bool initialized = false;
	BakeStatus status;
	BakeOptions options;
	int directionCount;
	int resolution;
	uint32_t lightmapWidth, lightmapHeight;
	std::vector<float> lightmapData;
	std::vector<float> denoisedLightmapData;
	uint32_t lightmapDataReadyFrameNo;
	std::thread *denoiseThread = nullptr;
	void *oidnLibrary = nullptr;
	bx::RngMwc rng;
	// shaders
	bgfx::ShaderHandle fs_atomicCounterClear;
	bgfx::ShaderHandle fs_lightmapClear;
	bgfx::ShaderHandle fs_lightmapAverage;
	bgfx::ShaderHandle fs_rayBundleClear;
	bgfx::ShaderHandle fs_rayBundleIntegrate;
	bgfx::ShaderHandle fs_rayBundleWrite;
	// programs
	bgfx::ProgramHandle atomicCounterClearProgram;
	bgfx::ProgramHandle lightmapClearProgram;
	bgfx::ProgramHandle lightmapAverageProgram;
	bgfx::ProgramHandle rayBundleClearProgram;
	bgfx::ProgramHandle rayBundleIntegrateProgram;
	bgfx::ProgramHandle rayBundleWriteProgram;
	// uniforms
	bgfx::UniformHandle u_clearLightmaps;
	bgfx::UniformHandle u_lightmapSize_dataSize;
	bgfx::UniformHandle u_rayNormal;
	bgfx::UniformHandle u_skyColor_enabled;
	bgfx::UniformHandle u_atomicCounterSampler;
	bgfx::UniformHandle u_rayBundleHeaderSampler;
	bgfx::UniformHandle u_rayBundleDataSampler;
	bgfx::UniformHandle u_rayBundleLightmapSampler;
	bgfx::UniformHandle u_lightmapSampler;
	// atomic counter
	bgfx::FrameBufferHandle atomicCounterFb;
	bgfx::TextureHandle atomicCounterTexture;
	// ray bundle write
	bgfx::FrameBufferHandle rayBundleFb;
	bgfx::TextureHandle rayBundleTarget;
	bgfx::TextureHandle rayBundleHeader;
	bgfx::TextureHandle rayBundleData;
	// ray bundle integrate
	bgfx::FrameBufferHandle rayBundleIntegrateFb;
	bgfx::TextureHandle rayBundleIntegrateTarget;
	bgfx::TextureHandle rayBundleLightmap; // Ray bundle is integrated into this. Cleared every ray bundle.
	// lightmap clear
	bgfx::FrameBufferHandle lightmapClearFb;
	bgfx::TextureHandle lightmapClearTarget;
	// lightmap average
	bgfx::FrameBufferHandle lightmapAverageFb;
	bgfx::TextureHandle lightmapAverageTarget;
	bgfx::TextureHandle lightmap; // Ray bundle lightmap is averaged (rgb / a). Only needs to be done at end of bake, but for visualization it's done every frame, erasing the previous frame result.
}
s_bake;

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

struct ScreenSpaceVertex
{
	float pos[2];
	float texcoord[2];
	static bgfx::VertexDecl decl;

	static void init()
	{
		decl.begin()
			.add(bgfx::Attrib::Position, 2, bgfx::AttribType::Float)
			.add(bgfx::Attrib::TexCoord0, 2, bgfx::AttribType::Float)
			.end();
		assert(decl.getStride() == sizeof(ScreenSpaceVertex));
	}
};

bgfx::VertexDecl ScreenSpaceVertex::decl;

void bakeInit()
{
	s_bake.enabled = (bgfx::getCaps()->supported & BGFX_CAPS_FRAMEBUFFER_RW) != 0;
	if (!s_bake.enabled) {
		printf("Read/Write frame buffer attachments are not supported. Baking is disabled.\n");
		return;
	}
	ScreenSpaceVertex::init();
	bx::mtxOrtho(s_bake.fsOrtho, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f, bgfx::getCaps()->homogeneousDepth);
}

void bakeShutdown()
{
	if (!s_bake.initialized)
		return;
	if (s_bake.oidnLibrary)
		bx::dlclose(s_bake.oidnLibrary);
	// shaders
	bgfx::destroy(s_bake.fs_atomicCounterClear);
	bgfx::destroy(s_bake.fs_lightmapClear);
	bgfx::destroy(s_bake.fs_lightmapAverage);
	bgfx::destroy(s_bake.fs_rayBundleClear);
	bgfx::destroy(s_bake.fs_rayBundleIntegrate);
	bgfx::destroy(s_bake.fs_rayBundleWrite);
	// programs
	bgfx::destroy(s_bake.atomicCounterClearProgram);
	bgfx::destroy(s_bake.rayBundleClearProgram);
	bgfx::destroy(s_bake.rayBundleWriteProgram);
	bgfx::destroy(s_bake.rayBundleIntegrateProgram);
	bgfx::destroy(s_bake.lightmapClearProgram);
	bgfx::destroy(s_bake.lightmapAverageProgram);
	// uniforms
	bgfx::destroy(s_bake.u_clearLightmaps);
	bgfx::destroy(s_bake.u_lightmapSize_dataSize);
	bgfx::destroy(s_bake.u_rayNormal);
	bgfx::destroy(s_bake.u_skyColor_enabled);
	bgfx::destroy(s_bake.u_atomicCounterSampler);
	bgfx::destroy(s_bake.u_rayBundleHeaderSampler);
	bgfx::destroy(s_bake.u_rayBundleDataSampler);
	bgfx::destroy(s_bake.u_rayBundleLightmapSampler);
	bgfx::destroy(s_bake.u_lightmapSampler);
	// framebuffers
	bgfx::destroy(s_bake.atomicCounterFb);
	bgfx::destroy(s_bake.rayBundleFb);
	bgfx::destroy(s_bake.rayBundleTarget);
	bgfx::destroy(s_bake.rayBundleHeader);
	bgfx::destroy(s_bake.rayBundleData);
	bgfx::destroy(s_bake.rayBundleLightmap);
	bgfx::destroy(s_bake.rayBundleIntegrateFb);
	bgfx::destroy(s_bake.rayBundleIntegrateTarget);
	bgfx::destroy(s_bake.lightmapClearFb);
	bgfx::destroy(s_bake.lightmapClearTarget);
	bgfx::destroy(s_bake.lightmapAverageFb);
	bgfx::destroy(s_bake.lightmapAverageTarget);
	bgfx::destroy(s_bake.lightmap);
}

static void setScreenSpaceQuadVertexBuffer()
{
	const uint32_t nVerts = 3;
	if (bgfx::getAvailTransientVertexBuffer(nVerts, ScreenSpaceVertex::decl) < nVerts)
		return;
	bgfx::TransientVertexBuffer vb;
	bgfx::allocTransientVertexBuffer(&vb, nVerts, ScreenSpaceVertex::decl);
	auto vertices = (ScreenSpaceVertex *)vb.data;
	vertices[0].pos[0] = -1.0f;
	vertices[0].pos[1] = 0.0f;
	vertices[1].pos[0] = 1.0f;
	vertices[1].pos[1] = 0.0f;
	vertices[2].pos[0] = 1.0f;
	vertices[2].pos[1] = 2.0f;
	bgfx::setVertexBuffer(0, &vb);
}

void bakeExecute()
{
	if (!(s_bake.status == BakeStatus::Idle || s_bake.status == BakeStatus::Finished))
		return;
	bakeClear();
	if (!s_bake.initialized) {
		// shaders
		s_bake.u_clearLightmaps = bgfx::createUniform("u_clearLightmaps", bgfx::UniformType::Vec4);
		s_bake.u_lightmapSize_dataSize = bgfx::createUniform("u_lightmapSize_dataSize", bgfx::UniformType::Vec4);
		s_bake.u_rayNormal = bgfx::createUniform("u_rayNormal", bgfx::UniformType::Vec4);
		s_bake.u_skyColor_enabled = bgfx::createUniform("u_skyColor_enabled", bgfx::UniformType::Vec4);
		s_bake.u_atomicCounterSampler = bgfx::createUniform("u_atomicCounterSampler", bgfx::UniformType::Sampler);
		s_bake.u_rayBundleHeaderSampler = bgfx::createUniform("u_rayBundleHeaderSampler", bgfx::UniformType::Sampler);
		s_bake.u_rayBundleDataSampler = bgfx::createUniform("u_rayBundleDataSampler", bgfx::UniformType::Sampler);
		s_bake.u_rayBundleLightmapSampler = bgfx::createUniform("u_rayBundleLightmapSampler", bgfx::UniformType::Sampler);
		s_bake.u_lightmapSampler = bgfx::createUniform("u_lightmapSampler", bgfx::UniformType::Sampler);
		s_bake.fs_atomicCounterClear = loadShader(ShaderId::fs_atomicCounterClear);
		s_bake.fs_lightmapClear = loadShader(ShaderId::fs_lightmapClear);
		s_bake.fs_lightmapAverage = loadShader(ShaderId::fs_lightmapAverage);
		s_bake.fs_rayBundleClear = loadShader(ShaderId::fs_rayBundleClear);
		s_bake.fs_rayBundleIntegrate = loadShader(ShaderId::fs_rayBundleIntegrate);
		s_bake.fs_rayBundleWrite = loadShader(ShaderId::fs_rayBundleWrite);
		s_bake.atomicCounterClearProgram = bgfx::createProgram(modelGet_vs_position(), s_bake.fs_atomicCounterClear);
		s_bake.lightmapClearProgram = bgfx::createProgram(modelGet_vs_position(), s_bake.fs_lightmapClear);
		s_bake.lightmapAverageProgram = bgfx::createProgram(modelGet_vs_position(), s_bake.fs_lightmapAverage);
		s_bake.rayBundleClearProgram = bgfx::createProgram(modelGet_vs_position(), s_bake.fs_rayBundleClear);
		s_bake.rayBundleIntegrateProgram = bgfx::createProgram(modelGet_vs_position(), s_bake.fs_rayBundleIntegrate);
		s_bake.rayBundleWriteProgram = bgfx::createProgram(modelGet_vs_model(), s_bake.fs_rayBundleWrite);
		// Atomic counter
		bgfx::TextureHandle target = bgfx::createTexture2D(1, 1, false, 1, bgfx::TextureFormat::RGBA8, BGFX_TEXTURE_RT);
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
			bgfx::destroy(s_bake.rayBundleFb);
			bgfx::destroy(s_bake.rayBundleTarget);
			bgfx::destroy(s_bake.rayBundleHeader);
			bgfx::destroy(s_bake.rayBundleData);
		}
		s_bake.rayBundleTarget = bgfx::createTexture2D((uint16_t)s_bake.resolution, (uint16_t)s_bake.resolution, false, 1, bgfx::TextureFormat::RGBA8, BGFX_TEXTURE_RT);
		s_bake.rayBundleHeader = bgfx::createTexture2D((uint16_t)s_bake.resolution, (uint16_t)s_bake.resolution, false, 1, bgfx::TextureFormat::R32U, BGFX_TEXTURE_COMPUTE_WRITE | BGFX_SAMPLER_POINT | BGFX_SAMPLER_UVW_CLAMP);
		s_bake.rayBundleData = bgfx::createTexture2D(s_bake.rbDataTextureSize, s_bake.rbDataTextureSize, false, 1, bgfx::TextureFormat::RGBA32U, BGFX_TEXTURE_COMPUTE_WRITE | BGFX_SAMPLER_POINT | BGFX_SAMPLER_UVW_CLAMP);
		bgfx::Attachment attachments[4];
		attachments[0].init(s_bake.rayBundleTarget);
		attachments[1].init(s_bake.atomicCounterTexture, bgfx::Access::ReadWrite);
		attachments[2].init(s_bake.rayBundleHeader, bgfx::Access::ReadWrite);
		attachments[3].init(s_bake.rayBundleData, bgfx::Access::ReadWrite);
		s_bake.rayBundleFb = bgfx::createFrameBuffer(BX_COUNTOF(attachments), attachments);
	}
	// Re-create lightmap if atlas resolution has changed.
	const bool lightmapResolutionChanged = s_bake.lightmapWidth != atlasGetWidth() || s_bake.lightmapHeight != atlasGetHeight();
	if (!s_bake.initialized || lightmapResolutionChanged) {
		if (s_bake.initialized) {
			bgfx::destroy(s_bake.lightmapClearFb);
			bgfx::destroy(s_bake.lightmapClearTarget);
			bgfx::destroy(s_bake.rayBundleLightmap);
			bgfx::destroy(s_bake.lightmapAverageFb);
			bgfx::destroy(s_bake.lightmapAverageTarget);
			bgfx::destroy(s_bake.lightmap);
		}
		s_bake.lightmapWidth = atlasGetWidth();
		s_bake.lightmapHeight = atlasGetHeight();
		s_bake.rayBundleLightmap = bgfx::createTexture2D((uint16_t)s_bake.lightmapWidth * 4, (uint16_t)s_bake.lightmapHeight, false, 1, bgfx::TextureFormat::R32U, BGFX_TEXTURE_COMPUTE_WRITE);
		s_bake.lightmap = bgfx::createTexture2D((uint16_t)s_bake.lightmapWidth, (uint16_t)s_bake.lightmapHeight, false, 1, bgfx::TextureFormat::RGBA32F, BGFX_TEXTURE_COMPUTE_WRITE);
		{
			s_bake.lightmapClearTarget = bgfx::createTexture2D((uint16_t)s_bake.lightmapWidth, (uint16_t)s_bake.lightmapHeight, false, 1, bgfx::TextureFormat::RGBA8, BGFX_TEXTURE_RT);
			bgfx::Attachment attachments[3];
			attachments[0].init(s_bake.lightmapClearTarget);
			attachments[1].init(s_bake.rayBundleLightmap, bgfx::Access::ReadWrite);
			attachments[2].init(s_bake.lightmap, bgfx::Access::ReadWrite);
			s_bake.lightmapClearFb = bgfx::createFrameBuffer(BX_COUNTOF(attachments), attachments);
		}
		{
			s_bake.lightmapAverageTarget = bgfx::createTexture2D((uint16_t)s_bake.lightmapWidth, (uint16_t)s_bake.lightmapHeight, false, 1, bgfx::TextureFormat::RGBA8, BGFX_TEXTURE_RT);
			bgfx::Attachment attachments[3];
			attachments[0].init(s_bake.lightmapAverageTarget);
			attachments[1].init(s_bake.rayBundleLightmap, bgfx::Access::ReadWrite);
			attachments[2].init(s_bake.lightmap, bgfx::Access::ReadWrite);
			s_bake.lightmapAverageFb = bgfx::createFrameBuffer(BX_COUNTOF(attachments), attachments);
		}
	}
	// Re-create ray bundle integrate fb if ray bundle or atlas resolution has changed.
	if (!s_bake.initialized || rayBundleResolutionChanged || lightmapResolutionChanged) {
		if (s_bake.initialized) {
			bgfx::destroy(s_bake.rayBundleIntegrateFb);
			bgfx::destroy(s_bake.rayBundleIntegrateTarget);
		}
		s_bake.rayBundleIntegrateTarget = bgfx::createTexture2D((uint16_t)s_bake.resolution, (uint16_t)s_bake.resolution, false, 1, bgfx::TextureFormat::RGBA8, BGFX_TEXTURE_RT);
		bgfx::Attachment attachments[4];
		attachments[0].init(s_bake.rayBundleIntegrateTarget);
		attachments[1].init(s_bake.rayBundleHeader, bgfx::Access::Read);
		attachments[2].init(s_bake.rayBundleData, bgfx::Access::Read);
		attachments[3].init(s_bake.rayBundleLightmap, bgfx::Access::ReadWrite);
		s_bake.rayBundleIntegrateFb = bgfx::createFrameBuffer(BX_COUNTOF(attachments), attachments);
	}
	s_bake.initialized = true;
	s_bake.status = BakeStatus::Executing;
	s_bake.directionCount = 0;
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

void bakeFrame(uint32_t bgfxFrame)
{
	bgfx::ViewId viewId = kFirstFreeView;
	if (s_bake.status == BakeStatus::Executing) {
		for (uint32_t i = 0; i < (uint32_t)s_bake.options.directionsPerFrame; i++) {
			if (s_bake.directionCount == 0) {
				// Lightmap clear integrate.
				bgfx::setViewFrameBuffer(viewId, s_bake.lightmapClearFb);
				bgfx::setViewRect(viewId, 0, 0, (uint16_t)s_bake.lightmapWidth * 4, (uint16_t)s_bake.lightmapHeight);
				bgfx::setViewTransform(viewId, nullptr, s_bake.fsOrtho);
				bgfx::setTexture(1, s_bake.u_rayBundleLightmapSampler, s_bake.rayBundleLightmap);
				bgfx::setTexture(2, s_bake.u_lightmapSampler, s_bake.lightmap);
				const float clear[] = { 1.0f, 0.0f, 0.0f, 0.0f };
				bgfx::setUniform(s_bake.u_clearLightmaps, clear);
				setScreenSpaceQuadVertexBuffer();
				bgfx::setState(0);
				bgfx::submit(viewId, s_bake.lightmapClearProgram);
				viewId++;
			}
			// Atomic counter clear.
			bgfx::setViewFrameBuffer(viewId, s_bake.atomicCounterFb);
			bgfx::setViewRect(viewId, 0, 0, 1, 1);
			bgfx::setViewTransform(viewId, nullptr, s_bake.fsOrtho);
			bgfx::setTexture(1, s_bake.u_atomicCounterSampler, s_bake.atomicCounterTexture);
			setScreenSpaceQuadVertexBuffer();
			bgfx::setState(0);
			bgfx::submit(viewId, s_bake.atomicCounterClearProgram);
			viewId++;
			// Ray bundle clear.
			bgfx::setViewFrameBuffer(viewId, s_bake.rayBundleFb);
			bgfx::setViewRect(viewId, 0, 0, (uint16_t)s_bake.resolution, (uint16_t)s_bake.resolution);
			bgfx::setViewTransform(viewId, nullptr, s_bake.fsOrtho);
			bgfx::setTexture(2, s_bake.u_rayBundleHeaderSampler, s_bake.rayBundleHeader);
			setScreenSpaceQuadVertexBuffer();
			bgfx::setState(0);
			bgfx::submit(viewId, s_bake.rayBundleClearProgram);
			viewId++;
			// Ray bundle write.
			const bx::Vec3 forward = bx::randUnitHemisphere(&s_bake.rng, bx::Vec3(0.0f, 0.0f, 1.0f));
			const bx::Vec3 right = bx::randUnitHemisphere(&s_bake.rng, bx::Vec3(1.0f, 0.0f, 0.0f));
			const bx::Vec3 up = bx::cross(forward, right);
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
				const float sizes[] = { (float)s_bake.lightmapWidth, (float)s_bake.lightmapHeight, (float)s_bake.rbDataTextureSize, 0.0f };
				const objzModel *model = modelGetData();
				for (uint32_t j = 0; j < model->numMeshes; j++) {
					const objzMesh &mesh = model->meshes[j];
					const objzMaterial *mat = mesh.materialIndex == -1 ? nullptr : &model->materials[mesh.materialIndex];
					bgfx::setIndexBuffer(atlasGetIb(), mesh.firstIndex, mesh.numIndices);
					bgfx::setVertexBuffer(0, atlasGetVb());
					bgfx::setState(0);
					bgfx::setTexture(1, s_bake.u_atomicCounterSampler, s_bake.atomicCounterTexture);
					bgfx::setTexture(2, s_bake.u_rayBundleHeaderSampler, s_bake.rayBundleHeader);
					bgfx::setTexture(3, s_bake.u_rayBundleDataSampler, s_bake.rayBundleData);
					bgfx::setUniform(s_bake.u_lightmapSize_dataSize, sizes);
					modelSetMaterialUniforms(mat);
					bgfx::submit(viewId, s_bake.rayBundleWriteProgram);
				}
			}
			viewId++;
			// Ray bundle integrate.
			bgfx::setViewFrameBuffer(viewId, s_bake.rayBundleIntegrateFb);
			bgfx::setViewRect(viewId, 0, 0, (uint16_t)s_bake.resolution, (uint16_t)s_bake.resolution);
			bgfx::setViewTransform(viewId, nullptr, s_bake.fsOrtho);
			bgfx::setTexture(1, s_bake.u_rayBundleHeaderSampler, s_bake.rayBundleHeader);
			bgfx::setTexture(2, s_bake.u_rayBundleDataSampler, s_bake.rayBundleData);
			bgfx::setTexture(3, s_bake.u_rayBundleLightmapSampler, s_bake.rayBundleLightmap);
			const float sizes[] = { (float)s_bake.lightmapWidth, (float)s_bake.lightmapHeight, (float)s_bake.rbDataTextureSize, 0.0f };
			bgfx::setUniform(s_bake.u_lightmapSize_dataSize, sizes);
			const float rayNormal[] = { -view[2], -view[6], -view[10], 0 };
			bgfx::setUniform(s_bake.u_rayNormal, rayNormal);
			const float sky[] = { s_bake.options.skyColor.x, s_bake.options.skyColor.y, s_bake.options.skyColor.z, s_bake.options.sky ? 1.0f : 0.0f };
			bgfx::setUniform(s_bake.u_skyColor_enabled, sky);
			setScreenSpaceQuadVertexBuffer();
			bgfx::setState(0);
			bgfx::submit(viewId, s_bake.rayBundleIntegrateProgram);
			viewId++;
			// Lightmap clear average.
			bgfx::setViewFrameBuffer(viewId, s_bake.lightmapClearFb);
			bgfx::setViewRect(viewId, 0, 0, (uint16_t)s_bake.lightmapWidth, (uint16_t)s_bake.lightmapHeight);
			bgfx::setViewTransform(viewId, nullptr, s_bake.fsOrtho);
			bgfx::setTexture(1, s_bake.u_rayBundleLightmapSampler, s_bake.rayBundleLightmap);
			bgfx::setTexture(2, s_bake.u_lightmapSampler, s_bake.lightmap);
			const float clear[] = { 0.0f, 1.0f, 0.0f, 0.0f };
			bgfx::setUniform(s_bake.u_clearLightmaps, clear);
			setScreenSpaceQuadVertexBuffer();
			bgfx::setState(0);
			bgfx::submit(viewId, s_bake.lightmapClearProgram);
			viewId++;
			// Lightmap average.
			bgfx::setViewFrameBuffer(viewId, s_bake.lightmapAverageFb);
			bgfx::setViewRect(viewId, 0, 0, (uint16_t)s_bake.lightmapWidth, (uint16_t)s_bake.lightmapHeight);
			bgfx::setViewTransform(viewId, nullptr, s_bake.fsOrtho);
			bgfx::setTexture(1, s_bake.u_rayBundleLightmapSampler, s_bake.rayBundleLightmap);
			bgfx::setTexture(2, s_bake.u_lightmapSampler, s_bake.lightmap);
			setScreenSpaceQuadVertexBuffer();
			bgfx::setState(0);
			bgfx::submit(viewId, s_bake.lightmapAverageProgram);
			viewId++;
			// Finished with this direction.
			s_bake.directionCount++;
			if (s_bake.directionCount >= s_bake.options.numDirections) {
				// Finished rendering.
				if (s_bake.options.denoise) {
					s_bake.status = BakeStatus::ReadingLightmap;
					s_bake.lightmapData.resize(s_bake.lightmapWidth * s_bake.lightmapHeight * 4 * sizeof(float));
					s_bake.lightmapDataReadyFrameNo = bgfx::readTexture(s_bake.lightmap, s_bake.lightmapData.data());
				} else {
					s_bake.status = BakeStatus::Finished;
				}
				return;
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
		bgfx::updateTexture2D(s_bake.lightmap, 0, 0, 0, 0, (uint16_t)s_bake.lightmapWidth, (uint16_t)s_bake.lightmapHeight, bgfx::makeRef(s_bake.denoisedLightmapData.data(), (uint32_t)s_bake.denoisedLightmapData.size() * sizeof(float)));
		s_bake.status = BakeStatus::Finished;
	}
}

void bakeClear()
{
	s_bake.status = BakeStatus::Idle;
	g_options.shadeMode = atlasIsReady() ? ShadeMode::Charts : ShadeMode::Flat;
}

void bakeShowGuiOptions()
{
	if (!s_bake.enabled)
		return;
	const ImVec2 buttonSize(ImVec2(ImGui::GetContentRegionAvailWidth() * 0.3f, 0.0f));
	ImGui::Text("Lightmap");
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
	ImGui::SliderInt("Ray bundle directions", &s_bake.options.numDirections, 300, 10000);
	ImGui::SliderInt("Directions per frame", &s_bake.options.directionsPerFrame, 1, 100);
	if (s_bake.status == BakeStatus::Idle || s_bake.status == BakeStatus::Finished) {
		if (ImGui::Button("Bake", buttonSize))
			bakeExecute();
	}
	else {
		if (s_bake.directionCount < s_bake.options.numDirections)
			ImGui::ProgressBar(s_bake.directionCount / (float)s_bake.options.numDirections);
		else {
			ImGui::AlignTextToFramePadding();
			ImGui::Text("Denoising...");
		}
	}
	if (s_bake.status != BakeStatus::Idle)
		ImGui::Checkbox("Show lightmap", &s_bake.options.showLightmap);
}

void bakeShowGuiWindow()
{
	if (s_bake.status == BakeStatus::Idle || !s_bake.options.showLightmap)
		return;
	ImGuiIO &io = ImGui::GetIO();
	const float size = 500;
	const float margin = 4.0f;
	ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x - size - margin, size + margin * 2.0f), ImGuiCond_FirstUseEver);
	ImGui::SetNextWindowSize(ImVec2(size, size), ImGuiCond_FirstUseEver);
	if (ImGui::Begin("Lightmap", &s_bake.options.showLightmap)) {
		const ImVec2 cursorPos = ImGui::GetCursorScreenPos();
		ImTextureID texture = (ImTextureID)(intptr_t)bakeGetLightmap().idx;
		ImGui::Image(texture, ImGui::GetContentRegionAvail());
		if (ImGui::IsItemHovered())
			guiImageMagnifierTooltip(texture, cursorPos, ImVec2((float)s_bake.lightmapWidth, (float)s_bake.lightmapHeight));
		ImGui::End();
	}
}

bgfx::TextureHandle bakeGetLightmap()
{
	return s_bake.lightmap;
}

bool bakeIsIdle()
{
	return s_bake.status == BakeStatus::Idle;
}
