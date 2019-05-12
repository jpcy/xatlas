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
#include <vector>
#include <assert.h>
#include <time.h>
#include <bx/os.h>
#include <bx/rng.h>
#include <embree3/rtcore.h>
#include <imgui/imgui.h>
#include <OpenImageDenoise/oidn.h>
#include "shaders/shared.h"
#include "viewer.h"

struct BakeStatus
{
	enum Enum
	{
		Idle,
		Rasterizing,
		Tracing,
		Denoising,
		ThreadFinished,
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
	int numSamples = 4;
	int numBounces = 1;
};

struct
{
	bool enabled;
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
	clock_t lastUpdateTime = 0;
	const double updateIntervalMs = 50;
	// worker thread
	std::thread *workerThread = nullptr;
	std::vector<float> lightmapData;
	std::vector<float> denoisedLightmapData;
	std::atomic<uint32_t> numTrianglesRasterized;
	bool denoiseSucceeded;
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

static void bakeRasterize()
{
	s_bake.status = BakeStatus::Rasterizing;
	s_bake.numTrianglesRasterized = 0;
	const objzModel *model = modelGetData();
	for (uint32_t tri = 0; tri < model->numIndices / 3; tri++) {
		s_bake.numTrianglesRasterized++;
	}
}

static void bakeTraceRays()
{
	s_bake.status = BakeStatus::Tracing;
	s_bake.lightmapData.resize(s_bake.lightmapWidth * s_bake.lightmapHeight * 4);
}

static void bakeDenoise()
{
	s_bake.denoiseSucceeded = false;
	if (!s_bake.options.denoise)
		return;
	s_bake.status = BakeStatus::Denoising;
	if (!s_bake.oidnLibrary) {
		s_bake.oidnLibrary = bx::dlopen("OpenImageDenoise.dll");
		if (!s_bake.oidnLibrary)
			return;
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
		return;
	}
	oidn::SetDevice1b(device, "setAffinity", false);
	oidn::CommitDevice(device);
	OIDNFilter filter = oidn::NewFilter(device, "RT");
	oidn::SetSharedFilterImage(filter, "color", input.data(), OIDN_FORMAT_FLOAT3, s_bake.lightmapWidth, s_bake.lightmapHeight, 0, 0, 0);
	oidn::SetSharedFilterImage(filter, "output", output.data(), OIDN_FORMAT_FLOAT3, s_bake.lightmapWidth, s_bake.lightmapHeight, 0, 0, 0);
	oidn::CommitFilter(filter);
	oidn::ExecuteFilter(filter);
	const char *errorMessage;
	if (oidn::GetDeviceError(device, &errorMessage) != OIDN_ERROR_NONE) {
		fprintf(stderr, "Denoiser error: %s\n", errorMessage);
		return;
	}
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
	s_bake.denoiseSucceeded = true;
}

static void bakeWorkerThread()
{
	bakeRasterize();
	bakeTraceRays();
	bakeDenoise();
	s_bake.status = BakeStatus::ThreadFinished;
}

static void bakeShutdownWorkerThread()
{
	if (s_bake.workerThread) {
		s_bake.workerThread->join();
		delete s_bake.workerThread;
		s_bake.workerThread = nullptr;
	}
}

static void bakeEmbreeError(void* /*userPtr*/, enum RTCError /*code*/, const char* str)
{
	fprintf(stderr, "Embree error: %s\n", str);
	exit(EXIT_FAILURE);
}

void bakeInit()
{
	s_bake.embreeLibrary = bx::dlopen(EMBREE_LIB);
	if (!s_bake.embreeLibrary) {
		printf("Embree not installed. Baking is disabled.\n");
		s_bake.enabled = false;
		return;
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
	s_bake.enabled = true;
}

void bakeShutdown()
{
	if (s_bake.embreeLibrary)
		bx::dlclose(s_bake.embreeLibrary);
	if (s_bake.oidnLibrary)
		bx::dlclose(s_bake.oidnLibrary);
	if (s_bake.initialized)
		bgfx::destroy(s_bake.lightmap);
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
	// Re-create lightmap if atlas resolution has changed.
	const bool lightmapResolutionChanged = s_bake.lightmapWidth != atlasGetWidth() || s_bake.lightmapHeight != atlasGetHeight();
	if (!s_bake.initialized || lightmapResolutionChanged) {
		s_bake.lightmapWidth = atlasGetWidth();
		s_bake.lightmapHeight = atlasGetHeight();
		if (s_bake.initialized)
			bgfx::destroy(s_bake.lightmap);
		s_bake.lightmap = bgfx::createTexture2D((uint16_t)s_bake.lightmapWidth, (uint16_t)s_bake.lightmapHeight, false, 1, bgfx::TextureFormat::RGBA32F, BGFX_TEXTURE_COMPUTE_WRITE | BGFX_SAMPLER_UVW_CLAMP);
	}
	s_bake.initialized = true;
	g_options.shadeMode = ShadeMode::Lightmap;
	s_bake.workerThread = new std::thread(bakeWorkerThread);
}

void bakeFrame()
{
	if (s_bake.status == BakeStatus::Tracing) {
		const double elapsedMs = (clock() - s_bake.lastUpdateTime) * 1000.0 / CLOCKS_PER_SEC;
		if (elapsedMs >= s_bake.updateIntervalMs) {
			// TODO
			s_bake.lastUpdateTime = clock();
		}
	} else if (s_bake.status == BakeStatus::ThreadFinished) {
		bakeShutdownWorkerThread();
		std::vector<float> *data = s_bake.denoiseSucceeded ? &s_bake.denoisedLightmapData : &s_bake.lightmapData;
		bgfx::updateTexture2D(s_bake.lightmap, 0, 0, 0, 0, (uint16_t)s_bake.lightmapWidth, (uint16_t)s_bake.lightmapHeight, bgfx::makeRef(data->data(), (uint32_t)data->size() * sizeof(float)));
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
		static const char * const sampleLabels[] = { "1", "4", "8", "16" };
		static const int samples[] = { 1, 4, 8, 16 };
		int sampleIndex = 0;
		for (int i = 0; i < (int)BX_COUNTOF(samples); i++) {
			if (samples[i] == s_bake.options.numSamples) {
				sampleIndex = i;
				break;
			}
		}
		ImGui::ListBox("Samples", &sampleIndex, sampleLabels, (int)BX_COUNTOF(sampleLabels));
		s_bake.options.numSamples = samples[sampleIndex];
		ImGui::SliderInt("Bounces", &s_bake.options.numBounces, 0, 4);
		const ImVec2 buttonSize(ImVec2(ImGui::GetContentRegionAvailWidth() * 0.3f, 0.0f));
		if (ImGui::Button("Bake", buttonSize))
			bakeExecute();
	} else if (s_bake.status == BakeStatus::Rasterizing) {
		ImGui::Text("Rasterizing...");
		ImGui::ProgressBar(s_bake.numTrianglesRasterized / float(modelGetData()->numIndices / 3));
		/*if (ImGui::Button("Cancel"))
			s_bake.status = BakeStatus::Finished;*/
	} else if (s_bake.status == BakeStatus::Tracing) {
		ImGui::Text("Tracing rays...");
	} else if (s_bake.status == BakeStatus::Denoising) {
		ImGui::AlignTextToFramePadding();
		ImGui::Text("Denoising...");
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

bool bakeIsIdle()
{
	return s_bake.status == BakeStatus::Idle;
}
