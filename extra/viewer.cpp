/*
xatlas
https://github.com/jpcy/xatlas
Copyright (c) 2018 Jonathan Young

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#include <mutex>
#include <thread>
#include <unordered_map>
#include <vector>
#include <stdio.h>
#include "bx/bx.h"
#include "bx/commandline.h"
#include "bx/math.h"
#include "bgfx/bgfx.h"
#include "bgfx/platform.h"
#include "GLFW/glfw3.h"
#if BX_PLATFORM_LINUX
#define GLFW_EXPOSE_NATIVE_X11
#elif BX_PLATFORM_WINDOWS
#define GLFW_EXPOSE_NATIVE_WIN32
#endif
#include "GLFW/glfw3native.h"
#undef Success
#include "imgui/imgui.h"
#include "nativefiledialog/nfd.h"
#include "objzero/objzero.h"
#include "../xatlas.h"

#define USE_LIBIGL 0

#if USE_LIBIGL
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4018)
#pragma warning(disable : 4127)
#pragma warning(disable : 4244)
#pragma warning(disable : 4267)
#pragma warning(disable : 4305)
#pragma warning(disable : 4427)
#pragma warning(disable : 4459)
#pragma warning(disable : 4702)
#pragma warning(disable : 4996)
#endif
#include "igl/arap.h"
#include "igl/boundary_loop.h"
#include "igl/harmonic.h"
#include "igl/lscm.h"
#include "igl/map_vertices_to_circle.h"
#ifdef _MSC_VER
#pragma warning(pop)
#endif
#endif

#include "shaders_bin/shaders.h"

#define WINDOW_TITLE "xatlas viewer"
#define WINDOW_DEFAULT_WIDTH 1920
#define WINDOW_DEFAULT_HEIGHT 1080

#ifdef _MSC_VER
#define STRNCPY(_dest, _destSize, _src) strncpy_s(_dest, _destSize, _src, (_destSize) - 1)
#else
#define STRNCPY(_dest, _destSize, _src) strncpy(_dest, _src, (_destSize) - 1)
#endif

struct AABB
{
	AABB() : min(FLT_MAX, FLT_MAX, FLT_MAX), max(-FLT_MAX, -FLT_MAX, -FLT_MAX) {}

	void addPoint(bx::Vec3 v)
	{
		min.x = bx::min(min.x, v.x);
		min.y = bx::min(min.y, v.y);
		min.z = bx::min(min.z, v.z);
		max.x = bx::max(max.x, v.x);
		max.y = bx::max(max.y, v.y);
		max.z = bx::max(max.z, v.z);
	}

	void getCorners(bx::Vec3 *corners) const
	{
		const bx::Vec3 aabb[] = { min, max };
		for (int i = 0; i < 8; i++) {
			corners[i].x = aabb[i & 1].x;
			corners[i].y = aabb[(i >> 1) & 1].y;
			corners[i].z = aabb[(i >> 2) & 1].z;
		}
	}

	bx::Vec3 min;
	bx::Vec3 max;
};

static const bgfx::ViewId kModelView = 0;
static const bgfx::ViewId kGuiView = 1;
static const bgfx::ViewId kAtomicCounterClearView = 2;
static const bgfx::ViewId kRayBundleClearView = 3;
static const bgfx::ViewId kRayBundleWriteView = 4;
static const bgfx::ViewId kRayBundleResolveView = 5;

static const bgfx::ViewId kNumBakeViews = 4;

static const uint8_t kPaletteBlack = 0;

static GLFWwindow *s_window;
static bool s_keyDown[GLFW_KEY_LAST + 1] = { 0 };
static bool s_showBgfxStats = false;

struct
{
	bgfx::VertexDecl vertexDecl;
	bgfx::TextureHandle font;
	bgfx::ProgramHandle program;
	bgfx::UniformHandle u_texture;
}
s_gui;

struct ModelStatus
{
	enum Enum
	{
		NotLoaded,
		Loading,
		Finalizing,
		Ready
	};

	Enum get()
	{
		m_lock.lock();
		Enum result = m_value;
		m_lock.unlock();
		return result;
	}

	void set(Enum value)
	{
		m_lock.lock();
		m_value = value;
		m_lock.unlock();
	}

private:
	std::mutex m_lock;
	Enum m_value = NotLoaded;
};

struct ModelVertex
{
	bx::Vec3 pos;
	bx::Vec3 normal;
	float texcoord[4];
	static bgfx::VertexDecl decl;

	static void init()
	{
		decl.begin()
			.add(bgfx::Attrib::Position, 3, bgfx::AttribType::Float)
			.add(bgfx::Attrib::Normal, 3, bgfx::AttribType::Float)
			.add(bgfx::Attrib::TexCoord0, 4, bgfx::AttribType::Float)
			.end();
		assert(decl.getStride() == sizeof(ModelVertex));
	}
};

bgfx::VertexDecl ModelVertex::decl;

struct
{
	ModelStatus status;
	std::thread *thread = nullptr;
	objzModel *data;
	AABB aabb;
	bx::Vec3 centroid = bx::Vec3(0.0f, 0.0f, 0.0f);
	bgfx::VertexBufferHandle vb = BGFX_INVALID_HANDLE;
	bgfx::IndexBufferHandle ib = BGFX_INVALID_HANDLE;
	bgfx::IndexBufferHandle wireframeIb = BGFX_INVALID_HANDLE;
	float scale = 1.0f;
	bgfx::ProgramHandle colorProgram;
	bgfx::ProgramHandle flatProgram;
	bgfx::ProgramHandle checkerboardProgram;
	bgfx::UniformHandle u_diffuse;
	bgfx::UniformHandle u_emission;
	bgfx::UniformHandle u_lightDir;
	bgfx::UniformHandle u_color;
	bgfx::UniformHandle u_textureSize_cellSize;
}
s_model;

struct AtlasStatus
{
	enum Enum
	{
		NotGenerated,
		AddingMeshes,
		Generating,
		Finalizing,
		Ready
	};

	Enum get()
	{
		m_lock.lock();
		Enum result = m_value;
		m_lock.unlock();
		return result;
	}

	void set(Enum value)
	{
		m_lock.lock();
		m_value = value;
		m_lock.unlock();
	}

	void getProgress(xatlas::ProgressCategory::Enum *category, int *progress)
	{
		m_lock.lock();
		*category = m_category;
		*progress = m_progress;
		m_lock.unlock();
	}

	void setProgress(xatlas::ProgressCategory::Enum category, int progress)
	{
		m_lock.lock();
		m_category = category;
		m_progress = progress;
		m_lock.unlock();
	}

private:
	std::mutex m_lock;
	Enum m_value = NotGenerated;
	xatlas::ProgressCategory::Enum m_category;
	int m_progress = 0;
};

enum class ParamMethod
{
	LSCM,
	libigl_Harmonic,
	libigl_LSCM,
	libigl_ARAP
};

struct
{
	const int chartColorSeed = 13;
	const int chartCellSize = 16;
	xatlas::Atlas *data = nullptr;
	std::thread *thread = nullptr;
	AtlasStatus status;
	bool verbose = false;
	bool showTexture = true;
	int currentTexture;
	std::vector<bgfx::TextureHandle> chartsTextures;
	std::vector<std::vector<uint8_t>> chartsImages;
	bgfx::VertexBufferHandle vb = BGFX_INVALID_HANDLE;
	bgfx::IndexBufferHandle ib = BGFX_INVALID_HANDLE;
	bgfx::IndexBufferHandle chartIb = BGFX_INVALID_HANDLE;
	bgfx::VertexBufferHandle chartBoundaryVb = BGFX_INVALID_HANDLE;
	xatlas::ChartOptions chartOptions;
	xatlas::PackOptions packOptions;
	ParamMethod paramMethod = ParamMethod::LSCM;
	bool paramMethodChanged = false;
	std::vector<ModelVertex> vertices;
	std::vector<uint32_t> indices;
	std::vector<uint32_t> chartIndices;
	std::vector<bx::Vec3> chartBoundaryVertices;
	bgfx::VertexDecl wireVertexDecl;
}
s_atlas;

enum class ShadeMode
{
	Flat,
	Charts
};

enum class WireframeMode
{
	Charts,
	Triangles
};

struct
{
	bool gui = true;
	bool wireframe = true;
	ShadeMode shadeMode = ShadeMode::Flat;
	WireframeMode wireframeMode = WireframeMode::Triangles;
}
s_options;

static void randomRGB(uint8_t *color)
{
	const int mix = 192;
	color[0] = uint8_t((rand() % 255 + mix) * 0.5f);
	color[1] = uint8_t((rand() % 255 + mix) * 0.5f);
	color[2] = uint8_t((rand() % 255 + mix) * 0.5f);
}

struct
{
	void set(const char *format, ...)
	{
		m_lock.lock();
		if (format) {
			va_list args;
			va_start(args, format);
			vsnprintf(m_text, sizeof(m_text), format, args);
			va_end(args);
		} else
			m_text[0] = 0;
		m_lock.unlock();
	}

	void get(char *buffer, size_t bufferSize)
	{
		m_lock.lock();
		STRNCPY(buffer, bufferSize, m_text);
		m_lock.unlock();
	}

private:
	char m_text[1024] = { 0 };
	std::mutex m_lock;
}
s_errorMessage;

static void axisFromEulerAngles(float pitch, float yaw, bx::Vec3 *forward, bx::Vec3 *right, bx::Vec3 *up)
{
	const float ryaw = bx::toRad(yaw);
	const float rpitch = bx::toRad(pitch);
	const bx::Vec3 f = bx::Vec3
	(
		bx::sin(ryaw) * bx::cos(rpitch),
		bx::sin(rpitch),
		bx::cos(ryaw) * bx::cos(rpitch)
	);
	const bx::Vec3 r = bx::Vec3
	(
		bx::sin(ryaw - float(bx::kPi * 0.5f)),
		0.0f,
		bx::cos(ryaw - float(bx::kPi * 0.5f))
	);
	if (forward)
		*forward = f;
	if (right)
		*right = r;
	if (up)
		*up = bx::mul(bx::cross(f, r), -1.0f);
}

static float cleanAngle(float degrees)
{
	while (degrees < 0.0f)
		degrees += 360.0f;
	while (degrees >= 360.0f)
		degrees -= 360.0f;
	return degrees;
}

struct FirstPersonCamera
{
	FirstPersonCamera() : position(bx::Vec3(0.0f, 0.0f, 0.0f)), pitch(0.0f), yaw(0.0f) {}

	void calculateViewMatrix(float *mat)
	{
		bx::Vec3 forward, up;
		axisFromEulerAngles(pitch, yaw, &forward, nullptr, &up);
		const bx::Vec3 at = bx::add(position, forward);
		bx::mtxLookAt(mat, position, at, up, bx::Handness::Right);
	}

	void move(float deltaForward, float deltaRight)
	{
		bx::Vec3 forward, right;
		axisFromEulerAngles(pitch, yaw, &forward, &right, nullptr);
		const bx::Vec3 velocity = bx::add(bx::mul(forward, deltaForward), bx::mul(right, deltaRight));
		position = bx::add(position, velocity);
	}

	void rotate(float deltaX, float deltaY)
	{
		yaw = cleanAngle(yaw + deltaX);
		pitch = bx::clamp(pitch + deltaY, -90.0f, 90.0f);
	}

	bx::Vec3 position;
	float pitch;
	float yaw;
};

struct OrbitCamera
{
	OrbitCamera() : distance(32.0f), pitch(0.0f), yaw(0.0f) {}

	void calculateViewMatrix(float *mat)
	{
		bx::Vec3 forward;
		axisFromEulerAngles(pitch, yaw, &forward, nullptr, nullptr);
		const bx::Vec3 center = bx::mul(s_model.centroid, s_model.scale);
		const bx::Vec3 eye = bx::add(bx::mul(forward, -distance), center);
		const bx::Vec3 up = bx::Vec3(0.0f, 1.0f, 0.0f);
		bx::mtxLookAt(mat, eye, center, up, bx::Handness::Right);
	}

	void rotate(float deltaX, float deltaY)
	{
		yaw = cleanAngle(yaw - deltaX);
		pitch = bx::clamp(pitch + deltaY, -75.0f, 75.0f);
	}

	void zoom(float delta)
	{
		distance = bx::clamp(distance + delta, 0.1f, 500.0f);
	}

	float distance;
	float pitch;
	float yaw;
};

enum class CameraMode
{
	FirstPerson,
	Orbit
};

struct
{
	CameraMode mode = CameraMode::Orbit;
	FirstPersonCamera firstPerson;
	OrbitCamera orbit;
	double lastCursorPos[2];
	float fov = 90.0f;
	float sensitivity = 0.25f;
}
s_camera;

static void glfw_errorCallback(int error, const char *description)
{
	fprintf(stderr, "GLFW error %d: %s\n", error, description);
}

static void glfw_charCallback(GLFWwindow * /*window*/, unsigned int c)
{
	if (!s_options.gui)
		return;
	if (c > 0 && c < 0x10000) {
		ImGuiIO &io = ImGui::GetIO();
		io.AddInputCharacter((unsigned short)c);
	}
}

static void glfw_cursorPosCallback(GLFWwindow * /*window*/, double xpos, double ypos)
{
	const float deltaX = float(xpos - s_camera.lastCursorPos[0]);
	const float deltaY = float(ypos - s_camera.lastCursorPos[1]);
	s_camera.lastCursorPos[0] = xpos;
	s_camera.lastCursorPos[1] = ypos;
	if (glfwGetInputMode(s_window, GLFW_CURSOR) == GLFW_CURSOR_DISABLED) {
		if (s_camera.mode == CameraMode::FirstPerson)
			s_camera.firstPerson.rotate(-deltaX * s_camera.sensitivity, -deltaY * s_camera.sensitivity);
		else if (s_camera.mode == CameraMode::Orbit)
			s_camera.orbit.rotate(deltaX * s_camera.sensitivity, -deltaY * s_camera.sensitivity);
	}
}

static void glfw_keyCallback(GLFWwindow * /*window*/, int key, int, int action, int /*mods*/)
{
	if (action == GLFW_REPEAT)
		return;
	if (key >= 0 && key <= GLFW_KEY_LAST)
		s_keyDown[key] = action == GLFW_PRESS;
	if (key == GLFW_KEY_F1 && action == GLFW_RELEASE)
		s_showBgfxStats = !s_showBgfxStats;
	if (key == GLFW_KEY_F2 && action == GLFW_RELEASE)
		s_options.gui = !s_options.gui;
	if (s_options.gui) {
		ImGuiIO &io = ImGui::GetIO();
		if (key >= 0 && key < 512)
			io.KeysDown[key] = action == GLFW_PRESS;
		io.KeyCtrl = io.KeysDown[GLFW_KEY_LEFT_CONTROL] || io.KeysDown[GLFW_KEY_RIGHT_CONTROL];
		io.KeyShift = io.KeysDown[GLFW_KEY_LEFT_SHIFT] || io.KeysDown[GLFW_KEY_RIGHT_SHIFT];
		io.KeyAlt = io.KeysDown[GLFW_KEY_LEFT_ALT] || io.KeysDown[GLFW_KEY_RIGHT_ALT];
		io.KeySuper = io.KeysDown[GLFW_KEY_LEFT_SUPER] || io.KeysDown[GLFW_KEY_RIGHT_SUPER];
	}
}

static void glfw_mouseButtonCallback(GLFWwindow *window, int button, int action, int /*mods*/)
{
	if (button != 0)
		return;
	if (glfwGetInputMode(s_window, GLFW_CURSOR) == GLFW_CURSOR_DISABLED && action == GLFW_RELEASE)
		glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
	if (s_options.gui) {
		ImGuiIO &io = ImGui::GetIO();
		io.MouseDown[button] = action == GLFW_PRESS;
		if (action == GLFW_PRESS && !io.WantCaptureMouse)
			glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
	} else if (action == GLFW_PRESS)
		glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
}

static void glfw_scrollCallback(GLFWwindow * /*window*/, double /*xoffset*/, double yoffset)
{
	ImGuiIO &io = ImGui::GetIO();
	if (s_options.gui && io.WantCaptureMouse)
		io.MouseWheel += (float)yoffset;
	else if (s_camera.mode == CameraMode::Orbit) 
		s_camera.orbit.zoom((float)-yoffset);
}

struct ProgramSource
{
	const uint8_t *vertexData;
	uint32_t vertexSize;
	const uint8_t *fragmentData;
	uint32_t fragmentSize;
};

struct ProgramSourceBundle
{
#if BX_PLATFORM_WINDOWS
	ProgramSource d3d11;
#endif
	ProgramSource gl;
};

#if BX_PLATFORM_WINDOWS
#define PROGRAM_SOURCE_BUNDLE(vname, fname) { \
	{ vname##_vertex_d3d11, sizeof(vname##_vertex_d3d11), fname##_fragment_d3d11, sizeof(fname##_fragment_d3d11) }, \
	{ vname##_vertex_gl, sizeof(vname##_vertex_gl), fname##_fragment_gl, sizeof(fname##_fragment_gl) }}
#else
#define PROGRAM_SOURCE_BUNDLE(vname, fname) { \
	{ vname##_vertex_gl, sizeof(vname##_vertex_gl), fname##_fragment_gl, sizeof(fname##_fragment_gl) }}
#endif

#define LOAD_PROGRAM(name) loadProgram(BX_STRINGIZE(name), PROGRAM_SOURCE_BUNDLE(name, name))

static bgfx::ProgramHandle loadProgram(const char *name, ProgramSourceBundle sourceBundle)
{
	ProgramSource source;
	if (bgfx::getRendererType() == bgfx::RendererType::OpenGL)
		source = sourceBundle.gl;
#if BX_PLATFORM_WINDOWS
	else if (bgfx::getRendererType() == bgfx::RendererType::Direct3D11)
		source = sourceBundle.d3d11;
#endif
	else {
		fprintf(stderr, "Unsupported renderer type.");
		return BGFX_INVALID_HANDLE;
	}
	bgfx::ShaderHandle vs = bgfx::createShader(bgfx::makeRef(source.vertexData, source.vertexSize));
	if (!bgfx::isValid(vs)) {
		fprintf(stderr, "Creating vertex shader '%s' failed.", name);
		return BGFX_INVALID_HANDLE;
	}
	bgfx::ShaderHandle fs = bgfx::createShader(bgfx::makeRef(source.fragmentData, source.fragmentSize));
	if (!bgfx::isValid(fs)) {
		fprintf(stderr, "Creating fragment shader '%s' failed.", name);
		bgfx::destroy(vs);
		return BGFX_INVALID_HANDLE;
	}
	bgfx::ProgramHandle program = bgfx::createProgram(vs, fs, true);
	if (!bgfx::isValid(program)) {
		bgfx::destroy(vs);
		bgfx::destroy(fs);
		fprintf(stderr, "Creating shader program '%s' failed.", name);
		return BGFX_INVALID_HANDLE;
	}
	return program;
}

static void guiInit()
{
	bgfx::setViewMode(kGuiView, bgfx::ViewMode::Sequential);
	bgfx::setViewRect(kGuiView, 0, 0, bgfx::BackbufferRatio::Equal);
	ImGui::CreateContext();
	int w, h;
	glfwGetWindowSize(s_window, &w, &h);
	if (w == 0 || h == 0) {
		w = WINDOW_DEFAULT_WIDTH;
		h = WINDOW_DEFAULT_HEIGHT;
	}
	ImGuiIO &io = ImGui::GetIO();
	io.DisplaySize.x = (float)w;
	io.DisplaySize.y = (float)h;
	io.IniFilename = nullptr;
	io.KeyMap[ImGuiKey_Tab] = GLFW_KEY_TAB;
	io.KeyMap[ImGuiKey_LeftArrow] = GLFW_KEY_LEFT;
	io.KeyMap[ImGuiKey_RightArrow] = GLFW_KEY_RIGHT;
	io.KeyMap[ImGuiKey_UpArrow] = GLFW_KEY_UP;
	io.KeyMap[ImGuiKey_DownArrow] = GLFW_KEY_DOWN;
	io.KeyMap[ImGuiKey_PageUp] = GLFW_KEY_PAGE_UP;
	io.KeyMap[ImGuiKey_PageDown] = GLFW_KEY_PAGE_DOWN;
	io.KeyMap[ImGuiKey_Home] = GLFW_KEY_HOME;
	io.KeyMap[ImGuiKey_End] = GLFW_KEY_END;
	io.KeyMap[ImGuiKey_Delete] = GLFW_KEY_DELETE;
	io.KeyMap[ImGuiKey_Backspace] = GLFW_KEY_BACKSPACE;
	io.KeyMap[ImGuiKey_Enter] = GLFW_KEY_ENTER;
	io.KeyMap[ImGuiKey_Escape] = GLFW_KEY_ESCAPE;
	io.KeyMap[ImGuiKey_A] = GLFW_KEY_A;
	io.KeyMap[ImGuiKey_C] = GLFW_KEY_C;
	io.KeyMap[ImGuiKey_V] = GLFW_KEY_V;
	io.KeyMap[ImGuiKey_X] = GLFW_KEY_X;
	io.KeyMap[ImGuiKey_Y] = GLFW_KEY_Y;
	io.KeyMap[ImGuiKey_Z] = GLFW_KEY_Z;
	// font
	int fontWidth, fontHeight;
	uint8_t *fontData;
	io.Fonts->GetTexDataAsRGBA32(&fontData, &fontWidth, &fontHeight);
	s_gui.font = bgfx::createTexture2D((uint16_t)fontWidth, (uint16_t)fontHeight, false, 0, bgfx::TextureFormat::RGBA8, 0, bgfx::makeRef(fontData, fontWidth * fontHeight * 4));
	io.Fonts->TexID = (ImTextureID)(intptr_t)s_gui.font.idx;
	// ImDrawVert vertex decl
	s_gui.vertexDecl
	.begin()
	.add(bgfx::Attrib::Position, 2, bgfx::AttribType::Float)
	.add(bgfx::Attrib::TexCoord0, 2, bgfx::AttribType::Float)
	.add(bgfx::Attrib::Color0, 4, bgfx::AttribType::Uint8, true)
	.end();
	// shader program
	s_gui.u_texture = bgfx::createUniform("u_texture", bgfx::UniformType::Sampler);
	s_gui.program = LOAD_PROGRAM(Gui);
}

static void guiShutdown()
{
	ImGui::DestroyContext();
	bgfx::destroy(s_gui.font);
	bgfx::destroy(s_gui.u_texture);
	bgfx::destroy(s_gui.program);
}

static void guiResize(int width, int height)
{
	ImGuiIO &io = ImGui::GetIO();
	io.DisplaySize.x = (float)width;
	io.DisplaySize.y = (float)height;
	bgfx::setViewRect(kGuiView, 0, 0, bgfx::BackbufferRatio::Equal);
}

static void guiRunFrame(float deltaTime)
{
	ImGuiIO &io = ImGui::GetIO();
	io.DeltaTime = deltaTime;
	if (glfwGetInputMode(s_window, GLFW_CURSOR) == GLFW_CURSOR_NORMAL) {
		double x, y;
		glfwGetCursorPos(s_window, &x, &y);
		io.MousePos.x = (float)x;
		io.MousePos.y = (float)y;
	}
	else
		io.MousePos.x = io.MousePos.y = -1.0f;
}

static void guiRender()
{
	ImGui::Render();
	ImDrawData *drawData = ImGui::GetDrawData();
	if (drawData->CmdListsCount == 0)
		return;
	ImGuiIO &io = ImGui::GetIO();
	float projection[16];
	bx::mtxOrtho(projection, 0, io.DisplaySize.x, io.DisplaySize.y, 0, 0, 1, 0, bgfx::getCaps()->homogeneousDepth);
	bgfx::setViewTransform(kGuiView, nullptr, projection);
	for (int n = 0; n < drawData->CmdListsCount; n++) {
		const ImDrawList* cmd_list = drawData->CmdLists[n];
		bgfx::TransientVertexBuffer tvb;
		bgfx::TransientIndexBuffer tib;
		if (!bgfx::allocTransientBuffers(&tvb, s_gui.vertexDecl, cmd_list->VtxBuffer.Size, &tib, cmd_list->IdxBuffer.Size))
			return;
		assert(sizeof(ImDrawVert) == s_gui.vertexDecl.getStride());
		memcpy(tvb.data, cmd_list->VtxBuffer.Data, tvb.size);
		assert(sizeof(ImDrawIdx) == sizeof(uint16_t));
		memcpy(tib.data, cmd_list->IdxBuffer.Data, tib.size);
		uint32_t firstIndex = 0;
		for (int cmd_i = 0; cmd_i < cmd_list->CmdBuffer.Size; cmd_i++) {
			const ImDrawCmd* pcmd = &cmd_list->CmdBuffer[cmd_i];
			if (pcmd->UserCallback)
				pcmd->UserCallback(cmd_list, pcmd);
			else {
				bgfx::setScissor((uint16_t)pcmd->ClipRect.x, (uint16_t)pcmd->ClipRect.y, uint16_t(pcmd->ClipRect.z - pcmd->ClipRect.x), uint16_t(pcmd->ClipRect.w - pcmd->ClipRect.y));
				bgfx::setState(BGFX_STATE_WRITE_RGB | BGFX_STATE_WRITE_A | BGFX_STATE_BLEND_ALPHA);
				bgfx::TextureHandle texture;
				texture.idx = (uint16_t)(intptr_t)pcmd->TextureId;
				bgfx::setTexture(0, s_gui.u_texture, texture);
				bgfx::setIndexBuffer(&tib, firstIndex, pcmd->ElemCount);
				bgfx::setVertexBuffer(0, &tvb);
				bgfx::submit(kGuiView, s_gui.program);
			}
			firstIndex += pcmd->ElemCount;
		}
	}
}

static void modelInit()
{
	s_model.u_color = bgfx::createUniform("u_color", bgfx::UniformType::Vec4);
	s_model.u_textureSize_cellSize = bgfx::createUniform("u_textureSize_cellSize", bgfx::UniformType::Vec4);
	s_model.u_diffuse = bgfx::createUniform("u_diffuse", bgfx::UniformType::Vec4);
	s_model.u_emission = bgfx::createUniform("u_emission", bgfx::UniformType::Vec4);
	s_model.u_lightDir = bgfx::createUniform("u_lightDir", bgfx::UniformType::Vec4);
	s_model.colorProgram = LOAD_PROGRAM(Color);
	s_model.checkerboardProgram = LOAD_PROGRAM(Checkerboard);
	s_model.flatProgram = LOAD_PROGRAM(Flat);
	ModelVertex::init();
	bgfx::setViewClear(kModelView, BGFX_CLEAR_COLOR | BGFX_CLEAR_DEPTH, 0x444444ff);
	bgfx::setViewRect(kModelView, 0, 0, bgfx::BackbufferRatio::Equal);
}

static void modelDestroy();

static void modelShutdown()
{
	modelDestroy();
	bgfx::destroy(s_model.u_color);
	bgfx::destroy(s_model.u_textureSize_cellSize);
	bgfx::destroy(s_model.u_diffuse);
	bgfx::destroy(s_model.u_emission);
	bgfx::destroy(s_model.u_lightDir);
	bgfx::destroy(s_model.colorProgram);
	bgfx::destroy(s_model.checkerboardProgram);
	bgfx::destroy(s_model.flatProgram);
}

struct ModelLoadThreadArgs
{
	char filename[256];
};

static void modelLoadThread(ModelLoadThreadArgs args)
{
	objz_setIndexFormat(OBJZ_INDEX_FORMAT_U32);
	objz_setVertexFormat(sizeof(ModelVertex), offsetof(ModelVertex, pos), offsetof(ModelVertex, texcoord), offsetof(ModelVertex, normal));
	objzModel *model = objz_load(args.filename);
	if (!model) {
		fprintf(stderr, "%s\n", objz_getError());
		s_errorMessage.set(objz_getError());
		s_model.data = nullptr;
		s_model.status.set(ModelStatus::NotLoaded);
		return;
	} else if (objz_getError())
		printf("%s\n", objz_getError());
	s_model.data = model;
	for (uint32_t i = 0; i < model->numVertices; i++) {
		auto v = &((ModelVertex *)model->vertices)[i];
		v->texcoord[1] = 1.0f - v->texcoord[1];
	}
	s_model.status.set(ModelStatus::Finalizing);
}

static void modelFinalize()
{
	if (s_model.thread) {
		s_model.thread->join();
		delete s_model.thread;
		s_model.thread = nullptr;
	}
	s_model.aabb = AABB();
	s_model.centroid = bx::Vec3(0.0f, 0.0f, 0.0f);
	for (uint32_t i = 0; i < s_model.data->numVertices; i++) {
		const bx::Vec3 &pos = ((const ModelVertex *)s_model.data->vertices)[i].pos;
		s_model.aabb.addPoint(pos);
		s_model.centroid = bx::add(s_model.centroid, pos);
	}
	s_model.centroid = bx::mul(s_model.centroid, 1.0f / s_model.data->numVertices);
	s_model.vb = bgfx::createVertexBuffer(bgfx::makeRef(s_model.data->vertices, s_model.data->numVertices * sizeof(ModelVertex)), ModelVertex::decl);
	s_model.ib = bgfx::createIndexBuffer(bgfx::makeRef(s_model.data->indices, s_model.data->numIndices * sizeof(uint32_t)), BGFX_BUFFER_INDEX32);
	const uint32_t numWireframeIndices = bgfx::topologyConvert(bgfx::TopologyConvert::TriListToLineList, nullptr, 0, s_model.data->indices, s_model.data->numIndices, true);
	const bgfx::Memory *wireframeIndices = bgfx::alloc(numWireframeIndices * sizeof(uint32_t));
	bgfx::topologyConvert(bgfx::TopologyConvert::TriListToLineList, wireframeIndices->data, wireframeIndices->size, s_model.data->indices, s_model.data->numIndices, true);
	s_model.wireframeIb = bgfx::createIndexBuffer(wireframeIndices, BGFX_BUFFER_INDEX32);
	s_camera.firstPerson = FirstPersonCamera();
	s_camera.orbit = OrbitCamera();
	s_options.shadeMode = ShadeMode::Flat;
	s_options.wireframeMode = WireframeMode::Triangles;
	s_model.status.set(ModelStatus::Ready);
}

static void modelOpenDialog()
{
	if (s_model.status.get() == ModelStatus::Loading || s_model.status.get() == ModelStatus::Finalizing)
		return;
	if (!(s_atlas.status.get() == AtlasStatus::NotGenerated || s_atlas.status.get() == AtlasStatus::Ready))
		return;
	nfdchar_t *filename = nullptr;
	nfdresult_t result = NFD_OpenDialog("obj", nullptr, &filename);
	if (result != NFD_OKAY)
		return;
	modelDestroy();
	s_model.status.set(ModelStatus::Loading);
	char windowTitle[256];
	snprintf(windowTitle, sizeof(windowTitle), "%s - %s\n", WINDOW_TITLE, filename);
	glfwSetWindowTitle(s_window, windowTitle);
	printf("Loading '%s'\n", filename);
	ModelLoadThreadArgs args;
	STRNCPY(args.filename, sizeof(args.filename), filename);
	s_model.thread = new std::thread(modelLoadThread, args);
	free(filename);
}

static void atlasDestroy();

static void modelDestroy()
{
	atlasDestroy();
	if (s_model.thread) {
		s_model.thread->join();
		delete s_model.thread;
		s_model.thread = nullptr;
	}
	if (s_model.data) {
		objz_destroy(s_model.data);
		s_model.data = nullptr;
	}
	if (bgfx::isValid(s_model.vb)) {
		bgfx::destroy(s_model.vb);
		bgfx::destroy(s_model.ib);
		bgfx::destroy(s_model.wireframeIb);
		s_model.vb = BGFX_INVALID_HANDLE;
		s_model.ib = BGFX_INVALID_HANDLE;
		s_model.wireframeIb = BGFX_INVALID_HANDLE;
	}
	glfwSetWindowTitle(s_window, WINDOW_TITLE);
	s_model.status.set(ModelStatus::NotLoaded);
}

static void modelSetMaterialUniforms(const objzMaterial *mat)
{
	if (!mat) {
		const float diffuse[] = { 0.75f, 0.75f, 0.75f, 1.0f };
		const float emission[] = { 0.0f, 0.0f, 0.0f, 0.0f };
		bgfx::setUniform(s_model.u_diffuse, diffuse);
		bgfx::setUniform(s_model.u_emission, emission);
	} else {
		const float diffuse[] = { mat->diffuse[0], mat->diffuse[1], mat->diffuse[2], 1.0f };
		const float emission[] = { mat->emission[0], mat->emission[1], mat->emission[2], 1.0f };
		bgfx::setUniform(s_model.u_diffuse, diffuse);
		bgfx::setUniform(s_model.u_emission, emission);
	}
}

static void modelRender(const float *view, const float *projection)
{
	if (s_model.status.get() != ModelStatus::Ready)
		return;
	float modelMatrix[16];
	bx::mtxScale(modelMatrix, s_model.scale);
	bgfx::setViewTransform(kModelView, view, projection);
	if (s_options.shadeMode == ShadeMode::Flat) {
		const float lightDir[] = { view[2], view[6], view[10], 0 };
		for (uint32_t i = 0; i < s_model.data->numMeshes; i++) {
			const objzMesh &mesh = s_model.data->meshes[i];
			const objzMaterial *mat = mesh.materialIndex == -1 ? nullptr : &s_model.data->materials[mesh.materialIndex];
			if (s_atlas.status.get() == AtlasStatus::Ready) {
				bgfx::setIndexBuffer(s_atlas.ib, mesh.firstIndex, mesh.numIndices);
				bgfx::setVertexBuffer(0, s_atlas.vb);
			} else {
				bgfx::setIndexBuffer(s_model.ib, mesh.firstIndex, mesh.numIndices);
				bgfx::setVertexBuffer(0, s_model.vb);
			}
			bgfx::setState(BGFX_STATE_DEFAULT);
			bgfx::setTransform(modelMatrix);
			bgfx::setUniform(s_model.u_lightDir, lightDir);
			modelSetMaterialUniforms(mat);
			bgfx::submit(kModelView, s_model.flatProgram);
		}
	} else if (s_options.shadeMode == ShadeMode::Charts && s_atlas.status.get() == AtlasStatus::Ready) {
		srand(s_atlas.chartColorSeed);
		uint32_t firstIndex = 0;
		for (uint32_t i = 0; i < s_atlas.data->meshCount; i++) {
			const xatlas::Mesh &mesh = s_atlas.data->meshes[i];
			for (uint32_t j = 0; j < mesh.chartCount; j++) {
				const xatlas::Chart &chart = mesh.chartArray[j];
				uint8_t bcolor[3];
				randomRGB(bcolor);
				float color[4];
				color[0] = bcolor[0] / 255.0f;
				color[1] = bcolor[1] / 255.0f;
				color[2] = bcolor[2] / 255.0f;
				color[3] = 1.0f;
				bgfx::setUniform(s_model.u_color, color);
				float textureSize_cellSize[4];
				textureSize_cellSize[0] = (float)s_atlas.data->width;
				textureSize_cellSize[1] = (float)s_atlas.data->height;
				textureSize_cellSize[2] = (float)s_atlas.chartCellSize;
				textureSize_cellSize[3] = (float)s_atlas.chartCellSize;
				bgfx::setUniform(s_model.u_textureSize_cellSize, textureSize_cellSize);
				bgfx::setState(BGFX_STATE_DEFAULT);
				bgfx::setTransform(modelMatrix);
				bgfx::setIndexBuffer(s_atlas.chartIb, firstIndex, chart.indexCount);
				bgfx::setVertexBuffer(0, s_atlas.vb);
				bgfx::submit(kModelView, s_model.checkerboardProgram);
				firstIndex += chart.indexCount;
			}
		}
	}
	if (s_options.wireframe) {
		const float color[] = { 1.0f, 1.0f, 1.0f, 0.5f };
		bgfx::setUniform(s_model.u_color, color);
		bgfx::setState(BGFX_STATE_WRITE_RGB | BGFX_STATE_WRITE_A | BGFX_STATE_PT_LINES | BGFX_STATE_BLEND_ALPHA);
		bgfx::setTransform(modelMatrix);
		if (s_options.wireframeMode == WireframeMode::Triangles) {
			bgfx::setIndexBuffer(s_model.wireframeIb);
			bgfx::setVertexBuffer(0, s_model.vb);
		} else if (s_atlas.status.get() == AtlasStatus::Ready && s_options.wireframeMode == WireframeMode::Charts) {
			bgfx::setVertexBuffer(0, s_atlas.chartBoundaryVb);
		}
		bgfx::submit(kModelView, s_model.colorProgram);
	}
}

static void atlasInit()
{
	s_atlas.wireVertexDecl
		.begin()
		.add(bgfx::Attrib::Position, 3, bgfx::AttribType::Float)
		.end();
	assert(s_atlas.wireVertexDecl.getStride() == sizeof(bx::Vec3));
	bgfx::setPaletteColor(kPaletteBlack, 0x000000ff);
}

static void bakeClear();

static void atlasDestroy()
{
	bakeClear();
	if (s_atlas.thread) {
		s_atlas.thread->join();
		delete s_atlas.thread;
		s_atlas.thread = nullptr;
	}
	if (s_atlas.data) {
		xatlas::Destroy(s_atlas.data);
		s_atlas.data = nullptr;
	}
	for (uint32_t i = 0; i < (uint32_t)s_atlas.chartsTextures.size(); i++) {
		if (bgfx::isValid(s_atlas.chartsTextures[i])) {
			bgfx::destroy(s_atlas.chartsTextures[i]);
			s_atlas.chartsTextures[i] = BGFX_INVALID_HANDLE;
		}
	}
	if (bgfx::isValid(s_atlas.vb)) {
		bgfx::destroy(s_atlas.vb);
		bgfx::destroy(s_atlas.ib);
		bgfx::destroy(s_atlas.chartIb);
		s_atlas.vb = BGFX_INVALID_HANDLE;
		s_atlas.ib = BGFX_INVALID_HANDLE;
		s_atlas.chartIb = BGFX_INVALID_HANDLE;
	}
	if (bgfx::isValid(s_atlas.chartBoundaryVb)) {
		bgfx::destroy(s_atlas.chartBoundaryVb);
		s_atlas.chartBoundaryVb = BGFX_INVALID_HANDLE;
	}
	s_atlas.status.set(AtlasStatus::NotGenerated);
}

static void atlasProgressCallback(xatlas::ProgressCategory::Enum category, int progress, void * /*userData*/)
{
	s_atlas.status.setProgress(category, progress);
}

static void atlasSetPixel(uint8_t *dest, int destWidth, int x, int y, const uint8_t *color, bool checkerboard)
{
	float scale = 1.0f;
	if (checkerboard)
		scale = (x / s_atlas.chartCellSize % 2) != (y / s_atlas.chartCellSize % 2) ? 0.75f : 1.0f;
	uint8_t *pixel = &dest[x * 3 + y * (destWidth * 3)];
	pixel[0] = uint8_t(color[0] * scale);
	pixel[1] = uint8_t(color[1] * scale);
	pixel[2] = uint8_t(color[2] * scale);
}

// https://github.com/miloyip/line/blob/master/line_bresenham.c
// License: public domain.
static void atlasRasterizeLine(uint8_t *dest, int destWidth, const int *p1, const int *p2, const uint8_t *color)
{
	const int dx = abs(p2[0] - p1[0]), sx = p1[0] < p2[0] ? 1 : -1;
	const int dy = abs(p2[1] - p1[1]), sy = p1[1] < p2[1] ? 1 : -1;
	int err = (dx > dy ? dx : -dy) / 2;
	int current[2];
	current[0] = p1[0];
	current[1] = p1[1];
	while (atlasSetPixel(dest, destWidth, current[0], current[1], color, false), current[0] != p2[0] || current[1] != p2[1]) {
		const int e2 = err;
		if (e2 > -dx) { err -= dy; current[0] += sx; }
		if (e2 < dy) { err += dx; current[1] += sy; }
	}
}

/*
https://github.com/ssloy/tinyrenderer/wiki/Lesson-2:-Triangle-rasterization-and-back-face-culling
Copyright Dmitry V. Sokolov

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
static void atlasRasterizeTriangle(uint8_t *dest, int destWidth, const int *t0, const int *t1, const int *t2, const uint8_t *color)
{
	if (t0[1] > t1[1]) std::swap(t0, t1);
	if (t0[1] > t2[1]) std::swap(t0, t2);
	if (t1[1] > t2[1]) std::swap(t1, t2);
	int total_height = t2[1] - t0[1];
	for (int i = 0; i < total_height; i++) {
		bool second_half = i > t1[1] - t0[1] || t1[1] == t0[1];
		int segment_height = second_half ? t2[1] - t1[1] : t1[1] - t0[1];
		float alpha = (float)i / total_height;
		float beta = (float)(i - (second_half ? t1[1] - t0[1] : 0)) / segment_height;
		int A[2], B[2];
		for (int j = 0; j < 2; j++) {
			A[j] = int(t0[j] + (t2[j] - t0[j]) * alpha);
			B[j] = int(second_half ? t1[j] + (t2[j] - t1[j]) * beta : t0[j] + (t1[j] - t0[j]) * beta);
		}
		if (A[0] > B[0])
			std::swap(A, B);
		for (int j = A[0]; j <= B[0]; j++)
			atlasSetPixel(dest, destWidth, j, t0[1] + i, color, true);
	}
}

static uint32_t sdbmHash(const void *data_in, uint32_t size, uint32_t h = 5381)
{
	const uint8_t *data = (const uint8_t *) data_in;
	uint32_t i = 0;
	while (i < size)
		h = (h << 16) + (h << 6) - h + (uint32_t ) data[i++];
	return h;
}

struct EdgeKey
{
	bx::Vec3 p0, p1;
};

struct EdgeKeyHash
{
	uint32_t operator()(const EdgeKey &key) const
	{
		int32_t data[6];
		data[0] = (int32_t)(key.p0.x * 100.0f);
		data[1] = (int32_t)(key.p0.y * 100.0f);
		data[2] = (int32_t)(key.p0.z * 100.0f);
		data[3] = (int32_t)(key.p1.x * 100.0f);
		data[4] = (int32_t)(key.p1.y * 100.0f);
		data[5] = (int32_t)(key.p1.z * 100.0f);
		return sdbmHash(data, sizeof(data));
	}
};

struct EdgeKeyEqual
{
	bool floatEqual(float f1, float f2) const
	{
		const float epsilon = 0.0001f;
		return bx::abs(f1 - f2) <= epsilon;
	}

	bool vec3Equal(const bx::Vec3 &v0, const bx::Vec3 &v1) const
	{
		return floatEqual(v0.x, v1.x) && floatEqual(v0.y, v1.y) && floatEqual(v0.z, v1.z);
	}

	bool operator()(const EdgeKey &k0, const EdgeKey &k1) const
	{
		return vec3Equal(k0.p0, k1.p0) && vec3Equal(k0.p1, k1.p1);
	}
};

#if USE_LIBIGL
static void atlasParameterizationCallback(const float *positions, float *texcoords, uint32_t vertexCount, const uint32_t *indices, uint32_t indexCount, bool isPlanar)
{
	if (isPlanar)
		return;
	Eigen::MatrixXd V(vertexCount, 3);
	for (uint32_t i = 0; i < vertexCount; i++) {
		V(i, 0) = positions[i * 3 + 0];
		V(i, 1) = positions[i * 3 + 1];
		V(i, 2) = positions[i * 3 + 2];
	}
	Eigen::MatrixXi F(indexCount / 3, 3);	
	for (uint32_t i = 0; i < indexCount / 3; i++) {
		F(i, 0) = indices[i * 3 + 0];
		F(i, 1) = indices[i * 3 + 1];
		F(i, 2) = indices[i * 3 + 2];
	}
	// Fix two points on the boundary
	Eigen::VectorXi bnd;
	igl::boundary_loop(F, bnd);
	Eigen::MatrixXd V_uv;
	if (s_atlas.paramMethod == ParamMethod::libigl_Harmonic) {
		Eigen::MatrixXd bnd_uv;
		igl::map_vertices_to_circle(V,bnd,bnd_uv);
		igl::harmonic(V,F,bnd,bnd_uv,1,V_uv);
	} else if (s_atlas.paramMethod == ParamMethod::libigl_LSCM) {
		Eigen::VectorXi b(2, 1);
		b(0) = bnd(0);
		b(1) = bnd((int)round(bnd.size()/2));
		// Use existing orthographic parameterization to scale the result.
		Eigen::MatrixXd bc(2,2);
		bc(0, 0) = texcoords[b(0) * 2 + 0];
		bc(0, 1) = texcoords[b(0) * 2 + 1];
		bc(1, 0) = texcoords[b(1) * 2 + 0];
		bc(1, 1) = texcoords[b(1) * 2 + 1];
		igl::lscm(V, F, b, bc, V_uv);
	} else if (s_atlas.paramMethod == ParamMethod::libigl_ARAP) {
		// Compute the initial solution for ARAP (harmonic parametrization)
		Eigen::MatrixXd initial_guess;
		Eigen::MatrixXd bnd_uv;
		igl::map_vertices_to_circle(V,bnd,bnd_uv);
		igl::harmonic(V,F,bnd,bnd_uv,1,initial_guess);
		// Add dynamic regularization to avoid to specify boundary conditions
		igl::ARAPData arap_data;
		arap_data.with_dynamics = true;
		Eigen::VectorXi b  = Eigen::VectorXi::Zero(0);
		Eigen::MatrixXd bc = Eigen::MatrixXd::Zero(0,0);
		// Initialize ARAP
		arap_data.max_iter = 100;
		// 2 means that we're going to *solve* in 2d
		arap_precomputation(V,F,2,b,arap_data);
		// Solve arap using the harmonic map as initial guess
		V_uv = initial_guess;
		arap_solve(bc,arap_data,V_uv);
	}
	for (uint32_t i = 0; i < vertexCount; i++) {
		texcoords[i * 2 + 0] = (float)V_uv(i, 0);
		texcoords[i * 2 + 1] = (float)V_uv(i, 1);
	}
}
#endif

static void atlasGenerateThread()
{
	int progress = 0;
	const bool firstRun = !s_atlas.data;
	if (firstRun) {
		// Create xatlas context and compute charts on first run only.
		s_atlas.data = xatlas::Create();
		for (uint32_t i = 0; i < s_model.data->numObjects; i++) {
			const objzObject &object = s_model.data->objects[i];
			auto v = &((const ModelVertex *)s_model.data->vertices)[object.firstVertex];
			xatlas::MeshDecl meshDecl;
			meshDecl.vertexCount = object.numVertices;
			meshDecl.vertexPositionData = &v->pos;
			meshDecl.vertexPositionStride = sizeof(ModelVertex);
			meshDecl.vertexNormalData = &v->normal;
			meshDecl.vertexNormalStride = sizeof(ModelVertex);
			meshDecl.vertexUvData = &v->texcoord;
			meshDecl.vertexUvStride = sizeof(ModelVertex);
			meshDecl.indexCount = object.numIndices;
			meshDecl.indexData = &((uint32_t *)s_model.data->indices)[object.firstIndex];
			meshDecl.indexFormat = xatlas::IndexFormat::UInt32;
			meshDecl.indexOffset = -(int32_t)object.firstVertex;
			xatlas::AddMeshError::Enum error = xatlas::AddMesh(s_atlas.data, meshDecl);
			if (error != xatlas::AddMeshError::Success) {
				fprintf(stderr, "Error adding mesh: %s\n", xatlas::StringForEnum(error));
				s_errorMessage.set("Error adding mesh: %s", xatlas::StringForEnum(error));
				xatlas::Destroy(s_atlas.data);
				s_atlas.data = nullptr;
				s_atlas.status.set(AtlasStatus::NotGenerated);
				return;
			}
			const int newProgress = int((i + 1) / (float)s_model.data->numObjects * 100.0f);
			if (newProgress != progress) {
				progress = newProgress;
				s_atlas.status.setProgress((xatlas::ProgressCategory::Enum)-1, progress);
			}
		}
		s_atlas.status.set(AtlasStatus::Generating);
		xatlas::ComputeCharts(s_atlas.data, s_atlas.chartOptions, atlasProgressCallback);
	} else
		s_atlas.status.set(AtlasStatus::Generating);
	if (firstRun || s_atlas.paramMethodChanged) {
		s_atlas.paramMethodChanged = false;
		xatlas::ParameterizeFunc paramFunc = nullptr;
#if USE_LIBIGL
		if (s_atlas.paramMethod != ParamMethod::LSCM)
			paramFunc = atlasParameterizationCallback;
#endif
		xatlas::ParameterizeCharts(s_atlas.data, paramFunc, atlasProgressCallback);
	}
	xatlas::PackCharts(s_atlas.data, s_atlas.packOptions, atlasProgressCallback);
	// Find chart boundary edges.
	uint32_t numEdges = 0;
	for (uint32_t i = 0; i < s_atlas.data->meshCount; i++) {
		const xatlas::Mesh &outputMesh = s_atlas.data->meshes[i];
		numEdges += outputMesh.indexCount;
	}
	std::vector<bool> boundaryEdges;
	boundaryEdges.resize(numEdges);
	std::unordered_map<EdgeKey, uint32_t, EdgeKeyHash, EdgeKeyEqual> edgeMap;
	numEdges = 0;
	for (uint32_t i = 0; i < s_atlas.data->meshCount; i++) {
		const xatlas::Mesh &mesh = s_atlas.data->meshes[i];
		const objzObject &object = s_model.data->objects[i];
		const ModelVertex *vertices = &((const ModelVertex *)s_model.data->vertices)[object.firstVertex];
		for (uint32_t j = 0; j < mesh.chartCount; j++) {
			const xatlas::Chart &chart = mesh.chartArray[j];
			// Hash edges for finding boundaries.
			edgeMap.clear();
			for (uint32_t k = 0; k < chart.indexCount; k += 3) {
				for (uint32_t l = 0; l < 3; l++) {
					EdgeKey key;
					key.p0 = vertices[mesh.vertexArray[chart.indexArray[k + l]].xref].pos;
					key.p1 = vertices[mesh.vertexArray[chart.indexArray[k + (l + 1) % 3]].xref].pos;
					edgeMap[key] = 0; // Don't care.
				}
			}
			for (uint32_t k = 0; k < chart.indexCount; k += 3) {
				for (uint32_t l = 0; l < 3; l++) {
					EdgeKey key;
					key.p0 = vertices[mesh.vertexArray[chart.indexArray[k + (l + 1) % 3]].xref].pos;
					key.p1 = vertices[mesh.vertexArray[chart.indexArray[k + l]].xref].pos;
					boundaryEdges[numEdges] = edgeMap.count(key) == 0;
					numEdges++;
				}
			}
		}
	}
	// Copy charts for rendering.
	s_atlas.vertices.clear();
	s_atlas.indices.clear();
	s_atlas.chartIndices.clear();
	s_atlas.chartBoundaryVertices.clear();
	uint32_t numIndices = 0, numVertices = 0;
	for (uint32_t i = 0; i < s_atlas.data->meshCount; i++) {
		const xatlas::Mesh &outputMesh = s_atlas.data->meshes[i];
		numIndices += outputMesh.indexCount;
		numVertices += outputMesh.vertexCount;
	}
	s_atlas.vertices.resize(numVertices);
	s_atlas.indices.resize(numIndices);
	s_atlas.chartIndices.resize(numIndices);
	uint32_t firstIndex = 0;
	uint32_t firstChartIndex = 0;
	uint32_t firstVertex = 0;
	numEdges = 0;
	for (uint32_t i = 0; i < s_atlas.data->meshCount; i++) {
		const xatlas::Mesh &mesh = s_atlas.data->meshes[i];
		const objzObject &object = s_model.data->objects[i];
		const ModelVertex *oldVertices = &((const ModelVertex *)s_model.data->vertices)[object.firstVertex];
		for (uint32_t j = 0; j < mesh.indexCount; j++)
			s_atlas.indices[firstIndex + j] = firstVertex + mesh.indexArray[j];
		firstIndex += mesh.indexCount;
		for (uint32_t j = 0; j < mesh.chartCount; j++) {
			const xatlas::Chart &chart = mesh.chartArray[j];
			for (uint32_t k = 0; k < chart.indexCount; k++)
				s_atlas.chartIndices[firstChartIndex + k] = firstVertex + chart.indexArray[k];
			firstChartIndex += chart.indexCount;
			for (uint32_t k = 0; k < chart.indexCount; k += 3) {
				for (int l = 0; l < 3; l++) {
					if (boundaryEdges[numEdges]) {
						const xatlas::Vertex &v0 = mesh.vertexArray[chart.indexArray[k + l]];
						const xatlas::Vertex &v1 = mesh.vertexArray[chart.indexArray[k + (l + 1) % 3]];
						s_atlas.chartBoundaryVertices.push_back(oldVertices[v0.xref].pos);
						s_atlas.chartBoundaryVertices.push_back(oldVertices[v1.xref].pos);
					}
					numEdges++;
				}
			}
		}
		for (uint32_t j = 0; j < mesh.vertexCount; j++) {
			const xatlas::Vertex &outputVertex = mesh.vertexArray[j];
			const ModelVertex &oldVertex = oldVertices[outputVertex.xref];
			ModelVertex &v = s_atlas.vertices[firstVertex + j];
			v.pos = oldVertex.pos;
			v.normal = oldVertex.normal;
			v.texcoord[0] = oldVertex.texcoord[0];
			v.texcoord[1] = oldVertex.texcoord[1];
			v.texcoord[2] = outputVertex.uv[0] / (float)s_atlas.data->width;
			v.texcoord[3] = outputVertex.uv[1] / (float)s_atlas.data->height;
		}
		firstVertex += mesh.vertexCount;
	}
	/*for (int i = 0; i < (int)s_atlas.atlas->meshCount; i++)
	{
		const xatlas::Mesh &outputMesh = s_atlas.atlas->meshes[i];
		const objzObject &object = s_model.data->objects[i];
		const int firstIndex = (int)indices.size();
		indices.resize(indices.size() + outputMesh.indexCount);
		for (int j = 0; j < (int)outputMesh.indexCount; j++)
			indices[firstIndex + j] = (uint32_t)vertices.size() + outputMesh.indexArray[j];
		const int firstVertex = (int)vertices.size();
		vertices.resize(vertices.size() + outputMesh.vertexCount);
		for (int j = 0; j < (int)outputMesh.vertexCount; j++)
		{
			const xatlas::Vertex &outputVertex = outputMesh.vertexArray[j];
			const ModelVertex &oldVertex = ((const ModelVertex *)s_model.data->vertices)[object.firstVertex + outputVertex.xref];
			ModelVertex &v = vertices[firstVertex + j];
			v.pos = oldVertex.pos;
			v.normal = oldVertex.normal;
			v.texcoord = vec4(oldVertex.texcoord.x, oldVertex.texcoord.y, outputVertex.uv[0] / (float)width, outputVertex.uv[1] / (float)height);
		}
	}*/
	// Rasterize charts to texture(s) for previewing UVs.
	s_atlas.chartsImages.resize(s_atlas.data->atlasCount);
	for (uint32_t i = 0; i < (uint32_t)s_atlas.chartsImages.size(); i++) {
		s_atlas.chartsImages[i].resize(s_atlas.data->width * s_atlas.data->height * 3);
		memset(s_atlas.chartsImages[i].data(), 0, s_atlas.chartsImages[i].size());
	}
	srand(s_atlas.chartColorSeed);
	numEdges = 0;
	for (uint32_t i = 0; i < s_atlas.data->meshCount; i++) {
		const xatlas::Mesh &mesh = s_atlas.data->meshes[i];
		for (uint32_t j = 0; j < mesh.chartCount; j++) {
			const xatlas::Chart &chart = mesh.chartArray[j];
			uint8_t color[3];
			randomRGB(color);
			for (uint32_t k = 0; k < chart.indexCount; k += 3) {
				int verts[3][2];
				for (int l = 0; l < 3; l++) {
					const xatlas::Vertex &v = mesh.vertexArray[chart.indexArray[k + l]];
					verts[l][0] = int(v.uv[0]);
					verts[l][1] = int(v.uv[1]);
				}
				uint8_t *imageData = s_atlas.chartsImages[chart.atlasIndex].data();
				atlasRasterizeTriangle(imageData, s_atlas.data->width, verts[0], verts[1], verts[2], color);
				for (int l = 0; l < 3; l++) {
					if (boundaryEdges[numEdges]) {
						const uint8_t white[] = { 255, 255, 255, 255 };
						atlasRasterizeLine(imageData, s_atlas.data->width, verts[l], verts[(l + 1) % 3], white);
					}
					numEdges++;
				}
			}
		}
	}
	s_atlas.status.set(AtlasStatus::Finalizing);
}

static void atlasGenerate()
{
	if (!(s_atlas.status.get() == AtlasStatus::NotGenerated || s_atlas.status.get() == AtlasStatus::Ready))
		return;
	bakeClear();
	xatlas::SetPrint(printf, s_atlas.verbose);
	s_atlas.status.set(AtlasStatus::AddingMeshes);
	s_atlas.thread = new std::thread(atlasGenerateThread);
}

static void atlasFinalize()
{
	if (s_atlas.thread) {
		s_atlas.thread->join();
		delete s_atlas.thread;
		s_atlas.thread = nullptr;
	}
	// Charts geometry.
	s_atlas.vb = bgfx::createVertexBuffer(bgfx::makeRef(s_atlas.vertices.data(), uint32_t(s_atlas.vertices.size() * sizeof(ModelVertex))), ModelVertex::decl);
	s_atlas.ib = bgfx::createIndexBuffer(bgfx::makeRef(s_atlas.indices.data(), uint32_t(s_atlas.indices.size() * sizeof(uint32_t))), BGFX_BUFFER_INDEX32);
	s_atlas.chartIb = bgfx::createIndexBuffer(bgfx::makeRef(s_atlas.chartIndices.data(), uint32_t(s_atlas.chartIndices.size() * sizeof(uint32_t))), BGFX_BUFFER_INDEX32);
	// Chart boundaries.
	s_atlas.chartBoundaryVb = bgfx::createVertexBuffer(bgfx::makeRef(s_atlas.chartBoundaryVertices.data(), uint32_t(s_atlas.chartBoundaryVertices.size() * sizeof(bx::Vec3))), s_atlas.wireVertexDecl);
	// Charts texture.
	if (s_atlas.chartsTextures.size() != s_atlas.data->atlasCount) {
		for (uint32_t i = 0; i < (uint32_t)s_atlas.chartsTextures.size(); i++)
			bgfx::destroy(s_atlas.chartsTextures[i]);
		s_atlas.chartsTextures.resize(s_atlas.data->atlasCount);
		for (uint32_t i = 0; i < (uint32_t)s_atlas.chartsTextures.size(); i++)
			s_atlas.chartsTextures[i] = BGFX_INVALID_HANDLE;
	}
	for (uint32_t i = 0; i < (uint32_t)s_atlas.chartsTextures.size(); i++) {
		bgfx::TextureHandle &tex = s_atlas.chartsTextures[i];
		const bgfx::Memory *mem = bgfx::makeRef(s_atlas.chartsImages[i].data(), (uint32_t)s_atlas.chartsImages[i].size());
		if (!bgfx::isValid(tex))
			tex = bgfx::createTexture2D((uint16_t)s_atlas.data->width, (uint16_t)s_atlas.data->height, false, 1, bgfx::TextureFormat::RGB8, BGFX_SAMPLER_UVW_BORDER | BGFX_SAMPLER_BORDER_COLOR(kPaletteBlack), mem);
		else
			bgfx::updateTexture2D(tex, 0, 0, 0, 0, (uint16_t)s_atlas.data->width, (uint16_t)s_atlas.data->height, mem);
	}
	s_atlas.currentTexture = 0;
	s_options.shadeMode = ShadeMode::Charts;
	s_options.wireframeMode = WireframeMode::Charts;
	s_atlas.status.set(AtlasStatus::Ready);
}

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

struct
{
	const uint16_t rbTextureSize = 512;
	const uint16_t rbDataTextureSize = 8192;
	bool enabled;
	bool initialized = false;
	bool executed = false;
	bool finished = false;
	int numDirections = 1000;
	int directionCount;
	uint32_t lightmapWidth, lightmapHeight;
	// programs
	bgfx::ProgramHandle atomicCounterClearProgram;
	bgfx::ProgramHandle rayBundleClearProgram;
	bgfx::ProgramHandle rayBundleWriteProgram;
	bgfx::ProgramHandle rayBundleIntegrateProgram;
	bgfx::ProgramHandle lightmapClearProgram;
	// uniforms
	bgfx::UniformHandle u_lightmapSize_dataSize;
	bgfx::UniformHandle u_rayNormal;
	bgfx::UniformHandle u_atomicCounterSampler;
	bgfx::UniformHandle u_rayBundleHeaderSampler;
	bgfx::UniformHandle u_rayBundleDataSampler;
	bgfx::UniformHandle u_lightmapSampler;
	// atomic counter
	bgfx::FrameBufferHandle atomicCounterFb;
	bgfx::TextureHandle atomicCounterTexture;
	// ray bundle write
	bgfx::FrameBufferHandle rayBundleFb;
	bgfx::TextureHandle rayBundleTarget;
	bgfx::TextureHandle rayBundleHeader;
	bgfx::TextureHandle rayBundleData;
	// lightmap clear
	bgfx::FrameBufferHandle lightmapClearFb;
	bgfx::TextureHandle lightmapClearTarget;
	// ray bundle resolve
	bgfx::FrameBufferHandle rayBundleIntegrateFb;
	bgfx::TextureHandle rayBundleIntegrateTarget;
	bgfx::TextureHandle lightmap;
}
s_bake;

static void bakeInit()
{
	ScreenSpaceVertex::init();
	s_bake.enabled = (bgfx::getCaps()->supported & BGFX_CAPS_FRAMEBUFFER_RW) != 0;
	if (!s_bake.enabled)
		printf("Read/Write frame buffer attachments are not supported. Baking is disabled.\n");
}

static void bakeShutdown()
{
	if (!s_bake.initialized)
		return;
	// shaders
	bgfx::destroy(s_bake.atomicCounterClearProgram);
	bgfx::destroy(s_bake.rayBundleClearProgram);
	bgfx::destroy(s_bake.rayBundleWriteProgram);
	bgfx::destroy(s_bake.rayBundleIntegrateProgram);
	bgfx::destroy(s_bake.lightmapClearProgram);
	bgfx::destroy(s_bake.u_lightmapSize_dataSize);
	bgfx::destroy(s_bake.u_rayNormal);
	bgfx::destroy(s_bake.u_atomicCounterSampler);
	bgfx::destroy(s_bake.u_rayBundleHeaderSampler);
	bgfx::destroy(s_bake.u_rayBundleDataSampler);
	bgfx::destroy(s_bake.u_lightmapSampler);
	// framebuffers
	bgfx::destroy(s_bake.atomicCounterFb);
	bgfx::destroy(s_bake.rayBundleFb);
	bgfx::destroy(s_bake.rayBundleTarget);
	bgfx::destroy(s_bake.rayBundleHeader);
	bgfx::destroy(s_bake.rayBundleData);
	bgfx::destroy(s_bake.lightmapClearFb);
	bgfx::destroy(s_bake.lightmapClearTarget);
	bgfx::destroy(s_bake.rayBundleIntegrateFb);
	bgfx::destroy(s_bake.rayBundleIntegrateTarget);
	bgfx::destroy(s_bake.lightmap);
}

static void bakeExecute()
{
	if (s_bake.executed)
		return;
	if (!s_bake.initialized) {
		// shaders
		s_bake.u_lightmapSize_dataSize = bgfx::createUniform("u_lightmapSize_dataSize", bgfx::UniformType::Vec4);
		s_bake.u_rayNormal = bgfx::createUniform("u_rayNormal", bgfx::UniformType::Vec4);
		s_bake.u_atomicCounterSampler = bgfx::createUniform("u_atomicCounterSampler", bgfx::UniformType::Sampler);
		s_bake.u_rayBundleHeaderSampler = bgfx::createUniform("u_rayBundleHeaderSampler", bgfx::UniformType::Sampler);
		s_bake.u_rayBundleDataSampler = bgfx::createUniform("u_rayBundleDataSampler", bgfx::UniformType::Sampler);
		s_bake.u_lightmapSampler = bgfx::createUniform("u_lightmapSampler", bgfx::UniformType::Sampler);
		s_bake.atomicCounterClearProgram = LOAD_PROGRAM(AtomicCounterClear);
		s_bake.rayBundleClearProgram = LOAD_PROGRAM(RayBundleClear);
		s_bake.rayBundleWriteProgram = LOAD_PROGRAM(RayBundleWrite);
		s_bake.rayBundleIntegrateProgram = LOAD_PROGRAM(RayBundleIntegrate);
		s_bake.lightmapClearProgram = LOAD_PROGRAM(LightmapClear);
		// framebuffers
		{
			bgfx::TextureHandle target = bgfx::createTexture2D(1, 1, false, 1, bgfx::TextureFormat::RGBA8, BGFX_TEXTURE_RT);
			s_bake.atomicCounterTexture = bgfx::createTexture2D(1, 1, false, 1, bgfx::TextureFormat::R32U, BGFX_TEXTURE_COMPUTE_WRITE | BGFX_SAMPLER_POINT | BGFX_SAMPLER_UVW_CLAMP);
			bgfx::Attachment attachments[2];
			attachments[0].init(target);
			attachments[1].init(s_bake.atomicCounterTexture, bgfx::Access::ReadWrite);
			s_bake.atomicCounterFb = bgfx::createFrameBuffer(BX_COUNTOF(attachments), attachments, true);
		}
		{
			s_bake.rayBundleTarget = bgfx::createTexture2D(s_bake.rbTextureSize, s_bake.rbTextureSize, false, 1, bgfx::TextureFormat::RGBA8, BGFX_TEXTURE_RT);
			s_bake.rayBundleHeader = bgfx::createTexture2D(s_bake.rbTextureSize, s_bake.rbTextureSize, false, 1, bgfx::TextureFormat::R32U, BGFX_TEXTURE_COMPUTE_WRITE | BGFX_SAMPLER_POINT | BGFX_SAMPLER_UVW_CLAMP);
			s_bake.rayBundleData = bgfx::createTexture2D(s_bake.rbDataTextureSize, s_bake.rbDataTextureSize, false, 1, bgfx::TextureFormat::RGBA32U, BGFX_TEXTURE_COMPUTE_WRITE | BGFX_SAMPLER_POINT | BGFX_SAMPLER_UVW_CLAMP);
			bgfx::Attachment attachments[4];
			attachments[0].init(s_bake.rayBundleTarget);
			attachments[1].init(s_bake.atomicCounterTexture, bgfx::Access::ReadWrite);
			attachments[2].init(s_bake.rayBundleHeader, bgfx::Access::ReadWrite);
			attachments[3].init(s_bake.rayBundleData, bgfx::Access::ReadWrite); // should be Write, but doesn't work
			s_bake.rayBundleFb = bgfx::createFrameBuffer(BX_COUNTOF(attachments), attachments);
		}
	}
	// Re-create lightmap if atlas resolution has changed.
	if (!s_bake.initialized || s_bake.lightmapWidth != s_atlas.data->width || s_bake.lightmapHeight != s_atlas.data->height) {
		if (s_bake.initialized) {
			bgfx::destroy(s_bake.lightmapClearFb);
			bgfx::destroy(s_bake.lightmapClearTarget);
			bgfx::destroy(s_bake.rayBundleIntegrateFb);
			bgfx::destroy(s_bake.rayBundleIntegrateTarget);
			bgfx::destroy(s_bake.lightmap);
		}
		s_bake.lightmapWidth = s_atlas.data->width;
		s_bake.lightmapHeight = s_atlas.data->height;
		{
			s_bake.rayBundleIntegrateTarget = bgfx::createTexture2D((uint16_t)s_bake.lightmapWidth, (uint16_t)s_bake.lightmapHeight, false, 1, bgfx::TextureFormat::RGBA8, BGFX_TEXTURE_RT);
			s_bake.lightmap = bgfx::createTexture2D((uint16_t)s_bake.lightmapWidth, (uint16_t)s_bake.lightmapHeight, false, 1, bgfx::TextureFormat::RGBA32F, BGFX_TEXTURE_COMPUTE_WRITE);
			bgfx::Attachment attachments[4];
			attachments[0].init(s_bake.rayBundleIntegrateTarget);
			attachments[1].init(s_bake.rayBundleHeader, bgfx::Access::Read);
			attachments[2].init(s_bake.rayBundleData, bgfx::Access::Read);
			attachments[3].init(s_bake.lightmap, bgfx::Access::ReadWrite);
			s_bake.rayBundleIntegrateFb = bgfx::createFrameBuffer(BX_COUNTOF(attachments), attachments);
		}
		{
			s_bake.lightmapClearTarget = bgfx::createTexture2D((uint16_t)s_bake.lightmapWidth, (uint16_t)s_bake.lightmapHeight, false, 1, bgfx::TextureFormat::RGBA8, BGFX_TEXTURE_RT);
			bgfx::Attachment attachments[2];
			attachments[0].init(s_bake.lightmapClearTarget);
			attachments[1].init(s_bake.lightmap, bgfx::Access::ReadWrite);
			s_bake.lightmapClearFb = bgfx::createFrameBuffer(BX_COUNTOF(attachments), attachments);
		}
	}
	s_bake.initialized = true;
	s_bake.executed = true;
	s_bake.directionCount = 0;
}

// https://en.wikipedia.org/wiki/Halton_sequence
static float haltonSequence(int index, int base)
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

static void bakeFrame()
{
	if (!s_bake.executed)
		return;
	bgfx::ViewId viewOffset = 0;
	for (uint32_t i = 0; i < 10; i++) {
		float fsOrtho[16];
		bx::mtxOrtho(fsOrtho, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f, bgfx::getCaps()->homogeneousDepth);
		// Atomic counter clear.
		bgfx::setViewFrameBuffer(viewOffset + kAtomicCounterClearView, s_bake.atomicCounterFb);
		bgfx::setViewRect(viewOffset + kAtomicCounterClearView, 0, 0, 1, 1);
		bgfx::setViewTransform(viewOffset + kAtomicCounterClearView, nullptr, fsOrtho);
		bgfx::setTexture(1, s_bake.u_atomicCounterSampler, s_bake.atomicCounterTexture);
		setScreenSpaceQuadVertexBuffer();
		bgfx::setState(0);
		bgfx::submit(viewOffset + kAtomicCounterClearView, s_bake.atomicCounterClearProgram);
		// Ray bundle clear.
		bgfx::setViewFrameBuffer(viewOffset + kRayBundleClearView, s_bake.rayBundleFb);
		bgfx::setViewRect(viewOffset + kRayBundleClearView, 0, 0, s_bake.rbTextureSize, s_bake.rbTextureSize);
		bgfx::setViewTransform(viewOffset + kRayBundleClearView, nullptr, fsOrtho);
		bgfx::setTexture(2, s_bake.u_rayBundleHeaderSampler, s_bake.rayBundleHeader);
		setScreenSpaceQuadVertexBuffer();
		bgfx::setState(0);
		bgfx::submit(viewOffset + kRayBundleClearView, s_bake.rayBundleClearProgram);
		// Ray bundle write.
		const float rx = -90.0f + haltonSequence(s_bake.directionCount, 2) * 180.0f;
		const float ry = -90.0f + haltonSequence(s_bake.directionCount, 3) * 180.0f;
		const float rz = -90.0f + haltonSequence(s_bake.directionCount, 5) * 180.0f;
		float rotation[16];
		bx::mtxRotateXYZ(rotation, rx, ry, rz);
		float view[16];
		bx::mtxTranspose(view, rotation);
		bx::Vec3 corners[8];
		s_model.aabb.getCorners(corners);
		AABB aabb;
		for (uint32_t j = 0; j < 8; j++) {
			const float in[4] = { corners[j].x, corners[j].y, corners[j].z, 1.0f };
			float out[4];
			bx::vec4MulMtx(out, in, view);
			aabb.addPoint(bx::Vec3(out[0], out[1], out[2]));
		}
		float projection[16];
		bx::mtxOrtho(projection, aabb.min.x, aabb.max.x, aabb.min.y, aabb.max.y, -aabb.max.z, -aabb.min.z, 0.0f, bgfx::getCaps()->homogeneousDepth, bx::Handness::Right);
		bgfx::setViewFrameBuffer(viewOffset + kRayBundleWriteView, s_bake.rayBundleFb);
		bgfx::setViewRect(viewOffset + kRayBundleWriteView, 0, 0, s_bake.rbTextureSize, s_bake.rbTextureSize);
		bgfx::setViewTransform(viewOffset + kRayBundleWriteView, view, projection);
		{
			const float sizes[] = { (float)s_bake.lightmapWidth, (float)s_bake.lightmapHeight, (float)s_bake.rbDataTextureSize, 0.0f };
			for (uint32_t j = 0; j < s_model.data->numMeshes; j++) {
				const objzMesh &mesh = s_model.data->meshes[j];
				const objzMaterial *mat = mesh.materialIndex == -1 ? nullptr : &s_model.data->materials[mesh.materialIndex];
				bgfx::setIndexBuffer(s_atlas.ib, mesh.firstIndex, mesh.numIndices);
				bgfx::setVertexBuffer(0, s_atlas.vb);
				bgfx::setState(0);
				bgfx::setTexture(1, s_bake.u_atomicCounterSampler, s_bake.atomicCounterTexture);
				bgfx::setTexture(2, s_bake.u_rayBundleHeaderSampler, s_bake.rayBundleHeader);
				bgfx::setTexture(3, s_bake.u_rayBundleDataSampler, s_bake.rayBundleData);
				bgfx::setUniform(s_bake.u_lightmapSize_dataSize, sizes);
				modelSetMaterialUniforms(mat);
				bgfx::submit(viewOffset + kRayBundleWriteView, s_bake.rayBundleWriteProgram);
			}
		}
		// Ray bundle resolve.
		bgfx::setViewFrameBuffer(viewOffset + kRayBundleResolveView, s_bake.rayBundleIntegrateFb);
		bgfx::setViewRect(viewOffset + kRayBundleResolveView, 0, 0, (uint16_t)s_bake.lightmapHeight, (uint16_t)s_bake.lightmapHeight);
		bgfx::setViewTransform(viewOffset + kRayBundleResolveView, nullptr, fsOrtho);
		bgfx::setTexture(1, s_bake.u_rayBundleHeaderSampler, s_bake.rayBundleHeader);
		bgfx::setTexture(2, s_bake.u_rayBundleDataSampler, s_bake.rayBundleData);
		bgfx::setTexture(3, s_bake.u_lightmapSampler, s_bake.lightmap);
		const float sizes[] = { (float)s_bake.lightmapWidth, (float)s_bake.lightmapHeight, (float)s_bake.rbDataTextureSize, 0.0f };
		bgfx::setUniform(s_bake.u_lightmapSize_dataSize, sizes);
		const float rayNormal[] = { view[2], view[6], view[10], 0 };
		bgfx::setUniform(s_bake.u_rayNormal, rayNormal);
		setScreenSpaceQuadVertexBuffer();
		bgfx::setState(0);
		bgfx::submit(viewOffset + kRayBundleResolveView, s_bake.rayBundleIntegrateProgram);
		s_bake.directionCount++;
		if (s_bake.directionCount >= s_bake.numDirections) {
			s_bake.executed = false;
			s_bake.finished = true;
			break;
		}
		viewOffset += kNumBakeViews;
	}
}

static void bakeClear()
{
	s_bake.executed = false;
	s_bake.finished = false;
}

int main(int argc, char **argv)
{
	bx::CommandLine commandLine(argc, argv);
	glfwSetErrorCallback(glfw_errorCallback);
	if (!glfwInit())
		return EXIT_FAILURE;
	glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
	s_window = glfwCreateWindow(WINDOW_DEFAULT_WIDTH, WINDOW_DEFAULT_HEIGHT, WINDOW_TITLE, nullptr, nullptr);
	if (!s_window)
		return EXIT_FAILURE;
	glfwMaximizeWindow(s_window);
	bgfx::renderFrame();
	bgfx::Init init;
	if (commandLine.hasArg("gl"))
		init.type = bgfx::RendererType::OpenGL;
#if BX_PLATFORM_LINUX || BX_PLATFORM_BSD
	init.platformData.ndt = glfwGetX11Display();
	init.platformData.nwh = (void*)(uintptr_t)glfwGetX11Window(s_window);
#elif BX_PLATFORM_OSX
	init.platformData.nwh = glfwGetCocoaWindow(s_window);
#elif BX_PLATFORM_WINDOWS
	init.platformData.nwh = glfwGetWin32Window(s_window);
#endif
	int width, height;
	glfwGetWindowSize(s_window, &width, &height);
	init.resolution.width = (uint32_t)width;
	init.resolution.height = (uint32_t)height;
	init.resolution.reset = BGFX_RESET_VSYNC;
	bgfx::init(init);
	guiInit();
	modelInit();
	atlasInit();
	bakeInit();
	glfwSetCharCallback(s_window, glfw_charCallback);
	glfwSetCursorPosCallback(s_window, glfw_cursorPosCallback);
	glfwSetKeyCallback(s_window, glfw_keyCallback);
	glfwSetMouseButtonCallback(s_window, glfw_mouseButtonCallback);
	glfwSetScrollCallback(s_window, glfw_scrollCallback);
	int frameCount = 0, progressDots = 0;
	double lastFrameTime = glfwGetTime();
	while (!glfwWindowShouldClose(s_window)) {
		glfwPollEvents();
		if (glfwGetKey(s_window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
			glfwSetWindowShouldClose(s_window, GLFW_TRUE);
			continue;
		}
		const double currentFrameTime = glfwGetTime();
		const float deltaTime = float(currentFrameTime - lastFrameTime);
		lastFrameTime = currentFrameTime;
		// Handle window resize.
		int oldWidth = width, oldHeight = height;
		glfwGetWindowSize(s_window, &width, &height);
		if (width != oldWidth || height != oldHeight) {
			bgfx::reset((uint32_t)width, (uint32_t)height, BGFX_RESET_VSYNC);
			guiResize(width, height);
			bgfx::setViewRect(kModelView, 0, 0, bgfx::BackbufferRatio::Equal);
		}
		// Update camera.
		if (s_camera.mode == CameraMode::FirstPerson && glfwGetInputMode(s_window, GLFW_CURSOR) == GLFW_CURSOR_DISABLED) {
			const float speed = (s_keyDown[GLFW_KEY_LEFT_SHIFT] ? 20.0f : 5.0f) * deltaTime;
			float deltaForward = 0.0f, deltaRight = 0.0f;
			if (s_keyDown[GLFW_KEY_W]) deltaForward += speed;
			if (s_keyDown[GLFW_KEY_S]) deltaForward -= speed;
			if (s_keyDown[GLFW_KEY_A]) deltaRight -= speed;
			if (s_keyDown[GLFW_KEY_D]) deltaRight += speed;
			s_camera.firstPerson.move(deltaForward, deltaRight);
			if (s_keyDown[GLFW_KEY_Q]) s_camera.firstPerson.position.y -= speed;
			if (s_keyDown[GLFW_KEY_E]) s_camera.firstPerson.position.y += speed;
		}
		float view[16];
		if (s_camera.mode == CameraMode::FirstPerson)
			s_camera.firstPerson.calculateViewMatrix(view);
		else if (s_camera.mode == CameraMode::Orbit)
			s_camera.orbit.calculateViewMatrix(view);
		float projection[16];
		const float ar = width / (float)height;
		bx::mtxProj(projection, s_camera.fov / ar, ar, 0.01f, 1000.0f, bgfx::getCaps()->homogeneousDepth, bx::Handness::Right);
		// GUI
		if (s_options.gui) {
			guiRunFrame(deltaTime);
			ImGui::NewFrame();
			ImGuiIO &io = ImGui::GetIO();
			char errorMessage[1024];
			s_errorMessage.get(errorMessage, sizeof(errorMessage));
			if (errorMessage[0])
				ImGui::OpenPopup("Error");
			if (ImGui::BeginPopupModal("Error", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
				ImGui::Text("%s", errorMessage);
				if (ImGui::Button("OK", ImVec2(120, 0))) {
					ImGui::CloseCurrentPopup();
					s_errorMessage.set(nullptr);
				}
				ImGui::EndPopup();
			}
			const float margin = 4.0f;
			ImGui::SetNextWindowPos(ImVec2(margin, margin), ImGuiCond_FirstUseEver);
			ImGui::SetNextWindowSize(ImVec2(400.0f, io.DisplaySize.y - margin * 2.0f), ImGuiCond_FirstUseEver);
			if (ImGui::Begin("##mainWindow", nullptr, ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse)) {
				ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.5f);
				ImGui::Separator();
				ImGui::Spacing();
				ImGui::Text("Model");
				ImGui::Spacing();
				if (ImGui::Button("Open model...", ImVec2(-1.0f, 0.0f)))
					modelOpenDialog();
				if (s_model.status.get() == ModelStatus::Ready) {
				ImGui::Text("%u objects", s_model.data->numObjects);
				ImGui::Text("%u vertices", s_model.data->numVertices);
				ImGui::Text("%u triangles", s_model.data->numIndices / 3);
				ImGui::InputFloat("Model scale", &s_model.scale, 0.01f, 0.1f);
				s_model.scale = bx::max(0.001f, s_model.scale);
				ImGui::Spacing();
				ImGui::Separator();
				ImGui::Spacing();
				ImGui::Text("View");
				ImGui::Spacing();
				if (s_atlas.status.get() == AtlasStatus::Ready) {
					ImGui::RadioButton("Flat shading", (int *)&s_options.shadeMode, (int)ShadeMode::Flat);
					ImGui::SameLine();
					ImGui::RadioButton("Charts shading", (int *)&s_options.shadeMode, (int)ShadeMode::Charts);
				}
				ImGui::Checkbox("Wireframe overlay", &s_options.wireframe);
				if (s_options.wireframe && s_atlas.status.get() == AtlasStatus::Ready) {
					ImGui::SameLine();
					ImGui::RadioButton("Charts", (int *)&s_options.wireframeMode, (int)WireframeMode::Charts);
					ImGui::SameLine();
					ImGui::RadioButton("Triangles", (int *)&s_options.wireframeMode, (int)WireframeMode::Triangles);
				}
				ImGui::RadioButton("First person camera", (int *)&s_camera.mode, (int)CameraMode::FirstPerson);
				ImGui::SameLine();
				ImGui::TextDisabled("(?)");
				if (ImGui::IsItemHovered()) {
					ImGui::BeginTooltip();
					ImGui::Text("Hold left mouse button on 3D view to enable camera\nW,A,S,D and Q,E to move\nHold SHIFT for faster movement");
					ImGui::EndTooltip();
				}
				ImGui::SameLine();
				ImGui::RadioButton("Orbit camera", (int *)&s_camera.mode, (int)CameraMode::Orbit);
				ImGui::SameLine();
				ImGui::TextDisabled("(?)");
				if (ImGui::IsItemHovered()) {
					ImGui::BeginTooltip();
					ImGui::Text("Hold left mouse button on 3D view to enable camera");
					ImGui::EndTooltip();
				}
				ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.25f);
				ImGui::DragFloat("FOV", &s_camera.fov, 1.0f, 45.0f, 150.0f, "%.0f");
				ImGui::SameLine();
				ImGui::Text(" ");
				ImGui::SameLine();
				ImGui::DragFloat("Sensitivity", &s_camera.sensitivity, 0.01f, 0.01f, 1.0f);
				ImGui::PopItemWidth();
					ImGui::Spacing();
					ImGui::Separator();
					ImGui::Spacing();
					ImGui::Text("Chart options");
					ImGui::Spacing();
					ImGui::InputFloat("Proxy fit metric weight", &s_atlas.chartOptions.proxyFitMetricWeight);
					ImGui::InputFloat("Roundness metric weight", &s_atlas.chartOptions.roundnessMetricWeight);
					ImGui::InputFloat("Straightness metric weight", &s_atlas.chartOptions.straightnessMetricWeight);
					ImGui::InputFloat("NormalSeam metric weight", &s_atlas.chartOptions.normalSeamMetricWeight);
					ImGui::InputFloat("Texture seam metric weight", &s_atlas.chartOptions.textureSeamMetricWeight);
					ImGui::InputFloat("Max chart area", &s_atlas.chartOptions.maxChartArea);
					ImGui::InputFloat("Max boundary length", &s_atlas.chartOptions.maxBoundaryLength);
					ImGui::InputFloat("Max threshold", &s_atlas.chartOptions.maxThreshold);
					ImGui::InputInt("Grow face count", (int *)&s_atlas.chartOptions.growFaceCount);
					ImGui::InputInt("Max iterations", (int *)&s_atlas.chartOptions.maxIterations);
					ImGui::Spacing();
					ImGui::Separator();
					ImGui::Spacing();
#if USE_LIBIGL
					ImGui::Text("Parameterization options");
					ImGui::Spacing();
					const ParamMethod oldParamMethod = s_atlas.paramMethod;
					ImGui::RadioButton("LSCM", (int *)&s_atlas.paramMethod, (int)ParamMethod::LSCM);
					ImGui::RadioButton("libigl Harmonic", (int *)&s_atlas.paramMethod, (int)ParamMethod::libigl_Harmonic);
					ImGui::RadioButton("libigl LSCM", (int *)&s_atlas.paramMethod, (int)ParamMethod::libigl_LSCM);
					ImGui::RadioButton("libigl ARAP", (int *)&s_atlas.paramMethod, (int)ParamMethod::libigl_ARAP);
					if (s_atlas.paramMethod != oldParamMethod)
						s_atlas.paramMethodChanged = true;
					ImGui::Spacing();
					ImGui::Separator();
					ImGui::Spacing();
#endif
					ImGui::Text("Pack options");
					ImGui::Spacing();
					ImGui::SliderInt("Attempts", &s_atlas.packOptions.attempts, 0, 4096);
					ImGui::InputFloat("Texels per unit", &s_atlas.packOptions.texelsPerUnit, 0.0f, 32.0f, 2);
					ImGui::InputInt("Resolution", (int *)&s_atlas.packOptions.resolution, 8);
					ImGui::InputInt("Max chart size", (int *)&s_atlas.packOptions.maxChartSize);
					ImGui::Checkbox("Block align", &s_atlas.packOptions.blockAlign);
					ImGui::SameLine();
					ImGui::Checkbox("Conservative", &s_atlas.packOptions.conservative);
					ImGui::SliderInt("Padding", &s_atlas.packOptions.padding, 0, 8);
					ImGui::Spacing();
					ImGui::Separator();
					ImGui::Spacing();
					ImGui::Text("Atlas");
					ImGui::Checkbox("Verbose output", &s_atlas.verbose);
					if (ImGui::Button("Generate atlas", ImVec2(-1.0f, 0.0f)))
						atlasGenerate();
					if (s_atlas.status.get() == AtlasStatus::Ready) {
						uint32_t numIndices = 0, numVertices = 0;
						for (uint32_t i = 0; i < s_atlas.data->meshCount; i++) {
							const xatlas::Mesh &outputMesh = s_atlas.data->meshes[i];
							numIndices += outputMesh.indexCount;
							numVertices += outputMesh.vertexCount;
						}
						ImGui::Text("%u atlases", s_atlas.data->atlasCount);
						ImGui::Text("%ux%u resolution", s_atlas.data->width, s_atlas.data->height);
						ImGui::Text("%u charts", s_atlas.data->chartCount);
						ImGui::Text("%u vertices", numVertices);
						ImGui::Text("%u triangles", numIndices / 3);
						ImGui::Checkbox("Show atlas", &s_atlas.showTexture);
						if (s_bake.enabled) {
							ImGui::Spacing();
							ImGui::Separator();
							ImGui::Spacing();
							ImGui::Text("Bake");
							ImGui::SliderInt("Num. directions", &s_bake.numDirections, 300, 10000);
							if (ImGui::Button("Bake", ImVec2(-1.0f, 0.0f)))
								bakeExecute();
						}
					}
				}
				ImGui::PopItemWidth();
				ImGui::End();
			}
			const ImGuiWindowFlags progressWindowFlags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings;
			if (s_model.status.get() == ModelStatus::Loading) {
				ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y * 0.5f), ImGuiCond_Always, ImVec2(0.5f, 0.5f));
				if (ImGui::Begin("##modelProgress", nullptr, progressWindowFlags)) {
					ImGui::Text("Loading model");
					for (int i = 0; i < 3; i++) {
						ImGui::SameLine();
						ImGui::Text(i < progressDots ? "." : " ");
					}
					ImGui::End();
				}
			}
			const AtlasStatus::Enum atlasStatus = s_atlas.status.get();
			if (atlasStatus == AtlasStatus::AddingMeshes || atlasStatus == AtlasStatus::Generating) {
				ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x - margin, margin), ImGuiCond_Always, ImVec2(1.0f, 0.0f));
				ImGui::SetNextWindowSize(ImVec2(300.0f, -1.0f), ImGuiCond_Always);
				if (ImGui::Begin("##atlasProgress", nullptr, progressWindowFlags)) {
					int progress;
					xatlas::ProgressCategory::Enum category;
					s_atlas.status.getProgress(&category, &progress);
					if (atlasStatus == AtlasStatus::AddingMeshes)
						ImGui::Text("Adding meshes");
					else
						ImGui::Text("%s", xatlas::StringForEnum(category));
					for (int i = 0; i < 3; i++) {
						ImGui::SameLine();
						ImGui::Text(i < progressDots ? "." : " ");
					}
					ImGui::ProgressBar(progress / 100.0f);
					ImGui::End();
				}
			} else if (atlasStatus == AtlasStatus::Ready && s_atlas.showTexture) {
				const float size = 500;
				ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x - size - margin, margin), ImGuiCond_FirstUseEver);
				ImGui::SetNextWindowSize(ImVec2(size, size), ImGuiCond_FirstUseEver);
				if (ImGui::Begin("Atlas", &s_atlas.showTexture)) {
					if (s_atlas.data->atlasCount > 1) {
						ImGui::Text("Atlas %d of %u", s_atlas.currentTexture + 1, s_atlas.data->atlasCount);
						ImGui::SameLine();
						if (ImGui::ArrowButton("##prevAtlas", ImGuiDir_Left)) {
							s_atlas.currentTexture--;
							if (s_atlas.currentTexture < 0)
								s_atlas.currentTexture = s_atlas.data->atlasCount - 1;
						}
						ImGui::SameLine();
						if (ImGui::ArrowButton("##nextAtlas", ImGuiDir_Right)) {
							s_atlas.currentTexture++;
							if (s_atlas.currentTexture > (int)s_atlas.data->atlasCount - 1)
								s_atlas.currentTexture = 0;
						}
					}
					const ImVec2 pos = ImGui::GetCursorScreenPos();
					ImTextureID texture = (ImTextureID)(intptr_t)s_atlas.chartsTextures[s_atlas.currentTexture].idx;
					ImGui::Image(texture, ImGui::GetContentRegionAvail());
					if (ImGui::IsItemHovered()) {
						const ImVec2 textureSize((float)s_atlas.data->width, (float)s_atlas.data->height);
						const ImVec2 imageSize(ImGui::GetItemRectSize());
						const ImVec2 imageToTex(textureSize.x / imageSize.x, textureSize.y / imageSize.y);
						const float magnifiedSize = 200.0f;
						const ImVec2 uv0 = ImVec2((io.MousePos.x - pos.x) * imageToTex.x - magnifiedSize * 0.5f, (io.MousePos.y - pos.y) * imageToTex.y - magnifiedSize * 0.5f);
						const ImVec2 uv1 = ImVec2(uv0.x + magnifiedSize, uv0.y + magnifiedSize);
						ImGui::BeginTooltip();
						ImGui::Image(texture, ImVec2(magnifiedSize, magnifiedSize), ImVec2(uv0.x / textureSize.x, uv0.y / textureSize.y), ImVec2(uv1.x / textureSize.x, uv1.y / textureSize.y));
						ImGui::EndTooltip();
					}
					ImGui::End();
				}
				if (s_bake.executed || s_bake.finished) {
					ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x - size - margin, size + margin * 2.0f), ImGuiCond_FirstUseEver);
					ImGui::SetNextWindowSize(ImVec2(size, size), ImGuiCond_FirstUseEver);
					if (ImGui::Begin("Bake", &s_atlas.showTexture)) {
						ImTextureID texture = (ImTextureID)(intptr_t)s_bake.lightmap.idx;
						ImGui::Image(texture, ImGui::GetContentRegionAvail(), ImVec2(0.0f, 1.0f), ImVec2(1.0f, 0.0f));
						ImGui::End();
					}
				}
			}
		}
		modelRender(view, projection);
		bakeFrame();
		if (s_options.gui)
			guiRender();
		bgfx::touch(kModelView);
		bgfx::setDebug(s_showBgfxStats ? BGFX_DEBUG_STATS : BGFX_DEBUG_NONE);
		bgfx::frame();
		frameCount++;
		if (frameCount % 20 == 0)
			progressDots = (progressDots + 1) % 4;
		if (s_model.status.get() == ModelStatus::Finalizing)
			modelFinalize();
		if (s_atlas.status.get() == AtlasStatus::Finalizing)
			atlasFinalize();
	}
	guiShutdown();
	bakeShutdown();
	atlasDestroy();
	modelShutdown();
	bgfx::shutdown();
	glfwTerminate();
	return 0;
}
