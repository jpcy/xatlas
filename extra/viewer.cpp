/*
xatlas
https://github.com/jpcy/xatlas
Copyright (c) 2018 Jonathan Young

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#include <mutex>
#include <stdio.h>
#include <bx/os.h>
#include <bgfx/platform.h>
#include <GLFW/glfw3.h>
#if BX_PLATFORM_LINUX
#define GLFW_EXPOSE_NATIVE_X11
#elif BX_PLATFORM_WINDOWS
#define GLFW_EXPOSE_NATIVE_WIN32
#endif
#include <GLFW/glfw3native.h>
#undef Success
#include <imgui/imgui.h>
#include <imgui/imgui_internal.h>
#include "shaders_bin/shaders.h"
#include "viewer.h"

// Progress spinner by zfedoran
// https://github.com/ocornut/imgui/issues/1901
namespace ImGui {
bool Spinner(const char* label, float radius, float thickness, const ImU32& color) {
	ImGuiContext& g = *GetCurrentContext();
	ImGuiWindow* window = g.CurrentWindow;
	if (window->SkipItems)
		return false;

	const ImGuiStyle& style = g.Style;
	const ImGuiID id = window->GetID(label);

	ImVec2 pos = window->DC.CursorPos;
	ImVec2 size((radius )*2, (radius + style.FramePadding.y)*2);

	const ImRect bb(pos, ImVec2(pos.x + size.x, pos.y + size.y));
	ItemSize(bb, style.FramePadding.y);
	if (!ItemAdd(bb, id))
		return false;

	// Render
	window->DrawList->PathClear();

	int num_segments = 30;
	int start = (int)fabsf(ImSin((float)g.Time*1.8f)*(num_segments-5));

	const float a_min = IM_PI*2.0f * ((float)start) / (float)num_segments;
	const float a_max = IM_PI*2.0f * ((float)num_segments-3) / (float)num_segments;

	const ImVec2 centre = ImVec2(pos.x+radius, pos.y+radius+style.FramePadding.y);

	for (int i = 0; i < num_segments; i++) {
		const float a = a_min + ((float)i / (float)num_segments) * (a_max - a_min);
		window->DrawList->PathLineTo(ImVec2(centre.x + ImCos(a+(float)g.Time*8) * radius,
			centre.y + ImSin(a+(float)g.Time*8) * radius));
	}

	window->DrawList->PathStroke(color, false, thickness);
	return true;
}

bool Spinner(const char* label)
{
	const ImU32 col = ImGui::GetColorU32(ImGuiCol_ButtonHovered);
	return Spinner(label, GetCurrentContext()->Font->FontSize * 0.5f, 2.0f, col);
}
} // namespace ImGui

#define WINDOW_DEFAULT_WIDTH 1920
#define WINDOW_DEFAULT_HEIGHT 1080

bgfx::VertexLayout WireframeVertex::layout;
Options g_options;
GLFWwindow *g_window;
int g_windowSize[2];
static bool s_keyDown[GLFW_KEY_LAST + 1] = { 0 };
static bool s_showBgfxStats = false;

void randomRGB(uint8_t *color)
{
	const int mix = 192;
	color[0] = uint8_t((rand() % 255 + mix) * 0.5f);
	color[1] = uint8_t((rand() % 255 + mix) * 0.5f);
	color[2] = uint8_t((rand() % 255 + mix) * 0.5f);
}

uint32_t encodeRGBA(const uint8_t *rgba)
{
	return rgba[3] << 24 | rgba[2] << 16 | rgba[1] << 8 | rgba[0];
}

void decodeRGBA(uint32_t rgbaIn, uint8_t *rgbaOut)
{
	rgbaOut[0] = rgbaIn & 0xff;
	rgbaOut[1] = (rgbaIn >> 8) & 0xff;
	rgbaOut[2] = (rgbaIn >> 16) & 0xff;
	rgbaOut[3] = (rgbaIn >> 24) & 0xff;
}

struct
{
	char text[1024] = { 0 };
	std::mutex mutex;
}
s_errorMessage;

void setErrorMessage(const char *format, ...)
{
	s_errorMessage.mutex.lock();
	if (format) {
		va_list args;
		va_start(args, format);
		vsnprintf(s_errorMessage.text, sizeof(s_errorMessage.text), format, args);
		va_end(args);
	} else
		s_errorMessage.text[0] = 0;
	s_errorMessage.mutex.unlock();
}

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
		const bx::Vec3 center = modelGetCentroid();
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

void resetCamera()
{
	s_camera.firstPerson = FirstPersonCamera();
	s_camera.orbit = OrbitCamera();
}

static void glfw_errorCallback(int error, const char *description)
{
	fprintf(stderr, "GLFW error %d: %s\n", error, description);
}

static void glfw_charCallback(GLFWwindow * /*window*/, unsigned int c)
{
	if (!g_options.gui)
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
	if (glfwGetInputMode(g_window, GLFW_CURSOR) == GLFW_CURSOR_DISABLED) {
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
		g_options.gui = !g_options.gui;
	if (g_options.gui) {
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
	if (glfwGetInputMode(g_window, GLFW_CURSOR) == GLFW_CURSOR_DISABLED && action == GLFW_RELEASE)
		glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
	if (g_options.gui) {
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
	if (g_options.gui && io.WantCaptureMouse)
		io.MouseWheel += (float)yoffset;
	else if (s_camera.mode == CameraMode::Orbit) 
		s_camera.orbit.zoom((float)-yoffset);
}

struct ShaderSource
{
	const uint8_t *data;
	uint32_t size;
};

struct ShaderSourceBundle
{
	const char *name;
#if BX_PLATFORM_WINDOWS
	ShaderSource d3d11;
#endif
	ShaderSource gl;
	ShaderSource vk;
};

#if BX_PLATFORM_WINDOWS
#define SHADER_SOURCE_BUNDLE(name) { BX_STRINGIZE(name), { name##_d3d11, sizeof(name##_d3d11) }, { name##_gl, sizeof(name##_gl) }, { name##_vk, sizeof(name##_vk) }}
#else
#define SHADER_SOURCE_BUNDLE(name) { BX_STRINGIZE(name), { name##_gl, sizeof(name##_gl) }, { name##_vk, sizeof(name##_vk) }}
#endif

// Sync with ShaderId
static ShaderSourceBundle s_shaders[] = 
{
	SHADER_SOURCE_BUNDLE(fs_blit),
	SHADER_SOURCE_BUNDLE(fs_chart),
	SHADER_SOURCE_BUNDLE(fs_color),
	SHADER_SOURCE_BUNDLE(fs_gui),
	SHADER_SOURCE_BUNDLE(fs_material),
	SHADER_SOURCE_BUNDLE(fs_wireframe),
	SHADER_SOURCE_BUNDLE(vs_blit),
	SHADER_SOURCE_BUNDLE(vs_chart),
	SHADER_SOURCE_BUNDLE(vs_color),
	SHADER_SOURCE_BUNDLE(vs_gui),
	SHADER_SOURCE_BUNDLE(vs_model),
	SHADER_SOURCE_BUNDLE(vs_wireframe)
};

bgfx::ShaderHandle loadShader(ShaderId id)
{
	const ShaderSourceBundle &sourceBundle = s_shaders[(int)id];
	ShaderSource source;
	if (bgfx::getRendererType() == bgfx::RendererType::OpenGL)
		source = sourceBundle.gl;
#if BX_PLATFORM_WINDOWS
	else if (bgfx::getRendererType() == bgfx::RendererType::Direct3D11)
		source = sourceBundle.d3d11;
#endif
	else if (bgfx::getRendererType() == bgfx::RendererType::Vulkan)
		source = sourceBundle.vk;
	else {
		fprintf(stderr, "Unsupported renderer type.");
		exit(EXIT_FAILURE);
	}
	bgfx::ShaderHandle shader = bgfx::createShader(bgfx::makeRef(source.data, source.size));
	if (!bgfx::isValid(shader)) {
		fprintf(stderr, "Creating shader '%s' failed.", sourceBundle.name);
		exit(EXIT_FAILURE);
	}
#if _DEBUG
	bgfx::setName(shader, sourceBundle.name);
#endif
	return shader;
}

struct
{
	bgfx::ShaderHandle vs_color, fs_color;
	bgfx::ProgramHandle colorProgram;
	bgfx::UniformHandle u_thickness;
	bgfx::ShaderHandle vs_wireframe, fs_wireframe;
	bgfx::ProgramHandle wireframeProgram;
}
s_commonShaders;

static void commonShadersInit()
{
	s_commonShaders.vs_color = loadShader(ShaderId::vs_color);
	s_commonShaders.fs_color = loadShader(ShaderId::fs_color);
	s_commonShaders.colorProgram = bgfx::createProgram(s_commonShaders.vs_color, s_commonShaders.fs_color);
	s_commonShaders.u_thickness = bgfx::createUniform("u_thickness", bgfx::UniformType::Vec4);
	s_commonShaders.vs_wireframe = loadShader(ShaderId::vs_wireframe);
	s_commonShaders.fs_wireframe = loadShader(ShaderId::fs_wireframe);
	s_commonShaders.wireframeProgram = bgfx::createProgram(s_commonShaders.vs_wireframe, s_commonShaders.fs_wireframe, true);
}

static void commonShadersShutdown()
{
	bgfx::destroy(s_commonShaders.vs_color);
	bgfx::destroy(s_commonShaders.fs_color);
	bgfx::destroy(s_commonShaders.colorProgram);
	bgfx::destroy(s_commonShaders.u_thickness);
	bgfx::destroy(s_commonShaders.wireframeProgram);
}

bgfx::ProgramHandle getColorProgram()
{
	return s_commonShaders.colorProgram;
}

void setWireframeThicknessUniform(float thickness)
{
	if (bgfx::getRendererType() == bgfx::RendererType::Direct3D11)
		thickness *= 0.75f;
	const float data[] = { thickness, 0.0f, 0.0f, 0.0f };
	bgfx::setUniform(s_commonShaders.u_thickness, data);
}

bgfx::ProgramHandle getWireframeProgram()
{
	return s_commonShaders.wireframeProgram;
}

struct BgfxCallback : public bgfx::CallbackI
{
	virtual ~BgfxCallback() {}

	void fatal(const char* /*_filePath*/, uint16_t /*_line*/, bgfx::Fatal::Enum /*_code*/, const char* _str) override
	{
		fprintf(stderr, "%s\n", _str);
		exit(EXIT_FAILURE);
	}

	void traceVargs(const char* /*_filePath*/, uint16_t /*_line*/, const char* _format, va_list _argList) override
	{
		bx::debugPrintfVargs(_format, _argList);
	}

	void profilerBegin(const char*, uint32_t, const char*, uint16_t) override {}
	void profilerBeginLiteral(const char*, uint32_t, const char*, uint16_t) override {}
	void profilerEnd() override {}
	uint32_t cacheReadSize(uint64_t) override { return 0; }
	bool cacheRead(uint64_t, void*, uint32_t) override { return false; }
	void cacheWrite(uint64_t, const void*, uint32_t) override {}
	void screenShot(const char*, uint32_t, uint32_t, uint32_t, const void*, uint32_t, bool) override {}
	void captureBegin(uint32_t, uint32_t, uint32_t, bgfx::TextureFormat::Enum, bool) override {}
	void captureEnd() override {}
	void captureFrame(const void*, uint32_t) override {}
};

int main(int argc, char **argv)
{
	glfwSetErrorCallback(glfw_errorCallback);
	if (!glfwInit())
		return EXIT_FAILURE;
	glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
	g_window = glfwCreateWindow(WINDOW_DEFAULT_WIDTH, WINDOW_DEFAULT_HEIGHT, WINDOW_TITLE, nullptr, nullptr);
	if (!g_window)
		return EXIT_FAILURE;
	glfwMaximizeWindow(g_window);
	bgfx::renderFrame();
	BgfxCallback bgfxCallback;
	bgfx::Init init;
	//init.type = bgfx::RendererType::OpenGL;
	init.callback = &bgfxCallback;
#if BX_PLATFORM_LINUX || BX_PLATFORM_BSD
	init.platformData.ndt = glfwGetX11Display();
	init.platformData.nwh = (void*)(uintptr_t)glfwGetX11Window(g_window);
#elif BX_PLATFORM_OSX
	init.platformData.nwh = glfwGetCocoaWindow(g_window);
#elif BX_PLATFORM_WINDOWS
	init.platformData.nwh = glfwGetWin32Window(g_window);
#endif
	int width, height;
	glfwGetWindowSize(g_window, &width, &height);
	g_windowSize[0] = width;
	g_windowSize[1] = height;
	init.resolution.width = (uint32_t)width;
	init.resolution.height = (uint32_t)height;
	init.resolution.reset = BGFX_RESET_VSYNC | BGFX_RESET_MSAA_X16;
	bgfx::init(init);
	WireframeVertex::init();
	commonShadersInit();
	guiInit();
	modelInit();
	atlasInit();
	bakeInit();
	glfwSetCharCallback(g_window, glfw_charCallback);
	glfwSetCursorPosCallback(g_window, glfw_cursorPosCallback);
	glfwSetKeyCallback(g_window, glfw_keyCallback);
	glfwSetMouseButtonCallback(g_window, glfw_mouseButtonCallback);
	glfwSetScrollCallback(g_window, glfw_scrollCallback);
	if (argc >= 2)
		modelOpen(argv[1]);
	double lastFrameTime = glfwGetTime();
	uint32_t bgfxFrameNo = 0;
	bool initDockLayout = true;
	while (!glfwWindowShouldClose(g_window)) {
		glfwPollEvents();
		while (glfwGetWindowAttrib(g_window, GLFW_ICONIFIED)) {
			glfwWaitEvents();
			lastFrameTime = glfwGetTime();
		}
		if (glfwGetKey(g_window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
			glfwSetWindowShouldClose(g_window, GLFW_TRUE);
			continue;
		}
		const double currentFrameTime = glfwGetTime();
		const float deltaTime = float(currentFrameTime - lastFrameTime);
		lastFrameTime = currentFrameTime;
		// Handle window resize.
		int oldWidth = width, oldHeight = height;
		glfwGetWindowSize(g_window, &width, &height);
		if (width != oldWidth || height != oldHeight) {
			if (width != 0 && height != 0) {
				g_windowSize[0] = width;
				g_windowSize[1] = height;
			}
			bgfx::reset((uint32_t)width, (uint32_t)height, BGFX_RESET_VSYNC | BGFX_RESET_MSAA_X16);
			guiResize(width, height);
			bgfx::setViewRect(kModelView, 0, 0, bgfx::BackbufferRatio::Equal);
			bgfx::setViewRect(kModelTransparentView, 0, 0, bgfx::BackbufferRatio::Equal);
		}
		// Update camera.
		if (s_camera.mode == CameraMode::FirstPerson && glfwGetInputMode(g_window, GLFW_CURSOR) == GLFW_CURSOR_DISABLED) {
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
		if (g_options.gui) {
			guiRunFrame(deltaTime);
			ImGui::NewFrame();
			ImGuiIO &io = ImGui::GetIO();
			char errorMessage[1024];
			s_errorMessage.mutex.lock();
			bx::strCopy(errorMessage, sizeof(errorMessage), s_errorMessage.text);
			s_errorMessage.mutex.unlock();
			if (errorMessage[0])
				ImGui::OpenPopup("Error");
			if (ImGui::BeginPopupModal("Error", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
				ImGui::Text("%s", errorMessage);
				if (ImGui::Button("OK", ImVec2(120, 0))) {
					ImGui::CloseCurrentPopup();
					setErrorMessage(nullptr);
				}
				ImGui::EndPopup();
			}
			ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
			ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
			ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
			ImGui::SetNextWindowPos(ImVec2(0.0f, 0.0f));
			ImGui::SetNextWindowSize(io.DisplaySize);
			ImGui::SetNextWindowBgAlpha(0.0f);
			ImGui::Begin("DockSpaceWindow", nullptr, ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus);
			ImGui::PopStyleVar(3);
			ImGuiID dockSpaceId = ImGui::GetID("DockSpace");
			ImGui::DockSpace(dockSpaceId, ImVec2(0.0f, 0.0f), ImGuiDockNodeFlags_PassthruCentralNode);
			const char *atlasOptionsWindowName = ICON_FA_GLOBE "Atlas Options";
			if (initDockLayout) {
				ImGuiID dockLeftId = ImGui::DockBuilderSplitNode(dockSpaceId, ImGuiDir_Left, 0.20f, nullptr, &dockSpaceId);
				ImGuiID dockRightId = ImGui::DockBuilderSplitNode(dockSpaceId, ImGuiDir_Right, 0.30f, nullptr, &dockSpaceId);
				ImGuiID dockRightBottomId = ImGui::DockBuilderSplitNode(dockRightId, ImGuiDir_Down, 0.50f, nullptr, &dockRightId);
				ImGui::DockBuilderDockWindow(atlasOptionsWindowName, dockLeftId);
				ImGui::DockBuilderDockWindow("Atlas", dockRightId);
				ImGui::DockBuilderDockWindow("Lightmap", dockRightBottomId);
				ImGui::DockBuilderFinish(dockSpaceId);
				initDockLayout = false;
			}
			ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(5.0f, 5.0f));
			if (ImGui::BeginMainMenuBar()) {
				ImGui::PopStyleVar();
				if (ImGui::BeginMenu(ICON_FA_CUBE " Model")) {
					if (ImGui::MenuItem(ICON_FA_FOLDER_OPEN " Open..."))
						modelOpenDialog();
					if (modelIsLoaded()) {
						ImGui::Spacing();
						ImGui::Separator();
						ImGui::Spacing();
						modelShowGuiMenu();
					}
					ImGui::EndMenu();
				}
				if (ImGui::BeginMenu(ICON_FA_CAMERA " Camera")) {
					ImGui::RadioButton("First person", (int *)&s_camera.mode, (int)CameraMode::FirstPerson);
					ImGui::SameLine();
					ImGui::TextDisabled("(?)");
					if (ImGui::IsItemHovered()) {
						ImGui::BeginTooltip();
						ImGui::Text("Hold left mouse button on 3D view to enable camera\nW,A,S,D and Q,E to move\nHold SHIFT for faster movement");
						ImGui::EndTooltip();
					}
					ImGui::RadioButton("Orbit", (int *)&s_camera.mode, (int)CameraMode::Orbit);
					ImGui::SameLine();
					ImGui::TextDisabled("(?)");
					if (ImGui::IsItemHovered()) {
						ImGui::BeginTooltip();
						ImGui::Text("Hold left mouse button on 3D view to enable camera");
						ImGui::EndTooltip();
					}
					ImGui::Spacing();
					ImGui::Separator();
					ImGui::Spacing();
					ImGui::PushItemWidth(100.0f);
					ImGui::DragFloat("Sensitivity", &s_camera.sensitivity, 0.01f, 0.01f, 1.0f);
					ImGui::PopItemWidth();
					ImGui::EndMenu();
				}
				if (ImGui::BeginMenu(ICON_FA_EYE " View")) {
					ImGui::Checkbox("Wireframe", &g_options.wireframe);
					if (g_options.wireframe && atlasIsReady()) {
						if (ImGui::BeginMenu("Wireframe mode")) {
							ImGui::RadioButton("Charts##wireframe", (int *)&g_options.wireframeMode, (int)WireframeMode::Charts);
							ImGui::RadioButton("Triangles", (int *)&g_options.wireframeMode, (int)WireframeMode::Triangles);
							ImGui::EndMenu();
						}
					}
					if (atlasIsReady()) {
						if (ImGui::BeginMenu("Shading")) {
							ImGui::RadioButton("Flat", (int *)&g_options.shadeMode, (int)ShadeMode::Flat);
							ImGui::RadioButton("Charts##shading", (int *)&g_options.shadeMode, (int)ShadeMode::Charts);
							if (bakeIsLightmapReady()) {
								ImGui::RadioButton("Lightmap", (int *)&g_options.shadeMode, (int)ShadeMode::Lightmap);
								ImGui::RadioButton("Lightmap only", (int *)&g_options.shadeMode, (int)ShadeMode::LightmapOnly);
							}
							ImGui::EndMenu();
						}
						if (g_options.shadeMode == ShadeMode::Charts) {
							if (ImGui::BeginMenu("Chart color")) {
								ImGui::RadioButton("Individual", (int *)&g_options.chartColorMode, (int)ChartColorMode::Individual);
								ImGui::RadioButton("Invalid", (int *)&g_options.chartColorMode, (int)ChartColorMode::Invalid);
								ImGui::EndMenu();
							}
							ImGui::PushItemWidth(100.0f);
							ImGui::SliderInt("Chart cell size", &g_options.chartCellSize, 1, 32);
							ImGui::PopItemWidth();
						}
						if (g_options.shadeMode == ShadeMode::Lightmap || g_options.shadeMode == ShadeMode::LightmapOnly) {
							ImGui::Checkbox("Lightmap point sampling", &g_options.lightmapPointSampling);
							if (bakeIsDenoised())
								ImGui::Checkbox("Use denoised lightmap", &g_options.useDenoisedLightmap);
						}
					}
					ImGui::EndMenu();
				}
				if (atlasIsReady() || bakeIsLightmapReady()) {
					if (ImGui::BeginMenu(ICON_FA_WINDOWS " Window")) {
						ImGui::MenuItem("Atlas Options", nullptr, &g_options.showAtlasOptionsWindow);
						if (atlasIsReady())
							ImGui::MenuItem("Atlas", nullptr, &g_options.showAtlasWindow);
						if (bakeIsLightmapReady())
							ImGui::MenuItem("Lightmap", nullptr, &g_options.showLightmapWindow);
						ImGui::EndMenu();
					}
				}
				ImGui::EndMainMenuBar();
			} else {
				ImGui::PopStyleVar();
			}
			ImGui::End(); // DockSpaceWindow
			if (ImGui::Begin(atlasOptionsWindowName, &g_options.showAtlasOptionsWindow)) {
				ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.5f);
				if (modelIsLoaded()) {
					atlasShowGuiOptions();
					if (atlasIsReady()) {
						ImGui::Spacing();
						ImGui::Separator();
						ImGui::Spacing();
						bakeShowGuiOptions();
					}
				}
				ImGui::PopItemWidth();
				ImGui::End();
			}
			modelShowGuiWindow();
			atlasShowGuiWindow();
			bakeShowGuiWindow();
		}
		modelRender(view, projection);
		bakeFrame(bgfxFrameNo);
		if (g_options.gui)
			guiRender();
		bgfx::touch(kModelView);
		bgfx::setDebug(s_showBgfxStats ? BGFX_DEBUG_STATS : BGFX_DEBUG_NONE);
		bgfxFrameNo = bgfx::frame();
		modelFinalize();
		atlasFinalize();
	}
	commonShadersShutdown();
	guiShutdown();
	bakeShutdown();
	atlasDestroy();
	atlasShutdown();
	modelShutdown();
	bgfx::shutdown();
	glfwTerminate();
	return 0;
}
