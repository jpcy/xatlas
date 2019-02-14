/*
https://github.com/jpcy/xatlas

Copyright (c) 2018 Jonathan Young
Copyright (c) 2013 Thekla, Inc
Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#include <mutex>
#include <thread>
#include <unordered_map>
#include <vector>
#include <stdio.h>
#include "flextGL.h"
#include "GLFW/glfw3.h"

#if defined(__clang__) || defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-braces"
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#include "HandmadeMath.h"
#pragma GCC diagnostic pop
#else
#include "HandmadeMath.h"
#endif

#include "imgui/imgui.h"
#include "nativefiledialog/nfd.h"
#include "objzero/objzero.h"
#include "../xatlas.h"

#define WINDOW_TITLE "xatlas viewer"
#define WINDOW_DEFAULT_WIDTH 1920
#define WINDOW_DEFAULT_HEIGHT 1080

#ifdef _MSC_VER
#define STRNCPY(_dest, _destSize, _src) strncpy_s(_dest, _destSize, _src, (_destSize) - 1)
#else
#define STRNCPY(_dest, _destSize, _src) strncpy(_dest, _src, (_destSize) - 1)
#endif

static GLFWwindow *s_window;
static bool s_keyDown[GLFW_KEY_LAST + 1] = { 0 };

struct
{
	GLuint id;
	GLuint u_color;
	GLint u_mvp;
}
s_colorShader;

struct
{
	GLuint id;
	GLuint u_color;
	GLint u_mvp;
	GLint u_textureSize_cellSize;
}
s_texcoordShader;

struct
{
	GLuint fontTexture;
	GLuint shaderProgram;
	GLint u_texture;
	GLint u_mvp;
	GLuint vao, vbo, ibo;
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
	hmm_vec3 pos;
	hmm_vec3 normal;
	hmm_vec4 texcoord;
};

struct
{
	ModelStatus status;
	std::thread *thread = nullptr;
	objzModel *data;
	hmm_vec3 centroid = HMM_Vec3(0.0f, 0.0f, 0.0f);
	GLuint vao = 0, vbo, ibo;
	float scale = 1.0f;
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
	std::vector<GLuint> chartsTextures;
	std::vector<std::vector<uint8_t>> chartsImages;
	GLuint chartVao = 0, chartVbo, chartIbo;
	GLuint chartBoundaryVao = 0, chartBoundaryVbo;
	xatlas::CharterOptions charterOptions;
	xatlas::PackerOptions packerOptions;
	std::vector<ModelVertex> chartVertices;
	std::vector<uint32_t> chartIndices;
	std::vector<hmm_vec3> chartBoundaryVertices;
}
s_atlas;

enum class WireframeMode
{
	Charts,
	Triangles
};

struct
{
	bool gui = true;
	bool wireframe = true;
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

static void axisFromEulerAngles(float pitch, float yaw, hmm_vec3 *forward, hmm_vec3 *right, hmm_vec3 *up)
{
	const float ryaw = HMM_ToRadians(yaw);
	const float rpitch = HMM_ToRadians(pitch);
	const hmm_vec3 f = HMM_Vec3
	(
		HMM_SinF(ryaw) * HMM_CosF(rpitch),
		HMM_SinF(rpitch),
		HMM_CosF(ryaw) * HMM_CosF(rpitch)
	);
	const hmm_vec3 r = HMM_Vec3
	(
		HMM_SinF(ryaw - float(HMM_PI * 0.5f)),
		0.0f,
		HMM_CosF(ryaw - float(HMM_PI * 0.5f))
	);
	if (forward)
		*forward = f;
	if (right)
		*right = r;
	if (up)
		*up = HMM_Multiply(HMM_Cross(f, r), -1.0f);
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
	FirstPersonCamera() : position(HMM_Vec3(0.0f, 0.0f, 0.0f)), pitch(0.0f), yaw(0.0f) {}

	hmm_mat4 calculateViewMatrix()
	{
		hmm_vec3 forward, up;
		axisFromEulerAngles(pitch, yaw, &forward, nullptr, &up);
		const hmm_vec3 at = HMM_Add(position, forward);
		return HMM_LookAt(position, at, up);
	}

	void move(float deltaForward, float deltaRight)
	{
		hmm_vec3 forward, right;
		axisFromEulerAngles(pitch, yaw, &forward, &right, nullptr);
		const hmm_vec3 velocity = HMM_Add(HMM_Multiply(forward, deltaForward), HMM_Multiply(right, deltaRight));
		position = HMM_Add(position, velocity);
	}

	void rotate(float deltaX, float deltaY)
	{
		yaw = cleanAngle(yaw + deltaX);
		pitch = HMM_Clamp(-90.0f, pitch + deltaY, 90.0f);
	}

	hmm_vec3 position;
	float pitch;
	float yaw;
};

struct OrbitCamera
{
	OrbitCamera() : distance(32.0f), pitch(0.0f), yaw(0.0f) {}

	hmm_mat4 calculateViewMatrix()
	{
		hmm_vec3 forward;
		axisFromEulerAngles(pitch, yaw, &forward, nullptr, nullptr);
		const hmm_vec3 center = HMM_Multiply(s_model.centroid, s_model.scale);
		const hmm_vec3 eye = HMM_Add(HMM_Multiply(forward, -distance), center);
		const hmm_vec3 up = HMM_Vec3(0.0f, 1.0f, 0.0f);
		return HMM_LookAt(eye, center, up);
	}

	void rotate(float deltaX, float deltaY)
	{
		yaw = cleanAngle(yaw - deltaX);
		pitch = HMM_Clamp(-75.0f, pitch + deltaY, 75.0f);
	}

	void zoom(float delta)
	{
		distance = HMM_Clamp(0.1f, distance + delta, 500.0f);
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

static GLuint createShader(GLenum type, const char *source)
{
	GLuint shader = glCreateShader(type);
	if (shader == 0) {
		fprintf(stderr, "glCreateShader failed\n");
		return 0;
	}
	glShaderSource(shader, 1, &source, nullptr);
	glCompileShader(shader);
	GLint compiled;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &compiled);
	if (!compiled) {
		if (type == GL_VERTEX_SHADER)
			fprintf(stderr, "Vertex shader compile failed\n");
		else
			fprintf(stderr, "Fragment shader compile failed\n");
		GLint infoLen = 0;
		glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infoLen);
		if (infoLen) {
			std::vector<char> infoLog;
			infoLog.resize(infoLen);
			glGetShaderInfoLog(shader, infoLen, nullptr, infoLog.data());
			fprintf(stderr, "%s\n", infoLog.data());
		}
		glDeleteShader(shader);
		return 0;
	}
	return shader;
}

static GLuint createShaderProgram(const char *vp, const char *fp, const char **attributes, int attributeCount)
{
	GLuint vertexShader = createShader(GL_VERTEX_SHADER, vp);
	if (!vertexShader)
		return 0;
	GLuint fragmentShader = createShader(GL_FRAGMENT_SHADER, fp);
	if (!fragmentShader) {
		glDeleteShader(vertexShader);
		return 0;
	}
	GLuint program = glCreateProgram();
	if (program == 0) {
		fprintf(stderr, "glCreateProgram failed\n");
		return 0;
	}
	glAttachShader(program, vertexShader);
	glAttachShader(program, fragmentShader);
	for (int i = 0; i < attributeCount; i++)
		glBindAttribLocation(program, i, attributes[i]);
	glLinkProgram(program);
	glDeleteShader(vertexShader);
	glDeleteShader(fragmentShader);
	GLint linked;
	glGetProgramiv(program, GL_LINK_STATUS, &linked);
	if (!linked) {
		fprintf(stderr, "Linking shader program failed\n");
		GLint infoLen = 0;
		glGetProgramiv(program, GL_INFO_LOG_LENGTH, &infoLen);
		if (infoLen) {
			char* infoLog = (char*)malloc(sizeof(char) * infoLen);
			glGetProgramInfoLog(program, infoLen, nullptr, infoLog);
			fprintf(stderr, "%s\n", infoLog);
			free(infoLog);
		}
		glDeleteProgram(program);
		return 0;
	}
	return program;
}

static void shadersInit()
{
	// create color shader
	const char *colorVertex = R"(
#version 330 core
in vec3 a_position;
uniform mat4 u_mvp;

void main()
{
	gl_Position = u_mvp * vec4(a_position, 1.0);
}
	)";
	const char *colorFragment = R"(
#version 330 core
uniform vec4 u_color;
out vec4 o_color;

void main()
{
	o_color = u_color;
}
	)";
	const char *colorAttribs[] = {
		"a_position"
	};
	s_colorShader.id = createShaderProgram(colorVertex, colorFragment, colorAttribs, 1);
	if (!s_colorShader.id) {
		fprintf(stderr, "Error creating color shader\n");
		exit(EXIT_FAILURE);
	}
	s_colorShader.u_color = glGetUniformLocation(s_colorShader.id, "u_color");
	s_colorShader.u_mvp = glGetUniformLocation(s_colorShader.id, "u_mvp");
	// create texcoord shader
	const char *texcoordVertex = R"(
#version 330 core
in vec3 a_position;
in vec3 a_normal;
in vec4 a_texcoord;
out vec4 v_texcoord;
uniform mat4 u_mvp;

void main()
{
	v_texcoord = a_texcoord;
	gl_Position = u_mvp * vec4(a_position, 1.0);
}
	)";
	const char *texcoordFragment = R"(
#version 330 core
uniform vec4 u_color;
uniform vec4 u_textureSize_cellSize;
in vec4 v_texcoord;
out vec4 o_color;

void main()
{
	int x = int(v_texcoord.z * u_textureSize_cellSize.x);
	int y = int(v_texcoord.w * u_textureSize_cellSize.y);
	int cellSize = int(u_textureSize_cellSize.z);
	float scale = (x / cellSize % 2) != (y / cellSize % 2) ? 0.75 : 1.0;
	o_color = vec4(u_color.rgb * scale, u_color.a);
}
	)";
	const char *texcoordAttribs[] = {
		"a_position",
		"a_normal",
		"a_texcoord"
	};
	s_texcoordShader.id = createShaderProgram(texcoordVertex, texcoordFragment, texcoordAttribs, 3);
	if (!s_texcoordShader.id) {
		fprintf(stderr, "Error creating texcoord shader\n");
		exit(EXIT_FAILURE);
	}
	s_texcoordShader.u_color = glGetUniformLocation(s_texcoordShader.id, "u_color");
	s_texcoordShader.u_mvp = glGetUniformLocation(s_texcoordShader.id, "u_mvp");
	s_texcoordShader.u_textureSize_cellSize = glGetUniformLocation(s_texcoordShader.id, "u_textureSize_cellSize");
	// create gui shader
	const char *guiVertex = R"(
#version 330 core
uniform mat4 u_mvp;
in vec2 a_position;
in vec2 a_texcoord;
in vec4 a_color;
out vec2 v_texcoord;
out vec4 v_color;

void main()
{
	v_texcoord = a_texcoord;
	v_color = a_color;
	gl_Position = u_mvp * vec4(a_position.xy, 0, 1);
}
	)";

	const char *guiFragment = R"(
#version 330 core
uniform sampler2D u_texture;
in vec2 v_texcoord;
in vec4 v_color;
out vec4 o_color;

void main()
{
	o_color = v_color * texture(u_texture, v_texcoord.st);
}
	)";
	const char *guiAttribs[] =
	{
		"a_position",
		"a_texcoord",
		"a_color"
	};
	s_gui.shaderProgram = createShaderProgram(guiVertex, guiFragment, guiAttribs, 3);
	if (!s_gui.shaderProgram) {
		fprintf(stderr, "Error creating GUI shader");
		exit(EXIT_FAILURE);
	}
	s_gui.u_mvp = glGetUniformLocation(s_gui.shaderProgram, "u_mvp");
	s_gui.u_texture = glGetUniformLocation(s_gui.shaderProgram, "u_texture");
}

static void shadersShutdown()
{
	glDeleteProgram(s_colorShader.id);
	glDeleteProgram(s_texcoordShader.id);
	glDeleteProgram(s_gui.shaderProgram);
}

static void guiInit()
{
	ImGui::CreateContext();
	int w, h;
	glfwGetFramebufferSize(s_window, &w, &h);
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
	glGenTextures(1, &s_gui.fontTexture);
	glBindTexture(GL_TEXTURE_2D, s_gui.fontTexture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, fontWidth, fontHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, fontData);
	io.Fonts->TexID = (ImTextureID)(size_t)s_gui.fontTexture;
	glBindTexture(GL_TEXTURE_2D, 0);
	// vertex buffer
	glGenBuffers(1, &s_gui.vbo);
	glGenBuffers(1, &s_gui.ibo);
	glGenVertexArrays(1, &s_gui.vao);
	glBindVertexArray(s_gui.vao);
	glBindBuffer(GL_ARRAY_BUFFER, s_gui.vbo);
	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	glEnableVertexAttribArray(2);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(ImDrawVert), (void *)offsetof(ImDrawVert, pos));
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(ImDrawVert), (void *)offsetof(ImDrawVert, uv));
	glVertexAttribPointer(2, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(ImDrawVert), (void *)offsetof(ImDrawVert, col));
	glBindVertexArray(0);
}

static void guiResize(int width, int height)
{
	ImGuiIO &io = ImGui::GetIO();
	io.DisplaySize.x = (float)width;
	io.DisplaySize.y = (float)height;
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
	// Setup render state: alpha-blending enabled, no face culling, no depth testing, scissor enabled
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glDisable(GL_CULL_FACE);
	glDisable(GL_DEPTH_TEST);
	glEnable(GL_SCISSOR_TEST);
	glActiveTexture(GL_TEXTURE0);
	// Setup viewport, orthographic projection matrix
	ImGuiIO &io = ImGui::GetIO();
	glViewport(0, 0, (int)io.DisplaySize.x, (int)io.DisplaySize.y);
	const hmm_mat4 projection = HMM_Orthographic(0.0f, io.DisplaySize.x, io.DisplaySize.y, 0.0f, 0.0f, 1.0f);
	glUseProgram(s_gui.shaderProgram);
	glUniform1i(s_gui.u_texture, 0);
	glUniformMatrix4fv(s_gui.u_mvp, 1, GL_FALSE, (const float *)&projection);
	glBindVertexArray(s_gui.vao);
	ImDrawData *draw_data = ImGui::GetDrawData();
	for (int n = 0; n < draw_data->CmdListsCount; n++) {
		const ImDrawList* cmd_list = draw_data->CmdLists[n];
		const ImDrawIdx* idx_buffer_offset = 0;
		glBindBuffer(GL_ARRAY_BUFFER, s_gui.vbo);
		glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)cmd_list->VtxBuffer.Size * sizeof(ImDrawVert), (const GLvoid*)cmd_list->VtxBuffer.Data, GL_STREAM_DRAW);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, s_gui.ibo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, (GLsizeiptr)cmd_list->IdxBuffer.Size * sizeof(ImDrawIdx), (const GLvoid*)cmd_list->IdxBuffer.Data, GL_STREAM_DRAW);
		for (int cmd_i = 0; cmd_i < cmd_list->CmdBuffer.Size; cmd_i++) {
			const ImDrawCmd* pcmd = &cmd_list->CmdBuffer[cmd_i];
			if (pcmd->UserCallback)
				pcmd->UserCallback(cmd_list, pcmd);
			else {
				glBindTexture(GL_TEXTURE_2D, (GLuint)(intptr_t)pcmd->TextureId);
				glScissor((int)pcmd->ClipRect.x, (int)(io.DisplaySize.y - pcmd->ClipRect.w), (int)(pcmd->ClipRect.z - pcmd->ClipRect.x), (int)(pcmd->ClipRect.w - pcmd->ClipRect.y));
				glDrawElements(GL_TRIANGLES, (GLsizei)pcmd->ElemCount, sizeof(ImDrawIdx) == 2 ? GL_UNSIGNED_SHORT : GL_UNSIGNED_INT, idx_buffer_offset);
			}
			idx_buffer_offset += pcmd->ElemCount;
		}
	}
	glBindTexture(GL_TEXTURE_2D, 0);
	glBindVertexArray(0);
	glDisable(GL_BLEND);
	glDisable(GL_SCISSOR_TEST);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
}

static void guiShutdown()
{
	ImGui::DestroyContext();
	glDeleteTextures(1, &s_gui.fontTexture);
}


static void atlasDestroy()
{
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
		glDeleteTextures(1, &s_atlas.chartsTextures[i]);
		s_atlas.chartsTextures[i] = 0;
	}
	if (s_atlas.chartVao > 0) {
		glDeleteVertexArrays(1, &s_atlas.chartVao);
		glDeleteBuffers(1, &s_atlas.chartVbo);
		glDeleteBuffers(1, &s_atlas.chartIbo);
		s_atlas.chartVao = 0;
	}
	if (s_atlas.chartBoundaryVao > 0) {
		glDeleteVertexArrays(1, &s_atlas.chartBoundaryVao);
		glDeleteBuffers(1, &s_atlas.chartBoundaryVbo);
		s_atlas.chartBoundaryVao = 0;
	}
	s_atlas.status.set(AtlasStatus::NotGenerated);
}

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
	if (s_model.vao > 0) {
		glDeleteVertexArrays(1, &s_model.vao);
		glDeleteBuffers(1, &s_model.vbo);
		glDeleteBuffers(1, &s_model.ibo);
		s_model.vao = 0;
	}
	glfwSetWindowTitle(s_window, WINDOW_TITLE);
	s_model.status.set(ModelStatus::NotLoaded);
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
		v->texcoord.Y = 1.0f - v->texcoord.Y;
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
	s_model.centroid = HMM_Vec3(0.0f, 0.0f, 0.0f);
	for (uint32_t i = 0; i < s_model.data->numVertices; i++)
		s_model.centroid = HMM_Add(s_model.centroid, ((const ModelVertex *)s_model.data->vertices)[i].pos);
	s_model.centroid = HMM_Multiply(s_model.centroid, 1.0f / s_model.data->numVertices);
	glGenBuffers(1, &s_model.vbo);
	glGenBuffers(1, &s_model.ibo);
	glGenVertexArrays(1, &s_model.vao);
	glBindVertexArray(s_model.vao);
	glBindBuffer(GL_ARRAY_BUFFER, s_model.vbo);
	glBufferData(GL_ARRAY_BUFFER, s_model.data->numVertices * sizeof(ModelVertex), s_model.data->vertices, GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, s_model.ibo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, s_model.data->numIndices * sizeof(uint32_t), s_model.data->indices, GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	glEnableVertexAttribArray(2);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(ModelVertex), (void *)offsetof(ModelVertex, pos));
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(ModelVertex), (void *)offsetof(ModelVertex, normal));
	glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(ModelVertex), (void *)offsetof(ModelVertex, texcoord));
	glBindVertexArray(0);
	s_camera.firstPerson = FirstPersonCamera();
	s_camera.orbit = OrbitCamera();
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
	ModelLoadThreadArgs args;
	STRNCPY(args.filename, sizeof(args.filename), filename);
	s_model.thread = new std::thread(modelLoadThread, args);
	free(filename);
}

static void modelRender(const hmm_mat4 &view, const hmm_mat4 &projection)
{
	if (s_model.status.get() != ModelStatus::Ready)
		return;
	const hmm_mat4 model = HMM_Scale(HMM_Vec3(s_model.scale, s_model.scale, s_model.scale));
	const hmm_mat4 mvp = HMM_Multiply(projection, HMM_Multiply(view, model));
	if (s_atlas.status.get() == AtlasStatus::Ready) {
		glBindVertexArray(s_atlas.chartVao);
		glUseProgram(s_texcoordShader.id);
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
				glUniform4fv(s_texcoordShader.u_color, 1, color);
				float textureSize_cellSize[4];
				textureSize_cellSize[0] = (float)s_atlas.data->width;
				textureSize_cellSize[1] = (float)s_atlas.data->height;
				textureSize_cellSize[2] = (float)s_atlas.chartCellSize;
				textureSize_cellSize[3] = (float)s_atlas.chartCellSize;
				glUniform4fv(s_texcoordShader.u_textureSize_cellSize, 1, textureSize_cellSize);
				glUniformMatrix4fv(s_texcoordShader.u_mvp, 1, GL_FALSE, (const float *)&mvp);
				glDrawElements(GL_TRIANGLES, chart.indexCount, GL_UNSIGNED_INT, (void *)(firstIndex * sizeof(uint32_t)));
				firstIndex += chart.indexCount;
			}
		}
	} else {
		glBindVertexArray(s_model.vao);
		glUseProgram(s_colorShader.id);
		const float color[] = { 0.75f, 0.75f, 0.75f, 1.0f };
		glUniform4fv(s_colorShader.u_color, 1, color);
		glUniformMatrix4fv(s_colorShader.u_mvp, 1, GL_FALSE, (const float *)&mvp);
		glDrawElements(GL_TRIANGLES, s_model.data->numIndices, GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
	}
	if (s_atlas.status.get() == AtlasStatus::Ready && s_options.wireframe && s_options.wireframeMode == WireframeMode::Charts) {
		// Chart boundary edges.
		glDisable(GL_DEPTH_TEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glBindVertexArray(s_atlas.chartBoundaryVao);
		glUseProgram(s_colorShader.id);
		const float wcolor[] = { 1.0f, 1.0f, 1.0f, 0.5f };
		glUniform4fv(s_colorShader.u_color, 1, wcolor);
		glUniformMatrix4fv(s_colorShader.u_mvp, 1, GL_FALSE, (const float *)&mvp);
		glDrawArrays(GL_LINES, 0, (GLsizei)s_atlas.chartBoundaryVertices.size());
		glBindVertexArray(0);
		glDisable(GL_BLEND);
		glEnable(GL_DEPTH_TEST);
	}
	else if (s_options.wireframe) {
		glBindVertexArray(s_model.vao);
		glUseProgram(s_colorShader.id);
		glDisable(GL_DEPTH_TEST);
		glDisable(GL_CULL_FACE);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		const float wcolor[] = { 1.0f, 1.0f, 1.0f, 0.5f };
		glUniform4fv(s_colorShader.u_color, 1, wcolor);
		glUniformMatrix4fv(s_colorShader.u_mvp, 1, GL_FALSE, (const float *)&mvp);
		glDrawElements(GL_TRIANGLES, s_model.data->numIndices, GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_CULL_FACE);
		glDisable(GL_BLEND);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
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

// https://github.com/ssloy/tinyrenderer/wiki/Lesson-2:-Triangle-rasterization-and-back-face-culling
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
	hmm_vec3 p0, p1;
};

struct EdgeKeyHash
{
	uint32_t operator()(const EdgeKey &key) const
	{
		int32_t data[6];
		data[0] = (int32_t)(key.p0.X * 100.0f);
		data[1] = (int32_t)(key.p0.Y * 100.0f);
		data[2] = (int32_t)(key.p0.Z * 100.0f);
		data[3] = (int32_t)(key.p1.X * 100.0f);
		data[4] = (int32_t)(key.p1.Y * 100.0f);
		data[5] = (int32_t)(key.p1.Z * 100.0f);
		return sdbmHash(data, sizeof(data));
	}
};

struct EdgeKeyEqual
{
	bool floatEqual(float f1, float f2) const
	{
		const float epsilon = 0.0001f;
		return fabs(f1 - f2) <= epsilon;
	}

	bool vec3Equal(const hmm_vec3 &v0, const hmm_vec3 &v1) const
	{
		return floatEqual(v0.X, v1.X) && floatEqual(v0.Y, v1.Y) && floatEqual(v0.Z, v1.Z);
	}

	bool operator()(const EdgeKey &k0, const EdgeKey &k1) const
	{
		return vec3Equal(k0.p0, k1.p0) && vec3Equal(k0.p1, k1.p1);
	}
};

static void atlasGenerateThread()
{
	int progress = 0;
	if (!s_atlas.data) {
		// Create xatlas context and generate charts on first run only.
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
		xatlas::GenerateCharts(s_atlas.data, s_atlas.charterOptions, atlasProgressCallback);
	} else
		s_atlas.status.set(AtlasStatus::Generating);
	xatlas::PackCharts(s_atlas.data, s_atlas.packerOptions, atlasProgressCallback);
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
	s_atlas.chartIndices.clear();
	s_atlas.chartVertices.clear();
	s_atlas.chartBoundaryVertices.clear();
	uint32_t numIndices = 0, numVertices = 0;
	for (uint32_t i = 0; i < s_atlas.data->meshCount; i++) {
		const xatlas::Mesh &outputMesh = s_atlas.data->meshes[i];
		numIndices += outputMesh.indexCount;
		numVertices += outputMesh.vertexCount;
	}
	s_atlas.chartIndices.resize(numIndices);
	s_atlas.chartVertices.resize(numVertices);
	uint32_t firstIndex = 0;
	uint32_t firstVertex = 0;
	numEdges = 0;
	for (uint32_t i = 0; i < s_atlas.data->meshCount; i++) {
		const xatlas::Mesh &mesh = s_atlas.data->meshes[i];
		const objzObject &object = s_model.data->objects[i];
		const ModelVertex *oldVertices = &((const ModelVertex *)s_model.data->vertices)[object.firstVertex];
		for (uint32_t j = 0; j < mesh.chartCount; j++) {
			const xatlas::Chart &chart = mesh.chartArray[j];
			for (uint32_t k = 0; k < chart.indexCount; k++)
				s_atlas.chartIndices[firstIndex + k] = firstVertex + chart.indexArray[k];
			firstIndex += chart.indexCount;
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
			ModelVertex &v = s_atlas.chartVertices[firstVertex + j];
			v.pos = oldVertex.pos;
			v.normal = oldVertex.normal;
			v.texcoord = HMM_Vec4(oldVertex.texcoord.X, oldVertex.texcoord.Y, outputVertex.uv[0] / (float)s_atlas.data->width, outputVertex.uv[1] / (float)s_atlas.data->height);
		}
		firstVertex += mesh.vertexCount;
	}
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
	xatlas::SetPrint(s_atlas.verbose ? printf : nullptr);
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
	glGenBuffers(1, &s_atlas.chartVbo);
	glGenBuffers(1, &s_atlas.chartIbo);
	glGenVertexArrays(1, &s_atlas.chartVao);
	glBindVertexArray(s_atlas.chartVao);
	glBindBuffer(GL_ARRAY_BUFFER, s_atlas.chartVbo);
	glBufferData(GL_ARRAY_BUFFER, s_atlas.chartVertices.size() * sizeof(ModelVertex), s_atlas.chartVertices.data(), GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, s_atlas.chartIbo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, s_atlas.chartIndices.size() * sizeof(uint32_t), s_atlas.chartIndices.data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	glEnableVertexAttribArray(2);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(ModelVertex), (void *)offsetof(ModelVertex, pos));
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(ModelVertex), (void *)offsetof(ModelVertex, normal));
	glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(ModelVertex), (void *)offsetof(ModelVertex, texcoord));
	glBindVertexArray(0);
	// Chart boundaries.
	glGenVertexArrays(1, &s_atlas.chartBoundaryVao);
	glBindVertexArray(s_atlas.chartBoundaryVao);
	glGenBuffers(1, &s_atlas.chartBoundaryVbo);
	glBindBuffer(GL_ARRAY_BUFFER, s_atlas.chartBoundaryVbo);
	glBufferData(GL_ARRAY_BUFFER, s_atlas.chartBoundaryVertices.size() * sizeof(hmm_vec3), s_atlas.chartBoundaryVertices.data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(hmm_vec3), 0);
	glBindVertexArray(0);
	// Charts texture.
	s_atlas.chartsTextures.resize(s_atlas.data->atlasCount);
	for (uint32_t i = 0; i < (uint32_t)s_atlas.chartsTextures.size(); i++) {
		if (s_atlas.chartsTextures[i] == 0)
			glGenTextures(1, &s_atlas.chartsTextures[i]);
		glBindTexture(GL_TEXTURE_2D, s_atlas.chartsTextures[i]);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
		const float color[] = { 0.0f, 0.0f, 0.0f, 0.0f };
		glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, color);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, s_atlas.data->width, s_atlas.data->height, 0, GL_RGB, GL_UNSIGNED_BYTE, s_atlas.chartsImages[i].data());
		glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
	}
	s_atlas.currentTexture = 0;
	s_options.wireframeMode = WireframeMode::Charts;
	s_atlas.status.set(AtlasStatus::Ready);
}

int main(int /*argc*/, char ** /*argv*/)
{
	glfwSetErrorCallback(glfw_errorCallback);
	if (!glfwInit())
		return EXIT_FAILURE;
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_SAMPLES, 0);
	s_window = glfwCreateWindow(WINDOW_DEFAULT_WIDTH, WINDOW_DEFAULT_HEIGHT, WINDOW_TITLE, nullptr, nullptr);
	if (!s_window)
		return EXIT_FAILURE;
	glfwMaximizeWindow(s_window);
	glfwMakeContextCurrent(s_window);
	glfwSwapInterval(1);
	flextInit();
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_CULL_FACE);
	shadersInit();
	guiInit();
	glfwSetCharCallback(s_window, glfw_charCallback);
	glfwSetCursorPosCallback(s_window, glfw_cursorPosCallback);
	glfwSetKeyCallback(s_window, glfw_keyCallback);
	glfwSetMouseButtonCallback(s_window, glfw_mouseButtonCallback);
	glfwSetScrollCallback(s_window, glfw_scrollCallback);
	int frameCount = 0, progressDots = 0;
	double lastFrameTime = glfwGetTime();
	while (!glfwWindowShouldClose(s_window)) {
		glfwPollEvents();
		if (glfwGetKey(s_window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
			glfwSetWindowShouldClose(s_window, GLFW_TRUE);
		const double currentFrameTime = glfwGetTime();
		const float deltaTime = float(currentFrameTime - lastFrameTime);
		lastFrameTime = currentFrameTime;
		int width, height;
		glfwGetFramebufferSize(s_window, &width, &height);
		guiResize(width, height);
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
			const float margin = 8.0f;
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
				if (s_model.data) {
					ImGui::Text("%u objects", s_model.data->numObjects);
					ImGui::Text("%u vertices", s_model.data->numVertices);
					ImGui::Text("%u triangles", s_model.data->numIndices / 3);
				}
				ImGui::InputFloat("Model scale", &s_model.scale, 0.01f, 0.1f);
				s_model.scale = HMM_MAX(0.001f, s_model.scale);
				ImGui::Spacing();
				ImGui::Separator();
				ImGui::Spacing();
				ImGui::Text("View");
				ImGui::Spacing();
				ImGui::Checkbox("Wireframe", &s_options.wireframe);
				if (s_atlas.status.get() == AtlasStatus::Ready) {
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
				ImGui::DragFloat("FOV", &s_camera.fov, 1.0f, 45.0f, 150.0f, "%.0f");
				ImGui::DragFloat("Sensitivity", &s_camera.sensitivity, 0.01f, 0.01f, 1.0f);
				if (s_model.status.get() == ModelStatus::Ready) {
					ImGui::Spacing();
					ImGui::Separator();
					ImGui::Spacing();
					ImGui::Text("Charter options");
					ImGui::Spacing();
					ImGui::InputFloat("Proxy fit metric weight", &s_atlas.charterOptions.proxyFitMetricWeight);
					ImGui::InputFloat("Roundness metric weight", &s_atlas.charterOptions.roundnessMetricWeight);
					ImGui::InputFloat("Straightness metric weight", &s_atlas.charterOptions.straightnessMetricWeight);
					ImGui::InputFloat("NormalSeam metric weight", &s_atlas.charterOptions.normalSeamMetricWeight);
					ImGui::InputFloat("Texture seam metric weight", &s_atlas.charterOptions.textureSeamMetricWeight);
					ImGui::InputFloat("Max chart area", &s_atlas.charterOptions.maxChartArea);
					ImGui::InputFloat("Max boundary length", &s_atlas.charterOptions.maxBoundaryLength);
					ImGui::InputFloat("Max threshold", &s_atlas.charterOptions.maxThreshold);
					ImGui::InputInt("Grow face count", (int *)&s_atlas.charterOptions.growFaceCount);
					ImGui::InputInt("Max iterations", (int *)&s_atlas.charterOptions.maxIterations);
					ImGui::Spacing();
					ImGui::Separator();
					ImGui::Spacing();
					ImGui::Text("Packer options");
					ImGui::Spacing();
					ImGui::SliderInt("Attempts", &s_atlas.packerOptions.attempts, 0, 4096);
					ImGui::InputFloat("Texels per unit", &s_atlas.packerOptions.texelsPerUnit, 0.0f, 32.0f, 2);
					ImGui::InputInt("Resolution", (int *)&s_atlas.packerOptions.resolution, 8);
					ImGui::InputInt("Max chart size", (int *)&s_atlas.packerOptions.maxChartSize);
					ImGui::Checkbox("Block align", &s_atlas.packerOptions.blockAlign);
					ImGui::SameLine();
					ImGui::Checkbox("Conservative", &s_atlas.packerOptions.conservative);
					ImGui::SliderInt("Padding", &s_atlas.packerOptions.padding, 0, 8);
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
						ImGui::Text("%u charts", s_atlas.data->chartCount);
						ImGui::Text("%u vertices", numVertices);
						ImGui::Text("%u triangles", numIndices / 3);
						ImGui::Checkbox("Show atlas", &s_atlas.showTexture);
					}
				}
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
					ImTextureID texture = (ImTextureID)(size_t)s_atlas.chartsTextures[s_atlas.currentTexture];
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
			}
		}
		if (s_camera.mode == CameraMode::FirstPerson && glfwGetInputMode(s_window, GLFW_CURSOR) == GLFW_CURSOR_DISABLED) {
			const float speed = (s_keyDown[GLFW_KEY_LEFT_SHIFT] ? 20.0f : 5.0f) * deltaTime;
			float deltaForward = 0.0f, deltaRight = 0.0f;
			if (s_keyDown[GLFW_KEY_W]) deltaForward += speed;
			if (s_keyDown[GLFW_KEY_S]) deltaForward -= speed;
			if (s_keyDown[GLFW_KEY_A]) deltaRight -= speed;
			if (s_keyDown[GLFW_KEY_D]) deltaRight += speed;
			s_camera.firstPerson.move(deltaForward, deltaRight);
			if (s_keyDown[GLFW_KEY_Q]) s_camera.firstPerson.position.Y -= speed;
			if (s_keyDown[GLFW_KEY_E]) s_camera.firstPerson.position.Y += speed;
		}
		glClearColor(0.25f, 0.25f, 0.25f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glViewport(0, 0, width, height);
		if (s_model.status.get() == ModelStatus::Ready) {
			hmm_mat4 view;
			if (s_camera.mode == CameraMode::FirstPerson)
				view = s_camera.firstPerson.calculateViewMatrix();
			else if (s_camera.mode == CameraMode::Orbit)
				view = s_camera.orbit.calculateViewMatrix();
			const hmm_mat4 projection = HMM_Perspective(s_camera.fov, width / (float)height, 0.01f, 1000.0f);
			modelRender(view, projection);
		}
		if (s_options.gui)
			guiRender();
		glfwSwapBuffers(s_window);
		frameCount++;
		if (frameCount % 20 == 0)
			progressDots = (progressDots + 1) % 4;
		if (s_model.status.get() == ModelStatus::Finalizing)
			modelFinalize();
		if (s_atlas.status.get() == AtlasStatus::Finalizing)
			atlasFinalize();
	}
	guiShutdown();
	atlasDestroy();
	modelDestroy();
	shadersShutdown();
	glfwTerminate();
	return 0;
}
