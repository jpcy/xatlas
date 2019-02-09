/*
Copyright (c) 2018 Jonathan Young
Copyright (c) 2013 Thekla, Inc
Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#include <mutex>
#include <thread>
#include <vector>
#include <stdio.h>
#include "flextGL.h"
#include "GLFW/glfw3.h"

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#include "HandmadeMath.h"
#pragma GCC diagnostic pop
#else
#include "HandmadeMath.h"
#endif

#include "imgui/imgui.h"
#include "nativefiledialog/nfd.h"
#include "objzero/objzero.h"

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
		Ready,
		Error
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
	hmm_vec2 texcoord;
};

struct
{
	ModelStatus status;
	std::thread *thread = nullptr;
	objzModel *data;
	hmm_vec3 centroid = HMM_Vec3(0.0f, 0.0f, 0.0f);
	GLuint vao, vbo, ibo;
	float scale = 1.0f;
}
s_model;

struct
{
	bool gui = true;
	hmm_vec3 clearColor;
}
s_options;

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

	hmm_vec3 position = HMM_Vec3(0.0f, 0.0f, 0.0f);
	float pitch = 0.0f;
	float yaw = 0.0f;
};

struct OrbitCamera
{
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

	float distance = 10.0f;
	float pitch = 0.0f;
	float yaw = 0.0f;
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
	float fov = 75.0f;
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
	if (!s_options.gui)
		return;
	ImGuiIO &io = ImGui::GetIO();
	if (io.WantCaptureMouse)
		io.MouseWheel += (float)yoffset;
	else {
		if (s_camera.mode == CameraMode::Orbit) 
			s_camera.orbit.zoom((float)-yoffset);
	}
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

static void modelDestroy()
{
	if (s_model.thread) {
		s_model.thread->join();
		delete s_model.thread;
		s_model.thread = nullptr;
	}
	if (s_model.data) {
		objz_destroy(s_model.data);
		s_model.data = nullptr;
	}
	if (s_model.vbo > 0) {
		glDeleteBuffers(1, &s_model.vbo);
		s_model.vbo = 0;
	}
	if (s_model.ibo > 0) {
		glDeleteBuffers(1, &s_model.ibo);
		s_model.ibo = 0;
	}
	if (s_model.vao > 0) {
		glDeleteVertexArrays(1, &s_model.vao);
		s_model.vao = 0;
	}
	glfwSetWindowTitle(s_window, WINDOW_TITLE);
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
		s_model.data = nullptr;
		s_model.status.set(ModelStatus::Error);
		return;
	}
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
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, s_model.ibo);
	glBufferData(GL_ARRAY_BUFFER, (int)s_model.data->numVertices * sizeof(ModelVertex), s_model.data->vertices, GL_STATIC_DRAW);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, (int)s_model.data->numIndices * sizeof(uint32_t), s_model.data->indices, GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	glEnableVertexAttribArray(2);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(ModelVertex), (void *)offsetof(ModelVertex, pos));
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(ModelVertex), (void *)offsetof(ModelVertex, normal));
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(ModelVertex), (void *)offsetof(ModelVertex, texcoord));
	glBindVertexArray(0);
	s_model.status.set(ModelStatus::Ready);
}

static void modelLoad(const char *filename)
{
	modelDestroy();
	s_model.status.set(ModelStatus::Loading);
	printf("Loading '%s'\n", filename);
	char windowTitle[256];
	snprintf(windowTitle, sizeof(windowTitle), "%s - %s\n", WINDOW_TITLE, filename);
	glfwSetWindowTitle(s_window, windowTitle);
	ModelLoadThreadArgs args;
	STRNCPY(args.filename, sizeof(args.filename), filename);
	s_model.thread = new std::thread(modelLoadThread, args);
}

static void modelOpenDialog()
{
	nfdchar_t *nfdPath = nullptr;
	nfdresult_t result = NFD_OpenDialog("obj", nullptr, &nfdPath);
	if (result != NFD_OKAY)
		return;
	modelLoad(nfdPath);
	free(nfdPath);
}

static void modelRender(const hmm_mat4 &view, const hmm_mat4 &projection)
{
	if (s_model.status.get() != ModelStatus::Ready)
		return;
	const float scale = HMM_MAX(s_model.scale, 0.001f);
	const hmm_mat4 model = HMM_Scale(HMM_Vec3(scale, scale, scale));
	const hmm_mat4 mvp = HMM_Multiply(projection, HMM_Multiply(view, model));
	glDisable(GL_CULL_FACE);
	glDisable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glBindVertexArray(s_model.vao);
	glUseProgram(s_colorShader.id);
	const float color[] = { 1, 1, 1, 0.5f };
	glUniform4fv(s_colorShader.u_color, 1, color);
	glUniformMatrix4fv(s_colorShader.u_mvp, 1, GL_FALSE, (const float *)&mvp);
	glDrawElements(GL_TRIANGLES, s_model.data->numIndices, GL_UNSIGNED_INT, 0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glDisable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
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
	glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
	s_window = glfwCreateWindow(WINDOW_DEFAULT_WIDTH, WINDOW_DEFAULT_HEIGHT, WINDOW_TITLE, nullptr, nullptr);
	if (!s_window)
		return EXIT_FAILURE;
	glfwMakeContextCurrent(s_window);
	glfwSwapInterval(1);
	flextInit();
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_CULL_FACE);
	s_options.clearColor.X = s_options.clearColor.Y = s_options.clearColor.Z = 0.5f;
	shadersInit();
	guiInit();
	glfwSetCharCallback(s_window, glfw_charCallback);
	glfwSetCursorPosCallback(s_window, glfw_cursorPosCallback);
	glfwSetKeyCallback(s_window, glfw_keyCallback);
	glfwSetMouseButtonCallback(s_window, glfw_mouseButtonCallback);
	glfwSetScrollCallback(s_window, glfw_scrollCallback);
	glfwMaximizeWindow(s_window);
	glfwShowWindow(s_window);
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
			const float margin = 8.0f;
			ImGui::SetNextWindowPos(ImVec2(margin, margin), ImGuiCond_FirstUseEver);
			ImGui::SetNextWindowSize(ImVec2(350.0f, 400.0f), ImGuiCond_FirstUseEver);
			if (ImGui::Begin("##mainWindow", nullptr, ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse)) {
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
				ImGui::ColorEdit3("Clear color", &s_options.clearColor.X);
				ImGui::RadioButton("First person camera", (int *)&s_camera.mode, (int)CameraMode::FirstPerson);
				ImGui::SameLine();
				ImGui::RadioButton("Orbit camera", (int *)&s_camera.mode, (int)CameraMode::Orbit);
				ImGui::DragFloat("FOV", &s_camera.fov, 1.0f, 45.0f, 150.0f, "%.0f");
				ImGui::DragFloat("Sensitivity", &s_camera.sensitivity, 0.01f, 0.01f, 1.0f);
				ImGui::End();
			}
			if (s_model.status.get() == ModelStatus::Loading) {
				ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y * 0.5f), ImGuiCond_Always, ImVec2(0.5f, 0.5f));
				if (ImGui::Begin("##progress", nullptr, ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings)) {
					ImGui::Text("Loading");
					for (int i = 0; i < 3; i++) {
						ImGui::SameLine();
						ImGui::Text(i < progressDots ? "." : " ");
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
		glClearColor(s_options.clearColor.X, s_options.clearColor.Y, s_options.clearColor.Z, 1.0f);
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
	}
	guiShutdown();
	modelDestroy();
	shadersShutdown();
	glfwTerminate();
	return 0;
}
