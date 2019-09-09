/*
xatlas
https://github.com/jpcy/xatlas
Copyright (c) 2018 Jonathan Young

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#include <GLFW/glfw3.h>
#include <imgui/imgui.h>
#include <fonts/fontawesome-webfont.h>
#include <fonts/Roboto-Regular.h>
#include "viewer.h"

struct
{
	GLFWcursor *cursors[ImGuiMouseCursor_COUNT];
	GLFWcursor *currentCursor = nullptr;
	bgfx::VertexLayout vertexFormat;
	bgfx::TextureHandle font;
	bgfx::ProgramHandle program;
	bgfx::UniformHandle s_texture;
}
s_gui;

void guiInit()
{
	bgfx::setViewMode(kGuiView, bgfx::ViewMode::Sequential);
	bgfx::setViewRect(kGuiView, 0, 0, bgfx::BackbufferRatio::Equal);
	ImGui::CreateContext();
	int w, h;
	glfwGetWindowSize(g_window, &w, &h);
	ImGuiIO &io = ImGui::GetIO();
	io.BackendFlags |= ImGuiBackendFlags_HasMouseCursors;
	io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
	io.ConfigWindowsResizeFromEdges = true;
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
	// cursors
	for (uint32_t i = 0; i < ImGuiMouseCursor_COUNT; i++)
		s_gui.cursors[i] = nullptr;
	s_gui.cursors[ImGuiMouseCursor_Arrow] = glfwCreateStandardCursor(GLFW_ARROW_CURSOR);
	s_gui.cursors[ImGuiMouseCursor_TextInput] = glfwCreateStandardCursor(GLFW_IBEAM_CURSOR);
	s_gui.cursors[ImGuiMouseCursor_ResizeNS] = glfwCreateStandardCursor(GLFW_VRESIZE_CURSOR);
	s_gui.cursors[ImGuiMouseCursor_ResizeEW] = glfwCreateStandardCursor(GLFW_HRESIZE_CURSOR);
	s_gui.cursors[ImGuiMouseCursor_Hand] = glfwCreateStandardCursor(GLFW_HAND_CURSOR);
	// merge in icons from Font Awesome
	const float fontSize = 16.0f;
	io.Fonts->AddFontFromMemoryCompressedTTF(s_robotoRegular_compressed_data, s_robotoRegular_compressed_size, fontSize);
	static const ImWchar icons_ranges[] = { ICON_MIN_FA, ICON_MAX_FA, 0 };
	ImFontConfig icons_config;
	icons_config.MergeMode = true;
	icons_config.PixelSnapH = true;
	io.Fonts->AddFontFromMemoryCompressedTTF(s_fontAwesome_compressed_data, s_fontAwesome_compressed_size, fontSize, &icons_config, icons_ranges);
	// font
	int fontWidth, fontHeight;
	uint8_t *fontData;
	io.Fonts->GetTexDataAsRGBA32(&fontData, &fontWidth, &fontHeight);
	s_gui.font = bgfx::createTexture2D((uint16_t)fontWidth, (uint16_t)fontHeight, false, 0, bgfx::TextureFormat::RGBA8, 0, bgfx::makeRef(fontData, fontWidth * fontHeight * 4));
	io.Fonts->TexID = (ImTextureID)(intptr_t)s_gui.font.idx;
	// ImDrawVert vertex decl
	s_gui.vertexFormat
		.begin()
		.add(bgfx::Attrib::Position, 2, bgfx::AttribType::Float)
		.add(bgfx::Attrib::TexCoord0, 2, bgfx::AttribType::Float)
		.add(bgfx::Attrib::Color0, 4, bgfx::AttribType::Uint8, true)
		.end();
	// shader program
	s_gui.s_texture = bgfx::createUniform("s_texture", bgfx::UniformType::Sampler);
	bgfx::ShaderHandle vertex = loadShader(ShaderId::vs_gui);
	bgfx::ShaderHandle fragment = loadShader(ShaderId::fs_gui);
	s_gui.program = bgfx::createProgram(vertex, fragment, true);
}

void guiShutdown()
{
	ImGui::DestroyContext();
	bgfx::destroy(s_gui.font);
	bgfx::destroy(s_gui.s_texture);
	bgfx::destroy(s_gui.program);
}

void guiResize(int width, int height)
{
	ImGuiIO &io = ImGui::GetIO();
	io.DisplaySize.x = (float)width;
	io.DisplaySize.y = (float)height;
	bgfx::setViewRect(kGuiView, 0, 0, bgfx::BackbufferRatio::Equal);
}

void guiRunFrame(float deltaTime)
{
	ImGuiIO &io = ImGui::GetIO();
	io.DeltaTime = deltaTime;
	if (glfwGetInputMode(g_window, GLFW_CURSOR) == GLFW_CURSOR_NORMAL) {
		double x, y;
		glfwGetCursorPos(g_window, &x, &y);
		io.MousePos.x = (float)x;
		io.MousePos.y = (float)y;
	}
	else
		io.MousePos.x = io.MousePos.y = -1.0f;
}

void guiRender()
{
	GLFWcursor *cursor = s_gui.cursors[ImGui::GetMouseCursor()];
	if (cursor != s_gui.currentCursor) {
		glfwSetCursor(g_window, cursor ? cursor : s_gui.cursors[ImGuiMouseCursor_Arrow]);
		s_gui.currentCursor = cursor;
	}
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
		if (!bgfx::allocTransientBuffers(&tvb, s_gui.vertexFormat, cmd_list->VtxBuffer.Size, &tib, cmd_list->IdxBuffer.Size))
			return;
		assert(sizeof(ImDrawVert) == s_gui.vertexFormat.getStride());
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
				GuiTexture texture;
				texture.imgui = pcmd->TextureId;
				uint32_t flags = UINT32_MAX;
				if (texture.bgfx.flags & GuiTextureFlags::PointSampler)
					flags = BGFX_SAMPLER_POINT | BGFX_SAMPLER_UVW_CLAMP;
				bgfx::setTexture(0, s_gui.s_texture, texture.bgfx.handle, flags);
				bgfx::setIndexBuffer(&tib, firstIndex, pcmd->ElemCount);
				bgfx::setVertexBuffer(0, &tvb);
				bgfx::submit(kGuiView, s_gui.program);
			}
			firstIndex += pcmd->ElemCount;
		}
	}
}

bool guiColumnCheckbox(const char *label, const char *id, bool *value)
{
	ImGui::AlignTextToFramePadding();
	ImGui::Text("%s", label);
	ImGui::NextColumn();
	ImGui::PushItemWidth(-1);
	bool result = ImGui::Checkbox(id, value);
	ImGui::PopItemWidth();
	ImGui::NextColumn();
	return result;
}

bool guiColumnColorEdit(const char *label, const char *id, float *color)
{
	ImGui::AlignTextToFramePadding();
	ImGui::Text("%s", label);
	ImGui::NextColumn();
	ImGui::PushItemWidth(-1);
	bool result = ImGui::ColorEdit3(id, color, ImGuiColorEditFlags_NoInputs);
	ImGui::PopItemWidth();
	ImGui::NextColumn();
	return result;
}

bool guiColumnInputFloat(const char *label, const char *id, float *value, float step, float stepFast, const char *format)
{
	ImGui::AlignTextToFramePadding();
	ImGui::Text("%s", label);
	ImGui::NextColumn();
	ImGui::PushItemWidth(-1);
	bool result = ImGui::InputFloat(id, value, step, stepFast, format);
	ImGui::PopItemWidth();
	ImGui::NextColumn();
	return result;
}

bool guiColumnInputInt(const char *label, const char *id, int *value, int step)
{
	ImGui::AlignTextToFramePadding();
	ImGui::Text("%s", label);
	ImGui::NextColumn();
	ImGui::PushItemWidth(-1);
	bool result = ImGui::InputInt(id, value, step);
	ImGui::PopItemWidth();
	ImGui::NextColumn();
	return result;
}

bool guiColumnSliderInt(const char *label, const char *id, int *value, int valueMin, int valueMax)
{
	ImGui::AlignTextToFramePadding();
	ImGui::Text("%s", label);
	ImGui::NextColumn();
	ImGui::PushItemWidth(-1);
	bool result = ImGui::SliderInt(id, value, valueMin, valueMax);
	ImGui::PopItemWidth();
	ImGui::NextColumn();
	return result;
}
