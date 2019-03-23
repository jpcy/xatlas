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
#include <GLFW/glfw3.h>
#include <imgui/imgui.h>
#include <nativefiledialog/nfd.h>
#include "viewer.h"

bgfx::VertexDecl ModelVertex::decl;

struct ModelStatus
{
	enum Enum
	{
		NotLoaded,
		Loading,
		Finalizing,
		Loaded
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
	bgfx::ShaderHandle vs_model;
	bgfx::ShaderHandle vs_position;
	bgfx::ShaderHandle fs_color;
	bgfx::ShaderHandle fs_material;
	bgfx::ProgramHandle colorProgram;
	bgfx::ProgramHandle materialProgram;
	bgfx::UniformHandle u_diffuse;
	bgfx::UniformHandle u_emission;
	bgfx::UniformHandle u_lightDir_shadeType;
	bgfx::UniformHandle u_lightmapSampler;
	bgfx::UniformHandle u_color;
}
s_model;

void modelInit()
{
	s_model.u_color = bgfx::createUniform("u_color", bgfx::UniformType::Vec4);
	s_model.u_diffuse = bgfx::createUniform("u_diffuse", bgfx::UniformType::Vec4);
	s_model.u_emission = bgfx::createUniform("u_emission", bgfx::UniformType::Vec4);
	s_model.u_lightDir_shadeType = bgfx::createUniform("u_lightDir_shadeType", bgfx::UniformType::Vec4);
	s_model.u_lightmapSampler = bgfx::createUniform("u_lightmapSampler", bgfx::UniformType::Sampler);
	s_model.vs_model = loadShader(ShaderId::vs_model);
	s_model.vs_position = loadShader(ShaderId::vs_position);
	s_model.fs_color = loadShader(ShaderId::fs_color);
	s_model.fs_material = loadShader(ShaderId::fs_material);
	s_model.colorProgram = bgfx::createProgram(s_model.fs_color, s_model.vs_position);
	s_model.materialProgram = bgfx::createProgram(s_model.fs_material, s_model.vs_model);
	ModelVertex::init();
	bgfx::setViewClear(kModelView, BGFX_CLEAR_COLOR | BGFX_CLEAR_DEPTH, 0x444444ff);
	bgfx::setViewRect(kModelView, 0, 0, bgfx::BackbufferRatio::Equal);
}

void modelShutdown()
{
	modelDestroy();
	bgfx::destroy(s_model.u_color);
	bgfx::destroy(s_model.u_diffuse);
	bgfx::destroy(s_model.u_emission);
	bgfx::destroy(s_model.u_lightDir_shadeType);
	bgfx::destroy(s_model.u_lightmapSampler);
	bgfx::destroy(s_model.vs_model);
	bgfx::destroy(s_model.vs_position);
	bgfx::destroy(s_model.fs_color);
	bgfx::destroy(s_model.fs_material);
	bgfx::destroy(s_model.colorProgram);
	bgfx::destroy(s_model.materialProgram);
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
		setErrorMessage(objz_getError());
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

void modelFinalize()
{
	if (s_model.status.get() != ModelStatus::Finalizing)
		return;
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
	resetCamera();
	g_options.shadeMode = ShadeMode::Flat;
	g_options.wireframeMode = WireframeMode::Triangles;
	s_model.status.set(ModelStatus::Loaded);
}

void modelOpenDialog()
{
	if (s_model.status.get() == ModelStatus::Loading || s_model.status.get() == ModelStatus::Finalizing)
		return;
	if (!(atlasIsNotGenerated() || atlasIsReady()))
		return;
	nfdchar_t *filename = nullptr;
	nfdresult_t result = NFD_OpenDialog("obj", nullptr, &filename);
	if (result != NFD_OKAY)
		return;
	modelDestroy();
	s_model.status.set(ModelStatus::Loading);
	char windowTitle[256];
	snprintf(windowTitle, sizeof(windowTitle), "%s - %s\n", WINDOW_TITLE, filename);
	glfwSetWindowTitle(g_window, windowTitle);
	printf("Loading '%s'\n", filename);
	ModelLoadThreadArgs args;
	STRNCPY(args.filename, sizeof(args.filename), filename);
	s_model.thread = new std::thread(modelLoadThread, args);
	free(filename);
}

void modelDestroy()
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
	glfwSetWindowTitle(g_window, WINDOW_TITLE);
	s_model.status.set(ModelStatus::NotLoaded);
}

void modelSetMaterialUniforms(const objzMaterial *mat)
{
	if (!mat) {
		const float diffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f };
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

void modelRender(const float *view, const float *projection)
{
	if (s_model.status.get() != ModelStatus::Loaded)
		return;
	float modelMatrix[16];
	bx::mtxScale(modelMatrix, s_model.scale);
	bgfx::setViewTransform(kModelView, view, projection);
	const bool renderCharts = g_options.shadeMode == ShadeMode::Charts && atlasIsReady();
	if (g_options.shadeMode == ShadeMode::Flat || g_options.shadeMode == ShadeMode::Lightmap || renderCharts) {
		const float lightDir[] = { view[2], view[6], view[10], g_options.shadeMode == ShadeMode::Lightmap ? 1.0f : 0.0f };
		for (uint32_t i = 0; i < s_model.data->numMeshes; i++) {
			const objzMesh &mesh = s_model.data->meshes[i];
			const objzMaterial *mat = mesh.materialIndex == -1 ? nullptr : &s_model.data->materials[mesh.materialIndex];
			// When rendering charts, emissive meshes won't be rendered, so do that here.
			if (renderCharts && (!mat || (mat->emission[0] <= 0.0f && mat->emission[1] <= 0.0f && mat->emission[2] <= 0.0f)))
				continue;
			if (atlasIsReady()) {
				bgfx::setIndexBuffer(atlasGetIb(), mesh.firstIndex, mesh.numIndices);
				bgfx::setVertexBuffer(0, atlasGetVb());
			} else {
				bgfx::setIndexBuffer(s_model.ib, mesh.firstIndex, mesh.numIndices);
				bgfx::setVertexBuffer(0, s_model.vb);
			}
			bgfx::setState(BGFX_STATE_DEFAULT);
			bgfx::setTransform(modelMatrix);
			bgfx::setUniform(s_model.u_lightDir_shadeType, lightDir);
			modelSetMaterialUniforms(mat);
			if (g_options.shadeMode == ShadeMode::Lightmap)
				bgfx::setTexture(0, s_model.u_lightmapSampler, bakeGetLightmap());
			bgfx::submit(kModelView, s_model.materialProgram);
		}
	}
	if (renderCharts)
		atlasRenderCharts(modelMatrix);
	if (g_options.wireframe) {
		const float color[] = { 0.5f, 0.5f, 0.5f, 0.5f };
		bgfx::setUniform(s_model.u_color, color);
		bgfx::setState(BGFX_STATE_WRITE_RGB | BGFX_STATE_WRITE_A | BGFX_STATE_PT_LINES | BGFX_STATE_BLEND_ALPHA);
		bgfx::setTransform(modelMatrix);
		if (g_options.wireframeMode == WireframeMode::Triangles) {
			bgfx::setIndexBuffer(s_model.wireframeIb);
			bgfx::setVertexBuffer(0, s_model.vb);
		} else if (atlasIsReady() && g_options.wireframeMode == WireframeMode::Charts) {
			bgfx::setVertexBuffer(0, atlasGetChartBoundaryVb());
		}
		bgfx::submit(kModelView, s_model.colorProgram);
	}
}

void modelShowGuiOptions()
{
	ImGui::Text("%u objects", s_model.data->numObjects);
	ImGui::Text("%u vertices", s_model.data->numVertices);
	ImGui::Text("%u triangles", s_model.data->numIndices / 3);
	ImGui::InputFloat("Model scale", &s_model.scale, 0.01f, 0.1f);
	s_model.scale = bx::max(0.001f, s_model.scale);
}

void modelShowGuiWindow(int progressDots)
{
	ImGuiIO &io = ImGui::GetIO();
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
}

AABB modelGetAABB()
{
	return s_model.aabb;
}

const objzModel *modelGetData()
{
	return s_model.data;
}

bx::Vec3 modelGetCentroid()
{
	return bx::mul(s_model.centroid, s_model.scale);
}

bgfx::ShaderHandle modelGet_vs_position()
{
	return s_model.vs_position;
}

bgfx::ShaderHandle modelGet_vs_model()
{
	return s_model.vs_model;
}

bool modelIsLoaded()
{
	return s_model.status.get() == ModelStatus::Loaded;
}
