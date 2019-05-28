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
#include <imgui/imgui.h>

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
#include <igl/arap.h>
#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/lscm.h>
#include <igl/map_vertices_to_circle.h>
#ifdef _MSC_VER
#pragma warning(pop)
#endif
#endif

#include "../xatlas.h"
#include "viewer.h"

static const uint8_t kPaletteBlack = 0;

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
	int chartTextureCellSize = 0;
	xatlas::Atlas *data = nullptr;
	std::thread *thread = nullptr;
	AtlasStatus status;
	int currentTexture;
	bool fitToWindow = true;
	bool wireframe = false;
	std::vector<bgfx::FrameBufferHandle> chartsFrameBuffers;
	bgfx::VertexBufferHandle vb = BGFX_INVALID_HANDLE;
	bgfx::IndexBufferHandle ib = BGFX_INVALID_HANDLE;
	bgfx::IndexBufferHandle chartIb = BGFX_INVALID_HANDLE;
	bgfx::VertexBufferHandle chartBoundaryVb = BGFX_INVALID_HANDLE;
	xatlas::ChartOptions chartOptions;
	bool chartOptionsChanged = false;
	xatlas::PackOptions packOptions;
	bool packOptionsChanged = false;
	ParamMethod paramMethod = ParamMethod::LSCM;
	bool paramMethodChanged = false;
	std::vector<ModelVertex> vertices;
	std::vector<uint32_t> indices;
	std::vector<uint32_t> chartIndices;
	std::vector<bool> boundaryEdges;
	std::vector<bx::Vec3> chartBoundaryVertices;
	bgfx::VertexDecl wireVertexDecl;
	// Chart rendering with checkerboard pattern.
	bgfx::UniformHandle u_color;
	bgfx::UniformHandle u_textureSize_cellSize;
	bgfx::ShaderHandle fs_chart;
	bgfx::ShaderHandle vs_chart;
	bgfx::ShaderHandle vs_chartTexcoordSpace;
	bgfx::ProgramHandle chartProgram;
	bgfx::ProgramHandle chartTexcoordSpaceProgram;
}
s_atlas;

static void clearPackOptions()
{
	s_atlas.packOptions = xatlas::PackOptions();
	// Baking needs 1 pixel padding for dilate filter.
	s_atlas.packOptions.conservative = true;
	s_atlas.packOptions.padding = 1;
}

void atlasInit()
{
	s_atlas.wireVertexDecl
		.begin()
		.add(bgfx::Attrib::Position, 3, bgfx::AttribType::Float)
		.end();
	assert(s_atlas.wireVertexDecl.getStride() == sizeof(bx::Vec3));
	bgfx::setPaletteColor(kPaletteBlack, 0x000000ff);
	s_atlas.u_color = bgfx::createUniform("u_color", bgfx::UniformType::Vec4);
	s_atlas.u_textureSize_cellSize = bgfx::createUniform("u_textureSize_cellSize", bgfx::UniformType::Vec4);
	s_atlas.fs_chart = loadShader(ShaderId::fs_chart);
	s_atlas.vs_chart = loadShader(ShaderId::vs_chart);
	s_atlas.vs_chartTexcoordSpace = loadShader(ShaderId::vs_chartTexcoordSpace);
	s_atlas.chartProgram = bgfx::createProgram(s_atlas.vs_chart, s_atlas.fs_chart);
	s_atlas.chartTexcoordSpaceProgram = bgfx::createProgram(s_atlas.vs_chartTexcoordSpace, s_atlas.fs_chart);
	clearPackOptions();
}

void atlasShutdown()
{
	bgfx::destroy(s_atlas.u_color);
	bgfx::destroy(s_atlas.u_textureSize_cellSize);
	bgfx::destroy(s_atlas.fs_chart);
	bgfx::destroy(s_atlas.vs_chart);
	bgfx::destroy(s_atlas.vs_chartTexcoordSpace);
	bgfx::destroy(s_atlas.chartProgram);
	bgfx::destroy(s_atlas.chartTexcoordSpaceProgram);
}

void atlasDestroy()
{
	bakeClear();
	if (s_atlas.thread) {
		if (s_atlas.thread->joinable())
			s_atlas.thread->join();
		delete s_atlas.thread;
		s_atlas.thread = nullptr;
	}
	if (s_atlas.data) {
		xatlas::Destroy(s_atlas.data);
		s_atlas.data = nullptr;
	}
	for (uint32_t i = 0; i < (uint32_t)s_atlas.chartsFrameBuffers.size(); i++) {
		if (bgfx::isValid(s_atlas.chartsFrameBuffers[i])) {
			bgfx::destroy(s_atlas.chartsFrameBuffers[i]);
			s_atlas.chartsFrameBuffers[i] = BGFX_INVALID_HANDLE;
		}
	}
	s_atlas.chartsFrameBuffers.clear();
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
	const objzModel *model = modelGetData();
	int progress = 0;
	const bool firstRun = !s_atlas.data;
	const clock_t startTime = clock();
	if (firstRun) {
		// Create xatlas context on first run only.
		s_atlas.data = xatlas::Create();
		std::vector<uint8_t> ignoreFaces; // Should be bool, workaround stupid C++ specialization.
		for (uint32_t i = 0; i < model->numObjects; i++) {
			const objzObject &object = model->objects[i];
			auto v = &((const ModelVertex *)model->vertices)[object.firstVertex];
			// Ignore faces with an emissive material.
			ignoreFaces.resize(object.numIndices / 3);
			memset(ignoreFaces.data(), 0, ignoreFaces.size() * sizeof(uint8_t));
			for (uint32_t j = 0; j < object.numMeshes; j++) {
				const objzMesh &mesh = model->meshes[object.firstMesh + j];
				const objzMaterial *mat = mesh.materialIndex == -1 ? nullptr : &model->materials[mesh.materialIndex];
				if (mat && (mat->emission[0] > 0.0f || mat->emission[1] > 0.0f || mat->emission[2] > 0.0f)) {
					for (uint32_t k = 0; k < mesh.numIndices / 3; k++)
						ignoreFaces[(mesh.firstIndex - object.firstIndex) / 3 + k] = true;
				}
			}
			xatlas::MeshDecl meshDecl;
			meshDecl.vertexCount = object.numVertices;
			meshDecl.vertexPositionData = &v->pos;
			meshDecl.vertexPositionStride = sizeof(ModelVertex);
			meshDecl.vertexNormalData = &v->normal;
			meshDecl.vertexNormalStride = sizeof(ModelVertex);
			meshDecl.vertexUvData = &v->texcoord;
			meshDecl.vertexUvStride = sizeof(ModelVertex);
			meshDecl.indexCount = object.numIndices;
			meshDecl.indexData = &((uint32_t *)model->indices)[object.firstIndex];
			meshDecl.indexFormat = xatlas::IndexFormat::UInt32;
			meshDecl.indexOffset = -(int32_t)object.firstVertex;
			meshDecl.faceIgnoreData = (const bool *)ignoreFaces.data();
			xatlas::AddMeshError::Enum error = xatlas::AddMesh(s_atlas.data, meshDecl);
			if (error != xatlas::AddMeshError::Success) {
				fprintf(stderr, "Error adding mesh: %s\n", xatlas::StringForEnum(error));
				setErrorMessage("Error adding mesh: %s", xatlas::StringForEnum(error));
				xatlas::Destroy(s_atlas.data);
				s_atlas.data = nullptr;
				s_atlas.status.set(AtlasStatus::NotGenerated);
				return;
			}
			const int newProgress = int((i + 1) / (float)model->numObjects * 100.0f);
			if (newProgress != progress) {
				progress = newProgress;
				s_atlas.status.setProgress((xatlas::ProgressCategory::Enum)-1, progress);
			}
		}
	}
	s_atlas.status.set(AtlasStatus::Generating);
	if (firstRun || s_atlas.chartOptionsChanged) {
		xatlas::ComputeCharts(s_atlas.data, s_atlas.chartOptions, atlasProgressCallback);
	}
	if (firstRun || s_atlas.chartOptionsChanged || s_atlas.paramMethodChanged) {
		xatlas::ParameterizeFunc paramFunc = nullptr;
#if USE_LIBIGL
		if (s_atlas.paramMethod != ParamMethod::LSCM)
			paramFunc = atlasParameterizationCallback;
#endif
		xatlas::ParameterizeCharts(s_atlas.data, paramFunc, atlasProgressCallback);
	}
	if (firstRun || s_atlas.chartOptionsChanged || s_atlas.paramMethodChanged || s_atlas.packOptionsChanged) {
		xatlas::PackCharts(s_atlas.data, s_atlas.packOptions, atlasProgressCallback);
	}
	const double elapsedTime = (clock() - startTime) * 1000.0 / CLOCKS_PER_SEC;
	printf("Generated atlas in %.2f seconds (%g ms)\n", elapsedTime / 1000.0, elapsedTime);
	s_atlas.chartOptionsChanged = false;
	s_atlas.paramMethodChanged = false;
	s_atlas.packOptionsChanged = false;
	// Find chart boundary edges.
	uint32_t numEdges = 0;
	for (uint32_t i = 0; i < s_atlas.data->meshCount; i++) {
		const xatlas::Mesh &outputMesh = s_atlas.data->meshes[i];
		numEdges += outputMesh.indexCount;
	}
	s_atlas.boundaryEdges.resize(numEdges);
	std::unordered_map<EdgeKey, uint32_t, EdgeKeyHash, EdgeKeyEqual> edgeMap;
	numEdges = 0;
	for (uint32_t i = 0; i < s_atlas.data->meshCount; i++) {
		const xatlas::Mesh &mesh = s_atlas.data->meshes[i];
		const objzObject &object = model->objects[i];
		const ModelVertex *vertices = &((const ModelVertex *)model->vertices)[object.firstVertex];
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
					s_atlas.boundaryEdges[numEdges] = edgeMap.count(key) == 0;
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
		const objzObject &object = model->objects[i];
		const ModelVertex *oldVertices = &((const ModelVertex *)model->vertices)[object.firstVertex];
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
					if (s_atlas.boundaryEdges[numEdges]) {
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
	s_atlas.status.set(AtlasStatus::Finalizing);
}

void atlasGenerate()
{
	if (!(s_atlas.status.get() == AtlasStatus::NotGenerated || s_atlas.status.get() == AtlasStatus::Ready))
		return;
	if (s_atlas.data && !s_atlas.chartOptionsChanged && !s_atlas.paramMethodChanged && !s_atlas.packOptionsChanged) {
		// Already have an atlas and none of the options that affect atlas creation have changed.
		return;
	}
	bakeClear();
	xatlas::SetPrint(printf, true);
	g_options.shadeMode = ShadeMode::Flat;
	g_options.wireframeMode = WireframeMode::Triangles;
	s_atlas.status.set(AtlasStatus::AddingMeshes);
	s_atlas.thread = new std::thread(atlasGenerateThread);
}

static void atlasRenderChartsTextures()
{
	bgfx::ViewId viewId = kFirstFreeView;
	float bottom = (float)s_atlas.data->height;
	float top = 0.0f;
	if (bgfx::getCaps()->originBottomLeft)
		bx::swap(bottom, top);
	float projection[16];
	bx::mtxOrtho(projection, 0.0f, (float)s_atlas.data->width, bottom, top, 0.0f, 1.0f, 0.0f, bgfx::getCaps()->homogeneousDepth);
	for (uint32_t i = 0; i < (uint32_t)s_atlas.chartsFrameBuffers.size(); i++) {
		// Setup view for rendering into atlas texture.
		bgfx::setViewClear(viewId, BGFX_CLEAR_COLOR);
		bgfx::setViewRect(viewId, 0, 0, (uint16_t)s_atlas.data->width, (uint16_t)s_atlas.data->height);
		bgfx::setViewFrameBuffer(viewId, s_atlas.chartsFrameBuffers[i]);
		bgfx::setViewTransform(viewId, nullptr, projection);
		// Render charts with checkerboard outline.
		srand(s_atlas.chartColorSeed);
		uint32_t firstIndex = 0;
		for (uint32_t mi = 0; mi < s_atlas.data->meshCount; mi++) {
			const xatlas::Mesh &mesh = s_atlas.data->meshes[mi];
			for (uint32_t ci = 0; ci < mesh.chartCount; ci++) {
				const xatlas::Chart &chart = mesh.chartArray[ci];
				if (chart.atlasIndex == i) {
					uint8_t bcolor[3];
					randomRGB(bcolor);
					float color[4];
					color[0] = bcolor[0] / 255.0f;
					color[1] = bcolor[1] / 255.0f;
					color[2] = bcolor[2] / 255.0f;
					color[3] = 1.0f;
					bgfx::setUniform(s_atlas.u_color, color);
					float textureSize_cellSize[4];
					textureSize_cellSize[0] = (float)s_atlas.data->width;
					textureSize_cellSize[1] = (float)s_atlas.data->height;
					textureSize_cellSize[2] = (float)s_atlas.chartTextureCellSize;
					textureSize_cellSize[3] = (float)s_atlas.chartTextureCellSize;
					bgfx::setUniform(s_atlas.u_textureSize_cellSize, textureSize_cellSize);
					bgfx::setState(BGFX_STATE_WRITE_RGB | BGFX_STATE_WRITE_A);
					bgfx::setIndexBuffer(s_atlas.chartIb, firstIndex, chart.indexCount);
					bgfx::setVertexBuffer(0, s_atlas.vb);
					bgfx::submit(viewId, s_atlas.chartTexcoordSpaceProgram);
				}
				firstIndex += chart.indexCount;
			}
		}
		if (s_atlas.wireframe) {
			// Render chart boundary lines.
			std::vector<PosVertex> boundaryVertices;
			int edge = 0;
			for (uint32_t mi = 0; mi < s_atlas.data->meshCount; mi++) {
				const xatlas::Mesh &mesh = s_atlas.data->meshes[mi];
				for (uint32_t ci = 0; ci < mesh.chartCount; ci++) {
					const xatlas::Chart &chart = mesh.chartArray[ci];
					for (uint32_t k = 0; k < chart.indexCount; k += 3) {
						for (int l = 0; l < 3; l++) {
							if (chart.atlasIndex == i && s_atlas.boundaryEdges[edge]) {
								const xatlas::Vertex &v0 = mesh.vertexArray[chart.indexArray[k + l]];
								const xatlas::Vertex &v1 = mesh.vertexArray[chart.indexArray[k + (l + 1) % 3]];
								PosVertex p;
								p.pos[0] = v0.uv[0];
								p.pos[1] = v0.uv[1];
								boundaryVertices.push_back(p);
								p.pos[0] = v1.uv[0];
								p.pos[1] = v1.uv[1];
								boundaryVertices.push_back(p);
							}
							edge++;
						}
					}
				}
			}
			if (!boundaryVertices.empty() && bgfx::getAvailTransientVertexBuffer((uint32_t)boundaryVertices.size(), PosVertex::decl) == (uint32_t)boundaryVertices.size()) {
				bgfx::TransientVertexBuffer vb;
				bgfx::allocTransientVertexBuffer(&vb, (uint32_t)boundaryVertices.size(), PosVertex::decl);
				memcpy(vb.data, boundaryVertices.data(), boundaryVertices.size() * sizeof(PosVertex));
				bgfx::setState(BGFX_STATE_WRITE_RGB | BGFX_STATE_WRITE_A | BGFX_STATE_PT_LINES | BGFX_STATE_BLEND_ALPHA);
				bgfx::setVertexBuffer(0, &vb);
				const float color[] = { 1.0f, 1.0f, 1.0f, 1.0f };
				bgfx::setUniform(s_atlas.u_color, color);
				bgfx::submit(viewId, getColorProgram(), 1);
			}
		}
		viewId++;
	}
}

void atlasFinalize()
{
	if (s_atlas.status.get() != AtlasStatus::Finalizing)
		return;
	if (s_atlas.thread) {
		if (s_atlas.thread->joinable())
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
	// Create framebuffer/texture for atlas.
	for (uint32_t i = 0; i < (uint32_t)s_atlas.chartsFrameBuffers.size(); i++) {
		bgfx::destroy(s_atlas.chartsFrameBuffers[i]);
		s_atlas.chartsFrameBuffers[i] = BGFX_INVALID_HANDLE;
	}
	s_atlas.chartsFrameBuffers.resize(s_atlas.data->atlasCount);
	for (uint32_t i = 0; i < (uint32_t)s_atlas.chartsFrameBuffers.size(); i++) {
		bgfx::TextureHandle texture = bgfx::createTexture2D((uint16_t)s_atlas.data->width, (uint16_t)s_atlas.data->height, false, 1, bgfx::TextureFormat::BGRA8, BGFX_TEXTURE_RT | BGFX_SAMPLER_UVW_BORDER | BGFX_SAMPLER_BORDER_COLOR(kPaletteBlack));
		s_atlas.chartsFrameBuffers[i] = bgfx::createFrameBuffer(1, &texture, true);
	}
	// Render charts to texture(s) for previewing UVs. Render in UV space.
	atlasRenderChartsTextures();
	s_atlas.currentTexture = 0;
	g_options.shadeMode = ShadeMode::Charts;
	g_options.wireframeMode = WireframeMode::Charts;
	s_atlas.status.set(AtlasStatus::Ready);
}

void atlasRenderCharts(const float *modelMatrix)
{
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
			bgfx::setUniform(s_atlas.u_color, color);
			float textureSize_cellSize[4];
			textureSize_cellSize[0] = (float)s_atlas.data->width;
			textureSize_cellSize[1] = (float)s_atlas.data->height;
			textureSize_cellSize[2] = (float)g_options.chartCellSize;
			textureSize_cellSize[3] = (float)g_options.chartCellSize;
			bgfx::setUniform(s_atlas.u_textureSize_cellSize, textureSize_cellSize);
			bgfx::setState(BGFX_STATE_DEFAULT);
			bgfx::setTransform(modelMatrix);
			bgfx::setIndexBuffer(s_atlas.chartIb, firstIndex, chart.indexCount);
			bgfx::setVertexBuffer(0, s_atlas.vb);
			bgfx::submit(kModelView, s_atlas.chartProgram);
			firstIndex += chart.indexCount;
		}
	}
}

void atlasRenderChartsWireframe(const float *modelMatrix)
{
	const float color[] = { 1.0f, 1.0f, 1.0f, 0.5f };
	bgfx::setUniform(s_atlas.u_color, color);
	bgfx::setState(BGFX_STATE_WRITE_RGB | BGFX_STATE_WRITE_A | BGFX_STATE_PT_LINES | BGFX_STATE_BLEND_ALPHA);
	bgfx::setTransform(modelMatrix);
	bgfx::setVertexBuffer(0, s_atlas.chartBoundaryVb);
	bgfx::submit(kModelView, getColorProgram());
}

void atlasShowGuiOptions()
{
	if (!(s_atlas.status.get() == AtlasStatus::NotGenerated || s_atlas.status.get() == AtlasStatus::Ready))
		return;
	const ImVec2 buttonSize(ImVec2(ImGui::GetContentRegionAvailWidth() * 0.35f, 0.0f));
	ImGui::Text("Atlas");
	ImGui::Spacing();
	if (ImGui::Button("Generate", buttonSize))
		atlasGenerate();
	if (s_atlas.status.get() == AtlasStatus::Ready) {
		uint32_t numIndices = 0, numVertices = 0;
		for (uint32_t i = 0; i < s_atlas.data->meshCount; i++) {
			const xatlas::Mesh &outputMesh = s_atlas.data->meshes[i];
			numIndices += outputMesh.indexCount;
			numVertices += outputMesh.vertexCount;
		}
		ImGui::Text("%u atlas%s", s_atlas.data->atlasCount, s_atlas.data->atlasCount > 1 ? "es" : "");
		ImGui::Text("%ux%u resolution", s_atlas.data->width, s_atlas.data->height);
		ImGui::Text("%u charts", s_atlas.data->chartCount);
		ImGui::Text("%u vertices", numVertices);
		ImGui::Text("%u triangles", numIndices / 3);
		ImGui::Text("%g texels per unit", s_atlas.data->texelsPerUnit);
	}
	if (ImGui::TreeNodeEx("Chart options", ImGuiTreeNodeFlags_FramePadding)) {
		bool changed = false;
		changed |= ImGui::InputFloat("Proxy fit metric weight", &s_atlas.chartOptions.proxyFitMetricWeight);
		changed |= ImGui::InputFloat("Roundness metric weight", &s_atlas.chartOptions.roundnessMetricWeight);
		changed |= ImGui::InputFloat("Straightness metric weight", &s_atlas.chartOptions.straightnessMetricWeight);
		changed |= ImGui::InputFloat("Normal seam metric weight", &s_atlas.chartOptions.normalSeamMetricWeight);
		changed |= ImGui::InputFloat("Texture seam metric weight", &s_atlas.chartOptions.textureSeamMetricWeight);
		changed |= ImGui::InputFloat("Max chart area", &s_atlas.chartOptions.maxChartArea);
		changed |= ImGui::InputFloat("Max boundary length", &s_atlas.chartOptions.maxBoundaryLength);
		changed |= ImGui::InputFloat("Max threshold", &s_atlas.chartOptions.maxThreshold);
		changed |= ImGui::InputInt("Grow face count", (int *)&s_atlas.chartOptions.growFaceCount);
		changed |= ImGui::InputInt("Max iterations", (int *)&s_atlas.chartOptions.maxIterations);
		if (ImGui::Button("Reset to default", buttonSize)) {
			s_atlas.chartOptions = xatlas::ChartOptions();
			changed = true;
		}
		ImGui::TreePop();
		if (changed)
			s_atlas.chartOptionsChanged = true;
	}
#if USE_LIBIGL
	if (ImGui::TreeNodeEx("Parameterization options", ImGuiTreeNodeFlags_FramePadding)) {
		const ParamMethod oldParamMethod = s_atlas.paramMethod;
		ImGui::RadioButton("LSCM", (int *)&s_atlas.paramMethod, (int)ParamMethod::LSCM);
		ImGui::RadioButton("libigl Harmonic", (int *)&s_atlas.paramMethod, (int)ParamMethod::libigl_Harmonic);
		ImGui::RadioButton("libigl LSCM", (int *)&s_atlas.paramMethod, (int)ParamMethod::libigl_LSCM);
		ImGui::RadioButton("libigl ARAP", (int *)&s_atlas.paramMethod, (int)ParamMethod::libigl_ARAP);
		ImGui::TreePop();
		if (s_atlas.paramMethod != oldParamMethod)
			s_atlas.paramMethodChanged = true;
	}
#endif
	if (ImGui::TreeNodeEx("Pack options", ImGuiTreeNodeFlags_FramePadding)) {
		bool changed = false;
		changed |= ImGui::SliderInt("Attempts", &s_atlas.packOptions.attempts, 0, 4096);
		changed |= ImGui::InputFloat("Texels per unit", &s_atlas.packOptions.texelsPerUnit, 0.0f, 32.0f, 2);
		changed |= ImGui::InputInt("Resolution", (int *)&s_atlas.packOptions.resolution, 8);
		changed |= ImGui::InputInt("Max chart size", (int *)&s_atlas.packOptions.maxChartSize);
		changed |= ImGui::Checkbox("Block align", &s_atlas.packOptions.blockAlign);
		ImGui::SameLine();
		changed |= ImGui::Checkbox("Conservative", &s_atlas.packOptions.conservative);
		changed |= ImGui::SliderInt("Padding", (int *)&s_atlas.packOptions.padding, 0, 8);
		if (ImGui::Button("Reset to default", buttonSize)) {
			clearPackOptions();
			changed = true;
		}
		ImGui::TreePop();
		if (changed)
			s_atlas.packOptionsChanged = true;
	}
}

void atlasShowGuiWindow(int progressDots)
{
	ImGuiIO &io = ImGui::GetIO();
	const float margin = 4.0f;
	const ImGuiWindowFlags progressWindowFlags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings;
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
	} else if (atlasStatus == AtlasStatus::Ready && g_options.showAtlasWindow) {
		const float size = 500;
		ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x - size - margin, margin), ImGuiCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(size, size), ImGuiCond_FirstUseEver);
		if (ImGui::Begin("Atlas", &g_options.showAtlasWindow, ImGuiWindowFlags_HorizontalScrollbar)) {
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
				ImGui::SameLine();
			}
			ImGui::Checkbox("Fit to window", &s_atlas.fitToWindow);
			ImGui::SameLine();
			if (ImGui::Checkbox("Wireframe", &s_atlas.wireframe))
				atlasRenderChartsTextures();
			GuiTexture texture;
			texture.bgfx.handle = bgfx::getTexture(s_atlas.chartsFrameBuffers[s_atlas.currentTexture]);
			texture.bgfx.flags = GuiTextureFlags::PointSampler;
			if (s_atlas.fitToWindow)
				ImGui::Image(texture.imgui, ImGui::GetContentRegionAvail());
			else 
				ImGui::Image(texture.imgui, ImVec2((float)s_atlas.data->width, (float)s_atlas.data->height));
			ImGui::End();
		}
	}
}

uint32_t atlasGetCount()
{
	return s_atlas.data->atlasCount;
}

uint32_t atlasGetWidth()
{
	return s_atlas.data->width;
}

uint32_t atlasGetHeight()
{
	return s_atlas.data->height;
}

std::vector<ModelVertex> *atlasGetVertices()
{
	return &s_atlas.vertices;
}

std::vector<uint32_t> *atlasGetIndices()
{
	return &s_atlas.indices;
}

bgfx::VertexBufferHandle atlasGetVb()
{
	return s_atlas.vb;
}

bgfx::IndexBufferHandle atlasGetIb()
{
	return s_atlas.ib;
}

bool atlasIsNotGenerated()
{
	return s_atlas.status.get() == AtlasStatus::NotGenerated;
}

bool atlasIsReady()
{
	return s_atlas.status.get() == AtlasStatus::Ready;
}
