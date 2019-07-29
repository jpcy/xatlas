/*
xatlas
https://github.com/jpcy/xatlas
Copyright (c) 2018 Jonathan Young

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
/*
*  Copyright (c) 2004-2010, Bruno Levy
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions are met:
*
*  * Redistributions of source code must retain the above copyright notice,
*  this list of conditions and the following disclaimer.
*  * Redistributions in binary form must reproduce the above copyright notice,
*  this list of conditions and the following disclaimer in the documentation
*  and/or other materials provided with the distribution.
*  * Neither the name of the ALICE Project-Team nor the names of its
*  contributors may be used to endorse or promote products derived from this
*  software without specific prior written permission.
* 
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
*  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
*  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
*  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
*  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
*  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
*  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
*  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
*  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
*  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
*
*  If you modify this software, you should include a notice giving the
*  name of the person performing the modification, the date of modification,
*  and the reason for such modification.
*
*  Contact: Bruno Levy
*
*     levy@loria.fr
*
*     ALICE Project
*     LORIA, INRIA Lorraine, 
*     Campus Scientifique, BP 239
*     54506 VANDOEUVRE LES NANCY CEDEX 
*     FRANCE
*
*/
#include <mutex>
#include <thread>
#include <unordered_map>
#include <imgui/imgui.h>

#define USE_LIBIGL 0
#define USE_OPENNL 1

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

#if USE_OPENNL
#define GEO_STATIC_LIBS
#include <OpenNL_psm.h>
#endif

#include "../xatlas.h"
#include "viewer.h"

namespace std { typedef std::lock_guard<std::mutex> mutex_lock; }

namespace bgfx
{
	template<typename T>
	const Memory* makeRef(const std::vector<T> &_data)
	{
		return bgfx::makeRef(_data.data(), uint32_t(_data.size() * sizeof(T)));
	}

	template<typename T>
	static void destroyAndClear(T &handle)
	{
		if (bgfx::isValid(handle)) {
			bgfx::destroy(handle);
			handle = BGFX_INVALID_HANDLE;
		}
	}
}

static const uint8_t kPaletteBlack = 0;

struct AtlasStatus
{
	enum Enum
	{
		NotGenerated,
		Generating,
		Finalizing,
		Ready
	};

	Enum get()
	{
		std::mutex_lock lock(m_lock);
		return m_value;
	}

	void set(Enum value)
	{
		std::mutex_lock lock(m_lock);
		m_value = value;
	}

	void getProgress(xatlas::ProgressCategory::Enum *category, int *progress)
	{
		std::mutex_lock lock(m_lock);
		*category = m_category;
		*progress = m_progress;
	}

	void setProgress(xatlas::ProgressCategory::Enum category, int progress)
	{
		std::mutex_lock lock(m_lock);
		m_category = category;
		m_progress = progress;
	}

	bool getCancel()
	{
		std::mutex_lock lock(m_lock);
		return m_cancel;
	}

	void setCancel(bool value)
	{
		std::mutex_lock lock(m_lock);
		m_cancel = value;
	}

private:
	std::mutex m_lock;
	Enum m_value = NotGenerated;
	bool m_cancel = false;
	xatlas::ProgressCategory::Enum m_category;
	int m_progress = 0;
};

enum class ParamMethod
{
	LSCM,
	libigl_Harmonic,
	libigl_LSCM,
	libigl_ARAP,
	OpenNL_LSCM
};

struct AtlasVertex
{
	float pos[2];
	uint32_t color;
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
	bool showBlockGrid = false;
	bgfx::FrameBufferHandle chartsFrameBuffer = BGFX_INVALID_HANDLE;
	bgfx::VertexBufferHandle vb = BGFX_INVALID_HANDLE;
	bgfx::IndexBufferHandle ib = BGFX_INVALID_HANDLE;
	bgfx::VertexBufferHandle chartColorVb = BGFX_INVALID_HANDLE;
	bgfx::VertexBufferHandle chartInvalidColorVb = BGFX_INVALID_HANDLE;
	bgfx::IndexBufferHandle chartIb = BGFX_INVALID_HANDLE;
	bgfx::VertexBufferHandle chartBoundaryVb = BGFX_INVALID_HANDLE;
	bgfx::IndexBufferHandle textureChartIb = BGFX_INVALID_HANDLE;
	bgfx::VertexBufferHandle textureChartBoundaryVb = BGFX_INVALID_HANDLE;
	xatlas::ChartOptions chartOptions;
	bool chartOptionsChanged = false;
	xatlas::PackOptions packOptions;
	bool packOptionsChanged = false;
	ParamMethod paramMethod = ParamMethod::LSCM;
	std::mutex paramMutex; // Used by OpenNL
	bool paramOptionsChanged = false;
	std::vector<ModelVertex> vertices;
	std::vector<uint32_t> indices;
	std::vector<uint32_t> chartColorVertices;
	std::vector<uint32_t> chartInvalidColorVertices;
	std::vector<uint32_t> chartIndices;
	std::vector<bool> boundaryEdges;
	std::vector<WireframeVertex> chartBoundaryVertices;
	std::vector<uint32_t> textureChartIndices;
	std::vector<AtlasVertex> textureChartBoundaryVertices;
	bgfx::VertexDecl wireVertexDecl;
	bgfx::VertexDecl chartColorDecl;
	bgfx::VertexDecl atlasVertexDecl;
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
	s_atlas.packOptions.padding = 1;
	s_atlas.packOptions.createImage = true;
}

void atlasInit()
{
	s_atlas.wireVertexDecl
		.begin()
		.add(bgfx::Attrib::Position, 3, bgfx::AttribType::Float)
		.end();
	s_atlas.chartColorDecl.begin()
		.add(bgfx::Attrib::Color0, 4, bgfx::AttribType::Uint8, true)
		.end();
	s_atlas.atlasVertexDecl.begin()
		.add(bgfx::Attrib::Position, 2, bgfx::AttribType::Float)
		.add(bgfx::Attrib::Color0, 4, bgfx::AttribType::Uint8, true)
		.end();
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
	bgfx::destroyAndClear(s_atlas.chartsFrameBuffer);
	bgfx::destroyAndClear(s_atlas.chartsFrameBuffer);
	bgfx::destroyAndClear(s_atlas.vb);
	bgfx::destroyAndClear(s_atlas.ib);
	bgfx::destroyAndClear(s_atlas.chartColorVb);
	bgfx::destroyAndClear(s_atlas.chartInvalidColorVb);
	bgfx::destroyAndClear(s_atlas.chartIb);
	bgfx::destroyAndClear(s_atlas.chartBoundaryVb);
	bgfx::destroyAndClear(s_atlas.textureChartIb);
	bgfx::destroyAndClear(s_atlas.textureChartBoundaryVb);
	s_atlas.status.set(AtlasStatus::NotGenerated);
}

static bool atlasProgressCallback(xatlas::ProgressCategory::Enum category, int progress, void * /*userData*/)
{
	s_atlas.status.setProgress(category, progress);
	if (s_atlas.status.getCancel())
		return false;
	return true;
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

#if USE_OPENNL
static bx::Vec3 operator-(bx::Vec3 a, bx::Vec3 b) { return bx::sub(a, b); }

struct vec2
{
	vec2() {}
	vec2(double x, double y) : x(x), y(y) {}
	double x, y;
};

/**
* \brief Computes the coordinates of the vertices of a triangle
* in a local 2D orthonormal basis of the triangle's plane.
* \param[in] p0 , p1 , p2 the 3D coordinates of the vertices of 
*   the triangle
* \param[out] z0 , z1 , z2 the 2D coordinates of the vertices of
*   the triangle
*/
static void project_triangle(
	const bx::Vec3& p0, 
	const bx::Vec3& p1, 
	const bx::Vec3& p2,
	vec2& z0,
	vec2& z1,
	vec2& z2
) {
	bx::Vec3 X = p1 - p0;
	X = bx::normalize(X);
	bx::Vec3 Z = cross(X,(p2 - p0));
	Z = bx::normalize(Z);
	bx::Vec3 Y = cross(Z,X);
	const bx::Vec3& O = p0;

	double x0 = 0;
	double y0 = 0;
	double x1 = bx::length(p1 - O);
	double y1 = 0;
	double x2 = dot((p2 - O),X);
	double y2 = dot((p2 - O),Y);        

	z0 = vec2(x0,y0);
	z1 = vec2(x1,y1);
	z2 = vec2(x2,y2);        
}
#endif

#if USE_LIBIGL || USE_OPENNL
static void atlasParameterizationCallback(const float *positions, float *texcoords, uint32_t vertexCount, const uint32_t *indices, uint32_t indexCount)
{
#if USE_LIBIGL
	if (s_atlas.paramMethod == ParamMethod::libigl_Harmonic || s_atlas.paramMethod == ParamMethod::libigl_LSCM || s_atlas.paramMethod == ParamMethod::libigl_ARAP) {
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
#if USE_OPENNL
	if (s_atlas.paramMethod == ParamMethod::OpenNL_LSCM) {
		std::mutex_lock lock(s_atlas.paramMutex); // Current context should really be TLS internally...
		/**
		* \brief Computes the least squares conformal map and stores it in
		*  the texture coordinates of the mesh.
		* \details Outline of the algorithm (steps 1,2,3 are not used 
		*   in spetral mode):
		*   - 1) Find an initial solution by projecting on a plane
		*   - 2) Lock two vertices of the mesh
		*   - 3) Copy the initial u,v coordinates to OpenNL
		*   - 4) Construct the LSCM equation with OpenNL
		*   - 5) Solve the equation with OpenNL
		*   - 6) Copy OpenNL solution to the u,v coordinates
		*/
		NLContext context = nlNewContext();
		uint32_t vxmin = UINT32_MAX, vxmax = UINT32_MAX;
		double umin = 1e30, umax = -1e30;
		for (uint32_t i = 0; i < vertexCount; i++) {
			double u = texcoords[i * 2 + 0];
			if (u < umin) {
				vxmin = i;
				umin = u;
			} 
			if (u > umax) {
				vxmax = i;
				umax = u;
			} 
		}
		NLuint nb_vertices = NLuint(vertexCount);
		nlSolverParameteri(NL_NB_VARIABLES, NLint(2*nb_vertices));
		nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
		nlSolverParameteri(NL_MAX_ITERATIONS, NLint(5*nb_vertices));
		nlSolverParameterd(NL_THRESHOLD, 1e-6);	  
		nlBegin(NL_SYSTEM);
		// Copies u,v coordinates from the mesh to OpenNL solver.
		for (uint32_t i = 0; i < vertexCount; i++) {
			nlSetVariable(2 * i    , texcoords[i * 2 + 0]);
			nlSetVariable(2 * i + 1, texcoords[i * 2 + 1]);
			if (i == vxmin || i == vxmax) {
				nlLockVariable(2 * i    );
				nlLockVariable(2 * i + 1);
			} 
		}
		nlBegin(NL_MATRIX);
		// Creates the LSCM equations in OpenNL.
		const uint32_t faceCount = indexCount / 3;
		for (uint32_t f = 0; f < faceCount; f++) {
			const uint32_t v0 = indices[f * 3 + 0];
			const uint32_t v1 = indices[f * 3 + 1];
			const uint32_t v2 = indices[f * 3 + 2];
			const bx::Vec3 p0(positions[v0 * 3 + 0], positions[v0 * 3 + 1], positions[v0 * 3 + 2]);
			const bx::Vec3 p1(positions[v1 * 3 + 0], positions[v1 * 3 + 1], positions[v1 * 3 + 2]);
			const bx::Vec3 p2(positions[v2 * 3 + 0], positions[v2 * 3 + 1], positions[v2 * 3 + 2]);

			vec2 z0,z1,z2;
			project_triangle(p0,p1,p2,z0,z1,z2);
			double a = z1.x - z0.x;
			double b = z1.y - z0.y;
			double c = z2.x - z0.x;
			double d = z2.y - z0.y;
			assert(b == 0.0);

			// Note  : 2*id + 0 --> u
			//         2*id + 1 --> v
			NLuint u0_id = 2*v0    ;
			NLuint v0_id = 2*v0 + 1;
			NLuint u1_id = 2*v1    ;
			NLuint v1_id = 2*v1 + 1;
			NLuint u2_id = 2*v2    ;
			NLuint v2_id = 2*v2 + 1;

			// Note : b = 0

			// Real part
			nlBegin(NL_ROW);
			nlCoefficient(u0_id, -a+c) ;
			nlCoefficient(v0_id,  b-d) ;
			nlCoefficient(u1_id,   -c) ;
			nlCoefficient(v1_id,    d) ;
			nlCoefficient(u2_id,    a);
			nlEnd(NL_ROW);

			// Imaginary part
			nlBegin(NL_ROW);
			nlCoefficient(u0_id, -b+d);
			nlCoefficient(v0_id, -a+c);
			nlCoefficient(u1_id,   -d);
			nlCoefficient(v1_id,   -c);
			nlCoefficient(v2_id,    a);
			nlEnd(NL_ROW);
		}
		nlEnd(NL_MATRIX);
		nlEnd(NL_SYSTEM);
		nlSolve();
		// Copies u,v coordinates from OpenNL solver to the mesh.
		for (uint32_t i = 0; i < vertexCount; i++) {
			double u = nlGetVariable(2 * i    );
			double v = nlGetVariable(2 * i + 1);
			texcoords[i * 2 + 0] = (float)u;
			texcoords[i * 2 + 1] = (float)v;
		}
		nlDeleteContext(context);
	}
#endif
}
#endif

static void atlasGenerateThread()
{
	const objzModel *model = modelGetData();
	const bool firstRun = !s_atlas.data;
	const clock_t startTime = clock();
	if (firstRun) {
		// Create xatlas context on first run only.
		s_atlas.data = xatlas::Create();
		xatlas::SetProgressCallback(s_atlas.data, atlasProgressCallback);
		std::vector<uint8_t> ignoreFaces; // Should be bool, workaround stupid C++ specialization.
		for (uint32_t i = 0; i < model->numObjects; i++) {
			const objzObject &object = model->objects[i];
			auto v = &((const ModelVertex *)model->vertices)[object.firstVertex];
			// Ignore faces with an emissive or transparent material.
			ignoreFaces.resize(object.numIndices / 3);
			memset(ignoreFaces.data(), 0, ignoreFaces.size() * sizeof(uint8_t));
			for (uint32_t j = 0; j < object.numMeshes; j++) {
				const objzMesh &mesh = model->meshes[object.firstMesh + j];
				const objzMaterial *mat = mesh.materialIndex == -1 ? nullptr : &model->materials[mesh.materialIndex];
				if (mat && (mat->emission[0] > 0.0f || mat->emission[1] > 0.0f || mat->emission[2] > 0.0f || mat->opacity < 1.0f)) {
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
			xatlas::AddMeshError::Enum error = xatlas::AddMesh(s_atlas.data, meshDecl, model->numObjects);
			if (error != xatlas::AddMeshError::Success) {
				fprintf(stderr, "Error adding mesh: %s\n", xatlas::StringForEnum(error));
				setErrorMessage("Error adding mesh: %s", xatlas::StringForEnum(error));
				xatlas::Destroy(s_atlas.data);
				s_atlas.data = nullptr;
				s_atlas.status.set(AtlasStatus::NotGenerated);
				return;
			}
			// Destroy context if cancelled while adding meshes.
			if (s_atlas.status.getCancel()) {
				xatlas::Destroy(s_atlas.data);
				s_atlas.data = nullptr;
				s_atlas.status.set(AtlasStatus::NotGenerated);
				s_atlas.status.setCancel(false);
				return;
			}
		}
	}
	if (firstRun || s_atlas.chartOptionsChanged) {
		xatlas::ComputeCharts(s_atlas.data, s_atlas.chartOptions);
		if (s_atlas.status.getCancel()) {
			s_atlas.chartOptionsChanged = true; // Force ComputeCharts to be called next time.
			s_atlas.status.set(AtlasStatus::NotGenerated);
			s_atlas.status.setCancel(false);
			return;
		}
	}
	if (firstRun || s_atlas.chartOptionsChanged || s_atlas.paramOptionsChanged) {
		xatlas::ParameterizeFunc paramFunc = nullptr;
#if USE_LIBIGL || USE_OPENNL
		if (s_atlas.paramMethod != ParamMethod::LSCM)
			paramFunc = atlasParameterizationCallback;
#endif
		xatlas::ParameterizeCharts(s_atlas.data, paramFunc);
		if (s_atlas.status.getCancel()) {
			s_atlas.paramOptionsChanged = true; // Force ParameterizeCharts to be called next time.
			s_atlas.status.set(AtlasStatus::NotGenerated);
			s_atlas.status.setCancel(false);
			return;
		}
	}
	if (firstRun || s_atlas.chartOptionsChanged || s_atlas.paramOptionsChanged || s_atlas.packOptionsChanged) {
		xatlas::PackCharts(s_atlas.data, s_atlas.packOptions);
		if (s_atlas.status.getCancel()) {
			s_atlas.packOptionsChanged = true; // Force PackCharts to be called next time.
			s_atlas.status.set(AtlasStatus::NotGenerated);
			s_atlas.status.setCancel(false);
			return;
		}
	}
	const double elapsedTime = (clock() - startTime) * 1000.0 / CLOCKS_PER_SEC;
	printf("Generated atlas in %.2f seconds (%g ms)\n", elapsedTime / 1000.0, elapsedTime);
	s_atlas.chartOptionsChanged = false;
	s_atlas.paramOptionsChanged = false;
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
			for (uint32_t k = 0; k < chart.faceCount; k++) {
				for (uint32_t l = 0; l < 3; l++) {
					EdgeKey key;
					key.p0 = vertices[mesh.vertexArray[mesh.indexArray[chart.faceArray[k] * 3 + l]].xref].pos;
					key.p1 = vertices[mesh.vertexArray[mesh.indexArray[chart.faceArray[k] * 3 + (l + 1) % 3]].xref].pos;
					edgeMap[key] = 0; // Don't care.
				}
			}
			for (uint32_t k = 0; k < chart.faceCount; k++) {
				for (uint32_t l = 0; l < 3; l++) {
					EdgeKey key;
					key.p0 = vertices[mesh.vertexArray[mesh.indexArray[chart.faceArray[k] * 3 + (l + 1) % 3]].xref].pos;
					key.p1 = vertices[mesh.vertexArray[mesh.indexArray[chart.faceArray[k] * 3 + l]].xref].pos;
					s_atlas.boundaryEdges[numEdges] = edgeMap.count(key) == 0;
					numEdges++;
				}
			}
		}
	}
	// Copy charts for rendering.
	s_atlas.vertices.clear();
	s_atlas.indices.clear();
	s_atlas.chartColorVertices.clear();
	s_atlas.chartInvalidColorVertices.clear();
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
	s_atlas.chartColorVertices.resize(numVertices);
	s_atlas.chartInvalidColorVertices.resize(numVertices);
	s_atlas.chartIndices.resize(numIndices);
	uint32_t firstIndex = 0;
	uint32_t firstChartIndex = 0;
	uint32_t firstVertex = 0;
	numEdges = 0;
	srand(s_atlas.chartColorSeed);
	for (uint32_t i = 0; i < s_atlas.data->meshCount; i++) {
		const xatlas::Mesh &mesh = s_atlas.data->meshes[i];
		const objzObject &object = model->objects[i];
		const ModelVertex *oldVertices = &((const ModelVertex *)model->vertices)[object.firstVertex];
		for (uint32_t j = 0; j < mesh.indexCount; j++)
			s_atlas.indices[firstIndex + j] = firstVertex + mesh.indexArray[j];
		firstIndex += mesh.indexCount;
		for (uint32_t j = 0; j < mesh.chartCount; j++) {
			const xatlas::Chart &chart = mesh.chartArray[j];
			uint8_t bcolor[4];
			randomRGB(bcolor);
			bcolor[3] = 255;
			const uint32_t color = encodeRGBA(bcolor);
			for (uint32_t k = 0; k < chart.faceCount; k++) {
				for (uint32_t l = 0; l < 3; l++) {
					uint32_t &index = s_atlas.chartIndices[firstChartIndex + k * 3 + l];
					index = firstVertex + mesh.indexArray[chart.faceArray[k] * 3 + l];
					s_atlas.chartColorVertices[index] = color;
					if (chart.flags & xatlas::ChartFlags::Invalid)
						s_atlas.chartInvalidColorVertices[index] = color;
					else
						s_atlas.chartInvalidColorVertices[index] = 0xffc0c0c0;
				}
			}
			firstChartIndex += chart.faceCount * 3;
			for (uint32_t k = 0; k < chart.faceCount; k++) {
				bool removeEdge[3];
				for (int l = 0; l < 3; l++) {
					removeEdge[l] = !s_atlas.boundaryEdges[numEdges];
					numEdges++;
				}
				if (removeEdge[0] && removeEdge[1] && removeEdge[2])
					continue;
				WireframeVertex verts[3];
				for (int l = 0; l < 3; l++)
					verts[l].pos = oldVertices[mesh.vertexArray[mesh.indexArray[chart.faceArray[k] * 3 + l]].xref].pos;
				verts[0].barycentric = bx::Vec3(1.0f, 0.0f, 0.0f);
				verts[1].barycentric = bx::Vec3(0.0f, 1.0f, 0.0f);
				verts[2].barycentric = bx::Vec3(0.0f, 0.0f, 1.0f);
				if (removeEdge[0])
					verts[0].barycentric.z = verts[1].barycentric.z = 1.0f;
				if (removeEdge[1])
					verts[1].barycentric.x = verts[2].barycentric.x = 1.0f;
				if (removeEdge[2])
					verts[2].barycentric.y = verts[0].barycentric.y = 1.0f;
				for (int l = 0; l < 3; l++)
					s_atlas.chartBoundaryVertices.push_back(verts[l]);
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
	if (s_atlas.data && !s_atlas.chartOptionsChanged && !s_atlas.paramOptionsChanged && !s_atlas.packOptionsChanged) {
		// Already have an atlas and none of the options that affect atlas creation have changed.
		return;
	}
	bakeClear();
	xatlas::SetPrint(printf, true);
	g_options.shadeMode = ShadeMode::Flat;
	g_options.wireframeMode = WireframeMode::Triangles;
	s_atlas.status.set(AtlasStatus::Generating);
	s_atlas.thread = new std::thread(atlasGenerateThread);
}

static void atlasRenderChartsTextures()
{
	const bgfx::ViewId viewId = kFirstFreeView;
	float bottom = (float)s_atlas.data->height;
	float top = 0.0f;
	if (bgfx::getCaps()->originBottomLeft)
		bx::swap(bottom, top);
	float projection[16];
	bx::mtxOrtho(projection, 0.0f, (float)s_atlas.data->width, bottom, top, 0.0f, 1.0f, 0.0f, bgfx::getCaps()->homogeneousDepth);
	// Setup view for rendering into atlas texture.
	bgfx::setViewClear(viewId, BGFX_CLEAR_COLOR);
	bgfx::setViewRect(viewId, 0, 0, (uint16_t)s_atlas.data->width, (uint16_t)s_atlas.data->height);
	bgfx::setViewFrameBuffer(viewId, s_atlas.chartsFrameBuffer);
	bgfx::setViewTransform(viewId, nullptr, projection);
	// Build chart and chart boundary geometry.
	s_atlas.textureChartIndices.clear();
	s_atlas.textureChartBoundaryVertices.clear();
	srand(s_atlas.chartColorSeed);
	int edge = 0;
	uint32_t firstSourceChartIndex = 0;
	for (uint32_t mi = 0; mi < s_atlas.data->meshCount; mi++) {
		const xatlas::Mesh &mesh = s_atlas.data->meshes[mi];
		for (uint32_t ci = 0; ci < mesh.chartCount; ci++) {
			const xatlas::Chart &chart = mesh.chartArray[ci];
			uint8_t bcolor[4];
			randomRGB(bcolor);
			bcolor[3] = 255;
			const uint32_t color = encodeRGBA(bcolor);
			if (chart.atlasIndex == (uint32_t)s_atlas.currentTexture) {
				for (uint32_t k = 0; k < chart.faceCount * 3; k++)
					s_atlas.textureChartIndices.push_back(s_atlas.chartIndices[firstSourceChartIndex + k]);
			}
			// Build boundary vertices for this chart.
			for (uint32_t k = 0; k < chart.faceCount; k++) {
				for (int l = 0; l < 3; l++) {
					if (chart.atlasIndex == (uint32_t)s_atlas.currentTexture && s_atlas.boundaryEdges[edge]) {
						const xatlas::Vertex &v0 = mesh.vertexArray[mesh.indexArray[chart.faceArray[k] * 3 + l]];
						const xatlas::Vertex &v1 = mesh.vertexArray[mesh.indexArray[chart.faceArray[k] * 3 + (l + 1) % 3]];
						AtlasVertex cbv;
						cbv.pos[0] = v0.uv[0];
						cbv.pos[1] = v0.uv[1];
						cbv.color = color;
						s_atlas.textureChartBoundaryVertices.push_back(cbv);
						cbv.pos[0] = v1.uv[0];
						cbv.pos[1] = v1.uv[1];
						s_atlas.textureChartBoundaryVertices.push_back(cbv);
					}
					edge++;
				}
			}
			firstSourceChartIndex += chart.faceCount * 3;
		}
	}
	// Render charts.
	bgfx::destroyAndClear(s_atlas.textureChartIb);
	s_atlas.textureChartIb = bgfx::createIndexBuffer(bgfx::makeRef(s_atlas.textureChartIndices), BGFX_BUFFER_INDEX32);
	float textureSize_cellSize[4];
	textureSize_cellSize[0] = (float)s_atlas.data->width;
	textureSize_cellSize[1] = (float)s_atlas.data->height;
	textureSize_cellSize[2] = (float)s_atlas.chartTextureCellSize;
	textureSize_cellSize[3] = (float)s_atlas.chartTextureCellSize;
	bgfx::setUniform(s_atlas.u_textureSize_cellSize, textureSize_cellSize);
	bgfx::setState(BGFX_STATE_WRITE_RGB | BGFX_STATE_WRITE_A);
	bgfx::setIndexBuffer(s_atlas.textureChartIb);
	bgfx::setVertexBuffer(0, s_atlas.vb);
	bgfx::setVertexBuffer(1, s_atlas.chartColorVb);
	bgfx::submit(viewId, s_atlas.chartTexcoordSpaceProgram);
	// Render chart boundary lines to emulate conservative rasterization.
	bgfx::destroyAndClear(s_atlas.textureChartBoundaryVb);
	s_atlas.textureChartBoundaryVb = bgfx::createVertexBuffer(bgfx::makeRef(s_atlas.textureChartBoundaryVertices), s_atlas.atlasVertexDecl);
	bgfx::setState(BGFX_STATE_WRITE_RGB | BGFX_STATE_WRITE_A | BGFX_STATE_PT_LINES);
	bgfx::setVertexBuffer(0, s_atlas.textureChartBoundaryVb);
	bgfx::submit(viewId, getColorProgram());
	// Render 4x4 block grid.
	if (s_atlas.showBlockGrid) {
		const uint32_t nVertices = ((s_atlas.data->width + 1) / 4 + (s_atlas.data->height + 1) / 4) * 2;
		if (bgfx::getAvailTransientVertexBuffer(nVertices, s_atlas.atlasVertexDecl) != nVertices)
			return;
		bgfx::TransientVertexBuffer tvb;
		bgfx::allocTransientVertexBuffer(&tvb, nVertices, s_atlas.atlasVertexDecl);
		auto vertices = (AtlasVertex *)tvb.data;
		uint32_t i = 0;
		for (uint32_t x = 0; x <= s_atlas.data->width; x += 4) {
			vertices[i].pos[0] = (float)x + 0.5f;
			vertices[i].pos[1] = 0.5f;
			i++;
			vertices[i].pos[0] = (float)x + 0.5f;
			vertices[i].pos[1] = (float)s_atlas.data->height - 0.5f;
			i++;
		}
		for (uint32_t y = 0; y <= s_atlas.data->height; y += 4) {
			vertices[i].pos[0] = 0.5f;
			vertices[i].pos[1] = (float)y + 0.5f;
			i++;
			vertices[i].pos[0] = (float)s_atlas.data->width - 0.5f;
			vertices[i].pos[1] = (float)y + 0.5f;
			i++;
		}
		for (i = 0; i < nVertices; i++)
			vertices[i].color = 0x77ffffff;
		bgfx::setState(BGFX_STATE_WRITE_RGB | BGFX_STATE_WRITE_A | BGFX_STATE_BLEND_ALPHA | BGFX_STATE_PT_LINES);
		bgfx::setVertexBuffer(0, &tvb);
		bgfx::submit(viewId, getColorProgram(), 1u); // Render last
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
	s_atlas.vb = bgfx::createVertexBuffer(bgfx::makeRef(s_atlas.vertices), ModelVertex::decl);
	s_atlas.ib = bgfx::createIndexBuffer(bgfx::makeRef(s_atlas.indices), BGFX_BUFFER_INDEX32);
	s_atlas.chartColorVb = bgfx::createVertexBuffer(bgfx::makeRef(s_atlas.chartColorVertices), s_atlas.chartColorDecl);
	s_atlas.chartInvalidColorVb = bgfx::createVertexBuffer(bgfx::makeRef(s_atlas.chartInvalidColorVertices), s_atlas.chartColorDecl);
	s_atlas.chartIb = bgfx::createIndexBuffer(bgfx::makeRef(s_atlas.chartIndices), BGFX_BUFFER_INDEX32);
	// Chart boundaries.
	s_atlas.chartBoundaryVb = bgfx::createVertexBuffer(bgfx::makeRef(s_atlas.chartBoundaryVertices), WireframeVertex::decl);
	// Create framebuffer/texture for atlas.
	bgfx::TextureHandle texture = bgfx::createTexture2D((uint16_t)s_atlas.data->width, (uint16_t)s_atlas.data->height, false, 1, bgfx::TextureFormat::BGRA8, BGFX_TEXTURE_RT | BGFX_SAMPLER_UVW_BORDER | BGFX_SAMPLER_BORDER_COLOR(kPaletteBlack));
	s_atlas.chartsFrameBuffer = bgfx::createFrameBuffer(1, &texture, true);
	// Render charts to texture(s) for previewing UVs. Render in UV space.
	s_atlas.currentTexture = 0;
	atlasRenderChartsTextures();
	g_options.shadeMode = ShadeMode::Charts;
	g_options.wireframeMode = WireframeMode::Charts;
	g_options.chartColorMode = ChartColorMode::Individual;
	s_atlas.status.set(AtlasStatus::Ready);
}

void atlasRenderCharts(const float *modelMatrix, uint64_t state)
{
	float textureSize_cellSize[4];
	textureSize_cellSize[0] = (float)s_atlas.data->width;
	textureSize_cellSize[1] = (float)s_atlas.data->height;
	textureSize_cellSize[2] = (float)g_options.chartCellSize;
	textureSize_cellSize[3] = (float)g_options.chartCellSize;
	bgfx::setUniform(s_atlas.u_textureSize_cellSize, textureSize_cellSize);
	bgfx::setState(state);
	bgfx::setTransform(modelMatrix);
	bgfx::setIndexBuffer(s_atlas.chartIb);
	bgfx::setVertexBuffer(0, s_atlas.vb);
	if (g_options.chartColorMode == ChartColorMode::Invalid)
		bgfx::setVertexBuffer(1, s_atlas.chartInvalidColorVb);
	else
		bgfx::setVertexBuffer(1, s_atlas.chartColorVb);
	bgfx::submit(kModelView, s_atlas.chartProgram);
}

void atlasRenderChartsWireframe(const float *modelMatrix)
{
	const float color[] = { 0.0f, 0.0f, 0.0f, 0.75f };
	bgfx::setUniform(s_atlas.u_color, color);
	setWireframeThicknessUniform(1.5f);
	bgfx::setState(BGFX_STATE_WRITE_RGB | BGFX_STATE_WRITE_A | BGFX_STATE_WRITE_Z | BGFX_STATE_DEPTH_TEST_LEQUAL | BGFX_STATE_CULL_CW | BGFX_STATE_BLEND_ALPHA | BGFX_STATE_MSAA);
	bgfx::setTransform(modelMatrix);
	bgfx::setVertexBuffer(0, s_atlas.chartBoundaryVb);
	bgfx::submit(kModelView, getWireframeProgram(), 1);
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
		ImGui::Columns(2, nullptr, false);
		if (s_atlas.data->atlasCount == 1)
			ImGui::Text("%ux%u atlas", s_atlas.data->width, s_atlas.data->height);
		else
			ImGui::Text("%ux%ux%u atlas", s_atlas.data->atlasCount, s_atlas.data->width, s_atlas.data->height);
		ImGui::NextColumn();
		ImGui::Text("%u charts", s_atlas.data->chartCount);
		ImGui::NextColumn();
		ImGui::Text("%u vertices", numVertices);
		ImGui::NextColumn();
		ImGui::Text("%u triangles", numIndices / 3);
		ImGui::NextColumn();
		ImGui::Text("%g texels per unit", s_atlas.data->texelsPerUnit);
		if (s_atlas.data->atlasCount == 1) {
			ImGui::NextColumn();
			ImGui::Text("%g%% utilization", s_atlas.data->utilization[0] * 100.0f);
			ImGui::Columns(1);
		}
		else {
			ImGui::Columns(1);
			for (uint32_t i = 0; i < s_atlas.data->atlasCount; i++)
				ImGui::Text("%u: %g%% utilization", i, s_atlas.data->utilization[i]);
		}
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
#if USE_LIBIGL || USE_OPENNL
	if (ImGui::TreeNodeEx("Parameterization options", ImGuiTreeNodeFlags_FramePadding)) {
		const ParamMethod oldParamMethod = s_atlas.paramMethod;
		ImGui::RadioButton("LSCM", (int *)&s_atlas.paramMethod, (int)ParamMethod::LSCM);
#if USE_LIBIGL
		ImGui::RadioButton("libigl Harmonic", (int *)&s_atlas.paramMethod, (int)ParamMethod::libigl_Harmonic);
		ImGui::RadioButton("libigl LSCM", (int *)&s_atlas.paramMethod, (int)ParamMethod::libigl_LSCM);
		ImGui::RadioButton("libigl ARAP", (int *)&s_atlas.paramMethod, (int)ParamMethod::libigl_ARAP);
#endif
#if USE_OPENNL
		ImGui::RadioButton("OpenNL LSCM", (int *)&s_atlas.paramMethod, (int)ParamMethod::OpenNL_LSCM);
#endif
		ImGui::TreePop();
		if (s_atlas.paramMethod != oldParamMethod)
			s_atlas.paramOptionsChanged = true;
	}
#endif
	if (ImGui::TreeNodeEx("Pack options", ImGuiTreeNodeFlags_FramePadding)) {
		bool changed = false;
		changed |= ImGui::Checkbox("Brute force", &s_atlas.packOptions.bruteForce);
		changed |= ImGui::Checkbox("Block align", &s_atlas.packOptions.blockAlign);
		changed |= ImGui::InputFloat("Texels per unit", &s_atlas.packOptions.texelsPerUnit, 0.0f, 32.0f, 2);
		changed |= ImGui::InputInt("Resolution", (int *)&s_atlas.packOptions.resolution, 8);
		changed |= ImGui::SliderInt("Padding", (int *)&s_atlas.packOptions.padding, 0, 8);
		changed |= ImGui::InputInt("Max chart size", (int *)&s_atlas.packOptions.maxChartSize);
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
	const float margin = 4.0f;
	const ImGuiWindowFlags progressWindowFlags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings;
	const AtlasStatus::Enum atlasStatus = s_atlas.status.get();
	if (atlasStatus == AtlasStatus::Generating) {
		ImGui::SetNextWindowPos(ImVec2(g_windowSize[0] - margin, margin), ImGuiCond_Always, ImVec2(1.0f, 0.0f));
		ImGui::SetNextWindowSize(ImVec2(300.0f, -1.0f), ImGuiCond_Always);
		if (ImGui::Begin("##atlasProgress", nullptr, progressWindowFlags)) {
			int progress;
			xatlas::ProgressCategory::Enum category;
			s_atlas.status.getProgress(&category, &progress);
			ImGui::Text("%s", xatlas::StringForEnum(category));
			for (int i = 0; i < 3; i++) {
				ImGui::SameLine();
				ImGui::Text(i < progressDots ? "." : " ");
			}
			ImGui::ProgressBar(progress / 100.0f);
			if (ImGui::Button("Cancel"))
				s_atlas.status.setCancel(true);
			ImGui::End();
		}
	} else if (atlasStatus == AtlasStatus::Ready && g_options.showAtlasWindow) {
		const float size = 500;
		ImGui::SetNextWindowPos(ImVec2(g_windowSize[0] - size - margin, margin), ImGuiCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(size, size), ImGuiCond_FirstUseEver);
		if (ImGui::Begin("Atlas", &g_options.showAtlasWindow, ImGuiWindowFlags_HorizontalScrollbar)) {
			if (s_atlas.data->atlasCount > 1) {
				ImGui::Text("Atlas %d of %u", s_atlas.currentTexture + 1, s_atlas.data->atlasCount);
				ImGui::SameLine();
				if (ImGui::ArrowButton("##prevAtlas", ImGuiDir_Left)) {
					s_atlas.currentTexture--;
					if (s_atlas.currentTexture < 0)
						s_atlas.currentTexture = s_atlas.data->atlasCount - 1;
					atlasRenderChartsTextures();
				}
				ImGui::SameLine();
				if (ImGui::ArrowButton("##nextAtlas", ImGuiDir_Right)) {
					s_atlas.currentTexture++;
					if (s_atlas.currentTexture > (int)s_atlas.data->atlasCount - 1)
						s_atlas.currentTexture = 0;
					atlasRenderChartsTextures();
				}
				ImGui::SameLine();
			}
			ImGui::Checkbox("Fit to window", &s_atlas.fitToWindow);
			ImGui::SameLine();
			if (ImGui::Checkbox("Show 4x4 grid", &s_atlas.showBlockGrid))
				atlasRenderChartsTextures();
			GuiTexture texture;
			texture.bgfx.handle = bgfx::getTexture(s_atlas.chartsFrameBuffer);
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

const uint32_t *atlasGetImage()
{
	return s_atlas.data->image;
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
