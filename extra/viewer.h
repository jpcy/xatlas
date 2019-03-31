/*
xatlas
https://github.com/jpcy/xatlas
Copyright (c) 2018 Jonathan Young

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#pragma once
#include <bx/bx.h>
#include <bx/math.h>
#include <bgfx/bgfx.h>
#include <objzero/objzero.h>

constexpr bgfx::ViewId kModelView = 0;
constexpr bgfx::ViewId kGuiView = 1;
constexpr bgfx::ViewId kFirstFreeView = 2;

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

void atlasInit();
void atlasShutdown();
void atlasDestroy();
void atlasGenerate();
void atlasFinalize();
void atlasRenderCharts(const float *modelMatrix);
void atlasRenderChartsWireframe(const float *modelMatrix);
void atlasShowGuiOptions();
void atlasShowGuiWindow(int progressDots);
uint32_t atlasGetCount();
uint32_t atlasGetWidth();
uint32_t atlasGetHeight();
bgfx::VertexBufferHandle atlasGetVb();
bgfx::IndexBufferHandle atlasGetIb();
bool atlasIsNotGenerated();
bool atlasIsReady();

void bakeInit();
void bakeShutdown();
void bakeExecute();
void bakeFrame(uint32_t bgfxFrame);
void bakeClear();
void bakeShowGuiOptions();
void bakeShowGuiWindow();
bgfx::TextureHandle bakeGetLightmap();
uint32_t bakeGetLightmapSamplerFlags();
bool bakeIsIdle();

struct GuiTextureFlags
{
	enum
	{
		PointSampler = 1 << 0,
	};
};

union GuiTexture
{
	ImTextureID imgui;
	struct
	{
		bgfx::TextureHandle handle;
		uint16_t flags;
	}
	bgfx;
};

void guiInit();
void guiShutdown();
void guiResize(int width, int height);
void guiRunFrame(float deltaTime);
void guiRender();

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

void modelInit();
void modelShutdown();
void modelFinalize();
void modelOpenDialog();
void modelDestroy();
void modelSetMaterialTexturesAndUniforms(int32_t materialIndex);
void modelRender(const float *view, const float *projection);
void modelShowGuiOptions();
void modelShowGuiWindow(int progressDots);
AABB modelGetAABB();
const objzModel *modelGetData();
bx::Vec3 modelGetCentroid();
bgfx::ShaderHandle modelGet_vs_model();
bool modelIsLoaded();

enum class ShadeMode
{
	Flat,
	Charts,
	Lightmap
};

enum class WireframeMode
{
	Charts,
	Triangles
};

struct Options
{
	bool gui = true;
	bool wireframe = true;
	ShadeMode shadeMode = ShadeMode::Flat;
	WireframeMode wireframeMode = WireframeMode::Triangles;
	int chartCellSize = 1;
	bool lightmapPointSampling = false;
	bool showAtlasWindow = true;
	bool showLightmapWindow = true;
};

extern Options g_options;
struct GLFWwindow;
extern GLFWwindow *g_window;

void randomRGB(uint8_t *color);
void setErrorMessage(const char *format, ...);
void resetCamera();

enum class ShaderId
{
	fs_atomicCounterClear,
	fs_chart,
	fs_color,
	fs_gui,
	fs_lightmapAverage,
	fs_lightmapOp,
	fs_material,
	fs_rayBundleClear,
	fs_rayBundleIntegrate,
	fs_rayBundleLightmapClear,
	fs_rayBundleWrite,
	vs_chart,
	vs_chartTexcoordSpace,
	vs_gui,
	vs_model,
	vs_position
};

bgfx::ShaderHandle loadShader(ShaderId id);
bgfx::ShaderHandle get_fs_color();
bgfx::ShaderHandle get_vs_position();
bgfx::ProgramHandle getColorProgram();

struct PosVertex
{
	float pos[2];
	static bgfx::VertexDecl decl;

	static void init()
	{
		decl.begin()
			.add(bgfx::Attrib::Position, 2, bgfx::AttribType::Float)
			.end();
		assert(decl.getStride() == sizeof(PosVertex));
	}
};

#define WINDOW_TITLE "xatlas viewer"
