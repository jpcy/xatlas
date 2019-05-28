/*
xatlas
https://github.com/jpcy/xatlas
Copyright (c) 2018 Jonathan Young

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#include <stdio.h>
#include <time.h>

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4201)
#endif
#include <tiny_obj_loader.h>
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#include "../xatlas.h"

#define MODEL_PATH "../../models/"
#define ASSERT(_condition) if (!(_condition)) printf("[FAILED] '%s' %s %d\n", #_condition, __FILE__, __LINE__);

struct AtlasResult {
	uint32_t chartCount;
};

bool generateAtlas(const char *name, AtlasResult *result)
{
	char filename[256];
	snprintf(filename, sizeof(filename), MODEL_PATH "%s.obj", name);
	printf("%s", name);
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	std::string err;
	if (!tinyobj::LoadObj(shapes, materials, err, filename, NULL, tinyobj::triangulation)) {
		printf("\n[FAILED]: %s\n", err.c_str());
		return false;
	}
	if (shapes.size() == 0) {
		printf("\n[FAILED]: no shapes in obj file\n");
		return false;
	}
	const clock_t start = clock();
	xatlas::Atlas *atlas = xatlas::Create();
	uint32_t totalVertices = 0, totalFaces = 0;
	for (int i = 0; i < (int)shapes.size(); i++) {
		const tinyobj::mesh_t &objMesh = shapes[i].mesh;
		xatlas::MeshDecl meshDecl;
		meshDecl.vertexCount = (int)objMesh.positions.size() / 3;
		meshDecl.vertexPositionData = objMesh.positions.data();
		meshDecl.vertexPositionStride = sizeof(float) * 3;
		if (!objMesh.normals.empty()) {
			meshDecl.vertexNormalData = objMesh.normals.data();
			meshDecl.vertexNormalStride = sizeof(float) * 3;
		}
		if (!objMesh.texcoords.empty()) {
			meshDecl.vertexUvData = objMesh.texcoords.data();
			meshDecl.vertexUvStride = sizeof(float) * 2;
		}
		meshDecl.indexCount = (int)objMesh.indices.size();
		meshDecl.indexData = objMesh.indices.data();
		meshDecl.indexFormat = xatlas::IndexFormat::UInt32;
		xatlas::AddMeshError::Enum error = xatlas::AddMesh(atlas, meshDecl);
		if (error != xatlas::AddMeshError::Success) {
			xatlas::Destroy(atlas);
			printf("\n[FAILED]: Error adding mesh %d '%s': %s\n", i, shapes[i].name.c_str(), xatlas::StringForEnum(error));
			return false;
		}
		totalVertices += meshDecl.vertexCount;
		totalFaces += meshDecl.indexCount / 3;
	}
	xatlas::Generate(atlas);
	const clock_t end = clock();
	printf(" [%g ms]\n", (end - start) * 1000.0 / (double)CLOCKS_PER_SEC);
	uint32_t atlasTotalVertices = 0, atlasTotalFaces = 0;
	for (uint32_t i = 0; i < atlas->meshCount; i++) {
		const xatlas::Mesh &mesh = atlas->meshes[i];
		atlasTotalVertices += mesh.vertexCount;
		atlasTotalFaces += mesh.indexCount / 3;
	}
	ASSERT(atlasTotalVertices == totalVertices);
	ASSERT(atlasTotalFaces == totalFaces);
	result->chartCount = atlas->chartCount;
	xatlas::Destroy(atlas);
	return true;
}

int main(int /*argc*/, char ** /*argv*/)
{
	AtlasResult result;
	if (generateAtlas("cube", &result)) {
		ASSERT(result.chartCount == 6);
	}
	if (generateAtlas("degenerate_edge", &result)) {
		ASSERT(result.chartCount == 1);
	}
	// double sided quad
	if (generateAtlas("double_sided", &result)) {
		ASSERT(result.chartCount == 2);
	}
	if (generateAtlas("duplicate_edge", &result)) {
		ASSERT(result.chartCount == 2);
	}
	if (generateAtlas("gazebo", &result)) {
		ASSERT(result.chartCount == 333);
	}
	if (generateAtlas("zero_area_face", &result)) {
		ASSERT(result.chartCount == 0);
	}
	if (generateAtlas("zero_length_edge", &result)) {
		ASSERT(result.chartCount == 1);
	}
	return 0;
}
