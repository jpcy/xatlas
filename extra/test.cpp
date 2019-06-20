/*
xatlas
https://github.com/jpcy/xatlas
Copyright (c) 2018 Jonathan Young

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#ifdef _MSC_VER
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#endif
#include <stdio.h>
#include <stdarg.h>
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
#define ASSERT(_condition) if (!(_condition)) logf("[FAILED] '%s' %s %d\n", #_condition, __FILE__, __LINE__);

static FILE *logFile = nullptr;

int logf(const char *format, ...)
{
	va_list args;
	va_start(args, format);
	char buffer[2048];
	vsnprintf(buffer, sizeof(buffer), format, args);
	va_end(args);
	if (logFile) {
		fprintf(logFile, "%s", buffer);
		fflush(logFile);
	}
	return printf("%s", buffer);
}

struct AtlasResult {
	uint32_t chartCount;
};

bool generateAtlas(const char *filename, AtlasResult *result)
{
	logf("%s\n", filename);
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	std::string err;
	if (!tinyobj::LoadObj(shapes, materials, err, filename, NULL, tinyobj::triangulation)) {
		logf("   [FAILED]: %s\n", err.c_str());
		return false;
	}
	if (shapes.size() == 0) {
		logf("   [FAILED]: no shapes in obj file\n");
		return false;
	}
	const clock_t start = clock();
	xatlas::Atlas *atlas = xatlas::Create();
	uint32_t totalFaces = 0;
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
			logf("   [FAILED]: Error adding mesh %d '%s': %s\n", i, shapes[i].name.c_str(), xatlas::StringForEnum(error));
			return false;
		}
		totalFaces += meshDecl.indexCount / 3;
	}
	xatlas::Generate(atlas);
	const clock_t end = clock();
	logf("   %g ms\n", (end - start) * 1000.0 / (double)CLOCKS_PER_SEC);
	uint32_t atlasTotalFaces = 0;
	for (uint32_t i = 0; i < atlas->meshCount; i++) {
		const xatlas::Mesh &mesh = atlas->meshes[i];
		atlasTotalFaces += mesh.indexCount / 3;
	}
	ASSERT(atlasTotalFaces == totalFaces);
	if (result)
		result->chartCount = atlas->chartCount;
	xatlas::Destroy(atlas);
	return true;
}

#ifdef _MSC_VER
void processFilesRecursive(const char *path)
{
	WIN32_FIND_DATAA ffd;
	const char lastChar = path[strlen(path) - 1];
	char cleanPath[256];
	sprintf_s(cleanPath, sizeof(cleanPath), "%s%s", path, lastChar == '/' || lastChar == '\\' ? "" : "/");
	char findPath[256];
	sprintf_s(findPath, sizeof(findPath), "%s*", cleanPath);
	HANDLE hFind = FindFirstFileA(findPath, &ffd);
	if (hFind != INVALID_HANDLE_VALUE) {
		do {
			if (strcmp(ffd.cFileName, ".") == 0 || strcmp(ffd.cFileName, "..") == 0)
				continue;
			if (ffd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
				char childPath[256];
				sprintf_s(childPath, sizeof(childPath), "%s%s/", cleanPath, ffd.cFileName);
				processFilesRecursive(childPath);
			} else {
				const char *dot = strrchr(ffd.cFileName, (int)'.');
				if (!dot)
					continue;
				if (_strcmpi(dot + 1, "obj") != 0)
					continue;
				char filename[256];
				sprintf_s(filename, sizeof(filename), "%s%s", cleanPath, ffd.cFileName);
				if (!generateAtlas(filename, nullptr))
					exit(1);
			}
		}
		while (FindNextFileA(hFind, &ffd) != 0);
		FindClose(hFind);
	}
}
#endif

int main(int argc, char **argv)
{
#ifdef _MSC_VER
	if (fopen_s(&logFile, "test.log", "w") != 0)
		logFile = nullptr;
#else
	logFile = fopen("test.log", "w");
#endif
	xatlas::SetPrint(logf, false);
	if (argc > 1) {
		logf("Search path is '%s'\n", argv[1]);
#ifdef _MSC_VER
		processFilesRecursive(argv[1]);
#else
		logf("not implemented\n");
#endif
	} else {
		AtlasResult result;
		if (generateAtlas(MODEL_PATH "cube.obj", &result)) {
			ASSERT(result.chartCount == 6);
		}
		if (generateAtlas(MODEL_PATH "degenerate_edge.obj", &result)) {
			ASSERT(result.chartCount == 1);
		}
		// double sided quad
		if (generateAtlas(MODEL_PATH "double_sided.obj", &result)) {
			ASSERT(result.chartCount == 2);
		}
		if (generateAtlas(MODEL_PATH "duplicate_edge.obj", &result)) {
			ASSERT(result.chartCount == 2);
		}
		if (generateAtlas(MODEL_PATH "gazebo.obj", &result)) {
			ASSERT(result.chartCount == 333);
		}
		if (generateAtlas(MODEL_PATH "zero_area_face.obj", &result)) {
			ASSERT(result.chartCount == 0);
		}
		if (generateAtlas(MODEL_PATH "zero_length_edge.obj", &result)) {
			ASSERT(result.chartCount == 1);
		}
	}
	return 0;
}
