/*
MIT License

Copyright (c) 2018-2020 Jonathan Young

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
/*
example

This example shows how to use the xatlas API to generate a unique set of texture coordinates.

Input: an .obj model file.

Output:
	* an .obj model file (example_output.obj). This is simplistic for example purposes, it doesn't copy materials from the input .obj file.
	* texture coordinates rasterized to images, colored by chart (example_charts*.tga) and by triangle (example_tris*.tga).
*/
#include <assert.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <xatlas_c.h>
#include "objzero/objzero.h"

#ifdef _MSC_VER
#define STRICMP _stricmp
#include <windows.h>
static CRITICAL_SECTION s_progressMutex;
#else
#include <strings.h>
#define STRICMP strcasecmp
#endif

static bool s_verbose = false;

static int Print(const char *format, ...)
{
	va_list arg;
	va_start(arg, format);
	printf("\r"); // Clear progress text.
	const int result = vprintf(format, arg);
	va_end(arg);
	return result;
}

// May be called from any thread.
static bool ProgressCallback(xatlasProgressCategory category, int progress, void *userData)
{
	(void)userData;
	// Don't interupt verbose printing.
	if (s_verbose)
		return true;
#ifdef _MSC_VER
	EnterCriticalSection(&s_progressMutex);
#endif
	printf("\r   %s [", xatlasProgressCategoryString(category));
	for (int i = 0; i < 10; i++)
		printf(progress / ((i + 1) * 10) ? "*" : " ");
	printf("] %d%%", progress);
	fflush(stdout);
	if (progress == 100)
		printf("\n");
#ifdef _MSC_VER
	LeaveCriticalSection(&s_progressMutex);
#endif
	return true;
}

int main(int argc, char *argv[])
{
#ifdef _MSC_VER
	InitializeCriticalSection(&s_progressMutex);
#endif
	if (argc < 2) {
	    printf("Usage: %s input_file.obj [options]\n", argv[0]);
		printf("  Options:\n");
		printf("    -verbose\n");  
	    return 1;
	}
	s_verbose = (argc >= 3 && STRICMP(argv[2], "-verbose") == 0);
	// Load object file.
	printf("Loading '%s'...\n", argv[1]);
	objz_setIndexFormat(OBJZ_INDEX_FORMAT_U32);
	objzModel *model = objz_load(argv[1]);
	if (!model) {
		printf("%s\n", objz_getError());
		return EXIT_FAILURE;
	}
	if (objz_getError()) // Print warnings.
		printf("%s\n", objz_getError());
	// Create empty atlas.
	xatlasSetPrint(Print, s_verbose);
	xatlasAtlas *atlas = xatlasCreate();
	// Set progress callback.
	xatlasSetProgressCallback(atlas, ProgressCallback, NULL);
	// Add mesh to atlas.
	const uint32_t vertexStride = (uint32_t)(sizeof(float) * 8);
	xatlasMeshDecl meshDecl;
	xatlasMeshDeclInit(&meshDecl);
	meshDecl.vertexCount = model->numVertices;
	meshDecl.vertexPositionData = model->vertices;
	meshDecl.vertexPositionStride = vertexStride;
	if (model->flags & OBJZ_FLAG_NORMALS) {
		meshDecl.vertexNormalData = (float *)model->vertices + 5;
		meshDecl.vertexNormalStride = vertexStride;
	}
	if (model->flags & OBJZ_FLAG_TEXCOORDS) {
		meshDecl.vertexUvData = (float *)model->vertices + 3;
		meshDecl.vertexUvStride = vertexStride;
	}
	meshDecl.indexCount = model->numIndices;
	meshDecl.indexData = model->indices;
	meshDecl.indexFormat = XATLAS_INDEX_FORMAT_UINT32;
	xatlasAddMeshError error = xatlasAddMesh(atlas, &meshDecl, 1);
	if (error != XATLAS_ADD_MESH_ERROR_SUCCESS) {
		xatlasDestroy(atlas);
		printf("\rError adding mesh: %s\n", xatlasAddMeshErrorString(error));
		return EXIT_FAILURE;
	}
	xatlasAddMeshJoin(atlas); // Not necessary. Only called here so geometry totals are printed after the AddMesh progress indicator.
	printf("   %u total vertices\n", meshDecl.vertexCount);
	printf("   %u total faces\n", meshDecl.indexCount / 3);
	// Generate atlas.
	printf("Generating atlas\n");
	xatlasGenerate(atlas, NULL, NULL);
	printf("   %d charts\n", atlas->chartCount);
	printf("   %d atlases\n", atlas->atlasCount);
	for (uint32_t i = 0; i < atlas->atlasCount; i++)
		printf("      %d: %0.2f%% utilization\n", i, atlas->utilization[i] * 100.0f);
	printf("   %ux%u resolution\n", atlas->width, atlas->height);
	printf("   %u total vertices\n", atlas->meshes[0].vertexCount);
	// Cleanup.
	xatlasDestroy(atlas);
	objz_destroy(model);
	printf("Done\n");
	return EXIT_SUCCESS;
}
