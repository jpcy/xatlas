#include <chrono>
#include <stdio.h>
#include <assert.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4996)
#endif
#include "stb_image_write.h"
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#define TINYOBJLOADER_IMPLEMENTATION
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4201)
#pragma warning(disable : 4706)
#endif
#include "tiny_obj_loader.h"
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#include "../xatlas.h"

#ifdef _MSC_VER
#define STRICMP _stricmp
#else
#include <strings.h>
#define STRICMP strcasecmp
#endif

class Stopwatch
{
	typedef std::chrono::high_resolution_clock Clock;
public:
	Stopwatch() { reset(); }
	void reset() { start_ = Clock::now(); }
	double elapsed() const { return (double)std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - start_).count(); }
private:
	std::chrono::time_point<Clock> start_;
};

static void Print(const char *format, ...)
{
	va_list arg;
	va_start(arg, format);
	vprintf(format, arg);
	va_end(arg);
}

static void PrintProgress(const char *name, const char *indent1, const char *indent2, int progress, Stopwatch *stopwatch)
{
	if (progress == 0)
		stopwatch->reset();
	printf("\r%s%s [", indent1, name);
	for (int i = 0; i < 10; i++)
		printf(progress / ((i + 1) * 10) ? "*" : " ");
	printf("] %d%%", progress);
	fflush(stdout);
	if (progress == 100)
		printf("\n%s%.2f seconds (%g ms) elapsed\n", indent2, stopwatch->elapsed() / 1000.0, stopwatch->elapsed());
}

static void ProgressCallback(xatlas::ProgressCategory::Enum category, int progress, void *userData)
{
	Stopwatch *stopwatch = (Stopwatch *)userData;
	PrintProgress(xatlas::StringForEnum(category), "   ", "      ", progress, stopwatch);
}

static void SetPixel(uint8_t *dest, int destWidth, int x, int y, const uint8_t *color)
{
	uint8_t *pixel = &dest[x * 3 + y * (destWidth * 3)];
	pixel[0] = color[0];
	pixel[1] = color[1];
	pixel[2] = color[2];
}

// https://github.com/miloyip/line/blob/master/line_bresenham.c
static void RasterizeLine(uint8_t *dest, int destWidth, const int *p1, const int *p2, const uint8_t *color)
{
	const int dx = abs(p2[0] - p1[0]), sx = p1[0] < p2[0] ? 1 : -1;
	const int dy = abs(p2[1] - p1[1]), sy = p1[1] < p2[1] ? 1 : -1;
	int err = (dx > dy ? dx : -dy) / 2;
	int current[2];
	current[0] = p1[0];
	current[1] = p1[1];
	while (SetPixel(dest, destWidth, current[0], current[1], color), current[0] != p2[0] || current[1] != p2[1])
	{
		const int e2 = err;
		if (e2 > -dx) { err -= dy; current[0] += sx; }
		if (e2 < dy) { err += dx; current[1] += sy; }
	}
}

// https://github.com/ssloy/tinyrenderer/wiki/Lesson-2:-Triangle-rasterization-and-back-face-culling
static void RasterizeTriangle(uint8_t *dest, int destWidth, const int *t0, const int *t1, const int *t2, const uint8_t *color)
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
		if (A[0] > B[0]) std::swap(A, B);
		for (int j = A[0]; j <= B[0]; j++)
			SetPixel(dest, destWidth, j, t0[1] + i, color);
	}
}

int main(int argc, char *argv[])
{
	if (argc < 2) {
	    printf("Usage: %s input_file.obj [options]\n", argv[0]);
		printf("  Options:\n");
		printf("    -verbose\n");  
	    return 1;
	}
	const bool verbose = (argc >= 3 && STRICMP(argv[2], "-verbose") == 0);
	// Load object file.
	printf("Loading '%s'...\n", argv[1]);
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	std::string err;
	if (!tinyobj::LoadObj(shapes, materials, err, argv[1], NULL, tinyobj::triangulation)) {
		printf("Error: %s\n", err.c_str());
		return 0;
	}
	if (shapes.size() == 0) {
		printf("Error: no shapes in obj file\n");
		return 0;
	}
	printf("   %d shapes\n", (int)shapes.size());
	// Create atlas.
	if (verbose)
		xatlas::SetPrint(Print);
	xatlas::Atlas *atlas = xatlas::Create();
	// Add meshes to atlas.
	Stopwatch stopwatch;
	int progress = 0;
	if (verbose)
		printf("Adding meshes...\n");
	else
		PrintProgress("Adding meshes", "", "   ", 0, &stopwatch);
	uint32_t totalVertices = 0, totalFaces = 0;
	for (int i = 0; i < (int)shapes.size(); i++) {
		const tinyobj::mesh_t &objMesh = shapes[i].mesh;
		xatlas::InputMesh mesh;
		mesh.vertexCount = (int)objMesh.positions.size() / 3;
		mesh.vertexPositionData = objMesh.positions.data();
		mesh.vertexPositionStride = sizeof(float) * 3;
		if (!objMesh.normals.empty()) {
			mesh.vertexNormalData = objMesh.normals.data();
			mesh.vertexNormalStride = sizeof(float) * 3;
		}
		if (!objMesh.texcoords.empty()) {
			mesh.vertexUvData = objMesh.texcoords.data();
			mesh.vertexUvStride = sizeof(float) * 2;
		}
		mesh.indexCount = (int)objMesh.indices.size();
		mesh.indexData = objMesh.indices.data();
		mesh.indexFormat = xatlas::IndexFormat::UInt32;
		if (verbose)
			printf("   mesh %d: %u vertices, %u triangles\n", i, mesh.vertexCount, mesh.indexCount / 3);
		xatlas::AddMeshError::Enum error = xatlas::AddMesh(atlas, mesh);
		if (error != xatlas::AddMeshError::Success) {
			if (!verbose)
				printf("\n");
			printf("Error adding mesh %d '%s': %s\n", i, shapes[i].name.c_str(), xatlas::StringForEnum(error));
			return EXIT_FAILURE;
		}
		totalVertices += mesh.vertexCount;
		totalFaces += mesh.indexCount / 3;
		if (!verbose)
		{
			const int newProgress = int((i + 1) / (float)shapes.size() * 100.0f);
			if (newProgress != progress)
			{
				progress = newProgress;
				PrintProgress("Adding meshes", "", "   ", progress, &stopwatch);
			}
		}
	}
	if (!verbose && progress != 100)
		PrintProgress("Adding meshes", "", "   ", 100, &stopwatch);
	printf("   %u total vertices\n", totalVertices);
	printf("   %u total triangles\n", totalFaces);
	// Generate output meshes.
	printf("Generating atlas\n");
	xatlas::PackerOptions packerOptions;
	packerOptions.resolution = 1024;
	packerOptions.conservative = true;
	packerOptions.padding = 1;
	xatlas::GenerateCharts(atlas, xatlas::CharterOptions(), verbose ? NULL : ProgressCallback, &stopwatch);
	xatlas::PackCharts(atlas, packerOptions, verbose ? NULL : ProgressCallback, &stopwatch);
	printf("   %d charts\n", xatlas::GetNumCharts(atlas));
	const uint32_t width = xatlas::GetWidth(atlas);
	const uint32_t height = xatlas::GetHeight(atlas);
	printf("   %ux%u resolution\n", width, height);
	totalVertices = 0;
	totalFaces = 0;
	for (int i = 0; i < (int)shapes.size(); i++) {
		const xatlas::OutputMesh *mesh = xatlas::GetOutputMeshes(atlas)[i];
		totalVertices += mesh->vertexCount;
		totalFaces += mesh->indexCount / 3;
	}
	printf("   %u total vertices\n", totalVertices);
	printf("   %u total triangles\n", totalFaces);
	if (verbose) {
		for (int i = 0; i < (int)shapes.size(); i++) {
			const xatlas::OutputMesh *mesh = xatlas::GetOutputMeshes(atlas)[i];
			printf("   output mesh %d: %u vertices, %u triangles, %u charts\n", i, mesh->vertexCount, mesh->indexCount / 3, mesh->chartCount);
		}
	}
	printf("Rasterizing result...\n");
	if (width > 0 && height > 0) {
		// Dump images.
		std::vector<uint8_t> outputTrisImage, outputChartsImage;
		outputTrisImage.resize(width * height * 3);
		outputChartsImage.resize(width * height * 3);
		for (int i = 0; i < (int)shapes.size(); i++) {
			const xatlas::OutputMesh *mesh = xatlas::GetOutputMeshes(atlas)[i];
			// Rasterize mesh triangles.
			const uint8_t white[] = { 255, 255, 255 };
			for (uint32_t j = 0; j < mesh->indexCount; j += 3) {
				int verts[3][2];
				uint8_t color[3];
				for (int k = 0; k < 3; k++) {
					const xatlas::OutputVertex &v = mesh->vertexArray[mesh->indexArray[j + k]];
					verts[k][0] = int(v.uv[0]);
					verts[k][1] = int(v.uv[1]);
					color[k] = rand() % 255;
				}
				if (!verts[0][0] && !verts[0][1] && !verts[1][0] && !verts[1][1] && !verts[2][0] && !verts[2][1])
					continue; // Skip triangles that weren't atlased.
				RasterizeTriangle(outputTrisImage.data(), width, verts[0], verts[1], verts[2], color);
				RasterizeLine(outputTrisImage.data(), width, verts[0], verts[1], white);
				RasterizeLine(outputTrisImage.data(), width, verts[1], verts[2], white);
				RasterizeLine(outputTrisImage.data(), width, verts[2], verts[0], white);
			}
			// Rasterize mesh charts.
			for (uint32_t j = 0; j < mesh->chartCount; j++) {
				const xatlas::OutputChart *chart = &mesh->chartArray[j];
				uint8_t color[3];
				color[0] = rand() % 255;
				color[1] = rand() % 255;
				color[2] = rand() % 255;
				for (uint32_t k = 0; k < chart->indexCount; k += 3) {
					int verts[3][2];
					for (int l = 0; l < 3; l++) {
						const xatlas::OutputVertex &v = mesh->vertexArray[chart->indexArray[k + l]];
						verts[l][0] = int(v.uv[0]);
						verts[l][1] = int(v.uv[1]);
					}
					RasterizeTriangle(outputChartsImage.data(), width, verts[0], verts[1], verts[2], color);
					RasterizeLine(outputChartsImage.data(), width, verts[0], verts[1], white);
					RasterizeLine(outputChartsImage.data(), width, verts[1], verts[2], white);
					RasterizeLine(outputChartsImage.data(), width, verts[2], verts[0], white);
				}
			}
		}
		const char *outputTrisFilename = "output_tris.tga";
		printf("Writing '%s'...\n", outputTrisFilename);
		stbi_write_tga(outputTrisFilename, width, height, 3, outputTrisImage.data());
		const char *outputChartsFilename = "output_charts.tga";
		printf("Writing '%s'...\n", outputChartsFilename);
		stbi_write_tga(outputChartsFilename, width, height, 3, outputChartsImage.data());
	}
	// Cleanup.
	xatlas::Destroy(atlas);
	printf("Done\n");
	return 0;
}
