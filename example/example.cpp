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
#include "tiny_obj_loader.h"

#include "../xatlas.h"
#include "../xatlas_raster.h"

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

struct RasterParam
{
	uint8_t color[3];
	uint8_t *imageData;
	uint32_t imageWidth;
};

static bool RasterCallback(void *param, int x, int y, xatlas::internal::Vector3::Arg bar, xatlas::internal::Vector3::Arg dx, xatlas::internal::Vector3::Arg dy, float coverage)
{
	RasterParam *data = (RasterParam *)param;
	uint8_t *rgba = &data->imageData[(x + y * data->imageWidth) * 4];
	rgba[0] = data->color[0];
	rgba[1] = data->color[1];
	rgba[2] = data->color[2];
	rgba[3] = 255;
	return true;
}

int main(int argc, char *argv[])
{
	if (argc < 2) {
	    printf("Usage: %s input_file.obj [options]\n", argv[0]);
		printf("  Options:\n");
		printf("    -verbose\n");  
	    return 1;
	}
	const bool verbose = (argc >= 3 && _stricmp(argv[2], "-verbose") == 0);
	// Load object file.
	printf("Loading '%s'...\n", argv[1]);
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	std::string err;
	bool ret = tinyobj::LoadObj(shapes, materials, err, argv[1]);
	if (!ret) {
		printf("Error: %s\n", err.c_str());
		return 0;
	}
	if (shapes.size() == 0) {
		printf("Error: no shapes in obj file\n");
		return 0;
	}
	printf("   %lu shapes\n", shapes.size());
	// Create atlas.
	if (verbose)
		xatlas::SetPrint(Print);
	xatlas::Options options;
	options.packer.resolution = 1024;
	options.packer.conservative = true;
	options.packer.padding = 1;
	xatlas::Atlas *atlas = xatlas::Create(options);
	// Add meshes to atlas.
	Stopwatch stopwatch;
	uint32_t totalVertices = 0, totalFaces = 0;
	for (int i = 0; i < (int)shapes.size(); i++) {
		const tinyobj::mesh_t &objMesh = shapes[i].mesh;
		xatlas::InputMesh mesh;
		mesh.vertexCount = (int)objMesh.positions.size() / 3;
		mesh.vertexPositionData = objMesh.positions.data();
		mesh.vertexPositionStride = sizeof(float) * 3;
		if (objMesh.normals.size() == 0) {
			mesh.vertexNormalData = NULL;
		} else {
			mesh.vertexNormalData = objMesh.normals.data();
			mesh.vertexNormalStride = sizeof(float) * 3;
		}
		if (objMesh.texcoords.size() == 0) {
			mesh.vertexUvData = NULL;
		} else {
			mesh.vertexUvData = objMesh.texcoords.data();
			mesh.vertexUvStride = sizeof(float) * 2;
		}
		mesh.indexCount = (int)objMesh.indices.size();
		mesh.indexData = objMesh.indices.data();
		mesh.indexFormat = xatlas::IndexFormat::Float;
		mesh.faceMaterialData = NULL;
		if (verbose)
			printf("      shape %d: %u vertices, %u triangles\n", i, mesh.vertexCount, mesh.indexCount / 3);
		xatlas::AddMeshError error = xatlas::AddMesh(atlas, mesh);
		if (error.code != xatlas::AddMeshErrorCode::Success) {
			printf("Error adding mesh %d '%s': %s\n", i, shapes[i].name.c_str(), xatlas::StringForEnum(error.code));
			switch (error.code) {
			case xatlas::AddMeshErrorCode::AlreadyAddedEdge:
			case xatlas::AddMeshErrorCode::DegenerateColocalEdge:
			case xatlas::AddMeshErrorCode::DegenerateEdge:
			case xatlas::AddMeshErrorCode::DuplicateEdge:
			case xatlas::AddMeshErrorCode::ZeroLengthEdge:
				printf("   indices: %u %u\n", error.index0, error.index1);
				for (int j = 0; j < 2; j++) {
					const float *pos = &objMesh.positions[(j == 0 ? error.index0 : error.index1) * 3];
					printf("   position %d: %g %g %g\n", j + 1, pos[0], pos[1], pos[2]);
				}
				break;
			case xatlas::AddMeshErrorCode::IndexOutOfRange:
				printf("   index: %u\n", error.index0);
				break;
			case xatlas::AddMeshErrorCode::ZeroAreaFace: {
					const uint32_t *indices = &objMesh.indices[error.face * 3];
					printf("   face: %u\n", error.face);
					printf("   indices: %u %u %u\n", indices[0], indices[1], indices[2]);
					for (int j = 0; j < 3; j++) {
						const float *pos = &objMesh.positions[indices[j] * 3];
						printf("   position %d: %g %g %g\n", j + 1, pos[0], pos[1], pos[2]);
					}
				}
				break;
			}
			return 0;
		}
		totalVertices += mesh.vertexCount;
		totalFaces += mesh.indexCount / 3;
	}
	printf("   %u vertices\n", totalVertices);
	printf("   %u triangles\n", totalFaces);
	double elapsedMs = stopwatch.elapsed();
	printf("   %.2f seconds elapsed (%g milliseconds)\n", elapsedMs / 1000.0, elapsedMs);
	// Generate output meshes.
	printf("Generating atlas...\n");
	stopwatch.reset();
	xatlas::Generate(atlas);
	elapsedMs = stopwatch.elapsed();
	printf("   %.2f seconds elapsed (%g milliseconds)\n", elapsedMs / 1000.0, elapsedMs);
	printf("   %d charts\n", xatlas::GetNumCharts(atlas));
	const uint32_t width = xatlas::GetWidth(atlas);
	const uint32_t height = xatlas::GetHeight(atlas);
	printf("   %ux%u resolution\n", width, height);
	// Dump images.
	std::vector<uint8_t> outputTrisImage, outputChartsImage;
	outputTrisImage.resize(width * height * 4);
	outputChartsImage.resize(width * height * 4);
	for (int i = 0; i < (int)shapes.size(); i++) {
		const xatlas::OutputMesh *mesh = xatlas::GetOutputMeshes(atlas)[i];
		if (verbose)
			printf("   output mesh %d: %u vertices, %u triangles, %u charts\n", i, mesh->vertexCount, mesh->indexCount / 3, mesh->chartCount);
		// Rasterize mesh triangles.
		for (uint32_t j = 0; j < mesh->indexCount; j += 3) {
			RasterParam raster;
			xatlas::internal::Vector2 verts[3];
			for (int k = 0; k < 3; k++) {
				const xatlas::OutputVertex &v = mesh->vertexArray[mesh->indexArray[j + k]];
				verts[k] = xatlas::internal::Vector2(v.uv[0], v.uv[1]);
				raster.color[k] = rand() % 255;
			}
			raster.imageData = outputTrisImage.data();
			raster.imageWidth = width;
			xatlas::internal::raster::drawTriangle(xatlas::internal::raster::Mode_Antialiased, xatlas::internal::Vector2((float)width, (float)height), true, verts, RasterCallback, &raster);
		}
		// Rasterize mesh charts.
		for (uint32_t j = 0; j < mesh->chartCount; j++) {
			const xatlas::OutputChart *chart = &mesh->chartArray[j];
			uint8_t color[3];
			color[0] = rand() % 255;
			color[1] = rand() % 255;
			color[2] = rand() % 255;
			for (uint32_t k = 0; k < chart->indexCount; k += 3) {
				RasterParam raster;
				xatlas::internal::Vector2 verts[3];
				for (int l = 0; l < 3; l++) {
					const xatlas::OutputVertex &v = mesh->vertexArray[chart->indexArray[k + l]];
					verts[l] = xatlas::internal::Vector2(v.uv[0], v.uv[1]);
					raster.color[l] = color[l];
				}
				raster.imageData = outputChartsImage.data();
				raster.imageWidth = width;
				xatlas::internal::raster::drawTriangle(xatlas::internal::raster::Mode_Antialiased, xatlas::internal::Vector2((float)width, (float)height), true, verts, RasterCallback, &raster);
			}
		}
	}
	const char *outputTrisFilename = "output_tris.tga";
	printf("Writing '%s'...\n", outputTrisFilename);
	stbi_write_tga(outputTrisFilename, width, height, 4, outputTrisImage.data());
	const char *outputChartsFilename = "output_charts.tga";
	printf("Writing '%s'...\n", outputChartsFilename);
	stbi_write_tga(outputChartsFilename, width, height, 4, outputChartsImage.data());
	// Cleanup.
	xatlas::Destroy(atlas);
	printf("Done\n");
	return 0;
}
