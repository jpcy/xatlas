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

#define xaPrint
#include "../xatlas.h"

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

struct RasterParam
{
	xatlas::Atlas *atlas;
	uint8_t color[3];
	uint8_t *output_image;
};

static bool RasterCallback(void *param, int x, int y, xatlas::internal::Vector3::Arg bar, xatlas::internal::Vector3::Arg dx, xatlas::internal::Vector3::Arg dy, float coverage)
{
	RasterParam *data = (RasterParam *)param;
	uint8_t *rgba = &data->output_image[(x + y * data->atlas->width) * 4];
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
	std::vector<xatlas::Input_Mesh> inputMeshes;
	inputMeshes.resize(shapes.size());
	int totalVertices = 0, totalFaces = 0;
	for (int i = 0; i < (int)shapes.size(); i++) {
		const tinyobj::mesh_t &objMesh = shapes[i].mesh;
		if (objMesh.normals.size() == 0) {
			printf("Shape %d has no normals\n", i);
			return 0;
		}
		xatlas::Input_Mesh &mesh = inputMeshes[i];
		mesh.vertex_count = objMesh.positions.size() / 3;
		mesh.vertex_array = new xatlas::Input_Vertex[mesh.vertex_count];
		for (int nvert = 0; nvert < mesh.vertex_count; nvert++) {
			for (int j = 0; j < 3; j++) {
				mesh.vertex_array[nvert].position[j] = objMesh.positions[nvert * 3 + j];
				mesh.vertex_array[nvert].normal[j] = objMesh.normals[nvert * 3 + j];
			}
			mesh.vertex_array[nvert].uv[0] = 0;
			mesh.vertex_array[nvert].uv[1] = 0;
			mesh.vertex_array[nvert].first_colocal = nvert;
		}
		mesh.face_count = objMesh.indices.size() / 3;
		mesh.face_array = new xatlas::Input_Face[mesh.face_count];
		for (int nface = 0; nface < mesh.face_count; nface++) {
			mesh.face_array[nface].material_index = 0;
			mesh.face_array[nface].vertex_index[0] = objMesh.indices[nface * 3];
			mesh.face_array[nface].vertex_index[1] = objMesh.indices[nface * 3 + 1];
			mesh.face_array[nface].vertex_index[2] = objMesh.indices[nface * 3 + 2];
		}
		if (verbose)
			printf("      shape %d: %d vertices, %d triangles\n", i, mesh.vertex_count, mesh.face_count);
		xatlas::add_mesh(&mesh);
		totalVertices += mesh.vertex_count;
		totalFaces += mesh.face_count;
	}
	printf("   %d vertices\n", totalVertices);
	printf("   %d triangles\n", totalFaces);
	// Generate Output_Mesh.
	xatlas::Options atlas_options;
	xatlas::set_default_options(&atlas_options);
	if (!verbose)
		atlas_options.Print = NULL;
	// Avoid brute force packing, since it can be unusably slow in some situations.
	atlas_options.packer.packing_quality = 1;
	atlas_options.packer.texel_area = 256;
	atlas_options.packer.conservative = true;
	atlas_options.packer.padding = 2;
	Stopwatch stopwatch;
	printf("Generating atlas...\n");
	xatlas::Atlas atlas = xatlas::atlas_generate(&atlas_options);
	const double elapsedMs = stopwatch.elapsed();
	if (atlas.error != xatlas::Error_Success)
	{
		printf("Error");
		if (atlas.error == xatlas::Error_Invalid_Args)
			printf(": invalid arguments");
		else if (atlas.error == xatlas::Error_Invalid_Options)
			printf(": invalid options");
		else if (atlas.error == xatlas::Error_Invalid_Mesh)
			printf(": invalid mesh, index %d", atlas.errorMeshIndex);
		else if (atlas.error == xatlas::Error_Invalid_Mesh_Non_Manifold)
			printf(": non-manifold mesh, index %d", atlas.errorMeshIndex);
		printf("\n");
		// Free meshes.
		for (int i = 0; i < (int)inputMeshes.size(); i++) {
			delete [] inputMeshes[i].face_array;
			delete [] inputMeshes[i].vertex_array;
		}
		return 0;
	}
	printf("   %.2f seconds elapsed (%g milliseconds)\n", elapsedMs / 1000.0, elapsedMs);
	printf("   %d charts\n", atlas.nCharts);
	printf("   %dx%d resolution\n", atlas.width, atlas.height);
	uint8_t *output_image = new uint8_t[atlas.width * atlas.height * 4];
	memset(output_image, 0, atlas.width * atlas.height * 4);
	for (int i = 0; i < atlas.nMeshes; i++) {
		xatlas::Output_Mesh *output_mesh = atlas.meshes[i];
		if (verbose)
			printf("   output mesh %d: %d vertices, %d triangles\n", i, output_mesh->vertex_count, output_mesh->index_count / 3);
		for (int j = 0; j < output_mesh->index_count; j += 3) {
			const xatlas::Output_Vertex *v[3];
			v[0] = &output_mesh->vertex_array[output_mesh->index_array[j + 0]];
			v[1] = &output_mesh->vertex_array[output_mesh->index_array[j + 1]];
			v[2] = &output_mesh->vertex_array[output_mesh->index_array[j + 2]];
			xatlas::internal::Vector2 verts[3];
			for (int k = 0; k < 3; k++) {
				verts[k] = xatlas::internal::Vector2(v[k]->uv[0], v[k]->uv[1]);
			}
			RasterParam raster;
			raster.atlas = &atlas;
			raster.color[0] = rand() % 255;
			raster.color[1] = rand() % 255;
			raster.color[2] = rand() % 255;
			raster.output_image = output_image;
			xatlas::internal::raster::drawTriangle(xatlas::internal::raster::Mode_Nearest, xatlas::internal::Vector2((float)atlas.width, (float)atlas.height), true, verts, RasterCallback, &raster);
		}
	}
	const char *outputFilename = "output.tga";
	printf("Writing '%s'...\n", outputFilename);
	stbi_write_tga(outputFilename, atlas.width, atlas.height, 4, output_image);
	delete [] output_image;
	// Free meshes.
	for (int i = 0; i < (int)inputMeshes.size(); i++) {
		delete [] inputMeshes[i].face_array;
		delete [] inputMeshes[i].vertex_array;
	}
	xatlas::atlas_free(atlas);
	printf("Done\n");
	return 0;
}
