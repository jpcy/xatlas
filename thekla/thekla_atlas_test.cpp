#include <stdio.h>
#include <assert.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../extern/stb_image_write.h"

#define TINYOBJLOADER_IMPLEMENTATION
#include "../extern/tinyobj/tiny_obj_loader.h"

#define XATLAS_IMPLEMENTATION
#include "../xatlas.h"

struct Obj_Vertex {
	float position[3];
	float normal[3];
	float uv[2];
	int first_colocal;
};

struct Obj_Face {
	int vertex_index[3];
	int material_index;
};

struct Obj_Material {
	// @@ Read obj mtl parameters as is.
};

struct Obj_Mesh {
	int vertex_count;
	Obj_Vertex *vertex_array;

	int face_count;
	Obj_Face *face_array;

	int material_count;
	Obj_Material *material_array;
};

enum Load_Flags {
	Load_Flag_Weld_Attributes,
};

struct Obj_Load_Options {
	int load_flags;
};

Obj_Mesh *obj_mesh_load(const char *filename, const Obj_Load_Options *options)
{
	using namespace std;
	vector<tinyobj::shape_t> shapes;
	vector<tinyobj::material_t> materials;
	string err;
	bool ret = tinyobj::LoadObj(shapes, materials, err, filename);
	if (!ret) {
		printf("%s\n", err.c_str());
		return NULL;
	}
	printf("%lu shapes\n", shapes.size());
	printf("%lu materials\n", materials.size());
	assert(shapes.size() > 0);
	Obj_Mesh *mesh = new Obj_Mesh();
	mesh->vertex_count = shapes[0].mesh.positions.size() / 3;
	mesh->vertex_array = new Obj_Vertex[mesh->vertex_count];
	for (int nvert = 0; nvert < mesh->vertex_count; nvert++) {
		mesh->vertex_array[nvert].position[0] = shapes[0].mesh.positions[nvert * 3];
		mesh->vertex_array[nvert].position[1] = shapes[0].mesh.positions[nvert * 3 + 1];
		mesh->vertex_array[nvert].position[2] = shapes[0].mesh.positions[nvert * 3 + 2];
		mesh->vertex_array[nvert].normal[0] = shapes[0].mesh.normals[nvert * 3];
		mesh->vertex_array[nvert].normal[1] = shapes[0].mesh.normals[nvert * 3 + 1];
		mesh->vertex_array[nvert].normal[2] = shapes[0].mesh.normals[nvert * 3 + 2];
		mesh->vertex_array[nvert].uv[0] = 0;
		mesh->vertex_array[nvert].uv[1] = 0;
		mesh->vertex_array[nvert].first_colocal = nvert;
	}
	mesh->face_count = shapes[0].mesh.indices.size() / 3;
	mesh->face_array = new Obj_Face[mesh->face_count];
	for (int nface = 0; nface < mesh->face_count; nface++) {
		mesh->face_array[nface].material_index = 0;
		mesh->face_array[nface].vertex_index[0] = shapes[0].mesh.indices[nface * 3];
		mesh->face_array[nface].vertex_index[1] = shapes[0].mesh.indices[nface * 3 + 1];
		mesh->face_array[nface].vertex_index[2] = shapes[0].mesh.indices[nface * 3 + 2];
	}
	printf("Reading %d verts\n", mesh->vertex_count);
	printf("Reading %d triangles\n", mesh->face_count);
	mesh->material_count = 0;
	mesh->material_array = 0;
	return mesh;
}

void obj_mesh_free(Obj_Mesh *mesh)
{
	if (mesh != NULL) {
		delete [] mesh->vertex_array;
		delete [] mesh->face_array;
		delete [] mesh->material_array;
		delete mesh;
	}
}

struct RasterParam
{
	xatlas::Atlas *atlas;
	uint8_t color[3];
	uint8_t *output_image;
};

bool RasterCallback(void *param, int x, int y, xatlas::internal::Vector3::Arg bar, xatlas::internal::Vector3::Arg dx, xatlas::internal::Vector3::Arg dy, float coverage)
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
	/*if (argc != 2) {
	    printf("Usage: %s input_file.obj\n", argv[0]);
		system("pause");
	    return EXIT_FAILURE;
	}*/
	// Load Obj_Mesh.
	Obj_Load_Options load_options = {0};
	Obj_Mesh *obj_mesh = obj_mesh_load("C:/Projects/atlas/gazebo.obj", &load_options);
	if (obj_mesh == NULL) {
		printf("Error loading obj file.\n");
		return 1;
	}
	// Convert Obj_Mesh to Atlast_Input_Mesh.
	assert(sizeof(xatlas::Input_Vertex) == sizeof(Obj_Vertex));
	assert(sizeof(xatlas::Input_Face) == sizeof(Obj_Face));
	xatlas::Input_Mesh input_mesh;
	input_mesh.vertex_count = obj_mesh->vertex_count;
	input_mesh.vertex_array = (xatlas::Input_Vertex *)obj_mesh->vertex_array;
	input_mesh.face_count = obj_mesh->face_count;
	input_mesh.face_array = (xatlas::Input_Face *)obj_mesh->face_array;
	// Generate Output_Mesh.
	xatlas::Options atlas_options;
	xatlas::set_default_options(&atlas_options);
	// Avoid brute force packing, since it can be unusably slow in some situations.
	atlas_options.packer.packing_quality = 1;
	atlas_options.packer.texel_area = 256;
	atlas_options.packer.conservative = true;
	atlas_options.packer.padding = 2;
	//xatlas::add_mesh(&input_mesh);
	xatlas::add_mesh(&input_mesh);
	xatlas::Atlas atlas = xatlas::atlas_generate(&atlas_options);
	printf("%d output meshes\n", atlas.nMeshes);
	uint8_t *output_image = new uint8_t[atlas.width * atlas.height * 4];
	memset(output_image, 0, atlas.width * atlas.height * 4);
	for (int i = 0; i < atlas.nMeshes; i++) {
		xatlas::Output_Mesh *output_mesh = atlas.meshes[i];
		printf("Output mesh %d has %d verts, %d triangles\n", i, output_mesh->vertex_count, output_mesh->index_count / 3);
		for (int j = 0; j < output_mesh->index_count; j += 3) {
			const xatlas::Output_Vertex *v[3];
			v[0] = &output_mesh->vertex_array[output_mesh->index_array[j + 0]];
			v[1] = &output_mesh->vertex_array[output_mesh->index_array[j + 1]];
			v[2] = &output_mesh->vertex_array[output_mesh->index_array[j + 2]];
			xatlas::internal::Vector2 verts[3];
			for (int k = 0; k < 3; k++) {
				//normalizedVerts[j] = nv::Vector2(v[j].uv[0] * atlas.width, v[j].uv[1] * atlas.height);
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
	stbi_write_tga("debug_packer_final.tga", atlas.width, atlas.height, 4, output_image);
	delete [] output_image;
	printf("Produced debug_packer_final.tga\n");
	// Free meshes.
	obj_mesh_free(obj_mesh);
	xatlas::atlas_free(atlas);
	return 0;
}
