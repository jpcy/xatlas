#include "nvmesh/raster/Raster.h"
#include "xatlas.h"

#include "thekla_atlas.h"
#include "thekla_mesh_load.h"

#include <stdio.h>
#include <assert.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../extern/stb_image_write.h"

using namespace Thekla;

Atlas_Output_Mesh * output_mesh;
static uint8_t *output_image;

bool RasterCallback(void * param, int x, int y, nv::Vector3::Arg bar, nv::Vector3::Arg dx, nv::Vector3::Arg dy, float coverage)
{
	uint8_t *color = (uint8_t *)param;
	uint8_t *rgba = &output_image[(x + y * output_mesh->atlas_width) * 4];
	rgba[0] = color[0];
	rgba[1] = color[1];
	rgba[2] = color[2];
	rgba[3] = 255;
	return true;
}

int main(int argc, char * argv[]) {

    /*if (argc != 2) {
        printf("Usage: %s input_file.obj\n", argv[0]);
		system("pause");
        return EXIT_FAILURE;
    }*/

    // Load Obj_Mesh.
    Obj_Load_Options load_options = {0};
    Obj_Mesh * obj_mesh = obj_mesh_load("C:/Projects/atlas/gazebo.obj", &load_options);

    if (obj_mesh == NULL) {
        printf("Error loading obj file.\n");
        return 1;
    }

    // Convert Obj_Mesh to Atlast_Input_Mesh.
    assert(sizeof(Atlas_Input_Vertex) == sizeof(Obj_Vertex));
    assert(sizeof(Atlas_Input_Face) == sizeof(Obj_Face));

    Atlas_Input_Mesh input_mesh;
    input_mesh.vertex_count = obj_mesh->vertex_count;
    input_mesh.vertex_array = (Atlas_Input_Vertex *)obj_mesh->vertex_array;
    input_mesh.face_count = obj_mesh->face_count;
    input_mesh.face_array = (Atlas_Input_Face *)obj_mesh->face_array;


    // Generate Atlas_Output_Mesh.
    Atlas_Options atlas_options;
    atlas_set_default_options(&atlas_options);

    // Avoid brute force packing, since it can be unusably slow in some situations.
    atlas_options.packer_options.witness.packing_quality = 1;
	atlas_options.packer_options.witness.texel_area = 256;

    Atlas_Error error = Atlas_Error_Success;
    output_mesh = atlas_generate(&input_mesh, &atlas_options, &error);

    printf("Atlas mesh has %d verts\n", output_mesh->vertex_count);
    printf("Atlas mesh has %d triangles\n", output_mesh->index_count / 3);

	output_image = new uint8_t[output_mesh->atlas_width * output_mesh->atlas_height * 4];
	memset(output_image, 0, output_mesh->atlas_width * output_mesh->atlas_height * 4);

	for (int i = 0; i < output_mesh->index_count; i += 3)
	{
		const Atlas_Output_Vertex *v[3];
		v[0] = &output_mesh->vertex_array[output_mesh->index_array[i + 0]];
		v[1] = &output_mesh->vertex_array[output_mesh->index_array[i + 1]];
		v[2] = &output_mesh->vertex_array[output_mesh->index_array[i + 2]];
		nv::Vector2 verts[3];

		for (int j = 0; j < 3; j++)
		{
			//normalizedVerts[j] = nv::Vector2(v[j].uv[0] * output_mesh->atlas_width, v[j].uv[1] * output_mesh->atlas_height);
			verts[j] = nv::Vector2(v[j]->uv[0], v[j]->uv[1]);
		}

		uint8_t color[3];
		color[0] = rand() % 255;
		color[1] = rand() % 255;
		color[2] = rand() % 255;
		nv::Raster::drawTriangle(nv::Raster::Mode_Nearest, nv::Vector2((float)output_mesh->atlas_width, (float)output_mesh->atlas_height), true, verts, RasterCallback, color);
	}

	stbi_write_tga("debug_packer_final.tga", output_mesh->atlas_width, output_mesh->atlas_height, 4, output_image);
    printf("Produced debug_packer_final.tga\n");

    // Free meshes.
    obj_mesh_free(obj_mesh);
    atlas_free(output_mesh);
 
    return 0;

}
