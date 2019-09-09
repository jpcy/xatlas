/*
xatlas
https://github.com/jpcy/xatlas
Copyright (c) 2018 Jonathan Young

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#include <atomic>
#include <cstddef>
#include <cmath>
#include <thread>
#include <unordered_map>
#include <bx/filepath.h>
#include <bx/string.h>
#include <cgltf.h>
#include <GLFW/glfw3.h>
#include <imgui/imgui.h>
#include <nativefiledialog/nfd.h>
#include <ofbx.h>
#include <stb_image.h>
#include <stb_image_resize.h>
#define STL_READER_NO_EXCEPTIONS
#include <stl_reader.h>
#include "shaders/shared.h"
#include "viewer.h"

#define ONE_GLTF_OBJECT_PER_MESH 0

bgfx::VertexLayout ModelVertex::layout;

enum class ModelStatus
{
	NotLoaded,
	Loading,
	Finalizing,
	Loaded
};

enum class ModelFormat
{
	Gltf,
	Obj,
	Stl
};

struct
{
	std::atomic<ModelStatus> status;
	std::thread *thread = nullptr;
	objzModel *data = nullptr;
	void (*destroyModelData)(objzModel *) = nullptr;
	std::vector<uint32_t> diffuseTextures;
	std::vector<uint32_t> emissionTextures;
	AABB aabb;
	bx::Vec3 centroid = bx::Vec3(0.0f, 0.0f, 0.0f);
	bgfx::VertexBufferHandle vb = BGFX_INVALID_HANDLE;
	bgfx::IndexBufferHandle ib = BGFX_INVALID_HANDLE;
	bgfx::VertexBufferHandle wireframeVb = BGFX_INVALID_HANDLE;
	std::vector<WireframeVertex> wireframeVertices;
	float scale = 1.0f;
	bool rightHandedAxis = false; // Default is z/-x/y, right handed is -z/x/y.
	bool clockwiseFaceWinding = true;
	bgfx::ShaderHandle vs_model;
	bgfx::ShaderHandle fs_material;
	bgfx::ProgramHandle materialProgram;
	bgfx::UniformHandle u_diffuse;
	bgfx::UniformHandle u_emission;
	bgfx::UniformHandle u_lightDir;
	bgfx::UniformHandle u_shade_diffuse_emission;
	bgfx::UniformHandle s_diffuse;
	bgfx::UniformHandle s_emission;
	bgfx::UniformHandle s_lightmap;
	bgfx::UniformHandle u_color;
	bgfx::TextureHandle u_dummyTexture;
}
s_model;

static const float s_rightHandedAxisMatrix[16] = {
	0.0f, 0.0f, -1.0f, 0.0f,
	1.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 1.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 1.0f
};

static bool readFileData(const char *filename, std::vector<uint8_t> *fileData)
{
#if _MSC_VER
	FILE *f;
	if (fopen_s(&f, filename, "rb") != 0)
		f = nullptr;
#else
	FILE *f = fopen(filename, "rb");
#endif
	if (!f) {
		fprintf(stderr, "Error opening '%s'\n", filename);
		return false;
	}
	fseek(f, 0, SEEK_END);
	const long length = ftell(f);
	fseek(f, 0, SEEK_SET);
	fileData->resize(length);
	if (fread(fileData->data(), 1, (size_t)length, f) < (size_t)length) {
		fclose(f);
		fprintf(stderr, "Error reading '%s'\n", filename);
		return false;
	}
	fclose(f);
	return true;
}

struct TextureData
{
	uint16_t width;
	uint16_t height;
	const bgfx::Memory *mem;
	int numComponents;
	// Data used by baking to sample the texture.
	uint8_t *sampleData;
	uint32_t sampleWidth, sampleHeight;
};

static TextureData textureLoad(const char *basePath, const char *filename)
{
	char fullFilename[256] = { 0 };
	bx::strCopy(fullFilename, sizeof(fullFilename), basePath);
	bx::strCat(fullFilename, sizeof(fullFilename), filename);
	TextureData td;
	td.mem = nullptr;
	td.sampleData = nullptr;
	std::vector<uint8_t> fileData;
	if (!readFileData(fullFilename, &fileData))
		return td;
	int width, height, numComponents;
	const uint8_t *imageData = stbi_load_from_memory(fileData.data(), (int)fileData.size(), &width, &height, &numComponents, 0);
	if (!imageData) {
		fprintf(stderr, "Error loading '%s': %s\n", fullFilename, stbi_failure_reason());
		return td;
	}
	printf("Texture '%s': %dx%d %d bpp\n", fullFilename, width, height, numComponents * 8);
	// Generate mipmaps.
	const int nMips = 1 + (int)bx::floor(bx::log2((float)bx::max(width, height)));
	int mipWidth = width, mipHeight = height;
	uint32_t memSize = 0;
	for (int i = 0; i < nMips; i++) {
		memSize += uint32_t(mipWidth * mipHeight * numComponents);
		mipWidth = bx::max(mipWidth >> 1, 1);
		mipHeight = bx::max(mipHeight >> 1, 1);
	}
	const bgfx::Memory *mem = bgfx::alloc(memSize);
	memcpy(mem->data, imageData, width * height * numComponents);
	stbi_image_free((void *)imageData);
	const uint8_t *src = mem->data;
	int srcWidth = width, srcHeight = height;
	uint8_t *dest = mem->data;
	mipWidth = width;
	mipHeight = height;
	for (int i = 0; i < nMips - 1; i++) {
		dest += mipWidth * mipHeight * numComponents;
		mipWidth = bx::max(mipWidth >> 1, 1);
		mipHeight = bx::max(mipHeight >> 1, 1);
		stbir_resize_uint8_srgb(src, srcWidth, srcHeight, srcWidth * numComponents, dest, mipWidth, mipHeight, mipWidth * numComponents, numComponents, numComponents == 4 ? 3 : STBIR_ALPHA_CHANNEL_NONE, 0);
		src = dest;
		srcWidth = mipWidth;
		srcHeight = mipHeight;
		// Copy a small mip for baking to use for sampling textures.
		if (!td.sampleData && (mipWidth <= 32 || mipHeight <= 32)) {
			const size_t size = mipWidth * mipHeight * numComponents;
			td.sampleData = new uint8_t[size];
			memcpy(td.sampleData, dest, size);
			td.sampleWidth = (uint32_t)mipWidth;
			td.sampleHeight = (uint32_t)mipHeight;
		}
	}
	td.mem = mem;
	td.width = (uint16_t)width;
	td.height = (uint16_t)height;
	td.numComponents = numComponents;
	return td;
}

struct CachedTexture
{
	char filename[256];
	TextureData data;
	bgfx::TextureHandle handle;
};

static std::vector<CachedTexture> s_textureCache;

static uint32_t textureLoadCached(const char *basePath, const char *filename)
{
	for (uint32_t i = 0; i < (uint32_t)s_textureCache.size(); i++) {
		if (bx::strCmpI(s_textureCache[i].filename, filename) == 0)
			return i;
	}
	CachedTexture texture;
	bx::strCopy(texture.filename, sizeof(texture.filename), filename);
	texture.data = textureLoad(basePath, filename);
	texture.handle = BGFX_INVALID_HANDLE;
	s_textureCache.push_back(texture);
	return (uint32_t)s_textureCache.size() - 1;
}

static void textureCreateCachedTextures()
{
	for (uint32_t i = 0; i < (uint32_t)s_textureCache.size(); i++) {
		CachedTexture &texture = s_textureCache[i];
		if (!texture.data.mem) {
			texture.handle = BGFX_INVALID_HANDLE;
			continue;
		}
		bgfx::TextureFormat::Enum format = bgfx::TextureFormat::RGBA8;
		if (texture.data.numComponents == 1)
			format = bgfx::TextureFormat::R8;
		else if (texture.data.numComponents == 3)
			format = bgfx::TextureFormat::RGB8;
		texture.handle = bgfx::createTexture2D(texture.data.width, texture.data.height, true, 1, format, BGFX_SAMPLER_MIN_ANISOTROPIC | BGFX_SAMPLER_MAG_ANISOTROPIC, texture.data.mem);
	}
}

static bgfx::TextureHandle textureGetHandle(uint32_t index)
{
	if (index == UINT32_MAX)
		return BGFX_INVALID_HANDLE;
	return s_textureCache[index].handle;
}

static void textureDestroyCache()
{
	for (int i = 0; i < (int)s_textureCache.size(); i++) {
		bgfx::destroy(s_textureCache[i].handle);
		delete s_textureCache[i].data.sampleData;
	}
	s_textureCache.clear();
}

void modelInit()
{
	s_model.status = ModelStatus::NotLoaded;
	s_model.u_color = bgfx::createUniform("u_color", bgfx::UniformType::Vec4);
	s_model.u_diffuse = bgfx::createUniform("u_diffuse", bgfx::UniformType::Vec4);
	s_model.u_emission = bgfx::createUniform("u_emission", bgfx::UniformType::Vec4);
	s_model.u_lightDir = bgfx::createUniform("u_lightDir", bgfx::UniformType::Vec4);
	s_model.u_shade_diffuse_emission = bgfx::createUniform("u_shade_diffuse_emission", bgfx::UniformType::Vec4);
	s_model.s_diffuse = bgfx::createUniform("s_diffuse", bgfx::UniformType::Sampler);
	s_model.s_emission = bgfx::createUniform("s_emission", bgfx::UniformType::Sampler);
	s_model.s_lightmap = bgfx::createUniform("s_lightmap", bgfx::UniformType::Sampler);
	s_model.vs_model = loadShader(ShaderId::vs_model);
	s_model.fs_material = loadShader(ShaderId::fs_material);
	s_model.materialProgram = bgfx::createProgram(s_model.vs_model, s_model.fs_material);
	s_model.u_dummyTexture = bgfx::createTexture2D(16, 16, false, 1, bgfx::TextureFormat::BGRA8);
	ModelVertex::init();
	bgfx::setViewClear(kModelView, BGFX_CLEAR_COLOR | BGFX_CLEAR_DEPTH, 0x444444ff);
	bgfx::setViewRect(kModelView, 0, 0, bgfx::BackbufferRatio::Equal);
	bgfx::setViewRect(kModelTransparentView, 0, 0, bgfx::BackbufferRatio::Equal);
	bgfx::setViewMode(kModelTransparentView, bgfx::ViewMode::DepthDescending);
}

void modelShutdown()
{
	modelDestroy();
	bgfx::destroy(s_model.u_color);
	bgfx::destroy(s_model.u_diffuse);
	bgfx::destroy(s_model.u_emission);
	bgfx::destroy(s_model.u_lightDir);
	bgfx::destroy(s_model.u_shade_diffuse_emission);
	bgfx::destroy(s_model.s_diffuse);
	bgfx::destroy(s_model.s_emission);
	bgfx::destroy(s_model.s_lightmap);
	bgfx::destroy(s_model.vs_model);
	bgfx::destroy(s_model.fs_material);
	bgfx::destroy(s_model.materialProgram);
	bgfx::destroy(s_model.u_dummyTexture);
}

static void fbxDestroy(objzModel *model)
{
	delete [] (uint32_t *)model->indices;
	delete [] model->meshes;
	delete [] model->objects;
	delete [] (ModelVertex *)model->vertices;
	if (model->materials)
		delete [] model->materials;
	delete model;
}

static objzModel *fbxLoad(const char *filename, const char * /*basePath*/) 
{
	std::vector<uint8_t> fileData;
	if (!readFileData(filename, &fileData))
		return nullptr;
	ofbx::IScene *scene = ofbx::load(fileData.data(), (int)fileData.size(), ofbx::LoadFlags::TRIANGULATE);
	if (!scene) {
		fprintf(stderr, "%s\n", ofbx::getError());
		return nullptr;
	}
	objzModel *model = new objzModel();
	model->numIndices = 0;
	model->numMaterials = 0;
	model->numMeshes = (uint32_t)scene->getMeshCount();
	model->numObjects = 1;
	model->numVertices = 0;
	// Count array lengths.
	std::unordered_map<const ofbx::Material *, uint32_t> materialToIndex;
	for (int i = 0; i < scene->getAllObjectCount(); i++) {
		const ofbx::Object *object = scene->getAllObjects()[i];
		if (object->getType() == ofbx::Object::Type::MATERIAL) {
			materialToIndex[(const ofbx::Material *)object] = model->numMaterials;
			model->numMaterials++;
		}
	}
	for (int i = 0; i < scene->getMeshCount(); i++) {
		const ofbx::Geometry *geo = scene->getMesh(i)->getGeometry();
		model->numIndices += (uint32_t)geo->getIndexCount();
		model->numVertices += (uint32_t)geo->getVertexCount();
	}
	// Alloc data.
	auto indices = new uint32_t[model->numIndices];
	auto vertices = new ModelVertex[model->numVertices];
	model->indices = indices;
	model->meshes = new objzMesh[model->numMeshes];
	model->objects = new objzObject[model->numObjects];
	model->vertices = vertices;
	if (model->numMaterials > 0)
		model->materials = new objzMaterial[model->numMaterials];
	else
		model->materials = nullptr;
	// Populate data.
	{
		objzObject &object = model->objects[0];
		object.name[0] = 0;
		object.firstMesh = 0;
		object.numMeshes = model->numMeshes;
		object.firstIndex = 0;
		object.numIndices = model->numIndices;
		object.firstVertex = 0;
		object.numVertices = model->numVertices;
	}
	uint32_t currentIndex = 0, currentVertex = 0;
	for (int i = 0; i < scene->getMeshCount(); i++) {
		const ofbx::Mesh *sourceMesh = scene->getMesh(i);
		ofbx::Matrix dtransform = sourceMesh->getGlobalTransform();
		float transform[16];
		for (uint32_t j = 0; j < 16; j++)
			transform[j] = (float)dtransform.m[j];
		const ofbx::Geometry *sourceGeo = scene->getMesh(i)->getGeometry();
		objzMesh &mesh = model->meshes[i];
		mesh.firstIndex = currentIndex;
		mesh.numIndices = (uint32_t)sourceGeo->getIndexCount();
		// ignoring all but the first material for now
		if (sourceMesh->getMaterialCount() > 0)
			mesh.materialIndex = materialToIndex[sourceMesh->getMaterial(0)];
		else
			mesh.materialIndex = -1;
		for (uint32_t j = 0; j < mesh.numIndices; j++) {
			int sourceIndex = sourceGeo->getFaceIndices()[j];
			if (sourceIndex < 0)
				sourceIndex = -sourceIndex - 1; // index is negative if last in face
			if (sourceIndex >= sourceGeo->getVertexCount()) {
				fprintf(stderr, "Index '%d' out of range of vertex count '%d'\n", sourceIndex, sourceGeo->getVertexCount());
				scene->destroy();
				fbxDestroy(model);
				return nullptr;
			}
			const uint32_t index = currentVertex + (uint32_t)sourceIndex;
			assert(index < model->numVertices);
			indices[mesh.firstIndex + j] = index;
		}
		for (uint32_t j = 0; j < (uint32_t)sourceGeo->getVertexCount(); j++) {
			ModelVertex &vertex = vertices[currentVertex + j];
			const ofbx::Vec3 &pos = sourceGeo->getVertices()[j];
			vertex.pos.x = (float)pos.x;
			vertex.pos.y = (float)pos.y;
			vertex.pos.z = (float)pos.z;
			vertex.pos = bx::mul(vertex.pos, transform);
			if (sourceGeo->getNormals()) {
				const ofbx::Vec3 &normal = sourceGeo->getNormals()[j];
				vertex.normal.x = (float)normal.x;
				vertex.normal.y = (float)normal.y;
				vertex.normal.z = (float)normal.z;
			} else {
				vertex.normal = bx::Vec3(0.0f);
			}
			if (sourceGeo->getUVs(0)) {
				const ofbx::Vec2 &uv = sourceGeo->getUVs(0)[j];
				vertex.texcoord[0] = (float)uv.x;
				vertex.texcoord[1] = (float)uv.y;
			} else {
				vertex.texcoord[0] = vertex.texcoord[1] = 0.0f;
			}
			vertex.texcoord[2] = vertex.texcoord[3] = 0.0f;
		}
		currentIndex += mesh.numIndices;
		currentVertex += (uint32_t)sourceGeo->getVertexCount();
	}
	uint32_t currentMaterial = 0;
	for (int i = 0; i < scene->getAllObjectCount(); i++) {
		const ofbx::Object *object = scene->getAllObjects()[i];
		if (object->getType() != ofbx::Object::Type::MATERIAL)
			continue;
		auto sourceMat = (const ofbx::Material *)object;
		objzMaterial &destMat = model->materials[currentMaterial];
		memset(&destMat, 0, sizeof(destMat));
		destMat.opacity = 1.0f;
		const ofbx::Color &diffuse = sourceMat->getDiffuseColor();
		destMat.diffuse[0] = diffuse.r;
		destMat.diffuse[1] = diffuse.g;
		destMat.diffuse[2] = diffuse.b;
		currentMaterial++;
	}
	scene->destroy();
	return model;
}

BX_PRAGMA_DIAGNOSTIC_PUSH();
BX_PRAGMA_DIAGNOSTIC_IGNORED_MSVC(4702) // 'unreachable code'
bool gltfAnyNodeInHierarchyHasMesh(const cgltf_node *node)
{
	if (node->mesh)
		return true;
	for (cgltf_size ci = 0; ci < node->children_count; ci++)
		return gltfAnyNodeInHierarchyHasMesh(node->children[ci]);
	return false;
}
BX_PRAGMA_DIAGNOSTIC_POP();

static void gltfCountMeshData(const cgltf_node *node, objzModel *model)
{
	if (node->mesh) {
		for (cgltf_size pi = 0; pi < node->mesh->primitives_count; pi++) {
			const cgltf_primitive &primitive = node->mesh->primitives[pi];
			const cgltf_accessor *apositions = nullptr;
			for (cgltf_size ai = 0; ai < primitive.attributes_count; ai++) {
				const cgltf_attribute &attrib = primitive.attributes[ai];
				if (attrib.type == cgltf_attribute_type_position) {
					apositions = attrib.data;
					break;
				}
			}
			const cgltf_accessor *aindices = primitive.indices;
			if (apositions && aindices) {
				model->numVertices += (uint32_t)apositions->count;
				model->numIndices += (uint32_t)aindices->count;
			}
			model->numMeshes++;
#if ONE_GLTF_OBJECT_PER_MESH
			model->numObjects++;
#endif
		}
	}
	for (cgltf_size ci = 0; ci < node->children_count; ci++)
		gltfCountMeshData(node->children[ci], model);
}

template<typename T>
static const T *gltfGetBufferData(const cgltf_accessor *accessor)
{
	auto buffer = (const uint8_t *)accessor->buffer_view->buffer->data;
	const cgltf_size offset = accessor->offset + accessor->buffer_view->offset;
	return (const T *)&buffer[offset];
}

static void gltfPopulateMeshData(const cgltf_node *node, const cgltf_material *firstMaterial, objzModel *model, uint32_t &currentObject, uint32_t &currentMesh, uint32_t &firstMeshIndex, uint32_t &firstMeshVertex)
{
	const cgltf_mesh *sourceMesh = node->mesh;
	if (sourceMesh) {
		float transform[16];
		cgltf_node_transform_world(node, transform);
		float rotation[16];
		if (node->has_rotation)
			bx::mtxQuat(rotation, *(bx::Quaternion *)node->rotation);
		else
			bx::mtxIdentity(rotation);
		for (cgltf_size pi = 0; pi < sourceMesh->primitives_count; pi++) {
			const cgltf_primitive &primitive = sourceMesh->primitives[pi];
			const cgltf_accessor *apositions = nullptr, *anormals = nullptr, *atexcoords = nullptr;
			for (cgltf_size ai = 0; ai < primitive.attributes_count; ai++) {
				const cgltf_attribute &attrib = primitive.attributes[ai];
				if (attrib.type == cgltf_attribute_type_position)
					apositions = attrib.data;
				else if (attrib.type == cgltf_attribute_type_normal)
					anormals = attrib.data;
				else if (attrib.type == cgltf_attribute_type_texcoord)
					atexcoords = attrib.data;
			}
			const cgltf_accessor *aindices = primitive.indices;
			if (!apositions || !aindices)
				continue;
			// Copy vertex data.
			const float *meshPosition = gltfGetBufferData<float>(apositions);
			const float *meshNormal = anormals && anormals->count == apositions->count ? gltfGetBufferData<float>(anormals) : nullptr;
			const float *meshTexcoord = atexcoords && atexcoords->count == apositions->count ? gltfGetBufferData<float>(atexcoords) : nullptr;
			for (cgltf_size vi = 0; vi < apositions->count; vi++) {
				assert(vi + firstMeshVertex < model->numVertices);
				ModelVertex &vertex = ((ModelVertex *)model->vertices)[vi + firstMeshVertex];
				vertex.pos = bx::mul(bx::Vec3(meshPosition[0], meshPosition[1], meshPosition[2]), transform);
				meshPosition += apositions->stride / sizeof(float);
				if (meshNormal) {
					vertex.normal = bx::Vec3(meshNormal[0], meshNormal[1], meshNormal[2]);
					if (node->has_rotation)
						vertex.normal = bx::mul(vertex.normal, rotation);
					meshNormal += anormals->stride / sizeof(float);
				}
				if (meshTexcoord) {
					vertex.texcoord[0] = meshTexcoord[0];
					vertex.texcoord[1] = meshTexcoord[1];
					meshTexcoord += atexcoords->stride / sizeof(float);
				}
			}
			// Copy indices.
			auto indices = (uint32_t *)model->indices;
			if (aindices->component_type == cgltf_component_type_r_16u) {
				const uint16_t *meshIndices = gltfGetBufferData<uint16_t>(aindices);
				for (uint32_t ii = 0; ii < (uint32_t)aindices->count; ii++) {
					assert(ii + firstMeshIndex < model->numIndices);
					indices[ii + firstMeshIndex] = firstMeshVertex + (uint32_t)meshIndices[ii];
				}
			} else  if (aindices->component_type == cgltf_component_type_r_32u) {
				const uint32_t *meshIndices = gltfGetBufferData<uint32_t>(aindices);
				for (uint32_t ii = 0; ii < (uint32_t)aindices->count; ii++) {
					assert(ii + firstMeshIndex < model->numIndices);
					indices[ii + firstMeshIndex] = firstMeshVertex + meshIndices[ii];
				}
			}
#if ONE_GLTF_OBJECT_PER_MESH
			// Create object.
			assert(currentObject < model->numObjects);
			objzObject &object = model->objects[currentObject];
			object.name[0] = 0;
			object.firstMesh = currentMesh;
			object.numMeshes = 1;
			object.firstIndex = firstMeshIndex;
			object.numIndices = (uint32_t)aindices->count;
			object.firstVertex = firstMeshVertex;
			object.numVertices = (uint32_t)apositions->count;
			currentObject++;
#endif
			// Create mesh.
			assert(currentMesh < model->numMeshes);
			objzMesh &mesh = model->meshes[currentMesh];
			mesh.materialIndex = primitive.material ? int32_t(primitive.material - firstMaterial) : -1;
			mesh.firstIndex = firstMeshIndex;
			mesh.numIndices = (uint32_t)aindices->count;
			currentMesh++;
#if !ONE_GLTF_OBJECT_PER_MESH
			// Update object.
			objzObject &object = model->objects[currentObject];
			object.numMeshes++;
			object.numIndices += (uint32_t)aindices->count;
			object.numVertices += (uint32_t)apositions->count;
#endif
			firstMeshVertex += (uint32_t)apositions->count;
			firstMeshIndex += (uint32_t)aindices->count;
		}
	}
	for (cgltf_size ci = 0; ci < node->children_count; ci++)
		gltfPopulateMeshData(node->children[ci], firstMaterial, model, currentObject, currentMesh, firstMeshIndex, firstMeshVertex);
}

static objzModel *gltfLoad(const char *filename, const char *basePath) 
{
	cgltf_data *gltfData = nullptr;
	cgltf_options options;
	bx::memSet(&options, 0, sizeof(options));
	cgltf_result result = cgltf_parse_file(&options, filename, &gltfData);
	if (result == cgltf_result_success) {
		result = cgltf_load_buffers(&options, gltfData, basePath);
		if (result == cgltf_result_success)
			result = cgltf_validate(gltfData);
	}
	if (result != cgltf_result_success) {
		if (gltfData)
			cgltf_free(gltfData);
		return nullptr;
	}
	objzModel *model = new objzModel();
	model->numIndices = 0;
	model->numMaterials = (uint32_t)gltfData->materials_count;
	model->numMeshes = 0;
	model->numObjects = 0;
	model->numVertices = 0;
	// Count array lengths.
	for (cgltf_size ni = 0; ni < gltfData->nodes_count; ni++) {
		// Objects are root nodes with a mesh, or any ancestor with a mesh.
		const cgltf_node &node = gltfData->nodes[ni];
		if (node.parent)
			continue;
		if (!gltfAnyNodeInHierarchyHasMesh(&node))
			continue;
		gltfCountMeshData(&node, model);
#if !ONE_GLTF_OBJECT_PER_MESH
		model->numObjects++;
#endif
	}
	// Alloc data.
	model->indices = new uint32_t[model->numIndices];
	model->meshes = new objzMesh[model->numMeshes];
	model->objects = new objzObject[model->numObjects];
	model->vertices = new ModelVertex[model->numVertices];
	// Populate data.
	uint32_t currentObject = 0, currentMesh = 0, firstMeshIndex = 0, firstMeshVertex = 0;
	for (cgltf_size ni = 0; ni < gltfData->nodes_count; ni++) {
		const cgltf_node &node = gltfData->nodes[ni];
		if (node.parent)
			continue;
		if (!gltfAnyNodeInHierarchyHasMesh(&node))
			continue;
#if !ONE_GLTF_OBJECT_PER_MESH
		// Create object.
		assert(currentObject < model->numObjects);
		objzObject &object = model->objects[currentObject];
		bx::strCopy(object.name, sizeof(object.name), node.name);
		object.firstMesh = currentMesh;
		object.numMeshes = 0;
		object.firstIndex = firstMeshIndex;
		object.numIndices = 0;
		object.firstVertex = firstMeshVertex;
		object.numVertices = 0;
#endif
		// Create mesh data.
		gltfPopulateMeshData(&node, gltfData->materials, model, currentObject, currentMesh, firstMeshIndex, firstMeshVertex);
#if !ONE_GLTF_OBJECT_PER_MESH
		currentObject++;
#endif
	}
	// Materials.
	model->materials = new objzMaterial[model->numMaterials];
	for (uint32_t i = 0; i < model->numMaterials; i++) {
		const cgltf_material &sourceMat = gltfData->materials[i];
		const cgltf_texture *diffuse = sourceMat.pbr_metallic_roughness.base_color_texture.texture;
		const cgltf_texture *emission = sourceMat.emissive_texture.texture;
		objzMaterial &destMat = model->materials[i];
		memset(&destMat, 0, sizeof(destMat));
		destMat.opacity = 1.0f;
		memcpy(destMat.diffuse, sourceMat.pbr_metallic_roughness.base_color_factor, sizeof(float) * 3);
		if (diffuse)
			bx::strCopy(destMat.diffuseTexture, sizeof(destMat.diffuseTexture), diffuse->image->uri);
		if (emission)
			bx::strCopy(destMat.emissionTexture, sizeof(destMat.emissionTexture), emission->image->uri);
	}
	cgltf_free(gltfData);
	return model;
}

static void gltfDestroy(objzModel *model)
{
	delete [] (uint32_t *)model->indices;
	delete [] model->meshes;
	delete [] model->objects;
	delete [] (ModelVertex *)model->vertices;
	delete [] model->materials;
	delete model;
}

static objzModel *stlLoad(const char *filename, const char * /*basePath*/) 
{
	std::vector<float> coords, normals;
	std::vector<unsigned int> tris, solids;
	if (!stl_reader::ReadStlFile(filename, coords, normals, tris, solids))
		return nullptr;
	objzModel *model = new objzModel();
	model->numIndices = (uint32_t)tris.size();
	model->numMaterials = 0;
	model->numMeshes = (uint32_t)solids.size() - 1;
	model->numObjects = (uint32_t)solids.size() - 1;
	model->numVertices = (uint32_t)coords.size() / 3;
	model->indices = new uint32_t[model->numIndices];
	model->materials = nullptr;
	model->meshes = new objzMesh[model->numMeshes];
	model->objects = new objzObject[model->numObjects];
	model->vertices = new ModelVertex[model->numVertices];
	for (uint32_t i = 0; i < model->numObjects; i++) {
		objzObject &object = model->objects[i];
		object.name[0] = 0;
		object.firstMesh = i;
		object.numMeshes = 1;
		object.firstIndex = solids[i] * 3;
		object.numIndices = solids[i + 1] * 3;
		object.firstVertex = 0;
		object.numVertices = model->numVertices;
		objzMesh &mesh = model->meshes[i];
		mesh.materialIndex = -1;
		mesh.firstIndex = object.firstIndex;
		mesh.numIndices = object.numIndices;
	}
	auto vertices = (ModelVertex *)model->vertices;
	for (uint32_t i = 0; i < model->numVertices; i++) {
		ModelVertex &v = vertices[i];
		v.pos.x = coords[i * 3 + 0];
		v.pos.y = coords[i * 3 + 1];
		v.pos.z = coords[i * 3 + 2];
		v.normal.x = normals[i * 3 + 0];
		v.normal.y = normals[i * 3 + 1];
		v.normal.z = normals[i * 3 + 2];
		v.texcoord[0] = v.texcoord[1] = v.texcoord[2] = v.texcoord[3] = 0.0f;
	}
	memcpy(model->indices, tris.data(), sizeof(uint32_t) * model->numIndices);
	return model;
}

static void stlDestroy(objzModel *model)
{
	delete [] (uint32_t *)model->indices;
	delete [] model->meshes;
	delete [] model->objects;
	delete [] (ModelVertex *)model->vertices;
	delete model;
}

struct ModelLoadThreadArgs
{
	char filename[256];
};

static void modelLoadThread(ModelLoadThreadArgs args)
{
	s_model.data = nullptr;
	char basePath[256] = { 0 };
	const char *lastSlash = strrchr(args.filename, '/');
	if (!lastSlash)
		lastSlash = strrchr(args.filename, '\\');
	if (lastSlash) {
		for (int i = 0;; i++) {
			basePath[i] = args.filename[i];
			if (&args.filename[i] == lastSlash)
				break;
		}
	}
	const bx::StringView ext = bx::FilePath(args.filename).getExt();
	if (bx::strCmpI(ext, ".fbx") == 0) {
		objzModel *model = fbxLoad(args.filename, basePath);
		if (!model) {
			fprintf(stderr, "Error loading '%s'\n", args.filename);
			setErrorMessage("Error loading '%s'\n", args.filename);
			s_model.status = ModelStatus::NotLoaded;
			return;
		}
		s_model.data = model;
		s_model.destroyModelData = fbxDestroy;
	} else if (bx::strCmpI(ext, ".glb") == 0 || bx::strCmpI(ext, ".gltf") == 0) {
		objzModel *model = gltfLoad(args.filename, basePath);
		if (!model) {
			fprintf(stderr, "Error loading '%s'\n", args.filename);
			setErrorMessage("Error loading '%s'\n", args.filename);
			s_model.status = ModelStatus::NotLoaded;
			return;
		}
		s_model.data = model;
		s_model.destroyModelData = gltfDestroy;
	} else if (bx::strCmpI(ext, ".obj") == 0) {
		objz_setIndexFormat(OBJZ_INDEX_FORMAT_U32);
		objz_setVertexFormat(sizeof(ModelVertex), offsetof(ModelVertex, pos), offsetof(ModelVertex, texcoord), offsetof(ModelVertex, normal));
		objzModel *model = objz_load(args.filename);
		if (!model) {
			fprintf(stderr, "%s\n", objz_getError());
			setErrorMessage("Error loading' %s'\n%s\n", args.filename, objz_getError());
			s_model.status = ModelStatus::NotLoaded;
			return;
		}
		if (objz_getError()) // Print warnings.
			printf("%s\n", objz_getError());
		s_model.data = model;
		s_model.destroyModelData = objz_destroy;
		for (uint32_t i = 0; i < model->numVertices; i++) {
			auto v = &((ModelVertex *)model->vertices)[i];
			v->texcoord[1] = 1.0f - v->texcoord[1];
		}
	} else if (bx::strCmpI(ext, ".stl") == 0) {
		objzModel *model = stlLoad(args.filename, basePath);
		if (!model) {
			fprintf(stderr, "Error loading '%s'\n", args.filename);
			setErrorMessage("Error loading '%s'\n", args.filename);
			s_model.status = ModelStatus::NotLoaded;
			return;
		}
		s_model.data = model;
		s_model.destroyModelData = stlDestroy;
	} else {
		abort();
	}
	uint32_t numWireframeVertices = 0;
	for (uint32_t i = 0; i < s_model.data->numMeshes; i++) {
		const objzMesh &mesh = s_model.data->meshes[i];
		const objzMaterial *mat = mesh.materialIndex == -1 ? nullptr : &s_model.data->materials[mesh.materialIndex];
		if (mat && mat->opacity < 1.0f)
			continue;
		numWireframeVertices += mesh.numIndices;
	}
	s_model.wireframeVertices.resize(numWireframeVertices);
	uint32_t currentWireframeVertex = 0;
	for (uint32_t i = 0; i < s_model.data->numMeshes; i++) {
		const objzMesh &mesh = s_model.data->meshes[i];
		const objzMaterial *mat = mesh.materialIndex == -1 ? nullptr : &s_model.data->materials[mesh.materialIndex];
		if (mat && mat->opacity < 1.0f)
			continue;
		for (uint32_t j = 0; j < mesh.numIndices / 3; j++) {
			WireframeVertex *dest = &s_model.wireframeVertices[currentWireframeVertex];
			for (uint32_t k = 0; k < 3; k++)
				dest[k].pos = ((const ModelVertex *)s_model.data->vertices)[((const uint32_t *)s_model.data->indices)[mesh.firstIndex + j * 3 + k]].pos;
			dest[0].barycentric = bx::Vec3(1.0f, 0.0f, 0.0f);
			dest[1].barycentric = bx::Vec3(0.0f, 1.0f, 0.0f);
			dest[2].barycentric = bx::Vec3(0.0f, 0.0f, 1.0f);
			currentWireframeVertex += 3;
		}
	}
	s_model.diffuseTextures.resize(s_model.data->numMaterials);
	s_model.emissionTextures.resize(s_model.data->numMaterials);
	for (uint32_t i = 0; i < s_model.data->numMaterials; i++) {
		const objzMaterial &mat = s_model.data->materials[i];
		s_model.diffuseTextures[i] = mat.diffuseTexture[0] ? textureLoadCached(basePath, mat.diffuseTexture) : UINT32_MAX;
		s_model.emissionTextures[i] = mat.emissionTexture[0] ? textureLoadCached(basePath, mat.emissionTexture) : UINT32_MAX;
	}
	s_model.status = ModelStatus::Finalizing;
}

void modelFinalize()
{
	if (s_model.status != ModelStatus::Finalizing)
		return;
	if (s_model.thread) {
		if (s_model.thread->joinable())
			s_model.thread->join();
		delete s_model.thread;
		s_model.thread = nullptr;
	}
	printf("   %u object%s\n", s_model.data->numObjects, s_model.data->numObjects > 1 ? "s" : "");
	printf("   %u mesh%s\n", s_model.data->numMeshes, s_model.data->numMeshes > 1 ? "es" : "");
	printf("   %u triangles\n", s_model.data->numIndices / 3);
	printf("   %u vertices\n", s_model.data->numVertices);
	textureCreateCachedTextures();
	s_model.aabb = AABB();
	s_model.centroid = bx::Vec3(0.0f, 0.0f, 0.0f);
	uint32_t centroidCount = 0;
	for (uint32_t i = 0; i < s_model.data->numVertices; i++) {
		const bx::Vec3 &pos = ((const ModelVertex *)s_model.data->vertices)[i].pos;
		s_model.aabb.addPoint(pos);
		if (!std::isnan(pos.x) && !std::isnan(pos.y) && !std::isnan(pos.z)) {
			s_model.centroid = bx::add(s_model.centroid, pos);
			centroidCount++;
		}
	}
	s_model.centroid = bx::mul(s_model.centroid, 1.0f / centroidCount);
	s_model.vb = bgfx::createVertexBuffer(bgfx::makeRef(s_model.data->vertices, s_model.data->numVertices * sizeof(ModelVertex)), ModelVertex::layout);
	s_model.ib = bgfx::createIndexBuffer(bgfx::makeRef(s_model.data->indices, s_model.data->numIndices * sizeof(uint32_t)), BGFX_BUFFER_INDEX32);
	s_model.wireframeVb = bgfx::createVertexBuffer(bgfx::makeRef(s_model.wireframeVertices.data(), uint32_t(s_model.wireframeVertices.size() * sizeof(WireframeVertex))), WireframeVertex::layout);
	resetCamera();
	g_options.shadeMode = ShadeMode::Flat;
	g_options.wireframeMode = WireframeMode::Triangles;
	s_model.status = ModelStatus::Loaded;
}

static bool modelCanOpen()
{
	if (s_model.status == ModelStatus::Loading || s_model.status == ModelStatus::Finalizing)
		return false;
	if (!(atlasIsNotGenerated() || atlasIsReady()))
		return false;
	return true;
}

void modelOpen(const char *filename)
{
	if (!modelCanOpen())
		return;
	modelDestroy();
	s_model.status = ModelStatus::Loading;
	char windowTitle[256];
	snprintf(windowTitle, sizeof(windowTitle), "%s - %s\n", WINDOW_TITLE, filename);
	glfwSetWindowTitle(g_window, windowTitle);
	printf("Loading '%s'\n", filename);
	ModelLoadThreadArgs args;
	bx::strCopy(args.filename, sizeof(args.filename), filename);
	s_model.thread = new std::thread(modelLoadThread, args);
}

void modelOpenDialog()
{
	if (!modelCanOpen())
		return;
	nfdchar_t *filename = nullptr;
	nfdresult_t result = NFD_OpenDialog("fbx,glb,gltf,obj,stl", nullptr, &filename);
	if (result != NFD_OKAY)
		return;
	modelOpen(filename);
	free(filename);
}

void modelDestroy()
{
	textureDestroyCache();
	atlasDestroy();
	if (s_model.thread) {
		if (s_model.thread->joinable())
			s_model.thread->join();
		delete s_model.thread;
		s_model.thread = nullptr;
	}
	if (s_model.data) {
		s_model.destroyModelData(s_model.data);
		s_model.data = nullptr;
	}
	if (bgfx::isValid(s_model.vb)) {
		bgfx::destroy(s_model.vb);
		bgfx::destroy(s_model.ib);
		bgfx::destroy(s_model.wireframeVb);
		s_model.vb = BGFX_INVALID_HANDLE;
		s_model.ib = BGFX_INVALID_HANDLE;
		s_model.wireframeVb = BGFX_INVALID_HANDLE;
	}
	glfwSetWindowTitle(g_window, WINDOW_TITLE);
	s_model.status = ModelStatus::NotLoaded;
}

void modelRender(const float *view, const float *projection)
{
	if (s_model.status != ModelStatus::Loaded)
		return;
	float transform[16];
	if (s_model.rightHandedAxis)
		memcpy(transform, s_rightHandedAxisMatrix, sizeof(float) * 16);
	else
		bx::mtxIdentity(transform);
	float scaleMatrix[16];
	bx::mtxScale(scaleMatrix, s_model.scale);
	float modelMatrix[16];
	bx::mtxMul(modelMatrix, transform, scaleMatrix);
	bgfx::setViewTransform(kModelView, view, projection);
	bgfx::setViewTransform(kModelTransparentView, view, projection);
	const bool renderCharts = g_options.shadeMode == ShadeMode::Charts && atlasIsReady();
	if (g_options.shadeMode != ShadeMode::Charts || renderCharts) {
		const float lightDir[] = { view[2], view[6], view[10], 0.0f };
		for (uint32_t i = 0; i < s_model.data->numMeshes; i++) {
			const objzMesh &mesh = s_model.data->meshes[i];
			const objzMaterial *mat = mesh.materialIndex == -1 ? nullptr : &s_model.data->materials[mesh.materialIndex];
			// When rendering charts, emissive and transparent meshes won't be rendered, so do that here.
			const bool emissive = mat ? mat->emission[0] > 0.0f || mat->emission[1] > 0.0f || mat->emission[2] > 0.0f : false;
			const bool transparent = mat ? mat->opacity < 1.0f : false;
			if (renderCharts && !emissive && !transparent)
				continue;
			if (atlasIsReady()) {
				bgfx::setIndexBuffer(atlasGetIb(), mesh.firstIndex, mesh.numIndices);
				bgfx::setVertexBuffer(0, atlasGetVb());
			} else {
				bgfx::setIndexBuffer(s_model.ib, mesh.firstIndex, mesh.numIndices);
				bgfx::setVertexBuffer(0, s_model.vb);
			}
			uint64_t state = BGFX_STATE_DEFAULT;
			if (!s_model.clockwiseFaceWinding)
				state = (state & ~BGFX_STATE_CULL_CW) | BGFX_STATE_CULL_CCW;
			if (transparent)
				state |= BGFX_STATE_BLEND_ALPHA;
			bgfx::setState(state);
			bgfx::setTransform(modelMatrix);
			bgfx::setUniform(s_model.u_lightDir, lightDir);
			if (!mat) {
				const float diffuse[] = { 0.5f, 0.5f, 0.5f, 1.0f };
				const float emission[] = { 0.0f, 0.0f, 0.0f, 0.0f };
				bgfx::setUniform(s_model.u_diffuse, diffuse);
				bgfx::setUniform(s_model.u_emission, emission);
			} else {
				const float diffuse[] = { mat->diffuse[0], mat->diffuse[1], mat->diffuse[2], mat->opacity };
				const float emission[] = { mat->emission[0], mat->emission[1], mat->emission[2], mat->opacity };
				bgfx::setUniform(s_model.u_diffuse, diffuse);
				bgfx::setUniform(s_model.u_emission, emission);
			}
			float shade_diffuse_emission[4];
			shade_diffuse_emission[1] = DIFFUSE_COLOR;
			shade_diffuse_emission[2] = EMISSION_COLOR;
			if (g_options.shadeMode == ShadeMode::Lightmap)
				shade_diffuse_emission[0] = (float)SHADE_LIGHTMAP;
			else if (g_options.shadeMode == ShadeMode::LightmapOnly)
				shade_diffuse_emission[0] = (float)SHADE_LIGHTMAP_ONLY;
			else
				shade_diffuse_emission[0] = (float)SHADE_FLAT;
			bgfx::TextureHandle diffuseTexture = BGFX_INVALID_HANDLE;
			bgfx::TextureHandle emissionTexture = BGFX_INVALID_HANDLE;
			if (mat) {
				diffuseTexture = textureGetHandle(s_model.diffuseTextures[mesh.materialIndex]);
				emissionTexture = textureGetHandle(s_model.emissionTextures[mesh.materialIndex]);
			}
			if (bgfx::isValid(diffuseTexture))
				shade_diffuse_emission[1] = DIFFUSE_TEXTURE;
			if (bgfx::isValid(emissionTexture))
				shade_diffuse_emission[2] = EMISSION_TEXTURE;
			bgfx::setUniform(s_model.u_shade_diffuse_emission, shade_diffuse_emission);
			bgfx::setTexture(0, s_model.s_diffuse, bgfx::isValid(diffuseTexture) ? diffuseTexture : s_model.u_dummyTexture);
			bgfx::setTexture(1, s_model.s_emission, bgfx::isValid(emissionTexture) ? emissionTexture : s_model.u_dummyTexture);
			if (g_options.shadeMode == ShadeMode::Lightmap || g_options.shadeMode == ShadeMode::LightmapOnly)
				bgfx::setTexture(2, s_model.s_lightmap, bakeGetLightmap(), bakeGetLightmapSamplerFlags());
			else
				bgfx::setTexture(2, s_model.s_lightmap, s_model.u_dummyTexture);
			bgfx::submit(transparent ? kModelTransparentView : kModelView, s_model.materialProgram);
		}
	}
	if (renderCharts) {
		uint64_t state = BGFX_STATE_DEFAULT;
		if (!s_model.clockwiseFaceWinding)
			state = (state & ~BGFX_STATE_CULL_CW) | BGFX_STATE_CULL_CCW;
		atlasRenderCharts(modelMatrix, state);
	}
	if (g_options.wireframe) {
		if (g_options.wireframeMode == WireframeMode::Triangles) {
			const float color[] = { 0.0f, 0.0f, 0.0f, 0.75f };
			bgfx::setUniform(s_model.u_color, color);
			setWireframeThicknessUniform(1.5f);
			bgfx::setState(BGFX_STATE_WRITE_RGB | BGFX_STATE_WRITE_A | BGFX_STATE_WRITE_Z | BGFX_STATE_DEPTH_TEST_LEQUAL | BGFX_STATE_CULL_CW | BGFX_STATE_BLEND_ALPHA | BGFX_STATE_MSAA);
			bgfx::setTransform(modelMatrix);
			bgfx::setVertexBuffer(0, s_model.wireframeVb);
			bgfx::submit(kModelView, getWireframeProgram(), 1);
		} else {
			atlasRenderChartsWireframe(modelMatrix);
		}
	}
}

void modelShowGuiMenu()
{
	ImGui::Checkbox("Right-handed axis", &s_model.rightHandedAxis);
	ImGui::Checkbox("Clockwise face winding", &s_model.clockwiseFaceWinding);
	ImGui::PushItemWidth(100.0f);
	ImGui::InputFloat("Scale", &s_model.scale, 0.01f, 0.1f);
	ImGui::PopItemWidth();
	s_model.scale = bx::max(0.001f, s_model.scale);
}

void modelShowGuiWindow()
{
	const ImGuiWindowFlags progressWindowFlags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings;
	if (s_model.status == ModelStatus::Loading) {
		ImGui::SetNextWindowPos(ImVec2(g_windowSize[0] * 0.5f, g_windowSize[1] * 0.5f), ImGuiCond_Always, ImVec2(0.5f, 0.5f));
		if (ImGui::Begin("##modelProgress", nullptr, progressWindowFlags)) {
			ImGui::AlignTextToFramePadding();
			ImGui::Text("Loading model");
			ImGui::SameLine();
			ImGui::Spinner("##modelSpinner");
			ImGui::End();
		}
	}
}

AABB modelGetAABB()
{
	return s_model.aabb;
}

const objzModel *modelGetData()
{
	return s_model.data;
}

bx::Vec3 modelGetCentroid()
{
	bx::Vec3 centroid(s_model.centroid);
	if (s_model.rightHandedAxis)
		centroid = bx::mul(centroid, s_rightHandedAxisMatrix);
	return bx::mul(centroid, s_model.scale);
}

float modelGetScale()
{
	return s_model.scale;
}

bgfx::ShaderHandle modelGet_vs_model()
{
	return s_model.vs_model;
}

bool modelIsLoaded()
{
	return s_model.status == ModelStatus::Loaded;
}

static bool modelSampleTexture(uint32_t textureIndex, const float *uv, bx::Vec3 *color)
{
	if (textureIndex == UINT32_MAX)
		return false;
	const CachedTexture &texture = s_textureCache[textureIndex];
	if (!texture.data.mem)
		return false;
	const uint32_t x = uint32_t(uv[0] * texture.data.sampleWidth) % texture.data.sampleWidth;
	const uint32_t y = uint32_t(uv[1] * texture.data.sampleHeight) % texture.data.sampleHeight;
	const uint8_t *rgb = &texture.data.sampleData[(x + y * texture.data.sampleWidth) * texture.data.numComponents];
	if (texture.data.numComponents == 1)
		*color = bx::Vec3(rgb[0] / 255.0f);
	else
		*color = bx::Vec3(rgb[0] / 255.0f, rgb[1] / 255.0f, rgb[2] / 255.0f);
	return true;
}

bool modelSampleMaterialDiffuse(const objzMaterial *mat, const float *uv, bx::Vec3 *color)
{
	return modelSampleTexture(s_model.diffuseTextures[mat - s_model.data->materials], uv, color);
}

bool modelSampleMaterialEmission(const objzMaterial *mat, const float *uv, bx::Vec3 *color)
{
	return modelSampleTexture(s_model.emissionTextures[mat - s_model.data->materials], uv, color);
}
