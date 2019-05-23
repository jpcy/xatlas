/*
xatlas
https://github.com/jpcy/xatlas
Copyright (c) 2018 Jonathan Young

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#include <mutex>
#include <thread>
#include <bx/filepath.h>
#include <bx/string.h>
#include <cgltf.h>
#include <GLFW/glfw3.h>
#include <imgui/imgui.h>
#include <nativefiledialog/nfd.h>
#include <stb_image.h>
#include <stb_image_resize.h>
#include "shaders/shared.h"
#include "viewer.h"

bgfx::VertexDecl ModelVertex::decl;

struct ModelStatus
{
	enum Enum
	{
		NotLoaded,
		Loading,
		Finalizing,
		Loaded
	};

	Enum get()
	{
		m_lock.lock();
		Enum result = m_value;
		m_lock.unlock();
		return result;
	}

	void set(Enum value)
	{
		m_lock.lock();
		m_value = value;
		m_lock.unlock();
	}

private:
	std::mutex m_lock;
	Enum m_value = NotLoaded;
};

struct
{
	ModelStatus status;
	std::thread *thread = nullptr;
	objzModel *data = nullptr;
	bool isGltf = false;
	std::vector<uint32_t> diffuseTextures;
	std::vector<uint32_t> emissionTextures;
	AABB aabb;
	bx::Vec3 centroid = bx::Vec3(0.0f, 0.0f, 0.0f);
	bgfx::VertexBufferHandle vb = BGFX_INVALID_HANDLE;
	bgfx::IndexBufferHandle ib = BGFX_INVALID_HANDLE;
	bgfx::IndexBufferHandle wireframeIb = BGFX_INVALID_HANDLE;
	float scale = 1.0f;
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
#if _MSC_VER
	FILE *f;
	if (fopen_s(&f, fullFilename, "rb") != 0)
		f = nullptr;
#else
	FILE *f = fopen(fullFilename, "rb");
#endif
	if (!f) {
		fprintf(stderr, "Error opening '%s'\n", fullFilename);
		return td;
	}
	fseek(f, 0, SEEK_END);
	const long length = ftell(f);
	fseek(f, 0, SEEK_SET);
	std::vector<uint8_t> fileData;
	fileData.resize(length);
	fread(fileData.data(), 1, (size_t)length, f);
	fclose(f);
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
		if (!texture.data.mem)
			texture.handle = BGFX_INVALID_HANDLE;
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

static bool gltfAnyNodeInHierarchyHasMesh(const cgltf_node *node)
{
	if (node->mesh)
		return true;
	for (cgltf_size ci = 0; ci < node->children_count; ci++)
		return gltfAnyNodeInHierarchyHasMesh(node->children[ci]);
	return false;
}

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

static void gltfPopulateMeshData(const cgltf_node *node, const cgltf_material *firstMaterial, objzModel *model, uint32_t currentObject, uint32_t &currentMesh, uint32_t &firstMeshIndex, uint32_t &firstMeshVertex)
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
			// Create mesh.
			assert(currentMesh < model->numMeshes);
			objzMesh &mesh = model->meshes[currentMesh];
			mesh.materialIndex = uint32_t(primitive.material - firstMaterial);
			mesh.firstIndex = firstMeshIndex;
			mesh.numIndices = (uint32_t)aindices->count;
			currentMesh++;
			objzObject &object = model->objects[currentObject];
			object.numMeshes++;
			object.numIndices += (uint32_t)aindices->count;
			object.numVertices += (uint32_t)apositions->count;
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
		model->numObjects++;
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
		// Create mesh data.
		gltfPopulateMeshData(&node, gltfData->materials, model, currentObject, currentMesh, firstMeshIndex, firstMeshVertex);
		currentObject++;
	}
	// Materials.
	model->materials = new objzMaterial[model->numMaterials];
	for (uint32_t i = 0; i < model->numMaterials; i++) {
		const cgltf_material &sourceMat = gltfData->materials[i];
		const cgltf_texture *diffuse = sourceMat.pbr_metallic_roughness.base_color_texture.texture;
		const cgltf_texture *emission = sourceMat.emissive_texture.texture;
		objzMaterial &destMat = model->materials[i];
		memset(&destMat, 0, sizeof(destMat));
		memcpy(destMat.diffuse, sourceMat.pbr_metallic_roughness.base_color_factor, sizeof(float) * 3);
		if (diffuse)
			bx::strCopy(destMat.diffuseTexture, sizeof(destMat.diffuseTexture), diffuse->image->uri);
		if (emission)
			bx::strCopy(destMat.emissionTexture, sizeof(destMat.emissionTexture), emission->image->uri);
	}
	cgltf_free(gltfData);
	return model;
}

void gltfDestroy(objzModel *model)
{
	delete [] (uint32_t *)model->indices;
	delete [] model->meshes;
	delete [] model->objects;
	delete [] (ModelVertex *)model->vertices;
	delete [] model->materials;
	delete model;
}

struct ModelLoadThreadArgs
{
	char filename[256];
};

static void modelLoadThread(ModelLoadThreadArgs args)
{
	s_model.data = nullptr;
	s_model.isGltf = false;
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
	if (bx::strCmpI(ext, ".glb") == 0 || bx::strCmpI(ext, ".gltf") == 0) {
		objzModel *model = gltfLoad(args.filename, basePath);
		if (!model) {
			fprintf(stderr, "Error loading %s\n", args.filename);
			setErrorMessage("Error loading %s\n", args.filename);
			s_model.status.set(ModelStatus::NotLoaded);
			return;
		}
		s_model.data = model;
		s_model.isGltf = true;
	} else {
		objz_setIndexFormat(OBJZ_INDEX_FORMAT_U32);
		objz_setVertexFormat(sizeof(ModelVertex), offsetof(ModelVertex, pos), offsetof(ModelVertex, texcoord), offsetof(ModelVertex, normal));
		objzModel *model = objz_load(args.filename);
		if (!model) {
			fprintf(stderr, "%s\n", objz_getError());
			setErrorMessage(objz_getError());
			s_model.status.set(ModelStatus::NotLoaded);
			return;
		}
		if (objz_getError()) // Print warnings.
			printf("%s\n", objz_getError());
		s_model.data = model;
		for (uint32_t i = 0; i < model->numVertices; i++) {
			auto v = &((ModelVertex *)model->vertices)[i];
			v->texcoord[1] = 1.0f - v->texcoord[1];
		}
	}
	s_model.diffuseTextures.resize(s_model.data->numMaterials);
	s_model.emissionTextures.resize(s_model.data->numMaterials);
	for (uint32_t i = 0; i < s_model.data->numMaterials; i++) {
		const objzMaterial &mat = s_model.data->materials[i];
		s_model.diffuseTextures[i] = mat.diffuseTexture[0] ? textureLoadCached(basePath, mat.diffuseTexture) : UINT32_MAX;
		s_model.emissionTextures[i] = mat.emissionTexture[0] ? textureLoadCached(basePath, mat.emissionTexture) : UINT32_MAX;
	}
	s_model.status.set(ModelStatus::Finalizing);
}

void modelFinalize()
{
	if (s_model.status.get() != ModelStatus::Finalizing)
		return;
	if (s_model.thread) {
		if (s_model.thread->joinable())
			s_model.thread->join();
		delete s_model.thread;
		s_model.thread = nullptr;
	}
	textureCreateCachedTextures();
	s_model.aabb = AABB();
	s_model.centroid = bx::Vec3(0.0f, 0.0f, 0.0f);
	for (uint32_t i = 0; i < s_model.data->numVertices; i++) {
		const bx::Vec3 &pos = ((const ModelVertex *)s_model.data->vertices)[i].pos;
		s_model.aabb.addPoint(pos);
		s_model.centroid = bx::add(s_model.centroid, pos);
	}
	s_model.centroid = bx::mul(s_model.centroid, 1.0f / s_model.data->numVertices);
	s_model.vb = bgfx::createVertexBuffer(bgfx::makeRef(s_model.data->vertices, s_model.data->numVertices * sizeof(ModelVertex)), ModelVertex::decl);
	s_model.ib = bgfx::createIndexBuffer(bgfx::makeRef(s_model.data->indices, s_model.data->numIndices * sizeof(uint32_t)), BGFX_BUFFER_INDEX32);
	const uint32_t numWireframeIndices = bgfx::topologyConvert(bgfx::TopologyConvert::TriListToLineList, nullptr, 0, s_model.data->indices, s_model.data->numIndices, true);
	const bgfx::Memory *wireframeIndices = bgfx::alloc(numWireframeIndices * sizeof(uint32_t));
	bgfx::topologyConvert(bgfx::TopologyConvert::TriListToLineList, wireframeIndices->data, wireframeIndices->size, s_model.data->indices, s_model.data->numIndices, true);
	s_model.wireframeIb = bgfx::createIndexBuffer(wireframeIndices, BGFX_BUFFER_INDEX32);
	resetCamera();
	g_options.shadeMode = ShadeMode::Flat;
	g_options.wireframeMode = WireframeMode::Triangles;
	s_model.status.set(ModelStatus::Loaded);
}

void modelOpenDialog()
{
	if (s_model.status.get() == ModelStatus::Loading || s_model.status.get() == ModelStatus::Finalizing)
		return;
	if (!(atlasIsNotGenerated() || atlasIsReady()))
		return;
	nfdchar_t *filename = nullptr;
	nfdresult_t result = NFD_OpenDialog("glb,gltf,obj", nullptr, &filename);
	if (result != NFD_OKAY)
		return;
	modelDestroy();
	s_model.status.set(ModelStatus::Loading);
	char windowTitle[256];
	snprintf(windowTitle, sizeof(windowTitle), "%s - %s\n", WINDOW_TITLE, filename);
	glfwSetWindowTitle(g_window, windowTitle);
	printf("Loading '%s'\n", filename);
	ModelLoadThreadArgs args;
	bx::strCopy(args.filename, sizeof(args.filename), filename);
	s_model.thread = new std::thread(modelLoadThread, args);
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
		if (s_model.isGltf)
			gltfDestroy(s_model.data);
		else
			objz_destroy(s_model.data);
		s_model.data = nullptr;
	}
	if (bgfx::isValid(s_model.vb)) {
		bgfx::destroy(s_model.vb);
		bgfx::destroy(s_model.ib);
		bgfx::destroy(s_model.wireframeIb);
		s_model.vb = BGFX_INVALID_HANDLE;
		s_model.ib = BGFX_INVALID_HANDLE;
		s_model.wireframeIb = BGFX_INVALID_HANDLE;
	}
	glfwSetWindowTitle(g_window, WINDOW_TITLE);
	s_model.status.set(ModelStatus::NotLoaded);
}

void modelSetMaterialTexturesAndUniforms(int32_t materialIndex, uint8_t stageOffset)
{
	const objzMaterial *mat = materialIndex == -1 ? nullptr : &s_model.data->materials[materialIndex];
	if (!mat) {
		const float diffuse[] = { 0.5f, 0.5f, 0.5f, 1.0f };
		const float emission[] = { 0.0f, 0.0f, 0.0f, 0.0f };
		bgfx::setUniform(s_model.u_diffuse, diffuse);
		bgfx::setUniform(s_model.u_emission, emission);
	} else {
		const float diffuse[] = { mat->diffuse[0], mat->diffuse[1], mat->diffuse[2], 1.0f };
		const float emission[] = { mat->emission[0], mat->emission[1], mat->emission[2], 1.0f };
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
		diffuseTexture = textureGetHandle(s_model.diffuseTextures[materialIndex]);
		emissionTexture = textureGetHandle(s_model.emissionTextures[materialIndex]);
	}
	if (bgfx::isValid(diffuseTexture))
		shade_diffuse_emission[1] = DIFFUSE_TEXTURE;
	if (bgfx::isValid(emissionTexture))
		shade_diffuse_emission[2] = EMISSION_TEXTURE;
	bgfx::setUniform(s_model.u_shade_diffuse_emission, shade_diffuse_emission);
	bgfx::setTexture(stageOffset + 0, s_model.s_diffuse, bgfx::isValid(diffuseTexture) ? diffuseTexture : s_model.u_dummyTexture);
	bgfx::setTexture(stageOffset + 1, s_model.s_emission, bgfx::isValid(emissionTexture) ? emissionTexture : s_model.u_dummyTexture);
}

void modelRender(const float *view, const float *projection)
{
	if (s_model.status.get() != ModelStatus::Loaded)
		return;
	float modelMatrix[16];
	bx::mtxScale(modelMatrix, s_model.scale);
	bgfx::setViewTransform(kModelView, view, projection);
	const bool renderCharts = g_options.shadeMode == ShadeMode::Charts && atlasIsReady();
	if (g_options.shadeMode != ShadeMode::Charts || renderCharts) {
		const float lightDir[] = { view[2], view[6], view[10], 0.0f };
		for (uint32_t i = 0; i < s_model.data->numMeshes; i++) {
			const objzMesh &mesh = s_model.data->meshes[i];
			const objzMaterial *mat = mesh.materialIndex == -1 ? nullptr : &s_model.data->materials[mesh.materialIndex];
			// When rendering charts, emissive meshes won't be rendered, so do that here.
			const bool emissive = mat ? mat->emission[0] > 0.0f || mat->emission[1] > 0.0f || mat->emission[2] > 0.0f : false;
			if (renderCharts && !emissive)
				continue;
			if (atlasIsReady()) {
				bgfx::setIndexBuffer(atlasGetIb(), mesh.firstIndex, mesh.numIndices);
				bgfx::setVertexBuffer(0, atlasGetVb());
			} else {
				bgfx::setIndexBuffer(s_model.ib, mesh.firstIndex, mesh.numIndices);
				bgfx::setVertexBuffer(0, s_model.vb);
			}
			bgfx::setState(BGFX_STATE_DEFAULT);
			bgfx::setTransform(modelMatrix);
			bgfx::setUniform(s_model.u_lightDir, lightDir);
			modelSetMaterialTexturesAndUniforms(mesh.materialIndex);
			if (g_options.shadeMode == ShadeMode::Lightmap || g_options.shadeMode == ShadeMode::LightmapOnly)
				bgfx::setTexture(2, s_model.s_lightmap, bakeGetLightmap(), bakeGetLightmapSamplerFlags());
			else
				bgfx::setTexture(2, s_model.s_lightmap, s_model.u_dummyTexture);
			bgfx::submit(kModelView, s_model.materialProgram);
		}
	}
	if (renderCharts)
		atlasRenderCharts(modelMatrix);
	if (g_options.wireframe) {
		if (g_options.wireframeMode == WireframeMode::Triangles) {
			const float color[] = { 1.0f, 1.0f, 1.0f, 0.5f };
			bgfx::setUniform(s_model.u_color, color);
			bgfx::setState(BGFX_STATE_WRITE_RGB | BGFX_STATE_WRITE_A | BGFX_STATE_PT_LINES | BGFX_STATE_BLEND_ALPHA);
			bgfx::setTransform(modelMatrix);
			bgfx::setIndexBuffer(s_model.wireframeIb);
			bgfx::setVertexBuffer(0, s_model.vb);
			bgfx::submit(kModelView, getColorProgram());
		} else {
			atlasRenderChartsWireframe(modelMatrix);
		}
	}
}

void modelShowGuiOptions()
{
	ImGui::Text("%u objects", s_model.data->numObjects);
	ImGui::Text("%u meshes", s_model.data->numMeshes);
	ImGui::Text("%u vertices", s_model.data->numVertices);
	ImGui::Text("%u triangles", s_model.data->numIndices / 3);
	ImGui::InputFloat("Model scale", &s_model.scale, 0.01f, 0.1f);
	s_model.scale = bx::max(0.001f, s_model.scale);
}

void modelShowGuiWindow(int progressDots)
{
	ImGuiIO &io = ImGui::GetIO();
	const ImGuiWindowFlags progressWindowFlags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings;
	if (s_model.status.get() == ModelStatus::Loading) {
		ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y * 0.5f), ImGuiCond_Always, ImVec2(0.5f, 0.5f));
		if (ImGui::Begin("##modelProgress", nullptr, progressWindowFlags)) {
			ImGui::Text("Loading model");
			for (int i = 0; i < 3; i++) {
				ImGui::SameLine();
				ImGui::Text(i < progressDots ? "." : " ");
			}
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
	return bx::mul(s_model.centroid, s_model.scale);
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
	return s_model.status.get() == ModelStatus::Loaded;
}

static bool modelSampleTexture(uint32_t textureIndex, const float *uv, bx::Vec3 *color)
{
	if (textureIndex == UINT32_MAX)
		return false;
	const CachedTexture &texture = s_textureCache[textureIndex];
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
