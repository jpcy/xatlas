/*
https://github.com/jpcy/objzero

MIT License

Copyright (c) 2018-2019 Jonathan Young

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
Ear clipping triangulation from tinyobjloader, also under MIT license.
https://github.com/syoyo/tinyobjloader
Copyright (c) 2012-2018 Syoyo Fujita and many contributors.
*/
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include "objzero.h"

#ifdef _MSC_VER
#define OBJZ_FOPEN(_file, _filename, _mode) { if (fopen_s(&_file, _filename, _mode) != 0) _file = NULL; }
#define OBJZ_STRICMP _stricmp
#define OBJZ_STRTOK(_str, _delim, _context) strtok_s(_str, _delim, _context)
#else
#include <strings.h>
#define OBJZ_FOPEN(_file, _filename, _mode) _file = fopen(_filename, _mode)
#define OBJZ_STRICMP strcasecmp
#define OBJZ_STRTOK(_str, _delim, _context) strtok(_str, _delim)
#endif

#define OBJZ_MAX_ERROR_LENGTH 1024
#define OBJZ_MAX_TOKEN_LENGTH 256
#define OBJZ_RAW_ARRAY_LEN(_x) (sizeof(_x) / sizeof((_x)[0]))
#define OBJZ_SMALLEST(_a, _b) ((_a) < (_b) ? (_a) : (_b))
#define OBJZ_LARGEST(_a, _b) ((_a) > (_b) ? (_a) : (_b))

static char s_error[OBJZ_MAX_ERROR_LENGTH] = { 0 };
static objzReallocFunc s_realloc = NULL;
static uint32_t s_indexFormat = OBJZ_INDEX_FORMAT_AUTO;

typedef struct {
	size_t stride;
	size_t positionOffset;
	size_t texcoordOffset;
	size_t normalOffset;
} VertexFormat;

static VertexFormat s_vertexDecl = {
	.stride = sizeof(float) * (3 + 2 + 3),
	.positionOffset = 0,
	.texcoordOffset = sizeof(float) * 3,
	.normalOffset = sizeof(float) * (3 + 2)
};

static void *objz_realloc(void *_ptr, size_t _size, char *_file, int _line) {
	if (!_ptr && !_size)
		return NULL;
	void *result;
	if (s_realloc)
		result = s_realloc(_ptr, _size);
	else
		result = realloc(_ptr, _size);
	if (_size > 0 && !result) {
		fprintf(stderr, "Memory allocation failed %s %d\n", _file, _line);
		abort();
	}
	return result;
}

#define OBJZ_MALLOC(_size) objz_realloc(NULL, (_size), __FILE__, __LINE__)
#define OBJZ_REALLOC(_ptr, _size) objz_realloc((_ptr), (_size), __FILE__, __LINE__)
#define OBJZ_FREE(_ptr) objz_realloc((_ptr), 0, __FILE__, __LINE__)

static size_t strLength(const char *_str, size_t _size)
{
	const char *c = _str;
	size_t len = 0;
	while (*c != 0 && len < _size) {
		c++;
		len++;
	}
	return len;
}

static void strCopy(char *_dest, size_t _destSize, const char *_src, size_t _count)
{
	const size_t n = OBJZ_SMALLEST(_destSize - 1, strLength(_src, _count));
	memcpy(_dest, _src, n);
	_dest[n] = 0;
}

static void strConcat(char *_dest, size_t _destSize, const char *_src, size_t _count)
{
	const size_t start = strLength(_dest, _destSize);
	strCopy(&_dest[start], _destSize - start, _src, _count);
}

typedef struct {
	float x, y, z;
} vec3;

#define OBJZ_VEC3_ABS(_out, _v) (_out).x = fabsf((_v).x); (_out).y = fabsf((_v).y); (_out).z = fabsf((_v).z);
#define OBJZ_VEC3_ADD(_out, _a, _b) (_out).x = (_a).x + (_b).x; (_out).y = (_a).y + (_b).y; (_out).z = (_a).z + (_b).z;

#define OBJZ_VEC3_CROSS(_out, _a, _b)             \
	(_out).x = (_a).y * (_b).z - (_a).z * (_b).y; \
	(_out).y = (_a).z * (_b).x - (_a).x * (_b).z; \
	(_out).z = (_a).x * (_b).y - (_a).y * (_b).x;

#define OBJZ_VEC3_COPY(_out, _in) (_out).x = (_in).x; (_out).y = (_in).y; (_out).z = (_in).z;
#define OBJZ_VEC3_DOT(_out, _a, _b) (_out) = (_a).x * (_b).x + (_a).y * (_b).y + (_a).z * (_b).z;
#define OBJZ_VEC3_MUL(_out, _v, _s) (_out).x = (_v).x * _s; (_out).y = (_v).y * _s; (_out).z = (_v).z * _s;
#define OBJZ_VEC3_SET(_out, _x, _y, _z) (_out).x = (_x); (_out).y = (_y); (_out).z = (_z);
#define OBJZ_VEC3_SUB(_out, _a, _b) (_out).x = (_a).x - (_b).x; (_out).y = (_a).y - (_b).y; (_out).z = (_a).z - (_b).z;

static bool vec3Equal(const vec3 *_a, const vec3 *_b, float epsilon) {
	return fabsf(_a->x - _b->x) <= epsilon && fabsf(_a->y - _b->y) <= epsilon && fabsf(_a->z - _b->z) <= epsilon;
}

static void vec3Normalize(vec3 *_out, const vec3 *_in) {
	float len;
	OBJZ_VEC3_DOT(len, *_in, *_in);
	if (len > 0) {
		len = 1.0f / sqrtf(len);
		OBJZ_VEC3_MUL(*_out, *_in, len);
	} else {
		OBJZ_VEC3_COPY(*_out, *_in);
	}
}

static void appendError(const char *_format, ...) {
	va_list args;
	va_start(args, _format);
	char buffer[OBJZ_MAX_ERROR_LENGTH];
	vsnprintf(buffer, sizeof(buffer), _format, args);
	va_end(args);
	if (s_error[0]) {
		const char *newline = "\n";
		strConcat(s_error, sizeof(s_error), newline, 1);
	}
	strConcat(s_error, sizeof(s_error), buffer, strLength(buffer, sizeof(buffer)));
}

typedef struct {
	uint8_t *data;
	uint32_t length;
	uint32_t capacity;
	uint32_t elementSize;
	uint32_t initialCapacity;
} Array;

static void arrayInit(Array *_array, size_t _elementSize, uint32_t _initialCapacity) {
	_array->data = NULL;
	_array->length = _array->capacity = 0;
	_array->elementSize = (uint32_t)_elementSize;
	_array->initialCapacity = _initialCapacity;
}

static void arrayDestroy(Array *_array) {
	OBJZ_FREE(_array->data);
}

static void arrayAppend(Array *_array, const void *_element) {
	if (!_array->data) {
		_array->data = OBJZ_MALLOC(_array->elementSize * _array->initialCapacity);
		_array->capacity = _array->initialCapacity;
	} else if (_array->length == _array->capacity) {
		_array->capacity *= 2;
		_array->data = OBJZ_REALLOC(_array->data, _array->capacity * _array->elementSize);
	}
	memcpy(&_array->data[_array->length * _array->elementSize], _element, _array->elementSize);
	_array->length++;
}

#define OBJZ_ARRAY_ELEMENT(_array, _index) (void *)&(_array).data[(_array).elementSize * (_index)]

// Array: reallocates the buffer when full. The buffer is a contiguous area of memory.
// ChunkedArray: allocates another chunk of memory when full. Buffer is a linked list of chunks, not contiguous.
typedef struct {
	Array chunks;
	uint32_t elementsPerChunk;
	size_t elementSize;
	uint32_t length;
} ChunkedArray;

static void chunkedArrayInit(ChunkedArray *_array, size_t _elementSize, uint32_t _chunkLength) {
	arrayInit(&_array->chunks, sizeof(void *), 32);
	_array->elementsPerChunk = _chunkLength;
	_array->elementSize = _elementSize;
	_array->length = 0;
}

static void chunkedArrayDestroy(ChunkedArray *_array) {
	for (uint32_t i = 0; i < _array->chunks.length; i++) {
		void **chunk = OBJZ_ARRAY_ELEMENT(_array->chunks, i);
		OBJZ_FREE(*chunk);
	}
	arrayDestroy(&_array->chunks);
}

static void chunkedArrayAppend(ChunkedArray *_array, const void *_element) {
	if (_array->length >= _array->chunks.length * _array->elementsPerChunk) {
		void *newChunk = OBJZ_MALLOC(_array->elementsPerChunk * _array->elementSize);
		arrayAppend(&_array->chunks, &newChunk);
	}
	uint8_t **chunk = OBJZ_ARRAY_ELEMENT(_array->chunks, _array->length / _array->elementsPerChunk);
	memcpy(&(*chunk)[_array->elementSize * (_array->length % _array->elementsPerChunk)], _element, _array->elementSize);
	_array->length++;
}

static void *chunkedArrayElement(const ChunkedArray *_array, uint32_t _index) {
	uint8_t **chunk = OBJZ_ARRAY_ELEMENT(_array->chunks, _index / _array->elementsPerChunk);
	return &(*chunk)[_array->elementSize * (_index % _array->elementsPerChunk)];
}

typedef struct {
	char *buf;
	uint32_t line, column;
} Lexer;

typedef struct {
	char text[OBJZ_MAX_TOKEN_LENGTH];
	uint32_t line, column;
} Token;

static void initLexer(Lexer *_lexer) {
	_lexer->buf = NULL;
	_lexer->column = 1;
	_lexer->line = 0;
}

static bool isEol(const Lexer *_lexer) {
	return (_lexer->buf[0] == 0);
}

static bool isWhitespace(const Lexer *_lexer) {
	return (_lexer->buf[0] == ' ' || _lexer->buf[0] == '\t' || _lexer->buf[0] == '\r');
}

static void skipWhitespace(Lexer *_lexer) {
	for (;;) {
		if (isEol(_lexer))
			break;
		if (!isWhitespace(_lexer))
			break;
		_lexer->buf++;
		_lexer->column++;
	}
}

static void lexerSetLine(Lexer *_lexer, char *_buf) {
	_lexer->column = 1;
	_lexer->line++;
	_lexer->buf = _buf;
}

static void tokenize(Lexer *_lexer, Token *_token, bool includeWhitespace) {
	uint32_t i = 0;
	skipWhitespace(_lexer);
	_token->line = _lexer->line;
	_token->column = _lexer->column;
	for (;;) {
		if (isEol(_lexer) || (!includeWhitespace && isWhitespace(_lexer)))
			break;
		_token->text[i++] = _lexer->buf[0];
		_lexer->buf++;
		_lexer->column++;
	}
	_token->text[i] = 0;
}

static bool parseFloats(Lexer *_lexer, float *_result, uint32_t n) {
	Token token;
	for (uint32_t i = 0; i < n; i++) {
		tokenize(_lexer, &token, false);
		if (strLength(token.text, sizeof(token.text)) == 0) {
			appendError("(%u:%u) Error parsing float", token.line, token.column);
			return false;
		}
		_result[i] = (float)atof(token.text);
	}
	return true;
}

static bool skipTokens(Lexer *_lexer, int _n) {
	Token token;
	for (int i = 0; i < _n; i++) {
		tokenize(_lexer, &token, false);
		if (strLength(token.text, sizeof(token.text)) == 0) {
			appendError("(%u:%u) Error skipping tokens", token.line, token.column);
			return false;
		}
	}
	return true;
}

typedef struct {
	char *buffer;
	size_t length;
	size_t pos;
} File;

bool fileOpen(File *_file, const char *_filename) {
	FILE *handle;
	OBJZ_FOPEN(handle, _filename, "rb");
	if (!handle)
		return false;
	fseek(handle, 0, SEEK_END);
	_file->length = (size_t)ftell(handle);
	fseek(handle, 0, SEEK_SET);
	if (_file->length == 0) {
		fclose(handle);
		return false;
	}
	_file->pos = 0;
	_file->buffer = OBJZ_MALLOC(_file->length + 1);
	const size_t bytesRead = fread(_file->buffer, 1, _file->length, handle);
	fclose(handle);
	if (bytesRead < _file->length) {
		OBJZ_FREE(_file->buffer);
		return false;
	}
	_file->buffer[_file->length] = 0;
	return true;
}

void fileClose(File *_file) {
	OBJZ_FREE(_file->buffer);
}

char *fileReadLine(File *_file) {
	if (_file->buffer[_file->pos] == 0)
		return NULL; // eof
	char *start = &_file->buffer[_file->pos];
	// Find eol. Newline is replaced with a null terminator. Position is set to the start of the next line.
	for (;;) {
		char *c = &_file->buffer[_file->pos];
		if (*c == 0) {
			break;
		}
		if (*c == '\r')
			*c = 0; // Null terminate here, but keep reading until newline.
		else if (*c == '\n') {
			*c = 0;
			_file->pos++;
			break;
		}
		_file->pos++;
	}
	return start;
}

#define OBJZ_MAT_TOKEN_STRING 0
#define OBJZ_MAT_TOKEN_FLOAT  1

typedef struct {
	const char *name;
	uint32_t type;
	size_t offset;
	uint32_t n;
} MaterialProperty;

static MaterialProperty s_materialProperties[] = {
	{ "d", OBJZ_MAT_TOKEN_FLOAT, offsetof(objzMaterial, opacity), 1 },
	{ "Ka", OBJZ_MAT_TOKEN_FLOAT, offsetof(objzMaterial, ambient), 3 },
	{ "Kd", OBJZ_MAT_TOKEN_FLOAT, offsetof(objzMaterial, diffuse), 3 },
	{ "Ke", OBJZ_MAT_TOKEN_FLOAT, offsetof(objzMaterial, emission), 3 },
	{ "Ks", OBJZ_MAT_TOKEN_FLOAT, offsetof(objzMaterial, specular), 3 },
	{ "Ns", OBJZ_MAT_TOKEN_FLOAT, offsetof(objzMaterial, specularExponent), 1 },
	{ "bump", OBJZ_MAT_TOKEN_STRING, offsetof(objzMaterial, bumpTexture), 1 },
	{ "map_Bump", OBJZ_MAT_TOKEN_STRING, offsetof(objzMaterial, bumpTexture), 1 },
	{ "map_Ka", OBJZ_MAT_TOKEN_STRING, offsetof(objzMaterial, ambientTexture), 1 },
	{ "map_Kd", OBJZ_MAT_TOKEN_STRING, offsetof(objzMaterial, diffuseTexture), 1 },
	{ "map_Ke", OBJZ_MAT_TOKEN_STRING, offsetof(objzMaterial, emissionTexture), 1 },
	{ "map_Ks", OBJZ_MAT_TOKEN_STRING, offsetof(objzMaterial, specularTexture), 1 },
	{ "map_Ns", OBJZ_MAT_TOKEN_STRING, offsetof(objzMaterial, specularExponentTexture), 1 },
	{ "map_d", OBJZ_MAT_TOKEN_STRING, offsetof(objzMaterial, opacityTexture), 1 }
};

typedef struct {
	const char *name;
	uint32_t n;
} MaterialMapArg;

static MaterialMapArg s_materialMapArgs[] = {
	{ "-blendu", 1 },
	{ "-blendv", 1 },
	{ "-bm", 1 },
	{ "-boost", 1 },
	{ "-clamp", 1 },
	{ "-imfchan", 1 },
	{ "-mm", 2 },
	{ "-o", 3 },
	{ "-s", 3 },
	{ "-t", 3 },
	{ "-texres", 1 },
	{ "-type", 1 }
};

static void materialInit(objzMaterial *_mat) {
	memset(_mat, 0, sizeof(*_mat));
	_mat->diffuse[0] = _mat->diffuse[1] = _mat->diffuse[2] = 1;
	_mat->opacity = 1;
}

static bool loadMaterialFile(const char *_objFilename, const char *_materialName, Array *_materials) {
	char filename[256] = { 0 };
	const char *lastSlash = strrchr(_objFilename, '/');
	if (!lastSlash)
		lastSlash = strrchr(_objFilename, '\\');
	if (lastSlash) {
		for (int i = 0;; i++) {
			filename[i] = _objFilename[i];
			if (&_objFilename[i] == lastSlash)
				break;
		}
		strConcat(filename, sizeof(filename), _materialName, strLength(_materialName, OBJZ_MAX_TOKEN_LENGTH));
	} else
		strCopy(filename, sizeof(filename), _materialName, strLength(_materialName, OBJZ_MAX_TOKEN_LENGTH));
	File file;
	if (!fileOpen(&file, filename)) {
		// Treat missing material file as a warning, not an error.
		appendError("Failed to read material file '%s'", filename);
		return true;
	}
	Lexer lexer;
	initLexer(&lexer);
	Token token;
	objzMaterial mat;
	materialInit(&mat);
	bool result = false;
	for (;;) {
		char *line = fileReadLine(&file);
		if (!line)
			break;
		lexerSetLine(&lexer, line);
		tokenize(&lexer, &token, false);
		if (OBJZ_STRICMP(token.text, "newmtl") == 0) {
			tokenize(&lexer, &token, false);
			if (token.text[0] == 0) {
				appendError("(%u:%u) Expected name after 'newmtl'", token.line, token.column);
				goto cleanup;
			}
			if (mat.name[0] != 0)
				arrayAppend(_materials, &mat);
			materialInit(&mat);
			strCopy(mat.name, sizeof(mat.name), token.text, strLength(token.text, sizeof(token.text)));
		} else {
			for (size_t i = 0; i < OBJZ_RAW_ARRAY_LEN(s_materialProperties); i++) {
				const MaterialProperty *prop = &s_materialProperties[i];
				uint8_t *dest = &((uint8_t *)&mat)[prop->offset];
				if (OBJZ_STRICMP(token.text, prop->name) == 0) {
					if (prop->type == OBJZ_MAT_TOKEN_STRING) {
						Token argToken;
						for (int j = 0;; j++) {
							tokenize(&lexer, &argToken, false);
							if (argToken.text[0] == 0) {
								if (j == 0) {
									appendError("(%u:%u) Expected token after '%s'", token.line, token.column, prop->name);
									goto cleanup;
								}
								break;
							}
							bool match = false;
							for (size_t k = 0; k < OBJZ_RAW_ARRAY_LEN(s_materialMapArgs); k++) {
								const MaterialMapArg *arg = &s_materialMapArgs[k];
								if (OBJZ_STRICMP(argToken.text, arg->name) == 0) {
									match = true;
									skipTokens(&lexer, arg->n);
									break;
								}
							}
							if (!match)
								strCopy((char *)dest, OBJZ_NAME_MAX, argToken.text, strLength(argToken.text, sizeof(argToken.text)));
						}
					} else if (prop->type == OBJZ_MAT_TOKEN_FLOAT) {
						if (!parseFloats(&lexer, (float *)dest, prop->n))
							goto cleanup;
					}
					break;
				}
			}
		}
	}
	if (mat.name[0] != 0)
		arrayAppend(_materials, &mat);
	result = true;
cleanup:
	fileClose(&file);
	return result;
}

static uint32_t sdbmHash(const uint8_t *_data, uint32_t _size)
{
	uint32_t hash = 0;
	for (uint32_t i = 0; i < _size; i++)
		hash = (uint32_t)_data[i] + (hash << 6) + (hash << 16) - hash;
	return hash;
}

typedef struct {
	uint32_t object;
	uint32_t pos;
	uint32_t texcoord;
	uint32_t normal;
	uint32_t hashNext; // For hash collisions: next HashedVertex with the same hash.
} HashedVertex;

typedef struct {
	uint32_t *slots;
	uint32_t numSlots;
	Array vertices;
} VertexHashMap;

static void vertexHashMapInit(VertexHashMap *_map, uint32_t _initialCapacity) {
	_map->numSlots = (uint32_t)(_initialCapacity * 1.3f);
	_map->slots = OBJZ_MALLOC(sizeof(uint32_t) * _map->numSlots);
	for (uint32_t i = 0; i < _map->numSlots; i++)
		_map->slots[i] = UINT32_MAX;
	arrayInit(&_map->vertices, sizeof(HashedVertex), _initialCapacity);
}

static void vertexHashMapDestroy(VertexHashMap *_map) {
	OBJZ_FREE(_map->slots);
	arrayDestroy(&_map->vertices);
}

static uint32_t vertexHashMapInsert(VertexHashMap *_map, uint32_t _object, uint32_t _pos, uint32_t _texcoord, uint32_t _normal) {
	uint32_t hashData[4] = { 0 };
	hashData[0] = _object;
	hashData[1] = _pos;
	if (_texcoord != UINT32_MAX)
		hashData[2] = _texcoord;
	if (_normal != UINT32_MAX)
		hashData[3] = _normal;
	const uint32_t hash = sdbmHash((const uint8_t *)hashData, sizeof(hashData)) % _map->numSlots;
	uint32_t i = _map->slots[hash];
	while (i != UINT32_MAX) {
		const HashedVertex *v = OBJZ_ARRAY_ELEMENT(_map->vertices, i);
		if (v->object == _object && v->pos == _pos && v->texcoord == _texcoord && v->normal == _normal)
			return i;
		i = v->hashNext;
	}
	HashedVertex v;
	v.object = _object;
	v.pos = _pos;
	v.texcoord = _texcoord;
	v.normal = _normal;
	v.hashNext = _map->slots[hash];
	_map->slots[hash] = _map->vertices.length;
	arrayAppend(&_map->vertices, &v);
	return _map->slots[hash];
}

typedef struct {
	uint32_t normalIndex;
	uint32_t hashNext; // For hash collisions: next HashedNormal with the same hash.
} HashedNormal;

typedef struct {
	uint32_t *slots;
	uint32_t numSlots;
	Array hashedNormals;
	ChunkedArray *normals;
} NormalHashMap;

static void normalHashMapClear(NormalHashMap *_map) {
	for (uint32_t i = 0; i < _map->numSlots; i++)
		_map->slots[i] = UINT32_MAX;
	_map->hashedNormals.length = 0;
}

static void normalHashMapInit(NormalHashMap *_map, uint32_t _initialCapacity, ChunkedArray *_normals) {
	_map->numSlots = (uint32_t)(_initialCapacity * 1.3f);
	_map->slots = OBJZ_MALLOC(sizeof(uint32_t) * _map->numSlots);
	_map->normals = _normals;
	arrayInit(&_map->hashedNormals, sizeof(HashedNormal), _initialCapacity);
	normalHashMapClear(_map);
}

static void normalHashMapDestroy(NormalHashMap *_map) {
	OBJZ_FREE(_map->slots);
	arrayDestroy(&_map->hashedNormals);
}

static uint32_t normalHashMapInsert(NormalHashMap *_map, const vec3 *_normal) {
	uint32_t hashData[3] = { 0 };
	hashData[0] = (uint32_t)(_normal->x * 0.5f + 0.5f * 255);
	hashData[1] = (uint32_t)(_normal->y * 0.5f + 0.5f * 255);
	hashData[2] = (uint32_t)(_normal->z * 0.5f + 0.5f * 255);
	const uint32_t hash = sdbmHash((const uint8_t *)hashData, sizeof(hashData)) % _map->numSlots;
	uint32_t i = _map->slots[hash];
	while (i != UINT32_MAX) {
		const HashedNormal *n = OBJZ_ARRAY_ELEMENT(_map->hashedNormals, i);
		if (vec3Equal(chunkedArrayElement(_map->normals, n->normalIndex), _normal, FLT_EPSILON))
			return n->normalIndex;
		i = n->hashNext;
	}
	HashedNormal n;
	n.normalIndex = _map->normals->length;
	n.hashNext = _map->slots[hash];
	_map->slots[hash] = _map->hashedNormals.length;
	arrayAppend(&_map->hashedNormals, &n);
	chunkedArrayAppend(_map->normals, _normal);
	return n.normalIndex;
}

typedef struct {
	char name[OBJZ_NAME_MAX];
	uint32_t firstFace;
	uint32_t numFaces;
} TempObject;

typedef struct {
	uint32_t v;
	uint32_t vt;
	uint32_t vn;
} IndexTriplet;

typedef struct {
	int16_t materialIndex;
	uint16_t smoothingGroup; // 0 is off
	IndexTriplet indices[3];
} Face;

static bool parseVertexAttribIndices(Token *_token, int32_t *_out) {
	int32_t *v = &_out[0];
	int32_t *vt = &_out[1];
	int32_t *vn = &_out[2];
	*v = *vt = *vn = INT_MAX;
	if (strLength(_token->text, sizeof(_token->text)) == 0)
		return false; // Empty token.
	const char *delim = "/";
	char *start = _token->text;
	bool eol = false;
	// v
	char *end = strstr(start, delim);
	if (!end) {
		end = &_token->text[strLength(_token->text, sizeof(_token->text))];
		eol = true;
	} else if (end == start)
		return false; // Token is just a delimiter.
	*end = 0;
	*v = atoi(start);
	// vt
	if (eol)
		return true;  // No vt or vn.
	start = end + 1;
	if (*start == 0)
		return true; // No vt or vn.
	end = strstr(start, delim);
	bool skipNormal = false;
	if (!end) {
		// No delimiter, must be no normal, i.e. "v/vt".
		skipNormal = true;
		end = &_token->text[strLength(_token->text, sizeof(_token->text)) - 1];
	}
	*end = 0;
	if (start != end)
		*vt = atoi(start);
	// vn
	if (skipNormal)
		return true;
	start = end + 1;
	if (*start != 0)
		*vn = atoi(start);
	return true;
}

// Convert 1-indexed relative vertex attrib index to a 0-indexed absolute index.
static uint32_t fixVertexAttribIndex(int32_t _index, uint32_t _n) {
	if (_index == INT_MAX)
		return UINT32_MAX;
	// Handle relative index.
	if (_index < 0)
		return (uint32_t)(_index + _n);
	// Convert from 1-indexed to 0-indexed.
	return (uint32_t)(_index - 1);
}

// code from https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html
static int pnpoly(int nvert, float *vertx, float *verty, float testx, float testy) {
	bool c = false;
	for (int i = 0, j = nvert - 1; i < nvert; j = i++) {
		if (((verty[i] > testy) != (verty[j] > testy)) && (testx < (vertx[j] - vertx[i]) * (testy - verty[i]) / (verty[j] - verty[i]) + vertx[i]))
			c = !c;
	}
	return c;
}

// Ear clipping triangulation from tinyobjloader
// https://github.com/syoyo/tinyobjloader
static void triangulate(const Array *_indices, const ChunkedArray *_positions, Array *_tempIndices, ChunkedArray *_faces, int32_t _materialIndex, uint16_t _smoothingGroup) {
	// find the two axes to work in
	uint32_t axes[2] = { 1, 2 };
	for (uint32_t i = 0; i < _indices->length; i++) {
		const IndexTriplet *indices = (const IndexTriplet *)_indices->data;
		vec3 v[3];
		for (int j = 0; j < 3; j++)
			v[j] = *(const vec3 *)chunkedArrayElement(_positions, indices[(i + j) % _indices->length].v);
		vec3 edges[2];
		OBJZ_VEC3_SUB(edges[0], v[1], v[0]);
		OBJZ_VEC3_SUB(edges[1], v[2], v[1]);
		vec3 corner;
		OBJZ_VEC3_CROSS(corner, edges[0], edges[1]);
		OBJZ_VEC3_ABS(corner, corner);
		if (corner.x > FLT_EPSILON || corner.y > FLT_EPSILON || corner.z > FLT_EPSILON) {
			// found a corner
			if (!(corner.x > corner.y && corner.x > corner.z)) {
				axes[0] = 0;
				if (corner.z > corner.x && corner.z > corner.y)
					axes[1] = 1;
			}
			break;
		}
	}
	float area = 0;
	for (uint32_t i = 0; i < _indices->length; i++) {
		const IndexTriplet *i0 = OBJZ_ARRAY_ELEMENT(*_indices, (i + 0) % _indices->length);
		const IndexTriplet *i1 = OBJZ_ARRAY_ELEMENT(*_indices, (i + 1) % _indices->length);
		const float *v0 = chunkedArrayElement(_positions, i0->v);
		const float *v1 = chunkedArrayElement(_positions, i1->v);
		area += (v0[axes[0]] * v1[axes[1]] - v0[axes[1]] * v1[axes[0]]) * 0.5f;
	}
	// Copy vertices.
	Array *remainingIndices = _tempIndices;
	remainingIndices->length = 0;
	for (uint32_t i = 0; i < _indices->length; i++)
		arrayAppend(remainingIndices, OBJZ_ARRAY_ELEMENT(*_indices, i));
	// How many iterations can we do without decreasing the remaining vertices.
	uint32_t remainingIterations = remainingIndices->length;
	uint32_t previousRemainingIndices = remainingIndices->length;
	uint32_t guess_vert = 0;
	while (remainingIndices->length > 3 && remainingIterations > 0) {
		if (guess_vert >= remainingIndices->length)
			guess_vert -= remainingIndices->length;
		if (previousRemainingIndices != remainingIndices->length) {
			// The number of remaining vertices decreased. Reset counters.
			previousRemainingIndices = remainingIndices->length;
			remainingIterations = remainingIndices->length;
		} else {
			// We didn't consume a vertex on previous iteration, reduce the
			// available iterations.
			remainingIterations--;
		}
		IndexTriplet *ind[3];
		float vx[3];
		float vy[3];
		for (uint32_t i = 0; i < 3; i++) {
			ind[i] = OBJZ_ARRAY_ELEMENT(*remainingIndices, (guess_vert + i) % remainingIndices->length);
			const float *pos = (float *)chunkedArrayElement(_positions, ind[i]->v);
			vx[i] = pos[axes[0]];
			vy[i] = pos[axes[1]];
		}
		float edge0[2], edge1[2];
		edge0[0] = vx[1] - vx[0];
		edge0[1] = vy[1] - vy[0];
		edge1[0] = vx[2] - vx[1];
		edge1[1] = vy[2] - vy[1];
		const float cross = edge0[0] * edge1[1] - edge0[1] * edge1[0];
		// if an internal angle
		if (cross * area < 0.0f) {
			guess_vert += 1;
			continue;
		}
		// check all other verts in case they are inside this triangle
		bool overlap = false;
		for (uint32_t otherVert = 3; otherVert < remainingIndices->length; ++otherVert) {
			uint32_t idx = (guess_vert + otherVert) % remainingIndices->length;
			if (idx >= remainingIndices->length)
				continue; // ???
			uint32_t ovi = ((IndexTriplet *)OBJZ_ARRAY_ELEMENT(*remainingIndices, idx))->v;
			float tx = ((float *)chunkedArrayElement(_positions, ovi))[axes[0]];
			float ty = ((float *)chunkedArrayElement(_positions, ovi))[axes[1]];
			if (pnpoly(3, vx, vy, tx, ty)) {
				overlap = true;
				break;
			}
		}
		if (overlap) {
			guess_vert += 1;
			continue;
		}
		// this triangle is an ear
		Face face;
		for (int i = 0; i < 3; i++)
			face.indices[i] = *ind[i];
		face.materialIndex = (int16_t)_materialIndex;
		face.smoothingGroup = _smoothingGroup;
		chunkedArrayAppend(_faces, &face);
		// remove v1 from the list
		uint32_t removed_vert_index = (guess_vert + 1) % remainingIndices->length;
		while (removed_vert_index + 1 < remainingIndices->length) {
			IndexTriplet *remainingIndicesData = (IndexTriplet *)remainingIndices->data;
			remainingIndicesData[removed_vert_index] = remainingIndicesData[removed_vert_index + 1];
			removed_vert_index += 1;
		}
		remainingIndices->length--;
	}
	if (remainingIndices->length == 3) {
		Face face;
		for (int i = 0; i < 3; i++)
			face.indices[i] = *(IndexTriplet *)OBJZ_ARRAY_ELEMENT(*remainingIndices, i);
		face.materialIndex = (int16_t)_materialIndex;
		face.smoothingGroup = _smoothingGroup;
		chunkedArrayAppend(_faces, &face);
	}
}

static vec3 calculateSmoothNormal(uint32_t _pos, ChunkedArray *_faces, Array *_faceNormals, uint16_t _smoothingGroup) {
	vec3 normal;
	OBJZ_VEC3_SET(normal, 0, 0, 0);
	int n = 0;
	for (uint32_t i = 0; i < _faces->length; i++) {
		const Face *face = chunkedArrayElement(_faces, i);
		if (face->smoothingGroup != _smoothingGroup)
			continue;
		for (int j = 0; j < 3; j++) {
			if (face->indices[j].v == _pos) {
				OBJZ_VEC3_ADD(normal, normal, *(vec3 *)OBJZ_ARRAY_ELEMENT(*_faceNormals, i));
				n++;
				break;
			}
		}
	}
	const float s = 1.0f / n;
	OBJZ_VEC3_MUL(normal, normal, s);
	vec3Normalize(&normal, &normal);
	return normal;
}

void objz_setRealloc(objzReallocFunc _realloc) {
	s_realloc = _realloc;
}

void objz_setIndexFormat(uint32_t _format) {
	s_indexFormat = _format;
}

void objz_setVertexFormat(size_t _stride, size_t _positionOffset, size_t _texcoordOffset, size_t _normalOffset) {
	s_vertexDecl.stride = _stride;
	s_vertexDecl.positionOffset = _positionOffset;
	s_vertexDecl.texcoordOffset = _texcoordOffset;
	s_vertexDecl.normalOffset = _normalOffset;
}

objzModel *objz_load(const char *_filename) {
	s_error[0] = 0;
	File file;
	if (!fileOpen(&file, _filename)) {
		appendError("Failed to read file '%s'", _filename);
		return NULL;
	}
	// Parse the obj file and any material files.
	// Faces are triangulated. Other than that, this is straight parsing.
	Array materialLibs, materials, tempObjects;
	ChunkedArray positions, texcoords, normals, faces;
	Array faceIndices, tempFaceIndices; // Re-used per face.
	arrayInit(&materialLibs, sizeof(char) * OBJZ_MAX_TOKEN_LENGTH, 1);
	arrayInit(&materials, sizeof(objzMaterial), 16);
	arrayInit(&tempObjects, sizeof(TempObject), 64);
	chunkedArrayInit(&positions, sizeof(float) * 3, 100000);
	chunkedArrayInit(&texcoords, sizeof(float) * 2, 100000);
	chunkedArrayInit(&normals, sizeof(float) * 3, 100000);
	chunkedArrayInit(&faces, sizeof(Face), 100000);
	arrayInit(&faceIndices, sizeof(IndexTriplet), 8);
	arrayInit(&tempFaceIndices, sizeof(IndexTriplet), 8);
	bool generateNormals = false;
	char currentGroupName[OBJZ_NAME_MAX] = { 0 };
	char currentObjectName[OBJZ_NAME_MAX] = { 0 };
	int32_t currentMaterialIndex = -1;
	uint16_t currentSmoothingGroup = 0;
	uint32_t flags = 0;
	Lexer lexer;
	initLexer(&lexer);
	Token token;
	for (;;) {
		char *line = fileReadLine(&file);
		if (!line)
			break;
		lexerSetLine(&lexer, line);
		tokenize(&lexer, &token, false);
		if (OBJZ_STRICMP(token.text, "f") == 0) {
			// Get current object.
			if (tempObjects.length == 0) {
				// No objects specifed, but there's a face, so create one.
				TempObject o;
				o.name[0] = 0;
				o.firstFace = o.numFaces = 0;
				arrayAppend(&tempObjects, &o);
			}
			TempObject *object = OBJZ_ARRAY_ELEMENT(tempObjects, tempObjects.length - 1);
			// Parse triplets.
			faceIndices.length = 0;
			for (;;) {
				Token tripletToken;
				tokenize(&lexer, &tripletToken, false);
				if (tripletToken.text[0] == 0) {
					if (isEol(&lexer))
						break;
					appendError("(%u:%u) Failed to parse face", tripletToken.line, tripletToken.column);
					goto error;
				}
				// Parse v/vt/vn triplet.
				int32_t rawTriplet[3];
				if (!parseVertexAttribIndices(&tripletToken, rawTriplet)) {
					appendError("(%u:%u) Failed to parse face", tripletToken.line, tripletToken.column);
					goto error;
				}
				IndexTriplet triplet;
				triplet.v = fixVertexAttribIndex(rawTriplet[0], positions.length);
				triplet.vt = fixVertexAttribIndex(rawTriplet[1], texcoords.length);
				triplet.vn = fixVertexAttribIndex(rawTriplet[2], normals.length);
				arrayAppend(&faceIndices, &triplet);
				if (triplet.vn == UINT32_MAX && s_vertexDecl.normalOffset != SIZE_MAX)
					generateNormals = true;
			}
			if (faceIndices.length < 3) {
				appendError("(%u:%u) Face needs at least 3 vertices", token.line, token.column);
				goto error;
			}
			// Triangulate.
			if (faceIndices.length == 3) {
				Face face;
				face.materialIndex = (int16_t)currentMaterialIndex;
				face.smoothingGroup = currentSmoothingGroup;
				for (int i = 0; i < 3; i++)
					face.indices[i] = *(IndexTriplet *)OBJZ_ARRAY_ELEMENT(faceIndices, i);
				chunkedArrayAppend(&faces, &face);
				object->numFaces++;
			} else {
				const uint32_t prevFacesLength = faces.length;
				triangulate(&faceIndices, &positions, &tempFaceIndices, &faces, currentMaterialIndex, currentSmoothingGroup);
				object->numFaces += faces.length - prevFacesLength;
			}
		} else if (OBJZ_STRICMP(token.text, "g") == 0 || OBJZ_STRICMP(token.text, "o") == 0) {
			const bool isGroup = OBJZ_STRICMP(token.text, "g") == 0;
			tokenize(&lexer, &token, true);
			if (isGroup) {
				// Empty group names are permitted.
				if (token.text[0] != 0)
					strCopy(currentGroupName, sizeof(currentGroupName), token.text, strLength(token.text, sizeof(token.text)));
			}
			else {
				if (token.text[0] == 0) {
					appendError("(%u:%u) Expected name after 'o'", token.line, token.column);
					goto error;
				}
				strCopy(currentObjectName, sizeof(currentObjectName), token.text, strLength(token.text, sizeof(token.text)));
			}
			TempObject o;
			o.name[0] = 0;
			if (currentGroupName[0] != 0)
				strCopy(o.name, sizeof(o.name), currentGroupName, strLength(currentGroupName, sizeof(currentGroupName)));
			if (currentObjectName[0] != 0) {
				if (strLength(o.name, sizeof(o.name)) > 0)
					strConcat(o.name, sizeof(o.name), " ", 1);
				strConcat(o.name, sizeof(o.name), currentObjectName, strLength(currentObjectName, sizeof(currentObjectName)));
			}
			o.firstFace = faces.length;
			o.numFaces = 0;
			arrayAppend(&tempObjects, &o);
		} else if (OBJZ_STRICMP(token.text, "mtllib") == 0) {
			tokenize(&lexer, &token, true);
			if (token.text[0] == 0) {
				appendError("(%u:%u) Expected name after 'mtllib'", token.line, token.column);
				goto error;
			}
			// Don't load the same material library twice.
			bool alreadyLoaded = false;
			for (uint32_t i = 0; i < materialLibs.length; i++) {
				if (OBJZ_STRICMP(token.text, (const char *)OBJZ_ARRAY_ELEMENT(materialLibs, i)) == 0) {
					alreadyLoaded = true;
					break;
				}
			}
			if (!alreadyLoaded) {
				if (!loadMaterialFile(_filename, token.text, &materials))
					goto error;
				arrayAppend(&materialLibs, token.text);
			}
		} else if (OBJZ_STRICMP(token.text, "s") == 0) {
			tokenize(&lexer, &token, false);
			if (token.text[0] == 0) {
				appendError("(%u:%u) Expected value after 's'", token.line, token.column);
				goto error;
			}
			if (OBJZ_STRICMP(token.text, "off") == 0)
				currentSmoothingGroup = 0;
			else
				currentSmoothingGroup = (uint16_t)atoi(token.text);
		} else if (OBJZ_STRICMP(token.text, "usemtl") == 0) {
			tokenize(&lexer, &token, false);
			if (token.text[0] == 0) {
				appendError("(%u:%u) Expected name after 'usemtl'", token.line, token.column);
				goto error;
			}
			currentMaterialIndex = -1;
			for (uint32_t i = 0; i < materials.length; i++) {
				const objzMaterial *mat = OBJZ_ARRAY_ELEMENT(materials, i);
				if (OBJZ_STRICMP(mat->name, token.text) == 0) {
					currentMaterialIndex = (int)i;
					break;
				}
			}
		} else if (OBJZ_STRICMP(token.text, "v") == 0) {
			float pos[3];
			if (!parseFloats(&lexer, pos, 3))
				goto error;
			chunkedArrayAppend(&positions, pos);
		} else if (OBJZ_STRICMP(token.text, "vn") == 0) {
			float normal[3];
			if (!parseFloats(&lexer, normal, 3))
				goto error;
			chunkedArrayAppend(&normals, normal);
			flags |= OBJZ_FLAG_NORMALS;
		} else if (OBJZ_STRICMP(token.text, "vt") == 0) {
			float texcoord[2];
			if (!parseFloats(&lexer, texcoord, 2))
				goto error;
			chunkedArrayAppend(&texcoords, texcoord);
			flags |= OBJZ_FLAG_TEXCOORDS;
		}
	}
	if (normals.length == 0)
		generateNormals = true;
	arrayDestroy(&materialLibs);
	arrayDestroy(&faceIndices);
	arrayDestroy(&tempFaceIndices);
	fileClose(&file);
	// Do some post-processing of parsed data:
	//   * generate normals
	//   * find unique vertices from separately index vertex attributes (pos, texcoord, normal).
	//   * build meshes by batching object faces by material
	Array faceNormals;
	arrayInit(&faceNormals, sizeof(vec3), faces.length); // Exact capacity
	if (generateNormals) {
		for (uint32_t i = 0; i < tempObjects.length; i++) {
			const TempObject *tempObject = OBJZ_ARRAY_ELEMENT(tempObjects, i);
			for (uint32_t j = 0; j < tempObject->numFaces; j++) {
				const Face *face = chunkedArrayElement(&faces, tempObject->firstFace + j);
				vec3 edge0, edge1;
				vec3 normal;
				const vec3 *p0 = chunkedArrayElement(&positions, face->indices[0].v);
				const vec3 *p1 =  chunkedArrayElement(&positions, face->indices[1].v);
				const vec3 *p2 = chunkedArrayElement(&positions, face->indices[2].v);
				OBJZ_VEC3_SUB(edge0, *p1, *p0);
				OBJZ_VEC3_SUB(edge1, *p2, *p0);
				OBJZ_VEC3_CROSS(normal, edge0, edge1);
				vec3Normalize(&normal, &normal);
				arrayAppend(&faceNormals, &normal);
			}
		}
	}
	Array meshes, objects, indices;
	arrayInit(&meshes, sizeof(objzMesh), tempObjects.length * 4); // Guess capacity: 4 meshes per object
	arrayInit(&objects, sizeof(objzObject), tempObjects.length); // Exact capacity
	arrayInit(&indices, sizeof(uint32_t), faces.length * 3); // Exact capacity
	VertexHashMap vertexHashMap;
	vertexHashMapInit(&vertexHashMap, positions.length * 2); // Guess capacity
	NormalHashMap normalHashMap; // Re-used for each object.
	if (generateNormals) {
		uint32_t maxObjectFaces = 0;
		for (uint32_t i = 0; i < tempObjects.length; i++) {
			const TempObject *tempObject = OBJZ_ARRAY_ELEMENT(tempObjects, i);
			maxObjectFaces = OBJZ_LARGEST(maxObjectFaces, tempObject->numFaces);
		}
		normalHashMapInit(&normalHashMap, OBJZ_LARGEST(maxObjectFaces, 32), &normals); // Guess capacity.
	}
	for (uint32_t i = 0; i < tempObjects.length; i++) {
		const TempObject *tempObject = OBJZ_ARRAY_ELEMENT(tempObjects, i);
		if (!tempObject->numFaces)
			continue;
		objzObject object;
		strCopy(object.name, sizeof(object.name), tempObject->name, strLength(tempObject->name, sizeof(tempObject->name)));
		if (generateNormals)
			normalHashMapClear(&normalHashMap);
		// Create one mesh per material. No material (-1) gets a mesh too.
		object.firstMesh = meshes.length;
		object.numMeshes = 0;
		for (int32_t material = -1; material < (int32_t)materials.length; material++) {
			objzMesh mesh;
			mesh.firstIndex = indices.length;
			mesh.numIndices = 0;
			mesh.materialIndex = material;
			for (uint32_t j = 0; j < tempObject->numFaces; j++) {
				const Face *face = chunkedArrayElement(&faces, tempObject->firstFace + j);
				if (face->materialIndex != (int16_t)material)
					continue;
				uint32_t faceNormalIndex = UINT32_MAX;
				if (generateNormals && face->smoothingGroup == 0) {
					for (int k = 0; k < 3; k++) {
						if (face->indices[k].vn >= normals.length) {
							faceNormalIndex = normalHashMapInsert(&normalHashMap, OBJZ_ARRAY_ELEMENT(faceNormals, tempObject->firstFace + j));
							break;
						}
					}
				}
				for (int k = 0; k < 3; k++) {
					const IndexTriplet *triplet = &face->indices[k];
					uint32_t vn = triplet->vn;
					if (generateNormals) {
						if (face->smoothingGroup > 0) {
							vec3 normal = calculateSmoothNormal(triplet->v, &faces, &faceNormals, face->smoothingGroup);
							vn = normalHashMapInsert(&normalHashMap, &normal);
						} else if (faceNormalIndex != UINT32_MAX)
							vn = faceNormalIndex;
					}
					const uint32_t index = vertexHashMapInsert(&vertexHashMap, i, triplet->v, triplet->vt, vn);
					if (index > UINT16_MAX)
						flags |= OBJZ_FLAG_INDEX32;
					arrayAppend(&indices, &index);
					mesh.numIndices++;
				}
			}
			if (mesh.numIndices > 0) {
				arrayAppend(&meshes, &mesh);
				object.numMeshes++;
			}
		}
		if (objects.length > 0) {
			const objzObject *prev = OBJZ_ARRAY_ELEMENT(objects, objects.length - 1);
			object.firstIndex = prev->firstIndex + prev->numIndices;
			object.firstVertex = prev->firstVertex + prev->numVertices;
		} else {
			object.firstIndex = 0;
			object.firstVertex = 0;
		}
		object.numIndices = indices.length - object.firstIndex;
		object.numVertices = vertexHashMap.vertices.length - object.firstVertex;
		arrayAppend(&objects, &object);
	}
	if (generateNormals)
		normalHashMapDestroy(&normalHashMap);
	arrayDestroy(&tempObjects);
	chunkedArrayDestroy(&faces);
	arrayDestroy(&faceNormals);
	// Build output data structure.
	objzModel *model = OBJZ_MALLOC(sizeof(objzModel));
	model->flags = flags;
	if (s_indexFormat == OBJZ_INDEX_FORMAT_U32 || (flags & OBJZ_FLAG_INDEX32))
		model->indices = indices.data;
	else {
		flags &= ~OBJZ_FLAG_INDEX32;
		model->indices = OBJZ_MALLOC(sizeof(uint16_t) * indices.length);
		for (uint32_t i = 0; i < indices.length; i++) {
			uint32_t *index = (uint32_t *)OBJZ_ARRAY_ELEMENT(indices, i);
			((uint16_t *)model->indices)[i] = (uint16_t)*index;
		}
		arrayDestroy(&indices);
	}
	model->numIndices = indices.length;
	model->materials = (objzMaterial *)materials.data;
	model->numMaterials = materials.length;
	model->meshes = (objzMesh *)meshes.data;
	model->numMeshes = meshes.length;
	model->objects = (objzObject *)objects.data;
	model->numObjects = objects.length;
	model->vertices = OBJZ_MALLOC(s_vertexDecl.stride * vertexHashMap.vertices.length);
	for (uint32_t i = 0; i < vertexHashMap.vertices.length; i++) {
		uint8_t *vOut = &((uint8_t *)model->vertices)[i * s_vertexDecl.stride];
		const HashedVertex *vIn = OBJZ_ARRAY_ELEMENT(vertexHashMap.vertices, i);
		if (s_vertexDecl.positionOffset != SIZE_MAX)
			memcpy(&vOut[s_vertexDecl.positionOffset], chunkedArrayElement(&positions, vIn->pos), sizeof(float) * 3);
		if (s_vertexDecl.texcoordOffset != SIZE_MAX) {
			if (vIn->texcoord == UINT32_MAX)
				memset(&vOut[s_vertexDecl.texcoordOffset], 0, sizeof(float) * 2);
			else
				memcpy(&vOut[s_vertexDecl.texcoordOffset], chunkedArrayElement(&texcoords, vIn->texcoord), sizeof(float) * 2);
		}
		if (s_vertexDecl.normalOffset != SIZE_MAX) {
			if (vIn->normal == UINT32_MAX)
				memset(&vOut[s_vertexDecl.normalOffset], 0, sizeof(float) * 3);
			else
			memcpy(&vOut[s_vertexDecl.normalOffset], chunkedArrayElement(&normals, vIn->normal), sizeof(float) * 3);
		}
	}
	model->numVertices = vertexHashMap.vertices.length;
	chunkedArrayDestroy(&positions);
	chunkedArrayDestroy(&texcoords);
	chunkedArrayDestroy(&normals);
	vertexHashMapDestroy(&vertexHashMap);
	return model;
error:
	fileClose(&file);
	arrayDestroy(&materialLibs);
	arrayDestroy(&materials);
	arrayDestroy(&tempObjects);
	chunkedArrayDestroy(&positions);
	chunkedArrayDestroy(&texcoords);
	chunkedArrayDestroy(&normals);
	chunkedArrayDestroy(&faces);
	arrayDestroy(&faceIndices);
	arrayDestroy(&tempFaceIndices);
	return NULL;
}

void objz_destroy(objzModel *_model) {
	if (!_model)
		return;
	OBJZ_FREE(_model->indices);
	OBJZ_FREE(_model->materials);
	OBJZ_FREE(_model->meshes);
	OBJZ_FREE(_model->objects);
	OBJZ_FREE(_model->vertices);
	OBJZ_FREE(_model);
}

const char *objz_getError() {
	if (s_error[0])
		return s_error;
	return NULL;
}
