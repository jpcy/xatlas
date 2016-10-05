// This code is in the public domain -- castano@gmail.com

#include <vector>
#include "OrthogonalProjectionMap.h"
#include "xatlas.h"

using namespace nv;

bool nv::computeOrthogonalProjectionMap(HalfEdge::Mesh *mesh)
{
	Vector3 axis[2];
	uint32_t vertexCount = mesh->vertexCount();
	std::vector<Vector3> points(vertexCount);
	points.resize(vertexCount);
	for (uint32_t i = 0; i < vertexCount; i++) {
		points[i] = mesh->vertexAt(i)->pos;
	}
	// Avoid redundant computations.
	float matrix[6];
	fit::computeCovariance(vertexCount, points.data(), matrix);
	if (matrix[0] == 0 && matrix[3] == 0 && matrix[5] == 0) {
		return false;
	}
	float eigenValues[3];
	Vector3 eigenVectors[3];
	if (!nv::fit::eigenSolveSymmetric3(matrix, eigenValues, eigenVectors)) {
		return false;
	}
	axis[0] = normalize(eigenVectors[0]);
	axis[1] = normalize(eigenVectors[1]);
	// Project vertices to plane.
	for (HalfEdge::Mesh::VertexIterator it(mesh->vertices()); !it.isDone(); it.advance()) {
		HalfEdge::Vertex *vertex = it.current();
		vertex->tex.x = dot(axis[0], vertex->pos);
		vertex->tex.y = dot(axis[1], vertex->pos);
	}
	return true;
}
