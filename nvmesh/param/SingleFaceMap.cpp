// Copyright NVIDIA Corporation 2008 -- Ignacio Castano <icastano@nvidia.com>

#include "xatlas.h"

namespace nv {
namespace param {

void computeSingleFaceMap(HalfEdge::Mesh *mesh)
{
	nvDebugCheck(mesh != NULL);
	nvDebugCheck(mesh->faceCount() == 1);
	HalfEdge::Face *face = mesh->faceAt(0);
	nvCheck(face != NULL);
	Vector3 p0 = face->edge->from()->pos;
	Vector3 p1 = face->edge->to()->pos;
	Vector3 X = normalizeSafe(p1 - p0, Vector3(0.0f), 0.0f);
	Vector3 Z = face->normal();
	Vector3 Y = normalizeSafe(cross(Z, X), Vector3(0.0f), 0.0f);
	uint32_t i = 0;
	for (HalfEdge::Face::EdgeIterator it(face->edges()); !it.isDone(); it.advance(), i++) {
		HalfEdge::Vertex *vertex = it.vertex();
		nvCheck(vertex != NULL);
		if (i == 0) {
			vertex->tex = Vector2(0);
		} else {
			Vector3 pn = vertex->pos;
			float xn = dot((pn - p0), X);
			float yn = dot((pn - p0), Y);
			vertex->tex = Vector2(xn, yn);
		}
	}
}

}
}
