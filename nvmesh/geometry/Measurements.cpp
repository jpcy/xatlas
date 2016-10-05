// This code is in the public domain -- castano@gmail.com

#include "xatlas.h"

namespace nv {
namespace HalfEdge {

float computeSurfaceArea(const HalfEdge::Mesh *mesh)
{
	float area = 0;
	for (HalfEdge::Mesh::ConstFaceIterator it(mesh->faces()); !it.isDone(); it.advance()) {
		const HalfEdge::Face *face = it.current();
		area += face->area();
	}
	nvDebugCheck(area >= 0);
	return area;
}

float computeParametricArea(const HalfEdge::Mesh *mesh)
{
	float area = 0;
	for (HalfEdge::Mesh::ConstFaceIterator it(mesh->faces()); !it.isDone(); it.advance()) {
		const HalfEdge::Face *face = it.current();
		area += face->parametricArea();
	}
	return area;
}

}
}
