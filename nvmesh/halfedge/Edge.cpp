// This code is in the public domain -- castanyo@yahoo.es

#include "xatlas.h"

namespace nv {
namespace HalfEdge {

bool Edge::isNormalSeam() const
{
	return (vertex->nor != pair->next->vertex->nor || next->vertex->nor != pair->vertex->nor);
}

bool Edge::isTextureSeam() const
{
	return (vertex->tex != pair->next->vertex->tex || next->vertex->tex != pair->vertex->tex);
}

Vector3 Edge::midPoint() const
{
	return (to()->pos + from()->pos) * 0.5f;
}

float Edge::length() const
{
	return nv::length(to()->pos - from()->pos);
}

float Edge::angle() const
{
	Vector3 p = vertex->pos;
	Vector3 a = prev->vertex->pos;
	Vector3 b = next->vertex->pos;
	Vector3 v0 = a - p;
	Vector3 v1 = b - p;
	return acosf(dot(v0, v1) / (nv::length(v0) * nv::length(v1)));
}

}
}
