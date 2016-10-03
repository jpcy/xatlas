// This code is in the public domain -- castano@gmail.com

#include "nvcore/nvcore.h"

namespace nv
{

namespace HalfEdge
{
class Mesh;
class Vertex;
}

uint32_t countMeshTriangles(const HalfEdge::Mesh *mesh);

HalfEdge::Mesh *unifyVertices(const HalfEdge::Mesh *inputMesh);
HalfEdge::Mesh *triangulate(const HalfEdge::Mesh *inputMesh);

} // nv namespace
