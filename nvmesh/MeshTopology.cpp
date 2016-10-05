// This code is in the public domain -- castanyo@yahoo.es
#include <vector>
#include "xatlas.h"

#include "nvmesh/MeshTopology.h"

using namespace nv;

void MeshTopology::buildTopologyInfo(const HalfEdge::Mesh *mesh)
{
	const uint32_t vertexCount = mesh->colocalVertexCount();
	const uint32_t faceCount = mesh->faceCount();
	const uint32_t edgeCount = mesh->edgeCount();
	nvDebug( "--- Building mesh topology:\n" );
	std::vector<uint32_t> stack(faceCount);
	BitArray bitFlags(faceCount);
	bitFlags.clearAll();
	// Compute connectivity.
	nvDebug( "---   Computing connectivity.\n" );
	m_connectedCount = 0;
	for (uint32_t f = 0; f < faceCount; f++ ) {
		if ( bitFlags.bitAt(f) == false ) {
			m_connectedCount++;
			stack.push_back( f );
			while ( !stack.empty() ) {
				const uint32_t top = stack.back();
				nvCheck(top != uint32_t(~0));
				stack.pop_back();
				if ( bitFlags.bitAt(top) == false ) {
					bitFlags.setBitAt(top);
					const HalfEdge::Face *face = mesh->faceAt(top);
					const HalfEdge::Edge *firstEdge = face->edge;
					const HalfEdge::Edge *edge = firstEdge;
					do {
						const HalfEdge::Face *neighborFace = edge->pair->face;
						if (neighborFace != NULL) {
							stack.push_back(neighborFace->id);
						}
						edge = edge->next;
					} while (edge != firstEdge);
				}
			}
		}
	}
	nvCheck(stack.empty());
	nvDebug( "---   %d connected components.\n", m_connectedCount );
	// Count boundary loops.
	nvDebug( "---   Counting boundary loops.\n" );
	m_boundaryCount = 0;
	bitFlags.resize(edgeCount);
	bitFlags.clearAll();
	// Don't forget to link the boundary otherwise this won't work.
	for (uint32_t e = 0; e < edgeCount; e++) {
		const HalfEdge::Edge *startEdge = mesh->edgeAt(e);
		if (startEdge != NULL && startEdge->isBoundary() && bitFlags.bitAt(e) == false) {
			nvDebugCheck(startEdge->face != NULL);
			nvDebugCheck(startEdge->pair->face == NULL);
			startEdge = startEdge->pair;
			m_boundaryCount++;
			const HalfEdge::Edge *edge = startEdge;
			do {
				bitFlags.setBitAt(edge->id / 2);
				edge = edge->next;
			} while (startEdge != edge);
		}
	}
	nvDebug("---   %d boundary loops found.\n", m_boundaryCount );
	// Compute euler number.
	m_eulerNumber = vertexCount - edgeCount + faceCount;
	nvDebug("---   Euler number: %d.\n", m_eulerNumber);
	// Compute genus. (only valid on closed connected surfaces)
	m_genus = -1;
	if ( isClosed() && isConnected() ) {
		m_genus = (2 - m_eulerNumber) / 2;
		nvDebug("---   Genus: %d.\n", m_genus);
	}
}
