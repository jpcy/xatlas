// Copyright NVIDIA Corporation 2008 -- Ignacio Castano <icastano@nvidia.com>

#include "ParameterizationQuality.h"
#include "nvmesh/halfedge/Mesh.h"
#include "nvmesh/halfedge/Face.h"
#include "xatlas.h"
#include <float.h>

using namespace nv;

ParameterizationQuality::ParameterizationQuality()
{
	m_totalTriangleCount = 0;
	m_flippedTriangleCount = 0;
	m_zeroAreaTriangleCount = 0;
	m_parametricArea = 0.0f;
	m_geometricArea = 0.0f;
	m_stretchMetric = 0.0f;
	m_maxStretchMetric = 0.0f;
	m_conformalMetric = 0.0f;
	m_authalicMetric = 0.0f;
}

ParameterizationQuality::ParameterizationQuality(const HalfEdge::Mesh *mesh)
{
	nvDebugCheck(mesh != NULL);
	m_totalTriangleCount = 0;
	m_flippedTriangleCount = 0;
	m_zeroAreaTriangleCount = 0;
	m_parametricArea = 0.0f;
	m_geometricArea = 0.0f;
	m_stretchMetric = 0.0f;
	m_maxStretchMetric = 0.0f;
	m_conformalMetric = 0.0f;
	m_authalicMetric = 0.0f;
	const uint32_t faceCount = mesh->faceCount();
	for (uint32_t f = 0; f < faceCount; f++) {
		const HalfEdge::Face *face = mesh->faceAt(f);
		const HalfEdge::Vertex *vertex0 = NULL;
		Vector3 p[3];
		Vector2 t[3];
		for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance()) {
			const HalfEdge::Edge *edge = it.current();
			if (vertex0 == NULL) {
				vertex0 = edge->vertex;
				p[0] = vertex0->pos;
				t[0] = vertex0->tex;
			} else if (edge->to() != vertex0) {
				p[1] = edge->from()->pos;
				p[2] = edge->to()->pos;
				t[1] = edge->from()->tex;
				t[2] = edge->to()->tex;
				processTriangle(p, t);
			}
		}
	}
	if (m_flippedTriangleCount + m_zeroAreaTriangleCount == faceCount) {
		// If all triangles are flipped, then none is.
		m_flippedTriangleCount = 0;
	}
	nvDebugCheck(std::isfinite(m_parametricArea) && m_parametricArea >= 0);
	nvDebugCheck(std::isfinite(m_geometricArea) && m_geometricArea >= 0);
	nvDebugCheck(std::isfinite(m_stretchMetric));
	nvDebugCheck(std::isfinite(m_maxStretchMetric));
	nvDebugCheck(std::isfinite(m_conformalMetric));
	nvDebugCheck(std::isfinite(m_authalicMetric));
}

bool ParameterizationQuality::isValid() const
{
	return m_flippedTriangleCount == 0; // @@ Does not test for self-overlaps.
}

float ParameterizationQuality::rmsStretchMetric() const
{
	if (m_geometricArea == 0) return 0.0f;
	float normFactor = sqrtf(m_parametricArea / m_geometricArea);
	return sqrtf(m_stretchMetric / m_geometricArea) * normFactor;
}

float ParameterizationQuality::maxStretchMetric() const
{
	if (m_geometricArea == 0) return 0.0f;
	float normFactor = sqrtf(m_parametricArea / m_geometricArea);
	return m_maxStretchMetric * normFactor;
}

float ParameterizationQuality::rmsConformalMetric() const
{
	if (m_geometricArea == 0) return 0.0f;
	return sqrtf(m_conformalMetric / m_geometricArea);
}

float ParameterizationQuality::maxAuthalicMetric() const
{
	if (m_geometricArea == 0) return 0.0f;
	return sqrtf(m_authalicMetric / m_geometricArea);
}

void ParameterizationQuality::operator += (const ParameterizationQuality &pq)
{
	m_totalTriangleCount += pq.m_totalTriangleCount;
	m_flippedTriangleCount += pq.m_flippedTriangleCount;
	m_zeroAreaTriangleCount += pq.m_zeroAreaTriangleCount;
	m_parametricArea += pq.m_parametricArea;
	m_geometricArea += pq.m_geometricArea;
	m_stretchMetric += pq.m_stretchMetric;
	m_maxStretchMetric = std::max(m_maxStretchMetric, pq.m_maxStretchMetric);
	m_conformalMetric += pq.m_conformalMetric;
	m_authalicMetric += pq.m_authalicMetric;
}


void ParameterizationQuality::processTriangle(Vector3 q[3], Vector2 p[3])
{
	m_totalTriangleCount++;
	// Evaluate texture stretch metric. See:
	// - "Texture Mapping Progressive Meshes", Sander, Snyder, Gortler & Hoppe
	// - "Mesh Parameterization: Theory and Practice", Siggraph'07 Course Notes, Hormann, Levy & Sheffer.
	float t1 = p[0].x;
	float s1 = p[0].y;
	float t2 = p[1].x;
	float s2 = p[1].y;
	float t3 = p[2].x;
	float s3 = p[2].y;
	float geometricArea = length(cross(q[1] - q[0], q[2] - q[0])) / 2;
	float parametricArea = ((s2 - s1) * (t3 - t1) - (s3 - s1) * (t2 - t1)) / 2;
	if (isZero(parametricArea)) {
		m_zeroAreaTriangleCount++;
		return;
	}
	Vector3 Ss = (q[0] * (t2 - t3) + q[1] * (t3 - t1) + q[2] * (t1 - t2)) / (2 * parametricArea);
	Vector3 St = (q[0] * (s3 - s2) + q[1] * (s1 - s3) + q[2] * (s2 - s1)) / (2 * parametricArea);
	float a = dot(Ss, Ss); // E
	float b = dot(Ss, St); // F
	float c = dot(St, St); // G
	// Compute eigen-values of the first fundamental form:
	float sigma1 = sqrtf(0.5f * std::max(0.0f, a + c - sqrtf(square(a - c) + 4 * square(b)))); // gamma uppercase, min eigenvalue.
	float sigma2 = sqrtf(0.5f * std::max(0.0f, a + c + sqrtf(square(a - c) + 4 * square(b)))); // gamma lowercase, max eigenvalue.
	nvCheck(sigma2 >= sigma1);
	// isometric: sigma1 = sigma2 = 1
	// conformal: sigma1 / sigma2 = 1
	// authalic: sigma1 * sigma2 = 1
	float rmsStretch = sqrtf((a + c) * 0.5f);
	float rmsStretch2 = sqrtf((square(sigma1) + square(sigma2)) * 0.5f);
	nvDebugCheck(equal(rmsStretch, rmsStretch2, 0.01f));
	if (parametricArea < 0.0f) {
		// Count flipped triangles.
		m_flippedTriangleCount++;
		parametricArea = fabsf(parametricArea);
	}
	m_stretchMetric += square(rmsStretch) * geometricArea;
	m_maxStretchMetric = std::max(m_maxStretchMetric, sigma2);
	if (!isZero(sigma1, 0.000001f)) {
		// sigma1 is zero when geometricArea is zero.
		m_conformalMetric += (sigma2 / sigma1) * geometricArea;
	}
	m_authalicMetric += (sigma1 * sigma2) * geometricArea;
	// Accumulate total areas.
	m_geometricArea += geometricArea;
	m_parametricArea += parametricArea;
	//triangleConformalEnergy(q, p);
}
