// Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>

#include "Atlas.h"
#include "xatlas.h"

namespace nv {
namespace param {


/// Ctor.
Atlas::Atlas()
{
}

// Dtor.
Atlas::~Atlas()
{
	for (size_t i = 0; i < m_meshChartsArray.size(); i++)
		delete m_meshChartsArray[i];
}

uint32_t Atlas::chartCount() const
{
	uint32_t count = 0;
	for (uint32_t c = 0; c < m_meshChartsArray.size(); c++) {
		count += m_meshChartsArray[c]->chartCount();
	}
	return count;
}

const Chart *Atlas::chartAt(uint32_t i) const
{
	for (uint32_t c = 0; c < m_meshChartsArray.size(); c++) {
		uint32_t count = m_meshChartsArray[c]->chartCount();
		if (i < count) {
			return m_meshChartsArray[c]->chartAt(i);
		}
		i -= count;
	}
	return NULL;
}

Chart *Atlas::chartAt(uint32_t i)
{
	for (uint32_t c = 0; c < m_meshChartsArray.size(); c++) {
		uint32_t count = m_meshChartsArray[c]->chartCount();
		if (i < count) {
			return m_meshChartsArray[c]->chartAt(i);
		}
		i -= count;
	}
	return NULL;
}

// Extract the charts and add to this atlas.
void Atlas::addMeshCharts(MeshCharts *meshCharts)
{
	m_meshChartsArray.push_back(meshCharts);
}

void Atlas::extractCharts(const HalfEdge::Mesh *mesh)
{
	MeshCharts *meshCharts = new MeshCharts(mesh);
	meshCharts->extractCharts();
	addMeshCharts(meshCharts);
}

void Atlas::computeCharts(const HalfEdge::Mesh *mesh, const SegmentationSettings &settings, const std::vector<uint32_t> &unchartedMaterialArray)
{
	MeshCharts *meshCharts = new MeshCharts(mesh);
	meshCharts->computeCharts(settings, unchartedMaterialArray);
	addMeshCharts(meshCharts);
}

void Atlas::parameterizeCharts()
{
	for (uint32_t i = 0; i < m_meshChartsArray.size(); i++) {
		m_meshChartsArray[i]->parameterizeCharts();
	}
}


float Atlas::packCharts(int quality, float texelsPerUnit, bool blockAlign, bool conservative)
{
	AtlasPacker packer(this);
	packer.packCharts(quality, texelsPerUnit, blockAlign, conservative);
	return packer.computeAtlasUtilization();
}




/// Ctor.
MeshCharts::MeshCharts(const HalfEdge::Mesh *mesh) : m_mesh(mesh)
{
}

// Dtor.
MeshCharts::~MeshCharts()
{
	for (size_t i = 0; i < m_chartArray.size(); i++)
		delete m_chartArray[i];
}

void MeshCharts::extractCharts()
{
	const uint32_t faceCount = m_mesh->faceCount();
	int first = 0;
	std::vector<uint32_t> queue;
	queue.reserve(faceCount);
	BitArray bitFlags(faceCount);
	bitFlags.clearAll();
	for (uint32_t f = 0; f < faceCount; f++) {
		if (bitFlags.bitAt(f) == false) {
			// Start new patch. Reset queue.
			first = 0;
			queue.clear();
			queue.push_back(f);
			bitFlags.setBitAt(f);
			while (first != queue.size()) {
				const HalfEdge::Face *face = m_mesh->faceAt(queue[first]);
				// Visit face neighbors of queue[first]
				for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance()) {
					const HalfEdge::Edge *edge = it.current();
					nvDebugCheck(edge->pair != NULL);
					if (!edge->isBoundary() && /*!edge->isSeam()*/
					        //!(edge->from()->tex() != edge->pair()->to()->tex() || edge->to()->tex() != edge->pair()->from()->tex()))
					        !(edge->from() != edge->pair->to() || edge->to() != edge->pair->from())) { // Preserve existing seams (not just texture seams).
						const HalfEdge::Face *neighborFace = edge->pair->face;
						nvDebugCheck(neighborFace != NULL);
						if (bitFlags.bitAt(neighborFace->id) == false) {
							queue.push_back(neighborFace->id);
							bitFlags.setBitAt(neighborFace->id);
						}
					}
				}
				first++;
			}
			Chart *chart = new Chart();
			chart->build(m_mesh, queue);
			m_chartArray.push_back(chart);
		}
	}
}


/*
LSCM:
- identify sharp features using local dihedral angles.
- identify seed faces farthest from sharp features.
- grow charts from these seeds.

MCGIM:
- phase 1: chart growth
  - grow all charts simultaneously using dijkstra search on the dual graph of the mesh.
  - graph edges are weighted based on planarity metric.
  - metric uses distance to global chart normal.
  - terminate when all faces have been assigned.
- phase 2: seed computation:
  - place new seed of the chart at the most interior face.
  - most interior is evaluated using distance metric only.

- method repeates the two phases, until the location of the seeds does not change.
  - cycles are detected by recording all the previous seeds and chartification terminates.

D-Charts:

- Uniaxial conic metric:
  - N_c = axis of the generalized cone that best fits the chart. (cone can a be cylinder or a plane).
  - omega_c = angle between the face normals and the axis.
  - Fitting error between chart C and tringle t: F(c,t) = (N_c*n_t - cos(omega_c))^2

- Compactness metrics:
  - Roundness:
    - C(c,t) = pi * D(S_c,t)^2 / A_c
    - S_c = chart seed.
    - D(S_c,t) = length of the shortest path inside the chart betwen S_c and t.
    - A_c = chart area.
  - Straightness:
    - P(c,t) = l_out(c,t) / l_in(c,t)
    - l_out(c,t) = lenght of the edges not shared between C and t.
    - l_in(c,t) = lenght of the edges shared between C and t.

- Combined metric:
  - Cost(c,t) = F(c,t)^alpha + C(c,t)^beta + P(c,t)^gamma
  - alpha = 1, beta = 0.7, gamma = 0.5




Our basic approach:
- Just one iteration of k-means?
- Avoid dijkstra by greedily growing charts until a threshold is met. Increase threshold and repeat until no faces left.
- If distortion metric is too high, split chart, add two seeds.
- If chart size is low, try removing chart.


Postprocess:
- If topology is not disk:
  - Fill holes, if new faces fit proxy.
  - Find best cut, otherwise.
- After parameterization:
  - If boundary self-intersects:
    - cut chart along the closest two diametral boundary vertices, repeat parametrization.
    - what if the overlap is on an appendix? How do we find that out and cut appropiately?
      - emphasize roundness metrics to prevent those cases.
  - If interior self-overlaps: preserve boundary parameterization and use mean-value map.

*/

void MeshCharts::computeCharts(const SegmentationSettings &settings, const std::vector<uint32_t> &unchartedMaterialArray)
{
	Chart *vertexMap = NULL;
	if (unchartedMaterialArray.size() != 0) {
		vertexMap = new Chart();
		vertexMap->buildVertexMap(m_mesh, unchartedMaterialArray);
		if (vertexMap->faceCount() == 0) {
			delete vertexMap;
			vertexMap = NULL;
		}
	}
	AtlasBuilder builder(m_mesh);
	if (vertexMap != NULL) {
		// Mark faces that do not need to be charted.
		builder.markUnchartedFaces(vertexMap->faceArray());
		m_chartArray.push_back(vertexMap);
	}
	if (builder.facesLeft != 0) {
		// Tweak these values:
		const float maxThreshold = 2;
		const uint32_t growFaceCount = 32;
		const uint32_t maxIterations = 4;
		builder.settings = settings;
		//builder.settings.proxyFitMetricWeight *= 0.75; // relax proxy fit weight during initial seed placement.
		//builder.settings.roundnessMetricWeight = 0;
		//builder.settings.straightnessMetricWeight = 0;
		// This seems a reasonable estimate.
		uint32_t maxSeedCount = std::max(6U, builder.facesLeft);
		// Create initial charts greedely.
		nvDebug("### Placing seeds\n");
		builder.placeSeeds(maxThreshold, maxSeedCount);
		nvDebug("###   Placed %d seeds (max = %d)\n", builder.chartCount(), maxSeedCount);
		builder.updateProxies();
		builder.mergeCharts();
#if 1
		nvDebug("### Relocating seeds\n");
		builder.relocateSeeds();
		nvDebug("### Reset charts\n");
		builder.resetCharts();
		if (vertexMap != NULL) {
			builder.markUnchartedFaces(vertexMap->faceArray());
		}
		builder.settings = settings;
		nvDebug("### Growing charts\n");
		// Restart process growing charts in parallel.
		uint32_t iteration = 0;
		while (true) {
			if (!builder.growCharts(maxThreshold, growFaceCount)) {
				nvDebug("### Can't grow anymore\n");
				// If charts cannot grow more: fill holes, merge charts, relocate seeds and start new iteration.
				nvDebug("### Filling holes\n");
				builder.fillHoles(maxThreshold);
				nvDebug("###   Using %d charts now\n", builder.chartCount());
				builder.updateProxies();
				nvDebug("### Merging charts\n");
				builder.mergeCharts();
				nvDebug("###   Using %d charts now\n", builder.chartCount());
				nvDebug("### Reseeding\n");
				if (!builder.relocateSeeds()) {
					nvDebug("### Cannot relocate seeds anymore\n");
					// Done!
					break;
				}
				if (iteration == maxIterations) {
					nvDebug("### Reached iteration limit\n");
					break;
				}
				iteration++;
				nvDebug("### Reset charts\n");
				builder.resetCharts();
				if (vertexMap != NULL) {
					builder.markUnchartedFaces(vertexMap->faceArray());
				}
				nvDebug("### Growing charts\n");
			}
		};
#endif
		// Make sure no holes are left!
		nvDebugCheck(builder.facesLeft == 0);
		const uint32_t chartCount = builder.chartArray.size();
		for (uint32_t i = 0; i < chartCount; i++) {
			Chart *chart = new Chart();
			m_chartArray.push_back(chart);
			chart->build(m_mesh, builder.chartFaces(i));
		}
	}
	const uint32_t chartCount = m_chartArray.size();
	// Build face indices.
	m_faceChart.resize(m_mesh->faceCount());
	m_faceIndex.resize(m_mesh->faceCount());
	for (uint32_t i = 0; i < chartCount; i++) {
		const Chart *chart = m_chartArray[i];
		const uint32_t faceCount = chart->faceCount();
		for (uint32_t f = 0; f < faceCount; f++) {
			uint32_t idx = chart->faceAt(f);
			m_faceChart[idx] = i;
			m_faceIndex[idx] = f;
		}
	}
	// Build an exclusive prefix sum of the chart vertex counts.
	m_chartVertexCountPrefixSum.resize(chartCount);
	if (chartCount > 0) {
		m_chartVertexCountPrefixSum[0] = 0;
		for (uint32_t i = 1; i < chartCount; i++) {
			const Chart *chart = m_chartArray[i - 1];
			m_chartVertexCountPrefixSum[i] = m_chartVertexCountPrefixSum[i - 1] + chart->vertexCount();
		}
		m_totalVertexCount = m_chartVertexCountPrefixSum[chartCount - 1] + m_chartArray[chartCount - 1]->vertexCount();
	} else {
		m_totalVertexCount = 0;
	}
}


void MeshCharts::parameterizeCharts()
{
	ParameterizationQuality globalParameterizationQuality;
	// Parameterize the charts.
	uint32_t diskCount = 0;
	const uint32_t chartCount = m_chartArray.size();
	for (uint32_t i = 0; i < chartCount; i++)
	{
		Chart *chart = m_chartArray[i];

		bool isValid = false;

		if (chart->isVertexMapped())
		{
			continue;
		}

		if (chart->isDisk())
		{
			diskCount++;
			ParameterizationQuality chartParameterizationQuality;
			if (chart->faceCount() == 1) {
				computeSingleFaceMap(chart->unifiedMesh());
				chartParameterizationQuality = ParameterizationQuality(chart->unifiedMesh());
			} else {
				computeOrthogonalProjectionMap(chart->unifiedMesh());
				ParameterizationQuality orthogonalQuality(chart->unifiedMesh());
				computeLeastSquaresConformalMap(chart->unifiedMesh());
				ParameterizationQuality lscmQuality(chart->unifiedMesh());
				chartParameterizationQuality = lscmQuality;
			}
			isValid = chartParameterizationQuality.isValid();
			if (!isValid) {
				nvDebug("*** Invalid parameterization.\n");
#if 0
				// Dump mesh to inspect problem:
				static int pieceCount = 0;
				StringBuilder fileName;
				fileName.format("invalid_chart_%d.obj", pieceCount++);
				exportMesh(chart->unifiedMesh(), fileName.str());
#endif
			}
			// @@ Check that parameterization quality is above a certain threshold.
			// @@ Detect boundary self-intersections.
			globalParameterizationQuality += chartParameterizationQuality;
		}

		// Transfer parameterization from unified mesh to chart mesh.
		chart->transferParameterization();

	}
	nvDebug("  Parameterized %d/%d charts.\n", diskCount, chartCount);
	nvDebug("  RMS stretch metric: %f\n", globalParameterizationQuality.rmsStretchMetric());
	nvDebug("  MAX stretch metric: %f\n", globalParameterizationQuality.maxStretchMetric());
	nvDebug("  RMS conformal metric: %f\n", globalParameterizationQuality.rmsConformalMetric());
	nvDebug("  RMS authalic metric: %f\n", globalParameterizationQuality.maxAuthalicMetric());
}

}
}
