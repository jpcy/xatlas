// Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>

#include "Atlas.h"
#include "Util.h"
#include "AtlasBuilder.h"
#include "AtlasPacker.h"
#include "SingleFaceMap.h"
#include "OrthogonalProjectionMap.h"
#include "LeastSquaresConformalMap.h"
#include "ParameterizationQuality.h"

#include "nvmesh/halfedge/Mesh.h"
#include "nvmesh/halfedge/Face.h"
#include "nvmesh/halfedge/Vertex.h"

#include "nvmesh/MeshTopology.h"
#include "nvmesh/param/Util.h"
#include "nvmesh/geometry/Measurements.h"

#include "nvmath/Vector.h"
#include "nvmath/Box.h"
#include "nvmath/Fitting.h"
#include "nvmath/ProximityGrid.h"
#include "nvmath/Morton.h"

using namespace nv;


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


SegmentationSettings::SegmentationSettings()
{
	// Charts have no area or boundary limits right now.
	maxChartArea = NV_FLOAT_MAX;
	maxBoundaryLength = NV_FLOAT_MAX;
	proxyFitMetricWeight = 1.0f;
	roundnessMetricWeight = 0.1f;
	straightnessMetricWeight = 0.25f;
	normalSeamMetricWeight = 1.0f;
	textureSeamMetricWeight = 0.1f;
}



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
		uint32_t maxSeedCount = max(6U, builder.facesLeft);
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



Chart::Chart() : m_isDisk(false), m_isVertexMapped(false)
{
}

void Chart::build(const HalfEdge::Mesh *originalMesh, const std::vector<uint32_t> &faceArray)
{
	// Copy face indices.
	m_faceArray = faceArray;
	const uint32_t meshVertexCount = originalMesh->vertexCount();
	m_chartMesh.reset(new HalfEdge::Mesh());
	m_unifiedMesh.reset(new HalfEdge::Mesh());
	std::vector<uint32_t> chartMeshIndices(meshVertexCount, ~0);
	std::vector<uint32_t> unifiedMeshIndices(meshVertexCount, ~0);
	// Add vertices.
	const uint32_t faceCount = faceArray.size();
	for (uint32_t f = 0; f < faceCount; f++) {
		const HalfEdge::Face *face = originalMesh->faceAt(faceArray[f]);
		nvDebugCheck(face != NULL);
		for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance()) {
			const HalfEdge::Vertex *vertex = it.current()->vertex;
			const HalfEdge::Vertex *unifiedVertex = vertex->firstColocal();
			if (unifiedMeshIndices[unifiedVertex->id] == ~0) {
				unifiedMeshIndices[unifiedVertex->id] = m_unifiedMesh->vertexCount();
				nvDebugCheck(vertex->pos == unifiedVertex->pos);
				m_unifiedMesh->addVertex(vertex->pos);
			}
			if (chartMeshIndices[vertex->id] == ~0) {
				chartMeshIndices[vertex->id] = m_chartMesh->vertexCount();
				m_chartToOriginalMap.push_back(vertex->id);
				m_chartToUnifiedMap.push_back(unifiedMeshIndices[unifiedVertex->id]);
				HalfEdge::Vertex *v = m_chartMesh->addVertex(vertex->pos);
				v->nor = vertex->nor;
				v->tex = vertex->tex;
			}
		}
	}
	// This is ignoring the canonical map:
	// - Is it really necessary to link colocals?
	m_chartMesh->linkColocals();
	//m_unifiedMesh->linkColocals();  // Not strictly necessary, no colocals in the unified mesh. # Wrong.
	// This check is not valid anymore, if the original mesh vertices were linked with a canonical map, then it might have
	// some colocal vertices that were unlinked. So, the unified mesh might have some duplicate vertices, because firstColocal()
	// is not guaranteed to return the same vertex for two colocal vertices.
	//nvCheck(m_chartMesh->colocalVertexCount() == m_unifiedMesh->vertexCount());
	// Is that OK? What happens in meshes were that happens? Does anything break? Apparently not...
	std::vector<uint32_t> faceIndices;
	faceIndices.reserve(7);
	// Add faces.
	for (uint32_t f = 0; f < faceCount; f++) {
		const HalfEdge::Face *face = originalMesh->faceAt(faceArray[f]);
		nvDebugCheck(face != NULL);
		faceIndices.clear();
		for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance()) {
			const HalfEdge::Vertex *vertex = it.current()->vertex;
			nvDebugCheck(vertex != NULL);
			faceIndices.push_back(chartMeshIndices[vertex->id]);
		}
		m_chartMesh->addFace(faceIndices);
		faceIndices.clear();
		for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance()) {
			const HalfEdge::Vertex *vertex = it.current()->vertex;
			nvDebugCheck(vertex != NULL);
			vertex = vertex->firstColocal();
			faceIndices.push_back(unifiedMeshIndices[vertex->id]);
		}
		m_unifiedMesh->addFace(faceIndices);
	}
	m_chartMesh->linkBoundary();
	m_unifiedMesh->linkBoundary();
	//exportMesh(m_unifiedMesh.ptr(), "debug_input.obj");
	if (m_unifiedMesh->splitBoundaryEdges()) {
		m_unifiedMesh.reset(unifyVertices(m_unifiedMesh.get()));
	}
	//exportMesh(m_unifiedMesh.ptr(), "debug_split.obj");
	// Closing the holes is not always the best solution and does not fix all the problems.
	// We need to do some analysis of the holes and the genus to:
	// - Find cuts that reduce genus.
	// - Find cuts to connect holes.
	// - Use minimal spanning trees or seamster.
	if (!closeHoles()) {
		/*static int pieceCount = 0;
		StringBuilder fileName;
		fileName.format("debug_hole_%d.obj", pieceCount++);
		exportMesh(m_unifiedMesh.ptr(), fileName.str());*/
	}
	m_unifiedMesh.reset(triangulate(m_unifiedMesh.get()));
	//exportMesh(m_unifiedMesh.ptr(), "debug_triangulated.obj");
	// Analyze chart topology.
	MeshTopology topology(m_unifiedMesh.get());
	m_isDisk = topology.isDisk();
}


void Chart::buildVertexMap(const HalfEdge::Mesh *originalMesh, const std::vector<uint32_t> &unchartedMaterialArray)
{
	nvCheck(m_chartMesh.get() == NULL && m_unifiedMesh.get() == NULL);
	m_isVertexMapped = true;
	// Build face indices.
	m_faceArray.clear();
	const uint32_t meshFaceCount = originalMesh->faceCount();
	for (uint32_t f = 0; f < meshFaceCount; f++) {
		const HalfEdge::Face *face = originalMesh->faceAt(f);
		if (std::find(unchartedMaterialArray.begin(), unchartedMaterialArray.end(), face->material) != unchartedMaterialArray.end()) {
			m_faceArray.push_back(f);
		}
	}
	const uint32_t faceCount = m_faceArray.size();
	if (faceCount == 0) {
		return;
	}
	// @@ The chartMesh construction is basically the same as with regular charts, don't duplicate!
	const uint32_t meshVertexCount = originalMesh->vertexCount();
	m_chartMesh.reset(new HalfEdge::Mesh());
	std::vector<uint32_t> chartMeshIndices(meshVertexCount, ~0);
	// Vertex map mesh only has disconnected vertices.
	for (uint32_t f = 0; f < faceCount; f++) {
		const HalfEdge::Face *face = originalMesh->faceAt(m_faceArray[f]);
		nvDebugCheck(face != NULL);
		for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance()) {
			const HalfEdge::Vertex *vertex = it.current()->vertex;
			if (chartMeshIndices[vertex->id] == ~0) {
				chartMeshIndices[vertex->id] = m_chartMesh->vertexCount();
				m_chartToOriginalMap.push_back(vertex->id);
				HalfEdge::Vertex *v = m_chartMesh->addVertex(vertex->pos);
				v->nor = vertex->nor;
				v->tex = vertex->tex; // @@ Not necessary.
			}
		}
	}
	// @@ Link colocals using the original mesh canonical map? Build canonical map on the fly? Do we need to link colocals at all for this?
	//m_chartMesh->linkColocals();
	std::vector<uint32_t> faceIndices;
	faceIndices.reserve(7);
	// Add faces.
	for (uint32_t f = 0; f < faceCount; f++) {
		const HalfEdge::Face *face = originalMesh->faceAt(m_faceArray[f]);
		nvDebugCheck(face != NULL);
		faceIndices.clear();
		for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance()) {
			const HalfEdge::Vertex *vertex = it.current()->vertex;
			nvDebugCheck(vertex != NULL);
			nvDebugCheck(chartMeshIndices[vertex->id] != ~0);
			faceIndices.push_back(chartMeshIndices[vertex->id]);
		}
		HalfEdge::Face *new_face = m_chartMesh->addFace(faceIndices);
		nvDebugCheck(new_face != NULL);
	}
	m_chartMesh->linkBoundary();
	const uint32_t chartVertexCount = m_chartMesh->vertexCount();
	Box bounds;
	bounds.clearBounds();
	for (uint32_t i = 0; i < chartVertexCount; i++) {
		HalfEdge::Vertex *vertex = m_chartMesh->vertexAt(i);
		bounds.addPointToBounds(vertex->pos);
	}
	ProximityGrid grid;
	grid.init(bounds, chartVertexCount);
	for (uint32_t i = 0; i < chartVertexCount; i++) {
		HalfEdge::Vertex *vertex = m_chartMesh->vertexAt(i);
		grid.add(vertex->pos, i);
	}
	uint32_t texelCount = 0;
	const float positionThreshold = 0.01f;
	const float normalThreshold = 0.01f;
	uint32_t verticesVisited = 0;
	uint32_t cellsVisited = 0;
	std::vector<int> vertexIndexArray(chartVertexCount, -1); // Init all indices to -1.
	// Traverse vertices in morton order. @@ It may be more interesting to sort them based on orientation.
	const uint32_t cellCodeCount = grid.mortonCount();
	for (uint32_t cellCode = 0; cellCode < cellCodeCount; cellCode++) {
		int cell = grid.mortonIndex(cellCode);
		if (cell < 0) continue;
		cellsVisited++;
		const std::vector<uint32_t> &indexArray = grid.cellArray[cell].indexArray;
		for (uint32_t i = 0; i < indexArray.size(); i++) {
			uint32_t idx = indexArray[i];
			HalfEdge::Vertex *vertex = m_chartMesh->vertexAt(idx);
			nvDebugCheck(vertexIndexArray[idx] == -1);
			std::vector<uint32_t> neighbors;
			grid.gather(vertex->pos, positionThreshold, /*ref*/neighbors);
			// Compare against all nearby vertices, cluster greedily.
			for (uint32_t j = 0; j < neighbors.size(); j++) {
				uint32_t otherIdx = neighbors[j];
				if (vertexIndexArray[otherIdx] != -1) {
					HalfEdge::Vertex *otherVertex = m_chartMesh->vertexAt(otherIdx);
					if (distance(vertex->pos, otherVertex->pos) < positionThreshold &&
					        distance(vertex->nor, otherVertex->nor) < normalThreshold) {
						vertexIndexArray[idx] = vertexIndexArray[otherIdx];
						break;
					}
				}
			}
			// If index not assigned, assign new one.
			if (vertexIndexArray[idx] == -1) {
				vertexIndexArray[idx] = texelCount++;
			}
			verticesVisited++;
		}
	}
	nvDebugCheck(cellsVisited == grid.cellArray.size());
	nvDebugCheck(verticesVisited == chartVertexCount);
	vertexMapWidth = ftoi_ceil(sqrtf(float(texelCount)));
	vertexMapWidth = (vertexMapWidth + 3) & ~3;                             // Width aligned to 4.
	vertexMapHeight = vertexMapWidth == 0 ? 0 : (texelCount + vertexMapWidth - 1) / vertexMapWidth;
	//vertexMapHeight = (vertexMapHeight + 3) & ~3;                           // Height aligned to 4.
	nvDebugCheck(vertexMapWidth >= vertexMapHeight);
	nvDebug("Reduced vertex count from %d to %d.\n", chartVertexCount, texelCount);
	// Lay down the clustered vertices in morton order.
	std::vector<uint32_t> texelCodes(texelCount);
	// For each texel, assign one morton code.
	uint32_t texelCode = 0;
	for (uint32_t i = 0; i < texelCount; i++) {
		uint32_t x, y;
		do {
			x = decodeMorton2X(texelCode);
			y = decodeMorton2Y(texelCode);
			texelCode++;
		} while (x >= uint32_t(vertexMapWidth) || y >= uint32_t(vertexMapHeight));
		texelCodes[i] = texelCode - 1;
	}
	for (uint32_t i = 0; i < chartVertexCount; i++) {
		HalfEdge::Vertex *vertex = m_chartMesh->vertexAt(i);
		int idx = vertexIndexArray[i];
		if (idx != -1) {
			uint32_t texelCode = texelCodes[idx];
			uint32_t x = decodeMorton2X(texelCode);
			uint32_t y = decodeMorton2Y(texelCode);
			vertex->tex.x = float(x);
			vertex->tex.y = float(y);
		}
	}
}



static void getBoundaryEdges(HalfEdge::Mesh *mesh, std::vector<HalfEdge::Edge *> &boundaryEdges)
{
	nvDebugCheck(mesh != NULL);
	const uint32_t edgeCount = mesh->edgeCount();
	BitArray bitFlags(edgeCount);
	bitFlags.clearAll();
	boundaryEdges.clear();
	// Search for boundary edges. Mark all the edges that belong to the same boundary.
	for (uint32_t e = 0; e < edgeCount; e++) {
		HalfEdge::Edge *startEdge = mesh->edgeAt(e);
		if (startEdge != NULL && startEdge->isBoundary() && bitFlags.bitAt(e) == false) {
			nvDebugCheck(startEdge->face != NULL);
			nvDebugCheck(startEdge->pair->face == NULL);
			startEdge = startEdge->pair;
			const HalfEdge::Edge *edge = startEdge;
			do {
				nvDebugCheck(edge->face == NULL);
				nvDebugCheck(bitFlags.bitAt(edge->id / 2) == false);
				bitFlags.setBitAt(edge->id / 2);
				edge = edge->next;
			} while (startEdge != edge);
			boundaryEdges.push_back(startEdge);
		}
	}
}


bool Chart::closeLoop(uint32_t start, const std::vector<HalfEdge::Edge *> &loop)
{
	const uint32_t vertexCount = loop.size() - start;
	nvDebugCheck(vertexCount >= 3);
	if (vertexCount < 3) return false;
	nvDebugCheck(loop[start]->vertex->isColocal(loop[start + vertexCount - 1]->to()));
	// If the hole is planar, then we add a single face that will be properly triangulated later.
	// If the hole is not planar, we add a triangle fan with a vertex at the hole centroid.
	// This is still a bit of a hack. There surely are better hole filling algorithms out there.
	std::vector<Vector3> points(vertexCount);
	for (uint32_t i = 0; i < vertexCount; i++) {
		points[i] = loop[start + i]->vertex->pos;
	}
	bool isPlanar = Fit::isPlanar(vertexCount, points.data());
	if (isPlanar) {
		// Add face and connect edges.
		HalfEdge::Face *face = m_unifiedMesh->addFace();
		for (uint32_t i = 0; i < vertexCount; i++) {
			HalfEdge::Edge *edge = loop[start + i];
			edge->face = face;
			edge->setNext(loop[start + (i + 1) % vertexCount]);
		}
		face->edge = loop[start];
		nvDebugCheck(face->isValid());
	} else {
		// If the polygon is not planar, we just cross our fingers, and hope this will work:
		// Compute boundary centroid:
		Vector3 centroidPos(0);
		for (uint32_t i = 0; i < vertexCount; i++) {
			centroidPos += points[i];
		}
		centroidPos *= (1.0f / vertexCount);
		HalfEdge::Vertex *centroid = m_unifiedMesh->addVertex(centroidPos);
		// Add one pair of edges for each boundary vertex.
		for (uint32_t j = vertexCount - 1, i = 0; i < vertexCount; j = i++) {
			HalfEdge::Face *face = m_unifiedMesh->addFace(centroid->id, loop[start + j]->vertex->id, loop[start + i]->vertex->id);
			nvDebugCheck(face != NULL);
		}
	}
	return true;
}


bool Chart::closeHoles()
{
	nvDebugCheck(!m_isVertexMapped);
	std::vector<HalfEdge::Edge *> boundaryEdges;
	getBoundaryEdges(m_unifiedMesh.get(), boundaryEdges);
	uint32_t boundaryCount = boundaryEdges.size();
	if (boundaryCount <= 1) {
		// Nothing to close.
		return true;
	}
	// Compute lengths and areas.
	std::vector<float> boundaryLengths;
	for (uint32_t i = 0; i < boundaryCount; i++) {
		const HalfEdge::Edge *startEdge = boundaryEdges[i];
		nvCheck(startEdge->face == NULL);
		//float boundaryEdgeCount = 0;
		float boundaryLength = 0.0f;
		//Vector3 boundaryCentroid(zero);
		const HalfEdge::Edge *edge = startEdge;
		do {
			Vector3 t0 = edge->from()->pos;
			Vector3 t1 = edge->to()->pos;
			//boundaryEdgeCount++;
			boundaryLength += length(t1 - t0);
			//boundaryCentroid += edge->vertex()->pos;
			edge = edge->next;
		} while (edge != startEdge);
		boundaryLengths.push_back(boundaryLength);
		//boundaryCentroids.append(boundaryCentroid / boundaryEdgeCount);
	}
	// Find disk boundary.
	uint32_t diskBoundary = 0;
	float maxLength = boundaryLengths[0];
	for (uint32_t i = 1; i < boundaryCount; i++) {
		if (boundaryLengths[i] > maxLength) {
			maxLength = boundaryLengths[i];
			diskBoundary = i;
		}
	}
	// Close holes.
	for (uint32_t i = 0; i < boundaryCount; i++) {
		if (diskBoundary == i) {
			// Skip disk boundary.
			continue;
		}
		HalfEdge::Edge *startEdge = boundaryEdges[i];
		nvDebugCheck(startEdge != NULL);
		nvDebugCheck(startEdge->face == NULL);
		std::vector<HalfEdge::Vertex *> vertexLoop;
		std::vector<HalfEdge::Edge *> edgeLoop;
		HalfEdge::Edge *edge = startEdge;
		do {
			HalfEdge::Vertex *vertex = edge->next->vertex;  // edge->to()
			uint32_t i;
			for (i = 0; i < vertexLoop.size(); i++) {
				if (vertex->isColocal(vertexLoop[i])) {
					break;
				}
			}
			bool isCrossing = (i != vertexLoop.size());
			if (isCrossing) {
				HalfEdge::Edge *prev = edgeLoop[i];     // Previous edge before the loop.
				HalfEdge::Edge *next = edge->next;    // Next edge after the loop.
				nvDebugCheck(prev->to()->isColocal(next->from()));
				// Close loop.
				edgeLoop.push_back(edge);
				closeLoop(i + 1, edgeLoop);
				// Link boundary loop.
				prev->setNext(next);
				vertex->setEdge(next);
				// Start over again.
				vertexLoop.clear();
				edgeLoop.clear();
				edge = startEdge;
				vertex = edge->to();
			}
			vertexLoop.push_back(vertex);
			edgeLoop.push_back(edge);
			edge = edge->next;
		} while (edge != startEdge);
		closeLoop(0, edgeLoop);
	}
	getBoundaryEdges(m_unifiedMesh.get(), boundaryEdges);
	boundaryCount = boundaryEdges.size();
	nvDebugCheck(boundaryCount == 1);
	return boundaryCount == 1;
}


// Transfer parameterization from unified mesh to chart mesh.
void Chart::transferParameterization()
{
	nvDebugCheck(!m_isVertexMapped);
	uint32_t vertexCount = m_chartMesh->vertexCount();
	for (uint32_t v = 0; v < vertexCount; v++) {
		HalfEdge::Vertex *vertex = m_chartMesh->vertexAt(v);
		HalfEdge::Vertex *unifiedVertex = m_unifiedMesh->vertexAt(mapChartVertexToUnifiedVertex(v));
		vertex->tex = unifiedVertex->tex;
	}
}

float Chart::computeSurfaceArea() const
{
	return nv::computeSurfaceArea(m_chartMesh.get()) * scale;
}

float Chart::computeParametricArea() const
{
	// This only makes sense in parameterized meshes.
	nvDebugCheck(m_isDisk);
	nvDebugCheck(!m_isVertexMapped);
	return nv::computeParametricArea(m_chartMesh.get());
}

Vector2 Chart::computeParametricBounds() const
{
	// This only makes sense in parameterized meshes.
	nvDebugCheck(m_isDisk);
	nvDebugCheck(!m_isVertexMapped);
	Box bounds;
	bounds.clearBounds();
	uint32_t vertexCount = m_chartMesh->vertexCount();
	for (uint32_t v = 0; v < vertexCount; v++) {
		HalfEdge::Vertex *vertex = m_chartMesh->vertexAt(v);
		bounds.addPointToBounds(Vector3(vertex->tex, 0));
	}
	return bounds.extents().xy();
}
