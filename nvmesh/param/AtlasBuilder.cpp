// This code is in the public domain -- castano@gmail.com

#include "AtlasBuilder.h"
#include "Util.h"

#include "xatlas.h"

//#include "nvcore/IntroSort.h"

#include <algorithm> // std::sort

#include <float.h> // FLT_MAX
#include <limits.h> // UINT_MAX

using namespace nv;

namespace
{

// Dummy implementation of a priority queue using sort at insertion.
// - Insertion is o(n)
// - Smallest element goes at the end, so that popping it is o(1).
// - Resorting is n*log(n)
// @@ Number of elements in the queue is usually small, and we'd have to rebalance often. I'm not sure it's worth implementing a heap.
// @@ Searcing at removal would remove the need for sorting when priorities change.
struct PriorityQueue {
	PriorityQueue(uint32_t size = UINT_MAX) : maxSize(size) {}

	void push(float priority, uint32_t face)
	{
		uint32_t i = 0;
		const uint32_t count = pairs.size();
		for (; i < count; i++) {
			if (pairs[i].priority > priority) break;
		}
		Pair p = { priority, face };
		pairs.insert(pairs.begin() + i, p);
		if (pairs.size() > maxSize) {
			pairs.erase(pairs.begin());
		}
	}

	// push face out of order, to be sorted later.
	void push(uint32_t face)
	{
		Pair p = { 0.0f, face };
		pairs.push_back(p);
	}

	uint32_t pop()
	{
		uint32_t f = pairs.back().face;
		pairs.pop_back();
		return f;
	}

	void sort()
	{
		//nv::sort(pairs); // @@ My intro sort appears to be much slower than it should!
		std::sort(pairs.begin(), pairs.end());
	}

	void clear()
	{
		pairs.clear();
	}

	uint32_t count() const
	{
		return pairs.size();
	}

	float firstPriority() const
	{
		return pairs.back().priority;
	}


	const uint32_t maxSize;

	struct Pair {
		bool operator <(const Pair &p) const
		{
			return priority > p.priority;    // !! Sort in inverse priority order!
		}
		float priority;
		uint32_t face;
	};


	std::vector<Pair> pairs;
};

static bool isNormalSeam(const HalfEdge::Edge *edge)
{
	return (edge->vertex->nor != edge->pair->next->vertex->nor || edge->next->vertex->nor != edge->pair->vertex->nor);
}

static bool isTextureSeam(const HalfEdge::Edge *edge)
{
	return (edge->vertex->tex != edge->pair->next->vertex->tex || edge->next->vertex->tex != edge->pair->vertex->tex);
}

} // namespace


struct nv::ChartBuildData {
	ChartBuildData(int id) : id(id)
	{
		planeNormal = Vector3(0);
		centroid = Vector3(0);
		coneAxis = Vector3(0);
		coneAngle = 0;
		area = 0;
		boundaryLength = 0;
		normalSum = Vector3(0);
		centroidSum = Vector3(0);
	}

	int id;

	// Proxy info:
	Vector3 planeNormal;
	Vector3 centroid;
	Vector3 coneAxis;
	float coneAngle;

	float area;
	float boundaryLength;
	Vector3 normalSum;
	Vector3 centroidSum;

	std::vector<uint32_t> seeds;  // @@ These could be a pointers to the HalfEdge faces directly.
	std::vector<uint32_t> faces;
	PriorityQueue candidates;
};



AtlasBuilder::AtlasBuilder(const HalfEdge::Mesh *m) : mesh(m), facesLeft(m->faceCount())
{
	const uint32_t faceCount = m->faceCount();
	faceChartArray.resize(faceCount, -1);
	faceCandidateArray.resize(faceCount, -1);
	// @@ Floyd for the whole mesh is too slow. We could compute floyd progressively per patch as the patch grows. We need a better solution to compute most central faces.
	//computeShortestPaths();
	// Precompute edge lengths and face areas.
	uint32_t edgeCount = m->edgeCount();
	edgeLengths.resize(edgeCount);
	for (uint32_t i = 0; i < edgeCount; i++) {
		uint32_t id = m->edgeAt(i)->id;
		nvDebugCheck(id / 2 == i);
		edgeLengths[i] = m->edgeAt(i)->length();
	}
	faceAreas.resize(faceCount);
	for (uint32_t i = 0; i < faceCount; i++) {
		faceAreas[i] = m->faceAt(i)->area();
	}
}

AtlasBuilder::~AtlasBuilder()
{
	const uint32_t chartCount = chartArray.size();
	for (uint32_t i = 0; i < chartCount; i++) {
		delete chartArray[i];
	}
}


void AtlasBuilder::markUnchartedFaces(const std::vector<uint32_t> &unchartedFaces)
{
	const uint32_t unchartedFaceCount = unchartedFaces.size();
	for (uint32_t i = 0; i < unchartedFaceCount; i++) {
		uint32_t f = unchartedFaces[i];
		faceChartArray[f] = -2;
		//faceCandidateArray[f] = -2; // @@ ?
		removeCandidate(f);
	}
	nvDebugCheck(facesLeft >= unchartedFaceCount);
	facesLeft -= unchartedFaceCount;
}


void AtlasBuilder::computeShortestPaths()
{
	const uint32_t faceCount = mesh->faceCount();
	shortestPaths.resize(faceCount * faceCount, FLT_MAX);
	// Fill edges:
	for (uint32_t i = 0; i < faceCount; i++) {
		shortestPaths[i * faceCount + i] = 0.0f;
		const HalfEdge::Face *face_i = mesh->faceAt(i);
		Vector3 centroid_i = face_i->centroid();
		for (HalfEdge::Face::ConstEdgeIterator it(face_i->edges()); !it.isDone(); it.advance()) {
			const HalfEdge::Edge *edge = it.current();
			if (!edge->isBoundary()) {
				const HalfEdge::Face *face_j = edge->pair->face;
				uint32_t j = face_j->id;
				Vector3 centroid_j = face_j->centroid();
				shortestPaths[i * faceCount + j] = shortestPaths[j * faceCount + i] = length(centroid_i - centroid_j);
			}
		}
	}
	// Use Floyd-Warshall algorithm to compute all paths:
	for (uint32_t k = 0; k < faceCount; k++) {
		for (uint32_t i = 0; i < faceCount; i++) {
			for (uint32_t j = 0; j < faceCount; j++) {
				shortestPaths[i * faceCount + j] = std::min(shortestPaths[i * faceCount + j], shortestPaths[i * faceCount + k] + shortestPaths[k * faceCount + j]);
			}
		}
	}
}


void AtlasBuilder::placeSeeds(float threshold, uint32_t maxSeedCount)
{
	// Instead of using a predefiened number of seeds:
	// - Add seeds one by one, growing chart until a certain treshold.
	// - Undo charts and restart growing process.
	// @@ How can we give preference to faces far from sharp features as in the LSCM paper?
	//   - those points can be found using a simple flood filling algorithm.
	//   - how do we weight the probabilities?
	for (uint32_t i = 0; i < maxSeedCount; i++) {
		if (facesLeft == 0) {
			// No faces left, stop creating seeds.
			break;
		}
		createRandomChart(threshold);
	}
}


void AtlasBuilder::createRandomChart(float threshold)
{
	ChartBuildData *chart = new ChartBuildData(chartArray.size());
	chartArray.push_back(chart);
	// Pick random face that is not used by any chart yet.
	uint32_t randomFaceIdx = rand.getRange(facesLeft - 1);
	uint32_t i = 0;
	for (uint32_t f = 0; f != randomFaceIdx; f++, i++) {
		while (faceChartArray[i] != -1) i++;
	}
	while (faceChartArray[i] != -1) i++;
	chart->seeds.push_back(i);
	addFaceToChart(chart, i, true);
	// Grow the chart as much as possible within the given threshold.
	growChart(chart, threshold * 0.5f, facesLeft);
	//growCharts(threshold - threshold * 0.75f / chartCount(), facesLeft);
}

void AtlasBuilder::addFaceToChart(ChartBuildData *chart, uint32_t f, bool recomputeProxy)
{
	// Add face to chart.
	chart->faces.push_back(f);
	nvDebugCheck(faceChartArray[f] == -1);
	faceChartArray[f] = chart->id;
	facesLeft--;
	// Update area and boundary length.
	chart->area = evaluateChartArea(chart, f);
	chart->boundaryLength = evaluateBoundaryLength(chart, f);
	chart->normalSum = evaluateChartNormalSum(chart, f);
	chart->centroidSum = evaluateChartCentroidSum(chart, f);
	if (recomputeProxy) {
		// Update proxy and candidate's priorities.
		updateProxy(chart);
	}
	// Update candidates.
	removeCandidate(f);
	updateCandidates(chart, f);
	updatePriorities(chart);
}

// @@ Get N best candidates in one pass.
const AtlasBuilder::Candidate &AtlasBuilder::getBestCandidate() const
{
	uint32_t best = 0;
	float bestCandidateMetric = FLT_MAX;
	const uint32_t candidateCount = candidateArray.size();
	nvCheck(candidateCount > 0);
	for (uint32_t i = 0; i < candidateCount; i++) {
		const Candidate &candidate = candidateArray[i];
		if (candidate.metric < bestCandidateMetric) {
			bestCandidateMetric = candidate.metric;
			best = i;
		}
	}
	return candidateArray[best];
}


// Returns true if any of the charts can grow more.
bool AtlasBuilder::growCharts(float threshold, uint32_t faceCount)
{
	// Using one global list.
	faceCount = std::min(faceCount, facesLeft);
	for (uint32_t i = 0; i < faceCount; i++) {
		const Candidate &candidate = getBestCandidate();
		if (candidate.metric > threshold) {
			return false; // Can't grow more.
		}
		addFaceToChart(candidate.chart, candidate.face);
	}
	return facesLeft != 0; // Can continue growing.
}

bool AtlasBuilder::growChart(ChartBuildData *chart, float threshold, uint32_t faceCount)
{
	// Try to add faceCount faces within threshold to chart.
	for (uint32_t i = 0; i < faceCount; ) {
		if (chart->candidates.count() == 0 || chart->candidates.firstPriority() > threshold) {
			return false;
		}
		uint32_t f = chart->candidates.pop();
		if (faceChartArray[f] == -1) {
			addFaceToChart(chart, f);
			i++;
		}
	}
	if (chart->candidates.count() == 0 || chart->candidates.firstPriority() > threshold) {
		return false;
	}
	return true;
}


void AtlasBuilder::resetCharts()
{
	const uint32_t faceCount = mesh->faceCount();
	for (uint32_t i = 0; i < faceCount; i++) {
		faceChartArray[i] = -1;
		faceCandidateArray[i] = -1;
	}
	facesLeft = faceCount;
	candidateArray.clear();
	const uint32_t chartCount = chartArray.size();
	for (uint32_t i = 0; i < chartCount; i++) {
		ChartBuildData *chart = chartArray[i];
		const uint32_t seed = chart->seeds.back();
		chart->area = 0.0f;
		chart->boundaryLength = 0.0f;
		chart->normalSum = Vector3(0);
		chart->centroidSum = Vector3(0);
		chart->faces.clear();
		chart->candidates.clear();
		addFaceToChart(chart, seed);
	}
}


void AtlasBuilder::updateCandidates(ChartBuildData *chart, uint32_t f)
{
	const HalfEdge::Face *face = mesh->faceAt(f);
	// Traverse neighboring faces, add the ones that do not belong to any chart yet.
	for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance()) {
		const HalfEdge::Edge *edge = it.current()->pair;
		if (!edge->isBoundary()) {
			uint32_t f = edge->face->id;
			if (faceChartArray[f] == -1) {
				chart->candidates.push(f);
			}
		}
	}
}


void AtlasBuilder::updateProxies()
{
	const uint32_t chartCount = chartArray.size();
	for (uint32_t i = 0; i < chartCount; i++) {
		updateProxy(chartArray[i]);
	}
}


namespace
{

// Unnormalized face normal assuming it's a triangle.
static Vector3 triangleNormal(const HalfEdge::Face *face)
{
	Vector3 p0 = face->edge->vertex->pos;
	Vector3 p1 = face->edge->next->vertex->pos;
	Vector3 p2 = face->edge->next->next->vertex->pos;
	Vector3 e0 = p2 - p0;
	Vector3 e1 = p1 - p0;
	return normalizeSafe(cross(e0, e1), Vector3(0), 0.0f);
}

static Vector3 triangleNormalAreaScaled(const HalfEdge::Face *face)
{
	Vector3 p0 = face->edge->vertex->pos;
	Vector3 p1 = face->edge->next->vertex->pos;
	Vector3 p2 = face->edge->next->next->vertex->pos;
	Vector3 e0 = p2 - p0;
	Vector3 e1 = p1 - p0;
	return cross(e0, e1);
}

// Average of the edge midpoints weighted by the edge length.
// I want a point inside the triangle, but closer to the cirumcenter.
static Vector3 triangleCenter(const HalfEdge::Face *face)
{
	Vector3 p0 = face->edge->vertex->pos;
	Vector3 p1 = face->edge->next->vertex->pos;
	Vector3 p2 = face->edge->next->next->vertex->pos;
	float l0 = length(p1 - p0);
	float l1 = length(p2 - p1);
	float l2 = length(p0 - p2);
	Vector3 m0 = (p0 + p1) * l0 / (l0 + l1 + l2);
	Vector3 m1 = (p1 + p2) * l1 / (l0 + l1 + l2);
	Vector3 m2 = (p2 + p0) * l2 / (l0 + l1 + l2);
	return m0 + m1 + m2;
}

} // namespace

void AtlasBuilder::updateProxy(ChartBuildData *chart)
{
	//#pragma message(NV_FILE_LINE "TODO: Use best fit plane instead of average normal.")
	chart->planeNormal = normalizeSafe(chart->normalSum, Vector3(0), 0.0f);
	chart->centroid = chart->centroidSum / float(chart->faces.size());
}



bool AtlasBuilder::relocateSeeds()
{
	bool anySeedChanged = false;
	const uint32_t chartCount = chartArray.size();
	for (uint32_t i = 0; i < chartCount; i++) {
		if (relocateSeed(chartArray[i])) {
			anySeedChanged = true;
		}
	}
	return anySeedChanged;
}


bool AtlasBuilder::relocateSeed(ChartBuildData *chart)
{
	Vector3 centroid = computeChartCentroid(chart);
	const uint32_t N = 10;  // @@ Hardcoded to 10?
	PriorityQueue bestTriangles(N);
	// Find the first N triangles that fit the proxy best.
	const uint32_t faceCount = chart->faces.size();
	for (uint32_t i = 0; i < faceCount; i++) {
		float priority = evaluateProxyFitMetric(chart, chart->faces[i]);
		bestTriangles.push(priority, chart->faces[i]);
	}
	// Of those, choose the most central triangle.
	uint32_t mostCentral;
	float maxDistance = -1;
	const uint32_t bestCount = bestTriangles.count();
	for (uint32_t i = 0; i < bestCount; i++) {
		const HalfEdge::Face *face = mesh->faceAt(bestTriangles.pairs[i].face);
		Vector3 faceCentroid = triangleCenter(face);
		float distance = length(centroid - faceCentroid);
		if (distance > maxDistance) {
			maxDistance = distance;
			mostCentral = bestTriangles.pairs[i].face;
		}
	}
	nvDebugCheck(maxDistance >= 0);
	// In order to prevent k-means cyles we record all the previously chosen seeds.
	uint32_t index = std::find(chart->seeds.begin(), chart->seeds.end(), mostCentral) - chart->seeds.begin();
	if (index < chart->seeds.size()) {
		// Move new seed to the end of the seed array.
		uint32_t last = chart->seeds.size() - 1;
		std::swap(chart->seeds[index], chart->seeds[last]);
		return false;
	} else {
		// Append new seed.
		chart->seeds.push_back(mostCentral);
		return true;
	}
}

void AtlasBuilder::removeCandidate(uint32_t f)
{
	int c = faceCandidateArray[f];
	if (c != -1) {
		faceCandidateArray[f] = -1;
		if (c == candidateArray.size() - 1) {
			candidateArray.pop_back();
		} else {
			// Replace with last.
			candidateArray[c] = candidateArray[candidateArray.size() - 1];
			candidateArray.pop_back();
			faceCandidateArray[candidateArray[c].face] = c;
		}
	}
}

void AtlasBuilder::updateCandidate(ChartBuildData *chart, uint32_t f, float metric)
{
	if (faceCandidateArray[f] == -1) {
		const uint32_t index = candidateArray.size();
		faceCandidateArray[f] = index;
		candidateArray.resize(index + 1);
		candidateArray[index].face = f;
		candidateArray[index].chart = chart;
		candidateArray[index].metric = metric;
	} else {
		int c = faceCandidateArray[f];
		nvDebugCheck(c != -1);
		Candidate &candidate = candidateArray[c];
		nvDebugCheck(candidate.face == f);
		if (metric < candidate.metric || chart == candidate.chart) {
			candidate.metric = metric;
			candidate.chart = chart;
		}
	}
}


void AtlasBuilder::updatePriorities(ChartBuildData *chart)
{
	// Re-evaluate candidate priorities.
	uint32_t candidateCount = chart->candidates.count();
	for (uint32_t i = 0; i < candidateCount; i++) {
		chart->candidates.pairs[i].priority = evaluatePriority(chart, chart->candidates.pairs[i].face);
		if (faceChartArray[chart->candidates.pairs[i].face] == -1) {
			updateCandidate(chart, chart->candidates.pairs[i].face, chart->candidates.pairs[i].priority);
		}
	}
	// Sort candidates.
	chart->candidates.sort();
}


// Evaluate combined metric.
float AtlasBuilder::evaluatePriority(ChartBuildData *chart, uint32_t face)
{
	// Estimate boundary length and area:
	float newBoundaryLength = evaluateBoundaryLength(chart, face);
	float newChartArea = evaluateChartArea(chart, face);
	float F = evaluateProxyFitMetric(chart, face);
	float C = evaluateRoundnessMetric(chart, face, newBoundaryLength, newChartArea);
	float P = evaluateStraightnessMetric(chart, face);
	// Penalize faces that cross seams, reward faces that close seams or reach boundaries.
	float N = evaluateNormalSeamMetric(chart, face);
	float T = evaluateTextureSeamMetric(chart, face);
	//float R = evaluateCompletenessMetric(chart, face);
	//float D = evaluateDihedralAngleMetric(chart, face);
	// @@ Add a metric based on local dihedral angle.
	// @@ Tweaking the normal and texture seam metrics.
	// - Cause more impedance. Never cross 90 degree edges.
	// -
	float cost = float(
	                 settings.proxyFitMetricWeight * F +
	                 settings.roundnessMetricWeight * C +
	                 settings.straightnessMetricWeight * P +
	                 settings.normalSeamMetricWeight * N +
	                 settings.textureSeamMetricWeight * T);
	// Enforce limits strictly:
	if (newChartArea > settings.maxChartArea) cost = FLT_MAX;
	if (newBoundaryLength > settings.maxBoundaryLength) cost = FLT_MAX;
	// Make sure normal seams are fully respected:
	if (settings.normalSeamMetricWeight >= 1000 && N != 0) cost = FLT_MAX;
	nvCheck(std::isfinite(cost));
	return cost;
}


// Returns a value in [0-1].
float AtlasBuilder::evaluateProxyFitMetric(ChartBuildData *chart, uint32_t f)
{
	const HalfEdge::Face *face = mesh->faceAt(f);
	Vector3 faceNormal = triangleNormal(face);
	// Use plane fitting metric for now:
	return 1 - dot(faceNormal, chart->planeNormal); // @@ normal deviations should be weighted by face area
}

float AtlasBuilder::evaluateDistanceToBoundary(ChartBuildData *chart, uint32_t face)
{
//#pragma message(NV_FILE_LINE "TODO: Evaluate distance to boundary metric.")
	// @@ This is needed for the seed relocation code.
	// @@ This could provide a better roundness metric.
	return 0.0f;
}

float AtlasBuilder::evaluateDistanceToSeed(ChartBuildData *chart, uint32_t f)
{
	const HalfEdge::Face *seed = mesh->faceAt(chart->seeds.back());
	const HalfEdge::Face *face = mesh->faceAt(f);
	return length(triangleCenter(seed) - triangleCenter(face));
}


float AtlasBuilder::evaluateRoundnessMetric(ChartBuildData *chart, uint32_t face, float newBoundaryLength, float newChartArea)
{
	float roundness = square(chart->boundaryLength) / chart->area;
	float newRoundness = square(newBoundaryLength) / newChartArea;
	if (newRoundness > roundness) {
		return square(newBoundaryLength) / (newChartArea * 4 * PI);
	} else {
		// Offer no impedance to faces that improve roundness.
		return 0;
	}
}

float AtlasBuilder::evaluateStraightnessMetric(ChartBuildData *chart, uint32_t f)
{
	float l_out = 0.0f;
	float l_in = 0.0f;
	const HalfEdge::Face *face = mesh->faceAt(f);
	for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance()) {
		const HalfEdge::Edge *edge = it.current();
		float l = edgeLengths[edge->id / 2];
		if (edge->isBoundary()) {
			l_out += l;
		} else {
			uint32_t neighborFaceId = edge->pair->face->id;
			if (faceChartArray[neighborFaceId] != chart->id) {
				l_out += l;
			} else {
				l_in += l;
			}
		}
	}
	nvDebugCheck(l_in != 0.0f); // Candidate face must be adjacent to chart. @@ This is not true if the input mesh has zero-length edges.
	float ratio = (l_out - l_in) / (l_out + l_in);
	return std::min(ratio, 0.0f); // Only use the straightness metric to close gaps.
}


float AtlasBuilder::evaluateNormalSeamMetric(ChartBuildData *chart, uint32_t f)
{
	float seamFactor = 0.0f;
	float totalLength = 0.0f;
	const HalfEdge::Face *face = mesh->faceAt(f);
	for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance()) {
		const HalfEdge::Edge *edge = it.current();
		if (edge->isBoundary()) {
			continue;
		}
		const uint32_t neighborFaceId = edge->pair->face->id;
		if (faceChartArray[neighborFaceId] != chart->id) {
			continue;
		}
		//float l = edge->length();
		float l = edgeLengths[edge->id / 2];
		totalLength += l;
		if (!edge->isSeam()) {
			continue;
		}
		// Make sure it's a normal seam.
		if (isNormalSeam(edge)) {
			float d0 = clamp(dot(edge->vertex->nor, edge->pair->next->vertex->nor), 0.0f, 1.0f);
			float d1 = clamp(dot(edge->next->vertex->nor, edge->pair->vertex->nor), 0.0f, 1.0f);
			l *= 1 - (d0 + d1) * 0.5f;
			seamFactor += l;
		}
	}
	if (seamFactor == 0) return 0.0f;
	return seamFactor / totalLength;
}


float AtlasBuilder::evaluateTextureSeamMetric(ChartBuildData *chart, uint32_t f)
{
	float seamLength = 0.0f;
	float totalLength = 0.0f;
	const HalfEdge::Face *face = mesh->faceAt(f);
	for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance()) {
		const HalfEdge::Edge *edge = it.current();
		if (edge->isBoundary()) {
			continue;
		}
		const uint32_t neighborFaceId = edge->pair->face->id;
		if (faceChartArray[neighborFaceId] != chart->id) {
			continue;
		}
		//float l = edge->length();
		float l = edgeLengths[edge->id / 2];
		totalLength += l;
		if (!edge->isSeam()) {
			continue;
		}
		// Make sure it's a texture seam.
		if (isTextureSeam(edge)) {
			seamLength += l;
		}
	}
	if (seamLength == 0.0f) {
		return 0.0f; // Avoid division by zero.
	}
	return seamLength / totalLength;
}


float AtlasBuilder::evaluateSeamMetric(ChartBuildData *chart, uint32_t f)
{
	float newSeamLength = 0.0f;
	float oldSeamLength = 0.0f;
	float totalLength = 0.0f;
	const HalfEdge::Face *face = mesh->faceAt(f);
	for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance()) {
		const HalfEdge::Edge *edge = it.current();
		//float l = edge->length();
		float l = edgeLengths[edge->id / 2];
		if (edge->isBoundary()) {
			newSeamLength += l;
		} else {
			if (edge->isSeam()) {
				uint32_t neighborFaceId = edge->pair->face->id;
				if (faceChartArray[neighborFaceId] != chart->id) {
					newSeamLength += l;
				} else {
					oldSeamLength += l;
				}
			}
		}
		totalLength += l;
	}
	return (newSeamLength - oldSeamLength) / totalLength;
}


float AtlasBuilder::evaluateChartArea(ChartBuildData *chart, uint32_t f)
{
	const HalfEdge::Face *face = mesh->faceAt(f);
	//return chart->area + face->area();
	return chart->area + faceAreas[face->id];
}


float AtlasBuilder::evaluateBoundaryLength(ChartBuildData *chart, uint32_t f)
{
	float boundaryLength = chart->boundaryLength;
	// Add new edges, subtract edges shared with the chart.
	const HalfEdge::Face *face = mesh->faceAt(f);
	for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance()) {
		const HalfEdge::Edge *edge = it.current();
		//float edgeLength = edge->length();
		float edgeLength = edgeLengths[edge->id / 2];
		if (edge->isBoundary()) {
			boundaryLength += edgeLength;
		} else {
			uint32_t neighborFaceId = edge->pair->face->id;
			if (faceChartArray[neighborFaceId] != chart->id) {
				boundaryLength += edgeLength;
			} else {
				boundaryLength -= edgeLength;
			}
		}
	}
	return std::max(0.0f, boundaryLength);  // @@ Hack!
}

Vector3 AtlasBuilder::evaluateChartNormalSum(ChartBuildData *chart, uint32_t f)
{
	const HalfEdge::Face *face = mesh->faceAt(f);
	return chart->normalSum + triangleNormalAreaScaled(face);
}

Vector3 AtlasBuilder::evaluateChartCentroidSum(ChartBuildData *chart, uint32_t f)
{
	const HalfEdge::Face *face = mesh->faceAt(f);
	return chart->centroidSum + face->centroid();
}


Vector3 AtlasBuilder::computeChartCentroid(const ChartBuildData *chart)
{
	Vector3 centroid(0);
	const uint32_t faceCount = chart->faces.size();
	for (uint32_t i = 0; i < faceCount; i++) {
		const HalfEdge::Face *face = mesh->faceAt(chart->faces[i]);
		centroid += triangleCenter(face);
	}
	return centroid / float(faceCount);
}


void AtlasBuilder::fillHoles(float threshold)
{
	while (facesLeft > 0) {
		createRandomChart(threshold);
	}
}


void AtlasBuilder::mergeChart(ChartBuildData *owner, ChartBuildData *chart, float sharedBoundaryLength)
{
	const uint32_t faceCount = chart->faces.size();
	for (uint32_t i = 0; i < faceCount; i++) {
		uint32_t f = chart->faces[i];
		nvDebugCheck(faceChartArray[f] == chart->id);
		faceChartArray[f] = owner->id;
		owner->faces.push_back(f);
	}
	// Update adjacencies?
	owner->area += chart->area;
	owner->boundaryLength += chart->boundaryLength - sharedBoundaryLength;
	owner->normalSum += chart->normalSum;
	owner->centroidSum += chart->centroidSum;
	updateProxy(owner);
}

void AtlasBuilder::mergeCharts()
{
	std::vector<float> sharedBoundaryLengths;
	const uint32_t chartCount = chartArray.size();
	for (int c = chartCount - 1; c >= 0; c--) {
		sharedBoundaryLengths.clear();
		sharedBoundaryLengths.resize(chartCount, 0.0f);
		ChartBuildData *chart = chartArray[c];
		float externalBoundary = 0.0f;
		const uint32_t faceCount = chart->faces.size();
		for (uint32_t i = 0; i < faceCount; i++) {
			uint32_t f = chart->faces[i];
			const HalfEdge::Face *face = mesh->faceAt(f);
			for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance()) {
				const HalfEdge::Edge *edge = it.current();
				//float l = edge->length();
				float l = edgeLengths[edge->id / 2];
				if (edge->isBoundary()) {
					externalBoundary += l;
				} else {
					uint32_t neighborFace = edge->pair->face->id;
					uint32_t neighborChart = faceChartArray[neighborFace];
					if (neighborChart != c) {
						if ((edge->isSeam() && (isNormalSeam(edge) || isTextureSeam(edge))) || neighborChart == -2) {
							externalBoundary += l;
						} else {
							sharedBoundaryLengths[neighborChart] += l;
						}
					}
				}
			}
		}
		for (int cc = chartCount - 1; cc >= 0; cc--) {
			if (cc == c)
				continue;
			ChartBuildData *chart2 = chartArray[cc];
			if (chart2 == NULL)
				continue;
			if (sharedBoundaryLengths[cc] > 0.8 * std::max(0.0f, chart->boundaryLength - externalBoundary)) {
				// Try to avoid degenerate configurations.
				if (chart2->boundaryLength > sharedBoundaryLengths[cc]) {
					if (dot(chart2->planeNormal, chart->planeNormal) > -0.25) {
						mergeChart(chart2, chart, sharedBoundaryLengths[cc]);
						delete chart;
						chartArray[c] = NULL;
						break;
					}
				}
			}
			if (sharedBoundaryLengths[cc] > 0.20 * std::max(0.0f, chart->boundaryLength - externalBoundary)) {
				// Compare proxies.
				if (dot(chart2->planeNormal, chart->planeNormal) > 0) {
					mergeChart(chart2, chart, sharedBoundaryLengths[cc]);
					delete chart;
					chartArray[c] = NULL;
					break;
				}
			}
		}
	}
	// Remove deleted charts.
	for (int c = 0; c < int32_t(chartArray.size()); /*do not increment if removed*/) {
		if (chartArray[c] == NULL) {
			chartArray.erase(chartArray.begin() + c);
			// Update faceChartArray.
			const uint32_t faceCount = faceChartArray.size();
			for (uint32_t i = 0; i < faceCount; i++) {
				nvDebugCheck (faceChartArray[i] != -1);
				nvDebugCheck (faceChartArray[i] != c);
				nvDebugCheck (faceChartArray[i] <= int32_t(chartArray.size()));
				if (faceChartArray[i] > c) {
					faceChartArray[i]--;
				}
			}
		} else {
			chartArray[c]->id = c;
			c++;
		}
	}
}

const std::vector<uint32_t> &AtlasBuilder::chartFaces(uint32_t i) const
{
	return chartArray[i]->faces;
}
