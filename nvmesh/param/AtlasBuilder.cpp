// This code is in the public domain -- castano@gmail.com

#include "xatlas.h"

//#include "nvcore/IntroSort.h"

#include <algorithm> // std::sort

#include <float.h> // FLT_MAX
#include <limits.h> // UINT_MAX

namespace nv {
namespace param {

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
		Vector3 faceCentroid = face->triangleCenter();
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
	Vector3 faceNormal = face->triangleNormal();
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
	return length(seed->triangleCenter() - face->triangleCenter());
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
		if (edge->isNormalSeam()) {
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
		if (edge->isTextureSeam()) {
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
	return chart->normalSum + face->triangleNormalAreaScaled();
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
		centroid += face->triangleCenter();
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
						if ((edge->isSeam() && (edge->isNormalSeam() || edge->isTextureSeam())) || neighborChart == -2) {
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

}
}
