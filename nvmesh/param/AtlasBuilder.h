// This code is in the public domain -- castano@gmail.com

#pragma once
#ifndef NV_MESH_ATLASBUILDER_H
#define NV_MESH_ATLASBUILDER_H

#include "Atlas.h"

#include "xatlas.h"
#include "xatlas.h"


namespace nv
{
namespace HalfEdge
{
class Mesh;
}

struct ChartBuildData;

struct AtlasBuilder {
	AtlasBuilder(const HalfEdge::Mesh *m);
	~AtlasBuilder();

	void markUnchartedFaces(const std::vector<uint32_t> &unchartedFaces);

	void computeShortestPaths();

	void placeSeeds(float threshold, uint32_t maxSeedCount);
	void createRandomChart(float threshold);

	void addFaceToChart(ChartBuildData *chart, uint32_t f, bool recomputeProxy = false);

	bool growCharts(float threshold, uint32_t faceCount);
	bool growChart(ChartBuildData *chart, float threshold, uint32_t faceCount);

	void resetCharts();

	void updateCandidates(ChartBuildData *chart, uint32_t face);

	void updateProxies();
	void updateProxy(ChartBuildData *chart);

	bool relocateSeeds();
	bool relocateSeed(ChartBuildData *chart);

	void updatePriorities(ChartBuildData *chart);

	float evaluatePriority(ChartBuildData *chart, uint32_t face);
	float evaluateProxyFitMetric(ChartBuildData *chart, uint32_t face);
	float evaluateDistanceToBoundary(ChartBuildData *chart, uint32_t face);
	float evaluateDistanceToSeed(ChartBuildData *chart, uint32_t face);
	float evaluateRoundnessMetric(ChartBuildData *chart, uint32_t face, float newBoundaryLength, float newChartArea);
	float evaluateStraightnessMetric(ChartBuildData *chart, uint32_t face);

	float evaluateNormalSeamMetric(ChartBuildData *chart, uint32_t f);
	float evaluateTextureSeamMetric(ChartBuildData *chart, uint32_t f);
	float evaluateSeamMetric(ChartBuildData *chart, uint32_t f);

	float evaluateChartArea(ChartBuildData *chart, uint32_t f);
	float evaluateBoundaryLength(ChartBuildData *chart, uint32_t f);
	Vector3 evaluateChartNormalSum(ChartBuildData *chart, uint32_t f);
	Vector3 evaluateChartCentroidSum(ChartBuildData *chart, uint32_t f);

	Vector3 computeChartCentroid(const ChartBuildData *chart);


	void fillHoles(float threshold);
	void mergeCharts();

	// @@ Cleanup.
	struct Candidate {
		uint32_t face;
		ChartBuildData *chart;
		float metric;
	};

	const Candidate &getBestCandidate() const;
	void removeCandidate(uint32_t f);
	void updateCandidate(ChartBuildData *chart, uint32_t f, float metric);

	void mergeChart(ChartBuildData *owner, ChartBuildData *chart, float sharedBoundaryLength);


	uint32_t chartCount() const
	{
		return chartArray.size();
	}
	const std::vector<uint32_t> &chartFaces(uint32_t i) const;

	const HalfEdge::Mesh *mesh;
	uint32_t facesLeft;
	std::vector<int> faceChartArray;
	std::vector<ChartBuildData *> chartArray;
	std::vector<float> shortestPaths;

	std::vector<float> edgeLengths;
	std::vector<float> faceAreas;

	std::vector<Candidate> candidateArray; //
	std::vector<uint32_t> faceCandidateArray; // Map face index to candidate index.

	MTRand rand;

	SegmentationSettings settings;
};

} // nv namespace

#endif // NV_MESH_ATLASBUILDER_H
