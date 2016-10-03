// Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>

#pragma once
#ifndef NV_MESH_ATLAS_H
#define NV_MESH_ATLAS_H

#include <memory>
#include "nvmath/Vector.h"
#include "nvmesh/halfedge/Mesh.h"


namespace nv
{
namespace HalfEdge
{
class Mesh;
}

class Chart;
class MeshCharts;
class VertexMap;

struct SegmentationSettings {
	SegmentationSettings();

	float maxChartArea;
	float maxBoundaryLength;

	float proxyFitMetricWeight;
	float roundnessMetricWeight;
	float straightnessMetricWeight;
	float normalSeamMetricWeight;
	float textureSeamMetricWeight;
};


/// An atlas is a set of charts.
class Atlas
{
public:

	Atlas();
	~Atlas();

	uint32_t meshCount() const
	{
		return m_meshChartsArray.size();
	}
	const MeshCharts *meshAt(uint32_t i) const
	{
		return m_meshChartsArray[i];
	}
	MeshCharts *meshAt(uint32_t i)
	{
		return m_meshChartsArray[i];
	}

	uint32_t chartCount() const;
	const Chart *chartAt(uint32_t i) const;
	Chart *chartAt(uint32_t i);

	// Add mesh charts and takes ownership.
	void addMeshCharts(MeshCharts *meshCharts);

	void extractCharts(const HalfEdge::Mesh *mesh);
	void computeCharts(const HalfEdge::Mesh *mesh, const SegmentationSettings &settings, const std::vector<uint32_t> &unchartedMaterialArray);


	// Compute a trivial seamless texture similar to ZBrush.
	//bool computeSeamlessTextureAtlas(bool groupFaces = true, bool scaleTiles = false, uint32_t w = 1024, uint32_t h = 1024);

	void parameterizeCharts();

	// Pack charts in the smallest possible rectangle.
	float packCharts(int quality, float texelArea, bool blockAlign, bool conservative);

private:

	std::vector<MeshCharts *> m_meshChartsArray;

};


// Set of charts corresponding to a single mesh.
class MeshCharts
{
public:
	MeshCharts(const HalfEdge::Mesh *mesh);
	~MeshCharts();

	uint32_t chartCount() const
	{
		return m_chartArray.size();
	}
	uint32_t vertexCount () const
	{
		return m_totalVertexCount;
	}

	const Chart *chartAt(uint32_t i) const
	{
		return m_chartArray[i];
	}
	Chart *chartAt(uint32_t i)
	{
		return m_chartArray[i];
	}

	void computeVertexMap(const std::vector<uint32_t> &unchartedMaterialArray);

	// Extract the charts of the input mesh.
	void extractCharts();

	// Compute charts using a simple segmentation algorithm.
	void computeCharts(const SegmentationSettings &settings, const std::vector<uint32_t> &unchartedMaterialArray);

	void parameterizeCharts();

	uint32_t faceChartAt(uint32_t i) const
	{
		return m_faceChart[i];
	}
	uint32_t faceIndexWithinChartAt(uint32_t i) const
	{
		return m_faceIndex[i];
	}

	uint32_t vertexCountBeforeChartAt(uint32_t i) const
	{
		return m_chartVertexCountPrefixSum[i];
	}

private:

	const HalfEdge::Mesh *m_mesh;

	std::vector<Chart *> m_chartArray;

	std::vector<uint32_t> m_chartVertexCountPrefixSum;
	uint32_t m_totalVertexCount;

	std::vector<uint32_t> m_faceChart; // the chart of every face of the input mesh.
	std::vector<uint32_t> m_faceIndex; // the index within the chart for every face of the input mesh.
};


/// A chart is a connected set of faces with a certain topology (usually a disk).
class Chart
{
public:

	Chart();

	void build(const HalfEdge::Mesh *originalMesh, const std::vector<uint32_t> &faceArray);
	void buildVertexMap(const HalfEdge::Mesh *originalMesh, const std::vector<uint32_t> &unchartedMaterialArray);

	bool closeHoles();

	bool isDisk() const
	{
		return m_isDisk;
	}
	bool isVertexMapped() const
	{
		return m_isVertexMapped;
	}

	uint32_t vertexCount() const
	{
		return m_chartMesh->vertexCount();
	}
	uint32_t colocalVertexCount() const
	{
		return m_unifiedMesh->vertexCount();
	}

	uint32_t faceCount() const
	{
		return m_faceArray.size();
	}
	uint32_t faceAt(uint32_t i) const
	{
		return m_faceArray[i];
	}

	const HalfEdge::Mesh *chartMesh() const
	{
		return m_chartMesh.get();
	}
	HalfEdge::Mesh *chartMesh()
	{
		return m_chartMesh.get();
	}
	const HalfEdge::Mesh *unifiedMesh() const
	{
		return m_unifiedMesh.get();
	}
	HalfEdge::Mesh *unifiedMesh()
	{
		return m_unifiedMesh.get();
	}

	//uint32_t vertexIndex(uint32_t i) const { return m_vertexIndexArray[i]; }

	uint32_t mapChartVertexToOriginalVertex(uint32_t i) const
	{
		return m_chartToOriginalMap[i];
	}
	uint32_t mapChartVertexToUnifiedVertex(uint32_t i) const
	{
		return m_chartToUnifiedMap[i];
	}

	const std::vector<uint32_t> &faceArray() const
	{
		return m_faceArray;
	}

	void transferParameterization();

	float computeSurfaceArea() const;
	float computeParametricArea() const;
	Vector2 computeParametricBounds() const;


	float scale = 1.0f;
	uint32_t vertexMapWidth;
	uint32_t vertexMapHeight;

private:

	bool closeLoop(uint32_t start, const std::vector<HalfEdge::Edge *> &loop);

	// Chart mesh.
	std::auto_ptr<HalfEdge::Mesh> m_chartMesh;
	std::auto_ptr<HalfEdge::Mesh> m_unifiedMesh;

	bool m_isDisk;
	bool m_isVertexMapped;

	// List of faces of the original mesh that belong to this chart.
	std::vector<uint32_t> m_faceArray;

	// Map vertices of the chart mesh to vertices of the original mesh.
	std::vector<uint32_t> m_chartToOriginalMap;

	std::vector<uint32_t> m_chartToUnifiedMap;
};

} // nv namespace

#endif // NV_MESH_ATLAS_H
