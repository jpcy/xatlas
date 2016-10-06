// Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>

#pragma once
#ifndef NV_MESH_ATLAS_H
#define NV_MESH_ATLAS_H

#include <memory>
#include "xatlas.h"


namespace nv
{
namespace HalfEdge
{
class Mesh;
}

namespace param {
class Chart;
class MeshCharts;
class VertexMap;

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

}
} // nv namespace

#endif // NV_MESH_ATLAS_H
