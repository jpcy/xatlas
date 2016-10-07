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

}
} // nv namespace

#endif // NV_MESH_ATLAS_H
