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

}
}
