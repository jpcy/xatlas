#pragma once
#ifndef NV_MATH_PROXIMITYGRID_H
#define NV_MATH_PROXIMITYGRID_H

#include <vector>
#include "Vector.h"
#include "ftoi.h"

#include "nvcore/Array.h"


// A simple, dynamic proximity grid based on Jon's code.
// Instead of storing pointers here I store indices.

namespace nv
{

class Box;

struct Cell {
	std::vector<uint> indexArray;
};

struct ProximityGrid {
	ProximityGrid();

	void init(const Box &box, uint count);

	int index_x(float x) const;
	int index_y(float y) const;
	int index_z(float z) const;
	int index(int x, int y, int z) const;

	uint32 mortonCount() const;
	int mortonIndex(uint32 code) const;

	void add(const Vector3 &pos, uint key);

	void gather(const Vector3 &pos, float radius, std::vector<uint> &indices);

	Array<Cell> cellArray;

	Vector3 corner;
	Vector3 invCellSize;
	int sx, sy, sz;
};

inline int ProximityGrid::index_x(float x) const
{
	return clamp(ftoi_floor((x - corner.x) * invCellSize.x),  0, sx - 1);
}

inline int ProximityGrid::index_y(float y) const
{
	return clamp(ftoi_floor((y - corner.y) * invCellSize.y),  0, sy - 1);
}

inline int ProximityGrid::index_z(float z) const
{
	return clamp(ftoi_floor((z - corner.z) * invCellSize.z),  0, sz - 1);
}

inline int ProximityGrid::index(int x, int y, int z) const
{
	nvDebugCheck(x >= 0 && x < sx);
	nvDebugCheck(y >= 0 && y < sy);
	nvDebugCheck(z >= 0 && z < sz);
	int idx = (z * sy + y) * sx + x;
	nvDebugCheck(idx >= 0 && uint(idx) < cellArray.count());
	return idx;
}

inline void ProximityGrid::add(const Vector3 &pos, uint key)
{
	int x = index_x(pos.x);
	int y = index_y(pos.y);
	int z = index_z(pos.z);
	uint idx = index(x, y, z);
	cellArray[idx].indexArray.push_back(key);
}

} // nv namespace

#endif // NV_MATH_PROXIMITYGRID_H
