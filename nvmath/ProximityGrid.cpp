#include "ProximityGrid.h"
#include "Morton.h"
#include "Box.h"

using namespace nv;

ProximityGrid::ProximityGrid()
{
}

void ProximityGrid::init(const Box &box, uint count)
{
	cellArray.clear();
	// Determine grid size.
	float cellWidth;
	Vector3 diagonal = box.extents() * 2.f;
	float volume = box.volume();
	if (equal(volume, 0)) {
		// Degenerate box, treat like a quad.
		Vector2 quad;
		if (diagonal.x < diagonal.y && diagonal.x < diagonal.z) {
			quad.x = diagonal.y;
			quad.y = diagonal.z;
		} else if (diagonal.y < diagonal.x && diagonal.y < diagonal.z) {
			quad.x = diagonal.x;
			quad.y = diagonal.z;
		} else {
			quad.x = diagonal.x;
			quad.y = diagonal.y;
		}
		float cellArea = quad.x * quad.y / count;
		cellWidth = sqrtf(cellArea); // pow(cellArea, 1.0f / 2.0f);
	} else {
		// Ideally we want one cell per point.
		float cellVolume = volume / count;
		cellWidth = powf(cellVolume, 1.0f / 3.0f);
	}
	nvDebugCheck(cellWidth != 0);
	sx = max(1, ftoi_ceil(diagonal.x / cellWidth));
	sy = max(1, ftoi_ceil(diagonal.y / cellWidth));
	sz = max(1, ftoi_ceil(diagonal.z / cellWidth));
	invCellSize.x = float(sx) / diagonal.x;
	invCellSize.y = float(sy) / diagonal.y;
	invCellSize.z = float(sz) / diagonal.z;
	cellArray.resize(sx * sy * sz);
	corner = box.minCorner; // @@ Align grid better?
}

// Gather all points inside the given sphere.
// Radius is assumed to be small, so we don't bother culling the cells.
void ProximityGrid::gather(const Vector3 &position, float radius, std::vector<uint> &indexArray)
{
	int x0 = index_x(position.x - radius);
	int x1 = index_x(position.x + radius);
	int y0 = index_y(position.y - radius);
	int y1 = index_y(position.y + radius);
	int z0 = index_z(position.z - radius);
	int z1 = index_z(position.z + radius);
	for (int z = z0; z <= z1; z++) {
		for (int y = y0; y <= y1; y++) {
			for (int x = x0; x <= x1; x++) {
				int idx = index(x, y, z);
				indexArray.insert(indexArray.begin(), cellArray[idx].indexArray.begin(), cellArray[idx].indexArray.end());
			}
		}
	}
}


uint32_t ProximityGrid::mortonCount() const
{
	uint64_t s = uint64_t(max3(sx, sy, sz));
	s = nextPowerOfTwo(s);
	if (s > 1024) {
		return uint32_t(s * s * min3(sx, sy, sz));
	}
	return uint32_t(s * s * s);
}

int ProximityGrid::mortonIndex(uint32_t code) const
{
	uint32_t x, y, z;
	uint s = uint32_t(max3(sx, sy, sz));
	if (s > 1024) {
		// Use layered two-dimensional morton order.
		s = nextPowerOfTwo(s);
		uint layer = code / (s * s);
		code = code % (s * s);
		uint layer_count = uint32_t(min3(sx, sy, sz));
		if (sx == layer_count) {
			x = layer;
			y = decodeMorton2X(code);
			z = decodeMorton2Y(code);
		} else if (sy == layer_count) {
			x = decodeMorton2Y(code);
			y = layer;
			z = decodeMorton2X(code);
		} else { /*if (sz == layer_count)*/
			x = decodeMorton2X(code);
			y = decodeMorton2Y(code);
			z = layer;
		}
	} else {
		x = decodeMorton3X(code);
		y = decodeMorton3Y(code);
		z = decodeMorton3Z(code);
	}
	if (x >= uint32_t(sx) || y >= uint32_t(sy) || z >= uint32_t(sz)) {
		return -1;
	}
	return index(x, y, z);
}
