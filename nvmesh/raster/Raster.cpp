// This code is in the public domain -- castanyo@yahoo.es

/** @file Raster.cpp
 * @brief Triangle rasterization library using affine interpolation. Not
 * specially optimized, but enough for my purposes.
**/

#include "xatlas.h"

namespace nv {
namespace raster {
/// Process the given triangle.
bool drawTriangle(Mode mode, Vector2::Arg extents, bool enableScissors, const Vector2 v[3], SamplingCallback cb, void *param)
{
	Triangle tri(v[0], v[1], v[2], Vector3(1, 0, 0), Vector3(0, 1, 0), Vector3(0, 0, 1));
	// @@ It would be nice to have a conservative drawing mode that enlarges the triangle extents by one texel and is able to handle degenerate triangles.
	// @@ Maybe the simplest thing to do would be raster triangle edges.
	if (tri.valid) {
		if (mode == Mode_Antialiased) {
			return tri.drawAA(extents, enableScissors, cb, param);
		}
		if (mode == Mode_Nearest) {
			return tri.draw(extents, enableScissors, cb, param);
		}
	}
	return true;
}

/// Process the given quad.
bool drawQuad(Mode mode, Vector2::Arg extents, bool enableScissors, const Vector2 v[4], SamplingCallback cb, void *param)
{
	bool sign0 = triangleArea2(v[0], v[1], v[2]) > 0.0f;
	bool sign1 = triangleArea2(v[0], v[2], v[3]) > 0.0f;
	// Divide the quad into two non overlapping triangles.
	if (sign0 == sign1) {
		Triangle tri0(v[0], v[1], v[2], Vector3(0, 0, 0), Vector3(1, 0, 0), Vector3(1, 1, 0));
		Triangle tri1(v[0], v[2], v[3], Vector3(0, 0, 0), Vector3(1, 1, 0), Vector3(0, 1, 0));
		if (tri0.valid && tri1.valid) {
			if (mode == Mode_Antialiased) {
				return tri0.drawAA(extents, enableScissors, cb, param) && tri1.drawAA(extents, enableScissors, cb, param);
			} else {
				return tri0.draw(extents, enableScissors, cb, param) && tri1.draw(extents, enableScissors, cb, param);
			}
		}
	} else {
		Triangle tri0(v[0], v[1], v[3], Vector3(0, 0, 0), Vector3(1, 0, 0), Vector3(0, 1, 0));
		Triangle tri1(v[1], v[2], v[3], Vector3(1, 0, 0), Vector3(1, 1, 0), Vector3(0, 1, 0));
		if (tri0.valid && tri1.valid) {
			if (mode == Mode_Antialiased) {
				return tri0.drawAA(extents, enableScissors, cb, param) && tri1.drawAA(extents, enableScissors, cb, param);
			} else {
				return tri0.draw(extents, enableScissors, cb, param) && tri1.draw(extents, enableScissors, cb, param);
			}
		}
	}
	return true;
}

} // namespace raster
} // namespace nv
