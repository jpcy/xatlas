// This code is in the public domain -- castanyo@yahoo.es

/** @file Raster.cpp
 * @brief Triangle rasterization library using affine interpolation. Not
 * specially optimized, but enough for my purposes.
**/

#include "Raster.h"

#include "nvcore/Utils.h" // min, max

#include "nvmath/Vector.h"
#include "nvmath/ftoi.h"


#define RA_EPSILON		0.00001f

using namespace nv;
using namespace nv::Raster;

namespace
{
class ClippedTriangle
{
public:
	ClippedTriangle(Vector2::Arg a, Vector2::Arg b, Vector2::Arg c)
	{
		m_numVertices = 3;
		m_activeVertexBuffer = 0;
		m_verticesA[0] = a;
		m_verticesA[1] = b;
		m_verticesA[2] = c;
		m_vertexBuffers[0] = m_verticesA;
		m_vertexBuffers[1] = m_verticesB;
	}

	uint32_t vertexCount()
	{
		return m_numVertices;
	}

	const Vector2 *vertices()
	{
		return m_vertexBuffers[m_activeVertexBuffer];
	}

	inline void clipHorizontalPlane(float offset, float clipdirection)
	{
		Vector2 *v  = m_vertexBuffers[m_activeVertexBuffer];
		m_activeVertexBuffer ^= 1;
		Vector2 *v2 = m_vertexBuffers[m_activeVertexBuffer];
		v[m_numVertices] = v[0];
		float dy2,   dy1 = offset - v[0].y;
		int   dy2in, dy1in = clipdirection * dy1 >= 0;
		uint32_t  p = 0;
		for (uint32_t k = 0; k < m_numVertices; k++) {
			dy2   = offset - v[k + 1].y;
			dy2in = clipdirection * dy2 >= 0;
			if (dy1in) v2[p++] = v[k];
			if ( dy1in + dy2in == 1 ) { // not both in/out
				float dx = v[k + 1].x - v[k].x;
				float dy = v[k + 1].y - v[k].y;
				v2[p++] = Vector2(v[k].x + dy1 * (dx / dy), offset);
			}
			dy1 = dy2;
			dy1in = dy2in;
		}
		m_numVertices = p;
		//for (uint32_t k=0; k<m_numVertices; k++) printf("(%f, %f)\n", v2[k].x, v2[k].y); printf("\n");
	}

	inline void clipVerticalPlane(float offset, float clipdirection )
	{
		Vector2 *v  = m_vertexBuffers[m_activeVertexBuffer];
		m_activeVertexBuffer ^= 1;
		Vector2 *v2 = m_vertexBuffers[m_activeVertexBuffer];
		v[m_numVertices] = v[0];
		float dx2,   dx1   = offset - v[0].x;
		int   dx2in, dx1in = clipdirection * dx1 >= 0;
		uint32_t  p = 0;
		for (uint32_t k = 0; k < m_numVertices; k++) {
			dx2 = offset - v[k + 1].x;
			dx2in = clipdirection * dx2 >= 0;
			if (dx1in) v2[p++] = v[k];
			if ( dx1in + dx2in == 1 ) { // not both in/out
				float dx = v[k + 1].x - v[k].x;
				float dy = v[k + 1].y - v[k].y;
				v2[p++] = Vector2(offset, v[k].y + dx1 * (dy / dx));
			}
			dx1 = dx2;
			dx1in = dx2in;
		}
		m_numVertices = p;
		//for (uint32_t k=0; k<m_numVertices; k++) printf("(%f, %f)\n", v2[k].x, v2[k].y); printf("\n");
	}

	void computeAreaCentroid()
	{
		Vector2 *v  = m_vertexBuffers[m_activeVertexBuffer];
		v[m_numVertices] = v[0];
		m_area = 0;
		float centroidx = 0, centroidy = 0;
		for (uint32_t k = 0; k < m_numVertices; k++) {
			// http://local.wasp.uwa.edu.au/~pbourke/geometry/polyarea/
			float f = v[k].x * v[k + 1].y - v[k + 1].x * v[k].y;
			m_area += f;
			centroidx += f * (v[k].x + v[k + 1].x);
			centroidy += f * (v[k].y + v[k + 1].y);
		}
		m_area = 0.5f * fabsf(m_area);
		if (m_area == 0) {
			m_centroid = Vector2(0.0f);
		} else {
			m_centroid = Vector2(centroidx / (6 * m_area), centroidy / (6 * m_area));
		}
	}

	void clipAABox(float x0, float y0, float x1, float y1)
	{
		clipVerticalPlane  ( x0, -1);
		clipHorizontalPlane( y0, -1);
		clipVerticalPlane  ( x1,  1);
		clipHorizontalPlane( y1,  1);
		computeAreaCentroid();
	}

	Vector2 centroid()
	{
		return m_centroid;
	}

	float area()
	{
		return m_area;
	}

private:
	Vector2 m_verticesA[7 + 1];
	Vector2 m_verticesB[7 + 1];
	Vector2 *m_vertexBuffers[2];
	uint32_t    m_numVertices;
	uint32_t    m_activeVertexBuffer;
	float   m_area;
	Vector2 m_centroid;
};

static inline float delta(float bot, float top, float ih)
{
	return (bot - top) * ih;
}

static inline Vector2 delta(Vector2::Arg bot, Vector2::Arg top, float ih)
{
	return (bot - top) * ih;
}

static inline Vector3 delta(Vector3::Arg bot, Vector3::Arg top, float ih)
{
	return (bot - top) * ih;
}

// @@ The implementation in nvmath.h should be equivalent.
static inline int iround(float f)
{
	// @@ Optimize this.
	return int(floorf(f + 0.5f));
	//return int(round(f));
	//return int(f);
}

/// A triangle vertex.
struct Vertex {
	Vector2 pos;	// Position.
	Vector3 tex;	// Texcoord. (Barycentric coordinate)
};


/// A triangle for rasterization.
struct Triangle {
	Triangle(Vector2::Arg v0, Vector2::Arg v1, Vector2::Arg v2, Vector3::Arg t0, Vector3::Arg t1, Vector3::Arg t2);

	bool computeDeltas();

	bool draw(const Vector2 &extents, bool enableScissors, SamplingCallback cb, void *param);
	bool drawAA(const Vector2 &extents, bool enableScissors, SamplingCallback cb, void *param);
	bool drawC(const Vector2 &extents, bool enableScissors, SamplingCallback cb, void *param);
	void flipBackface();
	void computeUnitInwardNormals();

	// Vertices.
	Vector2 v1, v2, v3;
	Vector2 n1, n2, n3; // unit inward normals
	Vector3 t1, t2, t3;

	// Deltas.
	Vector3 dx, dy;

	float sign;
	bool valid;
};


/// Triangle ctor.
Triangle::Triangle(Vector2::Arg v0, Vector2::Arg v1, Vector2::Arg v2,
                   Vector3::Arg t0, Vector3::Arg t1, Vector3::Arg t2)
{
	// Init vertices.
	this->v1 = v0;
	this->v2 = v2;
	this->v3 = v1;
	// Set barycentric coordinates.
	this->t1 = t0;
	this->t2 = t2;
	this->t3 = t1;
	// make sure every triangle is front facing.
	flipBackface();
	// Compute deltas.
	valid = computeDeltas();
	computeUnitInwardNormals();
}


/// Compute texture space deltas.
/// This method takes two edge vectors that form a basis, determines the
/// coordinates of the canonic vectors in that basis, and computes the
/// texture gradient that corresponds to those vectors.
bool Triangle::computeDeltas()
{
	Vector2 e0 = v3 - v1;
	Vector2 e1 = v2 - v1;
	Vector3 de0 = t3 - t1;
	Vector3 de1 = t2 - t1;
	float denom = 1.0f / (e0.y * e1.x - e1.y * e0.x);
	if (!std::isfinite(denom)) {
		return false;
	}
	float lambda1 = - e1.y * denom;
	float lambda2 = e0.y * denom;
	float lambda3 = e1.x * denom;
	float lambda4 = - e0.x * denom;
	dx = de0 * lambda1 + de1 * lambda2;
	dy = de0 * lambda3 + de1 * lambda4;
	return true;
}

// compute unit inward normals for each edge.
void Triangle::computeUnitInwardNormals()
{
	n1 = v1 - v2;
	n1 = Vector2(-n1.y, n1.x);
	n1 = n1 * (1.0f / sqrtf(n1.x * n1.x + n1.y * n1.y));
	n2 = v2 - v3;
	n2 = Vector2(-n2.y, n2.x);
	n2 = n2 * (1.0f / sqrtf(n2.x * n2.x + n2.y * n2.y));
	n3 = v3 - v1;
	n3 = Vector2(-n3.y, n3.x);
	n3 = n3 * (1.0f / sqrtf(n3.x * n3.x + n3.y * n3.y));
}

void Triangle::flipBackface()
{
	// check if triangle is backfacing, if so, swap two vertices
	if ( ((v3.x - v1.x) * (v2.y - v1.y) - (v3.y - v1.y) * (v2.x - v1.x)) < 0 ) {
		Vector2 hv = v1;
		v1 = v2;
		v2 = hv; // swap pos
		Vector3 ht = t1;
		t1 = t2;
		t2 = ht; // swap tex
	}
}

bool Triangle::draw(const Vector2 &extents, bool enableScissors, SamplingCallback cb, void *param)
{
	// 28.4 fixed-point coordinates
	const int Y1 = iround(16.0f * v1.y);
	const int Y2 = iround(16.0f * v2.y);
	const int Y3 = iround(16.0f * v3.y);
	const int X1 = iround(16.0f * v1.x);
	const int X2 = iround(16.0f * v2.x);
	const int X3 = iround(16.0f * v3.x);
	// Deltas
	const int DX12 = X1 - X2;
	const int DX23 = X2 - X3;
	const int DX31 = X3 - X1;
	const int DY12 = Y1 - Y2;
	const int DY23 = Y2 - Y3;
	const int DY31 = Y3 - Y1;
	// Fixed-point deltas
	const int FDX12 = DX12 << 4;
	const int FDX23 = DX23 << 4;
	const int FDX31 = DX31 << 4;
	const int FDY12 = DY12 << 4;
	const int FDY23 = DY23 << 4;
	const int FDY31 = DY31 << 4;
	int minx, miny, maxx, maxy;
	if (enableScissors) {
		int frustumX0 =  0 << 4;
		int frustumY0 =  0 << 4;
		int frustumX1 =  (int)extents.x << 4;
		int frustumY1 =  (int)extents.y << 4;
		// Bounding rectangle
		minx = (nv::max(min3(X1, X2, X3), frustumX0) + 0xF) >> 4;
		miny = (nv::max(min3(Y1, Y2, Y3), frustumY0) + 0xF) >> 4;
		maxx = (nv::min(max3(X1, X2, X3), frustumX1) + 0xF) >> 4;
		maxy = (nv::min(max3(Y1, Y2, Y3), frustumY1) + 0xF) >> 4;
	} else {
		// Bounding rectangle
		minx = (min3(X1, X2, X3) + 0xF) >> 4;
		miny = (min3(Y1, Y2, Y3) + 0xF) >> 4;
		maxx = (max3(X1, X2, X3) + 0xF) >> 4;
		maxy = (max3(Y1, Y2, Y3) + 0xF) >> 4;
	}
	// Block size, standard 8x8 (must be power of two)
	const int q = 8;
	// @@ This won't work when minx,miny are negative. This code path is not used. Leaving as is for now.
	nvCheck(minx >= 0);
	nvCheck(miny >= 0);
	// Start in corner of 8x8 block
	minx &= ~(q - 1);
	miny &= ~(q - 1);
	// Half-edge constants
	int C1 = DY12 * X1 - DX12 * Y1;
	int C2 = DY23 * X2 - DX23 * Y2;
	int C3 = DY31 * X3 - DX31 * Y3;
	// Correct for fill convention
	if (DY12 < 0 || (DY12 == 0 && DX12 > 0)) C1++;
	if (DY23 < 0 || (DY23 == 0 && DX23 > 0)) C2++;
	if (DY31 < 0 || (DY31 == 0 && DX31 > 0)) C3++;
	// Loop through blocks
	for (int y = miny; y < maxy; y += q) {
		for (int x = minx; x < maxx; x += q) {
			// Corners of block
			int x0 = x << 4;
			int x1 = (x + q - 1) << 4;
			int y0 = y << 4;
			int y1 = (y + q - 1) << 4;
			// Evaluate half-space functions
			bool a00 = C1 + DX12 * y0 - DY12 * x0 > 0;
			bool a10 = C1 + DX12 * y0 - DY12 * x1 > 0;
			bool a01 = C1 + DX12 * y1 - DY12 * x0 > 0;
			bool a11 = C1 + DX12 * y1 - DY12 * x1 > 0;
			int a = (a00 << 0) | (a10 << 1) | (a01 << 2) | (a11 << 3);
			bool b00 = C2 + DX23 * y0 - DY23 * x0 > 0;
			bool b10 = C2 + DX23 * y0 - DY23 * x1 > 0;
			bool b01 = C2 + DX23 * y1 - DY23 * x0 > 0;
			bool b11 = C2 + DX23 * y1 - DY23 * x1 > 0;
			int b = (b00 << 0) | (b10 << 1) | (b01 << 2) | (b11 << 3);
			bool c00 = C3 + DX31 * y0 - DY31 * x0 > 0;
			bool c10 = C3 + DX31 * y0 - DY31 * x1 > 0;
			bool c01 = C3 + DX31 * y1 - DY31 * x0 > 0;
			bool c11 = C3 + DX31 * y1 - DY31 * x1 > 0;
			int c = (c00 << 0) | (c10 << 1) | (c01 << 2) | (c11 << 3);
			// Skip block when outside an edge
			if (a == 0x0 || b == 0x0 || c == 0x0) continue;
			// Accept whole block when totally covered
			if (a == 0xF && b == 0xF && c == 0xF) {
				Vector3 texRow = t1 + dy * (y0 - v1.y) + dx * (x0 - v1.x);
				for (int iy = y; iy < y + q; iy++) {
					Vector3 tex = texRow;
					for (int ix = x; ix < x + q; ix++) {
						//Vector3 tex = t1 + dx * (ix - v1.x) + dy * (iy - v1.y);
						if (!cb(param, ix, iy, tex, dx, dy, 1.0)) {
							// early out.
							return false;
						}
						tex += dx;
					}
					texRow += dy;
				}
			} else { // Partially covered block
				int CY1 = C1 + DX12 * y0 - DY12 * x0;
				int CY2 = C2 + DX23 * y0 - DY23 * x0;
				int CY3 = C3 + DX31 * y0 - DY31 * x0;
				Vector3 texRow = t1 + dy * (y0 - v1.y) + dx * (x0 - v1.x);
				for (int iy = y; iy < y + q; iy++) {
					int CX1 = CY1;
					int CX2 = CY2;
					int CX3 = CY3;
					Vector3 tex = texRow;
					for (int ix = x; ix < x + q; ix++) {
						if (CX1 > 0 && CX2 > 0 && CX3 > 0) {
							if (!cb(param, ix, iy, tex, dx, dy, 1.0)) {
								// early out.
								return false;
							}
						}
						CX1 -= FDY12;
						CX2 -= FDY23;
						CX3 -= FDY31;
						tex += dx;
					}
					CY1 += FDX12;
					CY2 += FDX23;
					CY3 += FDX31;
					texRow += dy;
				}
			}
		}
	}
	return true;
}


#define PX_INSIDE    1.0f/sqrt(2.0f)
#define PX_OUTSIDE  -1.0f/sqrt(2.0f)

#define BK_SIZE 8
#define BK_INSIDE   sqrt(BK_SIZE*BK_SIZE/2.0f)
#define BK_OUTSIDE -sqrt(BK_SIZE*BK_SIZE/2.0f)

// extents has to be multiple of BK_SIZE!!
bool Triangle::drawAA(const Vector2 &extents, bool enableScissors, SamplingCallback cb, void *param)
{
	float minx, miny, maxx, maxy;
	if (enableScissors) {
		// Bounding rectangle
		minx = floorf(max(min3(v1.x, v2.x, v3.x), 0.0f));
		miny = floorf(max(min3(v1.y, v2.y, v3.y), 0.0f));
		maxx = ceilf( min(max3(v1.x, v2.x, v3.x), extents.x - 1.0f));
		maxy = ceilf( min(max3(v1.y, v2.y, v3.y), extents.y - 1.0f));
	} else {
		// Bounding rectangle
		minx = floorf(min3(v1.x, v2.x, v3.x));
		miny = floorf(min3(v1.y, v2.y, v3.y));
		maxx = ceilf( max3(v1.x, v2.x, v3.x));
		maxy = ceilf( max3(v1.y, v2.y, v3.y));
	}
	// There's no reason to align the blocks to the viewport, instead we align them to the origin of the triangle bounds.
	minx = floorf(minx);
	miny = floorf(miny);
	//minx = (float)(((int)minx) & (~((int)BK_SIZE - 1))); // align to blocksize (we don't need to worry about blocks partially out of viewport)
	//miny = (float)(((int)miny) & (~((int)BK_SIZE - 1)));
	minx += 0.5;
	miny += 0.5; // sampling at texel centers!
	maxx += 0.5;
	maxy += 0.5;
	// Half-edge constants
	float C1 = n1.x * (-v1.x) + n1.y * (-v1.y);
	float C2 = n2.x * (-v2.x) + n2.y * (-v2.y);
	float C3 = n3.x * (-v3.x) + n3.y * (-v3.y);
	// Loop through blocks
	for (float y0 = miny; y0 <= maxy; y0 += BK_SIZE) {
		for (float x0 = minx; x0 <= maxx; x0 += BK_SIZE) {
			// Corners of block
			float xc = (x0 + (BK_SIZE - 1) / 2.0f);
			float yc = (y0 + (BK_SIZE - 1) / 2.0f);
			// Evaluate half-space functions
			float aC = C1 + n1.x * xc + n1.y * yc;
			float bC = C2 + n2.x * xc + n2.y * yc;
			float cC = C3 + n3.x * xc + n3.y * yc;
			// Skip block when outside an edge
			if ( (aC <= BK_OUTSIDE) || (bC <= BK_OUTSIDE) || (cC <= BK_OUTSIDE) ) continue;
			// Accept whole block when totally covered
			if ( (aC >= BK_INSIDE) && (bC >= BK_INSIDE) && (cC >= BK_INSIDE) ) {
				Vector3 texRow = t1 + dy * (y0 - v1.y) + dx * (x0 - v1.x);
				for (float y = y0; y < y0 + BK_SIZE; y++) {
					Vector3 tex = texRow;
					for (float x = x0; x < x0 + BK_SIZE; x++) {
						if (!cb(param, (int)x, (int)y, tex, dx, dy, 1.0f)) {
							return false;
						}
						tex += dx;
					}
					texRow += dy;
				}
			} else { // Partially covered block
				float CY1 = C1 + n1.x * x0 + n1.y * y0;
				float CY2 = C2 + n2.x * x0 + n2.y * y0;
				float CY3 = C3 + n3.x * x0 + n3.y * y0;
				Vector3 texRow = t1 + dy * (y0 - v1.y) + dx * (x0 - v1.x);
				for (float y = y0; y < y0 + BK_SIZE; y++) { // @@ This is not clipping to scissor rectangle correctly.
					float CX1 = CY1;
					float CX2 = CY2;
					float CX3 = CY3;
					Vector3 tex = texRow;
					for (float x = x0; x < x0 + BK_SIZE; x++) { // @@ This is not clipping to scissor rectangle correctly.
						if (CX1 >= PX_INSIDE && CX2 >= PX_INSIDE && CX3 >= PX_INSIDE) {
							// pixel completely covered
							Vector3 tex = t1 + dx * (x - v1.x) + dy * (y - v1.y);
							if (!cb(param, (int)x, (int)y, tex, dx, dy, 1.0f)) {
								return false;
							}
						} else if ((CX1 >= PX_OUTSIDE) && (CX2 >= PX_OUTSIDE) && (CX3 >= PX_OUTSIDE)) {
							// triangle partially covers pixel. do clipping.
							ClippedTriangle ct(v1 - Vector2(x, y), v2 - Vector2(x, y), v3 - Vector2(x, y));
							ct.clipAABox(-0.5, -0.5, 0.5, 0.5);
							Vector2 centroid = ct.centroid();
							float area = ct.area();
							if (area > 0.0f) {
								Vector3 texCent = tex - dx * centroid.x - dy * centroid.y;
								//nvCheck(texCent.x >= -0.1f && texCent.x <= 1.1f); // @@ Centroid is not very exact...
								//nvCheck(texCent.y >= -0.1f && texCent.y <= 1.1f);
								//nvCheck(texCent.z >= -0.1f && texCent.z <= 1.1f);
								//Vector3 texCent2 = t1 + dx * (x - v1.x) + dy * (y - v1.y);
								if (!cb(param, (int)x, (int)y, texCent, dx, dy, area)) {
									return false;
								}
							}
						}
						CX1 += n1.x;
						CX2 += n2.x;
						CX3 += n3.x;
						tex += dx;
					}
					CY1 += n1.y;
					CY2 += n2.y;
					CY3 += n3.y;
					texRow += dy;
				}
			}
		}
	}
	return true;
}

} // namespace


/// Process the given triangle.
bool nv::Raster::drawTriangle(Mode mode, Vector2::Arg extents, bool enableScissors, const Vector2 v[3], SamplingCallback cb, void *param)
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

inline static float triangleArea(Vector2::Arg v1, Vector2::Arg v2, Vector2::Arg v3)
{
	return 0.5f * (v3.x * v1.y + v1.x * v2.y + v2.x * v3.y - v2.x * v1.y - v3.x * v2.y - v1.x * v3.y);
}

/// Process the given quad.
bool nv::Raster::drawQuad(Mode mode, Vector2::Arg extents, bool enableScissors, const Vector2 v[4], SamplingCallback cb, void *param)
{
	bool sign0 = triangleArea(v[0], v[1], v[2]) > 0.0f;
	bool sign1 = triangleArea(v[0], v[2], v[3]) > 0.0f;
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
