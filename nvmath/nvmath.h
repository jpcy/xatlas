// This code is in the public domain -- castanyo@yahoo.es

#pragma once
#ifndef NV_MATH_H
#define NV_MATH_H

#include "nvcore/nvcore.h"

#include <math.h>
#include <float.h>  // finite, isnan

#ifndef PI
#define PI                  float(3.1415926535897932384626433833)
#endif

#define NV_EPSILON          (0.0001f)
#define NV_NORMAL_EPSILON   (0.001f)

namespace nv
{
// Robust floating point comparisons:
// http://realtimecollisiondetection.net/blog/?p=89
inline bool equal(const float f0, const float f1, const float epsilon = NV_EPSILON)
{
	//return fabs(f0-f1) <= epsilon;
	return fabs(f0 - f1) <= epsilon * max3(1.0f, fabsf(f0), fabsf(f1));
}

inline bool isZero(const float f, const float epsilon = NV_EPSILON)
{
	return fabs(f) <= epsilon;
}

inline float lerp(float f0, float f1, float t)
{
	const float s = 1.0f - t;
	return f0 * s + f1 * t;
}

inline float square(float f)
{
	return f * f;
}
inline int square(int i)
{
	return i * i;
}

} // nv

#endif // NV_MATH_H
