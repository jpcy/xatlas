// This code is in the public domain -- Ignacio Castaño <castano@gmail.com>

#pragma once
#ifndef NV_MATH_FITTING_H
#define NV_MATH_FITTING_H

#include "Vector.h"

namespace nv
{
    namespace Fit
    {
        Vector3 computeCentroid(int n, const Vector3 * points);

        Vector3 computeCovariance(int n, const Vector3 * points, float * covariance);

        bool isPlanar(int n, const Vector3 * points, float epsilon = NV_EPSILON);

        bool eigenSolveSymmetric3(const float matrix[6], float eigenValues[3], Vector3 eigenVectors[3]);
    }

} // nv namespace

#endif // NV_MATH_FITTING_H
