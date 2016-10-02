// This code is in the public domain -- Ignacio Castaño <castano@gmail.com>

#pragma once
#ifndef NV_MATH_BASIS_H
#define NV_MATH_BASIS_H

#include "nvmath.h"
#include "Vector.inl"

namespace nv
{

    /// Basis class to compute tangent space basis, ortogonalizations and to
    /// transform vectors from one space to another.
    class Basis
    {
    public:

        /// Create a null basis.
        Basis() : tangent(0, 0, 0), bitangent(0, 0, 0), normal(0, 0, 0) {}

        void buildFrameForDirection(Vector3::Arg d, float angle = 0);

        Vector3 tangent;
        Vector3 bitangent;
        Vector3 normal;
    };

} // nv namespace

#endif // NV_MATH_BASIS_H
