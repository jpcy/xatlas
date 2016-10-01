// This code is in the public domain -- castanyo@yahoo.es

#pragma once
#ifndef NV_MATH_MATRIX_H
#define NV_MATH_MATRIX_H

#include "Vector.h"

// - Matrices are stored in memory in *column major* order.
// - Points are to be though of as column vectors.
// - Transformation of a point p by a matrix M is: p' = M * p

namespace nv
{
    enum identity_t { identity };
	
    // 4x4 matrix.
    class Matrix
    {
    public:
        typedef Matrix const & Arg;

        Matrix();
        explicit Matrix(float f);
        explicit Matrix(identity_t);
        Matrix(const Matrix & m);
        Matrix(Vector4::Arg v0, Vector4::Arg v1, Vector4::Arg v2, Vector4::Arg v3);
        //explicit Matrix(const float m[]);	// m is assumed to contain 16 elements

        float data(uint idx) const;
        float & data(uint idx);
        float get(uint row, uint col) const;
        float operator()(uint row, uint col) const;
        float & operator()(uint row, uint col);
        const float * ptr() const;

        Vector4 row(uint i) const;
        Vector4 column(uint i) const;

        void zero();
        void identity();

        void scale(float s);
        void scale(Vector3::Arg s);
        void translate(Vector3::Arg t);
        void rotate(float theta, float v0, float v1, float v2);
        float determinant() const;

        void operator+=(const Matrix & m);
        void operator-=(const Matrix & m);

        void apply(Matrix::Arg m);

    private:
        float m_data[16];
    };
} // nv namespace

#endif // NV_MATH_MATRIX_H
