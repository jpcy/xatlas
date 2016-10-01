// This code is in the public domain -- castanyo@yahoo.es

#pragma once
#ifndef NV_MATH_MATRIX_INL
#define NV_MATH_MATRIX_INL

#include "Matrix.h"

namespace nv
{
    inline Matrix::Matrix()
    {
    }

    inline Matrix::Matrix(float f)
    {
        for(int i = 0; i < 16; i++) {
            m_data[i] = 0.0f;
        }
    }

    inline Matrix::Matrix(identity_t)
    {
        for(int i = 0; i < 4; i++) {
            for(int j = 0; j < 4; j++) {
                m_data[4*j+i] = (i == j) ? 1.0f : 0.0f;
            }
        }
    }

    inline Matrix::Matrix(const Matrix & m)
    {
        for(int i = 0; i < 16; i++) {
            m_data[i] = m.m_data[i];
        }
    }

    inline Matrix::Matrix(Vector4::Arg v0, Vector4::Arg v1, Vector4::Arg v2, Vector4::Arg v3)
    {
        m_data[ 0] = v0.x; m_data[ 1] = v0.y; m_data[ 2] = v0.z; m_data[ 3] = v0.w;
        m_data[ 4] = v1.x; m_data[ 5] = v1.y; m_data[ 6] = v1.z; m_data[ 7] = v1.w;
        m_data[ 8] = v2.x; m_data[ 9] = v2.y; m_data[10] = v2.z; m_data[11] = v2.w;
        m_data[12] = v3.x; m_data[13] = v3.y; m_data[14] = v3.z; m_data[15] = v3.w;
    }

    // Accessors
    inline float Matrix::data(uint idx) const
    {
        nvDebugCheck(idx < 16);
        return m_data[idx];
    }
    inline float & Matrix::data(uint idx)
    {
        nvDebugCheck(idx < 16);
        return m_data[idx];
    }
    inline float Matrix::get(uint row, uint col) const
    {
        nvDebugCheck(row < 4 && col < 4);
        return m_data[col * 4 + row];
    }
    inline float Matrix::operator()(uint row, uint col) const
    {
        nvDebugCheck(row < 4 && col < 4);
        return m_data[col * 4 + row];
    }
    inline float & Matrix::operator()(uint row, uint col)
    {
        nvDebugCheck(row < 4 && col < 4);
        return m_data[col * 4 + row];
    }

    inline const float * Matrix::ptr() const
    {
        return m_data;
    }

    inline Vector4 Matrix::row(uint i) const
    {
        nvDebugCheck(i < 4);
        return Vector4(get(i, 0), get(i, 1), get(i, 2), get(i, 3));
    }

    inline Vector4 Matrix::column(uint i) const
    {
        nvDebugCheck(i < 4);
        return Vector4(get(0, i), get(1, i), get(2, i), get(3, i));
    }

    inline void Matrix::zero()
    {
        m_data[0] = 0; m_data[1] = 0; m_data[2] = 0; m_data[3] = 0;
        m_data[4] = 0; m_data[5] = 0; m_data[6] = 0; m_data[7] = 0;
        m_data[8] = 0; m_data[9] = 0; m_data[10] = 0; m_data[11] = 0;
        m_data[12] = 0; m_data[13] = 0; m_data[14] = 0; m_data[15] = 0;
    }

    inline void Matrix::identity()
    {
        m_data[0] = 1; m_data[1] = 0; m_data[2] = 0; m_data[3] = 0;
        m_data[4] = 0; m_data[5] = 1; m_data[6] = 0; m_data[7] = 0;
        m_data[8] = 0; m_data[9] = 0; m_data[10] = 1; m_data[11] = 0;
        m_data[12] = 0; m_data[13] = 0; m_data[14] = 0; m_data[15] = 1;
    }

    // Apply scale.
    inline void Matrix::scale(float s)
    {
        m_data[0] *= s; m_data[1] *= s; m_data[2] *= s; m_data[3] *= s;
        m_data[4] *= s; m_data[5] *= s; m_data[6] *= s; m_data[7] *= s;
        m_data[8] *= s; m_data[9] *= s; m_data[10] *= s; m_data[11] *= s;
        m_data[12] *= s; m_data[13] *= s; m_data[14] *= s; m_data[15] *= s;
    }

    // Apply scale.
    inline void Matrix::scale(Vector3::Arg s)
    {
        m_data[0] *= s.x; m_data[1] *= s.x; m_data[2] *= s.x; m_data[3] *= s.x;
        m_data[4] *= s.y; m_data[5] *= s.y; m_data[6] *= s.y; m_data[7] *= s.y;
        m_data[8] *= s.z; m_data[9] *= s.z; m_data[10] *= s.z; m_data[11] *= s.z;
    }

    // Apply translation.
    inline void Matrix::translate(Vector3::Arg t)
    {
        m_data[12] = m_data[0] * t.x + m_data[4] * t.y + m_data[8]  * t.z + m_data[12];
        m_data[13] = m_data[1] * t.x + m_data[5] * t.y + m_data[9]  * t.z + m_data[13];
        m_data[14] = m_data[2] * t.x + m_data[6] * t.y + m_data[10] * t.z + m_data[14];
        m_data[15] = m_data[3] * t.x + m_data[7] * t.y + m_data[11] * t.z + m_data[15];
    }

    Matrix rotation(float theta, float v0, float v1, float v2);

    // Apply rotation.
    inline void Matrix::rotate(float theta, float v0, float v1, float v2)
    {
        Matrix R(rotation(theta, v0, v1, v2));
        apply(R);
    }

    // Apply transform.
    inline void Matrix::apply(Matrix::Arg m)
    {
        nvDebugCheck(this != &m);

        for(int i = 0; i < 4; i++) {
            const float ai0 = get(i,0), ai1 = get(i,1), ai2 = get(i,2), ai3 = get(i,3);
            m_data[0 + i] = ai0 * m(0,0) + ai1 * m(1,0) + ai2 * m(2,0) + ai3 * m(3,0);
            m_data[4 + i] = ai0 * m(0,1) + ai1 * m(1,1) + ai2 * m(2,1) + ai3 * m(3,1);
            m_data[8 + i] = ai0 * m(0,2) + ai1 * m(1,2) + ai2 * m(2,2) + ai3 * m(3,2);
            m_data[12+ i] = ai0 * m(0,3) + ai1 * m(1,3) + ai2 * m(2,3) + ai3 * m(3,3);
        }
    }

    // Get scale matrix.
    inline Matrix scale(Vector3::Arg s)
    {
        Matrix m(identity);
        m(0,0) = s.x;
        m(1,1) = s.y;
        m(2,2) = s.z;
        return m;
    }

    // Get scale matrix.
    inline Matrix scale(float s)
    {
        Matrix m(identity);
        m(0,0) = m(1,1) = m(2,2) = s;
        return m;
    }

    // Get translation matrix.
    inline Matrix translation(Vector3::Arg t)
    {
        Matrix m(identity);
        m(0,3) = t.x;
        m(1,3) = t.y;
        m(2,3) = t.z;
        return m;
    }

    // Get rotation matrix.
    inline Matrix rotation(float theta, float v0, float v1, float v2)
    {
        float cost = cosf(theta);
        float sint = sinf(theta);

        Matrix m(identity);

        if( 1 == v0 && 0 == v1 && 0 == v2 ) {
            m(1,1) = cost; m(2,1) = -sint;
            m(1,2) = sint; m(2,2) = cost;
        }
        else if( 0 == v0  && 1 == v1 && 0 == v2 ) {
            m(0,0) = cost; m(2,0) = sint;
            m(1,2) = -sint; m(2,2) = cost;
        }
        else if( 0 == v0 && 0 == v1 && 1 == v2 ) {
            m(0,0) = cost; m(1,0) = -sint;
            m(0,1) = sint; m(1,1) = cost;
        } 
        else {
            float a2, b2, c2;
            a2 = v0 * v0;
            b2 = v1 * v1;
            c2 = v2 * v2;

            float iscale = 1.0f / sqrtf(a2 + b2 + c2);
            v0 *= iscale;
            v1 *= iscale;
            v2 *= iscale;

            float abm, acm, bcm;
            float mcos, asin, bsin, csin;
            mcos = 1.0f - cost;
            abm = v0 * v1 * mcos;
            acm = v0 * v2 * mcos;
            bcm = v1 * v2 * mcos;
            asin = v0 * sint;
            bsin = v1 * sint;
            csin = v2 * sint;
            m(0,0) = a2 * mcos + cost;
            m(1,0) = abm - csin;
            m(2,0) = acm + bsin;
            m(3,0) = abm + csin;
            m(1,1) = b2 * mcos + cost;
            m(2,1) = bcm - asin;
            m(3,1) = acm - bsin;
            m(1,2) = bcm + asin;
            m(2,2) = c2 * mcos + cost;
        }
        return m;
    }

    // Get frustum matrix.
    inline Matrix frustum(float xmin, float xmax, float ymin, float ymax, float zNear, float zFar)
    {
        Matrix m(0.0f);

        float doubleznear = 2.0f * zNear;
        float one_deltax = 1.0f / (xmax - xmin);
        float one_deltay = 1.0f / (ymax - ymin);
        float one_deltaz = 1.0f / (zFar - zNear);

        m(0,0) = doubleznear * one_deltax;
        m(1,1) = doubleznear * one_deltay;
        m(0,2) = (xmax + xmin) * one_deltax;
        m(1,2) = (ymax + ymin) * one_deltay;
        m(2,2) = -(zFar + zNear) * one_deltaz;
        m(3,2) = -1.0f;
        m(2,3) = -(zFar * doubleznear) * one_deltaz;

        return m;
    }

    // Get inverse frustum matrix.
    inline Matrix frustumInverse(float xmin, float xmax, float ymin, float ymax, float zNear, float zFar)
    {
        Matrix m(0.0f);

        float one_doubleznear = 1.0f / (2.0f * zNear);
        float one_doubleznearzfar = 1.0f / (2.0f * zNear * zFar);

        m(0,0) = (xmax - xmin) * one_doubleznear;
        m(0,3) = (xmax + xmin) * one_doubleznear;
        m(1,1) = (ymax - ymin) * one_doubleznear;
        m(1,3) = (ymax + ymin) * one_doubleznear;
        m(2,3) = -1;
        m(3,2) = -(zFar - zNear) * one_doubleznearzfar;
        m(3,3) = (zFar + zNear) * one_doubleznearzfar;

        return m;
    }

    // Get infinite frustum matrix.
    inline Matrix frustum(float xmin, float xmax, float ymin, float ymax, float zNear)
    {
        Matrix m(0.0f);

        float doubleznear = 2.0f * zNear;
        float one_deltax = 1.0f / (xmax - xmin);
        float one_deltay = 1.0f / (ymax - ymin);
        float nudge = 1.0; // 0.999;

        m(0,0) = doubleznear * one_deltax;
        m(1,1) = doubleznear * one_deltay;
        m(0,2) = (xmax + xmin) * one_deltax;
        m(1,2) = (ymax + ymin) * one_deltay;
        m(2,2) = -1.0f * nudge;
        m(3,2) = -1.0f;
        m(2,3) = -doubleznear * nudge;

        return m;
    }

    // Get perspective matrix.
    inline Matrix perspective(float fovy, float aspect, float zNear, float zFar)
    {
        float xmax = zNear * tan(fovy / 2);
        float xmin = -xmax;

        float ymax = xmax / aspect;
        float ymin = -ymax;

        return frustum(xmin, xmax, ymin, ymax, zNear, zFar);	
    }

    // Get inverse perspective matrix.
    inline Matrix perspectiveInverse(float fovy, float aspect, float zNear, float zFar)
    {
        float xmax = zNear * tan(fovy / 2);
        float xmin = -xmax;

        float ymax = xmax / aspect;
        float ymin = -ymax;

        return frustumInverse(xmin, xmax, ymin, ymax, zNear, zFar);	
    }

    // Get infinite perspective matrix.
    inline Matrix perspective(float fovy, float aspect, float zNear)
    {
        float x = zNear * tan(fovy / 2);
        float y = x / aspect;
        return frustum( -x, x, -y, y, zNear );	
    }

    // Get matrix determinant.
    inline float Matrix::determinant() const
    {
        return 
            m_data[3] * m_data[6] * m_data[ 9] * m_data[12] - m_data[2] * m_data[7] * m_data[ 9] * m_data[12] - m_data[3] * m_data[5] * m_data[10] * m_data[12] + m_data[1] * m_data[7] * m_data[10] * m_data[12] +
            m_data[2] * m_data[5] * m_data[11] * m_data[12] - m_data[1] * m_data[6] * m_data[11] * m_data[12] - m_data[3] * m_data[6] * m_data[ 8] * m_data[13] + m_data[2] * m_data[7] * m_data[ 8] * m_data[13] +
            m_data[3] * m_data[4] * m_data[10] * m_data[13] - m_data[0] * m_data[7] * m_data[10] * m_data[13] - m_data[2] * m_data[4] * m_data[11] * m_data[13] + m_data[0] * m_data[6] * m_data[11] * m_data[13] +
            m_data[3] * m_data[5] * m_data[ 8] * m_data[14] - m_data[1] * m_data[7] * m_data[ 8] * m_data[14] - m_data[3] * m_data[4] * m_data[ 9] * m_data[14] + m_data[0] * m_data[7] * m_data[ 9] * m_data[14] +
            m_data[1] * m_data[4] * m_data[11] * m_data[14] - m_data[0] * m_data[5] * m_data[11] * m_data[14] - m_data[2] * m_data[5] * m_data[ 8] * m_data[15] + m_data[1] * m_data[6] * m_data[ 8] * m_data[15] +
            m_data[2] * m_data[4] * m_data[ 9] * m_data[15] - m_data[0] * m_data[6] * m_data[ 9] * m_data[15] - m_data[1] * m_data[4] * m_data[10] * m_data[15] + m_data[0] * m_data[5] * m_data[10] * m_data[15];
    }

    inline Matrix transpose(Matrix::Arg m)
    {
        Matrix r;
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                r(i, j) = m(j, i);
            }
        }
        return r;
    }

    inline Matrix isometryInverse(Matrix::Arg m)
    {
        Matrix r(identity);

        // transposed 3x3 upper left matrix
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                r(i, j) = m(j, i);
            }
        }

        // translate by the negative offsets
        r.translate(-Vector3(m.data(12), m.data(13), m.data(14)));

        return r;
    }

    // Transform the given 3d point with the given matrix.
    inline Vector3 transformPoint(Matrix::Arg m, Vector3::Arg p)
    {
        return Vector3(
            p.x * m(0,0) + p.y * m(0,1) + p.z * m(0,2) + m(0,3),
            p.x * m(1,0) + p.y * m(1,1) + p.z * m(1,2) + m(1,3),
            p.x * m(2,0) + p.y * m(2,1) + p.z * m(2,2) + m(2,3));
    }

    // Transform the given 3d vector with the given matrix.
    inline Vector3 transformVector(Matrix::Arg m, Vector3::Arg p)
    {
        return Vector3(
            p.x * m(0,0) + p.y * m(0,1) + p.z * m(0,2),
            p.x * m(1,0) + p.y * m(1,1) + p.z * m(1,2),
            p.x * m(2,0) + p.y * m(2,1) + p.z * m(2,2));
    }

    // Transform the given 4d vector with the given matrix.
    inline Vector4 transform(Matrix::Arg m, Vector4::Arg p)
    {
        return Vector4(
            p.x * m(0,0) + p.y * m(0,1) + p.z * m(0,2) + p.w * m(0,3),
            p.x * m(1,0) + p.y * m(1,1) + p.z * m(1,2) + p.w * m(1,3),
            p.x * m(2,0) + p.y * m(2,1) + p.z * m(2,2) + p.w * m(2,3),
            p.x * m(3,0) + p.y * m(3,1) + p.z * m(3,2) + p.w * m(3,3));
    }

    inline Matrix mul(Matrix::Arg a, Matrix::Arg b)
    {
        // @@ Is this the right order? mul(a, b) = b * a
        Matrix m = a;
        m.apply(b);
        return m;
    }

    inline void Matrix::operator+=(const Matrix & m)
    {
        for(int i = 0; i < 16; i++) {
            m_data[i] += m.m_data[i];
        }
    }

    inline void Matrix::operator-=(const Matrix & m)
    {
        for(int i = 0; i < 16; i++) {
            m_data[i] -= m.m_data[i];
        }
    }

    inline Matrix operator+(const Matrix & a, const Matrix & b)
    {
        Matrix m = a;
        m += b;
        return m;
    }

    inline Matrix operator-(const Matrix & a, const Matrix & b)
    {
        Matrix m = a;
        m -= b;
        return m;
    }


} // nv namespace

#endif // NV_MATH_MATRIX_INL
