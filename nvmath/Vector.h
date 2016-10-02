// This code is in the public domain -- castanyo@yahoo.es

#pragma once
#ifndef NV_MATH_VECTOR_H
#define NV_MATH_VECTOR_H

#include "nvmath.h"

namespace nv
{
    class Vector2
    {
    public:
        typedef Vector2 const & Arg;

        Vector2();
        explicit Vector2(float f);
        Vector2(float x, float y);
        Vector2(Vector2::Arg v);

        //template <typename T> explicit Vector2(const T & v) : x(v.x), y(v.y) {}
        //template <typename T> operator T() const { return T(x, y); }

        const Vector2 & operator=(Vector2::Arg v);

        const float * ptr() const;

        void set(float x, float y);

        Vector2 operator-() const;
        void operator+=(Vector2::Arg v);
        void operator-=(Vector2::Arg v);
        void operator*=(float s);
        void operator*=(Vector2::Arg v);

        friend bool operator==(Vector2::Arg a, Vector2::Arg b);
        friend bool operator!=(Vector2::Arg a, Vector2::Arg b);

        union {
            struct {
                float x, y;
            };
            float component[2];
        };
    };

    class Vector3
    {
    public:
        typedef Vector3 const & Arg;

        Vector3();
        explicit Vector3(float x);
        //explicit Vector3(int x) : x(float(x)), y(float(x)), z(float(x)) {}
        Vector3(float x, float y, float z);
        Vector3(Vector2::Arg v, float z);
        Vector3(Vector3::Arg v);

        //template <typename T> explicit Vector3(const T & v) : x(v.x), y(v.y), z(v.z) {}
        //template <typename T> operator T() const { return T(x, y, z); }

        const Vector3 & operator=(Vector3::Arg v);

        Vector2 xy() const;

        const float * ptr() const;

        void set(float x, float y, float z);

        Vector3 operator-() const;
        void operator+=(Vector3::Arg v);
        void operator-=(Vector3::Arg v);
        void operator*=(float s);
        void operator/=(float s);
        void operator*=(Vector3::Arg v);
        void operator/=(Vector3::Arg v);

        friend bool operator==(Vector3::Arg a, Vector3::Arg b);
        friend bool operator!=(Vector3::Arg a, Vector3::Arg b);

        union {
            struct {
                float x, y, z;
            };
            float component[3];
        };
    };

} // nv namespace

#endif // NV_MATH_VECTOR_H
