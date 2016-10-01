// This code is in the public domain -- castanyo@yahoo.es

#pragma once
#ifndef NV_MATH_H
#define NV_MATH_H

#include "nvcore/nvcore.h"
#include "nvcore/Debug.h"   // nvDebugCheck
#include "nvcore/Utils.h"   // max, clamp

#include <math.h>

#if NV_OS_WIN32 || NV_OS_XBOX || NV_OS_DURANGO
#include <float.h>  // finite, isnan
#endif

#if NV_CPU_X86 || NV_CPU_X86_64
    //#include <intrin.h>
    #include <xmmintrin.h>
#endif

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
        return fabs(f0-f1) <= epsilon * max3(1.0f, fabsf(f0), fabsf(f1));
    }

    inline bool isZero(const float f, const float epsilon = NV_EPSILON)
    {
        return fabs(f) <= epsilon;
    }

    inline bool isFinite(const float f)
    {
#if NV_OS_WIN32 || NV_OS_XBOX || NV_OS_DURANGO
        return _finite(f) != 0;
#elif NV_OS_DARWIN || NV_OS_FREEBSD || NV_OS_OPENBSD || NV_OS_ORBIS
        return isfinite(f);
#elif NV_OS_LINUX
        return finitef(f);
#else
#   error "isFinite not supported"
#endif
        //return std::isfinite (f);
        //return finite (f);
    }

    inline float lerp(float f0, float f1, float t)
    {
        const float s = 1.0f - t;
        return f0 * s + f1 * t;
    }

    inline float square(float f) { return f * f; }
    inline int square(int i) { return i * i; }

    // Eliminates negative zeros from a float array.
    inline void floatCleanup(float * fp, int n)
    {
        for (int i = 0; i < n; i++) {
            //nvDebugCheck(isFinite(fp[i]));
            union { float f; uint32 i; } x = { fp[i] };
            if (x.i == 0x80000000) fp[i] = 0.0f;
        }
    }
} // nv

#endif // NV_MATH_H
