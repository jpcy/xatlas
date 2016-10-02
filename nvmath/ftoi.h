// This code is in the public domain -- castano@gmail.com

#pragma once
#ifndef NV_MATH_FTOI_H
#define NV_MATH_FTOI_H

#if NV_CPU_X86 || NV_CPU_X86_64
    //#include <intrin.h>
    #include <xmmintrin.h>
#endif

#include "nvmath/nvmath.h"

#include <math.h>

namespace nv
{
    // Optimized float to int conversions. See:
    // http://cbloomrants.blogspot.com/2009/01/01-17-09-float-to-int.html
    // http://www.stereopsis.com/sree/fpu2006.html
    // http://assemblyrequired.crashworks.org/2009/01/12/why-you-should-never-cast-floats-to-ints/
    // http://chrishecker.com/Miscellaneous_Technical_Articles#Floating_Point


    union DoubleAnd64 {
        uint64    i;
        double    d;
    };

    static const double floatutil_xs_doublemagic = (6755399441055744.0);                            // 2^52 * 1.5
    static const double floatutil_xs_doublemagicdelta = (1.5e-8);                                   // almost .5f = .5f + 1e^(number of exp bit)
    static const double floatutil_xs_doublemagicroundeps = (0.5f - floatutil_xs_doublemagicdelta);  // almost .5f = .5f - 1e^(number of exp bit)

    NV_FORCEINLINE int ftoi_round_xs(double val, double magic) {
#if 1
        DoubleAnd64 dunion;
        dunion.d = val + magic;
        return (int32) dunion.i; // just cast to grab the bottom bits
#else
        val += magic;
        return ((int*)&val)[0]; // @@ Assumes little endian.
#endif
    }

    NV_FORCEINLINE int ftoi_round_xs(float val) {
        return ftoi_round_xs(val, floatutil_xs_doublemagic);
    }

    NV_FORCEINLINE int ftoi_floor_xs(float val) {
        return ftoi_round_xs(val - floatutil_xs_doublemagicroundeps, floatutil_xs_doublemagic);
    }

    NV_FORCEINLINE int ftoi_ceil_xs(float val) {
        return ftoi_round_xs(val + floatutil_xs_doublemagicroundeps, floatutil_xs_doublemagic);
    }

    NV_FORCEINLINE int ftoi_trunc_xs(float val) {
        return (val<0) ? ftoi_ceil_xs(val) : ftoi_floor_xs(val);
    }

#if NV_CPU_X86 || NV_CPU_X86_64

    NV_FORCEINLINE int ftoi_round_sse(float f) {
        return _mm_cvt_ss2si(_mm_set_ss(f));
    }

    NV_FORCEINLINE int ftoi_trunc_sse(float f) {
      return _mm_cvtt_ss2si(_mm_set_ss(f));
    }

#endif



#if NV_USE_SSE

    NV_FORCEINLINE int ftoi_round(float val) {
        return ftoi_round_sse(val);
    }

    NV_FORCEINLINE int ftoi_trunc(float f) {
      return ftoi_trunc_sse(f);
    }

    // We can probably do better than this. See for example:
    // http://dss.stephanierct.com/DevBlog/?p=8
    NV_FORCEINLINE int ftoi_floor(float val) {
        return ftoi_round(floorf(val));
    }

    NV_FORCEINLINE int ftoi_ceil(float val) {
        return ftoi_round(ceilf(val));
    }

#else

    // In theory this should work with any double floating point math implementation, but it appears that MSVC produces incorrect code
    // when SSE2 is targeted and fast math is enabled (/arch:SSE2 & /fp:fast). These problems go away with /fp:precise, which is the default mode.

    NV_FORCEINLINE int ftoi_round(float val) {
        return ftoi_round_xs(val);
    }

    NV_FORCEINLINE int ftoi_floor(float val) {
        return ftoi_floor_xs(val);
    }

    NV_FORCEINLINE int ftoi_ceil(float val) {
        return ftoi_ceil_xs(val);
    }

    NV_FORCEINLINE int ftoi_trunc(float f) {
      return ftoi_trunc_xs(f);
    }

#endif
} // nv

#endif // NV_MATH_FTOI_H
