// This code is in the public domain -- castanyo@yahoo.es

#pragma once
#ifndef NV_MATH_COLOR_H
#define NV_MATH_COLOR_H

#include "nvmath.h"

namespace nv
{
    /// 32 bit color stored as BGRA.
    class Color32
    {
    public:
        Color32() { }
        Color32(const Color32 & c) : u(c.u) { }
        Color32(uint8 R, uint8 G, uint8 B) { setRGBA(R, G, B, 0xFF); }
        Color32(uint8 R, uint8 G, uint8 B, uint8 A) { setRGBA( R, G, B, A); }
        //Color32(uint8 c[4]) { setRGBA(c[0], c[1], c[2], c[3]); }
        //Color32(float R, float G, float B) { setRGBA(uint(R*255), uint(G*255), uint(B*255), 0xFF); }
        //Color32(float R, float G, float B, float A) { setRGBA(uint(R*255), uint(G*255), uint(B*255), uint(A*255)); }
        explicit Color32(uint32 U) : u(U) { }

        void setRGBA(uint8 R, uint8 G, uint8 B, uint8 A)
        {
            r = R;
            g = G;
            b = B;
            a = A;
        }

        void setBGRA(uint8 B, uint8 G, uint8 R, uint8 A = 0xFF)
        {
            r = R;
            g = G;
            b = B;
            a = A;
        }

        operator uint32 () const {
            return u;
        }

        union {
            struct {
#if NV_LITTLE_ENDIAN
                uint8 b, g, r, a;
#else
                uint8 a: 8;
                uint8 r: 8;
                uint8 g: 8;
                uint8 b: 8;
#endif
            };
            uint8 component[4];
            uint32 u;
        };
    };

} // nv namespace

#endif // NV_MATH_COLOR_H
