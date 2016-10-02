// This code is in the public domain -- castanyo@yahoo.es

#pragma once
#ifndef NV_MESH_RASTER_H
#define NV_MESH_RASTER_H

/** @file Raster.h
 * @brief Rasterization library.
 *
 * This is just a standard scanline rasterizer that I took from one of my old
 * projects. The perspective correction wasn't necessary so I just removed it.
**/

#include "nvmath/Vector.h"

namespace nv
{
    namespace Raster 
    {
        enum Mode {
            Mode_Nearest,
            Mode_Antialiased,
            //Mode_Conservative
        };

        /// A callback to sample the environment. Return false to terminate rasterization.
        typedef bool (* SamplingCallback)(void * param, int x, int y, Vector3::Arg bar, Vector3::Arg dx, Vector3::Arg dy, float coverage);

        // Process the given triangle. Returns false if rasterization was interrupted by the callback.
        bool drawTriangle(Mode mode, Vector2::Arg extents, bool enableScissors, const Vector2 v[3], SamplingCallback cb, void * param);

        // Process the given quad. Returns false if rasterization was interrupted by the callback.
        bool drawQuad(Mode mode, Vector2::Arg extents, bool enableScissors, const Vector2 v[4], SamplingCallback cb, void * param);
    }
}


#endif // NV_MESH_RASTER_H
