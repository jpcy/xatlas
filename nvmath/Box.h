// This code is in the public domain -- castanyo@yahoo.es

#pragma once
#ifndef NV_MATH_BOX_H
#define NV_MATH_BOX_H

#include "Vector.h"
#include "Vector.h"

#include <float.h> // FLT_MAX

namespace nv
{
    class Vector;

    // Axis Aligned Bounding Box.
    class Box
    {
    public:

        inline Box() {}
        inline Box(const Box & b) : minCorner(b.minCorner), maxCorner(b.maxCorner) {}
        inline Box(const Vector3 & mins, const Vector3 & maxs) : minCorner(mins), maxCorner(maxs) {}

        operator const float * () const { return reinterpret_cast<const float *>(this); }

        // Clear the bounds.
        void clearBounds()
		{
			minCorner.set(FLT_MAX, FLT_MAX, FLT_MAX);
			maxCorner.set(-FLT_MAX, -FLT_MAX, -FLT_MAX);
		}

        // Return extents of the box.
		Vector3 Box::extents() const
		{
			return (maxCorner - minCorner) * 0.5f;
		}

        // Add a point to this box.
        void addPointToBounds(const Vector3 & p)
		{
			minCorner = min(minCorner, p);
			maxCorner = max(maxCorner, p);
		}
 
        // Get the volume of the box.
        float volume() const
		{
			Vector3 d = extents();
			return 8.0f * (d.x * d.y * d.z);
		}

        const Vector3 & corner(int i) const { return (&minCorner)[i]; }

        Vector3 minCorner;
        Vector3 maxCorner;
    };
} // nv namespace


#endif // NV_MATH_BOX_H
