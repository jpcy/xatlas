// This code is in the public domain -- Ignacio Castaño <castano@gmail.com>

#include "Basis.h"

using namespace nv;

/// Build an arbitrary frame for the given direction.
void Basis::buildFrameForDirection(Vector3::Arg d, float angle/*= 0*/)
{
    nvCheck(isNormalized(d));
    normal = d;

    // Choose minimum axis.
    if (fabsf(normal.x) < fabsf(normal.y) && fabsf(normal.x) < fabsf(normal.z))
    {
        tangent = Vector3(1, 0, 0);
    }
    else if (fabsf(normal.y) < fabsf(normal.z))
    {
        tangent = Vector3(0, 1, 0);
    }
    else
    {
        tangent = Vector3(0, 0, 1);
    }

    // Ortogonalize
    tangent -= normal * dot(normal, tangent);
    tangent = ::normalize(tangent);

    bitangent = cross(normal, tangent);

    // Rotate frame around normal according to angle.
    if (angle != 0.0f) {
        float c = cosf(angle);
        float s = sinf(angle);
        Vector3 tmp = c * tangent - s * bitangent;
        bitangent = s * tangent + c * bitangent;
        tangent = tmp;
    }
}
