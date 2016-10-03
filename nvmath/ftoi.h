// This code is in the public domain -- castano@gmail.com

#pragma once
#ifndef NV_MATH_FTOI_H
#define NV_MATH_FTOI_H

#include "nvmath/nvmath.h"

#include <math.h>

namespace nv
{

NV_FORCEINLINE int ftoi_floor(float val)
{
	return (int)val;
}

NV_FORCEINLINE int ftoi_ceil(float val)
{
	return (int)ceilf(val);
}
} // nv

#endif // NV_MATH_FTOI_H
