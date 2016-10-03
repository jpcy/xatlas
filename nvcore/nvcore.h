// This code is in the public domain -- Ignacio Castaño <castano@gmail.com>

#pragma once
#ifndef NV_CORE_H
#define NV_CORE_H

// Platform definitions
#include <posh.h>

// Type definitions:
typedef posh_u8_t   uint8;
typedef posh_i8_t   int8;

typedef posh_u16_t  uint16;
typedef posh_i16_t  int16;

typedef posh_u32_t  uint32;
typedef posh_i32_t  int32;

//#if NV_OS_DARWIN
// OSX-64 is supposed to be LP64 (longs and pointers are 64 bits), thus uint64 is defined as
// unsigned long. However, some OSX headers define it as unsigned long long, producing errors,
// even though both types are 64 bit. Ideally posh should handle that, but it has not been
// updated in ages, so here I'm just falling back to the standard C99 types defined in inttypes.h
//#include <inttypes.h>
//typedef posh_u64_t  uint64_t;
//typedef posh_i64_t  int64_t;
//#else
typedef posh_u64_t  uint64;
typedef posh_i64_t  int64;
//#endif

// Aliases
typedef uint32      uint;

#ifdef _MSC_VER
// Ignore gcc attributes.
#define __attribute__(X)
#define restrict
#define NV_FORCEINLINE __forceinline
#else
#define restrict __restrict__
#define NV_FORCEINLINE  inline __attribute__((always_inline))
#endif

#endif // NV_CORE_H
