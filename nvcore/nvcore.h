// This code is in the public domain -- Ignacio Castaño <castano@gmail.com>

#pragma once
#ifndef NV_CORE_H
#define NV_CORE_H

// Platform definitions
#include <posh.h>

// OS:
// NV_OS_WIN32
// NV_OS_WIN64
// NV_OS_MINGW
// NV_OS_CYGWIN
// NV_OS_LINUX
// NV_OS_UNIX
// NV_OS_DARWIN
// NV_OS_XBOX
// NV_OS_ORBIS
// NV_OS_IOS

#if defined POSH_OS_LINUX
#   define NV_OS_LINUX 1
#   define NV_OS_UNIX 1
#elif defined POSH_OS_ORBIS
#   define NV_OS_ORBIS 1
#elif defined POSH_OS_FREEBSD
#   define NV_OS_FREEBSD 1
#   define NV_OS_UNIX 1
#elif defined POSH_OS_OPENBSD
#   define NV_OS_OPENBSD 1
#   define NV_OS_UNIX 1
#elif defined POSH_OS_CYGWIN32
#   define NV_OS_CYGWIN 1
#elif defined POSH_OS_MINGW
#   define NV_OS_MINGW 1
#   define NV_OS_WIN32 1
#elif defined POSH_OS_OSX
#   define NV_OS_OSX 1      // IC: Adding this, because iOS defines NV_OS_DARWIN too.
#   define NV_OS_DARWIN 1
#   define NV_OS_UNIX 1
#elif defined POSH_OS_IOS
#   define NV_OS_DARWIN 1 //ACS should we keep this on IOS?
#   define NV_OS_UNIX 1
#   define NV_OS_IOS 1
#elif defined POSH_OS_UNIX
#   define NV_OS_UNIX 1
#elif defined POSH_OS_WIN64
#   define NV_OS_WIN32 1
#   define NV_OS_WIN64 1
#elif defined POSH_OS_WIN32
#   define NV_OS_WIN32 1
#elif defined POSH_OS_XBOX
#   define NV_OS_XBOX 1
#elif defined POSH_OS_DURANGO
#   define NV_OS_DURANGO 1
#else
#   error "Unsupported OS"
#endif

// CPUs:
// NV_CPU_X86
// NV_CPU_X86_64
// NV_CPU_PPC
// NV_CPU_ARM

#define NV_CPU_STRING   POSH_CPU_STRING

#if defined POSH_CPU_X86_64
//#   define NV_CPU_X86 1
#   define NV_CPU_X86_64 1
#elif defined POSH_CPU_X86
#   define NV_CPU_X86 1
#elif defined POSH_CPU_PPC
#   define NV_CPU_PPC 1
#elif defined POSH_CPU_STRONGARM
#   define NV_CPU_ARM 1
#else
#   error "Unsupported CPU"
#endif


// Compiler:
// NV_CC_GNUC
// NV_CC_MSVC
// NV_CC_CLANG

#if defined POSH_COMPILER_CLANG
#   define NV_CC_CLANG  1
#   define NV_CC_GNUC   1    // Clang is compatible with GCC.
#elif defined POSH_COMPILER_GCC
#   define NV_CC_GNUC   1
#elif defined POSH_COMPILER_MSVC
#   define NV_CC_MSVC   1
#else
#   error "Unsupported compiler"
#endif

// Endiannes:
#define NV_LITTLE_ENDIAN    POSH_LITTLE_ENDIAN
#define NV_BIG_ENDIAN       POSH_BIG_ENDIAN

#if NV_OS_DARWIN
#include <stdint.h>
//#include <inttypes.h>

// Type definitions:
typedef uint8_t     uint8;
typedef int8_t      int8;

typedef uint16_t    uint16;
typedef int16_t     int16;

typedef uint32_t    uint32;
typedef int32_t     int32;

typedef uint64_t    uint64;
typedef int64_t     int64;

// POSH gets this wrong due to __LP64__
#undef POSH_I64_PRINTF_PREFIX
#define POSH_I64_PRINTF_PREFIX "ll"

#else

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
#endif

// Aliases
typedef uint32      uint;

// Null index. @@ Move this somewhere else... it's only used by nvmesh.
//const unsigned int NIL = unsigned int(~0);
#define NIL uint(~0)

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
