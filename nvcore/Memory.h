// This code is in the public domain -- Ignacio Castaño <castano@gmail.com>

#pragma once
#ifndef NV_CORE_MEMORY_H
#define NV_CORE_MEMORY_H

#include "nvcore.h"

#include <stdlib.h> // malloc(), realloc() and free()
#include <string.h> // memset
//#include <stddef.h> // size_t

//#include <new>	// new and delete

#define TRACK_MEMORY_LEAKS 0
#if TRACK_MEMORY_LEAKS
#include <vld.h>
#endif


#if NV_CC_GNUC
#   define NV_ALIGN_16 __attribute__ ((__aligned__ (16)))
#else
#   define NV_ALIGN_16 __declspec(align(16))
#endif


#define NV_OVERRIDE_ALLOC 0

#if NV_OVERRIDE_ALLOC

// Custom memory allocator
extern "C" {
    void * malloc(size_t size);
    void * debug_malloc(size_t size, const char * file, int line);
    void free(void * ptr);
    void * realloc(void * ptr, size_t size);
}

/*
#ifdef _DEBUG
#define new new(__FILE__, __LINE__)
#define malloc(i) debug_malloc(i, __FILE__, __LINE__)
#endif
*/

#endif

namespace nv {
    void * aligned_malloc(size_t size, size_t alignment);
    void aligned_free(void * );

    // C++ helpers.
    template <typename T> NV_FORCEINLINE T * malloc(size_t count) {
        return (T *)::malloc(sizeof(T) * count);
    }

    template <typename T> NV_FORCEINLINE T * realloc(T * ptr, size_t count) {
        return (T *)::realloc(ptr, sizeof(T) * count);
    }

    template <typename T> NV_FORCEINLINE void free(const T * ptr) {
        ::free((void *)ptr);
    }

    template <typename T> NV_FORCEINLINE void zero(T & data) {
        memset(&data, 0, sizeof(T));
    }

} // nv namespace

#endif // NV_CORE_MEMORY_H
