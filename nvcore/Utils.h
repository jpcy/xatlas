// This code is in the public domain -- Ignacio Castaño <castano@gmail.com>

#pragma once
#ifndef NV_CORE_UTILS_H
#define NV_CORE_UTILS_H

#include "Debug.h" // nvDebugCheck

#include <new> // for placement new


// Just in case. Grrr.
#undef min
#undef max

#define NV_UINT32_MAX   0xffffffff
#define NV_FLOAT_MAX    3.402823466e+38F

namespace nv
{
    /// Swap two values.
    template <typename T> 
    inline void swap(T & a, T & b)
    {
        T temp(a);
        a = b; 
        b = temp;
    }

    /// Return the maximum of the two arguments. For floating point values, it returns the second value if the first is NaN.
    template <typename T> 
    //inline const T & max(const T & a, const T & b)
    inline T max(const T & a, const T & b)
    {
        return (b < a) ? a : b;
    }

    /// Return the maximum of the three arguments.
    template <typename T> 
    //inline const T & max3(const T & a, const T & b, const T & c)
    inline T max3(const T & a, const T & b, const T & c)
    {
        return max(a, max(b, c));
    }

    /// Return the minimum of two values.
    template <typename T> 
    //inline const T & min(const T & a, const T & b)
    inline T min(const T & a, const T & b)
    {
        return (a < b) ? a : b;
    }

    /// Return the maximum of the three arguments.
    template <typename T> 
    //inline const T & min3(const T & a, const T & b, const T & c)
    inline T min3(const T & a, const T & b, const T & c)
    {
        return min(a, min(b, c));
    }

    /// Clamp between two values.
    template <typename T> 
    //inline const T & clamp(const T & x, const T & a, const T & b)
    inline T clamp(const T & x, const T & a, const T & b)
    {
        return min(max(x, a), b);
    }

    /** Return the next power of two. 
    * @see http://graphics.stanford.edu/~seander/bithacks.html
    * @warning Behaviour for 0 is undefined.
    * @note isPowerOfTwo(x) == true -> nextPowerOfTwo(x) == x
    * @note nextPowerOfTwo(x) = 2 << log2(x-1)
    */
    inline uint32 nextPowerOfTwo(uint32 x)
    {
        nvDebugCheck( x != 0 );
#if 1	// On modern CPUs this is supposed to be as fast as using the bsr instruction.
        x--;
        x |= x >> 1;
        x |= x >> 2;
        x |= x >> 4;
        x |= x >> 8;
        x |= x >> 16;
        return x+1;	
#else
        uint p = 1;
        while( x > p ) {
            p += p;
        }
        return p;
#endif
    }

    inline uint64 nextPowerOfTwo(uint64 x)
    {
        nvDebugCheck(x != 0);
        uint p = 1;
        while (x > p) {
            p += p;
        }
        return p;
    }

    // @@ Move this to utils?
    /// Delete all the elements of a container.
    template <typename T>
    void deleteAll(T & container)
    {
        for (typename T::PseudoIndex i = container.start(); !container.isDone(i); container.advance(i))
        {
            delete container[i];
        }
    }

    // @@ Specialize these methods for numeric, pointer, and pod types.

    template <typename T>
    void construct_range(T * restrict ptr, uint new_size, uint old_size) {
        for (uint i = old_size; i < new_size; i++) {
            new(ptr+i) T; // placement new
        }
    }

    template <typename T>
    void construct_range(T * restrict ptr, uint new_size, uint old_size, const T & elem) {
        for (uint i = old_size; i < new_size; i++) {
            new(ptr+i) T(elem); // placement new
        }
    }

    template <typename T>
    void construct_range(T * restrict ptr, uint new_size, uint old_size, const T * src) {
        for (uint i = old_size; i < new_size; i++) {
            new(ptr+i) T(src[i]); // placement new
        }
    }

    template <typename T>
    void destroy_range(T * restrict ptr, uint new_size, uint old_size) {
        for (uint i = new_size; i < old_size; i++) {
            (ptr+i)->~T(); // Explicit call to the destructor
        }
    }

    template <typename T>
    bool find(const T & element, const T * restrict ptr, uint begin, uint end, uint * index) {
        for (uint i = begin; i < end; i++) {
            if (ptr[i] == element) {
                if (index != NULL) *index = i;
                return true;
            }
        }
        return false;
    }

} // nv namespace

#endif // NV_CORE_UTILS_H
