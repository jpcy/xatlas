// This code is in the public domain -- Ignacio Castaño <castano@gmail.com>

#pragma once
#ifndef NV_CORE_BITARRAY_H
#define NV_CORE_BITARRAY_H

#include "nvcore.h"
#include "Array.h"

namespace nv
{
    // Simple bit array.
    class BitArray
    {
    public:

        BitArray() {}
        BitArray(uint sz) {
            resize(sz);
        }

        uint size() const { return m_size; }
        void clear() { resize(0); }

        void resize(uint new_size)
        {
            m_size = new_size;
            m_wordArray.resize( (m_size + 31) >> 5 );
        }

        /// Get bit.
        bool bitAt(uint b) const
        {
            nvDebugCheck( b < m_size );
            return (m_wordArray[b >> 5] & (1 << (b & 31))) != 0;
        }

        // Set a bit.
        void setBitAt(uint idx)
        {
            nvDebugCheck(idx < m_size);
            m_wordArray[idx >> 5] |=  (1 << (idx & 31));
        }

        // Toggle a bit.
        void toggleBitAt(uint idx)
        {
            nvDebugCheck(idx < m_size);
            m_wordArray[idx >> 5] ^= (1 << (idx & 31));
        }

        // Set a bit to the given value. @@ Rename modifyBitAt? 
        void setBitAt(uint idx, bool b)
        {
            nvDebugCheck(idx < m_size);
            m_wordArray[idx >> 5] = setBits(m_wordArray[idx >> 5], 1 << (idx & 31), b);
            nvDebugCheck(bitAt(idx) == b);
        }

        // Clear all the bits.
        void clearAll()
        {
            memset(m_wordArray.buffer(), 0, m_wordArray.size() * sizeof(uint));
        }

        // Set all the bits.
        void setAll()
        {
            memset(m_wordArray.buffer(), 0xFF, m_wordArray.size() * sizeof(uint));
        }

    private:
		// See "Conditionally set or clear bits without branching" at http://graphics.stanford.edu/~seander/bithacks.html
		inline uint setBits(uint w, uint m, bool b) {
			return (w & ~m) | (-int(b) & m);
		}

        // Number of bits stored.
        uint m_size;

        // Array of bits.
        Array<uint> m_wordArray;

    };

} // nv namespace

#endif // NV_CORE_BITARRAY_H

