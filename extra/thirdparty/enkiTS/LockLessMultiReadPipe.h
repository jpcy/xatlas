// Copyright (c) 2013 Doug Binks
// 
// This software is provided 'as-is', without any express or implied
// warranty. In no event will the authors be held liable for any damages
// arising from the use of this software.
// 
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
// 
// 1. The origin of this software must not be misrepresented; you must not
//    claim that you wrote the original software. If you use this software
//    in a product, an acknowledgement in the product documentation would be
//    appreciated but is not required.
// 2. Altered source versions must be plainly marked as such, and must not be
//    misrepresented as being the original software.
// 3. This notice may not be removed or altered from any source distribution.

#pragma once

#include <stdint.h>
#include <assert.h>

#include <atomic>
#include <string.h>


namespace enki
{
    // LockLessMultiReadPipe - Single writer, multiple reader thread safe pipe using (semi) lockless programming
    // Readers can only read from the back of the pipe
    // The single writer can write to the front of the pipe, and read from both ends (a writer can be a reader)
    // for many of the principles used here, see http://msdn.microsoft.com/en-us/library/windows/desktop/ee418650(v=vs.85).aspx
    // Note: using log2 sizes so we do not need to clamp (multi-operation)
    // T is the contained type
    // Note this is not true lockless as the use of flags as a form of lock state.
    template<uint8_t cSizeLog2, typename T> class LockLessMultiReadPipe
    {
    public:
        LockLessMultiReadPipe();
        ~LockLessMultiReadPipe() {}

        // ReaderTryReadBack returns false if we were unable to read
        // This is thread safe for both multiple readers and the writer
        bool ReaderTryReadBack(   T* pOut );

        // WriterTryReadFront returns false if we were unable to read
        // This is thread safe for the single writer, but should not be called by readers
        bool WriterTryReadFront(  T* pOut );

        // WriterTryWriteFront returns false if we were unable to write
        // This is thread safe for the single writer, but should not be called by readers
        bool WriterTryWriteFront( const T& in );

        // IsPipeEmpty() is a utility function, not intended for general use
        // Should only be used very prudently.
        bool IsPipeEmpty() const
        {
            return 0 == m_WriteIndex.load( std::memory_order_relaxed ) - m_ReadCount.load( std::memory_order_relaxed );
        }

        void Clear()
        {
            m_WriteIndex = 0;
            m_ReadIndex = 0;
            m_ReadCount = 0;
            memset( (void*)m_Flags, 0, sizeof( m_Flags ) );
        }

    private:
        const static uint32_t           ms_cSize        = ( 1 << cSizeLog2 );
        const static uint32_t           ms_cIndexMask   = ms_cSize - 1;
        const static uint32_t           FLAG_INVALID    = 0xFFFFFFFF; // 32bit for CAS
        const static uint32_t           FLAG_CAN_WRITE  = 0x00000000; // 32bit for CAS
        const static uint32_t           FLAG_CAN_READ   = 0x11111111; // 32bit for CAS

        T                               m_Buffer[ ms_cSize ];

        // read and write indexes allow fast access to the pipe, but actual access
        // controlled by the access flags. 
        std::atomic<uint32_t>            m_WriteIndex;
        std::atomic<uint32_t>            m_ReadCount;
        std::atomic<uint32_t>            m_Flags[  ms_cSize ];
        std::atomic<uint32_t>            m_ReadIndex;
    };

    template<uint8_t cSizeLog2, typename T> inline
        LockLessMultiReadPipe<cSizeLog2,T>::LockLessMultiReadPipe()
        : m_WriteIndex(0)
        , m_ReadCount(0)
        , m_ReadIndex(0)
    {
        assert( cSizeLog2 < 32 );
        memset( (void*)m_Flags, 0, sizeof( m_Flags ) );
    }

    template<uint8_t cSizeLog2, typename T> inline
        bool LockLessMultiReadPipe<cSizeLog2,T>::ReaderTryReadBack(   T* pOut )
    {

        uint32_t actualReadIndex;
        uint32_t readCount  = m_ReadCount.load( std::memory_order_relaxed );

        // We get hold of read index for consistency
        // and do first pass starting at read count
        uint32_t readIndexToUse  = readCount;
        while(true)
        {

            uint32_t writeIndex = m_WriteIndex.load( std::memory_order_relaxed );
            // power of two sizes ensures we can use a simple calc without modulus
            uint32_t numInPipe = writeIndex - readCount;
            if( 0 == numInPipe )
            {
                return false;
            }
            if( readIndexToUse >= writeIndex )
            {
                readIndexToUse = m_ReadIndex.load( std::memory_order_relaxed );
            }

            // power of two sizes ensures we can perform AND for a modulus
            actualReadIndex    = readIndexToUse & ms_cIndexMask;

            // Multiple potential readers mean we should check if the data is valid,
            // using an atomic compare exchange
            uint32_t previous = FLAG_CAN_READ;
            bool bSuccess = m_Flags[  actualReadIndex ].compare_exchange_strong( previous, FLAG_INVALID, std::memory_order_relaxed );
            if( bSuccess )
            {
                break;
            }
            ++readIndexToUse;

            // Update read count
            readCount  = m_ReadCount.load( std::memory_order_relaxed );
        }

        // we update the read index using an atomic add, as we've only read one piece of data.
        // this ensure consistency of the read index, and the above loop ensures readers
        // only read from unread data
        m_ReadCount.fetch_add(1, std::memory_order_relaxed );

        std::atomic_thread_fence( std::memory_order_acquire );
        // now read data, ensuring we do so after above reads & CAS
        *pOut = m_Buffer[ actualReadIndex ];

        m_Flags[  actualReadIndex ].store( FLAG_CAN_WRITE, std::memory_order_relaxed );

        return true;
    }

    template<uint8_t cSizeLog2, typename T> inline
        bool LockLessMultiReadPipe<cSizeLog2,T>::WriterTryReadFront(  T* pOut )
    {
        uint32_t writeIndex = m_WriteIndex.load( std::memory_order_relaxed );
        uint32_t frontReadIndex  = writeIndex;

        // Multiple potential readers mean we should check if the data is valid,
        // using an atomic compare exchange - which acts as a form of lock (so not quite lockless really).
        uint32_t actualReadIndex    = 0;
        while(true)
        {
            uint32_t readCount  = m_ReadCount.load( std::memory_order_relaxed );
            // power of two sizes ensures we can use a simple calc without modulus
            uint32_t numInPipe = writeIndex - readCount;
            if( 0 == numInPipe )
            {
                m_ReadIndex.store( readCount, std::memory_order_release );
                return false;
            }
            --frontReadIndex;
            actualReadIndex    = frontReadIndex & ms_cIndexMask;
            uint32_t previous = FLAG_CAN_READ;
            bool success = m_Flags[  actualReadIndex ].compare_exchange_strong( previous, FLAG_INVALID, std::memory_order_relaxed );
            if( success )
            {
                break;
            }
            else if( m_ReadIndex.load( std::memory_order_acquire ) >= frontReadIndex  )
            {
                return false;
            }
        }

        std::atomic_thread_fence( std::memory_order_acquire );
        // now read data, ensuring we do so after above reads & CAS
        *pOut = m_Buffer[ actualReadIndex ];

        m_Flags[  actualReadIndex ].store( FLAG_CAN_WRITE, std::memory_order_relaxed );

        m_WriteIndex.store(writeIndex-1, std::memory_order_relaxed);
        return true;
    }


    template<uint8_t cSizeLog2, typename T> inline
        bool LockLessMultiReadPipe<cSizeLog2,T>::WriterTryWriteFront( const T& in )
    {
        // The writer 'owns' the write index, and readers can only reduce
        // the amount of data in the pipe.
        // We get hold of both values for consistency and to reduce false sharing
        // impacting more than one access
        uint32_t writeIndex = m_WriteIndex;

        // power of two sizes ensures we can perform AND for a modulus
        uint32_t actualWriteIndex    = writeIndex & ms_cIndexMask;

        // a reader may still be reading this item, as there are multiple readers
        if( m_Flags[ actualWriteIndex ].load(std::memory_order_relaxed)  != FLAG_CAN_WRITE ) 
        {
            return false; // still being read, so have caught up with tail. 
        }

        // as we are the only writer we can update the data without atomics
        //  whilst the write index has not been updated
        m_Buffer[ actualWriteIndex ] = in;
        m_Flags[  actualWriteIndex ].store( FLAG_CAN_READ, std::memory_order_relaxed );

        // We need to ensure the above writes occur prior to updating the write index,
        // otherwise another thread might read before it's finished
        std::atomic_thread_fence( std::memory_order_release );

        m_WriteIndex.fetch_add(1, std::memory_order_relaxed);
        return true;
    }


    // Lockless multiwriter intrusive list
    // Type T must implement T* volatile pNext;
    template<typename T> class  LocklessMultiWriteIntrusiveList
    {

        std::atomic<T*> pHead;
        T               tail;
    public:
        LocklessMultiWriteIntrusiveList() : pHead( &tail )
        {
            tail.pNext = NULL;
        }

        bool IsListEmpty() const
        {
            return pHead == &tail;
        }

        // Add - safe to perform from any thread
        void WriterWriteFront( T* pNode_ )
        {
            assert( pNode_ );
            pNode_->pNext = NULL;
            T* pPrev = pHead.exchange( pNode_ );
            pPrev->pNext = pNode_;
        }

        // Remove - only thread safe for owner
        T* ReaderReadBack()
        {
            T* pTailPlus1 = tail.pNext;
            if( pTailPlus1 )
            {
                T* pTailPlus2 = pTailPlus1->pNext;
                if( pTailPlus2 )
                {
                    //not head
                    tail.pNext = pTailPlus2;
                }
                else
                {
                    tail.pNext = NULL;
                    T* pCompare = pTailPlus1; // we need preserve pTailPlus1 as compare will alter it on failure
                    // pTailPlus1 is the head, attempt swap with tail
                    if( !pHead.compare_exchange_strong( pCompare, &tail ) )
                    {
                        // pCompare receives the revised pHead on failure.
                        // pTailPlus1 is no longer the head, so pTailPlus1->pNext should be non NULL
                        while( (T*)NULL == pTailPlus1->pNext ) {;} // wait for pNext to be updated as head may have just changed.
                        tail.pNext = pTailPlus1->pNext.load();
                        pTailPlus1->pNext = NULL;
                    }
                }
            }
            return pTailPlus1;
        }
    };

}
