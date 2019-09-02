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

#include <assert.h>

#include "TaskScheduler.h"
#include "LockLessMultiReadPipe.h"

#include <algorithm>

#if defined __i386__ || defined __x86_64__
#include "x86intrin.h"
#elif defined _WIN32
#include <intrin.h>
#endif

using namespace enki;


static const uint32_t PIPESIZE_LOG2              = 8;
static const uint32_t SPIN_COUNT                 = 10;
static const uint32_t SPIN_BACKOFF_MULTIPLIER    = 100;
static const uint32_t MAX_NUM_INITIAL_PARTITIONS = 8;

// thread_local not well supported yet by C++11 compilers.
#ifdef _MSC_VER
    #if _MSC_VER <= 1800
        #define thread_local __declspec(thread)
    #endif
#elif __APPLE__
        // Apple thread_local currently not implemented despite it being in Clang.
        #define thread_local __thread
#endif


// each software thread gets it's own copy of gtl_threadNum, so this is safe to use as a static variable
static thread_local uint32_t                             gtl_threadNum       = 0;

namespace enki 
{
    struct SubTaskSet
    {
        ITaskSet*           pTask;
        TaskSetPartition    partition;
    };

    // we derive class TaskPipe rather than typedef to get forward declaration working easily
    class TaskPipe : public LockLessMultiReadPipe<PIPESIZE_LOG2,enki::SubTaskSet> {};

    struct ThreadArgs
    {
        uint32_t        threadNum;
        TaskScheduler*  pTaskScheduler;
    };

    class PinnedTaskList : public LocklessMultiWriteIntrusiveList<IPinnedTask> {};

    semaphoreid_t* SemaphoreCreate();
    void SemaphoreDelete( semaphoreid_t* pSemaphore_ );
    void SemaphoreWait(   semaphoreid_t& semaphoreid );
    void SemaphoreSignal( semaphoreid_t& semaphoreid, int32_t countWaiting );
}

namespace
{
    SubTaskSet       SplitTask( SubTaskSet& subTask_, uint32_t rangeToSplit_ )
    {
        SubTaskSet splitTask = subTask_;
        uint32_t rangeLeft = subTask_.partition.end - subTask_.partition.start;

        if( rangeToSplit_ > rangeLeft )
        {
            rangeToSplit_ = rangeLeft;
        }
        splitTask.partition.end = subTask_.partition.start + rangeToSplit_;
        subTask_.partition.start = splitTask.partition.end;
        return splitTask;
    }

    #if ( defined _WIN32 && ( defined _M_IX86  || defined _M_X64 ) ) || ( defined __i386__ || defined __x86_64__ )
    static void SpinWait( uint32_t spinCount_ )
    {
        uint64_t end = __rdtsc() + spinCount_;
        while( __rdtsc() < end )
        {
            _mm_pause();
        }        
    }
    #else
    static void SpinWait( uint32_t spinCount_ )
    {
        while( spinCount_ )
        {
            // TODO: may have NOP or yield equiv
            --spinCount_;
        }        
    }
    #endif
}

static void SafeCallback(ProfilerCallbackFunc func_, uint32_t threadnum_)
{
    if( func_ != nullptr )
    {
        func_(threadnum_);
    }
}

ProfilerCallbacks* TaskScheduler::GetProfilerCallbacks()
{
    return &m_ProfilerCallbacks;
}


void TaskScheduler::TaskingThreadFunction( const ThreadArgs& args_ )
{
    uint32_t threadNum  = args_.threadNum;
    TaskScheduler*  pTS = args_.pTaskScheduler;
    gtl_threadNum       = threadNum;

    SafeCallback( pTS->m_ProfilerCallbacks.threadStart, threadNum );

    uint32_t spinCount = 0;
    uint32_t hintPipeToCheck_io = threadNum + 1;    // does not need to be clamped.
    while( pTS->m_bRunning.load( std::memory_order_relaxed ) )
    {
        if(!pTS->TryRunTask( threadNum, hintPipeToCheck_io ) )
        {
            // no tasks, will spin then wait
            ++spinCount;
            if( spinCount > SPIN_COUNT )
            {
                pTS->WaitForTasks( threadNum );
                spinCount = 0;
            }
            else
            {
                // Note: see https://software.intel.com/en-us/articles/a-common-construct-to-avoid-the-contention-of-threads-architecture-agnostic-spin-wait-loops
                uint32_t spinBackoffCount = spinCount * SPIN_BACKOFF_MULTIPLIER;
                SpinWait( spinBackoffCount );
            }
        }
    }

    pTS->m_NumThreadsRunning.fetch_sub( 1, std::memory_order_release );
    SafeCallback( pTS->m_ProfilerCallbacks.threadStop, threadNum );
    return;

}


void TaskScheduler::StartThreads()
{
    if( m_bHaveThreads )
    {
        return;
    }

    for( int priority = 0; priority < TASK_PRIORITY_NUM; ++priority )
    {
        m_pPipesPerThread[ priority ]          = new TaskPipe[ m_NumThreads ];
        m_pPinnedTaskListPerThread[ priority ] = new PinnedTaskList[ m_NumThreads ];
    }

    m_pNewTaskSemaphore = SemaphoreCreate();

    // we create one less thread than m_NumThreads as the main thread counts as one
    m_pThreadArgStore   = new ThreadArgs[m_NumThreads];
    m_pThreads          = new std::thread*[m_NumThreads];
    m_pThreadArgStore[0].threadNum      = 0;
    m_pThreadArgStore[0].pTaskScheduler = this;
    m_NumThreadsRunning = 1; // account for main thread
    m_bRunning = 1;

   for( uint32_t thread = 1; thread < m_NumThreads; ++thread )
    {
        m_pThreadArgStore[thread].threadNum      = thread;
        m_pThreadArgStore[thread].pTaskScheduler = this;
        m_pThreads[thread] = new std::thread( TaskingThreadFunction, m_pThreadArgStore[thread] );
        ++m_NumThreadsRunning;
    }

    // ensure we have sufficient tasks to equally fill either all threads including main
    // or just the threads we've launched, this is outside the firstinit as we want to be able
    // to runtime change it
    if( 1 == m_NumThreads )
    {
        m_NumPartitions = 1;
        m_NumInitialPartitions = 1;
    }
    else
    {
        m_NumPartitions = m_NumThreads * (m_NumThreads - 1);
        m_NumInitialPartitions = m_NumThreads - 1;
        if( m_NumInitialPartitions > MAX_NUM_INITIAL_PARTITIONS )
        {
            m_NumInitialPartitions = MAX_NUM_INITIAL_PARTITIONS;
        }
    }

    m_bHaveThreads = true;
}

void TaskScheduler::StopThreads( bool bWait_ )
{
    if( m_bHaveThreads )
    {
        // wait for them threads quit before deleting data
        m_bRunning = 0;
        while( bWait_ && m_NumThreadsRunning > 1)
        {
            // keep firing event to ensure all threads pick up state of m_bRunning
           WakeThreads();
        }

        for( uint32_t thread = 1; thread < m_NumThreads; ++thread )
        {
            m_pThreads[thread]->detach();
            delete m_pThreads[thread];
        }

        m_NumThreads = 0;
        delete[] m_pThreadArgStore;
        delete[] m_pThreads;
        m_pThreadArgStore = 0;
        m_pThreads = 0;

        SemaphoreDelete( m_pNewTaskSemaphore );
        m_pNewTaskSemaphore = 0;

        m_bHaveThreads = false;
        m_NumThreadsWaiting = 0;
        m_NumThreadsRunning = 0;

        for( int priority = 0; priority < TASK_PRIORITY_NUM; ++priority )
        {
            delete[] m_pPipesPerThread[ priority ];
            m_pPipesPerThread[ priority ] = NULL;
            delete[] m_pPinnedTaskListPerThread[ priority ];
            m_pPinnedTaskListPerThread[ priority ] = NULL;
        }
    }
}

bool TaskScheduler::TryRunTask( uint32_t threadNum_, uint32_t& hintPipeToCheck_io_ )
{
    for( int priority = 0; priority < TASK_PRIORITY_NUM; ++priority )
    {
        if( TryRunTask( threadNum_, priority, hintPipeToCheck_io_ ) )
        {
            return true;
        }
    }
    return false;
}

bool TaskScheduler::TryRunTask( uint32_t threadNum, uint32_t priority_, uint32_t& hintPipeToCheck_io_ )
{
    // Run any tasks for this thread
    RunPinnedTasks( threadNum, priority_ );

    // check for tasks
    SubTaskSet subTask;
    bool bHaveTask = m_pPipesPerThread[ priority_ ][ threadNum ].WriterTryReadFront( &subTask );

    uint32_t threadToCheck = hintPipeToCheck_io_;
    uint32_t checkCount = 0;
    while( !bHaveTask && checkCount < m_NumThreads )
    {
        threadToCheck = ( hintPipeToCheck_io_ + checkCount ) % m_NumThreads;
        if( threadToCheck != threadNum )
        {
            bHaveTask = m_pPipesPerThread[ priority_ ][ threadToCheck ].ReaderTryReadBack( &subTask );
        }
        ++checkCount;
    }
        
    if( bHaveTask )
    {
        // update hint, will preserve value unless actually got task from another thread.
        hintPipeToCheck_io_ = threadToCheck;

        uint32_t partitionSize = subTask.partition.end - subTask.partition.start;
        if( subTask.pTask->m_RangeToRun < partitionSize )
        {
            SubTaskSet taskToRun = SplitTask( subTask, subTask.pTask->m_RangeToRun );
            SplitAndAddTask( threadNum, subTask, subTask.pTask->m_RangeToRun );
            taskToRun.pTask->ExecuteRange( taskToRun.partition, threadNum );
            taskToRun.pTask->m_RunningCount.fetch_sub(1,std::memory_order_release );

        }
        else
        {
            // the task has already been divided up by AddTaskSetToPipe, so just run it
            subTask.pTask->ExecuteRange( subTask.partition, threadNum );
            subTask.pTask->m_RunningCount.fetch_sub(1,std::memory_order_release );
        }
    }

    return bHaveTask;

}

bool TaskScheduler::HaveTasks(  uint32_t threadNum_ )
{
    for( int priority = 0; priority < TASK_PRIORITY_NUM; ++priority )
    {
        for( uint32_t thread = 0; thread < m_NumThreads; ++thread )
        {
            if( !m_pPipesPerThread[ priority ][ thread ].IsPipeEmpty() )
            {
                return true;
            }
        }
        if( !m_pPinnedTaskListPerThread[ priority ][ threadNum_ ].IsListEmpty() )
        {
            return true;
        }
    }
    return false;
}

void TaskScheduler::WaitForTasks( uint32_t threadNum )
{
    // We incrememt the number of threads waiting here in order
    // to ensure that the check for tasks occurs after the increment
    // to prevent a task being added after a check, then the thread waiting.
    // This will occasionally result in threads being mistakenly awoken,
    // but they will then go back to sleep.

    bool bHaveTasks = HaveTasks( threadNum );
    if( !bHaveTasks )
    {
        SafeCallback( m_ProfilerCallbacks.waitStart, threadNum );
        m_NumThreadsWaiting.fetch_add( 1, std::memory_order_acquire );
        SemaphoreWait( *m_pNewTaskSemaphore );
        SafeCallback( m_ProfilerCallbacks.waitStop, threadNum );
    }
}

void TaskScheduler::WakeThreads()
{
    int32_t waiting;
    do
    {
        waiting = m_NumThreadsWaiting;
    } while( waiting && !m_NumThreadsWaiting.compare_exchange_weak(waiting, 0, std::memory_order_relaxed ) );

    if( waiting )
    {
        SemaphoreSignal( *m_pNewTaskSemaphore, waiting );
    }
}

void TaskScheduler::SplitAndAddTask( uint32_t threadNum_, SubTaskSet subTask_, uint32_t rangeToSplit_ )
{
    int32_t numAdded = 0;
    while( subTask_.partition.start != subTask_.partition.end )
    {
        SubTaskSet taskToAdd = SplitTask( subTask_, rangeToSplit_ );

        // add the partition to the pipe
        ++numAdded;
        subTask_.pTask->m_RunningCount.fetch_add( 1, std::memory_order_acquire );
        if( !m_pPipesPerThread[ subTask_.pTask->m_Priority ][ threadNum_ ].WriterTryWriteFront( taskToAdd ) )
        {
            if( numAdded > 1 )
            {
                WakeThreads();
            }
            numAdded = 0;
            // alter range to run the appropriate fraction
            if( taskToAdd.pTask->m_RangeToRun < rangeToSplit_ )
            {
                taskToAdd.partition.end = taskToAdd.partition.start + taskToAdd.pTask->m_RangeToRun;
                subTask_.partition.start = taskToAdd.partition.end;
            }
            taskToAdd.pTask->ExecuteRange( taskToAdd.partition, threadNum_ );
            subTask_.pTask->m_RunningCount.fetch_sub( 1, std::memory_order_release );
        }
    }

    WakeThreads();
}

void    TaskScheduler::AddTaskSetToPipe( ITaskSet* pTaskSet_ )
{
    assert( pTaskSet_->m_RunningCount == 0 );
    pTaskSet_->m_RunningCount.store( 0, std::memory_order_relaxed );

    // divide task up and add to pipe
    pTaskSet_->m_RangeToRun = pTaskSet_->m_SetSize / m_NumPartitions;
    if( pTaskSet_->m_RangeToRun < pTaskSet_->m_MinRange ) { pTaskSet_->m_RangeToRun = pTaskSet_->m_MinRange; }

    uint32_t rangeToSplit = pTaskSet_->m_SetSize / m_NumInitialPartitions;
    if( rangeToSplit < pTaskSet_->m_MinRange ) { rangeToSplit = pTaskSet_->m_MinRange; }

    SubTaskSet subTask;
    subTask.pTask = pTaskSet_;
    subTask.partition.start = 0;
    subTask.partition.end = pTaskSet_->m_SetSize;
    SplitAndAddTask( gtl_threadNum, subTask, rangeToSplit );
}

void TaskScheduler::AddPinnedTask( IPinnedTask* pTask_ )
{
    assert( pTask_->m_RunningCount == 0 );

    pTask_->m_RunningCount = 1;
    m_pPinnedTaskListPerThread[ pTask_->m_Priority ][ pTask_->threadNum ].WriterWriteFront( pTask_ );
    WakeThreads();
}

void TaskScheduler::RunPinnedTasks()
{
    uint32_t threadNum = gtl_threadNum;
    for( int priority = 0; priority < TASK_PRIORITY_NUM; ++priority )
    {
        RunPinnedTasks( threadNum, priority );
    }
}

void TaskScheduler::RunPinnedTasks( uint32_t threadNum_, uint32_t priority_ )
{
    IPinnedTask* pPinnedTaskSet = NULL;
    do
    {
        pPinnedTaskSet = m_pPinnedTaskListPerThread[ priority_ ][ threadNum_ ].ReaderReadBack();
        if( pPinnedTaskSet )
        {
            pPinnedTaskSet->Execute();
            pPinnedTaskSet->m_RunningCount = 0;
        }
    } while( pPinnedTaskSet );
}

void    TaskScheduler::WaitforTask( const ICompletable* pCompletable_, enki::TaskPriority priorityOfLowestToRun_ )
{
    uint32_t hintPipeToCheck_io = gtl_threadNum + 1;    // does not need to be clamped.

    if( pCompletable_ )
    {
        // We need to ensure that the task we're waiting on can complete even if we're the only thread,
        // so we clamp the priorityOfLowestToRun_ to no smaller than the task we're waiting for
        priorityOfLowestToRun_ = std::max( priorityOfLowestToRun_, pCompletable_->m_Priority );
        while( !pCompletable_->GetIsComplete() )
        {
            for( int priority = 0; priority <= priorityOfLowestToRun_; ++priority )
            {
                if( TryRunTask( gtl_threadNum, priority, hintPipeToCheck_io ) )
                {
                    break;
                }
            }
        }
    }
    else
    {
            for( int priority = 0; priority <= priorityOfLowestToRun_; ++priority )
            {
                if( TryRunTask( gtl_threadNum, priority, hintPipeToCheck_io ) )
                {
                    break;
                }
            }
    }
}

void    TaskScheduler::WaitforAll()
{
    bool bHaveTasks = true;
    uint32_t hintPipeToCheck_io = gtl_threadNum  + 1;    // does not need to be clamped.
    int32_t numThreadsRunning = m_NumThreadsRunning.load( std::memory_order_relaxed ) - 1; // account for this thread
    while( bHaveTasks || m_NumThreadsWaiting.load( std::memory_order_relaxed ) < numThreadsRunning )
    {
        bHaveTasks = TryRunTask( gtl_threadNum, hintPipeToCheck_io );
     }
}

void    TaskScheduler::WaitforAllAndShutdown()
{
    WaitforAll();
    StopThreads(true);
}

uint32_t        TaskScheduler::GetNumTaskThreads() const
{
    return m_NumThreads;
}


uint32_t TaskScheduler::GetThreadNum() const
{
    return gtl_threadNum;
}


TaskScheduler::TaskScheduler()
        : m_pPipesPerThread()
        , m_pPinnedTaskListPerThread()
        , m_NumThreads(0)
        , m_pThreadArgStore(NULL)
        , m_pThreads(NULL)
        , m_bRunning(0)
        , m_NumThreadsRunning(0)
        , m_NumThreadsWaiting(0)
        , m_NumPartitions(0)
        , m_bHaveThreads(false)
{
}

TaskScheduler::~TaskScheduler()
{
    StopThreads( true ); // Stops threads, waiting for them.
}

void    TaskScheduler::Initialize( uint32_t numThreads_ )
{
    assert( numThreads_ );
    StopThreads( true ); // Stops threads, waiting for them.

    m_NumThreads = numThreads_;

    StartThreads();
}

void   TaskScheduler::Initialize()
{
    Initialize( std::thread::hardware_concurrency() );
}



// Semaphore implementation
#ifdef _WIN32

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <Windows.h>

namespace enki
{
    struct semaphoreid_t
    {
        HANDLE      sem;
    };
    
    inline void SemaphoreCreate( semaphoreid_t& semaphoreid )
    {
        semaphoreid.sem = CreateSemaphore(NULL, 0, MAXLONG, NULL );
    }

    inline void SemaphoreClose( semaphoreid_t& semaphoreid )
    {
        CloseHandle( semaphoreid.sem );
    }

    inline void SemaphoreWait( semaphoreid_t& semaphoreid  )
    {
        DWORD retval = WaitForSingleObject( semaphoreid.sem, INFINITE );

        assert( retval != WAIT_FAILED );
    }

    inline void SemaphoreSignal( semaphoreid_t& semaphoreid, int32_t countWaiting )
    {
        if( countWaiting )
        {
            ReleaseSemaphore( semaphoreid.sem, countWaiting, NULL );
        }
    }
}
#elif defined(__MACH__)

// OS X does not have POSIX semaphores
// see https://developer.apple.com/library/content/documentation/Darwin/Conceptual/KernelProgramming/synchronization/synchronization.html
#include <mach/mach.h>

namespace enki
{
    
    struct semaphoreid_t
    {
        semaphore_t   sem;
    };
    
    inline void SemaphoreCreate( semaphoreid_t& semaphoreid )
    {
        semaphore_create( mach_task_self(), &semaphoreid.sem, SYNC_POLICY_FIFO, 0 );
    }
    
    inline void SemaphoreClose( semaphoreid_t& semaphoreid )
    {
        semaphore_destroy( mach_task_self(), semaphoreid.sem );
    }
    
    inline void SemaphoreWait( semaphoreid_t& semaphoreid  )
    {
        semaphore_wait( semaphoreid.sem );
    }
    
    inline void SemaphoreSignal( semaphoreid_t& semaphoreid, int32_t countWaiting )
    {
        while( countWaiting-- > 0 )
        {
            semaphore_signal( semaphoreid.sem );
        }
    }
}

#else // POSIX

#include <semaphore.h>

namespace enki
{
    
    struct semaphoreid_t
    {
        sem_t   sem;
    };
    
    inline void SemaphoreCreate( semaphoreid_t& semaphoreid )
    {
        int err = sem_init( &semaphoreid.sem, 0, 0 );
        assert( err == 0 );
    }
    
    inline void SemaphoreClose( semaphoreid_t& semaphoreid )
    {
        sem_destroy( &semaphoreid.sem );
    }
    
    inline void SemaphoreWait( semaphoreid_t& semaphoreid  )
    {
        int err = sem_wait( &semaphoreid.sem );
        assert( err == 0 );
    }
    
    inline void SemaphoreSignal( semaphoreid_t& semaphoreid, int32_t countWaiting )
    {
        while( countWaiting-- > 0 )
        {
            sem_post( &semaphoreid.sem );
        }
    }
}
#endif

namespace enki
{
    semaphoreid_t* SemaphoreCreate()
    {
        semaphoreid_t* pSemaphore = new semaphoreid_t;
        SemaphoreCreate( *pSemaphore );
        return pSemaphore;
    }

    void SemaphoreDelete( semaphoreid_t* pSemaphore_ )
    {
        SemaphoreClose( *pSemaphore_ );
        delete pSemaphore_;
    }
}