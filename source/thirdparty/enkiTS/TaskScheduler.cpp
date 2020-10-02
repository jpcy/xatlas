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

#if defined(ENKI_CUSTOM_ALLOC_FILE_AND_LINE)
#define ENKI_FILE_AND_LINE __FILE__, __LINE__
#else
namespace
{
    const char* gc_File    = "";
    const uint32_t gc_Line = 0;
}
#define ENKI_FILE_AND_LINE  gc_File, gc_Line
#endif

namespace enki
{
    static const uint32_t gc_PipeSizeLog2            = 8;
    static const uint32_t gc_SpinCount               = 10;
    static const uint32_t gc_SpinBackOffMulitplier   = 100;
    static const uint32_t gc_MaxNumInitialPartitions = 8;
    static const uint32_t gc_CacheLineSize           = 64;
    // awaiting std::hardware_constructive_interference_size
};

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
    class TaskPipe : public LockLessMultiReadPipe<gc_PipeSizeLog2,enki::SubTaskSet> {};

    enum ThreadState : int32_t
    {
        THREAD_STATE_NONE,                  // shouldn't get this value
        THREAD_STATE_NOT_LAUNCHED,          // for debug purposes - indicates enki task thread not yet launched
        THREAD_STATE_RUNNING,
        THREAD_STATE_WAIT_TASK_COMPLETION,
        THREAD_STATE_EXTERNAL_REGISTERED,
        THREAD_STATE_EXTERNAL_UNREGISTERED,
        THREAD_STATE_WAIT_NEW_TASKS,
        THREAD_STATE_STOPPED,
    };

    struct ThreadArgs
    {
        uint32_t                 threadNum;
        TaskScheduler*           pTaskScheduler;
    };

    struct alignas(enki::gc_CacheLineSize) ThreadDataStore 
    {
        std::atomic<ThreadState> threadState = { THREAD_STATE_NONE };
        char prevent_false_Share[ enki::gc_CacheLineSize - sizeof(std::atomic<ThreadState>) ];
    };
    static_assert( sizeof( ThreadDataStore ) >= enki::gc_CacheLineSize, "ThreadDataStore may exhibit false sharing" );

    class PinnedTaskList : public LocklessMultiWriteIntrusiveList<IPinnedTask> {};

    semaphoreid_t* SemaphoreCreate();
    void SemaphoreDelete( semaphoreid_t* pSemaphore_ );
    void SemaphoreWait(   semaphoreid_t& semaphoreid );
    void SemaphoreSignal( semaphoreid_t& semaphoreid, int32_t countWaiting );
}

namespace
{
    SubTaskSet SplitTask( SubTaskSet& subTask_, uint32_t rangeToSplit_ )
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
    // Note: see https://software.intel.com/en-us/articles/a-common-construct-to-avoid-the-contention-of-threads-architecture-agnostic-spin-wait-loops
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

static void SafeCallback( ProfilerCallbackFunc func_, uint32_t threadnum_ )
{
    if( func_ != nullptr )
    {
        func_( threadnum_ );
    }
}

   
ENKITS_API void* enki::DefaultAllocFunc( size_t align_, size_t size_, void* userData_, const char* file_, int line_ )
{ 
    (void)userData_; (void)file_; (void)line_;
    void* pRet;
#ifdef _WIN32
    pRet = (void*)_aligned_malloc( size_, align_ );
#else
    int retval = posix_memalign( &pRet, align_, size_ );
    (void)retval;	//unused
#endif
    return pRet;
};

ENKITS_API void  enki::DefaultFreeFunc(  void* ptr_,   size_t size_, void* userData_, const char* file_, int line_ )
{
     (void)size_; (void)userData_; (void)file_; (void)line_;
#ifdef _WIN32
    _aligned_free( ptr_ );
#else
    free( ptr_ );
#endif
};

bool TaskScheduler::RegisterExternalTaskThread()
{
    bool bRegistered = false;
    while( !bRegistered && m_NumExternalTaskThreadsRegistered < (int32_t)m_Config.numExternalTaskThreads  )
    {
        for(uint32_t thread = 1; thread <= m_Config.numExternalTaskThreads; ++thread )
        {
            // ignore our thread
            ThreadState threadStateExpected = THREAD_STATE_EXTERNAL_UNREGISTERED;
            if( m_pThreadDataStore[thread].threadState.compare_exchange_strong(
                threadStateExpected, THREAD_STATE_EXTERNAL_REGISTERED ) )
            {
                ++m_NumExternalTaskThreadsRegistered;
                gtl_threadNum = thread;
                bRegistered = true;
                break;
            }
        }
    }
    return bRegistered;
}

void TaskScheduler::DeRegisterExternalTaskThread()
{
    assert( gtl_threadNum );
    ThreadState threadState = m_pThreadDataStore[gtl_threadNum].threadState.load( std::memory_order_acquire );
    assert( threadState == THREAD_STATE_EXTERNAL_REGISTERED );
    if( threadState == THREAD_STATE_EXTERNAL_REGISTERED )
    {
        --m_NumExternalTaskThreadsRegistered;
        m_pThreadDataStore[gtl_threadNum].threadState.store( THREAD_STATE_EXTERNAL_UNREGISTERED, std::memory_order_release );
        gtl_threadNum = 0;
    }
}

uint32_t TaskScheduler::GetNumRegisteredExternalTaskThreads()
{
    return m_NumExternalTaskThreadsRegistered;
}


void TaskScheduler::TaskingThreadFunction( const ThreadArgs& args_ )
{
    uint32_t threadNum  = args_.threadNum;
    TaskScheduler*  pTS = args_.pTaskScheduler;
    gtl_threadNum       = threadNum;

    pTS->m_pThreadDataStore[threadNum].threadState.store( THREAD_STATE_RUNNING, std::memory_order_release );
    SafeCallback( pTS->m_Config.profilerCallbacks.threadStart, threadNum );

    uint32_t spinCount = 0;
    uint32_t hintPipeToCheck_io = threadNum + 1;    // does not need to be clamped.
    while( pTS->m_bRunning.load( std::memory_order_relaxed ) )
    {
        if( !pTS->TryRunTask( threadNum, hintPipeToCheck_io ) )
        {
            // no tasks, will spin then wait
            ++spinCount;
            if( spinCount > gc_SpinCount )
            {
                pTS->WaitForNewTasks( threadNum );
            }
            else
            {
                uint32_t spinBackoffCount = spinCount * gc_SpinBackOffMulitplier;
                SpinWait( spinBackoffCount );
            }
        }
        else
        {
            spinCount = 0; // have run a task so reset spin count.
        }
    }

    pTS->m_NumInternalTaskThreadsRunning.fetch_sub( 1, std::memory_order_release );
    pTS->m_pThreadDataStore[threadNum].threadState.store( THREAD_STATE_STOPPED, std::memory_order_release );
    SafeCallback( pTS->m_Config.profilerCallbacks.threadStop, threadNum );
    return;

}


void TaskScheduler::StartThreads()
{
    if( m_bHaveThreads )
    {
        return;
    }

    m_NumThreads = m_Config.numTaskThreadsToCreate + m_Config.numExternalTaskThreads + 1;

    for( int priority = 0; priority < TASK_PRIORITY_NUM; ++priority )
    {
        m_pPipesPerThread[ priority ]          = NewArray<TaskPipe>( m_NumThreads, ENKI_FILE_AND_LINE );
        m_pPinnedTaskListPerThread[ priority ] = NewArray<PinnedTaskList>( m_NumThreads, ENKI_FILE_AND_LINE );
    }

    m_pNewTaskSemaphore      = SemaphoreNew();
    m_pTaskCompleteSemaphore = SemaphoreNew();

    // we create one less thread than m_NumThreads as the main thread counts as one
    m_pThreadDataStore   = NewArray<ThreadDataStore>( m_NumThreads, ENKI_FILE_AND_LINE );
    m_pThreads           = NewArray<std::thread>( m_NumThreads, ENKI_FILE_AND_LINE );
    m_bRunning = 1;

    for( uint32_t thread = 0; thread < m_Config.numExternalTaskThreads + 1; ++thread )
    {
        m_pThreadDataStore[thread].threadState   = THREAD_STATE_EXTERNAL_UNREGISTERED;
    }
    for( uint32_t thread = m_Config.numExternalTaskThreads + 1; thread < m_NumThreads; ++thread )
    {
        m_pThreadDataStore[thread].threadState   = THREAD_STATE_NOT_LAUNCHED;
    }
    // only launch threads once all thread states are set
    for( uint32_t thread = m_Config.numExternalTaskThreads + 1; thread < m_NumThreads; ++thread )
    {
        m_pThreads[thread]                       = std::thread( TaskingThreadFunction, ThreadArgs{ thread, this } );
        ++m_NumInternalTaskThreadsRunning;
    }

    // ensure we have sufficient tasks to equally fill either all threads including main
    // or just the threads we've launched, this is outside the firstinit as we want to be able
    // to runtime change it
    if( 1 == m_NumThreads )
    {
        m_NumPartitions        = 1;
        m_NumInitialPartitions = 1;
    }
    else
    {
        // There could be more threads than hardware threads if external threads are
        // being intended for blocking functionality such as io etc.
        // We only need to partition for a maximum of the available processor parallelism.
        uint32_t numThreadsToPartitionFor = std::min( m_NumThreads, GetNumHardwareThreads() );
        m_NumPartitions = numThreadsToPartitionFor * (numThreadsToPartitionFor - 1);
        m_NumInitialPartitions = numThreadsToPartitionFor - 1;
        if( m_NumInitialPartitions > gc_MaxNumInitialPartitions )
        {
            m_NumInitialPartitions = gc_MaxNumInitialPartitions;
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
        while( bWait_ && m_NumInternalTaskThreadsRunning )
        {
            // keep firing event to ensure all threads pick up state of m_bRunning
           WakeThreadsForNewTasks();
        }

        // detach threads starting with thread 1 (as 0 is initialization thread).
        for( uint32_t thread = m_Config.numExternalTaskThreads + 1; thread < m_NumThreads; ++thread )
        {
            assert( m_pThreads[thread].joinable() );
            m_pThreads[thread].join();
        }

        DeleteArray( m_pThreadDataStore, m_NumThreads, ENKI_FILE_AND_LINE );
        DeleteArray( m_pThreads, m_NumThreads, ENKI_FILE_AND_LINE );
        m_pThreadDataStore = 0;
        m_pThreads = 0;

        SemaphoreDelete( m_pNewTaskSemaphore );
        m_pNewTaskSemaphore = 0;
        SemaphoreDelete( m_pTaskCompleteSemaphore );
        m_pTaskCompleteSemaphore = 0;

        m_bHaveThreads = false;
        m_NumThreadsWaitingForNewTasks = 0;
        m_NumThreadsWaitingForTaskCompletion = 0;
        m_NumInternalTaskThreadsRunning = 0;
        m_NumExternalTaskThreadsRegistered = 0;

        for( int priority = 0; priority < TASK_PRIORITY_NUM; ++priority )
        {
            DeleteArray( m_pPipesPerThread[ priority ], m_NumThreads, ENKI_FILE_AND_LINE );
            m_pPipesPerThread[ priority ] = NULL;
            DeleteArray( m_pPinnedTaskListPerThread[ priority ], m_NumThreads, ENKI_FILE_AND_LINE );
            m_pPinnedTaskListPerThread[ priority ] = NULL;
        }
        m_NumThreads = 0;
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

bool TaskScheduler::TryRunTask( uint32_t threadNum_, uint32_t priority_, uint32_t& hintPipeToCheck_io_ )
{
    // Run any tasks for this thread
    RunPinnedTasks( threadNum_, priority_ );

    // check for tasks
    SubTaskSet subTask;
    bool bHaveTask = m_pPipesPerThread[ priority_ ][ threadNum_ ].WriterTryReadFront( &subTask );

    uint32_t threadToCheck = hintPipeToCheck_io_;
    uint32_t checkCount = 0;
    while( !bHaveTask && checkCount < m_NumThreads )
    {
        threadToCheck = ( hintPipeToCheck_io_ + checkCount ) % m_NumThreads;
        if( threadToCheck != threadNum_ )
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
            SplitAndAddTask( threadNum_, subTask, subTask.pTask->m_RangeToRun );
            taskToRun.pTask->ExecuteRange( taskToRun.partition, threadNum_ );
            int prevCount = taskToRun.pTask->m_RunningCount.fetch_sub(1,std::memory_order_release );
            if( 1 == prevCount )
            {
                TaskComplete( taskToRun.pTask, true, threadNum_ );
            }
        }
        else
        {
            // the task has already been divided up by AddTaskSetToPipe, so just run it
            subTask.pTask->ExecuteRange( subTask.partition, threadNum_ );
            int prevCount = subTask.pTask->m_RunningCount.fetch_sub(1,std::memory_order_release );
            if( 1 == prevCount )
            {
                TaskComplete( subTask.pTask, true, threadNum_ );
            }
        }
    }

    return bHaveTask;

}

void TaskScheduler::TaskComplete( ICompletable* pTask_, bool bWakeThreads_, uint32_t threadNum_ )
{
    assert( pTask_->GetIsComplete() );
    if( bWakeThreads_ && pTask_->m_WaitingForTaskCount.load( std::memory_order_acquire ) )
    {
        WakeThreadsForTaskCompletion();
    }

    Dependency* pDependent = pTask_->m_pDependents;
    while( pDependent )
    {
        int prevDeps = pDependent->pTaskToRunOnCompletion->m_DependenciesCompletedCount.fetch_add( 1, std::memory_order_release );
        assert( prevDeps < pDependent->pTaskToRunOnCompletion->m_DependenciesCount );
        if( pDependent->pTaskToRunOnCompletion->m_DependenciesCount == ( prevDeps + 1 ) )
        {
            // get temp copy of pDependent so OnDependenciesComplete can delete task if needed.
            Dependency* pDependentCurr = pDependent;
            pDependent = pDependent->pNext;
            // reset dependencies
            pDependentCurr->pTaskToRunOnCompletion->m_DependenciesCompletedCount.store(
                0,
                std::memory_order_release );
            pDependentCurr->pTaskToRunOnCompletion->OnDependenciesComplete( this, threadNum_ );
        }
        else
        {
            pDependent = pDependent->pNext;
        }
    }
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

void TaskScheduler::WaitForNewTasks( uint32_t threadNum_ )
{
    // We don't want to suspend this thread if there are task threads
    // with pinned tasks suspended, as it could result in this thread
    // being unsuspended and not the thread with pinned tasks
    if( WakeSuspendedThreadsWithPinnedTasks() )
    {
        return;
    }

    // We incrememt the number of threads waiting here in order
    // to ensure that the check for tasks occurs after the increment
    // to prevent a task being added after a check, then the thread waiting.
    // This will occasionally result in threads being mistakenly awoken,
    // but they will then go back to sleep.
    m_NumThreadsWaitingForNewTasks.fetch_add( 1, std::memory_order_acquire );
    ThreadState prevThreadState = m_pThreadDataStore[threadNum_].threadState.load( std::memory_order_relaxed );
    m_pThreadDataStore[threadNum_].threadState.store( THREAD_STATE_WAIT_NEW_TASKS, std::memory_order_seq_cst );

    if( HaveTasks( threadNum_ ) )
    {
        m_NumThreadsWaitingForNewTasks.fetch_sub( 1, std::memory_order_release );
    }
    else
    {
        SafeCallback( m_Config.profilerCallbacks.waitForNewTaskSuspendStart, threadNum_ );
        SemaphoreWait( *m_pNewTaskSemaphore );
        SafeCallback( m_Config.profilerCallbacks.waitForNewTaskSuspendStop, threadNum_ );
    }

    m_pThreadDataStore[threadNum_].threadState.store( prevThreadState, std::memory_order_release );
}

void TaskScheduler::WaitForTaskCompletion( const ICompletable* pCompletable_, uint32_t threadNum_ )
{
    // We don't want to suspend this thread if there are task threads
    // with pinned tasks suspended, as the completable could be a pinned task
    // or it could be waiting on one.
    if( WakeSuspendedThreadsWithPinnedTasks() )
    {
        return;
    }

    m_NumThreadsWaitingForTaskCompletion.fetch_add( 1, std::memory_order_acquire );
    pCompletable_->m_WaitingForTaskCount.fetch_add( 1, std::memory_order_acquire );
    ThreadState prevThreadState = m_pThreadDataStore[threadNum_].threadState.load( std::memory_order_relaxed );
    m_pThreadDataStore[threadNum_].threadState.store( THREAD_STATE_WAIT_TASK_COMPLETION, std::memory_order_seq_cst );

    if( pCompletable_->GetIsComplete() || HaveTasks( threadNum_ ) )
    {
        m_NumThreadsWaitingForTaskCompletion.fetch_sub( 1, std::memory_order_release );
    }
    else
    {
        SafeCallback( m_Config.profilerCallbacks.waitForTaskCompleteSuspendStart, threadNum_ );
        std::atomic_thread_fence(std::memory_order_acquire);

        SemaphoreWait( *m_pTaskCompleteSemaphore );
        if( !pCompletable_->GetIsComplete() )
        {
            // This thread which may not the one which was supposed to be awoken
            WakeThreadsForTaskCompletion();
        }
        SafeCallback( m_Config.profilerCallbacks.waitForTaskCompleteSuspendStop, threadNum_ );
    }

    m_pThreadDataStore[threadNum_].threadState.store( prevThreadState, std::memory_order_release );
    pCompletable_->m_WaitingForTaskCount.fetch_sub( 1, std::memory_order_release );
}

void TaskScheduler::WakeThreadsForNewTasks()
{
    int32_t waiting = m_NumThreadsWaitingForNewTasks.load( std::memory_order_relaxed );
    while( waiting > 0 && !m_NumThreadsWaitingForNewTasks.compare_exchange_weak(waiting, 0, std::memory_order_release, std::memory_order_relaxed ) ) {}

    if( waiting > 0 )
    {
        SemaphoreSignal( *m_pNewTaskSemaphore, waiting );
    }

    // We also wake tasks waiting for completion as they can run tasks
    WakeThreadsForTaskCompletion();
}

void TaskScheduler::WakeThreadsForTaskCompletion()
{
    // m_NumThreadsWaitingForTaskCompletion can go negative as this indicates that
    // we signalled more threads than the number which ended up waiting
    int32_t waiting = m_NumThreadsWaitingForTaskCompletion.load( std::memory_order_relaxed );
    while( waiting > 0 && !m_NumThreadsWaitingForTaskCompletion.compare_exchange_weak(waiting, 0, std::memory_order_release, std::memory_order_relaxed ) ) {}

    if( waiting > 0 )
    {
        SemaphoreSignal( *m_pTaskCompleteSemaphore, waiting );
    }
}

bool TaskScheduler::WakeSuspendedThreadsWithPinnedTasks()
{
    uint32_t threadNum = gtl_threadNum;
    for( uint32_t t = 1; t < m_NumThreads; ++t )
    {
        // distribute thread checks more evenly by starting at our thread number rather than 0.
        uint32_t thread = ( threadNum + t ) % m_NumThreads;

        ThreadState state = m_pThreadDataStore[ thread ].threadState.load( std::memory_order_acquire );
            
        assert( state != THREAD_STATE_NONE );

        if( state == THREAD_STATE_WAIT_NEW_TASKS || state == THREAD_STATE_WAIT_TASK_COMPLETION )
        {
            // thread is suspended, check if it has pinned tasks
            for( int priority = 0; priority < TASK_PRIORITY_NUM; ++priority )
            {
                if( !m_pPinnedTaskListPerThread[ priority ][ thread ].IsListEmpty() )
                {
                    WakeThreadsForNewTasks();
                    return true;
                }
            }
        }
    }
    return false;
}

void TaskScheduler::SplitAndAddTask( uint32_t threadNum_, SubTaskSet subTask_, uint32_t rangeToSplit_ )
{
    int32_t numAdded = 0;
    int32_t numRun   = 0;
    // ensure that an artificial completion is not registered whilst adding tasks by incrementing count
    subTask_.pTask->m_RunningCount.fetch_add( 1, std::memory_order_acquire );
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
                WakeThreadsForNewTasks();
            }
            numAdded = 0;
            // alter range to run the appropriate fraction
            if( taskToAdd.pTask->m_RangeToRun < rangeToSplit_ )
            {
                taskToAdd.partition.end = taskToAdd.partition.start + taskToAdd.pTask->m_RangeToRun;
                subTask_.partition.start = taskToAdd.partition.end;
            }
            taskToAdd.pTask->ExecuteRange( taskToAdd.partition, threadNum_ );
            ++numRun;
        }
    }
    int prevCount = subTask_.pTask->m_RunningCount.fetch_sub( numRun + 1, std::memory_order_release );
    if( numRun + 1 == prevCount )
    {
        TaskComplete( subTask_.pTask, false, threadNum_ );
    }

    // WakeThreadsForNewTasks also calls WakeThreadsForTaskCompletion() so do not need to do so above
    WakeThreadsForNewTasks();
}

TaskSchedulerConfig TaskScheduler::GetConfig() const
{
    return m_Config;
}

void TaskScheduler::AddTaskSetToPipeInt( ITaskSet* pTaskSet_, uint32_t threadNum_ )
{
    assert( pTaskSet_->m_RunningCount == 1 );
    ThreadState prevThreadState = m_pThreadDataStore[threadNum_].threadState.load( std::memory_order_relaxed );
    m_pThreadDataStore[threadNum_].threadState.store( THREAD_STATE_RUNNING, std::memory_order_relaxed );
    std::atomic_thread_fence(std::memory_order_acquire);


    // divide task up and add to pipe
    pTaskSet_->m_RangeToRun = pTaskSet_->m_SetSize / m_NumPartitions;
    if( pTaskSet_->m_RangeToRun < pTaskSet_->m_MinRange ) { pTaskSet_->m_RangeToRun = pTaskSet_->m_MinRange; }

    uint32_t rangeToSplit = pTaskSet_->m_SetSize / m_NumInitialPartitions;
    if( rangeToSplit < pTaskSet_->m_MinRange ) { rangeToSplit = pTaskSet_->m_MinRange; }

    SubTaskSet subTask;
    subTask.pTask = pTaskSet_;
    subTask.partition.start = 0;
    subTask.partition.end = pTaskSet_->m_SetSize;
    SplitAndAddTask( threadNum_, subTask, rangeToSplit );
    int prevCount = pTaskSet_->m_RunningCount.fetch_sub(1, std::memory_order_release );
    if( 1 == prevCount )
    {
        TaskComplete( pTaskSet_, true, threadNum_ );
    }

    m_pThreadDataStore[threadNum_].threadState.store( prevThreadState, std::memory_order_release );
}

void TaskScheduler::AddTaskSetToPipe( ITaskSet* pTaskSet_ )
{
    assert( pTaskSet_->m_RunningCount == 0 );
    InitDependencies( pTaskSet_ );
    pTaskSet_->m_RunningCount.store( 1, std::memory_order_relaxed );
    AddTaskSetToPipeInt( pTaskSet_, gtl_threadNum );
}

void  TaskScheduler::AddPinnedTaskInt( IPinnedTask* pTask_ )
{
    assert( pTask_->m_RunningCount == 1 );
    m_pPinnedTaskListPerThread[ pTask_->m_Priority ][ pTask_->threadNum ].WriterWriteFront( pTask_ );
    WakeThreadsForNewTasks();
}

void TaskScheduler::AddPinnedTask( IPinnedTask* pTask_ )
{
    assert( pTask_->m_RunningCount == 0 );
    InitDependencies( pTask_ );
    pTask_->m_RunningCount = 1;
    AddPinnedTaskInt( pTask_ );
}

void TaskScheduler::InitDependencies( ICompletable* pCompletable_ )
{
    // go through any dependencies and set thier running count so they show as not complete
    // and increment depedency count
    if( pCompletable_->m_RunningCount.load( std::memory_order_relaxed ) )
    {
        // already initialized
        return;
    }
    Dependency* pDependent = pCompletable_->m_pDependents;
    while( pDependent )
    {
        InitDependencies( pDependent->pTaskToRunOnCompletion );
        pDependent->pTaskToRunOnCompletion->m_RunningCount.store( 1, std::memory_order_relaxed );
        pDependent = pDependent->pNext;
    }
}


void TaskScheduler::RunPinnedTasks()
{
    uint32_t threadNum = gtl_threadNum;
    ThreadState prevThreadState = m_pThreadDataStore[threadNum].threadState.load( std::memory_order_relaxed );
    m_pThreadDataStore[threadNum].threadState.store( THREAD_STATE_RUNNING, std::memory_order_relaxed );
    std::atomic_thread_fence(std::memory_order_acquire);
    for( int priority = 0; priority < TASK_PRIORITY_NUM; ++priority )
    {
        RunPinnedTasks( threadNum, priority );
    }
    m_pThreadDataStore[threadNum].threadState.store( prevThreadState, std::memory_order_release );
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
            TaskComplete( pPinnedTaskSet, true, threadNum_ );
        }
    } while( pPinnedTaskSet );
}

void    TaskScheduler::WaitforTask( const ICompletable* pCompletable_, enki::TaskPriority priorityOfLowestToRun_ )
{
    uint32_t threadNum = gtl_threadNum;
    uint32_t hintPipeToCheck_io = threadNum + 1;    // does not need to be clamped.

    // waiting for a task is equivalent to 'running' for thread state purpose as we may run tasks whilst waiting
    ThreadState prevThreadState = m_pThreadDataStore[threadNum].threadState.load( std::memory_order_relaxed );
    m_pThreadDataStore[threadNum].threadState.store( THREAD_STATE_RUNNING, std::memory_order_relaxed );
    std::atomic_thread_fence(std::memory_order_acquire);


    if( pCompletable_ && !pCompletable_->GetIsComplete() )
    {
        SafeCallback( m_Config.profilerCallbacks.waitForTaskCompleteStart, threadNum );
        // We need to ensure that the task we're waiting on can complete even if we're the only thread,
        // so we clamp the priorityOfLowestToRun_ to no smaller than the task we're waiting for
        priorityOfLowestToRun_ = std::max( priorityOfLowestToRun_, pCompletable_->m_Priority );
        uint32_t spinCount = 0;
        while( !pCompletable_->GetIsComplete() )
        {
            ++spinCount;
            for( int priority = 0; priority <= priorityOfLowestToRun_; ++priority )
            {
                if( TryRunTask( threadNum, priority, hintPipeToCheck_io ) )
                {
                    spinCount = 0; // reset spin as ran a task
                    break;
                }
            }
            if( spinCount > gc_SpinCount )
            {
                WaitForTaskCompletion( pCompletable_, threadNum );
                spinCount = 0;
            }
            else
            {
                uint32_t spinBackoffCount = spinCount * gc_SpinBackOffMulitplier;
                SpinWait( spinBackoffCount );
            }
        }
        SafeCallback( m_Config.profilerCallbacks.waitForTaskCompleteStop, threadNum );
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

    m_pThreadDataStore[threadNum].threadState.store( prevThreadState, std::memory_order_release );

}

class TaskSchedulerWaitTask : public IPinnedTask
{
    void Execute() override
    {
        // do nothing
    }
};

void TaskScheduler::WaitforAll()
{
    bool bHaveTasks = true;
    uint32_t threadNum = gtl_threadNum;
    uint32_t hintPipeToCheck_io = threadNum  + 1;    // does not need to be clamped.
    int32_t numOtherThreadsRunning = 0; // account for this thread
    uint32_t spinCount = 0;
    TaskSchedulerWaitTask dummyWaitTask;
    dummyWaitTask.threadNum = 0;
    while( bHaveTasks || numOtherThreadsRunning )
    {
        bHaveTasks = TryRunTask( threadNum, hintPipeToCheck_io );
        ++spinCount;
        if( bHaveTasks )
        {
            spinCount = 0; // reset spin as ran a task
        }
        if( spinCount > gc_SpinCount )
        {
            // find a running thread and add a dummy wait task
            int32_t countThreadsToCheck = m_NumThreads - 1;
            bool bHaveThreadToWaitOn = false;
            do
            {
                --countThreadsToCheck;
                dummyWaitTask.threadNum = ( dummyWaitTask.threadNum + 1 ) % m_NumThreads;

                // We can only add a pinned task to wait on if we find an enki Task Thread which isn't this thread.
                // Otherwise we have to busy wait.
                if( dummyWaitTask.threadNum != threadNum && dummyWaitTask.threadNum > m_Config.numExternalTaskThreads )
                {
                    ThreadState state = m_pThreadDataStore[ dummyWaitTask.threadNum ].threadState.load( std::memory_order_acquire );
                    if( state == THREAD_STATE_RUNNING || state == THREAD_STATE_WAIT_TASK_COMPLETION )
                    {
                        bHaveThreadToWaitOn = true;
                        break;
                    }
                }
            } while( countThreadsToCheck );

            if( bHaveThreadToWaitOn )
            {
                assert( dummyWaitTask.threadNum != threadNum );
                AddPinnedTask( &dummyWaitTask );
                WaitforTask( &dummyWaitTask );
            }
            spinCount = 0;
        }
        else
        {
            uint32_t spinBackoffCount = spinCount * gc_SpinBackOffMulitplier;
            SpinWait( spinBackoffCount );
        }

        // count threads running
        numOtherThreadsRunning = 0;
        for(uint32_t thread = 0; thread < m_NumThreads; ++thread )
        {
            // ignore our thread
            if( thread != threadNum )
            {
                switch( m_pThreadDataStore[thread].threadState.load( std::memory_order_acquire ) )
                {
                case THREAD_STATE_NONE:
                    assert(false);
                    break;
                case THREAD_STATE_NOT_LAUNCHED:
                case THREAD_STATE_RUNNING:
                case THREAD_STATE_WAIT_TASK_COMPLETION:
                    ++numOtherThreadsRunning;
                    break;
                case THREAD_STATE_EXTERNAL_REGISTERED:
                case THREAD_STATE_EXTERNAL_UNREGISTERED:
                case THREAD_STATE_WAIT_NEW_TASKS:
                case THREAD_STATE_STOPPED:
                    break;
                 };
            }
        }
     }
}

void    TaskScheduler::WaitforAllAndShutdown()
{
    if( m_bHaveThreads )
    {
        WaitforAll();
        StopThreads(true);
    }
}

uint32_t        TaskScheduler::GetNumTaskThreads() const
{
    return m_NumThreads;
}


uint32_t TaskScheduler::GetThreadNum() const
{
    return gtl_threadNum;
}

template<typename T>
T* TaskScheduler::NewArray( size_t num_, const char* file_, int line_  )
{
    T* pRet = (T*)m_Config.customAllocator.alloc( alignof(T), num_*sizeof(T), m_Config.customAllocator.userData, file_, line_ );
    if( !std::is_pod<T>::value )
    {
		T* pCurr = pRet;
        for( size_t i = 0; i < num_; ++i )
        {
			void* pBuffer = pCurr;
            pCurr = new(pBuffer) T;
			++pCurr;
        }
    }
    return pRet;
}

template<typename T>
void TaskScheduler::DeleteArray( T* p_, size_t num_, const char* file_, int line_ )
{
    if( !std::is_pod<T>::value )
    {
        size_t i = num_;
        while(i)
        {
            p_[--i].~T();
        }
    }
    m_Config.customAllocator.free( p_, sizeof(T)*num_, m_Config.customAllocator.userData, file_, line_ );
}

template<class T, class... Args>
T* TaskScheduler::New( const char* file_, int line_, Args&&... args_ )
{
    T* pRet = (T*)m_Config.customAllocator.alloc( alignof(T), sizeof(T), m_Config.customAllocator.userData, file_, line_ );
    return new(pRet) T( std::forward<Args>(args_)... );
}

template< typename T >
void TaskScheduler::Delete( T* p_, const char* file_, int line_  )
{
    p_->~T(); 
    m_Config.customAllocator.free( p_, sizeof(T), m_Config.customAllocator.userData, file_, line_ );
}

TaskScheduler::TaskScheduler()
        : m_pPipesPerThread()
        , m_pPinnedTaskListPerThread()
        , m_NumThreads(0)
        , m_pThreadDataStore(NULL)
        , m_pThreads(NULL)
        , m_bRunning(0)
        , m_NumInternalTaskThreadsRunning(0)
        , m_NumThreadsWaitingForNewTasks(0)
        , m_NumThreadsWaitingForTaskCompletion(0)
        , m_NumPartitions(0)
        , m_pNewTaskSemaphore(NULL)
        , m_pTaskCompleteSemaphore(NULL)
        , m_NumInitialPartitions(0)
        , m_bHaveThreads(false)
        , m_NumExternalTaskThreadsRegistered(0)
{
}

TaskScheduler::~TaskScheduler()
{
    StopThreads( true ); // Stops threads, waiting for them.
}

void TaskScheduler::Initialize( uint32_t numThreadsTotal_ )
{
    assert( numThreadsTotal_ >= 1 );
    StopThreads( true ); // Stops threads, waiting for them.
    m_Config.numTaskThreadsToCreate = numThreadsTotal_ - 1;
    m_Config.numExternalTaskThreads = 0;
    StartThreads();}

void TaskScheduler::Initialize( TaskSchedulerConfig config_ )
{
    StopThreads( true ); // Stops threads, waiting for them.
    m_Config = config_;
    StartThreads();
}

void TaskScheduler::Initialize()
{
    Initialize( std::thread::hardware_concurrency() );
}

// Semaphore implementation
#ifdef _WIN32

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <Windows.h>

namespace enki
{
    struct semaphoreid_t
    {
        HANDLE      sem;
    };
    
    inline void SemaphoreCreate( semaphoreid_t& semaphoreid )
    {
#ifdef _XBOX_ONE
        semaphoreid.sem = CreateSemaphoreExW( NULL, 0, MAXLONG, NULL, 0, SEMAPHORE_ALL_ACCESS );
#else
        semaphoreid.sem = CreateSemaphore( NULL, 0, MAXLONG, NULL );
#endif
    }

    inline void SemaphoreClose( semaphoreid_t& semaphoreid )
    {
        CloseHandle( semaphoreid.sem );
    }

    inline void SemaphoreWait( semaphoreid_t& semaphoreid  )
    {
        DWORD retval = WaitForSingleObject( semaphoreid.sem, INFINITE );
        assert( retval != WAIT_FAILED );
        (void)retval; // only needed for assert
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
// Mach semaphores can now only be created by the kernel
// Named sempahores work, but would require unique name construction to ensure
// they are isolated to this process.
// Dispatch semaphores appear to be the way other developers use OSX Semaphores, e.g. Boost
// However the API could change
// OSX below 10.6 does not support dispatch, but I do not have an earlier OSX version
// to test alternatives
#include <dispatch/dispatch.h>

namespace enki
{
    
    struct semaphoreid_t
    {
        dispatch_semaphore_t   sem;
    };
    
    inline void SemaphoreCreate( semaphoreid_t& semaphoreid )
    {
        semaphoreid.sem = dispatch_semaphore_create(0);
    }
    
    inline void SemaphoreClose( semaphoreid_t& semaphoreid )
    {
        dispatch_release( semaphoreid.sem );
    }
    
    inline void SemaphoreWait( semaphoreid_t& semaphoreid  )
    {
        dispatch_semaphore_wait( semaphoreid.sem, DISPATCH_TIME_FOREVER );
    }
    
    inline void SemaphoreSignal( semaphoreid_t& semaphoreid, int32_t countWaiting )
    {
        while( countWaiting-- > 0 )
        {
            dispatch_semaphore_signal( semaphoreid.sem );
        }
    }
}

#else // POSIX

#include <semaphore.h>
#include <errno.h>

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
        while( sem_wait( &semaphoreid.sem ) == -1 && errno == EINTR ) {}
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

semaphoreid_t* TaskScheduler::SemaphoreNew()
{
    semaphoreid_t* pSemaphore = New<semaphoreid_t>( ENKI_FILE_AND_LINE );
    SemaphoreCreate( *pSemaphore );
    return pSemaphore;
}

void TaskScheduler::SemaphoreDelete( semaphoreid_t* pSemaphore_ )
{
    SemaphoreClose( *pSemaphore_ );
    Delete( pSemaphore_, ENKI_FILE_AND_LINE );
}

void TaskScheduler::SetCustomAllocator( CustomAllocator customAllocator_ )
{
    m_Config.customAllocator = customAllocator_;
}

Dependency::Dependency( const ICompletable* pDependencyTask_, ICompletable* pTaskToRunOnCompletion_ ) 
    : pDependencyTask( pDependencyTask_ )
    , pTaskToRunOnCompletion( pTaskToRunOnCompletion_ )
    , pNext( pDependencyTask->m_pDependents )
{
    assert( pDependencyTask->GetIsComplete() );
    assert( pTaskToRunOnCompletion->GetIsComplete() );
    pDependencyTask->m_pDependents = this;
    ++pTaskToRunOnCompletion->m_DependenciesCount;
}

Dependency::Dependency( Dependency&& rhs_ ) noexcept
{
    pDependencyTask   = rhs_.pDependencyTask;
    pTaskToRunOnCompletion = rhs_.pTaskToRunOnCompletion;
    pNext             = rhs_.pNext;
    if( rhs_.pDependencyTask )
    {
        assert( rhs_.pTaskToRunOnCompletion );
        assert( rhs_.pDependencyTask->GetIsComplete() );
        assert( rhs_.pTaskToRunOnCompletion->GetIsComplete() );
        Dependency** ppDependent = &(pDependencyTask->m_pDependents);
        while( *ppDependent )
        {
            if( &rhs_ == *ppDependent )
            {
                *ppDependent = this;
                break;
            }
            ppDependent = &((*ppDependent)->pNext);
        }
    }
}


Dependency::~Dependency()
{
    ClearDependency();
}

void Dependency::SetDependency( const ICompletable* pDependencyTask_, ICompletable* pTaskToRunOnCompletion_ )
{
    ClearDependency();
    assert( pDependencyTask_->GetIsComplete() );
    assert( pTaskToRunOnCompletion_->GetIsComplete() );
    pDependencyTask = pDependencyTask_;
    pTaskToRunOnCompletion = pTaskToRunOnCompletion_;
    pNext = pDependencyTask->m_pDependents;
    pDependencyTask->m_pDependents = this;
    ++pTaskToRunOnCompletion->m_DependenciesCount;
}

void Dependency::ClearDependency()
{
    if( pDependencyTask )
    {
        assert( pTaskToRunOnCompletion );
        assert( pDependencyTask->GetIsComplete() );
        assert( pTaskToRunOnCompletion->GetIsComplete() );
        Dependency* pDependent = pDependencyTask->m_pDependents;
        if( this == pDependent )
        {
            pDependencyTask->m_pDependents = pDependent->pNext;
        }
        else
        {
            while( pDependent )
            {
                Dependency* pPrev = pDependent;
                pDependent = pDependent->pNext;
                if( this == pDependent )
                {
                    pPrev->pNext = pDependent->pNext;
                    break;
                }
            }
        }
    }
    pDependencyTask = NULL;
    pDependencyTask =  NULL;
    pNext = NULL;
}
