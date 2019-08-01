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

#include <atomic>
#include <thread>
#include <condition_variable>
#include <stdint.h>
#include <functional>
#include <assert.h>

// ENKITS_TASK_PRIORITIES_NUM can be set from 1 to 5.
// 1 corresponds to effectively no priorities.
#ifndef ENKITS_TASK_PRIORITIES_NUM
    #define ENKITS_TASK_PRIORITIES_NUM 3
#endif

#if   defined(_WIN32) && defined(ENKITS_BUILD_DLL)
    // Building enkiTS as a DLL
    #define ENKITS_API __declspec(dllexport)
#elif defined(_WIN32) && defined(ENKITS_DLL)
    // Using enkiTS as a DLL
    #define ENKITS_API __declspec(dllimport)
#elif defined(__GNUC__) && defined(ENKITS_BUILD_DLL)
    // Building enkiTS as a shared library
    #define ENKITS_API __attribute__((visibility("default")))
#else
    #define ENKITS_API
#endif

namespace enki
{

    struct TaskSetPartition
    {
        uint32_t start;
        uint32_t end;
    };

    class  TaskScheduler;
    class  TaskPipe;
    class  PinnedTaskList;
    struct ThreadArgs;
    struct SubTaskSet;
    struct semaphoreid_t;


    enum TaskPriority
    {
        TASK_PRIORITY_HIGH   = 0,
#if ( ENKITS_TASK_PRIORITIES_NUM > 3 )
        TASK_PRIORITY_MED_HI,
#endif
#if ( ENKITS_TASK_PRIORITIES_NUM > 2 )
        TASK_PRIORITY_MED,
#endif
#if ( ENKITS_TASK_PRIORITIES_NUM > 4 )
        TASK_PRIORITY_MED_LO,
#endif 
#if ( ENKITS_TASK_PRIORITIES_NUM > 1 )
        TASK_PRIORITY_LOW,
#endif
        TASK_PRIORITY_NUM
    };

    // ICompletable is a base class used to check for completion.
    // Do not use this class directly, instead derive from ITaskSet or IPinnedTask.
    class ICompletable
    {
    public:
        ICompletable() :        m_Priority(TASK_PRIORITY_HIGH), m_RunningCount(0) {}
        bool                    GetIsComplete() const {
            return 0 == m_RunningCount.load( std::memory_order_acquire );
        }

        virtual                ~ICompletable() {}

        TaskPriority            m_Priority;
    private:
        friend class            TaskScheduler;
        std::atomic<int32_t>   m_RunningCount;
    };

    // Subclass ITaskSet to create tasks.
    // TaskSets can be re-used, but check completion first.
    class ITaskSet : public ICompletable
    {
    public:
        ITaskSet()
            : m_SetSize(1)
            , m_MinRange(1)
            , m_RangeToRun(1)
        {}

        ITaskSet( uint32_t setSize_ )
            : m_SetSize( setSize_ )
            , m_MinRange(1)
            , m_RangeToRun(1)
        {}

        ITaskSet( uint32_t setSize_, uint32_t minRange_ )
            : m_SetSize( setSize_ )
            , m_MinRange( minRange_ )
            , m_RangeToRun(minRange_)
        {}

        // Execute range should be overloaded to process tasks. It will be called with a
        // range_ where range.start >= 0; range.start < range.end; and range.end < m_SetSize;
        // The range values should be mapped so that linearly processing them in order is cache friendly
        // i.e. neighbouring values should be close together.
        // threadnum should not be used for changing processing of data, it's intended purpose
        // is to allow per-thread data buckets for output.
        virtual void            ExecuteRange( TaskSetPartition range, uint32_t threadnum  ) = 0;

        // Size of set - usually the number of data items to be processed, see ExecuteRange. Defaults to 1
        uint32_t                m_SetSize;

        // Minimum size of of TaskSetPartition range when splitting a task set into partitions.
        // This should be set to a value which results in computation effort of at least 10k
        // clock cycles to minimize tast scheduler overhead.
        // NOTE: The last partition will be smaller than m_MinRange if m_SetSize is not a multiple
        // of m_MinRange.
        // Also known as grain size in literature.
        uint32_t                m_MinRange;

    private:
        friend class           TaskScheduler;
        uint32_t               m_RangeToRun;
    };

    // Subclass IPinnedTask to create tasks which can be run on a given thread only.
    class IPinnedTask : public ICompletable
    {
    public:
        IPinnedTask()                      : threadNum(0), pNext(NULL) {}  // default is to run a task on main thread
        IPinnedTask( uint32_t threadNum_ ) : threadNum(threadNum_), pNext(NULL) {}  // default is to run a task on main thread


        // IPinnedTask needs to be non abstract for intrusive list functionality.
        // Should never be called as should be overridden.
        virtual void            Execute() { assert(false); }


        uint32_t                 threadNum; // thread to run this pinned task on
        std::atomic<IPinnedTask*> pNext;        // Do not use. For intrusive list only.
    };

    // A utility task set for creating tasks based on std::func.
    typedef std::function<void (TaskSetPartition range, uint32_t threadnum  )> TaskSetFunction;
    class TaskSet : public ITaskSet
    {
    public:
        TaskSet() = default;
        TaskSet( TaskSetFunction func_ ) : m_Function( func_ ) {}
        TaskSet( uint32_t setSize_, TaskSetFunction func_ ) : ITaskSet( setSize_ ), m_Function( func_ ) {}


        virtual void            ExecuteRange( TaskSetPartition range, uint32_t threadnum  )
        {
            m_Function( range, threadnum );
        }

        TaskSetFunction m_Function;
    };

    // TaskScheduler implements several callbacks intended for profilers
    typedef std::function<void ( uint32_t threadnum_ )> ProfilerCallbackFunc;
    struct ProfilerCallbacks
    {
        ProfilerCallbackFunc threadStart;
        ProfilerCallbackFunc threadStop;
        ProfilerCallbackFunc waitStart;
        ProfilerCallbackFunc waitStop;
    };

    class TaskScheduler
    {
    public:
        ENKITS_API TaskScheduler();
        ENKITS_API ~TaskScheduler();

        // Call either Initialize() or Initialize( numThreads_ ) before adding tasks.

        // Initialize() will create GetNumHardwareThreads()-1 threads, which is
        // sufficient to fill the system when including the main thread.
        // Initialize can be called multiple times - it will wait for completion
        // before re-initializing.
        ENKITS_API void            Initialize();

        // Initialize( numThreads_ ) - numThreads_ (must be > 0)
        // will create numThreads_-1 threads, as thread 0 is
        // the thread on which the initialize was called.
        ENKITS_API void            Initialize( uint32_t numThreads_ );


        // Adds the TaskSet to pipe and returns if the pipe is not full.
        // If the pipe is full, pTaskSet is run.
        // should only be called from main thread, or within a task
        ENKITS_API void            AddTaskSetToPipe( ITaskSet* pTaskSet_ );

        // Thread 0 is main thread, otherwise use threadNum
        ENKITS_API void            AddPinnedTask( IPinnedTask* pTask_ );

        // This function will run any IPinnedTask* for current thread, but not run other
        // Main thread should call this or use a wait to ensure it's tasks are run.
        ENKITS_API void            RunPinnedTasks();

       // Runs the TaskSets in pipe until true == pTaskSet->GetIsComplete();
        // should only be called from thread which created the taskscheduler , or within a task
        // if called with 0 it will try to run tasks, and return if none available.
        // To run only a subset of tasks, set priorityOfLowestToRun_ to a high priority.
        // Default is lowest priority available.
        ENKITS_API void            WaitforTask( const ICompletable* pCompletable_, enki::TaskPriority priorityOfLowestToRun_ = TaskPriority(TASK_PRIORITY_NUM - 1) );

        // WaitforTaskSet, deprecated interface use WaitforTask
        inline void     WaitforTaskSet( const ICompletable* pCompletable_ ) { WaitforTask( pCompletable_ ); }

        // Waits for all task sets to complete - not guaranteed to work unless we know we
        // are in a situation where tasks aren't being continuosly added.
        ENKITS_API void            WaitforAll();

        // Waits for all task sets to complete and shutdown threads - not guaranteed to work unless we know we
        // are in a situation where tasks aren't being continuosly added.
        ENKITS_API void            WaitforAllAndShutdown();

        // Returns the number of threads created for running tasks + 1
        // to account for the main thread.
        ENKITS_API uint32_t        GetNumTaskThreads() const;

        // Returns the current task threadNum
        // Will return 0 for main thread and all other non-enkiTS threads, and < GetNumTaskThreads()
        ENKITS_API uint32_t        GetThreadNum() const;

        // Returns the ProfilerCallbacks structure so that it can be modified to
        // set the callbacks.
        ENKITS_API ProfilerCallbacks* GetProfilerCallbacks();

    private:
        static void     TaskingThreadFunction( const ThreadArgs& args_ );
        bool            HaveTasks( uint32_t threadNum_ );
        void            WaitForTasks( uint32_t threadNum_ );
        void            RunPinnedTasks( uint32_t threadNum_, uint32_t priority_ );
        bool            TryRunTask( uint32_t threadNum_, uint32_t& hintPipeToCheck_io_ );
        bool            TryRunTask( uint32_t threadNum_, uint32_t priority_, uint32_t& hintPipeToCheck_io_ );
        void            StartThreads();
        void            StopThreads( bool bWait_ );
        void            SplitAndAddTask( uint32_t threadNum_, SubTaskSet subTask_, uint32_t rangeToSplit_ );
        void            WakeThreads();

        TaskPipe*                                                m_pPipesPerThread[ TASK_PRIORITY_NUM ];
        PinnedTaskList*                                          m_pPinnedTaskListPerThread[ TASK_PRIORITY_NUM ];

        uint32_t                                                 m_NumThreads;
        ThreadArgs*                                              m_pThreadArgStore;
        std::thread**                                            m_pThreads;
        std::atomic<int32_t>                                     m_bRunning;
        std::atomic<int32_t>                                     m_NumThreadsRunning;
        std::atomic<int32_t>                                     m_NumThreadsWaiting;
        uint32_t                                                 m_NumPartitions;
        semaphoreid_t*                                           m_pNewTaskSemaphore;
        uint32_t                                                 m_NumInitialPartitions;
        bool                                                     m_bHaveThreads;
        ProfilerCallbacks                                        m_ProfilerCallbacks;

        TaskScheduler( const TaskScheduler& nocopy );
        TaskScheduler& operator=( const TaskScheduler& nocopy );
    };

    inline uint32_t GetNumHardwareThreads()
    {
        return std::thread::hardware_concurrency();
    }
}