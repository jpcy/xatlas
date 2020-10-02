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

// Define ENKI_CUSTOM_ALLOC_FILE_AND_LINE (at project level) to get file and line report in custom allocators,
// this is default in Debug - to turn off define ENKI_CUSTOM_ALLOC_NO_FILE_AND_LINE
#ifndef ENKI_CUSTOM_ALLOC_FILE_AND_LINE
#if defined(_DEBUG ) && !defined(ENKI_CUSTOM_ALLOC_NO_FILE_AND_LINE)
#define ENKI_CUSTOM_ALLOC_FILE_AND_LINE
#endif
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
    class  Dependency;
    struct ThreadArgs;
    struct ThreadDataStore;
    struct SubTaskSet;
    struct semaphoreid_t;

    uint32_t GetNumHardwareThreads();

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
    // Can be used with dependencies to wait for their completion.
    // Derive from ITaskSet or IPinnedTask for running parallel tasks.
    class ICompletable
    {
    public:
        bool    GetIsComplete() const {
            return 0 == m_RunningCount.load( std::memory_order_acquire );
        }

        virtual ~ICompletable();

        // Dependency helpers, see Dependencies.cpp
        void SetDependency( Dependency& dependency_, const ICompletable* pDependencyTask_ );
        template<typename D, typename T, int SIZE> void SetDependenciesArr( D& dependencyArray_ , const T(&taskArray_)[SIZE] );
        template<typename D, typename T>           void SetDependenciesArr( D& dependencyArray_, std::initializer_list<T*> taskpList_ );
        template<typename D, typename T, int SIZE> void SetDependenciesArr( D(&dependencyArray_)[SIZE], const T(&taskArray_)[SIZE] );
        template<typename D, typename T, int SIZE> void SetDependenciesArr( D(&dependencyArray_)[SIZE], std::initializer_list<T*> taskpList_ );
        template<typename D, typename T, int SIZE> void SetDependenciesVec( D& dependencyVec_, const T(&taskArray_)[SIZE] );
        template<typename D, typename T>           void SetDependenciesVec( D& dependencyVec_, std::initializer_list<T*> taskpList_ );

        TaskPriority                   m_Priority            = TASK_PRIORITY_HIGH;
    protected:
        // Deriving from an ICompletable and overriding OnDependenciesComplete is advanced use.
        // If you do override OnDependenciesComplete() call:
        // ICompletable::OnDependenciesComplete( pTaskScheduler_, threadNum_ );
        // in your implementation.
        virtual void                   OnDependenciesComplete( TaskScheduler* pTaskScheduler_, uint32_t threadNum_ );
    private:
        friend class                   TaskScheduler;
        friend class                   Dependency;
        std::atomic<int32_t>           m_RunningCount               = {0};
        std::atomic<int32_t>           m_DependenciesCompletedCount = {0};
        int32_t                        m_DependenciesCount          = 0;
        mutable std::atomic<int32_t>   m_WaitingForTaskCount        = {0};
        mutable Dependency*            m_pDependents                = NULL;
    };

    // Subclass ITaskSet to create tasks.
    // TaskSets can be re-used, but check completion first.
    class ITaskSet : public ICompletable
    {
    public:
        ITaskSet() = default;
        ITaskSet( uint32_t setSize_ )
            : m_SetSize( setSize_ )
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
        virtual void ExecuteRange( TaskSetPartition range_, uint32_t threadnum_  ) = 0;

        // Set Size - usually the number of data items to be processed, see ExecuteRange. Defaults to 1
        uint32_t     m_SetSize  = 1;

        // Min Range - Minimum size of of TaskSetPartition range when splitting a task set into partitions.
        // Designed for reducing scheduling overhead by preventing set being
        // divided up too small. Ranges passed to ExecuteRange will *not* be a mulitple of this,
        // only attempts to deliver range sizes larger than this most of the time.
        // This should be set to a value which results in computation effort of at least 10k
        // clock cycles to minimize task scheduler overhead.
        // NOTE: The last partition will be smaller than m_MinRange if m_SetSize is not a multiple
        // of m_MinRange.
        // Also known as grain size in literature.
        uint32_t     m_MinRange  = 1;

    private:
        friend class TaskScheduler;
        void         OnDependenciesComplete( TaskScheduler* pTaskScheduler_, uint32_t threadNum_ ) override final;
        uint32_t     m_RangeToRun = 1;
    };

    // Subclass IPinnedTask to create tasks which can be run on a given thread only.
    class IPinnedTask : public ICompletable
    {
    public:
        IPinnedTask() = default;
        IPinnedTask( uint32_t threadNum_ ) : threadNum(threadNum_) {}  // default is to run a task on main thread

        // IPinnedTask needs to be non abstract for intrusive list functionality.
        // Should never be called as should be overridden.
        virtual void Execute() { assert(false); }

        uint32_t                  threadNum = 0; // thread to run this pinned task on
        std::atomic<IPinnedTask*> pNext = {NULL};
    private:
        void         OnDependenciesComplete( TaskScheduler* pTaskScheduler_, uint32_t threadNum_ ) override final;
    };

    // A utility task set for creating tasks based on std::func.
    typedef std::function<void (TaskSetPartition range, uint32_t threadnum  )> TaskSetFunction;
    class TaskSet : public ITaskSet
    {
    public:
        TaskSet() = default;
        TaskSet( TaskSetFunction func_ ) : m_Function( func_ ) {}
        TaskSet( uint32_t setSize_, TaskSetFunction func_ ) : ITaskSet( setSize_ ), m_Function( func_ ) {}


        virtual void ExecuteRange( TaskSetPartition range_, uint32_t threadnum_  )
        {
            m_Function( range_, threadnum_ );
        }

        TaskSetFunction m_Function;
    };

    class Dependency
    {
    public:
                        Dependency() = default; 
                        Dependency( const Dependency& ) = delete;
        ENKITS_API      Dependency( Dependency&& ) noexcept;	
        ENKITS_API      Dependency(    const ICompletable* pDependencyTask_, ICompletable* pTaskToRunOnCompletion_ );
        ENKITS_API      ~Dependency();

        ENKITS_API void SetDependency( const ICompletable* pDependencyTask_, ICompletable* pTaskToRunOnCompletion_ );
        ENKITS_API void ClearDependency();
              ICompletable* GetTaskToRunOnCompletion() { return pTaskToRunOnCompletion; }
        const ICompletable* GetDependencyTask()        { return pDependencyTask; }
    private:
        friend class TaskScheduler; friend class ICompletable;
        ICompletable*       pTaskToRunOnCompletion = NULL;
        const ICompletable* pDependencyTask        = NULL;
        Dependency*         pNext                  = NULL;
    };

    // TaskScheduler implements several callbacks intended for profilers
    typedef void (*ProfilerCallbackFunc)( uint32_t threadnum_ );
    struct ProfilerCallbacks
    {
        ProfilerCallbackFunc threadStart;
        ProfilerCallbackFunc threadStop;
        ProfilerCallbackFunc waitForNewTaskSuspendStart;      // thread suspended waiting for new tasks
        ProfilerCallbackFunc waitForNewTaskSuspendStop;       // thread unsuspended
        ProfilerCallbackFunc waitForTaskCompleteStart;        // thread waiting for task completion
        ProfilerCallbackFunc waitForTaskCompleteStop;         // thread stopped waiting
        ProfilerCallbackFunc waitForTaskCompleteSuspendStart; // thread suspended waiting task completion
        ProfilerCallbackFunc waitForTaskCompleteSuspendStop;  // thread unsuspended
    };

    // Custom allocator, set in TaskSchedulerConfig. Also see ENKI_CUSTOM_ALLOC_FILE_AND_LINE for file_ and line_
    typedef void* (*AllocFunc)( size_t align_, size_t size_, void* userData_, const char* file_, int line_ );
    typedef void  (*FreeFunc)(  void* ptr_,    size_t size_, void* userData_, const char* file_, int line_ );
    ENKITS_API void* DefaultAllocFunc(  size_t align_, size_t size_, void* userData_, const char* file_, int line_ );
    ENKITS_API void  DefaultFreeFunc(   void* ptr_,    size_t size_, void* userData_, const char* file_, int line_ );
    struct CustomAllocator
    {
        AllocFunc alloc    = DefaultAllocFunc;
        FreeFunc  free     = DefaultFreeFunc;
        void*     userData = nullptr;
    };

    // TaskSchedulerConfig - configuration struct for advanced Initialize
    struct TaskSchedulerConfig
    {
        // numTaskThreadsToCreate - Number of tasking threads the task scheduler will create. Must be > 0.
        // Defaults to GetNumHardwareThreads()-1 threads as thread which calls initialize is thread 0.
        uint32_t          numTaskThreadsToCreate = GetNumHardwareThreads()-1;

        // numExternalTaskThreads - Advanced use. Number of external threads which need to use TaskScheduler API.
        // See TaskScheduler::RegisterExternalTaskThread() for usage.
        // Defaults to 0. The thread used to initialize the TaskScheduler can also use the TaskScheduler API.
        // Thus there are (numTaskThreadsToCreate + numExternalTaskThreads + 1) able to use the API, with this
        // defaulting to the number of harware threads available to the system.
        uint32_t          numExternalTaskThreads = 0;

        ProfilerCallbacks profilerCallbacks = {};

        CustomAllocator   customAllocator;
    };

    class TaskScheduler
    {
    public:
        ENKITS_API TaskScheduler();
        ENKITS_API ~TaskScheduler();

        // Call an Initialize function before adding tasks.

        // Initialize() will create GetNumHardwareThreads()-1 tasking threads, which is
        // sufficient to fill the system when including the main thread.
        // Initialize can be called multiple times - it will wait for completion
        // before re-initializing.
        ENKITS_API void            Initialize();

        // Initialize( numThreadsTotal_ )
        // will create numThreadsTotal_-1 threads, as thread 0 is
        // the thread on which the initialize was called.
        // numThreadsTotal_ must be > 0
        ENKITS_API void            Initialize( uint32_t numThreadsTotal_ );

        // Initialize with advanced TaskSchedulerConfig settings. See TaskSchedulerConfig.
        ENKITS_API void            Initialize( TaskSchedulerConfig config_ );

        // Get config. Can be called before Initialize to get the defaults.
        ENKITS_API TaskSchedulerConfig GetConfig() const;

        // Adds the TaskSet to pipe and returns if the pipe is not full.
        // If the pipe is full, pTaskSet is run.
        // should only be called from main thread, or within a task
        ENKITS_API void            AddTaskSetToPipe( ITaskSet* pTaskSet_ );

        // Thread 0 is main thread, otherwise use threadNum
        // Pinned tasks can be added from any thread
        ENKITS_API void            AddPinnedTask( IPinnedTask* pTask_ );

        // This function will run any IPinnedTask* for current thread, but not run other
        // Main thread should call this or use a wait to ensure it's tasks are run.
        ENKITS_API void            RunPinnedTasks();

        // Runs the TaskSets in pipe until true == pTaskSet->GetIsComplete();
        // should only be called from thread which created the taskscheduler , or within a task
        // if called with 0 it will try to run tasks, and return if none available.
        // To run only a subset of tasks, set priorityOfLowestToRun_ to a high priority.
        // Default is lowest priority available.
        // Only wait for child tasks of the current task otherwise a deadlock could occur.
        ENKITS_API void            WaitforTask( const ICompletable* pCompletable_, enki::TaskPriority priorityOfLowestToRun_ = TaskPriority(TASK_PRIORITY_NUM - 1) );

        // Waits for all task sets to complete - not guaranteed to work unless we know we
        // are in a situation where tasks aren't being continuously added.
        ENKITS_API void            WaitforAll();

        // Waits for all task sets to complete and shutdown threads - not guaranteed to work unless we know we
        // are in a situation where tasks aren't being continuously added.
        // This function can be safely called even if TaskScheduler::Initialize() has not been called.
        ENKITS_API void            WaitforAllAndShutdown();

        // Returns the number of threads created for running tasks + number of external threads
        // plus 1 to account for the thread used to initialize the task scheduler.
        // Equivalent to config values: numTaskThreadsToCreate + numExternalTaskThreads + 1.
        // It is guaranteed that GetThreadNum() < GetNumTaskThreads()
        ENKITS_API uint32_t        GetNumTaskThreads() const;

        // Returns the current task threadNum
        // Will return 0 for thread which initialized the task scheduler,
        // and all other non-enkiTS threads which have not been registered ( see RegisterExternalTaskThread() ),
        // and < GetNumTaskThreads() for all threads.
        // It is guaranteed that GetThreadNum() < GetNumTaskThreads()
        ENKITS_API uint32_t        GetThreadNum() const;

         // Call on a thread to register the thread to use the TaskScheduling API.
        // This is implicitly done for the thread which initializes the TaskScheduler
        // Intended for developers who have threads who need to call the TaskScheduler API
        // Returns true if successfull, false if not.
        // Can only have numExternalTaskThreads registered at any one time, which must be set
        // at initialization time.
        ENKITS_API bool            RegisterExternalTaskThread();

        // Call on a thread on which RegisterExternalTaskThread has been called to deregister that thread.
        ENKITS_API void            DeRegisterExternalTaskThread();

        // Get the number of registered external task threads.
        ENKITS_API uint32_t        GetNumRegisteredExternalTaskThreads();


        // ------------- Start DEPRECATED Functions -------------
        // DEPRECATED - WaitforTaskSet, deprecated interface use WaitforTask
        inline void                WaitforTaskSet( const ICompletable* pCompletable_ ) { WaitforTask( pCompletable_ ); }

        // DEPRECATED - GetProfilerCallbacks.  Use TaskSchedulerConfig instead
        // Returns the ProfilerCallbacks structure so that it can be modified to
        // set the callbacks. Should be set prior to initialization.
        inline ProfilerCallbacks* GetProfilerCallbacks() { return &m_Config.profilerCallbacks; }
        // -------------  End DEPRECATED Functions  -------------

    private:
        friend class ICompletable; friend class ITaskSet; friend class IPinnedTask;
        static void TaskingThreadFunction( const ThreadArgs& args_ );
        bool        HaveTasks( uint32_t threadNum_ );
        void        WaitForNewTasks( uint32_t threadNum_ );
        void        WaitForTaskCompletion( const ICompletable* pCompletable_, uint32_t threadNum_ );
        void        RunPinnedTasks( uint32_t threadNum_, uint32_t priority_ );
        bool        TryRunTask( uint32_t threadNum_, uint32_t& hintPipeToCheck_io_ );
        bool        TryRunTask( uint32_t threadNum_, uint32_t priority_, uint32_t& hintPipeToCheck_io_ );
        void        StartThreads();
        void        StopThreads( bool bWait_ );
        void        SplitAndAddTask( uint32_t threadNum_, SubTaskSet subTask_, uint32_t rangeToSplit_ );
        void        WakeThreadsForNewTasks();
        void        WakeThreadsForTaskCompletion();
        bool        WakeSuspendedThreadsWithPinnedTasks();
        void        InitDependencies( ICompletable* pCompletable_  );
        ENKITS_API void TaskComplete( ICompletable* pTask_, bool bWakeThreads_, uint32_t threadNum_ );
        ENKITS_API void AddTaskSetToPipeInt( ITaskSet* pTaskSet_, uint32_t threadNum_ );
        ENKITS_API void AddPinnedTaskInt( IPinnedTask* pTask_ );

        template< typename T > T*   NewArray( size_t num_, const char* file_, int line_  );
        template< typename T > void DeleteArray( T* p_, size_t num_, const char* file_, int line_ );
        template<class T, class... Args> T* New( const char* file_, int line_,  Args&&... args_ );
        template< typename T > void Delete( T* p_, const char* file_, int line_ );
        semaphoreid_t* SemaphoreNew();
        void SemaphoreDelete( semaphoreid_t* pSemaphore_ );

        TaskPipe*              m_pPipesPerThread[ TASK_PRIORITY_NUM ];
        PinnedTaskList*        m_pPinnedTaskListPerThread[ TASK_PRIORITY_NUM ];

        uint32_t               m_NumThreads;
        ThreadDataStore*       m_pThreadDataStore;
        std::thread*           m_pThreads;
        std::atomic<int32_t>   m_bRunning;
        std::atomic<int32_t>   m_NumInternalTaskThreadsRunning;
        std::atomic<int32_t>   m_NumThreadsWaitingForNewTasks;
        std::atomic<int32_t>   m_NumThreadsWaitingForTaskCompletion;
        uint32_t               m_NumPartitions;
        semaphoreid_t*         m_pNewTaskSemaphore;
        semaphoreid_t*         m_pTaskCompleteSemaphore;
        uint32_t               m_NumInitialPartitions;
        bool                   m_bHaveThreads;
        TaskSchedulerConfig    m_Config;
        std::atomic<int32_t>   m_NumExternalTaskThreadsRegistered;

        TaskScheduler( const TaskScheduler& nocopy_ );
        TaskScheduler& operator=( const TaskScheduler& nocopy_ );

    protected:
        void SetCustomAllocator( CustomAllocator customAllocator_ ); // for C interface
    };

    inline uint32_t GetNumHardwareThreads()
    {
        return std::thread::hardware_concurrency();
    }

    inline void ICompletable::OnDependenciesComplete( TaskScheduler* pTaskScheduler_, uint32_t threadNum_ )
    {
        // mark complete
        m_RunningCount.store( 0, std::memory_order_seq_cst );
        pTaskScheduler_->TaskComplete( this, true, threadNum_ );
    }

    inline void ITaskSet::OnDependenciesComplete( TaskScheduler* pTaskScheduler_, uint32_t threadNum_ )
    {
        pTaskScheduler_->AddTaskSetToPipeInt( this, threadNum_ );
    }

    inline void IPinnedTask::OnDependenciesComplete( TaskScheduler* pTaskScheduler_, uint32_t threadNum_ )
    {
        (void)threadNum_;
        pTaskScheduler_->AddPinnedTaskInt( this );
    }

    inline ICompletable::~ICompletable()
    {
        Dependency* pDependency = m_pDependents;
        while( pDependency )
        {
            Dependency* pNext = pDependency->pNext;
            pDependency->pDependencyTask = NULL;
            pDependency->pNext = NULL;
            pDependency = pNext;
        }
    }

    inline void ICompletable::SetDependency( Dependency& dependency_, const ICompletable* pDependencyTask_ )
    {
        assert( pDependencyTask_ != this );
        dependency_.SetDependency( pDependencyTask_, this );
    }

    template<typename D, typename T, int SIZE>
    void ICompletable::SetDependenciesArr( D& dependencyArray_ , const T(&taskArray_)[SIZE] ) {
        static_assert( std::tuple_size<D>::value >= SIZE, "Size of dependency array too small" );
        for( int i = 0; i < SIZE; ++i )
        {
            dependencyArray_[i].SetDependency( &taskArray_[i], this );
        }
    }
    template<typename D, typename T>
    void ICompletable::SetDependenciesArr( D& dependencyArray_, std::initializer_list<T*> taskpList_ ) {
        assert( std::tuple_size<D>::value >= taskpList_.size() );
        int i = 0;
        for( auto pTask : taskpList_ )
        {
            dependencyArray_[i++].SetDependency( pTask, this );
        }
    }
    template<typename D, typename T, int SIZE>
    void ICompletable::SetDependenciesArr( D(&dependencyArray_)[SIZE], const T(&taskArray_)[SIZE] ) {
        for( int i = 0; i < SIZE; ++i )
        {
            dependencyArray_[i].SetDependency( &taskArray_[i], this );
        }
    }
    template<typename D, typename T, int SIZE>
    void ICompletable::SetDependenciesArr( D(&dependencyArray_)[SIZE], std::initializer_list<T*> taskpList_ ) {
        assert( SIZE >= taskpList_.size() );
        int i = 0;
        for( auto pTask : taskpList_ )
        {
            dependencyArray_[i++].SetDependency( pTask, this );
        }
    }
    template<typename D, typename T, int SIZE>
    void ICompletable::SetDependenciesVec( D& dependencyVec_, const T(&taskArray_)[SIZE] ) {
        dependencyVec_.resize( SIZE );
        for( int i = 0; i < SIZE; ++i )
        {
            dependencyVec_[i].SetDependency( &taskArray_[i], this );
        }
    }

    template<typename D, typename T>
    void ICompletable::SetDependenciesVec( D& dependencyVec_, std::initializer_list<T*> taskpList_ ) {
        dependencyVec_.resize( taskpList_.size() );
        int i = 0;
        for( auto pTask : taskpList_ )
        {
            dependencyVec_[i++].SetDependency( pTask, this );
        }
    }
}