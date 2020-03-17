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

#include "TaskScheduler_c.h"
#include "TaskScheduler.h"

#include <assert.h>

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

void* enkiDefaultAllocFunc(  size_t align_, size_t size_, void* userData_, const char* file_, int line_ )
{
    return enki::DefaultAllocFunc( align_, size_, userData_, file_, line_ );
}

void  enkiDefaultFreeFunc(   void* ptr_,    size_t size_, void* userData_, const char* file_, int line_ )
{
    return enki::DefaultFreeFunc( ptr_, size_, userData_, file_, line_ );
}

struct enkiTaskScheduler : TaskScheduler
{
    void enkiSetCustomAllocator( CustomAllocator customAllocator_ )
    {
        SetCustomAllocator( customAllocator_ );
    }
};

struct enkiCompletable : ICompletable {}; // empty struct which we will use for completables

struct enkiTaskSet : ITaskSet
{
    enkiTaskSet( enkiTaskExecuteRange taskFun_ ) : taskFun(taskFun_), pArgs(NULL) {}

    void ExecuteRange( TaskSetPartition range_, uint32_t threadnum_  ) override
    {
        taskFun( range_.start, range_.end, threadnum_, pArgs );
    }

    enkiTaskExecuteRange taskFun;
    void* pArgs = NULL;
};

struct enkiPinnedTask : IPinnedTask
{
    enkiPinnedTask( enkiPinnedTaskExecute taskFun_, uint32_t threadNum_ )
        : IPinnedTask( threadNum_ ), taskFun(taskFun_), pArgs(NULL) {}

    void Execute() override
    {
        taskFun( pArgs );
    }

    enkiPinnedTaskExecute taskFun;
    void* pArgs = NULL;
};

struct enkiCompletionAction : ICompletable
{
    void OnDependenciesComplete( TaskScheduler* pTaskScheduler_, uint32_t threadNum_ )
    {
        ICompletable::OnDependenciesComplete( pTaskScheduler_, threadNum_ );

        completionFunction( pArgs, threadNum_ );
    }

    enkiCompletionFunction completionFunction;
    Dependency             dependency;
    void* pArgs = NULL;
};

struct enkiDependency : Dependency {}; // empty struct which we will use for dependencies

enkiTaskScheduler* enkiNewTaskScheduler()
{
    CustomAllocator defaultAllocator;
    enkiTaskScheduler* pETS = (enkiTaskScheduler*)defaultAllocator.alloc(
        alignof(enkiTaskScheduler), sizeof(enkiTaskScheduler), defaultAllocator.userData, ENKI_FILE_AND_LINE );
    new(pETS) enkiTaskScheduler;
    return pETS;
}

enkiTaskScheduler* enkiNewTaskSchedulerWithCustomAllocator( struct enkiCustomAllocator customAllocator_ )
{
    enkiTaskScheduler* pETS = (enkiTaskScheduler*)customAllocator_.alloc(
        alignof(enkiTaskScheduler), sizeof(enkiTaskScheduler), customAllocator_.userData, ENKI_FILE_AND_LINE );

    CustomAllocator customAllocatorCpp;
    customAllocatorCpp.alloc      = customAllocator_.alloc;
    customAllocatorCpp.free       = customAllocator_.free;
    customAllocatorCpp.userData = customAllocator_.userData;

    new(pETS) enkiTaskScheduler;
    pETS->enkiSetCustomAllocator( customAllocatorCpp );
    return pETS;
}

struct enkiTaskSchedulerConfig enkiGetTaskSchedulerConfig( enkiTaskScheduler* pETS_ )
{
    TaskSchedulerConfig config = pETS_->GetConfig();
    enkiTaskSchedulerConfig configC;
    configC.numExternalTaskThreads                             = config.numExternalTaskThreads;
    configC.numTaskThreadsToCreate                             = config.numTaskThreadsToCreate;
    configC.profilerCallbacks.threadStart                      = config.profilerCallbacks.threadStart;                      
    configC.profilerCallbacks.threadStop                       = config.profilerCallbacks.threadStop;                       
    configC.profilerCallbacks.waitForNewTaskSuspendStart       = config.profilerCallbacks.waitForNewTaskSuspendStart;      
    configC.profilerCallbacks.waitForNewTaskSuspendStop        = config.profilerCallbacks.waitForNewTaskSuspendStop;        
    configC.profilerCallbacks.waitForTaskCompleteStart         = config.profilerCallbacks.waitForTaskCompleteStart;         
    configC.profilerCallbacks.waitForTaskCompleteStop          = config.profilerCallbacks.waitForTaskCompleteStop;          
    configC.profilerCallbacks.waitForTaskCompleteSuspendStart  = config.profilerCallbacks.waitForTaskCompleteSuspendStart;  
    configC.profilerCallbacks.waitForTaskCompleteSuspendStop   = config.profilerCallbacks.waitForTaskCompleteSuspendStop;
    configC.customAllocator.alloc                              = config.customAllocator.alloc;
    configC.customAllocator.free                               = config.customAllocator.free;
    configC.customAllocator.userData                         = config.customAllocator.userData;
    return configC;
}

void enkiInitTaskScheduler(  enkiTaskScheduler* pETS_ )
{
    pETS_->Initialize();
}

void enkiInitTaskSchedulerNumThreads(  enkiTaskScheduler* pETS_, uint32_t numThreads_ )
{
    pETS_->Initialize( numThreads_ );
}

void enkiInitTaskSchedulerWithConfig( enkiTaskScheduler* pETS_, struct enkiTaskSchedulerConfig config_ )
{
    TaskSchedulerConfig config;
    config.numExternalTaskThreads                             = config_.numExternalTaskThreads;
    config.numTaskThreadsToCreate                             = config_.numTaskThreadsToCreate;
    config.profilerCallbacks.threadStart                      = config_.profilerCallbacks.threadStart;                      
    config.profilerCallbacks.threadStop                       = config_.profilerCallbacks.threadStop;                       
    config.profilerCallbacks.waitForNewTaskSuspendStart       = config_.profilerCallbacks.waitForNewTaskSuspendStart;      
    config.profilerCallbacks.waitForNewTaskSuspendStop        = config_.profilerCallbacks.waitForNewTaskSuspendStop;        
    config.profilerCallbacks.waitForTaskCompleteStart         = config_.profilerCallbacks.waitForTaskCompleteStart;         
    config.profilerCallbacks.waitForTaskCompleteStop          = config_.profilerCallbacks.waitForTaskCompleteStop;          
    config.profilerCallbacks.waitForTaskCompleteSuspendStart  = config_.profilerCallbacks.waitForTaskCompleteSuspendStart;  
    config.profilerCallbacks.waitForTaskCompleteSuspendStop   = config_.profilerCallbacks.waitForTaskCompleteSuspendStop;   
    config.customAllocator.alloc                              = config_.customAllocator.alloc;
    config.customAllocator.free                               = config_.customAllocator.free;
    config.customAllocator.userData                         = config_.customAllocator.userData;
   pETS_->Initialize( config );
}

void enkiDeleteTaskScheduler( enkiTaskScheduler* pETS_ )
{
    CustomAllocator customAllocator = pETS_->GetConfig().customAllocator;
    pETS_->~enkiTaskScheduler();
    customAllocator.free( pETS_, sizeof(enkiTaskScheduler), customAllocator.userData, ENKI_FILE_AND_LINE );
}

enkiTaskSet* enkiCreateTaskSet( enkiTaskScheduler* pETS_, enkiTaskExecuteRange taskFunc_  )
{
    const CustomAllocator& customAllocator = pETS_->GetConfig().customAllocator;
    enkiTaskSet* pTask = (enkiTaskSet*)customAllocator.alloc(
        alignof(enkiTaskSet), sizeof(enkiTaskSet), customAllocator.userData, ENKI_FILE_AND_LINE );
    new(pTask) enkiTaskSet( taskFunc_ );

    return pTask;
}

void enkiDeleteTaskSet( enkiTaskScheduler* pETS_, enkiTaskSet* pTaskSet_ )
{
    const CustomAllocator& customAllocator = pETS_->GetConfig().customAllocator;

    pTaskSet_->~enkiTaskSet();
    customAllocator.free( pTaskSet_, sizeof(enkiTaskSet), customAllocator.userData, ENKI_FILE_AND_LINE );
}

enkiParamsTaskSet enkiGetParamsTaskSet( enkiTaskSet* pTaskSet_ )
{
    enkiParamsTaskSet paramsTaskSet;
    paramsTaskSet.pArgs    = pTaskSet_->pArgs;
    paramsTaskSet.setSize  = pTaskSet_->m_SetSize;
    paramsTaskSet.minRange = pTaskSet_->m_MinRange;
    paramsTaskSet.priority = pTaskSet_->m_Priority;
    return paramsTaskSet;
}

void enkiSetParamsTaskSet( enkiTaskSet* pTaskSet_, enkiParamsTaskSet params_ )
{
    assert( params_.priority < ENKITS_TASK_PRIORITIES_NUM );
    pTaskSet_->pArgs      = params_.pArgs;
    pTaskSet_->m_SetSize  = params_.setSize;
    pTaskSet_->m_MinRange = params_.minRange;
    pTaskSet_->m_Priority = TaskPriority( params_.priority );
}

void enkiSetPriorityTaskSet( enkiTaskSet* pTaskSet_, int priority_ )
{
    assert( priority_ < ENKITS_TASK_PRIORITIES_NUM );
    pTaskSet_->m_Priority = TaskPriority( priority_ );
}

void enkiSetArgsTaskSet( enkiTaskSet* pTaskSet_, void* pArgs_ )
{
    pTaskSet_->pArgs = pArgs_;
}

void enkiSetSetSizeTaskSet( enkiTaskSet* pTaskSet_, uint32_t setSize_ )
{
    pTaskSet_->m_SetSize = setSize_;
}

void enkiSetMinRangeTaskSet( enkiTaskSet* pTaskSet_, uint32_t minRange_ )
{
    pTaskSet_->m_MinRange = minRange_;
}

void enkiAddTaskSet( enkiTaskScheduler* pETS_, enkiTaskSet* pTaskSet_ )
{
    assert( pTaskSet_ );
    assert( pTaskSet_->taskFun );

    pETS_->AddTaskSetToPipe( pTaskSet_ );
}

void enkiAddTaskSetArgs( enkiTaskScheduler* pETS_, enkiTaskSet* pTaskSet_, void* pArgs_, uint32_t setSize_ )
{
    assert( pTaskSet_ );
    assert( pTaskSet_->taskFun );

    pTaskSet_->m_SetSize = setSize_;
    pTaskSet_->pArgs = pArgs_;
    pETS_->AddTaskSetToPipe( pTaskSet_ );
}

void enkiAddTaskSetMinRange(enkiTaskScheduler* pETS_, enkiTaskSet* pTaskSet_, void* pArgs_, uint32_t setSize_, uint32_t minRange_)
{
    assert( pTaskSet_ );
    assert( pTaskSet_->taskFun );

    pTaskSet_->m_SetSize = setSize_;
    pTaskSet_->m_MinRange = minRange_;
    pTaskSet_->pArgs = pArgs_;
    pETS_->AddTaskSetToPipe( pTaskSet_ );
}

int enkiIsTaskSetComplete( enkiTaskScheduler* pETS_, enkiTaskSet* pTaskSet_ )
{
    (void)pETS_;
    assert( pTaskSet_ );
    return ( pTaskSet_->GetIsComplete() ) ? 1 : 0;
}

enkiPinnedTask* enkiCreatePinnedTask(enkiTaskScheduler* pETS_, enkiPinnedTaskExecute taskFunc_, uint32_t threadNum_)
{
    const CustomAllocator& customAllocator = pETS_->GetConfig().customAllocator;
    enkiPinnedTask* pTask = (enkiPinnedTask*)customAllocator.alloc(
        alignof(enkiPinnedTask), sizeof(enkiPinnedTask), customAllocator.userData, ENKI_FILE_AND_LINE );
    new(pTask) enkiPinnedTask( taskFunc_, threadNum_ );
    return pTask;
}

void enkiDeletePinnedTask( enkiTaskScheduler* pETS_, enkiPinnedTask* pPinnedTask_ )
{
    const CustomAllocator& customAllocator = pETS_->GetConfig().customAllocator;

    pPinnedTask_->~enkiPinnedTask();
    customAllocator.free( pPinnedTask_, sizeof(enkiPinnedTask), customAllocator.userData, ENKI_FILE_AND_LINE );
}

enkiParamsPinnedTask enkiGetParamsPinnedTask( enkiPinnedTask* pTask_ )
{
    enkiParamsPinnedTask paramsPinnedTask;
    paramsPinnedTask.pArgs    = pTask_->pArgs;
    paramsPinnedTask.priority = pTask_->m_Priority;
    return paramsPinnedTask;
}

void enkiSetParamsPinnedTask( enkiPinnedTask* pTask_, enkiParamsPinnedTask params_ )
{
    assert( params_.priority < ENKITS_TASK_PRIORITIES_NUM );
    pTask_->pArgs      = params_.pArgs;
    pTask_->m_Priority = TaskPriority( params_.priority );
}

void enkiSetPriorityPinnedTask( enkiPinnedTask* pTask_, int priority_ )
{
    assert( priority_ < ENKITS_TASK_PRIORITIES_NUM );
    pTask_->m_Priority = TaskPriority( priority_ );
}

void enkiSetArgsPinnedTask( enkiPinnedTask* pTask_, void* pArgs_ )
{
    pTask_->pArgs = pArgs_;
}

void enkiAddPinnedTask(enkiTaskScheduler* pETS_, enkiPinnedTask* pTask_)
{
    assert( pTask_ );
    pETS_->AddPinnedTask( pTask_ );
}

void enkiAddPinnedTaskArgs(enkiTaskScheduler* pETS_, enkiPinnedTask* pTask_, void* pArgs_)
{
    assert( pTask_ );
    pTask_->pArgs = pArgs_;
    pETS_->AddPinnedTask( pTask_ );
}

void enkiRunPinnedTasks(enkiTaskScheduler* pETS_)
{
    pETS_->RunPinnedTasks();
}

int enkiIsPinnedTaskComplete(enkiTaskScheduler* pETS_, enkiPinnedTask* pTask_)
{
    assert( pTask_ );
    return ( pTask_->GetIsComplete() ) ? 1 : 0;
}

void enkiWaitForTaskSet( enkiTaskScheduler* pETS_, enkiTaskSet* pTaskSet_ )
{
    pETS_->WaitforTask( pTaskSet_ );
}

void enkiWaitForTaskSetPriority( enkiTaskScheduler * pETS_, enkiTaskSet * pTaskSet_, int maxPriority_ )
{
    pETS_->WaitforTask( pTaskSet_, TaskPriority( maxPriority_ ) );
}

void enkiWaitForPinnedTask( enkiTaskScheduler* pETS_, enkiPinnedTask* pTask_ )
{
    pETS_->WaitforTask( pTask_ );
}

void enkiWaitForPinnedTaskPriority( enkiTaskScheduler * pETS_, enkiPinnedTask * pTask_, int maxPriority_ )
{
    pETS_->WaitforTask( pTask_, TaskPriority( maxPriority_ ) );
}

void enkiWaitForAll( enkiTaskScheduler* pETS_ )
{
    pETS_->WaitforAll();
}

uint32_t enkiGetNumTaskThreads( enkiTaskScheduler* pETS_ )
{
    return pETS_->GetNumTaskThreads();
}

uint32_t enkiGetThreadNum( enkiTaskScheduler* pETS_ )
{
    return pETS_->GetThreadNum();
}

int enkiRegisterExternalTaskThread( enkiTaskScheduler* pETS_)
{
    return (int)pETS_->RegisterExternalTaskThread();
}

void enkiDeRegisterExternalTaskThread( enkiTaskScheduler* pETS_)
{
    return pETS_->DeRegisterExternalTaskThread();
}

uint32_t enkiGetNumRegisteredExternalTaskThreads( enkiTaskScheduler* pETS_)
{
    return pETS_->GetNumRegisteredExternalTaskThreads();
}

enkiCompletable* enkiGetCompletableFromTaskSet( enkiTaskSet* pTaskSet_ )
{
    return reinterpret_cast<enkiCompletable*>( static_cast<ICompletable*>( pTaskSet_ ) );
}

enkiCompletable* enkiGetCompletableFromPinnedTask( enkiPinnedTask* pPinnedTask_ )
{
    return reinterpret_cast<enkiCompletable*>( static_cast<ICompletable*>( pPinnedTask_ ) );
}

enkiCompletable* enkiGetCompletableFromCompletionAction( enkiCompletionAction* pCompletionAction_ )
{
    return reinterpret_cast<enkiCompletable*>( static_cast<ICompletable*>( pCompletionAction_ ) );
}

enkiCompletable* enkiCreateCompletable( enkiTaskScheduler* pETS_ )
{
    const CustomAllocator& customAllocator = pETS_->GetConfig().customAllocator;
    enkiCompletable* pTask = (enkiCompletable*)customAllocator.alloc(
        alignof(enkiCompletable), sizeof(enkiCompletable), customAllocator.userData, ENKI_FILE_AND_LINE );
    new(pTask) enkiCompletable();
    return pTask;
}

void enkiDeleteCompletable( enkiTaskScheduler* pETS_, enkiCompletable* pCompletable_ )
{
    const CustomAllocator& customAllocator = pETS_->GetConfig().customAllocator;

    pCompletable_->~enkiCompletable();
    customAllocator.free( pCompletable_, sizeof(enkiCompletable), customAllocator.userData, ENKI_FILE_AND_LINE );
}

void enkiWaitForCompletable( enkiTaskScheduler* pETS_, enkiCompletable* pTask_ )
{
    return pETS_->WaitforTask( reinterpret_cast<ICompletable*>( pTask_ ) );
}

void enkiWaitForCompletablePriority( enkiTaskScheduler* pETS_, enkiCompletable* pTask_, int maxPriority_ )
{
    return pETS_->WaitforTask( reinterpret_cast<ICompletable*>( pTask_ ), TaskPriority( maxPriority_ ) );
}

enkiDependency* enkiCreateDependency( enkiTaskScheduler* pETS_ )
{
    const CustomAllocator& customAllocator = pETS_->GetConfig().customAllocator;
    enkiDependency* pDep = (enkiDependency*)customAllocator.alloc(
        alignof(enkiDependency), sizeof(enkiDependency), customAllocator.userData, ENKI_FILE_AND_LINE );
    new(pDep) enkiDependency();
    return pDep;
}

void enkiDeleteDependency( enkiTaskScheduler* pETS_, enkiDependency* pDependency_ )
{
    const CustomAllocator& customAllocator = pETS_->GetConfig().customAllocator;

    pDependency_->~enkiDependency();
    customAllocator.free( pDependency_, sizeof(enkiDependency), customAllocator.userData, ENKI_FILE_AND_LINE );
}

void enkiSetDependency( enkiDependency* pDependency_, enkiCompletable* pDependencyTask_, enkiCompletable* pTaskToRunOnCompletion_ )
{
    pTaskToRunOnCompletion_->SetDependency( *pDependency_, pDependencyTask_ );
}

enkiCompletionAction* enkiCreateCompletionAction( enkiTaskScheduler* pETS_, enkiCompletionFunction completionFunc_ )
{
    const CustomAllocator& customAllocator = pETS_->GetConfig().customAllocator;
    enkiCompletionAction* pCA = (enkiCompletionAction*)customAllocator.alloc(
        alignof(enkiCompletionAction), sizeof(enkiCompletionAction), customAllocator.userData, ENKI_FILE_AND_LINE );
    new(pCA) enkiCompletionAction();
    pCA->completionFunction = completionFunc_;
    return pCA;
}

void enkiDeleteCompletionAction( enkiTaskScheduler* pETS_, enkiCompletionAction* pCompletionAction_ )
{
    const CustomAllocator& customAllocator = pETS_->GetConfig().customAllocator;

    pCompletionAction_->~enkiCompletionAction();
    customAllocator.free( pCompletionAction_, sizeof(enkiCompletionAction), customAllocator.userData, ENKI_FILE_AND_LINE );
}

enkiParamsCompletionAction enkiGetParamsCompletionAction( enkiCompletionAction* pCompletionAction_ )
{
    enkiParamsCompletionAction params;
    params.pArgs = pCompletionAction_->pArgs;
    params.pDependency = reinterpret_cast<const enkiCompletable*>(
        pCompletionAction_->dependency.GetDependencyTask() );
    return params;
}

void enkiSetParamsCompletionAction( enkiCompletionAction* pCompletionAction_, enkiParamsCompletionAction params_ )
{
    pCompletionAction_->pArgs = params_.pArgs;
    pCompletionAction_->SetDependency( pCompletionAction_->dependency, params_.pDependency );
}
