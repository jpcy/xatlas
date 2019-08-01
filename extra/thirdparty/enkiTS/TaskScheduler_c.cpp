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

struct enkiTaskScheduler : TaskScheduler
{
};

struct enkiTaskSet : ITaskSet
{
    enkiTaskSet( enkiTaskExecuteRange taskFun_ ) : taskFun(taskFun_), pArgs(NULL) {}

    virtual void ExecuteRange( TaskSetPartition range, uint32_t threadnum  )
    {
        taskFun( range.start, range.end, threadnum, pArgs );
    }

    enkiTaskExecuteRange taskFun;
    void* pArgs;
};

struct enkiPinnedTask : IPinnedTask
{
    enkiPinnedTask( enkiPinnedTaskExecute taskFun_, uint32_t threadNum_ )
        : IPinnedTask( threadNum_ ), taskFun(taskFun_), pArgs(NULL) {}

    virtual void Execute()
    {
        taskFun( pArgs );
    }

    enkiPinnedTaskExecute taskFun;
    void* pArgs;
};

enkiTaskScheduler* enkiNewTaskScheduler()
{
    enkiTaskScheduler* pETS = new enkiTaskScheduler();
    return pETS;
}

void enkiInitTaskScheduler(  enkiTaskScheduler* pETS_ )
{
    pETS_->Initialize();
}

void enkiInitTaskSchedulerNumThreads(  enkiTaskScheduler* pETS_, uint32_t numThreads_ )
{
    pETS_->Initialize( numThreads_ );
}

void enkiDeleteTaskScheduler( enkiTaskScheduler* pETS_ )
{
    delete pETS_;
}

enkiTaskSet* enkiCreateTaskSet( enkiTaskScheduler* pETS_, enkiTaskExecuteRange taskFunc_  )
{
    return new enkiTaskSet( taskFunc_ );
}

void enkiDeleteTaskSet( enkiTaskSet* pTaskSet_ )
{
    delete pTaskSet_;
}

void enkiSetPriorityTaskSet( enkiTaskSet* pTaskSet_, int priority_ )
{
    assert( priority_ < ENKITS_TASK_PRIORITIES_NUM );
    pTaskSet_->m_Priority = TaskPriority( priority_ );
}

void enkiAddTaskSetToPipe( enkiTaskScheduler* pETS_, enkiTaskSet* pTaskSet_, void* pArgs_, uint32_t setSize_ )
{
    assert( pTaskSet_ );
    assert( pTaskSet_->taskFun );

    pTaskSet_->m_SetSize = setSize_;
    pTaskSet_->pArgs = pArgs_;
    pETS_->AddTaskSetToPipe( pTaskSet_ );
}

void enkiAddTaskSetToPipeMinRange(enkiTaskScheduler* pETS_, enkiTaskSet* pTaskSet_, void* pArgs_, uint32_t setSize_, uint32_t minRange_)
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
    assert( pTaskSet_ );
    return ( pTaskSet_->GetIsComplete() ) ? 1 : 0;
}

enkiPinnedTask* enkiCreatePinnedTask(enkiTaskScheduler* pETS_, enkiPinnedTaskExecute taskFunc_, uint32_t threadNum_)
{
    return new enkiPinnedTask( taskFunc_, threadNum_ );
}

void enkiDeletePinnedTask(enkiPinnedTask* pTaskSet_)
{
    delete pTaskSet_;
}

void enkiSetPriorityPinnedTask( enkiPinnedTask* pTask_, int priority_ )
{
    assert( priority_ < ENKITS_TASK_PRIORITIES_NUM );
    pTask_->m_Priority = TaskPriority( priority_ );
}

void enkiAddPinnedTask(enkiTaskScheduler* pETS_, enkiPinnedTask* pTask_, void* pArgs_)
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

enkiProfilerCallbacks*    enkiGetProfilerCallbacks( enkiTaskScheduler* pETS_ )
{
    assert( sizeof(enkiProfilerCallbacks) == sizeof(enki::ProfilerCallbacks) );
    return (enkiProfilerCallbacks*)pETS_->GetProfilerCallbacks();
}

