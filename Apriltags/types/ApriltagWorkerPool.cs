using System.Collections;
using System.Collections.Generic;
using System.Threading;
using System.Threading.Tasks;
using System;

namespace Apriltags
{
    public class WorkerPool
    {
        public abstract class WorkTask
        {
            public abstract void DoTask();
        }

        public int NThreads;
        public int TaskPos;
        public Task[] Threads;
        public List<int> Status;
        public List<WorkTask> Tasks;
        public int EndCount;
        public Mutex GetTaskLock;
        public Mutex FinishTaskLock;
        public SemaphoreSlim FinishTasksSemaphore;

        public WorkerPool(int howManyThreads)
        {
            NThreads = howManyThreads;
            Tasks = new List<WorkTask>();

            if(NThreads > 1)
            {
                Threads = new Task[NThreads];
                GetTaskLock = new Mutex();
                FinishTaskLock = new Mutex();
                FinishTasksSemaphore = new SemaphoreSlim(0,1);
            }
        }

        public void Run()
        {
            EndCount = 0;
            if(NThreads > 1)
            {
                for (int i = 0; i < NThreads; i++)
                {
                    Threads[i] = new Task(() => completeTasks());
                    Threads[i].Start();
                }

                FinishTasksSemaphore.Wait();
                Tasks.Clear();
            }
            else
            {
                runSingle();
            }
        }

        private void runSingle()
        {
            for (int i = 0; i < Tasks.Count; i++) 
            {
                Tasks[i].DoTask();
            }

            Tasks.Clear();
        }

        private void completeTasks()
        {
            bool moreWork = true;
            while(moreWork == true)
            {
                WorkTask currentTask = null;
                GetTaskLock.WaitOne();
                if(Tasks.Count > 0)
                {
                    currentTask = Tasks[0];
                    Tasks.RemoveAt(0);
                }
                GetTaskLock.ReleaseMutex();

                if(currentTask != null)
                {
                    currentTask.DoTask();
                }
                else
                {
                    moreWork = false;
                }
            }

            FinishTaskLock.WaitOne();
            EndCount++;
            if(EndCount == NThreads)
            {
                FinishTasksSemaphore.Release(1);
            }
            FinishTaskLock.ReleaseMutex();
        }
    }
}
