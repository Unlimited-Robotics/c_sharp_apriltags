using System.Collections;
using System.Collections.Generic;
using System.Threading;
using System;
using UnityEngine;


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
        public Thread[] Threads;
        public List<int> Status;
        public List<WorkTask> Tasks;
        public int EndCount;
        public Mutex MutexLock;

        public WorkerPool(int howManyThreads)
        {
            NThreads = howManyThreads;
            Tasks = new List<WorkTask>();

            if(NThreads > 1)
            {
                Threads = new Thread[NThreads];
                MutexLock = new Mutex();
            }

        }

        public void Run()
        {
            if(NThreads > 1)
            {

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
    }
}
