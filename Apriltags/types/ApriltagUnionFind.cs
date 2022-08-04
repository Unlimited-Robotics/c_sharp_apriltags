using System.Collections;
using System.Collections.Generic;

namespace Apriltags
{
    public class UFrec
    {
        public uint Parent;
        public uint Size;

        public UFrec(uint parent, uint size)
        {
            Parent = parent;
            Size = size;
        }
    }

    public class UnionFind
    {
        public uint MaxID;
        public UFrec[] Data;

        public UnionFind(uint maxID)
        {
            MaxID = maxID;
            Data = new UFrec[MaxID+1];
            for (uint i = 0; i <= MaxID; i++)
            {
                Data[i] = new UFrec(i, 1);
            }
        }

        public UnionFind(Detector detector, Image threshImage, int w, int h) : this((uint)(w * h))
        {
            if(detector.HowManyThreadsToUse <= 1)
            {
                doUnionFindFirstLine(threshImage, w);
                for (int y = 1; y < h; y++) 
                {
                    doUnionfindLine2(threshImage, w, y);
                }
            }
            else
            {
                doUnionFindFirstLine(threshImage, w);

                int sz = h;
                int chunksize = 1 + sz / (Utils.Calculations.APRILTAG_TASKS_PER_THREAD_TARGET * detector.HowManyThreadsToUse);
                UnionFind.UnionFindTask[] tasks = new UnionFind.UnionFindTask[sz / chunksize + 1];

                int ntasks = 0;

                for (int i = 1; i < sz; i += chunksize) {
                    // each task will process [y0, y1). Note that this attaches
                    // each cell to the right and down, so row y1 *is* potentially modified.
                    //
                    // for parallelization, make sure that each task doesn't touch rows
                    // used by another thread.
                    tasks[ntasks] = new UnionFind.UnionFindTask();
                    tasks[ntasks].Y0 = i;
                    tasks[ntasks].Y1 = Utils.Calculations.IMin(sz, i + chunksize - 1);
                    tasks[ntasks].W = w;
                    tasks[ntasks].unionFind = this;
                    tasks[ntasks].ATImage = threshImage;

                    detector.WorkPool.Tasks.Add(tasks[ntasks]);
                    ntasks++;
                }

                detector.WorkPool.Run();

                // XXX stitch together the different chunks.
                for (int i = 1; i < ntasks; i++) 
                {
                    doUnionfindLine2(threshImage, w, tasks[i].Y0 - 1);
                }
            }
        }

        private void doUnionFindFirstLine(Image im, int w)
        {
            int y = 0;
            byte v;

            for (int x = 1; x < w - 1; x++) 
            {
                v = im.GetPixel(x, y);

                if (v == 127)
                {
                    continue;
                }

                int dx = -1, dy = 0;
                if (im.GetPixelCustom((y + dy)*im.Stride + x + dx) == v)
                {
                    unionfindConnect((uint)(y*w + x), (uint)((y + dy)*w + x + dx));
                } 
            }
        }

        private uint  unionfindConnect(uint aID, uint bID)
        {
            uint aRoot = GetRepresentative(aID);
            uint bRoot = GetRepresentative(bID);
            uint output = aRoot;

            if(aRoot != bRoot)
            {
                uint aSize = Data[aRoot].Size;
                uint bSize = Data[bRoot].Size;

                if(aSize > bSize)
                {
                    Data[bRoot].Parent = (uint)aRoot;
                    Data[aRoot].Size += bSize;
                    output = aRoot;
                }
                else
                {
                    Data[aRoot].Parent = (uint)bRoot;
                    Data[bRoot].Size += aSize;
                    output = bRoot;
                }
            }

            return output; 
        }

        public uint GetRepresentative(uint id)
        {
            uint root = id;

            while(Data[root].Parent != root)
            {
                root = Data[root].Parent;
            }

            while(Data[id].Parent != root)
            {
                uint tmp = Data[id].Parent;
                Data[id].Parent = root;
                id = tmp;
            }

            return root;
        }

        private void doUnionfindLine2(Image im, int w, int y)
        {
            byte v_m1_m1;
            byte v_0_m1 = im.GetPixelCustom((y - 1)*im.Stride);
            byte v_1_m1 = im.GetPixelCustom((y - 1)*im.Stride + 1);
            byte v_m1_0;
            byte v = im.GetPixelCustom(y*im.Stride);

            for (int x = 1; x < w - 1; x++) 
            {
                v_m1_m1 = v_0_m1;
                v_0_m1 = v_1_m1;
                v_1_m1 = im.GetPixelCustom((y - 1)*im.Stride + x + 1);
                v_m1_0 = v;
                v = im.GetPixel(x, y);

                if (v == 127)
                {
                    continue;
                }

                // (dx,dy) pairs for 8 connectivity:
                // (-1, -1)    (0, -1)    (1, -1)
                // (-1, 0)    (REFERENCE)
                int dx = -1, dy = 0;
                if (im.GetPixelCustom((y + dy)*im.Stride + x + dx) == v)
                {
                    unionfindConnect((uint)(y*w + x), (uint)((y + dy)*w + x + dx));
                } 

                if (x == 1 || !((v_m1_0 == v_m1_m1) && (v_m1_m1 == v_0_m1))) 
                {
                    int dx1 = 0, dy1 = -1;
                    if (im.GetPixelCustom((y + dy1)*im.Stride + x + dx1) == v)
                    {
                        unionfindConnect((uint)(y*w + x), (uint)((y + dy1)*w + x + dx1));
                    } 
                }

                if (v == 255) {
                    if (x == 1 || !(v_m1_0 == v_m1_m1 || v_0_m1 == v_m1_m1) ) 
                    {
                        int dx1 = -1, dy1 = -1;
                        if (im.GetPixelCustom((y + dy1)*im.Stride + x + dx1) == v)
                        {
                            unionfindConnect((uint)(y*w + x), (uint)((y + dy1)*w + x + dx1));
                        } 
                    }
                    if (!(v_0_m1 == v_1_m1)) 
                    {
                        int dx1 = 1, dy1 = -1;
                        if (im.GetPixelCustom((y + dy1)*im.Stride + x + dx1) == v)
                        {
                            unionfindConnect((uint)(y*w + x), (uint)((y + dy1)*w + x + dx1));
                        } 
                    }
                }
            }
        }

        public uint GetSetSize(uint id)
        {
            uint repid = GetRepresentative(id);
            return Data[repid].Size;
        }

        public class UnionFindTask : WorkerPool.WorkTask
        {
            public int Y0;
            public int Y1;
            public int W;
            public UnionFind unionFind;
            public Image ATImage;

            public override void DoTask()
            {
                for (int y = Y0; y < Y1; y++) 
                {
                    unionFind.doUnionfindLine2(ATImage, W, y);
                }
            }
        }
    }
}
