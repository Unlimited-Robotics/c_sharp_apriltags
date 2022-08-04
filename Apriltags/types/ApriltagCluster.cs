using System.Collections;
using System.Collections.Generic;

namespace Apriltags
{
    

    public class Cluster
    {
        public class ClusterPoint
        {
            public ushort X, Y;
            public short GX, GY;
            public double Slope;

            public ClusterPoint(ushort x, ushort y, short gx, short gy)
            {
                X = x;
                Y = y;
                GX = gx;
                GY = gy;
                Slope = 0;
            }

            public void Copy(ClusterPoint vals)
            {
                X = vals.X;
                Y = vals.Y;
                GX = vals.GX;
                GY = vals.GY;
                Slope = vals.Slope;
            }
        }

        public List<ClusterPoint> Points = new List<ClusterPoint>();



        public class ClusterHash
        {
            public uint Hash;
            public ulong ID;
            public Cluster Data;

            public ClusterHash(uint hash, ulong id, Cluster data)
            {
                Hash = hash;
                ID = id;
                Data = data;
            }

            public ClusterHash(ClusterHash clusterHash)
            {
                Hash = clusterHash.Hash;
                ID = clusterHash.ID;
                Data = new Cluster();
                for (int i = 0; i < clusterHash.Data.Points.Count; i++)
                {
                    ClusterPoint current = clusterHash.Data.Points[i];
                    Data.Points.Add(new ClusterPoint(current.X, current.Y, current.GX, current.GY));
                }
            }
        }

        public class ClusterTask : WorkerPool.WorkTask
        {
            public int Y0;
            public int Y1;
            public int W;
            public int S;
            public int NClusterMap;
            public UnionFind UF;
            public Image Im;
            public List<Cluster.ClusterHash> Clusters;

            public override void DoTask()
            {
                ClusterMap.DoGradientClusters(Im, S, Y0, Y1, W, NClusterMap, UF, Clusters);
            }
        }

        public Cluster()
        {

        }

        public Cluster(List<ClusterPoint> points, int start, int size)
        {
            for (int i = start; i < start + size; i++)            
            {
                Points.Add(points[i]);
            }
        }

        public void PtSort()
        {
            int sz = Points.Count;

            if (sz <= 1)
            {
                return;
            }

            if (sz == 2) 
            {
                // MAYBE_SWAP(pts, 0, 1);
                int ap = 0;
                int bp = 1;
                maybeSwap(ap, bp);
                return;
            }

            // NB: Using less-branch-intensive sorting networks here on the
            // hunch that it's better for performance.
            if (sz == 3)
            { // 3 element bubble sort is optimal
                int ap, bp;
                ap = 0;
                bp = 1;
                maybeSwap(ap, bp);
                ap = 1;
                bp = 2;
                maybeSwap(ap, bp);
                ap = 0;
                bp = 1;
                maybeSwap(ap, bp);
                return;
            }

            if (sz == 4)
            { // 4 element optimal sorting network.
                int ap, bp;
                ap = 0;
                bp = 1;
                maybeSwap(ap, bp);
                ap = 2;
                bp = 3;
                maybeSwap(ap, bp);
                ap = 0;
                bp = 2;
                maybeSwap(ap, bp);
                ap = 1;
                bp = 3;
                maybeSwap(ap, bp);
                ap = 1;
                bp = 2;
                maybeSwap(ap, bp);
                return;
            }
            if (sz == 5)
            {
                // this 9-step swap is optimal for a sorting network, but two
                // steps slower than a generic sort.
                int ap, bp;
                ap = 0;
                bp = 1;
                maybeSwap(ap, bp);
                ap = 3;
                bp = 4;
                maybeSwap(ap, bp);
                ap = 1;
                bp = 2;
                maybeSwap(ap, bp);
                ap = 0;
                bp = 1;
                maybeSwap(ap, bp);
                ap = 0;
                bp = 3;
                maybeSwap(ap, bp);
                ap = 2;
                bp = 4;
                maybeSwap(ap, bp);
                ap = 1;
                bp = 2;
                maybeSwap(ap, bp);
                ap = 2;
                bp = 3;
                maybeSwap(ap, bp);
                ap = 1;
                bp = 2;
                maybeSwap(ap, bp);
                return;
            }

            // a merge sort with temp storage.

            // struct pt *tmp = malloc(sizeof(struct pt) *sz);

            // memcpy(tmp, pts, sizeof(struct pt) *sz);

            int asz = sz / 2;
            int bsz = sz - asz;

            Cluster ac = new Cluster(Points, 0, asz);
            Cluster bc = new Cluster(Points, asz, bsz);

            ac.PtSort();
            bc.PtSort();

            int apos = 0, bpos = 0, outpos = 0;
            while (apos + 8 < asz && bpos + 8 < bsz)
            {
                if (ac.Points[apos].Slope - bc.Points[bpos].Slope < 0)   
                {
                    Points[outpos++] = ac.Points[apos++];             
                }     
                else
                {
                    Points[outpos++] = bc.Points[bpos++];
                }  
                if (ac.Points[apos].Slope - bc.Points[bpos].Slope < 0)   
                {
                    Points[outpos++] = ac.Points[apos++];             
                }     
                else
                {
                    Points[outpos++] = bc.Points[bpos++];
                }  
                if (ac.Points[apos].Slope - bc.Points[bpos].Slope < 0)   
                {
                    Points[outpos++] = ac.Points[apos++];             
                }     
                else
                {
                    Points[outpos++] = bc.Points[bpos++];
                }  
                if (ac.Points[apos].Slope - bc.Points[bpos].Slope < 0)   
                {
                    Points[outpos++] = ac.Points[apos++];             
                }     
                else
                {
                    Points[outpos++] = bc.Points[bpos++];
                }  
                if (ac.Points[apos].Slope - bc.Points[bpos].Slope < 0)   
                {
                    Points[outpos++] = ac.Points[apos++];             
                }     
                else
                {
                    Points[outpos++] = bc.Points[bpos++];
                }   
                if (ac.Points[apos].Slope - bc.Points[bpos].Slope < 0)   
                {
                    Points[outpos++] = ac.Points[apos++];             
                }     
                else
                {
                    Points[outpos++] = bc.Points[bpos++];
                }   
                if (ac.Points[apos].Slope - bc.Points[bpos].Slope < 0)   
                {
                    Points[outpos++] = ac.Points[apos++];             
                }     
                else
                {
                    Points[outpos++] = bc.Points[bpos++];
                }  
                if (ac.Points[apos].Slope - bc.Points[bpos].Slope < 0)   
                {
                    Points[outpos++] = ac.Points[apos++];             
                }     
                else
                {
                    Points[outpos++] = bc.Points[bpos++];
                }   
            }

            while (apos < asz && bpos < bsz)
            {
                // Debug.Log("apos " + apos + ", asz " + asz + ", a size " + ac.Points.Count);
                // Debug.Log("bpos " + bpos + ", bsz " + bsz + ", b size " + bc.Points.Count);
                if (ac.Points[apos].Slope - bc.Points[bpos].Slope < 0)   
                {
                    Points[outpos++] = ac.Points[apos++];             
                }     
                else
                {
                    Points[outpos++] = bc.Points[bpos++];
                }  
            }

            while(apos < asz)
            {
                Points[outpos++] = ac.Points[apos++];
            }

            while(bpos < bsz)
            {
                Points[outpos++] = bc.Points[bpos++];
            }
        }

        private void maybeSwap(int apos, int bpos)
        {
            if (Points[apos].Slope - Points[bpos].Slope > 0) 
            {                        
                ClusterPoint tmp = Points[apos]; 
                Points[apos] = Points[bpos]; 
                Points[bpos] = tmp;    
            }
        }
    }
}