using System.Collections;
using System.Collections.Generic;

namespace Apriltags
{
    public class ClusterMap
    {
        public class ClusterMapEntry
        {
            public ulong ID;
            public Cluster Data;
            public ClusterMapEntry Next;
        }
        public ClusterMapEntry[] Entries;

        public ClusterMap()
        {

        }

        public ClusterMap(int size, bool init)
        {
            Entries = new ClusterMapEntry[size];

            if(init == true)
            {
                for (int i = 0; i < size; i++)
                {
                    Entries[i] = new ClusterMapEntry();
                }
            }
        }

        public static List<Cluster> GradientClusters(Detector detector, Image threshImage,
            int w, int h, int ts, UnionFind uf)
        {
            List<Cluster> output = new List<Cluster>();

            int nclustermap = (int)(0.2*w*h);

            int sz = h - 1;
            int chunksize = 1 + sz / detector.HowManyThreadsToUse;
            Cluster.ClusterTask[] tasks = new Cluster.ClusterTask[sz / chunksize + 1];

            int ntasks = 0;

            for (int i = 1; i < sz; i += chunksize) 
            {
                // each task will process [y0, y1). Note that this processes
                // each cell to the right and down.
                // Debug.Log("ntasks " + ntasks);
                tasks[ntasks] = new Cluster.ClusterTask();
                tasks[ntasks].Y0 = i;
                // Debug.Log("y0 " + i);
                tasks[ntasks].Y1 = Utils.Calculations.IMin(sz, i + chunksize);
                // Debug.Log("y1 " + tasks[ntasks].Y1);
                tasks[ntasks].W = w;
                // Debug.Log("w " + tasks[ntasks].W);
                tasks[ntasks].S = ts;
                // Debug.Log("s " + tasks[ntasks].S);
                tasks[ntasks].UF = uf;
                tasks[ntasks].Im = threshImage;
                tasks[ntasks].NClusterMap = nclustermap/(sz / chunksize + 1);
                tasks[ntasks].Clusters = new List<Cluster.ClusterHash>();

                detector.WorkPool.Tasks.Add(tasks[ntasks]);
                ntasks++;
            }

            detector.WorkPool.Run();

            // Utils.PrintClustersData(tasks[0].Clusters);
            List<Cluster.ClusterHash>[] clusters_list = new List<Cluster.ClusterHash>[ntasks];
            for (int i = 0; i < ntasks; i++) 
            {
                clusters_list[i] = tasks[i].Clusters;
            }

            int length = ntasks;
            while (length > 1) 
            {
                int write = 0;
                for (int i = 0; i < length - 1; i += 2) 
                {
                    clusters_list[write] = mergeClusters(clusters_list[i], clusters_list[i + 1]);
                    write++;
                }

                if (length % 2 == 1) 
                {
                    clusters_list[write] = clusters_list[length - 1];
                }

                length = (length >> 1) + length % 2;
            }

            List<Cluster.ClusterHash> finalClustersHash = clusters_list[0];
            for (int i = 0; i < finalClustersHash.Count; i++)
            {
                output.Add(finalClustersHash[i].Data);  
            }

            return output;
        }

        public static List<Cluster.ClusterHash> DoGradientClusters(Image threshImage, int ts, int y0, int y1,
            int w, int nclustermap, UnionFind uf, List<Cluster.ClusterHash> clusters) 
        {
            ClusterMap clustermap = new ClusterMap(nclustermap, false);

            int mem_chunk_size = 2048;
            ClusterMap[] mem_pools = new ClusterMap[2*nclustermap/mem_chunk_size];
            int mem_pool_idx = 0;
            int mem_pool_loc = 0;
            mem_pools[mem_pool_idx] = new ClusterMap(mem_chunk_size, true);

            for (int y = y0; y < y1; y++) 
            {
                for (int x = 1; x < w-1; x++) 
                {

                    byte v0 = threshImage.GetPixelCustom(y*ts + x);
                    if (v0 == 127)
                        continue;

                    // XXX don't query this until we know we need it?
                    ulong rep0 = uf.GetRepresentative((uint)(y*w + x));
                    if (uf.GetSetSize((uint)rep0) < 25) 
                    {
                        continue;
                    }

                    // whenever we find two adjacent pixels such that one is
                    // white and the other black, we add the point half-way
                    // between them to a cluster associated with the unique
                    // ids of the white and black regions.
                    //
                    // We additionally compute the gradient direction (i.e., which
                    // direction was the white pixel?) Note: if (v1-v0) == 255, then
                    // (dx,dy) points towards the white pixel. if (v1-v0) == -255, then
                    // (dx,dy) points towards the black pixel. p.gx and p.gy will thus
                    // be -255, 0, or 255.
                    //
                    // Note that any given pixel might be added to multiple
                    // different clusters. But in the common case, a given
                    // pixel will be added multiple times to the same cluster,
                    // which increases the size of the cluster and thus the
                    // computational costs.
                    //
                    // A possible optimization would be to combine entries
                    // within the same cluster.

                    // do 4 connectivity. NB: Arguments must be [-1, 1] or we'll overflow .gx, .gy
                    int dx, dy;
                    byte v1;

                    dx = 1;
                    dy = 0;
                    v1 = threshImage.GetPixelCustom((y + dy)*ts + x + dx);       
                                                                                
                    if (v0 + v1 == 255) 
                    {                                   
                        ulong rep1 = (ulong)uf.GetRepresentative((uint)((y + dy)*w + x + dx)); 
                        if (uf.GetSetSize((uint)rep1) > 24) 
                        {        
                            ulong clusterid;                                 
                            if (rep0 < rep1)       
                            {
                                clusterid = (rep1 << 32) + rep0;                
                            }                             
                            else
                            {
                                clusterid = (rep0 << 32) + rep1;                
                            }                                                
                                
                            /* XXX lousy hash function */                       
                            uint clustermap_bucket = (uint)(Utils.Calculations.Hash2(clusterid) % nclustermap); 
                            ClusterMapEntry entry = clustermap.Entries[clustermap_bucket]; 
                            while (entry != null && entry.ID != clusterid) 
                            {           
                                entry = entry.Next;                            
                            }                                                   
                            
                            if (entry == null) 
                            {                                       
                                if (mem_pool_loc == mem_chunk_size) 
                                {           
                                    mem_pool_loc = 0;                           
                                    mem_pool_idx++;                             
                                    mem_pools[mem_pool_idx] = new ClusterMap(mem_chunk_size, true);
                                }         
                                
                                entry = mem_pools[mem_pool_idx].Entries[mem_pool_loc]; 
                                mem_pool_loc++;                                 
                                                                                            
                                entry.ID = clusterid;                          
                                entry.Data =  new Cluster(); 
                                entry.Next = clustermap.Entries[clustermap_bucket];    
                                clustermap.Entries[clustermap_bucket] = entry;          
                            }                                                  
                    
                            Cluster.ClusterPoint p = new Cluster.ClusterPoint((ushort)(2*x + dx),
                                (ushort)(2*y + dy), (short)(dx*((int) v1-v0)), (short)(dy*((int) v1-v0))); 
                            entry.Data.Points.Add(p);                     
                        }                                                   
                    }

                    dx = 0;
                    dy = 1;
                    v1 = threshImage.GetPixelCustom((y + dy)*ts + x + dx);       
                                                                                
                    if (v0 + v1 == 255) 
                    {                                   
                        ulong rep1 = (ulong)uf.GetRepresentative((uint)((y + dy)*w + x + dx)); 
                        if (uf.GetSetSize((uint)rep1) > 24) 
                        {        
                            ulong clusterid;                                 
                            if (rep0 < rep1)       
                            {
                                clusterid = (rep1 << 32) + rep0;                
                            }                             
                            else
                            {
                                clusterid = (rep0 << 32) + rep1;                
                            }                                                
                                
                            /* XXX lousy hash function */                       
                            uint clustermap_bucket = (uint)(Utils.Calculations.Hash2(clusterid) % nclustermap); 
                            ClusterMapEntry entry = clustermap.Entries[clustermap_bucket]; 
                            while (entry != null && entry.ID != clusterid) 
                            {           
                                entry = entry.Next;                            
                            }                                                   
                            
                            if (entry == null) 
                            {                                       
                                if (mem_pool_loc == mem_chunk_size) 
                                {           
                                    mem_pool_loc = 0;                           
                                    mem_pool_idx++;                             
                                    mem_pools[mem_pool_idx] = new ClusterMap(mem_chunk_size, true);
                                }         
                                
                                entry = mem_pools[mem_pool_idx].Entries[mem_pool_loc]; 
                                mem_pool_loc++;                                 
                                                                                            
                                entry.ID = clusterid;                          
                                entry.Data =  new Cluster(); 
                                entry.Next = clustermap.Entries[clustermap_bucket];    
                                clustermap.Entries[clustermap_bucket] = entry;          
                            }                                                  
                    
                            Cluster.ClusterPoint p = new Cluster.ClusterPoint((ushort)(2*x + dx),
                                (ushort)(2*y + dy), (short)(dx*((int) v1-v0)), (short)(dy*((int) v1-v0))); 
                            entry.Data.Points.Add(p);                     
                        }                                                   
                    }

                    // do 8 connectivity
                    dx = -1;
                    dy = 1;
                    v1 = threshImage.GetPixelCustom((y + dy)*ts + x + dx);       
                                                                                
                    if (v0 + v1 == 255) 
                    {                                   
                        ulong rep1 = (ulong)uf.GetRepresentative((uint)((y + dy)*w + x + dx)); 
                        if (uf.GetSetSize((uint)rep1) > 24) 
                        {        
                            ulong clusterid;                                 
                            if (rep0 < rep1)       
                            {
                                clusterid = (rep1 << 32) + rep0;                
                            }                             
                            else
                            {
                                clusterid = (rep0 << 32) + rep1;                
                            }                                                
                                
                            /* XXX lousy hash function */                       
                            uint clustermap_bucket = (uint)(Utils.Calculations.Hash2(clusterid) % nclustermap); 
                            ClusterMapEntry entry = clustermap.Entries[clustermap_bucket]; 
                            while (entry != null && entry.ID != clusterid) 
                            {           
                                entry = entry.Next;                            
                            }                                                   
                            
                            if (entry == null) 
                            {                                       
                                if (mem_pool_loc == mem_chunk_size) 
                                {           
                                    mem_pool_loc = 0;                           
                                    mem_pool_idx++;                             
                                    mem_pools[mem_pool_idx] = new ClusterMap(mem_chunk_size, true);
                                }         
                                
                                entry = mem_pools[mem_pool_idx].Entries[mem_pool_loc]; 
                                mem_pool_loc++;                                 
                                                                                            
                                entry.ID = clusterid;                          
                                entry.Data =  new Cluster(); 
                                entry.Next = clustermap.Entries[clustermap_bucket];    
                                clustermap.Entries[clustermap_bucket] = entry;          
                            }                                                  
                    
                            Cluster.ClusterPoint p = new Cluster.ClusterPoint((ushort)(2*x + dx),
                                (ushort)(2*y + dy), (short)(dx*((int) v1-v0)), (short)(dy*((int) v1-v0))); 
                            entry.Data.Points.Add(p);                     
                        }                                                   
                    }
                    
                    dx = 1;
                    dy = 1;
                    v1 = threshImage.GetPixelCustom((y + dy)*ts + x + dx);       
                                                                                
                    if (v0 + v1 == 255) 
                    {                                   
                        ulong rep1 = (ulong)uf.GetRepresentative((uint)((y + dy)*w + x + dx)); 
                        if (uf.GetSetSize((uint)rep1) > 24) 
                        {        
                            ulong clusterid;                                 
                            if (rep0 < rep1)       
                            {
                                clusterid = (rep1 << 32) + rep0;                
                            }                             
                            else
                            {
                                clusterid = (rep0 << 32) + rep1;                
                            }                                                
                                
                            /* XXX lousy hash function */                       
                            uint clustermap_bucket = (uint)(Utils.Calculations.Hash2(clusterid) % nclustermap); 
                            ClusterMapEntry entry = clustermap.Entries[clustermap_bucket]; 
                            while (entry != null && entry.ID != clusterid) 
                            {           
                                entry = entry.Next;                            
                            }                                                   
                            
                            if (entry == null) 
                            {                                       
                                if (mem_pool_loc == mem_chunk_size) 
                                {           
                                    mem_pool_loc = 0;                           
                                    mem_pool_idx++;                             
                                    mem_pools[mem_pool_idx] = new ClusterMap(mem_chunk_size, true);
                                }         
                                
                                entry = mem_pools[mem_pool_idx].Entries[mem_pool_loc]; 
                                mem_pool_loc++;                                 
                                                                                            
                                entry.ID = clusterid;                          
                                entry.Data =  new Cluster(); 
                                entry.Next = clustermap.Entries[clustermap_bucket];    
                                clustermap.Entries[clustermap_bucket] = entry;          
                            }                                                  
                    
                            Cluster.ClusterPoint p = new Cluster.ClusterPoint((ushort)(2*x + dx),
                                (ushort)(2*y + dy), (short)(dx*((int) v1-v0)), (short)(dy*((int) v1-v0))); 
                            entry.Data.Points.Add(p);                     
                        }                                                   
                    }
                }
            }

            for (int i = 0; i < nclustermap; i++) 
            {
                int start = clusters.Count;
                for (ClusterMapEntry entry = clustermap.Entries[i]; entry != null; entry = entry.Next) 
                {
                    uint hash = (uint)(Utils.Calculations.Hash2(entry.ID) % nclustermap);
                    Cluster.ClusterHash cluster_hash = new Cluster.ClusterHash(hash, entry.ID, entry.Data);
                    clusters.Add(cluster_hash);
                }
                int end = clusters.Count;

                // Do a quick bubblesort on the secondary key.
                int n = end - start;
                for (int j = 0; j < n - 1; j++) 
                {
                    for (int k = 0; k < n - j - 1; k++) 
                    {
                        Cluster.ClusterHash hash1 = clusters[start + k];
                        Cluster.ClusterHash hash2 = clusters[start + k + 1];
                        if (hash1.ID > hash2.ID) 
                        {
                            Cluster.ClusterHash tmp = hash2;
                            hash2 = hash1;
                            hash1 = tmp;
                        }
                    }
                }
            }

            return clusters;
        }

            private static List<Cluster.ClusterHash> mergeClusters(List<Cluster.ClusterHash> c1, List<Cluster.ClusterHash> c2) 
            {
                List<Cluster.ClusterHash> ret = new List<Cluster.ClusterHash>();

                int i1 = 0;
                int i2 = 0;
                int l1 = c1.Count;
                int l2 = c2.Count;

                while (i1 < l1 && i2 < l2) 
                {
                    Cluster.ClusterHash h1 = c1[i1];
                    Cluster.ClusterHash h2 = c2[i2];

                    if (h1.Hash == h2.Hash && h1.ID == h2.ID) 
                    {
                        for (int i = 0; i < h2.Data.Points.Count; i++)
                        {
                            h1.Data.Points.Add(h2.Data.Points[i]);
                        }
                        ret.Add(h1);
                        i1++;
                        i2++;
                    } 
                    else if (h2.Hash < h1.Hash || (h2.Hash == h1.Hash && h2.ID < h1.ID)) 
                    {
                        ret.Add(h2);
                        i2++;
                    } 
                    else 
                    {
                        ret.Add(h1);
                        i1++;
                    }
                }

                for (; i1 < l1; i1++) 
                {
                    ret.Add(c1[i1]);
                }

                for (; i2 < l2; i2++) 
                {
                    ret.Add(c2[i2]);
                }

                return ret;
            }
    }
}
