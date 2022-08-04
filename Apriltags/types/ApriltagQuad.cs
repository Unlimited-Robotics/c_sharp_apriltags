using System.Collections;
using System.Collections.Generic;
using System;
using System.Linq;

namespace Apriltags
{
    public class Quad
    {
        public Vector2[] Corners;
        public bool ReversedBorder;

        // H: tag coordinates ([-1,1] at the black corners) to pixels
        public Matd H;
        // Hinv: pixels to tag
        public Matd Hinv;

        public Quad()
        {
            Corners = new Vector2[4];
            Corners[0] = new Vector2();
            Corners[1] = new Vector2();
            Corners[2] = new Vector2();
            Corners[3] = new Vector2();
        }

        public Quad(Quad copy)
        {
            Corners = new Vector2[4];
            Corners[0] = new Vector2(copy.Corners[0].x, copy.Corners[0].y);
            Corners[1] = new Vector2(copy.Corners[1].x, copy.Corners[1].y);
            Corners[2] = new Vector2(copy.Corners[2].x, copy.Corners[2].y);
            Corners[3] = new Vector2(copy.Corners[3].x, copy.Corners[3].y);
            ReversedBorder = copy.ReversedBorder;
            H = copy.H;
            Hinv = copy.Hinv;
        }

        public class QuadTask : WorkerPool.WorkTask
        {
            public List<Cluster> Clusters;
            public int Cidx0;
            public int Cidx1;
            public List<Quad> Quads;
            public Detector ATDetector;
            public int Width;
            public int Height;
            public Image ATImage;
            public int TagWidth;
            public bool NormalBorder;
            public bool ReversedBorder;

            public override void DoTask()
            {
                for (int cidx = Cidx0; cidx < Cidx1; cidx++) 
                {
                    Cluster cluster = Clusters[cidx];

                    if (cluster.Points.Count < ATDetector.QuadThreshParams.MinClusterPixels)
                        continue;

                    // a cluster should contain only boundary points around the
                    // tag. it cannot be bigger than the whole screen. (Reject
                    // large connected blobs that will be prohibitively slow to
                    // fit quads to.) A typical point along an edge is added three
                    // times (because it has 3 neighbors). The maximum perimeter
                    // is 2w+2h.
                    if (cluster.Points.Count > 3*(2*Width+2*Height)) 
                    {
                        continue;
                    }

                    Quad quad = new Quad();
                    // memset(&quad, 0, sizeof(struct quad));

                    if (fitQuad(ATDetector, ATImage, cluster, quad, TagWidth, NormalBorder, ReversedBorder)) 
                    {
                        ATDetector.MutexLock.WaitOne();
                        Quads.Add(quad);
                        ATDetector.MutexLock.ReleaseMutex();
                    }
                }
            }

            private bool fitQuad(Detector detector, Image image, Cluster cluster, Quad quad,
                int tagWidth, bool normalBorder, bool reversedBorder) 
            {
                bool res = false;

                int sz = cluster.Points.Count;
                if (sz < 24) // Synchronize with later check.
                {
                    return false;
                }

                /////////////////////////////////////////////////////////////
                // Step 1. Sort points so they wrap around the center of the
                // quad. We will constrain our quad fit to simply partition this
                // ordered set into 4 groups.

                // compute a bounding box so that we can order the points
                // according to their angle WRT the center.
                Cluster.ClusterPoint p1 = cluster.Points[0];
                ushort xmax = p1.X;
                ushort xmin = p1.X;
                ushort ymax = p1.Y;
                ushort ymin = p1.Y;
                for (int pidx = 1; pidx < cluster.Points.Count; pidx++) 
                {
                    Cluster.ClusterPoint point = cluster.Points[pidx];

                    if (point.X > xmax) 
                    {
                        xmax = point.X;
                    } 
                    else if (point.X < xmin) 
                    {
                        xmin = point.X;
                    }

                    if (point.Y > ymax) 
                    {
                        ymax = point.Y;
                    } 
                    else if (point.Y < ymin) 
                    {
                        ymin = point.Y;
                    }
                }

                if ((xmax - xmin)*(ymax - ymin) < tagWidth) 
                {
                    return false;
                }

                // add some noise to (cx,cy) so that pixels get a more diverse set
                // of theta estimates. This will help us remove more points.
                // (Only helps a small amount. The actual noise values here don't
                // matter much at all, but we want them [-1, 1]. (XXX with
                // fixed-point, should range be bigger?)
                double cx = (xmin + xmax) * 0.5f + 0.05118f;
                double cy = (ymin + ymax) * 0.5f + -0.028581f;
                // if(sz == 577)
                // {
                //     Debug.Log("ymax " + ymax);
                //     Debug.Log("ymin " + ymin);
                //     Debug.Log("cx " + cx);
                //     Debug.Log("cy " + cy);
                // }

                double dot = 0;

                // float quadrants[2][2] = {{-1*(2 << 15), 0}, {2*(2 << 15), 2 << 15}};
                float[][] quadrants = new float[2][];
                quadrants[0] = new float[]{-1*(2 << 15), 0};
                quadrants[1] = new float[]{2*(2 << 15), 2 << 15};

                for (int pidx = 0; pidx < cluster.Points.Count; pidx++) 
                {
                    // if(sz == 577)
                    // {
                    //     Debug.Log("index " + pidx);
                    // }
                    Cluster.ClusterPoint point = cluster.Points[pidx];

                    double dx = point.X - cx;
                    double dy = point.Y - cy;

                    dot += dx*point.GX + dy*point.GY;
                    // if(sz == 577)
                    // {
                    //     Debug.Log("dx before " + dx);
                    //     Debug.Log("dy before " + dy);
                    //     Debug.Log("dot " + dot);
                    // }

                    // double quadrant = getQuadrant((int)dx, (int)dy);
                    double quadrant = quadrants[Convert.ToInt32(dy > 0)][Convert.ToInt32(dx > 0)];
                    if (dy < 0) 
                    {
                        dy = -dy;
                        dx = -dx;
                    }

                    if (dx < 0) 
                    {
                        double tmp = dx;
                        dx = dy;
                        dy = -tmp;
                    }
                    point.Slope = quadrant + dy/dx;
                    // if(sz == 577)
                    // {
                    //     Debug.Log("dx after " + dx);
                    //     Debug.Log("dy after " + dy);
                    //     Debug.Log("quadrant " + quadrant);
                    //     Debug.Log("slope " + point.Slope);
                    // }
                }

                               // Debug.Log("sz " + sz);
                // if(sz == 577)
                // {
                //     Utils.Log.SaveClustersDataToFile("debug_cluster", cluster);
                //     // Debug.Log("copy");
                // }

                // Ensure that the black border is inside the white border.
                quad.ReversedBorder = dot < 0;
                if (!reversedBorder && quad.ReversedBorder) 
                {
                    return false;
                }
                if (!normalBorder && !quad.ReversedBorder) 
                {
                    return false;
                }

                // we now sort the points according to theta. This is a prepatory
                // step for segmenting them into four lines.
                // Utils.Calculations.Plaster1(cluster);
                cluster.PtSort();

                // remove duplicate points. (A byproduct of our segmentation system.)
                int outpos = 1;

                Cluster.ClusterPoint last = cluster.Points[0];

                for (int i = 1; i < sz; i++) 
                {
                    Cluster.ClusterPoint point = cluster.Points[i];

                    if (point.X != last.X || point.Y != last.Y) 
                    {

                        if (i != outpos)  
                        {
                            cluster.Points[outpos].Copy(point);
                        }

                        outpos++;
                    }

                    last = point;
                }

                for (int i = cluster.Points.Count - 1; i >= outpos; i--)
                {
                    cluster.Points.RemoveAt(i);
                }
                sz = outpos;
                // Debug.Log("count " + cluster.Points.Count);
                // Debug.Log("sz " + sz);

                if (sz < 24)
                {
                    return false;
                }
                

                FitLine[] lfps = FitLine.ComputeLfps(cluster, image);
                // if(sz == 38)
                // {
                //     Utils.Log.SaveFitLinesDataToFile("lfps_debug", lfps.ToList());  
                // }

                int[] indices = new int[4];
                if (FitLine.QuadSegmentMaxima(detector, cluster, lfps, indices) == false)
                {
                    return res;
                }

                // Debug.Log("best indices " + indices[0] + ", " + indices[1] + ", " + indices[2] + ", " + indices[3]);

                double[][] lines = new double[4][];
                lines[0] = new double[4];
                lines[1] = new double[4];
                lines[2] = new double[4];
                lines[3] = new double[4];

                for (int i = 0; i < 4; i++) 
                {
                    int i0 = indices[i];
                    int i1 = indices[(i+1)&3];

                    double err;
                    double empty;
                    FitLine.FitLines(lfps, sz, i0, i1, lines[i], out empty, out err);

                    if (err > detector.QuadThreshParams.MaxLineFitMse) 
                    {
                        res = false;
                        return res;
                    }
                }

                for (int i = 0; i < 4; i++) 
                {
                    // solve for the intersection of lines (i) and (i+1)&3.
                    // p0 + lambda0*u0 = p1 + lambda1*u1, where u0 and u1
                    // are the line directions.
                    //
                    // lambda0*u0 - lambda1*u1 = (p1 - p0)
                    //
                    // rearrange (solve for lambdas)
                    //
                    // [u0_x   -u1_x ] [lambda0] = [ p1_x - p0_x ]
                    // [u0_y   -u1_y ] [lambda1]   [ p1_y - p0_y ]
                    //
                    // remember that lines[i][0,1] = p, lines[i][2,3] = NORMAL vector.
                    // We want the unit vector, so we need the perpendiculars. Thus, below
                    // we have swapped the x and y components and flipped the y components.

                    double A00 =  lines[i][3],  A01 = -lines[(i+1)&3][3];
                    double A10 =  -lines[i][2],  A11 = lines[(i+1)&3][2];
                    double B0 = -lines[i][0] + lines[(i+1)&3][0];
                    double B1 = -lines[i][1] + lines[(i+1)&3][1];

                    double det = A00 * A11 - A10 * A01;

                    // inverse.
                    double W00 = A11 / det, W01 = -A01 / det;
                    if ((float)Math.Abs(det) < 0.001) 
                    {
                        res = false;
                        return res;
                    }

                    // solve
                    double L0 = W00*B0 + W01*B1;

                    // compute intersection
                    quad.Corners[i].x = (float)(lines[i][0] + L0*A00);
                    quad.Corners[i].y = (float)(lines[i][1] + L0*A10);

                    res = true;
                }

                // reject quads that are too small
                double area = 0;

                // get area of triangle formed by points 0, 1, 2, 0
                double[] length = new double[3];
                double p;
                for (int i = 0; i < 3; i++) 
                {
                    int idxa = i; // 0, 1, 2,
                    int idxb = (i+1) % 3; // 1, 2, 0
                    float xVal = quad.Corners[idxb].x - quad.Corners[idxa].x;
                    float yVal = quad.Corners[idxb].y - quad.Corners[idxa].y;
                    length[i] = Math.Sqrt(xVal*xVal + yVal*yVal);
                }
                p = (length[0] + length[1] + length[2]) / 2;

                area += Math.Sqrt(p*(p-length[0])*(p-length[1])*(p-length[2]));

                // get area of triangle formed by points 2, 3, 0, 2
                for (int i = 0; i < 3; i++) 
                {
                    int[] idxs = new int[]{ 2, 3, 0, 2 };
                    int idxa = idxs[i];
                    int idxb = idxs[i+1];
                    float xVal = quad.Corners[idxb].x - quad.Corners[idxa].x;
                    float yVal = quad.Corners[idxb].y - quad.Corners[idxa].y;
                    length[i] = Math.Sqrt(xVal*xVal + yVal*yVal);
                }
                p = (length[0] + length[1] + length[2]) / 2;

                area += Math.Sqrt(p*(p-length[0])*(p-length[1])*(p-length[2]));

                if (area < tagWidth*tagWidth) 
                {
                    res = false;
                    return res;
                }

                // reject quads whose cumulative angle change isn't equal to 2PI
                for (int i = 0; i < 4; i++) 
                {
                    int i0 = i, i1 = (i+1)&3, i2 = (i+2)&3;

                    double dx1 = quad.Corners[i1].x - quad.Corners[i0].x;
                    double dy1 = quad.Corners[i1].y - quad.Corners[i0].y;
                    double dx2 = quad.Corners[i2].x - quad.Corners[i1].x;
                    double dy2 = quad.Corners[i2].y - quad.Corners[i1].y;
                    double cos_dtheta = (dx1*dx2 + dy1*dy2)/Math.Sqrt((dx1*dx1 + dy1*dy1)*(dx2*dx2 + dy2*dy2));

                    float cosCritRad = detector.QuadThreshParams.CosCriticalRad;
                    if ((cos_dtheta > cosCritRad || cos_dtheta < -cosCritRad) || dx1*dy2 < dy1*dx2) 
                    {
                        res = false;
                        return res;
                    }
                }

                return res;
            }

            private double getQuadrant(int dx, int dy)
            {
                double output;
                
                if(dy > 0)
                {
                    if(dx > 0)
                    {
                        output = 2 << 15;
                    }
                    else
                    {
                        output = 2 * (2 << 15);
                    }
                }
                else
                {
                    if(dx > 0)
                    {
                        output = 0;
                    }
                    else
                    {
                        output = -1 * (2 << 15);
                    }
                }

                return output;
            }
        }
    }
}
