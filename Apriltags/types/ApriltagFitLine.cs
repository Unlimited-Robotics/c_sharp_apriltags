using System.Collections;
using System.Collections.Generic;
using System;
using System.Linq;

namespace Apriltags
{
    public class FitLine
    {
        public double Mx;
        public double My;
        public double Mxx;
        public double Myy;
        public double Mxy;
        public double W;

        public FitLine()
        {

        }

        public FitLine(FitLine copy)
        {
            Mx = copy.Mx;
            My = copy.My;
            Mxx = copy.Mxx;
            Myy = copy.Myy;
            Mxy = copy.Mxy;
            W = copy.W;
        }

        public static FitLine[] ComputeLfps(Cluster cluster, Image im)
        {
            FitLine[] output = new FitLine[cluster.Points.Count];
            output[0] = new FitLine();
            
            for (int i = 0; i < cluster.Points.Count; i++) 
            {
                // if(cluster.Points.Count == 38)
                // {
                //     Debug.Log("point index " + i);
                // }
                Cluster.ClusterPoint p = cluster.Points[i];

                if (i > 0) 
                {
                    output[i] = new FitLine(output[i - 1]);
                }

                {
                    // we now undo our fixed-point arithmetic.
                    double delta = 0.5; // adjust for pixel center bias
                    double x = p.X * .5 + delta;
                    double y = p.Y * .5 + delta;
                    int ix = (int)x, iy = (int)y;
                    double W = 1;

                    if (ix > 0 && ix+1 < im.Width && iy > 0 && iy+1 < im.Height) 
                    {
                        int grad_x = im.GetPixelCustom(iy * im.Stride + ix + 1) -
                            im.GetPixelCustom(iy * im.Stride + ix - 1);

                        int grad_y = im.GetPixelCustom((iy+1) * im.Stride + ix) -
                            im.GetPixelCustom((iy-1) * im.Stride + ix);

                        // XXX Tunable. How to shape the gradient magnitude?
                        W = Math.Sqrt(grad_x*grad_x + grad_y*grad_y) + 1;
                        // if(cluster.Points.Count == 38)
                        // {
                        //     Debug.Log("x " + x.ToString("F2") + ", y " + y.ToString("F2") +
                        //         ", px " + p.X + ", py " + p.Y + 
                        //         ", grad x " + grad_x + ", grad y " + grad_y + ", w " + W.ToString("F2"));
                        // }
                    }

                    // if(cluster.Points.Count == 38)
                    // {
                    //     Debug.Log("w " + W.ToString("F2") + ", fx " + x + ", fy " + y);
                    // }

                    double fx = x, fy = y;
                    output[i].Mx  += W * fx;
                    output[i].My  += W * fy;
                    output[i].Mxx += W * fx * fx;
                    output[i].Mxy += W * fx * fy;
                    output[i].Myy += W * fy * fy;
                    output[i].W   += W;
                }
            }
            return output;
        }

        public static bool QuadSegmentMaxima(Detector detector, Cluster cluster, FitLine[] lfps, int[] indices)
        {
            int sz = cluster.Points.Count;

            // ksz: when fitting points, how many points on either side do we consider?
            // (actual "kernel" width is 2ksz).
            //
            // This value should be about: 0.5 * (points along shortest edge).
            //
            // If all edges were equally-sized, that would give a value of
            // sz/8. We make it somewhat smaller to account for tags at high
            // aspects.

            // XXX Tunable. Maybe make a multiple of JPEG block size to increase robustness
            // to JPEG compression artifacts?
            int ksz = Utils.Calculations.IMin(20, sz / 12);

            // can't fit a quad if there are too few points.
            if (ksz < 2)
            {
                return false;
            }

            double[] errs = new double[sz];

            for (int i = 0; i < sz; i++)
            {
                double empty = 0;
                FitLines(lfps, sz, (i + sz - ksz) % sz, (i + ksz) % sz, null, out errs[i], out empty);
            }
            // Utils.Log.SaveFitLinesDataToFile("lfps_debug", lfps.ToList());

            // if(sz == 73)
            // {
            //     for (int i = 0; i < errs.Length; i++)
            //     {
            //         Debug.Log("errs[" + i +"] " + errs[i]);
            //     }
            // }

            // apply a low-pass filter to errs
            double[] y = new double[sz];

            // how much filter to apply?

            // XXX Tunable
            double sigma = 1; // was 3

            // cutoff = exp(-j*j/(2*sigma*sigma));
            // log(cutoff) = -j*j / (2*sigma*sigma)
            // log(cutoff)*2*sigma*sigma = -j*j;

            // how big a filter should we use? We make our kernel big
            // enough such that we represent any values larger than
            // 'cutoff'.

            // XXX Tunable (though not super useful to change)
            double cutoff = 0.05;
            int fsz = (int)Math.Sqrt(-Math.Log(cutoff)*2*sigma*sigma) + 1;
            fsz = 2*fsz + 1;

            // For default values of cutoff = 0.05, sigma = 3,
            // we have fsz = 17.
            float[] f = new float[fsz];

            for (int i = 0; i < fsz; i++) 
            {
                int j = i - fsz / 2;
                f[i] = (float)Math.Exp(-j*j/(2*sigma*sigma));
            }

            for (int iy = 0; iy < sz; iy++) 
            {
                double acc = 0;

                for (int i = 0; i < fsz; i++) 
                {
                    acc += errs[(iy + i - fsz / 2 + sz) % sz] * f[i];
                }
                y[iy] = acc;
            }

            for (int i = 0; i < sz; i++)
            {
                errs[i] = y[i];
            }

            int[] maxima = new int[sz];
            double[] maxima_errs = new double[sz];
            int nmaxima = 0;

            for (int i = 0; i < sz; i++) 
            {
                if (errs[i] > errs[(i+1)%sz] && errs[i] > errs[(i+sz-1)%sz]) 
                {
                    maxima[nmaxima] = i;
                    maxima_errs[nmaxima] = errs[i];
                    nmaxima++;
                }
            }

            // if we didn't get at least 4 maxima, we can't fit a quad.
            if (nmaxima < 4)
                return false;

            // select only the best maxima if we have too many
            int max_nmaxima = detector.QuadThreshParams.MaxNMaxima;

            if (nmaxima > max_nmaxima) 
            {
                double[] maxima_errs_copy = new double[nmaxima];
                for (int i = 0; i < nmaxima; i++)
                {
                    maxima_errs_copy[i] = maxima_errs[i];
                }

                // throw out all but the best handful of maxima. Sorts descending.
                // qsort(maxima_errs_copy, nmaxima, sizeof(double), err_compare_descending);
                Utils.Calculations.MaximaComparer comparer = new Utils.Calculations.MaximaComparer();
                Array.Sort(maxima_errs_copy, comparer);

                double maxima_thresh = maxima_errs_copy[max_nmaxima];
                int outP = 0;
                for (int inP = 0; inP < nmaxima; inP++) 
                {
                    if (maxima_errs[inP] <= maxima_thresh)
                        continue;
                    maxima[outP++] = maxima[inP];
                }
                nmaxima = outP;
            }

            int[] best_indices = new int[4];
            double best_error = Utils.Calculations.HUGE_VALF;

            double err01, err12, err23, err30;
            double mse01, mse12, mse23, mse30;
            double[] params01 = new double[4];
            double[] params12 = new double[4]; 
            double[] params23 = new double[4];
            double[] params30 = new double[4];

            // disallow quads where the angle is less than a critical value.
            double max_dot = detector.QuadThreshParams.CosCriticalRad; //25*M_PI/180);

            for (int m0 = 0; m0 < nmaxima - 3; m0++) {
                int i0 = maxima[m0];

                for (int m1 = m0+1; m1 < nmaxima - 2; m1++) {
                    int i1 = maxima[m1];

                    FitLines(lfps, sz, i0, i1, params01, out err01, out mse01);

                    if (mse01 > detector.QuadThreshParams.MaxLineFitMse)
                    {
                        continue;
                    }

                    for (int m2 = m1+1; m2 < nmaxima - 1; m2++) 
                    {
                        int i2 = maxima[m2];

                        FitLines(lfps, sz, i1, i2, params12, out err12, out mse12);
                        if (mse12 > detector.QuadThreshParams.MaxLineFitMse)
                        {
                            continue;
                        }

                        double dot = params01[2]*params12[2] + params01[3]*params12[3];
                        if ((float)Math.Abs(dot) > max_dot)
                        {
                            continue;
                        }

                        for (int m3 = m2+1; m3 < nmaxima; m3++) {
                            int i3 = maxima[m3];

                            FitLines(lfps, sz, i2, i3, params23, out err23, out mse23);
                            if (mse23 > detector.QuadThreshParams.MaxLineFitMse)
                            {
                                continue;
                            }

                            FitLines(lfps, sz, i3, i0, params30, out err30, out mse30);
                            if (mse30 > detector.QuadThreshParams.MaxLineFitMse)
                            {
                                continue;
                            }

                            double err = err01 + err12 + err23 + err30;
                            if (err < best_error) {
                                best_error = err;
                                best_indices[0] = i0;
                                best_indices[1] = i1;
                                best_indices[2] = i2;
                                best_indices[3] = i3;
                            }
                        }
                    }
                }
            }

            if (best_error == Utils.Calculations.HUGE_VALF)
                return false;

            for (int i = 0; i < 4; i++)
                indices[i] = best_indices[i];

            if (best_error / sz < detector.QuadThreshParams.MaxLineFitMse)
                return true;
            return false;
        }

        public static void FitLines(FitLine[] lfps, int sz, int i0, int i1, double[] lineparm, out double err, out double mse)
        {
            // if(i0 == i1)
            // {
            //     err = 0;
            //     mse = 0;
            //     return;
            // }
            // if((i0 >= 0 && i1 >= 0 && i0 < sz && i1 < sz) == false)
            // {
            //     err = 0;
            //     mse = 0;
            //     return;
            // }

            double Mx, My, Mxx, Myy, Mxy, W;
            int N; // how many points are included in the set?
            // Debug.Log("i0 " + i0 + ", i1 " + i1);
            if (i0 < i1) {
                N = i1 - i0 + 1;

                Mx  = lfps[i1].Mx;
                My  = lfps[i1].My;
                Mxx = lfps[i1].Mxx;
                Mxy = lfps[i1].Mxy;
                Myy = lfps[i1].Myy;
                W   = lfps[i1].W;

                if (i0 > 0) {
                    Mx  -= lfps[i0-1].Mx;
                    My  -= lfps[i0-1].My;
                    Mxx -= lfps[i0-1].Mxx;
                    Mxy -= lfps[i0-1].Mxy;
                    Myy -= lfps[i0-1].Myy;
                    W   -= lfps[i0-1].W;
                }

            } else {
                // i0 > i1, e.g. [15, 2]. Wrap around.
                // if((i0>0) == false)
                // {
                //     err = 0;
                //     mse = 0;
                //     return;
                // }

                Mx  = lfps[sz-1].Mx   - lfps[i0-1].Mx;
                My  = lfps[sz-1].My   - lfps[i0-1].My;
                Mxx = lfps[sz-1].Mxx  - lfps[i0-1].Mxx;
                Mxy = lfps[sz-1].Mxy  - lfps[i0-1].Mxy;
                Myy = lfps[sz-1].Myy  - lfps[i0-1].Myy;
                W   = lfps[sz-1].W    - lfps[i0-1].W;

                Mx  += lfps[i1].Mx;
                My  += lfps[i1].My;
                Mxx += lfps[i1].Mxx;
                Mxy += lfps[i1].Mxy;
                Myy += lfps[i1].Myy;
                W   += lfps[i1].W;

                // Debug.Log("lfps[sz-1].Mx " + lfps[sz-1].Mx + ", lfps[i0-1].Mx " + lfps[i0-1].Mx + ", lfps[i1].Mx " + lfps[i1].Mx);
                N = sz - i0 + i1 + 1;
            }
            // Debug.Log("mx " + Mx.ToString("F2") + ", my " + My.ToString("F2") + 
            //     ", mxx " + Mxx.ToString("F2") + ", myy " + Myy.ToString("F2") +
            //     ", mxy " + Mxy.ToString("F2") + ", w " + W.ToString("F2") + ", n " + N);

            double Ex = Mx / W;
            double Ey = My / W;
            double Cxx = Mxx / W - Ex*Ex;
            double Cxy = Mxy / W - Ex*Ey;
            double Cyy = Myy / W - Ey*Ey;

            //if (1) {
            //    // on iOS about 5% of total CPU spent in these trig functions.
            //    // 85 ms per frame on 5S, example.pnm
            //    //
            //    // XXX this was using the double-precision atan2. Was there a case where
            //    // we needed that precision? Seems doubtful.
            //    double normal_theta = .5 * atan2f(-2*Cxy, (Cyy - Cxx));
            //    nx_old = cosf(normal_theta);
            //    ny_old = sinf(normal_theta);
            //}

            // Instead of using the above cos/sin method, pose it as an eigenvalue problem.
            double eig_small = 0.5*(Cxx + Cyy - (float)Math.Sqrt((Cxx - Cyy)*(Cxx - Cyy) + 4*Cxy*Cxy));

            if (lineparm != null) 
            {
                lineparm[0] = Ex;
                lineparm[1] = Ey;

                double eig = 0.5*(Cxx + Cyy + (float)Math.Sqrt((Cxx - Cyy)*(Cxx - Cyy) + 4*Cxy*Cxy));
                double nx1 = Cxx - eig;
                double ny1 = Cxy;
                double M1 = nx1*nx1 + ny1*ny1;
                double nx2 = Cxy;
                double ny2 = Cyy - eig;
                double M2 = nx2*nx2 + ny2*ny2;

                double nx, ny, M;
                if (M1 > M2) {
                    nx = nx1;
                    ny = ny1;
                    M = M1;
                } else {
                    nx = nx2;
                    ny = ny2;
                    M = M2;
                }

                double length = (float)Math.Sqrt(M);
                lineparm[2] = nx/length;
                lineparm[3] = ny/length;
            }

            // sum of squared errors =
            //
            // SUM_i ((p_x - ux)*nx + (p_y - uy)*ny)^2
            // SUM_i  nx*nx*(p_x - ux)^2 + 2nx*ny(p_x -ux)(p_y-uy) + ny*ny*(p_y-uy)*(p_y-uy)
            //  nx*nx*SUM_i((p_x -ux)^2) + 2nx*ny*SUM_i((p_x-ux)(p_y-uy)) + ny*ny*SUM_i((p_y-uy)^2)
            //
            //  nx*nx*N*Cxx + 2nx*ny*N*Cxy + ny*ny*N*Cyy

            // sum of squared errors
            err = N*eig_small;

            // mean squared error
            mse = eig_small;
        }
    }
}
