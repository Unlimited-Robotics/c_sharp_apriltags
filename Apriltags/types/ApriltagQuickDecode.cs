using System.Collections;
using System.Collections.Generic;
using System;

namespace Apriltags
{
    public class QuickDecode
    {
        public QuickDecodeEntry[] Entries;

        public QuickDecode(int entryAmount)
        {
            Entries = new QuickDecodeEntry[entryAmount];
            for (int i = 0; i < Entries.Length; i++)
            {
                Entries[i] = new QuickDecodeEntry();
            }
        }

        public void AddCode(ulong code, ushort id, byte hamming)
        {
            int bucket = (int)(code % (ulong)Entries.Length);
            while(Entries[bucket].RCode != ulong.MaxValue)
            {
                bucket = (bucket + 1) % Entries.Length;
            }

            Entries[bucket].RCode = code;
            Entries[bucket].ID = id;
            Entries[bucket].Hamming = hamming;
        }

        public class QuickDecodeTask : WorkerPool.WorkTask
        {
            public int I0;
            public int I1;

            public Image ATImage;
            public List<Quad> Quads;
            public Detector ATDetector;
            public List<Detection> Detections;

            public override void DoTask()
            {
                for (int quadidx = I0; quadidx < I1; quadidx++) 
                {
                    Quad quad_original = Quads[quadidx];

                    // refine edges is not dependent upon the tag family, thus
                    // apply this optimization BEFORE the other work.
                    //if (td->quad_decimate > 1 && td->refine_edges) {
                    if (ATDetector.RefineEdges) 
                    {
                        refineEdges(quad_original);
                    }

                    // make sure the homographies are computed...
                    if (quadUpdateHomographies(quad_original) != 0)
                    {
                        continue;
                    }

                    for (int famidx = 0; famidx < ATDetector.TagFamilies.Count; famidx++) 
                    {
                        ApriltagFamily family = ATDetector.TagFamilies[famidx];

                        if (family.ReversedBorder != quad_original.ReversedBorder) 
                        {
                            continue;
                        }

                        // since the geometry of tag families can vary, start any
                        // optimization process over with the original quad.
                        Quad quad = new Quad(quad_original);

                        QuickDecodeEntry entry;

                        float decision_margin = quadDecode(family, quad, out entry);

                        if (decision_margin >= 0 && entry.Hamming < 255) 
                        {
                            Detection det = new Detection();

                            det.Family = family;
                            det.ID = entry.ID;
                            det.Hamming = entry.Hamming;
                            det.DecisionMargin = decision_margin;

                            double theta = entry.Rotation * Math.PI / 2.0;
                            double c = Math.Cos(theta), s = Math.Sin(theta);

                            // Fix the rotation of our homography to properly orient the tag
                            Matd R = new Matd(3,3);
                            R.SetCell(0,0,c);
                            R.SetCell(0,1,-s);
                            R.SetCell(1,0,s);
                            R.SetCell(1,1,c);
                            R.SetCell(2,2,1);

                            det.H = Matd.MatdOp("M*M", new List<Matd>(){quad.H, R});

                            Utils.Homography.HomographyProject(det.H, 0, 0, out det.Center[0], out det.Center[1]);

                            // [-1, -1], [1, -1], [1, 1], [-1, 1], Desired points
                            // [-1, 1], [1, 1], [1, -1], [-1, -1], FLIP Y
                            // adjust the points in det->p so that they correspond to
                            // counter-clockwise around the quad, starting at -1,-1.
                            for (int i = 0; i < 4; i++) 
                            {
                                int tcx = (i == 1 || i == 2) ? 1 : -1;
                                int tcy = (i < 2) ? 1 : -1;

                                double[] p = new double[2];

                                Utils.Homography.HomographyProject(det.H, tcx, tcy, out p[0], out p[1]);

                                det.Corners[i][0] = p[0];
                                det.Corners[i][1] = p[1];
                            }

                            ATDetector.MutexLock.WaitOne();
                            Detections.Add(det);
                            ATDetector.MutexLock.ReleaseMutex();
                        }
                    }
                }
            }

            private void refineEdges(Quad quad)
            {
                double[][] lines = new double[4][]; // for each line, [Ex Ey nx ny]
                lines[0] = new double[4];
                lines[1] = new double[4];
                lines[2] = new double[4];
                lines[3] = new double[4];

                for (int edge = 0; edge < 4; edge++) 
                {
                    int a = edge, b = (edge + 1) & 3; // indices of the end points.

                    // compute the normal to the current line estimate
                    double nx = quad.Corners[b].y - quad.Corners[a].y;
                    double ny = -quad.Corners[b].x + quad.Corners[a].x;
                    double mag = Math.Sqrt(nx*nx + ny*ny);
                    nx /= mag;
                    ny /= mag;

                    if (quad.ReversedBorder) 
                    {
                        nx = -nx;
                        ny = -ny;
                    }

                    // we will now fit a NEW line by sampling points near
                    // our original line that have large gradients. On really big tags,
                    // we're willing to sample more to get an even better estimate.
                    int nsamples = Utils.Calculations.IMax(16, (int)(mag / 8)); // XXX tunable

                    // stats for fitting a line...
                    double Mx = 0, My = 0, Mxx = 0, Mxy = 0, Myy = 0, N = 0;

                    for (int s = 0; s < nsamples; s++) {
                        // compute a point along the line... Note, we're avoiding
                        // sampling *right* at the corners, since those points are
                        // the least reliable.
                        double alpha = (1.0 + s) / (nsamples + 1);
                        double x0 = alpha*quad.Corners[a].x + (1-alpha)*quad.Corners[b].x;
                        double y0 = alpha*quad.Corners[a].y + (1-alpha)*quad.Corners[b].y;

                        // search along the normal to this line, looking at the
                        // gradients along the way. We're looking for a strong
                        // response.
                        double Mn = 0;
                        double Mcount = 0;

                        // XXX tunable: how far to search?  We want to search far
                        // enough that we find the best edge, but not so far that
                        // we hit other edges that aren't part of the tag. We
                        // shouldn't ever have to search more than quad_decimate,
                        // since otherwise we would (ideally) have started our
                        // search on another pixel in the first place. Likewise,
                        // for very small tags, we don't want the range to be too
                        // big.
                        double range = ATDetector.QuadDecimate + 1;

                        // XXX tunable step size.
                        for (double n = -range; n <= range; n +=  0.25) {
                            // Because of the guaranteed winding order of the
                            // points in the quad, we will start inside the white
                            // portion of the quad and work our way outward.
                            //
                            // sample to points (x1,y1) and (x2,y2) XXX tunable:
                            // how far +/- to look? Small values compute the
                            // gradient more precisely, but are more sensitive to
                            // noise.
                            double grange = 1;
                            int x1 = (int)(x0 + (n + grange)*nx);
                            int y1 = (int)(y0 + (n + grange)*ny);
                            if (x1 < 0 || x1 >= ATImage.Width || y1 < 0 || y1 >= ATImage.Height)
                                continue;

                            int x2 = (int)(x0 + (n - grange)*nx);
                            int y2 = (int)(y0 + (n - grange)*ny);
                            if (x2 < 0 || x2 >= ATImage.Width || y2 < 0 || y2 >= ATImage.Height)
                                continue;

                            int g1 = ATImage.GetPixelCustom(y1*ATImage.Stride + x1);
                            int g2 = ATImage.GetPixelCustom(y2*ATImage.Stride + x2);

                            if (g1 < g2) // reject points whose gradient is "backwards". They can only hurt us.
                                continue;

                            double weight = (g2 - g1)*(g2 - g1); // XXX tunable. What shape for weight=f(g2-g1)?

                            // compute weighted average of the gradient at this point.
                            Mn += weight*n;
                            Mcount += weight;
                        }

                        // what was the average point along the line?
                        if (Mcount == 0)
                            continue;

                        double n0 = Mn / Mcount;

                        // where is the point along the line?
                        double bestx = x0 + n0*nx;
                        double besty = y0 + n0*ny;

                        // update our line fit statistics
                        Mx += bestx;
                        My += besty;
                        Mxx += bestx*bestx;
                        Mxy += bestx*besty;
                        Myy += besty*besty;
                        N++;
                    }

                    // fit a line
                    double Ex = Mx / N, Ey = My / N;
                    double Cxx = Mxx / N - Ex*Ex;
                    double Cxy = Mxy / N - Ex*Ey;
                    double Cyy = Myy / N - Ey*Ey;

                    // TODO: Can replace this with same code as in fit_line.
                    double normal_theta = .5 * Math.Atan2(-2*Cxy, (Cyy - Cxx));
                    nx = Math.Cos(normal_theta);
                    ny = Math.Sin(normal_theta);
                    lines[edge][0] = Ex;
                    lines[edge][1] = Ey;
                    lines[edge][2] = nx;
                    lines[edge][3] = ny;
                }

                // now refit the corners of the quad
                for (int i = 0; i < 4; i++) 
                {

                    // solve for the intersection of lines (i) and (i+1)&3.
                    double A00 =  lines[i][3],  A01 = -lines[(i+1)&3][3];
                    double A10 =  -lines[i][2],  A11 = lines[(i+1)&3][2];
                    double B0 = -lines[i][0] + lines[(i+1)&3][0];
                    double B1 = -lines[i][1] + lines[(i+1)&3][1];

                    double det = A00 * A11 - A10 * A01;

                    // inverse.
                    if (Math.Abs(det) > 0.001) {
                        // solve
                        double W00 = A11 / det, W01 = -A01 / det;

                        double L0 = W00*B0 + W01*B1;

                        // compute intersection
                        quad.Corners[i].x = (float)(lines[i][0] + L0*A00);
                        quad.Corners[i].y = (float)(lines[i][1] + L0*A10);
                    } else {
                        // this is a bad sign. We'll just keep the corner we had.
            //            printf("bad det: %15f %15f %15f %15f %15f\n", A00, A11, A10, A01, det);
                    }
                }
            }

            private int quadUpdateHomographies(Quad quad)
            {
                // Debug.Log("quad homo");
                //zarray_t *correspondences = zarray_create(sizeof(float[4]));

                double[][] corr_arr = new double[4][];
                corr_arr[0] = new double[4];
                corr_arr[1] = new double[4];
                corr_arr[2] = new double[4];
                corr_arr[3] = new double[4];

                for (int i = 0; i < 4; i++) 
                {
                    corr_arr[i][0] = (i==0 || i==3) ? -1 : 1;
                    corr_arr[i][1] = (i==0 || i==1) ? -1 : 1;
                    corr_arr[i][2] = quad.Corners[i].x;
                    corr_arr[i][3] = quad.Corners[i].y;
                }

                quad.H = null;
                quad.Hinv = null;

                // XXX Tunable
                quad.H = Utils.Homography.HomographyCompute2(corr_arr);

                quad.Hinv = quad.H.GetInverseMatrix();

                if (quad.H != null && quad.Hinv != null)
                {
                    return 0;
                }

                return -1;
            }

            private float quadDecode(ApriltagFamily family, Quad quad, out QuickDecodeEntry entry)
            {
                // decode the tag binary contents by sampling the pixel
                // closest to the center of each bit cell.

                // We will compute a threshold by sampling known white/black cells around this tag.
                // This sampling is achieved by considering a set of samples along lines.
                //
                // coordinates are given in bit coordinates. ([0, fam->border_width]).
                //
                // { initial x, initial y, delta x, delta y, WHITE=1 }
                float[] patterns = new float[] 
                {
                    // left white column
                    -0.5f, 0.5f,
                    0, 1,
                    1,

                    // left black column
                    0.5f, 0.5f,
                    0, 1,
                    0,

                    // right white column
                    family.WidthAtBorder + 0.5f, 0.5f,
                    0, 1,
                    1,

                    // right black column
                    family.WidthAtBorder - 0.5f, 0.5f,
                    0, 1,
                    0,

                    // top white row
                    0.5f, -0.5f,
                    1, 0,
                    1,

                    // top black row
                    0.5f, 0.5f,
                    1, 0,
                    0,

                    // bottom white row
                    0.5f, family.WidthAtBorder + 0.5f,
                    1, 0,
                    1,

                    // bottom black row
                    0.5f, family.WidthAtBorder - 0.5f,
                    1, 0,
                    0

                    // XXX double-counts the corners.
                };

                GrayModel whitemodel = new GrayModel(), blackmodel = new GrayModel();

                for (int pattern_idx = 0; pattern_idx < patterns.Length/5; pattern_idx ++) 
                {
                    int size = patterns.Length-pattern_idx * 5;
                    float[] pattern = new float[size];
                    Array.Copy(patterns, pattern_idx*5, pattern, 0, size);
                    // float *pattern = &patterns[pattern_idx * 5];

                    int is_white = (int)pattern[4];

                    for (int i = 0; i < family.WidthAtBorder; i++) 
                    {
                        double tagx01 = (pattern[0] + i*pattern[2]) / (family.WidthAtBorder);
                        double tagy01 = (pattern[1] + i*pattern[3]) / (family.WidthAtBorder);

                        double tagx = 2*(tagx01-0.5);
                        double tagy = 2*(tagy01-0.5);

                        double px, py;
                        Utils.Homography.HomographyProject(quad.H, tagx, tagy, out px, out py);

                        // don't round
                        int ix = (int)px;
                        int iy = (int)py;
                        if (ix < 0 || iy < 0 || ix >= ATImage.Width || iy >= ATImage.Height)
                        {
                            continue;
                        }

                        int v = ATImage.GetPixelCustom(iy*ATImage.Stride + ix);

                        if (is_white != 0)
                        {
                            whitemodel.Add(tagx, tagy, v);
                        }
                        else
                        {
                            blackmodel.Add(tagx, tagy, v);
                        }
                    }
                }

                whitemodel.Solve();
                blackmodel.Solve();

                // XXX Tunable
                if ((whitemodel.Interpolate(0,0) - blackmodel.Interpolate(0,0) < 0) != family.ReversedBorder) 
                {
                    entry = new QuickDecodeEntry();
                    entry.RCode = 0;
                    entry.ID = 65535;
                    entry.Hamming = 255;
                    entry.Rotation = 0;
                    return -1;
                }

                // compute the average decision margin (how far was each bit from
                // the decision boundary?
                //
                // we score this separately for white and black pixels and return
                // the minimum average threshold for black/white pixels. This is
                // to penalize thresholds that are too close to an extreme.
                float black_score = 0, white_score = 0;
                float black_score_count = 1, white_score_count = 1;

                double[] values = new double[family.TotalWidth*family.TotalWidth];

                int min_coord = (family.WidthAtBorder - family.TotalWidth)/2;
                for (int i = 0; i < family.BitX.Length; i++) 
                {
                    int bity = (int)family.BitY[i];
                    int bitx = (int)family.BitX[i];

                    double tagx01 = (bitx + 0.5) / (family.WidthAtBorder);
                    double tagy01 = (bity + 0.5) / (family.WidthAtBorder);

                    // scale to [-1, 1]
                    double tagx = 2*(tagx01-0.5);
                    double tagy = 2*(tagy01-0.5);

                    double px, py;
                    Utils.Homography.HomographyProject(quad.H, tagx, tagy, out px, out py);

                    double v = ATImage.ValueForPixel(px, py);

                    if (v == -1) 
                    {
                        continue;
                    }

                    double thresh = (blackmodel.Interpolate(tagx, tagy) + whitemodel.Interpolate(tagx, tagy)) / 2.0;
                    values[family.TotalWidth*(bity - min_coord) + bitx - min_coord] = v - thresh;
                }

                Utils.Calculations.Sharpen(ATDetector, values, family.TotalWidth);
                // for (int i = 0; i < values.Length; i++)
                // {
                //     Debug.Log("value " + i + " " + values[i]);
                // }


                ulong rcode = 0;
                for (int i = 0; i < family.BitX.Length; i++) 
                {
                    int bity = (int)family.BitY[i];
                    int bitx = (int)family.BitX[i];
                    rcode = (rcode << 1);
                    double v = values[(bity - min_coord)*family.TotalWidth + bitx - min_coord];

                    if (v > 0) 
                    {
                        white_score += (float)v;
                        white_score_count++;
                        rcode |= 1;
                    } 
                    else 
                    {
                        black_score -= (float)v;
                        black_score_count++;
                    }
                }

                // Debug.Log("white score " + white_score);
                // Debug.Log("white score count " + white_score_count);
                // Debug.Log("black score " + black_score);
                // Debug.Log("white score count " + black_score_count);
                family.QuickDecodeCodeword(rcode, out entry);
                return Mathf.Min(white_score / white_score_count, black_score / black_score_count);
            }
        }
    }

    public class QuickDecodeEntry
    {
        public ulong RCode;
        public ushort ID;
        public byte Hamming;
        public byte Rotation;

        public QuickDecodeEntry()
        {
            RCode = ulong.MaxValue;
        }
    }
}
