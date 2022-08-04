using System.Collections;
using System.Collections.Generic;
using System;

namespace Apriltags
{
    public class TagPose
    {
        public Matd R;
        public Matd T;
        public double Err;

        private TagPose()
        {

        }

        public TagPose(Detection detection, CameraParams cameraParams)
        {
            TagPose pose1 = new TagPose(), pose2 = new TagPose();
            estimateTagPoseOrthogonalIteration(detection, cameraParams, pose1, pose2, 50);
            if (pose1.Err <= pose2.Err) 
            {
                R = pose1.R;
                T = pose1.T;
                Err = pose1.Err;
            } 
            else 
            {
                R = pose2.R;
                T = pose2.T;
                Err = pose2.Err;
            }
        }

        void estimateTagPoseOrthogonalIteration(Detection det, CameraParams par, TagPose solution1, 
            TagPose solution2, int nIters) 
        {
            double scale = par.TagSize/2.0;
            Matd[] p = new Matd[]
            {
                new Matd(3, 1, new double[] {-scale, scale, 0}),
                new Matd(3, 1, new double[] {scale, scale, 0}),
                new Matd(3, 1, new double[] {scale, -scale, 0}),
                new Matd(3, 1, new double[] {-scale, -scale, 0})
            };

            Matd[] v = new Matd[4];
            for (int i = 0; i < 4; i++) 
            {
                v[i] = new Matd(3, 1, new double[] 
                {
                    (det.Corners[i][0] - par.Cx)/par.Fx, (det.Corners[i][1] - par.Cy)/par.Fy, 1
                });
            }

            estimatePoseForTagHomography(det, par, solution1);
            solution1.OrthogonalIteration(v, p, 4, nIters);
            solution2.R = fixPoseAmbiguities(v, p, solution1.T, solution1.R, 4);
            if (solution2.R != null) 
            {
                solution2.T = new Matd(3, 1);
                solution2.OrthogonalIteration(v, p, 4, nIters);
            } 
            else 
            {
                solution2.Err = double.MaxValue;
            }
        }

        void estimatePoseForTagHomography(Detection det, CameraParams par, TagPose solution) 
        {
            double scale = par.TagSize/2.0;

            Matd M_H = Utils.Homography.HomographyToPose(det.H, -par.Fx, par.Fy, par.Cx, par.Cy);
            M_H.SetCell(0,3,M_H.GetCell(0,3)*scale);
            M_H.SetCell(1,3,M_H.GetCell(0,3)*scale);
            M_H.SetCell(2,3,M_H.GetCell(0,3)*scale);

            Matd fix = new Matd(4,4);
            fix.SetCell(0,0,1);
            fix.SetCell(1,1,-1);
            fix.SetCell(2,2,-1);
            fix.SetCell(3,3,1);

            Matd initial_pose = Matd.MatdMultiply(fix, M_H);

            solution.R = new Matd(3, 3);
            for (int i = 0; i < 3; i++) 
            {
                for (int j = 0; j < 3; j++) 
                {
                    solution.R.SetCell(i,j,initial_pose.GetCell(i,j));
                }
            }

            solution.T = new Matd(3, 1);
            for (int i = 0; i < 3; i++) 
            {
                solution.T.SetCell(i,0,initial_pose.GetCell(i,3));
            }
        }

        private void OrthogonalIteration(Matd[] v, Matd[] p, int n_points, int n_steps) 
        {
            Matd p_mean = new Matd(3,1);
            for (int i = 0; i < n_points; i++) 
            {
                p_mean.AddInPlace(p[i]);
            }
            p_mean.ScaleInPlace(1.0/n_points);

            Matd[] p_res = new Matd[n_points];
            for (int i = 0; i < n_points; i++) 
            {
                p_res[i] = Matd.MatdOp("M-M", new List<Matd>() { p[i], p_mean });
            }

            // Compute M1_inv.
            Matd[] F = new Matd[n_points];
            Matd avg_F = new Matd(3,3);
            for (int i = 0; i < n_points; i++) 
            {
                F[i] = calculateF(v[i]);
                avg_F.AddInPlace(F[i]);
            }
            avg_F.ScaleInPlace(1.0/n_points);
            Matd I3 = Matd.GetIdentityMatrix(3);
            Matd M1 = Matd.MatdSubtract(I3, avg_F);
            Matd M1_inv = M1.GetInverseMatrix();

            double prev_error = double.MaxValue;
            // Iterate.
            for (int i = 0; i < n_steps; i++) 
            {
                // Calculate translation.
                Matd M2 = new Matd(3,1);
                for (int j = 0; j < n_points; j++) 
                {
                    Matd M2_update = Matd.MatdOp("(M - M)*M*M", new List<Matd>() { F[j], I3, R, p[j] });
                    M2.AddInPlace(M2_update);
                }
                M2.ScaleInPlace(1.0/n_points);
                T = Matd.MatdMultiply(M1_inv, M2);

                // Calculate rotation.
                Matd[] q = new Matd[n_points];
                Matd q_mean = new Matd(3,1);
                for (int j = 0; j < n_points; j++) 
                {
                    q[j] = Matd.MatdOp("M*(M*M+M)", new List<Matd>() { F[j], R, p[j], T });
                    q_mean.AddInPlace(q[j]);
                }
                q_mean.ScaleInPlace(1.0/n_points);

                Matd M3 = new Matd(3,3);
                for (int j = 0; j < n_points; j++) 
                {
                    Matd M3_update = Matd.MatdOp("(M-M)*M'", new List<Matd>() { q[j], q_mean, p_res[j] });
                    M3.AddInPlace(M3_update);
                }
                MatdSVD M3_svd = new MatdSVD(M3);

                R = Matd.MatdOp("M*M'", new List<Matd>() { M3_svd.U, M3_svd.V });

                double error = 0;
                for (int j = 0; j < 4; j++) 
                {
                    Matd err_vec = Matd.MatdOp("(M-M)(MM+M)", new List<Matd>() { I3, F[j], R, p[j], T});
                    error += (Matd.MatdOp("M'M", new List<Matd>() { err_vec, err_vec })).ToDouble();
                }
                prev_error = error;
            }

            Err = prev_error;
        }

        private static Matd calculateF(Matd v) 
        {
            Matd outer_product = Matd.MatdOp("MM'", new List<Matd>() { v, v, v, v });
            Matd inner_product = Matd.MatdOp("M'M", new List<Matd>() { v, v });
            outer_product.ScaleInPlace(1.0/inner_product.Data[0]);

            return outer_product;
        }

        private static Matd fixPoseAmbiguities(Matd[] v, Matd[] p, Matd t, Matd R, int n_points) 
        {
            Matd I3 = Matd.GetIdentityMatrix(3);

            // 1. Find R_t
            Matd R_t_3 = t.GetVecNormalizeMatrix();

            Matd e_x = new Matd(3,1);
            e_x.SetCell(0,0,1);
            Matd R_t_1_tmp = Matd.MatdOp("M-(M'*M)*M", new List<Matd>() { e_x, e_x, R_t_3, R_t_3 });
            Matd R_t_1 = R_t_1_tmp.GetVecNormalizeMatrix();

            Matd R_t_2 = Matd.MatdCrossProduct(R_t_3, R_t_1);

            Matd R_t = new Matd(3, 3, new double[] 
            {
                    R_t_1.GetCell(0,0), R_t_1.GetCell(0,1), R_t_1.GetCell(0,2),
                    R_t_2.GetCell(0,0), R_t_2.GetCell(0,1), R_t_2.GetCell(0,2),
                    R_t_3.GetCell(0,0), R_t_3.GetCell(0,1), R_t_3.GetCell(0,2)
            });

            // 2. Find R_z
            Matd R_1_prime = Matd.MatdMultiply(R_t, R);
            double r31 = R_1_prime.GetCell(2,0);
            double r32 = R_1_prime.GetCell(2,1);
            double hypotenuse = Math.Sqrt(r31*r31 + r32*r32);
            if (hypotenuse < 1e-100) 
            {
                r31 = 1;
                r32 = 0;
                hypotenuse = 1;
            }
            Matd R_z = new Matd(3, 3, new double[] 
            {
                    r31/hypotenuse, -r32/hypotenuse, 0,
                    r32/hypotenuse, r31/hypotenuse, 0,
                    0, 0, 1
            });

            // 3. Calculate parameters of Eos
            Matd R_trans = Matd.MatdMultiply(R_1_prime, R_z);
            double sin_gamma = -R_trans.GetCell(0,1);
            double cos_gamma = R_trans.GetCell(1,1);
            Matd R_gamma = new Matd(3, 3, new double[] 
            {
                    cos_gamma, -sin_gamma, 0,
                    sin_gamma, cos_gamma, 0,
                    0, 0, 1
            });

            double sin_beta = -R_trans.GetCell(2,0);
            double cos_beta = R_trans.GetCell(2,2);
            double t_initial = Math.Atan2(sin_beta, cos_beta);

            Matd[] v_trans = new Matd[n_points];
            Matd[] p_trans = new Matd[n_points];
            Matd[] F_trans = new Matd[n_points];
            Matd avg_F_trans = new Matd(3,3);
            for (int i = 0; i < n_points; i++) 
            {
                p_trans[i] = Matd.MatdOp("M'*M", new List<Matd>() { R_z, p[i] });
                v_trans[i] = Matd.MatdOp("M*M", new List<Matd>() { R_t, v[i] });
                F_trans[i] = calculateF(v_trans[i]);
                avg_F_trans.AddInPlace(F_trans[i]);
            }
            avg_F_trans.ScaleInPlace(1.0/n_points);

            Matd G = Matd.MatdOp("(M-M)^-1", new List<Matd>() { I3, avg_F_trans });
            G.ScaleInPlace(1.0/n_points);

            Matd M1 = new Matd(3, 3, new double[] 
            {
                    0, 0, 2,
                    0, 0, 0,
                    -2, 0, 0
            });
            Matd M2 = new Matd(3, 3, new double[] 
            {
                    -1, 0, 0,
                    0, 1, 0,
                    0, 0, -1
            });

            Matd b0 = new Matd(3,1);
            Matd b1 = new Matd(3,1);
            Matd b2 = new Matd(3,1);
            for (int i = 0; i < n_points; i++) 
            {
                Matd op_tmp1 = Matd.MatdOp("(M-M)MM", new List<Matd>() { F_trans[i], I3, R_gamma, p_trans[i]});
                Matd op_tmp2 = Matd.MatdOp("(M-M)MMM", new List<Matd>() { F_trans[i], I3, R_gamma, M1, p_trans[i] });
                Matd op_tmp3 = Matd.MatdOp("(M-M)MMM", new List<Matd>() { F_trans[i], I3, R_gamma, M2, p_trans[i] });

                b0.AddInPlace(op_tmp1);
                b1.AddInPlace(op_tmp2);
                b2.AddInPlace(op_tmp3);
            }
            Matd b0_ = Matd.MatdMultiply(G, b0);
            Matd b1_ = Matd.MatdMultiply(G, b1);
            Matd b2_ = Matd.MatdMultiply(G, b2);

            double a0 = 0;
            double a1 = 0;
            double a2 = 0;
            double a3 = 0;
            double a4 = 0;
            for (int i = 0; i < n_points; i++) 
            {
                Matd c0 = Matd.MatdOp("(M-M)(MM+M)", new List<Matd>() { I3, F_trans[i], R_gamma, p_trans[i], b0_ });
                Matd c1 = Matd.MatdOp("(M-M)(MMM+M)", new List<Matd>() { I3, F_trans[i], R_gamma, M1, p_trans[i], b1_ });
                Matd c2 = Matd.MatdOp("(M-M)(MMM+M)", new List<Matd>() { I3, F_trans[i], R_gamma, M2, p_trans[i], b2_ });

                a0 += (Matd.MatdOp("M'M", new List<Matd>() { c0, c0 })).ToDouble();
                a1 += (Matd.MatdOp("2M'M", new List<Matd>() { c0, c1 })).ToDouble();
                a2 += (Matd.MatdOp("M'M+2M'M", new List<Matd>() { c1, c1, c0, c2 })).ToDouble();
                a3 += (Matd.MatdOp("2M'M", new List<Matd>() { c1, c2 })).ToDouble();
                a4 += (Matd.MatdOp("M'M", new List<Matd>() { c2, c2 })).ToDouble();
            }

            // 4. Solve for minima of Eos.
            double p0 = a1;
            double p1 = 2*a2 - 4*a0;
            double p2 = 3*a3 - 3*a1;
            double p3 = 4*a4 - 2*a2;
            double p4 = -a3;

            double[] roots = new double[4];
            int n_roots;
            solvePolyApprox(new double [] {p0, p1, p2, p3, p4}, 4, roots, out n_roots);

            double[] minima = new double[4];
            int n_minima = 0;
            for (int i = 0; i < n_roots; i++) 
            {
                double t1 = roots[i];
                double t2 = t1*t1;
                double t3 = t1*t2;
                double t4 = t1*t3;
                double t5 = t1*t4;
                // Check extrema is a minima.
                if (a2 - 2*a0 + (3*a3 - 6*a1)*t1 + (6*a4 - 8*a2 + 10*a0)*t2 + (-8*a3 + 6*a1)*t3 + (-6*a4 + 3*a2)*t4 + a3*t5 >= 0) 
                {
                    // And that it corresponds to an angle different than the known minimum.
                    double tOther = 2*Math.Atan(roots[i]);
                    // We only care about finding a second local minima which is qualitatively
                    // different than the first.
                    if (Math.Abs(tOther - t_initial) > 0.1) 
                    {
                        minima[n_minima++] = roots[i];
                    }
                }
            }

            // 5. Get poses for minima.
            Matd ret = null;
            if (n_minima == 1) 
            {
                double tOther = minima[0];
                Matd R_beta = new Matd(M2);
                R_beta.ScaleInPlace(tOther);
                R_beta.AddInPlace(M1);
                R_beta.ScaleInPlace(tOther);
                R_beta.AddInPlace(I3);
                R_beta.ScaleInPlace(1/(1 + tOther*tOther));
                ret = Matd.MatdOp("M'MMM'", new List<Matd>() { R_t, R_gamma, R_beta, R_z });
            }

            return ret;
        }

        private static void solvePolyApprox(double[] p, int degree, double[] roots, out int n_roots) 
        {
            int MAX_ROOT = 1000;
            if (degree == 1) 
            {
                if (Math.Abs(p[0]) > MAX_ROOT*Math.Abs(p[1])) 
                {
                    n_roots = 0;
                } 
                else 
                {
                    roots[0] = -p[0]/p[1];
                    n_roots = 1;
                }
                return;
            }

            // Calculate roots of derivative.
            double[] p_der = new double[degree];
            for (int i = 0; i < degree; i++) 
            {
                p_der[i] = (i + 1) * p[i+1];
            }

            double[] der_roots = new double[degree - 1];
            int n_der_roots;
            solvePolyApprox(p_der, degree - 1, der_roots, out n_der_roots);


            // Go through all possibilities for roots of the polynomial.
            n_roots = 0;
            for (int i = 0; i <= n_der_roots; i++) 
            {
                double min;
                if (i == 0) 
                {
                    min = -MAX_ROOT;
                } 
                else 
                {
                    min = der_roots[i - 1];
                }

                double max;
                if (i == n_der_roots) 
                {
                    max = MAX_ROOT;
                } 
                else 
                {
                    max = der_roots[i];
                }

                if (polyval(p, degree, min)*polyval(p, degree, max) < 0) 
                {
                    // We have a zero-crossing in this interval, use a combination of Newton' and bisection.
                    // Some thanks to Numerical Recipes in C.

                    double lower;
                    double upper;
                    if (polyval(p, degree, min) < polyval(p, degree, max)) 
                    {
                        lower = min;
                        upper = max;
                    } 
                    else 
                    {
                        lower = max;
                        upper = min;
                    }
                    double root = 0.5*(lower + upper);
                    double dx_old = upper - lower;
                    double dx = dx_old;
                    double f = polyval(p, degree, root);
                    double df = polyval(p_der, degree - 1, root);

                    for (int j = 0; j < 100; j++) 
                    {
                        if (((f + df*(upper - root))*(f + df*(lower - root)) > 0)
                                || (Math.Abs(2*f) > Math.Abs(dx_old*df))) 
                        {
                            dx_old = dx;
                            dx = 0.5*(upper - lower);
                            root = lower + dx;
                        } 
                        else 
                        {
                            dx_old = dx;
                            dx = -f/df;
                            root += dx;
                        }

                        if (root == upper || root == lower) {
                            break;
                        }

                        f = polyval(p, degree, root);
                        df = polyval(p_der, degree - 1, root);

                        if (f > 0) 
                        {
                            upper = root;
                        } 
                        else 
                        {
                            lower = root;
                        }
                    }

                    roots[n_roots++] = root;
                } 
                else if(polyval(p, degree, max) == 0) 
                {
                    // Double/triple root.
                    roots[n_roots++] = max;
                }
            }
        }

        private static double polyval(double[] p, int degree, double x) 
        {
            double ret = 0;
            for (int i = 0; i <= degree; i++) 
            {
                ret += p[i]*Math.Pow(x, i);
            }
            return ret;
        }
    }
}
