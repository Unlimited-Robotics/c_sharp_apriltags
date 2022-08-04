using System.Collections;
using System.Collections.Generic;
using System;
using UnityEngine;

namespace Apriltags.Utils
{
    public static class Homography
    {
        public static Matd HomographyCompute2(double[][] c) 
        {
            double[] A = new double[]  
            {
                c[0][0], c[0][1], 1,       0,       0, 0, -c[0][0]*c[0][2], -c[0][1]*c[0][2], c[0][2],
                    0,       0, 0, c[0][0], c[0][1], 1, -c[0][0]*c[0][3], -c[0][1]*c[0][3], c[0][3],
                c[1][0], c[1][1], 1,       0,       0, 0, -c[1][0]*c[1][2], -c[1][1]*c[1][2], c[1][2],
                    0,       0, 0, c[1][0], c[1][1], 1, -c[1][0]*c[1][3], -c[1][1]*c[1][3], c[1][3],
                c[2][0], c[2][1], 1,       0,       0, 0, -c[2][0]*c[2][2], -c[2][1]*c[2][2], c[2][2],
                    0,       0, 0, c[2][0], c[2][1], 1, -c[2][0]*c[2][3], -c[2][1]*c[2][3], c[2][3],
                c[3][0], c[3][1], 1,       0,       0, 0, -c[3][0]*c[3][2], -c[3][1]*c[3][2], c[3][2],
                    0,       0, 0, c[3][0], c[3][1], 1, -c[3][0]*c[3][3], -c[3][1]*c[3][3], c[3][3],
            };

            // double epsilon = 1e-10;

            // Eliminate.
            for (int col = 0; col < 8; col++) {
                // Find best row to swap with.
                double max_val = 0;
                int max_val_idx = -1;
                for (int row = col; row < 8; row++) {
                    double val = Math.Abs(A[row*9 + col]);
                    if (val > max_val) {
                        max_val = val;
                        max_val_idx = row;
                    }
                }

                // if (max_val < epsilon) {
                //     fprintf(stderr, "WRN: Matrix is singular.\n");
                // }

                // Swap to get best row.
                if (max_val_idx != col) {
                    for (int i = col; i < 9; i++) {
                        double tmp = A[col*9 + i];
                        A[col*9 + i] = A[max_val_idx*9 + i];
                        A[max_val_idx*9 + i] = tmp;
                    }
                }

                // Do eliminate.
                for (int i = col + 1; i < 8; i++) {
                    double f = A[i*9 + col]/A[col*9 + col];
                    A[i*9 + col] = 0;
                    for (int j = col + 1; j < 9; j++) {
                        A[i*9 + j] -= f*A[col*9 + j];
                    }
                }
            }

            // Back solve.
            for (int col = 7; col >=0; col--) {
                double sum = 0;
                for (int i = col + 1; i < 8; i++) {
                    sum += A[col*9 + i]*A[i*9 + 8];
                }
                A[col*9 + 8] = (A[col*9 + 8] - sum)/A[col*9 + col];
            }

            return new  Matd(3, 3, new double[] { A[8], A[17], A[26], A[35], A[44], A[53], A[62], A[71], 1 });
        }

        public static void HomographyProject(Matd H, double x, double y, out double ox, out double oy)
        {
            double xx = H.GetCell(0,0)*x + H.GetCell(0,1)*y + H.GetCell(0,2);
            double yy = H.GetCell(1,0)*x + H.GetCell(1,1)*y + H.GetCell(1,2);
            double zz = H.GetCell(2,0)*x + H.GetCell(2,1)*y + H.GetCell(2,2);

            ox = xx / zz;
            oy = yy / zz;
        }

        public static Matd HomographyToPose(Matd H, double fx, double fy, double cx, double cy)
        {
            // Note that every variable that we compute is proportional to the scale factor of H.
            double R20 = H.GetCell(2,0);
            double R21 = H.GetCell(2,1);
            double TZ  = H.GetCell(2,2);
            double R00 = (H.GetCell(0,0) - cx*R20) / fx;
            double R01 = (H.GetCell(0,1) - cx*R21) / fx;
            double TX  = (H.GetCell(0,2) - cx*TZ)  / fx;
            double R10 = (H.GetCell(1,0) - cy*R20) / fy;
            double R11 = (H.GetCell(1,1) - cy*R21) / fy;
            double TY  = (H.GetCell(1,2) - cy*TZ)  / fy;

            // compute the scale by requiring that the rotation columns are unit length
            // (Use geometric average of the two length vectors we have)
            double length1 = Math.Sqrt(R00*R00 + R10*R10 + R20*R20);
            double length2 = Math.Sqrt(R01*R01 + R11*R11 + R21*R21);
            double s = 1.0 / Math.Sqrt(length1 * length2);

            // get sign of S by requiring the tag to be in front the camera;
            // we assume camera looks in the -Z direction.
            if (TZ > 0)
            {
                s *= -1;
            }

            R20 *= s;
            R21 *= s;
            TZ  *= s;
            R00 *= s;
            R01 *= s;
            TX  *= s;
            R10 *= s;
            R11 *= s;
            TY  *= s;

            // now recover [R02 R12 R22] by noting that it is the cross product of the other two columns.
            double R02 = R10*R21 - R20*R11;
            double R12 = R20*R01 - R00*R21;
            double R22 = R00*R11 - R10*R01;

            Matd R = new Matd(3, 3, new double[] { R00, R01, R02, R10, R11, R12, R20, R21, R22 });

            MatdSVD svd = new MatdSVD(R);

            R = Matd.MatdOp("M*M'", new List<Matd>() { svd.U, svd.V });

            R00 = R.GetCell(0,0);
            R01 = R.GetCell(0,1);
            R02 = R.GetCell(0,2);
            R10 = R.GetCell(1,0);
            R11 = R.GetCell(1,1);
            R12 = R.GetCell(1,2);
            R20 = R.GetCell(2,0);
            R21 = R.GetCell(2,1);
            R22 = R.GetCell(2,2);

            return new Matd(4, 4, new double[] 
            { 
                R00, R01, R02, TX,
                R10, R11, R12, TY,
                R20, R21, R22, TZ,
                0, 0, 0, 1 
            });
        }
    }
}
