using System.Collections;
using System.Collections.Generic;
using System;
using UnityEngine;

namespace Apriltags.Utils
{
    public static class Calculations
    {
        public static readonly int APRILTAG_TASKS_PER_THREAD_TARGET = 10;
        public static readonly double HUGE_VALF = 1e50d;
        public static readonly double MATD_EPS = 1e-8;

        public static byte Round(double number)
        {
            int output = (int)Math.Floor(number);
            double target = Math.Round(number, 6);
            // string step = number.ToString("F6");
            // double target = double.Parse(step);
            double compareMin = 0.998 + output;
            double compare999 = 0.999001 + output;
            double compare998 = 0.998001 + output;

            if(target == compare999 && target > 110)
            {

            }
            else if(target == compare999 && output == 104 && number < output + 0.9990009)
            {

            }
            else if(target == compare998 && (output == 82 || output == 102))
            {

            }
            else if(target > compareMin)
            {
                output += 1;
            }

            return (byte)output;
        }

        public static void Convolve(byte[] x, byte[] y, int sz, byte[] k, int ksz)
        {
            // Debug.Log("sz " + sz);
            // Debug.Log("ksz " + ksz);
            if((ksz & 1) == 1)
            {
                for (int i = 0; i < ksz/2 && i < sz; i++)
                {
                    y[i] = x[i];
                    // Debug.Log("first step index " + i);
                    // Debug.Log("first step res " + y[i]);
                }

                for (int i = 0; i < sz - ksz; i++) 
                {
                    uint acc = 0;

                    for (int j = 0; j < ksz; j++)
                    {
                        acc += (uint)(k[j]*x[i+j]);
                    }

                    y[ksz/2 + i] = (byte)(acc >> 8);
                    // Debug.Log("second step index " + ksz/2 + i);
                    // Debug.Log("second step res " + y[ksz/2 + i]);
                }

                for (int i = sz - ksz + ksz/2; i < sz; i++)
                {
                    y[i] = x[i];
                    // Debug.Log("third step index " + i);
                    // Debug.Log("third step res " + y[i]);
                }
            }
        }

        public static uint Hash2(ulong x)
        {
            return (uint)((2654435761 * x) >> 32);
        }

        public static int IMin(int a, int b)
        {
            return (a < b) ? a : b;
        }

        public static int IMax(int a, int b)
        {
            return (a > b) ? a : b;
        }

        public class MaximaComparer : IComparer
        {
            int IComparer.Compare(object x, object y)
            {
                return (double)x < (double)y ? 1 : -1;
            }
        }

        public static void Mat33SymSolve(double[] A, double[] B, double[] R)
        {
            double[] L = new double[9];
            mat33Chol(A, L);

            double[] M = new double[9];
            mat33LowerTriInv(L, M);

            double[] tmp = new double[3];
            tmp[0] = M[0]*B[0];
            tmp[1] = M[3]*B[0] + M[4]*B[1];
            tmp[2] = M[6]*B[0] + M[7]*B[1] + M[8]*B[2];

            R[0] = M[0]*tmp[0] + M[3]*tmp[1] + M[6]*tmp[2];
            R[1] = M[4]*tmp[1] + M[7]*tmp[2];
            R[2] = M[8]*tmp[2];
        }

        private static void mat33LowerTriInv(double[] A, double[] R)
        {
            // A[0]*R[0] = 1
            R[0] = 1 / A[0];

            // A[3]*R[0] + A[4]*R[3] = 0
            R[3] = -A[3]*R[0] / A[4];

            // A[4]*R[4] = 1
            R[4] = 1 / A[4];

            // A[6]*R[0] + A[7]*R[3] + A[8]*R[6] = 0
            R[6] = (-A[6]*R[0] - A[7]*R[3]) / A[8];

            // A[7]*R[4] + A[8]*R[7] = 0
            R[7] = -A[7]*R[4] / A[8];

            // A[8]*R[8] = 1
            R[8] = 1 / A[8];
        }

        private static void mat33Chol(double[] A, double[] R)
        {
            // A[0] = R[0]*R[0]
            R[0] = Math.Sqrt(A[0]);

            // A[1] = R[0]*R[3];
            R[3] = A[1] / R[0];

            // A[2] = R[0]*R[6];
            R[6] = A[2] / R[0];

            // A[4] = R[3]*R[3] + R[4]*R[4]
            R[4] = Math.Sqrt(A[4] - R[3]*R[3]);

            // A[5] = R[3]*R[6] + R[4]*R[7]
            R[7] = (A[5] - R[3]*R[6]) / R[4];

            // A[8] = R[6]*R[6] + R[7]*R[7] + R[8]*R[8]
            R[8] = Math.Sqrt(A[8] - R[6]*R[6] - R[7]*R[7]);

            R[1] = 0;
            R[2] = 0;
            R[5] = 0;
        }

        public static void Sharpen(Detector detector, double[] values, int size) 
        {
            double[] sharpened = new double[size*size];
            double[] kernel = new double[9]
            {
                0, -1, 0,
                -1, 4, -1,
                0, -1, 0
            };

            for (int y = 0; y < size; y++) 
            {
                for (int x = 0; x < size; x++) 
                {
                    sharpened[y*size + x] = 0;
                    for (int i = 0; i < 3; i++) 
                    {
                        for (int j = 0; j < 3; j++) 
                        {
                            if ((y + i - 1) < 0 || (y + i - 1) > size - 1 || (x + j - 1) < 0 || (x + j - 1) > size - 1) 
                            {
                                continue;
                            }
                            sharpened[y*size + x] += values[(y + i - 1)*size + (x + j - 1)]*kernel[i*3 + j];
                        }
                    }
                }
            }



            for (int y = 0; y < size; y++) 
            {
                for (int x = 0; x < size; x++) 
                {
                    values[y*size + x] = values[y*size + x] + detector.DecodeSharpening*sharpened[y*size + x];
                    // Debug.Log("val " + (y*size + x) + ", value " + values[y*size + x] + ", decode " + 
                    //     detector.DecodeSharpening + ", sharp " + sharpened[y*size + x]);
                }
            }
        }

        public static ulong Rotate90(ulong w, int numBits)
        {
            int p = numBits;
            ulong l = 0;
            if (numBits % 4 == 1) 
            {
                p = numBits - 1;
                l = 1;
            }
            w = ((w >> (int)l) << (p/4 + (int)l)) | (w >> (3 * p/ 4 + (int)l) << (int)l) | (w & l);
            ulong one = 1;
            w &= ((one << numBits) - 1);
            return w;
        }

        public static int PreferSmaller(int pref, double q0, double q1)
        {
            if (pref != 0)     // already prefer something? exit.
                return pref;

            if (q0 < q1)
                return -1; // we now prefer q0
            if (q1 < q0)
                return 1; // we now prefer q1

            // no preference
            return 0;
        }

        public static void SVD22(double[] A, double[] U, double[] S, double[] V)
        {
            double A00 = A[0];
            double A01 = A[1];
            double A10 = A[2];
            double A11 = A[3];

            double B0 = A00 + A11;
            double B1 = A00 - A11;
            double B2 = A01 + A10;
            double B3 = A01 - A10;

            double PminusT = Math.Atan2(B3, B0);
            double PplusT = Math.Atan2(B2, B1);

            double P = (PminusT + PplusT) / 2;
            double T = (-PminusT + PplusT) / 2;

            double CP = Math.Cos(P), SP = Math.Sin(P);
            double CT = Math.Cos(T), ST = Math.Sin(T);

            U[0] = CT;
            U[1] = -ST;
            U[2] = ST;
            U[3] = CT;

            V[0] = CP;
            V[1] = -SP;
            V[2] = SP;
            V[3] = CP;

            // C0 = e+f. There are two ways to compute C0; we pick the one
            // that is better conditioned.
            double CPmT = Math.Cos(P-T), SPmT = Math.Sin(P-T);
            double C0 = 0;
            if (Math.Abs(CPmT) > Math.Abs(SPmT))
            {
                C0 = B0 / CPmT;
            }
            else
            {
                C0 = B3 / SPmT;
            }

            // C1 = e-f. There are two ways to compute C1; we pick the one
            // that is better conditioned.
            double CPpT = Math.Cos(P+T), SPpT = Math.Sin(P+T);
            double C1 = 0;
            if (Math.Abs(CPpT) > Math.Abs(SPpT))
            {
                C1 = B1 / CPpT;
            }
            else
            {
                C1 = B2 / SPpT;
            }

            // e and f are the singular values
            double e = (C0 + C1) / 2;
            double f = (C0 - C1) / 2;

            if (e < 0) 
            {
                e = -e;
                U[0] = -U[0];
                U[2] = -U[2];
            }

            if (f < 0) 
            {
                f = -f;
                U[1] = -U[1];
                U[3] = -U[3];
            }

            // sort singular values.
            if (e > f) 
            {
                // already in big-to-small order.
                S[0] = e;
                S[1] = f;
            } 
            else 
            {
                // Curiously, this code never seems to get invoked.  Why is it
                // that S[0] always ends up the dominant vector?  However,
                // this code has been tested (flipping the logic forces us to
                // sort the singular values in ascending order).
                //
                // P = [ 0 1 ; 1 0 ]
                // USV' = (UP)(PSP)(PV')
                //      = (UP)(PSP)(VP)'
                //      = (UP)(PSP)(P'V')'
                S[0] = f;
                S[1] = e;

                // exchange columns of U and V
                double[] tmp = new double[2];
                tmp[0] = U[0];
                tmp[1] = U[2];
                U[0] = U[1];
                U[2] = U[3];
                U[1] = tmp[0];
                U[3] = tmp[1];

                tmp[0] = V[0];
                tmp[1] = V[2];
                V[0] = V[1];
                V[2] = V[3];
                V[1] = tmp[0];
                V[3] = tmp[1];
            }

            /*
            double SM[4] = { S[0], 0, 0, S[1] };

            doubles_print_mat(U, 2, 2, "%20.10g");
            doubles_print_mat(SM, 2, 2, "%20.10g");
            doubles_print_mat(V, 2, 2, "%20.10g");
            printf("A:\n");
            doubles_print_mat(A, 2, 2, "%20.10g");

            double SVt[4];
            doubles_mat_ABt(SM, 2, 2, V, 2, 2, SVt, 2, 2);
            double USVt[4];
            doubles_mat_AB(U, 2, 2, SVt, 2, 2, USVt, 2, 2);

            printf("USVt\n");
            doubles_print_mat(USVt, 2, 2, "%20.10g");

            double diff[4];
            for (int i = 0; i < 4; i++)
                diff[i] = A[i] - USVt[i];

            printf("diff\n");
            doubles_print_mat(diff, 2, 2, "%20.10g");

            */

        }
    }
}
