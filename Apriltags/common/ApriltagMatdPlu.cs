using System.Collections;
using System.Collections.Generic;
using System;
using UnityEngine;

namespace Apriltags
{
    public class MatdPlu
    {
        public int Singular;
        public uint[] Piv;
        public int PivSign;
        public Matd Lu; 

        public MatdPlu(Matd a)
        {
            uint[] piv = new uint[a.Rows];
            int pivsign = 1;
            Matd lu = new Matd(a);

            // matd_plu_t *mlu = calloc(1, sizeof(matd_plu_t));

            for (int i = 0; i < a.Rows; i++)
            {
                piv[i] = (uint)i;
            }

            for (int j = 0; j < a.Columns; j++) 
            {
                for (int i = 0; i < a.Rows; i++) 
                {
                    int kmax = i < j ? i : j; // min(i,j)

                    // compute dot product of row i with column j (up through element kmax)
                    double acc = 0;
                    for (int k = 0; k < kmax; k++)
                    {
                        acc += lu.GetCell(i,k) * lu.GetCell(k,j);
                    }

                    lu.SetCell(i,j,lu.GetCell(i,j)-acc);
                }

                // find pivot and exchange if necessary.
                int p = j;
                for (int i = j+1; i < lu.Rows; i++) 
                {
                    if (Math.Abs(lu.GetCell(i,j)) > Math.Abs(lu.GetCell(p,j))) 
                    {
                        p = i;
                    }
                }

                // swap rows p and j?
                if (p != j) 
                {
                    for (int i = 0; i < lu.Columns; i++)
                    {
                        double tmp = lu.GetCell(p,i);
                        lu.SetCell(p,i,lu.GetCell(j,i));
                        lu.SetCell(j,i,tmp);
                    }
                    int k = (int)piv[p];
                    piv[p] = piv[j];
                    piv[j] = (uint)k;
                    pivsign = -pivsign;
                }

                double LUjj = lu.GetCell(j,j);

                // If our pivot is very small (which means the matrix is
                // singular or nearly singular), replace with a new pivot of the
                // right sign.
                if (Math.Abs(LUjj) < Utils.Calculations.MATD_EPS) 
                {
        /*
                    if (LUjj < 0)
                        LUjj = -MATD_EPS;
                    else
                        LUjj = MATD_EPS;

                    MATD_EL(lu, j, j) = LUjj;
        */
                    Singular= 1;
                }

                if (j < lu.Columns && j < lu.Rows && LUjj != 0) {
                    LUjj = 1.0 / LUjj;
                    for (int i = j+1; i < lu.Rows; i++)
                    {
                        lu.SetCell(i,j,lu.GetCell(i,j)*LUjj);
                    }
                }
            }

            Lu = lu;
            Piv = piv;
            PivSign = pivsign;
        }

        public Matd Solve(Matd b)
        {
            Matd x = new Matd(b);

            // permute right hand side
            for (int i = 0; i < Lu.Rows; i++)
            {
                for (int j = 0; j < b.Columns; j++)
                {
                    x.SetCell(i,j,b.GetCell((int)Piv[i],j));
                }
            }

            // solve Ly = b
            for (int k = 0; k < Lu.Rows; k++) 
            {
                for (int i = k+1; i < Lu.Rows; i++) 
                {
                    double LUik = -Lu.GetCell(i,k);
                    for (int t = 0; t < b.Columns; t++)
                    {
                        x.SetCell(i,((byte)t),x.GetCell(i,t) + x.GetCell(k,t)* LUik);
                    }
                }
            }

            // solve Ux = y
            for (int k = (int)Lu.Columns-1; k >= 0; k--) 
            {
                double LUkk = 1.0 / Lu.GetCell(k,k);
                for (int t = 0; t < b.Columns; t++)
                {
                    x.SetCell(k,t,x.GetCell(k,t)*LUkk);
                }

                for (int i = 0; i < k; i++) 
                {
                    double LUik = -Lu.GetCell(i,k);
                    for (int t = 0; t < b.Columns; t++)
                    {
                        x.SetCell(i,t,x.GetCell(i,t) + x.GetCell(k,t)*LUik);
                    }
                }
            }

            return x;
        }
    }
}
