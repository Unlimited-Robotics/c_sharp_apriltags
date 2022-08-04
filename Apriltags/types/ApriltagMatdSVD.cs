using System.Collections;
using System.Collections.Generic;
using System;

namespace Apriltags
{
    public class MatdSVD
    {
        public Matd U;
        public Matd S;
        public Matd V;

        public MatdSVD(Matd A)
        {
            flags(A, 0);
        }

        private MatdSVD()
        {

        }

        private void flags(Matd A, int flags)
        {
            if (A.Columns <= A.Rows) 
            {
                MatdSVD res = matdSvdTall(A, flags);
                U = res.U;
                S = res.S;
                V = res.V;
            } 
            else 
            {
                Matd At = A.GetTransposeMatrix();

                // A =U  S  V'
                // A'=V  S' U'

                MatdSVD tmp = matdSvdTall(At, flags);


                U = tmp.V; //matd_transpose(tmp.V);
                S = tmp.S.GetTransposeMatrix();
                V = tmp.U; //matd_transpose(tmp.U);
            }

        /*
        matd_t *check = matd_op("M*M*M'-M", res.U, res.S, res.V, A);
        double maxerr = 0;

        for (int i = 0; i < check->nrows; i++)
        for (int j = 0; j < check->ncols; j++)
        maxerr = fmax(maxerr, fabs(MATD_EL(check, i, j)));

        matd_destroy(check);

        if (maxerr > 1e-7) {
        printf("bad maxerr: %15f\n", maxerr);
        }

        if (maxerr > 1e-5) {
        printf("bad maxerr: %15f\n", maxerr);
        matd_print(A, "%15f");
        assert(0);
        }

        */
        }

        private static MatdSVD matdSvdTall(Matd A, int flags)
        {
            Matd B = new Matd(A);

            // Apply householder reflections on each side to reduce A to
            // bidiagonal form. Specifically:
            //
            // A = LS*B*RS'
            //
            // Where B is bidiagonal, and LS/RS are unitary.
            //
            // Why are we doing this? Some sort of transformation is necessary
            // to reduce the matrix's nz elements to a square region. QR could
            // work too. We need nzs confined to a square region so that the
            // subsequent iterative process, which is based on rotations, can
            // work. (To zero out a term at (i,j), our rotations will also
            // affect (j,i).
            //
            // We prefer bidiagonalization over QR because it gets us "closer"
            // to the SVD, which should mean fewer iterations.

            // LS: cumulative left-handed transformations
            Matd LS = Matd.GetIdentityMatrix((int)A.Rows);

            // RS: cumulative right-handed transformations.
            Matd RS = Matd.GetIdentityMatrix((int)A.Columns);

            for (int hhidx = 0; hhidx < A.Rows; hhidx++)  
            {

                if (hhidx < A.Columns) 
                {
                    // We construct the normal of the reflection plane: let u
                    // be the vector to reflect, x =[ M 0 0 0 ] the target
                    // location for u (u') after reflection (with M = ||u||).
                    //
                    // The normal vector is then n = (u - x), but since we
                    // could equally have the target location be x = [-M 0 0 0
                    // ], we could use n = (u + x).
                    //
                    // We then normalize n. To ensure a reasonable magnitude,
                    // we select the sign of M so as to maximize the magnitude
                    // of the first element of (x +/- M). (Otherwise, we could
                    // end up with a divide-by-zero if u[0] and M cancel.)
                    //
                    // The householder reflection matrix is then H=(I - nn'), and
                    // u' = Hu.
                    //
                    //
                    int vlen = (int)A.Rows - hhidx;

                    double[] v = new double[vlen];

                    double mag2 = 0;
                    for (int i = 0; i < vlen; i++) 
                    {
                        v[i] = B.GetCell(hhidx+i, hhidx);
                        mag2 += v[i]*v[i];
                    }

                    double oldv0 = v[0];
                    if (oldv0 < 0)
                    {
                        v[0] -= Math.Sqrt(mag2);
                    }
                    else
                    {
                        v[0] += Math.Sqrt(mag2);
                    }

                    mag2 += -oldv0*oldv0 + v[0]*v[0];

                    // normalize v
                    double mag = Math.Sqrt(mag2);

                    // this case arises with matrices of all zeros, for example.
                    if (mag == 0)
                    {
                        continue;
                    }

                    for (int i = 0; i < vlen; i++)
                    {
                        v[i] /= mag;
                    }

                    // Q = I - 2vv'
                    //matd_t *Q = matd_identity(A->nrows);
                    //for (int i = 0; i < vlen; i++)
                    //  for (int j = 0; j < vlen; j++)
                    //    MATD_EL(Q, i+hhidx, j+hhidx) -= 2*v[i]*v[j];


                    // LS = matd_op("F*M", LS, Q);
                    // Implementation: take each row of LS, compute dot product with n,
                    // subtract n (scaled by dot product) from it.
                    for (int i = 0; i < LS.Rows; i++) 
                    {
                        double dot = 0;
                        for (int j = 0; j < vlen; j++)
                        {
                            dot += LS.GetCell(i,hhidx+j) * v[j];
                        }
                        for (int j = 0; j < vlen; j++)
                        {
                            LS.SetCell(i,hhidx+j,LS.GetCell(i,hhidx+j) - 2*dot*v[j]);
                        }
                    }

                    //  B = matd_op("M*F", Q, B); // should be Q', but Q is symmetric.
                    for (int i = 0; i < B.Columns; i++) 
                    {
                        double dot = 0;
                        for (int j = 0; j < vlen; j++)
                        {
                            dot += B.GetCell(hhidx+j,i) * v[j];
                        }
                        for (int j = 0; j < vlen; j++)
                        {
                            B.SetCell(hhidx+j,i,B.GetCell(hhidx+j,i) - 2*dot*v[j]);
                        }
                    }
                }

                if (hhidx+2 < A.Columns) 
                {
                    int vlen = (int)A.Columns - hhidx - 1;

                    double[] v = new double[vlen];

                    double mag2 = 0;
                    for (int i = 0; i < vlen; i++) 
                    {
                        v[i] = B.GetCell(hhidx,hhidx+i+1);
                        mag2 += v[i]*v[i];
                    }

                    double oldv0 = v[0];
                    if (oldv0 < 0)
                    {
                        v[0] -= Math.Sqrt(mag2);
                    }
                    else
                    {
                        v[0] += Math.Sqrt(mag2);
                    }

                    mag2 += -oldv0*oldv0 + v[0]*v[0];

                    // compute magnitude of ([1 0 0..]+v)
                    double mag = Math.Sqrt(mag2);

                    // this case can occur when the vectors are already perpendicular
                    if (mag == 0)
                    {
                        continue;
                    }

                    for (int i = 0; i < vlen; i++)
                    {
                        v[i] /= mag;
                    }

                    // TODO: optimize these multiplications
                    // matd_t *Q = matd_identity(A->ncols);
                    //  for (int i = 0; i < vlen; i++)
                    //    for (int j = 0; j < vlen; j++)
                    //       MATD_EL(Q, i+1+hhidx, j+1+hhidx) -= 2*v[i]*v[j];

                    //  RS = matd_op("F*M", RS, Q);
                    for (int i = 0; i < RS.Rows; i++) 
                    {
                        double dot = 0;
                        for (int j = 0; j < vlen; j++)
                        {
                            dot += RS.GetCell(i,hhidx+1+j) * v[j];
                        }
                        for (int j = 0; j < vlen; j++)
                        {
                            RS.SetCell(i,hhidx+1+j,RS.GetCell(i,hhidx+1+j) - 2*dot*v[j]);
                        }
                    }

                    //   B = matd_op("F*M", B, Q); // should be Q', but Q is symmetric.
                    for (int i = 0; i < B.Rows; i++) 
                    {
                        double dot = 0;
                        for (int j = 0; j < vlen; j++)
                        {
                            dot += B.GetCell(i,hhidx+1+j) * v[j];
                        }
                        for (int j = 0; j < vlen; j++)
                        {
                            B.SetCell(i,hhidx+1+j,B.GetCell(i,hhidx+1+j) - 2*dot*v[j]);
                        }
                    }
                }
            }

            // maxiters used to be smaller to prevent us from looping forever,
            // but this doesn't seem to happen any more with our more stable
            // svd22 implementation.
            ulong one = 1;
            int maxiters = (int)(one << 30);
            int iter;

            double maxv; // maximum non-zero value being reduced this iteration

            double tol = 1E-10;

            // which method will we use to find the largest off-diagonal
            // element of B?
            const int find_max_method = 1; //(B->ncols < 6) ? 2 : 1;

            // for each of the first B->ncols rows, which index has the
            // maximum absolute value? (used by method 1)
            int[] maxrowidx = new int[B.Columns];
            int lastmaxi, lastmaxj;

            if (find_max_method == 1) 
            {
                for (int i = 2; i < B.Columns; i++)
                {
                    maxrowidx[i] = B.MaxIndex(i,(int)B.Columns);
                }

                // note that we started the array at 2. That's because by setting
                // these values below, we'll recompute first two entries on the
                // first iteration!
                lastmaxi = 0;
                lastmaxj = 1;
            }

            for (iter = 0; iter < maxiters; iter++) 
            {

                // No diagonalization required for 0x0 and 1x1 matrices.
                if (B.Columns < 2)
                {
                    break;
                }

                // find the largest off-diagonal element of B, and put its
                // coordinates in maxi, maxj.
                int maxi, maxj;

                if (find_max_method == 1) 
                {
                    // method 1 is the "smarter" method which does at least
                    // 4*ncols work. More work might be needed (up to
                    // ncols*ncols), depending on data. Thus, this might be a
                    // bit slower than the default method for very small
                    // matrices.
                    maxi = -1;
                    maxv = -1;

                    // every iteration, we must deal with the fact that rows
                    // and columns lastmaxi and lastmaxj have been
                    // modified. Update maxrowidx accordingly.

                    // now, EVERY row also had columns lastmaxi and lastmaxj modified.
                    for (int rowi = 0; rowi < B.Columns; rowi++) 
                    {

                        // the magnitude of the largest off-diagonal element
                        // in this row.
                        double thismaxv;

                        // row 'lastmaxi' and 'lastmaxj' have been completely
                        // changed. compute from scratch.
                        if (rowi == lastmaxi || rowi == lastmaxj) 
                        {
                            maxrowidx[rowi] = B.MaxIndex(rowi,(int)B.Columns);
                            thismaxv = Math.Abs(B.GetCell(rowi,maxrowidx[rowi]));
                            goto endrowi;
                        }

                        // our maximum entry was just modified. We don't know
                        // if it went up or down, and so we don't know if it
                        // is still the maximum. We have to update from
                        // scratch.
                        if (maxrowidx[rowi] == lastmaxi || maxrowidx[rowi] == lastmaxj) 
                        {
                            maxrowidx[rowi] = B.MaxIndex(rowi,(int)B.Columns);
                            thismaxv = Math.Abs(B.GetCell(rowi,maxrowidx[rowi]));
                            goto endrowi;
                        }

                        // This row is unchanged, except for columns
                        // 'lastmaxi' and 'lastmaxj', and those columns were
                        // not previously the largest entry...  just check to
                        // see if they are now the maximum entry in their
                        // row. (Remembering to consider off-diagonal entries
                        // only!)
                        thismaxv = Math.Abs(B.GetCell(rowi,maxrowidx[rowi]));

                        // check column lastmaxi. Is it now the maximum?
                        if (lastmaxi != rowi) 
                        {
                            double v = Math.Abs(B.GetCell(rowi,lastmaxi));
                            if (v > thismaxv) 
                            {
                                thismaxv = v;
                                maxrowidx[rowi] = lastmaxi;
                            }
                        }

                        // check column lastmaxj
                        if (lastmaxj != rowi) 
                        {
                            double v = Math.Abs(B.GetCell(rowi,lastmaxj));
                            if (v > thismaxv) {
                                thismaxv = v;
                                maxrowidx[rowi] = lastmaxj;
                            }
                        }

                        // does this row have the largest value we've seen so far?
                    endrowi:
                        if (thismaxv > maxv) 
                        {
                            maxv = thismaxv;
                            maxi = rowi;
                        }
                    }

                    maxj = maxrowidx[maxi];

                    // save these for the next iteration.
                    lastmaxi = maxi;
                    lastmaxj = maxj;

                    if (maxv < tol)
                    {
                        break;
                    }
                } 
                else if (find_max_method == 2) 
                {
                    // brute-force (reference) version.
                    maxv = -1;

                    // only search top "square" portion
                    for (int i = 0; i < B.Columns; i++) 
                    {
                        for (int j = 0; j < B.Columns; j++) 
                        {
                            if (i == j)
                            {
                                continue;
                            }

                            double v = Math.Abs(B.GetCell(i,j));

                            if (v > maxv) 
                            {
                                maxi = i;
                                maxj = j;
                                maxv = v;
                            }
                        }
                    }

                    // termination condition.
                    if (maxv < tol)
                    {
                        break;
                    }
                }

        //        printf(">>> %5d %3d, %3d %15g\n", maxi, maxj, iter, maxv);

                // Now, solve the 2x2 SVD problem for the matrix
                // [ A0 A1 ]
                // [ A2 A3 ]
                double A0 = B.GetCell(maxi,maxi);
                double A1 = B.GetCell(maxi,maxj);
                double A2 = B.GetCell(maxj,maxi);
                double A3 = B.GetCell(maxj,maxj);

                double[] AQ = new double[4];
                AQ[0] = A0;
                AQ[1] = A1;
                AQ[2] = A2;
                AQ[3] = A3;

                double[] U = new double[4], S = new double [2], V = new double[4];
                Utils.Calculations.SVD22(AQ, U, S, V);

        /*  Reference (slow) implementation...

                    // LS = LS * ROT(theta) = LS * QL
                    matd_t *QL = matd_identity(A->nrows);
                    MATD_EL(QL, maxi, maxi) = U[0];
                    MATD_EL(QL, maxi, maxj) = U[1];
                    MATD_EL(QL, maxj, maxi) = U[2];
                    MATD_EL(QL, maxj, maxj) = U[3];

                    matd_t *QR = matd_identity(A->ncols);
                    MATD_EL(QR, maxi, maxi) = V[0];
                    MATD_EL(QR, maxi, maxj) = V[1];
                    MATD_EL(QR, maxj, maxi) = V[2];
                    MATD_EL(QR, maxj, maxj) = V[3];

                    LS = matd_op("F*M", LS, QL);
                    RS = matd_op("F*M", RS, QR); // remember we'll transpose RS.
                    B = matd_op("M'*F*M", QL, B, QR);

                    matd_destroy(QL);
                    matd_destroy(QR);
        */

                    //  LS = matd_op("F*M", LS, QL);
                for (int i = 0; i < LS.Rows; i++) 
                {
                    double vi = LS.GetCell(i,maxi);
                    double vj = LS.GetCell(i,maxj);

                    LS.SetCell(i,maxi,U[0]*vi + U[2]*vj);
                    LS.SetCell(i,maxj,U[1]*vi + U[3]*vj);
                }

                    //  RS = matd_op("F*M", RS, QR); // remember we'll transpose RS.
                for (int i = 0; i < RS.Rows; i++) 
                {
                    double vi = RS.GetCell(i,maxi);
                    double vj = RS.GetCell(i,maxj);

                    RS.SetCell(i,maxi,V[0]*vi + V[2]*vj);
                    RS.SetCell(i,maxj,V[1]*vi + V[3]*vj);
                }

                    // B = matd_op("M'*F*M", QL, B, QR);
                    // The QL matrix mixes rows of B.
                for (int i = 0; i < B.Columns; i++) 
                {
                    double vi = B.GetCell(maxi,i);
                    double vj = B.GetCell(maxj,i);

                    B.SetCell(maxi,i,U[0]*vi + U[2]*vj);
                    B.SetCell(maxj,i,U[1]*vi + U[3]*vj);
                }

                    // The QR matrix mixes columns of B.
                for (int i = 0; i < B.Rows; i++) 
                {
                    double vi = B.GetCell(i,maxi);
                    double vj = B.GetCell(i,maxj);

                    B.SetCell(i,maxi,V[0]*vi + V[2]*vj);
                    B.SetCell(i,maxj,V[1]*vi + V[3]*vj);
                }
            }

            // them all positive by flipping the corresponding columns of
            // U/LS.
            int[] idxs = new int[A.Columns];
            double[] vals = new double[A.Columns];
            for (int i = 0; i < A.Columns; i++) 
            {
                idxs[i] = i;
                vals[i] = B.GetCell(i,i);
            }

            // A bubble sort. Seriously.
            int changed;
            do {
                changed = 0;

                for (int i = 0; i + 1 < A.Columns; i++) 
                {
                    if (Math.Abs(vals[i+1]) > Math.Abs(vals[i])) 
                    {
                        int tmpi = idxs[i];
                        idxs[i] = idxs[i+1];
                        idxs[i+1] = tmpi;

                        double tmpv = vals[i];
                        vals[i] = vals[i+1];
                        vals[i+1] = tmpv;

                        changed = 1;
                    }
                }
            } while (changed != 0);

            Matd LP = Matd.GetIdentityMatrix((int)A.Rows);
            Matd RP = Matd.GetIdentityMatrix((int)A.Columns);

            for (int i = 0; i < A.Columns; i++) 
            {
                LP.SetCell(idxs[i],idxs[i],0);
                RP.SetCell(idxs[i],idxs[i],0);

                LP.SetCell(idxs[i],i,vals[i] < 0 ? -1 : 1);
                RP.SetCell(idxs[i],i,1);
            }

            // we've factored:
            // LP*(something)*RP'

            // solve for (something)
            B = Matd.MatdOp("M'*F*M", new List<Matd>() { LP, B, RP });

            // update LS and RS, remembering that RS will be transposed.
            LS = Matd.MatdOp("F*M", new List<Matd>() { LS, LP });
            RS = Matd.MatdOp("F*M", new List<Matd>() { RS, RP });

            MatdSVD res = new MatdSVD();

            // make B exactly diagonal

            for (int i = 0; i < B.Rows; i++) 
            {
                for (int j = 0; j < B.Columns; j++) 
                {
                    if (i != j)
                    {
                        B.SetCell(i,j,0);
                    }
                }
            }

            res.U = LS;
            res.S = B;
            res.V = RS;

            return res;
        }
    }
}
