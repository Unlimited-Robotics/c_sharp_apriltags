using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace Apriltags
{
    public class GrayModel
    {
        public double[][] A;
        public double[] B;
        public double[] C;

        public GrayModel()
        {
            A = new double[3][];
            A[0] = new double[3];
            A[1] = new double[3];
            A[2] = new double[3];
            B = new double[3];
            C = new double[3];
        }

        public void Add(double x, double y, double gray)
        {
            // update upper right entries of A = J'J
            A[0][0] += x*x;
            A[0][1] += x*y;
            A[0][2] += x;
            A[1][1] += y*y;
            A[1][2] += y;
            A[2][2] += 1;

            // update B = J'gray
            B[0] += x * gray;
            B[1] += y * gray;
            B[2] += gray;
        }

        public void Solve()
        {
            Utils.Calculations.Mat33SymSolve(A.SelectMany(a => a).ToArray(), B, C);
        }

        public double Interpolate(double x, double y)
        {
            return C[0]*x + C[1]*y + C[2];
        }
    }
}
