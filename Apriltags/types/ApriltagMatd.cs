using System.Collections;
using System.Collections.Generic;
using System.Text.RegularExpressions;
using System;

namespace Apriltags
{
    public class Matd
    {
        public uint Rows;
        public uint Columns;
        public double[] Data;

        public Matd(int rows, int cols)
        {
            Rows = (uint)rows;
            Columns = (uint)cols;
            Data = new double[Rows*Columns];
        }

        public Matd(int rows, int cols, double[] data)
        {
            Rows = (uint)rows;
            Columns = (uint)cols;
            Data = new double[Rows*Columns];
            for (int i = 0; i < Data.Length; i++)
            {
                Data[i] = data[i];
            }
        }

        public Matd(double value)
        {
            Rows = Columns = 0;
            Data = new double[] {value};
        }

        public Matd(Matd m)
        {
            Rows = m.Rows;
            Columns = m.Columns;
            Data = new double[m.Data.Length];
            for (int i = 0; i < m.Data.Length; i++)
            {
                Data[i] = m.Data[i];
            }
        }

        public bool IsScalar()
        {
            bool output = false;

            if(Columns <= 1 && Rows <= 1)
            {
                output = true;
            }

            return output;
        }

        public int MaxIndex(int row, int maxcol)
        {
            int maxi = 0;
            double maxv = -1;

            for (int i = 0; i < maxcol; i++) 
            {
                if (i == row)
                {
                    continue;
                }
                double v = Math.Abs(GetCell(row,i));
                if (v > maxv) 
                {
                    maxi = i;
                    maxv = v;
                }
            }

            return maxi;
        }

        public void AddInPlace(Matd b)
        {
            if (IsScalar() == true) 
            {
                Data[0] += b.Data[0];
                return;
            }

            for (int i = 0; i < Rows; i++) 
            {
                for (int j = 0; j < Columns; j++) 
                {
                    SetCell(i,j,GetCell(i,j) + b.GetCell(i,j));
                }
            }
        }

        public void ScaleInPlace(double s)
        {
            if (IsScalar() == true) 
            {
                Data[0] *= s;
                return;
            }

            for (int i = 0; i < Rows; i++) 
            {
                for (int j = 0; j < Columns; j++) 
                {
                    SetCell(i,j,GetCell(i,j)*s);
                }
            }
        }

        public double ToDouble()
        {
            double output = 0;

            if(IsScalar() == false)
            {
                // Debug.LogError("cant convert none scalat matrix to double");
            }
            else
            {
                output = Data[0];
            }

            return output;
        }

        public static Matd GetIdentityMatrix(int dim)
        {
            Matd output;

            if (dim == 0)
            {
                output = new Matd(1);
            }
            else
            {
                output = new Matd(dim, dim);
                for (int i = 0; i < dim; i++)
                {
                    output.SetCell(i, i, 1);
                }
            }

            return output;
        }

        public double GetCell(int row, int col)
        {
            return Data[row*Columns + col];
        }

        public void SetCell(int row, int col, double value)
        {
            Data[row*Columns + col] = value;
        }

        public Matd GetInverseMatrix()
        {

            if (IsScalar() == true) 
            {
                if (Data[0] == 0)
                {
                    return null;
                }

                return new Matd(1.0 / Data[0]);
            }

            switch(Rows) 
            {
                case 1: 
                {
                    double det = Data[0];
                    if (det == 0)
                    {
                        return null;
                    }

                    double invdet = 1.0 / det;

                    Matd m = new Matd((int)Rows, (int)Rows);
                    m.SetCell(0,0, 1.0 * invdet);
                    return m;
                }

                case 2: 
                {
                    double det = Data[0] * Data[3] - Data[1] * Data[2];
                    if (det == 0)
                    {
                        return null;
                    }

                    double invdet = 1.0 / det;

                    Matd m = new Matd((int)Rows, (int)Rows);
                    m.SetCell(0,0,GetCell(1,1)*invdet);
                    m.SetCell(0,1,-GetCell(0,1)*invdet);
                    m.SetCell(1,0,-GetCell(1,0)*invdet);
                    m.SetCell(1,1,GetCell(0,0)*invdet);
                    return m;
                }

                default: 
                {
                    MatdPlu plu = new MatdPlu(this);

                    Matd inv = null;
                    if (plu.Singular == 0) 
                    {
                        Matd ident = GetIdentityMatrix((int)Rows);
                        inv = plu.Solve(ident);
                    }

                    return inv;
                }
            }

            // return null; // unreachable
        }

        public Matd GetTransposeMatrix()
        {
            if (IsScalar() == true)
            {
                return new Matd(Data[0]);
            }

            Matd m = new Matd((int)Columns, (int)Rows);

            for (int i = 0; i < Rows; i++) 
            {
                for (int j = 0; j < Columns; j++) 
                {
                    m.SetCell(j,i,GetCell(i,j));
                }
            }

            return m;
        }

        public Matd GetScaledMatrix(double s)
        {
            if (IsScalar() == true)
            {
                return new Matd(Data[0]*s);
            }

            Matd m = new Matd((int)Rows, (int)Columns);

            for (int i = 0; i < m.Rows; i++) 
            {
                for (int j = 0; j < m.Columns; j++) 
                {
                    m.SetCell(i,j, s*GetCell(i,j));
                }
            }

            return m;
        }

        public Matd GetVecNormalizeMatrix()
        {
            double mag = GetVecMag();

            Matd b = new Matd((int)Rows, (int)Columns);

            int len = (int)(Rows*Columns);
            for(int i = 0; i < len; i++)
            {
                b.Data[i] = Data[i] / mag;
            }

            return b;
        }

        public double GetVecMag()
        {
            double mag = 0.0;
            int len = (int)(Rows*Columns);
            for (int i = 0; i < len; i++)
            {
                mag += Data[i]*Data[i];
            }
            return Math.Sqrt(mag);
        }

        public static Matd MatdOp(string expr, List<Matd> matricies)
        {
            int nargs = 0;
            int exprlen = 0;


            for (int i = 0; i < expr.Length; i++) 
            {
                if (expr[i] == 'M' || expr[i] == 'F')
                {
                    nargs++;
                }
                exprlen++;
            }

            if (expr.Length == 0) // expr = ""
            {
                return null;
            }

            int pos = 0;
            int argpos = 0;
            int garbpos = 0;

            // can't create more than 2 new result per character
            // one result, and possibly one argument to free
            // matd_t **garb = malloc(sizeof(matd_t*)*2*exprlen);
            Matd[] garb = new Matd[2*exprlen];

            Matd res = matdOpRecurse(expr, ref pos, null, matricies, ref argpos, garb, ref garbpos, 0);

            // 'res' may need to be freed as part of garbage collection (i.e. expr = "F")
            Matd res_copy = (res != null ? new Matd(res) : null);

            return res_copy;
        }

        private static Matd matdOpRecurse(string expr, ref int pos, Matd acc, List<Matd> args, ref int argpos,
            Matd[] garb, ref int garbpos, int oneterm)
        {
            while (pos < expr.Length) 
            {

                switch (expr[pos]) 
                {

                    case '(': 
                    {
                        if (oneterm != 0 && acc != null)
                        {
                            return acc;
                        }

                        pos++;
                        Matd rhs = matdOpRecurse(expr, ref pos, null, args, ref argpos, garb, ref garbpos, 0);
                        rhs = matdOpGobbleRight(expr, ref pos, rhs, garb, ref garbpos);

                        if (acc == null) 
                        {
                            acc = rhs;
                        }
                        else 
                        {
                            Matd res = MatdMultiply(acc, rhs);
                            garb[garbpos] = res;
                            garbpos++;
                            acc = res;
                        }

                        break;
                    }

                    case ')': 
                    {
                        if (oneterm != 0)
                        {
                            return acc;
                        }

                        pos++;
                        return acc;
                    }

                    case '*': 
                    {
                        pos++;

                        Matd rhs = matdOpRecurse(expr, ref pos, null, args, ref argpos, garb, ref garbpos, 1);
                        rhs = matdOpGobbleRight(expr, ref pos, rhs, garb, ref garbpos);

                        if (acc == null) 
                        {
                            acc = rhs;
                        }
                        else 
                        {
                            Matd res = MatdMultiply(acc, rhs);
                            garb[garbpos] = res;
                            garbpos++;
                            acc = res;
                        }

                        break;
                    }

                    case 'F': 
                    {
                        Matd rhs = args[argpos];
                        garb[garbpos] = rhs;
                        garbpos++;

                        pos++;
                        argpos++;

                        rhs = matdOpGobbleRight(expr, ref pos, rhs, garb, ref garbpos);

                        if (acc == null) 
                        {
                            acc = rhs;
                        } 
                        else 
                        {
                            Matd res = MatdMultiply(acc, rhs);
                            garb[garbpos] = res;
                            garbpos++;
                            acc = res;
                        }

                        break;
                    }

                    case 'M': 
                    {
                        Matd rhs = args[argpos];

                        pos++;
                        argpos++;

                        rhs = matdOpGobbleRight(expr, ref pos, rhs, garb, ref garbpos);

                        if (acc == null) 
                        {
                            acc = rhs;
                        } 
                        else 
                        {
                            Matd res = MatdMultiply(acc, rhs);
                            garb[garbpos] = res;
                            garbpos++;
                            acc = res;
                        }

                        break;
                    }

        /*
        case 'D': {
        int rows = expr[*pos+1]-'0';
        int cols = expr[*pos+2]-'0';

        matd_t *rhs = matd_create(rows, cols);

        break;
        }
        */
                        // a constant (SCALAR) defined inline. Treat just like M, creating a matd_t on the fly.
                    case '0':
                    case '1':
                    case '2':
                    case '3':
                    case '4':
                    case '5':
                    case '6':
                    case '7':
                    case '8':
                    case '9':
                    case '.': 
                    {
                        double number;
                        Utils.Print.Strtod(expr, ref pos, out number);
                        // string start = expr.Substring(pos);
                        // var matches = Regex.Matches(start, @"\d+.\d*");
                        // string number = "";
                        // number += matches[0];
                        // double s = double.Parse(number);
                        // pos += number.Length;
                        Matd rhs = new Matd(number);
                        garb[garbpos] = rhs;
                        garbpos++;

                        rhs = matdOpGobbleRight(expr, ref pos, rhs, garb, ref garbpos);

                        if (acc == null) 
                        {
                            acc = rhs;
                        } 
                        else 
                        {
                            Matd res = MatdMultiply(acc, rhs);
                            garb[garbpos] = res;
                            garbpos++;
                            acc = res;
                        }

                        break;
                    }

                    case '+': 
                    {
                        if (oneterm != 0 && acc != null)
                            return acc;

                        // don't support unary plus
                        pos++;
                        Matd rhs = matdOpRecurse(expr, ref pos, null, args, ref argpos, garb, ref garbpos, 1);
                        rhs = matdOpGobbleRight(expr, ref pos, rhs, garb, ref garbpos);

                        Matd res = matdAdd(acc, rhs);

                        garb[garbpos] = res;
                        garbpos++;
                        acc = res;
                        break;
                    }

                    case '-': 
                    {
                        if (oneterm != 0 && acc != null)
                        {
                            return acc;
                        }

                        if (acc == null) 
                        {
                            // unary minus
                            pos++;
                            Matd rhs = matdOpRecurse(expr, ref pos, null, args, ref argpos, garb, ref garbpos, 1);
                            rhs = matdOpGobbleRight(expr, ref pos, rhs, garb, ref garbpos);

                            Matd res = rhs.GetScaledMatrix(-1);
                            garb[garbpos] = res;
                            garbpos++;
                            acc = res;
                        } 
                        else 
                        {
                            // subtract
                            pos++;
                            Matd rhs = matdOpRecurse(expr, ref pos, null, args, ref argpos, garb, ref garbpos, 1);
                            rhs = matdOpGobbleRight(expr, ref pos, rhs, garb, ref garbpos);

                            Matd res = MatdSubtract(acc, rhs);
                            garb[garbpos] = res;
                            garbpos++;
                            acc = res;
                        }
                        break;
                    }

                    case ' ': 
                    {
                        // nothing to do. spaces are meaningless.
                        pos++;
                        break;
                    }

                    default: 
                    {
                        break;
                    }
                }
            }
            return acc;
        }

        public static Matd matdOpGobbleRight(string expr, ref int pos, Matd acc, Matd[] garb, ref int garbpos)
        {
            while (pos < expr.Length) 
            {

                switch (expr[pos]) 
                {

                    case '\'': 
                    {
                        Matd res = acc.GetTransposeMatrix();
                        garb[garbpos] = res;
                        garbpos++;
                        acc = res;

                        pos++;
                        break;
                    }

                        // handle inverse ^-1. No other exponents are allowed.
                    case '^': 
                    {
                        Matd res = acc.GetInverseMatrix();
                        garb[garbpos] = res;
                        garbpos++;
                        acc = res;

                        pos += 3;
                        break;
                    }

                    default:
                        return acc;
                }
            }

            return acc;
        }

        public static Matd MatdMultiply(Matd a, Matd b)
        {
            if (a.IsScalar() == true)
            {
                return b.GetScaledMatrix(a.Data[0]);
            }
            if (b.IsScalar() == true)
            {
                return a.GetScaledMatrix(b.Data[0]);
            }

            Matd m = new Matd((int)a.Rows, (int)b.Columns);

            for (int i = 0; i < m.Rows; i++) 
            {
                for (int j = 0; j < m.Columns; j++) 
                {
                    double acc = 0;
                    for (int k = 0; k < a.Columns; k++) 
                    {
                        acc += a.GetCell(i,k) * b.GetCell(k,j);
                    }
                    m.SetCell(i,j,acc);
                }
            }

            return m;
        }

        private static Matd matdAdd(Matd a, Matd b)
        {
            if (a.IsScalar() == true)
            {
                return new Matd(a.Data[0] + b.Data[0]);
            }

            Matd m = new Matd((int)a.Rows, (int)a.Columns);

            for (int i = 0; i < m.Rows; i++) 
            {
                for (int j = 0; j < m.Columns; j++) 
                {
                    m.SetCell(i,j,a.GetCell(i,j) + b.GetCell(i,j));
                }
            }

            return m;
        }

        public static Matd MatdSubtract(Matd a, Matd b)
        {
            if (a.IsScalar() == true)
            {
                return new Matd(a.Data[0] - b.Data[0]);
            }

            Matd m = new Matd((int)a.Rows, (int)a.Columns);

            for (int i = 0; i < m.Rows; i++) 
            {
                for (int j = 0; j < m.Columns; j++) 
                {
                    m.SetCell(i,j,a.GetCell(i,j) - b.GetCell(i,j));
                }
            }

            return m;
        }

        public static Matd MatdCrossProduct(Matd a, Matd b)
        { // only defined for vecs (col or row) of length 3
            Matd r = new Matd((int)a.Rows,(int)a.Columns);

            r.Data[0] = a.Data[1] * b.Data[2] - a.Data[2] * b.Data[1];
            r.Data[1] = a.Data[2] * b.Data[0] - a.Data[0] * b.Data[2];
            r.Data[2] = a.Data[0] * b.Data[1] - a.Data[1] * b.Data[0];

            return r;
        }
    }
}
