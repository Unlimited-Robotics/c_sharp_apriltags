using System.Collections;
using System.Collections.Generic;
using System;

namespace Apriltags.Utils
{
    public static class G2D 
    {
        public static int PolygonOverlapsPolygon(double[][] polya, double[][] polyb)
        {
            // do any of the line segments collide? If so, the answer is yes.
            if (g2dPolygonIntersectsPolygon(polya, polyb) != 0)
            {
                return 1;
            }

            // if none of the edges cross, then the polygon is either fully
            // contained or fully outside.
            double[] p = new double[2];
            polygonGetInteriorPoint(polyb, p);

            if (polygonContainsPoint(polya, p) != 0)
            {
                return 1;
            }

            polygonGetInteriorPoint(polya, p);

            if (polygonContainsPoint(polyb, p) != 0)
            {
                return 1;
            }

            return 0;
        }

        private static void polygonGetInteriorPoint(double[][] poly, double[] p)
        {
            // take the first three points, which form a triangle. Find the middle point
            double[] a = poly[0];
            double[] b = poly[1]; 
            double[] c = poly[2];

            p[0] = (a[0]+b[0]+c[0])/3;
            p[1] = (a[1]+b[1]+c[1])/3;
        }

        private static int polygonContainsPoint(double[][] poly, double[] q)
        {
            // use winding. If the point is inside the polygon, we'll wrap
            // around it (accumulating 6.28 radians). If we're outside the
            // polygon, we'll accumulate zero.
            int psz = poly.Length;

            int last_quadrant = 0;
            int quad_acc = 0;

            for (int i = 0; i <= psz; i++) 
            {
                double[] p = poly[i%psz];

                // p[0] < q[0]       p[1] < q[1]    quadrant
                //     0                 0              0
                //     0                 1              3
                //     1                 0              1
                //     1                 1              2

                // p[1] < q[1]       p[0] < q[0]    quadrant
                //     0                 0              0
                //     0                 1              1
                //     1                 0              3
                //     1                 1              2

                int quadrant;
                if (p[0] < q[0])
                {
                    quadrant = (p[1] < q[1]) ? 2 : 1;
                }
                else
                {
                    quadrant = (p[1] < q[1]) ? 3 : 0;
                }

                if (i > 0) 
                {
                    int dquadrant = quadrant - last_quadrant;

                    // encourage a jump table by mapping to small positive integers.
                    switch (dquadrant) {
                        case -3:
                        case 1:
                            quad_acc ++;
                            break;
                        case -1:
                        case 3:
                            quad_acc --;
                            break;
                        case 0:
                            break;
                        case -2:
                        case 2:
                        {
                            // get the previous point.
                            double[] p0 = poly[i-1];

                            // Consider the points p0 and p (the points around the
                            //polygon that we are tracing) and the query point q.
                            //
                            // If we've moved diagonally across quadrants, we want
                            // to measure whether we have rotated +PI radians or
                            // -PI radians. We can test this by computing the dot
                            // product of vector (p0-q) with the vector
                            // perpendicular to vector (p-q)
                            double nx = p[1] - q[1];
                            double ny = -p[0] + q[0];

                            double dot = nx*(p0[0]-q[0]) + ny*(p0[1]-q[1]);
                            if (dot < 0)
                                quad_acc -= 2;
                            else
                                quad_acc += 2;

                            break;
                        }
                    }
                }

                last_quadrant = quadrant;
            }

            int v = Convert.ToInt32((quad_acc >= 2) || (quad_acc <= -2));

            return v;
        }

        private static int g2dPolygonIntersectsPolygon(double[][] polya, double[][] polyb)
        {
            // do any of the line segments collide? If so, the answer is no.

            // dumb N^2 method.
            for (int ia = 0; ia <polya.Length; ia++) 
            {
                double[] pa0 = polya[ia];
                double[] pa1 = polya[(ia+1)%polya.Length];

                G2DLineSegment sega = new G2DLineSegment(pa0,pa1);

                for (int ib = 0; ib < polyb.Length; ib++) 
                {
                    double[] pb0 = polyb[ib];
                    double[] pb1 = polyb[(ib+1)%polyb.Length];

                    G2DLineSegment segb = new G2DLineSegment(pb0,pb1);

                    if (G2DLineSegment.G2DLineSegmentIntersectSegment(sega, segb, null) != 0)
                    {
                        return 1;
                    }
                }
            }

            return 0;
        }

        public class G2DLineSegment
        {
            public G2DLine Line;
            public double[] P1;

            public G2DLineSegment(double[] p0, double[] p1)
            {
                Line = new G2DLine(p0,p1);
                P1 = new double[2];
                P1[0] = p1[0];
                P1[1] = p1[1];
            }

            public static int G2DLineSegmentIntersectSegment(G2DLineSegment sega, G2DLineSegment segb, double[] p)
            {
                double[] tmp = new double[2];

                if (G2DLine.G2DLineIntersectLine(sega.Line, segb.Line, tmp) == 0)
                {
                    return 0;
                }

                double a = sega.Line.GetCoordinate(sega.Line.P);
                double b = sega.Line.GetCoordinate(sega.P1);
                double c = sega.Line.GetCoordinate(tmp);

                // does intersection lie on the first line?
                if ((c<a && c<b) || (c>a && c>b))
                {
                    return 0;
                }

                a = segb.Line.GetCoordinate(segb.Line.P);
                b = segb.Line.GetCoordinate(segb.P1);
                c = segb.Line.GetCoordinate(tmp);

                // does intersection lie on second line?
                if ((c<a && c<b) || (c>a && c>b))
                {
                    return 0;
                }

                if (p != null) 
                {
                    p[0] = tmp[0];
                    p[1] = tmp[1];
                }

                return 1;
            }
        }

        public class G2DLine
        {
            public double[] P;
            public double[] U;

            public G2DLine(double[] p0, double[] p1)
            {
                P = new double[2];
                U = new double[2];
                P[0] = p0[0];
                P[1] = p0[1];
                U[0] = p1[0]-p0[0];
                U[1] = p1[1]-p0[1];
                double mag = Math.Sqrt(U[0]*U[0] + U[1]*U[1]);

                U[0] /= mag;
                U[1] /= mag;
            }

            public static int G2DLineIntersectLine(G2DLine linea, G2DLine lineb, double[] p)
            {
                // this implementation is many times faster than the original,
                // mostly due to avoiding a general-purpose LU decomposition in
                // Matrix.inverse().
                double m00, m01, m10, m11;
                double i00, i01;
                double b00, b10;

                m00 = linea.U[0];
                m01= -lineb.U[0];
                m10 = linea.U[1];
                m11= -lineb.U[1];

                // determinant of m
                double det = m00*m11-m01*m10;

                // parallel lines?
                if (Math.Abs(det) < 0.00000001)
                {
                    return 0;
                }

                // inverse of m
                i00 = m11/det;
                i01 = -m01/det;

                b00 = lineb.P[0] - linea.P[0];
                b10 = lineb.P[1] - linea.P[1];

                double x00; //, x10;
                x00 = i00*b00+i01*b10;

                if (p != null) 
                {
                    p[0] = linea.U[0]*x00 + linea.P[0];
                    p[1] = linea.U[1]*x00 + linea.P[1];
                }

                return 1;
            }

            public double GetCoordinate(double[] q)
            {
                return (q[0]-P[0])*U[0] + (q[1]-P[1])*U[1];
            }
        }
    }
}
