using System.Collections;
using System.Collections.Generic;

namespace Apriltags
{
    public class CameraParams
    {
        public double Fx;
        public double Fy;
        public double Cx;
        public double Cy;
        public float TagSize;

        public CameraParams(double fx, double fy, double cx, double cy, float tagSize)
        {
            Fx = fx;
            Fy = fy;
            Cx = cx;
            Cy = cy;
            TagSize = tagSize;
        }
    }
}
