using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace Apriltags
{
    public class Detection
    {
        public ApriltagFamily Family;
        public int ID;
        public int Hamming;
        public float DecisionMargin;
        public Matd H;
        public double[] Center;
        public double[][] Corners;
        public TagPose Pose;

        public Detection()
        {
            Center = new double[2];
            Corners = new double[4][];
            Corners[0] = new double[2];
            Corners[1] = new double[2];
            Corners[2] = new double[2];
            Corners[3] = new double[2];
        }
    }
}
