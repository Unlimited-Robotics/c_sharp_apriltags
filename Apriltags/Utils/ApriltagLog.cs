using System.Collections;
using System.Collections.Generic;
using System;
using System.IO;
// unity content start 
using UnityEngine;
//unity content end

namespace Apriltags.Utils
{
    public static class Log
    {

        // unity content start
        public static int[,,] ConvertTexture2dToArray(Texture2D tex)
        {
            // Debug.Log("tex dims " + tex.width + ", " + tex.height);
            int[,,] output = new int[tex.height, tex.width, 3];

            Color[] pixels = tex.GetPixels();

            for (int i = 0; i < tex.width; i++)
            {
                for (int j = 0; j < tex.height; j++)
                {
                    Color pixel = pixels[i + j * tex.width];
                    Color32 pixelBig = pixel;
                    // output[i, j, 0] = Mathf.Round(pixel.r * 255);
                    // output[i, j, 1] = Mathf.Round(pixel.g * 255);
                    // output[i, j, 2] = Mathf.Round(pixel.b * 255);
                    output[tex.height - 1 - j, i, 0] = pixelBig.b;
                    output[tex.height - 1 - j, i, 1] = pixelBig.g;
                    output[tex.height - 1 - j, i, 2] = pixelBig.r;
                }
            }
        
            return output;
        }
        // unity contnet end

        public static void SaveImageDataToFile(string fileName, int[,,] imageData)
        {
            string path = "/home/ros2/_Alon/Compares/Apriltags/unity/";
            using StreamWriter file = new(path + fileName + ".txt");
            for (int i = 0; i < imageData.GetLength(0); i++)
            {
                for (int j = 0; j < imageData.GetLength(1); j++)
                {
                    file.WriteLine(imageData[i, j, 0] + " " + imageData[i, j, 1] + " " + imageData[i, j, 2]);
                }
            }
        }

        public static void SaveImageDataToFile(string fileName, Image imageData)
        {
            string path = "/home/ros2/_Alon/Compares/Apriltags/unity/";
            using StreamWriter file = new(path + fileName + ".txt");
            for (int y = 0; y < imageData.Height; y++)
            {
                for (int x = 0; x < imageData.Width; x++)
                {
                    file.WriteLine(imageData.GetPixel(x, y).ToString());
                }
            }
        }

        public static void SaveUnionfindDataToFile(string fileName, UnionFind unionfindData)
        {
            string path = "/home/ros2/_Alon/Compares/Apriltags/unity/";
            using StreamWriter file = new(path + fileName + ".txt");
            file.WriteLine("max id " + unionfindData.MaxID.ToString());
            for (int i = 0; i < unionfindData.Data.Length; i++)
            {
                file.WriteLine("parent " + unionfindData.Data[i].Parent + " size " + unionfindData.Data[i].Size);
            }
        }

        public static void SaveClustersDataToFile(string fileName, List<Cluster.ClusterHash> clustersData)
        {
            string path = "/home/ros2/_Alon/Compares/Apriltags/unity/";
            using StreamWriter file = new(path + fileName + ".txt");
            for (int i = 0; i < clustersData.Count; i++)
            {
                file.WriteLine("cluster hash " + clustersData[i].Hash);
                file.WriteLine("cluster id " + clustersData[i].ID);
                for (int j = 0; j < clustersData[i].Data.Points.Count; j++)
                {
                    Cluster.ClusterPoint p = clustersData[i].Data.Points[j];
                    file.WriteLine("cluster point index " + j);
                    file.WriteLine("x " + p.X + ", y " + p.Y + ", gx " + p.GX + ", gy " + p.GY + ", slope " + p.Slope);
                }
            }
        }

        public static void SaveClustersDataToFile(string fileName, List<Cluster> clustersData)
        {
            string path = "/home/ros2/_Alon/Compares/Apriltags/unity/";
            using StreamWriter file = new(path + fileName + ".txt");
            for (int i = 0; i < clustersData.Count; i++)
            {
                file.WriteLine("cluster size " + clustersData[i].Points.Count);
                for (int j = 0; j < clustersData[i].Points.Count; j++)
                {
                    Cluster.ClusterPoint p = clustersData[i].Points[j];
                    file.WriteLine("cluster point index " + j);
                    file.WriteLine("x " + p.X + ", y " + p.Y + ", gx " + p.GX + ", gy " + p.GY + ", slope " + p.Slope.ToString("F6"));
                }
            }
        }

        public static void SaveClustersDataToFile(string fileName, Cluster clustersData)
        {
            string path = "/home/ros2/_Alon/Compares/Apriltags/unity/";
            using StreamWriter file = new(path + fileName + ".txt");
            file.WriteLine("cluster size " + clustersData.Points.Count);
            for (int j = 0; j < clustersData.Points.Count; j++)
            {
                Cluster.ClusterPoint p = clustersData.Points[j];
                file.WriteLine("cluster point index " + j);
                file.WriteLine("x " + p.X + ", y " + p.Y + ", gx " + p.GX + ", gy " + p.GY + ", slope " + p.Slope.ToString("F2"));
            }
        }

        public static void SaveQuadsDataToFile(string fileName, List<Quad> quadsData)
        {
            string path = "/home/ros2/_Alon/Compares/Apriltags/unity/";
            using StreamWriter file = new(path + fileName + ".txt");
            for (int i = 0; i < quadsData.Count; i++)
            {
                file.WriteLine("quad index " + i);
                for (int j = 0; j < quadsData[i].Corners.Length; j++)
                {
                    double[] p = quadsData[i].Corners[j];
                    file.WriteLine("quad corner index " + j);
                    file.WriteLine("x " + p[0].ToString("F2") + ", y " + p[1].ToString("F2"));
                }
                file.WriteLine("quad reversed border " + Convert.ToInt32(quadsData[i].ReversedBorder));
            }
        }

        public static void SaveFitLinesDataToFile(string fileName, List<FitLine> linesData)
        {
            string path = "/home/ros2/_Alon/Compares/Apriltags/unity/";
            using StreamWriter file = new(path + fileName + ".txt");
            for (int i = 0; i < linesData.Count; i++)
            {
                file.WriteLine("fit line index " + i);
                file.WriteLine("mx " + linesData[i].Mx.ToString("F2") + ", my " + linesData[i].My.ToString("F2") + 
                ", mxx " + linesData[i].Mxx.ToString("F2") + ", myy " + linesData[i].Myy.ToString("F2") +
                ", mxy " + linesData[i].Mxy.ToString("F2") + ", w " + linesData[i].W.ToString("F2"));
            }
        }

        public static void SaveDetectionsDataToFile(string fileName, List<Detection> detsData)
        {
            string path = "/home/ros2/_Alon/Compares/Apriltags/unity/";
            using StreamWriter file = new(path + fileName + ".txt");
            for (int i = 0; i < detsData.Count; i++)
            {
                file.WriteLine("detection index " + i);
                file.WriteLine("id " + detsData[i].ID + ", hamming " + detsData[i].Hamming + 
                ", decision margine " + detsData[i].DecisionMargin.ToString("F2") +
                ", center " + detsData[i].Center[0].ToString("F2") + " " + detsData[i].Center[1].ToString("F2"));
                for (int j = 0; j < detsData[i].Corners.Length; j++)
                {
                    double[] p = detsData[i].Corners[j];
                    file.WriteLine("detection corner index " + j);
                    file.WriteLine("x " + p[0].ToString("F2") + ", y " + p[1].ToString("F2"));
                }

                file.WriteLine("homography");
                for (int j = 0; j < 3; j++)
                {
                    file.WriteLine("[" + detsData[i].H.GetCell(j,0).ToString("F2") + " " + 
                        detsData[i].H.GetCell(j,1).ToString("F2") + " " + detsData[i].H.GetCell(j,2).ToString("F2") + "]");
                }

                if(detsData[i].Pose != null)
                {
                    file.WriteLine("pose_R");
                    Matd m = detsData[i].Pose.R;
                    for (int j = 0; j < 3; j++)
                    {
                        file.WriteLine("[" + m.GetCell(j,0).ToString("F2") + " " + 
                            m.GetCell(j,1).ToString("F2") + " " + m.GetCell(j,2).ToString("F2") + "]");
                    }
                    file.WriteLine("pose_t");
                    m = detsData[i].Pose.T;
                    file.WriteLine("[" + m.GetCell(0,0).ToString("F2") + " " + 
                        m.GetCell(0,1).ToString("F2") + " " + m.GetCell(0,2).ToString("F2") + "]");
                    file.WriteLine("pose_err " + detsData[i].Pose.Err);
                }
            }
        }
    }
}
