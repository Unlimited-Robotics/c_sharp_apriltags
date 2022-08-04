using System.Collections;
using System.Collections.Generic;
using System;
using UnityEngine;

namespace Apriltags.Utils
{
    public static class Print
    {
        public static void PrintClustersData(List<Cluster.ClusterHash> clustersData)
        {
            for (int i = 0; i < clustersData.Count; i++)
            {
                Debug.Log("cluster hash " + clustersData[i].Hash);
                Debug.Log("cluster id " + clustersData[i].ID);
                Debug.Log("cluster size " + clustersData[i].Data.Points.Count);
                // for (int j = 0; j < clustersData[i].Data.Points.Count; j++)
                // {
                //     Cluster.ClusterPoint p = clustersData[i].Data.Points[j];
                //     file.WriteLine("cluster point index " + j);
                //     file.WriteLine("x " + p.X + ", y " + p.Y + ", gx " + p.GX + ", gy " + p.GY);
                // }
            }
        }

        public static void Strtod(string expr, ref int pos, out double num)
        {
            string number = "";
            while(expr[pos] == '.' || Char.IsNumber(expr[pos]) == true)
            {
                number += expr[pos];
                pos++;
            }

            num = double.Parse(number);
        }
    }
}
