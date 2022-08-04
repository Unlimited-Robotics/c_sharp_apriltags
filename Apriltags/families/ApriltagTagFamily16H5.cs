using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace Apriltags
{
    public class TagFamily16H5 : ApriltagFamily
    {
        public TagFamily16H5()
        {
            Name = "tag16h5";
            H = 5;

            Codes = new string[30];
            Codes[0] = "27c8"; //0
            Codes[1] = "31b6"; //1
            Codes[2] = "3859"; //2
            Codes[3] = "569c"; //3
            Codes[4] = "6c76"; //4
            Codes[5] = "7ddb"; //5
            Codes[6] = "af09"; //6
            Codes[7] = "f5a1"; //7
            Codes[8] = "fb8b"; //8
            Codes[9] = "1cb9"; //9
            Codes[10] = "28ca"; //10
            Codes[11] = "e8dc"; //11
            Codes[12] = "1426"; //12
            Codes[13] = "5770"; //13
            Codes[14] = "9253"; //14
            Codes[15] = "b702"; //15
            Codes[16] = "63a"; //16
            Codes[17] = "8f34"; //17
            Codes[18] = "b4c0"; //18
            Codes[19] = "51ec"; //19
            Codes[20] = "e6f0"; //20
            Codes[21] = "5fa4"; //21
            Codes[22] = "dd43"; //22
            Codes[23] = "1aaa"; //23
            Codes[24] = "e62f"; //24
            Codes[25] = "6dbc"; //25
            Codes[26] = "b6eb"; //26
            Codes[27] = "de10"; //27
            Codes[28] = "154d"; //28
            Codes[29] = "b57a"; //29

            BitX = new int[16];
            BitY = new int[16];
            BitX[0] = 1; //0
            BitY[0] = 1; //0
            BitX[1] = 2; //1
            BitY[1] = 1; //1
            BitX[2] = 3; //2
            BitY[2] = 1; //2
            BitX[3] = 2; //3
            BitY[3] = 2; //3
            BitX[4] = 4; //4
            BitY[4] = 1; //4
            BitX[5] = 4; //5
            BitY[5] = 2; //5
            BitX[6] = 4; //6
            BitY[6] = 3; //6
            BitX[7] = 3; //7
            BitY[7] = 2; //7
            BitX[8] = 4; //8
            BitY[8] = 4; //8
            BitX[9] = 3; //9
            BitY[9] = 4; //9
            BitX[10] = 2; //10
            BitY[10] = 4; //10
            BitX[11] = 3; //11
            BitY[11] = 3; //11
            BitX[12] = 1; //12
            BitY[12] = 4; //12
            BitX[13] = 1; //13
            BitY[13] = 3; //13
            BitX[14] = 1; //14
            BitY[14] = 2; //14
            BitX[15] = 2; //15
            BitY[15] = 3; //15

            WidthAtBorder = 6;
            TotalWidth = 8;
            ReversedBorder = false;
        }
    }
}