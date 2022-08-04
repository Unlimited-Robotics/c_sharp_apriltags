using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace Apriltags
{
    public class TagFamilyCircle21H7 : ApriltagFamily
    {
        public TagFamilyCircle21H7()
        {
            Name = "tagCircle21h7";
            H = 7;

            Codes = new string[38];
            Codes[0] = "0157863";
            Codes[1] = "0047e28";
            Codes[2] = "01383ed";
            Codes[3] = "000953c";
            Codes[4] = "00da68b";
            Codes[5] = "01cac50";
            Codes[6] = "00bb215";
            Codes[7] = "016ceee";
            Codes[8] = "005d4b3";
            Codes[9] = "01ff751";
            Codes[10] = "00efd16";
            Codes[11] = "0072b3e";
            Codes[12] = "0163103";
            Codes[13] = "0106e56";
            Codes[14] = "01996b9";
            Codes[15] = "00c0234";
            Codes[16] = "00624d2";
            Codes[17] = "01fa985";
            Codes[18] = "00344a5";
            Codes[19] = "00762fb";
            Codes[20] = "019e92b";
            Codes[21] = "0043755";
            Codes[22] = "001a4f4";
            Codes[23] = "010fad8";
            Codes[24] = "0001b52";
            Codes[25] = "017e59f";
            Codes[26] = "00e6f70";
            Codes[27] = "00ed47a";
            Codes[28] = "00c9931";
            Codes[29] = "0014df2";
            Codes[30] = "00a06f1";
            Codes[31] = "00e5041";
            Codes[32] = "012ec03";
            Codes[33] = "016724e";
            Codes[34] = "00af1a5";
            Codes[35] = "008a8ac";
            Codes[36] = "0015b39";
            Codes[37] = "01ec1e3";

            BitX = new int[21];
            BitY = new int[21];
            BitX[0] = 1;
            BitY[0] = -2;
            BitX[1] = 2;
            BitY[1] = -2;
            BitX[2] = 3;
            BitY[2] = -2;
            BitX[3] = 1;
            BitY[3] = 1;
            BitX[4] = 2;
            BitY[4] = 1;
            BitX[5] = 6;
            BitY[5] = 1;
            BitX[6] = 6;
            BitY[6] = 2;
            BitX[7] = 6;
            BitY[7] = 3;
            BitX[8] = 3;
            BitY[8] = 1;
            BitX[9] = 3;
            BitY[9] = 2;
            BitX[10] = 3;
            BitY[10] = 6;
            BitX[11] = 2;
            BitY[11] = 6;
            BitX[12] = 1;
            BitY[12] = 6;
            BitX[13] = 3;
            BitY[13] = 3;
            BitX[14] = 2;
            BitY[14] = 3;
            BitX[15] = -2;
            BitY[15] = 3;
            BitX[16] = -2;
            BitY[16] = 2;
            BitX[17] = -2;
            BitY[17] = 1;
            BitX[18] = 1;
            BitY[18] = 3;
            BitX[19] = 1;
            BitY[19] = 2;
            BitX[20] = 2;
            BitY[20] = 2;

            WidthAtBorder = 5;
            TotalWidth = 9;
            ReversedBorder = true;
        }
    }
}
