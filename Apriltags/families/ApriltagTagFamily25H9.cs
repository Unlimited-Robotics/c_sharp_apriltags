using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace Apriltags
{
    public class TagFamily25H9 : ApriltagFamily
    {
        public TagFamily25H9()
        {
            Name = "tag25h9";
            H = 9;

            Codes = new string[35];
            Codes[0] = "0156f1f4";
            Codes[1] = "01f28cd5";
            Codes[2] = "016ce32c";
            Codes[3] = "01ea379c";
            Codes[4] = "01390f89";
            Codes[5] = "0034fad0";
            Codes[6] = "007dcdb5";
            Codes[7] = "0119ba95";
            Codes[8] = "01ae9daa";
            Codes[9] = "00df02aa";
            Codes[10] = "0082fc15";
            Codes[11] = "00465123";
            Codes[12] = "00ceee98";
            Codes[13] = "01f17260";
            Codes[14] = "014429cd";
            Codes[15] = "017248a8";
            Codes[16] = "016ad452";
            Codes[17] = "009670ad";
            Codes[18] = "016f65b2";
            Codes[19] = "00b8322b";
            Codes[20] = "005d715b";
            Codes[21] = "01a1c7e7";
            Codes[22] = "00d7890d";
            Codes[23] = "01813522";
            Codes[24] = "01c9c611";
            Codes[25] = "0099e4a4";
            Codes[26] = "00855234";
            Codes[27] = "017b81c0";
            Codes[28] = "00c294bb";
            Codes[29] = "0089fae3";
            Codes[30] = "0044df5f";
            Codes[31] = "01360159";
            Codes[32] = "00ec31e8";
            Codes[33] = "01bcc0f6";
            Codes[34] = "00a64f8d";

            BitX = new int[25];
            BitY = new int[25];
            BitX[0] = 1;
            BitY[0] = 1;
            BitX[1] = 2;
            BitY[1] = 1;
            BitX[2] = 3;
            BitY[2] = 1;
            BitX[3] = 4;
            BitY[3] = 1;
            BitX[4] = 2;
            BitY[4] = 2;
            BitX[5] = 3;
            BitY[5] = 2;
            BitX[6] = 5;
            BitY[6] = 1;
            BitX[7] = 5;
            BitY[7] = 2;
            BitX[8] = 5;
            BitY[8] = 3;
            BitX[9] = 5;
            BitY[9] = 4;
            BitX[10] = 4;
            BitY[10] = 2;
            BitX[11] = 4;
            BitY[11] = 3;
            BitX[12] = 5;
            BitY[12] = 5;
            BitX[13] = 4;
            BitY[13] = 5;
            BitX[14] = 3;
            BitY[14] = 5;
            BitX[15] = 2;
            BitY[15] = 5;
            BitX[16] = 4;
            BitY[16] = 4;
            BitX[17] = 3;
            BitY[17] = 4;
            BitX[18] = 1;
            BitY[18] = 5;
            BitX[19] = 1;
            BitY[19] = 4;
            BitX[20] = 1;
            BitY[20] = 3;
            BitX[21] = 1;
            BitY[21] = 2;
            BitX[22] = 2;
            BitY[22] = 4;
            BitX[23] = 2;
            BitY[23] = 3;
            BitX[24] = 3;
            BitY[24] = 3;

            WidthAtBorder = 7;
            TotalWidth = 9;
            ReversedBorder = false;
        }
    }
}