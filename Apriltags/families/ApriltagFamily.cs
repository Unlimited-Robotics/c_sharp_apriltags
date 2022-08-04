using System.Collections;
using System.Collections.Generic;

namespace Apriltags
{
    public abstract class ApriltagFamily
    {
        public string[] Codes;
        public int WidthAtBorder;
        public int TotalWidth;
        public bool ReversedBorder;
        public int[] BitX;
        public int[] BitY;

        // minimum hamming distance between any two codes. (e.g. 36h11 => 11)
        public uint H; 
        public string Name;
        public QuickDecode Implementation;

        public ulong getCodeAsInt(int index)
        {
            return ulong.Parse(Codes[index], System.Globalization.NumberStyles.HexNumber);
        }

        public void QuickDecodeCodeword(ulong rcode, out QuickDecodeEntry entry)
        {
            for (int ridx = 0; ridx < 4; ridx++) 
            {

                for (int bucket = (int)(rcode % (ulong)Implementation.Entries.Length);
                    Implementation.Entries[bucket].RCode != ulong.MaxValue;
                    bucket = (bucket + 1) % Implementation.Entries.Length) {

                    if (Implementation.Entries[bucket].RCode == rcode) {
                        entry = Implementation.Entries[bucket];
                        entry.Rotation = (byte)ridx;
                        return;
                    }
                }

                rcode = Utils.Calculations.Rotate90(rcode, BitX.Length);
            }

            entry = new QuickDecodeEntry();
            entry.RCode = 0;
            entry.ID = 65535;
            entry.Hamming = 255;
            entry.Rotation = 0;
        }
    }
}

