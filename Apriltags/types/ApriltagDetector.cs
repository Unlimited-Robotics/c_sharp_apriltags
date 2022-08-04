using System.Collections;
using System.Collections.Generic;
using System.Threading;
using System;

namespace Apriltags
{
    public class Detector
    {
        public int HowManyThreadsToUse;

        // detection of quads can be done on a lower-resolution image,
        // improving speed at a cost of pose accuracy and a slight
        // decrease in detection rate. Decoding the binary payload is
        // still done at full resolution. .
        public float QuadDecimate;

        // What Gaussian blur should be applied to the segmented image
        // (used for quad detection?)  Parameter is the standard deviation
        // in pixels.  Very noisy images benefit from non-zero values
        // (e.g. 0.8).
        public float QuadSigma;

        // When non-zero, the edges of the each quad are adjusted to "snap
        // to" strong gradients nearby. This is useful when decimation is
        // employed, as it can increase the quality of the initial quad
        // estimate substantially. Generally recommended to be on (1).
        //
        // Very computationally inexpensive. Option is ignored if
        // quad_decimate = 1.
        public bool RefineEdges;

        // How much sharpening should be done to decoded images? This
        // can help decode small tags but may or may not help in odd
        // lighting conditions or low light conditions.
        //
        // The default value is 0.25.
        public double DecodeSharpening;

        // When non-zero, write a variety of debugging images to the
        // current working directory at various stages through the
        // detection process. (Somewhat slow).

        public QuadThreshParams QuadThreshParams;

        ///////////////////////////////////////////////////////////////
        // Statistics relating to last processed frame
        //timeprofile_t* tp;

        public uint EdgeAmount;
        public uint SegmentAmount;
        public uint QuadAmount;

        ///////////////////////////////////////////////////////////////
        // Internal variables below

        // Not freed on apriltag_destroy; a tag family can be shared
        // between multiple users. The user should ultimately destroy the
        // tag family passed into the constructor.
        public List<ApriltagFamily> TagFamilies;

        // Used to manage multi-threading.
        public WorkerPool WorkPool;

        // Used for thread safety.
        public Mutex MutexLock;

        public bool IsDebug;

        public Detector(List<ApriltagFamily> families, int howManyThreadsToUse=1, float quadDecimate=2, 
            float quadSigma=0, bool refineEdges=true, double decodeSharpening=0.25f, bool isDebug=false)
        {
            foreach(ApriltagFamily family in families)
            {
                loadFamily(family);
            }

            TagFamilies = families;
            HowManyThreadsToUse = howManyThreadsToUse;
            QuadDecimate = quadDecimate;
            QuadSigma = quadSigma;
            RefineEdges = refineEdges;
            DecodeSharpening = (int)decodeSharpening;
            IsDebug = isDebug;
            QuadThreshParams = new QuadThreshParams();
            MutexLock = new Mutex();
        }

        private void loadFamily(ApriltagFamily family, int maxHammingDistance=2)
        {
            uint codeAmount = (uint)family.Codes.Length;
            uint bitAmount = (uint)family.BitX.Length;
            uint capacity = codeAmount;

            if(maxHammingDistance >= 1)
            {
                capacity += codeAmount * bitAmount;
            }

            if(maxHammingDistance >= 2)
            {
                capacity += codeAmount * bitAmount * (bitAmount - 1);
            }

            if(maxHammingDistance >= 3)
            {
                capacity += codeAmount * bitAmount * (bitAmount - 1) * (bitAmount - 2);
            }

            QuickDecode quickDecode = new QuickDecode((int)capacity * 3);

            for (ushort i = 0; i < codeAmount; i++)
            {
                ulong code = family.getCodeAsInt(i);
                quickDecode.AddCode(code, i, 0);

                if(maxHammingDistance >= 1)
                {
                    ulong power = 1;
                    for (ushort j = 0; j < bitAmount; j++)
                    {
                        quickDecode.AddCode(code ^ power, i, 1);
                        power <<= 1;
                    }
                }

                if(maxHammingDistance >= 2)
                {
                    ulong power1 = 1;
                    for (ushort j = 0; j < bitAmount; j++)
                    {
                        ulong power2 = 1;
                        for (ushort k = 0; k < j; k++)
                        {
                            quickDecode.AddCode(code ^ power1 ^ power2, i, 2);
                            power2 <<= 1;
                        }
                        power1 <<= 1;
                    }
                }

                if(maxHammingDistance >= 3)
                {
                    ulong power1 = 1;
                    for (ushort j = 0; j < bitAmount; j++)
                    {
                        ulong power2 = 1;
                        for (ushort k = 0; k < j; k++)
                        {
                            ulong power3 = 1;
                            for (ushort m = 0; m < k; m++)
                            {
                                quickDecode.AddCode(code ^ power1 ^ power2 ^ power3, i, 3);
                                power3 <<= 1;
                            }
                            power2 <<= 1;
                        }
                        power1 <<= 1;
                    }
                }
            }

            family.Implementation = quickDecode;
        }

        public List<Detection> Detect(Image image, CameraParams cameraParams=null)
        {
            Image im_orig = new Image(image);
            List<Detection> output = new List<Detection>();

            WorkPool = new WorkerPool(HowManyThreadsToUse);

            if(QuadDecimate > 1)
            {
                image.Decimate(QuadDecimate);
                if(IsDebug == true)
                {
                    Utils.Log.SaveImageDataToFile("image_decimate", image);
                }
            }

            if(QuadSigma != 0)
            {
                float sigma = Math.Abs(QuadSigma);
                int ksz = (int)(4 * sigma);

                if((ksz & 1) == 0)
                {
                    ksz++;
                }

                if(ksz > 1)
                {
                    if(QuadSigma > 0)
                    {
                        image.GaussianBlur(sigma, ksz);
                    }
                    else
                    {
                        // Debug.Log("ned quad");
                        Image preBlur = new Image(image);
                        image.GaussianBlur(sigma, ksz);

                        for (int y = 0; y < preBlur.Height; y++)
                        {
                            for (int x = 0; x < preBlur.Width; x++) 
                            {
                                int vorig = preBlur.GetPixel(x, y);
                                int vblur = image.GetPixel(x, y);

                                int v = 2 * vorig - vblur;
                                if (v < 0)
                                {
                                    v = 0;
                                }
                                if (v > 255)
                                {
                                    v = 255;
                                }
                                // Debug.Log("vorig " + vorig);
                                // Debug.Log("vblur " + vblur);
                                // Debug.Log("v " + v);

                                image.SetPixel(x, y, (byte)v);
                            }
                        }
                    }
                }

                if(IsDebug == true)
                {
                    Utils.Log.SaveImageDataToFile("image_sigma", image);
                }
            }

            List<Quad> quads = QuadThresh.GetQuadThresh(this, image);
            if(IsDebug == true)
            {
                Utils.Log.SaveQuadsDataToFile("quads_start", quads);
            }

            if (QuadDecimate > 1) 
            {
                for (int i = 0; i < quads.Count; i++) 
                {
                    Quad q = quads[i];

                    for (int j = 0; j < 4; j++) 
                    {
                        if (QuadDecimate == 1.5) 
                        {
                            q.Corners[j][0] *= QuadDecimate;
                            q.Corners[j][1] *= QuadDecimate;
                        } 
                        else 
                        {
                            q.Corners[j][0] = (q.Corners[j][0] - 0.5f)*QuadDecimate + 0.5f;
                            q.Corners[j][1] = (q.Corners[j][1] - 0.5f)*QuadDecimate + 0.5f;
                        }
                    }
                }
            }

            // if(IsDebug == true)
            // {
            //     Utils.Log.SaveQuadsDataToFile("quads_decimate", quads);
            // }

            QuadAmount = (uint)quads.Count;

            List<Detection> detections = new List<Detection>();

            int chunksize = 1 + quads.Count / (Utils.Calculations.APRILTAG_TASKS_PER_THREAD_TARGET * HowManyThreadsToUse);

            QuickDecode.QuickDecodeTask[] tasks = new QuickDecode.QuickDecodeTask[quads.Count / chunksize + 1];

            int ntasks = 0;
            for (int i = 0; i < quads.Count; i+= chunksize) 
            {
                tasks[ntasks] = new QuickDecode.QuickDecodeTask();
                tasks[ntasks].I0 = i;
                tasks[ntasks].I1 = Utils.Calculations.IMin(quads.Count, i + chunksize);
                tasks[ntasks].Quads = quads;
                tasks[ntasks].ATDetector = this;
                tasks[ntasks].ATImage = im_orig;
                tasks[ntasks].Detections = detections;

                WorkPool.Tasks.Add(tasks[ntasks]);
                ntasks++;
            }

            WorkPool.Run();

            if(IsDebug == true)
            {
                Utils.Log.SaveDetectionsDataToFile("detections_start", detections);
            }

            double[][] poly0 = new double[4][];
            poly0[0] = new double[] {0,0};
            poly0[1] = new double[] {0,0};
            poly0[2] = new double[] {0,0};
            poly0[3] = new double[] {0,0};
            double[][] poly1 = new double[4][];
            poly1[0] = new double[] {0,0};
            poly1[1] = new double[] {0,0};
            poly1[2] = new double[] {0,0};
            poly1[3] = new double[] {0,0};

            for (int i0 = 0; i0 < detections.Count; i0++) 
            {

                Detection det0 = detections[i0];

                for (int k = 0; k < 4; k++)
                {
                    poly0[k][0] = det0.Corners[k][0];
                    poly0[k][1] = det0.Corners[k][1];
                }

                for (int i1 = i0+1; i1 < detections.Count; i1++) 
                {

                    Detection det1 = detections[i1];

                    if (det0.ID != det1.ID || det0.Family.Name != det1.Family.Name)
                    {
                        continue;
                    }

                    for (int k = 0; k < 4; k++)
                    {
                        poly1[k][0] = det1.Corners[k][0];
                        poly1[k][1] = det1.Corners[k][1];
                    }

                    if (Utils.G2D.PolygonOverlapsPolygon(poly0, poly1) != 0) 
                    {
                        // the tags overlap. Delete one, keep the other.

                        int pref = 0; // 0 means undecided which one we'll keep.
                        pref = Utils.Calculations.PreferSmaller(pref, det0.Hamming, det1.Hamming);     // want small hamming
                        pref = Utils.Calculations.PreferSmaller(pref, -det0.DecisionMargin, -det1.DecisionMargin);      // want bigger margins

                        // if we STILL don't prefer one detection over the other, then pick
                        // any deterministic criterion.
                        for (int i = 0; i < 4; i++) 
                        {
                            pref = Utils.Calculations.PreferSmaller(pref, det0.Corners[i][0], det1.Corners[i][0]);
                            pref = Utils.Calculations.PreferSmaller(pref, det0.Corners[i][1], det1.Corners[i][1]);
                        }

                        if (pref < 0) 
                        {
                            // keep det0, destroy det1
                            detections.RemoveAt(i1);
                            i1--; // retry the same index
                            goto retry1;
                        } 
                        else 
                        {
                            // keep det1, destroy det0
                            detections.RemoveAt(i0);
                            i0--; // retry the same index.
                            goto retry0;
                        }
                    }

                retry1: ;
                }

            retry0: ;
            }

            if(cameraParams != null)
            {
                foreach(Detection detection in detections)
                {
                    detection.Pose = new TagPose(detection, cameraParams); 
                }
            }

            if(IsDebug == true)
            {
                Utils.Log.SaveDetectionsDataToFile("detections_filtered", detections);
            }

            output = detections;
            return output;
        }
    }
}
