using System.Collections;
using System.Collections.Generic;
using System;
using System.IO;

namespace Apriltags.Utils
{
    public class Test
    {
        public static void ValidateImageDataFiles()
        {
            DirectoryInfo folder = new DirectoryInfo("/home/ros2/_Alon/Compares/Apriltags/unity");

            foreach(FileInfo file in folder.GetFiles())
            {
                FileInfo pupil = new FileInfo("/home/ros2/_Alon/Compares/Apriltags/pupil/" + file.Name);
                if(pupil == null)
                {
                    Debug.LogError("no pupil image data counterpart to: " + file.Name);
                }
                else
                {
                    bool areSameFile = true;
                    int lineCounter = 0;
                    List<int> disparityLineNumber = new List<int>();
                    List<string> lineUnityDisparity = new List<string>();
                    List<string> linePupilDisparity = new List<string>();
                    using(StreamReader srUnity = file.OpenText())
                    {
                        using(StreamReader srPupil = pupil.OpenText())
                        {
                            string lineUnity = srUnity.ReadLine(); 
                            string linePupil = srPupil.ReadLine();
                            while(lineUnity != null)
                            {
                                lineCounter++;
                                if(linePupil != lineUnity)
                                {
                                    areSameFile = false;
                                    disparityLineNumber.Add(lineCounter);
                                    lineUnityDisparity.Add(lineUnity);
                                    linePupilDisparity.Add(linePupil);
                                }
                                lineUnity = srUnity.ReadLine(); 
                                linePupil = srPupil.ReadLine();
                            }

                            if(linePupil != null)
                            {
                                areSameFile = false;
                            }
                        }
                    }

                    if(areSameFile == true)
                    {
                        Debug.Log(file.Name + " validated");
                    }
                    else
                    {
                        Debug.LogError(file.Name + " wrong data, disparity amount: " + disparityLineNumber.Count);
                        for (int i = 0; i < disparityLineNumber.Count; i++)
                        {
                            Debug.LogError("disparity at line " + disparityLineNumber[i].ToString());
                            if(lineUnityDisparity[i] != null)
                            {
                                Debug.LogError("unity line:\n" + lineUnityDisparity[i]);
                            }
                            else
                            {
                                Debug.LogError("unity EOF");
                            }
                            if(linePupilDisparity[i] != null)
                            {
                                Debug.LogError("pupil line:\n" + linePupilDisparity[i]);
                            }
                            else
                            {
                                Debug.LogError("pupil EOF");
                            }
                        }
                    }
                }
            }
        }
    }
}
