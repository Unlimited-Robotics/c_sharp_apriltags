using System.Collections;
using System.Collections.Generic;

namespace Apriltags
{
    public static class QuadThresh
    {
        public static List<Quad> GetQuadThresh(Detector detector, Image image, bool isDebug)
        {
            int w = image.Width, h = image.Height;
            Image threshImage = getThresholdImage(detector, image);
            if(isDebug == true)
            {
                Utils.Log.SaveImageDataToFile("image_thresh", threshImage);
            }

            int ts = threshImage.Stride;

        ////////////////////////////////////////////////////////
            // step 2. find connected components.
            UnionFind uf = new UnionFind(detector, threshImage, w, h);
            if(isDebug == true)
            {
                Utils.Log.SaveUnionfindDataToFile("unionfind_start", uf);
            }

            List<Cluster> clusters = ClusterMap.GradientClusters(detector, threshImage, w, h, ts, uf);
            if(isDebug == true)
            {
                Utils.Log.SaveClustersDataToFile("clusters_start", clusters);
            }

        
            List<Quad> quads = fitQuads(detector, w, h, clusters, image);

        // if (td->debug) {
        //     FILE *f = fopen("debug_lines.ps", "w");
        //     fprintf(f, "%%!PS\n\n");

        //     image_u8_t *im2 = image_u8_copy(im);
        //     image_u8_darken(im2);
        //     image_u8_darken(im2);

        //     // assume letter, which is 612x792 points.
        //     double scale = fmin(612.0/im->width, 792.0/im2->height);
        //     fprintf(f, "%.15f %.15f scale\n", scale, scale);
        //     fprintf(f, "0 %d translate\n", im2->height);
        //     fprintf(f, "1 -1 scale\n");

        //     postscript_image(f, im2);

        //     image_u8_destroy(im2);

        //     for (int i = 0; i < zarray_size(quads); i++) {
        //         struct quad *q;
        //         zarray_get_volatile(quads, i, &q);

        //         float rgb[3];
        //         int bias = 100;

        //         for (int i = 0; i < 3; i++)
        //             rgb[i] = bias + (random() % (255-bias));

        //         fprintf(f, "%f %f %f setrgbcolor\n", rgb[0]/255.0f, rgb[1]/255.0f, rgb[2]/255.0f);
        //         fprintf(f, "%.15f %.15f moveto %.15f %.15f lineto %.15f %.15f lineto %.15f %.15f lineto %.15f %.15f lineto stroke\n",
        //                 q->p[0][0], q->p[0][1],
        //                 q->p[1][0], q->p[1][1],
        //                 q->p[2][0], q->p[2][1],
        //                 q->p[3][0], q->p[3][1],
        //                 q->p[0][0], q->p[0][1]);
        //     }

        //     fclose(f);
        // }

            return quads;
        }

        private static Image getThresholdImage(Detector detector, Image image)
        {
            int w = image.Width, h = image.Height, s = image.Stride;
            Image output = new Image(w, h, s);

            const int tileSize = 4;

            int tw = w / tileSize;
            int th = h / tileSize;

            byte[] im_max = new byte[tw*th];
            byte[] im_min = new byte[tw*th];

            for (int ty = 0; ty < th; ty++)
            {
                for (int tx = 0; tx < tw; tx++)
                {
                    byte max = 0, min = 255;
                    for (int dy = 0; dy < tileSize; dy++)
                    {
                        for (int dx = 0; dx < tileSize; dx++)
                        {
                            byte v = image.GetPixelCustom((ty*tileSize+dy)*s + tx *tileSize + dx);
                            if(v < min)
                            {
                                min = v;
                            }
                            if(v > max)
                            {
                                max = v;
                            }
                        }
                    }

                    // Debug.Log("first step index " + (ty*tw + tx));
                    // Debug.Log("first step max " + max);
                    // Debug.Log("first step min " + min);
                    im_max[ty*tw + tx] = max;
                    im_min[ty*tw + tx] = min;
                }
            }

            byte[] im_max_tmp = new byte[tw*th];
            byte[] im_min_tmp = new byte[tw*th];

            for (int ty = 0; ty < th; ty++) 
            {
                for (int tx = 0; tx < tw; tx++) 
                {
                    byte max = 0, min = 255;

                    for (int dy = -1; dy <= 1; dy++) 
                    {
                        if (ty+dy < 0 || ty+dy >= th)
                        {
                            continue;
                        }
                        for (int dx = -1; dx <= 1; dx++) 
                        {
                            if (tx+dx < 0 || tx+dx >= tw)
                            {
                                continue;
                            }

                            byte m = im_max[(ty+dy)*tw+tx+dx];
                            if (m > max)
                            {
                                max = m;
                            }
                            m = im_min[(ty+dy)*tw+tx+dx];
                            if (m < min)
                            {
                                min = m;
                            }
                        }
                    }

                    im_max_tmp[ty*tw + tx] = max;
                    im_min_tmp[ty*tw + tx] = min;
                }
            }
            im_max = im_max_tmp;
            im_min = im_min_tmp;

            for (int ty = 0; ty < th; ty++) 
            {
                for (int tx = 0; tx < tw; tx++) 
                {
                    int min = im_min[ty*tw + tx];
                    int max = im_max[ty*tw + tx];

                    // Debug.Log("param " + detector.QuadThreshParams.MinWhiteBlackDiff);

                    // low contrast region? (no edges)
                    //     Debug.Log("second step outside index " + (ty*tw + tx));
                    // Debug.Log("max " + max);
                    // Debug.Log("min " + min);
                    // Debug.Log("diff " + (max - min));
                    // Debug.Log("tx " + tx);
                    // Debug.Log("ty " + ty);
                    if (max - min < detector.QuadThreshParams.MinWhiteBlackDiff) 
                    {
                        for (int dy = 0; dy < tileSize; dy++) 
                        {
                            int y = ty*tileSize + dy;

                            for (int dx = 0; dx < tileSize; dx++) 
                            {
                                int x = tx*tileSize + dx;

                                // Debug.Log("second step res 127");
                                output.SetPixel(x, y, 127, false);
                            }
                        }
                        // Debug.Log("cont");
                        continue;
                    }

                    // otherwise, actually threshold this tile.

                    // argument for biasing towards dark; specular highlights
                    // can be substantially brighter than white tag parts
                    byte thresh = (byte)(min + (max - min) / 2);

                    for (int dy = 0; dy < tileSize; dy++) 
                    {
                        int y = ty*tileSize + dy;

                        for (int dx = 0; dx < tileSize; dx++) 
                        {
                            int x = tx*tileSize + dx;

                            // Debug.Log("third step");
                            byte v = image.GetPixel(x, y, false);
                            if (v > thresh)
                            {
                                // Debug.Log("third step res 255");
                                output.SetPixel(x, y, 255);
                            }
                            else
                            {
                                // Debug.Log("third step res 0");
                                output.SetPixel(x, y, 0);
                            }
                        }
                    }
                }
            }

            for (int y = 0; y < h; y++) {

                // what is the first x coordinate we need to process in this row?

                int x0;

                if (y >= th*tileSize) 
                {
                    x0 = 0; // we're at the bottom; do the whole row.
                }
                else 
                {
                    x0 = tw*tileSize; // we only need to do the right most part.
                }

                // compute tile coordinates and clamp.
                int ty = y / tileSize;
                if (ty >= th)
                {
                    ty = th - 1;
                }

                for (int x = x0; x < w; x++) 
                {
                    int tx = x / tileSize;
                    if (tx >= tw)
                    {
                        tx = tw - 1;
                    }

                    int max = im_max[ty*tw + tx];
                    int min = im_min[ty*tw + tx];
                    int thresh = min + (max - min) / 2;

                    byte v = image.GetPixel(x,y);
                    if (v > thresh)
                    {
                        output.SetPixel(x, y, 255);
                    }
                    else
                    {
                        output.SetPixel(x, y, 0);
                    }
                }
            }

            if (detector.QuadThreshParams.Deglitch != 0) 
            {
                Image tmp = new Image(w, h);

                for (int y = 1; y + 1 < h; y++) 
                {
                    for (int x = 1; x + 1 < w; x++) 
                    {
                        byte max = 0;
                        for (int dy = -1; dy <= 1; dy++) 
                        {
                            for (int dx = -1; dx <= 1; dx++) 
                            {
                                byte v = output.GetPixelCustom((y+dy)*s + x + dx);
                                if (v > max)
                                {
                                    max = v;
                                }
                            }
                        }
                        tmp.SetPixel(x, y, max);
                    }
                }

                for (int y = 1; y + 1 < h; y++) 
                {
                    for (int x = 1; x + 1 < w; x++) 
                    {
                        byte min = 255;
                        for (int dy = -1; dy <= 1; dy++) 
                        {
                            for (int dx = -1; dx <= 1; dx++) 
                            {
                                byte v = tmp.GetPixelCustom((y+dy)*s + x + dx);
                                if (v < min)
                                {
                                    min = v;
                                }
                            }
                        }
                        output.SetPixel(x, y, min);
                    }
                }
            }

            return output;
        }

        private static List<Quad> fitQuads(Detector detector, int w, int h, List<Cluster> clusters, Image image) 
        {
            List<Quad> quads = new List<Quad>();

            bool normal_border = false;
            bool reversed_border = false;
            int min_tag_width = 1000000;
            for (int i = 0; i < detector.TagFamilies.Count; i++) {
                ApriltagFamily family = detector.TagFamilies[i];
                if (family.WidthAtBorder < min_tag_width) 
                {
                    min_tag_width = family.WidthAtBorder;
                }
                normal_border |= !family.ReversedBorder;
                reversed_border |= family.ReversedBorder;
            }
            // Debug.Log("min_tag_width " + min_tag_width);
            // Debug.Log("QuadDecimate " + detector.QuadDecimate);
            min_tag_width = (int)(min_tag_width / detector.QuadDecimate);
            if (min_tag_width < 3) 
            {
                min_tag_width = 3;
            }

            int sz = clusters.Count;
            int chunksize = 1 + sz / (Utils.Calculations.APRILTAG_TASKS_PER_THREAD_TARGET * detector.HowManyThreadsToUse);
            Quad.QuadTask[] tasks = new Quad.QuadTask[sz / chunksize + 1];

            // Debug.Log("normal border " + normal_border);
            // Debug.Log("reversed border " + reversed_border);
            // Debug.Log("min tag width " + min_tag_width);
            // Debug.Log("sz " + sz);
            // Debug.Log("chunksize " + chunksize);

            int ntasks = 0;
            for (int i = 0; i < sz; i += chunksize) 
            {
                // Debug.Log("ntasks " + ntasks);
                tasks[ntasks] = new Quad.QuadTask();
                tasks[ntasks].ATDetector = detector;
                tasks[ntasks].Cidx0 = i;
                // Debug.Log("Cidx0 " + tasks[ntasks].Cidx0);
                tasks[ntasks].Cidx1 = Utils.Calculations.IMin(sz, i + chunksize);
                // Debug.Log("Cidx1 " + tasks[ntasks].Cidx1);
                tasks[ntasks].Height = h;
                // Debug.Log("Height " + tasks[ntasks].Height);
                tasks[ntasks].Width = w;
                // Debug.Log("Width " + tasks[ntasks].Width);
                tasks[ntasks].Quads = quads;
                tasks[ntasks].Clusters = clusters;
                tasks[ntasks].ATImage = image;
                tasks[ntasks].TagWidth = min_tag_width;
                tasks[ntasks].NormalBorder = normal_border;
                tasks[ntasks].ReversedBorder = reversed_border;

                detector.WorkPool.Tasks.Add(tasks[ntasks]);
                ntasks++;
            }

            detector.WorkPool.Run();

            return quads;
        }
    }
}
