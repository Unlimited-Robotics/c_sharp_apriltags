using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System;
using UnityEngine;

namespace Apriltags
{
    public class Image
    {
        //only one channel in Image!!
        private byte[] _pixels;
        public int Stride { get; private set; }
        public int Height { get; private set; }
        public int Width { get; private set; }
        public byte GetPixel(int x, int y)
        {
            return _pixels[y*Stride + x];
        }
        public byte GetPixel(int x, int y, bool print)
        {
            if(print == true)
            {
                Debug.Log("index " + (y*Stride + x));
            }
            return _pixels[y*Stride + x];
        }
        public byte GetPixelCustom(int index)
        {
            return _pixels[index];
        }
        public void SetPixel(int x, int y, byte value)
        {
            _pixels[y*Stride + x] = value;
        }
        public void SetPixel(int x, int y, byte value, bool print)
        {
            _pixels[y*Stride + x] = value;
            if(print == true)
            {
                Debug.Log("index " + (y*Stride + x));
            }
        }

        private static int getStride(int width, int alignment)
        {
            int stride = width;
            int modulo = stride % alignment;
            if(modulo != 0)
            {
                stride += alignment - modulo;
            }

            return stride;
        }
        
        public Image(Texture2D texture)
        {
            Width = texture.width;
            Height = texture.height;

            Stride = getStride(Width, 96);

            _pixels = new byte[Height * Stride];

            Color[] pixels = texture.GetPixels();

            for (int i = 0; i < texture.width; i++)
            {
                for (int j = 0; j < texture.height; j++)
                {
                    Color pixel = pixels[i + j * texture.width];
                    Color32 pixelBig = pixel;
                    
                    double res = pixelBig.r * 0.299f + pixelBig.g * 0.587f + pixelBig.b * 0.114f;

                    if((texture.height - 1 - j) * texture.width + i == 0)
                    {
                        // Debug.Log("r " + pixelBig.r + ", g " + pixelBig.g + ", b " + pixelBig.b);
                        // Debug.Log("r " + pixelBig.r * 0.299f + ", g " + pixelBig.g * 0.587f + ", b " + pixelBig.b * 0.114f);
                    }

                    _pixels[(texture.height - 1 - j)*Stride + i] = Utils.Calculations.Round(res);
                }
            }
        }

        public Image(int width, int height, int alignment)
        {
            Width = width;
            Height = height;
            Stride = getStride(Width, alignment);
            _pixels = new byte[Height * Stride];
        }

        public Image(int width, int height)
        {
            Width = width;
            Height = height;
            Stride = getStride(Width, 96);
            _pixels = new byte[Height * Stride];
        }

        public Image(Image copy)
        {
            Width = copy.Width;
            Height = copy.Height;
            _pixels = copy._pixels.Clone() as byte[];
            Stride = copy.Stride;
        }

        public void Decimate(float ffactor)
        {

            if(ffactor == 1.5f)
            {
                int swidth = Width / 3 * 2, sheight = Height / 3 * 2;
                int newStride = getStride(swidth, 96);
                byte[]  newPixels = new byte[sheight * newStride];
                int y = 0, sy = 0;
                while(sy < sheight)
                {
                    int x = 0, sx = 0;
                    while(sx < swidth)
                    {
                        byte a = _pixels[y*Stride + x];
                        byte b = _pixels[y*Stride + x + 1];
                        byte c = _pixels[y*Stride + x + 2];

                        byte d = _pixels[(y+1)*Stride + x];
                        byte e = _pixels[(y+1)*Stride + x + 1];
                        byte f = _pixels[(y+1)*Stride + x + 2];
                        
                        byte g = _pixels[(y+2)*Stride + x];
                        byte h = _pixels[(y+2)*Stride + x + 1];
                        byte i = _pixels[(y+2)*Stride + x + 2];

                        // Debug.Log("a: " + a + ", b " + b + ", d " + d + ", e " + e);

                        newPixels[sy*newStride + sx] = (byte)((4*a+2*b+2*d+e)/9);
                        newPixels[sy*newStride + sx + 1] = (byte)((4*c+2*b+2*f+e)/9);

                        newPixels[(sy+1)*newStride + sx] = (byte)((4*g+2*d+2*h+e)/9);
                        newPixels[(sy+1)*newStride +  sx + 1] = (byte)((4*i+2*f+2*h+e)/9);

                        x += 3;
                        sx += 2;
                    }

                    y += 3;
                    sy += 2;
                }

                _pixels = newPixels;
                Width = swidth;
                Height = sheight;
                Stride = newStride;
            }
            else
            {
                int factor = (int)ffactor;
                int swidth = 1 + (Width - 1)/factor;
                int sheight = 1 + (Height - 1)/factor;
                int newStride = getStride(swidth, 96);
                byte[]  newPixels = new byte[sheight * newStride];
                int sy = 0;
                for (int y = 0; y < Height; y+=factor)
                {
                    int sx = 0;
                    for (int x = 0; x < Width; x+=factor)
                    {
                        newPixels[sy*newStride + sx] = _pixels[y*Stride + x];
                        sx++; 
                    }
                    sy++;
                }

                _pixels = newPixels;
                Width = swidth;
                Height = sheight;
                Stride = newStride;
            }
        }

        public void GaussianBlur(float sigma, int ksz)
        {
            // Debug.Log("ksz " + ksz);
            // Debug.Log("sigma " + sigma);

            if(sigma != 0 && (ksz & 1) == 1)
            {
                double[] dk = new double[ksz];

                for (int i = 0; i < ksz; i++) 
                {
                    int x = -ksz/2 + i;
                    double v = Math.Exp(-.5*((x / sigma)*(x / sigma)));
                    dk[i] = v;
                    // Debug.Log("dk before index " + i);
                    // Debug.Log("dk before res " + dk[i]);
                    // Debug.Log("dk before x " + x);
                }

                // normalize
                double acc = 0;
                for (int i = 0; i < ksz; i++)
                {
                    acc += dk[i];
                    // Debug.Log("acc index " + i);
                    // Debug.Log("acc res " + acc);
                }

                for (int i = 0; i < ksz; i++)
                {
                    dk[i] /= acc;
                    // Debug.Log("dk before index " + i);
                    // Debug.Log("dk before res " + dk[i]);
                }

                byte[] k = new byte[ksz];
                for (int i = 0; i < ksz; i++)
                {
                    k[i] = (byte)(dk[i]*255);
                    // Debug.Log("k index " + i);
                    // Debug.Log("k res " + k[i]);
                }

                convolve2D(k, ksz);
            }
        }

        private void convolve2D(byte[] k, int ksz)
        {
            if((ksz & 1) == 1)
            {
                for (int y = 0; y < Height; y++) 
                {
                    byte[] x = new byte[Stride];
                    byte[] tmp = new byte[Stride];
                    for (int i = 0; i < Stride; i++)
                    {
                        x[i] = _pixels[y*Stride + i];
                    }

                    Utils.Calculations.Convolve(x, tmp, Width, k, ksz);

                    for (int i = 0; i < Stride; i++)
                    {
                        _pixels[y*Stride + i] = tmp[i];
                    }
                }

                for (int x = 0; x < Width; x++) 
                {
                    byte[] xb = new byte[Height];
                    byte[] yb = new byte[Height];

                    for (int y = 0; y < Height; y++)
                    {
                        xb[y] = _pixels[y*Stride + x];
                    }

                    Utils.Calculations.Convolve(xb, yb, Height, k, ksz);

                    for (int y = 0; y < Height; y++)
                    {
                        _pixels[y*Stride + x] = yb[y];
                    }
                }
            }
        }

        public double ValueForPixel(double px, double py) 
        {
            int x1 = (int)Math.Floor(px - 0.5);
            int x2 = (int)Math.Ceiling(px - 0.5);
            double x = px - 0.5 - x1;
            int y1 = (int)Math.Floor(py - 0.5);
            int y2 = (int)Math.Ceiling(py - 0.5);
            double y = py - 0.5 - y1;
            if (x1 < 0 || x2 >= Width || y1 < 0 || y2 >= Height) 
            {
                return -1;
            }
            return _pixels[y1*Stride + x1]*(1-x)*(1-y) +
                    _pixels[y1*Stride + x2]*x*(1-y) +
                    _pixels[y2*Stride + x1]*(1-x)*y +
                    _pixels[y2*Stride + x2]*x*y;
        }
    }
}
