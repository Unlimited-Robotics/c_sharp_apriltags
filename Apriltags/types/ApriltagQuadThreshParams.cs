using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace Apriltags
{
    public class QuadThreshParams
    {
        // reject quads containing too few pixels
        public int MinClusterPixels;

        // how many corner candidates to consider when segmenting a group
        // of pixels into a quad.
        public int MaxNMaxima;

        // Reject quads where pairs of edges have angles that are close to
        // straight or close to 180 degrees. Zero means that no quads are
        // rejected. (In radians).
        public float CriticalRad;
        public float CosCriticalRad;

        // When fitting lines to the contours, what is the maximum mean
        // squared error allowed?  This is useful in rejecting contours
        // that are far from being quad shaped; rejecting these quads "early"
        // saves expensive decoding processing.
        public float MaxLineFitMse;

        // When we build our model of black & white pixels, we add an
        // extra check that the white model must be (overall) brighter
        // than the black model.  How much brighter? (in pixel values,
        // [0,255]). .
        public int MinWhiteBlackDiff;

        // should the thresholded image be deglitched? Only useful for
        // very noisy images
        public int Deglitch;

        public QuadThreshParams()
        {
            MaxNMaxima = 10;
            MinClusterPixels = 5;
            MaxLineFitMse = 10;
            CosCriticalRad = Mathf.Cos(10 * Mathf.PI / 180);
            Deglitch = 0;
            MinWhiteBlackDiff = 5;
        }
    }
}