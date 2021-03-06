#pragma once
////
//  MLtexture - N-dimensional Perlin-style fractal noise
//  D. P. Mitchell
//
//  Cost is exponential in dimension.  Only 2, 3 and 4 dimensional texture have known
//  applications.  At 12 dimensions, it takes one second to evaluate a single point on
//  a 400 MHz Pentium II system.
//
#define MAX_TEXTURE_DIMENSION 8
//
//  Basic Perlin-Noise functions
//
extern float BandNoise(const float rgf[], int nDim, int nPow = 1);
extern float FractalNoise(const float rgf[], int nDim, float fLevelDetail = 5,
                          float fFractalDimension = 1.0f, float fLacunarity = 2.1f,
                          int bAbsoluteValue = 0);
extern float GradientBandNoise(float rgfGrad[], const float rgf[], int nDim, int nPow = 1);
extern float GradientFractalNoise(float rgfGrad[], const float rgf[], int nDim, float fLevelDetail = 5,
                          float fFractalDimension = 1.0f, float fLacunarity = 2.1f,
                          int bAbsoluteValue = 0);
//
//  Periodic texture generation
//
typedef float ML_TextureType(const float [], int, float);

extern float ML_PeriodicTexture(const float rgf[], int nDim, float fBandlimit, ML_TextureType *pt);
//
//  Voronoi textures
//
extern float Voronoi(const float rgf[], int nPow = 1);
extern float GradientVoronoi(float rgfGradient[], const float rgf[], int nPow = 1);