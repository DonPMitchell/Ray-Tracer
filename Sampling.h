//
//  Build a digital image by sampling a function f(x,y,t)
//  D.P. Mitchell  2020/12/18
//
#pragma once

//
//  Simple uniform sampling for testing
//
struct TestSamplingTile : SamplingTile {
            TestSamplingTile() {}
    int     Initialize();
    int     Render(Image &im, Radiance (*f)(Vector2 v), double fSamplesPerPixel);
};
//
//  Cast one ray in the center of the scene
//
struct DebugSamplingTile : SamplingTile {
            DebugSamplingTile() {}
    int     Initialize();
    int     Render(Image &im, Radiance (*f)(Vector2 v), double fSamplesPerPixel);
};
//
//  Simple jittered super sampling tile
//
struct JitterSamplingTile : SamplingTile {
            JitterSamplingTile() {}
    int     Initialize();
    int     Render(Image &im, Radiance (*f)(Vector2 v), double fSamplesPerPixel);
};
//
//  Blue Noise sampling tile
//

extern Vector2 g_rgvBlueNoise[BLUE_NSQRT*BLUE_NSQRT];

struct BlueNoiseSamplingTile : SamplingTile {
            BlueNoiseSamplingTile() {}
    int     Initialize();
    int     Render(Image &im, Radiance (*f)(Vector2 v), double fSamplesPerPixel);
};

struct FastBlueSamplingTile : SamplingTile {
            FastBlueSamplingTile() { nVenera = 0;}
    int     Initialize();
    int     Render(Image &im, Radiance (*f)(Vector2 v), double fSamplesPerPixel);
    void    QIndex1(int nPrime);
};

struct AdaptiveBlueSamplingTile : SamplingTile {
            AdaptiveBlueSamplingTile() { nVenera = 0;}
    int     Initialize();
    int     Render(Image &im, Radiance (*f)(Vector2 v), double fSamplesPerPixel);
    void    QIndex1(int nPrime);
};

#define TESTIMAGESIZE   512

extern Radiance ZoneFunction(Vector2 v);
extern Radiance FactoryFunction(Vector2 v);
extern void     PerfectFactory(Image &im);
extern void     PerfectZone(Image &im);

extern Vector2 DiskSamples16[16];
extern Vector2 DiskSamples32[32];
extern Vector2 DiskSamples64[64];
extern Vector2 DiskSamples128[128];

extern int rgnQuasiIndex2[128];
extern int rgnQuasiIndex3[128];
extern int rgnQuasiIndex5[128];
extern int rgnQuasiIndex7[128];

extern Vector2  g_rgvBlueNoise[BLUE_NSQRT*BLUE_NSQRT];
extern int      g_BlueNoiseQuasiIndex2[BLUE_NSQRT*BLUE_NSQRT];
extern int      g_BlueNoiseQuasiIndex3[BLUE_NSQRT*BLUE_NSQRT];
extern int      g_BlueNoiseQuasiIndex5[BLUE_NSQRT*BLUE_NSQRT];
extern int      g_BlueNoiseQuasiIndex7[BLUE_NSQRT*BLUE_NSQRT];