// RayTracer2020.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "RayTracer2020.h"
#include "Cameras.h"
#include "Lights.h"
#include "Patterns.h"
#include "Sampling.h"
#include "Solids.h"
#include "RadiometryColor.h"
#include "Video.h"

int g_nVerbose = 0;
//
//  Ray tracer main routines
//
static RayTracer *s_prt;

static Radiance
CastFunction(Vector2 v)
{
    Ray r;
    Radiance R;

    r = s_prt->pcCamera->Emit(v);
    s_prt->nRayLevel = 0;
    s_prt->fAttenuation = 1.0;
    R = s_prt->VisualRay(r);
    return R;
}

int
RayTracer::Render(Image &im, double fSamplesPerPixel)
{
    s_prt = this;                                       // This is not a great way to create the image function
    return psSample->Render(im, CastFunction, fSamplesPerPixel);
}

extern void TestSampling();
extern void TestImageFunctions();
extern void TestPatterns();
extern void TestSolids();
extern void PatternTestImage(char *, Vector2 rgv[], int);
extern void TestCamera();
extern void TestSolidRender();
extern void BuildSpectrumRGB();
extern void TestMaterial();
extern void TestCubeRoot();
extern void TestRayTraceScene();
extern void TestFresnel();
extern void TestSpectrum();
extern void TestHitLists();
extern void TestRefraction();
extern void MakeBlueNoiseTile();
extern void TestExtents();
extern void TestDiskDiscrepancy();
extern void TestOptimize();
extern void TestLightRayBug();
extern void BuildRGBTables();
extern void TestBandNoise();
extern void TestSampleGradient();
extern void TestRGBSpectrum();
extern void BuildSphericalMaps();
extern void BuildVenusMaps();
extern void BuildMoonColorMap();
extern void Laboratory();
extern void TestClipPlanes();
extern void TranslateSpectra();
extern void Venera();
extern void TestObjFiles();
extern void TestRayleigh();
extern void TestVRGB();
extern void TestVoronoi();
extern void TestDiscrepancy();
extern void DistributedRayTracer();
extern void TestAreaMaps();
extern void GeneratePhiValues();
extern void TestQuadirandomPatterns();
extern void TestNewLocality();
extern void BuildSpectraXYZ();
extern void TestGeneralHilbert();
extern void Megacycles();


int _tmain(int argc, _TCHAR* argv[])
{

    // Megacycles(); return 0;
    //MakeBlueNoiseTile(); return 0;
    //TestGeneralHilbert(); return 0;
    //BuildSpectraXYZ(); return 0;
    // TestNewLocality(); return 0;
    // TestRayleigh();
    // Laboratory();
    Venera(); return 0;
    TestRayTraceScene();
    // DistributedRayTracer();
    return 0;
}

