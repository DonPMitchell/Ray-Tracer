#include "stdafx.h"
#include "RayTracer2020.h"
#include "Sampling.h"
#include "MonteCarlo.h"

#pragma intrinsic(sin)

#define NSQRT 16
#define NTILE (NSQRT*NSQRT)
#define ML_KAISER4_HALFWIDTH            4.0
#define D_PI    3.14159265358979323846264338327950288419716939937510

Vector2 rgvN[NTILE];                // nonuniform samples
Vector2 rgvU[NTILE];                // uniform samples
double rgfU2N[NTILE][NTILE];
double rgfN2U[NTILE][NTILE];

//
//  Kaiser Window and Sinc
//
static double gs_fI0Alpha = 1.0/11.3019219521363;
static double gs_fAlphaSquared = 4.0*4.0;

static double
I0(double x2)       // Modified Bessel function of x**2, used by Kaiser
{
	int n;
	double u, s, x;

	x = 0.25 * x2;
	u = 1.0;
	s = 1.0;
	for (n = 1; u/s > 1.0e-10; n++) {
		u *= x / (double)(n*n);
		s += u;
	}
	return s;
}

double
KaiserWindow(double x)
{
    double x2;

	x2 = x*x;
	if (x2 > 1.0) return 0.0;
	return I0(gs_fAlphaSquared*(1.0-x2))*gs_fI0Alpha;
}

static double
Sinc(double x)
{
    x *= D_PI;
    if (x < 0.01 && x > -0.01) {
        x = x*x;
        return 1.0 + x*(-1.0/6.0 + x*(1.0/120.0 + x*(-1.0/5040.0 + x/362880.0)));       // accurate out to about x == 0.2
    } else
        return sin(x)/x;
}

inline double
Kaiser4Filter(double x)
{
    return Sinc(x) * KaiserWindow(x/ML_KAISER4_HALFWIDTH);
}

//
//  Prepare perfectly filtered test images for comparison
//
//
//  caching the filter coefficients speeds up ML_ImageResample by orders of magnitude
//
#define ML_SAMPLE_TO_INDEX(X)   int(floor(X))
#define ML_INDEX_TO_SAMPLE(I)   (double(I) + 0.5)
#define ML_LOWEST_INDEX(X)      int(floor(X + 0.5))
#define ML_HIGHEST_INDEX(X)     int(floor(X - 0.5))

typedef double   (*ML_PFILTER)(double);

static int s_nLastSrc;
static int s_nLastDst;
static ML_PFILTER s_pLastFilter = 0;
static double *s_pfTap = 0;
static int s_nTaps = 0;

static void
ComputeTaps(int nDst, int nSrc, double xHalfWidth, ML_PFILTER pFilter)
{
    double xSrc, x, fRatio, fDecimation, fScale;
    int iDst, iSrc, iFirst, iLast, nTaps;
    double *pf;

    fRatio = double(nSrc)/double(nDst);
    if (nDst < nSrc)
        fDecimation = fRatio;
    else
        fDecimation = 1.0;
    xHalfWidth *= fDecimation;
    fScale = 1.0/fDecimation;
    nTaps = nDst * int(2.0*xHalfWidth + 1.0);
    if (s_pfTap == 0 || s_nTaps < nTaps) {
        delete [] s_pfTap;
        s_pfTap = new double[nTaps];
        s_nTaps = nTaps;
    }
    pf = s_pfTap;
    for (iDst = 0; iDst < nDst; iDst++) {
        xSrc = ML_INDEX_TO_SAMPLE(iDst)*fRatio;
        iFirst = ML_LOWEST_INDEX(xSrc - xHalfWidth);
        iLast = ML_HIGHEST_INDEX(xSrc + xHalfWidth);
        for (iSrc = iFirst; iSrc <= iLast; iSrc++) {
            x = ML_INDEX_TO_SAMPLE(iSrc);
            *pf++ = pFilter(float((x - xSrc) * fScale));

        }
    }
    s_nLastDst = nDst;
    s_nLastSrc = nSrc;
    s_pLastFilter = pFilter;
}

static inline int
ML_Modulus(int n, int m)
{
    n = n % m;
    if (n < 0)
        n += m;
    return n;
}

static int
Resample(float rgfDst[], int nDst, const float rgfSrc[], int nSrc, int nWrap,
            double xHalfWidth, ML_PFILTER pFilter)
{
    double xSrc, y, fRatio, fDecimation, fScale, f, fSum, fWeight;
    int iDst, i, iSrc, iFirst, iLast;
    double *pf;

    if (rgfDst == rgfSrc || nSrc == 0 || nDst == 0)
        return 0;
    if (nDst != s_nLastDst || nSrc != s_nLastSrc || pFilter != s_pLastFilter)
        ComputeTaps(nDst, nSrc, xHalfWidth, pFilter);
    pf = s_pfTap;
    fRatio = float(nSrc)/float(nDst);
    if (nDst < nSrc)
        fDecimation = fRatio;
    else
        fDecimation = 1.0;
    xHalfWidth *= fDecimation;
    fScale = 1.0/fDecimation;
    for (iDst = 0; iDst < nDst; iDst++) {
        xSrc = ML_INDEX_TO_SAMPLE(iDst)*fRatio;
        iFirst = ML_LOWEST_INDEX(xSrc - xHalfWidth);
        iLast = ML_HIGHEST_INDEX(xSrc + xHalfWidth);
        fSum = fWeight = 0.0;
        for (iSrc = iFirst; iSrc <= iLast; iSrc++) {
            //x = double(iSrc) + 0.5;
            if (unsigned(iSrc) >= unsigned(nSrc)) {
                if (nWrap == WRAP_AROUND) {
                    i = ML_Modulus(iSrc, nSrc);
                } else if (nWrap == WRAP_MIRROR) {
                    i = ML_Modulus(iSrc, 2*nSrc);
                    if (i >= nSrc)
                        i = 2*nSrc - i - 1;
                } else {
                    if (iSrc < 0)
                        i = 0;
                    else if (iSrc >= nSrc)
                        i = nSrc - 1;
                }
            } else
                i = iSrc;
            y = rgfSrc[i];
            //f = pFilter(float((x - xSrc) * fScale));
            f = *pf++;
            fSum += f*y;
            fWeight += f;
        }
        if (fWeight)
            fSum /= fWeight;
        rgfDst[iDst] = float(fSum);
    }
    return 1;
}

static int
ImageResample(Image &imgNew, int nNewWidth, int nNewHeight, const Image &imgSrc,
                 float xHalfWidth, ML_PFILTER pFilter)
{
    int i, j, k, nOldWidth, nOldHeight;
    float *pfDst, *pfRow, *pfCol, *pfTransposed;
    const float *pfSrc;

    nOldWidth = imgSrc.m_nWidth;
    nOldHeight = imgSrc.m_nHeight;
    pfRow = new float[nNewWidth];
    if (pfRow == 0)
        return 0;
    pfCol = new float[nNewHeight];
    if (pfCol == 0)
        return 0;
    pfTransposed = new float[nOldHeight * nNewWidth];
    if (pfTransposed == 0)
        return 0;
    //
    //  Note: important to transfer the boundary condition from src to dst.
    //
    if (imgNew.NewImage(nNewWidth, nNewHeight, imgSrc.m_nChannels, imgSrc.m_bWrapX, imgSrc.m_bWrapY) == 0)
        return 0;
    for (k = 0; k < imgSrc.m_nChannels; k++) {
        //
        //  Resample horizontally and transpose
        //
        pfSrc = imgSrc.m_pfImage + k*nOldWidth*nOldHeight;
        for (j = 0; j < nOldHeight; j++) {
            Resample(pfRow, nNewWidth, pfSrc + j*nOldWidth, nOldWidth, imgSrc.m_bWrapX, xHalfWidth, pFilter);
            for (i = 0; i < nNewWidth; i++)
                pfTransposed[j + i*nOldHeight] = pfRow[i];
        }
        //
        //  Resample vertically
        //
        pfDst = imgNew.m_pfImage + k*nNewWidth*nNewHeight;
        for (j = 0; j < nNewWidth; j++) {
            Resample(pfCol, nNewHeight, pfTransposed + j*nOldHeight, nOldHeight, imgSrc.m_bWrapY, xHalfWidth, pFilter);
            for (i = 0; i < nNewHeight; i++)
                pfDst[j + i*nNewWidth] = pfCol[i];
        }
    }
    delete [] pfRow;
    delete [] pfCol;
    delete [] pfTransposed;
    return 1;
}
//
//  Test image functions
//
double
ImageZone(double x, double y)
{
    x = x - 0.5*double(TESTIMAGESIZE);
    y = y - 0.5*double(TESTIMAGESIZE);
    x /= double(TESTIMAGESIZE)/D_SQRT2;
    y /= double(TESTIMAGESIZE)/D_SQRT2;
    return 0.5*cos(D_PI*TESTIMAGESIZE*(x*x + y*y)) + 0.5;
}

Radiance
ZoneFunction(Vector2 v)
{
    double x, y;
    Radiance rad;

    x = (v.x - 0.5)*double(TESTIMAGESIZE);
    y = (v.y - 0.5)*double(TESTIMAGESIZE);
    x /= double(TESTIMAGESIZE)/D_SQRT2;
    y /= double(TESTIMAGESIZE)/D_SQRT2;
    rad = RGBtoSpectrum(0.5*cos(D_PI*TESTIMAGESIZE*(x*x + y*y)) + 0.5);
    return rad;
}

Radiance
FactoryFunction(Vector2 v)
{
    static Image im;
    static int bLoaded;
    double x, y;
    Radiance rad;

    if (!bLoaded) {
        if (im.ReadBMP("ImageFactory.bmp") == 0) {
            printf("Cannot open image\n");
            exit(0);
        }
        bLoaded = 1;
    }
    x = im.m_nWidth*v.x;
    y = im.m_nHeight*v.y;
    rad = RGBtoSpectrum(im.Sample(x, y));
    return rad;
}

void
PerfectFactory(Image &im)
{
    Image im2;

    im2.ReadBMP("ImageFactory.bmp");
    ImageResample(im, TESTIMAGESIZE, TESTIMAGESIZE, im2, ML_KAISER4_HALFWIDTH, Kaiser4Filter);
}

void
PerfectZone(Image &im)
{
    Image im2;
    int i, j;
    double f;

    im2.NewImage(TESTIMAGESIZE*2, TESTIMAGESIZE*2, 1);
    for (j = 0; j < TESTIMAGESIZE*2; j++) {
        for (i = 0; i < TESTIMAGESIZE*2; i++) {
            f = ImageZone(0.5*(double(i) + 0.5), 0.5*(double(j) + 0.5));
            im2.Set(f, i, j);
        }
    }
    // im2.WriteBMP("BigZone.bmp");
    ImageResample(im, TESTIMAGESIZE, TESTIMAGESIZE, im2, ML_KAISER4_HALFWIDTH, Kaiser4Filter);
}

//
//  Render method
//
//  See work in Sampling(2018) and Antialiasing(Project 2013) and Antialiasing(Project)
//  This is where forking multiple threads can happen
//
extern int g_nTestRow, g_nTestColumn, g_nTestFrame;

int
TestSamplingTile::Initialize()
{
    return nSamples;
}

int     
TestSamplingTile::Render(Image &im, Radiance (*f)(Vector2 v), double fSamplesPerPixel)
{
    int i, j, n, nHi, nSamples;
    Vector2 v;
    Radiance rad;
    DisplayRGB rgb;

    n = im.m_nWidth;
    nHi = im.m_nHeight;
    nSamples = 0;
    for (j = 0; j < nHi; j++) {
        g_nTestRow = j;
        v.y = (double(j) + (n - nHi)/2 + 0.5)/double(n);
        for (i = 0; i < n; i++) {
            g_nTestColumn = i;
            v.x = (double(i) + 0.5)/double(n);
            rad = f(v);
            nSamples++;
            rgb = rad.sRGB();
            im.SetRGB(rgb, i, j);
        }
    }
    return nSamples;
}

int
DebugSamplingTile::Initialize()
{
    return nSamples;
}

int     
DebugSamplingTile::Render(Image &im, Radiance (*f)(Vector2 v), double fSamplesPerPixel)
{
    Vector2 v;
    Radiance rad;
    DisplayRGB rgb;

    printf("\nDebug Sample at 0.5, 0.5:\n");
    g_nVerbose = 1;
    v = Vector2(0.5, 0.5);
    rad = f(v);
    rgb = rad.sRGB();
    rgb.Print();
    im.FillRGB(rgb);
    printf("\n");

    return 1;

    v = Vector2(0.5 + 0.005, 0.5);
    rad = f(v);
    rgb = rad.sRGB();
    rgb.Print();
    printf("\n");

    v = Vector2(0.5 - 0.005, 0.5);
    rad = f(v);
    rgb = rad.sRGB();
    rgb.Print();
    printf("\n");
    return 1;
}

int
JitterSamplingTile::Initialize()
{
    return nSamples;
}

int     
JitterSamplingTile::Render(Image &im, Radiance (*f)(Vector2 v), double fSamplesPerPixel)
{
    int i, j, n, nHi, nSamples;
    Vector2 v;
    Radiance rad;
    Image im2;

    n = int(im.m_nWidth * sqrt(fSamplesPerPixel) + 0.5);
    nHi = int(im.m_nHeight * sqrt(fSamplesPerPixel) + 0.5);
    im2.NewImage(n, nHi, 3);        // rendered image
    im2.FillRGB(ML_Black);
    nSamples = 0;
    for (j = 0; j < nHi; j++) {
        for (i = 0; i < n; i++) {
            v.y = (double(j) + (n - nHi)/2 + RandomDouble())/double(n);
            v.x = (double(i) + RandomDouble())/double(n);
            rad = f(v);
            nSamples++;
            im2.SetRGB(rad.sRGB(), i, j);
        }
    }
    //im2.WriteBMP("RenderJitter.bmp");
    ImageResample(im, im.m_nWidth, im.m_nHeight, im2,  ML_KAISER4_HALFWIDTH, Kaiser4Filter);
    return nSamples;
}

int
BlueNoiseSamplingTile::Initialize()
{
    return nSamples;
};

int
BlueNoiseSamplingTile::Render(Image &im, Radiance (*f)(Vector2 v), double fSamplesPerPixel)
{
    Image imWeight;
    Vector2 vTile, vAdjust, v;
    Radiance rad;
    double fSqrt, fSize;
    int i, nSamples;

    imWeight.NewImage(im.m_nWidth, im.m_nHeight, 1);
    imWeight.Fill(0.0);
    im.FillRGB(ML_Black);
    fSqrt = double(BLUE_NSQRT)/sqrt(fSamplesPerPixel);
    if (im.m_nWidth > im.m_nHeight) {
        fSize = im.m_nWidth;
        vAdjust = Vector2(0.0, 0.5*(im.m_nWidth - im.m_nHeight));       // center rectangular window in unit square
    } else {
        fSize = im.m_nHeight;
        vAdjust = Vector2(0.5*(im.m_nHeight - im.m_nWidth), 0.0);
    }
    nSamples = 0;
    for (vTile.x = 0.0; vTile.x < im.m_nWidth; vTile.x += fSqrt) {
        for (vTile.y = 0.0; vTile.y < im.m_nHeight; vTile.y += fSqrt) {
            for (i  = 0; i < BLUE_NSQRT*BLUE_NSQRT; i++) {
                v = g_rgvBlueNoise[i]*fSqrt + vTile;
                g_uCamera = DiskSamples128[rgnQuasiIndex2[i & 127]].x;
                g_vCamera = DiskSamples128[rgnQuasiIndex2[i & 127]].y;
                g_uLight = DiskSamples128[rgnQuasiIndex3[i & 127]].x;
                g_vLight = DiskSamples128[rgnQuasiIndex3[i & 127]].y;
                g_uMicroFacet = DiskSamples128[rgnQuasiIndex5[i & 127]].x;
                g_vMicroFacet = DiskSamples128[rgnQuasiIndex5[i & 127]].y;
                rad = f((v + vAdjust) /fSize);
                im.SplatRGB(rad.sRGB(), v.x, v.y, &imWeight);
                nSamples++;
            }
        }
    }
    im.Normalize(&imWeight);
    return nSamples;
}

int
FastBlueSamplingTile::Initialize()
{
    return nSamples;
};

int
FastBlueSamplingTile::Render(Image &im, Radiance (*f)(Vector2 v), double fSamplesPerPixel)
{
    Image im2, im2Weight, imWeight;
    Vector2 vTile, vAdjust, v;
    Radiance rad;
    DisplayRGB rgb;
    double fSize;
    int i, nSamples, nW, nH;

    nW = int(im.m_nWidth * sqrt(fSamplesPerPixel) + 0.5);
    nH = int(im.m_nHeight * sqrt(fSamplesPerPixel) + 0.5);
    im2.NewImage(nW, nH, 3, WRAP_CLAMP, WRAP_CLAMP);
    im2.FillRGB(ML_Black);
    im2Weight.NewImage(nW, nH, 1, WRAP_CLAMP, WRAP_CLAMP);
    im2Weight.Fill(0.0);
    if (nW > nH) {
        fSize = nW;
        vAdjust = Vector2(0.0, 0.5*(nW - nH));       // center rectangular window in unit square
    } else {
        fSize = nH;
        vAdjust = Vector2(0.5*(nH - nW), 0.0);
    }
    nSamples = 0;
    for (vTile.x = 0.0; vTile.x < nW; vTile.x += BLUE_NSQRT) {
        for (vTile.y = 0.0; vTile.y < nH; vTile.y += BLUE_NSQRT) {
            for (i  = 0; i < BLUE_NSQRT*BLUE_NSQRT; i++) {
                v = g_rgvBlueNoise[i]*BLUE_NSQRT + vTile;
                /*
                g_uCamera = DiskSamples128[rgnQuasiIndex3[i & 127]].x;
                g_vCamera = DiskSamples128[rgnQuasiIndex3[i & 127]].y;
                g_uLight = DiskSamples128[rgnQuasiIndex2[i & 127]].x;
                g_vLight = DiskSamples128[rgnQuasiIndex2[i & 127]].y;
                g_uMicroFacet = DiskSamples128[rgnQuasiIndex5[i & 127]].x;
                g_vMicroFacet = DiskSamples128[rgnQuasiIndex5[i & 127]].y;
                */
                g_uCamera = 2.0*g_rgvBlueNoise[g_BlueNoiseQuasiIndex3[i % BLUE_N]].x - 1.0;
                g_vCamera = 2.0*g_rgvBlueNoise[g_BlueNoiseQuasiIndex3[i % BLUE_N]].y - 1.0;
                g_uLight = 2.0*g_rgvBlueNoise[g_BlueNoiseQuasiIndex2[i % BLUE_N]].x - 1.0;
                g_vLight = 2.0*g_rgvBlueNoise[g_BlueNoiseQuasiIndex2[i % BLUE_N]].y - 1.0;
                g_uMicroFacet = 2.0*g_rgvBlueNoise[g_BlueNoiseQuasiIndex2[i % BLUE_N]].x - 1.0;
                g_vMicroFacet = 2.0*g_rgvBlueNoise[g_BlueNoiseQuasiIndex2[i % BLUE_N]].y - 1.0;

                rad = f((v + vAdjust) /fSize);
                if (nVenera)
                    rgb = rad.vRGB();
                else
                    rgb = rad.sRGB();
                im2.FastSplatRGB(rgb, v.x, v.y, &im2Weight);
                nSamples++;
            }
        }
    }
    ImageResample(im, im.m_nWidth, im.m_nHeight, im2,  ML_KAISER4_HALFWIDTH, Kaiser4Filter);
    ImageResample(imWeight, im.m_nWidth, im.m_nHeight, im2Weight,  ML_KAISER4_HALFWIDTH, Kaiser4Filter);
    im.Normalize(&imWeight);
    // im2.Normalize(&im2Weight);
    // im2.WriteBMP("FastBlue.bmp");
    return nSamples;
}

void
FastBlueSamplingTile::QIndex1(int nPrime)
{
    int i;

    if (nPrime > 0) {
        BuildQuasiIndex(g_BlueNoiseQuasiIndex2, 2, BLUE_NSQRT*BLUE_NSQRT);
    } else {
        for (i = 0; i < BLUE_NSQRT*BLUE_NSQRT; i++)
            g_BlueNoiseQuasiIndex2[i] = i;
        if (nPrime < 0)
            RandomShuffle(g_BlueNoiseQuasiIndex2, BLUE_NSQRT*BLUE_NSQRT);
    }
}

int
AdaptiveBlueSamplingTile::Initialize()
{
    return nSamples;
};

#define THRESHOLD1  0.30
#define THRESHOLD2  0.30
#define HYSTERESIS  2
#define NSQRT1      2
#define HSQRT2      3

int
AdaptiveBlueSamplingTile::Render(Image &im, Radiance (*f)(Vector2 v), double fSamplesPerPixel)
{
    Image im2, im2Weight, imWeight, imDetect;
    Vector2 vTile, vAdjust, v;
    Radiance rad;
    DisplayRGB rgb;
    double fSize, f00, f01, f10, f11, g1, g2, g, sum, maxg;
    int i, j, k, is, js, nSamples, nW, nH, nSubSample;

    nSubSample = 1;
    fSamplesPerPixel /= double(nSubSample);
    nW = int(im.m_nWidth * sqrt(fSamplesPerPixel) + 0.5);
    nH = int(im.m_nHeight * sqrt(fSamplesPerPixel) + 0.5);
    im2.NewImage(nW, nH, 3, WRAP_CLAMP, WRAP_CLAMP);
    im2.FillRGB(ML_Black);
    im2Weight.NewImage(nW, nH, 1, WRAP_CLAMP, WRAP_CLAMP);
    im2Weight.Fill(0.0);
    imDetect.NewImage(nW, nH, 1, WRAP_MIRROR, WRAP_MIRROR);
    if (nW > nH) {
        fSize = nW;
        vAdjust = Vector2(0.0, 0.5*(nW - nH));       // center rectangular window in unit square
    } else {
        fSize = nH;
        vAdjust = Vector2(0.5*(nH - nW), 0.0);
    }
    //
    //  Base sampling rate
    //
    nSamples = 0;
    for (vTile.x = 0.0; vTile.x < nW; vTile.x += BLUE_NSQRT) {
        for (vTile.y = 0.0; vTile.y < nH; vTile.y += BLUE_NSQRT) {
            for (k  = 0; k < BLUE_NSQRT*BLUE_NSQRT; k += nSubSample) {
                v = g_rgvBlueNoise[k]*BLUE_NSQRT + vTile;
                //
                //  Auxiliary sampling dimensions
                //
                g_uCamera = 2.0*g_rgvBlueNoise[g_BlueNoiseQuasiIndex3[k % BLUE_N]].x - 1.0;
                g_vCamera = 2.0*g_rgvBlueNoise[g_BlueNoiseQuasiIndex3[k % BLUE_N]].y - 1.0;
                g_uLight = 2.0*g_rgvBlueNoise[g_BlueNoiseQuasiIndex2[k % BLUE_N]].x - 1.0;
                g_vLight = 2.0*g_rgvBlueNoise[g_BlueNoiseQuasiIndex2[k % BLUE_N]].y - 1.0;
                g_uMicroFacet = 2.0*g_rgvBlueNoise[g_BlueNoiseQuasiIndex2[k % BLUE_N]].x - 1.0;
                g_vMicroFacet = 2.0*g_rgvBlueNoise[g_BlueNoiseQuasiIndex2[k % BLUE_N]].y - 1.0;
                //
                //  Sampling the image function
                //
                rad = f((v + vAdjust) /fSize);
                if (nVenera)
                    rgb = rad.vRGB();
                else
                    rgb = rad.sRGB();
                im2.FastSplatRGB(rgb, v.x, v.y, &im2Weight);
                nSamples++;
            }
        }
    }
    //
    //  normalized luminance, and then Robert's gradient magnitude to measure high frequency features
    //
    for (j = 0; j < nH; j++) {
        for (i = 0; i < nW; i++) {
            f00 = im2.GetRGB(i+0, j+0).Luminance();
            imDetect.Set(f00, i, j);
        }
    }
    imDetect.Normalize(&im2Weight);
    imDetect.WriteBMP("Luminance.bmp");
    for (j = 0; j < nH; j++) {
        for (i = 0; i < nW; i++) {
            f00 = imDetect.Get(i+0, j+0);
            f01 = imDetect.Get(i+0, j+1);
            f10 = imDetect.Get(i+1, j+0);
            f11 = imDetect.Get(i+1, j+1);
            g1 = f00 - f11;
            g2 = f01 - f10;
            sum = f00+f01+f10+f11;
            if (sum > 0.0)
                g = sqrt(g1*g1 + g2*g2)/sum * 1.0;
            else
                g = 0.0;
            imDetect.Set(g, i, j);
        }
    }
    imDetect.WriteBMP("Roberts.bmp", 1.0);
    //
    //  Analyse and supersample regions with gradients higher than the thresholds
    //
    for (j = 0; j < nH; j += BLUE_NSQRT/NSQRT1) {
        for (i = 0; i < nW; i += BLUE_NSQRT/NSQRT1) {
            maxg = 0;
            for (is = i - HYSTERESIS; is < i + BLUE_NSQRT/NSQRT1 + HYSTERESIS; is++) {
                for (js = j - HYSTERESIS; js < j + BLUE_NSQRT/NSQRT1 + HYSTERESIS; js++) {
                    g = imDetect.Get(is, js);
                    if (g > maxg)
                        maxg = g;
                }
            }
            if (maxg > THRESHOLD1) {
                vTile = Vector2(double(i), double(j));
                for (k  = 0; k < BLUE_NSQRT*BLUE_NSQRT; k++) {
                    g_uCamera = 2.0*g_rgvBlueNoise[g_BlueNoiseQuasiIndex3[k % BLUE_N]].x - 1.0;
                    g_vCamera = 2.0*g_rgvBlueNoise[g_BlueNoiseQuasiIndex3[k % BLUE_N]].y - 1.0;
                    g_uLight = 2.0*g_rgvBlueNoise[g_BlueNoiseQuasiIndex2[k % BLUE_N]].x - 1.0;
                    g_vLight = 2.0*g_rgvBlueNoise[g_BlueNoiseQuasiIndex2[k % BLUE_N]].y - 1.0;
                    g_uMicroFacet = 2.0*g_rgvBlueNoise[g_BlueNoiseQuasiIndex2[k % BLUE_N]].x - 1.0;
                    g_vMicroFacet = 2.0*g_rgvBlueNoise[g_BlueNoiseQuasiIndex2[k % BLUE_N]].y - 1.0;
                    v = g_rgvBlueNoise[k]*BLUE_NSQRT/NSQRT1 + vTile;
                    rad = f((v + vAdjust) /fSize);
                    if (nVenera)
                        rgb = rad.vRGB();
                    else
                        rgb = rad.sRGB();
                    //  rgb.red = 1.0;
                    im2.FastSplatRGB(rgb, v.x, v.y, &im2Weight);
                }
            }
        }
    }
    //
    //  Resample the image to output size
    //
    if (fSamplesPerPixel != 1.0) {
        ImageResample(im, im.m_nWidth, im.m_nHeight, im2,  ML_KAISER4_HALFWIDTH, Kaiser4Filter);
        ImageResample(imWeight, im.m_nWidth, im.m_nHeight, im2Weight,  ML_KAISER4_HALFWIDTH, Kaiser4Filter);
        im.Normalize(&imWeight);
    } else {
        delete im.m_pfImage;
        im.m_pfImage = im2.m_pfImage;
        im2.m_pfImage = 0;
        im.Normalize(&im2Weight);
    }
    return nSamples;
}

void
AdaptiveBlueSamplingTile::QIndex1(int nPrime)
{
    int i;

    if (nPrime > 0) {
        BuildQuasiIndex(g_BlueNoiseQuasiIndex2, 2, BLUE_N);
    } else {
        for (i = 0; i < BLUE_N; i++)
            g_BlueNoiseQuasiIndex2[i] = i;
        if (nPrime < 0)
            RandomShuffle(g_BlueNoiseQuasiIndex2, BLUE_N);
    }
}