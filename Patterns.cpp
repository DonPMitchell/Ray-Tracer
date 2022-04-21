#include "stdafx.h"
#include "Patterns.h"
#include "Montecarlo.h"
#include "Image.h"
#include "Utilities.h"
#include "MLsorting.h"
#pragma intrinsic(exp, sin, cos, sqrt, fabs, sqrt)
#define D_PI            3.14159265358979323846264338327950288419716939937510
#define D_PHI           ((sqrt(5.0) + 1.0)/2.0)

//
//  Data structures for locality and spectral optimization
//
struct Sample {
    Vector2  v;
    Sample   *psNext;       // grid-cell contents list
    short    j, k;
};

static Sample   *s_rgpsGrid[32][32], *s_psSamples = 0;
static int      s_nSqrt;


//
//  Uniform grid of samples
//
int
UniformPattern(Vector2 rgv[], int nSqrt)
{
    int k, n;
    double x, y, fScale;
    unsigned i, j;

    n = nSqrt*nSqrt;
    fScale = 1.0/double(nSqrt);
    for (k = 0; k < n; k++) {
        i = k % nSqrt;
        j = k / nSqrt;
        x = fScale*(double(i) + 0.5);
        y = fScale*(double(j) + 0.5);
        *rgv++ = Vector2(x, y);
    }
    return n;
}
//
//  Stratified sampling, one sample per grid area
//
int
JitterPattern(Vector2 rgv[], int nSqrt, double fJitter)
{
    int k, n;
    double x, y, fScale;
    unsigned i, j;

    n = nSqrt*nSqrt;
    fScale = 1.0/double(nSqrt);
    for (k = 0; k < n; k++) {
        i = k % nSqrt;
        j = k / nSqrt;
        x = fScale*(double(i) + fJitter*RandomDouble());
        y = fScale*(double(j) + fJitter*RandomDouble());
        *rgv++ = Vector2(x, y);
    }
    return n;
}

static void
ML_DeInterLeave(const unsigned rgnLocationCode[2], unsigned &ix, unsigned &iy)
{
    unsigned n;

    n = rgnLocationCode[0];
    n = ((n & 0x44444444) >> 1) | ((n & 0x22222222) << 1) | (n & 0x99999999);
    n = ((n & 0x30303030) >> 2) | ((n & 0x0C0C0C0C) << 2) | (n & 0xC3C3C3C3);
    n = ((n & 0x0F000F00) >> 4) | ((n & 0x00F000F0) << 4) | (n & 0xF00FF00F);
    n = ((n & 0x00FF0000) >> 8) | ((n & 0x0000FF00) << 8) | (n & 0xFF0000FF);
    ix = n & 0xFFFF;
    iy = n >> 16;
    n = rgnLocationCode[1];
    n = ((n & 0x44444444) >> 1) | ((n & 0x22222222) << 1) | (n & 0x99999999);
    n = ((n & 0x30303030) >> 2) | ((n & 0x0C0C0C0C) << 2) | (n & 0xC3C3C3C3);
    n = ((n & 0x0F000F00) >> 4) | ((n & 0x00F000F0) << 4) | (n & 0xF00FF00F);
    n = ((n & 0x00FF0000) >> 8) | ((n & 0x0000FF00) << 8) | (n & 0xFF0000FF);
    ix |= n << 16;
    iy |= n & 0xFFFF0000;
}

static void
ML_HilbertDecode(const unsigned rgnLocationCode[2], unsigned &ix, unsigned &iy)     // only works for size == 2**N
{
    unsigned nOdd, nEven, nV0, nV1, nTemp0, nTemp1;
    int k;

    ML_DeInterLeave(rgnLocationCode, nEven, nOdd);
    nV0 = nV1 = 0;
    nTemp1 = ~(nEven | nOdd);
    nTemp0 = ~(nEven ^ nOdd);
    for (k = 1; k < 32; k++) {
        nV1 = (nV1 ^ nTemp1) >> 1;
        nV0 = (nV0 ^ nTemp0) >> 1;
    }
    iy = (nV0 & ~nEven) ^ nV1 ^ nOdd;
    ix = (nV0 | nEven) ^ nV1 ^ nOdd;
}

int
JitterHilbert(Vector2 rgv[], int nSqrt, double fJitter)
{
    int k, n;
    double fScale;

    n = nSqrt*nSqrt;
    fScale = 1.0/double(nSqrt);
    GeneralHilbert(rgv, nSqrt, nSqrt);
    for (k = 0; k < n; k++) {
        rgv[k].x += fScale*fJitter*RandomDouble();
        rgv[k].y += fScale*fJitter*RandomDouble();
    }
    return n;
}
//
//  Generalized hilbert curve for arbitrary size
//
static Vector2 *s_pvHilbert;

static inline double
sign(double x)
{
    if (x < 0.0)
        return -1.0;
    else if (x > 0.0)
        return +1.0;
    else
        return 0.0;
}

static inline Vector2
floor(const Vector2 &v)
{
    return Vector2(floor(v.x + 0.1), floor(v.y + 0.1));
}

static void
DoHilbert(Vector2 v, Vector2 a, Vector2 b)
{
    int w, h, w2, h2;
    Vector2 da, db, a2, b2;
    int i;

    w = abs(int(a.x + a.y));
    h = abs(int(b.x + b.y));
    da = Vector2(sign(a.x), sign(a.y));
    db = Vector2(sign(b.x), sign(b.y));
    //
    //  fill simple columns and rows
    //
    if (h == 1) {
        for (i = 0; i < w; i++) {               
            *s_pvHilbert++ = v;
            v = v + da;
        }
        return;
    }
    if (w == 1) {
        for (i = 0; i < h; i++) {
            *s_pvHilbert++ = v;
            v = v + db;
        }
        return;
    }
    a2 = floor(a * 0.5);
    b2 = floor(b * 0.5);
    w2 = abs(int(a2.x + a2.y));
    h2 = abs(int(b2.x + b2.y));
    if (2*w > 3*h) {
        if ((w2 % 2) && (w > 2))    // prefer even steps
            a2 = a2 + da;
        DoHilbert(v, a2, b);
        DoHilbert(v + a2, a - a2, b);
    } else {
        if ((h2 % 2) && (h > 2))
            b2 = b2 + db;
        DoHilbert(v, b2, a2);
        DoHilbert(v+b2, a, b - b2);
        DoHilbert(v + (a-da) + (b2-db), -b2, -(a-a2));
    }
}

void
GeneralHilbert(Vector2 rgv[], int nW, int nH)       // rgv[] must have nW*nH elements
{
    int i;

    s_pvHilbert = rgv;
    if (nW >= nH)
        DoHilbert(Vector2(0.0, 0.0), Vector2(nW, 0.0), Vector2(0.0, nH));
    else
        DoHilbert(Vector2(0.0, 0.0), Vector2(0.0, nH), Vector2(nW, 0.0));
    for (i = 0; i < nW*nH; i++) {
        rgv[i].x /= nW;
        rgv[i].y /= nH;
    }
}

static void
DrawOnTorus(Image &im, double x1, double y1, double x2, double y2, int nW, int nH, int nB = 0)
{
    if (x2 - x1 >  nW/2) x2 -= nW;
    if (x2 - x1 < -nW/2) x2 += nW;
    if (y2 - y1 >  nH/2) y2 -= nH;
    if (y2 - y1 < -nH/2) y2 += nH;
    im.DrawLineRGB(ML_White, x1+nB, y1+nB, x2+nB, y2+nB, 1.5);
}

#define NGH_WIDE 40
#define NGH_HIGH 40
#define NGH     (NGH_WIDE*NGH_HIGH)

void
TestGeneralHilbert()
{
    Image im;
    static Vector2 rgv[NGH];
    double x1, x2, y1, y2;
    int i, nW, nH, nB;

    nW = 1280;
    nH = (nW*NGH_HIGH) / NGH_WIDE;
    nB = nW/NGH_WIDE;
    im.NewImage(nW, nH, 3);
    im.FillRGB(ML_Black);
    //GeneralHilbert(rgv, NGH_WIDE, NGH_HIGH);
    BlueJitterPattern(rgv, 40);
    for (i = 0; i < NGH; i++) {
        //  rgv[i].Print(); printf("\n");
        if (i == NGH-1) 
            break;
        
        x1 = rgv[i].x;  // / NGH_WIDE;
        y1 = rgv[i].y;  // / NGH_HIGH;
        x2 = rgv[i+1].x; // / NGH_WIDE;
        y2 = rgv[i+1].y; // / NGH_HIGH;
        
        DrawOnTorus(im, x1*double(nW-nB), y1*double(nH-nB), x2*double(nW-nB), y2*double(nH-nB), nW, nH, nB);
    }
    im.WriteBMP("BlueJitterHilbert.bmp");
}

//
//  Uniform random points on unit square.  Only used in some testing.
//
int
RandomPattern(Vector2 rgv[], int nSqrt)
{
    int k, n;
    double x, y;

    n = nSqrt*nSqrt;
    for (k = 0; k < n; k++) {
        x = RandomDouble();
        y = RandomDouble();
        *rgv++ = Vector2(x, y);
    }
    return n;
}

//
//  Find samples (on the torus) within fRadius of a point
//REVIEW: rgn not used?
//
static int
Neighborhood(Vector2 rgv[], int nSqrt, const Vector2 &vCenter, double fRadius)
{
    int jLo, jHi, kLo, kHi, j, k, dj, dk, i;
    Vector2 dv;
    Sample *ps;

    jLo = int(floor(32.0*(vCenter.x - fRadius)));
    jHi = int(floor(32.0*(vCenter.x + fRadius)));
    kLo = int(floor(32.0*(vCenter.y - fRadius)));
    kHi = int(floor(32.0*(vCenter.y + fRadius)));
    i = 0;
    for (j = jLo; j <= jHi; j++) {
        if (j < 0) {
            dj = 32;
            dv.x = -1.0;
        } else if (j >= 32) {
            dj = -32;
            dv.x = 1.0;
        } else {
            dj = 0;
            dv.x = 0.0;
        }
        for (k = kLo; k <= kHi; k++) {
            if (k < 0) {
                dk = 32;
                dv.y = -1.0;
            } else if (k >= 32) {
                dk = -32;
                dv.y = 1.0;
            } else {
                dk = 0;
                dv.y = 0.0;
            }
            for (ps = s_rgpsGrid[j+dj][k+dk]; ps; ps = ps->psNext) {
                if (Norm(ps->v + dv - vCenter) < fRadius) {
                    rgv[i++] = ps->v + dv;
                }
            }
        }
    }
    return i;
}
//
//  Change the coordinates of the i-th sample
//
static void
Move(int i, const Vector2 &vNew)
{
    int j, k;
    Sample *psBin, *ps;
    double x, y;

    ps = s_psSamples + i;
    x = vNew.x;
    y = vNew.y;
    while (x < 0.0)
        x += 1.0;
    while (x >= 1.0)
        x -= 1.0;
    while (y < 0.0)
        y += 1.0;
    while (y >= 1.0)
        y -= 1.0;
    j = int(floor(32.0*x));
    k = int(floor(32.0*y));
    if (j != ps->j || k != ps->k) {
        //
        //  Delete sample from old grid's bin and add it to the new grid.
        //
        psBin = s_rgpsGrid[ps->j][ps->k];
        if (psBin == ps) {
            s_rgpsGrid[ps->j][ps->k] = psBin->psNext;
        } else {
            for (; psBin; psBin = psBin->psNext)
                if (psBin->psNext == ps) {
                    psBin->psNext = ps->psNext;
                    break;
                }
        }
        ps->j = j;
        ps->k = k;
        ps->psNext = s_rgpsGrid[j][k];
        s_rgpsGrid[j][k] = ps;
    }
    ps->v = Vector2(x, y);
}
//
//  Load/unload a sample pattern in the search-grid structure
//
static int
LoadPattern(Vector2 rgv[], int nSqrt)
{
    int i, j, k, n;

    n = nSqrt*nSqrt;
    s_nSqrt = nSqrt;
    delete [] s_psSamples;
    s_psSamples = new Sample[n];
    if (s_psSamples == 0)
        return 0;
    for (i = 0; i < 32; i++)
        for (j = 0; j < 32; j++)
            s_rgpsGrid[i][j] = 0;
    for (i = 0; i < n; i++) {
        s_psSamples[i].v = rgv[i];
        j = int(floor(32.0*s_psSamples[i].v.x));
        k = int(floor(32.0*s_psSamples[i].v.y));
        s_psSamples[i].j = j;
        s_psSamples[i].k = k;
        s_psSamples[i].psNext = s_rgpsGrid[j][k];
        s_rgpsGrid[j][k] = s_psSamples + i;
    }
    return 1;
}

static int
UnloadPattern(Vector2 rgv[], int nSqrt)
{
    int i, n;

    if (s_psSamples == 0)
        return 0;
    n = nSqrt*nSqrt;
    for (i = 0; i < n; i++)
        rgv[i] = s_psSamples[i].v;
    delete [] s_psSamples;
    s_psSamples = 0;
    return n;
}

//
//  Generate a poisson-disk sampling pattern
//
int
PoissonDiskPattern(Vector2 rgv[], int nSqrt)
{
    int nIteration, i, j, nNHood, nSamples, nResult;
    int *pnPermute;
    Vector2 vClosest, vNew, vOld;
    Vector2 *pvNHood;
    double fDist, fClosest, fRadius, fJitter, fAveClosest, fN;

    nResult = 0;
    pnPermute = 0;
    pvNHood = 0;
    nSamples = nSqrt*nSqrt;
    pnPermute = new int[nSamples];
    pvNHood = new Vector2[nSamples];
    if (pnPermute == 0 || pvNHood == 0)
        goto Fail;
    for (i = 0; i < nSamples; i++)
        pnPermute[i] = i;
    fRadius = 2.0f/float(nSqrt);
    fJitter = 0.02f;
    JitterPattern(rgv, nSqrt);
    LoadPattern(rgv, nSqrt);
    for (nIteration = 0; nIteration < 130; nIteration++) {
        //
        //  Keep permuting the order that points are visited to avoid bias
        //
        RandomShuffle(pnPermute, nSamples);
        if ((nIteration & 31) == 31)
            fJitter = 0.1f*fJitter;
        fAveClosest = fN = 0.0;
        for (i = 0; i < nSamples; i++) {
            vOld = s_psSamples[(pnPermute[i])].v;
            nNHood = Neighborhood(pvNHood, nSqrt, vOld, fRadius);
            if (nNHood < 2)
                continue;
            fClosest = 10.0;
            for (j = 0; j < nNHood; j++) {
                if (pvNHood[j] == vOld)
                    continue;
                fDist = Norm(pvNHood[j] - vOld);
                if (fDist < fClosest) {
                    fClosest = fDist;
                    vClosest = pvNHood[j];
                }
            }
            fAveClosest += fClosest;
            fN += 1.0;
            vNew = vOld - Normalize(vClosest - vOld)*RandomDouble()*fJitter;
            Move(pnPermute[i], vNew);
        }
    }
    //printf("Average Closest Distance: %g\n", sqrt(SP_NTILE.0)*fAveClosest/fN);
    nResult = nSamples;
    UnloadPattern(rgv, nSqrt);
Fail:
    delete [] pnPermute;
    delete [] pvNHood;
    return nResult;
}
//
//  Generate a more optimal blue-noise pattern
//
static double
Sinc(double x)
{
    x *= D_PI;
    if (x < 0.01 && x > -0.01) {
        x = x*x;
        return 1.0 + x*(-1.0/6.0 + x*(1.0/120.0 + x*(-1.0/5040.0 + x/362880.0)));
    } else
        return sin(x)/x;
}

double
Jinc(double r)
{
    r *= D_PI;
    if (r < 0.01 && r > -0.01) {
        r = r*r;
        return 0.5 + r*(-1.0/16.0 + r*1.0/384.0);
    } else
        return _j1(r)/r;
}

static double
Blueness(Vector2 &vCenter, Vector2 *pvNHood)
{
    int i, n;
    double fRadius, fBlueness;

    fRadius = 1.0/float(s_nSqrt);
    n = Neighborhood(pvNHood, s_nSqrt, vCenter, 4.0*fRadius);
    if (n < 2)
        return 0.0;
    fBlueness = 0.0;
    for (i = 0; i < n; i++) {
        if (pvNHood[i] == vCenter)
            continue;
        //fBlueness += 1.0 - Jinc((5.52007/D_PI) * Norm(pvNHood[i] - vCenter)/fRadius);
        fBlueness += 1.0 - Sinc(sqrt(2.0)*(pvNHood[i].x - vCenter.x)/fRadius)                       // 1.4
                         * Sinc(sqrt(2.0)*(pvNHood[i].y - vCenter.y)/fRadius);
    }
    return fBlueness/double(n - 1);
}

int
BlueNoisePattern(Vector2 rgv[], int nSqrt)
{
    int nIteration, i, nResult, nSamples, nActual, *pnPermute;
    Vector2 vClosest, vNew, vOld, *pvNHood;
    double fRadius, fJitter;
    double fOldBlue, fNewBlue;

    nResult = 0;
    pnPermute = 0;
    pvNHood = 0;
    nSamples = nSqrt*nSqrt;
    pnPermute = new int[nSamples];
    pvNHood = new Vector2[nSamples];
    if (pnPermute == 0 || pvNHood == 0)
        goto Fail;
    for (i = 0; i < nSamples; i++)
        pnPermute[i] = i;
    fRadius = 2.0f/float(nSqrt);
    fJitter = 0.01f;
    if ((nActual = PoissonDiskPattern(rgv, nSqrt)) != nSamples)
        goto Fail;
    //
    //  Perturb poisson-disk pattern into higher-quality blue noise
    //
    LoadPattern(rgv, nSqrt);
    fJitter = 0.01 * fRadius;
    for (nIteration = 0; nIteration < 200; nIteration++) {
        RandomShuffle(pnPermute, nSamples);
        if (nIteration == 100)
            fJitter *= 0.1;
        for (i = 0; i < nSamples; i++) {
            vOld = s_psSamples[pnPermute[i]].v;
            fOldBlue = Blueness(vOld, pvNHood);
            vNew = vOld + SquareToUnitDisk(Vector2(RandomDouble(), RandomDouble()))*fJitter;
            Move(pnPermute[i], vNew);
            fNewBlue = Blueness(vNew, pvNHood);
            if (fOldBlue >= fNewBlue)
                Move(pnPermute[i], vOld);
            //else
            //    printf("Blueness %f to %f\n", fOldBlue, fNewBlue);
        }
    }
    UnloadPattern(rgv, nSqrt);
    nResult = nSamples;
    goto NoFail;
Fail:
    printf("failed N samples %d vs %d\n", nActual, nSamples);
NoFail:
    delete [] pvNHood;
    delete [] pnPermute;
    return nResult;
}

int
BlueJitterPattern(Vector2 rgv[], int nSqrt)
{
    int nIteration, i, nResult, nSamples, nActual, *pnPermute, iox, ioy, inx, iny;
    Vector2 vClosest, vNew, vOld, *pvNHood;
    double fRadius, fJitter;
    double fOldBlue, fNewBlue;
    extern void PatternTestImage(char *szFileName, Vector2 rgv[], int nSamples);

    nResult = 0;
    pnPermute = 0;
    pvNHood = 0;
    nSamples = nSqrt*nSqrt;
    pnPermute = new int[nSamples];
    pvNHood = new Vector2[nSamples];
    if (pnPermute == 0 || pvNHood == 0)
        goto Fail;
    for (i = 0; i < nSamples; i++)
        pnPermute[i] = i;
    fRadius = 1.0/float(nSqrt);
    fJitter = 0.01f;
    if ((nActual = JitterHilbert(rgv, nSqrt)) != nSamples)
        goto Fail;
    PatternTestImage("JitterBefore.bmp", rgv, nSamples);
    //
    //  Perturb poisson-disk pattern into higher-quality blue noise
    //
    LoadPattern(rgv, nSqrt);
    fJitter = 0.01 * fRadius;
    for (nIteration = 1; nIteration < 2000; nIteration++) {
        RandomShuffle(pnPermute, nSamples);
        if (nIteration % 100 == 0)
            fJitter *= sqrt(0.5);
        for (i = 0; i < nSamples; i++) {
            vOld = s_psSamples[pnPermute[i]].v;
            iox = int(floor(vOld.x * nSqrt));
            ioy = int(floor(vOld.y * nSqrt));
            fOldBlue = Blueness(vOld, pvNHood);
            vNew = vOld + SquareToUnitDisk(Vector2(RandomDouble(), RandomDouble()))*fJitter;
            inx = int(floor(vNew.x * nSqrt));
            iny = int(floor(vNew.y * nSqrt));
            if (inx < iox) vNew.x += fRadius;       // keep point within its cell
            if (inx > iox) vNew.x -= fRadius;
            if (iny < ioy) vNew.y += fRadius;
            if (iny > ioy) vNew.y -= fRadius;
            Move(pnPermute[i], vNew);
            fNewBlue = Blueness(vNew, pvNHood);
            if (fOldBlue >= fNewBlue || RandomDouble() < 0.01)
                Move(pnPermute[i], vOld);
            //else
            //    printf("Blueness %f to %f\n", fOldBlue, fNewBlue);
        }
    }
    UnloadPattern(rgv, nSqrt);
    PatternTestImage("JitterAfter.bmp", rgv, nSamples);
    nResult = nSamples;
    goto NoFail;
Fail:
    printf("failed N samples %d vs %d\n", nActual, nSamples);
NoFail:
    delete [] pvNHood;
    delete [] pnPermute;
    return nResult;
}

//
//  Reorder the 2D samples to maximize locality.  Uses simulated annealing.
//
static int s_nWrapLocality = 1;

static void
Locality(double *pfAve, double *pfMax, Vector2 rgv[], int nSamples)
{
    double dx, dy, di, fMax, fSum, fDist;
    int i, j;

    fMax = 0.0;
    fSum = 0.0;
    //
    //  2D distance is measured on the torus, but 1D index distance is
    //  measured on the segment (not the circle).
    //
    for (i = 0; i < nSamples; i++) {
        for (j = 0; j < i; j++) {
            dx = fabs(rgv[i].x - rgv[j].x);
            dy = fabs(rgv[i].y - rgv[j].y);
            if (s_nWrapLocality) {
                if (dx > 0.5)
                    dx = 1.0 - dx;
                if (dy > 0.5)
                    dy = 1.0 - dy;
            }
            di = fabs(double(i - j));
            fDist = (dx*dx + dy*dy)/di;
            if (fDist > fMax)
                fMax = fDist;
            fSum += fDist;
        }
    }
    *pfAve = fSum;
    *pfMax = nSamples*fMax;
}

static void
PartialLocality(double *pfAve, double *pfMax, Vector2 rgv[], int nSamples, int i)
{
    double fSum, dx, dy, di, fMax, fDist;
    int j;

    fSum = 0.0;
    fMax = 0.0;
    for (j = 0; j < nSamples; j++) {
        if (j == i)
            continue;
        dx = fabs(rgv[i].x - rgv[j].x);
        dy = fabs(rgv[i].y - rgv[j].y);
        if (s_nWrapLocality) {
            if (dx > 0.5)
                dx = 1.0 - dx;
            if (dy > 0.5)
                dy = 1.0 - dy;
        }
        di = fabs(double(i - j));
        fDist = (dx*dx + dy*dy)/di;
        if (fDist > fMax)
            fMax = fDist;
        fSum += fDist;
    }
    *pfAve = -fSum;
    *pfMax = -fMax;
}

static void
DeltaEnergy(double *pfDeltaAve, double *pfDeltaMax, Vector2 rgv[], int nSamples, int i1, int i2)
{
    double dE1, dE2, dM1, dM2, fAve1, fMax1, fAve2, fMax2;
    Vector2 vTmp;

    PartialLocality(&fAve1, &fMax1, rgv, nSamples, i1);
    PartialLocality(&fAve2, &fMax2, rgv, nSamples, i2);
    dE1 = fAve1 + fAve2;
    dM1 = fMax1;
    if (fMax2 < dM1)
        dM1 = fMax2;
    vTmp = rgv[i1];
    rgv[i1] = rgv[i2];
    rgv[i2] = vTmp;
    PartialLocality(&fAve1, &fMax1, rgv, nSamples, i1);
    PartialLocality(&fAve2, &fMax2, rgv, nSamples, i2);
    dE2 = fAve1 + fAve2;
    dM2 = fMax1;
    if (fMax2 < dM2)
        dM2 = fMax2;
    vTmp = rgv[i1];
    rgv[i1] = rgv[i2];
    rgv[i2] = vTmp;
    *pfDeltaAve = dE1 - dE2;
    *pfDeltaMax = dM1 - dM2;
}

static double
Distance(const Vector2 &v1, const Vector2 &v2)
{
    double dx, dy;

    dx = v1.x - v2.x;
    dy = v1.y - v2.y;
    if (dx > 0.5)
        dx = 1.0 - dx;
    if (dy > 0.5)
        dy = 1.0 - dy;
    return sqrt(dx*dx + dy*dy);
}

static int inline
Intersect(const Vector2 &vStart1, const Vector2 &vEnd1,
          const Vector2 &vStart2, const Vector2 &vEnd2)
{
    int n1, n2;

    n1 = CCW(vStart1, vEnd1, vStart2) * CCW(vStart1, vEnd1, vEnd2);
    n2 = CCW(vStart2, vEnd2, vStart1) * CCW(vStart2, vEnd2, vEnd1);

    return (n1 < 0 && n2 <= 0) || (n1 <= 0 && n2 < 0);
}

static void
Swap_2_Opt(Vector2 rgv[], int i, int k, int n)
{
    int m;
    Vector2 vTemp;

    m = (i + k + 1)/2;
    for (; i < m; i++, --k) {
        vTemp = rgv[i];
        rgv[i] = rgv[k];
        rgv[k] = vTemp;
    }
}

int
OptimizeLocality(Vector2 rgv[], int nSamples, int nWrap)
{
    double fTemperature, fInitialTemp, dEnergy, dEnergyAve, dEnergyMax, fAve, fMax;
    int i, j, i1, i2, nCooling, nTurns, maxTurns, nSucc, maxSucc;
    int *pnPermute1, *pnPermute2, nResult;
    Vector2 vS1, vE1, vS2, vE2, vVert, vTmp;

    s_nWrapLocality = nWrap;
    nResult = 0;
    pnPermute1 = pnPermute2 = 0;
    pnPermute1 = new int[nSamples];
    pnPermute2 = new int[nSamples];
    if (pnPermute1 == 0 || pnPermute2 == 0)
        goto Fail;
    for (i = 0; i < nSamples; i++)
        pnPermute1[i] = pnPermute2[i] = i;
    RandomShuffle(rgv, nSamples);
    maxTurns = 100*nSamples;
    maxSucc = 10*nSamples;
    Locality(&fAve, &fMax, rgv, nSamples);
    printf("starting Locality = %f max, %f ave\n", fMax, fAve);
    //
    //  Find a good starting temperature for average locality
    //
    fInitialTemp = 0.0;
    for (nTurns = 0; nTurns < 100; nTurns++) {
        i = RandomLessThanN(nSamples);
        do {
            j = RandomLessThanN(nSamples);
        } while (j == i);
        DeltaEnergy(&dEnergyAve, &dEnergyMax, rgv, nSamples, i, j);
        dEnergy = fabs(dEnergyAve);
        if (dEnergy > fInitialTemp)
            fInitialTemp = dEnergy;
    }
    fTemperature = 5.0*fInitialTemp;
    printf(" initial temp %f\n", fInitialTemp);
    //
    //  Anneal average locality at first, max locality won't converge starting
    //  from random order.
    //
    for (nCooling = 0; fTemperature > 1.0e-6*fInitialTemp; nCooling++) {
        nSucc = 0;
        RandomShuffle(pnPermute1, nSamples);
        RandomShuffle(pnPermute2, nSamples);
        for (nTurns = 0; nTurns < nSamples; nTurns++) {
            i = pnPermute1[nTurns];
            j = pnPermute2[nTurns];
            while (i == j)
                j = RandomLessThanN(nSamples);
            DeltaEnergy(&dEnergyAve, &dEnergyMax, rgv, nSamples, i, j);
            if (dEnergyAve < 0.0 || RandomDouble() < exp(-dEnergyAve/fTemperature)) {
                vTmp = rgv[i];
                rgv[i] = rgv[j];
                rgv[j] = vTmp;
                nSucc++;
            }
            if (nSucc > maxSucc)
                break;
        }
        fTemperature *= 0.999;   // Lower temperature 1 percent
    }
    printf("%d cooling turns, ", nCooling);
    printf(" %f final temp\n", fTemperature);
    Locality(&fAve, &fMax, rgv, nSamples);
    printf("after annealing Locality = %f max, %f ave\n", fMax, fAve);
    //
    //  Switch to annealing maximum locality
    //
    fInitialTemp = 0.0;
    for (nTurns = 0; nTurns < 100; nTurns++) {
        i = RandomLessThanN(nSamples);
        do {
            j = RandomLessThanN(nSamples);
        } while (j == i);
        DeltaEnergy(&dEnergyAve, &dEnergyMax, rgv, nSamples, i, j);
        dEnergy = fabs(dEnergyMax);
        if (dEnergy > fInitialTemp)
            fInitialTemp = dEnergy;
    }
    fTemperature = 0.0001*fInitialTemp;
    printf(" %f initial temp\n", fInitialTemp);
    for (nCooling = 0; fTemperature > 1.0e-8*fInitialTemp; nCooling++) {
        nSucc = 0;
        RandomShuffle(pnPermute1, nSamples);
        RandomShuffle(pnPermute2, nSamples);
        for (nTurns = 0; nTurns < nSamples; nTurns++) {
            i = pnPermute1[nTurns];
            j = pnPermute2[nTurns];
            while (i == j)
                j = RandomLessThanN(nSamples);
            DeltaEnergy(&dEnergyAve, &dEnergyMax, rgv, nSamples, i, j);
            if (dEnergyMax < 0.0 || RandomDouble() < exp(-dEnergyMax/fTemperature)) {
                vTmp = rgv[i];
                rgv[i] = rgv[j];
                rgv[j] = vTmp;
                nSucc++;
                //ML_Printf("%4d ", nCooling);
                //ML_PrintFloat(dEnergyMax, " dMax ");
                //ML_PrintFloat(dEnergyAve, " dAve\n");
            }
            if (nSucc > maxSucc)
                break;
        }
        fTemperature *= 0.995;   // Lower temperature 5 percent
    }
    printf("%d cooling turns, ", nCooling);
    printf(" %f final temp\n", fTemperature);
    Locality(&fAve, &fMax, rgv, nSamples);
    printf("after annealing Locality = ");
    printf("%f max %f ave\n", fMax, fAve);

    //
    //  Fix cross-overs
    //
    for (j = nSamples; j >= 2; --j) {
        for (i = 0; i < nSamples-j-1; i++) {
            vS1 = rgv[i                ];       // i doesn't wrap, 1D distance is not toroidal
            vE1 = rgv[(i+1)  % nSamples];
            vS2 = rgv[(i+j)  % nSamples];
            vE2 = rgv[(i+j+1)% nSamples];
            if (vS1.x - vE1.x > 0.5)               // segment wraps around torus
                vE1.x += 1.0;
            else if (vS1.x - vE1.x < -0.5)
                vS1.x += 1.0;
            if (vS1.y - vE1.y > 0.5)
                vE1.y += 1.0;
            else if (vS1.y - vE1.y < -0.5)
                vS1.y += 1.0;
            if (vS2.x - vE2.x > 0.5)
                vE2.x += 1.0;
            else if (vS2.x - vE2.x < -0.5)
                vS2.x += 1.0;
            if (vS2.y - vE2.y > 0.5)
                vE2.y += 1.0;
            else if (vS2.y - vE2.y < -0.5)
                vS2.y += 1.0;
            if (Intersect(vS1, vE1, vS2, vE2)) {    // segments cross
                //printf("segements cross j == %d\n", j);
                for (i1 = i + 1, i2 = i + j; i1 < i2; i1++, --i2) {
                    vTmp = rgv[i1%nSamples];
                    rgv[i1%nSamples] = rgv[i2%nSamples];
                    rgv[i2%nSamples] = vTmp;
                }
            }
        }
    }
    Locality(&fAve, &fMax, rgv, nSamples);
    printf("after fixing crossovers, Locality = ");
    printf("%f max %f ave\n", fMax, fAve);
    //
    //  One pass of greedy optimization
    //
    for (i = 0; i < nSamples; i++) {
        for (j = 0; j < i; j++) {
            DeltaEnergy(&dEnergyAve, &dEnergyMax, rgv, nSamples, i, j);
            if (dEnergyMax < 0.0 || (dEnergyMax == 0.0 && dEnergyAve < 0.0)) {
                vTmp = rgv[i];
                rgv[i] = rgv[j];
                rgv[j] = vTmp;
                nSucc++;
            }
        }
    }
    Locality(&fAve, &fMax, rgv, nSamples);
    printf("after greedy Locality = ");
    printf("%f max %f ave\n", fMax, fAve);

    nResult = 1;
Fail:
    delete [] pnPermute1;
    delete [] pnPermute2;
    return nResult;
}

static double
VertexPathLength(Vector2 rgv[], int n)
{
    int i;
    double f;

    f = 0.0;
    for (i = 0; i < n; i++) {
        f += Distance(rgv[i], rgv[(i+1)%n]);
    }
    return f;
}

//
//  New Experiment in localizing sequence of sampling pattern
//
// #define BLUE_NSQRT  32
#define NPOINTS     (BLUE_NSQRT*BLUE_NSQRT)
#define NEDGES      ((NPOINTS-1)*(NPOINTS)/2)
#define NIMAGE      512

struct Edge {
    float   fDist;
    short   iFirst, iSecond;

        Edge() {}
        Edge(double r, int i, int j) : fDist(float(r)), iFirst(i), iSecond(j) {}
    void    Print() { printf("<%8.6f, %4d %4d >\n", fDist, iFirst, iSecond); }
};

struct Vertex {
    Vertex      **rgpvOutEdges;
    Vertex      **rgpvInEdges;
    int         iPnt;
    int         nInDegree;
    int         nOutDegree;
    int         iIn;
    int         iOut;
    int         bTreeified;
};

static inline int
ML_IsLess(Edge e1, Edge e2)
{
    return e1.fDist < e2.fDist;
}

static Edge     *s_rge;
static Edge     *s_rgeMST;
static int       s_nMST;
static short    *s_rgiID;                   // Union-Find structure for MST algorithm
static Vector2   s_rgv[NPOINTS];
static Vertex    s_rgvVertices[NPOINTS];


static void
DrawTorus(Image &im, double x1, double y1, double x2, double y2, DisplayRGB rgb = ML_Red)
{
    x1 *= NIMAGE;
    y1 *= NIMAGE;
    x2 *= NIMAGE;
    y2 *= NIMAGE;
    if (x2 - x1 >  NIMAGE/2) x2 -= NIMAGE;
    if (x2 - x1 < -NIMAGE/2) x2 += NIMAGE;
    if (y2 - y1 >  NIMAGE/2) y2 -= NIMAGE;
    if (y2 - y1 < -NIMAGE/2) y2 += NIMAGE;
    im.DrawLineRGB(rgb, x1, y1, x2, y2, 1.5);
}

static void
DrawMST()
{
    Image im;
    Edge e;
    int i;

    im.NewImage(NIMAGE, NIMAGE, 3);
    im.m_bWrapX = WRAP_AROUND;
    im.m_bWrapY = WRAP_AROUND;
    for (i = 0; i < NPOINTS; i++) {
        im.FastSplatRGB(ML_White, NIMAGE * s_rgv[i].x, NIMAGE * s_rgv[i].y);
    }
    for (i = 0; i < s_nMST; i++) {
        e = s_rgeMST[i];
        DrawTorus(im, s_rgv[e.iFirst].x, s_rgv[e.iFirst].y, s_rgv[e.iSecond].x, s_rgv[e.iSecond].y);
    }
    im.WriteBMP("MST.bmp");
}

static int
MST_Initialize(Vector2 rgv[], int nPoints)
{
    int i, j, n;
    double r;

    s_rge = new Edge[NEDGES];               // Minimum spanning tree in fully-connected graph of sampling pattern
    s_rgeMST = new Edge[NEDGES];
    s_nMST = 0;
    s_rgiID = new short[NEDGES];
    n = 0;
    for (i = 0; i < NPOINTS; i++) {
        for (j = i+1; j < NPOINTS; j++) {
            r = Distance(s_rgv[i], s_rgv[j]);
            s_rge[n] = Edge(r, i, j);
            s_rgiID[n] = n;
            n++;
        }
    }
    ML_QuickSort(s_rge, n);
    printf("NEDGES = %d, n = %d\n", NEDGES, n);
    return n;
}

static int
Find(int i)
{
    while (s_rgiID[i] != i) {
        s_rgiID[i] = s_rgiID[s_rgiID[i]];
        i = s_rgiID[i];
    }
    return i;
}

static void
Union1(int i, int j)
{
    int p, q;

    p = Find(i);
    q = Find(j);
    s_rgiID[p] = s_rgiID[q];
}

static double
MST_Kruskal(Vector2 rgv[], int nPoints)
{
    int i, ex, ey;
    double fDist, minDist = 0.0;

    MST_Initialize(rgv, nPoints);
    s_nMST = 0;
    for (i = 0; i < NEDGES; i++) {
        ex = s_rge[i].iFirst;
        ey = s_rge[i].iSecond;
        fDist = s_rge[i].fDist;
        if (Find(ex) != Find(ey)) {
            minDist += sqrt(fDist);
            s_rgeMST[s_nMST++] = s_rge[i];
            Union1(ex, ey);
        }
    }
    return minDist;
}

static Image s_im;
static int s_bIm = 0;

static void
DrawTree(Vertex *pv, double a)
{
    DisplayRGB rgb;
    int i;
    Vector2 v, v2;

    if (a < 0.7) return;
    if (s_bIm == 0) {
        s_im.NewImage(1024, 1024, 3);
        s_im.m_bWrapX = WRAP_AROUND;
        s_im.m_bWrapY = WRAP_AROUND;
        s_bIm++;
        printf("Init treeify s_im\n");
    }
    rgb = ML_Red + ML_Cyan*a;
    v = s_rgv[pv->iPnt];
    for (i = 0; i < pv->nOutDegree; i++) {
        v2 = s_rgv[pv->rgpvOutEdges[i]->iPnt];
        DrawTorus(s_im, v.x*1024.0, v.y*1024.0, v2.x*1024, v2.y*1024.0, rgb);
        DrawTree(pv->rgpvOutEdges[i], a*0.9);
    }
}

static Vertex*
Treeify(Vertex *vRoot, Vertex *vParent)
{
    int i, j;
    Vertex *pvIn, *pvOut;

    if (vRoot->bTreeified)
        return vRoot;
    vRoot->bTreeified = 1;
    for (i = 0; i < vRoot->nInDegree; i++) {          // move unique incoming edges to outgoing list
        pvIn = vRoot->rgpvInEdges[i];
        if (pvIn == vParent)
            continue;
        for (j = 0; j < vRoot->nOutDegree; j++) {
            pvOut = vRoot->rgpvOutEdges[j];
            if (pvIn == pvOut)
                goto BreakContinue;
        }
        vRoot->rgpvOutEdges[vRoot->nOutDegree++] = pvIn;
    BreakContinue: continue;
    }
    vRoot->rgpvInEdges[0] = vParent;
    if (vParent)
        vRoot->nInDegree = 1;
    else
        vRoot->nInDegree = 0;
    for (j = 0; j < vRoot->nOutDegree; j++) {
        vRoot->rgpvOutEdges[j] = Treeify(vRoot->rgpvOutEdges[j], vRoot);
    }
    return vRoot;
}

static Vertex*
BuildVertex()      
{
    int i, i1, i2, nIn, nOut;

    for (i = 0; i < NPOINTS; i++) {
        s_rgvVertices[i].nOutDegree = 0;
        s_rgvVertices[i].nInDegree = 0;
        s_rgvVertices[i].iIn = 0;
        s_rgvVertices[i].iOut = 0;
        s_rgvVertices[i].iPnt = i;          // the vertext represents the point at s_rgv[i]
        s_rgvVertices[i].bTreeified = 0;
    }
    for (i = 0; i < s_nMST; i++) {
        s_rgvVertices[s_rgeMST[i].iFirst].nOutDegree++;    // set the outdegree of edges, and allocate edge list
        s_rgvVertices[s_rgeMST[i].iSecond].nInDegree++;
    }
    nIn = nOut = 0;
    for (i = 0; i < NPOINTS; i++) {
        s_rgvVertices[i].rgpvOutEdges = new Vertex*[s_rgvVertices[i].nOutDegree + s_rgvVertices[i].nInDegree];
        s_rgvVertices[i].rgpvInEdges = new Vertex*[s_rgvVertices[i].nInDegree + 1];  
    }
    //
    //  undirected MST from Kruskal --> bi-directional graph
    //
    for (i = 0; i < s_nMST; i++) {
        i1 = s_rgeMST[i].iFirst;
        i2 = s_rgeMST[i].iSecond;
        s_rgvVertices[i1].rgpvOutEdges[s_rgvVertices[i1].iOut++] = &s_rgvVertices[i2];
        s_rgvVertices[i2].rgpvInEdges[s_rgvVertices[i2].iIn++] = &s_rgvVertices[i1];
    }
    //
    //  Turn bidirectional graph into an oriented tree
    //
    printf("Call Treeify\n");
    return Treeify(&s_rgvVertices[0], 0);
}

//
//  fBlue = 84.325130       // total path length
//  fOpt = 60.136788
//  fMST = 175.061970
//
void
TestNewLocality()
{
    double fMST, fBlue, fOpt;
    Vertex *pv;

    BlueNoisePattern(s_rgv, BLUE_NSQRT);
    fBlue = VertexPathLength(s_rgv, BLUE_NSQRT*BLUE_NSQRT);
    // OptimizeLocality(s_rgv, BLUE_NSQRT*BLUE_NSQRT);
    fOpt = VertexPathLength(s_rgv, BLUE_NSQRT*BLUE_NSQRT);
    fMST = MST_Kruskal(s_rgv, NPOINTS);
    s_rge[0].Print();
    s_rge[1].Print();
    s_rge[2].Print();
    s_rge[3].Print();
    s_rge[500000].Print();
    s_rge[500001].Print();
    s_rge[500002].Print();
    s_rge[500003].Print();
    DrawMST();
    printf("fBlue = %f\nfOpt = %f\nfMST = %f\n", fBlue, fOpt, fMST);

    pv = BuildVertex();
    printf("call DrawTree\n");
    DrawTree(pv, 1.0);
    s_im.WriteBMP("Tree.bmp");
}
//
//  1D patterns
//
int
UniformPattern(double rgf[], int nSamples)
{
    double fScale;
    int i;

    fScale = 1.0/double(nSamples);
    for (i = 0; i < nSamples; i++) {
        rgf[i] = fScale*(double(i) + 0.5);
    }
    return nSamples;
}

int
JitterPattern(double rgf[], int nSamples)
{
    double fScale;
    int i;

    fScale = 1.0/double(nSamples);
    for (i = 0; i < nSamples; i++) {
        rgf[i] = fScale*(double(i) + RandomDouble());
    }
    return nSamples;
}
//
//  Disc pattern for area lights
//
#define EPSILON 0.0000000001

static double
SegmentArea(Vector2 v1, Vector2 v2)
{
    double h, fArea;
    Vector2 vN;

    vN = Normalize(Vector2(v2.y - v1.y, v1.x - v2.x));
    h = -Dot(vN, v1);
    if (h < 0.0) {
        h = -h;
        fArea = D_PI - (acos(h) - h*sqrt(1 - h*h));
    } else
        fArea = acos(h) - h*sqrt(1 - h*h);
    return fArea;
}

static double
SegmentCount(Vector2 v1, Vector2 v2, Vector2 rgv[], int n)
{
    Vector2 vN;
    int i;
    double fCount;

    vN = Vector2(v2.y - v1.y, v1.x - v2.x);     // need not be normalized for counting
    fCount = 0.0;
    for (i = 0; i < n; i++) {
        if(Dot(vN, v1 - rgv[i]) > EPSILON)
            fCount += 1.0;
    }
    return fCount;
}

static double
DiskDiscrepancy(Vector2 rgv[], int n)
{
    int i, j;
    double fArea, fCount, fDisc, fMaxDisc;

    fMaxDisc = 0.0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < i; j++) {
            fArea = SegmentArea(rgv[i], rgv[j]);
            fCount = SegmentCount(rgv[i], rgv[j], rgv, n);
            fDisc = fabs(fArea/D_PI - (fCount + 1)/n);
            if (fDisc > fMaxDisc)
                fMaxDisc = fDisc;
        }
    }
    return fMaxDisc;
}

struct PointAngle {
    Vector2 v;
    double  fAng;
};

static inline int
ML_IsLess(const PointAngle &a1, const PointAngle &a2)
{
    return a1.fAng < a2.fAng;
}

double
DiamondAngle(double y, double x)
{
    if (y >= 0)
        return (x >= 0 ? y/(x+y) : 1-x/(-x+y)); 
    else
        return (x < 0 ? 2-y/(-x-y) : 3+x/(x-y)); 
}

void
QuickSort(PointAngle rgtItem[], int nItems)
{
    int i, j;
    PointAngle tTemp, tPivot;

    nItems--;   // iLast
    while (nItems >= 30) {
        i = nItems/2;
        //
        //  Sort first, middle and last elements.  This provides median-of-3
        //  partitioning, limits partitioning to only N - 3 remaining items,
        //  and creates sentinals to simplify the inner loop.
        //
        if (ML_IsLess(rgtItem[i], rgtItem[0])) {      // 2.48 compares on average
            tTemp = rgtItem[0];
            rgtItem[0] = rgtItem[i];
            rgtItem[i] = tTemp;
        }
        if (ML_IsLess(rgtItem[nItems], rgtItem[i])) {
            tTemp = rgtItem[nItems];
            rgtItem[nItems] = rgtItem[i];
            if (ML_IsLess(tTemp, rgtItem[0])) {
                rgtItem[i] = rgtItem[0];
                rgtItem[0] = tTemp;
            } else
                rgtItem[i] = tTemp;
        }
        j = nItems - 1;
        tPivot = rgtItem[i];
        rgtItem[i] = rgtItem[j];
        rgtItem[j] = tPivot;
        i = 0;
        //
        //  Partition, using Sedgewick's "j < i" suggestion.  Oddly, it is
        //  faster to loop on i before looping on j (on the Pentium 4).
        //
        for(;;) {
            while(ML_IsLess(rgtItem[++i], tPivot))
                ;
            while(ML_IsLess(tPivot, rgtItem[--j]))
                ;
            if (j < i)
                break;
            tTemp = rgtItem[i];
            rgtItem[i] = rgtItem[j];
            rgtItem[j] = tTemp;
        }
        tTemp = rgtItem[nItems - 1];
        rgtItem[nItems - 1] = rgtItem[i];
        rgtItem[i] = tTemp;
        //
        //  Recursing on smaller partition yields O(log N) stack growth.
        //
        if (j < nItems - i - 1) {
            QuickSort(rgtItem ,j + 1);
            rgtItem += i + 1;
            nItems -= i + 1;
        } else {
            QuickSort(rgtItem + i + 1, nItems - i);
            nItems = j;
        }
    }
    //
    //  Small partitions are insertion sorted.  Distribution is wedge
    //  shaped, with only about 3.8 comparisons done on average, and
    //  benefit gained from structuring the loop for quick out.
    //
    for (i = 1; i <= nItems; i++) {
        j = i;
        tTemp = rgtItem[j];
        if (ML_IsLess(tTemp,rgtItem[j - 1])) {
            do {
                rgtItem[j] = rgtItem[j - 1];
                j = j - 1;
            } while (j > 0 && ML_IsLess(tTemp, rgtItem[j - 1]));    
            rgtItem[j] = tTemp;
        }
    }
}

static void
AddPoint(Vector2 rgv[], int n, int nTry)
{
    double fDisc, fMinDisc;
    Vector2 v, vMin;
    int j;

    fMinDisc = 1000000.0;
    for (j = 0; j < nTry; j++) {
        RandomInsideDisk(v.rgf);
        rgv[n - 1] = v;
        fDisc = DiskDiscrepancy(rgv, n);
        if (fDisc < fMinDisc) {
            fMinDisc = fDisc;
            vMin = v;
        }
    }
    rgv[n-1] = vMin;
}

static double
PerturbPoints(Vector2 rgv[], int n, float fDelta)
{
    double fDisc, fDiscNew;
    Vector2 vDelta, vNew, vOld;
    int i, j;

    fDisc = DiskDiscrepancy(rgv, n);
    for (j = 0; j < 500; j++) {
        for (i = 0; i < n; i++) {
            vOld = rgv[i];
            RandomInsideDisk(vDelta.rgf);
            vDelta = vDelta * fDelta;
            rgv[i] = rgv[i] + vDelta;
            if (Dot(rgv[i], rgv[i]) >= 1.0) {
                rgv[i] = vOld;
                continue;
            }
            fDiscNew = DiskDiscrepancy(rgv, n);
            if (fDiscNew < fDisc) {
                fDisc = fDiscNew;
                continue;
            }
            rgv[i] = vOld;
        }
    }
    return fDisc;
}

static void
Annealing(Vector2 rgv[], int n)
{
    int i;

    for (i = 0; i < 20; i++) {
        printf("Anneal: %f\n", PerturbPoints(rgv, n, 0.1f));
    }
    for (i = 0; i < 20; i++) {
        printf("Anneal: %f\n", PerturbPoints(rgv, n, 0.03f));
    }
    for (i = 0; i < 20; i++) {
        printf("Anneal: %f\n", PerturbPoints(rgv, n, 0.01f));
    }
    for (i = 0; i < 20; i++) {
        printf("Anneal: %f\n", PerturbPoints(rgv, n, 0.003f));
    }
    for (i = 0; i < 20; i++) {
        printf("Anneal: %f\n", PerturbPoints(rgv, n, 0.001f));
    }
}

#define ISIZE 720

static void
DrawLine(Image &im, DisplayRGB rgb, Vector2 v1, Vector2 v2)
{
    im.DrawLineRGB(rgb, double(ISIZE/2)*v1.x + double(ISIZE/2), double(ISIZE/2)*v1.y + double(ISIZE/2),
                        double(ISIZE/2)*v2.x + double(ISIZE/2), double(ISIZE/2)*v2.y + double(ISIZE/2), 1.0);
}

static void
Splat(Image &im, DisplayRGB rgb, Vector2 v)
{
    im.SplatRGB(rgb, double(ISIZE/2)*v.x + double(ISIZE/2), double(ISIZE/2)*v.y + double(ISIZE/2), 0, 2.0);
}

static void
DrawDiskPattern(Vector2 rgv[], int n)
{
    Image im;
    Vector2 v1, v2;
    double fAng;
    int i;

    im.NewImage(ISIZE, ISIZE, 3);
    for (fAng = 0.0; fAng <= 360.0; fAng += 1.0) {
        v2 = Vector2(sin(fAng * D_PI/180.0), cos(fAng * D_PI/180.0));
        if (fAng > 0.0)
            DrawLine(im, ML_White, v1, v2);
        v1 = v2;
    }
    for (i = 0; i < n; i++)
        Splat(im, ML_White, rgv[i]);
    im.WriteBMP("disk.bmp");
}

static void
DrawDiskPattern2(Vector2 rgv[], int n, Vector2 u1, Vector2 u2)
{
    Image im;
    Vector2 v1, v2, vN;
    double fAng;
    int i;

    im.NewImage(ISIZE, ISIZE, 3);
    for (fAng = 0.0; fAng <= 360.0; fAng += 1.0) {
        v2 = Vector2(sin(fAng * D_PI/180.0), cos(fAng * D_PI/180.0));
        if (fAng > 0.0)
            DrawLine(im, ML_White, v1, v2);
        v1 = v2;
    }
    DrawLine(im, ML_Green, u1, u2);
    vN = Vector2(u2.y - u1.y, u1.x - u2.x); 
    for (i = 0; i < n; i++) {
        if(Dot(vN, u1 - rgv[i]) > EPSILON)
            Splat(im, ML_Red, rgv[i]);
        else
            Splat(im, ML_White, rgv[i]);
    }
    im.WriteBMP("disk.bmp");
}

static void
DrawDiskPattern3(Vector2 rgv[], int n)
{
    Image im;
    Vector2 v1, v2;
    double fAng;
    int i;

    im.NewImage(ISIZE, ISIZE, 3);
    for (fAng = 0.0; fAng <= 360.0; fAng += 1.0) {
        v2 = Vector2(sin(fAng * D_PI/180.0), cos(fAng * D_PI/180.0));
        if (fAng > 0.0)
            DrawLine(im, ML_White, v1, v2);
        v1 = v2;
    }
    for (i = 0; i < n-1; i++)
        DrawLine(im, ML_Red, rgv[i], rgv[i+1]);
    for (i = 0; i < n; i++)
        Splat(im, ML_White, rgv[i]);
    im.WriteBMP("disk.bmp");
}

void
GenerateDiskPattern(int n)
{
    static Vector2 rgv[10000], v1, v2;
    int i;
    FILE *pf;

    printf("generate %d\n", n);
    fopen_s(&pf, "DiskSamples.cpp", "w");
    fprintf(pf, "Vector2\nDiskSamples%d[%d] = {\n", n, n);
    rgv[0] = Vector2(0.0, 0.0);
    for (i = 2; i < n; i++)
        AddPoint(rgv, i, 10*(i+2));
    Annealing(rgv, n);
    OptimizeLocality(rgv, n, 0);    // no wrap-around
    for (i = 0; i < n; i++) {
        fprintf(pf, "    Vector2(%f, %f),\n", rgv[i].x, rgv[i].y);
    }
    fprintf(pf, "};\n");
    fclose(pf);
    DrawDiskPattern3(rgv, n);
}

void
TestDiskDiscrepancy()
{
    static Vector2 rgv[10000], v1, v2;
    int i;
    FILE *pf;
    //
    //  Generate the tables of low discrepancy disk pattersn
    //
    GenerateDiskPattern(128);
    return;
    fopen_s(&pf, "DiskSamples.cpp", "w");
    rgv[0] = Vector2(0.0, 0.0);
    for (i = 2; i < 32; i++)
        AddPoint(rgv, i, 10*(i+2));
    Annealing(rgv, 32);
    for (i = 0; i < 32; i++) {
        rgv[i].Print();
        printf(" - %d\n", i);
        fprintf(pf, "Vector2(%f, %f),\n", rgv[i].x, rgv[i].y);
    }
    fclose(pf);
    //DrawDiskPattern2(rgv, i, rgv[2], rgv[5]);
    return;

    for (i = 0; i < 1000000; i++)
        RandomInsideDisk(rgv[i].rgf);
    for (i = 0; i < 20; i++) {
        RandomInsideDisk(v1.rgf);
        RandomInsideDisk(v2.rgf);
        //v1 = Vector2(0.0, 0.0);
        printf("Area = %f, Counted = %f\n", SegmentArea(v1, v2)/D_PI, SegmentCount(v1, v2, rgv, 1000000) /1000000.0);
    }
    printf("\n");
    DrawDiskPattern(rgv, 100);
    for (i = 2; i <= 1000000; i *= 2)
        printf("Max Discrepancy = %f\n", DiskDiscrepancy(rgv, i));
}
//
//  Quasi-Random Utility Functions
//
double
ML_RadicalInverse(unsigned n, unsigned nPrime)
{
	double invRadix, invPower, value;
	unsigned nDigit, nNext;

	invRadix = 1.0/double(nPrime);
	value = 0.0;
	invPower = 1.0;
    while (n) {
        invPower *= invRadix;
        nNext = n/nPrime;
        nDigit = n - nNext*nPrime;
        n = nNext;
        value += invPower*double(nDigit);
    }
	return value;
}

double
ML_FoldedRadicalInverse(unsigned n, unsigned nPrime)
{
	double invRadix, invPower, value;
	int i;
    unsigned nDigit;

	invRadix = 1.0/double(nPrime);
	value = 0;
	invPower = 1.0;
	for (i = 0; value + invPower != value; i++) {
		invPower *= invRadix;
		nDigit = (n + i) % nPrime;
		n /= nPrime;
		value += invPower*nDigit;
	}
	return value;
}

#define NPERMUTE 256

struct NumDigit {
    int N, nDigit;
};

static NumDigit rgnd[NPERMUTE];

static int
Digit(int n, int ith, int nBase)
{
    int i;

    for (i = 0; i < ith; i++)
        n = n / nBase;
    return n % nBase;
}

static void
InsertionSort(NumDigit rgtItem[], int nStart, int nLimit)
{
    int i, j;
    NumDigit tKey;

    for (j = nStart + 1; j < nLimit; j++) {
        tKey = rgtItem[j];
        for (i = j; i > nStart && (tKey.nDigit < rgtItem[i - 1].nDigit); i = i - 1)
            rgtItem[i] = rgtItem[i - 1];    // faster than [i+1] = [i]
        rgtItem[i] = tKey;
    }
}
//
//  Permute an array of numbers [0, n-1] into a radical-inverse order
//
static void
QuasiShuffle(NumDigit rgnd[], int ith, int nBase, int nStart, int nLimit)
{
    int i, iStart, iLimit, nDigit;

    if (nLimit - nStart < 2)
        return;
    for (i = nStart; i < nLimit; i++)
        rgnd[i].nDigit = Digit(rgnd[i].N, ith, nBase);
    InsertionSort(rgnd, nStart, nLimit);
    iStart = nStart;
    nDigit = rgnd[iStart].nDigit;
    for (iLimit = iStart + 1; iLimit <= nLimit; ) {
        if (rgnd[iLimit].nDigit != nDigit || iLimit == nLimit) {
            QuasiShuffle(rgnd, ith+1, nBase, iStart, iLimit);
            iStart = iLimit;
            iLimit = iStart+1;    
            nDigit = rgnd[iStart].nDigit;
            continue;
        }
        iLimit++;
    }
}

void
BuildQuasiIndex(int rgn[], int nPrime, int N)
{
    NumDigit *pnd;
    int i;

    pnd = new NumDigit[N];
    for (i = 0; i < N; i++)
        pnd[i].N = i;
    QuasiShuffle(pnd, 0, nPrime, 0, N);
    for (i = 0; i < N; i++)
        rgn[i] = pnd[i].N;
    delete [] pnd;
}

#define NB 5
#define NN 125

static void
DigitPrint(int N)
{
    printf("%d %d %d %d %d   ", Digit(N, 4, NB), Digit(N, 3, NB), Digit(N, 2, NB), Digit(N, 1, NB), Digit(N, 0, NB));
}

static void
PrintQuasiIndex(int nPrime, int N)
{
    int i, *rgn;

    rgn = new int[N];
    BuildQuasiIndex(rgn, nPrime, N);
    printf("\nint rgnQuasiIndex%d[%d] = {\n", nPrime, N);
    for (i = 0; i < N; i++) {
        printf("%4d,", rgn[i]);
        if (i%8 == 7)
            printf("\n");
    }
    printf("};\n");
    delete [] rgn;
}

static void
FPrintQuasiIndex(FILE *pf, int nPrime, int N)
{
    int i, *rgn;

    rgn = new int[N];
    BuildQuasiIndex(rgn, nPrime, N);
    fprintf(pf, "\nint g_BlueNoiseQuasiIndex%d[%d] = {\n", nPrime, N);
    for (i = 0; i < N; i++) {
        fprintf(pf, "%4d,", rgn[i]);
        if (i%8 == 7)
            fprintf(pf, "\n");
    }
    fprintf(pf, "};\n");
    delete [] rgn;
}
//
//  Write out samples
//
int
WritePattern(char *szFileName, char *szPatternName, Vector2 rgv[], int nSqrt)
{
    FILE *pf;
    int i;

    fopen_s(&pf, szFileName, "w");
    fprintf(pf, "#include \"stdafx.h\"\n");
    fprintf(pf, "#include \"Affine.h\"\n\n");
    fprintf(pf, "Vector2 %s[%d] = {\n", szPatternName, nSqrt*nSqrt);
    for (i = 0; i < nSqrt*nSqrt; i++)
        fprintf(pf, "    Vector2(%0.6f, %0.6f),\n", rgv[i].x, rgv[i].y);
    fprintf(pf, "};\n");
  

    FPrintQuasiIndex(pf, 2, BLUE_NSQRT*BLUE_NSQRT);
    FPrintQuasiIndex(pf, 3, BLUE_NSQRT*BLUE_NSQRT);
    FPrintQuasiIndex(pf, 5, BLUE_NSQRT*BLUE_NSQRT);
    FPrintQuasiIndex(pf, 7, BLUE_NSQRT*BLUE_NSQRT);
    fclose(pf);
    return 1;
}

void
MakeBlueNoiseTile()
{
    Vector2 *rgv;

    rgv = new Vector2[BLUE_N];
    printf("Make Blue Noise Tile: %d x %d\n", BLUE_NSQRT, BLUE_NSQRT);
    // BlueNoisePattern(rgv, BLUE_NSQRT);
    // OptimizeLocality(rgv, BLUE_NSQRT*BLUE_NSQRT);
    BlueJitterPattern(rgv, BLUE_NSQRT);
    WritePattern("BlueNoiseTile.cpp", "g_rgvBlueNoise", rgv, BLUE_NSQRT);
    delete [] rgv;
}


void
TestDiscrepancy()
{
    NumDigit rgnd[256];
    int i;

    //PrintQuasiIndex(7, 1024); return;

    PrintQuasiIndex(2, 128);
    PrintQuasiIndex(3, 128);
    PrintQuasiIndex(5, 128);
    PrintQuasiIndex(7, 128);
    return;

    for (i = 0; i < NN; i++)
        rgnd[i].N = i;
    QuasiShuffle(rgnd, 0, NB, 0, NN);
    for (i = 0; i < NN; i++) {
        printf("%3d   ", i);
        DigitPrint(i);
        DigitPrint(rgnd[i].N),
        printf("  %3d  %6.3f  %6.3f\n", rgnd[i].N, ML_RadicalInverse(i, NB)*NN, ML_FoldedRadicalInverse(i, NB)*NN);
    }
}
//
//  Fibonacci spiral on the sphere
//
double g_fPhi = D_PHI;

void
FibonacciSphere(double rgf[], int i, int N)
{
    double u, v;
    double x, y, z, r;

    u = fmod(double(i)/g_fPhi, 1.0);
    v = double(i)/double(N);
    z = 2.0*v - 1.0;
    r = sqrt(1.0 - z*z);
    x = r*sin(2.0*D_PI*u);
    y = r*cos(2.0*D_PI*u);
    rgf[0] = (x);
    rgf[1] = (y);
    rgf[2] = (z);
}
//
//  Equal-area mappings  (
//
static void
PolarToSquare(double fRadius, double fAngle, Vector2 *ppntI2)
{
    double x, y;

    fAngle *= 2.0 / D_PI;
    fAngle += 1.5;
    if (fAngle > 2.0)
        fAngle -= 4.0;
    y = fabs(fAngle) - 1.0;
    if (y < 0.0f)
        x = fAngle;
    else if (fAngle > 0.0)
        x = 2.0 - fAngle;
    else
        x = -2.0 - fAngle;
    ppntI2->x = (fRadius * 0.5 * (y + x) + 0.5);
    ppntI2->y = (fRadius * 0.5 * (y - x) + 0.5);
}

Vector2
DiskToUnitSquare(const Vector2 &pntD)
{
    Vector2 pntI2;
    double fRadius, fAngle;

    fRadius = Norm(pntD);
    fAngle = atan2(pntD.y, pntD.x);
    PolarToSquare(fRadius, fAngle, &pntI2);
    return pntI2;
}

static inline void
SquareToPolar(const Vector2 &pntI2, double *pfRadius, double *pfAngle)
{
    double x, y, r, a;

    y = pntI2.x + pntI2.y - 1.0;
    x = pntI2.x - pntI2.y;
    *pfRadius = r = fabs(y) + fabs(x);    // 0.0 to 1.0
    if (r == 0.0f)
        a = 0.0;
    else if (x >= 0.0)
        a = 0.5 * (y/r + 1.0);            // 0.0 to pi
    else
        a = 0.5 * (-1.0 - y/r);           // 0.0 to -pi
    a -= 0.75;                             // undo the implicit rotation
    if (a < -1.0)
        a += 2.f;
    *pfAngle = D_PI * a;
}

static Vector2
SquareToUnitDisk_Old(const Vector2 &pntI2)
{
    double fRadius, fAngle;

    SquareToPolar(pntI2, &fRadius, &fAngle);
    return Vector2(fRadius * cos(fAngle), fRadius * sin(fAngle));
}

Vector3
DiskToUnitHemisphere(const Vector2 &pntD)
{
    double c;
    double s;

    c = 1.0 - pntD.x*pntD.x - pntD.y*pntD.y;
    s = (sqrt(1.0 + c));
    return Vector3(s*pntD.x, s*pntD.y, (c));
}

Vector3
SquareToUnitSphere(const Vector2 &pntI2)
{
    double fRadius, fAngle, zSign;
    Vector2 pntI2rolled, pntD;
    Vector3 pntS2;

    //
    //  roll the edge of the square to the center, so the inverse image of the north
    //  and south poles will touch.
    //
    pntI2rolled = pntI2;
    pntI2rolled.y += 0.5f;
    if (pntI2rolled.y >= 1.0f)
        pntI2rolled.y -= 1.0f;
    //
    //  Map half the disk to the north hemisphere, half to the south
    //
    SquareToPolar(pntI2rolled, &fRadius, &fAngle);
    if (fAngle < 0.0) {
        zSign = -1.0f;
        fAngle = +2.0f * fAngle;
    } else {
        zSign = 1.0f;
        fAngle = -2.0f * fAngle;    // equator matches up on the square
    }
    pntD = Vector2(fRadius * cos(fAngle), fRadius * sin(fAngle)); //REVIEW: r*sin(x/r) for small r
    pntS2 = DiskToUnitHemisphere(pntD);
    pntS2.z *= zSign;
    return pntS2;
}

Vector2
HemisphereToUnitDisk(const Vector3 &pntH)
{
    double s;

    s = 1.0 / (sqrt(1.0 + double(pntH.z)));
    return Vector2(s * pntH.x, s * pntH.y);
}

Vector2
SphereToUnitSquare(const Vector3 &pntS2)
{
    double zSign, fRadius, fAngle;
    Vector2 pntI2, pntD;
    Vector3 pntH;

    pntH = pntS2;
    if (pntH.z < 0.0) {
        zSign = -1.0;
        pntH.z = -pntH.z;
    } else
        zSign = 1.0;
    pntD = HemisphereToUnitDisk(pntH);
    fRadius = Norm(pntD);
    fAngle = -atan2(pntD.y, pntD.x);
    if (fAngle < 0.0)
        fAngle = 2.0*D_PI + fAngle;    // (-pi, pi) -> (0, 2pi)
    if (zSign < 0.0)
        fAngle = - 0.5*fAngle;
    else
        fAngle = 0.5*fAngle;

    PolarToSquare(fRadius, fAngle, &pntI2);
    pntI2.y -= 0.5f;
    if (pntI2.y < 0.0)
        pntI2.y += 1.0;
    return pntI2;
}

void 
Shirley_to_unit_disk( double seedx, double seedy, double *x, double *y )
{
    double phi, r;
    double a = 2*seedx - 1;         /* (a,b) is now on [-1,1]^2 */
    double b = 2*seedy - 1;

    if (a > -b) {                   /* region 1 or 2 */
        if (a > b) {                    /* region 1, also |a| > |b| */
            r = a;
            phi = (D_PI/4 ) * (b/a);
        }
        else {                          /* region 2, also |b| > |a| */
            r = b;  // ???
            phi = (D_PI/4) * (2 - (a/b));
        }
    } else {                          /* region 3 or 4 */
        if (a < b)  {                    /* region 3, also |a| >= |b|, a != 0 */
            r = -a;
            phi = (D_PI/4) * (4 + (b/a));
        } else {                          /* region 4, |b| >= |a|, but a==0 and b==0 could occur. */
            r = -b;
            if (b != 0)
                phi = (D_PI/4) * (6 - (a/b));
            else
                phi = 0;
        }
    }
    *x = r * cos(phi);
    *y = r * sin(phi);
}

Vector2 
SquareToUnitDisk(const Vector2 &O)     // Dave Cline's code is much faster than mine or Shirley's
{
    double phi,r;
    double a = 2*O.x - 1;
    double b = 2*O.y - 1;

    if (a*a> b*b) {             // use squares instead of absolute values
        r = a;
        phi = (D_PI/4)*(b/a);
    } else {
        r = b;   
        phi = (D_PI/2) - (D_PI/4)*(a/b);             // Franz' fix
    }
    return Vector2( r*cos(phi), r*sin(phi) );
}

extern void sncndn ( double u, double m, double& sn, double& cn, double& dn, double& err);

void
TestAreaMaps()
{
    Image imSqr, imDsk;
    int i,  j, W, H;
    Vector2 vSqr, vDsk;
    double f, sn, cn, dn, err;
    ML_TimingInfo ti;

    imSqr.ReadBMP("ImageLady.bmp");
    imDsk.NewImage(1536, 1536, 1);
    W = imSqr.m_nWidth;
    H = imSqr.m_nHeight;

    ML_StartTiming(ti);
    for (i = 0; i < W; i++) {
        vSqr.x = (double(i) + 0.5)/W;
        for (j = 0; j < H; j++) {
            f = imSqr.Get(i, j) * 0.25;         // 4 times smaller area than input image
            vSqr.y = (double(j) + 0.5)/H;
            vDsk = vSqr;
            vDsk = vDsk*768.0 + Vector2(768.0, 768.0);
            imDsk.FastSplat(f, vDsk.x, vDsk.y);
        }
    }
    ML_StopTiming(ti);
    ML_ReportTiming(ti);

    imDsk.Fill(0.0);
    ML_StartTiming(ti);
    for (i = 0; i < W; i++) {
        vSqr.x = (double(i) + 0.5)/W;
        for (j = 0; j < H; j++) {
            f = imSqr.Get(i, j) * 0.25;         // 4 times smaller area than input image
            vSqr.y = (double(j) + 0.5)/H;
            vDsk = SquareToUnitDisk_Old(vSqr);
            vDsk = vDsk*768.0 + Vector2(768.0, 768.0);
            imDsk.FastSplat(f, vDsk.x, vDsk.y);
        }
    }
    ML_StopTiming(ti);
    ML_ReportTiming(ti);
    imDsk.WriteBMP("ML_SquareDisk.bmp");

    imDsk.Fill(0.0);
    ML_StartTiming(ti);
    for (i = 0; i < W; i++) {
        vSqr.x = (double(i) + 0.5)/W;
        for (j = 0; j < H; j++) {
            f = imSqr.Get(i, j) * 0.25;         // 4 times smaller area than input image
            vSqr.y = (double(j) + 0.5)/H;
            Shirley_to_unit_disk(vSqr.x, vSqr.y, &vDsk.x, &vDsk.y);
            vDsk = vDsk*768.0 + Vector2(768.0, 768.0);
            imDsk.FastSplat(f, vDsk.x, vDsk.y);
        }
    }
    ML_StopTiming(ti);
    ML_ReportTiming(ti);
    imDsk.WriteBMP("SH_SquareDisk.bmp");

    imDsk.Fill(0.0);
    ML_StartTiming(ti);
    for (i = 0; i < W; i++) {
        vSqr.x = (double(i) + 0.5)/W;
        for (j = 0; j < H; j++) {
            f = imSqr.Get(i, j) * 0.25;         // 4 times smaller area than input image
            vSqr.y = (double(j) + 0.5)/H;
            vDsk = SquareToUnitDisk(vSqr);
            vDsk = vDsk*768.0 + Vector2(768.0, 768.0);
            imDsk.FastSplat(f, vDsk.x, vDsk.y);
        }
    }
    ML_StopTiming(ti);
    ML_ReportTiming(ti);
    imDsk.WriteBMP("CL_SquareDisk.bmp");

    sncndn(0.5, 0.5, sn, cn, dn, err);
    printf("sncndn = %f %f %f, %f\n", sn, cn, dn, err);
    imDsk.Fill(0.0);
    ML_StartTiming(ti);
    for (i = 0; i < W; i++) {
        vSqr.x = (double(i) + 0.5)/W;
        for (j = 0; j < H; j++) {
            f = imSqr.Get(i, j) * 0.25;         // 4 times smaller area than input image
            vSqr.y = (double(j) + 0.5)/H;
            vDsk = SchwarzChristoffel(vSqr);
            vDsk = vDsk*768.0 + Vector2(768.0, 768.0);

            imDsk.FastSplat(f, vDsk.x, vDsk.y);
        }
    }
    ML_StopTiming(ti);
    ML_ReportTiming(ti);
    imDsk.WriteBMP("SC_SquareDisk.bmp");
}
//
//  Jacobi elliptic function (note m = k**2)
//
#define KE      1.854074677301371918434
const double  eps = 2.220446049e-16;
const int  Nmax = 16;

 void
sncndn ( double u, double m, double& sn, double& cn, double& dn, double& err) {
    double sqrt_eps, m1, t=0., si_u, co_u, se_u, ta_u, b, c[Nmax], a[Nmax], phi;
    int n, Nn, ii;

  if (m < 0. || m > 1.) {
     printf ("ellipj: expecting 0. <= m <= 1."); 
     return;
	}
  sqrt_eps = sqrt(eps);
  if (m < sqrt_eps) {
    /*  # For small m, ( Abramowitz and Stegun, Section 16.13 ) */
        si_u = sin(u);
        co_u = cos(u);
        t = 0.25*m*(u-si_u*co_u);
        sn = si_u - t * co_u;
        cn = co_u + t * si_u;
        dn = 1.0 - 0.5*m*si_u*si_u;
  } else if ( (1.0 - m) < sqrt_eps ) {
    /*  For m1 = (1-m) small ( Abramowitz and Stegun, Section 16.15 ) */
        m1 = 1.0-m;
        si_u = sinh(u);
        co_u = cosh(u);
        ta_u = tanh(u);
        se_u = 1.0/co_u;
        sn = ta_u + 0.25*m1*(si_u*co_u-u)*se_u*se_u;
        cn = se_u - 0.25*m1*(si_u*co_u-u)*ta_u*se_u;
        dn = se_u + 0.25*m1*(si_u*co_u+u)*ta_u*se_u;
  } else {
        /*
        //  Arithmetic-Geometric Mean (AGM) algorithm
        //    ( Abramowitz and Stegun, Section 16.4 )
        */
       
        a[0] = 1.0;
        b    = sqrt(1.0-m);
        c[0] = sqrt(m);
        for (n = 1; n<Nmax; ++n) {
          a[n] = (a[n-1]+b)/2;
          c[n] = (a[n-1]-b)/2;
          b = sqrt(a[n-1]*b);
          if ( c[n]/a[n] < eps) break; 
				}
        if ( n >= Nmax-1) {
           // fprintf(stderr, "Not enough workspace\n");
           err = 1.;
           return;
        }
        Nn = n;
        for ( ii = 1;  n>0;	ii = ii*2, --n) ; // pow(2, Nn)
        phi = ii*a[Nn]*u;
        for ( n = Nn; n > 0; --n) {
          t = phi;
          phi = (asin((c[n]/a[n])* sin(phi))+phi)/2.;
        }
        sn = sin(phi);
        cn = cos(phi);
        dn = cn/cos(t-phi);
  }
 return;
}

struct Complex {
    double re, im;
        
        Complex() {}
        Complex(double r, double i) : re(r), im(i) {}
};

inline Complex
operator *(const Complex &a, const Complex &b)
{
    return Complex(a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re);
}


 void
sncndn ( Complex& u, double m, Complex& sn, Complex& cn, Complex& dn, double& err) {
    double m1 = 1.0-m, ss1, cc1, dd1;
    double ss, cc, dd, ddd;

  sncndn( u.im, m1, ss1, cc1, dd1, err);
  if ( u.re == 0.) { 
    sn = Complex (0. , ss1/cc1);
    cn = Complex(1.0/cc1, 0.0);        
    dn = Complex(dd1/cc1, 0.0);       
  } else {
    sncndn( u.re, m, ss, cc, dd, err);
      ddd = cc1*cc1 + m*ss*ss*ss1*ss1;
      sn = Complex (ss*dd1/ddd, cc*dd*ss1*cc1/ddd); 
      cn = Complex (cc*cc1/ddd, -ss*dd*ss1*dd1/ddd);
      dn = Complex (dd*cc1*dd1/ddd, -m*ss*cc*ss1/ddd);
  }
 return;
}

Vector2
SchwarzChristoffel(const Vector2 &vSqr)
{
    Vector2 vDsk;
    double x, y, err;
    Complex xy, sn, cn, dn;

    x = 2.0*vSqr.x - 1.0;
    y = 2.0*vSqr.y - 1.0;
    xy = Complex(KE/2.0, KE/2.0) * Complex(x, y);
    xy.re = xy.re - KE;
    sncndn(xy, 0.5, sn, cn, dn, err);
    xy = Complex(sqrt(0.5), -sqrt(0.5)) * cn;
    vDsk.x = xy.re;
    vDsk.y = xy.im;
    return vDsk;
}