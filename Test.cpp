//
//  Testing routines, including the generation of vidual report images for sampling patterns and reconstruction
//
#include "stdafx.h"
#include "RayTracer2020.h"
#include "Cameras.h"
#include "Solids.h"
#include "Lights.h"
#include "Sampling.h"
#include "Video.h"
#include "Patterns.h"
#include "MonteCarlo.h"
#include "TextureMaps.h"
#include "Utilities.h"
#pragma intrinsic(fabs)

static void
Print(Ray r)
{
    printf("[%f %f %f -> %f %f %f]\n", r.pOrg.x, r.pOrg.y, r.pOrg.z,  r.vDir.x, r.vDir.y, r.vDir.z);
}

static Camera *s_pc;

Radiance
CamTestFunction(Vector2 v)
{
    Ray r;
    double t;
    Point3 p;
    Radiance rad;

    r = s_pc->Emit(v);
    t = -r.vDir.z/r.pOrg.z;
    p = r(t);
    rad = RGBtoSpectrum((int(8.0*p.x + 1000.0) ^ int(8.0*p.y + 1000.0)) & 1);
    return rad;
}

void
TestCamera()
{
    Ray r;
    PerspectiveCamera cp;
    StereographicCamera cs;
    TestSamplingTile st;
    Image im;

    printf("Test cameras:\n");
    r = cp.Emit(Vector2(0.5, 0.5));
    Print(r);
    r = cs.Emit(Vector2(0.5, 0.5));
    Print(r);

    r = cp.Emit(Vector2(0.25, 0.5));
    Print(r);
    r = cs.Emit(Vector2(0.25, 0.5));
    Print(r);

    r = cp.Emit(Vector2(0.0, 0.5));
    Print(r);
    r = cs.Emit(Vector2(0.0, 0.5));
    Print(r);
    s_pc = &cp;
    im.NewImage(512, 512, 1);
    st.Render(im, CamTestFunction, 1);
    im.WriteBMP("CamTestPerspective.bmp");
    s_pc = &cs;
    im.NewImage(512, 512, 1);
    st.Render(im, CamTestFunction, 1);
    im.WriteBMP("CamTestStereographic.bmp");
}

Solid   *pSolid;
Camera  *pc;
Vector3 vLight(sqrt(3.0)/2.0, 0.0, 1.0/2.0);

Radiance
SolidTestFunction(Vector2 v)
{
    Ray r;
    Hit *ph;
    Radiance rad;
    double f;
    int i, j;

    rad = g_sBlack;
    i = j = 0;
    //for (i = -3; i <= 3; i += 3) {
    //    for (j = -3; j <= 3; j += 3) {       // backup 8.0
            r = pc->Emit(v);
            r.pOrg.z += 16.0;
            r.pOrg.x += double(i);
            r.pOrg.y += double(j);
            ph = pSolid->Intersect(r);
            if (ph) {
                f = Normalize(ph->cvNormal) * vLight;
                if (f < 0.0)
                    f = 0.0;
                f = 0.05 + 0.95*f;
            } else
                f = 0.05/D_PI;         // irradiance from ambient light is Pi * background radiance
            rad = rad + RGBtoSpectrum(f);
            DeleteHitList(ph);
   //     }
   // }
    return rad;
}

void
TestSolidRender()
{
    TestSamplingTile st;
    Image im;
    Video vi;
    double fAng, fRad;
    Solid *ps1, *ps2, *psTorus, *rgps[9];
    int i;

    printf("TestSolidRender:\n");
    im.NewImage(512, 512, 3);
    pSolid = new Torus(0.5);
    pSolid = new Scale(0.5, 1.0, 1.0, pSolid);
    pSolid = new Rotate(0.0, 0.0, 45.0, pSolid);
    pc = new StereographicCamera;
    st.Render(im, SolidTestFunction, 1.0);
    im.WriteBMP("SolidTestStereographic.bmp");
    pc = new PerspectiveCamera;
    st.Render(im, SolidTestFunction, 1.0);
    im.WriteBMP("SolidTestPerspective.bmp");

    ps1 = new Translate(2.0, 0.0, 0.0, new Dodecahedron);
    ps2 = new Translate(-2.0, 0.0, 0.0, new Icosahedron);
    pSolid = new Union(ps1, ps2);
    st.Render(im, SolidTestFunction, 1.0);
    im.WriteBMP("SolidTestUnion.bmp");
    //ps1 = new Translate(0.5, 0.0, 0.0, new Sphere);
    //ps2 = new Translate(-0.5, 0.0, 0.0, new Cylinder);
    ps1 = new Hyperboloid1(0.05);
    ps2 = new Cube;
    pSolid = new Rotate(0.0, D_PI/4.0, 0.0, new Intersection(ps1, ps2));
    st.Render(im, SolidTestFunction, 1.0);
    im.WriteBMP("SolidTestIntersection.bmp");
    ps1 = new Translate(0.5, 0.0, 0.0, new Sphere);
    ps2 = new Translate(-0.5, 0.0, 0.0, new Cylinder);
    pSolid = new Difference(ps1, ps2);
    st.Render(im, SolidTestFunction, 1.0);
    im.WriteBMP("SolidTestDifference.bmp");

    

    printf("Video:\n");
    pc = new StereographicCamera;
    im.NewImage(1920, 1080, 3);
    vi.NewVideo("Torus.mpeg", im);
    psTorus = new Torus(0.5);
    for (fAng = 0.0; fAng < 360.0; fAng += 1.0) {
        fRad = fAng*D_PI/180.0;
        rgps[0] = new Translate(-3.0, -3.0, 0.0, new Rotate(0.0, fRad, 0.0, new Cylinder));
        rgps[1] = new Translate(-3.0, 0.0, 0.0, new Rotate(0.0, fRad, 0.0, new Torus(0.5)));
        rgps[2] = new Translate(-3.0, 3.0, 0.0, new Rotate(0.0, fRad, 0.0, new Scale(0.8, 0.8, 0.8, new Octahedron)));

        rgps[3] = new Translate(0.0, -3.0, 0.0, new Rotate(0.0, fRad, 0.0, new Scale(0.5, 0.5, 0.5, new Tetrahedron)));
        rgps[4] = new Translate(0.0, 0.0, 0.0, new Rotate(0.0, fRad, 0.0, new Icosahedron));
        rgps[5] = new Translate(0.0, 3.0, 0.0, new Rotate(0.0, fRad, 0.0, new Dodecahedron));

        rgps[6] = new Translate(3.0, -3.0, 0.0, new Rotate(0.0, fRad, 0.0, new Scale(0.2, 1.0, 1.0, new Sphere)));
        rgps[7] = new Translate(3.0, 0.0, 0.0, new Rotate(0.0, fRad, 0.0, new Torus(0.2)));
        rgps[8] = new Translate(3.0, 3.0, 0.0, new Rotate(0.0, fRad, 0.0, new Difference(new Cube,
                                                                                      new Scale(0.5, 0.5, 1.1, new Cylinder))));
        ps1 = rgps[0];
        for (i = 1; i < 9; i++)
            ps1 = new Union(rgps[i], ps1);
        pSolid = ps1;
        //pSolid = new RotateSolid(fAng*D_PI/180.0, 0.0, 0.0, psTorus);
        st.Render(im, SolidTestFunction, 1.0);
        vi.WriteFrame(im);
    }
    vi.Close();
}

void
TestImageFunctions()
{
    Image im;
    TestSamplingTile st;

    im.NewImage(TESTIMAGESIZE, TESTIMAGESIZE, 1);
    st.Render(im, ZoneFunction, 1.0);
    im.WriteBMP("NewTestZone.bmp");
    st.Render(im, FactoryFunction, 1.0);
    im.WriteBMP("NewTestFactory.bmp");
    PerfectFactory(im);
    im.WriteBMP("TestPerfectFactory.bmp");
    PerfectZone(im);
    im.WriteBMP("TestPerfectZone.bmp");
}

extern double CubeRoot(double);

void
TestSolids()
{
    Sphere s;

    CubeRoot(70000.0);
    CubeRoot(0.00000000000001);
    CubeRoot(7.0);
}

void
TestSampling()
{
    double x, x2, y;

    printf("/nTestSampling:\n");
    for (x = -5.0; x <= 5.0; x += 0.25) {
        //printf("%6.1f %8f\n", x, Kaiser4Filter(x));
    }
    for (x = 0.0; x < 0.5; x += 1.0/64.0) {
        x2 = x*x;
        y = 1.0 + x2*(-1.0/6.0 + x2*(1.0/120.0 + x2*(-1.0/5040.0 + x2/362880.0)));
        printf("%10.6f %16.14f %16.14f\n", x, y, sin(x)/x);
    }
}

//
//
//  Make an image of the sampling pattern, the locality of sample order, and spectrum of pattern
//
//

#define NSIZE 512
#define SQ 0.70710678118654752440084436210485f
#define SK 0.38268343236508977172845998403040f
#define CK 0.92387953251128675612818318939680f
#define TP float(2.0*D_PI)
#define BASE 0.09f

static void reorder(float a[], float b[], int n);

static void
FFT(float a[], float b[], int n)
{
	int span, j, jj, k, kb, m, mm, mk;
	int k0, k2, k1, k3;
	float rad, c1, c2, c3, s1, s2, s3;
	float a0, a1, a2, a3, b0, b1, b2, b3;
	int c[32];

	for (m = 0, c[0] = 1; c[m] < n; m++)
		c[m + 1] = c[m] + c[m];
	mm = m & (~1);
	rad = TP / (float)n;
	mk = m - 5;
	kb = 0;
	if (mm != m)
		for (k0 = 0, k2 = c[mm]; k2 < n; k0++, k2++) {
			a0 = a[k2];
			b0 = b[k2];
			a[k2] = a[k0] - a0;
			a[k0] += a0;
			b[k2] = b[k0] - b0;
			b[k0] += b0;
		}
	c1 = 1.0;
	s1 = 0.0;
	jj = 0;
	k = mm - 2;
	j = 3;
	if (k >= 0)
		goto L4;
	else
		goto end;
L3:
	if (c[j] <= jj) {
		jj -= c[j];
		j--;
		if (c[j] <= jj) {
			jj -= c[j];
			j--;
			k += 2;
			goto L3;
		}
	}
	jj += c[j];
	j = 3;
L4:
	span = c[k];
	if (jj) {
		c2 = jj*span*rad;
		c1 = cosf(c2);
		s1 = sinf(c2);
	L5:
		c2 = c1*c1 - s1*s1;
		s2 = 2.0f*c1*s1;
		c3 = c2*c1 - s2*s1;
		s3 = c2*s1 + s2*c1;
	}
	for (k0 = kb + span - 1; k0 >= kb; --k0) {
		k1 = k0 + span;
		k2 = k1 + span;
		k3 = k2 + span;
		a0 = a[k0];	b0 = b[k0];
		if (s1 == 0.0) {
			a1 = a[k1];	b1 = b[k1];
			a2 = a[k2];	b2 = b[k2];
			a3 = a[k3];	b3 = b[k3];
		} else {
			a1 = a[k1]*c1 - b[k1]*s1;
			b1 = a[k1]*s1 + b[k1]*c1;
			a2 = a[k2]*c2 - b[k2]*s2;
			b2 = a[k2]*s2 + b[k2]*c2;
			a3 = a[k3]*c3 - b[k3]*s3;
			b3 = a[k3]*s3 + b[k3]*c3;
		}
		a[k0] = a0 + a2 + a1 + a3;	b[k0] = b0 + b2 + b1 + b3;
		a[k1] = a0 + a2 - a1 - a3;	b[k1] = b0 + b2 - b1 - b3;
		a[k2] = a0 - a2 - b1 + b3;	b[k2] = b0 - b2 + a1 - a3;
		a[k3] = a0 - a2 + b1 - b3;	b[k3] = b0 - b2 - a1 + a3;
	}
	if (k > 0) {
		k -= 2;
		goto L4;
	}
	kb = k3 + span;
	if (kb < n) {
		if (j == 0) {
			k = 2;
			j = mk;
			goto L3;
		}
		j--;
		c2 = c1;
		if (j == 1) {
			c1 = c1*CK + s1*SK;
			s1 = s1*CK - c2*SK;
		} else {
			c1 = (c1 - s1)*SQ;
			s1 = (c2 + s1)*SQ;
		}
		goto L5;
	}
end:
	reorder(a, b, n);
}

static void
reorder(float a[], float b[], int n)
{
	int c[32], lst[32];
	int m, i, j, k, jj;
	int kb, ku, lim, p;
	int kk, k2;
	float t;

	for (m = 0, c[0] = 1; c[m] < n; m++)
		c[m + 1] = c[m] + c[m];
	p = j = m - 1;
	i = kb = 0;
	--m;
	lim = (m + 2) >> 1;
	if (p <= 0)
		return;
	for (;;) {
		ku = k2 = c[j] + kb;
		jj = c[m - j];
		kk = kb + jj;
		do {
			k = kk + jj;
			do {
				t = a[kk];
				a[kk] = a[k2];
				a[k2] = t;
				t = b[kk];
				b[kk] = b[k2];
				b[k2] = t;
				kk++;
				k2++;
			} while (kk < k);
			kk += jj;
			k2 += jj;
		} while (kk < ku);
		if (j > lim) {
			--j;
			i++;
			lst[i] = j;
			continue;
		}
		kb = k2;
		if (i > 0) {
			j = lst[i];
			--i;
			continue;
		}
		if (kb < n) {
			j = p;
			continue;
		}
		break;
	}
	a++;
	b++;
	for (k = n - 2; k > 0; k -= 2) {
		t = a[k];
		a[k] = *a;
		*a++ = t;
		t = b[k];
		b[k] = *b;
		*b++ = t;
	}
}

static void
ML_FFT(float a[], float b[], int n)
{
    float fTmp;
    int i, nHalf;

    FFT(a, b, n);
    //
    //  Shift zero frequency (DC) to center
    //
    nHalf = n/2;
    for (i = 0; i < nHalf; i++) {
        fTmp = a[i];
        a[i] = a[i + nHalf];
        a[i + nHalf] = fTmp;
        fTmp = b[i];
        b[i] = b[i + nHalf];
        b[i + nHalf] = fTmp;
    }
}

Complex
cis(double f)
{
    return Complex(cos(f), sin(f));
}

void
SlowTransform(Complex x[20][127], Complex S[20][127])
{
    int tau, nu, p, G_rbo;
    Complex c, rgc[20][127];

    G_rbo = 0;
    for (tau = 0; tau < 127; tau++) {
        for (nu = 0; nu < 20; nu++) {
            c = 0;
            for (p = 0; p < 20; p++) {
                c = c + x[p][tau]*cis(-2.0*D_PI*(nu - 9 + G_rbo)*(tau + 127*p)/2540.0);
            }
            rgc[nu][tau] = c;
        }
    }
    for (tau = 0; tau < 127; tau++)
        for (nu = 0; nu < 20; nu++)
            S[nu][tau] = rgc[nu][tau];
}


int
ImageToSpectrum(Image &imgSpec, const Image &img, int bLog)
{
    float *pfReal, *pfImag, *pfReChan, *pfImChan, fReal, fImag, fMag, fPhase, fScale;
    int i, j, k, nImage, nWide, nHigh, nResult;

    nWide = img.m_nWidth;
    nHigh = img.m_nHeight;
    nImage = nWide*nHigh;
    pfReal = pfImag = 0;
    pfReal = new float[nHigh];
    pfImag = new float[nHigh];
    nResult = 0;
    if (pfReal == 0 || pfImag == 0 || 
        imgSpec.NewImage(nWide, nHigh, 2*img.m_nChannels, img.m_bWrapX, img.m_bWrapY) == 0)
        goto Fail;
    fScale = 1.0f/float(nImage);
    for (k = 0; k < img.m_nChannels; k++) {
        pfReChan = imgSpec.m_pfImage + k*nImage;
        pfImChan = pfReChan + nImage*img.m_nChannels;
        for (i = 0; i < nWide; i++) {
            for (j = 0; j < nHigh; j++) {
                pfReal[j] = img.m_pfImage[k*nImage + j*nWide + i];
                pfImag[j] = 0.0;
            }
            ML_FFT(pfReal, pfImag, nHigh);  // columns
            for (j = 0; j < nHigh; j++) {
                pfReChan[j*nWide + i] = pfReal[j];  // zero frequency at center
                pfImChan[j*nWide + i] = pfImag[j];
            }
        }
        for (j = 0; j < nHigh; j++) {       // rows
            ML_FFT(pfReChan + j*nWide, pfImChan + j*nWide, nWide);
        }
        for (i = 0; i < nWide; i++) {
            for (j = 0; j < nHigh; j++) {
                fReal = pfReChan[j*nWide + i];
                fImag = pfImChan[j*nWide + i];
                fPhase = atan2f(fReal, fImag);
                if (bLog)
                    fMag = BASE*float(log(fScale*sqrt(fReal*fReal + fImag*fImag) + 1.0e-35)) + 1.0f;
                else
                    fMag = sqrt(fReal*fReal + fImag*fImag);
                pfReChan[j*nWide + i] = fMag;
                pfImChan[j*nWide + i] = fPhase;
            }
        }
    }
    nResult = 1;
Fail:
    delete [] pfReal;
    delete [] pfImag;
    return nResult;
}

static void
DrawLocality(Image &im, double x1, double y1, double x2, double y2)
{
    if (x2 - x1 >  NSIZE/2) x2 -= NSIZE;
    if (x2 - x1 < -NSIZE/2) x2 += NSIZE;
    if (y2 - y1 >  NSIZE/2) y2 -= NSIZE;
    if (y2 - y1 < -NSIZE/2) y2 += NSIZE;
    im.DrawLineRGB(ML_White, x1, y1, x2, y2, 1.5);
}

void
PatternTestImage(char *szFileName, Vector2 rgv[], int nSamples)
{
    Image im, imPattern, imLocality, imSpectrum;
    int i, j;

    im.NewImage(3*NSIZE, NSIZE, 3);
    imPattern.NewImage(NSIZE, NSIZE, 3);
    imLocality.NewImage(NSIZE, NSIZE, 3, WRAP_AROUND, WRAP_AROUND);
    im.FillRGB(ML_Black);
    imPattern.FillRGB(ML_Black);
    imLocality.FillRGB(ML_Black);
    for (i = 0; i < nSamples; i++) {
        im.SplatRGB(ML_White, NSIZE * rgv[i].x, NSIZE * rgv[i].y);
        imPattern.SplatRGB(ML_White, NSIZE * rgv[i].x, NSIZE * rgv[i].y);
    }
    for (i = 0; i < nSamples - 1; i++) {
        DrawLocality(imLocality, NSIZE * rgv[i].x + NSIZE, NSIZE * rgv[i].y, 
                                 NSIZE * rgv[i+1].x + NSIZE, NSIZE * rgv[i+1].y);
    }
    ImageToSpectrum(imSpectrum, imPattern, 1);
    for (i = 0; i < NSIZE; i++) {
        for (j = 0; j < NSIZE; j++) {
            im.SetRGB(imSpectrum.Get(i, j), i + 2*NSIZE, j);
            im.SetRGB(imLocality.Get(i, j), i + 1*NSIZE, j);
        }
    }
    /*
    imLocality.FillRGB(ML_Black);
    for (i = 0; i < nSamples; i++) {
        rgv[i].x *= 512.0;          // 512 == image width, not NPTEST
        rgv[i].y *= 512.0;
    }
    for (i = 0; i < nSamples - 1; i++) {
        imLocality.DrawLineRGB(ML_White, rgv[i].x, rgv[i].y, rgv[i+1].x, rgv[i+1].y, 1.5);
    }
    imLocality.WriteBMP("Locality.bmp");
    */
    im.WriteBMP("PatternPaths.bmp");
    im.WriteBMP(szFileName);
}

#define NPSQRT 32
#define NPTEST (NPSQRT*NPSQRT)

void
TestPatterns()
{
    Image im;
    static Vector2 rgv[NPTEST];
    int i;

    BlueJitterPattern(rgv, NPSQRT);
    PatternTestImage("BlueJitter32Test.bmp", rgv, NPTEST);
    return;

    BlueNoisePattern(rgv, NPSQRT);
    OptimizeLocality(rgv, NPTEST);
    PatternTestImage("BlueTest.bmp", rgv, NPTEST);
    return;

    //BlueNoisePattern(rgv, NPSQRT);
    for (i = 0; i < NPTEST; i++) {
        rgv[i].x *= 512.0;          // 512 == image width, not NPTEST
        rgv[i].y *= 512.0;
    }
    im.NewImage(512, 512, 3);
    im.FillRGB(ML_Black);
    for (i = 0; i < NPTEST; i++) {
        im.SplatRGB(ML_White, rgv[i].x, rgv[i].y);
    }
    im.WriteBMP("PatternDots.bmp");

    im.FillRGB(ML_Black);
    for (i = 0; i < NPTEST - 1; i++) {
        im.DrawLineRGB(ML_White, rgv[i].x, rgv[i].y, rgv[i+1].x, rgv[i+1].y, 1.5);
    }
    im.WriteBMP("PatternPaths.bmp");

    //OptimizeLocality(rgv, NPTEST);
    im.FillRGB(ML_Black);
    for (i = 0; i < NPTEST - 1; i++) {
        im.DrawLineRGB(ML_White, rgv[i].x, rgv[i].y, rgv[i+1].x, rgv[i+1].y, 1.5);
    }
    im.WriteBMP("PatternOptimalPaths.bmp");
}
//
//  Build the RGB to Spectrum tables
//

static double
Smoothness(Spectrum &s)
{
    double f;
    int i;

    f = 0;
    for (i = 0; i < NSPEC-1; i++)
        f += (s.rgf[i] - s.rgf[i+1]) * (s.rgf[i] - s.rgf[i+1]);
    for (i = 0; i < NSPEC; i++)
        if (s.rgf[i] < 0.0)
            f += 100.0 * s.rgf[i]*s.rgf[i];
    return f;
}

static float mat[3][3] ={ 0.412410914897918f,   0.357584565877914f, 0.180453807115554f,
                          0.212649390101432f,   0.715169131755828f, 0.0721815153956413f,
                          0.0193317625671625f,  0.119194857776165f, 0.950390040874481f };

float   s_rgfWeight[NSPEC];
int     s_iDominant;

static inline float
NormalizedDiff(DisplayRGB rgb, DisplayRGB rgb2)                 // this compares closeness of chroma
{
    float X, Y, Z, fNorm, X2, Y2, Z2, fNorm2;

    X = rgb.red*mat[0][0] + rgb.grn*mat[0][1] + rgb.blu*mat[0][2];
    Y = rgb.red*mat[1][0] + rgb.grn*mat[1][1] + rgb.blu*mat[1][2];
    Z = rgb.red*mat[2][0] + rgb.grn*mat[2][1] + rgb.blu*mat[2][2];
    fNorm = fabs(X) + fabs(Y) + fabs(Z);
    X /= fNorm; Y /= fNorm; Z /= fNorm;
    X2 = rgb2.red*mat[0][0] + rgb2.grn*mat[0][1] + rgb2.blu*mat[0][2];
    Y2 = rgb2.red*mat[1][0] + rgb2.grn*mat[1][1] + rgb2.blu*mat[1][2];
    Z2 = rgb2.red*mat[2][0] + rgb2.grn*mat[2][1] + rgb2.blu*mat[2][2];
    fNorm2 = fabs(X2) + fabs(Y2) + fabs(Z2);
    X2 /= fNorm2; Y2 /= fNorm2; Z2 /= fNorm2;
    return (X - X2)*(X - X2) + (Y - Y2)*(Y - Y2) + (Z - Z2)*(Z - Z2);

    if (rgb.Norm())
        rgb = rgb / rgb.Norm();
    rgb2 = rgb2 / rgb2.Norm();
    return (rgb.red - rgb2.red)*(rgb.red - rgb2.red) 
         + (rgb.grn - rgb2.grn)*(rgb.grn - rgb2.grn) 
         + (rgb.blu - rgb2.blu)*(rgb.blu - rgb2.blu);
}

static int
DominantWavelength(DisplayRGB rgb)
{
    float d, mind;
    int i, mini;

    mini = 0;
    mind = NormalizedDiff(rgb, SpectrumToRGB[0]);
    for (i = 1; i < NSPEC; i++) {
        d = NormalizedDiff(rgb, SpectrumToRGB[i]);
        if (d < mind) {
            mini = i;
            mind = d;
        }
    }
    for (i = 0; i < NSPEC; i++)
        s_rgfWeight[i] = 1.0f;
    s_rgfWeight[mini] = 2.0f;               // weight dominant wavelength
    s_iDominant = mini;
    return mini;
}

static double
DifferenceRGB(DisplayRGB rgb, Spectrum &s)
{
    DisplayRGB rgb2;

    rgb2 = s.sRGB();
    return (rgb.red - rgb2.red)*(rgb.red - rgb2.red) 
         + (rgb.grn - rgb2.grn)*(rgb.grn - rgb2.grn) 
         + (rgb.blu - rgb2.blu)*(rgb.blu - rgb2.blu);
}

static void
Positive1(Spectrum &s)
{
    int i;

    for (i = 0; i < NSPEC; i++)
        if (s.rgf[i] < 0.0)
            s.rgf[i] = 0.0;
}

static void
Positive(Spectrum &s)
{
    float minf;
    int i;

    minf = s.rgf[0];
    for (i = 1; i < NSPEC; i++)
        if (s.rgf[i] < minf)
            minf = s.rgf[i];
    if (minf < 0.0)
        for (i = 0; i < NSPEC; i++)
            s.rgf[i] -= minf;
}

static Spectrum
RGBtoSpectrumSearch(DisplayRGB rgb)
{
    Spectrum s;
    double f, minf, fDelta, fD, fLambda;
    DisplayRGB rgb2;
    int i, ii, j, iDominant, rgi[NSPEC];

    fLambda = 1.0;                      // regularize with a smoothness, or spectra will be random
    iDominant = DominantWavelength(rgb);
    for (i = 0; i < NSPEC; i++) {
        rgi[i] = i;
        s.rgf[i] = float(0.5/SCALESPEC);
    }
    s.rgf[iDominant] = 1.0/SCALESPEC;
    minf = DifferenceRGB(rgb, s) + fLambda*Smoothness(s);
    for (fDelta = 0.5/SCALESPEC; fDelta >= 0.0009765625*0.0009765625*0.0009765625; fDelta *= sqrt(0.5)) {
        for (j = 0; j < 50; j++) {
            i = iDominant;                                              // first adjust the dominant wavelength
            s.rgf[i] += float(fDelta);
            f = DifferenceRGB(rgb, s) + fLambda*Smoothness(s);
            if (f < minf) {
                minf = f;
                continue;
            }
            s.rgf[i] -= float(2.0*fDelta);
            if (s.rgf[i] >= 0.0) {                                      // don't let any values go negative
                f = DifferenceRGB(rgb, s) + fLambda*Smoothness(s);
                if (f < minf) {
                    minf = f;
                    continue;
                }
            }
            s.rgf[i] += float(fDelta);

            RandomShuffle(rgi, NSPEC);
            fD = 1.0*fDelta;
            for (ii = 0; ii < NSPEC; ii++) {
                i = rgi[ii];
                s.rgf[i] += float(fD);
                f = DifferenceRGB(rgb, s) + fLambda*Smoothness(s);
                if (f < minf) {
                    minf = f;
                    continue;
                }
                s.rgf[i] -= float(2.0*fD);
                if (s.rgf[i] >= 0.0) {                                      // don't let any values go negative
                    f = DifferenceRGB(rgb, s) + fLambda*Smoothness(s);
                    if (f < minf) {
                        minf = f;
                        continue;
                    }
                }
                s.rgf[i] += float(fD);
            }
        }
        fLambda *= 0.8;
    }
    //Positive(s);
    return s;
}

static void
ReportSpectrum(DisplayRGB rgb, Spectrum &s, char *sz)
{
    DisplayRGB rgb2;

    // rgb.Print();
    rgb2 = s.sRGB();
    printf("//  ");
    rgb2.Print();
    printf("Spectrum g_s%s =", sz);
    s.Print();

}

static void
PlotSpectrum(char *szName, Spectrum &s)
{
    DisplayRGB rgb;
    Image im;
    int i, j, k, kLo, kHi, mini;
    char sz[128];

    rgb = s.sRGB();
    mini = DominantWavelength(rgb);
    im.NewImage(640, 480, 3);
    im.FillRGB(ML_Black); 
    sprintf_s(sz, 128, "(%0.3f, %0.3f, %0.3f)", rgb.red, rgb.grn, rgb.blu);
    im.DrawTextRGB(rgb, 320.0, 440.0, sz, 2.0, 0.0, JUSTIFIED_CENTER);
    for (i = 0; i < NSPEC; i++) {
        rgb = SpectrumToRGB[i] * 0.5;
        kLo = 0;
        kHi = int(floor(200.0*s.rgf[i] + 0.5));
        if (kHi < kLo) {
            kLo = kHi;
            kHi = 0;
        }
        for (j = 16; j < 48; j++) {
            for (k = kLo + 100; k <= kHi + 100; k++)
                    im.SetRGB(rgb, j + 64*i, 480 - k);
        }
    }
    i = mini;
        for (j = 16; j < 48; j++) {
            for (k = 0; k < 32; k++)
                    im.SetRGB(ML_White, j + 64*i, 480 - k);
        }
    im.WriteBMP(szName);
}

void
BuildSpectrumRGB()
{
    Spectrum sRed, sGreen, sBlue, sYellow, sCyan, sMagenta, sWhite, sOrange, s;
    DisplayRGB rgbRed, rgbGreen, rgbBlue, rgbYellow, rgbCyan, rgbMagenta, rgbWhite, rgbOrange, rgb1, rgb2;
    int i, j;

    //
    //  Test the tables
    //
    for (i = 0; i < 0; i++) {
        rgb1 = DisplayRGB(RandomFloat(), RandomFloat(), RandomFloat());
        s = RGBtoSpectrum(rgb1);
        rgb2 = s.sRGB();
        rgb1.Print();
        rgb2.Print();
        printf("\n");
    }
    
    //
    //  Compute the tables
    //
    sRed = RGBtoSpectrumSearch(ML_Red);
    sGreen = RGBtoSpectrumSearch(ML_Green);
    sBlue = RGBtoSpectrumSearch(ML_Blue);
    sYellow = RGBtoSpectrumSearch(ML_Red + ML_Green);
    sCyan = RGBtoSpectrumSearch(ML_Green + ML_Blue);
    sMagenta = RGBtoSpectrumSearch(ML_Red + ML_Blue);
    sWhite = RGBtoSpectrumSearch(ML_White);
    sOrange = RGBtoSpectrumSearch(ML_Red + ML_Green*0.5);

    rgbRed = sRed.sRGB();
    rgbGreen = sGreen.sRGB();
    rgbBlue = sBlue.sRGB();
    rgbYellow = sYellow.sRGB();
    rgbCyan = sCyan.sRGB();
    rgbMagenta = sMagenta.sRGB();
    rgbWhite = sWhite.sRGB();
    rgbOrange = sOrange.sRGB();
    /*
    ReportSpectrum(ML_Red, sRed);
    ReportSpectrum(ML_Green, sGreen);
    ReportSpectrum(ML_Blue, sBlue);
    ReportSpectrum(ML_Red+ML_Green, sYellow);
    ReportSpectrum(ML_Green+ML_Blue, sCyan);
    ReportSpectrum(ML_Red+ML_Blue, sMagenta);
    ReportSpectrum(ML_White, sWhite);
    */
    for (j = 0; j < 3; j++) {
        for (i = 0; i < 3; i++) {
            sRed = sRed - sGreen * rgbRed.grn/rgbGreen.grn;
            rgbRed = sRed.sRGB();
            sRed = sRed - sBlue * rgbRed.blu/rgbBlue.blu;
            rgbRed = sRed.sRGB();
            sRed = sRed / rgbRed.red;
            rgbRed = sRed.sRGB();
            //ReportSpectrum(ML_Red, sRed);

            sGreen = sGreen - sRed * rgbGreen.red/rgbRed.red;
            rgbGreen = sGreen.sRGB();
            sGreen = sGreen - sBlue * rgbGreen.blu/rgbBlue.blu;
            rgbGreen = sGreen.sRGB();
            sGreen = sGreen / rgbGreen.grn;
            rgbGreen = sGreen.sRGB();
            //ReportSpectrum(ML_Green, sGreen);

            sBlue = sBlue - sRed * rgbBlue.red/rgbRed.red;
            rgbBlue = sBlue.sRGB();
            sBlue = sBlue - sGreen * rgbBlue.grn/rgbGreen.grn;
            rgbBlue = sBlue.sRGB();
            sBlue = sBlue / rgbBlue.blu;
            rgbBlue = sBlue.sRGB();
            //ReportSpectrum(ML_Blue, sBlue);
        }
        sYellow = sYellow - sBlue * rgbYellow.blu;
        rgbYellow = sYellow.sRGB();
        sYellow = sYellow + sGreen * (rgbYellow.red - rgbYellow.grn);
        rgbYellow = sYellow.sRGB();
        sYellow = sYellow / rgbYellow.grn;
        rgbYellow = sYellow.sRGB();
        //ReportSpectrum(ML_Green + ML_Red, sYellow);

        sCyan = sCyan - sRed * rgbCyan.red;
        rgbCyan = sCyan.sRGB();
        sCyan = sCyan + sGreen * (rgbCyan.blu - rgbCyan.grn);
        rgbCyan = sCyan.sRGB();
        sCyan = sCyan / rgbCyan.grn;
        rgbCyan = sCyan.sRGB();
        //ReportSpectrum(ML_Green + ML_Blue, sCyan);

        sMagenta = sMagenta - sGreen * rgbMagenta.grn;
        rgbMagenta = sMagenta.sRGB();
        sMagenta = sMagenta + sRed * (rgbMagenta.blu - rgbMagenta.red);
        rgbMagenta = sMagenta.sRGB();
        sMagenta = sMagenta / rgbMagenta.red;
        rgbMagenta = sMagenta.sRGB();
        //ReportSpectrum(ML_Red + ML_Blue, sMagenta);

        sWhite = sWhite - sGreen * (rgbWhite.grn - rgbWhite.red);
        rgbWhite = sWhite.sRGB();
        sWhite = sWhite - sBlue * (rgbWhite.blu - rgbWhite.red);
        rgbWhite = sWhite.sRGB();
        sWhite = sWhite / rgbWhite.red;
        rgbWhite = sWhite.sRGB();
        
        Positive(sRed);
        Positive(sGreen);
        Positive(sBlue);
        Positive(sYellow);
        Positive(sCyan);
        Positive(sMagenta);
        Positive(sWhite);
        /*
        ReportSpectrum(ML_Red, sRed);
        ReportSpectrum(ML_Green, sGreen);
        ReportSpectrum(ML_Blue, sBlue);
        ReportSpectrum(ML_Red+ML_Green, sYellow);
        ReportSpectrum(ML_Green+ML_Blue, sCyan);
        ReportSpectrum(ML_Red+ML_Blue, sMagenta);
        ReportSpectrum(ML_White, sWhite);
        */
    }
    ReportSpectrum(ML_Red, sRed, "Red");
    ReportSpectrum(ML_Green, sGreen,"Green");
    ReportSpectrum(ML_Blue, sBlue, "Blue");
    ReportSpectrum(ML_Red+ML_Green, sYellow, "Yellow");
    ReportSpectrum(ML_Green+ML_Blue, sCyan, "Cyan");
    ReportSpectrum(ML_Red+ML_Blue, sMagenta, "Magenta");
    ReportSpectrum(ML_White, sWhite, "White");
    ReportSpectrum(rgbOrange, sOrange, "Orange");

    PlotSpectrum("Red.bmp", sRed);
    PlotSpectrum("Green.bmp", sGreen);
    PlotSpectrum("Blue.bmp", sBlue);
    PlotSpectrum("Yellow.bmp", sYellow);
    PlotSpectrum("Cyan.bmp", sCyan);
    PlotSpectrum("Magenta.bmp", sMagenta);
    PlotSpectrum("White.bmp", sWhite);
    PlotSpectrum("Orange.bmp", sOrange);
}

void
TestSpectrum()
{
    Spectrum s;

    printf("Equal White: ");
    g_sEqualWhite.sRGB().Print();
    printf("\n");

    s = RGBtoSpectrum(g_rgbCopper);
    ReportSpectrum(s.sRGB(), s, "Copper");
    PlotSpectrum("Copper.bmp", s);

    s = g_sCopper / 1.3;
    ReportSpectrum(s.sRGB(), s, "RealCopper");
    PlotSpectrum("RealCopper.bmp", s);

    s = BlackBody(6500.0);
    ReportSpectrum(s.sRGB(), s, "D65");
    PlotSpectrum("D6500.bmp", s);
    //(s * 2.63).Print();
    //BuildSpectrumRGB();
}
//
//  Test and compare cube-root routines
//
//  CubeRoot interative:
//  11571845.886 cube roots/sec
//  CubeRoot2 halley/newton:
//  26062212.464 cube roots/sec
//  Error1 = 0.00121602, Error2 = 2.67112e-010
//
extern double CubeRoot(double);

double
CubeRoot2(double f)         // conclusion: this is not as fast or accurate as CubeRoot
{
    double x;
    int i, nExp;
    
    frexp(f, &nExp);
    x = ldexp(1.0, nExp/3);
    for (i = 0; i < 5; i++) {
        x = x - (x*x*x - f)/(3.0*x*x);
        //printf("%d. x = %0.10f, x*x*x = %0.10f, ERR=%g\n", i, x, x*x*x, f - x*x*x);
    }
    return x;
}

double
CubeRoot_BlinnHack(double x)
{
    double t, t2;
    unsigned *pnT;

    if (x == 0.0)
        return 0.0;
    pnT = (unsigned *)&t;
    t = fabs(x);
    pnT[1] = pnT[1]/3 + 715094163;  // plus 2/3 of one as integer
    if (x < 0.0)
        t = -t;
    t2 = x/(t*t);
    t = t*(t + t2 + t2)/(t + t + t2);               // Halley-method step
    t = (t + t + x/(t*t))*0.33333333333333333333;   // Newton-method steps
    t = (t + t + x/(t*t))*0.33333333333333333333;
    return t;
}

double
CubeRoot_frexp(double x)
{
    double t, t2;
    int nExp;

    t = frexp(x, &nExp);
    t = ldexp(t, nExp/3);
    t2 = x/(t*t);
    t = t*(t + t2 + t2)/(t + t + t2);               // Halley-method step
    t = (t + t + x/(t*t))*0.33333333333333333333;   // Newton-method steps
    t = (t + t + x/(t*t))*0.33333333333333333333;
    return t;
}

#define CBRT2   1.259921049894873164767210607278228350570251
#define CBRT4   1.587401051968199474751705639272308260391493
#define CBRT2I  0.7937005259840997373758528196361541301957467
#define CBRT4I  0.6299605249474365823836053036391141752851257

double
CubeRoot_woboq(double x)
{
    int nExp, nRemainder, nSign;
    double z;

    if (x == 0)
        return (x);
    if (x > 0)
        nSign = 1;
    else {
        nSign = -1;
        x = -x;
    }
    z =x;
    x = frexp(x, &nExp);
    x = (((((1.3584464340920900529734e-1) * x
         - (6.3986917220457538402318e-1)) * x
          + (1.2875551670318751538055e0)) * x
          - (1.4897083391357284957891e0)) * x
          + (1.3304961236013647092521e0)) * x + (3.7568280825958912391243e-1);
    if (nExp >= 0) {
        nRemainder = nExp;
        nExp /= 3;
        nRemainder -= 3 * nExp;
        if (nRemainder == 1)
            x *= CBRT2;
        else if (nRemainder == 2)
            x *= CBRT4;
    } else {                                /* argument less than 1 */
        nExp = -nExp;
        nRemainder = nExp;
        nExp /= 3;
        nRemainder -= 3 * nExp;
        if (nRemainder == 1)
            x *= CBRT2I;
        else if (nRemainder == 2)
            x *= CBRT4I;
        nExp = -nExp;
    }
    x = ldexpl (x, nExp);
    x -= (x - (z / (x * x))) / 3.0;
    x -= (x - (z / (x * x))) / 3.0;
    x -= (x - (z / (x * x))) / 3.0;
    if (nSign < 0)
        x = -x;
    return x;
}

void
TestCubeRoot()                                          // woboq most accurate by far
{
    double x, x3, x33, xError1, xError2;
    int i;
    ML_TimingInfo ti;

    printf("CubeRoot interative:\n");
    xError1 = xError2 = 0.0;
    ML_StartTiming(ti);
    for (i = 0; i < 10000000; i++) {
        x = 2000000.0 * RandomDouble() - 1000000.0;
        x3 = x*x*x;
        x33 = CubeRoot_BlinnHack(x3);
        xError1 += fabs(x - x33);
    }
    ML_StopTiming(ti);
    ML_ReportEventsPerSecond(ti, i, "cube roots");
    printf("CubeRoot2 halley/newton:\n");
    ML_StartTiming(ti);
    for (i = 0; i < 10000000; i++) {
        x = 2000000.0 * RandomDouble() - 1000000.0;
        x3 = x*x*x;
        x33 = CubeRoot_woboq(x3);
        xError2 += fabs(x - x33);
    }
    ML_StopTiming(ti);
    ML_ReportEventsPerSecond(ti, i, "cube roots");
    printf("Error1 = %g, Error2 = %g\n", xError1, xError2);
}
//
//  Test stable quartic solver, 
//
extern int SolveQuartic(double rgfRoots[], double a, double b, double c, double d, double e);
extern int SolveCubic(double rgfRoots[], double B, double C, double D, double E);
extern int SolveQuadratic(double rgfRoots[], double C, double D, double E);

int
SolveQuartic2(double rgfRoots[], double A, double B, double C, double D, double E)
{
    double a, b, c, d, p, q, r, y, g, G, h, H, g1, g2, h1, h2, x;
    double rgfCubic[3];
    int nRoots, i, j;

    if (A == 0.0)
        return SolveCubic(rgfRoots, B, C, D, E);
    a = B/A;
    b = C/A;
    c = D/A;
    d = E/A;
    p = -2.0*b;
    q = b*b + a*c - 4.0*d;
    r = c*(c - a*b) + a*a*d;
    SolveCubic(rgfCubic, 1.0, p, q, r);
    y = rgfCubic[0];
    g1 = 0.5*a;
    g2 = 0.5*sqrt(a*a - 4.0*y);
    h1 = 0.5*(b - y);
    h2 = 0.5*(a*h1 - c)/g2;
    if (g1*g2 > 0.0) {                  // avoid cancellation error
        G = g1 + g2;
        g = y/G;
    } else {
        g = g1 - g2;
        G = y/g;
    }
    if (h1*h2 > 0.0) {
        H = h1 + h2;
        h = d/H;
    } else {
        h = h1 - h2;
        H = d/h;
    }
    nRoots  = SolveQuadratic(rgfRoots,        1.0, G, H);
    nRoots += SolveQuadratic(rgfRoots+nRoots, 1.0, g, h);
    for (j = 1; j < nRoots; j++) {
        x = rgfRoots[j];
        for (i = j - 1; i >= 0 && x < rgfRoots[i]; --i)
            rgfRoots[i + 1] = rgfRoots[i];
        rgfRoots[i + 1] = x;
    }
    return nRoots;
}

int
SolveQuarticHybrid(double rgfRoots[], double a, double b, double c, double d, double e)
{
    int k;

    if (b < 0.0) k = 1; else k = 0;
    if (c < 0.0) k += k+1; else k += k;
    if (d < 0.0) k += k+1; else k += k;
    if (e < 0.0) k += k+1; else k += k;
    switch (k) {
        case 0 : return SolveQuartic(rgfRoots, a, b, c, d, e);       // Ferrari
        case 1 : return SolveQuartic2(rgfRoots, a, b, c, d, e);      // Neumark
        case 2 : return SolveQuartic2(rgfRoots, a, b, c, d, e);
        case 3 : return SolveQuartic(rgfRoots, a, b, c, d, e);
        case 4 : return SolveQuartic(rgfRoots, a, b, c, d, e);
        case 5 : return SolveQuartic2(rgfRoots, a, b, c, d, e);
        case 6 : return SolveQuartic(rgfRoots, a, b, c, d, e);
        case 7 : return SolveQuartic(rgfRoots, a, b, c, d, e);
        case 8 : return SolveQuartic2(rgfRoots, a, b, c, d, e);
        case 9 : return SolveQuartic(rgfRoots, a, b, c, d, e);
        case 10: return SolveQuartic(rgfRoots, a, b, c, d, e);
        case 11: return SolveQuartic2(rgfRoots, a, b, c, d, e);
        case 12: return SolveQuartic(rgfRoots, a, b, c, d, e);
        case 13: return SolveQuartic(rgfRoots, a, b, c, d, e);
        case 14: return SolveQuartic(rgfRoots, a, b, c, d, e);
        case 15: return SolveQuartic(rgfRoots, a, b, c, d, e);
    }
    return 0;
}
//
//  quick test of RayTrace class and camera/model setup
//
static Radiance
TestShade1(Hit *ph)
{
    if (ph)
        return Radiance(g_sWhite);
    else
        return Radiance(g_sBlack);
}

//
//  Test material nodes
//
void
TestMaterial()
{
    Solid *ps;
    Ray ray;
    Hit *ph;

    ps = new Union(new Material(0.9, g_sYellow, MAT_PLASTIC, new Rotate(0.5*D_PI, 0.0, 0.0, new Torus(0.5))),
                   new Material(0.5, g_sBlue, MAT_PLASTIC, new Sphere));
    ray = Ray(Point3(0.0, 0.0, -5.0), Vector3(0.0, 0.0, 1.0));
    ph = ps->Intersect(ray);
    for (; ph; ph = ph->phNext)
        printf("t = %f, density = %f\n", ph->t, ph->pm->fDensity);
    printf("\n");

    ps = new Intersection(new Material(0.9, g_sYellow, MAT_PLASTIC, new Rotate(0.5*D_PI, 0.0, 0.0, new Torus(0.5))),
                          new Material(0.5, g_sBlue, MAT_PLASTIC, new Sphere));
    ray = Ray(Point3(0.0, 0.0, -5.0), Vector3(0.0, 0.0, 1.0));
    ph = ps->Intersect(ray);
    for (; ph; ph = ph->phNext)
        printf("t = %f, density = %f\n", ph->t, ph->pm->fDensity);
    printf("\n");

    ps = new Difference(new Material(0.9, g_sYellow, MAT_PLASTIC, new Rotate(0.5*D_PI, 0.0, 0.0, new Torus(0.5))),
                        new Material(0.5, g_sBlue, MAT_PLASTIC, new Sphere));
    ray = Ray(Point3(0.0, 0.0, -5.0), Vector3(0.0, 0.0, 1.0));
    ph = ps->Intersect(ray);
    for (; ph; ph = ph->phNext)
        printf("t = %f, density = %f\n", ph->t, ph->pm->fDensity);
    printf("\n");
}
//
//  Test Fresnel functions
//
static double
Fresnel2(double fNE, int nInsideDenser)
{
    double x, fCritical, fFresnel;

    return 1.0 - fNE*fNE*fNE*fNE*fNE;

    x = 1.0 + fNE;
    if (x < 0.0)
        x = 0.0;
    if (nInsideDenser) {
        fCritical = 0.7453559925;           // sqrt(1 - 1/n**2)
        x = x / fCritical;
        if (x > 1.0)
            return 1.0;                     // total internal reflection
    }
    x = x*x*x*x;
    if (x > 1.0)
        x = 1.0;
    fFresnel = 0.04;                        // ((n - 1)/(n + 1))**2 
    x = fFresnel + (1.0 - fFresnel)*x;
    return x;
}

#define ETA     0.617
#define KAPPA   2.63

static double
FresnelMetal(double fNE)
{
    double r1, r2, fEK;

    fEK = (ETA*ETA + KAPPA*KAPPA)*fNE*fNE;
    r1 = (fEK - 2.0*ETA*fNE + 1.0) / (fEK + 2.0*ETA*fNE + 1.0);
    fEK = ETA*ETA + KAPPA*KAPPA;
    r2 = (fEK - 2.0*ETA*fNE + fNE*fNE) / (fEK + 2.0*ETA*fNE + fNE*fNE);
    return 0.5*(r1 + r2);
}

static double
FresnelSchlick(double fNE)
{
    double R0, f1;

    R0 = 0.04;
    f1 = 1.0 - fNE;
    f1 = f1*f1*f1*f1*f1;
    return R0 + (1.0 - R0)*f1;
}

static double
Fresnel2(double fNE, double fN1, double fN2)
{
    double R0, fNEcritical, fNT, fEta, Rs, Rp;

    R0 = (fN1 - fN2)/(fN1 + fN2) * (fN1 - fN2)/(fN1 + fN2);
    fEta = fN1/fN2;
    fNEcritical = 1.0 - 1.0/(fEta*fEta);
    if (fNEcritical < 0.0)
        fNEcritical = 0.0;
    else
        fNEcritical = sqrt(fNEcritical);
    if (fNE < fNEcritical)
        return 1.0;                                 // total internal reflection
    fNT = sqrt(1.0 - fEta*fEta*(1.0 - fNE*fNE));
    Rs = (fN1*fNE - fN2*fNT)/(fN1*fNE + fN2*fNT);
    Rp = (fN1*fNT - fN2*fNE)/(fN1*fNT + fN2*fNE);
    return 0.5*(Rs*Rs + Rp*Rp);
}

static double
FresnelReflectance(double VN, double index1, double index2)
{
    double tau1, n;
    double tau;

    if(VN > 0.0)
	    VN = -VN;
    n = index2 / index1;		        // relative refractive index
    tau1 = (n - 1.0)/(n + 1.0);
    tau1 = tau1*tau1;
    if (n < 1.0)
	    tau = (1.0 + VN) / (1.0 - sqrt(1.0 - n*n));
    else
	    tau = (1.0 + VN);
    tau = tau*tau;
    tau = tau*tau;
    if (tau > 1.0)
	    tau = 1.0;		                // total internal reflection 
    return tau1 + (1.0 - tau1)*tau;
}

static double
FresnelReflectance2(double fNE, double fN1, double fN2)
{
    double fEta, fNEcritical, fNT, Rs, Rp, fFresnel;

    fEta = fN1/fN2;  
    fNEcritical = 1.0 - 1.0/(fEta*fEta);
    if (fNEcritical < 0.0)
        fNEcritical = 0.0;
    else
        fNEcritical = sqrt(fNEcritical);
    if (fNE < fNEcritical) {
        fFresnel = 1.0;                                 // total internal reflection
    } else {
        fNT = sqrt(1.0 - fEta*fEta*(1.0 - fNE*fNE));
        Rs = (fN1*fNE - fN2*fNT)/(fN1*fNE + fN2*fNT);
        Rp = (fN1*fNT - fN2*fNE)/(fN1*fNT + fN2*fNE);
        fFresnel = 0.5*(Rs*Rs + Rp*Rp);
        if (fN1 == fN2)
            fFresnel = 0.0;
    }
    return fFresnel;
}

static inline void
FresnelLaw_Real(double fNE, double fN1, double fN2, double &fFresnel, double &fNT)
{
    double fEta, Rs, Rp, fNEcritical;

    fEta = fN1/fN2;  
    fNEcritical = 1.0 - 1.0/(fEta*fEta);
    if (fNEcritical < 0.0)
        fNEcritical = 0.0;
    else
        fNEcritical = sqrt(fNEcritical);
    if (fNE < fNEcritical) {
        fFresnel = 1.0;                                     // total internal reflection
    } else {
        fNT = sqrt(1.0 - fEta*fEta*(1.0 - fNE*fNE));        // this may be used later to calculate refraction ray.vDir
        if (fN1 == fN2) {
            fFresnel = 0.0;
        } else {
            Rs = (fN1*fNE - fN2*fNT)/(fN1*fNE + fN2*fNT);
            Rp = (fN1*fNT - fN2*fNE)/(fN1*fNT + fN2*fNE);
            fFresnel = 0.5*(Rs*Rs + Rp*Rp);
        }
    }
}
//
//  Complex refractive indices for metals
//
static void
FresnelLaw_Complex(double fNE, Complex fN1, Complex fN2, double &fFresnel, double &fNT)
{
    Complex fEta, Rs, Rp;
    double fNEcritical;

    fEta = fN1/fN2;  
    fNEcritical = 1.0 - 1.0/fEta.AbsSqr();
    if (fNEcritical < 0.0)
        fNEcritical = 0.0;
    else
        fNEcritical = sqrt(fNEcritical);
    if (fNE < fNEcritical) {
        fFresnel = 1.0;                                     // total internal reflection
    } else {
        fNT = sqrt(1.0 - fEta.AbsSqr()*(1.0 - fNE*fNE));        // this may be used later to calculate refraction ray.vDir
        Rs = (fN1*fNE - fN2*fNT)/(fN1*fNE + fN2*fNT);
        Rp = (fN1*fNT - fN2*fNE)/(fN1*fNT + fN2*fNE);
        fFresnel = 0.5*(Rs.AbsSqr() + Rp.AbsSqr());
        if (fN1 == fN2)
            fFresnel = 0.0;
    }
}

static Spectrum
Fresnel_Spectral(double fNE, double fN1)
{
    Spectrum s;
    double fFresnel, fNT;
    int i;

    for (i = 0; i < NSPEC; i++) {
        FresnelLaw_Complex(fNE, Complex(fN1, 0.0), Complex(g_sCopperRefraction.rgf[i], g_sCopperExtinction.rgf[i]), fFresnel, fNT);
        s.rgf[i] = float(fFresnel);
    }
    return s;
}

static double
NormalRefl(Complex n2)
{
    double f1, f2;

    f1 = (n2.re - 1.0)*(n2.re - 1.0) + n2.im*n2.im;
    f2 = (n2.re + 1.0)*(n2.re + 1.0) + n2.im*n2.im;
    return f1/f2;
}

void
TestFresnel()
{
    double fAng, fCos, fFresnel, fNT;
    Complex n1, n2;
    Spectrum s;

    for (fAng = 0.0; fAng <= 180.0; fAng += 10.0) {
        fCos = cos(D_PI*fAng/180.0);
        s = Fresnel_Spectral(fCos, 1.0);
        if (fAng == 0.0) s.Print();
        printf("%f  ", fAng); s.sRGB().Print();
    }
    s.Print();

    n1 = Complex(1.0, 0.0);
    n2 = REFRACT_COPPER;
    printf("Normal incidence R = %f\n", NormalRefl(n2));
    FresnelLaw_Complex(1.0, n1, n2, fFresnel, fNT);
    printf("Complex R(1) = %f\n", fFresnel);
    FresnelLaw_Complex(1.0, n2, n1, fFresnel, fNT);
    printf("Complex R(1) = %f\n", fFresnel);
    FresnelLaw_Complex(0.0, n1, n2, fFresnel, fNT);
    printf("Complex R(0) = %f\n", fFresnel);
    FresnelLaw_Complex(0.0, n2, n1, fFresnel, fNT);
    printf("Complex R(0) = %f\n", fFresnel);
    return;

    for (fAng = 0.0; fAng <= 180.0; fAng += 10.0) {
        fCos = cos(D_PI*fAng/180.0);
        printf("%6.2f %8.3f %8.3f %8.3f\n", fAng, FresnelReflectance(fCos, 1.0, 1.5), FresnelReflectance(fCos, 1.5, 1.0), FresnelReflectance(fCos, 1.0, 1.0));
    }
    printf("\n");
    for (fAng = 0.0; fAng <= 180.0; fAng += 10.0) {
        fCos = cos(D_PI*fAng/180.0);
        printf("%6.2f %8.3f %8.3f %8.3f\n", fAng, FresnelReflectance2(fCos, 1.0, 1.5), FresnelReflectance2(fCos, 1.5, 1.0), FresnelReflectance2(fCos, 1.0, 1.0));
    }
    return;

    Fresnel2(0.0, 1.0, 1.5);
    Fresnel2(0.0, 1.5, 1.0);
    Fresnel2(0.0, 1.0, 1.8);
    Fresnel2(0.0, 1.0, 1.33);
    Fresnel2(0.0, 1.33, 1.0);
    return;
    for (fAng = 0.0; fAng <= 90.0; fAng += 10.0) {
        fCos = cos(D_PI*fAng/180.0);
        printf("%6.2f %8.3f %8.3f %8.3f\n", fAng, Fresnel2(fCos, 1), FresnelMetal(fCos), FresnelSchlick(fCos));
    }
}
//
//  Test reflection and refraction rays
//
static void
PrintHitList(char *sz, Hit *ph)
{
    printf("%s { ", sz);
    for (; ph; ph = ph->phNext)
        printf("t=%0.1f d=%0.1f, ", ph->t, ph->pm->fDensity);
    printf("}\n");
}

void
TestRefraction()
{
    RayTracer rt;
    Hit *ph;
    Radiance R;
    Ray ray;
    double fAng;
    Spectrum s6500;

    fAng = 20.0 * D_PI/180.0;
    s6500 = BlackBody(6500.0);
    rt.plAmbient = new AmbientLight(s6500 * 0.02);
    rt.plPointLights = new PointLight(Point3(10.0, -10.0, 10.0), s6500 * 10000.0);
    rt.psModel = new Surface(SURF_POLISHED, new Material(0.5, g_sWhite, MAT_TRANSPARENT, new Rotate(0.0, fAng, 0.0, new Cube)));
    ray = Ray(Point3(0.0, 0.0, 0.0), Vector3(0.0, 0.0, -1.0));
    ph = rt.psModel->Intersect(ray);
    PrintHitList("cast", ph);
    ray.Print();
    //R = rt.Shade(ph, ray);
    R.Print();
}



//
//  extents
//
void
TestExtents()
{
    Solid *ps, *ps2;
    Ray ray;
    Hit *ph;
    int n;

    ps = new Icosahedron;
    ps->GetExtent().Print();
    ps = new Torus(0.5);
    ps->GetExtent().Print();
    ps = new Translate(0.0, 3.0, 0.0, new Scale(1.0, 5.0, 1.0, new Sphere));
    ps->GetExtent().Print();
    ps = new Sphere;
    ps->GetExtent().Print();
    ps = new Sphere;
    ps2 = new Translate(0.0, 1.0, 0.0, new Sphere);
    ps = new Intersection(ps, ps2);
    ps->GetExtent().Print();
    ps = new Sphere;
    ps2 = new Translate(0.0, 3.0, 0.0, new Sphere);
    ps = new Union(ps, ps2);
    ps->GetExtent().Print();
    ps = new BoundingBox(ps);
    ps->GetExtent().Print();
    printf("\n");

    ray = Ray(Point3(0.0, 5.0, 0.0), Vector3(0.0, -1.0, 0.0));
    ps = new BoundingBox(new Sphere);
    //ps = new Rotate(0.0, 0.5*D_PI, 0.0, ps);

    n = ps->GetExtent().TestIntersection(ray);
    printf("intersect %d\n", n);
    ph = ps->Intersect(ray);
    PrintHitList("TestExtents", ph);
    return;

    printf("\nray test:\n");
    ray = Ray(Point3(0.0, 0.0, 5.0), Vector3(0.0, 0.0, -1.0));
    ph = ps->Intersect(ray);
    PrintHitList("TestExtents", ph);
    //printf("Test Extent: %d\n", ps->GetExtent().TestIntersection(ray));
    printf("ray miss:\n");
    ray.vDir = Vector3(0.0, 1.0, 0.0);
    ph = ps->Intersect(ray);
    PrintHitList("TestExtents", ph);
    //printf("Test Extent: %d\n", ps->GetExtent().TestIntersection(ray));


}
//
//  Model optimization
//
void
TestOptimize()
{
    Solid *ps, *psUnion;

    ps = new Sphere();
    ps = new Union(new Sphere, ps);
    ps = new Translate(1.0, 0.0, 0.0, ps);
    ps = new Rotate(0.0, D_PI/2.0, 0.0, ps);
    ps->Optimize();

    printf("\n");
    psUnion = new Translate(0.0, 3.0, -10.0, new Scale(12.0, 1.0, 12.0, new Cube));
    ps = new Material(0.5, REFRACT_PYREX, g_sBlack, MAT_TRANSPARENT, 0, new Surface(SURF_POLISHED, 0, new Cube));
    psUnion = new Union(new Translate(2.0, 0.0, 0.0, ps), psUnion);
    ps = new Material(0.5, REFRACT_WATER, g_sWater, MAT_TRANSPARENT, ColoredGlass, new Surface(SURF_POLISHED, 0, new Cube));
    psUnion = new Union(new Translate(-2.0, 0.0, 0.0, ps), psUnion);
    ps->Optimize();
    printf("\n");
    psUnion->Optimize();
}
//
//  Test textures
//
void
TestBandNoise()
{
    double fB, fF, fAveB, fAveF;
    float rgf[3];
    int i;

    fAveB = fAveF = 0.0;
    for (i = 0; i < 50; i++) {
        rgf[0] = float(10*RandomDouble());
        rgf[1] = float(10*RandomDouble());
        rgf[2] = float(10*RandomDouble());
        fF = FractalNoise(rgf, 3, 5, 1.0f, 1.9f, 0);
        fB = BandNoise(rgf, 3);
        printf("B = %f,  F = %f\n", fB, fF);
        fAveB += fB;
        fAveF += fF;
    }
    printf("Averages %f %f\n", fAveB/i, fAveF/i);
}
//
//  Test SampleGradient
//
static double gs_fI0Alpha = 1.0/11.3019219521363;
static double gs_fAlphaSquared = 4.0*4.0;

static float
Kaiser4Filter(float x)
{
    float x1, x2, xTerm, fN, fSinc, fKaiser;

    x1 = x * 1.0f/KAISER4_HALFWIDTH;
    x2 = x1*x1;
    if (x2 > 1.0)
        return 0.0;
    x2 = 0.25f * float(gs_fAlphaSquared) * (1.0f - x2);
    fKaiser = 1.0f + x2;
    fN = 2.0f;
    xTerm = x2;
    while (xTerm > 1.0e-4) {
        xTerm *= x2/(fN*fN);
        fN += 1.0;
        fKaiser += xTerm;
        xTerm *= x2/(fN*fN);
        fN += 1.0;
        fKaiser += xTerm;
    }
    x1 = float(D_PI*x);
    x2 = x1*x1;
    if (x2 < 0.0001f)
        fSinc = 1.0f + x2*(-1.0f/6.0f + x2*1.0f/120.0f);
    else
        fSinc = sinf(x1)/x1;
    return fSinc * fKaiser * float(gs_fI0Alpha); // * 1.0033f;    // correct normalization
}

static float
Sinc(float x)
{
    x *= float(D_PI);
    if (x < 0.01f && x > -0.01f) {
        x = x*x;
        return float(1.0 + x*(-1.0/6.0 + x*(1.0/120.0 + x*(-1.0/5040.0 + x/362880.0))));
    } else
        return float(sin(x)/x);
}

static float
DSinc(float x)
{
    if (x < 0.001 && x > -0.001) {
        x *= float(D_PI);
        return float(D_PI * x*(-2.0/6.0 + x*x*(4.0/120.0 + x*x*(-6.0/5040.0))));
    }
    return float(cos(D_PI*x) - Sinc(x))/x;
}

static float
DerivativeFilter(float x)
{
    float x1, x2, xTerm, fN, fDSinc, fKaiser;

    x1 = x * 1.0f/KAISER4_HALFWIDTH;
    x2 = x1*x1;
    if (x2 > 1.0)
        return 0.0;
    x2 = 0.25f * float(gs_fAlphaSquared) * (1.0f - x2);
    fKaiser = 1.0f + x2;
    fN = 2.0f;
    xTerm = x2;
    while (xTerm > 1.0e-4) {
        xTerm *= x2/(fN*fN);
        fN += 1.0;
        fKaiser += xTerm;
        xTerm *= x2/(fN*fN);
        fN += 1.0;
        fKaiser += xTerm;
    }
    fDSinc = DSinc(x);
    return fDSinc * fKaiser * float(gs_fI0Alpha); 
}

void
TestSampleGradient()
{
    Image im, imX, imY, imF, imDiff;
    double x, y, xmin, xmax, fDiff, xmin2, xmax2, fDeriv;
    float rgf[2], f;
    int i,  j;

    xmin  = +100000.0;
    xmax  = -100000.0;
    xmin2 = +100000.0;
    xmax2 = -100000.0;
    im.ReadBMP("RayTraceScene_Retort.bmp");
    imX.NewImage(1920, 1080);
    imY.NewImage(1920, 1080);
    imF.NewImage(1920, 1080);
    imDiff.NewImage(1920, 1080);
    for (x = 0.5; x < 1920.5; x += 1.0) {
        for (y = 0.5; y < 1080.5; y += 1.0) {
            i = int(x);
            j = int(y);
            im.SampleGradient(rgf, x, y);
            imX.Set(0.5f + 0.5f*rgf[0], int(x), int(y));
            imY.Set(0.5f + 0.5f*rgf[1], int(x), int(y));
            f = im.Sample(x, y);
            imF.Set(f, int(x), int(y));
            fDiff = (im.Sample(x+0.01, y) - im.Sample(x-0.01, y)) / 0.02;
            fDeriv = rgf[0];
            if (fDiff < xmin)
                xmin = fDiff;
            if (fDiff > xmax)
                xmax = fDiff;
            if (fDeriv < xmin2)
                xmin2 = fDeriv;
            if (fDeriv > xmax2)
                xmax2 = fDeriv;
            imDiff.Set(fDiff, int(x), int(y));
            if (i % 100 == 0 && j % 100 == 0)
                printf("diff %f   deriv %f\n", fDiff, fDeriv);
        }
    }
    imX.WriteBMP("gradx.bmp");
    imY.WriteBMP("grady.bmp");
    imF.WriteBMP("value.bmp");
    imDiff.WriteBMP("diffx.bmp");
    printf("diff: min max %f %f\n", xmin, xmax);
    printf("derv: min max %f %f\n", xmin2, xmax2);
return;
/*
    for (x = -5.0; x <= 5.0; x += 0.25) {
        printf("x = %f ds = %f deltas = %f\n", x, DerivativeFilter(x), (Kaiser4Filter(x+0.001f) - Kaiser4Filter(x-0.001f))/0.002);
    }
    printf("sinc\n");
    for (x = -5.0; x <= 5.0; x += 0.25) {
        printf("x = %f ds = %f deltas = %f\n", x, DSinc(x), (Sinc(x+0.001) - Sinc(x-0.001))/0.002);
    }
*/
}

//
//  Test RGB Spectrum
//
void
TestRGBSpectrum()
{
    DisplayRGB rgb, rgb2;
    Spectrum s;
    int i;

    printf("D6500 ");
    s = BlackBody(6500.0);
    s.sRGB().Print();
    s.Print();
    printf("\n");
    for (i = 0; i < 10; i++) {
        rgb = DisplayRGB(RandomDouble(), RandomDouble(), RandomDouble());
        s = RGBtoSpectrum(rgb);
        rgb2 = s.sRGB();
        rgb.Print();
        rgb2.Print();
        printf("\n");
    }
}
//
//  Test CLipPlanes to make extent
//
void
TestClipPlanes()
{
    Extent e;

    e = Extent(3.0);
    e = ClipPlanes(e, rgvCube, 8, 1.0);
    e.Print();
}

//
//  Generate spectra from old format
//
#define ML_SPECTRUM_START    360
#define ML_SPECTRUM_END      830
#define ML_SPECTRUM_STEP       5
#define ML_SPECTRUM_NDATA    ((ML_SPECTRUM_END - ML_SPECTRUM_START)/ML_SPECTRUM_STEP + 1)

struct ML_SpectralDensity {
    double rgf[ML_SPECTRUM_NDATA];
};

static Spectrum
MLtoSpectrum(const ML_SpectralDensity &sd)
{
    Spectrum s;
    int i, isd;

    for (i = 0; i < NSPEC; i++) {
        isd = i*NDELTASPEC/ML_SPECTRUM_STEP + (380 - ML_SPECTRUM_START)/ML_SPECTRUM_STEP;
        s.rgf[i] = float(sd.rgf[isd]);
    }
    return s;
}

static void
PrintSpectrum(char *szName, const Spectrum &s)
{
    int i;

    printf("Spectrum %s = {\n", szName);
    for (i = 0; i < NSPEC; i++) {
        printf("    %0.8ff,    // %d\n", s.rgf[i], 380 + NDELTASPEC*i);
    }
    printf("};\n");
}

ML_SpectralDensity ML_Venera14_000 = {
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 4.2754406929016f,  4.7265758514404f,  5.2252173423767f,  5.7762556076050f, 
 6.3857541084290f,  7.0595622062683f,  7.8039236068726f,  8.4556789398193f, 
 9.0632972717285f, 10.0854539871216f, 11.1388387680054f, 12.3256616592407f, 
13.6272487640381f, 15.0372591018677f, 16.5918292999268f, 18.2113342285156f, 
20.0180130004883f, 21.8910293579102f, 23.9336261749268f, 26.1448841094971f, 
28.4565067291260f, 30.5991649627686f, 32.6827468872070f, 34.7127799987793f, 
37.1206626892090f, 39.5640258789063f, 42.0539283752441f, 44.7331848144531f, 
47.4939231872559f, 50.6179580688477f, 53.5910148620605f, 56.8915634155273f, 
59.7307281494141f, 62.6434860229492f, 65.4312286376953f, 67.3309936523438f, 
68.5714645385742f, 69.5750656127930f, 70.0902023315430f, 70.2403488159180f, 
70.0932998657227f, 69.5118179321289f, 68.9794006347656f, 68.1121826171875f, 
67.2705001831055f, 66.1285629272461f, 64.9302673339844f, 63.4994850158691f, 
61.1402702331543f, 57.7050170898438f, 54.4045181274414f, 51.5560798645020f, 
48.9190101623535f, 46.5636215209961f, 44.8165092468262f, 43.9801406860352f, 
44.0369110107422f, 44.7193489074707f, 46.0993003845215f, 48.0129737854004f, 
49.7849617004395f, 51.5939331054688f, 53.1084022521973f, 53.6891021728516f, 
53.2677078247070f, 52.0145111083984f, 50.0021781921387f, 48.2184028625488f, 
46.8142776489258f, 46.1816444396973f, 46.1850547790527f, 46.6093521118164f, 
47.4579315185547f, 48.6802635192871f, 50.4559822082520f, 52.2673301696777f, 
54.1361694335938f, 55.9414825439453f, 57.9172592163086f, 
};

ML_SpectralDensity ML_Venera14_ground = {
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.7146394252777f, 
 1.0287981033325f,  1.1131254434586f,  1.2159054279327f,  1.3355289697647f, 
 1.4641534090042f,  1.5981974601746f,  1.7191938161850f,  1.8665820360184f, 
 2.0099160671234f,  2.1566274166107f,  2.3255722522736f,  2.4822158813477f, 
 2.6429715156555f,  2.8165218830109f,  3.0065689086914f,  3.1857826709747f, 
 3.3597571849823f,  3.5166771411896f,  3.6863250732422f,  3.8372313976288f, 
 3.9635088443756f,  4.0609378814697f,  4.1328783035278f,  4.1797447204590f, 
 4.1832523345947f,  4.1606764793396f,  4.1110396385193f,  4.0160708427429f, 
 3.9138987064362f,  3.7656757831573f,  3.6214213371277f,  3.4614870548248f, 
 3.3196225166321f,  3.2153573036194f,  3.1296482086182f,  3.0887014865875f, 
 3.0718123912811f,  3.0938925743103f,  3.1561410427094f,  3.2719984054565f, 
 3.4615058898926f,  3.6341335773468f,  3.7684698104858f,  3.8399972915649f, 
 3.8398959636688f,  3.7584984302521f,  3.6036117076874f,  3.4640460014343f, 
 3.3575990200043f,  3.3307313919067f,  3.3380537033081f,  3.3587317466736f, 
 3.4293169975281f,  3.5611724853516f,  3.7251620292664f,  3.9415683746338f, 
 4.2622385025024f,  4.5615463256836f,  4.8975844383240f, 
};

ML_SpectralDensity ML_Venera11_000 = {
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f, 10.1761684417725f, 11.3461761474609f, 
13.0054445266724f, 14.4933013916016f, 15.9710569381714f, 17.5292530059814f, 
19.3625621795654f, 21.5083560943604f, 24.0771884918213f, 26.7785320281982f, 
29.7019824981689f, 32.7629661560059f, 35.8059806823730f, 39.3028259277344f, 
43.5647354125977f, 48.7843132019043f, 54.7558479309082f, 59.4658927917480f, 
62.5749359130859f, 64.3997650146484f, 66.2776336669922f, 68.0810852050781f, 
69.7926483154297f, 71.6311111450195f, 73.4144668579102f, 75.1607971191406f, 
77.0463790893555f, 79.1236572265625f, 80.9176788330078f, 82.7832565307617f, 
84.7465515136719f, 86.0251998901367f, 86.5220489501953f, 86.5485076904297f, 
86.4257812500000f, 86.1011505126953f, 85.6617584228516f, 85.0110473632813f, 
84.0065307617188f, 82.8132553100586f, 81.6675186157227f, 80.3923263549805f, 
79.0146560668945f, 76.6344375610352f, 73.2253570556641f, 69.3763504028320f, 
65.0117187500000f, 60.3875160217285f, 56.0095977783203f, 52.8779869079590f, 
53.4052925109863f, 55.2187995910645f, 58.9158897399902f, 63.9598846435547f, 
67.3519287109375f, 69.3126220703125f, 69.5337295532227f, 68.5373382568359f, 
67.1025772094727f, 63.9992294311523f, 59.3493041992188f, 54.8607368469238f, 
50.8408317565918f, 47.7862358093262f, 46.1126403808594f, 45.4479217529297f, 
45.4176445007324f, 45.7428092956543f, 46.4800949096680f, 47.2450218200684f, 
47.8922386169434f, 48.6508483886719f, 49.6700134277344f, 51.2002143859863f, 
53.1546325683594f, 54.5670547485352f, 55.0512580871582f, 
};

ML_SpectralDensity ML_Venera11_072 = {
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f, 11.3268709182739f, 13.2905998229980f, 
15.5952358245850f, 18.2988185882568f, 21.4717178344727f, 25.1945667266846f, 
29.5625247955322f, 34.6887664794922f, 40.5789413452148f, 44.8406181335449f, 
52.6450805664063f, 61.1192893981934f, 70.2478103637695f, 80.1135253906250f, 
90.3623962402344f, 101.2378463745117f, 113.4931716918945f, 126.8272094726563f, 
138.0222930908203f, 148.1359100341797f, 158.0055999755859f, 167.2462921142578f, 
175.0063781738281f, 181.3072814941406f, 187.2035675048828f, 191.7482757568359f, 
197.9413299560547f, 201.9744110107422f, 204.6788635253906f, 206.5121154785156f, 
208.3563079833984f, 209.0952453613281f, 209.3571777343750f, 208.9805603027344f, 
208.1353759765625f, 207.3684844970703f, 206.7562103271484f, 206.2508850097656f, 
205.9447479248047f, 205.8604583740234f, 204.6592712402344f, 202.6399688720703f, 
199.7153625488281f, 195.7463378906250f, 191.1268768310547f, 186.2112579345703f, 
181.6002349853516f, 178.9270019531250f, 178.0906219482422f, 178.3910217285156f, 
177.1901092529297f, 173.9102478027344f, 168.2518615722656f, 160.5976104736328f, 
152.1742706298828f, 142.4298095703125f, 136.4926300048828f, 131.6634063720703f, 
128.9111938476563f, 128.5836334228516f, 129.2254943847656f, 133.2222900390625f, 
137.0889739990234f, 138.8344116210938f, 138.8866729736328f, 138.1311187744141f, 
135.8213195800781f, 131.6275329589844f, 124.7729187011719f, 115.6852722167969f, 
106.3275604248047f, 97.7201461791992f, 92.5079421997070f, 90.2174072265625f, 
89.1273803710938f, 89.0563735961914f, 89.1893768310547f, 89.2572555541992f, 
89.2951354980469f, 89.8281402587891f, 90.8873672485352f, 92.4143371582031f, 
94.6329498291016f, 96.2755203247070f, 96.6841583251953f, 
};

ML_SpectralDensity ML_Venera11_162 = {
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
69.3289337158203f, 74.6118545532227f, 80.3257904052734f, 86.4768066406250f, 
93.1001052856445f, 100.2289505004883f, 107.9051132202148f, 116.1687164306641f, 
125.0644378662109f, 134.6431732177734f, 144.9530487060547f, 152.1889648437500f, 
163.3925781250000f, 175.9825439453125f, 189.9207763671875f, 204.8109588623047f, 
217.8246612548828f, 231.7831573486328f, 247.2216796875000f, 262.1100769042969f, 
275.7065124511719f, 286.6281738281250f, 298.2355346679688f, 308.6993713378906f, 
317.2845153808594f, 324.2751159667969f, 331.5311584472656f, 335.6611022949219f, 
338.4302978515625f, 340.6479492187500f, 341.3361816406250f, 340.7544860839844f, 
340.0108947753906f, 338.2261352539063f, 335.8829040527344f, 333.0346374511719f, 
329.9282226562500f, 327.0925292968750f, 323.9711914062500f, 320.9389038085938f, 
318.3675231933594f, 315.9929199218750f, 312.9908447265625f, 308.3563232421875f, 
302.9342651367188f, 296.3717346191406f, 288.7480163574219f, 281.3021240234375f, 
274.2943725585938f, 268.1767272949219f, 265.1970825195313f, 263.9938049316406f, 
262.2646789550781f, 257.8024597167969f, 250.8009643554688f, 240.7787170410156f, 
230.2847595214844f, 219.3279724121094f, 209.8757171630859f, 203.2767333984375f, 
198.9181365966797f, 197.2313232421875f, 197.0564575195313f, 197.4641113281250f, 
197.8442230224609f, 198.1833953857422f, 198.0050354003906f, 196.5275268554688f, 
193.9667816162109f, 188.0091705322266f, 179.9958648681641f, 168.8838806152344f, 
156.5214691162109f, 145.3936614990234f, 137.7867431640625f, 133.5764465332031f, 
131.4213104248047f, 130.3873443603516f, 129.6771392822266f, 128.9649047851563f, 
128.3854522705078f, 128.3448791503906f, 128.8155822753906f, 129.8925170898438f, 
131.4554595947266f, 132.4899597167969f, 132.3951416015625f, 
};

ML_SpectralDensity ML_Venera11_237 = {
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
163.0815429687500f, 170.8244934082031f, 178.9752349853516f, 187.5141906738281f, 
196.4622344970703f, 205.8350219726563f, 215.6567840576172f, 225.9466552734375f, 
236.7266387939453f, 248.0230560302734f, 259.8556823730469f, 267.8827514648438f, 
280.7392883300781f, 294.6671752929688f, 308.6408081054688f, 323.6221313476563f, 
336.8989868164063f, 349.7339477539063f, 362.9396057128906f, 375.1243286132813f, 
387.2769775390625f, 396.5361633300781f, 405.4095153808594f, 413.6255798339844f, 
420.4838867187500f, 425.6582336425781f, 428.5494079589844f, 430.1887512207031f, 
429.6323547363281f, 428.6672058105469f, 425.5140991210938f, 422.0432739257813f, 
417.3143005371094f, 412.2404174804688f, 407.0355224609375f, 402.7599792480469f, 
397.8481445312500f, 394.0824279785156f, 390.2658081054688f, 385.9298706054688f, 
382.6383056640625f, 379.3744201660156f, 375.0840148925781f, 369.3187866210938f, 
361.4179687500000f, 352.7572937011719f, 342.7475280761719f, 333.2415771484375f, 
324.9688720703125f, 318.2100219726563f, 313.5085449218750f, 310.4644775390625f, 
307.6968994140625f, 302.7357788085938f, 295.2323608398438f, 284.9907531738281f, 
275.0024414062500f, 263.7935791015625f, 254.1057739257813f, 246.9177856445313f, 
240.9966430664063f, 237.7115325927734f, 234.9036560058594f, 232.7840270996094f, 
231.0810852050781f, 229.9534301757813f, 228.7232513427734f, 226.8190002441406f, 
223.9866027832031f, 218.4021759033203f, 210.2440948486328f, 199.0089416503906f, 
186.1522521972656f, 174.1404571533203f, 164.9839935302734f, 159.3445129394531f, 
156.2708892822266f, 154.8769836425781f, 153.5985565185547f, 152.4906158447266f, 
151.4277496337891f, 150.9887695312500f, 150.8276214599609f, 151.0440826416016f, 
151.4446716308594f, 151.6188354492188f, 151.0125732421875f, 
};

ML_SpectralDensity ML_Venera11_377 = {
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
336.1604003906250f, 346.4671325683594f, 357.1422119140625f, 368.1452941894531f, 
379.4895019531250f, 391.1804809570313f, 403.2338562011719f, 415.6579589843750f, 
428.4638366699219f, 441.6667175292969f, 455.2731933593750f, 465.4540710449219f, 
479.3243713378906f, 493.7886352539063f, 506.3296813964844f, 517.0447387695313f, 
527.1596069335938f, 536.1258544921875f, 543.5343627929688f, 550.0176391601563f, 
555.1928710937500f, 558.7474975585938f, 560.7023315429688f, 561.2365112304688f, 
560.7464599609375f, 559.7979736328125f, 556.8595581054688f, 551.5872802734375f, 
544.2988891601563f, 535.5794677734375f, 525.6268920898438f, 515.8616333007813f, 
506.1788024902344f, 494.8086242675781f, 485.3645935058594f, 477.2481384277344f, 
471.3355102539063f, 465.9440002441406f, 460.5997924804688f, 455.3612976074219f, 
450.6788635253906f, 446.0133972167969f, 439.7545776367188f, 432.0183105468750f, 
422.2624206542969f, 411.3922119140625f, 399.9425964355469f, 388.6145324707031f, 
378.7023925781250f, 371.3183593750000f, 365.4731445312500f, 360.3965454101563f, 
355.3630981445313f, 349.0037841796875f, 341.5583190917969f, 332.7865600585938f, 
323.6680603027344f, 313.1502990722656f, 304.5351562500000f, 297.8777770996094f, 
291.7662963867188f, 286.1412048339844f, 280.0963439941406f, 274.2136230468750f, 
269.0977783203125f, 264.3115234375000f, 260.3864746093750f, 256.4470214843750f, 
252.4308166503906f, 246.5681304931641f, 240.2569885253906f, 231.4870758056641f, 
218.6643524169922f, 206.9500732421875f, 197.5857391357422f, 190.5561828613281f, 
185.8835144042969f, 183.5403747558594f, 181.6700592041016f, 179.7577209472656f, 
178.1116027832031f, 176.9338684082031f, 175.9772644042969f, 175.1673431396484f, 
174.1922302246094f, 172.9401702880859f, 171.2913055419922f, 
};

ML_SpectralDensity ML_Venera11_486 = {
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
361.0508728027344f, 373.0881958007813f, 385.5881958007813f, 398.5059814453125f, 
411.8589782714844f, 425.6561584472656f, 439.9181518554688f, 454.6571960449219f, 
469.8888854980469f, 485.6337890625000f, 501.9023437500000f, 514.7336425781250f, 
531.4433593750000f, 548.7018432617188f, 563.1466674804688f, 573.6550292968750f, 
583.3555297851563f, 592.5153808593750f, 601.4141235351563f, 608.6943969726563f, 
614.1123657226563f, 617.7197875976563f, 619.5156250000000f, 618.2125854492188f, 
614.5865478515625f, 609.6989135742188f, 602.5556030273438f, 594.6294555664063f, 
585.2559204101563f, 575.4600219726563f, 564.0883789062500f, 551.9868164062500f, 
541.5497436523438f, 530.2617187500000f, 522.4237060546875f, 514.2289428710938f, 
506.6693420410156f, 500.2748413085938f, 493.1835327148438f, 485.9837951660156f, 
479.0921936035156f, 470.8888854980469f, 462.8778991699219f, 453.4944763183594f, 
443.5351867675781f, 432.8656005859375f, 421.9091491699219f, 412.0095214843750f, 
402.4611816406250f, 393.4325561523438f, 387.5496215820313f, 383.8561096191406f, 
380.4463500976563f, 375.4335327148438f, 368.5727233886719f, 359.2666015625000f, 
348.4097900390625f, 336.4002685546875f, 324.6949462890625f, 316.0705566406250f, 
309.0790405273438f, 304.6963806152344f, 300.6466064453125f, 296.2043457031250f, 
291.3041076660156f, 285.5432434082031f, 279.3084411621094f, 273.4436035156250f, 
266.8871765136719f, 258.4973144531250f, 248.4845886230469f, 238.8181610107422f, 
230.0202484130859f, 221.9810485839844f, 213.8431091308594f, 207.6132659912109f, 
202.5420989990234f, 198.9353790283203f, 196.1158905029297f, 194.0046539306641f, 
192.0865936279297f, 190.6033325195313f, 189.6044464111328f, 189.1246337890625f, 
188.7989959716797f, 188.0345306396484f, 186.4326934814453f, 
};

ML_SpectralDensity ML_Venera11_511 = {
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
512.7476196289063f, 522.7415771484375f, 532.9802246093750f, 543.4186401367188f, 
554.0634155273438f, 564.9141845703125f, 575.9794311523438f, 587.2608642578125f, 
598.7623901367188f, 610.4912719726563f, 622.4470825195313f, 634.0526123046875f, 
646.3865966796875f, 659.5075683593750f, 671.0207519531250f, 679.1824340820313f, 
688.6021728515625f, 698.3911132812500f, 707.8794555664063f, 714.3681030273438f, 
718.5143432617188f, 720.9843139648438f, 721.2152709960938f, 719.3210449218750f, 
716.8303833007813f, 711.1995239257813f, 704.9725341796875f, 694.9288940429688f, 
683.9080200195313f, 671.8973999023438f, 657.0526123046875f, 643.0708007812500f, 
631.3302001953125f, 622.7111816406250f, 615.6292724609375f, 609.6785278320313f, 
602.9391479492188f, 596.6362915039063f, 589.9099731445313f, 581.6760253906250f, 
574.6509399414063f, 567.4446411132813f, 559.1452026367188f, 549.6550903320313f, 
539.9683227539063f, 530.9935302734375f, 522.3281860351563f, 512.9453735351563f, 
504.0101928710938f, 496.3333740234375f, 489.4014587402344f, 483.1847839355469f, 
477.2365722656250f, 470.4863281250000f, 462.7319030761719f, 454.6384582519531f, 
446.3745727539063f, 438.4006958007813f, 430.6026306152344f, 423.6518249511719f, 
417.8429260253906f, 412.9228210449219f, 408.5455627441406f, 403.6853637695313f, 
399.0003967285156f, 394.5638732910156f, 389.6445312500000f, 384.5600585937500f, 
379.0624694824219f, 371.8325805664063f, 363.0304565429688f, 352.2024536132813f, 
340.6506042480469f, 329.3891296386719f, 319.7914123535156f, 312.1740417480469f, 
306.0195922851563f, 301.8602905273438f, 298.9320983886719f, 296.4274902343750f, 
294.2778015136719f, 291.7517700195313f, 289.0052185058594f, 286.3060607910156f, 
284.2008972167969f, 281.5361938476563f, 279.3271789550781f, 
};

ML_SpectralDensity ML_Venera11_561 = {
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
541.6168823242188f, 554.7609863281250f, 568.2902221679688f, 582.1482543945313f, 
596.3468627929688f, 610.8883056640625f, 625.7871093750000f, 641.0484619140625f, 
656.6807861328125f, 672.6972045898438f, 689.1004638671875f, 704.3121337890625f, 
718.6840820312500f, 736.2269287109375f, 752.4652099609375f, 766.5521240234375f, 
778.9356079101563f, 788.6569213867188f, 795.1631469726563f, 800.2462768554688f, 
804.4552612304688f, 807.3598022460938f, 807.8225097656250f, 805.8179321289063f, 
802.3788452148438f, 796.0286865234375f, 787.6839599609375f, 776.8837280273438f, 
763.1860351562500f, 748.8508300781250f, 733.8710937500000f, 718.8386840820313f, 
705.3757934570313f, 693.5076293945313f, 685.0839843750000f, 677.3989868164063f, 
670.6698608398438f, 663.4845581054688f, 656.7745361328125f, 650.0254516601563f, 
642.7100830078125f, 634.1230468750000f, 625.7033691406250f, 617.2557983398438f, 
608.2528686523438f, 599.8612670898438f, 592.0787353515625f, 585.4542236328125f, 
579.6932373046875f, 574.4101562500000f, 568.8999633789063f, 564.0980834960938f, 
558.7735595703125f, 552.9705200195313f, 545.4848022460938f, 538.2055664062500f, 
530.7058715820313f, 523.9788208007813f, 518.6316528320313f, 513.6226806640625f, 
509.4734802246094f, 505.7113037109375f, 502.8280029296875f, 500.0353698730469f, 
496.6470642089844f, 492.3317260742188f, 487.5890197753906f, 481.9772033691406f, 
474.5484924316406f, 466.2550659179688f, 456.7453613281250f, 447.1746215820313f, 
437.0695495605469f, 427.9923706054688f, 420.0506591796875f, 413.1208496093750f, 
407.3994750976563f, 402.2223510742188f, 397.9039611816406f, 393.5854797363281f, 
389.6205139160156f, 386.7946166992188f, 384.1600646972656f, 382.2314758300781f, 
379.8375549316406f, 377.8297424316406f, 374.6924743652344f, 
};

ML_SpectralDensity ML_Venera11_621 = {
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
934.7803344726563f, 942.3978881835938f, 950.1150512695313f, 957.8947753906250f, 
965.7396850585938f, 973.6469116210938f, 981.6203613281250f, 989.6586303710938f, 
997.7620849609375f, 1005.9334106445312f, 1014.1697998046875f, 1021.9335937500000f, 
1030.1210937500000f, 1039.3040771484375f, 1044.8166503906250f, 1046.8864746093750f, 
1048.5610351562500f, 1052.9532470703125f, 1055.1226806640625f, 1055.7518310546875f, 
1053.5324707031250f, 1050.8411865234375f, 1047.7344970703125f, 1042.5377197265625f, 
1035.2701416015625f, 1029.9384765625000f, 1024.2741699218750f, 1019.2406616210937f, 
1014.0550537109375f, 1010.8923950195312f, 1004.2780761718750f, 996.7069091796875f, 
987.7683715820313f, 978.0752563476563f, 969.2519531250000f, 960.5021362304688f, 
952.4762573242188f, 944.1843261718750f, 936.9082641601563f, 929.6488037109375f, 
922.5667724609375f, 916.1752319335938f, 908.7924194335938f, 901.7980957031250f, 
894.5284423828125f, 888.1619873046875f, 881.3390502929688f, 875.3092651367188f, 
868.3131713867188f, 860.9246826171875f, 854.4833984375000f, 846.9782104492188f, 
840.3634643554688f, 833.4407958984375f, 825.8519287109375f, 819.3421630859375f, 
813.1521606445313f, 807.0564575195313f, 801.7816162109375f, 796.1551513671875f, 
790.7193603515625f, 785.2715454101563f, 779.7985229492188f, 773.2409057617188f, 
765.2770385742188f, 756.8049926757813f, 747.4616088867188f, 738.5473632812500f, 
729.5880737304688f, 720.0014648437500f, 709.4889526367188f, 697.2947387695313f, 
682.9015502929688f, 668.4013061523438f, 655.8574829101563f, 644.0608520507813f, 
633.6894531250000f, 626.1138305664063f, 621.0673217773438f, 617.5806274414063f, 
614.7899780273438f, 611.9105224609375f, 608.2891845703125f, 604.2615356445313f, 
600.0896606445313f, 594.7891845703125f, 588.2678222656250f, 
};

ML_SpectralDensity ML_Venera11_Sun = {
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
978.3179321289063f, 994.1447143554688f, 1010.3063964843750f, 1026.7294921875000f, 
1043.4226074218750f, 1060.3831787109375f, 1077.6225585937500f, 1095.1413574218750f, 
1112.9434814453125f, 1131.0383300781250f, 1149.4230957031250f, 1166.6159667968750f, 
1193.9080810546875f, 1215.4260253906250f, 1218.8944091796875f, 1205.4306640625000f, 
1192.7509765625000f, 1184.1165771484375f, 1172.7536621093750f, 1162.8565673828125f, 
1153.7912597656250f, 1143.4482421875000f, 1135.5789794921875f, 1131.6756591796875f, 
1126.7154541015625f, 1122.1539306640625f, 1116.9654541015625f, 1110.1583251953125f, 
1105.1032714843750f, 1101.2801513671875f, 1096.3297119140625f, 1091.0190429687500f, 
1085.0891113281250f, 1079.4855957031250f, 1074.5708007812500f, 1070.3306884765625f, 
1065.0505371093750f, 1059.0333251953125f, 1054.5787353515625f, 1049.8293457031250f, 
1045.2545166015625f, 1041.0120849609375f, 1032.2506103515625f, 1023.2916259765625f, 
1012.6386718750000f, 1002.1814575195312f, 992.1271972656250f, 981.6460571289063f, 
971.8814086914063f, 962.5303955078125f, 952.6730957031250f, 943.2692871093750f, 
933.6368408203125f, 923.4637451171875f, 913.0635986328125f, 903.6446533203125f, 
892.6660156250000f, 883.0415039062500f, 873.3912963867188f, 863.4489135742188f, 
852.8961181640625f, 843.6724853515625f, 834.0812988281250f, 825.2114257812500f, 
817.0713500976563f, 807.0634765625000f, 798.5275878906250f, 789.6691284179688f, 
780.1351928710938f, 770.9586791992188f, 763.1088256835938f, 754.0242919921875f, 
745.6994628906250f, 737.3438110351563f, 727.8187255859375f, 719.2243652343750f, 
710.4271850585938f, 701.5769653320313f, 693.1331787109375f, 684.4362792968750f, 
675.0478515625000f, 666.6995849609375f, 659.5302124023438f, 652.9783325195313f, 
647.4387817382813f, 640.6473999023438f, 633.8024902343750f, 
};

ML_SpectralDensity ML_Venera13_000 = {
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f, 
 0.0000000000000f,  0.0000000000000f,  3.4950106143951f,  4.4060382843018f, 
 4.8104062080383f,  5.2515540122986f,  5.7332673072815f,  6.2593512535095f, 
 6.8331832885742f,  7.4603042602539f,  8.1444664001465f,  8.8915395736694f, 
 9.7074260711670f, 10.3655605316162f, 11.1077995300293f, 12.2575922012329f, 
13.4115695953369f, 14.6272811889648f, 15.8955869674683f, 17.0919628143311f, 
18.4019432067871f, 19.6144371032715f, 20.7819480895996f, 22.2115917205811f, 
23.4218540191650f, 24.6335506439209f, 25.8852233886719f, 26.9293193817139f, 
27.9331703186035f, 28.8875293731689f, 29.6829013824463f, 30.2966213226318f, 
30.8392238616943f, 31.1434288024902f, 31.3822994232178f, 31.4634742736816f, 
31.2861671447754f, 31.0386753082275f, 30.6403923034668f, 29.9882736206055f, 
29.4057159423828f, 28.5112934112549f, 27.4797763824463f, 26.0864334106445f, 
24.6803417205811f, 23.4389095306396f, 22.5525474548340f, 21.9879989624023f, 
21.7640323638916f, 21.8337898254395f, 22.3733634948730f, 23.2993831634521f, 
24.5243453979492f, 26.0819816589355f, 27.8464298248291f, 29.3135871887207f, 
30.1867179870605f, 30.3464508056641f, 29.8531074523926f, 28.0867424011230f, 
25.6331367492676f, 23.7873210906982f, 22.8946971893311f, 22.5610466003418f, 
22.5190048217773f, 22.7417545318604f, 23.1762046813965f, 23.6999473571777f, 
24.3467369079590f, 25.1190052032471f, 25.7956695556641f, 
};

ML_SpectralDensity ML_Illuminant_D65 = {
46.8617591857910f, 49.3637008666992f, 51.6352233886719f, 51.0323028564453f, 
50.3826789855957f, 52.3118400573730f, 56.0542373657227f, 68.7015457153320f, 
81.5923614501953f, 87.1204376220703f, 91.0789031982422f, 92.4589080810547f, 
92.9100646972656f, 90.0570220947266f, 88.1782226562500f, 95.7736358642578f, 
104.5026016235352f, 110.9364013671875f, 116.3275985717773f, 117.4099960327148f, 
117.5868072509766f, 116.3364028930664f, 115.1016006469727f, 115.3919906616211f, 
115.4328002929688f, 112.3670043945313f, 109.2702026367188f, 109.0823974609375f, 
109.2283935546875f, 108.5780029296875f, 107.7144012451172f, 106.2959976196289f, 
105.1446075439453f, 106.2393951416016f, 107.3179931640625f, 106.0469970703125f, 
104.5803985595703f, 104.2254028320313f, 103.8248062133789f, 102.0229949951172f, 
100.0228424072266f, 98.1670989990234f, 96.5214080810547f, 96.0611038208008f, 
95.3946151733398f, 92.2368011474609f, 89.1909790039063f, 89.3459014892578f, 
89.9025421142578f, 89.8026428222656f, 89.5094985961914f, 88.6489028930664f, 
87.5481185913086f, 85.4936370849609f, 83.5778427124023f, 83.4938964843750f, 
83.4542160034180f, 81.8629989624023f, 80.2584152221680f, 80.1206970214844f, 
80.3271026611328f, 81.2462005615234f, 81.9143981933594f, 80.2809982299805f, 
78.0100326538086f, 74.0027389526367f, 70.3483657836914f, 70.6651992797852f, 
71.6602172851563f, 72.9790420532227f, 73.4198989868164f, 67.9765014648438f, 
62.8656005859375f, 65.7447967529297f, 69.7007827758789f, 72.4863052368164f, 
74.0852584838867f, 69.3398361206055f, 63.2518806457520f, 55.0054206848145f, 
48.6718826293945f, 56.6117973327637f, 65.3768157958984f, 65.0941009521484f, 
63.6434211730957f, 63.8434028625488f, 63.9576187133789f, 61.8779411315918f, 
59.2934379577637f, 55.7054443359375f, 52.7374801635742f, 54.6998062133789f, 
57.2840194702148f, 58.8765373229980f, 60.0252990722656f, 
};

ML_SpectralDensity ML_SolarSpectrum = {
882.5198974609375f, 1099.1999511718750f, 1186.4000244140625f, 1023.4801025390625f, 
1214.8000488281250f, 907.1599731445313f, 1172.0000000000000f, 977.6599731445313f, 
1537.5200195312500f, 1706.5999755859375f, 1716.8000488281250f, 1801.8000488281250f, 
1750.4000244140625f, 1719.5999755859375f, 1516.1999511718750f, 1769.1999511718750f, 
1798.8000488281250f, 1948.0000000000000f, 2101.1999511718750f, 2034.4000244140625f, 
2069.8000488281250f, 2051.1999511718750f, 2012.4000244140625f, 2045.4000244140625f, 
2092.0000000000000f, 1925.8000488281250f, 1954.5999755859375f, 1990.5999755859375f, 
1935.5000000000000f, 1965.8000488281250f, 1978.0000000000000f, 1865.5999755859375f, 
1819.8000488281250f, 1887.8000488281250f, 1944.8000488281250f, 1914.8000488281250f, 
1885.5999755859375f, 1904.8000488281250f, 1897.0000000000000f, 1901.1999511718750f, 
1853.1999511718750f, 1870.1999511718750f, 1861.5999755859375f, 1894.4000244140625f, 
1863.5999755859375f, 1872.4000244140625f, 1792.4000244140625f, 1823.8000488281250f, 
1782.0000000000000f, 1784.0000000000000f, 1737.8000488281250f, 1684.4000244140625f, 
1707.8000488281250f, 1683.0000000000000f, 1688.5999755859375f, 1670.0000000000000f, 
1654.8000488281250f, 1631.6666259765625f, 1619.4000244140625f, 1498.0000000000000f, 
1578.4000244140625f, 1564.8000488281250f, 1553.8000488281250f, 1530.0000000000000f, 
1514.0000000000000f, 1502.8000488281250f, 1503.5999755859375f, 1490.1999511718750f, 
1470.8000488281250f, 1465.1666259765625f, 1426.4000244140625f, 1405.8000488281250f, 
1387.0000000000000f, 1379.5000000000000f, 1352.1999511718750f, 1347.0000000000000f, 
1303.0000000000000f, 1317.0000000000000f, 1300.2500000000000f, 1282.5000000000000f, 
1262.6666259765625f, 1241.0000000000000f, 1215.0000000000000f, 1195.5000000000000f, 
1208.3333740234375f, 1204.5000000000000f, 1171.3333740234375f, 1163.5000000000000f, 
1144.6666259765625f, 1133.5000000000000f, 1098.6666259765625f, 1095.0000000000000f, 
1066.0000000000000f, 1070.0000000000000f, 1046.6666259765625f, 
};

ML_SpectralDensity ML_BlueSky = {
 0.0956052988768f,  0.2541843354702f,  0.3252299726009f,  0.4483328461647f, 
 0.4355113506317f,  0.4509140849113f,  0.4650315940380f,  0.5206413269043f, 
 0.6212207078934f,  0.7067369222641f,  0.7453636527061f,  0.7602809071541f, 
 0.7711473703384f,  0.7408665418625f,  0.7251814603806f,  0.7459014654160f, 
 0.8027479648590f,  0.8498325943947f,  0.8784006237984f,  0.8852601051331f, 
 0.8830512762070f,  0.8864428997040f,  0.8831475377083f,  0.8782640695572f, 
 0.8658743500710f,  0.8426309227943f,  0.8253563642502f,  0.8186481595039f, 
 0.8156746625900f,  0.7987827658653f,  0.7759097814560f,  0.7455062270164f, 
 0.7288818955421f,  0.7330437302589f,  0.7461776733398f,  0.7425632476807f, 
 0.7210410833359f,  0.7023914456367f,  0.6890467405319f,  0.6757579445839f, 
 0.6662161946297f,  0.6570122241974f,  0.6512525081635f,  0.6400557756424f, 
 0.6289314627647f,  0.6142150163651f,  0.5898039340973f,  0.5705909729004f, 
 0.5671805739403f,  0.5647823810577f,  0.5580387711525f,  0.5457829236984f, 
 0.5380914211273f,  0.5244552493095f,  0.5093741416931f,  0.4927855134010f, 
 0.4820905923843f,  0.4736015200615f,  0.4626024067402f,  0.4551009833813f, 
 0.4522045254707f,  0.4522112607956f,  0.4486874043941f,  0.4407587647438f, 
 0.4283075630665f,  0.4000492691994f,  0.3832249939442f,  0.3718059360981f, 
 0.3708406388760f,  0.3673564791679f,  0.3608773648739f,  0.3401287794113f, 
 0.3167705535889f,  0.3021306097507f,  0.3012179434299f,  0.3067002296448f, 
 0.3157644271851f,  0.3134954571724f,  0.3084856867790f,  0.2793076038361f, 
 0.2398619055748f,  0.2284522354603f,  0.2527685165405f,  0.2742228507996f, 
 0.2752481997013f,  0.2729110717773f,  0.2642267048359f,  0.2577526569366f, 
 0.2510012090206f,  0.2403025329113f,  0.2266324609518f,  0.2087972909212f, 
 0.2019419074059f,  0.1963204145432f,  0.1996869891882f, 
};

ML_SpectralDensity ML_AlbedoMars = {
 0.0444121435285f,  0.0455990731716f,  0.0463592968881f,  0.0473471283913f,
 0.0473360083997f,  0.0497060418129f,  0.0507811307907f,  0.0519857443869f,
 0.0540833622217f,  0.0557661317289f,  0.0582254305482f,  0.0594133809209f,
 0.0615872777998f,  0.0636164993048f,  0.0661344081163f,  0.0685027688742f,
 0.0710719749331f,  0.0735214427114f,  0.0768102705479f,  0.0793889313936f,
 0.0823791995645f,  0.0855076387525f,  0.0883682817221f,  0.0916211158037f,
 0.0947303920984f,  0.0979286804795f,  0.1008570343256f,  0.1041895672679f,
 0.1074512600899f,  0.1111317649484f,  0.1146685183048f,  0.1188438907266f,
 0.1228942573071f,  0.1270616948605f,  0.1314250081778f,  0.1363826394081f,
 0.1412972509861f,  0.1468490660191f,  0.1537520736456f,  0.1616970151663f,
 0.1716659963131f,  0.1805673092604f,  0.1892928481102f,  0.1965759992599f,
 0.2041985094547f,  0.2105647921562f,  0.2169668972492f,  0.2222914397717f,
 0.2274606376886f,  0.2320930510759f,  0.2368943244219f,  0.2415092885494f,
 0.2452704459429f,  0.2483697384596f,  0.2522217035294f,  0.2557595372200f,
 0.2589770853519f,  0.2617872953415f,  0.2646272480488f,  0.2671374976635f,
 0.2695948183537f,  0.2720769643784f,  0.2743667066097f,  0.2758504450321f,
 0.2776886522770f,  0.2801114320755f,  0.2816069424152f,  0.2836902737617f,
 0.2853116393089f,  0.2868651449680f,  0.2882736623287f,  0.2896606922150f,
 0.2908309102058f,  0.2918810248375f,  0.2931383550167f,  0.2941007018089f,
 0.2950029969215f,  0.2958471477032f,  0.2967388331890f,  0.2973304986954f,
 0.2980076670647f,  0.2986100316048f,  0.2989281713963f,  0.2993769943714f,
 0.2996523976326f,  0.2999357581139f,  0.3001045584679f,  0.3002387285233f,
 0.3002394139767f,  0.3002384603024f,  0.3001035749912f,  0.2999349832535f,
 0.3002389669418f,  0.2999349534512f,  0.2999349832535f,
};

ML_SpectralDensity ML_AlbedoVenus = {
 0.1137460619211f,  0.1133410260081f,  0.1114660054445f,  0.1095454394817f,
 0.1106775701046f,  0.1129597947001f,  0.1180993765593f,  0.1258796155453f,
 0.1347882598639f,  0.1428026407957f,  0.1506072133780f,  0.1569005846977f,
 0.1613186299801f,  0.1651804745197f,  0.1688438355923f,  0.1730578392744f,
 0.1760085225105f,  0.1784681975842f,  0.1825992166996f,  0.1859342902899f,
 0.1881722062826f,  0.1893689632416f,  0.1908621042967f,  0.1922747194767f,
 0.1928874701262f,  0.1945480704308f,  0.1962348967791f,  0.1962979733944f,
 0.1972590684891f,  0.1983463913202f,  0.1988587677479f,  0.1994663476944f,
 0.1993396133184f,  0.1982317417860f,  0.1977133899927f,  0.1976707428694f,
 0.1977879256010f,  0.1978515982628f,  0.1974700987339f,  0.1968590766191f,
 0.1962061524391f,  0.1959082782269f,  0.1953814923763f,  0.1930866092443f,
 0.1922950595617f,  0.1915534436703f,  0.1904397755861f,  0.1891153156757f,
 0.1882750540972f,  0.1880159527063f,  0.1883077621460f,  0.1882649809122f,
 0.1874666810036f,  0.1874866336584f,  0.1869596838951f,  0.1854841113091f,
 0.1851366758347f,  0.1844804286957f,  0.1829353570938f,  0.1819037646055f,
 0.1818179786205f,  0.1824352443218f,  0.1828412711620f,  0.1828113645315f,
 0.1825690865517f,  0.1818974316120f,  0.1814941167831f,  0.1811762452126f,
 0.1808367222548f,  0.1799072623253f,  0.1786752641201f,  0.1773234903812f,
 0.1767576932907f,  0.1764613091946f,  0.1765931099653f,  0.1770588606596f,
 0.1775839328766f,  0.1772814393044f,  0.1768646836281f,  0.1766058206558f,
 0.1767296493053f,  0.1764741539955f,  0.1760075539351f,  0.1752029359341f,
 0.1735697686672f,  0.1728819906712f,  0.1732082813978f,  0.1739709079266f,
 0.1746332049370f,  0.1749852299690f,  0.1758794337511f,  0.1747659593821f,
 0.1736490726471f,  0.1737406849861f,  0.1743721514940f,
};

ML_SpectralDensity ML_AlbedoMercury = {
 0.6353623867035f,  0.6196553111076f,  0.6195921897888f,  0.6194387078285f,
 0.6195932030678f,  0.6618647575378f,  0.7289669513702f,  0.7320570945740f,
 0.7324585914612f,  0.7320570945740f,  0.7296822071075f,  0.7238752841949f,
 0.7180691361427f,  0.7135618925095f,  0.7135591506958f,  0.7135590314865f,
 0.7135589122772f,  0.7505147457123f,  0.7957915067673f,  0.8410646319389f,
 0.8556934595108f,  0.8556934595108f,  0.8556934595108f,  0.8575626015663f,
 0.8726295828819f,  0.8877282738686f,  0.9021537303925f,  0.9023917317390f,
 0.9023897647858f,  0.9023917317390f,  0.9153583645821f,  0.9397447109222f,
 0.9641277790070f,  0.9793297648430f,  0.9793304800987f,  0.9793298244476f,
 0.9794343709946f,  0.9910805225372f,  1.0038315057755f,  1.0165807008743f,
 1.0205419063568f,  1.0205419063568f,  1.0205419063568f,  1.0257303714752f,
 1.0479670763016f,  1.0702066421509f,  1.0906473398209f,  1.0908449888229f,
 1.0908449888229f,  1.0908449888229f,  1.1010464429855f,  1.1175454854965f,
 1.1340466737747f,  1.1441782712936f,  1.1441782712936f,  1.1441782712936f,
 1.1444970369339f,  1.1750508546829f,  1.2083213329315f,  1.2415959835052f,
 1.2493183612823f,  1.2490353584290f,  1.2493370771408f,  1.2512943744659f,
 1.2480773925781f,  1.2448602914810f,  1.2364937067032f,  1.2377885580063f,
 1.2378284931183f,  1.2378324270248f,  1.2333401441574f,  1.2143507003784f,
 1.1953594684601f,  1.1853903532028f,  1.1853903532028f,  1.1853903532028f,
 1.1748974323273f,  1.2328470945358f,  1.2951987981796f,  1.3574513196945f,
 1.3642581701279f,  1.3655560016632f,  1.3651579618454f,  1.3668328523636f,
 1.3820129632950f,  1.3971910476685f,  1.4090573787689f,  1.4090579748154f,
 1.4090582132339f,  1.4090584516525f,  1.4172235727310f,  1.4280521869659f,
 1.4388792514801f,  1.4438338279724f,  1.4438368082047f,
};

ML_SpectralDensity ML_VenusUVcontrast = {
 0.2557919025421f,  0.2371650338173f,  0.2178729325533f,  0.1966274678707f,
 0.1707612574100f,  0.1446800678968f,  0.1161368861794f,  0.0991007015109f,
 0.0909952297807f,  0.0807544142008f,  0.0638795346022f,  0.0642958730459f,
 0.0484165996313f,  0.0474280230701f,  0.0385190546513f,  0.0379100851715f,
 0.0304844286293f,  0.0299910902977f,  0.0243423990905f,  0.0240692030638f,
 0.0176406316459f,  0.0159956496209f,  0.0124749643728f,  0.0109324669465f,
 0.0087981484830f,  0.0073923571035f,  0.0063837291673f,  0.0000000000000f,
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f,
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f,
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f,
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f,
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f,
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f,
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f,
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f,
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f,
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f,
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f,
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f,
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f,
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f,
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f,
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,  0.0000000000000f,
 0.0000000000000f,  0.0000000000000f,  0.0000000000000f,
};

ML_SpectralDensity sdBluChannel = {   // x,y = 0.150559782981872, 0.0327584072947502
    0,   // 360 nm
    0,   // 365 nm
    0,   // 370 nm
    0,   // 375 nm
    0,   // 380 nm
    0.00294117652811110,   // 385 nm
    0.0247549023479223,   // 390 nm
    0.0617647059261798,   // 395 nm
    0.120361991226673,   // 400 nm
    0.212009802460670,   // 405 nm
    0.368656724691390,   // 410 nm
    0.547295928001403,   // 415 nm
    0.715449452400207,   // 420 nm
    0.848774552345275,   // 425 nm
    0.931862771511077,   // 430 nm
    0.966559708118438,   // 435 nm
    0.986029386520385,   // 440 nm
    0.988480389118194,   // 445 nm
    0.977205872535705,   // 450 nm
    0.944852948188781,   // 455 nm
    0.881107032299041,   // 460 nm
    0.764406383037567,   // 465 nm
    0.579518079757690,   // 470 nm
    0.445965498685836,   // 475 nm
    0.347549021244049,   // 480 nm
    0.270343154668807,   // 485 nm
    0.204112216830253,   // 490 nm
    0.152941182255744,   // 495 nm
    0.103431373834609,   // 500 nm
    0.0625000000000000,   // 505 nm
    0.0286764707416296,   // 510 nm
    0.00416666688397526,   // 515 nm
    0,   // 520 nm
    0,   // 525 nm
    0,   // 530 nm
    0,   // 535 nm
    0,   // 540 nm
    0,   // 545 nm
    0,   // 550 nm
    0,   // 555 nm
    0,   // 560 nm
    0,   // 565 nm
    0,   // 570 nm
    0,   // 575 nm
    0,   // 580 nm
    0,   // 585 nm
    0,   // 590 nm
    0,   // 595 nm
    0,   // 600 nm
    0,   // 605 nm
    0,   // 610 nm
    0,   // 615 nm
    0,   // 620 nm
    0,   // 625 nm
    0,   // 630 nm
    0,   // 635 nm
    0,   // 640 nm
    0,   // 645 nm
    0,   // 650 nm
    0,   // 655 nm
    0,   // 660 nm
    0,   // 665 nm
    0,   // 670 nm
    0,   // 675 nm
    0,   // 680 nm
    0,   // 685 nm
    0,   // 690 nm
    0,   // 695 nm
    0,   // 700 nm
    0,   // 705 nm
    0,   // 710 nm
    0,   // 715 nm
    0,   // 720 nm
    0,   // 725 nm
    0,   // 730 nm
    0,   // 735 nm
    0,   // 740 nm
    0,   // 745 nm
    0,   // 750 nm
    0,   // 755 nm
    0,   // 760 nm
    0,   // 765 nm
    0,   // 770 nm
    0,   // 775 nm
    0,   // 780 nm
    0,   // 785 nm
    0,   // 790 nm
    0,   // 795 nm
    0,   // 800 nm
    0,   // 805 nm
    0,   // 810 nm
    0,   // 815 nm
    0,   // 820 nm
    0,   // 825 nm
    0,   // 830 nm
};
ML_SpectralDensity sdGrnChannel = {   // x,y = 0.297232359647750, 0.667712807655334
    0,   // 360 nm
    0,   // 365 nm
    0,   // 370 nm
    0,   // 375 nm
    0,   // 380 nm
    0,   // 385 nm
    0,   // 390 nm
    0,   // 395 nm
    0,   // 400 nm
    0,   // 405 nm
    0,   // 410 nm
    0,   // 415 nm
    0,   // 420 nm
    0,   // 425 nm
    0,   // 430 nm
    0,   // 435 nm
    0,   // 440 nm
    0,   // 445 nm
    0,   // 450 nm
    0,   // 455 nm
    0,   // 460 nm
    0,   // 465 nm
    0,   // 470 nm
    0,   // 475 nm
    0,   // 480 nm
    0.0147058824077248,   // 485 nm
    0.0470588244497776,   // 490 nm
    0.0874999985098838,   // 495 nm
    0.148039221763610,   // 500 nm
    0.245588228106498,   // 505 nm
    0.538540184497833,   // 510 nm
    0.815045237541198,   // 515 nm
    0.930147051811218,   // 520 nm
    0.967647075653076,   // 525 nm
    0.983088254928588,   // 530 nm
    0.985294103622436,   // 535 nm
    0.972794115543365,   // 540 nm
    0.950735270977020,   // 545 nm
    0.909852981567382,   // 550 nm
    0.857352912425994,   // 555 nm
    0.797860980033874,   // 560 nm
    0.713903725147247,   // 565 nm
    0.597794115543365,   // 570 nm
    0.471995383501052,   // 575 nm
    0.366318792104721,   // 580 nm
    0.276470601558685,   // 585 nm
    0.213970586657524,   // 590 nm
    0.135294124484062,   // 595 nm
    0.0801470652222633,   // 600 nm
    0.0345588251948356,   // 605 nm
    0.00441176490858197,   // 610 nm
    0,   // 615 nm
    0,   // 620 nm
    0,   // 625 nm
    0,   // 630 nm
    0,   // 635 nm
    0,   // 640 nm
    0,   // 645 nm
    0,   // 650 nm
    0,   // 655 nm
    0,   // 660 nm
    0,   // 665 nm
    0,   // 670 nm
    0,   // 675 nm
    0,   // 680 nm
    0,   // 685 nm
    0,   // 690 nm
    0,   // 695 nm
    0,   // 700 nm
    0,   // 705 nm
    0,   // 710 nm
    0,   // 715 nm
    0,   // 720 nm
    0,   // 725 nm
    0,   // 730 nm
    0,   // 735 nm
    0,   // 740 nm
    0,   // 745 nm
    0,   // 750 nm
    0,   // 755 nm
    0,   // 760 nm
    0,   // 765 nm
    0,   // 770 nm
    0,   // 775 nm
    0,   // 780 nm
    0,   // 785 nm
    0,   // 790 nm
    0,   // 795 nm
    0,   // 800 nm
    0,   // 805 nm
    0,   // 810 nm
    0,   // 815 nm
    0,   // 820 nm
    0,   // 825 nm
    0,   // 830 nm
};
ML_SpectralDensity sdRedChannel = {   // x,y = 0.671043872833251, 0.328718066215515
    0,   // 360 nm
    0,   // 365 nm
    0,   // 370 nm
    0,   // 375 nm
    0,   // 380 nm
    0,   // 385 nm
    0,   // 390 nm
    0,   // 395 nm
    0,   // 400 nm
    0,   // 405 nm
    0,   // 410 nm
    0,   // 415 nm
    0,   // 420 nm
    0,   // 425 nm
    0,   // 430 nm
    0,   // 435 nm
    0,   // 440 nm
    0,   // 445 nm
    0,   // 450 nm
    0,   // 455 nm
    0,   // 460 nm
    0,   // 465 nm
    0,   // 470 nm
    0,   // 475 nm
    0,   // 480 nm
    0,   // 485 nm
    0,   // 490 nm
    0,   // 495 nm
    0,   // 500 nm
    0,   // 505 nm
    0,   // 510 nm
    0,   // 515 nm
    0,   // 520 nm
    0,   // 525 nm
    0,   // 530 nm
    0,   // 535 nm
    0,   // 540 nm
    0,   // 545 nm
    0,   // 550 nm
    0,   // 555 nm
    0,   // 560 nm
    0,   // 565 nm
    0,   // 570 nm
    0,   // 575 nm
    0.0294117648154497,   // 580 nm
    0.103676468133926,   // 585 nm
    0.219411760568618,   // 590 nm
    0.407481431961059,   // 595 nm
    0.706849873065948,   // 600 nm
    0.854831933975219,   // 605 nm
    0.935294091701507,   // 610 nm
    0.972058832645416,   // 615 nm
    0.986029386520385,   // 620 nm
    0.966062068939208,   // 625 nm
    0.936764717102050,   // 630 nm
    0.882352948188781,   // 635 nm
    0.810294091701507,   // 640 nm
    0.724554359912872,   // 645 nm
    0.629954993724822,   // 650 nm
    0.538970589637756,   // 655 nm
    0.455245107412338,   // 660 nm
    0.382352948188781,   // 665 nm
    0.322058826684951,   // 670 nm
    0.276960790157318,   // 675 nm
    0.238970592617988,   // 680 nm
    0.201470583677291,   // 685 nm
    0.172058820724487,   // 690 nm
    0.144852936267852,   // 695 nm
    0.121078431606292,   // 700 nm
    0.0985294133424758,   // 705 nm
    0.0779411792755126,   // 710 nm
    0.0595588237047195,   // 715 nm
    0.0433823540806770,   // 720 nm
    0.0291666667908430,   // 725 nm
    0.0161764714866876,   // 730 nm
    0.00735294120386242,   // 735 nm
    0,   // 740 nm
    0,   // 745 nm
    0,   // 750 nm
    0,   // 755 nm
    0,   // 760 nm
    0,   // 765 nm
    0,   // 770 nm
    0,   // 775 nm
    0,   // 780 nm
    0,   // 785 nm
    0,   // 790 nm
    0,   // 795 nm
    0,   // 800 nm
    0,   // 805 nm
    0,   // 810 nm
    0,   // 815 nm
    0,   // 820 nm
    0,   // 825 nm
    0,   // 830 nm
};
ML_SpectralDensity sdClrChannel = {   // x,y = 0.358680069446563, 0.390770018100738
    0.132450327277183,   // 360 nm
    0.145695358514785,   // 365 nm
    0.157100811600685,   // 370 nm
    0.169977918267250,   // 375 nm
    0.184694632887840,   // 380 nm
    0.199779242277145,   // 385 nm
    0.216335535049438,   // 390 nm
    0.231788083910942,   // 395 nm
    0.256845206022262,   // 400 nm
    0.271891117095947,   // 405 nm
    0.295805752277374,   // 410 nm
    0.319720387458801,   // 415 nm
    0.348785877227783,   // 420 nm
    0.381898462772369,   // 425 nm
    0.401766002178192,   // 430 nm
    0.460264891386032,   // 435 nm
    0.501839637756347,   // 440 nm
    0.541574716567993,   // 445 nm
    0.576526880264282,   // 450 nm
    0.628035306930541,   // 455 nm
    0.659676253795623,   // 460 nm
    0.694260478019714,   // 465 nm
    0.728108942508697,   // 470 nm
    0.761221528053283,   // 475 nm
    0.813097894191741,   // 480 nm
    0.828918337821960,   // 485 nm
    0.857615888118743,   // 490 nm
    0.882634341716766,   // 495 nm
    0.903973519802093,   // 500 nm
    0.922737300395965,   // 505 nm
    0.934878587722778,   // 510 nm
    0.952632188796997,   // 515 nm
    0.957321584224700,   // 520 nm
    0.963576138019561,   // 525 nm
    0.969094932079315,   // 530 nm
    0.973509907722473,   // 535 nm
    0.976821184158325,   // 540 nm
    0.977924942970275,   // 545 nm
    0.979028701782226,   // 550 nm
    0.980132460594177,   // 555 nm
    0.980973422527313,   // 560 nm
    0.978660821914672,   // 565 nm
    0.974613666534423,   // 570 nm
    0.970198690891265,   // 575 nm
    0.965783655643463,   // 580 nm
    0.959161162376403,   // 585 nm
    0.952538609504699,   // 590 nm
    0.944812357425689,   // 595 nm
    0.935614466667175,   // 600 nm
    0.930463552474975,   // 605 nm
    0.916114807128906,   // 610 nm
    0.906181037425994,   // 615 nm
    0.895143508911132,   // 620 nm
    0.883002221584320,   // 625 nm
    0.869757175445556,   // 630 nm
    0.854304611682891,   // 635 nm
    0.838852107524871,   // 640 nm
    0.823399543762207,   // 645 nm
    0.796173691749572,   // 650 nm
    0.785871982574462,   // 655 nm
    0.768211901187896,   // 660 nm
    0.746136844158172,   // 665 nm
    0.721486449241638,   // 670 nm
    0.694628417491912,   // 675 nm
    0.666465997695922,   // 680 nm
    0.642752051353454,   // 685 nm
    0.615894019603729,   // 690 nm
    0.582781434059143,   // 695 nm
    0.551508486270904,   // 700 nm
    0.518763780593872,   // 705 nm
    0.473509937524795,   // 710 nm
    0.445916116237640,   // 715 nm
    0.412435621023178,   // 720 nm
    0.380794703960418,   // 725 nm
    0.350993394851684,   // 730 nm
    0.310154527425765,   // 735 nm
    0.294701993465423,   // 740 nm
    0.269315659999847,   // 745 nm
    0.242825612425804,   // 750 nm
    0.218175143003463,   // 755 nm
    0.194260492920875,   // 760 nm
    0.178594321012496,   // 765 nm
    0.151949971914291,   // 770 nm
    0.132450327277183,   // 775 nm
    0.112582780420780,   // 780 nm
    0.0971302464604377,   // 785 nm
    0.0816777050495147,   // 790 nm
    0.0676968395709991,   // 795 nm
    0.0625459924340248,   // 800 nm
    0.0452538616955280,   // 805 nm
    0.0353200882673263,   // 810 nm
    0.0264900662004947,   // 815 nm
    0.0198675505816936,   // 820 nm
    0.0132450349628925,   // 825 nm
    0.00662251934409141,   // 830 nm
};
ML_SpectralDensity sdVenus = {   // x,y = 0.340176731348037, 0.416712403297424
    0.192337170243263,   // 360 nm
    0.205605223774909,   // 365 nm
    0.205363988876342,   // 370 nm
    0.200893998146057,   // 375 nm
    0.200766280293464,   // 380 nm
    0.185823753476142,   // 385 nm
    0.179412528872489,   // 390 nm
    0.187739461660385,   // 395 nm
    0.260836243629455,   // 400 nm
    0.311749696731567,   // 405 nm
    0.333278596401214,   // 410 nm
    0.347151964902877,   // 415 nm
    0.344389706850051,   // 420 nm
    0.321455925703048,   // 425 nm
    0.287611752748489,   // 430 nm
    0.281992346048355,   // 435 nm
    0.297318011522293,   // 440 nm
    0.316985964775085,   // 445 nm
    0.328309923410415,   // 450 nm
    0.332694768905639,   // 455 nm
    0.333547890186309,   // 460 nm
    0.330217123031616,   // 465 nm
    0.329118788242340,   // 470 nm
    0.345210731029510,   // 475 nm
    0.384674340486526,   // 480 nm
    0.430906802415847,   // 485 nm
    0.563210189342498,   // 490 nm
    0.713409960269927,   // 495 nm
    0.791347980499267,   // 500 nm
    0.778780937194824,   // 505 nm
    0.767916560173034,   // 510 nm
    0.745976984500885,   // 515 nm
    0.712643682956695,   // 520 nm
    0.752873599529266,   // 525 nm
    0.791187763214111,   // 530 nm
    0.816475093364715,   // 535 nm
    0.770625829696655,   // 540 nm
    0.737164735794067,   // 545 nm
    0.708045959472656,   // 550 nm
    0.660919547080993,   // 555 nm
    0.616475105285644,   // 560 nm
    0.582120060920715,   // 565 nm
    0.561302661895751,   // 570 nm
    0.542911887168884,   // 575 nm
    0.529369115829467,   // 580 nm
    0.518177926540374,   // 585 nm
    0.503448247909545,   // 590 nm
    0.455172419548034,   // 595 nm
    0.346729099750518,   // 600 nm
    0.432694762945175,   // 605 nm
    0.488122612237930,   // 610 nm
    0.505566835403442,   // 615 nm
    0.512909293174743,   // 620 nm
    0.518915116786956,   // 625 nm
    0.521300613880157,   // 630 nm
    0.520977020263671,   // 635 nm
    0.526767671108245,   // 640 nm
    0.534982025623321,   // 645 nm
    0.542188167572021,   // 650 nm
    0.541981399059295,   // 655 nm
    0.546615600585937,   // 660 nm
    0.593486607074737,   // 665 nm
    0.628735661506652,   // 670 nm
    0.655172407627105,   // 675 nm
    0.675862073898315,   // 680 nm
    0.685808718204498,   // 685 nm
    0.689673840999603,   // 690 nm
    0.692337155342102,   // 695 nm
    0.677190065383911,   // 700 nm
    0.656704962253570,   // 705 nm
    0.630651354789733,   // 710 nm
    0.604214549064636,   // 715 nm
    0.585002720355987,   // 720 nm
    0.576158106327056,   // 725 nm
    0.569035172462463,   // 730 nm
    0.559163451194763,   // 735 nm
    0.548834443092346,   // 740 nm
    0.545407474040985,   // 745 nm
    0.540977001190185,   // 750 nm
    0.538553237915039,   // 755 nm
    0.533765614032745,   // 760 nm
    0.526436805725097,   // 765 nm
    0.518506646156311,   // 770 nm
    0.512586236000061,   // 775 nm
    0.505546450614929,   // 780 nm
    0.498488724231719,   // 785 nm
    0.495855093002319,   // 790 nm
    0.493933588266372,   // 795 nm
    0.495543450117111,   // 800 nm
    0.505605995655059,   // 805 nm
    0.510427653789520,   // 810 nm
    0.509310364723205,   // 815 nm
    0.502461433410644,   // 820 nm
    0.501232266426086,   // 825 nm
    0.494609296321868,   // 830 nm
};
ML_SpectralDensity sdHotRed = {   // x,y = 0.633847773075103, 0.307286202907562
    1.68540864251554e-4,   // 360 nm
    4.49439743533730e-4,   // 365 nm
    7.30338622815907e-4,   // 370 nm
    0.00101123750209808,   // 375 nm
    0.00129213638138025,   // 380 nm
    0.00157303526066243,   // 385 nm
    0.00185393413994461,   // 390 nm
    0.00213483301922678,   // 395 nm
    0.00241573178209364,   // 400 nm
    0.00241573178209364,   // 405 nm
    0.00297752954065799,   // 410 nm
    0.00269663077779114,   // 415 nm
    0.00325842853635549,   // 420 nm
    0.00297752954065799,   // 425 nm
    0.00325842853635549,   // 430 nm
    0.00261441641487181,   // 435 nm
    0.00382022629491984,   // 440 nm
    0.00372659345157444,   // 445 nm
    0.00410112505778670,   // 450 nm
    0.00410112505778670,   // 455 nm
    0.00438202405348420,   // 460 nm
    0.00438202405348420,   // 465 nm
    0.00367348804138600,   // 470 nm
    0.00499257259070873,   // 475 nm
    0.00494382157921791,   // 480 nm
    0.00494382157921791,   // 485 nm
    0.00466292304918169,   // 490 nm
    0.00494382157921791,   // 495 nm
    0.00438202405348420,   // 500 nm
    0.00438202405348420,   // 505 nm
    0.00409339414909482,   // 510 nm
    0.00411038566380739,   // 515 nm
    0.00353932729922235,   // 520 nm
    0.00325842853635549,   // 525 nm
    0.00269663077779114,   // 530 nm
    0.00213483301922678,   // 535 nm
    0.00213483301922678,   // 540 nm
    0.00213483301922678,   // 545 nm
    0.00241573178209364,   // 550 nm
    0.00297752954065799,   // 555 nm
    0.00382022629491984,   // 560 nm
    0.00494382157921791,   // 565 nm
    0.00634831609204411,   // 570 nm
    0.00775281060487031,   // 575 nm
    0.0102809006348252,   // 580 nm
    0.0125280916690826,   // 585 nm
    0.0161797776818275,   // 590 nm
    0.0206536967307329,   // 595 nm
    0.0276966318488121,   // 600 nm
    0.0333146080374717,   // 605 nm
    0.0442696660757064,   // 610 nm
    0.0549438223242759,   // 615 nm
    0.0678651705384254,   // 620 nm
    0.0829400792717933,   // 625 nm
    0.0982958897948265,   // 630 nm
    0.118801504373550,   // 635 nm
    0.138651683926582,   // 640 nm
    0.161123603582382,   // 645 nm
    0.185280904173851,   // 650 nm
    0.210655435919761,   // 655 nm
    0.237340837717056,   // 660 nm
    0.264213502407073,   // 665 nm
    0.288277179002761,   // 670 nm
    0.308595508337020,   // 675 nm
    0.325730353593826,   // 680 nm
    0.340524375438690,   // 685 nm
    0.352829515933990,   // 690 nm
    0.362809002399444,   // 695 nm
    0.370205998420715,   // 700 nm
    0.376292139291763,   // 705 nm
    0.380786538124084,   // 710 nm
    0.382752805948257,   // 715 nm
    0.383595526218414,   // 720 nm
    0.383595526218414,   // 725 nm
    0.383595526218414,   // 730 nm
    0.383595526218414,   // 735 nm
    0.383595526218414,   // 740 nm
    0.383595526218414,   // 745 nm
    0.383595526218414,   // 750 nm
    0.383595526218414,   // 755 nm
    0.383595526218414,   // 760 nm
    0.383595526218414,   // 765 nm
    0.383595526218414,   // 770 nm
    0.383595526218414,   // 775 nm
    0.383595526218414,   // 780 nm
    0.383595526218414,   // 785 nm
    0.383595526218414,   // 790 nm
    0.383595526218414,   // 795 nm
    0.383595526218414,   // 800 nm
    0.383595526218414,   // 805 nm
    0.383595526218414,   // 810 nm
    0.383595526218414,   // 815 nm
    0.383595526218414,   // 820 nm
    0.383595526218414,   // 825 nm
    0.383595526218414,   // 830 nm
};
ML_SpectralDensity sdHotGrn = {   // x,y = 0.412260174751281, 0.456324040889739
    0.00578651577234268,   // 360 nm
    0.00662921275943517,   // 365 nm
    0.00747190974652767,   // 370 nm
    0.00831460673362016,   // 375 nm
    0.00915730372071266,   // 380 nm
    0.0100000007078051,   // 385 nm
    0.0108426976948976,   // 390 nm
    0.0116853946819901,   // 395 nm
    0.0125280916690826,   // 400 nm
    0.0133707886561751,   // 405 nm
    0.0143071180209517,   // 410 nm
    0.0150561816990375,   // 415 nm
    0.0153370806947350,   // 420 nm
    0.0158988777548074,   // 425 nm
    0.0158988777548074,   // 430 nm
    0.0161797776818275,   // 435 nm
    0.0167415756732225,   // 440 nm
    0.0167415756732225,   // 445 nm
    0.0164606757462024,   // 450 nm
    0.0167415756732225,   // 455 nm
    0.0164606757462024,   // 460 nm
    0.0164606757462024,   // 465 nm
    0.0161797776818275,   // 470 nm
    0.0167415756732225,   // 475 nm
    0.0167415756732225,   // 480 nm
    0.0173033718019723,   // 485 nm
    0.0181460697203874,   // 490 nm
    0.0195505637675523,   // 495 nm
    0.0217041224241256,   // 500 nm
    0.0251685418188571,   // 505 nm
    0.0257814116775989,   // 510 nm
    0.0366853959858417,   // 515 nm
    0.0431460700929164,   // 520 nm
    0.0496067441999912,   // 525 nm
    0.0560674183070659,   // 530 nm
    0.0622471943497657,   // 535 nm
    0.0678651705384254,   // 540 nm
    0.0729213505983352,   // 545 nm
    0.0774157345294952,   // 550 nm
    0.0813483148813247,   // 555 nm
    0.0833146125078201,   // 560 nm
    0.0838764086365699,   // 565 nm
    0.0819101184606552,   // 570 nm
    0.0791011229157447,   // 575 nm
    0.0754494443535804,   // 580 nm
    0.0713295936584472,   // 585 nm
    0.0667415782809257,   // 590 nm
    0.0625280886888504,   // 595 nm
    0.0578464455902576,   // 600 nm
    0.0532584302127361,   // 605 nm
    0.0498876422643661,   // 610 nm
    0.0476404502987861,   // 615 nm
    0.0459550581872463,   // 620 nm
    0.0453932620584964,   // 625 nm
    0.0453932620584964,   // 630 nm
    0.0453932620584964,   // 635 nm
    0.0462359562516212,   // 640 nm
    0.0470786541700363,   // 645 nm
    0.0487640462815761,   // 650 nm
    0.0504494421184062,   // 655 nm
    0.0526029989123344,   // 660 nm
    0.0543820261955261,   // 665 nm
    0.0566292144358158,   // 670 nm
    0.0600000023841857,   // 675 nm
    0.0636516883969306,   // 680 nm
    0.0658988803625106,   // 685 nm
    0.0688951388001441,   // 690 nm
    0.0720786526799201,   // 695 nm
    0.0757303386926651,   // 700 nm
    0.0793820247054100,   // 705 nm
    0.0830337107181549,   // 710 nm
    0.0869662985205650,   // 715 nm
    0.0908988788723945,   // 720 nm
    0.0947378352284431,   // 725 nm
    0.0987640470266342,   // 730 nm
    0.102977529168128,   // 735 nm
    0.102977529168128,   // 740 nm
    0.102977529168128,   // 745 nm
    0.102977529168128,   // 750 nm
    0.102977529168128,   // 755 nm
    0.102977529168128,   // 760 nm
    0.102977529168128,   // 765 nm
    0.102977529168128,   // 770 nm
    0.102977529168128,   // 775 nm
    0.102977529168128,   // 780 nm
    0.102977529168128,   // 785 nm
    0.102977529168128,   // 790 nm
    0.102977529168128,   // 795 nm
    0.102977529168128,   // 800 nm
    0.102977529168128,   // 805 nm
    0.102977529168128,   // 810 nm
    0.102977529168128,   // 815 nm
    0.102977529168128,   // 820 nm
    0.102977529168128,   // 825 nm
    0.102977529168128,   // 830 nm
};
ML_SpectralDensity sdHotBlu = {   // x,y = 0.253991484642028, 0.275550991296768
    0,   // 360 nm
    0,   // 365 nm
    0,   // 370 nm
    0,   // 375 nm
    0.00110486708581447,   // 380 nm
    0.00831460580229759,   // 385 nm
    0.0155243445187807,   // 390 nm
    0.0227340832352638,   // 395 nm
    0.0299438219517469,   // 400 nm
    0.0371535606682300,   // 405 nm
    0.0444569326937198,   // 410 nm
    0.0511985048651695,   // 415 nm
    0.0577528104186058,   // 420 nm
    0.0639325901865959,   // 425 nm
    0.0695505663752555,   // 430 nm
    0.0751319080591201,   // 435 nm
    0.0802247226238250,   // 440 nm
    0.0839700400829315,   // 445 nm
    0.0864044949412345,   // 450 nm
    0.0875280946493148,   // 455 nm
    0.0880898907780647,   // 460 nm
    0.0880898907780647,   // 465 nm
    0.0872471928596496,   // 470 nm
    0.0864044949412345,   // 475 nm
    0.0850000008940696,   // 480 nm
    0.0833146125078201,   // 485 nm
    0.0801183506846427,   // 490 nm
    0.0788202285766601,   // 495 nm
    0.0760112404823303,   // 500 nm
    0.0737640485167503,   // 505 nm
    0.0709550604224205,   // 510 nm
    0.0675842687487602,   // 515 nm
    0.0642715245485305,   // 520 nm
    0.0608427003026008,   // 525 nm
    0.0571910142898559,   // 530 nm
    0.0552502572536468,   // 535 nm
    0.0507303401827812,   // 540 nm
    0.0479213520884513,   // 545 nm
    0.0451123602688312,   // 550 nm
    0.0428651720285415,   // 555 nm
    0.0407116152346134,   // 560 nm
    0.0380898900330066,   // 565 nm
    0.0366853959858417,   // 570 nm
    0.0347919277846813,   // 575 nm
    0.0335955061018466,   // 580 nm
    0.0324719138443470,   // 585 nm
    0.0313483178615570,   // 590 nm
    0.0304119884967803,   // 595 nm
    0.0296629238873720,   // 600 nm
    0.0293820239603519,   // 605 nm
    0.0296629238873720,   // 610 nm
    0.0287376102060079,   // 615 nm
    0.0302247218787670,   // 620 nm
    0.0305056199431419,   // 625 nm
    0.0316292159259319,   // 630 nm
    0.0327528119087219,   // 635 nm
    0.0344382040202617,   // 640 nm
    0.0364606752991676,   // 645 nm
    0.0389325879514217,   // 650 nm
    0.0417415760457515,   // 655 nm
    0.0448314622044563,   // 660 nm
    0.0493258461356163,   // 665 nm
    0.0551310889422893,   // 670 nm
    0.0600000023841857,   // 675 nm
    0.0670224726200103,   // 680 nm
    0.0742322131991386,   // 685 nm
    0.0829829499125480,   // 690 nm
    0.0875280946493148,   // 695 nm
    0.0938675403594970,   // 700 nm
    0.100636705756187,   // 705 nm
    0.106629215180873,   // 710 nm
    0.112715363502502,   // 715 nm
    0.118707865476608,   // 720 nm
    0.123483151197433,   // 725 nm
    0.127415731549263,   // 730 nm
    0.129475668072700,   // 735 nm
    0.129475668072700,   // 740 nm
    0.129475668072700,   // 745 nm
    0.129475668072700,   // 750 nm
    0.129475668072700,   // 755 nm
    0.129475668072700,   // 760 nm
    0.129475668072700,   // 765 nm
    0.129475668072700,   // 770 nm
    0.129475668072700,   // 775 nm
    0.129475668072700,   // 780 nm
    0.129475668072700,   // 785 nm
    0.129475668072700,   // 790 nm
    0.129475668072700,   // 795 nm
    0.129475668072700,   // 800 nm
    0.129475668072700,   // 805 nm
    0.129475668072700,   // 810 nm
    0.129475668072700,   // 815 nm
    0.129475668072700,   // 820 nm
    0.129475668072700,   // 825 nm
    0.129475668072700,   // 830 nm
};
ML_SpectralDensity sdHotGray = {   // x,y = 0.336299538612365, 0.347512692213058
    0,   // 360 nm
    0,   // 365 nm
    0,   // 370 nm
    0,   // 375 nm
    0,   // 380 nm
    0,   // 385 nm
    0,   // 390 nm
    0.0203932598233222,   // 395 nm
    0.0425842702388763,   // 400 nm
    0.0647752806544303,   // 405 nm
    0.0834082439541816,   // 410 nm
    0.101385772228240,   // 415 nm
    0.114400751888751,   // 420 nm
    0.131070971488952,   // 425 nm
    0.142865166068077,   // 430 nm
    0.152696639299392,   // 435 nm
    0.162059932947158,   // 440 nm
    0.170674160122871,   // 445 nm
    0.177696630358695,   // 450 nm
    0.183307766914367,   // 455 nm
    0.189426288008689,   // 460 nm
    0.193841412663459,   // 465 nm
    0.197359561920166,   // 470 nm
    0.199101135134696,   // 475 nm
    0.199606746435165,   // 480 nm
    0.199606746435165,   // 485 nm
    0.199232220649719,   // 490 nm
    0.198951318860054,   // 495 nm
    0.198483154177665,   // 500 nm
    0.198075234889984,   // 505 nm
    0.197188183665275,   // 510 nm
    0.195769131183624,   // 515 nm
    0.195112362504005,   // 520 nm
    0.194269672036170,   // 525 nm
    0.193146064877510,   // 530 nm
    0.192022472620010,   // 535 nm
    0.191179782152175,   // 540 nm
    0.190056189894676,   // 545 nm
    0.188932582736015,   // 550 nm
    0.187948644161224,   // 555 nm
    0.186953082680702,   // 560 nm
    0.185655444860458,   // 565 nm
    0.184719100594520,   // 570 nm
    0.184157311916351,   // 575 nm
    0.183595508337020,   // 580 nm
    0.183127343654632,   // 585 nm
    0.182846456766128,   // 590 nm
    0.182752817869186,   // 595 nm
    0.183033704757690,   // 600 nm
    0.183180510997772,   // 605 nm
    0.183314606547355,   // 610 nm
    0.183595508337020,   // 615 nm
    0.184157311916351,   // 620 nm
    0.185561805963516,   // 625 nm
    0.186123594641685,   // 630 nm
    0.187247201800346,   // 635 nm
    0.188089892268180,   // 640 nm
    0.189775288105010,   // 645 nm
    0.190898880362510,   // 650 nm
    0.192584276199340,   // 655 nm
    0.194331437349319,   // 660 nm
    0.195674166083335,   // 665 nm
    0.198202252388000,   // 670 nm
    0.199606746435165,   // 675 nm
    0.201573044061660,   // 680 nm
    0.203820228576660,   // 685 nm
    0.205786526203155,   // 690 nm
    0.208033710718154,   // 695 nm
    0.210280910134315,   // 700 nm
    0.212626218795776,   // 705 nm
    0.215051293373107,   // 710 nm
    0.216830283403396,   // 715 nm
    0.220112368464469,   // 720 nm
    0.222640454769134,   // 725 nm
    0.225449442863464,   // 730 nm
    0.228258430957794,   // 735 nm
    0.231067419052124,   // 740 nm
    0.231637328863143,   // 745 nm
    0.230786517262458,   // 750 nm
    0.229935705661773,   // 755 nm
    0.229084894061088,   // 760 nm
    0.228234082460403,   // 765 nm
    0.227383270859718,   // 770 nm
    0.226532459259033,   // 775 nm
    0.225681647658348,   // 780 nm
    0.224830836057662,   // 785 nm
    0.223980024456977,   // 790 nm
    0.223129212856292,   // 795 nm
    0.222278401255607,   // 800 nm
    0.221427589654922,   // 805 nm
    0.220576778054237,   // 810 nm
    0.219725966453552,   // 815 nm
    0.218875154852867,   // 820 nm
    0.218024343252182,   // 825 nm
    0.217173531651496,   // 830 nm
};
ML_SpectralDensity sdCoolRed = {   // x,y = 0.612012922763824, 0.345946282148361
    0.0211413055658340,   // 360 nm
    0.0211413055658340,   // 365 nm
    0.0211413055658340,   // 370 nm
    0.0211413055658340,   // 375 nm
    0.0208695661276578,   // 380 nm
    0.0195108707994222,   // 385 nm
    0.0186956524848937,   // 390 nm
    0.0195108707994222,   // 395 nm
    0.0178804360330104,   // 400 nm
    0.0178804360330104,   // 405 nm
    0.0178804360330104,   // 410 nm
    0.0162500012665987,   // 415 nm
    0.0162500012665987,   // 420 nm
    0.0151630453765392,   // 425 nm
    0.0154347838833928,   // 430 nm
    0.0154347838833928,   // 435 nm
    0.0146195665001869,   // 440 nm
    0.0138043491169810,   // 445 nm
    0.0129891317337751,   // 450 nm
    0.0129891317337751,   // 455 nm
    0.0129891317337751,   // 460 nm
    0.0121739143505692,   // 465 nm
    0.0121739143505692,   // 470 nm
    0.0118536502122879,   // 475 nm
    0.0113586969673633,   // 480 nm
    0.0113586969673633,   // 485 nm
    0.0113586969673633,   // 490 nm
    0.0113586969673633,   // 495 nm
    0.0105434795841574,   // 500 nm
    0.0105434795841574,   // 505 nm
    0.0105434795841574,   // 510 nm
    0.0105434795841574,   // 515 nm
    0.0113586969673633,   // 520 nm
    0.0113586969673633,   // 525 nm
    0.0113586969673633,   // 530 nm
    0.0129891317337751,   // 535 nm
    0.0140760885551571,   // 540 nm
    0.0162500012665987,   // 545 nm
    0.0213043484836816,   // 550 nm
    0.0263405814766883,   // 555 nm
    0.0339130461215972,   // 560 nm
    0.0469565242528915,   // 565 nm
    0.0648913085460662,   // 570 nm
    0.0871739163994789,   // 575 nm
    0.123858697712421,   // 580 nm
    0.164619565010070,   // 585 nm
    0.203749999403953,   // 590 nm
    0.239891305565834,   // 595 nm
    0.270597815513610,   // 600 nm
    0.302119582891464,   // 605 nm
    0.330108702182769,   // 610 nm
    0.350652188062667,   // 615 nm
    0.376576095819473,   // 620 nm
    0.390706539154052,   // 625 nm
    0.404293477535247,   // 630 nm
    0.417336970567703,   // 635 nm
    0.428750008344650,   // 640 nm
    0.438532620668411,   // 645 nm
    0.447445660829544,   // 650 nm
    0.454565227031707,   // 655 nm
    0.461358696222305,   // 660 nm
    0.467880427837371,   // 665 nm
    0.472771733999252,   // 670 nm
    0.479293465614318,   // 675 nm
    0.485639810562133,   // 680 nm
    0.490835249423980,   // 685 nm
    0.495597839355468,   // 690 nm
    0.499945640563964,   // 695 nm
    0.503750026226043,   // 700 nm
    0.508641302585601,   // 705 nm
    0.512717366218566,   // 710 nm
    0.516793489456176,   // 715 nm
    0.516793489456176,   // 720 nm
    0.516793489456176,   // 725 nm
    0.516793489456176,   // 730 nm
    0.516793489456176,   // 735 nm
    0.516793489456176,   // 740 nm
    0.516793489456176,   // 745 nm
    0.516793489456176,   // 750 nm
    0.516793489456176,   // 755 nm
    0.516793489456176,   // 760 nm
    0.516793489456176,   // 765 nm
    0.516793489456176,   // 770 nm
    0.516793489456176,   // 775 nm
    0.516793489456176,   // 780 nm
    0.516793489456176,   // 785 nm
    0.516793489456176,   // 790 nm
    0.516793489456176,   // 795 nm
    0.516793489456176,   // 800 nm
    0.516793489456176,   // 805 nm
    0.516793489456176,   // 810 nm
    0.516793489456176,   // 815 nm
    0.516793489456176,   // 820 nm
    0.516793489456176,   // 825 nm
    0.516793489456176,   // 830 nm
};
ML_SpectralDensity sdCoolGrn = {   // x,y = 0.332766711711883, 0.450448364019393
    0,   // 360 nm
    0,   // 365 nm
    0.00864130631089210,   // 370 nm
    0.0222282633185386,   // 375 nm
    0.0358152203261852,   // 380 nm
    0.0494021773338317,   // 385 nm
    0.0738587006926536,   // 390 nm
    0.103478260338306,   // 395 nm
    0.127934783697128,   // 400 nm
    0.139347821474075,   // 405 nm
    0.135543480515480,   // 410 nm
    0.120597824454307,   // 415 nm
    0.0969565212726593,   // 420 nm
    0.0811956524848937,   // 425 nm
    0.0676086992025375,   // 430 nm
    0.0578260868787765,   // 435 nm
    0.0504891313612461,   // 440 nm
    0.0464130453765392,   // 445 nm
    0.0431521758437156,   // 450 nm
    0.0415217392146587,   // 455 nm
    0.0423369593918323,   // 460 nm
    0.0431521758437156,   // 465 nm
    0.0447826087474822,   // 470 nm
    0.0480434782803058,   // 475 nm
    0.0518478304147720,   // 480 nm
    0.0578260868787765,   // 485 nm
    0.0651630461215972,   // 490 nm
    0.0744021758437156,   // 495 nm
    0.0869021788239479,   // 500 nm
    0.101032607257366,   // 505 nm
    0.123858697712421,   // 510 nm
    0.157554358243942,   // 515 nm
    0.187445655465126,   // 520 nm
    0.207010865211486,   // 525 nm
    0.215978264808654,   // 530 nm
    0.211902171373367,   // 535 nm
    0.201032608747482,   // 540 nm
    0.183097824454307,   // 545 nm
    0.162173911929130,   // 550 nm
    0.144239127635955,   // 555 nm
    0.127391308546066,   // 560 nm
    0.111902177333831,   // 565 nm
    0.0985869541764259,   // 570 nm
    0.0879891291260719,   // 575 nm
    0.0798369571566581,   // 580 nm
    0.0741304382681846,   // 585 nm
    0.0708695650100708,   // 590 nm
    0.0692391321063041,   // 595 nm
    0.0692391321063041,   // 600 nm
    0.0700543522834777,   // 605 nm
    0.0719565227627754,   // 610 nm
    0.0741304382681846,   // 615 nm
    0.0790217369794845,   // 620 nm
    0.0833695679903030,   // 625 nm
    0.0888043493032455,   // 630 nm
    0.0945108681917190,   // 635 nm
    0.100217394530773,   // 640 nm
    0.107826091349124,   // 645 nm
    0.115434788167476,   // 650 nm
    0.123043477535247,   // 655 nm
    0.132010877132415,   // 660 nm
    0.140978261828422,   // 665 nm
    0.149945646524429,   // 670 nm
    0.160543486475944,   // 675 nm
    0.171684786677360,   // 680 nm
    0.180923908948898,   // 685 nm
    0.192336961627006,   // 690 nm
    0.203749999403953,   // 695 nm
    0.215163052082061,   // 700 nm
    0.227391302585601,   // 705 nm
    0.238804355263710,   // 710 nm
    0.252119570970535,   // 715 nm
    0.252119570970535,   // 720 nm
    0.252119570970535,   // 725 nm
    0.252119570970535,   // 730 nm
    0.252119570970535,   // 735 nm
    0.252119570970535,   // 740 nm
    0.252119570970535,   // 745 nm
    0.252119570970535,   // 750 nm
    0.252119570970535,   // 755 nm
    0.252119570970535,   // 760 nm
    0.252119570970535,   // 765 nm
    0.252119570970535,   // 770 nm
    0.252119570970535,   // 775 nm
    0.252119570970535,   // 780 nm
    0.252119570970535,   // 785 nm
    0.252119570970535,   // 790 nm
    0.252119570970535,   // 795 nm
    0.252119570970535,   // 800 nm
    0.252119570970535,   // 805 nm
    0.252119570970535,   // 810 nm
    0.252119570970535,   // 815 nm
    0.252119570970535,   // 820 nm
    0.252119570970535,   // 825 nm
    0.252119570970535,   // 830 nm
};
ML_SpectralDensity sdCoolBlu = {   // x,y = 0.239894628524780, 0.240181058645248
    0.0972282588481903,   // 360 nm
    0.101032607257366,   // 365 nm
    0.104836955666542,   // 370 nm
    0.108641304075717,   // 375 nm
    0.112445652484893,   // 380 nm
    0.116250000894069,   // 385 nm
    0.119375005364418,   // 390 nm
    0.121413044631481,   // 395 nm
    0.125489130616188,   // 400 nm
    0.124673917889595,   // 405 nm
    0.124673917889595,   // 410 nm
    0.124673917889595,   // 415 nm
    0.127119570970535,   // 420 nm
    0.127934783697128,   // 425 nm
    0.128583803772926,   // 430 nm
    0.129565224051475,   // 435 nm
    0.129565224051475,   // 440 nm
    0.129565224051475,   // 445 nm
    0.129565224051475,   // 450 nm
    0.129565224051475,   // 455 nm
    0.129565224051475,   // 460 nm
    0.130380436778068,   // 465 nm
    0.129565224051475,   // 470 nm
    0.129565224051475,   // 475 nm
    0.129565224051475,   // 480 nm
    0.131195649504661,   // 485 nm
    0.131195649504661,   // 490 nm
    0.130819395184516,   // 495 nm
    0.127119570970535,   // 500 nm
    0.119782611727714,   // 505 nm
    0.110815219581127,   // 510 nm
    0.0977717414498329,   // 515 nm
    0.0874456539750099,   // 520 nm
    0.0764855146408081,   // 525 nm
    0.0673369616270065,   // 530 nm
    0.0611224025487899,   // 535 nm
    0.0556521788239479,   // 540 nm
    0.0526630468666553,   // 545 nm
    0.0487549416720867,   // 550 nm
    0.0472282618284225,   // 555 nm
    0.0447826087474822,   // 560 nm
    0.0415217392146587,   // 565 nm
    0.0398913063108921,   // 570 nm
    0.0374456532299518,   // 575 nm
    0.0355341099202632,   // 580 nm
    0.0350000001490116,   // 585 nm
    0.0343206524848937,   // 590 nm
    0.0350000001490116,   // 595 nm
    0.0358152203261852,   // 600 nm
    0.0363586992025375,   // 605 nm
    0.0366304367780685,   // 610 nm
    0.0382608696818351,   // 615 nm
    0.0382608696818351,   // 620 nm
    0.0398913063108921,   // 625 nm
    0.0415217392146587,   // 630 nm
    0.0447826087474822,   // 635 nm
    0.0521195679903030,   // 640 nm
    0.0635326132178306,   // 645 nm
    0.0773913040757179,   // 650 nm
    0.0896195694804191,   // 655 nm
    0.101736664772033,   // 660 nm
    0.113260872662067,   // 665 nm
    0.125038936734199,   // 670 nm
    0.133641317486763,   // 675 nm
    0.140335619449615,   // 680 nm
    0.146750256419181,   // 685 nm
    0.148858696222305,   // 690 nm
    0.150371521711349,   // 695 nm
    0.150760874152183,   // 700 nm
    0.150760874152183,   // 705 nm
    0.150760874152183,   // 710 nm
    0.149130433797836,   // 715 nm
    0.147499993443489,   // 720 nm
    0.145869553089141,   // 725 nm
    0.144239112734794,   // 730 nm
    0.142608672380447,   // 735 nm
    0.140978232026100,   // 740 nm
    0.139347791671752,   // 745 nm
    0.137717351317405,   // 750 nm
    0.136086910963058,   // 755 nm
    0.134456470608711,   // 760 nm
    0.132826030254364,   // 765 nm
    0.131195589900016,   // 770 nm
    0.129565149545669,   // 775 nm
    0.127934709191322,   // 780 nm
    0.126304268836975,   // 785 nm
    0.124673828482627,   // 790 nm
    0.123043388128280,   // 795 nm
    0.121412947773933,   // 800 nm
    0.119782507419586,   // 805 nm
    0.118152067065238,   // 810 nm
    0.116521626710891,   // 815 nm
    0.114891186356544,   // 820 nm
    0.113260746002197,   // 825 nm
    0.111630305647850,   // 830 nm
};
ML_SpectralDensity sdCoolGray = {   // x,y = 0.330785274505615, 0.332112520933151
    0,   // 360 nm
    0,   // 365 nm
    0,   // 370 nm
    0.0132608711719512,   // 375 nm
    0.0458695665001869,   // 380 nm
    0.0784782618284225,   // 385 nm
    0.111630439758300,   // 390 nm
    0.139347821474075,   // 395 nm
    0.162173911929130,   // 400 nm
    0.179293483495712,   // 405 nm
    0.190434783697128,   // 410 nm
    0.196141302585601,   // 415 nm
    0.198858693242073,   // 420 nm
    0.202119573950767,   // 425 nm
    0.203749999403953,   // 430 nm
    0.205380439758300,   // 435 nm
    0.207010865211486,   // 440 nm
    0.208641305565834,   // 445 nm
    0.210271745920181,   // 450 nm
    0.211243063211441,   // 455 nm
    0.212647527456283,   // 460 nm
    0.213295459747314,   // 465 nm
    0.214347824454307,   // 470 nm
    0.215163052082061,   // 475 nm
    0.216793477535247,   // 480 nm
    0.216793477535247,   // 485 nm
    0.216793477535247,   // 490 nm
    0.216793477535247,   // 495 nm
    0.215978264808654,   // 500 nm
    0.216793477535247,   // 505 nm
    0.215148225426673,   // 510 nm
    0.213532611727714,   // 515 nm
    0.211902171373367,   // 520 nm
    0.210271745920181,   // 525 nm
    0.211086958646774,   // 530 nm
    0.209456518292427,   // 535 nm
    0.206467390060424,   // 540 nm
    0.206195652484893,   // 545 nm
    0.203913047909736,   // 550 nm
    0.202934786677360,   // 555 nm
    0.202317193150520,   // 560 nm
    0.200760871171951,   // 565 nm
    0.200489133596420,   // 570 nm
    0.200489133596420,   // 575 nm
    0.198858693242073,   // 580 nm
    0.198858693242073,   // 585 nm
    0.199673920869827,   // 590 nm
    0.198586955666542,   // 595 nm
    0.198722824454307,   // 600 nm
    0.199673920869827,   // 605 nm
    0.199673920869827,   // 610 nm
    0.200489133596420,   // 615 nm
    0.201576098799705,   // 620 nm
    0.203749999403953,   // 625 nm
    0.206195652484893,   // 630 nm
    0.210271745920181,   // 635 nm
    0.216793477535247,   // 640 nm
    0.222116380929946,   // 645 nm
    0.227934792637825,   // 650 nm
    0.232282608747482,   // 655 nm
    0.236358702182769,   // 660 nm
    0.240434780716896,   // 665 nm
    0.244510874152183,   // 670 nm
    0.247771739959716,   // 675 nm
    0.252391308546066,   // 680 nm
    0.255499988794326,   // 685 nm
    0.259999990463256,   // 690 nm
    0.263285577297210,   // 695 nm
    0.267336964607238,   // 700 nm
    0.272228270769119,   // 705 nm
    0.277119576930999,   // 710 nm
    0.281195640563964,   // 715 nm
    0.285271733999252,   // 720 nm
    0.285271733999252,   // 725 nm
    0.285271733999252,   // 730 nm
    0.285271733999252,   // 735 nm
    0.285271733999252,   // 740 nm
    0.285271733999252,   // 745 nm
    0.285271733999252,   // 750 nm
    0.285271733999252,   // 755 nm
    0.285271733999252,   // 760 nm
    0.285271733999252,   // 765 nm
    0.285271733999252,   // 770 nm
    0.285271733999252,   // 775 nm
    0.285271733999252,   // 780 nm
    0.285271733999252,   // 785 nm
    0.285271733999252,   // 790 nm
    0.285271733999252,   // 795 nm
    0.285271733999252,   // 800 nm
    0.285271733999252,   // 805 nm
    0.285271733999252,   // 810 nm
    0.285271733999252,   // 815 nm
    0.285271733999252,   // 820 nm
    0.285271733999252,   // 825 nm
    0.285271733999252,   // 830 nm
};



void
TranslateSpectra()
{
    Spectrum s;

    s = MLtoSpectrum(sdCoolGray);
    PrintSpectrum("g_sCoolGray", s);
    s = MLtoSpectrum(sdCoolRed);
    PrintSpectrum("g_sCoolRed", s);
    s = MLtoSpectrum(sdCoolGrn);
    PrintSpectrum("g_sCoolGrn", s);
    s = MLtoSpectrum(sdCoolBlu);
    PrintSpectrum("g_sCoolBlu", s);
    s = MLtoSpectrum(sdHotGray);
    PrintSpectrum("g_sHotGray", s);
    s = MLtoSpectrum(sdHotRed);
    PrintSpectrum("g_sHotRed", s);
    s = MLtoSpectrum(sdHotGrn);
    PrintSpectrum("g_sHotGrn", s);
    s = MLtoSpectrum(sdHotBlu);
    PrintSpectrum("g_sHotBlu", s);
    /*
    s = MLtoSpectrum(sdRedChannel);
    PrintSpectrum("g_sRedChannel", s);
    s = MLtoSpectrum(sdGrnChannel);
    PrintSpectrum("g_sGrnChannel", s);
    s = MLtoSpectrum(sdBluChannel);
    PrintSpectrum("g_sBluChannel", s);
    s = MLtoSpectrum(sdClrChannel);
    PrintSpectrum("g_sClrChannel", s);
    */
    /*
    s = MLtoSpectrum(ML_Illuminant_D65);
    PrintSpectrum("g_sIllumD65", s);
    s = MLtoSpectrum(ML_SolarSpectrum);
    PrintSpectrum("g_sSolar", s);
    s = MLtoSpectrum(ML_BlueSky);
    PrintSpectrum("g_sEarthSky", s);
    s = MLtoSpectrum(ML_AlbedoMars);
    PrintSpectrum("g_sMarsAlbedo", s);
    s = MLtoSpectrum(ML_AlbedoVenus);
    PrintSpectrum("g_sVenusAlbedo", s);
    s = MLtoSpectrum(ML_VenusUVcontrast);
    PrintSpectrum("g_sVenusUVcontrast", s);
    s = MLtoSpectrum(ML_AlbedoMercury);
    PrintSpectrum("g_sMercuryAlbedo", s);
    */
    /*
    s = MLtoSpectrum(ML_Venera14_000);
    PrintSpectrum("g_sVenera14Sky", s);
    s = MLtoSpectrum(ML_Venera14_ground);
    PrintSpectrum("g_sVenera14Ground", s);
    s = MLtoSpectrum(ML_Venera13_000);
    PrintSpectrum("g_sVenera13Sky", s);
    s = MLtoSpectrum(ML_Venera11_000);
    PrintSpectrum("g_sVenera11Sky", s);
    s = MLtoSpectrum(ML_Venera11_072);
    PrintSpectrum("g_sVenera11_07200m", s);
    s = MLtoSpectrum(ML_Venera11_162);
    PrintSpectrum("g_sVenera11_16200m", s);
    s = MLtoSpectrum(ML_Venera11_237);
    PrintSpectrum("g_sVenera11_23700m", s);
    s = MLtoSpectrum(ML_Venera11_377);
    PrintSpectrum("g_sVenera11_37700m", s);
    s = MLtoSpectrum(ML_Venera11_486);
    PrintSpectrum("g_sVenera11_48600m", s);
    s = MLtoSpectrum(ML_Venera11_511);
    PrintSpectrum("g_sVenera11_51100m", s);
    s = MLtoSpectrum(ML_Venera11_561);
    PrintSpectrum("g_sVenera11_56100m", s);
    s = MLtoSpectrum(ML_Venera11_621);
    PrintSpectrum("g_sVenera11_62100m", s);
    s = MLtoSpectrum(ML_Venera11_Sun);
    PrintSpectrum("g_sVenera11_Sun", s);
    */
}

void
TestVRGB()
{
    DisplayRGB rgb;

    printf("sRGB and vRGB:\n");
    rgb = g_sEqualWhite.sRGB();
    rgb.Print();
    rgb = g_sEqualWhite.vRGB();
    rgb.Print();
}

//
//  3x3x3 cube  - 24.76 sec
//  if (ijk)    - 25.10 sec
//  iStart/iEnd - 24.5
//  dynamic S/E - 39.10
//  i != nx...  - 21.75
//
//  fmin = 0.0, fmax = 1.7958 (two hours) 
//      fmax can be as big as 3.0, but probability greater than 1.7958 is less than 1 in a billion
//
#define SCALE 32.0f

static void
IntegrateCubeSection()
{
    int i, nInside;
    float x, y, z, f;
    double p;

    nInside = 0;
    for (i = 0; i < 1000000000; i++) {
        x = RandomFloat();
        y = RandomFloat();
        z = RandomFloat();
        f = x*x + y*y + z* z;
        if (f < 1.7958)
            nInside++;
    }
    p = double(nInside)/double(i);
    printf("%d/%d = %f\n", nInside, i, p);
    p = 1.0 - p;
    p = p*p*p*p * p*p*p*p;
    printf("p = %g\n", p);
}

void
TestVoronoi()
{
    Image im;
    Video vi;
    ML_TimingInfo ti;
    int i, j, k, nFrames;
    float f, rgf[3], fmin, fmax;

    IntegrateCubeSection();
    return;
    printf("Test Voronoi\n");
    im.NewImage(1920, 1080, 3);
    vi.NewVideo("Voronoi.mpg", im);
    nFrames = 600;
    fmin = 1000.0;
    fmax = -1000.0;
    ML_StartTiming(ti);
    for (i = 0; i < nFrames; i += 1) {
        rgf[0] = float(i)/SCALE;
        for (j = 0; j < im.m_nWidth; j++) {
            rgf[1] = float(j)/SCALE;
            for (k = 0; k < im.m_nHeight; k++) {
                rgf[2] = float(k)/SCALE;
                f =  BandNoise(rgf, 3, 2);
                // f = Voronoi(rgf);    // SCALE = 64
                im.SetRGB(f*0.5f, j, k);
                if (f < fmin) fmin = f;
                if (f > fmax) fmax = f;
            }
        }
        vi.WriteFrame(im);
        break;
    }
    ML_StopTiming(ti);
    ML_ReportTiming(ti);
    printf("fmin = %f, fmax = %f, maxdist = %f\n", fmin, fmax, sqrt(fmax));
    im.WriteBMP("Veronoi.bmp");
    vi.Close();
}
//
//  Generalized golden ratio, roots of x**(d+1) = x + 1
//
//  d = 2, x = (1 + sqrt(5))/2
//  d = 3, 
//
//  g = 1.22074408460575947536
//  a1 = 1.0/g
//  a2 = 1.0/(g*g)
//  a3 = 1.0/(g*g*g)
//  x[n] = (0.5+a1*n) %1
//  y[n] = (0.5+a2*n) %1
//  z[n] = (0.5+a3*n) %1
//
extern double ML_RadicalInverse(unsigned n, unsigned nPrime);
extern double ML_FoldedRadicalInverse(unsigned n, unsigned nPrime);
static Vector2 s_rgv[1024];

void
TestQuadirandomPatterns()
{
    int i;
    double a1, a2;
    Vector2 v;

    //
    //  Hamersley points
    //
    for (i = 0; i < 1024; i++) {
        v.x = (double(i) + 0.5)/1024.0;
        v.y = ML_RadicalInverse(i, 2);
        s_rgv[i] = v;
    }
    OptimizeLocality(s_rgv, 1024);
    PatternTestImage("QuasiHammersley.bmp", s_rgv, 1024);
    for (i = 0; i < 1024; i++) {
        v.x = (double(i) + 0.5)/1024.0;
        v.y = ML_FoldedRadicalInverse(i, 2);
        s_rgv[i] = v;
    }
    OptimizeLocality(s_rgv, 1024);
    PatternTestImage("QuasiFoldedHammersley.bmp", s_rgv, 1024);
    //
    //  golden-ratio patterns
    //
    a1 = 1.0/PHI1;
    for (i = 0; i < 1024; i++) {
        v.x = (double(i) + 0.5)/1024.0;
        v.y = fmod(0.5 + a1*double(i+1), 1.0);
        s_rgv[i] = v;
    }
    OptimizeLocality(s_rgv, 1024);
    PatternTestImage("Golden1.bmp", s_rgv, 1024);
    a1 = 1.0/PHI2;
    a2 = 1.0/(PHI2*PHI2);
    for (i = 0; i < 1024; i++) {
        v.x = fmod(0.5 + a1*double(i+1), 1.0);
        v.y = fmod(0.5 + a2*double(i+1), 1.0);
        s_rgv[i] = v;
    }
    OptimizeLocality(s_rgv, 1024);
    PatternTestImage("Golden2.bmp", s_rgv, 1024);
}

extern double CubeRoot(double x);

static double
Golden2()
{
    return 0.5*(1.0 + sqrt(5.0));
}

static double
Golden3()
{
    return CubeRoot(27.0/2.0 - 3.0*sqrt(69.0)/2.0)/3.0 + CubeRoot((9.0 + sqrt(69.0))/2.0)/pow(3, 2.0/3.0);
}

static void
MakePhi(double d)
{
    double x;
    int i;

    x = 2.0;
    for (i = 0; i < 40; i++) {
        x = pow(x + 1.0, 1.0/(d+1.0));
    }
    printf("#define PHI%d  %0.15f\n", int(d), x);
}

void
GeneratePhiValues()
{
    printf("#define PHI  %0.15f\n\n", (sqrt(5.0) + 1.0)/2.0);
    MakePhi(1.0);
    MakePhi(2.0);
    MakePhi(3.0);
    MakePhi(4.0);
    MakePhi(5.0);
    printf("\ng2 = %0.15f\ng3 = %0.15f\n", Golden2(), Golden3());
}
//
//  Make X, Y, Z functions
//
struct ML_VisualXYZ {
    double X, Y, Z;

    ML_VisualXYZ(double x, double y, double z) : X(x), Y(y), Z(z) {}
};

const ML_VisualXYZ ML_CIE[471] = {
    ML_VisualXYZ(0.000129900000,  0.000003917000,  0.000606100000),   // 360 nm
    ML_VisualXYZ(0.000145847000,  0.000004393581,  0.000680879200),
    ML_VisualXYZ(0.000163802100,  0.000004929604,  0.000765145600),
    ML_VisualXYZ(0.000184003700,  0.000005532136,  0.000860012400),
    ML_VisualXYZ(0.000206690200,  0.000006208245,  0.000966592800),
    ML_VisualXYZ(0.000232100000,  0.000006965000,  0.001086000000),
    ML_VisualXYZ(0.000260728000,  0.000007813219,  0.001220586000),
    ML_VisualXYZ(0.000293075000,  0.000008767336,  0.001372729000),
    ML_VisualXYZ(0.000329388000,  0.000009839844,  0.001543579000),
    ML_VisualXYZ(0.000369914000,  0.000011043230,  0.001734286000),
    ML_VisualXYZ(0.000414900000,  0.000012390000,  0.001946000000),   // 370 nm
    ML_VisualXYZ(0.000464158700,  0.000013886410,  0.002177777000),
    ML_VisualXYZ(0.000518986000,  0.000015557280,  0.002435809000),
    ML_VisualXYZ(0.000581854000,  0.000017442960,  0.002731953000),
    ML_VisualXYZ(0.000655234700,  0.000019583750,  0.003078064000),
    ML_VisualXYZ(0.000741600000,  0.000022020000,  0.003486000000),
    ML_VisualXYZ(0.000845029600,  0.000024839650,  0.003975227000),
    ML_VisualXYZ(0.000964526800,  0.000028041260,  0.004540880000),
    ML_VisualXYZ(0.001094949000,  0.000031531040,  0.005158320000),
    ML_VisualXYZ(0.001231154000,  0.000035215210,  0.005802907000),
    ML_VisualXYZ(0.001368000000,  0.000039000000,  0.006450001000),   // 380 nm
    ML_VisualXYZ(0.001502050000,  0.000042826400,  0.007083216000),
    ML_VisualXYZ(0.001642328000,  0.000046914600,  0.007745488000),
    ML_VisualXYZ(0.001802382000,  0.000051589600,  0.008501152000),
    ML_VisualXYZ(0.001995757000,  0.000057176400,  0.009414544000),
    ML_VisualXYZ(0.002236000000,  0.000064000000,  0.010549990000),
    ML_VisualXYZ(0.002535385000,  0.000072344210,  0.011965800000),
    ML_VisualXYZ(0.002892603000,  0.000082212240,  0.013655870000),
    ML_VisualXYZ(0.003300829000,  0.000093508160,  0.015588050000),
    ML_VisualXYZ(0.003753236000,  0.000106136100,  0.017730150000),
    ML_VisualXYZ(0.004243000000,  0.000120000000,  0.020050010000),   // 390 nm
    ML_VisualXYZ(0.004762389000,  0.000134984000,  0.022511360000),
    ML_VisualXYZ(0.005330048000,  0.000151492000,  0.025202880000),
    ML_VisualXYZ(0.005978712000,  0.000170208000,  0.028279720000),
    ML_VisualXYZ(0.006741117000,  0.000191816000,  0.031897040000),
    ML_VisualXYZ(0.007650000000,  0.000217000000,  0.036210000000),
    ML_VisualXYZ(0.008751373000,  0.000246906700,  0.041437710000),
    ML_VisualXYZ(0.010028880000,  0.000281240000,  0.047503720000),
    ML_VisualXYZ(0.011421700000,  0.000318520000,  0.054119880000),
    ML_VisualXYZ(0.012869010000,  0.000357266700,  0.060998030000),
    ML_VisualXYZ(0.014310000000,  0.000396000000,  0.067850010000),   // 400 nm
    ML_VisualXYZ(0.015704430000,  0.000433714700,  0.074486320000),
    ML_VisualXYZ(0.017147440000,  0.000473024000,  0.081361560000),
    ML_VisualXYZ(0.018781220000,  0.000517876000,  0.089153640000),
    ML_VisualXYZ(0.020748010000,  0.000572218700,  0.098540480000),
    ML_VisualXYZ(0.023190000000,  0.000640000000,  0.110200000000),
    ML_VisualXYZ(0.026207360000,  0.000724560000,  0.124613300000),
    ML_VisualXYZ(0.029782480000,  0.000825500000,  0.141701700000),
    ML_VisualXYZ(0.033880920000,  0.000941160000,  0.161303500000),
    ML_VisualXYZ(0.038468240000,  0.001069880000,  0.183256800000),
    ML_VisualXYZ(0.043510000000,  0.001210000000,  0.207400000000),   // 410 nm
    ML_VisualXYZ(0.048995600000,  0.001362091000,  0.233692100000),
    ML_VisualXYZ(0.055022600000,  0.001530752000,  0.262611400000),
    ML_VisualXYZ(0.061718800000,  0.001720368000,  0.294774600000),
    ML_VisualXYZ(0.069212000000,  0.001935323000,  0.330798500000),
    ML_VisualXYZ(0.077630000000,  0.002180000000,  0.371300000000),
    ML_VisualXYZ(0.086958110000,  0.002454800000,  0.416209100000),
    ML_VisualXYZ(0.097176720000,  0.002764000000,  0.465464200000),
    ML_VisualXYZ(0.108406300000,  0.003117800000,  0.519694800000),
    ML_VisualXYZ(0.120767200000,  0.003526400000,  0.579530300000),
    ML_VisualXYZ(0.134380000000,  0.004000000000,  0.645600000000),   // 420 nm
    ML_VisualXYZ(0.149358200000,  0.004546240000,  0.718483800000),
    ML_VisualXYZ(0.165395700000,  0.005159320000,  0.796713300000),
    ML_VisualXYZ(0.181983100000,  0.005829280000,  0.877845900000),
    ML_VisualXYZ(0.198611000000,  0.006546160000,  0.959439000000),
    ML_VisualXYZ(0.214770000000,  0.007300000000,  1.039050100000),
    ML_VisualXYZ(0.230186800000,  0.008086507000,  1.115367300000),
    ML_VisualXYZ(0.244879700000,  0.008908720000,  1.188497100000),
    ML_VisualXYZ(0.258777300000,  0.009767680000,  1.258123300000),
    ML_VisualXYZ(0.271807900000,  0.010664430000,  1.323929600000),
    ML_VisualXYZ(0.283900000000,  0.011600000000,  1.385600000000),   // 430 nm
    ML_VisualXYZ(0.294943800000,  0.012573170000,  1.442635200000),
    ML_VisualXYZ(0.304896500000,  0.013582720000,  1.494803500000),
    ML_VisualXYZ(0.313787300000,  0.014629680000,  1.542190300000),
    ML_VisualXYZ(0.321645400000,  0.015715090000,  1.584880700000),
    ML_VisualXYZ(0.328500000000,  0.016840000000,  1.622960000000),
    ML_VisualXYZ(0.334351300000,  0.018007360000,  1.656404800000),
    ML_VisualXYZ(0.339210100000,  0.019214480000,  1.685295900000),
    ML_VisualXYZ(0.343121300000,  0.020453920000,  1.709874500000),
    ML_VisualXYZ(0.346129600000,  0.021718240000,  1.730382100000),
    ML_VisualXYZ(0.348280000000,  0.023000000000,  1.747060000000),   // 440 nm
    ML_VisualXYZ(0.349599900000,  0.024294610000,  1.760044600000),
    ML_VisualXYZ(0.350147400000,  0.025610240000,  1.769623300000),
    ML_VisualXYZ(0.350013000000,  0.026958570000,  1.776263700000),
    ML_VisualXYZ(0.349287000000,  0.028351250000,  1.780433400000),
    ML_VisualXYZ(0.348060000000,  0.029800000000,  1.782600000000),
    ML_VisualXYZ(0.346373300000,  0.031310830000,  1.782968200000),
    ML_VisualXYZ(0.344262400000,  0.032883680000,  1.781699800000),
    ML_VisualXYZ(0.341808800000,  0.034521120000,  1.779198200000),
    ML_VisualXYZ(0.339094100000,  0.036225710000,  1.775867100000),
    ML_VisualXYZ(0.336200000000,  0.038000000000,  1.772110000000),   // 450 nm
    ML_VisualXYZ(0.333197700000,  0.039846670000,  1.768258900000),
    ML_VisualXYZ(0.330041100000,  0.041768000000,  1.764039000000),
    ML_VisualXYZ(0.326635700000,  0.043766000000,  1.758943800000),
    ML_VisualXYZ(0.322886800000,  0.045842670000,  1.752466300000),
    ML_VisualXYZ(0.318700000000,  0.048000000000,  1.744100000000),
    ML_VisualXYZ(0.314025100000,  0.050243680000,  1.733559500000),
    ML_VisualXYZ(0.308884000000,  0.052573040000,  1.720858100000),
    ML_VisualXYZ(0.303290400000,  0.054980560000,  1.705936900000),
    ML_VisualXYZ(0.297257900000,  0.057458720000,  1.688737200000),
    ML_VisualXYZ(0.290800000000,  0.060000000000,  1.669200000000),   // 460 nm
    ML_VisualXYZ(0.283970100000,  0.062601970000,  1.647528700000),
    ML_VisualXYZ(0.276721400000,  0.065277520000,  1.623412700000),
    ML_VisualXYZ(0.268917800000,  0.068042080000,  1.596022300000),
    ML_VisualXYZ(0.260422700000,  0.070911090000,  1.564528000000),
    ML_VisualXYZ(0.251100000000,  0.073900000000,  1.528100000000),
    ML_VisualXYZ(0.240847500000,  0.077016000000,  1.486111400000),
    ML_VisualXYZ(0.229851200000,  0.080266400000,  1.439521500000),
    ML_VisualXYZ(0.218407200000,  0.083666800000,  1.389879900000),
    ML_VisualXYZ(0.206811500000,  0.087232800000,  1.338736200000),
    ML_VisualXYZ(0.195360000000,  0.090980000000,  1.287640000000),   // 470 nm
    ML_VisualXYZ(0.184213600000,  0.094917550000,  1.237422300000),
    ML_VisualXYZ(0.173327300000,  0.099045840000,  1.187824300000),
    ML_VisualXYZ(0.162688100000,  0.103367400000,  1.138761100000),
    ML_VisualXYZ(0.152283300000,  0.107884600000,  1.090148000000),
    ML_VisualXYZ(0.142100000000,  0.112600000000,  1.041900000000),
    ML_VisualXYZ(0.132178600000,  0.117532000000,  0.994197600000),
    ML_VisualXYZ(0.122569600000,  0.122674400000,  0.947347300000),
    ML_VisualXYZ(0.113275200000,  0.127992800000,  0.901453100000),
    ML_VisualXYZ(0.104297900000,  0.133452800000,  0.856619300000),
    ML_VisualXYZ(0.095640000000,  0.139020000000,  0.812950100000),   // 480 nm
    ML_VisualXYZ(0.087299550000,  0.144676400000,  0.770517300000),
    ML_VisualXYZ(0.079308040000,  0.150469300000,  0.729444800000),
    ML_VisualXYZ(0.071717760000,  0.156461900000,  0.689913600000),
    ML_VisualXYZ(0.064580990000,  0.162717700000,  0.652104900000),
    ML_VisualXYZ(0.057950010000,  0.169300000000,  0.616200000000),
    ML_VisualXYZ(0.051862110000,  0.176243100000,  0.582328600000),
    ML_VisualXYZ(0.046281520000,  0.183558100000,  0.550416200000),
    ML_VisualXYZ(0.041150880000,  0.191273500000,  0.520337600000),
    ML_VisualXYZ(0.036412830000,  0.199418000000,  0.491967300000),
    ML_VisualXYZ(0.032010000000,  0.208020000000,  0.465180000000),   // 490 nm
    ML_VisualXYZ(0.027917200000,  0.217119900000,  0.439924600000),
    ML_VisualXYZ(0.024144400000,  0.226734500000,  0.416183600000),
    ML_VisualXYZ(0.020687000000,  0.236857100000,  0.393882200000),
    ML_VisualXYZ(0.017540400000,  0.247481200000,  0.372945900000),
    ML_VisualXYZ(0.014700000000,  0.258600000000,  0.353300000000),
    ML_VisualXYZ(0.012161790000,  0.270184900000,  0.334857800000),
    ML_VisualXYZ(0.009919960000,  0.282293900000,  0.317552100000),
    ML_VisualXYZ(0.007967240000,  0.295050500000,  0.301337500000),
    ML_VisualXYZ(0.006296346000,  0.308578000000,  0.286168600000),
    ML_VisualXYZ(0.004900000000,  0.323000000000,  0.272000000000),   // 500 nm
    ML_VisualXYZ(0.003777173000,  0.338402100000,  0.258817100000),
    ML_VisualXYZ(0.002945320000,  0.354685800000,  0.246483800000),
    ML_VisualXYZ(0.002424880000,  0.371698600000,  0.234771800000),
    ML_VisualXYZ(0.002236293000,  0.389287500000,  0.223453300000),
    ML_VisualXYZ(0.002400000000,  0.407300000000,  0.212300000000),
    ML_VisualXYZ(0.002925520000,  0.425629900000,  0.201169200000),
    ML_VisualXYZ(0.003836560000,  0.444309600000,  0.190119600000),
    ML_VisualXYZ(0.005174840000,  0.463394400000,  0.179225400000),
    ML_VisualXYZ(0.006982080000,  0.482939500000,  0.168560800000),
    ML_VisualXYZ(0.009300000000,  0.503000000000,  0.158200000000),   // 510 nm
    ML_VisualXYZ(0.012149490000,  0.523569300000,  0.148138300000),
    ML_VisualXYZ(0.015535880000,  0.544512000000,  0.138375800000),
    ML_VisualXYZ(0.019477520000,  0.565690000000,  0.128994200000),
    ML_VisualXYZ(0.023992770000,  0.586965300000,  0.120075100000),
    ML_VisualXYZ(0.029100000000,  0.608200000000,  0.111700000000),
    ML_VisualXYZ(0.034814850000,  0.629345600000,  0.103904800000),
    ML_VisualXYZ(0.041120160000,  0.650306800000,  0.096667480000),
    ML_VisualXYZ(0.047985040000,  0.670875200000,  0.089982720000),
    ML_VisualXYZ(0.055378610000,  0.690842400000,  0.083845310000),
    ML_VisualXYZ(0.063270000000,  0.710000000000,  0.078249990000),   // 520 nm
    ML_VisualXYZ(0.071635010000,  0.728185200000,  0.073208990000),
    ML_VisualXYZ(0.080462240000,  0.745463600000,  0.068678160000),
    ML_VisualXYZ(0.089739960000,  0.761969400000,  0.064567840000),
    ML_VisualXYZ(0.099456450000,  0.777836800000,  0.060788350000),
    ML_VisualXYZ(0.109600000000,  0.793200000000,  0.057250010000),
    ML_VisualXYZ(0.120167400000,  0.808110400000,  0.053904350000),
    ML_VisualXYZ(0.131114500000,  0.822496200000,  0.050746640000),
    ML_VisualXYZ(0.142367900000,  0.836306800000,  0.047752760000),
    ML_VisualXYZ(0.153854200000,  0.849491600000,  0.044898590000),
    ML_VisualXYZ(0.165500000000,  0.862000000000,  0.042160000000),   // 530 nm
    ML_VisualXYZ(0.177257100000,  0.873810800000,  0.039507280000),
    ML_VisualXYZ(0.189140000000,  0.884962400000,  0.036935640000),
    ML_VisualXYZ(0.201169400000,  0.895493600000,  0.034458360000),
    ML_VisualXYZ(0.213365800000,  0.905443200000,  0.032088720000),
    ML_VisualXYZ(0.225749900000,  0.914850100000,  0.029840000000),
    ML_VisualXYZ(0.238320900000,  0.923734800000,  0.027711810000),
    ML_VisualXYZ(0.251066800000,  0.932092400000,  0.025694440000),
    ML_VisualXYZ(0.263992200000,  0.939922600000,  0.023787160000),
    ML_VisualXYZ(0.277101700000,  0.947225200000,  0.021989250000),
    ML_VisualXYZ(0.290400000000,  0.954000000000,  0.020300000000),   // 540 nm
    ML_VisualXYZ(0.303891200000,  0.960256100000,  0.018718050000),
    ML_VisualXYZ(0.317572600000,  0.966007400000,  0.017240360000),
    ML_VisualXYZ(0.331438400000,  0.971260600000,  0.015863640000),
    ML_VisualXYZ(0.345482800000,  0.976022500000,  0.014584610000),
    ML_VisualXYZ(0.359700000000,  0.980300000000,  0.013400000000),
    ML_VisualXYZ(0.374083900000,  0.984092400000,  0.012307230000),
    ML_VisualXYZ(0.388639600000,  0.987418200000,  0.011301880000),
    ML_VisualXYZ(0.403378400000,  0.990312800000,  0.010377920000),
    ML_VisualXYZ(0.418311500000,  0.992811600000,  0.009529306000),
    ML_VisualXYZ(0.433449900000,  0.994950100000,  0.008749999000),   // 550 nm
    ML_VisualXYZ(0.448795300000,  0.996710800000,  0.008035200000),
    ML_VisualXYZ(0.464336000000,  0.998098300000,  0.007381600000),
    ML_VisualXYZ(0.480064000000,  0.999112000000,  0.006785400000),
    ML_VisualXYZ(0.495971300000,  0.999748200000,  0.006242800000),
    ML_VisualXYZ(0.512050100000,  1.000000000000,  0.005749999000),
    ML_VisualXYZ(0.528295900000,  0.999856700000,  0.005303600000),
    ML_VisualXYZ(0.544691600000,  0.999304600000,  0.004899800000),
    ML_VisualXYZ(0.561209400000,  0.998325500000,  0.004534200000),
    ML_VisualXYZ(0.577821500000,  0.996898700000,  0.004202400000),
    ML_VisualXYZ(0.594500000000,  0.995000000000,  0.003900000000),   // 560 nm
    ML_VisualXYZ(0.611220900000,  0.992600500000,  0.003623200000),
    ML_VisualXYZ(0.627975800000,  0.989742600000,  0.003370600000),
    ML_VisualXYZ(0.644760200000,  0.986444400000,  0.003141400000),
    ML_VisualXYZ(0.661569700000,  0.982724100000,  0.002934800000),
    ML_VisualXYZ(0.678400000000,  0.978600000000,  0.002749999000),
    ML_VisualXYZ(0.695239200000,  0.974083700000,  0.002585200000),
    ML_VisualXYZ(0.712058600000,  0.969171200000,  0.002438600000),
    ML_VisualXYZ(0.728828400000,  0.963856800000,  0.002309400000),
    ML_VisualXYZ(0.745518800000,  0.958134900000,  0.002196800000),
    ML_VisualXYZ(0.762100000000,  0.952000000000,  0.002100000000),   // 570 nm
    ML_VisualXYZ(0.778543200000,  0.945450400000,  0.002017733000),
    ML_VisualXYZ(0.794825600000,  0.938499200000,  0.001948200000),
    ML_VisualXYZ(0.810926400000,  0.931162800000,  0.001889800000),
    ML_VisualXYZ(0.826824800000,  0.923457600000,  0.001840933000),
    ML_VisualXYZ(0.842500000000,  0.915400000000,  0.001800000000),
    ML_VisualXYZ(0.857932500000,  0.907006400000,  0.001766267000),
    ML_VisualXYZ(0.873081600000,  0.898277200000,  0.001737800000),
    ML_VisualXYZ(0.887894400000,  0.889204800000,  0.001711200000),
    ML_VisualXYZ(0.902318100000,  0.879781600000,  0.001683067000),
    ML_VisualXYZ(0.916300000000,  0.870000000000,  0.001650001000),   // 580 nm
    ML_VisualXYZ(0.929799500000,  0.859861300000,  0.001610133000),
    ML_VisualXYZ(0.942798400000,  0.849392000000,  0.001564400000),
    ML_VisualXYZ(0.955277600000,  0.838622000000,  0.001513600000),
    ML_VisualXYZ(0.967217900000,  0.827581300000,  0.001458533000),
    ML_VisualXYZ(0.978600000000,  0.816300000000,  0.001400000000),
    ML_VisualXYZ(0.989385600000,  0.804794700000,  0.001336667000),
    ML_VisualXYZ(0.999548800000,  0.793082000000,  0.001270000000),
    ML_VisualXYZ(1.009089200000,  0.781192000000,  0.001205000000),
    ML_VisualXYZ(1.018006400000,  0.769154700000,  0.001146667000),
    ML_VisualXYZ(1.026300000000,  0.757000000000,  0.001100000000),   // 590 nm
    ML_VisualXYZ(1.033982700000,  0.744754100000,  0.001068800000),
    ML_VisualXYZ(1.040986000000,  0.732422400000,  0.001049400000),
    ML_VisualXYZ(1.047188000000,  0.720003600000,  0.001035600000),
    ML_VisualXYZ(1.052466700000,  0.707496500000,  0.001021200000),
    ML_VisualXYZ(1.056700000000,  0.694900000000,  0.001000000000),
    ML_VisualXYZ(1.059794400000,  0.682219200000,  0.000968640000),
    ML_VisualXYZ(1.061799200000,  0.669471600000,  0.000929920000),
    ML_VisualXYZ(1.062806800000,  0.656674400000,  0.000886880000),
    ML_VisualXYZ(1.062909600000,  0.643844800000,  0.000842560000),
    ML_VisualXYZ(1.062200000000,  0.631000000000,  0.000800000000),   // 600 nm
    ML_VisualXYZ(1.060735200000,  0.618155500000,  0.000760960000),
    ML_VisualXYZ(1.058443600000,  0.605314400000,  0.000723680000),
    ML_VisualXYZ(1.055224400000,  0.592475600000,  0.000685920000),
    ML_VisualXYZ(1.050976800000,  0.579637900000,  0.000645440000),
    ML_VisualXYZ(1.045600000000,  0.566800000000,  0.000600000000),
    ML_VisualXYZ(1.039036900000,  0.553961100000,  0.000547866700),
    ML_VisualXYZ(1.031360800000,  0.541137200000,  0.000491600000),
    ML_VisualXYZ(1.022666200000,  0.528352800000,  0.000435400000),
    ML_VisualXYZ(1.013047700000,  0.515632300000,  0.000383466700),
    ML_VisualXYZ(1.002600000000,  0.503000000000,  0.000340000000),   // 610 nm
    ML_VisualXYZ(0.991367500000,  0.490468800000,  0.000307253300),
    ML_VisualXYZ(0.979331400000,  0.478030400000,  0.000283160000),
    ML_VisualXYZ(0.966491600000,  0.465677600000,  0.000265440000),
    ML_VisualXYZ(0.952847900000,  0.453403200000,  0.000251813300),
    ML_VisualXYZ(0.938400000000,  0.441200000000,  0.000240000000),
    ML_VisualXYZ(0.923194000000,  0.429080000000,  0.000229546700),
    ML_VisualXYZ(0.907244000000,  0.417036000000,  0.000220640000),
    ML_VisualXYZ(0.890502000000,  0.405032000000,  0.000211960000),
    ML_VisualXYZ(0.872920000000,  0.393032000000,  0.000202186700),
    ML_VisualXYZ(0.854449900000,  0.381000000000,  0.000190000000),   // 620 nm
    ML_VisualXYZ(0.835084000000,  0.368918400000,  0.000174213300),
    ML_VisualXYZ(0.814946000000,  0.356827200000,  0.000155640000),
    ML_VisualXYZ(0.794186000000,  0.344776800000,  0.000135960000),
    ML_VisualXYZ(0.772954000000,  0.332817600000,  0.000116853300),
    ML_VisualXYZ(0.751400000000,  0.321000000000,  0.000100000000),
    ML_VisualXYZ(0.729583600000,  0.309338100000,  0.000086133330),
    ML_VisualXYZ(0.707588800000,  0.297850400000,  0.000074600000),
    ML_VisualXYZ(0.685602200000,  0.286593600000,  0.000065000000),
    ML_VisualXYZ(0.663810400000,  0.275624500000,  0.000056933330),
    ML_VisualXYZ(0.642400000000,  0.265000000000,  0.000049999990),   // 630 nm
    ML_VisualXYZ(0.621514900000,  0.254763200000,  0.000044160000),
    ML_VisualXYZ(0.601113800000,  0.244889600000,  0.000039480000),
    ML_VisualXYZ(0.581105200000,  0.235334400000,  0.000035720000),
    ML_VisualXYZ(0.561397700000,  0.226052800000,  0.000032640000),
    ML_VisualXYZ(0.541900000000,  0.217000000000,  0.000030000000),
    ML_VisualXYZ(0.522599500000,  0.208161600000,  0.000027653330),
    ML_VisualXYZ(0.503546400000,  0.199548800000,  0.000025560000),
    ML_VisualXYZ(0.484743600000,  0.191155200000,  0.000023640000),
    ML_VisualXYZ(0.466193900000,  0.182974400000,  0.000021813330),
    ML_VisualXYZ(0.447900000000,  0.175000000000,  0.000020000000),   // 640 nm
    ML_VisualXYZ(0.429861300000,  0.167223500000,  0.000018133330),
    ML_VisualXYZ(0.412098000000,  0.159646400000,  0.000016200000),
    ML_VisualXYZ(0.394644000000,  0.152277600000,  0.000014200000),
    ML_VisualXYZ(0.377533300000,  0.145125900000,  0.000012133330),
    ML_VisualXYZ(0.360800000000,  0.138200000000,  0.000010000000),
    ML_VisualXYZ(0.344456300000,  0.131500300000,  0.000007733333),
    ML_VisualXYZ(0.328516800000,  0.125024800000,  0.000005400000),
    ML_VisualXYZ(0.313019200000,  0.118779200000,  0.000003200000),
    ML_VisualXYZ(0.298001100000,  0.112769100000,  0.000001333333),
    ML_VisualXYZ(0.283500000000,  0.107000000000,  0.000000000000),   // 650 nm
    ML_VisualXYZ(0.269544800000,  0.101476200000,  0.000000000000),
    ML_VisualXYZ(0.256118400000,  0.096188640000,  0.000000000000),
    ML_VisualXYZ(0.243189600000,  0.091122960000,  0.000000000000),
    ML_VisualXYZ(0.230727200000,  0.086264850000,  0.000000000000),
    ML_VisualXYZ(0.218700000000,  0.081600000000,  0.000000000000),
    ML_VisualXYZ(0.207097100000,  0.077120640000,  0.000000000000),
    ML_VisualXYZ(0.195923200000,  0.072825520000,  0.000000000000),
    ML_VisualXYZ(0.185170800000,  0.068710080000,  0.000000000000),
    ML_VisualXYZ(0.174832300000,  0.064769760000,  0.000000000000),
    ML_VisualXYZ(0.164900000000,  0.061000000000,  0.000000000000),   // 660 nm
    ML_VisualXYZ(0.155366700000,  0.057396210000,  0.000000000000),
    ML_VisualXYZ(0.146230000000,  0.053955040000,  0.000000000000),
    ML_VisualXYZ(0.137490000000,  0.050673760000,  0.000000000000),
    ML_VisualXYZ(0.129146700000,  0.047549650000,  0.000000000000),
    ML_VisualXYZ(0.121200000000,  0.044580000000,  0.000000000000),
    ML_VisualXYZ(0.113639700000,  0.041758720000,  0.000000000000),
    ML_VisualXYZ(0.106465000000,  0.039084960000,  0.000000000000),
    ML_VisualXYZ(0.099690440000,  0.036563840000,  0.000000000000),
    ML_VisualXYZ(0.093330610000,  0.034200480000,  0.000000000000),
    ML_VisualXYZ(0.087400000000,  0.032000000000,  0.000000000000),   // 670 nm
    ML_VisualXYZ(0.081900960000,  0.029962610000,  0.000000000000),
    ML_VisualXYZ(0.076804280000,  0.028076640000,  0.000000000000),
    ML_VisualXYZ(0.072077120000,  0.026329360000,  0.000000000000),
    ML_VisualXYZ(0.067686640000,  0.024708050000,  0.000000000000),
    ML_VisualXYZ(0.063600000000,  0.023200000000,  0.000000000000),
    ML_VisualXYZ(0.059806850000,  0.021800770000,  0.000000000000),
    ML_VisualXYZ(0.056282160000,  0.020501120000,  0.000000000000),
    ML_VisualXYZ(0.052971040000,  0.019281080000,  0.000000000000),
    ML_VisualXYZ(0.049818610000,  0.018120690000,  0.000000000000),
    ML_VisualXYZ(0.046770000000,  0.017000000000,  0.000000000000),   // 680 nm
    ML_VisualXYZ(0.043784050000,  0.015903790000,  0.000000000000),
    ML_VisualXYZ(0.040875360000,  0.014837180000,  0.000000000000),
    ML_VisualXYZ(0.038072640000,  0.013810680000,  0.000000000000),
    ML_VisualXYZ(0.035404610000,  0.012834780000,  0.000000000000),
    ML_VisualXYZ(0.032900000000,  0.011920000000,  0.000000000000),
    ML_VisualXYZ(0.030564190000,  0.011068310000,  0.000000000000),
    ML_VisualXYZ(0.028380560000,  0.010273390000,  0.000000000000),
    ML_VisualXYZ(0.026344840000,  0.009533311000,  0.000000000000),
    ML_VisualXYZ(0.024452750000,  0.008846157000,  0.000000000000),
    ML_VisualXYZ(0.022700000000,  0.008210000000,  0.000000000000),   // 690 nm
    ML_VisualXYZ(0.021084290000,  0.007623781000,  0.000000000000),
    ML_VisualXYZ(0.019599880000,  0.007085424000,  0.000000000000),
    ML_VisualXYZ(0.018237320000,  0.006591476000,  0.000000000000),
    ML_VisualXYZ(0.016987170000,  0.006138485000,  0.000000000000),
    ML_VisualXYZ(0.015840000000,  0.005723000000,  0.000000000000),
    ML_VisualXYZ(0.014790640000,  0.005343059000,  0.000000000000),
    ML_VisualXYZ(0.013831320000,  0.004995796000,  0.000000000000),
    ML_VisualXYZ(0.012948680000,  0.004676404000,  0.000000000000),
    ML_VisualXYZ(0.012129200000,  0.004380075000,  0.000000000000),
    ML_VisualXYZ(0.011359160000,  0.004102000000,  0.000000000000),   // 700 nm
    ML_VisualXYZ(0.010629350000,  0.003838453000,  0.000000000000),
    ML_VisualXYZ(0.009938846000,  0.003589099000,  0.000000000000),
    ML_VisualXYZ(0.009288422000,  0.003354219000,  0.000000000000),
    ML_VisualXYZ(0.008678854000,  0.003134093000,  0.000000000000),
    ML_VisualXYZ(0.008110916000,  0.002929000000,  0.000000000000),
    ML_VisualXYZ(0.007582388000,  0.002738139000,  0.000000000000),
    ML_VisualXYZ(0.007088746000,  0.002559876000,  0.000000000000),
    ML_VisualXYZ(0.006627313000,  0.002393244000,  0.000000000000),
    ML_VisualXYZ(0.006195408000,  0.002237275000,  0.000000000000),
    ML_VisualXYZ(0.005790346000,  0.002091000000,  0.000000000000),   // 710 nm
    ML_VisualXYZ(0.005409826000,  0.001953587000,  0.000000000000),
    ML_VisualXYZ(0.005052583000,  0.001824580000,  0.000000000000),
    ML_VisualXYZ(0.004717512000,  0.001703580000,  0.000000000000),
    ML_VisualXYZ(0.004403507000,  0.001590187000,  0.000000000000),
    ML_VisualXYZ(0.004109457000,  0.001484000000,  0.000000000000),
    ML_VisualXYZ(0.003833913000,  0.001384496000,  0.000000000000),
    ML_VisualXYZ(0.003575748000,  0.001291268000,  0.000000000000),
    ML_VisualXYZ(0.003334342000,  0.001204092000,  0.000000000000),
    ML_VisualXYZ(0.003109075000,  0.001122744000,  0.000000000000),
    ML_VisualXYZ(0.002899327000,  0.001047000000,  0.000000000000),   // 720 nm
    ML_VisualXYZ(0.002704348000,  0.000976589600,  0.000000000000),
    ML_VisualXYZ(0.002523020000,  0.000911108800,  0.000000000000),
    ML_VisualXYZ(0.002354168000,  0.000850133200,  0.000000000000),
    ML_VisualXYZ(0.002196616000,  0.000793238400,  0.000000000000),
    ML_VisualXYZ(0.002049190000,  0.000740000000,  0.000000000000),
    ML_VisualXYZ(0.001910960000,  0.000690082700,  0.000000000000),
    ML_VisualXYZ(0.001781438000,  0.000643310000,  0.000000000000),
    ML_VisualXYZ(0.001660110000,  0.000599496000,  0.000000000000),
    ML_VisualXYZ(0.001546459000,  0.000558454700,  0.000000000000),
    ML_VisualXYZ(0.001439971000,  0.000520000000,  0.000000000000),   // 730 nm
    ML_VisualXYZ(0.001340042000,  0.000483913600,  0.000000000000),
    ML_VisualXYZ(0.001246275000,  0.000450052800,  0.000000000000),
    ML_VisualXYZ(0.001158471000,  0.000418345200,  0.000000000000),
    ML_VisualXYZ(0.001076430000,  0.000388718400,  0.000000000000),
    ML_VisualXYZ(0.000999949300,  0.000361100000,  0.000000000000),
    ML_VisualXYZ(0.000928735800,  0.000335383500,  0.000000000000),
    ML_VisualXYZ(0.000862433200,  0.000311440400,  0.000000000000),
    ML_VisualXYZ(0.000800750300,  0.000289165600,  0.000000000000),
    ML_VisualXYZ(0.000743396000,  0.000268453900,  0.000000000000),
    ML_VisualXYZ(0.000690078600,  0.000249200000,  0.000000000000),   // 740 nm
    ML_VisualXYZ(0.000640515600,  0.000231301900,  0.000000000000),
    ML_VisualXYZ(0.000594502100,  0.000214685600,  0.000000000000),
    ML_VisualXYZ(0.000551864600,  0.000199288400,  0.000000000000),
    ML_VisualXYZ(0.000512429000,  0.000185047500,  0.000000000000),
    ML_VisualXYZ(0.000476021300,  0.000171900000,  0.000000000000),
    ML_VisualXYZ(0.000442453600,  0.000159778100,  0.000000000000),
    ML_VisualXYZ(0.000411511700,  0.000148604400,  0.000000000000),
    ML_VisualXYZ(0.000382981400,  0.000138301600,  0.000000000000),
    ML_VisualXYZ(0.000356649100,  0.000128792500,  0.000000000000),
    ML_VisualXYZ(0.000332301100,  0.000120000000,  0.000000000000),   // 750 nm
    ML_VisualXYZ(0.000309758600,  0.000111859500,  0.000000000000),
    ML_VisualXYZ(0.000288887100,  0.000104322400,  0.000000000000),
    ML_VisualXYZ(0.000269539400,  0.000097335600,  0.000000000000),
    ML_VisualXYZ(0.000251568200,  0.000090845870,  0.000000000000),
    ML_VisualXYZ(0.000234826100,  0.000084800000,  0.000000000000),
    ML_VisualXYZ(0.000219171000,  0.000079146670,  0.000000000000),
    ML_VisualXYZ(0.000204525800,  0.000073858000,  0.000000000000),
    ML_VisualXYZ(0.000190840500,  0.000068916000,  0.000000000000),
    ML_VisualXYZ(0.000178065400,  0.000064302670,  0.000000000000),
    ML_VisualXYZ(0.000166150500,  0.000060000000,  0.000000000000),   // 760 nm
    ML_VisualXYZ(0.000155023600,  0.000055981870,  0.000000000000),
    ML_VisualXYZ(0.000144621900,  0.000052225600,  0.000000000000),
    ML_VisualXYZ(0.000134909800,  0.000048718400,  0.000000000000),
    ML_VisualXYZ(0.000125852000,  0.000045447470,  0.000000000000),
    ML_VisualXYZ(0.000117413000,  0.000042400000,  0.000000000000),
    ML_VisualXYZ(0.000109551500,  0.000039561040,  0.000000000000),
    ML_VisualXYZ(0.000102224500,  0.000036915120,  0.000000000000),
    ML_VisualXYZ(0.000095394450,  0.000034448680,  0.000000000000),
    ML_VisualXYZ(0.000089023900,  0.000032148160,  0.000000000000),
    ML_VisualXYZ(0.000083075270,  0.000030000000,  0.000000000000),   // 770 nm
    ML_VisualXYZ(0.000077512690,  0.000027991250,  0.000000000000),
    ML_VisualXYZ(0.000072313040,  0.000026113560,  0.000000000000),
    ML_VisualXYZ(0.000067457780,  0.000024360240,  0.000000000000),
    ML_VisualXYZ(0.000062928440,  0.000022724610,  0.000000000000),
    ML_VisualXYZ(0.000058706520,  0.000021200000,  0.000000000000),
    ML_VisualXYZ(0.000054770280,  0.000019778550,  0.000000000000),
    ML_VisualXYZ(0.000051099180,  0.000018452850,  0.000000000000),
    ML_VisualXYZ(0.000047676540,  0.000017216870,  0.000000000000),
    ML_VisualXYZ(0.000044485670,  0.000016064590,  0.000000000000),
    ML_VisualXYZ(0.000041509940,  0.000014990000,  0.000000000000),   // 780 nm
    ML_VisualXYZ(0.000038733240,  0.000013987280,  0.000000000000),
    ML_VisualXYZ(0.000036142030,  0.000013051550,  0.000000000000),
    ML_VisualXYZ(0.000033723520,  0.000012178180,  0.000000000000),
    ML_VisualXYZ(0.000031464870,  0.000011362540,  0.000000000000),
    ML_VisualXYZ(0.000029353260,  0.000010600000,  0.000000000000),
    ML_VisualXYZ(0.000027375730,  0.000009885877,  0.000000000000),
    ML_VisualXYZ(0.000025524330,  0.000009217304,  0.000000000000),
    ML_VisualXYZ(0.000023793760,  0.000008592362,  0.000000000000),
    ML_VisualXYZ(0.000022178700,  0.000008009133,  0.000000000000),
    ML_VisualXYZ(0.000020673830,  0.000007465700,  0.000000000000),   // 790 nm
    ML_VisualXYZ(0.000019272260,  0.000006959567,  0.000000000000),
    ML_VisualXYZ(0.000017966400,  0.000006487995,  0.000000000000),
    ML_VisualXYZ(0.000016749910,  0.000006048699,  0.000000000000),
    ML_VisualXYZ(0.000015616480,  0.000005639396,  0.000000000000),
    ML_VisualXYZ(0.000014559770,  0.000005257800,  0.000000000000),
    ML_VisualXYZ(0.000013573870,  0.000004901771,  0.000000000000),
    ML_VisualXYZ(0.000012654360,  0.000004569720,  0.000000000000),
    ML_VisualXYZ(0.000011797230,  0.000004260194,  0.000000000000),
    ML_VisualXYZ(0.000010998440,  0.000003971739,  0.000000000000),
    ML_VisualXYZ(0.000010253980,  0.000003702900,  0.000000000000),   // 800 nm
    ML_VisualXYZ(0.000009559646,  0.000003452163,  0.000000000000),
    ML_VisualXYZ(0.000008912044,  0.000003218302,  0.000000000000),
    ML_VisualXYZ(0.000008308358,  0.000003000300,  0.000000000000),
    ML_VisualXYZ(0.000007745769,  0.000002797139,  0.000000000000),
    ML_VisualXYZ(0.000007221456,  0.000002607800,  0.000000000000),
    ML_VisualXYZ(0.000006732475,  0.000002431220,  0.000000000000),
    ML_VisualXYZ(0.000006276423,  0.000002266531,  0.000000000000),
    ML_VisualXYZ(0.000005851304,  0.000002113013,  0.000000000000),
    ML_VisualXYZ(0.000005455118,  0.000001969943,  0.000000000000),
    ML_VisualXYZ(0.000005085868,  0.000001836600,  0.000000000000),   // 810 nm
    ML_VisualXYZ(0.000004741466,  0.000001712230,  0.000000000000),
    ML_VisualXYZ(0.000004420236,  0.000001596228,  0.000000000000),
    ML_VisualXYZ(0.000004120783,  0.000001488090,  0.000000000000),
    ML_VisualXYZ(0.000003841716,  0.000001387314,  0.000000000000),
    ML_VisualXYZ(0.000003581652,  0.000001293400,  0.000000000000),
    ML_VisualXYZ(0.000003339127,  0.000001205820,  0.000000000000),
    ML_VisualXYZ(0.000003112949,  0.000001124143,  0.000000000000),
    ML_VisualXYZ(0.000002902121,  0.000001048009,  0.000000000000),
    ML_VisualXYZ(0.000002705645,  0.000000977058,  0.000000000000),
    ML_VisualXYZ(0.000002522525,  0.000000910930,  0.000000000000),   // 820 nm
    ML_VisualXYZ(0.000002351726,  0.000000849251,  0.000000000000),
    ML_VisualXYZ(0.000002192415,  0.000000791721,  0.000000000000),
    ML_VisualXYZ(0.000002043902,  0.000000738090,  0.000000000000),
    ML_VisualXYZ(0.000001905497,  0.000000688110,  0.000000000000),
    ML_VisualXYZ(0.000001776509,  0.000000641530,  0.000000000000),
    ML_VisualXYZ(0.000001656215,  0.000000598090,  0.000000000000),
    ML_VisualXYZ(0.000001544022,  0.000000557575,  0.000000000000),
    ML_VisualXYZ(0.000001439440,  0.000000519808,  0.000000000000),
    ML_VisualXYZ(0.000001341977,  0.000000484612,  0.000000000000),
    ML_VisualXYZ(0.000001251141,  0.000000451810,  0.000000000000),   // 830 nm
};

void
BuildSpectraXYZ()
{
    static ML_SpectralDensity sd;
    Spectrum s;
    int i;

    for (i = 0; i <= 830-360; i += 5)
        sd.rgf[i/5] = ML_CIE[i].X;
    s = MLtoSpectrum(sd);
    PrintSpectrum("g_sX", s);

    for (i = 0; i <= 830-360; i += 5)
        sd.rgf[i/5] = ML_CIE[i].Y;
    s = MLtoSpectrum(sd);
    PrintSpectrum("g_sY", s);

    for (i = 0; i <= 830-360; i += 5)
        sd.rgf[i/5] = ML_CIE[i].Z;
    s = MLtoSpectrum(sd);
    PrintSpectrum("g_sZ", s);
}