#include "stdafx.h"
#include "Sampling.h"
#include "Image.h"

#define ML_SAMPLE_TO_INDEX(X)   int(floor(X))
#define ML_INDEX_TO_SAMPLE(I)   (double(I) + 0.5)
#define ML_LOWEST_INDEX(X)      int(floor(X + 0.5))
#define ML_HIGHEST_INDEX(X)     int(floor(X - 0.5))
#define ML_KAISER4_HALFWIDTH    4.0f

#define IMAGESIZE   512
#define D_SQRT2     1.41421356237309504880168872420969807856967187537694
#define D_PI        3.14159265358979323846264338327950288419716939937510

typedef float   (*ML_PFILTER)(float);

double
ImageZone(double x, double y)
{
    x = x - 0.5*double(IMAGESIZE);
    y = y - 0.5*double(IMAGESIZE);
    x /= double(IMAGESIZE)/D_SQRT2;
    y /= double(IMAGESIZE)/D_SQRT2;
    return 0.5*cos(D_PI*IMAGESIZE*(x*x + y*y)) + 0.5;
}

double
ImageFactory(double x, double y)
{
    static Image im;
    static int bLoaded;

    if (!bLoaded) {
        if (im.ReadBMP("ImageFactory.bmp") == 0) {
            printf("Cannot open image\n");
            exit(0);
        }
        bLoaded = 1;
    }
    x = im.m_nWidth*x/double(IMAGESIZE);
    y = im.m_nHeight*y/double(IMAGESIZE);
    return im.Sample(x, y);
}
//
//  Filter down images into perfect subsampled forms
//
static int s_nLastSrc;
static int s_nLastDst;
static ML_PFILTER s_pLastFilter = 0;
static float *s_pfTap = 0;
static int s_nTaps = 0;
static double gs_fI0Alpha = 1.0/11.3019219521363;
static double gs_fAlphaSquared = 4.0*4.0;

static float
ML_Kaiser4Filter(float x)
{
    float x1, x2, xTerm, fN, fSinc, fKaiser;

    x1 = x * 1.0f/ML_KAISER4_HALFWIDTH;
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
    x1 = D_PI*x;
    x2 = x1*x1;
    if (x2 < 0.0001f)
        fSinc = 1.0f + x2*(-1.0f/6.0f + x2*1.0f/120.0f);
    else
        fSinc = sinf(x1)/x1;
    return fSinc * fKaiser * float(gs_fI0Alpha); // * 1.0033f;    // correct normalization
}
//
//  caching the filter coefficients speeds up ML_ImageResample by orders of magnitude
//
static void
ComputeTaps(int nDst, int nSrc, double xHalfWidth, ML_PFILTER pFilter)
{
    double xSrc, x, fRatio, fDecimation, fScale;
    int iDst, iSrc, iFirst, iLast, nTaps;
    float *pf;

    fRatio = float(nSrc)/float(nDst);
    if (nDst < nSrc)
        fDecimation = fRatio;
    else
        fDecimation = 1.0;
    xHalfWidth *= fDecimation;
    fScale = 1.0/fDecimation;
    nTaps = nDst * int(2.0*xHalfWidth + 1.0);
    if (s_pfTap == 0 || s_nTaps < nTaps) {
        delete [] s_pfTap;
        s_pfTap = new float[nTaps];
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
ML_Resample(float rgfDst[], int nDst, const float rgfSrc[], int nSrc, int nWrap,
            double xHalfWidth, ML_PFILTER pFilter)
{
    double xSrc, y, fRatio, fDecimation, fScale, f, fSum, fWeight;
    int iDst, i, iSrc, iFirst, iLast;
    float *pf;

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
ML_ImageResample(Image &imgNew, int nNewWidth, int nNewHeight, const Image &imgSrc,
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
            ML_Resample(pfRow, nNewWidth, pfSrc + j*nOldWidth, nOldWidth, imgSrc.m_bWrapX, xHalfWidth, pFilter);
            for (i = 0; i < nNewWidth; i++)
                pfTransposed[j + i*nOldHeight] = pfRow[i];
        }
        //
        //  Resample vertically
        //
        pfDst = imgNew.m_pfImage + k*nNewWidth*nNewHeight;
        for (j = 0; j < nNewWidth; j++) {
            ML_Resample(pfCol, nNewHeight, pfTransposed + j*nOldHeight, nOldHeight, imgSrc.m_bWrapY, xHalfWidth, pFilter);
            for (i = 0; i < nNewHeight; i++)
                pfDst[j + i*nNewWidth] = pfCol[i];
        }
    }
    delete [] pfRow;
    delete [] pfCol;
    delete [] pfTransposed;
    return 1;
}

void
PerfectFactory(Image &im)
{
    Image im2;

    im2.ReadBMP("ImageFactory.bmp");
    ML_ImageResample(im, IMAGESIZE, IMAGESIZE, im2, ML_KAISER4_HALFWIDTH, ML_Kaiser4Filter);
}

void
PerfectZone(Image &im)
{
    Image im2;
    int i, j;
    double f;

    im2.NewImage(IMAGESIZE*2, IMAGESIZE*2, 1);
    for (j = 0; j < IMAGESIZE*2; j++) {
        for (i = 0; i < IMAGESIZE*2; i++) {
            f = ImageZone(0.5*(double(i) + 0.5), 0.5*(double(j) + 0.5));
            im2.Set(f, i, j);
        }
    }
    // im2.WriteBMP("BigZone.bmp");
    ML_ImageResample(im, IMAGESIZE, IMAGESIZE, im2, ML_KAISER4_HALFWIDTH, ML_Kaiser4Filter);
}
//
//  Spectrum of image
//
//
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
