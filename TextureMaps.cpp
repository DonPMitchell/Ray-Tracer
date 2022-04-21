#include "stdafx.h"
#include "TextureMaps.h"
#include "MonteCarlo.h"
#pragma intrinsic(fabs,sqrt,sin,cos,atan2,asin,acos,log,tan,pow)
#define D_DTR (D_PI/180.0)
#define D_RTD (180.0/D_PI)
#define D_2PI (2.0*D_PI)
#define EQUIRECTANGULAR 1

//
//  Propagators
//
Propagation
ColoredGlass(const Point3 &p0, const Point3 &p1, Material *pm, RayTracer *scene)
{
    Propagation p;
    double t, fTransmittance, fMuT;
    int i;

    VERBOSE printf("ColoredGlass\n");
    VERBOSE p0.Print(); 
    VERBOSE p1.Print();
    p.cSterisent = g_sBlack;
    t = Norm(p1 - p0) * 100.0;          // centimeters, absorption coefficients usually given in units of 1/cm
    for (i = 0; i < NSPEC; i++) {
        fMuT = pm->color.rgf[i] * t;
        if (fMuT > 87.0)
            fTransmittance = 0.0;       // would be less than minimum float value
        else
            fTransmittance = exp(-fMuT);
        p.cTransmittance.rgf[i] = float(fTransmittance);
    }
    VERBOSE printf("End of ColoredGlass\n");
    return p;
}

#define BACKGROUND 0.1

Propagation
Foggy(const Point3 &p0, const Point3 &p1, Material *pm, RayTracer *scene)
{
    Propagation p;
    double t, fScatter, fTransmittance;
    int i;

    t = Norm(p1 - p0);
    for (i = 0; i < NSPEC; i++) {
        fScatter = pm->color.rgf[i] * t;
        if (fScatter > 87.0)
            fTransmittance = 0.0;
        else
            fTransmittance = exp(-fScatter);
         p.cTransmittance.rgf[i] = float(fTransmittance);
         p.cSterisent.rgf[i] = float((1.0-fTransmittance)*BACKGROUND);
    }
    return p;
}
//
//  Rayleigh scattering propagator
//
#define AVOGADRO            6.02214199e+23      // per Mol
#define BOLTZMANN           1.3806503e-23       // Joules per Kelvin
#define MOLARGAS_R          8.31446261815324    // Joules per Kelvin per Mol  (known more accurately than k*A)
#define VACUUM_PERMITTIVITY 8.854187812e-12     // Farad per Meter

#define ALPHA_N2            1.710               // molecular polarizability in Aangstrom3 = 1.0e-30 m3
#define ALPHA_O2            1.562
#define ALPHA_CO2           2.507
#define ALPHA_AR            1.664
#define MOLWEIGHT_N2        28.0134             // g/Mol
#define MOLWEIGHT_O2        31.9988
#define MOLWEIGHT_CO2       44.01
#define MOLWEIGHT_AR        39.96
#define DENSITY_O2          1.429               // kg/m3
#define DENSITY_N2          1.2506
#define DENSITY_CO2         1.977
#define DENSITY_AR          1.784
#define CROSSSECTION_O2     5.08e-28            // 532 nm, 273 K 1.013e5 Pa (for testing)
#define CROSSSECTION_N2     6.13e-28
#define CROSSSECTION_CO2    13.8e-28

#define EARTH_N2            0.78084
#define EARTH_O2            0.20946
#define EARTH_AR            0.00934
#define EARTH_CO2           0.00041
#define MARS_CO2            0.95
#define MARS_N2             0.028
#define MARS_AR             0.02
#define MARS_O2             0.002
#define VENUS_CO2           0.965
#define VENUS_N2            0.035

#define MERCURYRADIUS       2.4397e+6           // Meters
#define MOONRADIUS          1.73753e+6          // Meters

#define VENUS_SCALE_HEIGHT  15.9
#define VENUS_DENSITY       64.79               // Surface atmosphere density, Kg/M**3
#define VENUS_MOLWEIGHT     43.45               // Mean Mol Weight, grams/mole
#define VENUS_REFRACT       1.014716            // Surface refractive index
#define VENUS_RADIUS        6.05184e+6          // Meters

#define EARTH_SCALE_HEIGHT  8.5
#define EARTH_HAZE_SCALE    1.2
#define EARTH_DENSITY       1.217               // Kg/M**3  
#define EARTH_MOLWEIGHT     28.97               // Grams/Mole
#define REFRACTION_AIR      1.000271375         // dimensionless
#define POLARIZABILITY_AIR  2.133e-29           // m3
#define EARTH_RADIUS        6378136.55          // IAU 1994

#define MARS_SCALE_HEIGHT   11.1
#define MARS_DENSITY        0.020               // Kg/M**3
#define MARS_MOLWEIGHT      43.34               // grams/mole
#define MARS_RADIUS         3.397e+6            // Meters

#define LAYERS  300

static inline double
RefractiveIndex(double fPolarizabilityAlpha, double fNumberDensity)
{
    return sqrt(1.0 + fPolarizabilityAlpha * fNumberDensity);       // often approximated by 1 + alpha*N/2
}

static void
Rayleigh(double fScaleHeight, double fSurfDensity, double fMoleWeight, double fRefractive)
{
    int i;
    double fKm, fDensity, fFraction, fK;

    fRefractive = fRefractive*fRefractive;
    fRefractive = (fRefractive - 1.0)/(fRefractive + 2.0);
    fRefractive = fRefractive*fRefractive;
    fDensity = fSurfDensity;                            // Kg/M**3
    fDensity = 1000.0*fDensity/fMoleWeight;             // Moles/M**3
    fDensity *= AVOGADRO;                             // Molecules/M**3
    fK = 1.0e36 * 1.0e3 * (24.0*D_PI*D_PI*D_PI)*fRefractive/fDensity;
  //  s_fK = fK;
    for (i = 0; i < LAYERS; i++) {
     //   fKm = Kilometers(i);
        fFraction = exp(-fKm/fScaleHeight);
        //fDensity = fSurfDensity*fFraction;                  // Kg/M**3
        //fDensity = 1000.0*fDensity/fMoleWeight;             // Moles/M**3
        //fDensity *= C_AVOGADRO;                             // Molecules/M**3
        //
        //  This gets divided by wavelength**4 in nanometers and we want
        //  the resulting scattering coefficient in reciprocal kilometers
        //
       // s_rgfRayleigh[i] = fFraction * fK;
    }
}

void
TestRayleigh()
{
    double fN;

    fN = DENSITY_O2;        // kg/m3
    fN /= MOLWEIGHT_O2;     // moles/m3
    fN *= AVOGADRO;         // particles/m3
    printf("N = %g\n", fN);
}

Propagation
Atmosphere(const Point3 &p0, const Point3 &p1, Material *pm, RayTracer *scene)
{
    Propagation p;

    return p;
}

//
//  Texture Maps
//
extern CoVector3 rgvDodecahedron[20];
extern CoVector3 rgvIcosahedron[12];
//
//  Bump mapping textures
//
double
BumpNoise(Point3 p)
{
    float rgf[3];

    rgf[0] = float(p.x);
    rgf[1] = float(p.y);
    rgf[2] = float(p.z);
    return BandNoise(rgf, 3);
}

inline void
PerturbNormal(Hit *ph, CoVector3 cvP)
{
    if (ph->nPerturbed) {
        ph->cvPerturbed = ph->cvPerturbed + cvP;
    } else {
        ph->cvPerturbed = Normalize(ph->cvNormal) + cvP;
        ph->nPerturbed++;
    }
}

inline CoVector3
GetNormal(Hit *ph)
{
    if (ph->nPerturbed)
        return ph->cvPerturbed;
    else
        return ph->cvNormal;
}

void
Bumpy(const Point3 &p, Hit *ph, RayTracer *scene)
{
    float rgf[3], rgfGrad[3];

    rgf[0] = float(p.x * 50.0);
    rgf[1] = float(p.y * 50.0);
    rgf[2] = float(p.z * 50.0);
    GradientBandNoise(rgfGrad, rgf, 3);
    PerturbNormal(ph, CoVector3(rgfGrad[0], rgfGrad[1], rgfGrad[2]) * 0.3);
    //ph->cvPerturbed = ph->cvNormal + CoVector3(rgfGrad[0], rgfGrad[1], rgfGrad[2]) * 0.3;
    //ph->nPerturbed = 1;
}

void
MicroFacets(const Point3 &p, Hit *ph, RayTracer *scene)
{
    float rgf[3], rgfGrad[3];

    rgf[0] = float(p.x * 500.0);
    rgf[1] = float(p.y * 500.0);
    rgf[2] = float(p.z * 500.0);
    GradientBandNoise(rgfGrad, rgf, 3);
    PerturbNormal(ph, CoVector3(rgfGrad[0], rgfGrad[1], rgfGrad[2]) * 0.3);
    //ph->cvPerturbed = ph->cvNormal + CoVector3(rgfGrad[0], rgfGrad[1], rgfGrad[2]) * 0.3;
    //ph->nPerturbed = 1;
}

void
MicroFacets2(const Point3 &p, Hit *ph, RayTracer *scene)
{
    float rgf[3], rgfGrad[3];

    rgf[0] = float(p.x * 5000.0);
    rgf[1] = float(p.y * 5000.0);
    rgf[2] = float(p.z * 5000.0);
    GradientBandNoise(rgfGrad, rgf, 3);
    PerturbNormal(ph, CoVector3(rgfGrad[0], rgfGrad[1], rgfGrad[2]) * 0.3);
    //ph->cvPerturbed = ph->cvNormal + CoVector3(rgfGrad[0], rgfGrad[1], rgfGrad[2]) * 0.3;
    //ph->nPerturbed = 1;
}

void
QuasiMicroFacet(const Point3 &p, Hit *ph, RayTracer *scene)
{
    CoVector3 cvN, cvU, cvV;
    double u, v;

    cvN = GetNormal(ph);
    cvU = Normalize(Cross(cvN, CoVector3(0.0, 1.0, 0.0)));
    cvV = Normalize(Cross(cvN, cvU));
    u = g_uMicroFacet;
    v = g_vMicroFacet;
    u = u*u*u*u*u;
    v = v*v*v*v*v;
    PerturbNormal(ph, (cvU*u + cvV*v) * 0.25);
}

void
Wrinkled(const Point3 &p, Hit *ph, RayTracer *scene)
{
    float rgf[3], rgfGrad[3];

    rgf[0] = float(p.x * 50.0);
    rgf[1] = float(p.y * 50.0);
    rgf[2] = float(p.z * 50.0);
    GradientFractalNoise(rgfGrad, rgf, 3, 6);
    PerturbNormal(ph, CoVector3(rgfGrad[0], rgfGrad[1], rgfGrad[2]) * 0.3);
    //ph->cvPerturbed = ph->cvNormal + CoVector3(rgfGrad[0], rgfGrad[1], rgfGrad[2]) * 0.3;
    //ph->nPerturbed = 1;
}

void
Faceted20(const Point3 &p, Hit *ph, RayTracer *scene)
{
    float rgf[3], rgfGrad[3];
    CoVector3 cv, maxcv;
    double fDot, maxDot;
    int i;

    rgf[0] = float(p.x * 50.0);
    rgf[1] = float(p.y * 50.0);
    rgf[2] = float(p.z * 50.0);
    GradientBandNoise(rgfGrad, rgf, 3);
    cv = CoVector3(rgfGrad[0], rgfGrad[1], rgfGrad[2]);
    maxDot = -1000000.0;
    for (i = 0; i < 20; i++) {
        fDot = Dot(cv, rgvDodecahedron[i]);
        if (fDot > maxDot) {
            maxDot = fDot;
            maxcv = rgvDodecahedron[i];
        }
    }
    PerturbNormal(ph, maxcv * 0.15);
    //ph->cvPerturbed = ph->cvNormal + maxcv * 0.15;
    //ph->nPerturbed = 1;
}

void
Dented(const Point3 &p, Hit *ph, RayTracer *scene)
{
    float rgf[3], rgfGrad[3], f;
    Point3 p2;
    CoVector3 cv;

    p2 = p + Vector3(13.0, 17.0, 31.0);                     // uncorrelated with wrinkle function
    rgf[0] = float(p2.x * 40.0);
    rgf[1] = float(p2.y * 40.0);
    rgf[2] = float(p2.z * 40.0);
    f = GradientFractalNoise(rgfGrad, rgf, 3);
    if (f < 0.5f) {
        ph->cvPerturbed = ph->cvNormal;
        return;
    }
    cv = CoVector3(rgfGrad[0], rgfGrad[1], rgfGrad[2]);
    f = 3.0f*f*f;
    PerturbNormal(ph, cv * (1.0*f));
    //ph->cvPerturbed = ph->cvNormal +  cv * (0.6*f);
    //ph->nPerturbed = 1;
}

void
Wavy(const Point3 &p, Hit *ph, RayTracer *scene)
{
    PerturbNormal(ph, 0.3*cos(p.z/D_2PI));
   // ph->cvPerturbed = ph->cvNormal +  CoVector3(0.0, 0.0, 0.3*cos(p.z/D_2PI));
   // ph->nPerturbed = 1;
}

void
Fibers(const Point3 &p, Hit *ph, RayTracer *scene)
{
    float rgf[3], rgfGrad[3];

    rgf[0] = float(p.x * 250.0);
    rgf[1] = float(p.y * 250.0);
    GradientBandNoise(rgfGrad, rgf, 2);
    PerturbNormal(ph, CoVector3(rgfGrad[0], rgfGrad[1], 0.0) * 0.3);
    //ph->cvPerturbed = ph->cvNormal +  CoVector3(rgfGrad[0], rgfGrad[1], 0.0) * 0.3;
    //ph->nPerturbed = 1;
}

void
Blistered(const Point3 &p, Hit *ph, RayTracer *scene)
{
    float rgf[3], rgfGrad[3], f;
    CoVector3 cv;

    rgf[0] = float(p.x * 25.0);
    rgf[1] = float(p.y * 25.0);
    rgf[2] = float(p.z * 25.0);
    f = GradientBandNoise(rgfGrad, rgf, 3);
    cv = CoVector3(rgfGrad[0], rgfGrad[1], rgfGrad[2]);
    if (f > 0.5)
        f = -0.5;
    else
        f = 0.0;
    PerturbNormal(ph, cv * f);
   // ph->cvPerturbed = ph->cvNormal + cv * f;
   // ph->nPerturbed = 1;
}

void
Hammered(const Point3 &p, Hit *ph, RayTracer *scene)
{
    float rgf[3], rgfGrad[3], f;
    CoVector3 cv, cv2;

    rgf[0] = float(p.x * 8.0);
    rgf[1] = float(p.y * 8.0);
    rgf[2] = float(p.z * 8.0);
    f = GradientBandNoise(rgfGrad, rgf, 3);
    cv = CoVector3(rgfGrad[0], rgfGrad[1], rgfGrad[2]);
    f = GradientFractalNoise(rgfGrad, rgf, 3, 5);
    cv2 =  CoVector3(rgfGrad[0], rgfGrad[1], rgfGrad[2]);
    PerturbNormal(ph, cv*0.1 + cv2*0.1);
    //ph->cvPerturbed = ph->cvNormal +  cv*0.1 + cv2*0.1;
    //ph->nPerturbed = 1;
}

void
Tectonic(const Point3 &p, Hit *ph, RayTracer *scene)
{
    float rgf[3], rgfGrad[3], f;

    rgf[0] = float(p.x * 3.0);
    rgf[1] = float(p.y * 3.0);
    rgf[2] = float(p.z * 3.0);
    f = GradientFractalNoise(rgfGrad, rgf, 3, 6);
    if (f < 0.5f)
        f = 0.0f;
    else
        f = -0.4f;
    PerturbNormal(ph, CoVector3(rgfGrad[0], rgfGrad[1], rgfGrad[2]) * f);
    //ph->cvPerturbed = ph->cvNormal +  CoVector3(rgfGrad[0], rgfGrad[1], rgfGrad[2]) * f;
   // ph->nPerturbed = 1;
}

void
VoronoiBump(const Point3 &p, Hit *ph, RayTracer *scene)
{
    float rgf[3], rgfGrad[3];

    rgf[0] = float(p.x * 7.0);
    rgf[1] = float(p.y * 7.0);
    rgf[2] = float(p.z * 7.0);
    GradientVoronoi(rgfGrad, rgf);
    PerturbNormal(ph, CoVector3(rgfGrad[0], rgfGrad[1], rgfGrad[2]) * 0.3);
}

void
VoronoiAntiBump(const Point3 &p, Hit *ph, RayTracer *scene)
{
    float rgf[3], rgfGrad[3];

    rgf[0] = float(p.x * 7.0);
    rgf[1] = float(p.y * 7.0);
    rgf[2] = float(p.z * 7.0);
    GradientVoronoi(rgfGrad, rgf);
    PerturbNormal(ph, CoVector3(-rgfGrad[0], -rgfGrad[1], -rgfGrad[2]) * 0.3);
}

void
VoronoiCrinkled(const Point3 &p, Hit *ph, RayTracer *scene)
{
    float rgf[3], rgfGrad[3];

    rgf[0] = float(p.x * 7.0);
    rgf[1] = float(p.y * 7.0);
    rgf[2] = float(p.z * 7.0);
    GradientVoronoi(rgfGrad, rgf, 4);
    PerturbNormal(ph, CoVector3(-rgfGrad[0], -rgfGrad[1], -rgfGrad[2]) * 0.3);
}

double
Harmonic(Point3 p, double fChoppy)
{
    Point3 pS, pC;
    double f;

    f = BumpNoise(p);
    p.x += f;
    pS = Point3(1.0, 1.0, 0.0) - Vector3(fabs(sin(p.x)), fabs(sin(p.y)), 0.0);
    pC = Point3(fabs(cos(p.x)), fabs(cos(p.y)), 0.0);
    pS.x = pS.x + (pC.x - pS.x)*pS.x;
    pS.y = pS.y + (pC.y - pS.y)*pS.y;
    return pow(1.0 - pow(pS.x*pS.y, 0.65), fChoppy);
}

#define ITER    5
#define SEA_HEIGHT  0.6
#define SEA_CHOPPY  4.0;
#define SEA_SPEED   0.8
#define SEA_FREQ    0.16
#define SEA_TIME    (1.0 + g_tTime*SEA_SPEED)
// octave_m = mat2(1.6,1.2,-1.2,1.6)

double
SeaMap(const Point3 &p, Hit *ph, RayTracer *scene)
{
    double freq = SEA_FREQ;
    double amp = SEA_HEIGHT;
    double fChoppy = SEA_CHOPPY;
    double d, h;
    Point3 pUV, pt;
    int i;

    pUV = p;
    pUV.x *= 0.75;
    h = 0.0;
    for (i = 0; i < ITER; i++) {
        d = Harmonic(Point3((pUV.x + SEA_TIME)*freq, pUV.y, pUV.z), fChoppy);
        d += Harmonic(Point3((pUV.x - SEA_TIME)*freq, pUV.y, pUV.z), fChoppy);
        h += d * amp;
        pUV = Point3(pUV.x * 1.6 + pUV.y * 1.2, pUV.x * -1.2 + pUV.y * 1.6, p.z);
        freq *= 1.9;
        amp *= 0.22;
        fChoppy = fChoppy + (1.0 - fChoppy)*0.2;
    }
    return pUV.y - h;
}

static Vector3
Fibonacci(int i, int N)
{
    double rgf[3];

    FibonacciSphere(rgf, i, N);
    return Vector3(rgf[0], rgf[1], rgf[2]);
}

void
SeaScape(const Point3 &p, Hit *ph, RayTracer *scene)
{
    Point3 pUV;
    Vector3 v;
    int i;

    for (i = 0; i < ITER; i++) {
        v = Fibonacci(i, 5);
        pUV = p + v*BumpNoise(p);


    }
}
//
//  Solid color textures
//
void
Voronoi(const Point3 &p, Hit *ph, RayTracer *scene)
{
    double fBlend;
    float rgf[3];

    rgf[0] = float(p.x * 7.0);
    rgf[1] = float(p.y * 7.0);
    rgf[2] = float(p.z * 7.0);
    fBlend = Voronoi(rgf) / 1.6;
    if (fBlend > 1.0)
        fBlend = 1.0;
    ph->pm->color = g_sGreen*fBlend + g_sEqualWhite*(1.0-fBlend);
}

void
Agate(const Point3 &p, Hit *ph, RayTracer *scene)
{
    double fBlend;
    float rgf[3];

    rgf[0] = float(p.x * 7.0);
    rgf[1] = float(p.y * 7.0);
    rgf[2] = float(p.z * 7.0);
    fBlend = BandNoise(rgf, 3);
    fBlend = 0.5*(sin(100.0*fBlend) + 1.0);
    ph->pm->color = g_sOrange*fBlend + g_sEqualWhite*(1.0-fBlend);
}

static double
SawToothWave(double x)
{
	x = (0.5/D_PI)*x + 0.25;
	x = 2.0*(x - floor(x));
	if (x < 1.0)
		return 2.0*x - 1.0;
	else
		return 3.0 - 2.0*x;
}

void
Marble(const Point3 &p, Hit *ph, RayTracer *scene)
{
    Spectrum sPink;
    double fDensity;
    float rgf[3];

    rgf[0] = float(p.x * 7.0);
    rgf[1] = float(p.y * 7.0);
    rgf[2] = float(p.z * 7.0);
    sPink = g_sEqualWhite * 0.67 + g_sRed * 0.33;
    fDensity = SawToothWave(rgf[2] + 12.0*FractalNoise(rgf, 3, 5, 1.0f, 1.9f, 1));
    fDensity = fDensity*fDensity;
    ph->pm->color = sPink*(1.0 - fDensity);
}

void
Clouds(const Point3 &p, Hit *ph, RayTracer *scene)
{
    double fBlend;
    float rgf[3];

    rgf[0] = float(p.x * 7.0);
    rgf[1] = float(p.y * 7.0);
    rgf[2] = float(p.z * 7.0);
    fBlend = 2.0*FractalNoise(rgf, 3, 5, 1.0f, 1.9f, 0) - 1.0;
    if (fBlend < 0.0)
        fBlend = 0.0;
    else if (fBlend > 1.0)
        fBlend = 1.0;
    ph->pm->color = g_sEqualWhite*fBlend + g_sBlue*(1.0-fBlend);
}

static double
WoodWave(double x)
{
	x = (0.5/D_PI)*x + 0.25;
	x = 2.0*(x - floor(x));
	if (x < 1.0)                    // early wood
        return 0.0; 
    else if (x < 1.5) {             // late wood
        x = 2.0*(x - 1.0);
        return 0.5*x*x;
    } else {
        x = 2.0*(x - 2.0);
        return -0.5*x*x + 1.0;
    }
}

static float s_rgf[3];
double s_fBlend, s_fWarp;

void
Wood(const Point3 &p, Hit *ph, RayTracer *scene)
{
    double fBlend, fWarp;
    float rgf[3];

    rgf[0] = float(p.x * 3.0);
    rgf[1] = float(p.y * 3.0);
    rgf[2] = float(p.z * 3.0);
    fWarp = 3.0f*FractalNoise(rgf, 3, 2, 1.0f, 1.9f, 0);

    rgf[0] = float(p.x * 3.0);
    rgf[1] = float(p.y * 3.0);
    rgf[2] = float(p.z * 3.0);
    rgf[0] += float(fWarp);
    fBlend = 0.5 + 0.5*WoodWave(50.0*sqrt(rgf[0]*rgf[0] + rgf[1]*rgf[1]));
    fBlend = 0.9*fBlend + 0.1*SawToothWave(250.0*(1.0+fBlend)*sqrt(rgf[0]*rgf[0] + rgf[1]*rgf[1]));
    fBlend = 1.0 - fBlend;
    ph->pm->color = (g_sOrange + g_sRed) * (0.25*fBlend) + g_sEqualWhite * 0.03;
}

CoVector3
WoodFiber(Point3 const &p)
{
    double fWarp = 0.0;
    float rgf[3], rgfGrad[3];

    rgf[0] = float(p.x * 3.0);
    rgf[1] = float(p.y * 3.0);
    rgf[2] = float(p.z * 3.0);
    fWarp = 3.0f*FractalNoise(rgf, 3, 2, 1.0f, 1.9f, 0);

    rgf[0] = float(p.x * 3.0);
    rgf[1] = float(p.y * 3.0);
    rgf[2] = float(p.z * 3.0);
    rgf[0] += float(fWarp);
    rgf[0] *= 250.0f/3.0f;
    rgf[1] *= 250.0f/3.0f;
    rgf[2] *= 250.0f/3.0f;
    //rgf[0] = float(p.x * 250.0 );
    //rgf[1] = float(p.y * 250.0 );
    //rgf[2] = float(p.z * 250.0 );
    GradientBandNoise(rgfGrad, rgf, 2);
    return CoVector3(rgfGrad[0], rgfGrad[1], 0.0) * 1.0;
}

void
Spots(const Point3 &p, Hit *ph, RayTracer *scene)
{
    float rgf[3];
    double fBlend;

    rgf[0] = float(p.x * 7.0);
    rgf[1] = float(p.y * 7.0);
    rgf[2] = float(p.z * 7.0);
    fBlend = BandNoise(rgf, 3);
    if (fBlend > 0.5) {
        ph->pm->nLuster = MAT_COPPER;
    } else {
        ph->pm->nLuster = MAT_PLASTIC;
        ph->pm->color = g_sBlack;
    }
}

void 
Veins(const Point3 &p, Hit *ph, RayTracer *scene)
{
    double fBlend;
    float rgf[3];

    rgf[0] = float(p.x * 7.0);
    rgf[1] = float(p.y * 7.0);
    rgf[2] = float(p.z * 7.0);
    fBlend = FractalNoise(rgf, 3, 5, 1.0f, 1.9f, 0) - 0.5;
    if (fabs(fBlend) > 0.02)
        ph->pm->color = g_sGreen * 0.125;
    else
        ph->pm->color = g_sEqualWhite;
}
//
//  Spherical conformal image maps
//
#define MAPSELECT 32.70             // maximum scale is 1.188

struct Cartesian {
    double  xPolar;             // 0.5 x 1.0 cylindrical projection
    double  yPolar;
    double  xCylindrical;       // 0.5 x 0.5 azimuthal projection
    double  yCylindrical;
    double  fScalePolar;        // isotropic scaling of comformal mappings,
    double  fScaleCylindrical;  // zero if invalid (polar or cylindrical)
};

static void
SphericalMapping(Cartesian &cart, double fLatitude, double fLongitude)
{
    double x, y, z, r;
    static Cartesian cartLast;
    static double fLatLast, fLongLast;

    if (fLatitude == fLatLast && fLongitude == fLongLast) {
        cart = cartLast;
        return;
    }
    fLatLast = fLatitude;
    fLongLast = fLongitude;
    fLatitude  *= D_DTR;
    fLongitude *= D_DTR;
    z = sin(fLatitude);
    r = cos(fLatitude);
    x = cos(fLongitude) * r;                        // x axis, longitude == 0.0
    y = sin(fLongitude) * r;
    //
    //  Stereographic projection of north and south polar regions
    //
    if (fabs(fLatitude) < 0.5*D_PI - 0.001)
        cart.fScalePolar = (0.5*D_PI - fabs(fLatitude))/r;
    else
        cart.fScalePolar = 1/fabs(z);               // l'Hospital rule at latitude = pi/2
    if (fLatitude > 0.0) {
        cart.xPolar = -x/(1.0 + z)/D_PI + 0.25;
        cart.yPolar =  y/(1.0 + z)/D_PI + 0.25;
        if (cart.xPolar < 0.0 || cart.xPolar >= 0.5)
            cart.fScalePolar = 0.0;
    } else {
        cart.xPolar = -x/(1.0 - z)/D_PI + 0.75;
        cart.yPolar = -y/(1.0 - z)/D_PI + 0.25;
        if (cart.xPolar < 0.5 || cart.xPolar >= 1.0)
            cart.fScalePolar = 0.0;
    }
    if (cart.yPolar < 0.0 || cart.yPolar >= 0.5)
        cart.fScalePolar = 0.0;
    //
    //  Mercator projection of equatorial region
    //
    if (fLatitude <= 0.5*D_PI && fLatitude >= -0.5*D_PI) {
        cart.xCylindrical = fLongitude/(D_2PI);
        if (cart.xCylindrical >= 1.0)
            cart.xCylindrical -= 1.0;
        else if (cart.xCylindrical < 0.0)
            cart.xCylindrical += 1.0;
        cart.yCylindrical = -log(tan(0.25*D_PI + 0.5*fLatitude))/(D_2PI) + 0.75;
    }
    if (cart.yCylindrical >=0.5 && cart.yCylindrical < 1.0)
        cart.fScaleCylindrical = 1.0/r;
    else
        cart.fScaleCylindrical = 0.0;
    cartLast = cart;
}

DisplayRGB
SphericalSplatRGB(Image &im, DisplayRGB rgb, double fLatitude, double fLongitude,
                     Image *pimWeight)
{
    Cartesian cart;
    double w;

    im.m_bWrapY = WRAP_AROUND;
    w = double(im.m_nWidth);
    if (fLongitude < 0.0)
        fLongitude += 360.0;
    if (fLongitude >= 360.0)
        fLongitude -= 360.0;
    SphericalMapping(cart, fLatitude, fLongitude);
    if (fabs(fLatitude) > 16.0) {
        im.SplatRGB(rgb, w*cart.xPolar, w*cart.yPolar, pimWeight,
                        cart.fScalePolar);
    }
    if (fabs(fLatitude) < 65.0) {
        im.SplatRGB(rgb, w*cart.xCylindrical, w*cart.yCylindrical, pimWeight,
                        cart.fScaleCylindrical);
    }
    return rgb;
}

double
SphericalSplat(Image &im, double f, double fLatitude, double fLongitude,
                     Image *pimWeight)
{
    Cartesian cart;
    double w;

    im.m_bWrapY = WRAP_AROUND;
    w = double(im.m_nWidth);
    if (fLongitude < 0.0)
        fLongitude += 360.0;
    if (fLongitude >= 360.0)
        fLongitude -= 360.0;
    SphericalMapping(cart, fLatitude, fLongitude);
    if (fabs(fLatitude) > 16.0) {
        im.GaussianSplat(f, w*cart.xPolar, w*cart.yPolar, pimWeight,
                        cart.fScalePolar);
    }
    if (fabs(fLatitude) < 65.0) {
        im.GaussianSplat(f, w*cart.xCylindrical, w*cart.yCylindrical, pimWeight,
                        cart.fScaleCylindrical);
    }
    return f;
}

double
FastSphericalSplat(Image &im, double f, double fLatitude, double fLongitude,
                     Image *pimWeight)
{
    Cartesian cart;
    double w;

    im.m_bWrapY = WRAP_AROUND;
    w = double(im.m_nWidth);
    if (fLongitude < 0.0)
        fLongitude += 360.0;
    if (fLongitude >= 360.0)
        fLongitude -= 360.0;
    SphericalMapping(cart, fLatitude, fLongitude);
    if (fabs(fLatitude) > 16.0) {
        im.FastSplat(f, w*cart.xPolar, w*cart.yPolar, pimWeight,
                        0);
    }
    if (fabs(fLatitude) < 65.0) {
        im.FastSplat(f, w*cart.xCylindrical, w*cart.yCylindrical, pimWeight,
                        0);
    }
    return f;
}

DisplayRGB
SphericalSampleRGB(const Image &im,  double fLatitude, double fLongitude)
{
    Cartesian cart;
    double x, y, w;

    w = double(im.m_nWidth);
    if (fLongitude < 0.0)
        fLongitude += 360.0;
    if (fLongitude >= 360.0)
        fLongitude -= 360.0;
    SphericalMapping(cart, fLatitude, fLongitude);
    if (fabs(fLatitude) > MAPSELECT) {
        x = cart.xPolar;
        y = cart.yPolar;
    } else {
        x = cart.xCylindrical;
        y = cart.yCylindrical;
    }
    return im.SampleLerpRGB(w*x, w*y);
}

double
SphericalSample(const Image &im,  double fLatitude, double fLongitude)
{
    Cartesian cart;
    double x, y, w;

    w = double(im.m_nWidth);
    if (fLongitude < 0.0)
        fLongitude += 360.0;
    if (fLongitude >= 360.0)
        fLongitude -= 360.0;
    SphericalMapping(cart, fLatitude, fLongitude);
    if (fabs(fLatitude) > MAPSELECT) {
        x = cart.xPolar;
        y = cart.yPolar;
    } else {
        x = cart.xCylindrical;
        y = cart.yCylindrical;
    }
    return im.SampleLerp(w*x, w*y);
}

CoVector3
SampleNormal(const Image &im,  double fLatitude, double fLongitude, double fRadius)
{
    Cartesian cart;
    double x, y, z, w, xNormal, yNormal, zNormal, fKmPerPixel, fScale;
    CoVector3 cvN, cvN2, cvN3;
    float rgfGrad[2];
    int bPolar;

    w = double(im.m_nWidth);
    fKmPerPixel = D_2PI*fRadius/w;
    if (fLongitude < 0.0)
        fLongitude += 360.0;
    if (fLongitude >= 360.0)
        fLongitude -= 360.0;
    SphericalMapping(cart, fLatitude, fLongitude);
    if (fabs(fLatitude) > MAPSELECT) {
        x = cart.xPolar;
        y = cart.yPolar;
        fScale = cart.fScalePolar;
        bPolar = 1;
    } else {
        x = cart.xCylindrical;
        y = cart.yCylindrical;
        fScale = cart.fScaleCylindrical;
        bPolar = 0;
    }
    im.SampleGradient(rgfGrad, w*x, w*y, 0);
    rgfGrad[0] *= float(fScale/fKmPerPixel);
    rgfGrad[1] *= float(fScale/fKmPerPixel);
    cvN = CoVector3(rgfGrad[0], rgfGrad[1], 1.0);
    fLatitude *= D_DTR;
    fLongitude *= D_DTR;
    if (bPolar) {
        x = cos(-fLongitude)*cvN.x - sin(-fLongitude)*cvN.y;
        y = sin(-fLongitude)*cvN.x + cos(-fLongitude)*cvN.y;
        z = cvN.z;
        xNormal = cos(fLatitude)*z - sin(fLatitude)*x;
        yNormal = y;
        zNormal = sin(fLatitude)*z + cos(fLatitude)*x;
        x = cos(fLongitude)*xNormal - sin(fLongitude)*yNormal;
        y = sin(fLongitude)*xNormal + cos(fLongitude)*yNormal;
        z = zNormal;
    } else {
        y = cvN.x;
        x = cos(fLatitude)*cvN.z - sin(fLatitude)*cvN.y;
        z = sin(fLatitude)*cvN.z + cos(fLatitude)*cvN.y;
        zNormal = z;
        xNormal = cos(fLongitude)*x - sin(fLongitude)*y;
        yNormal = sin(fLongitude)*x + cos(fLongitude)*y;
        z = zNormal;
        x = xNormal;
        y = yNormal;
    }
    cvN3 = Normalize(CoVector3(x, y, z));
    /*
    z = sin(fLatitude);
    r = cos(fLatitude);
    x = cos(fLongitude) * r;                        // x axis, longitude == 0.0
    y = sin(fLongitude) * r;
    cvN2 = CoVector3(x, y, z);
    */
    return cvN3; // - cvN2;
}

#define SPHEREMAPWIDTH  (8*1024)

Image s_imSky;
Image s_imEarth, s_imEarthBump, s_imEarthWater;
Image s_imMars, s_imMarsBump, s_imMoon, s_imMoonBump, s_imVenus, s_imVenusBump, s_imMercury, s_imMercuryBump;
int   s_nSky, s_nEarth, s_nMars, s_nMoon, s_nVenus, s_nMercury;
Spectrum s_sOcean;

//
//  5.0857 sec.  - first version, xyz -> lat/long -> xyz
//  4.6957 sec.  - skip lat/lang -> xyz
//
void
PlanetMapping(const Point3 &p, Hit *ph, Image *pimColor, Image *pimBump, Image *pimWater)
{
    double fLat, fLong, w, x, y, z, r, fScale, fKmPerPixel, fSinLat, fCosLat, fSinLong, fCosLong;
    double xGrad, yGrad, zGrad, xNormal, yNormal, zNormal;
    CoVector3 cvN, cvPerturbation;
    struct Cartesian {
        double  xPolar;             // 0.5 x 1.0 cylindrical projection
        double  yPolar;
        double  xCylindrical;       // 0.5 x 0.5 azimuthal projection
        double  yCylindrical;
        double  fScalePolar;        // isotropic scaling of comformal mappings,
        double  fScaleCylindrical;  // zero if invalid (polar or cylindrical)
    } cart;
    float rgfGrad[2];
    int bPolar;

    x = p.x;
    y = p.y;
    z = p.z;
    cvN = Normalize(ph->cvNormal);
    x = cvN.x;
    y = cvN.y;
    z = cvN.z;
    fLat = asin(z);               // -pi/2 to pi/2
    fLong = D_PI - atan2(p.y, p.x);        // -pi to pi
    w = double(SPHEREMAPWIDTH);
    fKmPerPixel = D_2PI*6378.137/w;
    if (fLong < 0.0)
        fLong += 2.0*D_PI;
    fSinLat = z;
    fCosLat = sqrt(1.0 - z*z);
    r = sqrt(x*x + y*y);
    fSinLong = y/r;                 // r == 0.0 case
    fCosLong = x/r;

    if (fabs(fLat*180.0/D_PI) > MAPSELECT) {        // Only compute the polar or cylindrical mappings
        //
        //  Stereographic projection of north and south polar regions
        //
        if (fabs(fLat) < 0.5*D_PI - 0.001)
            cart.fScalePolar = (0.5*D_PI - fabs(fLat))/r;
        else
            cart.fScalePolar = 1/fabs(z);               // l'Hospital rule at latitude = pi/2
        if (fLat > 0.0) {
            cart.xPolar = +x/(1.0 + z)/D_PI + 0.25;
            cart.yPolar = +y/(1.0 + z)/D_PI + 0.25;
            if (cart.xPolar < 0.0 || cart.xPolar >= 0.5)
                cart.fScalePolar = 0.0;
        } else {
            cart.xPolar = +x/(1.0 - z)/D_PI + 0.75;
            cart.yPolar = -y/(1.0 - z)/D_PI + 0.25;
            if (cart.xPolar < 0.5 || cart.xPolar >= 1.0)
                cart.fScalePolar = 0.0;
        }
        if (cart.yPolar < 0.0 || cart.yPolar >= 0.5)
            cart.fScalePolar = 0.0;

        x = cart.xPolar;
        y = cart.yPolar;
        fScale = cart.fScalePolar;
        bPolar = 1;
    } else {
        //
        //  Mercator projection of equatorial region
        //
        if (fLat <= 0.5*D_PI && fLat >= -0.5*D_PI) {
            cart.xCylindrical = fLong/(D_2PI);
            if (cart.xCylindrical >= 1.0)
                cart.xCylindrical -= 1.0;
            else if (cart.xCylindrical < 0.0)
                cart.xCylindrical += 1.0;
            cart.yCylindrical = -log(tan(0.25*D_PI + 0.5*fLat))/(D_2PI) + 0.75;
        }
        if (cart.yCylindrical >=0.5 && cart.yCylindrical < 1.0)
            cart.fScaleCylindrical = 1.0/r;
        else
            cart.fScaleCylindrical = 0.0;

        x = cart.xCylindrical;
        y = cart.yCylindrical;
        fScale = cart.fScaleCylindrical;
        bPolar = 0;
    }
    if (pimWater && pimWater->SampleLerp(w*x, w*y) > 0.0) {
        ph->pm->color = s_sOcean;
        ph->ps->nRoughness = SURF_SHINY;  // SURF_WET;
        pimWater->SampleGradient(rgfGrad, w*x, w*y);
        xGrad = rgfGrad[0] * 256.0*fScale/fKmPerPixel;          
        yGrad = rgfGrad[1] * 256.0*fScale/fKmPerPixel;
        zGrad = 1.0;
    } else {
        ph->pm->color = RGBtoSpectrum(pimColor->SampleLerpRGB(w*x, w*y), 1);
        if (pimBump) {
            ph->ps->nRoughness = SURF_SHINY;
            pimBump->SampleGradient(rgfGrad, w*x, w*y);
            xGrad = rgfGrad[0] * 256.0*fScale/fKmPerPixel;
            yGrad = rgfGrad[1] * 256.0*fScale/fKmPerPixel;
          //  zGrad = 1.0;
              zGrad = 0.0;
        } else {
            ph->pm->nLuster = MAT_NOSHADE;      // sky dome textures
        }
    }
    //
    //  Perturbed surface normal
    //
    if (bPolar) {
        x = +fCosLong*xGrad - -fSinLong*yGrad;
        y = -fSinLong*xGrad + +fCosLong*yGrad;
        z = zGrad;
        xNormal = fCosLat*z - fSinLat*x;
        yNormal = y;
        zNormal = fSinLat*z + fCosLat*x;
        x = fCosLong*xNormal - fSinLong*yNormal;
        y = fSinLong*xNormal + fCosLong*yNormal;
        z = zNormal;
    } else {
        y = xGrad;
        x = fCosLat*zGrad - fSinLat*yGrad;
        z = fSinLat*zGrad + fCosLat*yGrad;
        zNormal = z;
        xNormal = fCosLong*x - fSinLong*y;
        yNormal = fSinLong*x + fCosLong*y;
        z = zNormal;
        x = xNormal;
        y = yNormal;
    }
    PerturbNormal(ph, CoVector3(x, y, z));
}

static void
LoadSkyMaps()
{
    if (s_nSky)
        return;
    VERBOSE printf("Loading Sky Maps\n");
    s_nSky = 1;
    s_imSky.ReadBMP("SkyColorMap.bmp");
    s_imSky.m_bWrapX = WRAP_AROUND;
    s_imSky.m_bWrapY = WRAP_MIRROR;
}

void
SkyTexture(const Point3 &p, Hit *ph, RayTracer *scene)
{
    if (!s_nSky) {
        LoadSkyMaps();
    }
    PlanetMapping(p, ph, &s_imSky, 0, 0);
}

static void
LoadEarthMaps()
{
    if (s_nEarth)
        return;
    VERBOSE printf("Loading Earth Maps\n");
    s_nEarth = 1;
    s_imEarth.ReadBMP("EarthColorMap.bmp");
    s_imEarth.m_bWrapX = WRAP_AROUND;
    s_imEarth.m_bWrapY = WRAP_MIRROR;
    s_imEarthBump.ReadBMP("EarthBump.bmp", 1.0);
    s_imEarthBump.m_bWrapX = WRAP_AROUND;
    s_imEarthBump.m_bWrapY = WRAP_MIRROR;
    s_imEarthWater.ReadBMP("EarthWater.bmp", 1.0);
    s_imEarthWater.m_bWrapX = WRAP_AROUND;
    s_imEarthWater.m_bWrapY = WRAP_MIRROR;
    s_sOcean = g_sEqualWhite * 0.20 + g_sCyan * 0.20 + g_sBlue * 0.30;  // Ocean color
    s_nEarth++;
}

void
EarthTexture(const Point3 &p, Hit *ph, RayTracer *scene)
{
    if (!s_nEarth) {
        LoadEarthMaps();
    }
    PlanetMapping(p, ph, &s_imEarth, &s_imEarthBump, &s_imEarthWater);
}

void
Earth2Texture(const Point3 &p, Hit *ph, RayTracer *scene)
{
    if (!s_nEarth) {
        LoadEarthMaps();
    }
    PlanetMapping(p, ph, &s_imEarth, &s_imEarthBump, 0);
}

static void
LoadMarsMaps()
{
    if (s_nMars)
        return;
    VERBOSE printf("Loading Mars Maps\n");
    s_nMars = 1;
    s_imMars.ReadBMP("MarsColorMap.bmp");
    s_imMars.m_bWrapX = WRAP_AROUND;
    s_imMars.m_bWrapY = WRAP_MIRROR;
    s_imMarsBump.ReadBMP("MarsBump.bmp", 1.0);
    s_imMarsBump.m_bWrapX = WRAP_AROUND;
    s_imMarsBump.m_bWrapY = WRAP_MIRROR;
    s_nMars++;
}

void
MarsTexture(const Point3 &p, Hit *ph, RayTracer *scene)
{
    if (!s_nMars) {
        LoadMarsMaps();
    }
    PlanetMapping(p, ph, &s_imMars, &s_imMarsBump, 0);
}

static void
LoadVenusMaps()
{
    if (s_nVenus)
        return;
    VERBOSE printf("Loading Venus Maps\n");
    s_nVenus = 1;
    s_imVenus.ReadBMP("VenusColorMap.bmp");
    s_imVenus.m_bWrapX = WRAP_AROUND;
    s_imVenus.m_bWrapY = WRAP_MIRROR;
    s_imVenusBump.ReadBMP("VenusBump.bmp", 1.0);
    s_imVenusBump.m_bWrapX = WRAP_AROUND;
    s_imVenusBump.m_bWrapY = WRAP_MIRROR;
    s_nVenus++;
}

void
VenusTexture(const Point3 &p, Hit *ph, RayTracer *scene)
{
    if (!s_nVenus) {
        LoadVenusMaps();
    }
    PlanetMapping(p, ph, &s_imVenus, &s_imVenusBump, 0);
}

static void
LoadMercuryMaps()
{
    if (s_nMercury)
        return;
    VERBOSE printf("Loading Mercury Maps\n");
    s_nMercury = 1;
    s_imMercury.ReadBMP("MercuryColorMap.bmp");
    s_imMercury.m_bWrapX = WRAP_AROUND;
    s_imMercury.m_bWrapY = WRAP_MIRROR;
    s_imMercuryBump.ReadBMP("MercuryBump.bmp", 1.0);
    s_imMercuryBump.m_bWrapX = WRAP_AROUND;
    s_imMercuryBump.m_bWrapY = WRAP_MIRROR;
    s_nMercury++;
}

void
MercuryTexture(const Point3 &p, Hit *ph, RayTracer *scene)
{
    if (!s_nMercury) {
        LoadMercuryMaps();
    }
    PlanetMapping(p, ph, &s_imMercury, &s_imMercuryBump, 0);
}

static void
LoadMoonMaps()
{
    if (s_nMoon)
        return;
    VERBOSE printf("Loading Moon Maps\n");
    s_nMoon = 1;
    s_imMoon.ReadBMP("MoonColorMap2.bmp");
    s_imMoon.m_bWrapX = WRAP_AROUND;
    s_imMoon.m_bWrapY = WRAP_MIRROR;
    s_imMoonBump.ReadBMP("MoonBump.bmp", 1.0);
    s_imMoonBump.m_bWrapX = WRAP_AROUND;
    s_imMoonBump.m_bWrapY = WRAP_MIRROR;
    s_nMoon++;
}

void
MoonTexture(const Point3 &p, Hit *ph, RayTracer *scene)
{
    if (!s_nMoon) {
        LoadMoonMaps();
    }
    PlanetMapping(p, ph, &s_imMoon, &s_imMoonBump, 0);
}

//
//  Generate the stereographic/mercator texture maps
//  Note: color maps use sRGB gamma, bump maps use gamma = 1.0, storyiing height linearly
//
#define MAPSTEP 1

int
MakeSphericalMapRGB(char *szLatLong, char *szSpherical, int nWide, double fLongOffset = 0.0)
{
    Image imInput, imOutput, imWeight;
    double fLatitude, fLongitude;
    DisplayRGB rgb;
    int i, j;

    if (imInput.ReadBMP(szLatLong) == 0) {
        printf("Failed to open %s\n", szLatLong);
        return 0;
    }
    imOutput.NewImage(nWide, nWide, 3, WRAP_AROUND, WRAP_CLAMP);
    imWeight.NewImage(nWide, nWide, 1, WRAP_AROUND, WRAP_CLAMP);
    imOutput.FillRGB(ML_Black);
    imWeight.FillRGB(0.0);
    for (j = 0; j < imInput.m_nHeight; j += MAPSTEP) {
        fLatitude = 90.0 - 180.0*INDEX_TO_SAMPLE(j)/double(imInput.m_nHeight);
        if (j%256 == 0) printf("row %d\n", j);
        for (i = 0; i < imInput.m_nWidth; i += MAPSTEP) {
            fLongitude = 360.0*INDEX_TO_SAMPLE(i)/double(imInput.m_nWidth) - 180.0;
            fLongitude = fmod(fLongitude + fLongOffset, 360.0);
            if (imInput.m_nChannels == 3)
                rgb = imInput.GetRGB(i, j);
            else
                rgb = DisplayRGB(imInput.Get(i, j));
            SphericalSplatRGB(imOutput, rgb, fLatitude, fLongitude, &imWeight);
        }
    }
    imOutput.Normalize(&imWeight);
    imOutput.WriteBMP(szSpherical);
    return 1;
}

#define BORDER_TOP      710
#define BORDER_BOTTOM   708
#define BORDER_RIGHT    300
#define BORDER_LEFT     301
#define SKY_WIDE        8600
#define SKY_HIGH        5416

int
MakeSphericalMapRGBfromElliptical(char *szElliptical, char *szSpherical, int nWide, double fLongOffset = 0.0)
{
    Image imInput, imOutput, imWeight;
    double fLatitude, fLongitude, x, y, theta, Rsqrt2;
    DisplayRGB rgb;
    int i, j;

    if (imInput.ReadBMP(szElliptical) == 0) {
        printf("Failed to open %s\n", szElliptical);
        return 0;
    }
    imOutput.NewImage(nWide, nWide, 3, WRAP_AROUND, WRAP_CLAMP);
    imWeight.NewImage(nWide, nWide, 1, WRAP_AROUND, WRAP_CLAMP);
    imOutput.FillRGB(ML_Black);
    imWeight.FillRGB(0.0);
    Rsqrt2 = double(SKY_HIGH - BORDER_TOP - BORDER_BOTTOM) / 2.0;
    for (j = 0; j < imInput.m_nHeight; j += MAPSTEP) {
        if (j < BORDER_TOP)
            continue;
        if (j >= SKY_HIGH - BORDER_BOTTOM)
            break;
        y = double(j-BORDER_TOP) - Rsqrt2 + 0.5;
        theta = asin(y/Rsqrt2);
        fLatitude = asin((2.0*theta + sin(2.0*theta))/D_PI);
        fLatitude = -fLatitude * 180.0/D_PI;
        if (fLatitude < -90.0 || fLatitude > 90.0)
            continue;
        if (j%256 == 0) printf("row %d\n", j);
        for (i = 0; i < imInput.m_nWidth; i += MAPSTEP) {
            if (i < BORDER_LEFT)
                continue;
            if (i >= SKY_WIDE - BORDER_RIGHT)
                break;
            x = double(i - BORDER_LEFT) - 2.0*Rsqrt2 + 0.5;
            fLongitude = D_PI * x /(2.0*Rsqrt2*cos(theta));
            fLongitude = 180.0 + fLongitude * 180.0/D_PI;
            if (fLongitude < 0.0 || fLongitude > 360.0)
                continue;
            if (fLongitude > 360.0)
                fLongitude -= 360.0;
            if (imInput.m_nChannels == 3)
                rgb = imInput.GetRGB(i, j);
            else
                rgb = DisplayRGB(imInput.Get(i, j));
            SphericalSplatRGB(imOutput, rgb, fLatitude, fLongitude, &imWeight);
        }
    }
    imOutput.Normalize(&imWeight);
    imOutput.WriteBMP(szSpherical);
    return 1;
}

int
MakeSphericalMap(char *szLatLong, char *szSpherical, int nWide, double fLongOffset = 0.0, int nEquiRectangular = 0)
{
    Image imInput, imOutput, imWeight;
    double fLatitude, fLongitude, f, y;
    int i, j;

    if (imInput.ReadBMP(szLatLong, 1.0) == 0) {         // Gamma 1.0
        printf("Failed to open %s\n", szLatLong);
        return 0;
    }
    imOutput.NewImage(nWide, nWide, 1, WRAP_AROUND, WRAP_CLAMP);
    imWeight.NewImage(nWide, nWide, 1, WRAP_AROUND, WRAP_CLAMP);
    imOutput.Fill(0.0);
    imWeight.Fill(0.0);
    for (j = 0; j < imInput.m_nHeight; j += MAPSTEP) {
        if (nEquiRectangular) {
            y = 1.0 - 2.0*INDEX_TO_SAMPLE(j)/double(imInput.m_nHeight);
            fLatitude = 90.0 - 180.0*acos(y)/D_PI;
        } else {
            fLatitude = 90.0 - 180.0*INDEX_TO_SAMPLE(j)/double(imInput.m_nHeight);
        }
        if (j%256 == 0) printf("row %d\n", j);
        for (i = 0; i < imInput.m_nWidth; i += MAPSTEP) {
            fLongitude = 360.0*INDEX_TO_SAMPLE(i)/double(imInput.m_nWidth) - 180.0;
            fLongitude = fmod(fLongitude + fLongOffset, 360.0);
            f = imInput.Get(i, j);
            SphericalSplat(imOutput, f, fLatitude, fLongitude, &imWeight);
        }
    }
    imOutput.Normalize(&imWeight);
    imOutput.WriteBMP(szSpherical, 1.0);
    return 1;
}

void
BuildSphericalMaps()
{
   // MakeSphericalMapRGB("EarthColor.bmp", "EarthColorMap.bmp", 4096*2);
   // MakeSphericalMap("EarthDepth.bmp", "EarthBump.bmp", 4096*2);
   // MakeSphericalMap("world.bathy.21601x10801.bmp", "EarthWater.bmp", 4096*2);

   // MakeSphericalMapRGB("MoonColor.bmp", "MoonColorMap.bmp", 4096*2);
   // MakeSphericalMap("MoonDepth.bmp", "MoonBump.bmp", 4096*2, 180.0);

   // MakeSphericalMapRGB("Mars.bmp", "MarsColorMap.bmp", 4096*2);
   // MakeSphericalMap("marsbump6k.bmp", "MarsBump.bmp", 4096*2);
   // MakeSphericalMapRGB("Magellan_LRS.bmp", "VenusColorMap.bmp", 4096*2);
   // MakeSphericalMapRGB("venus00n180_Composite.bmp", "VenusColorMap.bmp", 4096*2, 180.0);
   // MakeSphericalMap("Mercury_Messenger_DEM_Global_665m_v2_max.bmp", "MercuryBump.bmp", 4096*2, 0.0, 0);
   // MakeSphericalMapRGB("Mercury_MESSENGER_MDIS_Basemap_LOI_Mosaic_Global_32ppd.bmp", "MercuryColorMap.bmp", 4096*2, 180.0);
   // MakeSphericalMapRGB("Lunar_Clementine_12km.bmp", "MoonColorMap2.bmp", 4096*2);
   // MakeSphericalMapRGB("A1_stitch.bmp", "SkyColorMap1.bmp", 4096*2);
   MakeSphericalMapRGBfromElliptical("eso1908e.bmp", "SkyColorMap.bmp", 4096*2);
}

//
//  Build Moon color map from spectral images
//  nm = 360, 415, 566, 604, 643, 689
//  E:\Pictures\Maps\Moon\MoonWavelengths   \  WAC_HAPKE_360NM_E350N0450.bmp
//  6840 x 5321 panels  - 27,360 x 10,642 = 76*360 x 76*140  (-70 to +70 latitude)
//  Divide by solar spectrum in RGB conversion
//
static int s_rgnm[6] = {360, 415, 566, 604, 643, 689};
Image rg_imMoon[6];

static int
LoadMoonImage(Image &im, int nm, char cEast, int nNorth)
{
    char szFile[128];
    char szDir[64] = "E:\\Pictures\\Maps\\Moon\\MoonWavelengths\\";

    sprintf_s(szFile, sizeof(szFile), "%sWAC_HAPKE_%03dNM_E350%c%03d0.bmp", szDir, nm, cEast, nNorth);
    printf("|%s|\n", szFile);
    return im.ReadBMP(szFile);

}

#define LERP(v1, v2, vS, f1, f2) float(f1 + (f2 - f1)*(vS - v1)/(v2 - v1))

static Spectrum
MoonColor(double f360, double f415, double f566, double f604, double f643, double f689)
{
    Spectrum s;

    s.rgf[0] = LERP(360.0, 415.0, 380.0, f360, f415);       // 380 nm
    s.rgf[1] = LERP(415.0, 566.0, 420.0, f415, f566);       // 420 nm
    s.rgf[2] = LERP(415.0, 566.0, 460.0, f415, f566);       // 460 nm
    s.rgf[3] = LERP(415.0, 566.0, 500.0, f415, f566);       // 500 nm
    s.rgf[4] = LERP(415.0, 566.0, 540.0, f415, f566);       // 540 nm
    s.rgf[5] = LERP(566.0, 604.0, 580.0, f566, f604);       // 580 nm
    s.rgf[6] = LERP(604.0, 643.0, 620.0, f604, f643);       // 620 nm
    s.rgf[7] = LERP(643.0, 689.0, 660.0, f643, f689);       // 660 nm
    s.rgf[8] = LERP(643.0, 689.0, 700.0, f643, f689);       // 700 nm
    s.rgf[9] = LERP(643.0, 689.0, 740.0, f643, f689);       // 740 nm
    return s;
}

void
BuildMoonColorMap()
{
    Image im, imW, imMoon;
    Spectrum s;
    int nm, inm, iEast, cEast, nNorth, i, j;
    double fLat, fLong, f;
    double f360,  f415,  f566,  f604,  f643,  f689;
    char szFile[32];

    imW.NewImage(8*1024, 8*1024);
    for (inm = 0; inm < 6; inm++) {
        nm = s_rgnm[inm];
        rg_imMoon[inm].NewImage(8*1024, 8*1024);
        rg_imMoon[inm].Fill(0.0);
        imW.Fill(0.0);
        for (nNorth = 45; nNorth <= 315; nNorth += 90) {
            for (iEast = 0; iEast < 2; iEast++) {
                cEast = "NS"[iEast];
                LoadMoonImage(im, nm, cEast, nNorth);
                for (j = 0; j < 5321; j++) {
                    fLat = -70.0*double(cEast == 'S') + double(5320 - j)/76.0;
                    for (i = 0; i < 6840; i++) {
                        fLong = double(i)/76.0 + nNorth - 45;
                        f = im.Get(i, j);
                        FastSphericalSplat(rg_imMoon[inm], f, fLat, fLong, &imW);
                    }
                }
            }
        }
        sprintf_s(szFile, sizeof(szFile), "Moon%dnm.bmp", nm);
        rg_imMoon[inm].Normalize(&imW);
        rg_imMoon[inm].WriteBMP(szFile);
    }
    printf("Build color Moon map\n");
    imMoon.NewImage(8*1024, 8*1024, 3);
    for (j = 0; j < 8*1024; j++) {
        for (i = 0; i < 8*1024; i++) {
            f360 = rg_imMoon[0].Get(i, j);
            f415 = rg_imMoon[1].Get(i, j);
            f566 = rg_imMoon[2].Get(i, j);
            f604 = rg_imMoon[3].Get(i, j);
            f643 = rg_imMoon[4].Get(i, j);
            f689 = rg_imMoon[5].Get(i, j);
            s = MoonColor(f360, f415, f566, f604, f643, f689);
            imMoon.SetRGB(s.sRGB(), i, j);
        }
    }
    imMoon.WriteBMP("MoonColorMap_LRO.bmp");
}
//
//  Venus altimetry data from Pioneer Venus, Venera-15/16 and Magellan
//
//  Venera-15:  7 to 17 degree SAR beam angle
//  Magellan:   15 to 45 degree beam angle
//
struct AltimetryRecord {
    float   fLatitude;      // 90 to -90 degrees
    float   fLongitude;     // 0 to 360 degrees
    float   fRadius;        // R - 6000 kilometers
    float   fRadiusError;
    float   fReflectivity;  // Fresnel reflectivity
    int     nOrbit;
};

static AltimetryRecord s_rgar[512];
static double fmin = +1000000000.0;
static double fmax = -1000000000.0;
static double fminMagellan, fmaxMagellan, fminVenera15, fmaxVenera15, fminPioneer12, fmaxPioneer12;
static double fM, fB;

static Image s_imMagellan,  s_imMagellanWeight;
static Image s_imVenera15,  s_imVenera15Weight;
static Image s_imPioneer12, s_imPioneer12Weight;
static Image s_imAllVenus, s_imAllVenusWeight;
static Image s_imGradient, s_imVenera15_SAR, s_imVenera15map;
static Image s_imApprox;

#define VENUSRADIUS   6051.84  
#define VENUSHIGH       12.0
#define VENUSLOW        -3.0
#define VENUS_SNOW       2.6
#define VENUS_BUMP_LO   48.840431 
#define VENUS_BUMP_HI   63.789631
#define VENUS_BUMP_ALT(F) (VENUS_BUMP_LO + (F)*(VENUS_BUMP_HI - VENUS_BUMP_LO) + 6000.0 - VENUSRADIUS)

static void
BuildVenera15()
{
    int i, j, nWide, n8K;
    double x, y, fLong, fLat, f;
    Image imW;

    n8K = 8*1024;
    printf("BuildVenera15\n");
    nWide = s_imVenera15_SAR.m_nWidth;
    s_imVenera15map.NewImage(n8K, n8K, 1, WRAP_AROUND, WRAP_CLAMP);
    imW.NewImage(n8K, n8K, 1, WRAP_AROUND, WRAP_CLAMP);
    for (j = 0; j < nWide; j++) {
        for (i = 0; i < nWide; i++) {
            x = double(nWide - i - 1) + 0.5;
            y = double(nWide - j - 1) + 0.5;
            fLong = D_RTD*atan2(x - nWide/2, y - nWide/2); // + 180.0;
            fLat = hypot(x - nWide/2, y - nWide/2);
            fLat = 90.0 - 140.0*fLat/(nWide);             // Fit 20 to 90 degrees in square
            f = s_imVenera15_SAR.Sample(x, y);
            FastSphericalSplat(s_imVenera15map, f, fLat, fLong, &imW);
        }
    }
    s_imVenera15map.Normalize(&imW);
    s_imVenera15map.WriteBMP("Venera15map.bmp", 1.0);
}

static void
TestLSE()
{
    double fLong, fLat, x, y, z, a, b, c;
    Matrix3 m(0.0);
    Vector3 vB(0.0), vX;

    a = 1.7;
    b = 3.5;
    c = -0.4;
    printf("Test LSE1 %f %f %f\n", a, b, c);
    for (fLong = 0.0; fLong < 360.0; fLong += 1.0) {
        for (fLat = -90.0; fLat <= 90.0; fLat += 1.0) {
            x = RandomDouble();
            y = RandomDouble();
            z = a*x + b*y + c + 0.1*(RandomDouble() - 0.5);
            m.m[0][0] += x*x;  m.m[0][1] += x*y;  m.m[0][2] += x;
            m.m[1][0] += x*y;  m.m[1][1] += y*y;  m.m[1][2] += y;
            m.m[2][0] += x;    m.m[2][1] += y;    m.m[2][2] += 1.0;
            vB.rgf[0] += x*z;
            vB.rgf[1] += y*z;
            vB.rgf[2] += z;
        }
    }
    vX = InverseMatrix(m) * vB;
    a = vX.x;
    b = vX.y;
    c = vX.z;
    printf("Test LSE2 %f %f %f\n", a, b, c);
}

static void
HillShade()
{
    int i, j;
    CoVector3 cvN;
    Vector3 vL, vB(0.0), vX;
    double fLambert, a, b, c, x, y, z, fLat, fLong;
    float rgf[2];
    Matrix3 m(0.0);

    printf("HillShade\n");
    //
    //  Load height map
    //
    LoadVenusMaps();
    //
    //  Build a height-gradient map
    //
    s_imGradient.NewImage(SPHEREMAPWIDTH, SPHEREMAPWIDTH, 1, WRAP_AROUND, WRAP_CLAMP);
    s_imGradient.Fill(0.0);
    vL = Normalize(Vector3(0.0, 0.0, 1.0));
    a = 50.0;
    for (j = 0; j < SPHEREMAPWIDTH; j++) {
        for (i = 0; i < SPHEREMAPWIDTH; i++) {
            if (s_imVenusBump.Get(i, j) > 0.0) {
                s_imVenusBump.SampleGradient(rgf, double(i)+0.5, double(j)+0.5);
                cvN = Normalize(CoVector3(a*rgf[0], a*rgf[1], 1.0));
                fLambert = cvN * vL;
                // if (fLambert < 0.0)
                //     fLambert = 0.0;
                s_imGradient.Set(1.55 - 1.5*fLambert, i, j);
            }
        }
    }
    s_imGradient.WriteBMP("HillShadeVenus.bmp", 1.0);
    //
    //  Load Venera-15 SAR image
    //
    // s_imVenera15_SAR.ReadBMP("Inverse.bmp", 1.0);
    // BuildVenera15();
    s_imVenera15map.ReadBMP("Venera15map.bmp", 1.0);
    //
    //  Least-square-error plane fit - Magellan = a*Depth + b*Gradient + c;
    //
    for (fLong = 0.0; fLong < 360.0; fLong += 1.0) {
        for (fLat = -90.0; fLat <= 90.0; fLat += 1.0) {
            z = SphericalSample(s_imVenus, fLat, fLong);
            x = SphericalSample(s_imVenusBump, fLat, fLong);
            y = SphericalSample(s_imGradient, fLat, fLong);
            if (z == 0.0 || x == 0.0 || y == 0.0)
                continue;
            m.m[0][0] += x*x;  m.m[0][1] += x*y;  m.m[0][2] += x;
            m.m[1][0] += x*y;  m.m[1][1] += y*y;  m.m[1][2] += y;
            m.m[2][0] += x;    m.m[2][1] += y;    m.m[2][2] += 1.0;
            vB.rgf[0] += x*z;
            vB.rgf[1] += y*z;
            vB.rgf[2] += z;
        }
    }
    vX = InverseMatrix(m) * vB;
    printf("LSE vector: ");
    vX.Print();
    a = vX.x;
    b = vX.y;
    c = vX.z;
    s_imApprox.NewImage(SPHEREMAPWIDTH, SPHEREMAPWIDTH, 1, WRAP_AROUND, WRAP_CLAMP);
    for (j = 0; j < SPHEREMAPWIDTH; j++) {
        for (i = 0; i < SPHEREMAPWIDTH; i++) {
            x = s_imVenusBump.Get(i, j);
            y = s_imGradient.Get(i, j);
            z = a*x + b*y + c;
            s_imApprox.Set(z, i, j);
        }
    }
    s_imApprox.WriteBMP("VenusApprox.bmp");
    for (j = 0; j < SPHEREMAPWIDTH; j++) {
        for (i = 0; i < SPHEREMAPWIDTH; i++) {
            x = s_imApprox.Get(i, j);
            y = s_imVenera15map.Get(i, j);
            if (x == 0.0 || y == 0.0)
                continue;
            z = x + y - 0.5;
            s_imApprox.Set(z, i, j);
        }
    }
    s_imApprox.WriteBMP("ApproxV15.bmp");
}



static void
RenderAltimetry(const char *szAltFile, Image &imOutput, Image &imWeight, int nCheck = 0, double A = 1.0, double B = 0.0)
{
    HANDLE hFile;
    DWORD nBytesRead;
    double f, fMagellan;
    double x, y, xS, yS, xxS, xyS, fN, fD;
    int i, n;

    printf("Rendering %s\n", szAltFile);
    hFile = CreateFile(szAltFile, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING,
                                   FILE_ATTRIBUTE_NORMAL | FILE_FLAG_SEQUENTIAL_SCAN,
                                   NULL);
    fmin = +1000000000.0;
    fmax = -1000000000.0;
    xS = yS = xxS = xyS = fN = 0.0;
    for (;;) {
        ReadFile(hFile, &s_rgar, sizeof(s_rgar), &nBytesRead, NULL);
        n = nBytesRead/sizeof(AltimetryRecord);
        if (n == 0)
            break;
        for (i = 0; i < n; i++) {
            f = s_rgar[i].fRadius;
            f = (f - B)/A;
            if (f < VENUSRADIUS + VENUSLOW - 6000.0)
                continue;
            if (f > VENUSRADIUS + VENUSHIGH - 6000.0)
                continue;
            if (f < fmin) fmin = f;
            if (f > fmax) fmax = f;
            if (nCheck) {
                fMagellan = SphericalSample(s_imMagellan,  s_rgar[i].fLatitude, s_rgar[i].fLongitude);
                fMagellan = (fmaxMagellan - fminMagellan)*fMagellan + fminMagellan;
                if (fMagellan > fminMagellan && fMagellan < fmaxMagellan) {
                    y = f;
                    x = fMagellan;
                    fN += 1.0;
                    xS += x;
                    yS += y;
                    xxS += x*x;
                    xyS += x*y;
                }
            }
            if (nCheck == 2 && fMagellan > fminMagellan+EPSILON)
                continue;
            SphericalSplat(imOutput, f, s_rgar[i].fLatitude, s_rgar[i].fLongitude, &imWeight);
           // SphericalSplatGaussian(s_rgar[i].fRadius + 6000.0, 1.0,
           //                        s_rgar[i].fLatitude, s_rgar[i].fLongitude,
           //                        10.0);
        }
    }
    if (nCheck) {
        fD = fN*xxS - xS*xS;
        fB = (xxS*yS - xS*xyS)/fD;
        fM = (fN*xyS - xS*yS)/fD;
        printf(" %s = %f * Magellan + %f\n", szAltFile, fM, fB);
    }
    CloseHandle(hFile);
}

static void
LSE1_Adjust(Image &imAdjusted, const Image &im, const Image &imStandard)
{
    double fLat, fLong;
    double x, y, xS, yS, xxS, xyS, fN, fD, fM, fB;
    int i, j;

    xS = yS = xxS = xyS = fN = 0.0;
    for (fLong = 0.0; fLong < 360.0; fLong += 0.5) {
        for (fLat = -90.0; fLat <= 90.0; fLat += 0.5) {
            y = SphericalSample(imStandard, fLat, fLong);
            x = SphericalSample(im, fLat, fLong);
            if (x <= 0.0 || y <= 0.0)
                continue;
            fN += 1.0;
            xS += x;
            yS += y;
            xxS += x*x;
            xyS += x*y;
        }
    }
    if (&im != &imAdjusted)
        imAdjusted.NewImage(im);
    fD = fN*xxS - xS*xS;
    fB = (xxS*yS - xS*xyS)/fD;
    fM = (fN*xyS - xS*yS)/fD;
    printf("LSE1(%d):  %f x + %f\n", int(fN), fM, fB);
    for (j = 0; j < im.m_nHeight; j++) {
        for (i = 0; i < im.m_nWidth; i++) {
            x = im.Get(i, j);
            if (x == 0.0)
                continue;
            y = fM*x + fB;
            imAdjusted.Set(y, i, j);
        }
    }
}

static void
LSE2_Adjust(Image &imAdjusted, const Image &im1, const Image &im2, const Image &imStandard)
{
    double fLong, fLat, x, y, z, a, b, c;
    Matrix3 m(0.0);
    Vector3 vB(0.0), vX;
    int i, j;

    for (fLong = 0.0; fLong < 360.0; fLong += 0.5) {
        for (fLat = -90.0; fLat <= 90.0; fLat += 0.5) {
            z = SphericalSample(imStandard, fLat, fLong);
            x = SphericalSample(im1, fLat, fLong);
            y = SphericalSample(im2, fLat, fLong);
            if (x == 0.0 || y == 0.0 || z == 0.0)
                continue;
            m.m[0][0] += x*x;  m.m[0][1] += x*y;  m.m[0][2] += x;
            m.m[1][0] += x*y;  m.m[1][1] += y*y;  m.m[1][2] += y;
            m.m[2][0] += x;    m.m[2][1] += y;    m.m[2][2] += 1.0;
            vB.rgf[0] += x*z;
            vB.rgf[1] += y*z;
            vB.rgf[2] += z;
        }
    }
    imAdjusted.NewImage(im1);
    vX = InverseMatrix(m) * vB;
    a = vX.x;
    b = vX.y;
    c = vX.z;
    for (j = 0; j < im1.m_nHeight; j++) {
        for (i = 0; i < im1.m_nWidth; i++) {
            x = im1.Get(i, j);
            y = im2.Get(i, j);
            if (x == 0.0 || y == 0.0)
                continue;
            z = a*x + b*y + c;
            imAdjusted.Set(z, i, j);
        }
    }
}

static void
ApproximateMagellan(Image &imApprox)
{
    Image imVenera15;
    Vector3 vL;
    CoVector3 cvN;
    double fBeamAngle, f, fH, fAlt, fVenera, a;
    float rgf[2];
    int i, j, w;

    printf("Approx image\n");
    w = s_imVenusBump.m_nWidth;
    imApprox.NewImage(w, w);
    fBeamAngle = 30.0 * D_PI/180.0;
    vL = Normalize(Vector3(-sin(fBeamAngle), 0.0, cos(fBeamAngle)));
    a = 50.0;
    for (j = 0; j < w; j++) {
        for (i = 0; i < w; i++) {
            fAlt = VENUS_BUMP_ALT(fH = s_imVenusBump.SampleLerp(double(i)+0.5, double(j)+0.5));
            if (fH == 0.0)
                continue;
            s_imVenusBump.SampleGradient(rgf, double(i)+0.5, double(j)+0.5);
            cvN = Normalize(CoVector3(a*rgf[0], a*rgf[1], 1.0));
            f = 1.0 - cvN * vL;
            if (fAlt > VENUS_SNOW || fAlt < 0.0)
                f += fAlt/20.0;
            imApprox.Set(f, i, j);
        }
    }
    imApprox.WriteBMP("Approx1_NL.bmp");
   // LSE1_Adjust(imApprox, imApprox, s_imVenus);
   // imApprox.WriteBMP("Approx2_LSE.bmp");
   imVenera15.ReadBMP("Venera15map.bmp");
    for (j = 0; j < w; j++) {
        for (i = 0; i < w; i++) {
            fVenera = imVenera15.Get(i, j);
            if (fVenera == 0)
                continue;
            f = imApprox.Get(i, j);
            f = f + 0.5*(fVenera - 0.25);
            imApprox.Set(f, i, j);
        }
    }
    imApprox.WriteBMP("Approx3_V15.bmp");
}

void
BuildVenusMaps()
{
    Image imOutput, imWeight;
    int i, j, nWide;
    double f;
    
    printf("BuildVenusMaps:\n");
    LoadVenusMaps();
    ApproximateMagellan(imOutput);
    return;

    TestLSE();
    HillShade();
    return;
    
    //
    //  Build height map from Magellan, Venera-15 and Pioneer-12 altimeter data
    //
    nWide = 4*1024;
    s_imAllVenus.NewImage(nWide, nWide, 1, WRAP_AROUND, WRAP_CLAMP);
    s_imAllVenusWeight.NewImage(nWide, nWide, 1, WRAP_AROUND, WRAP_CLAMP);
    s_imAllVenus.Fill(0.0);
    s_imAllVenusWeight.Fill(0.0);
    //
    //  Magellan data
    //
    s_imMagellan.NewImage(nWide, nWide, 1, WRAP_AROUND, WRAP_CLAMP);
    s_imMagellanWeight.NewImage(nWide, nWide, 1, WRAP_AROUND, WRAP_CLAMP);
    s_imMagellan.Fill(0.0);
    s_imMagellanWeight.Fill(0.0);
    RenderAltimetry("Magellan.alt", s_imMagellan, s_imMagellanWeight);
    fminMagellan = fmin;
    fmaxMagellan = fmax;
    printf("Magellan altimetry: %f - %f\n", fmin, fmax);
    s_imMagellan.Normalize(&s_imMagellanWeight);
    RenderAltimetry("Magellan.alt", s_imAllVenus, s_imAllVenusWeight);
    for (i = 0; i < nWide; i++) {
        for (j = 0; j < nWide; j++) {
            f = s_imMagellan.Get(i, j);
            s_imMagellan.Set((f - fmin)/(fmax - fmin), i, j);
        }
    }
    s_imMagellan.WriteBMP("VenusMagellan.bmp", 1.0);
    //
    //  Venear-15 data
    //
    s_imVenera15.NewImage(nWide, nWide, 1, WRAP_AROUND, WRAP_CLAMP);
    s_imVenera15Weight.NewImage(nWide, nWide, 1, WRAP_AROUND, WRAP_CLAMP);
    s_imVenera15.Fill(0.0);
    s_imVenera15Weight.Fill(0.0);
    RenderAltimetry("Venera.alt", s_imVenera15, s_imVenera15Weight, 1);
    fminVenera15 = fmin;
    fmaxVenera15 = fmax;
    printf("Venera15 altimetry: %f - %f\n", fmin, fmax);
    s_imVenera15.Normalize(&s_imVenera15Weight);
    for (i = 0; i < nWide; i++) {
        for (j = 0; j < nWide; j++) {
            f = s_imVenera15.Get(i, j);
            s_imVenera15.Set((f - fmin)/(fmax - fmin), i, j);
        }
    }
    s_imVenera15.WriteBMP("VenusVenera15.bmp", 1.0);

    s_imVenera15.NewImage(nWide, nWide, 1, WRAP_AROUND, WRAP_CLAMP);
    s_imVenera15Weight.NewImage(nWide, nWide, 1, WRAP_AROUND, WRAP_CLAMP);
    s_imVenera15.Fill(0.0);
    s_imVenera15Weight.Fill(0.0);
    RenderAltimetry("Venera.alt", s_imAllVenus, s_imAllVenusWeight, 2, fM, fB);     // corrected Venera data into AllVenus
    //
    //  Pioneer Data
    //
    s_imPioneer12.NewImage(nWide, nWide, 1, WRAP_AROUND, WRAP_CLAMP);
    s_imPioneer12Weight.NewImage(nWide, nWide, 1, WRAP_AROUND, WRAP_CLAMP);
    s_imPioneer12.Fill(0.0);
    s_imPioneer12Weight.Fill(0.0);
    RenderAltimetry("Pioneer.alt", s_imPioneer12, s_imPioneer12Weight, 1);
    fminPioneer12 = fmin;
    fmaxPioneer12 = fmax;
    printf("Pioneer12 altimetry: %f - %f\n", fmin, fmax);
    s_imPioneer12.Normalize(&s_imPioneer12Weight);
    for (i = 0; i < nWide; i++) {
        for (j = 0; j < nWide; j++) {
            f = s_imPioneer12.Get(i, j);
            s_imPioneer12.Set((f - fmin)/(fmax - fmin), i, j);
        }
    }
    s_imPioneer12.WriteBMP("VenusPioneer12.bmp", 1.0);

    s_imPioneer12.NewImage(nWide, nWide, 1, WRAP_AROUND, WRAP_CLAMP);
    s_imPioneer12Weight.NewImage(nWide, nWide, 1, WRAP_AROUND, WRAP_CLAMP);
    s_imPioneer12.Fill(0.0);
    s_imPioneer12Weight.Fill(0.0);
    RenderAltimetry("Pioneer.alt", s_imAllVenus, s_imAllVenusWeight, 2, fM, fB);     // corrected Venera data into AllVenus
    //
    //  Composite Data
    //
    fmin = Min(fminMagellan, Min(fminVenera15, fminPioneer12));
    fmax = Max(fmaxMagellan, Max(fmaxVenera15, fmaxPioneer12));
    s_imAllVenus.Normalize(&s_imAllVenusWeight);
    for (i = 0; i < nWide; i++) {
        for (j = 0; j < nWide; j++) {
            f = s_imAllVenus.Get(i, j);
            s_imAllVenus.Set((f - fmin)/(fmax - fmin), i, j);
        }
    }
    printf("Global min max of VenusAll: %f %f\n", fmin, fmax);
    // s_imAllVenus.WriteBMP("AllVenus.bmp", 1.0);



    return;
    imOutput.NewImage(nWide, nWide, 1, WRAP_AROUND, WRAP_CLAMP);
    imWeight.NewImage(nWide, nWide, 1, WRAP_AROUND, WRAP_CLAMP);
    imOutput.Fill(0.0);
    imWeight.Fill(0.0);

    RenderAltimetry("Pioneer.alt", imOutput, imWeight);
    RenderAltimetry("Venera.alt", imOutput, imWeight);
    RenderAltimetry("Magellan.alt", imOutput, imWeight);

    imOutput.Normalize(&imWeight);
    for (i = 0; i < nWide; i++) {
        for (j = 0; j < nWide; j++) {
            f = imOutput.Get(i, j);
            imOutput.Set((f - fmin)/(fmax - fmin), i, j);
            f = imWeight.Get(i, j);
            imWeight.Set((f)/(fmax - fmin), i, j);
        }
    }
    printf("fmin, fmax %f %f\n", fmin, fmax);
    imOutput.WriteBMP("VenusBump.bmp", 1.0);
}