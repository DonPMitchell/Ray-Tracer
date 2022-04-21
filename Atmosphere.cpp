//
//  Data describing Venus and its atmosphere
//
#include "stdafx.h"
#include "BookGraphics.h"
#include "MLspecialfunctions.h"
#include "MLvector.h"
//
//  Data every 5 km from 0 to 100 (Seiff)
//
double g_rgfPressure[21] = {        // P, bars
    92.10,     66.65,    47.39,     33.04,
    22.52,     14.93,     9.581,     5.917,
     3.501,     1.979,    1.066,     0.5314,
     0.2357,    0.09599,  0.03690,   0.01335,
     0.004476,  0.001351, 0.0003736, 0.00009814,
     0.00002660
};

double g_rgfDensity[21] = {         // rho, kg/m**3
    64.79,  49.87,  37.72,  27.95,
    20.39,  14.57,  10.15,  6.831,
    4.404,  2.693,  1.594,  0.9207,
    0.4694, 0.2055, 0.08393, 0.03236,
    0.01186, 0.003900, 0.001151, 0.0003040,
    0.00007890
};

double g_rgfTemperature[21] = {     // T, Kelvin
    735.3,  696.8,  658.2,  620.8,
    580.7,  539.2,  496.9,  455.5,
    417.6,  385.4,  350.5,  302.3,
    262.8,  243.2,  229.8,  215.4,
    197.1,  181.0,  169.4,  168.2,
    175.4
};

double g_rgfRefractiveIndex[21] = {
    1.014716, 1.011321, 1.008559, 1.006339,
    1.004623, 1.003303, 1.002301, 1.001548,
    1.000998, 1.000610, 1.000361, 1.000209,
    1.000106, 1.000047, 1.000019, 1.000007,
    1.000003, 1.000001, 1.000000, 1.000000,
    1.000000
};

double g_rgfRayleighScattering[21] = {
    9.095310e-023, 7.000820e-023, 5.295186e-023, 3.923660e-023,
    2.862376e-023, 2.045357e-023, 1.424871e-023, 9.589452e-024,
    6.182396e-024, 3.780471e-024, 2.237679e-024, 1.292491e-024,
    6.589502e-025, 2.884837e-025, 1.178221e-025, 4.542742e-026,
    1.664923e-026, 5.474874e-027, 1.615790e-027, 4.267594e-028,
    1.107609e-028
};

inline double
ML_GeometricLerp(double f1, double f2, double a)
{
    return f1 * pow(f2/f1, a);
}

double
VenusPressure(double fKm)
{
    int iKm;
    double dKm;

    if (fKm >= 100.0)
        return g_rgfPressure[20];
    fKm /= 5.0;
    iKm = int(floor(fKm));
    if (iKm < 0)
        iKm = 0;
    dKm = fKm - double(iKm);

    return ML_GeometricLerp(g_rgfPressure[iKm], g_rgfPressure[iKm+1], dKm);
}

double
VenusDensity(double fKm)
{
    int iKm;
    double dKm;

    if (fKm >= 100.0)
        return g_rgfDensity[20];
    fKm /= 5.0;
    iKm = int(floor(fKm));
    if (iKm < 0)
        iKm = 0;
    dKm = fKm - double(iKm);

    return ML_GeometricLerp(g_rgfDensity[iKm], g_rgfDensity[iKm+1], dKm);
}

double
VenusTemperature(double fKm)
{
    int iKm;
    double dKm;

    if (fKm >= 100.0)
        return g_rgfTemperature[20];
    fKm /= 5.0;
    iKm = int(floor(fKm));
    if (iKm < 0)
        iKm = 0;
    dKm = fKm - double(iKm);

    return ML_Lerp(g_rgfTemperature[iKm], g_rgfTemperature[iKm+1], dKm);
}

double
VenusRefractiveIndex(double fKm)
{
    int iKm;
    double dKm;

    if (fKm >= 100.0)
        return g_rgfRefractiveIndex[20];
    fKm /= 5.0;
    iKm = int(floor(fKm));
    if (iKm < 0)
        iKm = 0;
    dKm = fKm - double(iKm);

    return ML_Lerp(g_rgfRefractiveIndex[iKm], g_rgfRefractiveIndex[iKm+1], dKm);
}

double
RefractionAir(double fNanometers)
{
    fNanometers /= 1000.0;   // nanometers to micrometers
    return 1.0 + 0.0472326/(173.3 - 1.0/(fNanometers*fNanometers));
}

static double
ComputeVenusRefraction(double fKm)
{
    double fDensity;

    fDensity = 1000.0 * VenusDensity(fKm) / MOLE_WEIGHT_VENUS;  // molar density
    return sqrt((2.0*LORENTZ_LORENZ_VENUS*fDensity + 1.0)
                 / (1.0 - LORENTZ_LORENZ_VENUS*fDensity));
}

void
TabulateVenusRefraction()
{
    int iKm;
    double fKm, fRefract;

    for (iKm = 0; iKm <= 100; iKm += 5) {
        fKm = iKm;
        fRefract = ComputeVenusRefraction(fKm);
        printf("%f, ", fRefract);
        if (iKm % 20 == 15)
            printf("\n");
    }
    printf("\n");
}

#define RAYLEIGH (24.0*D_PI*D_PI*D_PI*KING_FACTOR_VENUS)

void
TabulateVenusRayleigh()
{
    int iKm;
    double fKm, fDensity, fN, fRefract, fRayleigh, fLL;

    for (iKm = 0; iKm <= 100; iKm += 5) {
        fKm = iKm;
        fDensity = 1000.0 * VenusDensity(fKm) / MOLE_WEIGHT_VENUS;  // moles/meter**3
        fN = C_AVOGADRO * fDensity / 1000000.0;                     // molecules/cm**3
        fRefract = ComputeVenusRefraction(fKm);
        fLL = (fRefract*fRefract - 1.0)/(fRefract*fRefract + 2.0);
        //REVIEW: applying Lorentz-Lorenz effect twice?
        fRayleigh = RAYLEIGH * fLL * fLL / (fN);                    // volume scattering coeff.
        printf("%0.7g, ", fRayleigh);
        if (iKm % 20 == 15)
            printf("\n");
    }
    printf("\n");
}