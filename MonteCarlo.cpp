#include "stdafx.h"
#include "MonteCarlo.h"

#pragma intrinsic (sqrt, log, exp, sin, cos, pow)

#define D_PI            3.14159265358979323846264338327950288419716939937510
#define D_E             2.71828182845904523536028747135266249775724709369995

unsigned ML_nMarsagliaX  = 886459;
unsigned ML_nMarsagliaC  = 361290869;

//
//  Uniform distribution from 0 to N-1
//
unsigned
RandomLessThanN(unsigned nModulus)
{
    unsigned nRemainder, n;

    nRemainder = 0xFFFFFFFF % nModulus;
    do {
        n = RandomUnsigned();
    } while (n <= nRemainder);
    return n % nModulus;
}

//
//  True random numbers generated from jitter in clock interrupt
//
unsigned
TrueRandomUnsigned()
{
    static LARGE_INTEGER nPerformanceCount;
    static unsigned n;
    BOOL bResult;
    int i;


    for (i = 0; i < 32; i += 8) {
        bResult = QueryPerformanceCounter(&nPerformanceCount);
        Sleep(1L);
        n = (n << 8) | (n >> 24);
        n ^= nPerformanceCount.LowPart;
    }
    return n;
}
//
//  Non-uniform variates
//
double
RandomNormal()
{
    double s, x, y;
    static double fNormal;
    static int nSecondVariate = 1;

    nSecondVariate = !nSecondVariate;
    if (nSecondVariate) {
        return fNormal;
    } else {
        do {
            x = 1.0 - 2.0*RandomDouble();
            y = 1.0 - 2.0*RandomDouble();
            s = x*x + y*y;
        } while (s >= 1.0);
        s = sqrt(-2.0*log(s)/s);
        fNormal = x*s;
        return y*s;
    }
}
//
//  random point on disk
//
void
RandomInsideDisk(double rgf[2])
{
    do {
        rgf[0] = 2.0*RandomDouble() - 1.0;
        rgf[1] = 2.0*RandomDouble() - 1.0;
    } while (rgf[0]*rgf[0] + rgf[1]*rgf[1] >= 1.0);
}