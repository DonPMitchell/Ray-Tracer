#pragma once
//
//  MLpatterns - sampling patterns in the unit square
//  D.P. Mitchell  12/15/2003.
//
#pragma once
#include "Affine.h"

#define PHI     ((sqrt(5.0) + 1.0)/2.0)
#define PHI1    1.618033988749895
#define PHI2    1.324717957244746
#define PHI3    1.220744084605760
#define PHI4    1.167303978261419
#define PHI5    1.134724138401519

#define BLUE_NSQRT  40
#define BLUE_N      (BLUE_NSQRT*BLUE_NSQRT)

//
//  2D sampling-pattern tiles
//
extern void     GeneralHilbert(Vector2 rgv[], int nW, int nH);
extern int      UniformPattern(Vector2 rgv[], int nSqrt);
extern int      JitterPattern(Vector2 rgv[], int nSqrt, double fJitter = 1.0);
extern int      JitterHilbert(Vector2 rgv[], int nSqrt, double fJitter = 1.0);
extern int      PoissonDiskPattern(Vector2 rgv[], int nSqrt);
extern int      BlueNoisePattern(Vector2 rgv[], int nSqrt);
extern int      BlueJitterPattern(Vector2 rgv[], int nSqrt);
extern int      OptimizeLocality(Vector2 rgv[], int nSamples, int nWrap = 1);

extern int      WritePattern(char *szFileName, char *szPatternName, Vector2 rgv[], int nSqrt);
//
//  1D sampling patterns (e.g., for time samples)
//
extern int      UniformPattern(double rgf[], int nSamples);
extern int      JitterPattern(double rgf[], int nSamples);
//
//  Low arbitrary-edge discrepancy on unit radius disk
//
extern double DiamondAngle(double y, double x);

extern void BuildQuasiIndex(int rgn[], int nPrime, int N);
//
//  Fibonacci spiral on the unit sphere
//
extern void FibonacciSphere(double rgf[], int i, int N);
//
//  Equiarea maps
//
extern Vector2 SquareToUnitDisk(const Vector2 &pntI2);
extern Vector2 DiskToUnitSquare(const Vector2 &pntD);
extern Vector3 DiskToUnitHemisphere(const Vector2 &pntD);
extern Vector3 SquareToUnitSphere(const Vector2 &pntI2);
extern Vector2 HemisphereToUnitDisk(const Vector3 &pntH);
extern Vector2 SphereToUnitSquare(const Vector3 &pntS2);
extern Vector2 SchwarzChristoffel(const Vector2 &vSqr);         // conformal map square -> disk