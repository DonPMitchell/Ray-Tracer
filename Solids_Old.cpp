#include "stdafx.h"
#include "RayTracer2020.h"
#include "Solids.h"
#pragma intrinsic(sqrt)

#define EPSILON 1.0e-8
#define DENSITY(PH) ((PH)->pm->fDensity)

static HitWork s_hw;

void
DeleteHitList(Hit *ph)
{
    if (ph) {
        DeleteHitList(ph->phNext);
        delete ph;
    }
}

inline Hit *
DeleteFirstHit(Hit *ph)
{
    Hit *phNext = 0;

    if (ph) {
        phNext = ph->phNext;
        delete ph;
    }
    return phNext;
}

static void
PrintHitList(char *sz, Hit *ph)
{
    printf("%s { ", sz);
    for (; ph; ph = ph->phNext)
        printf("t=%0.1f d=%0.1f, ", ph->t, ph->fDensity);
    printf("}\n");
}
//
//  Solve quadratic equation: c + x*(b + x*a) == 0
//
int
SolveLinear(double rgfRoots[], double a, double b)
{
    if (a == 0.0) {
        if (b == 0.0) {             // handle the zero-order equation b = 0
            rgfRoots[0] = 0.0;
            return 1;
        } else {
            return 0;
        }
    } else {
        rgfRoots[0] = -b/a;
        return 1;
    }
}

int
SolveQuadratic(double rgfRoots[], double a, double b, double c)
{
    double fDiscriminant;

    if (a == 0.0)
        return SolveLinear(rgfRoots, b, c);
    fDiscriminant = b*b - 4.0*a*c;
    if (fDiscriminant <= 0.0) {
        if (fDiscriminant == 0.0) {
            rgfRoots[0] = rgfRoots[1] = -b/(a + a);     // double root is counted twice
            return 2;
        } else
            return 0;
    }
    fDiscriminant = sqrt(fDiscriminant);
    if (a < 0.0) {
        a = -a;
        b = -b;
        c = -c;
    }
    if (b > 0.0) {
		rgfRoots[0] = (-b - fDiscriminant)/(a + a);     // smaller root first
		rgfRoots[1] = ( c + c)/(-b - fDiscriminant);    // avoid cancellation error

	} else {
		rgfRoots[1] = (-b + fDiscriminant)/(a + a);
		rgfRoots[0] = ( c + c)/(-b + fDiscriminant);    // smaller root first
	}
    return 2;
}

static Material s_mVacuum(0.0, g_sBlack, 0);
static Material s_mNeutronium(1000000.0, g_sEqualWhite, 0);

Hit*
HitWork::HitList()
{
    int i, iFirst;
    Hit *ph;

    //
    //  Ray is semi-infinite, so trim off negative roots and determine if ray
    //  origin is outside solid (ray is entering) or inside.
    //
    for (iFirst = 0; iFirst < nHits; iFirst++)
        if (rgfRoots[iFirst] >= 0.0)
            break;
    //
    //  Regularization: (0, eps, ...) -> (...) or (eps, ...) -> (0.0, eps)
    //
    if (!nRayEnteringSolid && nHits >= 2 && rgfRoots[0] == 0.0 && rgfRoots[1] < EPSILON)
        iFirst += 2;
    else if (nRayEnteringSolid && nHits && rgfRoots[0] < EPSILON)
        rgfRoots[0] = 0.0;
    //
    //  Build the hit list
    //
    ph = 0;
    for (i = nHits - 1; i >= iFirst; --i) {
        ph = new Hit(ph);
        ph->t = rgfRoots[i];
        ph->cvNormal = rgvNormals[i];
        if (nRayInsideAtInfinity) {
            ph->pm = &s_mNeutronium;
            ph->fDensity = DENSITY(ph);
        } else {
            ph->pm = &s_mVacuum;
            ph->fDensity = DENSITY(ph);
        }
        nRayInsideAtInfinity = !nRayInsideAtInfinity;
    }
    //
    //  If the origin is embedded inside the solid, add a node at t == 0.0
    //
    if (nRayInsideAtInfinity) {
        ph = new Hit(ph);
        ph->t = 0.0;
        ph->cvNormal = CoVector3(0.0, 0.0, 0.0);
        ph->pm = &s_mNeutronium;
        ph->fDensity = DENSITY(ph);
    }
    return ph;
}
//
//  Liang-Barsky clipping of parametric line against implicit planes
//
static int
ClipPlane(Ray &ray, const CoVector3 &cvNormal, HitWork &hw, int nDoClip = 1)
{
    double B, C, t;

    if (nDoClip && hw.nHits == 0)
        return 0;
    B = cvNormal * ray.vDir;
    C = cvNormal * (ray.pOrg - Point3(0.0, 0.0, 0.0)) - 1.0;
    if (B)
        t = -C/B;
    if (!nDoClip) {
        hw.nHits = 2;
        hw.rgfRoots[0] = -HUGE;
        hw.rgfRoots[1] = HUGE;
        hw.nRayInsideAtInfinity = 1;
    }
    hw.nRayInsideAtInfinity &= B < 0.0;
    if (B < 0.0) {
        if (hw.nHits == 2 && hw.rgfRoots[1] <= t) {
            hw.nHits = 0;
        } else if (hw.rgfRoots[0] < t) {
            hw.rgfRoots[0] = t;
            hw.rgvNormals[0] = cvNormal;
        }
    } else if (B > 0.0) {
        if (t <= hw.rgfRoots[0]) {
            hw.nHits = 0;
        } else if (hw.nHits == 1 || t < hw.rgfRoots[1]) {
            hw.nHits = 2;
            hw.rgfRoots[1] = t;
            hw.rgvNormals[1] = cvNormal;
        }
    } else if (C > 0.0)
        hw.nHits = 0;
    return hw.nHits;
}

Hit*
Quadric::Intersect(Ray &ray) 
{
    double A, B, C;
    HitWork hw;
    Point3 p;
    int i;

    hw.nRayEnteringSolid = ray.nRayEnteringSolid;
    A = ray.vDir.x*ray.vDir.x +  ray.vDir.y*ray.vDir.y + ray.vDir.z*ray.vDir.z*P;
    B = 2.0*(ray.vDir.x*ray.pOrg.x +
             ray.vDir.y*ray.pOrg.y +
             ray.vDir.z*ray.pOrg.z*P) + ray.vDir.z*Q;
    C = ray.pOrg.x*ray.pOrg.x + ray.pOrg.y*ray.pOrg.y +
        ray.pOrg.z*ray.pOrg.z*P + ray.pOrg.z*Q + R;
    hw.nRayInsideAtInfinity = A < 0.0;
    hw.nHits = SolveQuadratic(hw.rgfRoots, A, B, C);
    for (i = 0; i < hw.nHits; i++) {
        p = ray(hw.rgfRoots[i]);
        hw.rgvNormals[i].x = p.x;
        hw.rgvNormals[i].y = p.y;
        hw.rgvNormals[i].z = p.z*P + 0.5*Q;
    }
    if (nClip) {
        ClipPlane(ray, CoVector3(0.0, 0.0,  1.0), hw);
        ClipPlane(ray, CoVector3(0.0, 0.0, -1.0), hw);
    }
    return hw.HitList();
}
//
//  Quartic solids
//
//
//  Solve cubic equation a*x**3 + b*x**2 + c*x + d == 0
//
//  Cardano's solution to the cubic is about 10 times faster than using isolation
//  and refinement (Fourier-Budan or Sturm sequences).  But great care must be
//  exercised to avoid unnecessary cancellation errors.
//
double
CubeRoot(double f)
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
CubeRoot2(double x)             // This one is much faster and much more accurate
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

int
SolveCubic(double rgfRoots[], double a, double b, double c, double d)
{
    double A, B, C, S, S3, T, U, V;
    double yInflect2, yPrimeInflect3, fDiscriminant;
    double fTheta, fCos, fSin;

    if (a == 0.0)
        return SolveQuadratic(rgfRoots, b, c, d);
    A = b/a;
    B = c/a;
    C = d/a;
    //
    //  Don't need to divide by 27 and 9 here, if we are careful, and thus
    //  more precision is preserved.
    //
    yInflect2 = 0.5*(A*(2.0*A*A - 9.0*B) + 27.0*C);
    yPrimeInflect3 = 3.0*B - A*A;
    fDiscriminant = yPrimeInflect3*yPrimeInflect3*yPrimeInflect3 + yInflect2*yInflect2;
    if (fDiscriminant > 0.0) {
        //
        //  One real root
        //
        U = CubeRoot2(fabs(yInflect2) + sqrt(fDiscriminant));
        if (yInflect2 > 0.0)
            U = -U;
        V = -yPrimeInflect3/U;
        rgfRoots[0] = (U + V - A)*1.0/3.0;
        return 1;
   } else {
        //
        //  Three real roots.  Handle the double-root case here too.
        //
        S = sqrt(-yPrimeInflect3);
        S3 = S*S*S;
        if (yInflect2 >= S3) {          // -2sqrt(Q), double sqrt(Q)
            if (yInflect2 > S3)
                goto OneRoot;           // no double root
            fCos = 1.0;
            fSin = 0.0;
        } else if (yInflect2 <= -S3) {  // double -sqrt(Q), 2sqrt(Q)
            if (yInflect2 < -S3)
                goto OneRoot;
            fCos = 0.5;
            fSin = 1.5;
        } else {
            T = yInflect2/S3;
            fTheta = acos(T)/3.0;
            fCos = cos(fTheta);
            fSin = D_SQRT3*sin(fTheta);
        }
        rgfRoots[0] =  (-2.0 * S * fCos  - A)*1.0/3.0;
        rgfRoots[1] = (S * (fCos - fSin) - A)*1.0/3.0;
        rgfRoots[2] = (S * (fCos + fSin) - A)*1.0/3.0;
        return 3;
    }
OneRoot:
    T = yInflect2/S3;
    U = CubeRoot2(T + sqrt(T*T - 1.0));
    rgfRoots[0] = (-S*(U + 1.0/U) - A)*1.0/3.0;
    return 1;
}
//
//  Solve quartic equation a*x**4 + b*x**3 + c*x**2 + d*x + e
//
int
SolveQuartic(double rgfRoots[], double a, double b, double c, double d, double e)
{
    double A, B, C, D, P, Q, R, y, EE, FF, EF, E, F, G, g, H, h;
    double rgfCubicRoots[3], x;
    int i, j, nRoots;

    if (a == 0.0)
        return SolveCubic(rgfRoots, b, c, d, e);
    A = b/a;    // ML_Normalize
    B = c/a;
    C = d/a;
    D = e/a;
    //
    //  Ferrari's resolvant cubic equation.
    //
    P = B; 
    Q = A*C - 4.0*D;
    R = D*(A*A - 4.0*B) + C*C;
    (void) SolveCubic(rgfCubicRoots, 1.0, P, Q, R);
    y = rgfCubicRoots[0];
    EE = 0.25*A*A - B - y;
    FF = 0.25*y*y - D;
    //EF = -(0.25*A*y + 0.5*C);
    EF = -0.5*A*y - C;
    if (EE < 0.0)
        E = 0.0;
    else
        E = sqrt(EE);
    if (FF < 0.0)
        F = 0.0;
    else
        F = sqrt(FF);
    if (EF < 0.0) {
            F = -F;
    }
    //
    //  The quadratic factors x**2 + G*x + H and x**2 + g*x + h
    //
    G = 0.5*A + E;
    g = 0.5*A - E;
    if (y*F >= 0.0) {   // This cancellation avoidance does improve accuracy
        h = -0.5*y - F;
        if (h)
            H = D/h;
        else
            H = -0.5*y + F;
    } else {
        H = -0.5*y + F;
        h = D/H;
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

Hit *
Torus::Intersect(Ray &ray) 
{
	double f1, f2, zOrg, zDir, zOrgDir, fR2;
    double a, b, c, d, e, x, y, z, fAxis;
    int i;
    Vector3 vDir;
    HitWork hw;

    hw.nRayEnteringSolid = ray.nRayEnteringSolid;
    vDir = Normalize(ray.vDir);
	zOrg = ray.pOrg.z * ray.pOrg.z;
	zDir = vDir.z * vDir.z;
	zOrgDir  = ray.pOrg.z * vDir.z;
    fR2 = fMinorRadius*fMinorRadius;
	f1 = ray.pOrg.x * ray.pOrg.x  +  ray.pOrg.y * ray.pOrg.y  +  zOrg  -  fR2 - 1.0;
	f2 = ray.pOrg.x * vDir.x  +  ray.pOrg.y * vDir.y  +  zOrgDir;
	a =   1.0;
	b =   4.0 * f2;
	c =   2.0 * (f1  +  2.0 * f2 * f2  +  2.0 * zDir);
	d =   4.0 * (f2*f1  +  2.0*zOrgDir);
	e =   f1 * f1  +  4.0 * (zOrg - fR2);
    hw.nRayInsideAtInfinity = a < 0.0;
    hw.nHits = SolveQuartic(hw.rgfRoots, a, b, c, d, e);
    for (i = 0; i < hw.nHits; i++) {
        x = ray.pOrg.x + vDir.x*hw.rgfRoots[i];
        y = ray.pOrg.y + vDir.y*hw.rgfRoots[i];
        z = ray.pOrg.z + vDir.z*hw.rgfRoots[i];
        fAxis = sqrt(x*x + y*y);
        hw.rgvNormals[i] = CoVector3(x - x/fAxis, y - y/fAxis, z); // from central ring to hit

    }
    return hw.HitList();
}
//
//  Polyhedra defined efficiently by the ClipPlane function
//
#define T   1.0/D_SQRT3
#define P1  D_PHI/D_SQRT3
#define P2  (D_PHI - 1.0)/D_SQRT3
#define I1  0.52573111211913
#define I2  0.85065080835204

static CoVector3 rgvTetrahedron[4] = {
    CoVector3( T,  T,  T),
    CoVector3( T, -T, -T),
    CoVector3(-T,  T, -T),
    CoVector3(-T, -T,  T)
};

static CoVector3 rgvCube[8] = {
    CoVector3( T,  T,  T),
    CoVector3( T,  T, -T),
    CoVector3( T, -T,  T),
    CoVector3( T, -T, -T),
    CoVector3(-T,  T,  T),
    CoVector3(-T,  T, -T),
    CoVector3(-T, -T,  T),
    CoVector3(-T, -T, -T)
};

static CoVector3 rgvOctahedron[6] = {
    CoVector3( 1.0,  0.0,  0.0),
    CoVector3(-1.0,  0.0,  0.0),
    CoVector3( 0.0,  1.0,  0.0),
    CoVector3( 0.0, -1.0,  0.0),
    CoVector3( 0.0,  0.0,  1.0),
    CoVector3( 0.0,  0.0, -1.0)
};

static CoVector3 rgvDodecahedron[20] = {
    CoVector3( T,  T,  T),
    CoVector3( T,  T, -T),
    CoVector3( T, -T,  T),
    CoVector3( T, -T, -T),
    CoVector3(-T,  T,  T),
    CoVector3(-T,  T, -T),
    CoVector3(-T, -T,  T),
    CoVector3(-T, -T, -T),

    CoVector3(0.0,  P2,  P1),
    CoVector3(0.0,  P2, -P1),
    CoVector3(0.0, -P2,  P1),
    CoVector3(0.0, -P2, -P1),

    CoVector3( P1, 0.0,  P2),
    CoVector3( P1, 0.0, -P2),
    CoVector3(-P1, 0.0,  P2),
    CoVector3(-P1, 0.0, -P2),

    CoVector3( P2,  P1, 0.0),
    CoVector3( P2, -P1, 0.0),
    CoVector3(-P2,  P1, 0.0),
    CoVector3(-P2, -P1, 0.0)
};

static CoVector3 rgvIcosahedron[12] = {
    CoVector3( I1, 0.0,  I2),
    CoVector3( I1, 0.0, -I2),
    CoVector3(-I1, 0.0,  I2),
    CoVector3(-I1, 0.0, -I2),

    CoVector3(0.0,  I2,  I1),
    CoVector3(0.0,  I2, -I1),
    CoVector3(0.0, -I2,  I1),
    CoVector3(0.0, -I2, -I1),

    CoVector3( I2,  I1, 0.0),
    CoVector3(-I2,  I1, 0.0),
    CoVector3( I2, -I1, 0.0),
    CoVector3(-I2, -I1, 0.0)
};

static Hit *
IntersectPolygon(Ray &ray, const CoVector3 rgv[], int nSides)
{
    int i, nHits;
    HitWork hw;

    hw.nRayEnteringSolid = ray.nRayEnteringSolid;
    if (nSides)
        nHits = ClipPlane(ray, rgv[0], hw, 0);
    for (i = 1; nHits && i < nSides; i++)
        nHits = ClipPlane(ray, rgv[i], hw);
    return hw.HitList();
}

Hit *
Tetrahedron::Intersect(Ray &ray) 
{
    return IntersectPolygon(ray, rgvTetrahedron, 4);
}

Hit *
Cube::Intersect(Ray &ray) 
{
    return IntersectPolygon(ray, rgvOctahedron, 6);
}

Hit *
Octahedron::Intersect(Ray &ray) 
{
    return IntersectPolygon(ray, rgvCube, 8);
}

Hit *
Dodecahedron::Intersect(Ray &ray) 
{
    return IntersectPolygon(ray, rgvIcosahedron, 12);
}

Hit *
Icosahedron::Intersect(Ray &ray) 
{
    return IntersectPolygon(ray, rgvDodecahedron, 20);
}

//
//  CSG methods: affine transformations and lattice set operations
//
Hit *
AffineSolid::Intersect(Ray &ray)
{
    Ray rayNew;
    Hit *phHit, *ph;

    rayNew.vDir = mInverse * ray.vDir;
    rayNew.pOrg = Point3(0.0, 0.0, 0.0) + mInverse * ((ray.pOrg - Point3(0.0, 0.0, 0.0)) - vTranslate);
    rayNew.nShadowRay = ray.nShadowRay;
    rayNew.nRayEnteringSolid = ray.nRayEnteringSolid;
    phHit = psChild->Intersect(rayNew);
    for (ph = phHit; ph; ph = ph->phNext) {
        ph->cvNormal = ph->cvNormal * mInverse;
    }
    return phHit;
}

static inline double
Max(double f1, double f2)
{
    return (f1 > f2) ? f1 : f2;
}

static inline double
Min(double f1, double f2)
{
    return (f1 < f2) ? f1 : f2;
}

static Hit *
Join(Hit *phLeft, Hit *phRight)
{
    Hit phHead, *phTail;
    double fLeft, fRight, fOld, fNew;

    //PrintHitList("Left ", phLeft);
    //PrintHitList("Right", phRight);
    if (phLeft == 0)
        return phRight;
    if (phRight == 0)
        return phLeft;
    fLeft = fRight = fOld = 0.0;
    phHead.phNext = 0;
    phTail = &phHead;
    while (phLeft && phRight) {
        PrintHitList("Left:  ", phLeft);
        PrintHitList("Right: ", phRight);
        PrintHitList("Head:  ", &phHead);
        if (phLeft->t < phRight->t) {
            fLeft = phLeft->fDensity;
            fNew = Max(fLeft, fRight);
            if (fNew != fOld) {
                phTail->phNext = phLeft;
                phTail = phTail->phNext;
                if (fLeft < fRight) {     // fix 2021 Jan 22
                    phLeft->fDensity = phRight->fDensity;      
                    phLeft->pm = phRight->pm;                   // The material property, but not surface
                }
                phLeft = phLeft->phNext;
                fOld = fNew;
            } else {
                phLeft = DeleteFirstHit(phLeft);
            }
        } else {
            fRight = phRight->fDensity;
            fNew = Max(fLeft, fRight);
            if (fNew != fOld) {
                phTail->phNext = phRight;
                phTail = phTail->phNext;
                if (fRight < fLeft) {
                    phRight->fDensity = phLeft->fDensity;
                    phRight->pm = phLeft->pm;
                }
                phRight = phRight->phNext;
                fOld = fNew;
            } else {
                phRight = DeleteFirstHit(phRight);
            }
        }
    }
    if (phLeft == 0) {
        phLeft = phRight;
        fRight = fLeft;
    }
    while (phLeft) {
        PrintHitList("Left1: ", phLeft);
        PrintHitList("Head:  ", &phHead);
        fLeft = phLeft->fDensity;
        fNew = Max(fLeft, fRight);
        if (fNew != fOld) {
            phTail->phNext = phLeft;
            phTail = phTail->phNext;
            phLeft = phLeft->phNext;
            fOld = fNew;
        } else {
            phLeft = DeleteFirstHit(phLeft);
        }
    }
    phTail->phNext = 0;
    PrintHitList("HeadN: ", &phHead);
    return phHead.phNext;
}

Material *
Max(Material *pmL, Material *pmR)
{
    if (pmL->fDensity > pmR->fDensity)
        return pmL;
    else
        return pmR;
}

Material *
Min(Material *pmL, Material *pmR)
{
    if (pmL->fDensity < pmR->fDensity)
        return pmL;
    else
        return pmR;
}

Material *
Diff(Material *pmL, Material *pmR)
{
    if (pmR->fDensity <= 0.0)
        return pmL;
    else
        return &s_mVacuum;
}

static Hit *
Merge(Hit *phLeft, Hit *phRight, Material *(*op)(Material *, Material *))
{
    Hit *ph, **pph;
    Material *pmLeft, *pmRight, *pmPrev;

    ph = 0;
    pph = &ph;
    pmLeft = pmRight = pmPrev = &s_mVacuum;
    while (phLeft || phRight) {
        if (phRight == 0 || (phLeft && phLeft->t < phRight->t)) {
            pmLeft = phLeft->pm;
            phLeft->pm  = op(pmLeft, pmRight);
            if (phLeft->pm->fDensity == pmPrev->fDensity) {
                phLeft = DeleteFirstHit(phLeft);
            } else {
                pmPrev = phLeft->pm;
                *pph = phLeft;
                pph = &phLeft->phNext;
                phLeft = phLeft->phNext;
            }
        } else {
            pmRight = phRight->pm;
            phRight->pm = op(pmLeft, pmRight);
            if (phRight->pm->fDensity == pmPrev->fDensity) {
                phRight = DeleteFirstHit(phRight);
            } else {
                pmPrev = phRight->pm;
                *pph = phRight;
                pph = &phRight->phNext;
                phRight = phRight->phNext;
            }
        }
    }
    // if (*pph != 0) printf("last pph %f\n", (*pph)->t);  // debug
    *pph = 0;
    return ph;
}

Hit *
Union::Intersect(Ray &ray) 
{
    return Merge(psLeft->Intersect(ray), psRight->Intersect(ray), Max);
}

static Hit *
Meet(Hit *phLeft, Hit *phRight)
{
    Hit phHead, *phTail;
    double fLeft, fRight, fOld, fNew;

    if (phLeft == 0)
        return 0;
    if (phRight == 0)
        return 0;
    fLeft = fRight = fOld = 0.0;
    phHead.phNext = 0;
    phTail = &phHead;
    while (phLeft && phRight) {
        if (phLeft->t < phRight->t) {
            fLeft = phLeft->fDensity;
            fNew = Min(fLeft, fRight);
            if (fNew != fOld) {
                phTail->phNext = phLeft;
                phTail = phTail->phNext;
                phLeft = phLeft->phNext;
                fOld = fNew;
            } else {
                phLeft = DeleteFirstHit(phLeft);
            }
        } else {
            fRight = phRight->fDensity;
            fNew = Min(fLeft, fRight);
            if (fNew != fOld) {
                phTail->phNext = phRight;
                phTail = phTail->phNext;
                phRight = phRight->phNext;
                fOld = fNew;
            } else {
                phRight = DeleteFirstHit(phRight);
            }
        }
    }
    if (phLeft == 0) {
        phLeft = phRight;
        fRight = fLeft;
    }
    while (phLeft) {
        fLeft = phLeft->fDensity;
        fNew = Min(fLeft, fRight);
        if (fNew != fOld) {
            phTail->phNext = phLeft;
            phTail = phTail->phNext;
            phLeft = phLeft->phNext;
            fOld = fNew;
        } else {
            phLeft = DeleteFirstHit(phLeft);
        }
    }
    phTail->phNext = 0;
    return phHead.phNext;
}

Hit *
Intersection::Intersect(Ray &ray) 
{
    return Merge(psLeft->Intersect(ray), psRight->Intersect(ray), Min);
}

static Hit *
MeetComplement(Hit *phLeft, Hit *phRight)
{
    Hit phHead, *phTail;
    double fLeft, fRight, fOld, fNew;

    if (phLeft == 0)
        return 0;
    if (phRight == 0)
        return phLeft;
    fLeft = fOld = 0.0;
    fRight = 1.0;
    phHead.phNext = 0;
    phTail = &phHead;
    while (phLeft && phRight) {
        if (phLeft->t <= phRight->t) {
            fLeft = phLeft->fDensity;
            fNew = Min(fLeft, fRight);
            if (fNew != fOld) {
                phTail->phNext = phLeft;
                phTail = phTail->phNext;
                phLeft = phLeft->phNext;
                fOld = fNew;
            } else {
                phLeft = DeleteFirstHit(phLeft);
            }
        } else {
            fRight = 1.0f - phRight->fDensity;
            fNew = Min(fLeft, fRight);
            if (fNew != fOld) {
                phRight->fDensity = fRight;                 //REVIEW: should be Left material
                phRight->cvNormal = -phRight->cvNormal;
                phTail->phNext = phRight;
                phTail = phTail->phNext;
                phRight = phRight->phNext;
                fOld = fNew;
            } else {
                phRight = DeleteFirstHit(phRight);
            }
        }
    }
    while (phRight) {
        fRight = 1.0f - phRight->fDensity;
        fNew = Min(fLeft, fRight);
        if (fNew != fOld) {
            phTail->phNext = phRight;
            phTail = phTail->phNext;
            phRight = phRight->phNext;
            fOld = fNew;
        } else {
            phRight = DeleteFirstHit(phRight);
        }
    }
    while (phLeft) {
        fLeft = phLeft->fDensity;
        fNew = Min(fLeft, fRight);
        if (fNew != fOld) {
            phTail->phNext = phLeft;
            phTail = phTail->phNext;
            phLeft = phLeft->phNext;
            fOld = fNew;
        } else {
            phLeft = DeleteFirstHit(phLeft);
        }
    }
    phTail->phNext = 0;
    return phHead.phNext;
}

Hit *
Difference::Intersect(Ray &ray) 
{
    return Merge(psLeft->Intersect(ray), psRight->Intersect(ray), Diff);
}
//
// Material and surface properties.  Much better to have these be nodes in the model than properties of solids
//
Hit *
Material::Intersect(Ray &ray)
{
    Hit *phHit, *ph;

    phHit = ps->Intersect(ray);
    for (ph = phHit; ph; ph = ph->phNext) {
        if (ph->fDensity) {
            ph->fDensity = fDensity;
            ph->pm = this;
        }
    }
    return phHit;
}

Hit *
Surface::Intersect(Ray &ray)
{
    Hit *ph;

    ph = ps->Intersect(ray);
    return ph;
}