#include "stdafx.h"
#include "RayTracer2020.h"
#pragma intrinsic(sqrt, sin, cos, tan, pow, exp)

Parent Parent::rgp[PNEWSIZE];
Parent *Parent::pNew = Parent::rgp;

int g_nTextureListCount = 0;
int g_nTextureLists = 0;
TextureList *TextureList::ptFree = 0;

inline static TextureList *
NewTextureList(TextureMap pText, TextureList *ptN)
{
    TextureList *pt;

    g_nTextureListCount++;
    if (pt = TextureList::ptFree) {
        TextureList::ptFree = pt->ptNext;
        pt->ptNext = ptN;
    } else {
        g_nTextureLists++;
        pt = new TextureList(ptN);
    }
    pt->Texture = pText;
    return pt;
}

static void
DeleteTextureList(TextureList *pt)
{
    TextureList *ptHead;

    if (ptHead = pt) {
        while (pt->ptNext) {
            --g_nTextureListCount;
            pt = pt->ptNext;
        }
        --g_nTextureListCount;
        pt->ptNext = TextureList::ptFree;
        TextureList::ptFree = ptHead;
    }
}

inline static TextureList*
Append(TextureList *ptOld, TextureList *ptNew)
{
    TextureList *pt;

    ptNew->ptNext = 0;
    if (ptOld == 0)
        return ptNew;
    for (pt = ptOld; pt->ptNext; pt = pt->ptNext)
        ;
    pt->ptNext = ptNew;
    return ptOld;
}

void
PrintTextureList(char *sz, TextureList *pt)
{
    printf("%s {", sz);
    for (; pt; pt = pt->ptNext)
        printf("  Texture %d pa %d\n", pt->Texture, pt->pa);
    printf("}\n");
}

static HitWork s_hw;
int g_nHitCount = 0;
int g_nNewHits = 0;
Hit *Hit::phFree = 0;
//
//  new/delete  - 15.8498 sec
//  Hit::phFree - 13.8529 sec
//
inline static Hit *
NewHit(Hit *phN)
{
    Hit *ph;

    g_nHitCount++;
    if (ph = Hit::phFree) {
        Hit::phFree = ph->phNext;
        ph->phNext = phN;
    } else {
        g_nNewHits++;
        ph = new Hit(phN);
    }
    ph->ps = 0;
    ph->pm = 0;
    ph->pt = 0;
    ph->paAffineTransforms = 0;
    ph->nPerturbed = 0;
    return ph;
}

void
DeleteHitList(Hit *ph)
{
    Hit *phHead;

    if (phHead = ph) {
        while (ph->phNext) {
            --g_nHitCount;
            if (ph->pt)
                DeleteTextureList(ph->pt);
            ph = ph->phNext;
        }
        if (ph->pt)
            DeleteTextureList(ph->pt);
        --g_nHitCount;
        ph->phNext = Hit::phFree;
        Hit::phFree = phHead;
    }
}

inline Hit *
DeleteFirstHit(Hit *ph)
{
    Hit *phN = 0;

    if (ph) {
        if (ph->pt)
            DeleteTextureList(ph->pt);
        phN = ph->phNext;
        ph->phNext = Hit::phFree;
        Hit::phFree = ph;
        --g_nHitCount;
    }
    return phN;
}

void
TestHitLists()
{
    Hit *ph = 0;
    int i;

    printf("Test Hit List - %d\n", g_nHitCount);
    for (i = 0; i < 100; i++)
        ph = NewHit(ph);
    printf("  Hit List - %d\n", g_nHitCount);
    ph = DeleteFirstHit(ph);
    ph = DeleteFirstHit(ph);
    ph = DeleteFirstHit(ph);
    printf("  Hit List - %d\n", g_nHitCount);
    DeleteHitList(ph);
    printf("  Hit List - %d\n", g_nHitCount);
}

int
LengthHitList(Hit *ph)
{
    if (ph)
        return 1 + LengthHitList(ph->phNext);
    else
        return 0;
}

void
PrintHitList(char *sz, Hit *ph)
{
    printf("%s { ", sz);
    for (; ph; ph = ph->phNext)
        printf("t=%0.1f d=%0.1f, ", ph->t, ph->pm->fDensity);
    printf("}\n");
}

static Surface  s_sMatte(SURF_MATTE, 0);
static AffineSolid  *s_paLastHitAffine = 0;

Hit*
HitWork::HitList()
{
    int i, iFirst, nOriginInside;
    Hit *ph;

    VERBOSE printf("Hitlist (%d)\n", nHits);
    //
    //  Ray is semi-infinite, so trim off negative roots and determine if ray
    //  origin is inside or outside the solid.
    //
VERBOSE for (i = 0; i < nHits; i++) printf("HitList[%d] = %g\n", i, rgfRoots[i]);
    for (iFirst = 0; iFirst < nHits; iFirst++)
        if (rgfRoots[iFirst] >= 0.0)
            break;
    nOriginInside = (iFirst & 1) ^ nRayInsideAtInfinity;       
    if (nOriginInside && (nHits == 0 || rgfRoots[iFirst] > fEPSILON)) {                                            // intersection starting inside ( 0, t1, t2, ... )
        --iFirst;
        rgfRoots[iFirst] = 0.0;
        rgvNormals[iFirst] = CoVector3(0.0, 0.0, 0.0);
    }
VERBOSE for (i = 0; i < nHits; i++) printf("HitList2[%d, %d, %d] = %g\n", iFirst, nOriginInside, i, rgfRoots[i]);
    //
    //  Regularization: (0, eps, ...) -> (...) or (eps, ...) -> (0.0, ...)
    //  if (!nRayEnteringSolid && nHits >= 2 && rgfRoots[iFirst] == 0.0 && rgfRoots[iFirst + 1] < fEPSILON) {
    //
    if ((nHits-iFirst) >= 2 && rgfRoots[iFirst] == 0.0 && rgfRoots[iFirst + 1] < fEPSILON) {
        iFirst += 2;
    } else if ((nHits-iFirst) && rgfRoots[iFirst] < fEPSILON) {
        rgfRoots[iFirst] = 0.0;
        rgvNormals[iFirst] = CoVector3(0.0, 0.0, 0.0);
    }
VERBOSE for (i = 0; i < nHits; i++) printf("HitList3[%d, %d] = %g\n", iFirst, i, rgfRoots[i]);
    //
    //  Build the hit list
    //
    ph = 0;
    for (i = nHits - 1; i >= iFirst; --i) {
        ph = NewHit(ph);                           // This is the only place where new Hits are allocated
        ph->t = rgfRoots[i];
        ph->cvNormal = rgvNormals[i];
        ph->paAffineTransforms = s_paLastHitAffine;
        if (nRayInsideAtInfinity) {
            ph->pm = &g_mNeutronium;
        } else {
            ph->pm = &g_mVacuum;
        }
        if (ph->t > 0.0) {
            ph->ps = &s_sMatte;
        } else {
            ph->ps = 0;             // no surface attribute, must material. ray starts inside a solid
        }
        nRayInsideAtInfinity = !nRayInsideAtInfinity;
    }
    return ph;
}
//
//  Liang-Barsky clipping of parametric line against implicit planes.
//  This can only be done on convex objects.
//
static int
ClipPlane(Ray &ray, const CoVector3 &cvNormal, HitWork &hw, int nDoClip = 1, double fPN = 1.0)
{
    double B, C, t;

    if (nDoClip && hw.nHits == 0)
        return 0;
    B = cvNormal * ray.vDir;
    C = cvNormal * ray.pOrg - fPN;          // faster form of cvNOrmal * (pRayOrg - pOnPlane)
 VERBOSE printf("ClipPlane %f t + %f:\n", B, C);
 //ray.Print();
 //cvNormal.Print();
    if (B)
        t = -C/B;
    if (!nDoClip) {
        hw.nHits = 2;
        hw.rgfRoots[0] = -INFINITY;
        hw.rgfRoots[1] = INFINITY;
    hw.nHits = 1;
    hw.rgfRoots[0] = 0.0;
        hw.nRayInsideAtInfinity = 1;        // should be sufficient to initialize to (0) with InsideInfinity == 1
    }
    hw.nRayInsideAtInfinity &= (B < 0.0);
    if (B < 0.0) {                                          // ----------t=============>
        if (hw.nHits == 2 && hw.rgfRoots[1] <= t) {         // -----r0=============r1--------t=============>
            hw.nHits = 0;
        } else if (hw.rgfRoots[0] < t) {                    // --------r0======t=========
            hw.rgfRoots[0] = t;
            hw.rgvNormals[0] = cvNormal;
        }
    } else if (B > 0.0) {                                   // =============t--------
        if (t <= hw.rgfRoots[0]) {                          // =======t-------r0=============
            hw.nHits = 0;
        } else if (hw.nHits == 1 || t < hw.rgfRoots[1]) {   // ========r0====t==========r1---------
            hw.nHits = 2;
            hw.rgfRoots[1] = t;
            hw.rgvNormals[1] = cvNormal;
        }
    } else if (C > 0.0)
        hw.nHits = 0;
    return hw.nHits;
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

Hit*
Quadric::Intersect(Ray &ray) 
{
    double A, B, C;
    HitWork hw;
    Point3 p;
    int i;

    A = ray.vDir.x*ray.vDir.x +  ray.vDir.y*ray.vDir.y + ray.vDir.z*ray.vDir.z*P;
    B = 2.0*(ray.vDir.x*ray.pOrg.x +
             ray.vDir.y*ray.pOrg.y +
             ray.vDir.z*ray.pOrg.z*P) + ray.vDir.z*Q;
    C = ray.pOrg.x*ray.pOrg.x + ray.pOrg.y*ray.pOrg.y +
        ray.pOrg.z*ray.pOrg.z*P + ray.pOrg.z*Q + R;
VERBOSE printf("Quadric Intersect A,B,C = %g %g %g\n", A, B, C);
    hw.nHits = SolveQuadratic(hw.rgfRoots, A, B, C);
    for (i = 0; i < hw.nHits; i++) {
        p = ray(hw.rgfRoots[i]);
        hw.rgvNormals[i].x = p.x;
        hw.rgvNormals[i].y = p.y;
        hw.rgvNormals[i].z = p.z*P + 0.5*Q;
    }
    if (nClip == 1) {
        ClipPlane(ray, CoVector3(0.0, 0.0,  1.0), hw);
        ClipPlane(ray, CoVector3(0.0, 0.0, -1.0), hw);
    }
    if (nClip == -1) {
        A = -A;
        B = -B;
        C = -C;
    }
    hw.nRayInsideAtInfinity = (A < 0.0);        // cones, paraboloids and hyperboloids can go to infinity in some ray directions
    return hw.HitList();
}

Extent
ConeHyperCylinder::GetExtent()
{
    double r, z;

    r = Max(1.0, R);
    if (R <= 1.0)
        z = 1.0;
    else
        z = sqrt(R*R/(R*R - 1.0));
    return Extent(Point3(-r, -r, -z), Point3(r, r, +z));     // For r outside [0, 1], it expands beyong unit cube
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
#define CBRT2   1.259921049894873164767210607278228350570251
#define CBRT4   1.587401051968199474751705639272308260391493
#define CBRT2I  0.7937005259840997373758528196361541301957467
#define CBRT4I  0.6299605249474365823836053036391141752851257

double
CubeRoot(double x)
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
        U = CubeRoot(fabs(yInflect2) + sqrt(fDiscriminant));
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
    U = CubeRoot(T + sqrt(T*T - 1.0));
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

extern int SolveQuartic2(double rgfRoots[], double A, double B, double C, double D, double E);
extern int SolveQuarticHybrid(double rgfRoots[], double A, double B, double C, double D, double E);

Hit *
Torus::Intersect(Ray &ray) 
{
	double f0, f1, f2, zOrg, zDir, zOrgDir, fR2;
    double a, b, c, d, e, x, y, z, fAxis, fT, fS;
    double db, dc, dd, de;  // derivative
    int i;
    Vector3 vDir;
    Point3  pOrg;
    HitWork hw;

    if (!GetExtent().TestIntersection(ray))
        return 0;
    fT = Norm(ray.pOrg - Point3(0.0, 0.0, 0.0)) - 1.0 - fMinorRadius - 0.1;     // this helps a lot
    fS = Norm(ray.vDir);                                                        // this doesn't really help
VERBOSE printf("Torus Intersect: fT, fS %g %g\n", fT, fS);
VERBOSE ray.pOrg.Print();
VERBOSE ray.vDir.Print();
VERBOSE printf("\n");
    vDir = (ray.vDir);
    pOrg = (ray.pOrg);
    if (fT >= 0.0) {
        vDir = vDir / fS;
        pOrg = pOrg + vDir*fT;
    } else {
        fT = 0.0;
        vDir = vDir / fS;
    }
VERBOSE pOrg.Print();
VERBOSE vDir.Print();
VERBOSE printf("scaled\n");
	zOrg = pOrg.z * pOrg.z;
	zDir = vDir.z * vDir.z;
	zOrgDir  = pOrg.z * vDir.z;
    fR2 = fMinorRadius*fMinorRadius;
    f0 = Dot(vDir, vDir);
    f1 = pOrg.x * vDir.x + pOrg.y * vDir.y + pOrg.z * vDir.z;
    f2 = pOrg.x * pOrg.x  +  pOrg.y * pOrg.y + pOrg.z * pOrg.z;
    a = f0*f0;
    b = 4.0*f0*f1;
    c = 2.0*f0*(f2 - fR2 - 1.0) + 4.0*f1*f1 + 4.0*vDir.z*vDir.z;
    d = 4.0*(f2 - fR2 - 1.0)*f1 + 8.0*pOrg.z*vDir.z;
    e = (f2 - fR2 - 1)*(f2 - fR2 - 1.0) - 4.0*(fR2 - pOrg.z*pOrg.z);
VERBOSE printf("a, b, c, d, e = %g %g %g %g %g\n", a, b, c, d, e);
    hw.nRayInsideAtInfinity = a < 0.0;
    hw.nHits = SolveQuartic(hw.rgfRoots, a, b, c, d, e);
    db = 4.0*a;
    dc = 3.0*b;
    dd = 2.0*c;
    de = d;
    
    for (i = 0; i < hw.nHits; i++) {        // refine the accuracy of roots with two steps of Newton's method
        x = hw.rgfRoots[i];
        x = x - (e + x*(d + x*(c + x*(b + x*a)))) / (de + x*(dd + x*(dc + x*db))); 
        x = x - (e + x*(d + x*(c + x*(b + x*a)))) / (de + x*(dd + x*(dc + x*db))); 
        hw.rgfRoots[i] = x;
    }
    
    for (i = 0; i < hw.nHits; i++) {
        x = pOrg.x + vDir.x*hw.rgfRoots[i];
        y = pOrg.y + vDir.y*hw.rgfRoots[i];
        z = pOrg.z + vDir.z*hw.rgfRoots[i];
        fAxis = sqrt(x*x + y*y);
        hw.rgvNormals[i] = CoVector3(x - x/fAxis, y - y/fAxis, z); // from central ring to hit
    }
    
    if (fT >= 0.0) {
        for (i = 0; i < hw.nHits; i++)
            hw.rgfRoots[i] = (hw.rgfRoots[i] + fT)/fS;
    }
    hw.nRayInsideAtInfinity = 0;            // torus is compact, not a general quartic surface
    hw.fEPSILON *= 10.0;
    return hw.HitList();
}

Extent
Torus::GetExtent()
{
    return Extent(Point3(-1.0-fMinorRadius, -1.0-fMinorRadius, -fMinorRadius), Point3(1.0+fMinorRadius, 1.0+fMinorRadius, fMinorRadius));
}
//
//  Polyhedra defined efficiently by the ClipPlane function
//
#define T   1.0/D_SQRT3
#define P1  D_PHI/D_SQRT3
#define P2  (D_PHI - 1.0)/D_SQRT3
#define I1  0.52573111211913
#define I2  0.85065080835204

CoVector3 rgvSlab[2] = {
    CoVector3( 0.0,  0.0,  1.0),
    CoVector3( 0.0,  0.0, -1.0)
};

CoVector3 rgvTetrahedron[4] = {
    CoVector3( T,  T,  T),
    CoVector3( T, -T, -T),
    CoVector3(-T,  T, -T),
    CoVector3(-T, -T,  T)
};

CoVector3 rgvCube[8] = {
    CoVector3( T,  T,  T),
    CoVector3( T,  T, -T),
    CoVector3( T, -T,  T),
    CoVector3( T, -T, -T),
    CoVector3(-T,  T,  T),
    CoVector3(-T,  T, -T),
    CoVector3(-T, -T,  T),
    CoVector3(-T, -T, -T)
};

CoVector3 rgvOctahedron[6] = {
    CoVector3( 1.0,  0.0,  0.0),
    CoVector3(-1.0,  0.0,  0.0),
    CoVector3( 0.0,  1.0,  0.0),
    CoVector3( 0.0, -1.0,  0.0),
    CoVector3( 0.0,  0.0,  1.0),
    CoVector3( 0.0,  0.0, -1.0)
};

CoVector3 rgvDodecahedron[20] = {
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

CoVector3 rgvIcosahedron[12] = {        //  I2 : I1 is golden ratio
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
/*
static Hit *
IntersectPolygon(Ray &ray, const CoVector3 rgv[], int nSides)
{
    int i, nHits;
    HitWork hw;

    hw.nRayInsideAtInfinity = 0; // (nSides == 2 && ray.vDir.z == 0.0 && ray.pOrg.z > -1.0 && ray.pOrg.z < 1.0);  // Slab
    if (nSides)
        nHits = ClipPlane(ray, rgv[0], hw, 0);
    for (i = 1; nHits && i < nSides; i++)
        nHits = ClipPlane(ray, rgv[i], hw);
    return hw.HitList();
}
*/
inline Point3
Max(const Point3 &p1, const CoVector3 &p2)
{
    return Point3(Max(p1.x, p2.x), Max(p1.y, p2.y), Max(p1.z, p2.z));
}

inline Point3
Min(const Point3 &p1, const CoVector3 &p2)
{
    return Point3(Min(p1.x, p2.x), Min(p1.y, p2.y), Min(p1.z, p2.z));
}

static Extent
ExtentPolygon(const CoVector3 rgv[], int nSides)
{
    Extent e(Point3(0.0, 0.0, 0.0), Point3(0.0, 0.0, 0.0));     //REVIEW: all wrong
    int i;

    for (i = 0; i < nSides; i++) {
        e.pMin = Min(e.pMin, rgv[i]);
        e.pMax = Max(e.pMax, rgv[i]);
    }
    return e;
}

Hit*
ConvexPolyhedron::Intersect(Ray &ray)
{
    int i, nHits;
    HitWork hw;

    hw.nRayInsideAtInfinity = 0; // (nSides == 2 && ray.vDir.z == 0.0 && ray.pOrg.z > -1.0 && ray.pOrg.z < 1.0);  // Slab
    if (nSides)
        nHits = ClipPlane(ray, rgv[0], hw, 0);
    for (i = 1; nHits && i < nSides; i++)
        nHits = ClipPlane(ray, rgv[i], hw);
    return hw.HitList();
}

Extent
ConvexPolyhedron::GetExtent()
{
    Extent e(Point3(0.0, 0.0, 0.0), Point3(0.0, 0.0, 0.0));
    int i;

    return Extent(1.0);

    for (i = 0; i < nSides; i++) {
        e.pMin = Min(e.pMin, rgv[i]);       //REVIEW: wrong.  Must use verticies of the dual polyhedron.
        e.pMax = Max(e.pMax, rgv[i]);
    }
    return e;
}

static Point3 s_rgpIn[128];
static Point3 s_rgpOut[128];

Extent
ClipPlanes(const Extent &eInitial, const CoVector3 rgcv[], int nv, double fFromOrigin)
{
    Extent e;
    Point3 rgp[2], p, p1, p2;
    double f, t;
    int iv, ip, ipIn, ipOut, i, j, nIn;

    printf("ClipPlanes %d\n", nv);
    e = eInitial;
    nIn = 0;
    rgp[0] = e.pMin;
    rgp[1] = e.pMax;
    e.Print();
    for (ip = 0; ip < 8; ip++) {                        // 8 points at corners of extent
        p.x = rgp[(ip >> 0) & 1].x;
        p.y = rgp[(ip >> 1) & 1].y;
        p.z = rgp[(ip >> 2) & 1].z;
        s_rgpIn[nIn++] = p;
        p.Print(); printf(" - %d\n", nIn);
    }
    for (iv = 0; iv < nv; iv++) {                       // for each plane cv * p - f = 0
        ipIn = ipOut = 0;
        for (ip = 0; ip < nIn; ip++) {                  // separate points by inside or outside half space
            f = rgcv[iv] * s_rgpIn[ip] - fFromOrigin;
            if (f < 0.0) {
                s_rgpIn[ipIn++] = s_rgpIn[ip];
                printf("IN: ");
            } else {
                s_rgpOut[ipOut++] = s_rgpIn[ip];
                printf("OUT: ");
            }
            s_rgpIn[ip].Print(); printf("ip %d\n", ip);
        }
        nIn = ipIn;
        for (i = 0; i < ipIn; i++) {                    // process ipIn, not the new nIn points
            p1 = s_rgpIn[i];
            for (j = 0; j < ipOut; j++) {
                p2 = s_rgpOut[j];
                t = (fFromOrigin - rgcv[iv]*p1) / (rgcv[iv] * (p2 - p1));   // plane intersect line from p1 to p2
                p = p1 + (p2 - p1)*t;
                s_rgpIn[nIn++] = p;
            }
        }
break;
    }
    printf("%d total points processed\n", nIn);
    if (nIn < 2)
        return eInitial;
    e.pMin = e.pMax = s_rgpIn[0];
    for (i = 1; i < nIn; i++) {
        e.pMin = Min(e.pMin, s_rgpIn[i]);   
        e.pMax = Max(e.pMax, s_rgpIn[i]);
    }

    return e;
}

Hit *
Slab::Intersect(Ray &ray)
{
    HitWork hw;
    int nHits;

    nHits = ClipPlane(ray, rgvSlab[0], hw, 0);
    if (nHits)
        nHits = ClipPlane(ray, rgvSlab[1], hw);
    return hw.HitList();
}

Prism::Prism(int n)
{
    int i;

    nSides = n + 2;
    rgv = new CoVector3[n + 2];        // sides and end faces
    for (i = 0; i < n; i++)
        rgv[i] = CoVector3(sin(2.0*i*D_PI/n), cos(2.0*i*D_PI/n), 0.0);
    rgv[i++] = CoVector3(0.0, 0.0, 1.0);
    rgv[i++] = CoVector3(0.0, 0.0, -1.0);
}

Fibonacci::Fibonacci(int n)
{
    int i;

    nSides = n;
    rgv = new CoVector3[n];
    for (i = 0; i < n; i++) {
        FibonacciSphere(rgv[i].rgf, i, n);
    }
}
//
//  CSG methods: affine transformations
//
//  14.8045 sec - ModelM old affine (every hit transformed in intersect routine)
//  14.4459 sec - ModelM new affine (visible hit transformed in shader)
//
Hit *
AffineSolid::Intersect(Ray &ray)
{
    Ray rayNew;
    Hit *phHit, *ph;

    rayNew.vDir = mInverse * ray.vDir; 
    rayNew.pOrg = mInverse * (ray.pOrg - vTranslate);
    rayNew.tMin = ray.tMin;
    rayNew.bPrimary = ray.bPrimary;
    s_paLastHitAffine = this;
    phHit = psChild->Intersect(rayNew);
    ray.tMin = rayNew.tMin;                     // update tMin
    return phHit;

    for (ph = phHit; ph; ph = ph->phNext) {                     //REVIEW: evaluate late, only as needed in shader (done above)
        ph->cvNormal = ph->cvNormal * mInverse;
        if (ph->nPerturbed)
            ph->cvPerturbed = ph->cvPerturbed * mInverse;
    }
    return phHit;
}

Extent
AffineSolid::GetExtent()
{
    Extent e;
    Point3 p, pTran, rgp[2];
    int i;

    if (psChild == 0)
        return Extent(0.0);
    e = psChild->GetExtent();
    rgp[0] = e.pMin;
    rgp[1] = e.pMax;
    for (i = 0; i < 8; i++) {
        p.x = rgp[(i >> 0) & 1].x;
        p.y = rgp[(i >> 1) & 1].y;
        p.z = rgp[(i >> 2) & 1].z;
        pTran = mTransform * p + vTranslate;
        if (i == 0) {
            e.pMin = e.pMax = pTran;
        } else {
            e.pMin = Min(e.pMin, pTran);
            e.pMax = Max(e.pMax, pTran);
        }
    }
    return e;
}
//
//  22.40 sec   -  ModelG uncombined affine nodes
//  19.86 sec   -  ModelG affine optimized.
//  19.9994 - ModelM no optimize
//  16.8794 - ModelM old optimize
//  16.8950 - ModelM optimize while-loop
//
static AffineSolid *s_paPreviousOptimizeAffine = 0;

Solid *
AffineSolid::Optimize()
{
    Matrix3 m, mi;
    Vector3 v;
    AffineSolid *pa, *paSavePreviousAffine;

    if (psChild == 0)
        return 0;
    while (psChild->nIsAffine) {
        pa = (AffineSolid *)psChild;
        m = mTransform * pa->mTransform;                    // merge chains of affine transforms into one transform
        mi = pa->mInverse * mInverse;
        v = vTranslate + mTransform * pa->vTranslate;
        mTransform = m;
        mInverse = mi;
        vTranslate = v;
        psChild = pa->psChild;
    }
    paParent = s_paPreviousOptimizeAffine;                  // build parent lines of Affine nodes by passing points down-call
    paSavePreviousAffine = s_paPreviousOptimizeAffine;
    VERBOSE printf("Optimize Affine %d -> affine %d\n", this, paParent);
    s_paPreviousOptimizeAffine = this;
    psChild = psChild->Optimize();
    s_paPreviousOptimizeAffine = paSavePreviousAffine;
    return this;
}
//
//  CSG - set operations (density lattice algebra rather than boolean set ops)
//
Material *
Max(Material *pmL, Material *pmR)               // join operation (union)
{
    if (pmL->fDensity > pmR->fDensity)
        return pmL;
    else
        return pmR;
}

Material *
Min(Material *pmL, Material *pmR)               // meet operation (intersection)
{
    if (pmL->fDensity < pmR->fDensity)
        return pmL;
    else
        return pmR;
}

Material *
Diff(Material *pmL, Material *pmR)              // different operation (not part of lattice algebra)
{
    if (pmR->fDensity <= 0.0)
        return pmL;
    else
        return &g_mVacuum;
}

static Hit *
Merge(Hit *phLeft, Hit *phRight, Material *(*op)(Material *, Material *))
{
    Hit *ph, **pph;
    Material *pmLeft, *pmRight, *pmPrev;

    ph = 0;
    pph = &ph;
    pmLeft = pmRight = pmPrev = &g_mVacuum;
    while (phLeft || phRight) {
        if (phRight == 0 || (phLeft && phLeft->t <= phRight->t)) {
            if (phRight && phLeft->t == phRight->t) {
                pmRight = phRight->pm;
                phRight = DeleteFirstHit(phRight);
            }
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
    *pph = 0;
    if (ph && ph->phNext == 0 && ph->t == 0.0 && ph->pm->fDensity == 0.0) {
        printf("(t=0, d=0)\n");
        ph = DeleteFirstHit(ph);
    }
    return ph;
}

Hit *
Union::Intersect(Ray &ray) 
{
    Hit *phLeft = 0, *phRight = 0;

    if (psLeft)
        phLeft = psLeft->Intersect(ray);
    if (psRight)
        phRight = psRight->Intersect(ray);
    if (phLeft) {
        if (phRight)
            return Merge(phLeft, phRight, Max);
        return phLeft;
    }
    return phRight;
}

Extent
Union::GetExtent()
{
    Extent e, eLeft, eRight;

    if (psLeft)
        eLeft = psLeft->GetExtent();
    if (psRight)
        eRight = psRight->GetExtent();
    if (psLeft) {
        if (psRight) {
            e.pMin = Min(eLeft.pMin, eRight.pMin);
            e.pMax = Max(eLeft.pMax, eRight.pMax);
            return e;
        }
        return eLeft;
    }
    if (psRight)
        return eRight;
    return Extent(0.0);
}
//
//  Top-level Union is allowed to update tMin.   
//  Only unions in the top of the CSG tree with no other set ops between them and the root
//  Does not work  if the material priorities are creating a bubble (a union that really acts like an intersect or diff)
//
Hit *
TopUnion::Intersect(Ray &ray) 
{
    Hit *phLeft = 0, *phRight = 0;

    if (psLeft)
        phLeft = psLeft->Intersect(ray);
    if (phLeft && phLeft->t < ray.tMin) {
        VERBOSE printf("TopUnion L %g to %g\n", ray.tMin, phLeft->t);
        ray.tMin = phLeft->t;
    }
    if (psRight)
        phRight = psRight->Intersect(ray);
    if (phRight && phRight->t < ray.tMin) {
        VERBOSE printf("TopUnion R %g to %g\n", ray.tMin, phLeft->t);
        ray.tMin = phRight->t;
        if (ray.bPrimary) {
            Solid *ps = psLeft;     // for primary rays, try to bring left branch closer to improve tMin optimization
            psLeft = psRight;
            psRight = ps;
        }
    }
    if (phLeft) {
        if (phRight)
            return Merge(phLeft, phRight, Max);
        return phLeft;
    }
    return phRight;
}

Extent
TopUnion::GetExtent()
{
    return Union::GetExtent();
}

Hit *
Intersection::Intersect(Ray &ray) 
{
    Hit *phLeft, *phRight;

    if (psLeft)
        phLeft = psLeft->Intersect(ray);
    else
        return 0;
    if (psRight)
        phRight = psRight->Intersect(ray);
    else
        return 0;
    if (phLeft) {
        if (phRight)
            return Merge(phLeft, phRight, Min);
        return 0;
    }
    return 0;
}

Extent
Intersection::GetExtent()
{
   Extent e, eLeft, eRight;

    if (psLeft)
        eLeft = psLeft->GetExtent();
    if (psRight)
        eRight = psRight->GetExtent();
    if (psLeft && psRight) {
        e.pMin = Max(eLeft.pMin, eRight.pMin);
        e.pMax = Min(eLeft.pMax, eRight.pMax);
    } else
        e = Extent(0.0);
    return e;
}

Hit *
Difference::Intersect(Ray &ray) 
{
    Hit *phLeft, *phRight;

    if (psLeft)
        phLeft = psLeft->Intersect(ray);
    else
        return 0;
    if (psRight)
        phRight = psRight->Intersect(ray);
    else
        return phLeft;
    if (phLeft) {
        if (phRight)
            return Merge(phLeft, phRight, Diff);
        return phLeft;
    }
    return 0;
}

Extent
Difference::GetExtent()
{
    if (psLeft)
        return psLeft->GetExtent();         // Always try to use intersect instead, smaller extents
    return Extent(0.0);
}

AngularClip::AngularClip(double f1, double f2, Solid *p)
{
    ps = p;
    if (f2 - f1 <= D_PI) {
        nIntersect = 1;
        fAng1 = f1;
        fAng2 = f2;
        cv1 = CoVector3(-cos(fAng1), +sin(fAng1), 0.0);
        cv2 = CoVector3(+cos(fAng2), -sin(fAng2), 0.0);
    } else {
        nIntersect = 0;
        fAng1 = f2;
        fAng2 = f1 + 2.0*D_PI;
        cv1 = CoVector3(-cos(fAng1), +sin(fAng1), 0.0);
        cv2 = CoVector3(+cos(fAng2), -sin(fAng2), 0.0);
    }
}

Hit *
AngularClip::Intersect(Ray &ray)
{
    HitWork hw;
    Hit *ph, *ph2;

    ph2 = ps->Intersect(ray);
    if (ph2) {
        ClipPlane(ray, cv1, hw, 0, 0.0);
        ClipPlane(ray, cv2, hw, 1, 0.0);
        ph  = hw.HitList();
        if (nIntersect)
            ph2 = Merge(ph2, ph, Min);
        else
            ph2 = Merge(ph2, ph, Diff);
    }
    return ph2;
}

Extent
AngularClip::GetExtent()
{
    Extent e;

    return ps->GetExtent();     // for now
}

Solid*
MakeUnion(Solid *ps1)
{
    return ps1;
}
//
//  Convenient union and intersect balanced-tree constructors
//
Solid* 
MakeUnion(Solid *ps1, Solid *ps2)
{
    if (ps1 == 0)
        return ps2;
    return new Union(ps1, ps2);
}

Solid* 
MakeUnion(Solid *ps1, Solid *ps2, Solid *ps3)
{
    if (ps1 == 0)
        return MakeUnion(ps2, ps3);
    return new Union(ps1, MakeUnion(ps2, ps3));
}

Solid* 
MakeUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4)
{
    if (ps1 == 0)
        return MakeUnion(ps2, ps3, ps4);
    return new Union(MakeUnion(ps1, ps2), MakeUnion(ps3, ps4));
}

Solid* 
MakeUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5)
{
    if (ps1 == 0)
        return MakeUnion(ps2, ps3, ps4, ps5);
    return new Union(MakeUnion(ps1, ps2, ps3), MakeUnion(ps4, ps5));
}

Solid* 
MakeUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6)
{
    return new Union(MakeUnion(ps1, ps2, ps3), MakeUnion(ps4, ps5, ps6));
}

Solid* 
MakeUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7)
{
    return new Union(MakeUnion(ps1, ps2, ps3), MakeUnion(ps4, ps5, ps6, ps7));
}

Solid* 
MakeUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8)
{
    return new Union(MakeUnion(ps1, ps2, ps3, ps4), MakeUnion(ps5, ps6, ps7, ps8));
}

Solid* 
MakeUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9)
{
    return new Union(MakeUnion(ps1, ps2, ps3, ps4), MakeUnion(ps5, ps6, ps7, ps8, ps9));
}

Solid* 
MakeUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9,
          Solid *ps10)
{
    return new Union(MakeUnion(ps1, ps2, ps3, ps4, ps5), MakeUnion(ps6, ps7, ps8, ps9, ps10));
}

Solid* 
MakeUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9,
          Solid *ps10, Solid *ps11)
{
    return new Union(MakeUnion(ps1, ps2, ps3, ps4, ps5), MakeUnion(ps6, ps7, ps8, ps9, ps10, ps11));
}

Solid* 
MakeUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9,
          Solid *ps10, Solid *ps11, Solid *ps12)
{
    return new Union(MakeUnion(ps1, ps2, ps3, ps4, ps5, ps6), MakeUnion(ps7, ps8, ps9, ps10, ps11, ps12));
}

Solid* 
MakeTopUnion(Solid *ps1, Solid *ps2)
{
    if (ps1 == 0)
        return ps2;
    return new TopUnion(ps1, ps2);
}

Solid* 
MakeTopUnion(Solid *ps1, Solid *ps2, Solid *ps3)
{
    if (ps1 == 0)
        return MakeTopUnion(ps2, ps3);
    return new TopUnion(ps1, MakeTopUnion(ps2, ps3));
}

Solid* 
MakeTopUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4)
{
    if (ps1 == 0)
        return MakeTopUnion(ps2, ps3, ps4);
    return new TopUnion(MakeTopUnion(ps1, ps2), MakeTopUnion(ps3, ps4));
}

Solid* 
MakeTopUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5)
{
    if (ps1 == 0)
        return MakeTopUnion(ps2, ps3, ps4, ps5);
    return new TopUnion(MakeTopUnion(ps1, ps2, ps3), MakeTopUnion(ps4, ps5));
}

Solid* 
MakeTopUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6)
{
    return new TopUnion(MakeTopUnion(ps1, ps2, ps3), MakeTopUnion(ps4, ps5, ps6));
}

Solid* 
MakeTopUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7)
{
    return new TopUnion(MakeTopUnion(ps1, ps2, ps3), MakeTopUnion(ps4, ps5, ps6, ps7));
}

Solid* 
MakeTopUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8)
{
    return new TopUnion(MakeTopUnion(ps1, ps2, ps3, ps4), MakeTopUnion(ps5, ps6, ps7, ps8));
}

Solid* 
MakeTopUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9)
{
    return new TopUnion(MakeTopUnion(ps1, ps2, ps3, ps4), MakeTopUnion(ps5, ps6, ps7, ps8, ps9));
}

Solid* 
MakeTopUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9,
          Solid *ps10)
{
    return new TopUnion(MakeTopUnion(ps1, ps2, ps3, ps4, ps5), MakeTopUnion(ps6, ps7, ps8, ps9, ps10));
}

Solid* 
MakeTopUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9,
          Solid *ps10, Solid *ps11)
{
    return new TopUnion(MakeTopUnion(ps1, ps2, ps3, ps4, ps5), MakeTopUnion(ps6, ps7, ps8, ps9, ps10, ps11));
}

Solid* 
MakeTopUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9,
          Solid *ps10, Solid *ps11, Solid *ps12)
{
    return new TopUnion(MakeTopUnion(ps1, ps2, ps3, ps4, ps5, ps6), MakeTopUnion(ps7, ps8, ps9, ps10, ps11, ps12));
}

Solid*
MakeIntersection(Solid *ps1)
{
    return ps1;
}

Solid* 
MakeIntersection(Solid *ps1, Solid *ps2)
{
    return new Intersection(ps1, ps2);
}

Solid* 
MakeIntersection(Solid *ps1, Solid *ps2, Solid *ps3)
{
    return new Intersection(ps1, MakeIntersection(ps2, ps3));
}

Solid* 
MakeIntersection(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4)
{
    return new Intersection(MakeIntersection(ps1, ps2), MakeIntersection(ps3, ps4));
}

Solid* 
MakeIntersection(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5)
{
    return new Intersection(MakeIntersection(ps1, ps2, ps3), MakeIntersection(ps4, ps5));
}

Solid* 
MakeIntersection(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6)
{
    return new Intersection(MakeIntersection(ps1, ps2, ps3), MakeIntersection(ps4, ps5, ps6));
}

Solid* 
MakeIntersection(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7)
{
    return new Intersection(MakeIntersection(ps1, ps2, ps3), MakeIntersection(ps4, ps5, ps6, ps7));
}

Solid* 
MakeIntersection(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8)
{
    return new Intersection(MakeIntersection(ps1, ps2, ps3, ps4), MakeIntersection(ps5, ps6, ps7, ps8));
}
//
// Material and surface properties.  Better to be nodes in the model than properties of solids
//
Material g_mVacuum(0.0, g_sBlack, MAT_TRANSPARENT, 0);
Material g_mNeutronium(1000000.0, g_sEqualWhite, MAT_PLASTIC, 0);

Hit *
Material::Intersect(Ray &ray)
{
    Hit *phHit, *ph;
    TextureList *pt;

    phHit = psChild->Intersect(ray);
    for (ph = phHit; ph; ph = ph->phNext) {
        if (ph->pm->fDensity > 0.0) {
            ph->pm = this;
            if (Texture) {
                pt = NewTextureList(Texture, 0);
                ph->pt = Append(ph->pt, pt);
                pt->p = ray(ph->t);
                pt->Texture = Texture;
                pt->pa = paAffineTransforms;
                VERBOSE printf("  Hit Material pt %d texture %d pa %d\n", ph->pt, Texture, ph->pt->pa);
            }
        }
    }
    return phHit;
}

Solid*
Material::Optimize()
{ 
    paAffineTransforms = s_paPreviousOptimizeAffine;  
    VERBOSE printf("Optimise Material %d Texture %d -> affine %d\n", this, Texture, paAffineTransforms);
    psChild = psChild->Optimize(); 
    return this; 
}

Hit *
Surface::Intersect(Ray &ray)
{
    Hit *phHit, *ph;
    TextureList *pt;

    phHit = psChild->Intersect(ray);
    for (ph = phHit; ph; ph = ph->phNext) {
        if (ph->t == 0.0)
            continue;                               // No surface property if ray starts inside a solid
        ph->ps = this;
        if (Texture) {
            pt = NewTextureList(Texture, 0);
            ph->pt = Append(ph->pt, pt);
            pt->p = ray(ph->t);
            pt->Texture = Texture;
            pt->pa = paAffineTransforms;
            VERBOSE printf("  Hit Surface pt %d texture %d pa %d\n", ph->pt, Texture, ph->pt->pa);
        }
    }
    return phHit;
}

Solid*
Surface::Optimize()
{ 
    paAffineTransforms = s_paPreviousOptimizeAffine;
    VERBOSE printf("Optimize Surface %d Texture %d -> affine %d\n", this, Texture, paAffineTransforms);
    psChild = psChild->Optimize(); 
    return this; 
}

//
//  Bounding-box extents
//
int
Extent::TestIntersection(Ray &ray)
{
    double B, CMin, CMax, tMin, tMax, t0, t1;
    int i;

    t0 = 0.0;
    t1 = INFINITY;  // ray.tMin;                              // formerly INFIINITY
    VERBOSE printf("Extent [%g, %g]\n", t0, t1);
    for (i = 0; i < 3; i++) {
        B = ray.vDir.rgf[i];                    // B == BMax == -BMin
        CMin = +pMin.rgf[i] - ray.pOrg.rgf[i];
        CMax = -pMax.rgf[i] + ray.pOrg.rgf[i];
        if (B) {
            tMin = +CMin/B;
            tMax = -CMax/B;
        }
        if (B < 0.0) {
            if (t1 < tMax || CMin > 0.0 || t0 >= tMin)
                return 0;
            if (t0 < tMax)
                t0 = tMax;
            if (t1 > tMin)
                t1 = tMin;
        } else if (B > 0.0) {
            if (t1 < tMin || CMax > 0.0 || t0 >= tMax)
                return 0;
            if (t0 < tMin)
                t0 = tMin;
            if (t1 > tMax)
                t1 = tMax;
        } else if (CMin > 0.0 || CMax > 0.0) {
            VERBOSE printf("Extent return [%g, %g]\n", t0, t1);
            return 0;
        }
    }
   // if (ray.tMin < t0)      // ray has already seen something closer than the bounding box
   //     return 0;
    return 1;
}

int
Extent::NullExtent()
{
    return (pMax.x <= pMin.x) || (pMax.y <= pMin.y) || (pMax.z <= pMin.z);
}

Extent
Solid::GetExtent()
{
    return Extent(Point3(0.0, 0.0, 0.0), Point3(0.0, 0.0, 0.0));
}

//
//  Compound shape library
//
Solid*
Cuboid2(double r)
{
    Solid *ps, *psCubes, *psCyl, *psSph;

    psCubes =           new Scale(1.0-r, 1.0-r, 1.0, new Cube);
    psCubes = new Union(new Scale(1.0-r, 1.0, 1.0-r, new Cube), psCubes);
    psCubes = new Union(new Scale(1.0, 1.0-r, 1.0-r, new Cube), psCubes);
    psCyl =           new Translate(+(1.0-r), +(1.0-r), 0.0, new Scale(r, r, 1.0-r, new Cylinder));
    psCyl = new Union(new Translate(+(1.0-r), -(1.0-r), 0.0, new Scale(r, r, 1.0-r, new Cylinder)), psCyl);
    psCyl = new Union(new Translate(-(1.0-r), +(1.0-r), 0.0, new Scale(r, r, 1.0-r, new Cylinder)), psCyl);
    psCyl = new Union(new Translate(-(1.0-r), -(1.0-r), 0.0, new Scale(r, r, 1.0-r, new Cylinder)), psCyl);
    ps = new Union(psCyl, psCubes);
    psCyl =           new Translate(+(1.0-r), +(1.0-r), 0.0, new Scale(r, r, 1.0-r, new Cylinder));
    psCyl = new Union(new Translate(+(1.0-r), -(1.0-r), 0.0, new Scale(r, r, 1.0-r, new Cylinder)), psCyl);
    psCyl = new Union(new Translate(-(1.0-r), +(1.0-r), 0.0, new Scale(r, r, 1.0-r, new Cylinder)), psCyl);
    psCyl = new Union(new Translate(-(1.0-r), -(1.0-r), 0.0, new Scale(r, r, 1.0-r, new Cylinder)), psCyl);
    ps = new Union(new Rotate(0.5*D_PI, 0.0, 0.0, psCyl), ps);
    psCyl =           new Translate(+(1.0-r), +(1.0-r), 0.0, new Scale(r, r, 1.0-r, new Cylinder));
    psCyl = new Union(new Translate(+(1.0-r), -(1.0-r), 0.0, new Scale(r, r, 1.0-r, new Cylinder)), psCyl);
    psCyl = new Union(new Translate(-(1.0-r), +(1.0-r), 0.0, new Scale(r, r, 1.0-r, new Cylinder)), psCyl);
    psCyl = new Union(new Translate(-(1.0-r), -(1.0-r), 0.0, new Scale(r, r, 1.0-r, new Cylinder)), psCyl);
    ps = new Union(new Rotate(0.0, 0.5*D_PI, 0.0, psCyl), ps);
    psSph =           new Translate(+(1.0-r), +(1.0-r), +(1.0-r), new Scale(r, r, r, new Sphere));
    psSph = new Union(new Translate(+(1.0-r), +(1.0-r), -(1.0-r), new Scale(r, r, r, new Sphere)), psSph);
    psSph = new Union(new Translate(+(1.0-r), -(1.0-r), +(1.0-r), new Scale(r, r, r, new Sphere)), psSph);
    psSph = new Union(new Translate(+(1.0-r), -(1.0-r), -(1.0-r), new Scale(r, r, r, new Sphere)), psSph);
    psSph = new Union(new Translate(-(1.0-r), +(1.0-r), +(1.0-r), new Scale(r, r, r, new Sphere)), psSph);
    psSph = new Union(new Translate(-(1.0-r), +(1.0-r), -(1.0-r), new Scale(r, r, r, new Sphere)), psSph);
    psSph = new Union(new Translate(-(1.0-r), -(1.0-r), +(1.0-r), new Scale(r, r, r, new Sphere)), psSph);
    psSph = new Union(new Translate(-(1.0-r), -(1.0-r), -(1.0-r), new Scale(r, r, r, new Sphere)), psSph);
    ps = new BoundingBox(new Union(psSph, ps));
    return ps;
}

Solid*
Cuboid(double r, double x, double y, double z)
{
    Solid *ps, *psCubes, *psCyl, *psSph;

    psCubes =           new Scale(x-r, y-r, z, new Cube);
    psCubes = new Union(new Scale(x-r, y, z-r, new Cube), psCubes);
    psCubes = new Union(new Scale(x, y-r, z-r, new Cube), psCubes);
    psCyl =           new Translate(+(x-r), +(y-r), 0.0, new Scale(r, r, z-r, new Cylinder));
    psCyl = new Union(new Translate(+(x-r), -(y-r), 0.0, new Scale(r, r, z-r, new Cylinder)), psCyl);
    psCyl = new Union(new Translate(-(x-r), +(y-r), 0.0, new Scale(r, r, z-r, new Cylinder)), psCyl);
    psCyl = new Union(new Translate(-(x-r), -(y-r), 0.0, new Scale(r, r, z-r, new Cylinder)), psCyl);
    ps = new Union(psCyl, psCubes);
    psCyl =           new Translate(+(x-r), +(z-r), 0.0, new Scale(r, r, y-r, new Cylinder));
    psCyl = new Union(new Translate(+(x-r), -(z-r), 0.0, new Scale(r, r, y-r, new Cylinder)), psCyl);
    psCyl = new Union(new Translate(-(x-r), +(z-r), 0.0, new Scale(r, r, y-r, new Cylinder)), psCyl);
    psCyl = new Union(new Translate(-(x-r), -(z-r), 0.0, new Scale(r, r, y-r, new Cylinder)), psCyl);
    ps = new Union(new Rotate(0.5*D_PI, 0.0, 0.0, psCyl), ps);
    psCyl =           new Translate(+(z-r), +(y-r), 0.0, new Scale(r, r, x-r, new Cylinder));
    psCyl = new Union(new Translate(+(z-r), -(y-r), 0.0, new Scale(r, r, x-r, new Cylinder)), psCyl);
    psCyl = new Union(new Translate(-(z-r), +(y-r), 0.0, new Scale(r, r, x-r, new Cylinder)), psCyl);
    psCyl = new Union(new Translate(-(z-r), -(y-r), 0.0, new Scale(r, r, x-r, new Cylinder)), psCyl);
    ps = new Union(new Rotate(0.0, 0.5*D_PI, 0.0, psCyl), ps);
    psSph =           new Translate(+(x-r), +(y-r), +(z-r), new Scale(r, r, r, new Sphere));
    psSph = new Union(new Translate(+(x-r), +(y-r), -(z-r), new Scale(r, r, r, new Sphere)), psSph);
    psSph = new Union(new Translate(+(x-r), -(y-r), +(z-r), new Scale(r, r, r, new Sphere)), psSph);
    psSph = new Union(new Translate(+(x-r), -(y-r), -(z-r), new Scale(r, r, r, new Sphere)), psSph);
    psSph = new Union(new Translate(-(x-r), +(y-r), +(z-r), new Scale(r, r, r, new Sphere)), psSph);
    psSph = new Union(new Translate(-(x-r), +(y-r), -(z-r), new Scale(r, r, r, new Sphere)), psSph);
    psSph = new Union(new Translate(-(x-r), -(y-r), +(z-r), new Scale(r, r, r, new Sphere)), psSph);
    psSph = new Union(new Translate(-(x-r), -(y-r), -(z-r), new Scale(r, r, r, new Sphere)), psSph);
    ps = new BoundingBox(new Union(psSph, ps));
    return ps;
}

Solid*
CylinderFromTo(const Point3 &pFrom, const Point3 &pTo, double r)
{
    double L;
    Solid *ps;

    L = Norm(pTo - pFrom);
    if (L == 0.0) {
        printf("Zero length cylinder\n");
        0;
    }
    // printf("L = %f\n", L);
    ps = new Scale(r, r, 0.5*L, new Translate(0.0, 0.0, -1.0, new Cylinder));
    ps = new Rotate(Vector3(0.0, 0.0, -1.0), pTo - pFrom, ps);
    ps = new Translate(pFrom, ps);
    return ps;
}

static inline Point3
Icosa(int i)
{
    return Point3(rgvIcosahedron[i].x, rgvIcosahedron[i].y, rgvIcosahedron[i].z);
}

static inline Point3
Dodeca(int i)
{
    return Point3(rgvDodecahedron[i].x, rgvDodecahedron[i].y, rgvDodecahedron[i].z);
}

static inline Point3
Octa(int i)
{
    return Point3(rgvOctahedron[i].x, rgvOctahedron[i].y, rgvOctahedron[i].z);
}

static inline Point3
Hexa(int i)
{
    return Point3(rgvCube[i].x, rgvCube[i].y, rgvCube[i].z);
}

static inline Point3
Tetra(int i)
{
    return Point3(rgvTetrahedron[i].x, rgvTetrahedron[i].y, rgvTetrahedron[i].z);
}

Solid*
TetrahedronCage(double r, double rs)
{
    int i, j;
    Solid *ps, *psCyl;

    ps = 0;
    if (r) {
        for (i = 0; i < 4; i++) {
            for (j = i + 1; j < 4; j++) {
                psCyl = CylinderFromTo(Tetra(i), Tetra(j), r);
                if (ps)
                    ps = new Union(ps, psCyl);
                else
                    ps = psCyl;
            }
        }
    }
    if (rs) {
        for (i = 0; i < 4; i++) {
            if (ps)
                ps = new Union(ps, new Translate(Tetra(i), new Scale(rs, rs, rs, new Sphere)));
            else
                ps = new Translate(Tetra(i), new Scale(rs, rs, rs, new Sphere));
        }
    }
    return new BoundingBox(ps);
}

Solid*
IcosahedronCage(double r, double rs)
{
    int i, j;
    double L;
    Point3 p1, p2;
    Solid *ps, *psCyl;

    ps = 0;
    if (r) {
        for (i = 0; i < 12; i++) {
            for (j = i + 1; j < 12; j++) {
                p1 = Icosa(i);
                p2 = Icosa(j);
                L = Norm(p1 - p2);
                if (L > 1.1)
                    continue;
                psCyl = CylinderFromTo(p1, p2, r);
                if (ps)
                    ps = new Union(ps, psCyl);
                else
                    ps = psCyl;
            }
        }
    }
    if (rs) {
        for (i = 0; i < 12; i++) {
            if (ps)
                ps = new Union(ps, new Translate(Icosa(i), new Scale(rs, rs, rs, new Sphere)));
            else
                ps = new Translate(Tetra(i), new Scale(rs, rs, rs, new Sphere));
        }
    }
    return new BoundingBox(ps);
}

Solid*
DodecahedronCage(double r, double rs)
{
    int i, j;
    double L;
    Point3 p1, p2;
    Solid *ps, *psCyl;

    ps = 0;
    if (r) {
        for (i = 0; i < 20; i++) {
            for (j = i + 1; j < 20; j++) {
                p1 = Dodeca(i);
                p2 = Dodeca(j);
                L = Norm(p1 - p2);
                if (L > 0.8)
                    continue;
                psCyl = CylinderFromTo(p1, p2, r);
                if (ps)
                    ps = new Union(ps, psCyl);
                else
                    ps = psCyl;
            }
        }
    }
    if (rs) {
        for (i = 0; i < 20; i++) {
            if (ps)
                ps = new Union(ps, new Translate(Dodeca(i), new Scale(rs, rs, rs, new Sphere)));
            else
                ps = new Translate(Dodeca(i), new Scale(rs, rs, rs, new Sphere));
        }
    }
    return new BoundingBox(ps);
}

Solid*
HexahedronCage(double r, double rs)         // Cube
{
    int i, j;
    double L;
    Point3 p1, p2;
    Solid *ps, *psCyl;

    ps = 0;
    if (r) {
        for (i = 0; i < 8; i++) {
            for (j = i + 1; j < 8; j++) {
                p1 = Hexa(i);
                p2 = Hexa(j);
                L = Norm(p1 - p2);
                if (L > 1.2)
                    continue;
                psCyl = CylinderFromTo(p1, p2, r);
                if (ps)
                    ps = new Union(ps, psCyl);
                else
                    ps = psCyl;
            }
        }
    }
    if (rs) {
        for (i = 0; i < 8; i++) {
            if (ps)
                ps = new Union(ps, new Translate(Hexa(i), new Scale(rs, rs, rs, new Sphere)));
            else
                ps = new Translate(Hexa(i), new Scale(rs, rs, rs, new Sphere));
        }
    }
    return new BoundingBox(ps);
}

Solid*
OctahedronCage(double r, double rs)      
{
    int i, j;
    double L;
    Point3 p1, p2;
    Solid *ps, *psCyl;

    ps = 0;
    if (r) {
        for (i = 0; i < 6; i++) {
            for (j = i + 1; j < 6; j++) {
                p1 = Octa(i);
                p2 = Octa(j);
                L = Norm(p1 - p2);
                if (L > 1.9)
                    continue;
                psCyl = CylinderFromTo(p1, p2, r);
                if (ps)
                    ps = new Union(ps, psCyl);
                else
                    ps = psCyl;
            }
        }
    }
    if (rs) {
        for (i = 0; i < 6; i++) {
            if (r)
                ps = new Union(ps, new Translate(Octa(i), new Scale(rs, rs, rs, new Sphere)));
            else
                ps = new Translate(Octa(i), new Scale(rs, rs, rs, new Sphere));
        }
    }
    return new BoundingBox(ps);
}

Solid*
PseudoHelix(double r, double zRepeat, double fAng)
{
    Solid *ps1, *ps2, *ps;
    double fRad, fCos;
    int n;

    fRad = asin(zRepeat/4.0);
    fCos = sqrt(1.0 - zRepeat*zRepeat/16.0);
    ps = 0;
    for (n = 0; fAng >= 0.0*D_PI; fAng -= 2.0*D_PI, n++) {
        if (fAng < D_PI) {
            ps1 = new AngularClip(0.0-0.00001, fAng+0.00001, new Torus(r));
            ps2 = 0;
        } else {
            ps1 = new AngularClip(0.0-0.00001, D_PI+0.00001, new Torus(r));
            if (fAng < 2.0*D_PI)
                ps2 = new AngularClip(D_PI-0.00001, fAng+0.00001, new Torus(r));
            else
                ps2 = new AngularClip(D_PI-0.00001, 2.0*D_PI+0.00001, new Torus(r));
            ps2 = new Translate(0.0, 1.0, 0.0, ps2);
            ps2 = new Rotate(-fRad, 0.0, 0.0, ps2);
            ps2 = new Translate(0.0, -fCos, 0.0, ps2);
        }
        ps1 = new Translate(0.0, 1.0, 0.0, ps1);
        ps1 = new Rotate(+fRad, 0.0, 0.0, ps1);
        ps1 = new Translate(0.0, -fCos, 0.0, ps1);
        if (ps2)
            ps1 = new Union(ps1, ps2);
        ps1 = new BoundingBox(ps1);
        if (n == 0)
            ps = ps1;
        else
            ps = new Union(ps, new Translate(0.0, 0.0, -zRepeat*n, ps1));
    }
    ps = new Translate(0.0, 0.0, -zRepeat/2.0, new Scale(1.0, 1.0/fCos, 1.0, ps));    // fCos to center and around the helix
    return new BoundingBox(ps);
}

Solid*
Diskoid(double r)
{
    return new BoundingBox(new Union(new Torus(r), new Scale(1.0, 1.0, r, new Cylinder)));
}

Solid*
Rodoid(double r)
{
    return new BoundingBox(MakeUnion(
        new Scale(r, r, 1.0, new Cylinder),
        new Translate(0.0, 0.0, +1.0, new Scale(r, r, r, new Sphere)),
        new Translate(0.0, 0.0, -1.0, new Scale(r, r, r, new Sphere))
    ));
}

Solid*
Squareoid(double r)
{
    double f90 = 90.0 * D_PI/180.0;

    return new BoundingBox(MakeUnion(
        new Scale(1.0, 1.0, r, new Cube),
        new Translate(+1.0, +1.0, 0.0, new Scale(r, r, r, new Sphere)),
        new Translate(+1.0, -1.0, 0.0, new Scale(r, r, r, new Sphere)),
        new Translate(-1.0, +1.0, 0.0, new Scale(r, r, r, new Sphere)),
        new Translate(-1.0, -1.0, 0.0, new Scale(r, r, r, new Sphere)),
        new Translate(+1.0, 0.0, 0.0, new Rotate(f90, 0.0, 0.0, new Scale(r, r, 1.0, new Cylinder))),
        new Translate(-1.0, 0.0, 0.0, new Rotate(f90, 0.0, 0.0, new Scale(r, r, 1.0, new Cylinder))),
        new Translate(0.0, +1.0, 0.0, new Rotate(0.0, f90, 0.0, new Scale(r, r, 1.0, new Cylinder))),
        new Translate(0.0, -1.0, 0.0, new Rotate(0.0, f90, 0.0, new Scale(r, r, 1.0, new Cylinder)))
    ));
}

Solid*
Knoboid(double r, double z)                     // centered on throat of hyperboloid
{
    double fTan, fCos, fSin, rs, zDelta;

    fTan = (1.0 - SSS(r)) / z;                  // tangent at base of ConeHyperCylinder
    fCos = 1.0/sqrt(1.0 + fTan*fTan);
    fSin = fTan * fCos;
    rs = 1.0/fCos;                              // necessary radius of sphere, so cosine is 1.0
    zDelta = z + fTan;
    return new BoundingBox(new Union(new Translate(0.0, 0.0, zDelta, new Scale(rs, rs, rs, new Sphere)), 
                new Scale(1.0, 1.0, z, new Intersection( new Scale(1.0, 1.0, 1.0, new Slab), new ConeHyperCylinder(r)))));
}

Solid*
Knoboid2(double r, double z)                    // centered on sphere
{
    double fTan, fCos, fSin, rs, zDelta;

    fTan = (1.0 - SSS(r)) / z;                  // tangent at base of ConeHyperCylinder
    fCos = 1.0/sqrt(1.0 + fTan*fTan);
    fSin = fTan * fCos;
    rs = 1.0/fCos;                              // necessary radius of sphere, so cosine is 1.0
    zDelta = z + fTan;
    return new BoundingBox(new Union(new Scale(rs, rs, rs, new Sphere), 
               new Translate(0.0, 0.0, -zDelta, new Scale(1.0, 1.0, z, 
                    new Intersection( new Scale(1.0, 1.0, 1.0, new Slab), new ConeHyperCylinder(r))))));
}

Solid*
Hyperboloidoid(double r, double rt, double z)
{
    double fTan, fCos, fSin, zDelta, fS;

    fTan = (1.0 - r*r) / z;                 // tangent at base of ConeHyperCylinder
    fCos = 1.0/sqrt(1.0 + fTan*fTan);
    fSin = fTan * fCos;
    fCos *= rt;
    fSin *= rt;
    zDelta = z + fSin;
    fS = 1.0 - fCos;
    return new BoundingBox(new Union(new Translate(0.0, 0.0, zDelta, new Scale(fS, fS, fS, Diskoid(rt/fS))), 
                new Scale(1.0, 1.0, z, new Intersection( new Cube, new ConeHyperCylinder(r)))));
}

Solid*
Slaboid(double x, double y, double z, double r)
{
    return new BoundingBox(MakeUnion(
        new Scale(x-r, y, z, new Cube),
        new Scale(x, y-r, z, new Cube),
        new Translate(+x-r, +y-r, 0.0, new Scale(r, r, z, new Cylinder)),
        new Translate(+x-r, -y+r, 0.0, new Scale(r, r, z, new Cylinder)),
        new Translate(-x+r, +y-r, 0.0, new Scale(r, r, z, new Cylinder)),
        new Translate(-x+r, -y+r, 0.0, new Scale(r, r, z, new Cylinder))
    ));
}

Solid*
Slaboidoid(double x, double y, double z, double r, double r2)
{
    Solid *ps, *ps2, *ps3;
    double f90 = 90.0 * D_PI/180.0;

    ps = MakeUnion(
        new Scale(x-r, y, z-r2, new Cube),
        new Scale(x, y-r, z-r2, new Cube),
        new Scale(x-r, y-r2, z, new Cube),
        new Scale(x-r2, y-r, z, new Cube),
        new Translate(+x-r, +y-r, 0.0, new Scale(r, r, z-r2, new Cylinder)),
        new Translate(+x-r, -y+r, 0.0, new Scale(r, r, z-r2, new Cylinder)),
        new Translate(-x+r, +y-r, 0.0, new Scale(r, r, z-r2, new Cylinder)),
        new Translate(-x+r, -y+r, 0.0, new Scale(r, r, z-r2, new Cylinder))
    );
    ps2 = MakeUnion(
        new Translate(+x-r, +y-r, +z-r2, new Scale(r-r2, r-r2, r-r2, Diskoid(r2/(r-r2)))),
        new Translate(+x-r, +y-r, -z+r2, new Scale(r-r2, r-r2, r-r2, Diskoid(r2/(r-r2)))),
        new Translate(+x-r, -y+r, +z-r2, new Scale(r-r2, r-r2, r-r2, Diskoid(r2/(r-r2)))),
        new Translate(+x-r, -y+r, -z+r2, new Scale(r-r2, r-r2, r-r2, Diskoid(r2/(r-r2)))),
        new Translate(-x+r, +y-r, +z-r2, new Scale(r-r2, r-r2, r-r2, Diskoid(r2/(r-r2)))),
        new Translate(-x+r, +y-r, -z+r2, new Scale(r-r2, r-r2, r-r2, Diskoid(r2/(r-r2)))),
        new Translate(-x+r, -y+r, +z-r2, new Scale(r-r2, r-r2, r-r2, Diskoid(r2/(r-r2)))),
        new Translate(-x+r, -y+r, -z+r2, new Scale(r-r2, r-r2, r-r2, Diskoid(r2/(r-r2))))
    );
    ps3 = MakeUnion(
        new Translate(+x-r2, 0.0, +z-r2, new Scale(r2, y-r, r2, new Rotate(f90, 0.0, 0.0, new Cylinder))),
        new Translate(+x-r2, 0.0, -z+r2, new Scale(r2, y-r, r2, new Rotate(f90, 0.0, 0.0, new Cylinder))),
        new Translate(-x+r2, 0.0, +z-r2, new Scale(r2, y-r, r2, new Rotate(f90, 0.0, 0.0, new Cylinder))),
        new Translate(-x+r2, 0.0, -z+r2, new Scale(r2, y-r, r2, new Rotate(f90, 0.0, 0.0, new Cylinder))),

        new Translate(0.0, +y-r2, +z-r2, new Scale(x-r, r2, r2, new Rotate(0.0, f90, 0.0, new Cylinder))),
        new Translate(0.0, +y-r2, -z+r2, new Scale(x-r, r2, r2, new Rotate(0.0, f90, 0.0, new Cylinder))),
        new Translate(0.0, -y+r2, +z-r2, new Scale(x-r, r2, r2, new Rotate(0.0, f90, 0.0, new Cylinder))),
        new Translate(0.0, -y+r2, -z+r2, new Scale(x-r, r2, r2, new Rotate(0.0, f90, 0.0, new Cylinder)))
    );
    return new BoundingBox(MakeUnion(ps, ps2, ps3));
}
//
//  Draw text with cylinders
//
#define FONT_UP 0xFE
#define FONT_LAST 0xFF
#define P(x,y)	((((x) & 0xF) << 4) | (((y) & 0xF) << 0))

typedef struct
{
	unsigned char vertex[8]; 
} VectorFont;

static VectorFont rgFont[ML_NFONT] = {
    { FONT_LAST },                                                                                      // ' '
    { P(4,0), P(3,2), P(5,2), P(4,0), FONT_UP, P(4,4), P(4,12), FONT_LAST },                            // '!'
    { P(2,10), P(2,6), FONT_UP, P(6,10), P(6,6), FONT_LAST },                                           // '"'
    { P(0,4), P(8,4), P(6,2), P(6,10), P(8,8), P(0,8), P(2,10), P(2,2) },                               // '#'
    { P(2,2), P(6,4), P(2,8), P(6,10), FONT_UP, P(4,12), P(4,0), FONT_LAST },                           // '$'
    { P(0,1), P(8,11), FONT_UP, P(2,10), P(2,8), FONT_UP, P(6,4), P(6,2) },                             // '%'
    { P(8,0), P(3,12), P(8,9), P(0,4), P(4,0), P(8,4), FONT_LAST },                                     // '&'
    { P(2,6), P(6,10), FONT_LAST },                                                                     // '''
    { P(6,0), P(2,4), P(2,8), P(6,12), FONT_LAST },                                                     // '('
    { P(2,0), P(6,4), P(6,8), P(2,12), FONT_LAST },                                                     // ')'
    { P(1,3), P(4,11), P(7,3), P(0,8), P(8,8), P(1,3), FONT_LAST },                                     // '*'
    { P(1,6), P(7,6), FONT_UP, P(4,9), P(4,3), FONT_LAST },                                             // '+'
    { P(2,0), P(4,2), FONT_LAST },                                                                      // ','
    { P(1,6), P(7,6), FONT_LAST },                                                                      // '-'
    { P(3,0), P(4,0), FONT_LAST },                                                                      // '.'
    { P(0,0), P(8,12), FONT_LAST },                                                                     // '/'
    { P(0,0), P(8,0), P(8,12), P(0,12), P(0,0), P(8,12), FONT_LAST },                                   // '0' to '9'   
    { P(4,0), P(4,12), P(3,10), FONT_LAST },
    { P(0,12), P(8,12), P(8,7), P(0,5), P(0,0), P(8,0), FONT_LAST },
    { P(0,12), P(8,12), P(8,0), P(0,0), FONT_UP, P(0,6), P(8,6), FONT_LAST },
    { P(0,12), P(0,6), P(8,6), FONT_UP, P(8,12), P(8,0), FONT_LAST },
    { P(0,0), P(8,0), P(8,6), P(0,7), P(0,12), P(8,12), FONT_LAST },
    { P(0,12), P(0,0), P(8,0), P(8,5), P(0,7), FONT_LAST },
    { P(0,12), P(8,12), P(8,6), P(4,0), FONT_LAST },
    { P(0,0), P(8,0), P(8,12), P(0,12), P(0,0), FONT_UP, P(0,6), P(8,6), },
    { P(8,0), P(8,12), P(0,12), P(0,7), P(8,5), FONT_LAST },
    { P(4,9), P(4,7), FONT_UP, P(4,5), P(4,3), FONT_LAST },                                             // ':'
    { P(4,9), P(4,7), FONT_UP, P(4,5), P(1,2), FONT_LAST },                                             // ';'
    { P(6,0), P(2,6), P(6,12), FONT_LAST },                                                             // '<'
    { P(1,4), P(7,4), FONT_UP, P(1,8), P(7,8), FONT_LAST },                                             // '='
    { P(2,0), P(6,6), P(2,12), FONT_LAST },                                                             // '>'
    { P(0,8), P(4,12), P(8,8), P(4,4), FONT_UP, P(4,1), P(4,0), FONT_LAST },                            // '?'
    { P(8,4), P(4,0), P(0,4), P(0,8), P(4,12), P(8,8), P(4,4), P(3,6) },                                // '@'
    { P(0,0), P(0,8), P(4,12), P(8,8), P(8,0), FONT_UP, P(0,4), P(8,4) },                               // 'A' to 'Z'
    { P(0,0), P(0,12), P(4,12), P(8,10), P(4,6), P(8,2), P(4,0), P(0,0) },
    { P(8,0), P(0,0), P(0,12), P(8,12), FONT_LAST },
    { P(0,0), P(0,12), P(4,12), P(8,8), P(8,4), P(4,0), P(0,0), FONT_LAST },
    { P(8,0), P(0,0), P(0,12), P(8,12), FONT_UP, P(0,6), P(6,6), FONT_LAST },
    { P(0,0), P(0,12), P(8,12), FONT_UP, P(0,6), P(6,6), FONT_LAST },
    { P(6,6), P(8,4), P(8,0), P(0,0), P(0,12), P(8,12), FONT_LAST },
    { P(0,0), P(0,12), FONT_UP, P(0,6), P(8,6), FONT_UP, P(8,12), P(8,0) },
    { P(0,0), P(8,0), FONT_UP, P(4,0), P(4,12), FONT_UP, P(0,12), P(8,12) },
    { P(0,4), P(4,0), P(8,0), P(8,12), FONT_LAST },
    { P(0,0), P(0,12), FONT_UP, P(8,12), P(0,6), P(6,0), FONT_LAST },
    { P(8,0), P(0,0), P(0,12), FONT_LAST },
    { P(0,0), P(0,12), P(4,8), P(8,12), P(8,0), FONT_LAST },
    { P(0,0), P(0,12), P(8,0), P(8,12), FONT_LAST },
    { P(0,0), P(0,12), P(8,12), P(8,0), P(0,0), FONT_LAST },
    { P(0,0), P(0,12), P(8,12), P(8,6), P(0,5), FONT_LAST },
    { P(0,0), P(0,12), P(8,12), P(8,4), P(0,0), FONT_UP, P(4,4), P(8,0) },
    { P(0,0), P(0,12), P(8,12), P(8,6), P(0,5), FONT_UP, P(4,5), P(8,0) },
    { P(0,2), P(2,0), P(8,0), P(8,5), P(0,7), P(0,12), P(6,12), P(8,10) },
    { P(0,12), P(8,12), FONT_UP, P(4,12), P(4,0), FONT_LAST },
    { P(0,12), P(0,2), P(4,0), P(8,2), P(8,12), FONT_LAST },
    { P(0,12), P(4,0), P(8,12), FONT_LAST },
    { P(0,12), P(2,0), P(4,4), P(6,0), P(8,12), FONT_LAST },
    { P(0,0), P(8,12), FONT_UP, P(0,12), P(8,0), FONT_LAST },
    { P(0,12), P(4,6), P(8,12), FONT_UP, P(4,6), P(4,0), FONT_LAST },
    { P(0,12), P(8,12), P(0,0), P(8,0), FONT_UP, P(2,6), P(6,6), FONT_LAST },
    { P(6,0), P(2,0), P(2,12), P(6,12), FONT_LAST },                                                    // '['
    { P(0,12), P(8,0), FONT_LAST },                                                                     // '\'
    { P(2,0), P(6,0), P(6,12), P(2,12), FONT_LAST },                                                    // ']'
    { P(2,6), P(4,12), P(6,6), FONT_LAST },                                                             // '^'
    { P(0,0), P(8,0), FONT_LAST },                                                                      // '_'
    { P(2,10), P(6,6), FONT_LAST },                                                                     // '`' grave accent
    { P(0,0), P(0,8), P(4,12), P(8,8), P(8,0), FONT_UP, P(0,4), P(8,4) },                               // lower case drawn as 'A' to 'Z'
    { P(0,0), P(0,12), P(4,12), P(8,10), P(4,6), P(8,2), P(4,0), P(0,0) },
    { P(8,0), P(0,0), P(0,12), P(8,12), FONT_LAST },
    { P(0,0), P(0,12), P(4,12), P(8,8), P(8,4), P(4,0), P(0,0), FONT_LAST },
    { P(8,0), P(0,0), P(0,12), P(8,12), FONT_UP, P(0,6), P(6,6), FONT_LAST },
    { P(0,0), P(0,12), P(8,12), FONT_UP, P(0,6), P(6,6), FONT_LAST },
    { P(6,6), P(8,4), P(8,0), P(0,0), P(0,12), P(8,12), FONT_LAST },
    { P(0,0), P(0,12), FONT_UP, P(0,6), P(8,6), FONT_UP, P(8,12), P(8,0) },
    { P(0,0), P(8,0), FONT_UP, P(4,0), P(4,12), FONT_UP, P(0,12), P(8,12) },
    { P(0,4), P(4,0), P(8,0), P(8,12), FONT_LAST },
    { P(0,0), P(0,12), FONT_UP, P(8,12), P(0,6), P(6,0), FONT_LAST },
    { P(8,0), P(0,0), P(0,12), FONT_LAST },
    { P(0,0), P(0,12), P(4,8), P(8,12), P(8,0), FONT_LAST },
    { P(0,0), P(0,12), P(8,0), P(8,12), FONT_LAST },
    { P(0,0), P(0,12), P(8,12), P(8,0), P(0,0), FONT_LAST },
    { P(0,0), P(0,12), P(8,12), P(8,6), P(0,5), FONT_LAST },
    { P(0,0), P(0,12), P(8,12), P(8,4), P(0,0), FONT_UP, P(4,4), P(8,0) },
    { P(0,0), P(0,12), P(8,12), P(8,6), P(0,5), FONT_UP, P(4,5), P(8,0) },
    { P(0,2), P(2,0), P(8,0), P(8,5), P(0,7), P(0,12), P(6,12), P(8,10) },
    { P(0,12), P(8,12), FONT_UP, P(4,12), P(4,0), FONT_LAST },
    { P(0,12), P(0,2), P(4,0), P(8,2), P(8,12), FONT_LAST },
    { P(0,12), P(4,0), P(8,12), FONT_LAST },
    { P(0,12), P(2,0), P(4,4), P(6,0), P(8,12), FONT_LAST },
    { P(0,0), P(8,12), FONT_UP, P(0,12), P(8,0), FONT_LAST },
    { P(0,12), P(4,6), P(8,12), FONT_UP, P(4,6), P(4,0), FONT_LAST },
    { P(0,12), P(8,12), P(0,0), P(8,0), FONT_UP, P(2,6), P(6,6), FONT_LAST },
    { P(6,0), P(4,2), P(4,10), P(6,12), FONT_UP, P(2,6), P(4,6), FONT_LAST },                           // '{'
    { P(4,0), P(4,5), FONT_UP, P(4,6), P(4,12), FONT_LAST },                                            // '|'
    { P(4,0), P(6,2), P(6,10), P(4,12), FONT_UP, P(6,6), P(8,6), FONT_LAST },                           // '}'
    { P(0,4), P(2,8), P(6,4), P(8,8), FONT_LAST },                                                      // '~'

    //{ P(3,12), P(5,11), P(6,8), P(4,6), P(2,6), P(0,8), P(1,11), P(3,12) },                           // degree sign
    { P(2,12), P(4,11), P(4,10), P(3,8), P(1,8), P(0,10), P(0,11), P(2,12) },
    { P(1,6), P(7,6), FONT_UP, P(4,9), P(4,3), FONT_UP, P(1,2), P(7,2) },                               // plus/minus sign
    { P(1,6), P(2,8), P(6,6), P(7,8), FONT_UP, P(1,4), P(7, 4), FONT_LAST },                            // aproximate equal
    { P(1,3), P(7,9), FONT_UP, P(7,3), P(1,9), FONT_LAST },                                             // times sign
    { P(0,2), P(8,2), P(8,10), P(0,10), P(0,2), FONT_LAST },                                            // square
    { P(4,10), P(7,8), P(8,5), P(6,2), P(2,2), P(0,5), P(1,8), P(4,10) },                               // circle
    { P(0,2), P(8,2), P(4,9), P(0,2), FONT_LAST },                                                      // triangle
    { P(0,9), P(8,9), P(4,2), P(0,9), FONT_LAST },                                                      // del
    { P(4,2), P(8,6), P(4,10), P(0,6), P(4,2), FONT_LAST },                                             // diamond
    { P(7,6), P(4,1), P(2,1), P(0,2), P(0,5), P(2, 6), P(4,6), P(7,1) },                                // alpha
    { P(0,0), P(0,10), P(3,12), P(7,10), P(3,6), P(7,2), P(4,0), P(0,2) },                              // beta
};

Solid *
SolidText(double xStart, double yStart, char *sz, int nJust, double r)
{
    double dx, dy, x, y;
    int i, j, iChar, nVertex, nPenUp;
    static double xPen, yPen, xOrg, yOrg, xJust, yJust;
    Point3 pFrom, pTo;
    Solid *ps, *psLetter;

    if (nJust == JUSTIFIED_LEFT) {
        xJust = 0.0;
        yJust = 0.0;
    } else if (nJust == JUSTIFIED_RIGHT) {
        xJust = -12.0*strlen(sz);
        yJust = 0.0;
    } else {
        xJust = -12.0*strlen(sz)/2.0;
        yJust = 0.0;
    }
    ps = 0;
    for (i = 0; sz[i]; i++) {
        xOrg = xPen = xStart + 12*i + xJust;      
        yOrg = yPen = yStart + yJust;
        nPenUp = 1;
        iChar = unsigned char(sz[i]) - ' ';
        if (iChar < 0 || iChar >= ML_NFONT)
            iChar = 0;
        psLetter = 0;
        for (j = 0; j < 8; j++) {
            nVertex = rgFont[iChar].vertex[j];
            if (nVertex == FONT_LAST)
                break;
            if (nVertex == FONT_UP) {
                nPenUp = 1;
                continue;
            }
            x = ((nVertex >> 4) & 0xF);
            y = ((nVertex >> 0) & 0xF);
            if (sz[i] >= 'a' && sz[i] <= 'z')               // smallcaps
                y *= 0.75;
            dx = x;
            dy = y;
            if (!nPenUp) {
                //  DrawLineRGB(rgb, xPen, m_nHeight - yPen, xOrg + dx, m_nHeight - yOrg - dy, 1.0*fSize);
                pFrom = Point3(xPen/12.0, -yPen/12.0, 0.0);
                pTo = Point3((xOrg+dx)/12.0, -(yOrg + dy)/12.0, 0.0);
                psLetter = MakeUnion(psLetter, CylinderFromTo(pFrom, pTo, r), new Translate(pTo, new Scale(r, r, r, new Sphere)));
            }
            xPen = xOrg + dx;
            yPen = yOrg + dy;
            psLetter = MakeUnion(psLetter, new Translate(xPen/12.0, -yPen/12.0, 0.0, new Scale(r, r, r, new Sphere)));
		    nPenUp = 0;
        }
        if (psLetter)
            ps = MakeUnion(ps, new BoundingBox(psLetter));
    }
    return new BoundingBox(ps);
}