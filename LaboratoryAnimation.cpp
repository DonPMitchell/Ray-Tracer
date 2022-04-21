//
//  Shape grammer experiment, laboratory glassware animation
//  D.P. Mitchell  2021/03/30.
//
#include "stdafx.h"
#include "RayTracer2020.h"
#include "TextureMaps.h"
#include "Utilities.h"
#include "Video.h"

#define THICK 0.05
#define TUBE  0.25

static Solid*
Connect(double fAng)
{
    Solid *ps, *ps2, *psHyp, *psTor1, *psTor2;
    double f90, w, r, R, R1, R2, theta, phi, fHyp;

    w = 3.0*TUBE;
    r = 0.5*TUBE;
    R = w - 2.0*r;
    f90 = 90.0 * D_PI/180.0;
    R1 = TUBE;
    R2 = 3.0/4.0 * TUBE;
    phi = atan(3.0/4.0 / 1.0);
    theta = 0.5*D_PI - phi;
    // printf("theta = %f\n", theta * 180.0/D_PI);

    fHyp = 1.0;
    psHyp = MakeIntersection(
        new Scale(1.0, 1.0, 1.0*fHyp, new ConeHyperCylinder(r/TUBE)),
        new Scale(1.0, 1.0, 0.5*fHyp, new Translate(0.0, 0.0, 1.0, new Cube))
    );
    psHyp = new Translate(0.0, 0.0, w - 3.0*TUBE - 0.00001, new Scale(TUBE, TUBE, 2.0*TUBE, psHyp));
    // psHyp = new Intersection(psHyp, new Scale(w, w, w, new Sphere));

    psTor1 = new Scale(R1, R1, R1, new AngularClip(D_PI-0.00001, 1.5*D_PI+0.00001, new Torus(r/R1)));
    psTor1 = new Translate(w - 2.0*TUBE, 0.0, w - 3.0*TUBE, new Rotate(f90, 0.0, 0.0, psTor1));
    psTor2 = new Scale(R2, R2, R2, new AngularClip(theta, 1.0*D_PI, new Torus(r/R2)));
    psTor2 = new Translate(TUBE, R2, w - 3.0*TUBE - TUBE, psTor2);
    ps = MakeUnion(
       psHyp, psTor1, psTor2
    );

    psHyp = MakeIntersection(
        new Scale(1.0, 1.0, 1.0*fHyp, new ConeHyperCylinder((r-THICK)/(TUBE-THICK))), 
        new Scale(1.0, 1.0, 0.5*fHyp+0.0001, new Translate(0.0, 0.0, 1.0-0.00001, new Cube))
    );
    psHyp = new Translate(0.0, 0.0, w - 3.0*TUBE - 0.0000, new Scale(TUBE-THICK, TUBE-THICK, 2.0*TUBE+0.0001, psHyp));
    // psHyp = new Intersection(psHyp, new Scale(w, w, w, new Sphere));

    psTor1 = new Scale(R1, R1, R1, new AngularClip(D_PI-0.001, 1.5*D_PI+0.001, new Torus((r-THICK)/R1)));
    psTor1 = new Translate(w - 2.0*TUBE, 0.0, w - 3.0*TUBE, new Rotate(f90, 0.0, 0.0, psTor1));

    psTor2 = new Scale(R2, R2, R2, new AngularClip(theta-0.0001, 1.0*D_PI+0.0001, new Torus((r-THICK)/R2)));
    psTor2 = new Translate(TUBE, R2, w - 3.0*TUBE - TUBE, psTor2);
    ps2 = MakeUnion(
        psHyp, psTor1, psTor2
    );

    ps = new Difference(ps, ps2);
    ps = new Translate(0.0, 0.0, 2.0*TUBE - w - TUBE, ps);
    ps = new Rotate(0.0, 0.0, fAng + theta, ps);
    return new BoundingBox(ps);
}

static Solid*
TestHyperBug(double fAng)
{
    Solid *ps;
    double r, f90;

    r = 0.5*TUBE;
    f90 = 0.5*D_PI;
    ps = MakeIntersection(
        new ConeHyperCylinder(r/TUBE),
        new Scale(1.0, 1.0, 1.0, new Translate(0.0, 0.0, 0.0, new Cube))
    );
    ps = new Rotate(-f90, 0.0, 0.0, ps);
    return ps;
}

static Solid*
SimpleTube(double L)
{
    Solid *ps, *ps2;
    double f90;

    f90 = 90.0 * D_PI/180.0;
    ps  = new Scale(TUBE, TUBE, 0.5*L, new Cylinder);
    ps2 = new Scale(TUBE-THICK, TUBE-THICK, 0.5001*L, new Cylinder);
    return new BoundingBox(new Rotate(-f90, 0.0, 0.0, new Difference(ps, ps2)));
}

static Solid*
SpiralTube(double L)
{
    Solid *ps, *ps2, *psBead;
    double f90, w, R, r, r1, z, zSphereCap1, lMainTube, zSphereCap2;
    double fAng, fRad;

    f90 = 90.0 * D_PI/180.0;
    w = 3.0*TUBE;
    zSphereCap1 = -w - 1.5*TUBE;
    lMainTube = L + 2.0*zSphereCap1;
    zSphereCap2 = -L - zSphereCap1;
    r = 0.5*TUBE;
    R = w - 2.0*r;      // - 2.0*THICK;
    z = 2.0*TUBE;
    fAng = 2.0*D_PI * lMainTube/z;
    fRad = fAng + 0.5*D_PI;
    while (fRad > 2.0*D_PI)
        fRad -= 2.0*D_PI;
    ps = new Scale(R, R, R, PseudoHelix(r/R, z/R, fAng));
    r1 = r*1.05;
    psBead = new Difference(new Scale(r1, r1, r1, new Sphere), new Rotate(0.5*D_PI, 0.0, 0.0,
                            new Scale(r-THICK, r-THICK, r1+0.001, new Cylinder)));
    ps = new Union(ps, new Rotate(0.0, 0.0, 0.5*D_PI,  new Translate(R, 0.0, 0.0, psBead)));

    psBead = new Difference(new Scale(r1, r1, r1, new Sphere), new Rotate(0.5*D_PI, 0.0, 0.0,
                            new Scale(r-THICK, r-THICK, r1+0.001, new Cylinder)));
    ps = new Union(ps, new Rotate(0.0, 0.0, -fRad+1.0*D_PI, new Translate(R, 0.0, -lMainTube, psBead)));

    ps2 = new Scale(R, R, R, PseudoHelix((r-THICK)/R, z/R, fAng));
    ps = new Translate(0.0, 0.0, zSphereCap1, new Difference(ps, ps2));

   // ps = new Scale(0.1, 0.1, 0.1, new Sphere);

    ps = new Union(ps, new Rotate(0.0, 0.0, 1.0*D_PI, new Translate(0.0, 0.0, zSphereCap1 + w, new Scale(1.0, -1.0, +1.0, Connect(0.0)))));
    ps = new Union(ps, new Rotate(0.0, 0.0, -fRad+0.5*D_PI,  new Translate(0.0, 0.0, zSphereCap2 - w, new Scale(1.0, +1.0, -1.0, Connect(0.0)))));
    //ps = new Union(ps, new Translate(0.0, 0.0, -2.0*TUBE - 1.5*TUBE, new Scale(1.0, -1.0, 1.0, Connect(D_PI))));
    //ps = new Union(ps, new Translate(0.0, 0.0, -L + 2.0*TUBE + 1.5*TUBE, new Scale(1.0, 1.0, -1.0, Connect(D_PI))));
    ps = new Rotate(-f90, 0.0, 0.0, ps);
    ps = new Rotate(0.0, 0.5*D_PI, 0.0, ps);
    return new BoundingBox(ps);
}

static Solid*
CondensorEnvelope(double L)
{
    Solid *ps, *ps2;
    double f90, w, zEndTube1, zSphereCap1, zMainTube, lMainTube, zSphereCap2, zEndTube2;

    f90 = 90.0 * D_PI/180.0;
    w = 3.0*TUBE;
    zEndTube1 = 0.0;
    zSphereCap1 = -w - 1.5*TUBE;
    zMainTube = zSphereCap1;            // start main wall at equator of end cap sphere
    lMainTube = L + 2.0*zSphereCap1;
    zSphereCap2 = -L - zSphereCap1;     // zSphereCap1 - lMainTube
    zEndTube2 = -L - zEndTube1;
    ps = MakeUnion(
        new Translate(0.0, 0.0, zEndTube1, new Scale(TUBE+THICK, TUBE+THICK, 2.0*THICK, new Sphere)),
        new Translate(0.0, 0.0, zEndTube1, new Scale(TUBE, TUBE, 2.0*TUBE, new Translate(0.0, 0.0, -1.0, new Cylinder))),
        new Translate(0.0, 0.0, zSphereCap1, new Scale(w, w, w, new Sphere)),
        new Translate(w+0.5*TUBE, 0.0, zSphereCap1, new Rotate(0.0, f90, 0.0, new Scale(0.5*TUBE, 0.5*TUBE, TUBE, new Cylinder))),

        new Translate(0.0, 0.0, zMainTube, new Scale(w, w, 0.5*lMainTube, new Translate(0.0, 0.0, -1.0, new Cylinder))),

        new Translate(w+0.5*TUBE, 0.0, zSphereCap2, new Rotate(0.0, f90, 0.0, new Scale(0.5*TUBE, 0.5*TUBE, TUBE, new Cylinder))),
        new Translate(0.0, 0.0, zSphereCap2, new Scale(w, w, w, new Sphere)),
        new Translate(0.0, 0.0, zEndTube2, new Scale(TUBE, TUBE, 2.0*TUBE, new Translate(0.0, 0.0, +1.0, new Cylinder))),
        new Translate(0.0, 0.0, -L, new Scale(TUBE+THICK, TUBE+THICK, 2.0*THICK, new Sphere))
    );
    ps2 = MakeUnion(
        new Translate(0.0, 0.0, zEndTube1, new Scale(TUBE-THICK, TUBE-THICK, 2.5*TUBE, new Translate(0.0, 0.0, -0.75, new Cylinder))),
        new Translate(0.0, 0.0, zSphereCap1, new Scale(w-THICK, w-THICK, w-THICK, new Sphere)),
        new Translate(w+0.5*TUBE, 0.0, zSphereCap1, new Rotate(0.0, f90, 0.0, new Scale(0.5*TUBE-THICK, 0.5*TUBE-THICK, 1.01*TUBE, new Cylinder))),

        new Translate(0.0, 0.0, zMainTube, new Scale(w-THICK, w-THICK, 0.5*lMainTube, new Translate(0.0, 0.0, -1.0, new Cylinder))),

        new Translate(w+0.5*TUBE, 0.0, zSphereCap2, new Rotate(0.0, f90, 0.0, new Scale(0.5*TUBE-THICK, 0.5*TUBE-THICK, 1.01*TUBE, new Cylinder))),
        new Translate(0.0, 0.0, zSphereCap2, new Scale(w-THICK, w-THICK, w-THICK, new Sphere)),
        new Translate(0.0, 0.0, zEndTube2, new Scale(TUBE-THICK, TUBE-THICK, 2.5*TUBE, new Translate(0.0, 0.0, +0.75, new Cylinder)))
    );
    ps = new Difference(ps, ps2);
    return new BoundingBox(new Rotate(-f90, 0.0, 0.0, ps));
}

static Solid*
StopCock(double fAng)
{
    Solid *ps, *psShell, *psShell2, *psCyl;
    double f90;

    f90 = 90.0 * D_PI/180.0;
    fAng = fAng * D_PI/180.0;
    psShell = new Union(new Scale(TUBE, TUBE, 2.0*TUBE, new Cylinder),
                    new Union(new Translate(0.0, 0.0, +2.0*TUBE, new Scale(TUBE, TUBE, 2.0*THICK, new Sphere)),
                              new Translate(0.0, 0.0, -2.0*TUBE, new Scale(TUBE+THICK, TUBE+THICK, 2.0*THICK, new Sphere))
                        ));
    psShell2 = new Union(new Scale(TUBE, TUBE, 1.5*TUBE, new Cylinder),
                    new Union(new Translate(0.0, 0.0, +1.5*TUBE, new Scale(TUBE, TUBE, 2.0*THICK, new Sphere)),
                              new Translate(0.0, 0.0, -1.5*TUBE, new Scale(TUBE+THICK, TUBE+THICK, 2.0*THICK, new Sphere))
                        ));
    psShell = new Union(psShell, new Rotate(-f90, 0.0, 0.0, psShell2));
    psShell = new Difference(psShell, new Scale(TUBE-THICK, TUBE-THICK, 2.5*(TUBE+THICK), new Cylinder));
    psShell = new Difference(psShell, new Rotate(-f90, 0.0, 0.0, new Scale(TUBE-THICK, TUBE-THICK, 2.0*(TUBE+THICK), new Cylinder)));
    psCyl = new Union(new Union(new Scale(TUBE-THICK-0.01, TUBE-THICK-0.01, 2.0*TUBE, new Cylinder),
                        new Translate(0.0, 0.0, -2.0*TUBE, new Scale(TUBE-THICK-0.01, TUBE-THICK-0.01, TUBE-THICK, new Sphere))),
                      new Translate(0.0, 0.0, -2.5*(TUBE+THICK), new Scale(0.3*TUBE, 1.5*TUBE, 0.5*TUBE, new Sphere)));
    psCyl = new Difference(psCyl, new Rotate(-f90*0.8, 0.0, 0.0, new Scale(0.25*TUBE, 0.25*TUBE, TUBE+THICK, new Cylinder)));
    ps = new Union(psShell, new Rotate(0.0, 0.0, fAng, psCyl));
    return new BoundingBox(new Rotate(0.0, f90, 0.0, ps));
}

static Solid*
SeparatorFunnel(double r, double z, double fNeck, int nStopCock = 1)
{
    Solid *ps, *ps2;
    double f90, fTan, rs, fOff;

    f90 = 90.0 * D_PI/180.0;
    fTan = (1.0 - TUBE*TUBE) / z;    
    rs = sqrt(1.0 + fTan*fTan);
    fOff = z + fTan + rs;
    ps = new Union(new Difference(Knoboid(TUBE/r, z), new Translate(0.0, 0.0, -1.0, new Cube)), 
            new Union(new Translate(0.0, 0.0, fOff, new Scale(TUBE, TUBE, fNeck, new Cylinder)),
            new Translate(0.0, 0.0, fOff + fNeck, new Scale(TUBE+THICK, TUBE+THICK, 2.0*THICK, new Sphere))));
    ps2 = new Union(new Scale(1.0-THICK/rs, 1.0-THICK/rs, 1.0-THICK/fOff, Knoboid((TUBE-THICK)/(1.0-THICK/rs)/r, z)), 
                new Translate(0.0, 0.0, fOff, new Scale(TUBE-THICK, TUBE-THICK, fNeck+0.1, new Cylinder))
            );
    ps = new Scale(r, r, 1.0, new Difference(ps, ps2));
    ps = new Rotate(f90, 0.0, 0.0, ps);
    if (nStopCock) {
        ps = new Union(ps, new Translate(0.0, 1.5*TUBE, 0.0, StopCock(0.0)));
    }
    return new BoundingBox(ps);
}

static Solid*
ErlenmeyerFlask(double r, double fNeck)
{
    Solid *ps, *ps2, *psNeck, *psBody;
    double f20, f90, rt, z;

    f20 = 20.0 * D_PI/180.0;
    f90 = 90.0 * D_PI/180.0;
    rt = 0.15;
    z = 2.0;

    psBody = new Scale(r, r, r, 
            new Intersection(new Scale(1.01+rt, 1.01+rt, 1.01+rt, new Translate(0.0, 0.0, 1.0, new Cube)), Hyperboloidoid(TUBE/r, rt, z)));
    psNeck = new Union(new Translate(0.0, 0.0, 0.0, new Scale(TUBE, TUBE, fNeck, new Cylinder)),
                 new Translate(0.0, 0.0, 0.0 - fNeck, new Scale(TUBE+THICK, TUBE+THICK, 2.0*THICK,  new Sphere)));
    ps = new Rotate(-f90, 0.0, 0.0, new Union(psNeck, psBody));

    psBody = new Scale(r-2.0*THICK, r-2.0*THICK, r-THICK, 
            new Intersection(new Scale(1.01+rt, 1.01+rt, 1.01+rt, new Translate(0.0, 0.0, 1.0, new Cube)), 
                                        Hyperboloidoid((TUBE-THICK)/(r-2.0*THICK), rt, z)));
    psNeck = new Translate(0.0, 0.0, 0.0, new Scale(TUBE-THICK, TUBE-THICK, fNeck + 2.1*THICK, new Cylinder));
    ps2 = new Rotate(-f90, 0.0, 0.0, new Union(psNeck, psBody));

    return new BoundingBox(new Difference(ps, ps2));
}

static Solid*
FlorenceFlask(double r, double fNeck, int nNecks = 1)
{
    Solid *ps, *ps2, *psNeck;
    double f20, f90, fRad;
    int i;

    f20 = 20.0 * D_PI/180.0;
    f90 = 90.0 * D_PI/180.0;
    fRad = 0.0;
    if (nNecks > 1)
        fRad = 2.0*D_PI/double(nNecks - 1);

    ps = new Union(new Scale(r, r, r, new Sphere),
             new Union(new Rotate(f90, 0.0, 0.0, new Translate(0.0, 0.0, r, new Scale(TUBE, TUBE, fNeck, new Cylinder))),
               new Rotate(f90, 0.0, 0.0, new Translate(0.0, 0.0, r + fNeck, new Scale(TUBE+THICK, TUBE+THICK, 2.0*THICK, new Sphere)))));
    for (i = 0; i < nNecks - 1; i++) {
        psNeck = new Union(new Rotate(f90, 0.0, 0.0, new Translate(0.0, 0.0, r, new Scale(TUBE, TUBE, fNeck, new Cylinder))),
               new Rotate(f90, 0.0, 0.0, new Translate(0.0, 0.0, r + fNeck, new Scale(TUBE+THICK, TUBE+THICK, 2.0*THICK, new Sphere))));
        ps = new Union(ps, new Rotate(2.0*f20, fRad*i + 0.5*D_PI, 0.0, psNeck));
    }

    ps2 = new Union(new Scale(r-THICK, r-THICK, r-THICK, new Sphere),
             new Rotate(f90, 0.0, 0.0, new Translate(0.0, 0.0, r, new Scale(TUBE-THICK, TUBE-THICK, fNeck+2.1*THICK, new Cylinder))));
    for (i = 0; i < nNecks - 1; i++) {
        psNeck = new Rotate(f90, 0.0, 0.0, new Translate(0.0, 0.0, r, new Scale(TUBE-THICK, TUBE-THICK, fNeck+2.1*THICK, new Cylinder)));
        ps2 = new Union(ps2, new Rotate(2.0*f20, fRad*i + 0.5*D_PI, 0.0, psNeck));
    }

    return new BoundingBox(new Translate(0.0, r, 0.0, new Difference(ps, ps2)));
}

static Solid*
Retort()
{
    Solid *ps, *ps2, *psCone;
    double f20, f90;

    f20 = 20.0 * D_PI/180.0;
    f90 = 90.0 * D_PI/180.0;
    psCone = new Intersection(new ConeFX, new Translate(0.0, 0.0, 1.0, new Slab));  // only use this once
    ps = new Union(new Sphere,
           new Union(new Rotate(f20, 0.0, 0.0, new Translate(0.0, -0.6, -4.0, new Scale(0.4, 0.4, 2.0, psCone))),
             new Union(new Rotate(f90, 0.0, 0.0, new Translate(0.0, 0.0, 1.0, new Scale(0.3, 0.3, 0.4, new Cylinder))),
               new Rotate(f90, 0.0, 0.0, new Translate(0.0, 0.0, 1.4, new Scale(0.35, 0.35, 0.1, new Sphere))))));

    psCone = new Intersection(new ConeFX, new Translate(0.0, 0.0, 1.0, new Slab));  // only use this once
    ps2 = new Union(new Scale(0.9, 0.9, 0.9, new Sphere),
           new Union(new Rotate(f20, 0.0, 0.0, new Translate(0.0, -0.6, -4.0, new Scale(0.3, 0.3, 2.0, psCone))),
             new Union(new Rotate(f90, 0.0, 0.0, new Translate(0.0, 0.0, 1.0, new Scale(0.2, 0.2, 0.6, new Cylinder))),
               new Rotate(f20, 0.0, 0.0, new Translate(0.0, -0.6, -6.0, new Scale(1.0, 1.0, 3.0, new Cube))))));
    return new BoundingBox(new Rotate(0.0, -0.5*D_PI, 0.0, new Difference(ps, ps2)));
}

static Solid*
Stand()
{
    Solid *psStand;
    double f;

    f = 1.05/sqrt(2.0);
    psStand = new Torus(0.05/f);
    psStand = new Union(psStand, new Translate(1.0, 0.0, -0.5, new Scale(0.05, 0.05, 0.5, new Cylinder)));
    psStand = new Union(psStand, new Translate(-0.5, +0.866, -0.5, new Scale(0.05, 0.05, 0.5, new Cylinder)));
    psStand = new Union(psStand, new Translate(-0.5, -0.866, -0.5, new Scale(0.05, 0.05, 0.5, new Cylinder)));
    psStand = new Translate(0.0, f, 0.0, new Rotate(0.5*D_PI, 0.0, 0.0, new Scale(f, f, f, psStand)));
    psStand = new Material(0.9, g_sCopper, MAT_COPPER,
               new Surface(SURF_POLISHED, 0, psStand));
    psStand = new BoundingBox(psStand);
    return psStand;
}

static Solid*
Stand2(double r, double L)
{
    Solid *ps;
    double fArm;

    fArm = r + 2.0*TUBE;
    ps = MakeUnion(
        new Translate(fArm, 0.0, 0.0, new Scale(TUBE, TUBE, L, new Cylinder)),
        new Translate(fArm, 0.0, L, new Scale(TUBE+THICK, TUBE+THICK, 2.0*THICK, new Sphere)),
        new Translate(0.0, 0.0, 0.0, new Scale(TUBE+THICK, TUBE+THICK, 3.0*THICK, new Cylinder))
    );
    return ps;
}

static Solid*
Glass(Solid *ps)
{
    return new Material(0.33, REFRACT_PYREX, g_sWater, MAT_TRANSPARENT, ColoredGlass,
               new Surface(SURF_POLISHED, 0, ps));
}

static Solid*
Plastic(Solid *ps)
{
    Spectrum sBlue;

    sBlue = (g_sBlue + g_sEqualWhite) * 0.5;
    return new Material(0.33, REFRACT_WATER, sBlue, MAT_PLASTIC, ColoredGlass,
               new Surface(SURF_POLISHED, 0, ps));
}

static Solid*
Slice(Solid *ps)
{
    return new Difference(ps, new Scale(5.0, 5.0, 5.0, new Translate(0.0, 0.0, 1.0, new Cube)));
}

static Solid*
Environment(Solid *ps)
{
    Spectrum sBlue;
    double f20;

    sBlue = (g_sBlue + g_sEqualWhite) * 0.5;
    f20 = 20.0 * D_PI/180.0;
    if (ps)
        ps = new Union(ps, new Translate(0.0, 5.0, -10.0, new Scale(12.0, 1.0, 12.0, 
                    new Material(0.9, g_sEqualWhite, MAT_PLASTIC, new Cube))));
    else
        ps = new Translate(0.0, 5.0, -10.0, new Scale(15.0, 1.0, 15.0, 
                    new Material(0.9, g_sEqualWhite, MAT_PLASTIC, new Cube)));
    ps = new Union(ps, new Scale (100.0, 100.0, 100.0, new Material(0.9, sBlue, MAT_PLASTIC, new HollowSphere)));
    ps = new Rotate(-f20, 0.0, 0.0, ps);
    return ps;
}

void
Laboratory()
{
    RayTracer scene;
    Image im;
    Video vid;
    ML_TimingInfo ti;
    Spectrum sBlue, s6500;
    Solid *ps;
    double t, fSqrt2;

    sBlue = (g_sBlue + g_sEqualWhite) * 0.5;
    s6500 = BlackBody(6500.0);
    fSqrt2 = 1.05/sqrt(2.0);
    scene.plAmbient = new HavercosineLight(sBlue * 0.02, s6500 * 0.02);
    scene.plPointLights = new PointLight(Point3(9.0, -6.5, 12.0), s6500 * 10000.0);
//  scene.plPointLights = new PointLight(Point3(-10.0, -10.0, -10.0), s6500 * 500.0, scene.plPointLights);
    scene.psSample = new FastBlueSamplingTile;
    scene.pcCamera = new StereographicCamera(Point3(0.0, 0.0, 12.0)); 

    im.NewImage(1920, 1080, 3);
    vid.NewVideo("TestRayTracer.mpg", im);
    ML_StartTiming(ti);
   for (t = 3.0; t <= 6.0; t += 3.0/120.0) {
   // for (t = 60.0; t < 360.0; t += 5.0) {
        ps = MakeUnion(
                //new Translate(-2.0, 0.0, 0.0, new Rotate(0.0, 0.0, 0.0, Connect(30.0))),
                    //new Translate(+2.0, 0.0, 0.0, new Rotate(0.6*D_PI, 0.0, 0.0, Connect(30.0))),
                new Translate(-6.5, 2.5, -8.0, new Material(0.9, sBlue, MAT_GOLD, new Surface(SURF_POLISHED, Bumpy, new Sphere))),
                new Translate(-3.1, 2.5, -9.3, new Material(0.9, sBlue, MAT_IRON, new Surface(SURF_POLISHED, Bumpy, new Sphere))),
                new Translate(-0.9, 2.5, -10.7, new Material(0.9, sBlue, MAT_ALUMINUM, new Surface(SURF_POLISHED, Bumpy, new Sphere))),
                new Translate(+3.3, 2.5, -9.0, new Material(0.9, sBlue, MAT_OSMIUM, new Surface(SURF_POLISHED, Bumpy, new Sphere))),
                new Translate(+6.4, 2.5, -10.9, new Material(0.9, sBlue, MAT_COPPER, new Surface(SURF_POLISHED, Bumpy, new Sphere))),
                Glass(new Translate(+2.5, 3.5,     -0.9, new Rotate(0.5*D_PI, 0.05*D_PI, 0.0, StopCock(-20.0)))),
                Glass(new Translate(-4.5, 2.8-1.3, -2.0, (FlorenceFlask(1.2, 0.8, 5)))),
                      new Translate(-4.5, 2.8,     -2.0, Stand()),
                Glass(new Translate(+4.0, 0.5,     -4.5, (ErlenmeyerFlask(1.5, 0.5)))),
                Glass(new Translate(+0.0, 0.5*t,   -0.5, new Rotate(0.0, -30.0 * D_PI/180.0, 0.0,CondensorEnvelope(t)))),
                Glass(new Translate(+0.0, 0.5*t,   -0.5, SpiralTube(t)))
                    // new Translate(+2.5, 2.0, 0.0, SeparatorFunnel(0.8, 2.0, 0.5))
            );
       // ps = Glass(ps);
        //ps = Slice(ps);
        ps = Environment(ps);
        scene.psModel = ps;
        scene.psModel = scene.psModel->Optimize();
        scene.Render(im, 4.0);
        vid.WriteFrame(im);
        printf("t = %f\n", t);

    }
    ML_StopTiming(ti);
    vid.Close();
    im.WriteBMP("RayTraceScene.bmp");
    ML_ReportTiming(ti);
}