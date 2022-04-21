#include "stdafx.h"
#include "RayTracer2020.h"
#include "Video.h"
#include "Utilities.h"

#define RAD(X)  ((X)*D_PI/180.0)

static Solid*
Robot(double fLeftShoulder, double fLeftElbow, double fRightShoulder, double fRightElbow,
      double fLeftHip, double fLeftKnee, double fRightHip, double fRightKnee)
{
    Solid *ps, *psTrunk, *psLeftArm, *psRightArm, *psLeftLeg, *psRightLeg;
    Spectrum sRobotOrange;

    psTrunk = new BoundingBox(MakeUnion(
        new Translate(0.0,  2.1,  0.00, new Rotate(0.0, 0.0, RAD(20.0), new Scale(0.32, 0.48, 0.32, new Sphere))),  // neck
        new Translate(0.0,  1.6,  0.00, new Scale(0.20, 0.40, 0.20, new Sphere)),                                   // chest, shoulders
        new Translate(0.0,  1.2,  0.15, new Scale(0.30, 0.40, 0.55, new Sphere)),
        new Translate(0.0,  1.2, -0.15, new Scale(0.30, 0.40, 0.55, new Sphere)),
        new Translate(0.0,  0.8,  0.00, new Scale(0.30, 0.80, 0.50, new Sphere)),                                   // abdomen, hip
        new Translate(0.0, -0.5,  0.00, new Scale(0.35, 0.35, 0.55, new Sphere)),
        new Translate(0.0,  0.34, 0.00, new Scale(0.25, 0.25, 0.35, new Sphere)),
        new Translate(0.0,  0.17, 0.00, new Scale(0.25, 0.25, 0.35, new Sphere)),
        new Translate(0.0,  0.00, 0.00, new Scale(0.25, 0.25, 0.35, new Sphere)),
        new Translate(0.0, -0.17, 0.00, new Scale(0.25, 0.25, 0.35, new Sphere))
    ));
    psLeftArm = new BoundingBox(new Translate(0.0, 1.3, 0.7, new Rotate(0.0, 0.0, RAD(fLeftShoulder), MakeUnion(
        new Translate(0.0, -1.1, 0.0, new Rotate(0.0, 0.0, RAD(fLeftElbow), MakeUnion(
            new Translate(0.0,  0.0, 0.0, new Scale(0.13, 0.20, 0.13, new Sphere)),                 // forearm, hand
            new Translate(0.0, -0.5, 0.0, new Scale(0.12, 0.50, 0.12, new Sphere)),
            new Translate(0.0, -1.0, 0.0, new Scale(0.10, 0.10, 0.10, new Sphere)),
            new Translate(0.0, -1.2, 0.0, new Scale(0.12, 0.24, 0.06, new Sphere))
        ))),
        new Translate(0.0, -0.55,  0.0, new Scale(0.23, 0.55, 0.16, new Sphere)),
        new Translate(0.0, -0.10, -0.1, new Rotate(-RAD(35.0), 0.0, 0.0, new Scale(0.22, 0.33, 0.22, new Sphere)))
    ))));
    psRightArm = new BoundingBox(new Translate(0.0, 1.3, -0.7, new Scale(1.0, 1.0, -1.0, new Rotate(0.0, 0.0, RAD(fRightShoulder), MakeUnion(
        new Translate(0.0, -1.1, 0.0, new Rotate(0.0, 0.0, RAD(fRightElbow), MakeUnion(
            new Translate(0.0,  0.0, 0.0, new Scale(0.13, 0.20, 0.13, new Sphere)),                 // forearm, hand
            new Translate(0.0, -0.5, 0.0, new Scale(0.12, 0.50, 0.12, new Sphere)),
            new Translate(0.0, -1.0, 0.0, new Scale(0.10, 0.10, 0.10, new Sphere)),
            new Translate(0.0, -1.2, 0.0, new Scale(0.12, 0.24, 0.06, new Sphere))
        ))),
        new Translate(0.0, -0.55,  0.0, new Scale(0.23, 0.55, 0.16, new Sphere)),
        new Translate(0.0, -0.10, -0.1, new Rotate(-RAD(35.0), 0.0, 0.0, new Scale(0.22, 0.33, 0.22, new Sphere)))
    )))));
    psLeftLeg = new BoundingBox(new Translate(0.0, -0.5, 0.3, new Rotate(0.0, 0.0, RAD(fLeftHip), MakeUnion(
        new Translate(0.0, -1.6, 0.0, new Rotate(0.0, 0.0, RAD(fLeftKnee), MakeUnion(               // calf, foot
            new Translate(0.0,  0.00, 0.0, new Scale(0.20, 0.20, 0.20, new Sphere)),
            new Translate(0.0, -0.75, 0.0, new Scale(0.20, 0.75, 0.20, new Sphere)),
            new Translate(0.0, -1.50, 0.0, new Scale(0.20, 0.12, 0.12, new Sphere)),
            new Translate(0.2, -1.60, 0.0, new Scale(0.45, 0.10, 0.20, new Sphere))
        ))),
        new Translate(0.0, -0.8, 0.0, new Scale(0.30, 0.80, 0.27, new Sphere))                      // thigh
    ))));
    psRightLeg = new BoundingBox(new Translate(0.0, -0.5, -0.3, new Scale(1.0, 1.0, -1.0, new Rotate(0.0, 0.0, RAD(fRightHip), MakeUnion(
        new Translate(0.0, -1.6, 0.0, new Rotate(0.0, 0.0, RAD(fRightKnee), MakeUnion(               // calf, foot
            new Translate(0.0,  0.00, 0.0, new Scale(0.20, 0.20, 0.20, new Sphere)),
            new Translate(0.0, -0.75, 0.0, new Scale(0.20, 0.75, 0.20, new Sphere)),
            new Translate(0.0, -1.50, 0.0, new Scale(0.20, 0.12, 0.12, new Sphere)),
            new Translate(0.2, -1.60, 0.0, new Scale(0.45, 0.10, 0.20, new Sphere))
        ))),
        new Translate(0.0, -0.8, 0.0, new Scale(0.30, 0.80, 0.27, new Sphere))                      // thigh
    )))));
    ps = new BoundingBox(MakeTopUnion(
        psTrunk,
        psLeftArm,
        psRightArm,
        psLeftLeg,
        psRightLeg
    ));
    sRobotOrange = g_sRed + g_sGreen*0.43 + g_sBlue*0.07;
    ps = new Material(0.9, sRobotOrange, MAT_PLASTIC, ps);
    ps = new Surface(SURF_SHINY, ps);
    return ps;
    return new Rotate(0.0, -0.5*D_PI, 0.0, new Scale(1.0, -1.0, 1.0, ps));
}

static inline double sqr(double x) { return x*x; }

static Solid*
Unicyclist(double t)
{
    double s, f, g, h, l, d, v, a, c;
    double fLeftKnee, fLeftHip, fRightKnee, fRightHip;
    double fLeftElbow, fLeftShoulder, fRightElbow, fRightShoulder;
    Solid *ps;
    //
    //  Angles for left and right legs derived kinematically from pedal angle on tire
    //
    s = D_PI * (2.0 - t);               // pedal angle (0 to 2pi)
    f = 0.8 * sin(s);                   // foot z
    g = 0.8 * cos(s);                   // foot y
    h = 16/3.0 - 1.0;                   // height of hip joint
    l = 8.0/3.0;                        // length of thigh and calf
    d = sqrt(sqr(h - g) + sqr(f));      // distance foot to hip
    v = sqrt(sqr(l) - sqr(d)/4.0);      // knee-chord distance
    a = atan(2.0*v/d);                  // internal angle
    fLeftKnee = -360.0 * a/D_PI;
    c = atan(f/(g-h));                  // external angle
    fLeftHip = 180.0 * (a+c)/D_PI;

    s = s + D_PI;                       // 180 degrees out of phase with left pedal
    f = 0.8 * sin(s);                   // foot z
    g = 0.8 * cos(s);                   // foot y
    h = 16/3.0 - 1.0;                   // height of hip joint
    l = 8.0/3.0;                        // length of thigh and calf
    d = sqrt(sqr(h - g) + sqr(f));      // distance foot to hip
    v = sqrt(sqr(l) - sqr(d)/4.0);      // knee-chord distance
    a = atan(2.0*v/d);                  // internal angle
    fRightKnee = -360.0 * a/D_PI;
    c = atan(f/(g-h));                  // external angle
    fRightHip = 180.0 * (a+c)/D_PI;

    fLeftElbow     = 15.0 + 15.0*cos(D_PI*(t + 1.0));
    fLeftShoulder  = 20.0*cos(D_PI*(t + 1.0));
    fRightElbow    = 15.0 + 15.0*cos(D_PI*t);
    fRightShoulder = 20.0*cos(D_PI*t);

    ps = Robot(fLeftShoulder, fLeftElbow, fRightShoulder, fRightElbow, fLeftHip, fLeftKnee, fRightHip, fRightKnee);
    ps = new Translate(0.0, -18.5/3, 0.0, new Scale(5.0/3.0, -5.0/3.0, -5.0/3.0, new Rotate(0.0, 0.5*D_PI, 0.0, ps)));
    ps = new TopUnion(ps, new Translate(0.0, -1.0, 0.0, new Rotate(0.0, 0.5*D_PI, 0.0, new Scale(1.27324, 1.27324, 0.1, new Cylinder))));
    ps = new Translate(0.0, -0.13662, 0.0, new Scale(0.5, 0.5, 0.5, ps));
    return ps;
}

static Solid*
eF(double t)
{
    return new BoundingBox(MakeTopUnion(
        new Translate(0.0, 32.0, 6.0, new Scale(2.0, 32.0, 8.0, new Cube)),
        new Translate(0.0, 0.0,  2.0+2.0*t, Unicyclist(t)),
        new Translate(0.0, 0.0,  6.0+2.0*t, Unicyclist(t)),
        new Translate(0.0, 0.0, 10.0+2.0*t, Unicyclist(t))
    ));
}

static Solid*
eR(double t)
{
    return new BoundingBox(MakeTopUnion(
        new Translate(0.0, 32.0, 6.0, new Scale(2.0, 32.0, 8.0, new Cube)),
        new Translate(0.0, 0.0,  2.0+2.0*t, Unicyclist(t)),
        new Translate(0.0, 0.0,  6.0+2.0*t, Unicyclist(t)),
        new Translate(2.0, 0.0, 10.0, new Rotate(0.0, 0.25*D_PI*t, 0.0, new Translate(-2.0, 0.0, 0.0, Unicyclist(t))))
    ));
}

static Solid*
eL(double t)
{
    return new BoundingBox(MakeTopUnion(
        new Translate(0.0, 32.0, 6.0, new Scale(2.0, 32.0, 8.0, new Cube)),
        new Translate(0.0, 0.0,  2.0+2.0*t, Unicyclist(t)),
        new Translate(0.0, 0.0,  6.0+2.0*t, Unicyclist(t)),
        new Translate(-2.0, 0.0, 10.0, new Rotate(0.0, -0.25*D_PI*t, 0.0, new Translate(2.0, 0.0, 0.0, Unicyclist(t))))
    ));
}

static Solid*
Clock1(double t, Solid *ps)
{
    return new BoundingBox(MakeTopUnion(
        eR(t),
        new Translate( 0.0, 0.0, 12.0, new Rotate(0.0, 0.5*D_PI, 0.0, eR(t))),
        new Translate(12.0, 0.0, 12.0, new Rotate(0.0, D_PI, 0.0, ps))
    ));
}

static Solid*
Clock1R(double t)
{
    return Clock1(t, eR(t));
}

static Solid*
Clock1F(double t)
{
    return Clock1(t, eF(t));
}

static Solid*
Clock1L(double t)
{
    return Clock1(t, eL(t));
}

static Solid*
Counter1(double t, Solid *ps)
{
    return new BoundingBox(MakeTopUnion(
        new Translate(12.0, 0.0,  0.0, eL(t)),
        new Translate(12.0, 0.0, 12.0, new Rotate(0.0, -0.5*D_PI, 0.0, eL(t))),
        new Translate( 0.0, 0.0, 12.0, new Rotate(0.0, -D_PI, 0.0, ps))
    ));
}

static Solid*
Counter1R(double t)
{
    return Counter1(t, eR(t));
}

static Solid*
Counter1F(double t)
{
    return Counter1(t, eF(t));
}

static Solid*
Counter1L(double t)
{
    return Counter1(t, eL(t));
}

static Solid*
Clockwise(double size, Solid* quad1, Solid* edge1, Solid* quad2, Solid* edge2, Solid* quad3, Solid* edge3, Solid* quad4)
{
    return new BoundingBox(MakeTopUnion(
        new Translate(0.0, 0.0, size, new Rotate(0.0, 0.5*D_PI, 0.0, quad1)),
        new Translate(0.0, 0.0, size, edge1),
        new Translate(0.0, 0.0, 12.0+size, quad2),
        new Translate(size, 0.0, 12.0+size, new Rotate(0.0, 0.5*D_PI, 0.0, edge2)),
        new Translate(12.0+size, 0.0, 12.0+size, quad3),
        new Translate(12.0+size+size, 0.0, 12.0+size, new Rotate(0.0, D_PI, 0.0, edge3)),
        new Translate(12.0+size+size, 0.0, 0.0, new Rotate(0.0, -0.5*D_PI, 0.0, quad4))
    ));
}

static Solid*
Counterwise(double size, Solid* quad1, Solid* edge1, Solid* quad2, Solid* edge2, Solid* quad3, Solid* edge3, Solid* quad4)
{
    return new BoundingBox(MakeTopUnion(
        new Translate(0.0, 0.0, size, new Rotate(0.0, 0.5*D_PI, 0.0, quad1)),
        new Translate(0.0, 0.0, 12.0+size, new Rotate(0.0, D_PI, 0.0, edge1)),
        new Translate(0.0, 0.0, 12.0+size, quad2),
        new Translate(12.0+size, 0.0, 12.0+size, new Rotate(0.0, -0.5*D_PI, 0.0, edge2)),
        new Translate(12.0+size, 0.0, 12.0+size, quad3),
        new Translate(12.0+size+size, 0.0, size, edge3),
        new Translate(12.0+size+size, 0.0, 0.0, new Rotate(0.0, -0.5*D_PI, 0.0, quad4))
    ));
}

typedef Solid *(*Prim)(double t);

static void
OddPeanos(Solid *rgs[], double t, double size, Prim cR, Prim cF, Prim cL, Prim ccR, Prim ccF, Prim ccL)
{
    rgs[0] = Clockwise(size, ccF(t), eR(t), cF(t), eF(t), cR(t), eF(t), ccR(t));
    rgs[1] = Clockwise(size, ccF(t), eR(t), cF(t), eF(t), cR(t), eF(t), ccF(t));
    rgs[2] = Clockwise(size, ccF(t), eR(t), cF(t), eF(t), cR(t), eF(t), ccL(t));
    rgs[3] = Counterwise(size, cR(t), eF(t), ccL(t), eF(t), ccF(t), eL(t), cF(t));
    rgs[4] = Counterwise(size, cF(t), eF(t), ccL(t), eF(t), ccF(t), eL(t), cF(t));
    rgs[5] = Counterwise(size, cL(t), eF(t), ccL(t), eF(t), ccF(t), eL(t), cF(t));
    rgs[6] = Counterwise(size, cF(t), eF(t), ccL(t), eF(t), ccF(t), eL(t), cF(t));  // another copy of [4]
}

static void
EvenPeanos(Solid *rgs[], double t, double size, Prim cR, Prim cF, Prim cL, Prim ccR, Prim ccF, Prim ccL)
{
    rgs[0] = Clockwise(size, ccR(t), eF(t), cL(t), eL(t), cF(t), eR(t), ccR(t));
    rgs[1] = Clockwise(size, ccR(t), eF(t), cL(t), eL(t), cF(t), eR(t), ccF(t));
    rgs[2] = Clockwise(size, ccR(t), eF(t), cL(t), eL(t), cF(t), eR(t), ccL(t));
    rgs[3] = Counterwise(size, cR(t), eL(t), ccF(t), eR(t), ccR(t), eF(t), cL(t));
    rgs[4] = Counterwise(size, cF(t), eL(t), ccF(t), eR(t), ccR(t), eF(t), cL(t));
    rgs[5] = Counterwise(size, cL(t), eL(t), ccF(t), eR(t), ccR(t), eF(t), cL(t));
    rgs[6] = Counterwise(size, cF(t), eL(t), ccF(t), eR(t), ccR(t), eF(t), cL(t));  // Another copy of [4]
}

typedef Solid *(Peanos)(Solid *rgs[6], double t, double size, Prim cR, Prim cF, Prim cL, Prim ccR, Prim ccF, Prim ccL);



static Solid*
Peano1(double t)
{
    return Clockwise(12.0, Counter1R(t), eF(t), Clock1L(t), eL(t), Clock1F(t), eR(t), Counter1F(t));
}

static Solid*
Peano2(double t)
{
    Solid *rgs[7];

    EvenPeanos(rgs, t, 12.0, Clock1R, Clock1F, Clock1L, Counter1R, Counter1F, Counter1L);
    // (clockWise size ccR eF cL eL cF eR ccF)
    //return Clockwise(12.0 + 12.0+12.0, rgs[3], eF(t), rgs[2], eL(t), rgs[1], eR(t), rgs[4]);
    // (clockWise size ccF eR cF eF cR eF ccF)
    return Clockwise(12.0 + 12.0+12.0, rgs[4], eR(t), rgs[1], eF(t), rgs[0], eF(t), rgs[6]);
}

static Solid*
Environment(Solid *ps, RayTracer &scene)
{
    Spectrum s6500;
    DisplayRGB rgb;

    s6500 = BlackBody(6500.0);
    scene.plAmbient = new HavercosineLight(g_sEarthSky*0.025, g_sEqualWhite*0.00625);
    scene.plPointLights = new PointLight(Point3(90.0, -65.0, 120.0), s6500 * 800000.0 * 1.0);
    // scene.plPointLights = new PointLight(Point3(-10.0, -10.0, -10.0), s6500 * 100.0, scene.plPointLights);
    // ps = new Union(ps, new Translate(0.0, 4.0, -10.0, new Scale(12.0, 0.1, 12.0, 
    //            new Material(0.9, g_sEqualWhite, MAT_PLASTIC, new Cube))));
    ps = new TopUnion(ps, new Scale (1000.0, 1000.0, 1000.0, new Material(0.9, g_sEarthSky*2.0, MAT_NOSHADE, new HollowSphere)));

    return ps;
}
//
//  Union    - 34.9599 sec
//  TopUnion - 20.9978
//  Swap L/R - 19.4222
//  SwapFix  - 19.1726
//  MoreTop  - 19.1882
//  (0,tMin) - 18.8450
//  Reset    - 28.4390
//
void
Megacycles()
{
    RayTracer scene;
    Image im;
    Video vid;
    ML_TimingInfo ti;
    Solid *ps;
    double t;

    scene.psSample = new DebugSamplingTile;
    scene.pcCamera = new StereographicCamera(Point3(17.0-0.5, -7.0, 172.0), 5000.0);      // 75 mm

    im.NewImage(1920, 1080, 3);
    vid.NewVideo("Megacycles.mpg", im);
    ML_StartTiming(ti);
    for (t = 0.0; t < 2.0; t += 1.0/60.0) {
        ps = Peano2(t);
        ps = new Rotate(0.0, (90 + 30.0)*D_PI/180.0, 0.0, ps);
        ps = new Rotate(-15.0*D_PI/180.0, 0.0, 0.0, ps);
        ps = Environment(ps, scene);
        scene.psModel = ps;
        scene.psModel = scene.psModel->Optimize();
        scene.Render(im, 4.0);
        vid.WriteFrame(im);
        printf("t = %f\n", t);
        break;
    }
    ML_StopTiming(ti);
    vid.Close();
    im.WriteBMP("Megacycles.bmp");
    ML_ReportTiming(ti);
}