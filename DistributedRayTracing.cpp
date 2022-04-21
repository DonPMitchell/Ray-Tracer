#include "stdafx.h"
#include "RayTracer2020.h"
#include "Patterns.h"
#include "Video.h"
#include "Utilities.h"
#include "TextureMaps.h"

extern void PerturbNormal(Hit *ph, CoVector3 cvP);

double  g_tTime;
double  g_uCamera, g_vCamera;
double  g_uLight, g_vLight;
double  g_uMicroFacet, g_vMicroFacet;
/*
struct ApertureCamera : StereographicCamera {
    Vector3 vX, vY;
    double  fApertureRadius;
    double  fFocusDistance;

    ApertureCamera(Point3 p, double focal, double fstop, double fdist) : StereographicCamera(p, focal), fFocusDistance(fdist)
        {   fApertureRadius = 0.001 * 0.5 * focal/fstop;    // 0.001 meters per millimeter
            vX = Vector3(1.0, 0.0, 0.0);
            vY = Vector3(0.0, 1.0, 0.0);
        }

    Ray Emit(Vector2 v) {
        Ray r, r2;
        double a, b;
        
        r = StereographicCamera::Emit(v);
        r.vDir = Normalize(r.vDir);
        a = g_uCamera * fApertureRadius;
        b = g_vCamera * fApertureRadius;
        r2.pOrg = r.pOrg + vX*a + vY*b;
        r2.vDir = r(fFocusDistance) - r2.pOrg;
        return r2;
    }
};
*/
CoVector3
SampleNormalPerturbation(double u, double v, int kappa)
{
    double r2;

    r2 = u*u + v*v;
    return CoVector3(u, v, sqrt(1.0 - r2));
}

void
StochasticBump(const Point3 &p, Hit *ph, RayTracer *scene)
{
}
//
//  Test distributed ray tracing functions
//
static Solid*
EnvironmentLab(Solid *ps, RayTracer &scene)
{
    double f20;
    Spectrum s6500;
    DisplayRGB rgb;

    s6500 = BlackBody(6500.0);
    f20 = 0.0 * D_PI/180.0;
    scene.plAmbient = new HavercosineLight(s6500*0.002, g_sBlack);
    scene.plPointLights = new AreaLight(1.0, Point3(9.0, -6.5, 12.0), s6500 * 8000.0 * 1.2);
    //scene.plPointLights = new PointLight(Point3(9.0, -6.5, 12.0), s6500 * 8000.0 * 1.2);
    scene.plPointLights = new PointLight(Point3(10.0, -10.0, -25.0), s6500 * 1000.0, scene.plPointLights);
    ps = MakeUnion(
            ps,
            new Translate(0.0, 1.5, -10.0, new Rotate(f20, 0.0, 0.0, new Scale(16.0, 0.1, 16.0, 
                    new Material(0.9, g_sEqualWhite, MAT_PLASTIC, new Cube)))),
            new Scale (1000.0, 1000.0, 1000.0, new Material(0.9, g_sEarthSky*2.0, MAT_NOSHADE, new HollowSphere))
        );

    return ps;
}

static Solid*
Plastic(Solid *ps)
{
    Spectrum sBlue;

    sBlue = (g_sEqualWhite + g_sBlue) * 0.5;
    return new Surface(SURF_POLISHED, 0, new Material(0.7, REFRACT_PYREX, sBlue, MAT_PLASTIC, 0, ps));
}

static Solid*
Quartz(Solid *ps)
{
    return new Material(1.0, REFRACT_WATER, g_sWater, MAT_TRANSPARENT, ColoredGlass,
               new Surface(SURF_POLISHED, 0, ps));
}

static Solid*
MakeModel(double fAng)
{
    Solid *ps;

    fAng = fAng * D_PI/180.0;
    ps = MakeUnion(
      //  new Translate(+15.0, 0.0, -30.0, new Rotate(0.0, fAng, 0.0, IcosahedronCage(0.05, 0.1))),
      //  new Translate(+7.0, 0.0, -15.0, new Rotate(0.0, fAng, 0.0, IcosahedronCage(0.05, 0.1))),
        new Translate(+0.0, 0.0, 0.0, new Rotate(0.0, fAng, 0.0, IcosahedronCage(0.05, 0.1)))
       // new Translate(-1.0, 0.0, +4.0, new Rotate(0.0, fAng, 0.0, IcosahedronCage(0.05, 0.1)))
    );
    ps = Plastic(ps);
    return ps;
}

void
DistributedRayTracer()
{
    RayTracer scene;
    Image im;
    Video vid;
    ML_TimingInfo ti;
    Spectrum s6500;
    Solid *ps;
    double t;
    int nQindex[128];

    s6500 = BlackBody(6500.0);
    scene.plAmbient = new AmbientLight(s6500 * 0.02);
    scene.plPointLights = new AreaLight(1.0, Point3(9.0, -6.5, 12.0), s6500 * 10000.0);
    scene.psSample = new FastBlueSamplingTile;
    scene.psSample->QIndex1(0);
    //scene.pcCamera = new ApertureCamera(1.4, 9.0+15.0, Point3(0.0, +0.0, 9.0), 35.0); 
    scene.pcCamera = new StereographicCamera(Point3(0.0, +0.0, 9.0), 50.0); 

    im.NewImage(1920, 1080, 3);
    vid.NewVideo("Distributed.mpg", im);
    ML_StartTiming(ti);
    BuildQuasiIndex(nQindex, 2, 16);
    for (t = 0.0; t < 180.0; t += 2.0) {
        ps = MakeModel(t);
        //ps = Plastic(new Sphere);
        ps = EnvironmentLab(ps, scene);
        scene.psModel = ps;
        scene.psModel = scene.psModel->Optimize();
        scene.Render(im, 4.0);
        vid.WriteFrame(im);
        printf("t = %f\n", t);
        break;
    }
    ML_StopTiming(ti);
    vid.Close();
    im.WriteBMP("Distributed.bmp");
    ML_ReportTiming(ti);
}