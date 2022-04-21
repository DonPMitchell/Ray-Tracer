//
//  Fundamental classes for ray tracer
//  D.P. Mitchell  2021/01/06.
//
#pragma once
#include "Affine.h"
#include "Image.h"
#include "RadiometryColor.h"
#include "Patterns.h"
#include "Texture.h"
#include "Complex.h"

#define D_PI        3.14159265358979323846264338327950288419716939937510
#define D_SQRT2     1.41421356237309504880168872420969807856967187537694
#define D_SQRT3     1.73205080756887729352744634150587236694280525381038
#define D_PHI       1.61803398874989484820458683436563811772030917980576
#define EPSILON     2.0e-7
#define INFINITY    1.0e+37

#define VERBOSE if(g_nVerbose)
extern int  g_nVerbose;

//
//  Ray is a parametric curve with semi-infinite domain (0, inf] or line segment (0, tMin]
//
struct Ray {
    Point3  pOrg;
    Vector3 vDir;
    double  tMin;
    int     bPrimary;

            Ray() {}
            Ray(Point3 &p, Vector3 &v, int nCameraRay=0) : pOrg(p), vDir(v), tMin(INFINITY), bPrimary(nCameraRay) {}
    Point3  operator()(double t) const { return pOrg + vDir*t; }
    void    Print() const { printf("{"); pOrg.Print(); printf("->"); vDir.Print(); printf(" [0,%g]}\n", tMin); }
};
//
//  Camera is the first stage of the image function, mapping the unit square to ray: [0,1]x[0,1] -> (Org, Dir)
//      Possible mappings: Perspective, Stereographic, equiangular, equal-area, orthograhic
//
struct Camera {
    Point3  pOrg;
    double  focal;

            Camera() {}
            Camera(Point3 p, double f = 35.0) : pOrg(p), focal(f) {}
    virtual Ray Emit(Vector2 v) = 0;            // normalized vector
};
//
//  Axis-aligned bounding box with fast ray-hit test
//
struct Extent {
    Point3  pMin, pMax;

            Extent() {}
            Extent(double w) { pMin = Point3(-w, -w, -w); pMax = Point3(w, w, w); }
            Extent(const Point3 &p1, const Point3 &p2) : pMin(p1), pMax(p2) {}
    int     NullExtent();
    int     TestIntersection(Ray &ray);
    void    Print() { printf("[(%f %f %f), (%f %f %f)]\n", pMin.x, pMin.y, pMin.z, pMax.x, pMax.y, pMax.z); }
};
//
//  CSG nodes, representing primitive solids, transformed solids, generalized set operations, volumetric and surface shading
//
struct Hit;
struct RayTracer;

struct Solid {
    int     nIsAffine;

            Solid() : nIsAffine(0) {}
    virtual Hit    *Intersect(Ray &r) = 0;
    virtual Extent  GetExtent();                    // bounding box for fast ray test
    virtual Solid  *Optimize() { return this; }   // traverse model, perform optimizations
};
//
//  Hit list returned by ray/solid intersection.  Surface hit points, normal vectors, shading properties
//
struct AffineSolid;
struct Surface;
struct Material;
struct Instance;
struct TextureList;

#define PNEWSIZE            (8192 * 10)
#define PARENT_AFFINE       0
#define PARENT_MATERIAL     1
#define PARENT_SURFACE      2
#define PARENT_INSTTANCE    3

struct Parent {
    union {
        AffineSolid *pa;
        Material    *pm;
        Surface     *ps;
    };
    Parent  *pParent;
    int     nType;
    static  Parent rgp[];
    static  Parent *pNew;

            Parent() {}
            Parent(AffineSolid *p, Parent *pp) : pa(p), pParent(pp), nType(PARENT_AFFINE)   {}
            Parent(Material *p, Parent *pp)    : pm(p), pParent(pp), nType(PARENT_MATERIAL) {}
            Parent(Surface *p, Parent *pp)     : ps(p), pParent(pp), nType(PARENT_SURFACE)  {}
};

struct Hit {
    double      t;
    CoVector3   cvNormal;
    CoVector3   cvPerturbed;        // normal vector perturbed by bump map during intersect
    Point3      pTexture;           // coordinates for Texture call-out computed during intersect.  Replace with TexturePoint

    Surface     *ps;
    Material    *pm;
    TextureList *pt;
    AffineSolid *paAffineTransforms;
    Hit         *phNext;
    static Hit  *phFree;
    int         nPerturbed;

            Hit() {}
            Hit(Material *pMat) : t(0.0), cvNormal(0.0, 0.0, 0.0), ps(0), pm(pMat), pt(0), paAffineTransforms(0), nPerturbed(0), phNext(0) {}
            Hit(Hit *ph) : ps(0), pm(0), pt(0), phNext(ph) {}
    void    Print() { 
                Hit *ph;
                printf("{ ");
                for (ph = this; ph; ph = ph->phNext)
                   printf("t=%g, (%0.2f %0.2f %0.2f)", ph->t, ph->cvNormal.x, ph->cvNormal.y, ph->cvNormal.z);
                printf("}\n");
            }
};

extern void DeleteHitList(Hit *ph);
extern int  LengthHitList(Hit *ph);
extern void PrintHitList(char *sz, Hit *ph);
extern int  g_nHitCount, g_nNewHits;
extern int  g_nTextureListCount, g_nTextureLists;

//
//  Light sources (point, area, ambient)
//
struct Light {
    Flux    flux;

            Light() {}
            Light(DisplayRGB rgb) { flux = RGBtoSpectrum(rgb); }
            Light(Spectrum &s) { flux = s; }
    virtual Radiance Cast(Ray &ray) = 0;
    virtual Point3 Origin() = 0;
};

struct PointLight : Light {
    Point3      pOrg;
    PointLight  *plNext;
    int         nNoShadow;

                PointLight() {}
                PointLight(Point3 &p, double fPower, PointLight *pl = 0) : plNext(pl), pOrg(p), Light(fPower), nNoShadow(0) {}
                PointLight(Point3 &p, Spectrum sPower, PointLight *pl = 0) : plNext(pl), pOrg(p), Light(sPower), nNoShadow(0) {} 
    Radiance    Cast(Ray &ray);
    Point3      Origin() { return pOrg; }
};

struct AmbientLight : Light {

                AmbientLight(double fBackground) : Light(fBackground) {}
                AmbientLight(Spectrum sBackground) : Light(sBackground) {}
    Irradiance  Illuminate(Point3 &p, CoVector3 &cvNormal, int nRoughness);
    Radiance    Cast(Ray &ray);
    Point3      Origin() { return Point3(0.0, 0.0, 0.0); }
};
//
//  Sample an image function (e.g. ray tracer) to generate a digital image
//
struct SamplingTile {
    int     nSamples;
    int     nVenera;

            SamplingTile() {}
    int     virtual Initialize() = 0;
    int     virtual Render(Image &im, Radiance (*f)(Vector2 v), double fSamplesPerPixel) = 0;   // parallelize here
    void    virtual QIndex1(int nPime)  {}
};
//
//  Scene structure groups model, lights, camera and tiling
//
struct RayTracer {
    Solid           *psModel;
    AmbientLight    *plAmbient;
    PointLight      *plPointLights;         // A unniversial interface for all light sources is too problematic
    Camera          *pcCamera;
    SamplingTile    *psSample;
    int             nRayLevel;
    double          fAttenuation;
    
                RayTracer() : psModel(0), plAmbient(0), plPointLights(0), pcCamera(0), psSample(0), nRayLevel(0), fAttenuation(1.0) {}
    //Radiance    Shade(Hit *ph, Ray &r);
    Irradiance  LightRay(Ray &ray, const PointLight *pl);
    Radiance    VisualRay(Ray &ray);
    int         Render(Image &im, double fSamplesPerPixel = 1.0);
};
//
//  distributed-ray-tracing parameters
//
extern double  g_tTime;
extern double  g_uCamera, g_vCamera;            // [-1, 1]x[-1, 1]
extern double  g_uLight, g_vLight;              // [-1, 1]x[-1, 1]
extern double  g_uMicroFacet, g_vMicroFacet;
//
//  debugging
//
extern int g_nTestRow, g_nTestColumn;
extern int IsNiceNumber(double f);
extern int BadRay(Ray &ray);
extern void Abort();

#include "Solids.h"
#include "Lights.h"
#include "Sampling.h"
#include "Cameras.h"

