#include "stdafx.h"
#include "RayTracer2020.h"
#pragma intrinsic(sqrt)

Radiance
PointLight::Cast(Ray &ray)
{
    Vector3 vL, vR;
    double fNL;

    vL = pOrg - ray.pOrg;
    vR = ray.vDir;
    fNL = Dot(vL, vR);
    if (fNL < 0.0) {
        return g_sBlack;       // for now, point lights are invisible
    } else {
        fNL /= sqrt(Dot(vL, vL) * Dot(vR, vR));
        if (fNL < 0.999)
            return g_sBlack;
        fNL *= fNL;
        fNL *= fNL;
        fNL *= fNL;
        fNL *= fNL;
        fNL *= fNL;
        fNL *= fNL;
        fNL *= fNL;
        fNL *= fNL;
        fNL *= fNL;
        return flux * fNL;
    }
}

 Irradiance
 AmbientLight::Illuminate(Point3 &p, CoVector3 &cvNormal, int nRoughness)
 {
    return Irradiance(flux * D_PI);      // The ambient radience integrated over the hemisphere
 }

 Radiance
 AmbientLight::Cast(Ray &ray)
{
    return Radiance(flux);
}

static inline double
root(double f)
{
    if (f < 0)
        return -sqrt(-f);
    else
        return sqrt(f);
}

 Irradiance
 HavercosineLight::Illuminate(Point3 &p, CoVector3 &cvNormal, int nRoughness)
 {
    double f, x, S;
    
    switch (nRoughness) {
    case SURF_SHINY:        
    case SURF_SMOOTH:       
    case SURF_BURNISHED:   
        S = 1.5;
        break;
    default:
        S = 1.0;
    }
    x = cvNormal.y;
    f = ((1.0 - S)*x*x + S)*x;
    f = 0.5*(1.0 - f);
    return (flux * (f*D_PI)) + (flux2 * ((1.0 - f)*D_PI));
 }