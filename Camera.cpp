#include "stdafx.h"
#include "RayTracer2020.h"
#pragma intrinsic(sin, cos, exp, atan)

static Vector3
InverseStereographic(Vector2 v)
{
    double d;

    d = 1.0 + Dot(v, v);
    return Vector3(2.0*v.x/d, 2.0*v.y/d, (Dot(v,v) - 1.0)/d);
}

static Vector3
InverseMercator(Vector2 v)
{
    double fLong, fLat, x, y, z, r;

    fLong = 2.0*D_PI*v.x;
    fLat  = 2.0*atan(exp(v.y)) - 0.5*D_PI;
    z = sin(fLat);
    r = cos(fLat);
    x = cos(fLong) * r;                        // x axis, longitude == 0.0
    y = sin(fLong) * r;
    return Vector3(x, y, z);
}

Ray
SphericalCamera::Emit(Vector2 v)
{
    Ray r;

    if (v.y > 0.5) {
        if (v.x > 0.5) {
            
        } else {
        }
    } else {
    }
    r = Ray(Point3(0.0, 0.0, 0.0), Vector3(0.0, 0.0, 1.0));     // place holder
    return r;
}