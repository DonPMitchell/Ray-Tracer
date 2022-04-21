//
//  Rays and cameras
//  D.P. Mitchell  2020/12/22.
//
#pragma once
#include "RayTracer2020.h"

struct PerspectiveCamera : Camera {

        PerspectiveCamera(Point3 p = Point3(0.0, 0.0, 1.0), double f = 35.0) : Camera(p, f) {}
    Ray Emit(Vector2 v) { return Ray(pOrg, (Vector3((v.x-0.5)*35.0/focal, (v.y-0.5)*35.0/focal, -1.0)), 1 ); }
};

struct StereographicCamera : Camera {

        StereographicCamera(Point3 p = Point3(0.0, 0.0, 1.0), double f = 35.0) : Camera(p, f) {}
    Ray Emit(Vector2 v) {
        double d, x, y;
        x = (v.x - 0.5)*35.0/focal; 
        y = (v.y - 0.5)*35.0/focal;
        d = x*x + y*y + 4.0;
        return Ray(pOrg, (Vector3(4.0*x/d, 4.0*y/d, 1.0 - 8.0/d)), 1 );
    }
};

struct OrthographicCamera : Camera {

        OrthographicCamera(Point3 p = Point3(0.0, 0.0, 1.0), double f = 35.0) : Camera(p, f) {}
    Ray Emit(Vector2 v) {
        double x, y;
        x = (v.x - 0.5)*35.0/focal; 
        y = (v.y - 0.5)*35.0/focal;
        return Ray(Point3(x, y, pOrg.z), Vector3(0.0, 0.0, -1.0), 1 );
    }
};

struct SphericalCamera : Camera {

        SphericalCamera(Point3 p) : Camera(p, 0.0) {}
    Ray Emit(Vector2 v);
};

//
//  Finite aperture camera
//
struct ApertureCamera : StereographicCamera {
    Vector3 vX, vY;
    double  fApertureRadius;
    double  fFocusDistance;

    ApertureCamera(double fstop, double fdist, Point3 p, double focal) : StereographicCamera(p, focal), fFocusDistance(fdist)
        {   fApertureRadius = 0.001 * 0.5 * focal/fstop;    // 0.001 meters per millimeter
            vX = Vector3(fApertureRadius, 0.0, 0.0);
            vY = Vector3(0.0, fApertureRadius, 0.0);
        }

    Ray Emit(Vector2 v) {
        Ray r, r2;
        
        r = StereographicCamera::Emit(v);
        r.vDir = Normalize(r.vDir);
        r2.pOrg = r.pOrg + vX*g_uCamera + vY*g_vCamera;
        r2.vDir = r(fFocusDistance) - r2.pOrg;
        r2.bPrimary = 1;
        return r2;
    }
};
