//
//  Light sources
//  D.P.Mitchell  2021/01/06.
//
#pragma once

//
//  Haversine ambient illumination from two colored hemispheres.  0.5*(1 - sine)
//
struct HavercosineLight : AmbientLight {
    Flux flux2;

                HavercosineLight(const Spectrum &sAbove, const Spectrum &sBelow) : AmbientLight(sAbove), flux2(sBelow) {}
    Irradiance  Illuminate(Point3 &p, CoVector3 &cvNormal, int nRoughness);
};
//
//  Area light source
//PointLight(Point3 &p, Spectrum sPower, PointLight *pl = 0) : plNext(pl), pOrg(p), Light(sPower) {}
//
struct AreaLight : PointLight {
    Vector3     vA, vB;

            AreaLight(double fDiam, Point3 &p, Spectrum sPower, PointLight *pl = 0) : PointLight(p, sPower, pl) {
                Vector3 vN;
                vN = Point3(0.0 , 0.0, 0.0) - p;
                vA = Normalize(Cross(vN, Vector3(0.0, 1.0, 0.0))) * fDiam;
                vB = Normalize(Cross(vN, vA)) * fDiam;
            }
    Point3  Origin() { return pOrg + vA*g_uLight + vB*g_vLight; }
};