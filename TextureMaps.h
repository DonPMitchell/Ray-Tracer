//
//  Procedural texture and light propagators
//  D.P. Mitchell 2021/03/16.
//
#pragma once
#include "RayTracer2020.h"

Propagation ColoredGlass(const Point3 &p0, const Point3 &p1, Material *pm, RayTracer *scene);

extern void Bumpy(const Point3 &p, Hit *ph, RayTracer *scene);
extern void MicroFacets(const Point3 &p, Hit *ph, RayTracer *scene);
extern void MicroFacets2(const Point3 &p, Hit *ph, RayTracer *scene);
extern void Wrinkled(const Point3 &p, Hit *ph, RayTracer *scene);
extern void Faceted20(const Point3 &p, Hit *ph, RayTracer *scene);
extern void Dented(const Point3 &p, Hit *ph, RayTracer *scene);
extern void Wavy(const Point3 &p, Hit *ph, RayTracer *scene);
extern void Fibers(const Point3 &p, Hit *ph, RayTracer *scene);
extern void Blistered(const Point3 &p, Hit *ph, RayTracer *scene);
extern void Hammered(const Point3 &p, Hit *ph, RayTracer *scene);
extern void Tectonic(const Point3 &p, Hit *ph, RayTracer *scene);
extern void SeaScape(const Point3 &p, Hit *ph, RayTracer *scene);
extern void VoronoiBump(const Point3 &p, Hit *ph, RayTracer *scene);
extern void VoronoiAntiBump(const Point3 &p, Hit *ph, RayTracer *scene);
extern void VoronoiCrinkled(const Point3 &p, Hit *ph, RayTracer *scene);
extern void QuasiMicroFacet(const Point3 &p, Hit *ph, RayTracer *scene);

extern void Voronoi(const Point3 &p, Hit *ph, RayTracer *scene);
extern void Agate(const Point3 &p, Hit *ph, RayTracer *scene);
extern void Marble(const Point3 &p, Hit *ph, RayTracer *scene);
extern void Clouds(const Point3 &p, Hit *ph, RayTracer *scene);
extern void Wood(const Point3 &p, Hit *ph, RayTracer *scene);
extern void Spots(const Point3 &p, Hit *ph, RayTracer *scene);
extern void Veins(const Point3 &p, Hit *ph, RayTracer *scene);

extern void EarthTexture(const Point3 &p, Hit *ph, RayTracer *scene);
extern void Earth2Texture(const Point3 &p, Hit *ph, RayTracer *scene);
extern void MarsTexture(const Point3 &p, Hit *ph, RayTracer *scene);
extern void MoonTexture(const Point3 &p, Hit *ph, RayTracer *scene);
extern void VenusTexture(const Point3 &p, Hit *ph, RayTracer *scene);
extern void MercuryTexture(const Point3 &p, Hit *ph, RayTracer *scene);
extern void SkyTexture(const Point3 &p, Hit *ph, RayTracer *scene);
