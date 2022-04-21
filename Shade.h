//
//  Shade - calculate illumination of surfaces and propagation through media
//  D.P. Mitchell 2021/01/04
//
#pragma once
#include "Solids.h"
#include "RadiometryColor.h"
#include "Sampling.h"

struct Scene {
    Solid           *psModel;
    Light           *plLights;
    Camera          *pcCamera;
    SamplingTile    *psSample;

    Radiance    Shade(Hit *ph);
    Radiance    Cast(Ray &ray)              { return Shade(psModel->Intersect(ray)); }
    Radiance    ImageFunction(Vector2 &v)   { return Cast(pcCamera->Emit(v)); }
};
