#include "stdafx.h"
#include "RayTracer2020.h"
#include "TextureMaps.h"
#include "Video.h"
#include "MonteCarlo.h"
#include "Utilities.h"

//
//  Test RayTracer object
//
#define RNDSPEC (RGBtoSpectrum(DisplayRGB(RandomDouble(), RandomDouble(), RandomDouble())))

static Spectrum s_rgs[20];
static double s_fRad;

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
    return new Material(0.33, REFRACT_WATER, sBlue, MAT_PLASTIC, 0,
               new Surface(SURF_POLISHED, 0, ps));
}

static Solid*
Plastic2(Solid *ps)
{
    Spectrum sBlue;

    sBlue = (g_sBlue + g_sEqualWhite) * 0.5;
    return new Material(0.33, REFRACT_WATER, sBlue, MAT_PLASTIC, 0,
               new Surface(SURF_SHINY, 0, ps));
}


static Solid*
MakeModel(double fAng)                              // geometric primitive types
{
    double fRad;
    Solid *rgps[9], *ps;
    int i;

    ps = 0;

        fRad = fAng*D_PI/180.0;
        rgps[0] = new Translate(-3.0, -3.0, 0.0, new Rotate(0.0, fRad, 0.0, new Cylinder));
        rgps[1] = new Translate(-3.0, 0.0, 0.0, new Rotate(0.0, fRad, 0.0, new Torus(0.5)));
        rgps[2] = new Translate(-3.0, 3.0, 0.0, new Rotate(0.0, fRad, 0.0, new Intersection(new Sphere, new Hyperboloid2(0.125))));

        rgps[3] = new Translate(0.0, -3.0, 0.0, new Rotate(0.0, fRad, 0.0, new Intersection(new Sphere, new Hyperboloid1(0.125))));
        rgps[4] = new Translate(0.0, 0.0, 0.0, new Rotate(0.0, fRad, 0.0, new Difference(new Icosahedron, new Cone)));
        rgps[5] = new Translate(0.0, 3.0, 0.0, new Rotate(0.0, fRad, 0.0, new Dodecahedron));

        rgps[6] = new Translate(3.0, -3.0, 0.0, new Rotate(0.0, fRad, 0.0, new Scale(1.2, 0.8, 0.8, new Sphere)));
        rgps[7] = new Translate(3.0, 0.0, 0.0, new Rotate(0.0, fRad, 0.0, new Union(new Torus(0.2), new Sphere)));
        rgps[8] = new Translate(3.0, 3.0, 0.0, new Rotate(0.0, fRad, 0.0, new Difference(new Cube,
                                                                               new Scale(0.5, 0.5, 1.2, new Cylinder))));
        for (i = 0; i < 9; i++)
            rgps[i] = new Material(0.9, s_rgs[i], MAT_METAL, rgps[i]);
        ps = rgps[0];
        for (i = 1; i < 9; i++) {
            ps = new Union(rgps[i], ps);
        }
        ps = new Surface(SURF_BURNISHED, ps);
        //ps = new Translate(3.0, 0.0, 0.0, new Rotate(0.0, fRad, 0.0, new Torus(0.2)));
        ps = new Union(ps, new Translate(0.0, 5.25, -10.0, new Scale(12.0, 1.0, 12.0, new Cube)));

    return ps;
}

static Solid*
MakeModel2(double fAng)                             // tori with various materials and surfaces
{
    double fRad, x, y;
    int i, n, iSurface, iMaterial;
    Solid *rgps[10], *ps;

    fRad = fAng*D_PI/180.0;
    ps = 0;
    i = 0;
    x = -6.0;
    for (iSurface = SURF_MATTE; iSurface <= SURF_POLISHED; iSurface++) {
        y = -2.0;
        for (iMaterial = MAT_PLASTIC; iMaterial <= MAT_METAL; iMaterial++) {
            rgps[i] = new Surface(iSurface, new Material(0.9, s_rgs[i], iMaterial, 
                new Translate(x, y, 0.0, new Rotate(0.0, fRad, 0.0, new Torus(0.5))) ));
            i++;
            y += 4.0;
        }
        x += 3.0;
    }
    n = i;
    ps = rgps[0];
    for (i = 1; i < n; i++) {
        ps = new Union(rgps[i], ps);
    }
    ps = new Union(ps, new Translate(0.0, 5.25, -10.0, new Scale(12.0, 1.0, 12.0, new Cube)));
    //return rgps[9];
    return ps;
}

static Solid*
MakeModel3(double fAng)                             // rotating textured spheres, camera 9
{
    double fRad;
    Spectrum sBlue;
    Solid *ps, *psUnion;

    fRad = fAng*D_PI/180.0;
    sBlue = g_sBlue * 0.5 + g_sWhite * 0.5;
    ps = new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC,
               new Surface(SURF_SHINY, Bumpy, new Sphere)));
    psUnion = new Translate(-3.0, 0.0, 0.0, ps);
    ps = new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC,
               new Surface(SURF_SHINY, Wrinkled, new Sphere)));
    psUnion = new Union(new Translate(0.0, 0.0, 0.0, ps), psUnion);
    ps = new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC,
               new Surface(SURF_SHINY, Dented, new Sphere)));
    psUnion = new Union(new Translate(3.0, 0.0, 0.0, ps), psUnion);
    return psUnion;
}

#define BUMPMAP Bumpy
#define SURF SURF_BURNISHED 

static Solid*
MakeModel4(double x)                             // differences of objects with materials and surfaces
{
    Spectrum sBlue, sOrange;
    Solid *ps, *psDiff, *psUnion;

    sBlue = g_sBlue * 0.5 + g_sWhite * 0.5;
    sOrange = g_sRed * 0.5 + g_sGreen * 0.25 + g_sWhite * 0.25;
    if (x > 60.0) {
        x = -1.0 + 2.5*(x - 60.0)/60.0;
        ps =  new Material(0.9, sOrange, MAT_PLASTIC,
                   new Surface(SURF_MATTE, 0, new Sphere));
        psDiff = new Translate(-1.5, 0.0, 0.0, ps);
        ps = new Material(0.9, sBlue, MAT_PLASTIC,
                   new Surface(SURF, BUMPMAP, new Sphere));                                      // bumpy
        psDiff = new Difference(psDiff, new Translate(-1.0, 0.0, 0.0, ps));
        psUnion = new Union(psDiff, new Translate(x, 0.0, 0.0, new Material(0.9, sBlue, MAT_PLASTIC,
                                        new Surface(SURF, BUMPMAP, new Sphere))));               // bumpy
    } else {
        x = 1.5 - 2.5*x/60.0;
        ps =  new Material(0.9, sOrange, MAT_PLASTIC,
                   new Surface(SURF_MATTE, 0, new Sphere));
        psDiff = new Translate(-1.5, 0.0, 0.0, ps);
        psUnion = new Union(psDiff, new Translate(x, 0.0, 0.0, new Material(0.9, sBlue, MAT_PLASTIC,
                                        new Surface(SURF, BUMPMAP, new Sphere))));               // bumpy
    }
    return psUnion;
}

static Solid*
MakeModel5(double fAng)                             // rotating textured spheres, camera 9
{
    double fRad, fRad2;
    Spectrum sBlue;
    Solid *ps, *psUnion;

    fRad2 = fAng*D_PI/180.0;
    fRad = 0.0;
    sBlue = g_sBlue * 0.5 + g_sWhite * 0.5;
    ps = new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, g_sGold, MAT_METAL,
               new Surface(SURF_POLISHED, 0, new Sphere)));
    psUnion = new Translate(-3.0, 0.0, -3.0, ps);
    ps = new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, g_sCopper, MAT_METAL,
               new Surface(SURF_POLISHED, Bumpy, new Sphere)));
    psUnion = new Union(new Translate(3.0, 0.0, -3.0, ps), psUnion);
    ps = new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_TRANSPARENT,
               new Surface(SURF_POLISHED, 0, new Sphere)));
    psUnion = new Union(new Translate(0.0, 0.0, 0.0, ps), psUnion);
    ps = new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC,
               new Surface(SURF_POLISHED, 0, new Sphere)));
    psUnion = new Rotate(0.0, fRad2, 0.0, new Union(new Translate(0.0, 0.0, -3.0, ps), psUnion));
    psUnion = new Union(psUnion, new Translate(0.0, 3.0, -10.0, new Scale(12.0, 1.0, 12.0, new Cube)));
    return psUnion;
}

static Solid*
MakeModel6(double fAng)             // water cubes with glass, copper, air spheres
{
    double fRad, fPitch;
    Solid *ps, *psUnion;

    fRad = fAng*D_PI/180.0;
    fPitch = -20.0 * D_PI/180.0;
    ps = new Rotate(0.0, fRad, 0.0, 
            new Material(0.33, REFRACT_WATER, g_sWater, MAT_TRANSPARENT, ColoredGlass,
               new Surface(SURF_POLISHED, 0, new Cube)));
    psUnion = new Translate(-3.0, 0.0, 0.0, ps);
    ps = new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, g_sCopper, MAT_METAL,
               new Surface(SURF_POLISHED, 0, new Scale(0.5, 0.5, 0.5, new Sphere))));
    psUnion = new Union(psUnion, new Translate(0.0, 0.0, 0.0, ps));

    ps = new Rotate(0.0, fRad, 0.0, 
            new Material(0.33, REFRACT_WATER, g_sWater, MAT_TRANSPARENT, ColoredGlass,
               new Surface(SURF_POLISHED, 0, new Cube)));
    psUnion = new Union (psUnion, new Translate(0.0, 0.0, 0.0, ps));
    ps = new Rotate(0.0, fRad, 0.0, 
            new Material(0.5, REFRACT_PYREX, g_sEqualWhite, MAT_TRANSPARENT, 0,
               new Surface(SURF_POLISHED, 0, new Scale(0.5, 0.5, 0.5, new Sphere))));
    psUnion = new Union(new Translate(-3.0, 0.0, 0.0, ps), psUnion);
    
    ps = new Rotate(0.0, fRad, 0.0, 
            new Material(0.33, REFRACT_WATER, g_sWater, MAT_TRANSPARENT, ColoredGlass,
               new Surface(SURF_POLISHED, 0, new Cube)));
    psUnion = new Union (psUnion, new Translate(3.0, 0.0, 0.0, ps));
    ps = new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, REFRACT_AIR, g_sEqualWhite, MAT_TRANSPARENT, 0,
               new Surface(SURF_POLISHED, 0, new Scale(0.5, 0.5, 0.5, new Sphere))));
    psUnion = new Union(psUnion, new Translate(3.0, 0.0, 0.0, ps));
   
    psUnion = new Rotate(fPitch, 0.0, 0.0, new Union(psUnion, new Translate(0.0, 3.0, -10.0, new Scale(12.0, 1.0, 12.0, new Cube))));
    return psUnion;
}

static Solid*
MakeModel7(double fAng)
{
    Solid *ps;

    fAng = D_PI * fAng / 180.0;
    ps = new Surface(SURF_POLISHED, new Material(0.5, g_sWhite, MAT_TRANSPARENT, new Rotate(0.0, fAng, 0.0, new Cube)));
    ps = new Union(ps, new Translate(0.0, 3.0, -10.0, new Scale(12.0, 1.0, 12.0, new Cube)));
    return ps;
}

static Solid*
MakeModel8(double fAng)                             // 3 spheres, camera 9  glass test
{
    double fRad, fPitch;
    Spectrum sBlue;
    Solid *ps, *psUnion;

    fRad = fAng*D_PI/180.0;
    fPitch = -10.0 * D_PI/180.0;
    sBlue = g_sBlue * 0.5 + g_sEqualWhite * 0.5;
    ps = new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC,
               new Surface(SURF_SHINY, 0, new Sphere)));
    psUnion = new Translate(-3.0, 0.0, 0.0, ps);

    ps = new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, REFRACT_WATER, g_sWater, MAT_TRANSPARENT, ColoredGlass,
               new Surface(SURF_POLISHED, Bumpy, new Scale(1.0, 1.0, 1.0, new Sphere))));

    psUnion = new Union(new Translate(0.0, 0.0, 0.0, ps), psUnion);
    ps = new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC,
               new Surface(SURF_SHINY, 0, new Sphere)));
    psUnion = new Union(new Translate(3.0, 0.0, 0.0, ps), psUnion);
    ps = new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC,
               new Surface(SURF_SHINY, 0, new Sphere)));
    psUnion = new Union(new Translate(-3.0, 0.0, -6.0, ps), psUnion);
    ps = new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC,
               new Surface(SURF_SHINY, 0, new Sphere)));
    psUnion = new Union(new Translate(0.0, 0.0, -6.0, ps), psUnion);
    ps = new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC,
               new Surface(SURF_SHINY, 0, new Sphere)));
    psUnion = new Union(new Translate(3.0, 0.0, -6.0, ps), psUnion);
    psUnion = new Rotate(fPitch, 0.0, 0.0, new Union(psUnion, new Translate(0.0, 3.0, -10.0, new Scale(12.0, 1.0, 12.0, new Cube))));
    return psUnion;
}

static Solid*
MakeModel9(double fAng)             // glass of whiskey & ice
{
    double fRad, fPitch;
    Solid *psGlass, *psWhiskey, *psIce, *psUnion;

    fRad = fAng*D_PI/180.0;
    fPitch = 90.0 * D_PI/180.0;
    psGlass = new Material(DENSE_PYREX, g_sBlack, MAT_TRANSPARENT, new Surface(SURF_POLISHED, 0, 
                new Union (new Difference(new Cylinder, new Translate(0.0, 0.0, 0.2, new Scale(0.9, 0.9, 1.0, new Cylinder))),
                           new Translate(0.0, 0.0, 1.0, new Torus(0.06)) )
               // new Union(new Translate(0.0, 0.0, 1.0, new Scale(0.625, 0.625, 0.625, new Torus(0.07))),  
               // new Intersection(new Cube, new Scale(1.0, 1.0, 4.0, new Difference(new Cone, new Scale(0.8, 0.8, 1.0, new Cone)))) )
                ));
    psWhiskey = new Material(DENSE_WATER, REFRACT_WATER, g_sWhiskey, MAT_TRANSPARENT, ColoredGlass, new Surface(SURF_POLISHED, 0, 
                  new Translate(0.0, 0.0, 0.0, new Scale(0.95, 0.95, 0.95, new Cylinder)) ));
    psIce = new Material(DENSE_ICE, REFRACT_ICE, g_sWater, MAT_TRANSPARENT, ColoredGlass, new Surface(SURF_POLISHED, Bumpy, 
               new Union(new Translate(0.1, 0.1, 0.4, new Scale(0.5, 0.5, 0.5, new Sphere)),
                         new Translate(-0.1, 0., -0.4, new Scale(0.5, 0.5, 0.5, new Sphere)) ) ));
    psUnion = psGlass;
    //psUnion = new Union(psUnion, new Translate(0.0, 0.0, 0.0, psWhiskey));
    psUnion = new Union(psUnion, new Translate(0.0, 0.0, 0.0, psIce));
    psUnion = new Rotate(fPitch, 0.0, 0.0, psUnion);
    psUnion = new Scale(0.1, 0.1, 0.1, psUnion);
    psUnion = new Difference(psUnion, new Translate(0.0, 0.0, 1.0, new Surface(SURF_POLISHED, 0, new Cube)));
    psUnion = new Union(psUnion, new Translate(0.0, 3.0, -10.0, new Scale(12.0, 1.0, 12.0, new Cube)));
    //psUnion = new Rotate(fPitch, 0.0, 0.0, psUnion);
    return psUnion;
}

static Solid*
MakeModelA(double fAng)             // various spheres moving through glass & water cubes
{
    double fPitch, fYaw;
    Solid *ps, *psUnion;

    psUnion = new Translate(0.0, 3.0, -10.0, new Scale(12.0, 1.0, 12.0, new Cube));

    ps = new Material(0.5, REFRACT_PYREX, g_sBlack, MAT_TRANSPARENT, 0, new Surface(SURF_POLISHED, 0, new Cube));
    psUnion = new Union(new Translate(2.0, 0.0, 0.0, ps), psUnion);
    ps = new Material(0.5, REFRACT_WATER, g_sWater, MAT_TRANSPARENT, ColoredGlass, new Surface(SURF_POLISHED, 0, new Cube));
    psUnion = new Union(new Translate(-2.0, 0.0, 0.0, ps), psUnion);

    ps =  new Material(0.9, g_sCopper, MAT_METAL, new Surface(SURF_POLISHED, 0, new Scale(0.5, 0.5, 0.5, new Sphere)));
    psUnion = new Union(new Translate(fmod(fAng/10 + 0.0, 8.0) - 4.0,  0.0, 0.0, ps), psUnion);
    ps =  new Material(0.9, REFRACT_PYREX, g_sBlack, MAT_TRANSPARENT, 0, new Surface(SURF_POLISHED, 0, new Scale(0.5, 0.5, 0.5, new Sphere)));
    psUnion = new Union(new Translate(fmod(fAng/10 + 2.0, 8.0) - 4.0,  0.0, 0.0, ps), psUnion);
    ps =  new Material(0.9, REFRACT_WATER, g_sWater, MAT_TRANSPARENT, ColoredGlass, new Surface(SURF_POLISHED, 0, new Scale(0.5, 0.5, 0.5, new Sphere)));
    psUnion = new Union(new Translate(fmod(fAng/10 + 4.0, 8.0) - 4.0,  0.0, 0.0, ps), psUnion);
    ps =  new Material(0.9, REFRACT_AIR, g_sBlack, MAT_TRANSPARENT, 0, new Surface(SURF_POLISHED, 0, new Scale(0.5, 0.5, 0.5, new Sphere)));
    psUnion = new Union(new Translate(fmod(fAng/10 + 6.0, 8.0) - 4.0,  0.0, 0.0, ps), psUnion);

    fPitch = -40.0 * D_PI/180.0;
    fYaw = 20.0 * D_PI/180.0;
    psUnion = new Rotate(fPitch, fYaw, 0.0, psUnion);
    return psUnion;
}

static Solid*
MakeModelB(double fAng)
{
    Solid *ps;

    //ps =  new Rotate(D_PI/2.0, 0.0, 0.0, new Intersection(new ConeHyperCylinder(fAng), new Slab));
    //ps = new Union(ps, new Translate(3.0, 0.0, 0.0, new Sphere));
    ps = new Translate(3.0, 0.0, 0.0, new Sphere);
    return ps;
}
//
//  19.48 - nobounds
//  14.49 - bounding boxes
//
static Solid*                               // Retort
MakeModelC(double fAng)
{
    Solid *ps, *ps2, *psDiff, *psCone;
    Solid *psStand;
    double f20, f90, fRad, f;

    f20 = 20.0 * D_PI/180.0;
    f90 = 90.0 * D_PI/180.0;
    fRad = (fAng - 90.0) * D_PI/180.0;
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

    psDiff = new Rotate(0.0, fRad, 0.0, new Difference(ps, ps2));
    psDiff = new Material(0.33, REFRACT_WATER, g_sWater, MAT_TRANSPARENT, ColoredGlass,
               new Surface(SURF_POLISHED, 0, psDiff));
    psDiff = new BoundingBox(psDiff);

    f = 1.05/sqrt(2.0);
    psStand = new Torus(0.05/f);
    psStand = new Union(psStand, new Translate(1.0, 0.0, -0.5, new Scale(0.05, 0.05, 0.5, new Cylinder)));
    psStand = new Union(psStand, new Translate(-0.5, +0.866, -0.5, new Scale(0.05, 0.05, 0.5, new Cylinder)));
    psStand = new Union(psStand, new Translate(-0.5, -0.866, -0.5, new Scale(0.05, 0.05, 0.5, new Cylinder)));
    psStand = new Translate(0.0, f, 0.0, new Rotate(0.5*D_PI, 0.0, 0.0, new Scale(f, f, f, psStand)));
    psStand = new Material(0.9, g_sCopper, MAT_METAL,
               new Surface(SURF_POLISHED, 0, psStand));
    psStand = new BoundingBox(psStand);
    //psStand->GetExtent().Print();

    ps = new Union(psStand, psDiff);
    ps = new Union(ps, new Translate(0.0, 2.5, -10.0, new Scale(12.0, 1.0, 12.0, new Cube)));
    ps = new Rotate(-f20*0.7, 0.0, 0.0, ps);
    ps = new Union(ps, new Scale (100.0, 100.0, 100.0, new HollowSphere));
    return ps;
}

#define PRISM  Rotate(0.5*D_PI, 0.0, 0.0, new Prism(17))

static Solid*
MakeModelD(double fAng)             // water prisms  with glass, copper, air spheres
{
    double fRad, fPitch;
    Solid *ps, *psUnion;

    fRad = fAng*D_PI/180.0;
    fPitch = -20.0 * D_PI/180.0;
    ps = new Rotate(0.0, fRad, 0.0, 
            new Material(0.33, REFRACT_WATER, g_sWater, MAT_TRANSPARENT, ColoredGlass,
               new Surface(SURF_POLISHED, 0, new PRISM)));
    psUnion = new Translate(-3.0, 0.0, 0.0, ps);
    ps = new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, g_sCopper, MAT_METAL,
               new Surface(SURF_POLISHED, 0, new Scale(0.5, 0.5, 0.5, new Sphere))));
    psUnion = new Union(psUnion, new Translate(0.0, 0.0, 0.0, ps));

    ps = new Rotate(0.0, fRad, 0.0, 
            new Material(0.33, REFRACT_WATER, g_sWater, MAT_TRANSPARENT, ColoredGlass,
               new Surface(SURF_POLISHED, 0, new PRISM)));
    psUnion = new Union (psUnion, new Translate(0.0, 0.0, 0.0, ps));
    ps = new Rotate(0.0, fRad, 0.0, 
            new Material(0.5, REFRACT_PYREX, g_sEqualWhite, MAT_TRANSPARENT, 0,
               new Surface(SURF_POLISHED, 0, new Scale(0.5, 0.5, 0.5, new Sphere))));
    psUnion = new Union(new Translate(-3.0, 0.0, 0.0, ps), psUnion);
    
    ps = new Rotate(0.0, fRad, 0.0, 
            new Material(0.33, REFRACT_WATER, g_sWater, MAT_TRANSPARENT, ColoredGlass,
               new Surface(SURF_POLISHED, 0, new PRISM)));
    psUnion = new Union (psUnion, new Translate(3.0, 0.0, 0.0, ps));
    ps = new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, REFRACT_AIR, g_sEqualWhite, MAT_TRANSPARENT, 0,
               new Surface(SURF_POLISHED, 0, new Scale(0.5, 0.5, 0.5, new Sphere))));
    psUnion = new Union(psUnion, new Translate(3.0, 0.0, 0.0, ps));
   
    psUnion = new Rotate(fPitch, 0.0, 0.0, new Union(psUnion, new Translate(0.0, 3.0, -10.0, new Scale(12.0, 1.0, 12.0, new Cube))));
    return psUnion;
}

static Solid*
MakeModelE(double fAng)             // water prisms  with glass, copper, air spheres
{
    Spectrum sBlue;
    double fRad, fPitch;
    Solid *ps, *ps2;

    fRad = fAng*D_PI/180.0;
    fPitch = -20.0 * D_PI/180.0;
    sBlue = g_sBlue * 0.5 + g_sEqualWhite * 0.5;
    ps = new Material(0.5, REFRACT_PYREX, g_sBlack, MAT_TRANSPARENT, 0, new Surface(SURF_POLISHED, 0, new Cube));
    ps = new Translate(0.0, 1.5, 0.0, new Scale(1.9, 1.0, 1.9, ps));
    ps2 = new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC,
               new Surface(SURF_SHINY, 0, new Sphere)));
    ps = new Union(ps, new Translate(3.0, 0.0, -3.0, ps2));
    ps = new Union(ps, new Translate(0.0, 1.9, -10.0, new Scale(12.0, 1.0, 12.0, new Cube)));
    ps = new Rotate(fPitch, 0.0, 0.0, ps);
    return ps;
}

static Solid*
ChainLink(double fAng)
{
    Solid *psT1, *psT2, *ps;
    double r, g, fRad;

    fRad = fAng * D_PI/180.0;
    r = 0.5;
    g = 0.1;
    psT1 = new Difference(new Torus(r), new Translate(+0.5*(1.0+r), 0.0, 0.0, new Scale(0.5*(1.0+r), (1.01+r), r+0.01, new Cube)));
    psT2 = new Difference(new Torus(r), new Translate(-0.5*(1.0+r), 0.0, 0.0, new Scale(0.5*(1.0+r), (1.01+r), r+0.01, new Cube)));
    ps = new Union(new Translate(-r-g, 0.0, 0.0, psT1), new Translate(+r+g, 0.0, 0.0, psT2));
    ps = new Union(new Translate(0.0, +1.0, 0.0, new Scale(r+0.000001+g, r, r, new Rotate(0.0, 0.5*D_PI, 0.0, new Cylinder))), ps);
    ps = new Union(new Translate(0.0, -1.0, 0.0, new Scale(r+0.000001+g, r, r, new Rotate(0.0, 0.5*D_PI, 0.0, new Cylinder))), ps);
    return new Rotate(fRad, 0.0, 0.0, new BoundingBox(ps));
}

static Solid*
MakeModelF(double fAng)
{
    Solid *ps;
    double fRad, fPitch, r, g, d;

    fRad = fAng*D_PI/180.0;
    fPitch = -20.0 * D_PI/180.0;
    r = 0.5;
    g = 0.1;
    d = 1.0 + r + r + g + g;
    // ps = new Rotate(0.0, fRad, 0.0, Cuboid(0.1));
    ps = ChainLink(0.0);
    ps = new Union(ps, new Translate(d, 0.0, 0.0, ChainLink(90.0)));
    ps = new Union(ps, new Translate(2.0*d, 0.0, 0.0, ChainLink(0.0)));
    ps = new Union(ps, new Translate(-d, 0.0, 0.0, ChainLink(90.0)));
    ps = new Union(ps, new Translate(0.0, 3, -10.0, new Scale(12.0, 1.0, 12.0, new Cube)));
    ps = new Rotate(fPitch, 0.0, 0.0, ps);
    return ps;
}

static Solid*
MakeModelG(double fAng)             // various spheres moving through glass & water cuboids
{
    double fPitch, fYaw, r;
    Solid *ps, *psUnion;
    Spectrum sEmerald, sBlue, sCobalt;

    r = 0.075;
    sEmerald = g_sEmerald * 0.01;
    sCobalt = g_sCobaltGlass * 0.001;
    sBlue = (g_sBlue + g_sEqualWhite) * 0.5;

    psUnion = new Translate(0.0, 3.0, -10.0, new Scale(12.0, 1.0, 12.0, new Cube));

    ps = new Material(0.5, REFRACT_PYREX, g_sBlack, MAT_TRANSPARENT, 0, new Surface(SURF_POLISHED, 0, Cuboid(r)));
    psUnion = new Union(new Translate(2.0, 0.0, 0.0, ps), psUnion);
    ps = new Material(0.5, REFRACT_WATER, g_sWater, MAT_TRANSPARENT, ColoredGlass, new Surface(SURF_POLISHED, 0, Cuboid(r)));
    psUnion = new Union(new Translate(-2.0, 0.0, 0.0, ps), psUnion);

    ps =  new Material(0.9, REFRACT_PYREX, g_sCopper, MAT_GOLD, 0, new Surface(SURF_POLISHED, 0, new Scale(0.5, 0.5, 0.5, new Sphere)));
    psUnion = new Union(new Translate(fmod(fAng/10 + 0.0, 8.0) - 4.0,  0.0, 0.0, ps), psUnion);
    ps =  new Material(0.9, REFRACT_PYREX, sEmerald, MAT_TRANSPARENT, ColoredGlass, new Surface(SURF_POLISHED, Bumpy, new Scale(0.5, 0.5, 0.5, new Sphere)));
    psUnion = new Union(new Translate(fmod(fAng/10 + 2.0, 8.0) - 4.0,  0.0, 0.0, ps), psUnion);
    ps =  new Material(0.9, REFRACT_WATER, sBlue, MAT_PLASTIC, ColoredGlass, new Surface(SURF_POLISHED, 0, new Scale(0.5, 0.5, 0.5, new Sphere)));
    psUnion = new Union(new Translate(fmod(fAng/10 + 4.0, 8.0) - 4.0,  0.0, 0.0, ps), psUnion);
    ps =  new Material(0.9, REFRACT_AIR, g_sBlack, MAT_TRANSPARENT, 0, new Surface(SURF_POLISHED, 0, new Scale(0.5, 0.5, 0.5, new Sphere)));
    psUnion = new Union(new Translate(fmod(fAng/10 + 6.0, 8.0) - 4.0,  0.0, 0.0, ps), psUnion);

    fPitch = -40.0 * D_PI/180.0;
    fYaw = 20.0 * D_PI/180.0;
    psUnion = new Rotate(fPitch, fYaw, 0.0, psUnion);
    return psUnion;
}

static Solid*
MakeModelG2(double fAng)             // various spheres moving through glass & water cuboids
{
    double fPitch, fYaw, r;
    Solid *ps, *psUnion;
    Spectrum sEmerald, sBlue;

    r = 0.075;
    sEmerald = g_sEmerald * 0.01;
    sBlue = (g_sBlue + g_sEqualWhite) * 0.5;

    psUnion = new Translate(0.0, 3.0, -10.0, new Scale(12.0, 1.0, 12.0, new Cube));

    ps = new Material(0.5, REFRACT_WATER, g_sWater, MAT_TRANSPARENT, ColoredGlass, new Surface(SURF_POLISHED, 0, Cuboid(r)));
    psUnion = new Union(new Translate(-2.0, 0.0, 0.0, ps), psUnion);

    fPitch = -40.0 * D_PI/180.0;
    fYaw = 20.0 * D_PI/180.0;
    psUnion = new Rotate(fPitch, fYaw, 0.0, psUnion);
    return psUnion;
}

static Solid*
MakeModelH(double fAng)             // chain links glass & water cuboids
{
    double fPitch, fYaw, r, rc, g, d;
    Solid *ps, *psUnion, *psCubes, *psChain;
    Spectrum sBlue;

    r = 0.075;
    rc = 0.5;
    g = 0.1;
    d = 1.0 + rc + rc + g + g;                  // 2.2
    sBlue = (g_sEqualWhite + g_sBlue)*0.5;

    psUnion = new Translate(0.0, 3.0, -10.0, new Scale(12.0, 1.0, 12.0, new Cube));

    ps = new Material(0.5, REFRACT_PYREX, g_sBlack, MAT_TRANSPARENT, 0, new Surface(SURF_POLISHED, 0, Cuboid(r)));
    psUnion = new Union(new Translate(2.0, 0.0, 0.0, ps), psUnion);
    ps = new Material(0.5, REFRACT_WATER, g_sWater, MAT_TRANSPARENT, ColoredGlass, new Surface(SURF_POLISHED, 0, Cuboid(r)));
    psCubes = new Union(new Translate(-2.0, 0.0, 0.0, ps), psUnion);

    ps =  new Material(0.9, g_sCopper, MAT_METAL, new Surface(SURF_POLISHED, 0, new Scale(0.5, 0.5, 0.5, ChainLink(0.0))));
    psUnion = (new Translate(fmod(fAng + 0.0*d, 5.0*d) - 2.5*d,  0.0, 0.0, ps) );
    ps =  new Material(0.9, REFRACT_PYREX, g_sBlack, MAT_TRANSPARENT, 0, new Surface(SURF_POLISHED, 0, 
                                                                                new Scale(0.5, 0.5, 0.5, ChainLink(90.0))));
    psUnion = new Union(new Translate(fmod(fAng + 0.5*d, 5.0*d) - 2.5*d,  0.0, 0.0, ps), psUnion);
    ps =  new Material(0.9, REFRACT_WATER, g_sWater, MAT_TRANSPARENT, ColoredGlass, new Surface(SURF_POLISHED, 0, 
                                                                                new Scale(0.5, 0.5, 0.5, ChainLink(0.0))));
    psUnion = new Union(new Translate(fmod(fAng + 1.0*d, 5.0*d) - 2.5*d,  0.0, 0.0, ps), psUnion);
    ps =  new Material(0.9, REFRACT_AIR, g_sBlack, MAT_TRANSPARENT, 0, new Surface(SURF_POLISHED, 0, 
                                                                                new Scale(0.5, 0.5, 0.5, ChainLink(0.0))));
    psUnion = new Union(new Translate(fmod(fAng + 2.0*d, 5.0*d) - 2.5*d,  0.0, 0.0, ps), psUnion);
    ps =  new Material(0.9, REFRACT_PYREX, sBlue, MAT_PLASTIC, 0, new Surface(SURF_SHINY, 0, 
                                                                                new Scale(0.5, 0.5, 0.5, ChainLink(90.0))));
    psUnion = new Union(new Translate(fmod(fAng + 1.5*d, 5.0*d) - 2.5*d,  0.0, 0.0, ps), psUnion);

        ps =  new Material(0.9, g_sCopper, MAT_METAL, new Surface(SURF_POLISHED, 0, new Scale(0.5, 0.5, 0.5, ChainLink(90.0))));
    psUnion = new Union(new Translate(fmod(fAng + 2.5*d, 5.0*d) - 2.5*d,  0.0, 0.0, ps), psUnion);
    ps =  new Material(0.9, REFRACT_PYREX, g_sBlack, MAT_TRANSPARENT, 0, new Surface(SURF_POLISHED, 0, 
                                                                                new Scale(0.5, 0.5, 0.5, ChainLink(0.0))));
    psUnion = new Union(new Translate(fmod(fAng + 3.0*d, 5.0*d) - 2.5*d,  0.0, 0.0, ps), psUnion);
    ps =  new Material(0.9, REFRACT_WATER, g_sWater, MAT_TRANSPARENT, ColoredGlass, new Surface(SURF_POLISHED, 0, 
                                                                                new Scale(0.5, 0.5, 0.5, ChainLink(90.0))));
    psUnion = new Union(new Translate(fmod(fAng + 3.5*d, 5.0*d) - 2.5*d,  0.0, 0.0, ps), psUnion);
    ps =  new Material(0.9, REFRACT_AIR, g_sBlack, MAT_TRANSPARENT, 0, new Surface(SURF_POLISHED, 0, 
                                                                                new Scale(0.5, 0.5, 0.5, ChainLink(90.0))));
    psUnion = new Union(new Translate(fmod(fAng + 4.5*d, 5.0*d) - 2.5*d,  0.0, 0.0, ps), psUnion);
    ps =  new Material(0.9, REFRACT_PYREX, sBlue, MAT_PLASTIC, 0, new Surface(SURF_SHINY, 0, 
                                                                                new Scale(0.5, 0.5, 0.5, ChainLink(0.0))));
    psChain = new Union(new Translate(fmod(fAng + 4.0*d, 5.0*d) - 2.5*d,  0.0, 0.0, ps), psUnion);
    psChain = new Scale(1.0/1.5, 1.0/1.5, 1.0/1.5, psChain);

    fPitch = -40.0 * D_PI/180.0;
    fYaw = 20.0 * D_PI/180.0;
    psUnion = new Rotate(fPitch, fYaw, 0.0, new Union(psCubes, psChain));
    return psUnion;
}

static Solid*
MakeModelI(double fAng)         // Just a sphere
{
    Spectrum sBlue;
    Solid *ps, *psUnion;

    sBlue = (g_sEqualWhite + g_sBlue)*0.5;
    ps = new Translate(-3.0, 0.0, 0.0, new Material(0.5, REFRACT_PYREX, g_sWater, MAT_TRANSPARENT, 0, new Surface(SURF_POLISHED, 0, new Sphere)));
    psUnion = ps;
    ps = new Translate(+0.0, 0.0, 0.0, new Material(0.5, REFRACT_PYREX, sBlue, MAT_PLASTIC, 0, new Surface(SURF_POLISHED, Bumpy, new Sphere)));
    psUnion = new Union(ps, psUnion);
    ps = new Translate(+3.0, 0.0, 0.0, new Material(0.5, REFRACT_PYREX, g_sCopper, MAT_METAL, 0, new Surface(SURF_POLISHED, 0, new Sphere)));
    psUnion = new Union(ps, psUnion);

    psUnion = new Union(psUnion, new Translate(0.0, 3.0, -10.0, new Scale(12.0, 1.0, 12.0, new Cube)));
    return psUnion;
}

static Solid*
MakeModelJ(double fAng)                             // rotating metal spheres, camera 9
{
    double fRad, fRadX, fPitch;
    Spectrum sBlue;
    Solid *ps, *psUnion;

    fRad  = 0.0 * D_PI/180.0;
    fRadX = -fAng * D_PI/180.0;
    fPitch = -30.0 * D_PI/180.0;
    sBlue = g_sBlue * 0.5 + g_sWhite * 0.5;         // not used by metal shaders

    ps = new Rotate(fRadX, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_COPPER,
               new Surface(SURF_POLISHED, VoronoiCrinkled, new Sphere)));
    psUnion = new Translate(-2.5, 1.0, 0.0, ps);
    ps = new Rotate(fRadX, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_GOLD,
               new Surface(SURF_POLISHED, Dented, new Sphere)));
    psUnion = new Union(new Translate(0.0, 1.0, 0.0, ps), psUnion);
    ps = new Rotate(fRadX, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_IRON,
               new Surface(SURF_POLISHED, Bumpy, new Sphere)));
    psUnion = new Union(new Translate(2.5, 1.0, 0.0, ps), psUnion);
    ps = new Rotate(fRadX, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_ALUMINUM,
               new Surface(SURF_POLISHED, MicroFacets2, new Sphere)));
    psUnion = new Union(new Translate(-1.25, -0.5, -1.5, ps), psUnion);
    ps = new Rotate(fRadX, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_OSMIUM,
               new Surface(SURF_POLISHED, Tectonic, new Sphere)));
    psUnion = new Union(new Translate(+1.25, -0.5, -1.5, ps), psUnion);

    psUnion = new Union(psUnion, new Translate(0.0, 3.0, 0.0, new Scale(32.0, 1.0, 32.0, new Cube)));
    psUnion = new Rotate(fPitch, 0.0, 0.0, psUnion);
    return psUnion;
}

static Solid*
MakeModelK(double x)             // various spheres moving through glass & water cuboids
{
    double fPitch = 0.0, fYaw = 0.0, r, fAng;
    Solid *ps, *psUnion;
    Spectrum sEmerald, sBlue, sCobalt;

    r = 0.075;
    sEmerald = g_sEmerald * 0.01;
    sCobalt = g_sCobaltGlass * 0.001;
    sBlue = (g_sBlue + g_sEqualWhite) * 0.5;

    fAng = 20.0;

    psUnion = new Translate(0.0, 3.0, -10.0, new Scale(12.0, 1.0, 12.0, new Cube));

    ps = new Translate(x, 0.0, -0.05, new Scale(0.015, 0.015, 0.015, new Sphere));
    psUnion = new Union(ps, psUnion);

    ps = new Material(0.5, REFRACT_PYREX, g_sBlack, MAT_TRANSPARENT, 0, new Surface(SURF_POLISHED, 0, Cuboid(r)));
    psUnion = new Union(new Translate(2.0, 0.0, 0.0, ps), psUnion);
    ps = new Material(0.5, REFRACT_WATER, g_sWater, MAT_TRANSPARENT, ColoredGlass, new Surface(SURF_POLISHED, 0, Cuboid(r)));
    psUnion = new Union(new Translate(-2.0, 0.0, 0.0, ps), psUnion);

    ps =  new Material(0.9, REFRACT_PYREX, g_sCopper, MAT_GOLD, 0, new Surface(SURF_POLISHED, 0, new Scale(0.5, 0.5, 0.5, new Sphere)));
    psUnion = new Union(new Translate(fmod(fAng/10 + 0.0, 8.0) - 4.0,  0.0, 0.0, ps), psUnion);
    ps =  new Material(0.9, REFRACT_PYREX, sEmerald, MAT_TRANSPARENT, ColoredGlass, new Surface(SURF_POLISHED, Bumpy, new Scale(0.5, 0.5, 0.5, new Sphere)));
    psUnion = new Union(new Translate(fmod(fAng/10 + 2.0, 8.0) - 4.0,  0.0, 0.0, ps), psUnion);
    ps =  new Material(0.9, REFRACT_WATER, g_sWater, MAT_TRANSPARENT, ColoredGlass, new Surface(SURF_POLISHED, 0, new Scale(0.5, 0.5, 0.5, new Sphere)));
    psUnion = new Union(new Translate(fmod(fAng/10 + 4.0, 8.0) - 4.0,  0.0, 0.0, ps), psUnion);
    ps =  new Material(0.9, REFRACT_AIR, g_sBlack, MAT_TRANSPARENT, 0, new Surface(SURF_POLISHED, 0, new Scale(0.5, 0.5, 0.5, new Sphere)));
    psUnion = new Union(new Translate(fmod(fAng/10 + 6.0, 8.0) - 4.0,  0.0, 0.0, ps), psUnion);

    fPitch = -40.0 * D_PI/180.0;
    //  fYaw = 20.0 * D_PI/180.0;
    psUnion = new Rotate(fPitch, fYaw, 0.0, psUnion);
    return psUnion;
}

static Solid*
MakeModelL(double fAng)                             // solid textures
{
    double fRad, fRadX, fPitch;
    Spectrum sBlue;
    Solid *ps, *psUnion;

    fRad  = 0.0 * D_PI/180.0;
    fRadX = -fAng * D_PI/180.0;
    fPitch = -30.0 * D_PI/180.0;
    sBlue = g_sBlue * 0.5 + g_sWhite * 0.5;         // not used by metal shaders

    ps = new Rotate(fRadX, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC, Agate,
               new Surface(SURF_BURNISHED, 0, new Sphere)));
    psUnion = new Translate(-2.5, 1.0, 0.0, ps);
    ps = new Rotate(fRadX, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC, Marble,
               new Surface(SURF_BURNISHED, 0, new Sphere)));
    psUnion = new Union(new Translate(0.0, 1.0, 0.0, ps), psUnion);
    ps = new Rotate(fRadX, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC, Clouds,
               new Surface(SURF_BURNISHED, 0, new Sphere)));
    psUnion = new Union(new Translate(2.5, 1.0, 0.0, ps), psUnion);
    ps = new Rotate(fRadX, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC, Voronoi,
               new Surface(SURF_BURNISHED, 0, new Sphere)));
    psUnion = new Union(new Translate(+1.25, -0.5, -1.5, ps), psUnion);
    ps = new Rotate(fRadX, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC, Wood,
               new Surface(SURF_BURNISHED, 0, new Cylinder)));
    psUnion = new Union(new Translate(-1.25, -0.5, -1.5, ps), psUnion);

    psUnion = new Union(psUnion, new Translate(0.0, 3.0, 0.0, new Scale(32.0, 1.0, 32.0, new Cube)));
    psUnion = new Rotate(fPitch, 0.0, 0.0, psUnion);
    return psUnion;
}

static Solid*
MakeModelL2(double fAng)                             // wood sphere
{
    double fRad, fRadX, fPitch, yDelta;
    Spectrum sBlue;
    Solid *ps, *psUnion;

    fRad  = 25.0 * D_PI/180.0;
    fRadX = -fAng * D_PI/180.0;
    fPitch = -30.0 * D_PI/180.0;
    sBlue = g_sBlue * 0.5 + g_sWhite * 0.5;         // not used by metal shaders
    yDelta = 0.0;

    ps = new Rotate(fRadX, fRad, 0.0, new Translate(0.0, -yDelta, 0.0, 
            new Surface(SURF_BURNISHED, 0, 
            new Material(0.9, sBlue, MAT_PLASTIC, Wood,
               new Translate(0.0, yDelta, 0.0, new Cylinder)))));
    psUnion = new Translate(0.0, 0.0, 0.0, ps);

    psUnion = new Union(psUnion, new Translate(0.0, 3.0, 0.0, new Scale(32.0, 1.0, 32.0, new Cube)));
    psUnion = new Rotate(fPitch, 0.0, 0.0, psUnion);
    return psUnion;
}

static Solid*
MakeModelL3(double fAng)                             // new bump maps
{
    double fRad, fRadX, fPitch;
    Spectrum sBlue, sGreen, sPink, sOrange, sLavender;
    Solid *ps, *psUnion;

    fRad  = 0.0 * D_PI/180.0;
    fRadX = -fAng * D_PI/180.0;
    fPitch = -30.0 * D_PI/180.0;
    sBlue = g_sBlue * 0.5 + g_sWhite * 0.5;         // not used by metal shaders
    sGreen = g_sGreen * 0.5 + g_sWhite * 0.5;
    sPink = g_sRed * 0.5 + g_sWhite * 0.5;
    sOrange = g_sOrange * 0.5 + g_sWhite * 0.5;
    sLavender = g_sRed * 0.25 + g_sBlue * 0.25 + g_sWhite * 0.5;

    ps = new Rotate(fRadX, fRad, 0.0, 
            new Material(0.9, sPink, MAT_PLASTIC,
               new Surface(SURF_SMOOTH, Bumpy, new Sphere)));
    psUnion = new Translate(+2.5, 1.0, 0.0, ps);
    ps = new Rotate(fRadX, fRad, 0.0, 
            new Material(0.9, sGreen, MAT_PLASTIC,
               new Surface(SURF_SMOOTH, VoronoiCrinkled, new Sphere)));
    psUnion = new Union(new Translate(-2.5, 1.0, 0.0, ps), psUnion);
    ps = new Rotate(fRadX, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC,
               new Surface(SURF_POLISHED, Blistered, new Sphere)));
    psUnion = new Union(new Translate(0.0, 1.0, 0.0, ps), psUnion);
    ps = new Rotate(fRadX, fRad, 0.0, 
            new Material(0.9, sLavender, MAT_PLASTIC,
               new Surface(SURF_SHINY, Dented, new Sphere)));
    psUnion = new Union(new Translate(+1.25, -0.5, -1.5, ps), psUnion);
    ps = new Rotate(fRadX, fRad, 0.0, 
            new Material(0.9, sOrange, MAT_PLASTIC,
               new Surface(SURF_SHINY, Tectonic, new Sphere)));
    psUnion = new Union(new Translate(-1.25, -0.5, -1.5, ps), psUnion);

    psUnion = new Union(psUnion, new Translate(0.0, 3.0, 0.0, new Scale(32.0, 1.0, 32.0, new Cube)));
    psUnion = new Rotate(fPitch, 0.0, 0.0, psUnion);
    return psUnion;
}

static Solid*
MakeModelL4(double fAng)                             // solid textures, planet.  Camera z = 5.0.  (12, 105 is good)
{
    double fRad, fRadX;
    Spectrum sBlue;
    Solid *ps, *psUnion;

    fRad  = fAng * D_PI/180.0;
    fRadX = 90.0 * D_PI/180.0;
    sBlue = g_sBlue * 0.5 + g_sWhite * 0.5;         // not used by metal shaders

    ps = new Rotate(fRadX, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC, EarthTexture,
               new Surface(SURF_SHINY, 0, new Rotate(0.0, 0.0, 0.0, new Sphere))));
    psUnion = new Translate(0.0, 0.0, 0.0, ps);
    return psUnion;
}

static Solid*
MakeModelL5(double fAng)                             // solid textures, planet.  Camera z = 5.0
{
    double fRad, fRadX;
    Spectrum sBlue;
    Solid *ps, *psUnion;

    fRad  = fAng * D_PI/180.0;
    fRadX = -90.0 * D_PI/180.0;
    sBlue = g_sBlue * 0.5 + g_sWhite * 0.5;         // not used by metal shaders

    ps = new Rotate(fRadX, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC, EarthTexture,
               new Surface(SURF_SHINY, 0, new Rotate(-fRadX, 0.0, 0.0, new Torus(0.3)))));
    psUnion = new Translate(0.0, 0.0, 0.0, ps);
    return psUnion;
}

static Solid*
MakeModelM(double fAng)                         // test image with retort and bump and material maps
{
    Solid *ps, *ps2, *psDiff, *psCone, *psBall;
    Solid *psStand;
    double f20, f90, fRad, f;
    Spectrum sBlue;

    f20 = 20.0 * D_PI/180.0;
    f90 = 90.0 * D_PI/180.0;
    fRad = (fAng - 120.0) * D_PI/180.0;

    sBlue = g_sBlue * 0.5 + g_sWhite * 0.5;

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

    psDiff = new Rotate(0.0, fRad, 0.0, new Difference(ps, ps2));
    psDiff = new Material(0.33, REFRACT_WATER, g_sWater, MAT_TRANSPARENT, ColoredGlass,
               new Surface(SURF_POLISHED, 0, psDiff));
    psDiff = new BoundingBox(psDiff);
    f = 1.05/sqrt(2.0);
    psStand = new Torus(0.05/f);
    psStand = new Union(psStand, new Translate(1.0, 0.0, -0.5, new Scale(0.05, 0.05, 0.5, new Cylinder)));
    psStand = new Union(psStand, new Translate(-0.5, +0.866, -0.5, new Scale(0.05, 0.05, 0.5, new Cylinder)));
    psStand = new Union(psStand, new Translate(-0.5, -0.866, -0.5, new Scale(0.05, 0.05, 0.5, new Cylinder)));
    psStand = new Translate(0.0, f, 0.0, new Rotate(0.5*D_PI, 0.0, 0.0, new Scale(f, f, f, psStand)));
    psStand = new Material(0.9, g_sCopper, MAT_METAL,
               new Surface(SURF_POLISHED, 0, psStand));
    psStand = new BoundingBox(psStand);
    ps = new Union(psStand, psDiff);

 

    psBall = new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC, Marble,
               new Surface(SURF_SMOOTH, 0, new Scale(0.7, 1.2, 0.7, new Dodecahedron))));
    psBall = new BoundingBox(psBall);
    ps = new Union(new Translate(3.0, 0.0, 0.0, psBall), ps);
    //ps = new Translate(3.0, 0.0, 0.0, psBall);

    psBall = new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC,
               new Surface(SURF_SMOOTH, Bumpy, new Scale(0.7, 1.2, 0.7, new Sphere))));
    psBall = new BoundingBox(psBall);
    ps = new Union(new Translate(-3.0, 0.0, 0.0, psBall), ps);

    ps = new Union(ps, new Translate(0.0, 2.5, -10.0, new Scale(12.0, 1.0, 12.0, new Cube)));
    ps = new Rotate(-f20*0.7, 0.0, 0.0, ps);
    return ps;
}

static Solid*
MakeModelN(double fAng)                                 // 4 spheres, for tuning shader KS values
{
    Spectrum sBlue;
    Solid *ps, *psBall;
    double fRad;

    sBlue = g_sBlue * 0.5 + g_sWhite * 0.5;
    fRad = (fAng - 120.0) * D_PI/180.0;

    psBall = new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC, 0,
               new Surface(SURF_BURNISHED, 0, new Scale(1.0, 1.0, 1.0, new Sphere))));
    ps = new Translate(3.3, 0.0, 0.0, psBall);
    psBall = new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC, 0,
               new Surface(SURF_SMOOTH, 0, new Scale(1.0, 1.0, 1.0, new Sphere))));
    ps = new Union(ps, new Translate(1.1, 0.0, 0.0, psBall));
    psBall = new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC, 0,
               new Surface(SURF_SHINY, 0, new Scale(1.0, 1.0, 1.0, new Sphere))));
    ps = new Union(ps, new Translate(-1.1, 0.0, 0.0, psBall));
    psBall = new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC, 0,
               new Surface(SURF_POLISHED, 0, new Scale(1.0, 1.0, 1.0, new Sphere))));
    ps = new Union(ps, new Translate(-3.3, 0.0, 0.0, psBall));

    return ps;
}

static Solid*
MakeModelO(double fAng)                                 // Test textures and coordinate systems
{
    Spectrum sBlue;
    Solid *ps, *psBall;
    double fRad;

    sBlue = g_sBlue * 0.5 + g_sWhite * 0.5;
    fRad = (fAng) * D_PI/180.0;

    psBall = new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC, 0,
               new Surface(SURF_SMOOTH, Bumpy, new Scale(1.5, 0.5, 0.5, new Sphere))));
    ps =               new Translate(-3.0, 1.0, 0.0, psBall);
    psBall = new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC, 0,
               new Scale(1.5, 0.5, 0.5, new Surface(SURF_SMOOTH, Bumpy, new Sphere))));
    ps = new Union(ps, new Translate(0.0, 1.0, 0.0, psBall));
    psBall = new Surface(SURF_SMOOTH, Bumpy, new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC, 0,
                new Scale(1.5, 0.5, 0.5, new Sphere))));
    ps = new Union(ps, new Translate(3.0, 1.0, 0.0, psBall));

    psBall = new Rotate(0.0, fRad, 0.0, 
            new Material(0.9, sBlue, MAT_PLASTIC, Marble,
               new Surface(SURF_SMOOTH, 0, new Scale(1.5, 0.5, 0.5, new Sphere))));
    ps = new Union(ps, new Translate(-3.0, -1.0, 0.0, psBall));
    psBall = new Rotate(0.0, fRad, 0.0, 
               new Scale(1.5, 0.5, 0.5, new Surface(SURF_SMOOTH, 0, 
               new Material(0.9, sBlue, MAT_PLASTIC, Marble, new Sphere))));
    ps = new Union(ps, new Translate(0.0, -1.0, 0.0, psBall));
    psBall = new Material(0.9, sBlue, MAT_PLASTIC, Marble, new Surface(SURF_SMOOTH, 0, new Rotate(0.0, fRad, 0.0, 
               new Scale(1.5, 0.5, 0.5, new Sphere))));
    ps = new Union(ps, new Translate(3.0, -1.0, 0.0, psBall));

    return ps;
}

static Solid*
MakeModelP(double fAng)                                 // DodecahedronCage
{
    Spectrum sBlue;
    Solid *ps, *ps1, *ps2, *psBall;
    double fRad, f20;

    sBlue = g_sBlue * 0.5 + g_sWhite * 0.5;
    fRad = (fAng) * D_PI/180.0;
    f20 = 10.0 * D_PI/180.0;

    //psBall = new Sphere;
    // psBall = new Surface(SURF_POLISHED, new Material(0.9, sBlue, MAT_PLASTIC, DodecahedronCage(0.05, 0.1)));
    psBall = new Surface(SURF_POLISHED, new Material(0.9, sBlue, MAT_PLASTIC, PseudoHelix(0.1, 0.1)));
    psBall = new Union(psBall, new Translate(0.0, 0.0, 0.1, 
                        new Surface(SURF_POLISHED, new Material(0.9, sBlue, MAT_PLASTIC, PseudoHelix(0.1, 0.1)))));
    psBall = new Union(psBall, new Translate(0.0, 0.0, 0.2, 
                        new Surface(SURF_POLISHED, new Material(0.9, sBlue, MAT_PLASTIC, PseudoHelix(0.1, 0.1)))));
    psBall = new Union(psBall, new Translate(0.0, 0.0, 0.3, 
                        new Surface(SURF_POLISHED, new Material(0.9, sBlue, MAT_PLASTIC, PseudoHelix(0.1, 0.1)))));
    psBall = new Union(psBall, new Translate(0.0, 0.0, 0.4, 
                        new Surface(SURF_POLISHED, new Material(0.9, sBlue, MAT_PLASTIC, PseudoHelix(0.1, 0.1)))));
    psBall = new Rotate(0.0, +fRad, 0.0, psBall);
    ps1 =    new Translate(+1.25, 0.0, 0.0, psBall);

    psBall = new Surface(SURF_POLISHED, new Material(0.9, sBlue, MAT_PLASTIC, PseudoHelix(0.1, 0.3)));
    psBall = new Union(psBall, new Translate(0.0, 0.0, 0.3, 
                        new Surface(SURF_POLISHED, new Material(0.9, sBlue, MAT_PLASTIC, PseudoHelix(0.1, 0.3)))));
    psBall = new Union(psBall, new Translate(0.0, 0.0, 0.6, 
                        new Surface(SURF_POLISHED, new Material(0.9, sBlue, MAT_PLASTIC, PseudoHelix(0.1, 0.3)))));
    psBall = new Union(psBall, new Translate(0.0, 0.0, 0.9, 
                        new Surface(SURF_POLISHED, new Material(0.9, sBlue, MAT_PLASTIC, PseudoHelix(0.1, 0.3)))));
    psBall = new Union(psBall, new Translate(0.0, 0.0, 1.2, 
                        new Surface(SURF_POLISHED, new Material(0.9, sBlue, MAT_PLASTIC, PseudoHelix(0.1, 0.3)))));
    psBall = new Rotate(0.0, +fRad, 0.0, psBall);
    psBall = new Rotate(0.0, +fRad, 0.0, psBall);
    ps2 =    new Translate(-1.25, 0.0, 0.0, psBall);

    ps = new Rotate(0.5*D_PI, 0.0, 0.0, new Union(ps1, ps2));
    ps = new Union(ps, new Translate(0.0, 2.5, -10.0, new Scale(12.0, 1.0, 12.0, new Material(0.9, g_sEqualWhite, MAT_PLASTIC, new Cube))));
    ps = new Rotate(-f20*0.7, 0.0, 0.0, ps);
    ps = new Union(ps, new Scale (100.0, 100.0, 100.0, new Material(0.9, sBlue, MAT_PLASTIC, new HollowSphere)));
    return ps;

}

static Solid*
MakeModelQ(double fAng)             // test array of objects
{
    Spectrum sBlue;
    Solid *ps;
    double fRad, f90, r, x, y;

    sBlue = (g_sBlue + g_sWhite) * 0.5;
    fRad = (fAng) * D_PI/180.0;
    f90 = 0.5*D_PI;

    x = y = 0.0;
    r = 1.0;

    r = (fAng+10)/360.0;
    x = 1.05*fAng/360.0;
    y = 0.57*fAng/360.0;
 

    ps = MakeUnion(
        /*
        new Translate(-3.0, +1.5, 0.0, new Rotate(0.0, fRad, 0.0, new Scale(0.5, 0.5, 0.5, new Tetrahedron))),
        new Translate(+0.0, +1.5, 0.0, new Rotate(0.0, fRad, 0.0, new Cube)),
        new Translate(+3.0, +1.5, 0.0, new Rotate(0.0, fRad, 0.0, new Octahedron)),
        new Translate(-3.0, -1.5, 0.0, new Rotate(0.0, fRad, 0.0, new Dodecahedron)),
        new Translate(+0.0, -1.5, 0.0, new Rotate(0.0, fRad, 0.0, new Icosahedron)),
        new Translate(+3.0, -1.5, 0.0, new Rotate(f90, fRad, 0.0, new Prism(7)))
        */

        new Translate(-3.0, +1.5, 0.0, new Rotate(0.0, fRad, 0.0, Cuboid(0.1, 0.5, 0.7) )),
        new Translate(+0.0, +1.5, 0.0, new Rotate(0.0, fRad, 0.0, IcosahedronCage(0.1, 0.15) )),
        new Translate(+3.0, +1.5, 0.0, new Rotate(0.0, fRad, 0.0, PseudoHelix(0.1, 0.3, 7.0*D_PI) )),
        new Translate(-3.0, -1.5, 0.0, new Rotate(0.0, fRad, 0.0, new Union(Diskoid(0.1), Rodoid(0.1)) )),
        new Translate(+0.0, -1.1, 0.0, new Rotate(0.1, fRad, 0.1, Slaboidoid(0.5, 1.1, 0.2, 0.1, 0.025) )),
        new Translate(+3.0, -1.5, 0.0, new Rotate(0.0, fRad, 0.0, new Scale(0.5, 0.5, 0.5, Knoboid(0.4)) ))
    );
    ps = // new Translate(x, y, 0.0, new Rotate(0.0, fRad, 0.0, new Scale(r, r, r, new Torus(0.1))));
       new Surface(SURF_POLISHED, new Material(0.9, sBlue, MAT_PLASTIC, ps));
      //  new Rotate(0.1, fRad, 0.1, new Scale(r, r, r, new Translate(0.0, 1.0, 0.0, new Torus(0.1))));
    return ps;
}

static Solid*
MakeModelR(double fAng)
{
    Solid *ps;
    double fAng1, fAng2;

    ps = PseudoHelix(0.2, 0.5, fAng * D_PI/180.0);
    return new Rotate(0.25*D_PI, 0.0, 0.0, ps);

    fAng1 = (0.0) * D_PI/180.0;
    fAng2 = (0.0+fAng) * D_PI/180.0;
    ps = new Rotate(0.0, 20.0*D_PI/180.0, 0.0, new AngularClip(fAng1, fAng2, new Torus(0.2)));
    return new Translate(0.0, +0.0, 0.0, ps);
}

static Solid*
MakeModelS(double r)            //  extended Knoboid test
{
    Solid *ps;

   // ps = new Intersection(new Scale(3.0, 3.0, 1.0, new Cube), new ConeHyperCylinder(r));
    ps = Knoboid(r, 1.0);
    ps = Plastic(ps);
    ps = new Rotate(0.0, 0.5*D_PI, 0.0, ps);
    return ps;
}

static Solid*
MakeModelT(double fAng)         //  test field of view, camera at z = 5.0,
{
    Solid *ps;
    double fRad;
    int i;

    ps = new Rotate(0.0, 0.0, 0.0, new Translate(0.0, 0.0, 5.0, new Sphere));
    for (i = 10; i < 360; i += 10) {
        fRad = double(i) * D_PI/180.0;
        ps = new Union(ps, new Rotate(0.0, fRad, 0.0, new Translate(0.0, 0.0, 25.0, new Sphere)));
    }
    ps = new Translate(0.0, 0.0, 5.0, ps);
    return ps;
}

static Solid*
MakeModelU(double fAng)
{
    Solid *ps, *psSky;
    double fRad, fGX, fGY, fGZ;

    fRad = fAng * D_PI/180.0;
    ps = new Rotate(0.5*D_PI, fRad, 0.0, new Material(0.9, g_sEqualWhite, MAT_PLASTIC, Earth2Texture,
               new Surface(SURF_SMOOTH, 0, new Sphere)),
               "xyz");
    psSky = new Scale (100.0, 100.0, 100.0, new Material(0.9, g_sEarthSky*2.0, MAT_NOSHADE, SkyTexture, new HollowSphere));
    fGX = fGY = fGZ = 0.0;
    fGX = 90.0 * D_PI/180.0;
    psSky = new Rotate(fGX, fGY, fGZ, psSky);
    fGX = fGY = fGZ = 0.0;
    fGZ = 33.0 * D_PI/180.0;
    psSky = new Rotate(fGX, fGY, fGZ, psSky);
    return new Union(ps, psSky);
}

static Solid*
MakeModelV(double fAng)
{
    Solid *ps;
    Spectrum sBlue;
    double fRad;

    fRad = fAng * D_PI/180.0;
    sBlue = (g_sBlue + g_sWhite) * 0.5;
    ps = new Rotate(0.0, fRad, 0.0, new Material(0.9, sBlue, MAT_PLASTIC, SeaScape,
               new Surface(SURF_POLISHED, 0, new Scale(1.2, 1.2, 1.2, new Sphere))));
    ps = new Union(ps, new Translate(0.0, 2.5, -10.0, new Scale(12.0, 1.0, 12.0, 
            new Material(0.9, g_sEqualWhite, MAT_PLASTIC, new Cube))));
    return ps;
}

static Solid*
MakeModelFib(int N, double fAng)           // Fibonacci sphere
{
    Solid *ps, *psN;
    Spectrum sBlue;
    double fRad;
    char szN[16], szPhi[16];
    extern double g_fPhi;

    fRad = fAng * D_PI/180.0;
    sBlue = (g_sBlue + g_sWhite) * 0.5;
    sprintf_s(szN, sizeof(szN), "N = %4d", N);
    sprintf_s(szPhi, sizeof(szPhi), "PHI = %6.3f", g_fPhi);
    /*
    ps = new Scale(0.1, 0.1, 0.1, new Sphere);
    for (i = 0; i < N; i++) {
        FibonacciSphere(rgf, i, N);
        ps = new Union(ps, new Translate(rgf[0], rgf[1], rgf[2], new Scale(0.05, 0.05, 0.05, Plastic(new Sphere))));
    }
    */
    ps = new Fibonacci(N);
    ps = new Rotate(0.0, fRad, 0.0, Plastic2(ps));
    psN = new Translate(1.5-0.6, 1.5, 0.0, new Scale(0.2, 0.2, 0.2, SolidText(0.0, 0.0, szN, JUSTIFIED_LEFT, 0.1)));
    ps = new Union(ps, Plastic2(psN));
    psN = new Translate(-1.5-0.9, 1.5, 0.0, new Scale(0.2, 0.2, 0.2, SolidText(0.0, 0.0, szPhi, JUSTIFIED_LEFT, 0.1)));
    ps = new Union(ps, Plastic2(psN));
    return ps;
}

#define DANGLE (((180.0 - DODECDIHEDRAL)/2.0) * D_PI/180.0)
#define GLASSSCALE 50.0
#define GLASS (0.005 * GLASSSCALE)

static Solid*
BumpyMetal(int nLuster, Solid *ps)
{
    return new Rotate(0.0, s_fRad, 0.0, new Rotate(DANGLE, 0.0, 0.0,
    new Scale(0.5, 0.5, 0.5, new Surface(SURF_POLISHED, 0, new Material(0.9, g_sWhite, nLuster, new Scale(2.0, 2.0, 2.0, ps))))));
}

static Solid*
GemStone(Spectrum sAbsorb, Solid *ps)
{
    return new Rotate(0.0, s_fRad, 0.0, new Rotate(DANGLE, 0.0, 0.0,
        new Surface(SURF_POLISHED, 0, new Material(0.8, REFRACT_SAPPHIRE, sAbsorb, MAT_TRANSPARENT, ColoredGlass, ps))));
}

static Solid*
MakeModelRuby(double fAng)
{
    Solid *ps;
    double fRad, fD, f20;

    fRad = fAng * D_PI/180.0;
    f20 = -20.0 * D_PI/180.0;
    s_fRad = fRad;
    fD = 0.82;              // icosa incribed/circumscribed 0.79465447229176612295553092832759
    fD = 0.79465447229176612295553092832759;
    ps = MakeUnion(
        new Translate(+0.0, +0.0, +0.0, BumpyMetal(MAT_IRON, DodecahedronCage(0.05, 0.1))),
        new Translate(+0.0, +0.0, +0.0, GemStone(g_sRuby*0.4/GLASSSCALE, 
            new Difference(new Scale(fD+GLASS, fD+GLASS, fD+GLASS, new Dodecahedron), new Scale(fD-GLASS, fD-GLASS, fD-GLASS, new Dodecahedron)))),

        new Translate(+2.5, +0.0, +0.0, BumpyMetal(MAT_COPPER, DodecahedronCage(0.05, 0.1))),
        new Translate(+2.5, +0.0, +0.00, GemStone(g_sSapphire*1.5/GLASSSCALE, 
            new Difference(new Scale(fD+GLASS, fD+GLASS, fD+GLASS, new Dodecahedron), new Scale(fD-GLASS, fD-GLASS, fD-GLASS, new Dodecahedron)))),

        new Translate(-2.5, +0.0, +0.0, BumpyMetal(MAT_GOLD, DodecahedronCage(0.05, 0.1))),
        new Translate(-2.5, +0.0, +0.0, GemStone(g_sEmerald*0.4/GLASSSCALE, 
            new Difference(new Scale(fD+GLASS, fD+GLASS, fD+GLASS, new Dodecahedron), new Scale(fD-GLASS, fD-GLASS, fD-GLASS, new Dodecahedron))))

    );
    ps = new Union(ps, new Translate(0.0, 2.5, -10.0, new Scale(12.0, 1.0, 12.0, 
            new Material(0.9, g_sEqualWhite, MAT_PLASTIC, new Cube))));
    ps = new Rotate(f20, 0.0, 0.0, ps);
    ps = new Union(ps, new Scale (1000.0, 1000.0, 1000.0, new Material(0.9, g_sEarthSky*2.0, MAT_NOSHADE, new HollowSphere)));
    return ps;
}

static Solid*
MakeModelLetters(double fAng)
{
    Solid *ps;
    double fRad;

    fRad = fAng * D_PI/180.0;
    ps = SolidText(0.0, 0.0, "Hello World!", JUSTIFIED_CENTER, 0.05);
    return ps;
}

static Solid*
MakeModelSetOperations(double fAng)
{
    Solid *ps, *psU, *psI, *psD;
    double fRad;
    Extent e;

    fRad = fAng * D_PI/180.0;
    psU = new Union(new Translate(+0.5, 0.0, 0.0, new Sphere), new Translate(-0.5, 0.0, 0.0, new Sphere));
    psI = new Intersection(new Translate(+0.5, 0.0, 0.0, new Sphere), new Translate(-0.5, 0.0, 0.0, new Sphere));
    psD = new Difference(new Translate(+0.5, 0.0, 0.0, new Sphere), new Translate(-0.5, 0.0, 0.0, new Sphere));
    ps = MakeUnion(
        new Translate(-3.0, 0.0, 0.0, psU),
        new Translate(+0.0, 0.0, 0.0, psI),
        new Translate(+3.0, 0.0, 0.0, psD)
    );
    psU = new Union(new Translate(+0.5, 0.0, 0.0, 0), new Translate(-0.5, 0.0, 0.0, new Sphere));
    psI = new Intersection(new Translate(+0.5, 0.0, 0.0, 0), new Translate(-0.5, 0.0, 0.0, new Sphere));
    psD = new Difference(new Translate(+0.5, 0.0, 0.0, 0), new Translate(-0.5, 0.0, 0.0, new Sphere));
    e = psU->GetExtent();
    e = psI->GetExtent();
    e = psD->GetExtent();
    ps = MakeUnion(
        ps,
        new Translate(-3.0, +2.5, 0.0, psU),
        new Translate(+0.0, +2.5, 0.0, psI),
        new Translate(+3.0, +2.5, 0.0, psD)
    );
    psU = new Union(new Translate(+0.5, 0.0, 0.0, new Sphere), new Translate(-0.5, 0.0, 0.0, 0));
    psI = new Intersection(new Translate(+0.5, 0.0, 0.0, new Sphere), new Translate(-0.5, 0.0, 0.0, 0));
    psD = new Difference(new Translate(+0.5, 0.0, 0.0, new Sphere), new Translate(-0.5, 0.0, 0.0, 0));
    ps = MakeUnion(
        ps,
        new Translate(-3.0, -2.5, 0.0, psU),
        new Translate(+0.0, -2.5, 0.0, psI),
        new Translate(+3.0, -2.5, 0.0, psD)
    );
    psU = new Union(new Sphere, 0);
    e = psU->GetExtent();
    ps = new BoundingBox(ps);
    return Plastic(ps);
}

static Solid*
MakeModelHyperTest(double fAng)
{
    double fRad;
    Solid *ps;

    fRad = fAng * D_PI/180.0;
    ps = new Rotate(0.0, fRad, 0.0, new Intersection(new Slab, new ConeHyperCylinder(0.4)));
    return ps;
}

#define RHO     28.0
#define SIGMA   10.0
#define BETA    (8.0/3.0)
#define NLORENZ 256

static Vector3 vCrit1, vCrit2;
static Vector3 rgvLorenz[NLORENZ];
static Spectrum rgsLorenz[NLORENZ];

static void
InitializeLorenz()
{
    int i;
    Spectrum sBlue;

    vCrit1 = Vector3(+sqrt(BETA*(RHO - 1.0)), +sqrt(BETA*(RHO - 1.0)), RHO - 1.0);
    vCrit2 = Vector3(-sqrt(BETA*(RHO - 1.0)), -sqrt(BETA*(RHO - 1.0)), RHO - 1.0);
    sBlue = g_sEqualWhite*0.5 + g_sBlue*0.5;
    for (i = 0; i < NLORENZ; i++) {
        rgvLorenz[i] = Vector3(RandomDouble(), RandomDouble(), RandomDouble());
        //rgsLorenz[i] = g_sEqualWhite*0.25 + g_sRed*(0.75*RandomDouble()) + g_sGreen*(0.75*RandomDouble()) + g_sBlue*(0.75*RandomDouble());
        rgsLorenz[i] = sBlue;
        if (i == 0) 
            rgsLorenz[i] = g_sRed;
    }
}

static Solid*
MakeModelLorenz(double fAng, int n)
{
    double dx, dy, dz, dt, f45, f35, fRad;
    Vector3 v;
    Solid *ps, *psUnion = 0;
    Spectrum sBlue;
    int i;

    sBlue = g_sEqualWhite*0.5 + g_sBlue*0.5;
    f45 = -45.0 * D_PI/180.0;
    f35 = +35.3 * D_PI/180.0;
    fRad = fAng * D_PI/180.0;
    dt = 0.01;
    while (--n >= 0) {
        for (i = 0; i < NLORENZ; i++) {
            v = rgvLorenz[i];
            dx = SIGMA*(v.y - v.x);
            dy = v.x*(RHO - v.z) - v.y;
            dz = v.x*v.y - BETA*v.z;
            rgvLorenz[i] = v + Vector3(dx, dy, dz)*dt;
        }
    }
    psUnion = new Union(
        new Translate(vCrit1, new Scale(0.45, 0.45,  0.45, new Sphere)),
        new Translate(vCrit2, new Scale(0.45, 0.45,  0.45, new Sphere)) );
    for (i = 0; i < NLORENZ; i++) {
        ps = new Translate(rgvLorenz[i], new Scale(0.15, 0.15,  0.15, new Sphere));
        if (i == 0)
            ps = new Translate(rgvLorenz[i], new Scale(0.45, 0.45,  0.45, new Sphere));
        ps = new Material(0.9, rgsLorenz[i], MAT_PLASTIC, ps);
        psUnion = new Union(ps, psUnion);
    }

    psUnion = new Translate(0.0, 0.0, -25.0, psUnion);          // center the system at (0,0,0)
    // psUnion = new Rotate(0.0, 0.5*D_PI, 0.0, psUnion);          // y-z plane
    psUnion = new Rotate(0.5*D_PI, 0.0, 0.0, psUnion);              // x-z plane
    psUnion = new Rotate(0.0, fRad, 0.0, psUnion);

    return psUnion;
}

static Vector5 rgvCSystem[NLORENZ];

static void
InitializeCSystem()
{
    int i;
    Spectrum sBlue;
    Vector5 v;

    sBlue = g_sEqualWhite*0.5 + g_sBlue*0.5;
    for (i = 0; i < NLORENZ; i++) {
        v = Vector5(RandomDouble(), RandomDouble(), RandomDouble(), RandomDouble(), RandomDouble());
        rgvCSystem[i] = v/v.Norm();
        //rgsLorenz[i] = g_sEqualWhite*0.25 + g_sRed*(0.75*RandomDouble()) + g_sGreen*(0.75*RandomDouble()) + g_sBlue*(0.75*RandomDouble());
        rgsLorenz[i] = sBlue;
        if (i == 0) 
            rgsLorenz[i] = g_sRed;
    }
}

static Solid*
MakeModelCSystem(double fAng, int n)
{
    double dx, dy, dz, da, db, dt, f45, f35, fRad;
    Vector5 v, dv;
    Vector3 v3;
    Solid *ps, *psUnion = 0;
    Spectrum sBlue;
    int i, j;

    sBlue = g_sEqualWhite*0.5 + g_sBlue*0.5;
    f45 = -45.0 * D_PI/180.0;
    f35 = +35.3 * D_PI/180.0;
    fRad = fAng * D_PI/180.0;
    dt = 0.01;
    while (--n >= 0) {
        for (i = 0; i < NLORENZ; i++) {
            v = rgvCSystem[i];
		    for (j = 0; j < 5; j++)
			    dv.rgf[j] = v.rgf[(j+1)%5]*v.rgf[(j+2)%5]
			          + v.rgf[(j+4)%5]*v.rgf[(j+3)%5]
			      - 2.0*v.rgf[(j+1)%5]*v.rgf[(j+4)%5];
		    for (j = 0; j < 5; j++)
			    v.rgf[j] += dv.rgf[j] * dt;
            v = v / v.Norm();
            rgvCSystem[i] = v;
        }
    }
    psUnion = new Translate(0.0, 0.0, 0.0, new Scale(0.5, 0.5,  0.5, new Sphere));
    for (i = 0; i < NLORENZ; i++) {
        v3 = Vector3(rgvCSystem[i].x, rgvCSystem[i].y, rgvCSystem[i].z);
        ps = new Translate(v3/Norm(v3), new Scale(0.05, 0.05,  0.05, new Sphere));
        if (i == 0)
            ps = new Translate(rgvCSystem[i].x, rgvCSystem[i].y, rgvCSystem[i].z, new Scale(0.1, 0.1,  0.1, new Sphere));
        ps = new Material(0.9, rgsLorenz[i], MAT_PLASTIC, ps);
        psUnion = new Union(ps, psUnion);
    }

   // psUnion = new Translate(0.0, 0.0, -25.0, psUnion);          // center the system at (0,0,0)
    // psUnion = new Rotate(0.0, 0.5*D_PI, 0.0, psUnion);          // y-z plane
   // psUnion = new Rotate(0.5*D_PI, 0.0, 0.0, psUnion);              // x-z plane
   // psUnion = new Rotate(0.0, fRad, 0.0, psUnion);

    return psUnion;
}

#define RND (10.0*RandomDouble() - 5.0)

static Solid*
MakeModel_TestTop(double fAng)
{
    Solid *ps;

    ps = new BoundingBox (MakeTopUnion(
        new BoundingBox(MakeTopUnion(
            new Translate(RND, -RND, RND, new Scale(0.5, 1.0, 1.0, new Sphere)),
            new Translate(RND, -RND, RND, new Scale(0.5, 1.0, 1.0, new Sphere)),
            new Translate(RND, -RND, RND, new Scale(0.5, 1.0, 1.0, new Sphere)),
            new Translate(RND, -RND, RND, new Sphere)
        )),
        new BoundingBox(MakeTopUnion(
            new Translate(RND, -RND, RND, new Scale(0.5, 1.0, 1.0, new Sphere)),
            new Translate(RND, -RND, RND, new Scale(0.5, 1.0, 1.0, new Sphere)),
            new Translate(RND, -RND, RND, new Scale(0.5, 1.0, 1.0, new Sphere)),
            new Translate(RND, -RND, RND, new Sphere)
        )),
        new BoundingBox(MakeTopUnion(
            new Translate(RND, -RND, RND, new Scale(0.5, 1.0, 1.0, new Sphere)),
            new Translate(RND, -RND, RND, new Scale(0.5, 1.0, 1.0, new Sphere)),
            new Translate(RND, -RND, RND, new Scale(0.5, 1.0, 1.0, new Sphere)),
            new Translate(RND, -RND, RND, new Sphere)
        ))
    ));
    return ps;
}

static void
MakeLights(RayTracer *ps)
{
    Spectrum s6500, sBlue;

    sBlue = (g_sBlue + g_sEqualWhite) * 0.5;
    s6500 = BlackBody(6500.0);
    printf("6500 K light: ");
    s6500.sRGB().Print();
     ps->plAmbient = new AmbientLight(s6500 * 0.02);
    //ps->plAmbient = new HavercosineLight(sBlue * 0.02, s6500 * 0.02);
   // ps->plPointLights = new PointLight(Point3(10.0, -10.0, 10.0), s6500 * 10000.0);
    ps->plPointLights = new AreaLight(1.0, Point3(10.0, -10.0, 10.0), s6500 * 10000.0);
    ps->plPointLights = new PointLight(Point3(-10.0, -10.0, -10.0), s6500 * 500.0, ps->plPointLights);
}

static void
MakeLightsLorenz(RayTracer *ps)
{
    Spectrum s6500, sBlue;

    sBlue = (g_sBlue + g_sEqualWhite) * 0.5;
    s6500 = BlackBody(6500.0);
    printf("6500 K light: ");
    s6500.sRGB().Print();
     ps->plAmbient = new AmbientLight(s6500 * 0.02);
    //ps->plAmbient = new HavercosineLight(sBlue * 0.02, s6500 * 0.02);
   // ps->plPointLights = new PointLight(Point3(10.0, -10.0, 10.0), s6500 * 10000.0);
    ps->plPointLights = new AreaLight(1.0, Point3(100.0, -100.0, 100.0), s6500 * 1000000.0);
    ps->plPointLights->nNoShadow = 1;
    //ps->plPointLights = new PointLight(Point3(-10.0, -10.0, -10.0), s6500 * 500.0, ps->plPointLights);
}

static void
MakeLights3(RayTracer *ps, double x)
{
    Spectrum s6500;

    s6500 = BlackBody(6500.0);
    ps->plAmbient = new AmbientLight(s6500 * 0.02);
    ps->plPointLights = new PointLight(Point3(10.0, -10.0, 10.0), s6500 * 10000.0);
    ps->plPointLights = new PointLight(Point3(x, 0.0, 0.0), s6500 * 100.0, ps->plPointLights);
}

void
TestLightRayBug()
{
    RayTracer scene;
    Ray ray;
    Radiance rad;

    g_nVerbose = 1;
    printf("Test LightRay bug\n");
    ray = Ray(Point3(0.0, 0.0, 5.0), Vector3(0.0, 0.0, -1.0));
    MakeLights3(&scene, 29.0/10 - 4.0);
    scene.psModel = new Material(0.9, REFRACT_PYREX, g_sEmerald * 0.01, MAT_TRANSPARENT, ColoredGlass, new Surface(SURF_POLISHED, new Sphere));
    rad = scene.VisualRay(ray);
    printf("\n");
    MakeLights3(&scene, 31.0/10 - 4.0);
    //scene.psModel = new Material(0.9, REFRACT_PYREX, g_sEmerald * 0.01, MAT_TRANSPARENT, ColoredGlass, new Surface(SURF_POLISHED, new Sphere));
    rad = scene.VisualRay(ray);
}

static void
MakeLights2(RayTracer *ps, double fAng)
{
    Spectrum s6500;
    double fRad, fSin, fCos;

    fRad = fAng * D_PI/180.0;
    fSin = sin(fRad);
    fCos = cos(fRad);
    s6500 = BlackBody(6500.0);
    ps->plAmbient = new AmbientLight(s6500 * 0.005);
    ps->plPointLights = new PointLight(Point3(10.0*fSin, -10.0*fCos, 0.0), s6500 * 4000.0);
}

int g_nTestFrame, g_nTestRow, g_nTestColumn;
extern double g_fPhi;
//
// 19.5  	- baseline Model6
// 18.47   - phLeft == 0 set ops
//
void
TestRayTraceScene()
{
    RayTracer scene;
    Image im;
    Video v;
    ML_TimingInfo ti;
    Ray ray;
    Solid *ps;
    double fAng, rc, g, d, fCycle, ff, f, fAve;
    int i, n;
    float rgf[3];
    char sz[32];

    rc = 0.5;
    g = 0.1;
    d = 1.0 + rc + rc + g + g;    
    printf("Start, hits = %d\n", g_nHitCount);
    fAve = 0;
    for (i = 0; i < 0; i++) {
        rgf[0] = RandomFloat()* 10.0f;
        rgf[1] = RandomFloat()* 10.0f;
        rgf[2] = RandomFloat()* 10.0f;
        printf("fractal %f\n", f = FractalNoise(rgf, 3));
        fAve += f;
    }
    if (i) printf("ave fractal= %f\n", fAve/double(i));
    for (i = 0; i < 1; i++)
        s_rgs[i] = RNDSPEC;
    for (i = 0; i < 20; i++)
        s_rgs[i] = RNDSPEC;
    for (i = 0; i < 20; i++)
        s_rgs[i] = RNDSPEC;
    s_rgs[9] = RGBtoSpectrum(g_rgbCopper);
    s_rgs[7] = RGBtoSpectrum(g_rgbCopper);
    s_rgs[5] = RGBtoSpectrum(g_rgbCopper);
    s_rgs[3] = RGBtoSpectrum(g_rgbCopper);
    s_rgs[1] = RGBtoSpectrum(g_rgbCopper);
    s_rgs[9] = g_sCopper;
    s_rgs[7] = g_sCopper;
    s_rgs[5] = g_sCopper;
    s_rgs[3] = g_sCopper;
    s_rgs[1] = g_sCopper;

    //
    //  Jitter 4.0 - 19.17  (Model8)
    //  Blue   4.0 - 29.34  (SplatRGB)
    //  Blue   4.0 - 23.88  (no calls to KaiserFilter)
    //  Blue   4.0 - 19.0   (FastSplatRGB)
    //  Fast   4.0 - 20.0   (FastSplatRGB & ImageResample)
    //
    //MakeLights(&scene);
    MakeLightsLorenz(&scene);
    scene.psSample = new TestSamplingTile;        // g_nVerbose = 1;
    scene.pcCamera = new StereographicCamera(Point3(+0.0, 0.0, 8.0), 70.0);
    im.NewImage(1920, 1080, 3);
    v.NewVideo("TestRayTracer.mpg", im);
    fCycle = 5.0*d;
    ray = Ray(Point3(0.0, 0.0, 5.0), Vector3(0.0, 0.0, -1.0));  // cvPerturbed doesn't get transformed properly
//TestLightRayBug();
    InitializeCSystem();
    ML_StartTiming(ti);
    for (fAng = 0.0; fAng < 360; fAng += 1.0/3.0) {
        if (fAng < 0.0)
            f = 0.0;
        else
            f = fAng;
        ff = 1.0 + fAng/180.0;
        ff = 1.0;
        n = (fAng == 0.0) ? 1000 : 1;
        g_nTestFrame = int(fAng);
        // MakeLights3(&scene, fAng);       // fAng = -4.0; fAng <= 4.0; fAng += 0.1
        // scene.psModel = MakeModelK(fAng);
        ps = MakeModelCSystem(fAng, n);
        scene.psModel = ps;
        scene.psModel = scene.psModel->Optimize();
        scene.Render(im, 4.0);
        v.WriteFrame(im);
        sprintf_s(sz, sizeof(sz), "Frame%d.bmp", int(fAng));
        //im.WriteBMP(sz);
        printf("fAng = %f\n", fAng);
       
    }
    ML_StopTiming(ti);
    v.Close();
    im.WriteBMP("RayTraceScene.bmp");
    ML_ReportTiming(ti);
    printf("End, hits = %d\nNew Hits = %d\n", g_nHitCount, g_nNewHits);
    printf("End, texturelists = %d\nNew Lists = %d\n", g_nTextureListCount, g_nTextureLists);
}