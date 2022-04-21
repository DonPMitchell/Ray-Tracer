//
//  Venera landers and Venus
//  D.P. Mitchell  2021/04/09.
//
#include "stdafx.h"
#include "RayTracer2020.h"
#include "TextureMaps.h"
#include "Utilities.h"
#include "Video.h"

static void
ColorizeVenusMap(const Spectrum &s)
{
    Image im;
    Spectrum sGround;
    DisplayRGB rgb, rgbGround, rgbWhite, rgbBlack;
    double f, fGround, fAverageVenus, fN;
    int i, j;

    sGround = g_sVenera14Ground / g_sVenera14Sky;       // should be called with sGround as the argument
    f = sGround.rgf[5] - sGround.rgf[4];
    sGround.rgf[3] = float(sGround.rgf[4] - f);
    sGround.rgf[2] = float(sGround.rgf[4] - 2.0*f);
    sGround.rgf[1] = float(sGround.rgf[4] - 3.0*f);
    sGround.rgf[0] = float(sGround.rgf[4] - 4.0*f);

    printf("Colorize Venus Map:\n");
    rgbGround = s.sRGB();
    fGround = rgbGround.Luminance();
    printf("VenusGround luminance = %g\n", fGround);
    rgbWhite = ML_White*0.75 + rgbGround*0.25;
    rgbBlack = ML_Black;
    fAverageVenus = fN = 0.0;
    im.ReadBMP("VenusColorMap_BW.bmp");
    for (i = 0; i < im.m_nWidth; i++) {
        for (j = 0; j < im.m_nHeight; j++) {
            rgb = im.GetRGB(i, j);
            f = rgb.Luminance();
            if (f == 0.0)
                continue;
            fAverageVenus += f;
            fN += 1.0;
        }
    }
    fAverageVenus /= fN;
    //rgbGround = rgbGround * fAverageVenus/fGround;
    printf("Average VenusMap = %g\n", fAverageVenus);

return;


    for (i = 0; i < im.m_nWidth; i++) {
        for (j = 0; j < im.m_nHeight; j++) {
            rgb = im.GetRGB(i, j);
            f = rgb.Luminance();
            if (f < fAverageVenus)
                rgb = rgbBlack + (rgbGround - rgbBlack)*(f/(fAverageVenus - 0.0));
            else
                rgb = rgbGround + (rgbWhite - rgbGround)*((f - fAverageVenus)/(1.0 - fAverageVenus));
            im.SetRGB(rgb, i, j);
        }
    }
    im.WriteBMP("VenusColorMap2.bmp");
}

static void
ColorizeMercuryMap(const Spectrum &s)
{
    Image im;
    DisplayRGB rgb, rgbGround, rgbWhite, rgbBlack;
    double f, fGround, fAverage, fN;
    int i, j;

    printf("Colorize Venus Map:\n");
    rgbGround = s.sRGB();
    fGround = rgbGround.Luminance();
    printf("MercuryGround luminance = %g\n", fGround);
    rgbWhite = ML_White*0.75 + rgbGround*0.25;
    rgbBlack = ML_Black;
    fAverage = fN = 0.0;
    im.ReadBMP("MercuryColorMap.bmp");
    for (i = 0; i < im.m_nWidth; i++) {
        for (j = 0; j < im.m_nHeight; j++) {
            rgb = im.GetRGB(i, j);
            f = rgb.Luminance();
            if (f == 0.0)
                continue;
            fAverage += f;
            fN += 1.0;
        }
    }
    fAverage /= fN;
    //rgbGround = rgbGround * fAverageVenus/fGround;
    printf("Average MercuryMap = %g\n", fAverage);
    for (i = 0; i < im.m_nWidth; i++) {
        for (j = 0; j < im.m_nHeight; j++) {
            rgb = im.GetRGB(i, j);
            f = rgb.Luminance();
            if (f < fAverage)
                rgb = rgbBlack + (rgbGround - rgbBlack)*(f/(fAverage - 0.0));
            else
                rgb = rgbGround + (rgbWhite - rgbGround)*((f - fAverage)/(1.0 - fAverage));
            im.SetRGB(rgb, i, j);
        }
    }
    im.WriteBMP("MercuryColorMap2.bmp");
}

static Solid*
ColorChips()
{
    Solid *ps, *psRow1;
    Spectrum sVenusSky, sGround, s6500;
    double f;

    sVenusSky = g_sEqualWhite;
    s6500 = BlackBody(6500.0);
    sGround = g_sVenera14Ground / g_sVenera14Sky;
    //sGround = sGround * s6500;
    f = sGround.rgf[5] - sGround.rgf[4];
    sGround.rgf[3] = float(sGround.rgf[4] - f);
    sGround.rgf[2] = float(sGround.rgf[4] - 2.0*f);
    sGround.rgf[1] = float(sGround.rgf[4] - 3.0*f);
    sGround.rgf[0] = float(sGround.rgf[4] - 4.0*f);
    sGround.Print();
    sGround.sRGB().Print();
    ColorizeVenusMap(sGround);
    //ColorizeMercuryMap(g_sMercuryAlbedo * 0.1);

    psRow1 = MakeUnion(
        new Translate(-8.75, 1.25, 0.0, new Material(0.9, g_sVenera14Ground/75.0, MAT_PLASTIC, new Cube)),
        new Translate(-6.25, 1.25, 0.0, new Material(0.9, g_sVenera14Sky/75.0, MAT_PLASTIC, new Cube)),
        new Translate(-3.75, 1.25, 0.0, new Material(0.9, g_sVenera13Sky/75.0, MAT_PLASTIC, new Cube)),
        new Translate(-1.25, 1.25, 0.0, new Material(0.9, g_sVenera11Sky/75.0, MAT_PLASTIC, new Cube)),
        new Translate(+1.25, 1.25, 0.0, new Material(0.9, sGround, MAT_PLASTIC, new Cube)),
        new Translate(+3.75, 1.25, 0.0, new Material(0.9, g_sEarthSky, MAT_PLASTIC, new Cube)),
        new Translate(+6.25, 1.25, 0.0, new Material(0.9, g_sMarsAlbedo, MAT_PLASTIC, new Cube)),
        new Translate(+8.75, 1.25, 0.0, new Material(0.9, g_sMercuryAlbedo*0.1, MAT_PLASTIC, new Cube))
    );
    ps = MakeUnion(psRow1);
    return ps;
}

static Solid*
Titanium(Solid *ps)
{
    return new Material(0.8, REFRACT_WATER, g_sEqualWhite, MAT_TITANIUM, 0,
               new Surface(SURF_BURNISHED, 0, ps));
}

static Solid*
Titanium2(Solid *ps)
{
    return new Material(0.8, REFRACT_WATER, g_sEqualWhite, MAT_OSMIUM, 0,
               new Surface(SURF_SMOOTH, 0, ps));
}

static Solid*
Titanium3(Solid *ps)
{
    return new Material(0.8, REFRACT_WATER, g_sEqualWhite, MAT_TITANIUM, 0,
               new Surface(SURF_SMOOTH, 0, ps));
}

static Solid*
Plastic(Solid *ps)
{
    Spectrum sBlue;

    sBlue = (g_sEqualWhite + g_sBlue) * 0.5;
    return new Material(0.7, REFRACT_WATER, sBlue, MAT_PLASTIC, 0,
               new Surface(SURF_POLISHED, 0, ps));
}

static Solid*
Plastic2(Solid *ps)
{
    Spectrum sBlue;

    sBlue = (g_sEqualWhite + g_sOrange) * 0.5;
    return new Material(0.7, REFRACT_WATER, sBlue, MAT_PLASTIC, 0,
               new Surface(SURF_POLISHED, 0, ps));
}

static Solid*
Quartz(Solid *ps)
{
    Spectrum sBlue;

    sBlue = (g_sEqualWhite + g_sOrange) * 0.5;
    return new Material(1.0, REFRACT_QUARTZ, g_sBlack, MAT_TRANSPARENT, ColoredGlass,
               new Surface(SURF_POLISHED, 0, ps));
}

static Solid*
KG_25(Solid *ps)
{
    Spectrum sYellow;

    sYellow = (g_sEqualWhite + g_sYellow) * 0.5;
    return new Material(0.9, REFRACT_WATER, sYellow, MAT_PLASTIC, 0,
               new Surface(SURF_MATTE, 0, ps));
}

static Solid*
PTKV_260(Solid *ps)
{
    Spectrum sBlue;

    sBlue = (g_sEqualWhite*0.7 + g_sBlue*0.3);
    return new Material(0.7, REFRACT_WATER, sBlue, MAT_PLASTIC, 0,
               new Surface(SURF_MATTE, 0, ps));
}

static Solid*
Slice(Solid *ps)
{
    return new Difference(ps, new Scale(25.0, 25.0, 25.0, new Translate(0.0, 0.0, 1.0, new Cube)));
}

static Solid*
VeneraColorChips()
{
    Solid *ps;

    ps = new BoundingBox(MakeUnion(
        new Translate(-5.0, 0.0, 0.0, new Material(0.9, g_sEqualWhite, MAT_PLASTIC, new Scale(1.0, 1.0, 0.1, new Cube))),
        new Translate(-3.0, 0.0, 0.0, new Material(0.9, g_sCoolGray, MAT_PLASTIC, new Scale(1.0, 1.0, 0.1, new Cube))),
        new Translate(-1.0, 0.0, 0.0, new Material(0.9, g_sCoolRed,  MAT_PLASTIC, new Scale(1.0, 1.0, 0.1, new Cube))),
        new Translate(+1.0, 0.0, 0.0, new Material(0.9, g_sCoolGrn,  MAT_PLASTIC, new Scale(1.0, 1.0, 0.1, new Cube))),
        new Translate(+3.0, 0.0, 0.0, new Material(0.9, g_sCoolBlu,  MAT_PLASTIC, new Scale(1.0, 1.0, 0.1, new Cube)))
    ));
    ps = new Scale(0.25, 0.25, 0.25, ps);
    return ps;
}

static Solid*
EnvironmentVenus(Solid *ps, RayTracer &scene)
{
    double f20, f;
    Spectrum sGround, s6500;
    DisplayRGB rgb;

    s6500 = BlackBody(6500.0);
    sGround = g_sVenera14Ground / g_sVenera14Sky;
    f = sGround.rgf[5] - sGround.rgf[4];
    sGround.rgf[3] = float(sGround.rgf[4] - f);
    sGround.rgf[2] = float(sGround.rgf[4] - 2.0*f);
    sGround.rgf[1] = float(sGround.rgf[4] - 3.0*f);
    sGround.rgf[0] = float(sGround.rgf[4] - 4.0*f);
    scene.plAmbient = new HavercosineLight(g_sVenera13Sky * 0.15 /250.0, g_sVenera14Ground * 0.15 /250.0);
    scene.plPointLights = new PointLight(Point3(-10.0, -10.0, -10.0), s6500 * 0.0);
    f20 = 10.0 * D_PI/180.0;
    ps = new Union(ps, new Translate(0.0, 3.0, -10.0, new Scale(20.0, 1.0, 20.0, 
                new Material(0.9, sGround*4.0, MAT_PLASTIC, new Cube))));
    ps = new Union(ps, new Scale (1000.0, 1000.0, 1000.0, new Material(0.9, g_sVenera13Sky*1.5, MAT_NOSHADE, new HollowSphere)));
    ps = new Rotate(-f20, 0.0, 0.0, ps);
    return ps;
}

static Solid*
EnvironmentLab(Solid *ps, RayTracer &scene)
{
    double f20, f;
    Spectrum sGround, s6500;
    DisplayRGB rgb;

    s6500 = BlackBody(6500.0);
    sGround = g_sVenera14Ground / g_sVenera14Sky;
    f = sGround.rgf[5] - sGround.rgf[4];
    sGround.rgf[3] = float(sGround.rgf[4] - f);
    sGround.rgf[2] = float(sGround.rgf[4] - 2.0*f);
    sGround.rgf[1] = float(sGround.rgf[4] - 3.0*f);
    sGround.rgf[0] = float(sGround.rgf[4] - 4.0*f);
    scene.plAmbient = new HavercosineLight(g_sEarthSky*0.025, g_sEqualWhite*0.00625);
    scene.plPointLights = new PointLight(Point3(9.0, -6.5, 12.0), s6500 * 8000.0 * 1.0);
    scene.plPointLights = new PointLight(Point3(-10.0, -10.0, -10.0), s6500 * 100.0, scene.plPointLights);
    ps = new Union(ps, new Translate(0.0, 1.0, -10.0, new Scale(12.0, 0.1, 12.0, 
                new Material(0.9, g_sEqualWhite, MAT_PLASTIC, new Cube))));
    ps = new Union(ps, new Scale (1000.0, 1000.0, 1000.0, new Material(0.9, g_sEarthSky*2.0, MAT_NOSHADE, new HollowSphere)));

    f20 = 10.0 * D_PI/180.0;
    ps = new Rotate(-f20, 0.0, 0.0, ps);
    return ps;
}

#define RINGRADIUS  1.0
#define RINGWIDE    0.17
#define RINGHIGH    (RINGWIDE*0.75)
#define RINGBOLTS   120.0 
#define RINGHOLES   60.0

#define FLANGE  0.075
#define FLANGEX 0.022
#define FTHICK  0.0025
#define FRAD    (FLANGEX*0.5)
#define FBOLT   0.5*0.012

#define PHULL   0.41
#define PHIGH   0.12
#define ITHICK  0.11
#define PF3     -atan2(FLANGEX, PHULL)

#define RSTRUT  (0.02)

static Solid*
HalfSlaboid(double x, double y, double z, double r)
{
    Solid *ps;

    ps = new Translate(0.0, y-r, 0.0, Slaboid(x, y+r, z, r));
    ps = new Difference(ps, new Scale(1.0, 1.01*y, 1.0, new Translate(0.0, -1.0, 0.0, new Cube)));
    ps = new Union(ps, new Rotate(0.0, 0.5*D_PI, 0.0, new Scale(z, z, x, new Cylinder)));
    return ps;
}

static Solid*
Flange1(double fAng)
{
    Solid *ps, *ps2, *ps3;

    ps2 = new Translate(0.0, FLANGE-FTHICK, -(FLANGEX-FTHICK), new Rotate(0.0, 0.5*D_PI, 0.0, Slaboid(FLANGEX, FLANGE, FTHICK, FRAD)));
    ps2 = new Rotate(+0.25*D_PI, 0.0, 0.0, new Scale(1.0, 1.0, 1.6, new Rotate(-0.25*D_PI, 0.0, 0.0, ps2)));
    ps3 = new Scale(2.0*FTHICK, 2.0*FLANGE, 2.0*FLANGEX, new Translate(0.0, 1.0, -1.0, new Cube));
    ps = MakeUnion(
        HalfSlaboid(FLANGEX, FLANGE, FTHICK, FRAD),
        new Translate(FLANGEX-FRAD, FRAD, 0.0, new Scale(FBOLT, FBOLT, FBOLT, new Prism(6))),
        new Translate(FLANGEX-FRAD, 2.0*FLANGE-FRAD, 0.0, new Scale(FBOLT, FBOLT, FBOLT, new Prism(6))),
        new Translate(-FLANGEX+FRAD, FRAD, 0.0, new Scale(FBOLT, FBOLT, FBOLT, new Prism(6))),
        new Translate(-FLANGEX+FRAD, 2.0*FLANGE-FRAD, 0.0, new Scale(FBOLT, FBOLT, FBOLT, new Prism(6))),

        new Rotate(-0.5*D_PI, 0.0, 0.0, MakeUnion(
            HalfSlaboid(FLANGEX, 0.7*FLANGE, FTHICK, FRAD),
            new Translate(FLANGEX-FRAD, FRAD, 0.0, new Scale(FBOLT, FBOLT, FBOLT, new Prism(6))),
            new Translate(FLANGEX-FRAD, 1.4*FLANGE-FRAD, 0.0, new Scale(FBOLT, FBOLT, FBOLT, new Prism(6))),
            new Translate(-FLANGEX+FRAD, FRAD, 0.0, new Scale(FBOLT, FBOLT, FBOLT, new Prism(6))),
            new Translate(-FLANGEX+FRAD, 1.4*FLANGE-FRAD, 0.0, new Scale(FBOLT, FBOLT, FBOLT, new Prism(6)))
        )),

        new Difference(new Translate(0.0, -0.4*FLANGE, 0.55*FLANGE, ps2), ps3)
    );
    fAng = fAng * D_PI/180.0;
   // ps = Titanium(ps);
    ps = Titanium3(ps);
    ps = new Rotate(0.0, 0.0, fAng, new Translate(0.0, 1.0-RINGWIDE-2.0*FTHICK, RINGHIGH + 2.0*FTHICK, ps));
    return new BoundingBox(ps);
}

static Solid*
Flange2(double fAng)
{
    Solid *ps, *ps2, *ps3;

    ps2 = new Translate(0.0, FLANGE-FTHICK, -(FLANGEX-FTHICK), new Rotate(0.0, 0.5*D_PI, 0.0, Slaboid(FLANGEX, FLANGE, FTHICK, FRAD)));
    ps2 = new Rotate(+0.25*D_PI, 0.0, 0.0, new Scale(1.0, 1.0, 1.6, new Rotate(-0.25*D_PI, 0.0, 0.0, ps2)));
    ps3 = new Scale(2.0*FTHICK, 2.0*FLANGE, 2.0*FLANGEX, new Translate(0.0, 1.0, -1.0, new Cube));
    ps = MakeUnion(
        HalfSlaboid(FLANGEX, FLANGE, FTHICK, FRAD),
        new Translate(FLANGEX-FRAD, FRAD, 0.0, new Scale(FBOLT, FBOLT, FBOLT, new Prism(6))),
        new Translate(FLANGEX-FRAD, 2.0*FLANGE-FRAD, 0.0, new Scale(FBOLT, FBOLT, FBOLT, new Prism(6))),
        new Translate(-FLANGEX+FRAD, FRAD, 0.0, new Scale(FBOLT, FBOLT, FBOLT, new Prism(6))),
        new Translate(-FLANGEX+FRAD, 2.0*FLANGE-FRAD, 0.0, new Scale(FBOLT, FBOLT, FBOLT, new Prism(6))),

        new Rotate(-0.5*D_PI, 0.0, 0.0, MakeUnion(
            HalfSlaboid(FLANGEX, 0.7*FLANGE, FTHICK, FRAD),
            new Translate(FLANGEX-FRAD, FRAD, 0.0, new Scale(FBOLT, FBOLT, FBOLT, new Prism(6))),
            new Translate(FLANGEX-FRAD, 1.4*FLANGE-FRAD, 0.0, new Scale(FBOLT, FBOLT, FBOLT, new Prism(6))),
            new Translate(-FLANGEX+FRAD, FRAD, 0.0, new Scale(FBOLT, FBOLT, FBOLT, new Prism(6))),
            new Translate(-FLANGEX+FRAD, 1.4*FLANGE-FRAD, 0.0, new Scale(FBOLT, FBOLT, FBOLT, new Prism(6)))
        )),
        new Rotate(0.75*D_PI, 0.0, 0.0, HalfSlaboid(1.5*FLANGEX, 0.3*FLANGE, FTHICK, FRAD)),
        new Difference(new Translate(0.0, -0.4*FLANGE, 0.55*FLANGE, ps2), ps3)
    );
    fAng = fAng * D_PI/180.0;
   // ps = Titanium(ps);
    ps = Titanium3(ps);
    ps = new Rotate(0.0, 0.0, fAng, new Translate(0.0, 1.0-RINGWIDE-2.0*FTHICK, RINGHIGH + 2.0*FTHICK, ps));
    return new BoundingBox(ps);
}

static Solid*
Flange3(double fAng)
{
    Solid *ps, *ps2;

    ps = new Scale(PHULL+FTHICK, PHULL+FTHICK, PHULL, new Cylinder);
    ps2 = new Translate(0.0, +PHULL, 0.0, new Rotate(-0.25*D_PI, 0.0, 0.0,
            new Translate(0.0, -PHULL, 0.0, new Scale(PHULL+0.3*FLANGE, PHULL+0.3*FLANGE, FTHICK, new Cylinder))));
    ps = new Difference(new Union(ps, ps2), new Scale(PHULL, PHULL, PHULL+FTHICK, new Cylinder));
    ps = new Intersection(ps, new Translate(0.0, +PHULL, 0.0, new Rotate(0.5*D_PI, 0.0, 0.0, Slaboid(1.5*FLANGEX, FLANGEX, FLANGEX, FRAD))));
    ps = new BoundingBox(ps);
    ps = new Rotate(PF3, 0.0, 0.0, ps);
    ps = new Rotate(0.0, 0.0, fAng * D_PI/180.0, ps);
    return Titanium3(ps);
}

static Solid*
LandingRing()
{
    Solid *ps, *psDiff, *psCrushPad, *psFlanges, *psBolts;
    double fRad, f;

    ps = new Translate(0.0, 0.0, 1.0, new Difference(new Scale(1.0+0.25*RINGHIGH, 1.0+0.25*RINGHIGH, 1.0, new Cylinder), 
                                                     new Scale(1.0-RINGWIDE, 1.0-RINGWIDE, 1.01, new Cylinder)) );
    //
    //  Flat top of landing ring is z = RINGHIGH
    //
    ps = new Scale(1.0, 1.0, 0.5*RINGHIGH, ps);
    ps = new Intersection(ps, new Scale(1.0+0.25*RINGHIGH, 1.0+0.25*RINGHIGH, 4.0, new Translate(0.0, 0.0, 1.0, new ConeUnit)));
    psDiff = new Difference(new Slab, new Scale(1.0-RINGWIDE+0.01, 1.0-RINGWIDE+0.01, 1.01, new Cylinder));
    ps = new Difference(ps, new Scale(1.0, 1.0, 0.5*RINGHIGH, psDiff));
    fRad = 2.0*D_PI/RINGBOLTS;
    
    for (f = 0.0; f < 2.0*D_PI; f += 6.0*fRad) {
        psBolts = new BoundingBox(MakeUnion(
            new Rotate(0.0, 0.0, 0.0*fRad+f, new Translate(1.0, 0.0, 0.0, new Rotate(0.0, 0.5*D_PI, 0.0, 
                new Scale(FBOLT, FBOLT, FBOLT, new Prism(6))))),
            new Rotate(0.0, 0.0, 1.0*fRad+f, new Translate(1.0, 0.0, 0.0, new Rotate(0.0, 0.5*D_PI, 0.0, 
                new Scale(FBOLT, FBOLT, FBOLT, new Prism(6))))),
            new Rotate(0.0, 0.0, 2.0*fRad+f, new Translate(1.0, 0.0, 0.0, new Rotate(0.0, 0.5*D_PI, 0.0, 
                new Scale(FBOLT, FBOLT, FBOLT, new Prism(6))))),
            new Rotate(0.0, 0.0, 3.0*fRad+f, new Translate(1.0, 0.0, 0.0, new Rotate(0.0, 0.5*D_PI, 0.0, 
                new Scale(FBOLT, FBOLT, FBOLT, new Prism(6))))),
            new Rotate(0.0, 0.0, 4.0*fRad+f, new Translate(1.0, 0.0, 0.0, new Rotate(0.0, 0.5*D_PI, 0.0, 
                new Scale(FBOLT, FBOLT, FBOLT, new Prism(6))))),
            new Rotate(0.0, 0.0, 5.0*fRad+f, new Translate(1.0, 0.0, 0.0, new Rotate(0.0, 0.5*D_PI, 0.0, 
                new Scale(FBOLT, FBOLT, FBOLT, new Prism(6)))))
        ));
        ps = Titanium(new Union(ps, new Translate(0.0, 0.0, RINGHIGH-3.0*FBOLT, psBolts)));
    }
    
    //
    //  z = 0 plane thru center of crush pad torus
    //
    psCrushPad = Titanium(new Difference(new Scale(1.0-0.4*RINGWIDE, 1.0-0.4*RINGWIDE, 1.25, new Torus(0.6*RINGWIDE)),
                                         new Scale(1.0-0.4*RINGWIDE, 1.0-0.4*RINGWIDE, 1.25, new Torus(0.5*RINGWIDE))));
                           
    psFlanges = MakeUnion(
        new BoundingBox(MakeUnion(
            Flange1(0.0),
            Flange2(30.0),
            Flange1(60.0),
            Flange2(90.0))),
        new BoundingBox(MakeUnion(
            Flange1(120.0),
            Flange2(150.0),
            Flange1(180.0),
            Flange2(210.0))),
        new BoundingBox(MakeUnion(
            Flange1(240.0),
            Flange2(270.0),
            Flange1(300.0),
            Flange2(330.0)))
    );
    
    ps = MakeUnion(ps, psCrushPad, psFlanges);
    return new BoundingBox(ps);
}

static Solid*
Body()
{
    Solid *psHull, *psInsulate, *psFlanges;
    double fRad, f50;

    psHull = new Union(
            new Scale(PHULL, PHULL, PHULL, new Sphere),
            new Scale(PHULL+0.02, PHULL+0.02, 0.01, new Cylinder)
            );
    for (fRad = 0; fRad < 1.9999*D_PI; fRad += 10.0 *D_PI/180.0)
        psHull = new Union(psHull, 
            new Rotate(0.0, 0.0, fRad, new Translate(0.0, PHULL+0.01, 0.0, new Scale(0.005, 0.005, 0.015, new Prism(6)))));
    psHull = new Union(Titanium(psHull), KG_25(new Scale(PHULL-0.01, PHULL-0.01, PHULL-0.01, new Sphere)));
    psHull = new Difference(psHull, new Scale(PHULL-0.04, PHULL-0.04, PHULL-0.04, new Sphere));

    f50 = (90.0 - 50.0)*D_PI/180.0;
    psHull = new Difference(psHull, new Rotate(0.0, +f50, 0.0, new Translate(0.0, 0.0, PHULL, new Scale(0.02, 0.02, 0.05, new Cylinder))));
    psHull = new Difference(psHull, new Rotate(0.0, -f50, 0.0, new Translate(0.0, 0.0, PHULL, new Scale(0.02, 0.02, 0.05, new Cylinder))));

    psInsulate = new Scale(PHULL + ITHICK, PHULL + ITHICK, -(PHULL + ITHICK), Knoboid2(0.666));
    psInsulate = new Difference (psInsulate, new Scale(PHULL, PHULL, PHULL, new Sphere));
    psInsulate = new Difference (psInsulate, new Translate(0.0, 0.0, 1.0 + 1.17 - PHIGH - PHULL, new Cube));
    psInsulate = PTKV_260(psInsulate);

    psFlanges = new BoundingBox(MakeUnion(
        Flange3(0.0),
        Flange3(60.0),
        Flange3(120.0),
        Flange3(180.0),
        Flange3(240.0),
        Flange3(300.0)
    ));
    return new BoundingBox(new Translate(0.0, 0.0, PHIGH + PHULL, 
                MakeUnion(
                    psHull,
                    psFlanges
                    // psInsulate
                )));
}

static Solid*
ShockAbsorber(const Point3 &p1, const Point3 &p2)
{
    Solid *ps;

    ps = MakeUnion(
        CylinderFromTo(p1, p2, RSTRUT),
        new Translate(p1, new Scale(RSTRUT, RSTRUT, RSTRUT, new Sphere)),
        new Translate(p2, new Scale(RSTRUT, RSTRUT, RSTRUT, new Sphere))
    );
    return ps;
}

static Solid*
Strut(double fAng)
{
    Solid *ps;
    Point3 p1, p2;
    Vector3 vRot1, vRot2, vRot1b, vRot2b;
    double fRad1, fRad2;

    p1 = Point3(1.0-RINGWIDE-0.5*RSTRUT, 0.0, RINGHIGH+1.0*RSTRUT);
    p2 = Point3(PHULL, 0.0, PHULL+PHIGH-FLANGEX);
    fRad1 = (30.0-1.5) * D_PI/180.0;
    fRad2 = 5.0 * D_PI/180.0;
    vRot1  = Point3(p1.x*cos(fRad1), p1.x*sin(fRad1), p1.z) - p1;
    vRot2  = Point3(p2.x*cos(fRad2), p2.x*sin(fRad2), p2.z) - p2;
    vRot1b = Point3(p1.x*cos(-fRad1), p1.x*sin(-fRad1), p1.z) - p1;
    vRot2b = Point3(p2.x*cos(-fRad2), p2.x*sin(-fRad2), p2.z) - p2;
    ps = MakeUnion(
        ShockAbsorber(p1, p2),
        ShockAbsorber(p1+vRot1,  p2+vRot2),
        ShockAbsorber(p1+vRot1b, p2+vRot2b)
    );
    ps = Titanium2(ps);
    return new BoundingBox(new Rotate(0.0, 0.0, fAng * D_PI/180.0 + 0.5*D_PI, ps)); 
}

static Solid*
Struts()
{
    return MakeUnion(
        Strut(0.0),
        Strut(60.0),
        Strut(120.0),
        Strut(180.0),
        Strut(240.0),
        Strut(300.0)
    );
}

#define CRAD    (10.5*0.5)
#define CQUARTZ (12*0.5)
#define CTRIM   ((12.0+1.0)*0.5*sqrt(2.0) - 0.7)
#define CDELTA  ((12.0+1.0)*0.5)
#define CLONG   16.0
#define CSTEM   (CLONG-7.0)
#define CSTDEL  (CLONG-6.0)/2.0

static Solid*
Camera(double fAng)
{
    Solid *ps, *ps2;

    ps = MakeUnion(
        new Translate(0.0, 0.0, 0.5, new Scale(CRAD, CRAD, 0.5, new Cylinder)),
        new Translate(0.0, 0.0, CSTDEL, new Scale(5.5/2.0, 5.5/2.0, CSTEM/2, new Cylinder)),
        new Rotate(0.0, 0.0, 0.0*D_PI/1.5, new Translate(5.5/2.0, 0.0, CSTDEL, new Scale(1.0, 0.25, CSTEM/2.0, new Cylinder))),
        new Rotate(0.0, 0.0, 1.0*D_PI/1.5, new Translate(5.5/2.0, 0.0, CSTDEL, new Scale(1.0, 0.25, CSTEM/2.0, new Cylinder))),
        new Rotate(0.0, 0.0, 2.0*D_PI/1.5, new Translate(5.5/2.0, 0.0, CSTDEL, new Scale(1.0, 0.25, CSTEM/2.0, new Cylinder))),
        new Translate(0.0, 0.0, CLONG+3.5/2.0, new Scale(CQUARTZ+0.5, CQUARTZ+0.5, 3.5/2.0, new Cylinder)),
        new Translate(0.0, 0.0, CLONG-3.5/2.0-6.0, new Scale(CQUARTZ+0.5, CQUARTZ+0.5, 3.5/2.0, new Cylinder))
    );
    ps = new Difference(ps, new Translate(0.0, 0.0, CSTDEL+1.0, new Scale(4.0/2.0, 4.0/2.0, CSTDEL+1.1, new Cylinder)));
    ps = new Difference(ps, new Translate(0.0, 0.0, CLONG-2.5, new Scale(4.9, 4.9, 4.9, new Sphere)));
    ps = new Difference(ps, new Translate(0.0, 0.0, CLONG-3.5, new Scale(4.9, 4.9, 4.9, new Sphere)));
    ps = new Intersection(ps, MakeUnion(
            new Translate(0.0, 0.0, CSTDEL, new Scale(CRAD+0.1, CRAD+0.1, CSTDEL+0.01, new Cylinder)),
            new Translate(0.0, 0.0, CLONG+3.5-CDELTA, new Scale(CTRIM, CTRIM, CTRIM, new Sphere)),
            new Translate(0.0, 0.0, CLONG-6.0-3.5+CDELTA, new Scale(CTRIM, CTRIM, CTRIM, new Sphere))
        ));
    ps = MakeUnion(
            new Rotate(0.0, 0.0, 0.0*D_PI/3.0, new Translate(0.0, 4.0, 0.0, new Scale(0.5, 0.5, 1.5, new Prism(6)))),
            new Rotate(0.0, 0.0, 1.0*D_PI/3.0, new Translate(0.0, 4.0, 0.0, new Scale(0.5, 0.5, 1.5, new Prism(6)))),
            new Rotate(0.0, 0.0, 2.0*D_PI/3.0, new Translate(0.0, 4.0, 0.0, new Scale(0.5, 0.5, 1.5, new Prism(6)))),
            new Rotate(0.0, 0.0, 3.0*D_PI/3.0, new Translate(0.0, 4.0, 0.0, new Scale(0.5, 0.5, 1.5, new Prism(6)))),
            new Rotate(0.0, 0.0, 4.0*D_PI/3.0, new Translate(0.0, 4.0, 0.0, new Scale(0.5, 0.5, 1.5, new Prism(6)))),
            new Rotate(0.0, 0.0, 5.0*D_PI/3.0, new Translate(0.0, 4.0, 0.0, new Scale(0.5, 0.5, 1.5, new Prism(6)))),
            ps
        );
    ps = new Scale(0.01, 0.01, 0.01, ps);       // modeled in centimeters
    ps = Titanium(ps);
    ps2 = (new Difference(new Scale(CQUARTZ, CQUARTZ, 6.5/2.0, new Cylinder), 
                                new Scale(CQUARTZ-1.0, CQUARTZ-1.0, 12.0/2.0+0.01, new Cylinder)));
    ps2 = new Translate(0.0, 0.0, CLONG-3.0, ps2);
    ps2 = Quartz(new Scale(0.01, 0.01, 0.01, ps2));     // modeled in centimeters
    ps = new Union(ps, ps2);
    ps = new Translate(0.0, 0.0 , PHULL+PHIGH, new Rotate(0.0, (90.0-50.0) *D_PI/180.0, 0.0, new Translate(0.0, 0.0, PHULL, ps)));
    ps = new Rotate(0.0, 0.0, fAng *D_PI/180.0, ps);
    return new BoundingBox(ps);
}

static Solid*
VeneraLander()
{
    return new Rotate(90.0 * D_PI/180.0, 0.0, 0.0, MakeUnion(
        new Rotate(+0.0 * D_PI/180.0, 0.0, 0.0, new Translate(0.0, 1.5, -0.2, VeneraColorChips())),
        LandingRing(),
        Camera(0.0),
        Camera(180.0),
        Struts(),
        Body()
    ));
}

static Solid*
TestTorus()
{
    return MakeUnion(
        Titanium(new Translate(-3.0, 0.0, 0.0, new Sphere)),
        Titanium(new Rotate(0.5*D_PI, 0.0, 0.0, new Scale(0.2, 0.2, 1.0, new Cylinder))),
        Titanium(new Translate(+3.0, 0.0, 0.0, new Torus(0.2)))
    );
}

void
Venera()
{
    RayTracer scene;
    Image im;
    Video vid;
    ML_TimingInfo ti;
    Spectrum s6500, sIllumE;
    Solid *ps;
    double t;

    ColorChips();
    return;

    s6500 = BlackBody(6500.0);
    sIllumE = g_sEqualWhite;
   // scene.plAmbient = new HavercosineLight(g_sVenera11Sky * 0.15 /70.0, g_sVenera14Ground * 0.15 / 70.0);
    scene.plAmbient = new AmbientLight(sIllumE * 0.02);
    scene.plPointLights = new PointLight(Point3(9.0, -6.5, 12.0), s6500 * 10000.0);
    //scene.plPointLights = new PointLight(Point3(3.0, -3.0, 100.0), s6500 * (10000.0 * 20.0));
    scene.plPointLights = new PointLight(Point3(-10.0, -10.0, -10.0), s6500 * 100.0, scene.plPointLights);

    scene.psSample = new FastBlueSamplingTile;
    //scene.psSample->nVenera = 1;
    scene.pcCamera = new StereographicCamera(Point3(0.0, -0.4, 9.0), 105.0); 

    im.NewImage(1920, 1080, 3);
    vid.NewVideo("Venera.mpg", im);
    ML_StartTiming(ti);
    for (t = 0.0; t < 360.0; t += 1.0) {
        ps = new Rotate(0.0, t*D_PI/180.0, 0.0, VeneraLander());
        //ps = Slice(ps);
        ps = EnvironmentLab(ps, scene);
        scene.psModel = ps;
        scene.psModel = scene.psModel->Optimize();
        scene.Render(im, 4.0);
        vid.WriteFrame(im);
        printf("t = %f\n", t);
      
    }
    ML_StopTiming(ti);
    vid.Close();
    im.WriteBMP("Venera.bmp");
    ML_ReportTiming(ti);
}