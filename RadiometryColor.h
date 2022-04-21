//
//  Colormetric and radiometric types, image sample values (spectral radiance)
//  D.P. Mitchell  2020/12/19
//
//  separate surface reflectance image and illumination image.  smoothing soft shadows.
//
#pragma once
#include "Image.h"

#define NDELTASPEC  40                          // this must be divisor of (740-380) and divisible by 5 (10,15,20,30,40,45...)
#define NSPEC       ((740-380)/NDELTASPEC + 1)  // 10
#define SCALESPEC   (5.0/10.0)

extern DisplayRGB SpectrumToRGB[NSPEC];

struct Spectrum {
    float  rgf[NSPEC];      // 380, 420, ... 740 nm

    Spectrum    operator+(Spectrum s) const {   Spectrum sNew; int i;
                                            for (i = 0; i < NSPEC; i++)
                                                sNew.rgf[i] = rgf[i] + s.rgf[i];
                                            return sNew;  
                                        }
    Spectrum    operator-(Spectrum s) const {   Spectrum sNew; int i;
                                            for (i = 0; i < NSPEC; i++)
                                                sNew.rgf[i] = rgf[i] - s.rgf[i];
                                            return sNew;  
                                        }
    Spectrum    operator*(Spectrum s) const {   Spectrum sNew; int i;
                                            for (i = 0; i < NSPEC; i++)
                                                sNew.rgf[i] = rgf[i] * s.rgf[i];
                                            return sNew;  
                                        }
    Spectrum    operator/(Spectrum s) const {   Spectrum sNew; int i;
                                            for (i = 0; i < NSPEC; i++) {
                                                if (s.rgf[i])
                                                    sNew.rgf[i] = rgf[i] / s.rgf[i];
                                                else
                                                    sNew.rgf[i] = 0.0;
                                            }
                                            return sNew;  
                                        }
    Spectrum    operator*(double f) const {   Spectrum sNew; int i;
                                            for (i = 0; i < NSPEC; i++)
                                                sNew.rgf[i] = float(rgf[i] * f);
                                            return sNew;  
                                        }
    Spectrum    operator*=(double f)    {   int i;
                                            for (i = 0; i < NSPEC; i++)
                                                rgf[i] *= float(f);
                                            return *this;
                                        }
    Spectrum    operator/(double f) const {   Spectrum sNew; int i;
                                            if (f) f = 1.0/f;
                                            for (i = 0; i < NSPEC; i++)
                                                sNew.rgf[i] = float(rgf[i] * f);
                                            return sNew;  
                                        }
    int         operator==(Spectrum const &s) const {   int i; for (i = 0; i < NSPEC; i++) if (rgf[i] != s.rgf[i]) return 0; return 1;
                                        }
    float&      operator[](int i)  { return rgf[i]; }
    DisplayRGB  sRGB() const {   DisplayRGB rgb(0.0, 0.0, 0.0); int i;
                                for (i = 0; i < NSPEC; i++)
                                    rgb = rgb + SpectrumToRGB[i] * (rgf[i] * SCALESPEC);
                                return rgb; 
                            }
    DisplayRGB  vRGB() const;
    void Print() const { int i; printf("{\n"); for (i = 0; i < NSPEC; i++) printf("    %0.8ff,\n", rgf[i]); printf("};\n"); }
};

extern Spectrum g_sX;
extern Spectrum g_sY;
extern Spectrum g_sZ;

inline double
Luminance(const Spectrum &s)
{
    int i;
    double f = 0.0;

    for (i = 0; i < NSPEC; i++)
        f += s.rgf[i] * g_sY.rgf[i];
    return f;
}

extern Spectrum RGBtoSpectrum(DisplayRGB rgb, int n6500 = 0);
extern Spectrum RGBtoSpectrum(double f);
extern Spectrum BlackBody(double fKelvin);
extern Spectrum g_sRed, g_sGreen, g_sBlue, g_sYellow, g_sCyan, g_sMagenta, g_sWhite, g_sOrange;
extern Spectrum g_sBlack, g_sEqualWhite, g_sD6500;
extern Spectrum g_sWater, g_sCobaltGlass, g_sCopper, g_sGold, g_sRuby, g_sEmerald, g_sSapphire, g_sWhiskey;
extern Spectrum g_sCopperExtinction, g_sCopperRefraction, g_sGoldExtinction, g_sGoldRefraction, g_sIronExtinction, g_sIronRefraction;
extern Spectrum g_sAluminumExtinction, g_sAluminumRefraction, g_sOsmiumExtinction, g_sOsmiumRefraction;
extern Spectrum g_sSilverExtinction, g_sSilverRefraction, g_sTitaniumExtinction, g_sTitaniumRefraction;

extern Spectrum g_sVenera14Ground;
extern Spectrum g_sVenera14Sky;
extern Spectrum g_sVenera13Sky;
extern Spectrum g_sVenera11Sky;
extern Spectrum g_sVenera11_07200m;
extern Spectrum g_sVenera11_16200m;
extern Spectrum g_sVenera11_23700m;
extern Spectrum g_sVenera11_37700m;
extern Spectrum g_sVenera11_48600m;
extern Spectrum g_sVenera11_51100m;
extern Spectrum g_sVenera11_56100m;
extern Spectrum g_sVenera11_62100m;
extern Spectrum g_sVenera11_Sun;

extern Spectrum g_sRedChannel;
extern Spectrum g_sGrnChannel;
extern Spectrum g_sBluChannel;
extern Spectrum g_sClrChannel;
extern Spectrum g_sCoolGray;
extern Spectrum g_sCoolRed;
extern Spectrum g_sCoolGrn;
extern Spectrum g_sCoolBlu;
extern Spectrum g_sHotGray;
extern Spectrum g_sHotRed;
extern Spectrum g_sHotGrn;
extern Spectrum g_sHotBlu;

extern Spectrum g_sIllumD65;
extern Spectrum g_sSolar;
extern Spectrum g_sEarthSky;
extern Spectrum g_sMarsAlbedo;
extern Spectrum g_sVenusAlbedo;
extern Spectrum g_sVenusUVcontrast;
extern Spectrum g_sMercuryAlbedo;

extern DisplayRGB g_rgbIron, g_rgbCopper, g_rgbGold, g_rgbZinc, g_rgbZilver, g_rgbAluminum, g_rgbSilicon;

//
//  Radiometry quantities
//

//  typedef double Color;           // Simple color-less flux for testing.  Later will use Spectrum
typedef Spectrum Color;             // Full multi-spectral color
typedef Color Flux;                 // watts
typedef Color Irradiance;           // watts per meter**2
typedef Color Emittance;            // watts per meter**2
typedef Color Radiance;             // watts per meter**2 per steradian
typedef Color Intensity;            // watts per steradian
typedef Color Reflectance;          // unitless ratio
typedef Color Transmittance;
typedef Color Absorptance;