#include "stdafx.h"
#include "RayTracer2020.h"
#include "Solids.h"
//
//  Shading, surfaces, materials
//  - reflect bumpmaps when subtracting solids or internal rays
//

int
IsNiceNumber(double f)
{
    switch (_fpclass(f)) {
        case _FPCLASS_SNAN:
        case _FPCLASS_QNAN:
        case _FPCLASS_NINF:
        case _FPCLASS_PINF:
            return 0;
        default:
            return 1;
    }
}

static int
BadSpectrum(const Spectrum &s)
{
    int i;

    for (i = 0; i < NSPEC; i++)
        if (!IsNiceNumber(s.rgf[i]))
            return 1;
    return 0;
}

int
BadRay(Ray &ray)
{
    int i;

    for (i = 0; i < 3; i++)
        if (!IsNiceNumber(ray.pOrg.rgf[i]))
            return 1;
    for (i = 0; i < 3; i++)
        if (!IsNiceNumber(ray.vDir.rgf[i]))
            return 1;
    return 0;
}

void
Abort()
{
    static int *pi = 0, n;

    n = *pi;
}


#define BADPIXEL (g_nTestRow == 226 && g_nTestColumn == 418)

//
//  Calculate a unit surface normal, oriented toward the origin of the ray, not messed up by bump mapping
//
static inline void
ConsistentSurfaceNormal(Ray &r, const Hit *ph, Vector3 &vE, CoVector3 &cvN, double &fNE)
{
    CoVector3 cvN2;
    double fNE2;

    vE = -Normalize(r.vDir);                            // normalized vector to camera/eye.
    if (ph->nPerturbed) {                    // Unit normal vector with no bumpmap artifacts
        cvN = Normalize(ph->cvPerturbed);               //REVIEW: SurfaceTexture could be called here
        fNE = cvN * vE;   
        cvN2 = ph->cvNormal;
        fNE2 = cvN2 * vE; 
        if (fNE * fNE2 < 0.0) {                         // Perturbed normal inconsistent with geometry?
            cvN = Normalize(ph->cvNormal);
            fNE = cvN * vE;
        }
    } else {
        cvN = Normalize(ph->cvNormal);
        fNE = cvN * vE;
    }
    if (fNE < 0.0) {                                    // orient surface toward camera
        cvN = -cvN;
        fNE = -fNE;
    }
}
//
//  Fresnel law for reflection coefficient at a boundary between two refractive indices
//
static inline void
FresnelLaw_Real(double fNE, double fN1, double fN2, double &fFresnel, double &fNT)
{
    double fEta, Rs, Rp, fNEcritical;

    fEta = fN1/fN2;  
    fNEcritical = 1.0 - 1.0/(fEta*fEta);
    if (fNEcritical < 0.0)
        fNEcritical = 0.0;
    else
        fNEcritical = sqrt(fNEcritical);
    if (fNE < fNEcritical) {
        fFresnel = 1.0;                                     // total internal reflection
    } else {
        fNT = sqrt(1.0 - fEta*fEta*(1.0 - fNE*fNE));        // this may be used later to calculate refraction ray.vDir
        if (fN1 == fN2) {
            fFresnel = 0.0;
        } else {
            Rs = (fN1*fNE - fN2*fNT)/(fN1*fNE + fN2*fNT);
            Rp = (fN1*fNT - fN2*fNE)/(fN1*fNT + fN2*fNE);
            fFresnel = 0.5*(Rs*Rs + Rp*Rp);
        }
    }
}
//
//  Complex refractive indices for metals
//
static void
FresnelLaw_Complex(double fNE, Complex fN1, Complex fN2, double &fFresnel, double &fNT)
{
    Complex fEta, Rs, Rp;
    double fNEcritical;

    fEta = fN1/fN2;  
    fNEcritical = 1.0 - 1.0/fEta.AbsSqr();
    if (fNEcritical < 0.0)
        fNEcritical = 0.0;
    else
        fNEcritical = sqrt(fNEcritical);
    if (fNE < fNEcritical) {
        fFresnel = 1.0;                                     // total internal reflection
    } else {
        fNT = sqrt(1.0 - fEta.AbsSqr()*(1.0 - fNE*fNE)); 
        Rs = (fN1*fNE - fN2*fNT)/(fN1*fNE + fN2*fNT);
        Rp = (fN1*fNT - fN2*fNE)/(fN1*fNT + fN2*fNE);
        fFresnel = 0.5*(Rs.AbsSqr() + Rp.AbsSqr());
        if (fN1 == fN2)
            fFresnel = 0.0;
    }
}

static Spectrum
Fresnel_Spectral(double fNE, double fN1, const Spectrum &sRefraction, const Spectrum &sExtinction)
{
    Spectrum s;
    double fFresnel, fNT;
    int i;

    for (i = 0; i < NSPEC; i++) {
        FresnelLaw_Complex(fNE, Complex(fN1, 0.0), Complex(sRefraction.rgf[i], sExtinction.rgf[i]), fFresnel, fNT);
        s.rgf[i] = float(fFresnel);
    }
    return s;
}

static inline Reflectance
BidirectionalReflectance(const Hit *ph, const Vector3 &vL, const Vector3 &vE, const CoVector3 &cvN, double fFresnel)
{
    double fNH, fBlinn, fLambert, fGeometry, fBRDF, f1, f2, fEH, fNE, fNL, fNLSpec, fKS;
    int logKappa, kappa;
    Vector3 vH;
    CoVector3 cvNSpec;

    fLambert = 1.0 / D_PI;
    fKS = 1.0;
    cvNSpec = cvN;
    switch (ph->ps->nRoughness) { 
        case SURF_MATTE: 
            return ph->pm->color * fLambert;    //REVIEW: what to do with transparent material?
            logKappa = 0;
            break;
        case SURF_BURNISHED:
            logKappa = 3;
            break;
        case SURF_SMOOTH:
            logKappa = 5;
            fKS = 1.0;
            break;
        case SURF_WET:
            cvNSpec = Normalize(ph->cvNormal);  //  Plastic with clear smooth resin over bump-mapped pigment
        case SURF_SHINY:
            logKappa = 7;
            fKS = 0.60;
            break;
        case SURF_POLISHED:   
            logKappa = 9;
            fKS = 0.35;
            break;
    }
    if (logKappa) {
        //
        //  Blinn style microfacet distribution
        //
        fNL = cvN * vL;
        fNLSpec = cvNSpec * vL;
        fNE = cvNSpec * vE;
        vH = Normalize(vE + vL);
        fNH = cvNSpec * vH;
        kappa = 1;
        fBlinn = fNH;
        while (logKappa > 0) {
            fBlinn = fBlinn * fBlinn;
            kappa  =  kappa + kappa;
            logKappa--;
        };
        fBlinn = fKS*(double(kappa) + 2.0)*fBlinn / (2.0*D_PI);
        //
        //  Geometric microfacet occlusion factor
        //
        fEH = Dot(vE, vH);
        f1 = fabs(2.0*fNH*fNE/fEH);
        f2 = fabs(2.0*fNH*fNLSpec/fEH);
        fGeometry = Min(1.0, Min(f1, f2));
        fBRDF = fBlinn * fGeometry / (4.0*fNL*fNE);      
    } else
        fBRDF = 0.0;
    switch (ph->pm->nLuster) {
        case MAT_PLASTIC:
            return ph->pm->color * fLambert + g_sEqualWhite * (fFresnel*fBRDF);
        case MAT_METAL:
            return ph->pm->color*((1.0-fFresnel)*fBRDF) + g_sEqualWhite*(fFresnel*fBRDF);
        case MAT_TRANSPARENT:
            return g_sEqualWhite*(fFresnel*fBRDF);
        case MAT_COPPER:
        case MAT_GOLD:
        case MAT_IRON:
        case MAT_ALUMINUM:
        case MAT_OSMIUM:
        case MAT_SILVER:
        case MAT_TITANIUM:
            return ph->pm->color*fBRDF;                     // spectral Fresnel law calculated and written into pm->color

    }
    return g_sBlack;
}

static inline Reflectance
RayReflectance(Hit *ph, double fFresnel)
{
    switch (ph->pm->nLuster) {
        case MAT_PLASTIC:
            return g_sEqualWhite * fFresnel;
        case MAT_METAL:
            return ph->pm->color*(1.0-fFresnel) + g_sEqualWhite*fFresnel;
        case MAT_TRANSPARENT:
            return g_sEqualWhite*fFresnel;
        case MAT_COPPER:
        case MAT_GOLD:
        case MAT_IRON:
        case MAT_ALUMINUM:
        case MAT_OSMIUM:
        case MAT_SILVER:
        case MAT_TITANIUM:
            return ph->pm->color;
    }
    return g_sBlack;
}
//
//  Late evalutation of the surface-normal affine transformations, so work only done on visible surfaces
//  Apply surface and material texture call-outs at intermediate coordinate systems
//
//  This is where instancing needs multiple parent pointers (FX)
//  FX built an ipList during descent through instance nodes.  During ascent:  parents, ip++, parents, ip++, ...
//
static inline void
TransformNormals(Hit *ph, RayTracer *scene)
{
    AffineSolid *pa;
    TextureList *pt;

    pt = ph->pt;
    for (pa = ph->paAffineTransforms; pa; pa = pa->paParent) {
        while (pt && pt->pa == pa) {
            (*pt->Texture)(pt->p, ph, scene);
            pt = pt->ptNext;
        }
        ph->cvNormal = ph->cvNormal * pa->mInverse;
        if (ph->nPerturbed)
            ph->cvPerturbed = ph->cvPerturbed * pa->mInverse;
    }
    while (pt && pt->pa == pa) {                                            // don't forget the pa == 0 case
        (*pt->Texture)(pt->p, ph, scene);
        pt = pt->ptNext;
    }
}
//
//  Process a hit list to calculate radiance from a visual ray or irradiance from a shadow ray.
//
//         pm         ps   pmNext
//  O------------------|-----------|-----|------------------>
// ph             ph->phNext
//
static void
Dummy()
{
    static int x = 0;

    VERBOSE printf("Dummy\n");
    x += 1;
}

Irradiance  
RayTracer::LightRay(Ray &ray, const PointLight *pl)
{
    Hit *ph, *phAllocated, hVacuum(&g_mVacuum);
    double t, r2, fNE, fFresnel, fNT;                   // fNT is used in VisualRay but not in LightRay
    Vector3 vE;
    CoVector3 cvN;
    Material *pm;
    Transmittance tr;
    Irradiance E;
    Propagation prop;

    VERBOSE printf("LightRay ");
    VERBOSE ray.Print();
    r2 = Dot(ray.vDir, ray.vDir);                       
    E = pl->flux / (4.0*D_PI * r2);                     // unimpeded irradicance from the light
    if (pl->nNoShadow)
        return E;
    ph = phAllocated = psModel->Intersect(ray);
    if (ph == 0 || ph->t > 0.0) {                       // regularize hit list into ( 0, t, ... )
        hVacuum.phNext = ph;
        ph = &hVacuum;
    }
    VERBOSE PrintHitList("  Light Intersect:", ph);
    t = 0.0;
    pm = ph->pm;
    tr = g_sEqualWhite;
    while ((ph = ph->phNext) && (ph->t < 1.0-EPSILON)) {
        VERBOSE printf("pm %d ps %d pmNext %d\n", pm, ph->ps, ph->pm);
        TransformNormals(ph, this);
        if (pm->nLuster != MAT_TRANSPARENT || ph->pm->nLuster != MAT_TRANSPARENT || ph->ps->nRoughness != SURF_POLISHED)
            goto Shadow;
        if (pm->Propagate) {                                                                // propagate light thru volume
            prop = pm->Propagate(ray(t), ray(ph->t), pm, this);
            tr = tr*prop.cTransmittance;                                                    // Transmittance accumulates
        }
        if (pm->fRefractiveIndex != ph->pm->fRefractiveIndex) {                             // surface boundary
            ConsistentSurfaceNormal(ray, ph, vE, cvN, fNE);
            FresnelLaw_Real(fNE, ph->pm->fRefractiveIndex, pm->fRefractiveIndex, fFresnel, fNT);
                if (fFresnel == 1.0)
                    goto Shadow;                                // light blocked by total internal reflection
            tr = tr * (1.0 - fFresnel);
        }
        t = ph->t;
        pm = ph->pm;
    }
    if (pm->Propagate && (pm->nLuster == MAT_TRANSPARENT) ) {
        VERBOSE printf("calling propagator\n");
        prop = pm->Propagate(ray(t), ray(1.0), pm, this);      // propgation in material before light
        // VERBOSE printf("returned from propagator\n");
        Dummy();
        tr = tr * prop.cTransmittance;
    }
    DeleteHitList(phAllocated);
    return E * tr;
Shadow:
    DeleteHitList(phAllocated);
    return g_sBlack;
}

static inline Vector3
ReflectionDirection(const CoVector3 &cvN, const Vector3 &vE, const Hit *ph, double fNE)
{
    Vector3 vR;
    CoVector3 cvN2;
    double fNE2;

    vR = Vector3(cvN.x, cvN.y, cvN.z) * 2.0*fNE - vE;
    if (ph->nPerturbed && ph->cvNormal * vR < EPSILON) {                // don't let bumpmap reflect into surface (necessary?)
        cvN2 = Normalize(ph->cvNormal);
        fNE2 = cvN2 * vE;
        vR = Vector3(cvN2.x, cvN2.y, cvN2.z) * 2.0*fNE2 - vE;   // almost tangential to surface
    }
    return vR;
}

static inline Vector3
RefractionDirection(const CoVector3 &cvN, const Vector3 &vE, const Hit *ph, double fEta, double fNE, double fNT)
{
    Vector3 vT;
    CoVector3 cvN2;
    double fNE2;

    vT = Vector3(cvN.x, cvN.y, cvN.z) * (fEta*fNE - fNT) - vE*fEta;
    if (ph->nPerturbed && ph->cvNormal * vT < EPSILON) {                // don't let bumpmap reflect into surface (necessary?)
        cvN2 = Normalize(ph->cvNormal);
        fNE2 = cvN2 * vE;
        vT = Vector3(cvN2.x, cvN2.y, cvN2.z) * (fEta*fNE2 - fNT) - vE*fEta;
    }
    return vT;
}
//
//  Cast a visual ray into the scane to sample radiance (image brightness)
//
static inline void
Reset_tMin(Ray &ray)
{
    ray.tMin = INFINITY;            // Be careful to do this if a ray is being reused
}

Radiance
RayTracer::VisualRay(Ray &ray)
{
    Material *pm;
    Hit *ph, *phAllocated, hVacuum(&g_mVacuum);
    PointLight *pl;
    double t, fNE, fNL, fFresnel, fEta, fNT;
    Propagation prop;
    Transmittance tr;
    Reflectance refl;
    Radiance R;
    Irradiance E;
    Vector3 vE, vL;
    CoVector3 cvN;
    Ray rSecondary;
    static int nRayDepth = 0;

    VERBOSE printf("VisualRay ");
    VERBOSE ray.Print();
    if (nRayDepth >= 11)
        return g_sBlack;
    ph = phAllocated = psModel->Intersect(ray);
    nRayDepth++;
    if (ph == 0 || ph->t > 0.0) {                       // regularize hit list into ( 0, t, ... )
        hVacuum.phNext = ph;
        ph = &hVacuum;
    }
    VERBOSE ph->Print();
    t = 0;
    pm = ph->pm;
    tr = g_sEqualWhite;
    R = g_sBlack;
    while (ph = ph->phNext) {
        VERBOSE printf("pm %d ps %d pmNext %d\n", pm, ph->ps, ph->pm);
        TransformNormals(ph, this);
        if (pm->nLuster == MAT_NOSHADE)                 // sky domes and such
            return pm->color;
        if (pm->Propagate) {                                                                // propagate light thru volume
            prop = pm->Propagate(ray(t), ray(ph->t), pm, this);
            tr = tr*prop.cTransmittance;                                                    // Transmittance accumulates
            R = R + prop.cSterisent;                                                        // Sterisent just contributes once
        }
        if ((pm->fRefractiveIndex != ph->pm->fRefractiveIndex) || (ph->pm->nLuster != MAT_TRANSPARENT)) {                             // surface boundary
            ConsistentSurfaceNormal(ray, ph, vE, cvN, fNE);
            switch (ph->pm->nLuster) {
                case MAT_COPPER:
                   ph->pm->color = Fresnel_Spectral(fNE, pm->fRefractiveIndex, g_sCopperRefraction, g_sCopperExtinction);
                    break;
                case MAT_GOLD:
                    ph->pm->color = Fresnel_Spectral(fNE, pm->fRefractiveIndex, g_sGoldRefraction, g_sGoldExtinction);
                    break;
                case MAT_IRON:
                    ph->pm->color = Fresnel_Spectral(fNE, pm->fRefractiveIndex, g_sIronRefraction, g_sIronExtinction);
                    break;
                case MAT_ALUMINUM:
                    ph->pm->color = Fresnel_Spectral(fNE, pm->fRefractiveIndex, g_sAluminumRefraction, g_sAluminumExtinction);
                    break;
                case MAT_OSMIUM:
                    ph->pm->color = Fresnel_Spectral(fNE, pm->fRefractiveIndex, g_sOsmiumRefraction, g_sOsmiumExtinction);
                    break;
                case MAT_SILVER:
                    ph->pm->color = Fresnel_Spectral(fNE, pm->fRefractiveIndex, g_sSilverRefraction, g_sSilverExtinction);
                    break;
                case MAT_TITANIUM:
                    ph->pm->color = Fresnel_Spectral(fNE, pm->fRefractiveIndex, g_sTitaniumRefraction, g_sTitaniumExtinction);
                    break;
                default:
                    FresnelLaw_Real(fNE, pm->fRefractiveIndex, ph->pm->fRefractiveIndex, fFresnel, fNT);
                    break;
            }
            if (ph->ps == 0) ph->Print();
            if (plAmbient && (ph->ps->nRoughness != SURF_POLISHED)) {
                E = plAmbient->Illuminate(rSecondary.pOrg, cvN, ph->ps->nRoughness);                   // ambient illumination if not reflective
                R = R + ph->pm->color * E;
            }
            rSecondary.pOrg = ray(ph->t);                                                   // Shadows, Reflections, Refractions start here
            for (pl = plPointLights; pl; pl = pl->plNext) {                                 // local shading per light
                rSecondary.vDir = (pl->Origin() - rSecondary.pOrg);
                Reset_tMin(rSecondary);
                if (cvN*rSecondary.vDir <= 0.0 && ph->pm->nLuster != MAT_TRANSPARENT)
                    continue;
                vL = Normalize(rSecondary.vDir);
                fNL = cvN * vL;
                E = LightRay(rSecondary, pl) * fNL;
                if (E == g_sBlack)
                    continue;
                refl = BidirectionalReflectance(ph, vL, vE, cvN, fFresnel);
                R = R + refl*E;
            }
            if (ph->ps->nRoughness == SURF_POLISHED) {                                      // global shading, reflection/refraction
                if (fFresnel && nRayDepth < 6) {
                    rSecondary.vDir = ReflectionDirection(cvN, vE, ph, fNE);
                    Reset_tMin(rSecondary);
                    R = R + VisualRay(rSecondary)*RayReflectance(ph, fFresnel);
                }
                if (ph->pm->nLuster == MAT_TRANSPARENT && fFresnel < 1.0) {
                    fEta = pm->fRefractiveIndex / ph->pm->fRefractiveIndex;
                    rSecondary.vDir = RefractionDirection(cvN, vE, ph, fEta, fNE, fNT);
                    Reset_tMin(rSecondary);
                    R = R + VisualRay(rSecondary)*(1.0 - fFresnel);
                }
            }
            R = R * tr;
            t = ph->t;
            pm = ph->pm;
            break;      // Visual rays stop a boundary and may split into reflect/refract rays
        } else {
            t = ph->t;
            pm = ph->pm;
        }
    }
    if (pm->Propagate) {
        prop = pm->Propagate(ray(t), ray(INFINITY), pm, this);
        R = R + tr*prop.cSterisent;
        tr = tr * prop.cTransmittance;
    }
    DeleteHitList(phAllocated);
    --nRayDepth;
    return R; 
}
