//
//  Solid geometry and ray intersection
//  D.P. Mitchell 2020/12/23.
//
// subclasses of Solid:
//
// Surface - roughness, bump maps, etc.
// Material - color, propagator
// NullSolid
// Quadric
// ConvexPolyhedron
// Torus
// AffineSolid
// Union
// Intersection
// Difference
// BoundingBox
//
#pragma once
#define MAXHITS         16

#define TETRADIHEDRAL   70.528779365509308630754000660038       // dihedral angle acos(1/3)
#define TETRAINSCRIBE   0.20412414523193150818310700622549      // inscribed radius/edge
#define TETRACIRCUM     0.61237243569579452454932101867647      // circumscribed radius/edge  (3x inscribed)

#define CUBEDIHEDRAL    90.0
#define CUBEINSCRIBE    0.5
#define CUBECIRCUM      0.86602540378443864676372317075294   

#define OCTADIHEDRAL    109.47122063449069136924599933996       // acos(-1/3)
#define OCTAINSCRIBE    0.40824829046386301636621401245098
#define OCTACIRCUM      0.70710678118654752440084436210485

#define DODECDIHEDRAL   116.56505117707798935157219372045       // dihedral angle acos(-sqrt(5)/5)
#define DODECINSCRIBE   0.95105651629515357211643933337938      // insribed radius / edge
#define DODECCIRCUM     1.4012585384440735446766779353221       // circumsribed radius / edge

#define ICOSADIHEDRAL   138.18968510422140193414208326942       // dihedral angle acos(-sqrt(5)/3)
#define ICOSAINSCRIBE   0.755761314076170730480133702025        // inscribed sphere radius / edge length
#define ICOSACIRCUM     0.95105651629515357211643933337938      // circumscribed radius / edge length

//
//  Primitive Solids
//
//  HitWork is a workspace for primitive intersection numerical calculations.
//
struct HitWork {
    double      rgfRootsStorage[MAXHITS+1];
    CoVector3   rgvNormalsStorage[MAXHITS+1];
    double      *rgfRoots;                      // This makes room for an rgfRoots[-1]
    double      fEPSILON;
    CoVector3   *rgvNormals;
    short       nHits;
    char        nRayInsideAtInfinity;

                HitWork() : fEPSILON(EPSILON) { rgfRoots = rgfRootsStorage + 1; rgvNormals = rgvNormalsStorage + 1; } 
    Hit         *HitList();                 // build a hit list
};

struct NullSolid : Solid {
            NullSolid() {}
    Hit*    Intersect(Ray &ray) { return 0; }
};

struct Quadric : Solid {
            Quadric(double p, double q, double r, int clip) :
                P(p), Q(q), R(r), nClip(clip) {}
    Hit*    Intersect(Ray &ray);
    Extent  GetExtent() { return Extent(1.0); }

    double  P;          // surface: x**2 + y**2 + P*z**2 + Q*z + R = 0
    double  Q;
    double  R;
    int     nClip;      // intersect with flat faces (clip == 1) or negate solid (clip == -1)
};

struct Sphere : Quadric {
    Sphere() : Quadric(1.0, 0.0, -1.0, 0) {}
};

struct HollowSphere : Quadric {
    HollowSphere() : Quadric(1.0, 0.0, -1.0, -1) {}
};

struct Cylinder : Quadric {
    Cylinder() : Quadric(0.0, 0.0, -1.0, 1) {}
};

struct Cone : Quadric {
    Cone() : Quadric(-0.25, -0.5, -0.25, 0) {}
};

struct ConeFX : Quadric {                           // end caps radius 0.5 at z = +/- 1.0
    ConeFX() : Quadric(-0.25, 0.0, 0.0, 0) {}
};

struct ConeUnit : Quadric {                           // end caps radius 1.0 at z = +/- 1.0
    ConeUnit() : Quadric(-1.0, 0.0, 0.0, 0) {}
};

struct Hyperboloid1 : Quadric {
    Hyperboloid1(double R) : Quadric(-1.0, 0.0, -R*R, 0) {}   // one-sheet  Quadric(-P, 0.0, -1.0, 0)
};

struct Hyperboloid2 : Quadric {
    Hyperboloid2(double R) : Quadric(-1.0, 0.0, +R*R, 0) {}   // two-sheet  Quadric(-P, 0.0, +1.0, 0)
};

struct Paraboloid : Quadric {
    Paraboloid() : Quadric(0.0, -0.5, 0.0, 1) {}
};

inline double SSS(double x) { if (x < 0.0) return -x*x; else return x*x; }
inline double USS(double x) { return x*x; }

struct ConeHyperCylinder : Quadric { 
    double R;
                                           // neck radius r, end caps radius 1.0 at z = +/- 1.0
        ConeHyperCylinder(double r) : Quadric(SSS(r) - 1.0, 0.0, -SSS(r), 0), R(r) {}       // r = 0 to 1, cone to hyperboloid I to cylinder
    Extent GetExtent();
};

extern CoVector3 rgvTetrahedron[4];
extern CoVector3 rgvCube[8];
extern CoVector3 rgvOctahedron[6];
extern CoVector3 rgvDodecahedron[20];
extern CoVector3 rgvIcosahedron[12];
extern Extent ClipPlanes(const Extent &eInitial, const CoVector3 rgcv[], int nv, double fFromOrigin);

struct ConvexPolyhedron : Solid {
    CoVector3   *rgv;
    int         nSides;

                ConvexPolyhedron() {}
                ConvexPolyhedron(CoVector3 *pcv, int n) : rgv(pcv), nSides(n) {}
    Hit*        Intersect(Ray &ray);
    Extent      GetExtent();
};

struct Slab : Solid {                       // infinite slab from z = -1 to 1, for intersecting other solids
            Slab() {}
    Hit*    Intersect(Ray &ray);
    Extent  GetExtent() { return Extent(Point3(-INFINITY, -INFINITY, -1.0), Point3(INFINITY, INFINITY, 1.0)); }
};

struct Tetrahedron : ConvexPolyhedron {
            Tetrahedron() : ConvexPolyhedron(rgvTetrahedron, 4) {}
};

struct Cube : ConvexPolyhedron {
            Cube() : ConvexPolyhedron(rgvOctahedron, 6) {}
};

struct Octahedron : ConvexPolyhedron {
            Octahedron() : ConvexPolyhedron(rgvCube, 8) {}
};

struct Dodecahedron : ConvexPolyhedron {
            Dodecahedron() : ConvexPolyhedron(rgvIcosahedron, 12) {}
};

struct Icosahedron : ConvexPolyhedron {
            Icosahedron() : ConvexPolyhedron(rgvDodecahedron, 20) {}
};

struct Prism : ConvexPolyhedron {
            Prism(int n);
};

struct Fibonacci : ConvexPolyhedron {
            Fibonacci(int n);
};

//
//  Quartic implicit solid, Torus
//
struct Torus : Solid {
            Torus(double rMinor) : fMinorRadius(rMinor) {}   // major radius is 1.0
    Hit*    Intersect(Ray &ray);
    Extent  GetExtent();

    double  fMinorRadius;
};
//
//  Surface and Material defintions
//
#define SURF_MATTE      1
#define SURF_BURNISHED  2
#define SURF_SMOOTH     3
#define SURF_SHINY      4
#define SURF_POLISHED   5
#define SURF_WET        6

struct AffineSolid;
struct Surface;
struct Material;
struct Instance;

typedef void (*TextureMap)(const Point3 &p, Hit *ph, RayTracer *scene);

struct Surface : Solid {
    TextureMap  Texture;
    AffineSolid *paAffineTransforms;        // Old parent link system
    Solid       *psChild;
    int         nRoughness;

            Surface(int nRough, TextureMap pText, Solid *ps) : Texture(pText), paAffineTransforms(0), psChild(ps), nRoughness(nRough) {}
            Surface(int nRough, Solid *ps) : Texture(0), paAffineTransforms(0), psChild(ps), nRoughness(nRough) {}
    Hit    *Intersect(Ray &r);
    Extent  GetExtent() { return psChild->GetExtent(); }
    Solid  *Optimize();
};

#define MAT_NOSHADE     0
#define MAT_PLASTIC     1
#define MAT_METAL       2
#define MAT_TRANSPARENT 3
#define MAT_COPPER      4
#define MAT_GOLD        5
#define MAT_IRON        6
#define MAT_ALUMINUM    7
#define MAT_OSMIUM      8
#define MAT_SILVER      9
#define MAT_TITANIUM   10

#define DENSE_AIR       0.000293
#define DENSE_VENUS     0.014716
#define DENSE_WATER     0.333
#define DENSE_ICE       0.36
#define DENSE_PYREX     0.474
#define DENSE_GLASS     0.52
#define DENSE_LEADGLASS 0.71
#define DENSE_SAPPHIRE  0.77
#define DENSE_DIAMOND   1.42

#define REFRACT_AIR       1.000293
#define REFRACT_VENUS     1.014716
#define REFRACT_WATER     1.333
#define REFRACT_ETHANOL   1.3571
#define REFRACT_ICE       1.36
#define REFRACT_QUARTZ    1.4298
#define REFRACT_CELLULOSE 1.4637
#define REFRACT_PYREX     1.474
#define REFRACT_CORNING   1.5078
#define REFRACT_LEADGLASS 1.71
#define REFRACT_SAPPHIRE  1.77
#define REFRACT_DIAMOND   2.42

#define REFRACT_COPPER    Complex(1.0697, 2.5866)
#define REFRACT_SILVER    Complex(0.056895, 3.5047)
#define REFRACT_ALUMINUM  Complex(1.7795, 17.029))
#define REFRACT_BRASS     Complex(0.32186, 3.3439)
#define REFRACT_GOLD      Complex(0.48899, 2.3389)
#define REFRACT_IRON      Complex(3.1228, 5.5463)
#define REFRACT_MERCURY   Complex(6.4621, 8.5349)
#define REFRACT_LEAD      Complex(1.8627, 9.3253)
#define REFRACT_SILICON   Complex(3.6730, 0.0050000)
#define REFRACT_ZINC      Complex(1.8491, 7.3489)

struct Propagation {
    Transmittance   cTransmittance;
    Radiance        cSterisent;
};

struct Material;

typedef Propagation (*Propagator)(const Point3 &p0, const Point3 &p1, Material *pm, RayTracer *scene);

struct Material : Solid {
    double      fDensity;                       // used in lattice-algebra set operations
    double      fRefractiveIndex;               
    Color       color;
    TextureMap  Texture;
    Propagator  Propagate;
    AffineSolid *paAffineTransforms;
    Solid       *psChild;
    int         nLuster;                

            Material(double f, Color &c, int nL, Solid *ps) :
                fDensity(f), fRefractiveIndex(1.0+f), color(c), 
                Texture(0), Propagate(0), paAffineTransforms(0), psChild(ps), nLuster(nL) {}
            Material(double f, double fN, Color &c, int nL,  Propagator pProp, Solid *ps) :
                fDensity(f), fRefractiveIndex(fN), color(c), 
                Texture(0), Propagate(pProp), paAffineTransforms(0), psChild(ps), nLuster(nL) {}
            Material(double f, Color &c, int nL, TextureMap pText, Solid *ps) : 
                fDensity(f), fRefractiveIndex(1.0+f), color(c), 
                Texture(pText), Propagate(0), paAffineTransforms(0), psChild(ps), nLuster(nL) {}
    Hit    *Intersect(Ray &r);
    Extent  GetExtent() { return psChild->GetExtent(); }
    Solid  *Optimize();
};

extern Material g_mVacuum, g_mNeutronium;

struct TextureList {
    Point3          p;
    TextureMap      Texture;
    AffineSolid     *pa;
    TextureList    *ptNext;
    static TextureList *ptFree;

            TextureList(TextureList *ptN) : ptNext(ptN) {}
};

extern void PrintTextureList(char *sz, TextureList *pt);
//
//  Affine transformations of solids
//
struct AffineSolid : Solid {
    Matrix3     mTransform;
    Matrix3     mInverse;
    Vector3     vTranslate;
    Solid       *psChild;
    AffineSolid *paParent;

            AffineSolid(Solid *ps) : psChild(ps), paParent(0) { Solid::nIsAffine = 1; }
    Hit*    Intersect(Ray &ray);
    Extent  GetExtent();
    Solid   *Optimize();
};

struct Scale : AffineSolid {
    Scale(double x, double y, double z, Solid *ps) : AffineSolid(ps) {
        if (x == 0.0 || y == 0.0 || z == 0.0)
            printf("Afine Scale singular %g %g %g\n", x, y, z);
        mTransform = ScaleMatrix(x, y, z);
        mInverse   = ScaleMatrix(1.0/x, 1.0/y, 1.0/z);
        vTranslate = Vector3(0.0, 0.0, 0.0);
    }
};

struct Rotate : AffineSolid {
public:
    Rotate(double x, double y, double z, Solid *ps, char *sz = "xyz") : AffineSolid(ps) {
        mTransform = EulerMatrix(x, y, z, sz);
        mInverse   = mTransform.Transpose();    //  InverseMatrix(mTransform);
        vTranslate = Vector3(0.0, 0.0, 0.0);
    }
    Rotate(const Vector3 &vAxis, double fRad, Solid *ps) : AffineSolid(ps) {
        mTransform = RotateAboutAxis(vAxis, fRad);
        mInverse   = mTransform.Transpose();    //  InverseMatrix(mTransform);
        vTranslate = Vector3(0.0, 0.0, 0.0);
    }
    Rotate(const Vector3 &vFrom, const Vector3 &vTo, Solid *ps) : AffineSolid(ps) {
        mTransform = RotateFromTo(vFrom, vTo);
        mInverse   = mTransform.Transpose();    //  InverseMatrix(mTransform);
        vTranslate = Vector3(0.0, 0.0, 0.0);
    }
};

struct Translate : AffineSolid {
    Translate(double x, double y, double z, Solid *ps) : AffineSolid(ps) {
        mTransform = Matrix3(1.0);
        mInverse   = Matrix3(1.0);
        vTranslate = Vector3(x, y, z);
    }
    Translate(const Vector3 &v, Solid *ps) : AffineSolid(ps) {
        mTransform = Matrix3(1.0);
        mInverse   = Matrix3(1.0);
        vTranslate = v; 
    }
    Translate(const Point3 &p, Solid *ps) : AffineSolid(ps) {
        mTransform = Matrix3(1.0);
        mInverse   = Matrix3(1.0);
        vTranslate = p - Point3(0.0, 0.0, 0.0);
    }
};
//
//  Set operations on solids, based on a lattice algebra of material density
//
struct Union : Solid {
    Solid   *psLeft;
    Solid   *psRight;

            Union(Solid *ps1, Solid *ps2) : psLeft(ps1), psRight(ps2) {}
    Hit*    Intersect(Ray &ray);
    Extent  GetExtent();
    Solid   *Optimize() { if(psLeft) psLeft = psLeft->Optimize(); if(psRight) psRight = psRight->Optimize();  return this; }
};

struct TopUnion : Union {

            TopUnion(Solid *ps1, Solid *ps2) : Union(ps1, ps2) {}
    Hit*    Intersect(Ray &ray);
    Extent  GetExtent();
    Solid   *Optimize() { if(psLeft) psLeft = psLeft->Optimize(); if(psRight) psRight = psRight->Optimize();  return this; }
};

struct Intersection : Solid {
    Solid   *psLeft;
    Solid   *psRight;

            Intersection(Solid *ps1, Solid *ps2) : psLeft(ps1), psRight(ps2) {}
    Hit*    Intersect(Ray &ray);
    Extent  GetExtent();
    Solid   *Optimize() { if(psLeft) psLeft = psLeft->Optimize(); if(psRight) psRight = psRight->Optimize();  return this; }
};

struct Difference : Solid {
    Solid   *psLeft;
    Solid   *psRight;

            Difference(Solid *ps1, Solid *ps2) : psLeft(ps1), psRight(ps2) {}
    Hit*    Intersect(Ray &ray);
    Extent  GetExtent();
    Solid   *Optimize() { if(psLeft) psLeft = psLeft->Optimize(); if(psRight) psRight = psRight->Optimize();  return this; }
};

struct AngularClip : Solid {
    Solid       *ps;
    CoVector3   cv1, cv2;
    double      fAng1, fAng2;
    int         nIntersect;

            AngularClip(double f1, double f2, Solid *p);
    Hit*    Intersect(Ray &ray);
    Extent  GetExtent();
    Solid   *Optimize() { if(ps) ps = ps->Optimize(); return this; }
};

extern Solid* MakeUnion(Solid *ps1);
extern Solid* MakeUnion(Solid *ps1, Solid *ps2);
extern Solid* MakeUnion(Solid *ps1, Solid *ps2, Solid *ps3);
extern Solid* MakeUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4);
extern Solid* MakeUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5);
extern Solid* MakeUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6);
extern Solid* MakeUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7);
extern Solid* MakeUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8);
extern Solid* MakeUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9);
extern Solid* MakeUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9,
                        Solid *ps10);
extern Solid* MakeUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9,
                        Solid *ps10, Solid *ps11);
extern Solid* MakeUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9,
                        Solid *ps10, Solid *ps11, Solid *ps12);
extern Solid* MakeUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9,
                        Solid *ps10, Solid *ps11, Solid *ps12, Solid *ps13);
extern Solid* MakeUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9,
                        Solid *ps10, Solid *ps11, Solid *ps12, Solid *ps13, Solid *ps14);
extern Solid* MakeUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9,
                        Solid *ps10, Solid *ps11, Solid *ps12, Solid *ps13, Solid *ps14, Solid *ps15);
extern Solid* MakeUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9,
                        Solid *ps10, Solid *ps11, Solid *ps12, Solid *ps13, Solid *ps14, Solid *ps15, Solid *ps16);
extern Solid* MakeUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9,
                        Solid *ps10, Solid *ps11, Solid *ps12, Solid *ps13, Solid *ps14, Solid *ps15, Solid *ps16, Solid *ps17);
extern Solid* MakeUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9,
                        Solid *ps10, Solid *ps11, Solid *ps12, Solid *ps13, Solid *ps14, Solid *ps15, Solid *ps16, Solid *ps17,
                        Solid *ps18);

extern Solid* MakeTopUnion(Solid *ps1);
extern Solid* MakeTopUnion(Solid *ps1, Solid *ps2);
extern Solid* MakeTopUnion(Solid *ps1, Solid *ps2, Solid *ps3);
extern Solid* MakeTopUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4);
extern Solid* MakeTopUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5);
extern Solid* MakeTopUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6);
extern Solid* MakeTopUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7);
extern Solid* MakeTopUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8);
extern Solid* MakeTopUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9);
extern Solid* MakeTopUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9,
                        Solid *ps10);
extern Solid* MakeTopUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9,
                        Solid *ps10, Solid *ps11);
extern Solid* MakeTopUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9,
                        Solid *ps10, Solid *ps11, Solid *ps12);
extern Solid* MakeTopUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9,
                        Solid *ps10, Solid *ps11, Solid *ps12, Solid *ps13);
extern Solid* MakeTopUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9,
                        Solid *ps10, Solid *ps11, Solid *ps12, Solid *ps13, Solid *ps14);
extern Solid* MakeTopUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9,
                        Solid *ps10, Solid *ps11, Solid *ps12, Solid *ps13, Solid *ps14, Solid *ps15);
extern Solid* MakeTopUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9,
                        Solid *ps10, Solid *ps11, Solid *ps12, Solid *ps13, Solid *ps14, Solid *ps15, Solid *ps16);
extern Solid* MakeTopUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9,
                        Solid *ps10, Solid *ps11, Solid *ps12, Solid *ps13, Solid *ps14, Solid *ps15, Solid *ps16, Solid *ps17);
extern Solid* MakeTopUnion(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8, Solid *ps9,
                        Solid *ps10, Solid *ps11, Solid *ps12, Solid *ps13, Solid *ps14, Solid *ps15, Solid *ps16, Solid *ps17,
                        Solid *ps18);

extern Solid* MakeIntersection(Solid *ps1);
extern Solid* MakeIntersection(Solid *ps1, Solid *ps2);
extern Solid* MakeIntersection(Solid *ps1, Solid *ps2, Solid *ps3);
extern Solid* MakeIntersection(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4);
extern Solid* MakeIntersection(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5);
extern Solid* MakeIntersection(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6);
extern Solid* MakeIntersection(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7);
extern Solid* MakeIntersection(Solid *ps1, Solid *ps2, Solid *ps3, Solid *ps4, Solid *ps5, Solid*ps6, Solid *ps7, Solid *ps8);
//
//  Hierarchical bounding boxes
//
struct BoundingBox : Solid {
    Extent  ex;
    Solid   *ps;

            BoundingBox(Solid *p) : ps(p) { ex = ps->GetExtent(); if(ex.NullExtent()) ps = new NullSolid; }
    Hit     *Intersect(Ray &r) { if (ex.TestIntersection(r)) return ps->Intersect(r); else return 0; }
    Extent  GetExtent() { return ex; }
    Solid   *Optimize() { ps = ps->Optimize(); return this; }
};
//
//  Shape Library
//
extern Solid *Cuboid(double r, double x = 1.0, double y = 1.0, double z = 1.0);
extern Solid* CylinderFromTo(const Point3 &pFrom, const Point3 &pTo, double r);
extern Solid* TetrahedronCage(double r, double rs);
extern Solid* HexahedronCage(double r, double rs);
extern Solid* OctahedronCage(double r, double rs);
extern Solid* IcosahedronCage(double r, double rs);
extern Solid* DodecahedronCage(double r, double rs);
extern Solid* PseudoHelix(double r, double zRepeat, double fAng = 2.0*D_PI);
extern Solid* Diskoid(double r);
extern Solid* Rodoid(double r);
extern Solid* Squareoid(double r);
extern Solid* Knoboid(double r, double z = 1.0);
extern Solid* Knoboid2(double r, double z = 1.0);
extern Solid* Hyperboloidoid(double r, double rt, double z = 1.0);
extern Solid* Slaboid(double x, double y, double z, double r);
extern Solid* Slaboidoid(double x, double y, double z, double r, double r2);
extern Solid* SolidText(double xStart, double yStart, char *sz, int nJust, double r);