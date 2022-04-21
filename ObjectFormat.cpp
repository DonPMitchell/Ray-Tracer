#include "stdafx.h"
#include "RayTracer2020.h"
#include "Affine.h"
#include "Image.h"
#define MAXVERT 200000
#define WIDE    1920
#define HIGH    1080
#define S0      0
#define S1      1
#define S2      2

struct Triangle {
    int iV0, iV1, iV2;

        Triangle() {}
        Triangle(int i0, int i1, int i2) : iV0(i0), iV1(i1), iV2(i2) {}
};

struct Mesh {
    Extent e;
    int iFirst, iEnd;
};

static Point3       pMin, pMax;
static Point3       s_rgpVertices[MAXVERT];
static Triangle     s_rgt[3*MAXVERT];
static Mesh         s_rgm[10000];
static int          nVertices, nTriangles, nMeshes;
static float        zBuffer[HIGH][WIDE];

int
LoadObjFile(char *szFile)
{
    FILE *pf;
    char szLine[64];
    int iV0, iV1, iV2, iVT, nState;
    double x, y, z;
    Point3 p;

    if (fopen_s(&pf, szFile, "r")) {
        printf("Failed to open %s\n");
        return 0;
    }
    nState = S0;
    while (fgets(szLine, sizeof(szLine), pf)) {
        if (szLine[0] == 'v' && szLine[1] == ' ') {
            if (nState == S2)
                s_rgm[nMeshes++].iEnd = nTriangles;
            sscanf_s(szLine,"v %lf%lf%lf", &x, &y, &z);
            p = Point3(x, y, z);
            if (nVertices == 0) {
                pMin = pMax = p;
            } else {
                pMin = Min(pMin, p);
                pMax = Max(pMax, p);
            }
            s_rgpVertices[nVertices++]  = p;
            nState = S1;
        }
        if (szLine[0] == 'f' && szLine[1] == ' ') {
            if (nState = S1)
                s_rgm[nMeshes].iFirst = nTriangles;
            sscanf_s(szLine, "f %d/%d %d/%d %d/%d", &iV0, &iVT, &iV1, &iVT, &iV2, &iVT);
            s_rgt[nTriangles++] = Triangle(iV0, iV1, iV2);
            nState = S2;
        }
            
    }
    if (nState == S2)
        s_rgm[nMeshes++].iEnd = nTriangles;
    printf("%d vertices, %d triangles, %d meshes\n", nVertices, nTriangles, nMeshes);
    pMin.Print();
    pMax.Print();
    return nVertices;
}
//
//  Moeller-Trumbore algorithm
//
double
TriangleIntersect(const Triangle &tri, const Ray &ray)
{
    Vector3 e1, e2, h, s, q;
    double a, f, u, v, t;

    e1 = s_rgpVertices[tri.iV1] - s_rgpVertices[tri.iV0];
    e2 = s_rgpVertices[tri.iV2] - s_rgpVertices[tri.iV0];
    h = Cross(ray.vDir, e2);
    a = Dot(e1, h);
    if (fabs(a) < EPSILON)
        return -1000000.0;
    f = 1.0/a;
    s = ray.pOrg - s_rgpVertices[tri.iV0];
    u = f * Dot(s, h);
    if (u < 0.0 || u > 1.0)
        return -1000000.0;
    q = Cross(s, e1);
    v = f * Dot(ray.vDir, q);
    if (v < 0.0 || v > 1.0)
        return -1000000.0;
    t = f * Dot(f, q);
    return t;
}

void
FastRender(int iMesh, char *szFile)
{
    Image im;
    Vector3 e1, e2, vN, vL, vert;
    Point3 p, puv;
    Triangle tri;
    double x, y, z, f, u, v;
    int i, j, k;
    Mesh m;

    im.NewImage(WIDE, HIGH);
    for (i = 0; i < WIDE; i++) {
        for (j = 0; j < HIGH; j++) {
            zBuffer[j][i] = -1000000.0;
        }
    }
    vL = Normalize(Vector3(1.0, 1.0, 1.0));
    m = s_rgm[iMesh];
    for (k = m.iFirst; k < m.iEnd; k++) {
        tri = s_rgt[k];
        p = s_rgpVertices[tri.iV0];
        e1 = s_rgpVertices[tri.iV1] - s_rgpVertices[tri.iV0];
        e2 = s_rgpVertices[tri.iV2] - s_rgpVertices[tri.iV0];
        vN = Normalize(Cross(e1, e2));
        f = fabs(Dot(vN, vL));
        for (u = 0.0; u <= 1.0; u += 0.003) {
            for (v = 0.0; v <= 1.0-u; v += 0.003) {
                puv = p + e1*u + e2*v;
                vert = (puv - pMin) / 3.0;
                z = vert.x;
                x = vert.y;
                y = vert.z;
                i = int(x * HIGH) + (WIDE-HIGH)/2;
                j = HIGH - int(y * HIGH) - 1;
                if (i < 0 || j < 0 || i > WIDE-1 || j > HIGH-1)
                    continue;
                if (z > zBuffer[j][i]) {
                    zBuffer[j][i] = float(z);
                    im.Set(f, i, j);
                }
            }
        }
    }
    im.WriteBMP(szFile);
}

void
TestObjFiles()
{
    int i;
    char sz[32];

    LoadObjFile("exported_from_VENERA13_v3_1_.obj");
    for (i = 0; i < nMeshes; i++) {
        sprintf_s(sz, sizeof(sz), "Meshes\\Mesh_%03d.bmp", i);
        FastRender(i, sz);
    }
}