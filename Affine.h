//
//  3D affine space: points, vectors, co-vectors, transformations
//  D.P. Mitchell  2018/09/29.
//
#pragma once
#include "stdafx.h"
#pragma intrinsic (sqrt, sin, cos)
//
//  Vector space
//
struct Vector3 {
    union {
        struct { double x, y, z; };
        double rgf[3];
    };

            Vector3() {};
            Vector3(double f) : x(f), y(f), z(f) {};
            Vector3(double f1, double f2, double f3) : x(f1), y(f2), z(f3) {};

    Vector3  operator +(const Vector3 &v) const { return Vector3(x+v.x, y+v.y, z+v.z); }
    Vector3  operator -(const Vector3 &v) const { return Vector3(x-v.x, y-v.y, z-v.z); }
    Vector3  operator -()         const { return Vector3(-x, -y, -z); }
    Vector3  operator *(double f) const { return Vector3(x*f, y*f, z*f); }
    Vector3  operator /(double f) const { return Vector3(x/f, y/f, z/f); }      // the compiler only computes 1/f once

    void    Print() const { printf("(%f %f %f)", x, y, z); }
};


inline double
Dot(const Vector3 &u, const Vector3 &v)
{
    return u.x*v.x + u.y*v.y + u.z*v.z;
}

inline Vector3
Cross(const Vector3 &u, const Vector3 &v)
{
    return Vector3(u.y*v.z-u.z*v.y, u.z*v.x-u.x*v.z, u.x*v.y-u.y*v.x);
}

inline double
Norm(const Vector3 &v)
{
    return sqrt(Dot(v, v));
}

inline Vector3
Normalize(const Vector3 &v)
{
    return v / Norm(v);
}
//
//  Affine point space
//
struct Point3 {
    union {
        struct { double x, y, z; };
        double rgf[3];
    };
            Point3() : x(0.0), y(0.0), z(0.0) {};
            Point3(double a, double b, double c) : x(a), y(b), z(c) {};

    Point3   operator +(const Vector3 &v) const { return Point3(x+v.x, y+v.y, z+v.z); }
    Point3   operator -(const Vector3 &v) const { return Point3(x-v.x, y-v.y, z-v.z); }
    Vector3  operator -(const Point3 &p)  const { return Vector3 (x-p.x, y-p.y, z-p.z); }
    int      operator <(const Point3 &p)  const { return x < p.x && y < p.y && z < p.z; }
    int      operator >(const Point3 &p)  const { return x > p.x && y > p.y && z > p.z; }
    void    Print() const { printf("(%f %f %f)", x, y, z); }
};

inline Point3
Blend(double u, const Point3 &p1, const Point3 &p2)
{
    return p1 + (p2 - p1)*u;
}

inline double
Max(double f1, double f2)
{
    if (f1 > f2)
        return f1;
    else
        return f2;
}

inline Point3
Max(const Point3 &p1, const Point3 &p2)
{
    return Point3(Max(p1.x, p2.x), Max(p1.y, p2.y), Max(p1.z, p2.z));
}

inline double
Min(double f1, double f2)
{
    if (f1 < f2)
        return f1;
    else
        return f2;
}

inline Point3
Min(const Point3 &p1, const Point3 &p2)
{
    return Point3(Min(p1.x, p2.x), Min(p1.y, p2.y), Min(p1.z, p2.z));
}

//
//  Matrix, linear transformations
//
struct Matrix3 {
    union {
        double   m[3][3];            // m[iRow][iColumn]
        struct {
            double f00, f01, f02;
            double f10, f11, f12;
            double f20, f21, f22;
        };
    };
                Matrix3() {};
                Matrix3(double x) :
                                f00(x),   f01(0.0), f02(0.0),
                                f10(0.0), f11(x),   f12(0.0),
                                f20(0.0), f21(0.0), f22(x)      {};
                Matrix3(double x00, double x01, double x02,
                       double x10, double x11, double x12,
                       double x20, double x21, double x22) :
                                f00(x00), f01(x01), f02(x02),
                                f10(x10), f11(x11), f12(x12),
                                f20(x20), f21(x21), f22(x22)    {};

    Matrix3  operator +(const Matrix3 &m2) const;
    Matrix3  operator -(const Matrix3 &m2) const;
    Matrix3  operator *(double f) const;
    Matrix3  operator /(double f) const;
    Matrix3  operator *(const Matrix3 &m2) const;
    Vector3  operator *(const Vector3 &v) const;        // v' = M*v
    Point3   operator *(const Point3 &p) const;         // not affine, but useful

    double  Determinant() const;
    Matrix3 Transpose() const;
    void    Print() const { printf("[%9f %9f %9f]\n[%9f %9f %9f]\n[%9f %9f %9f]\n", 
                                f00, f01, f02, f10, f11, f12, f20, f21, f22); }
};

inline Matrix3
Matrix3::operator +(const Matrix3 &m2) const {
    return Matrix3 ( f00 + m2.f00, f01 + m2.f01, f02 + m2.f02,
                     f10 + m2.f10, f11 + m2.f11, f12 + m2.f12,
                     f20 + m2.f20, f21 + m2.f21, f22 + m2.f22);
}

inline Matrix3
Matrix3::operator -(const Matrix3 &m2) const {
    return Matrix3( f00 - m2.f00, f01 - m2.f01, f02 - m2.f02,
                    f10 - m2.f10, f11 - m2.f11, f12 - m2.f12,
                    f20 - m2.f20, f21 - m2.f21, f22 - m2.f22);
}

inline Matrix3
Matrix3::operator *(double f) const {
    return Matrix3( f00 * f, f01 * f, f02 * f,
                    f10 * f, f11 * f, f12 * f,
                    f20 * f, f21 * f, f22 * f);
}

inline Matrix3
Matrix3::operator /(double f) const {
    return *this * (1.0/f);
}

inline Matrix3
Matrix3::operator *(const Matrix3 &m2) const
{
    return Matrix3(
        m[0][0]*m2.m[0][0] + m[0][1]*m2.m[1][0] + m[0][2]*m2.m[2][0],
        m[0][0]*m2.m[0][1] + m[0][1]*m2.m[1][1] + m[0][2]*m2.m[2][1],
        m[0][0]*m2.m[0][2] + m[0][1]*m2.m[1][2] + m[0][2]*m2.m[2][2],

        m[1][0]*m2.m[0][0] + m[1][1]*m2.m[1][0] + m[1][2]*m2.m[2][0],
        m[1][0]*m2.m[0][1] + m[1][1]*m2.m[1][1] + m[1][2]*m2.m[2][1],
        m[1][0]*m2.m[0][2] + m[1][1]*m2.m[1][2] + m[1][2]*m2.m[2][2],

        m[2][0]*m2.m[0][0] + m[2][1]*m2.m[1][0] + m[2][2]*m2.m[2][0],
        m[2][0]*m2.m[0][1] + m[2][1]*m2.m[1][1] + m[2][2]*m2.m[2][1],
        m[2][0]*m2.m[0][2] + m[2][1]*m2.m[1][2] + m[2][2]*m2.m[2][2]
    );
}

inline Vector3
Matrix3::operator *(const Vector3 &vec) const
{
return Vector3(vec.x*m[0][0] + vec.y*m[0][1] + vec.z*m[0][2],
               vec.x*m[1][0] + vec.y*m[1][1] + vec.z*m[1][2],
               vec.x*m[2][0] + vec.y*m[2][1] + vec.z*m[2][2]);
}

inline Point3
Matrix3::operator *(const Point3 &pnt) const
{
return Point3(pnt.x*m[0][0] + pnt.y*m[0][1] + pnt.z*m[0][2],
              pnt.x*m[1][0] + pnt.y*m[1][1] + pnt.z*m[1][2],
              pnt.x*m[2][0] + pnt.y*m[2][1] + pnt.z*m[2][2]);
}

inline double
Matrix3::Determinant() const
{
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
         + m[0][1] * (m[1][2] * m[2][0] - m[1][0] * m[2][2])
         + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}
inline Matrix3
Matrix3::Transpose() const
{
    return Matrix3 (f00, f10, f20,
                    f01, f11, f21,
                    f02, f12, f22);
}

extern Matrix3 InverseMatrix(const Matrix3 &mat);
extern double OrthogonalizeMatrix(const Matrix3 &mA, Matrix3 &mQ, Matrix3 &mR);

inline Matrix3
ScaleMatrix(double x, double y, double z)
{
    return Matrix3(  x, 0.0, 0.0,
                   0.0,   y, 0.0,
                   0.0, 0.0,   z);
}

extern Matrix3 EulerMatrix(double fRadians0, double fRadians1, double fRadians2, char *szConvention = "xyz");
extern Matrix3 RotateAboutAxis(const Vector3 &vAxis, double fRadians);
extern Matrix3 RotateFromTo(const Vector3 &vFrom, const Vector3 &vTo);
//
//  Dual vector space
//
struct CoVector3 {
    union {
        struct { double x, y, z; };
        double rgf[3];
    };

            CoVector3() {};
            CoVector3(double f) : x(f), y(f), z(f) {};
            CoVector3(double f1, double f2, double f3) : x(f1), y(f2), z(f3) {};

    CoVector3  operator +(const CoVector3 &v) const { return CoVector3(x+v.x, y+v.y, z+v.z); }
    CoVector3  operator -(const CoVector3 &v) const { return CoVector3(x-v.x, y-v.y, z-v.z); }
    CoVector3  operator -()          const { return CoVector3(-x, -y, -z); }
    CoVector3  operator *(double f)  const { return CoVector3(x*f, y*f, z*f); }
    CoVector3  operator /(double f)  const { return CoVector3(x/f, y/f, z/f); }
    CoVector3  operator *(const Matrix3 &m2) const;       // v' = v*M
    double     operator *(const Vector3 &v) const { return x*v.x + y*v.y + z*v.z; }
    double     operator *(const Point3 &p) const { return x*p.x + y*p.y + z*p.z; }      // not affine

    void    Print() const { printf("(%f %f %f)", x, y, z); }
};


inline double
Dot(const CoVector3 &u, const CoVector3 &v)
{
    return u.x*v.x + u.y*v.y + u.z*v.z;
}

inline CoVector3
Cross(const CoVector3 &u, const CoVector3 &v)
{
    return CoVector3(u.y*v.z-u.z*v.y, u.z*v.x-u.x*v.z, u.x*v.y-u.y*v.x);
}

inline double
Norm(const CoVector3 &v)
{
    return sqrt(Dot(v, v));
}

inline CoVector3
Normalize(const CoVector3 &v)
{
    return v / Norm(v);
}

inline CoVector3
CoVector3::operator *(const Matrix3 &m) const
{
return CoVector3( x*m.m[0][0] + y*m.m[1][0] + z*m.m[2][0],
                  x*m.m[0][1] + y*m.m[1][1] + z*m.m[2][1],
                  x*m.m[0][2] + y*m.m[1][2] + z*m.m[2][2] );
}
//
//  2D vector
//
struct Vector2 {
    union {
        struct {
            double x, y;
        };
        double rgf[2];
    };

        Vector2() {}
        Vector2(double a, double b) : x(a), y(b) {}
    Vector2  operator +(const Vector2 &v) const { return Vector2(x+v.x, y+v.y); }
    Vector2  operator -(const Vector2 &v) const { return Vector2(x-v.x, y-v.y); }
    Vector2  operator -()  { return Vector2(-x, -y); }
    Vector2  operator *(double f) const { return Vector2(x*f, y*f); }
    Vector2  operator /(double f) const { return Vector2(x/f, y/f); }
    int      operator ==(const Vector2 &v) const { return x == v.x && y == v.y; }
    void    Print() { printf("(%f %f)", x, y); }
};

inline double
Dot(const Vector2 &u, const Vector2 &v)
{
    return u.x*v.x + u.y*v.y;
}

inline double
Norm(const Vector2 &v)
{
    return sqrt(v.x*v.x + v.y*v.y);
}

inline Vector2
Normalize(const Vector2 &v)
{
    double f;

    f = Norm(v);
    return Vector2(v.x/f, v.y/f);
}

inline int
CCW(const Vector2 &vA, const Vector2 &vB, const Vector2 &vC)
{
    double fPos, fNeg, fEps;

    fPos = (vB.x - vA.x) * (vC.y - vA.y);
    fNeg = (vC.x - vA.x) * (vB.y - vA.y);
    fEps = (fabs(fPos) + fabs(fNeg)) * 5.0e-016;    // 2*machine epsilone
    if (fPos > fNeg + fEps)
        return 1;
    if (fPos < fNeg - fEps)
        return -1;
    return 0;
}
//
//  Unit bases
//
extern Vector3     g_vx, g_vy, g_vz;
extern CoVector3   g_cvx, g_cvy, g_cvz;
//
//  Simple 4-D and 5-D vector classes
//
struct Vector4
{
    union {
        double rgf[4];
        struct { double x, y, z, a; };
    };

            Vector4() {}
            Vector4(double x1, double x2, double x3, double x4) : x(x1), y(x2), z(x3), a(x4) {}
    Vector4 operator +(const Vector4 &v) const { return Vector4(x+v.x, y+v.y, z+v.z, a+v.a); }
    Vector4 operator -(const Vector4 &v) const { return Vector4(x-v.x, y-v.y, z-v.z, a-v.a); }
    Vector4 operator -()         const { return Vector4(-x, -y, -z, -a); }
    Vector4 operator *(double f) const { return Vector4(x*f, y*f, z*f, a*f); }
    Vector4 operator /(double f) const { return Vector4(x/f, y/f, z/f, a/f); }
    double  Norm() { return sqrt(x*x + y*y + z*z + a*a); }
};

struct Vector5
{
    union {
        double rgf[5];
        struct { double x, y, z, a, b; };
    };

            Vector5() {}
            Vector5(double x1, double x2, double x3, double x4, double x5) : x(x1), y(x2), z(x3), a(x4), b(x5) {}
    Vector5 operator +(const Vector5 &v) const { return Vector5(x+v.x, y+v.y, z+v.z, a+v.a, b+v.b); }
    Vector5 operator -(const Vector5 &v) const { return Vector5(x-v.x, y-v.y, z-v.z, a-v.a, b-v.b); }
    Vector5 operator -()         const { return Vector5(-x, -y, -z, -a, -b); }
    Vector5 operator *(double f) const { return Vector5(x*f, y*f, z*f, a*f, b*f); }
    Vector5 operator /(double f) const { return Vector5(x/f, y/f, z/f, a/f, b/f); }
    double  Norm() { return sqrt(x*x + y*y + z*z + a*a + b*b); }
};