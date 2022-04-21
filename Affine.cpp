#include "stdafx.h"
#include "Affine.h"
#pragma intrinsic(sin, cos)
#define D_PI        3.14159265358979323846264338327950288419716939937510

Matrix3
InverseMatrix(const Matrix3 &mat)
{
    double fDet, fDetInv, f1, f2, f3;

    fDet = mat.m[0][0] * ( f1 = mat.m[1][1] * mat.m[2][2] - mat.m[1][2] * mat.m[2][1] )
         + mat.m[0][1] * ( f2 = mat.m[1][2] * mat.m[2][0] - mat.m[1][0] * mat.m[2][2] )
         + mat.m[0][2] * ( f3 = mat.m[1][0] * mat.m[2][1] - mat.m[1][1] * mat.m[2][0] );
    if (fDet)
        fDetInv = (1.0/fDet);
    else
        fDetInv = 0.0;
    return Matrix3(
          f1 * fDetInv,
        - ( mat.m[0][1] * mat.m[2][2] - mat.m[0][2] * mat.m[2][1] ) * fDetInv,
          ( mat.m[0][1] * mat.m[1][2] - mat.m[0][2] * mat.m[1][1] ) * fDetInv,

          f2 * fDetInv,
          ( mat.m[0][0] * mat.m[2][2] - mat.m[0][2] * mat.m[2][0] ) * fDetInv,
        - ( mat.m[0][0] * mat.m[1][2] - mat.m[0][2] * mat.m[1][0] ) * fDetInv,

          f3 * fDetInv,
        - ( mat.m[0][0] * mat.m[2][1] - mat.m[0][1] * mat.m[2][0] ) * fDetInv,
          ( mat.m[0][0] * mat.m[1][1] - mat.m[0][1] * mat.m[1][0] ) * fDetInv
    );
}

//
//  Rotation by Euler angles
//
static void
RotateBasisAxis(int i1, int i2, double fRadians, Matrix3 &mat)
{
    double s, c, x;
    int j;

    s = sin(fRadians);
    c = cos(fRadians);
    //
    //  efficient multiplication by an axis-rotation matrix
    //
    for (j = 0; j < 3; j++) {
        x =            c * mat.m[i1][j] - s * mat.m[i2][j];
        mat.m[i2][j] = s * mat.m[i1][j] + c * mat.m[i2][j];
        mat.m[i1][j] = x;
    }    
}

static void
InitialBasisAxis(int i1, int i2, double fRadians, Matrix3 &mat)
{
    double s, c;

    s = sin(fRadians);
    c = cos(fRadians);
    mat.m[i1][i1] = c;
    mat.m[i1][i2] = -s;
    mat.m[i2][i1] = s;
    mat.m[i2][i2] = c;
}

Matrix3
EulerMatrix(double fRadians0, double fRadians1, double fRadians2, char *szConvention)
{
    Matrix3 mat(1.0);

    switch(szConvention[0]) {
        case 'x':   InitialBasisAxis(1, 2, fRadians0, mat);
                    break;
        case 'y':   InitialBasisAxis(2, 0, fRadians0, mat);
                    break;
        case 'z':   InitialBasisAxis(0, 1, fRadians0, mat);
                    break;
    }
    switch(szConvention[1]) {
        case 'x':   RotateBasisAxis(1, 2, fRadians1, mat);
                    break;
        case 'y':   RotateBasisAxis(2, 0, fRadians1, mat);
                    break;
        case 'z':   RotateBasisAxis(0, 1, fRadians1, mat);
                    break;
    }
    switch(szConvention[2]) {
        case 'x':   RotateBasisAxis(1, 2, fRadians2, mat);
                    break;
        case 'y':   RotateBasisAxis(2, 0, fRadians2, mat);
                    break;
        case 'z':   RotateBasisAxis(0, 1, fRadians2, mat);
                    break;
    }
    return mat;
}

Matrix3
RotateAboutAxis(const Vector3 &vAxis, double fRadians)
{
    double s, c, t;
    Matrix3 mat;
    Vector3 vNorm;

    vNorm = Normalize(vAxis);
    s = sin(-fRadians);
    c = cos(-fRadians);
    t = 1.0 - c;
    mat.m[0][0] = t * vNorm.x * vNorm.x + c;    //REVIEW: 9 redundant multiplies
    mat.m[1][1] = t * vNorm.y * vNorm.y + c;
    mat.m[2][2] = t * vNorm.z * vNorm.z + c;
    mat.m[0][1] = t * vNorm.x * vNorm.y + s * vNorm.z;
    mat.m[1][0] = t * vNorm.x * vNorm.y - s * vNorm.z;
    mat.m[0][2] = t * vNorm.x * vNorm.z - s * vNorm.y;
    mat.m[2][0] = t * vNorm.x * vNorm.z + s * vNorm.y;
    mat.m[1][2] = t * vNorm.y * vNorm.z + s * vNorm.x;
    mat.m[2][1] = t * vNorm.y * vNorm.z - s * vNorm.x;
    return mat;
}

Matrix3
RotateFromTo(const Vector3 &vFrom, const Vector3 &vTo)
{
    double s, c, t;
    Matrix3 mat;
    Vector3 vSNorm, vDNorm, vNorm;

    vSNorm = Normalize(vFrom);
    vDNorm = Normalize(vTo);
    vNorm = Cross(vSNorm, vDNorm);
    s = Norm(vNorm);                        // sine
    c = Dot(vSNorm, vDNorm);                // cosine
    if (s < 1.0e-8) {
        if (c > 0.0)
            return Matrix3(1.0);
        else
            return RotateAboutAxis(Vector3(1.0, 0.0, 0.0), D_PI);
    }
    vNorm = vNorm/(s);
    s = -s;
    t = 1.0 - c;
    mat.m[0][0] = (t * vNorm.x * vNorm.x + c);
    mat.m[1][1] = (t * vNorm.y * vNorm.y + c);
    mat.m[2][2] = (t * vNorm.z * vNorm.z + c);
    mat.m[0][1] = (t * vNorm.x * vNorm.y + s * vNorm.z);
    mat.m[1][0] = (t * vNorm.x * vNorm.y - s * vNorm.z);
    mat.m[0][2] = (t * vNorm.x * vNorm.z - s * vNorm.y);
    mat.m[2][0] = (t * vNorm.x * vNorm.z + s * vNorm.y);
    mat.m[1][2] = (t * vNorm.y * vNorm.z + s * vNorm.x);
    mat.m[2][1] = (t * vNorm.y * vNorm.z - s * vNorm.x);
    return mat;
}
//
//  Modified Gram-Schmidt method used to create a QR decomposition.  Numerically
//  stable and fast for small matrices.
//
//  Returns det(R) == fabs(det(A)).  Note that det(Q) = +/- 1.0
//
double
OrthogonalizeMatrix(const Matrix3 &mA, Matrix3 &mQ, Matrix3 &mR)
{
    int j, k;
    double fRkk, fRkj, fDet;

    mQ = mA;
    mR = Matrix3(1.0);
    fDet = 1.0;
    for (k = 0; k < 3; k++) {
        fRkk = sqrt(  double(mQ.m[0][k])*double(mQ.m[0][k])
                    + double(mQ.m[1][k])*double(mQ.m[1][k])
                    + double(mQ.m[2][k])*double(mQ.m[2][k]));
        mR.m[k][k] = (fRkk);
        fDet *= fRkk;
        if (fRkk)
            fRkk = 1.0/fRkk;
        else
            return 0.0;
        mQ.m[0][k] = (mQ.m[0][k] * fRkk);
        mQ.m[1][k] = (mQ.m[1][k] * fRkk);
        mQ.m[2][k] = (mQ.m[2][k] * fRkk);
        for (j = k + 1; j < 3; j++) {
            fRkj = double(mQ.m[0][j])*double(mQ.m[0][k])
                 + double(mQ.m[1][j])*double(mQ.m[1][k])
                 + double(mQ.m[2][j])*double(mQ.m[2][k]);
            mR.m[k][j] = (fRkj);
            mQ.m[0][j] = (mQ.m[0][j] - mQ.m[0][k] * fRkj);
            mQ.m[1][j] = (mQ.m[1][j] - mQ.m[1][k] * fRkj);
            mQ.m[2][j] = (mQ.m[2][j] - mQ.m[2][k] * fRkj);
        }
    }
    return (fDet);
}

//
//  Unit bases
//
Vector3     g_vx(1.0, 0.0, 0.0);
Vector3     g_vy(0.0, 1.0, 0.0);
Vector3     g_vz(0.0, 0.0, 1.0);

CoVector3   g_cvx(1.0, 0.0, 0.0);
CoVector3   g_cvy(0.0, 1.0, 0.0);
CoVector3   g_cvz(0.0, 0.0, 1.0);