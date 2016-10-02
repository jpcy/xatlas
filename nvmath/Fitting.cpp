// This code is in the public domain -- Ignacio Castaño <castano@gmail.com>

#include "Fitting.h"
#include "Vector.h"

#include "nvcore/Array.h"
#include "nvcore/Utils.h" // max, swap

#include <float.h> // FLT_MAX
//#include <vector>
#include <string.h>

using namespace nv;

// @@ Move to EigenSolver.h

Vector3 nv::Fit::computeCentroid(int n, const Vector3 *__restrict points)
{
    Vector3 centroid(0.0f);

    for (int i = 0; i < n; i++)
    {
        centroid += points[i];
    }
    centroid /= float(n);

    return centroid;
}

Vector3 nv::Fit::computeCovariance(int n, const Vector3 *__restrict points, float *__restrict covariance)
{
    // compute the centroid
    Vector3 centroid = computeCentroid(n, points);

    // compute covariance matrix
    for (int i = 0; i < 6; i++)
    {
        covariance[i] = 0.0f;
    }

    for (int i = 0; i < n; i++)
    {
        Vector3 v = points[i] - centroid;

        covariance[0] += v.x * v.x;
        covariance[1] += v.x * v.y;
        covariance[2] += v.x * v.z;
        covariance[3] += v.y * v.y;
        covariance[4] += v.y * v.z;
        covariance[5] += v.z * v.z;
    }

    return centroid;
}

bool nv::Fit::isPlanar(int n, const Vector3 * points, float epsilon/*=NV_EPSILON*/)
{
    // compute the centroid and covariance
    float matrix[6];
    computeCovariance(n, points, matrix);

    float eigenValues[3];
    Vector3 eigenVectors[3];
    if (!eigenSolveSymmetric3(matrix, eigenValues, eigenVectors)) {
        return false;
    }

    return eigenValues[2] < epsilon;
}

// Tridiagonal solver from Charles Bloom. 
// Householder transforms followed by QL decomposition. 
// Seems to be based on the code from Numerical Recipes in C.

static void EigenSolver3_Tridiagonal(float mat[3][3], float * diag, float * subd);
static bool EigenSolver3_QLAlgorithm(float mat[3][3], float * diag, float * subd);

bool nv::Fit::eigenSolveSymmetric3(const float matrix[6], float eigenValues[3], Vector3 eigenVectors[3])
{
    nvDebugCheck(matrix != NULL && eigenValues != NULL && eigenVectors != NULL);

    float subd[3];
    float diag[3];
    float work[3][3];

    work[0][0] = matrix[0];
    work[0][1] = work[1][0] = matrix[1];
    work[0][2] = work[2][0] = matrix[2];
    work[1][1] = matrix[3];
    work[1][2] = work[2][1] = matrix[4];
    work[2][2] = matrix[5];

    EigenSolver3_Tridiagonal(work, diag, subd);
    if (!EigenSolver3_QLAlgorithm(work, diag, subd))
    {
        for (int i = 0; i < 3; i++) {
            eigenValues[i] = 0;
            eigenVectors[i] = Vector3(0);
        }
        return false;
    }

    for (int i = 0; i < 3; i++) {
        eigenValues[i] = (float)diag[i];
    }

    // eigenvectors are the columns; make them the rows :

    for (int i=0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            eigenVectors[j].component[i] = (float) work[i][j];
        }
    }

    // shuffle to sort by singular value :
    if (eigenValues[2] > eigenValues[0] && eigenValues[2] > eigenValues[1])
    {
        swap(eigenValues[0], eigenValues[2]);
        swap(eigenVectors[0], eigenVectors[2]);
    }
    if (eigenValues[1] > eigenValues[0])
    {
        swap(eigenValues[0], eigenValues[1]);
        swap(eigenVectors[0], eigenVectors[1]);
    }
    if (eigenValues[2] > eigenValues[1])
    {
        swap(eigenValues[1], eigenValues[2]);
        swap(eigenVectors[1], eigenVectors[2]);
    }

    nvDebugCheck(eigenValues[0] >= eigenValues[1] && eigenValues[0] >= eigenValues[2]);
    nvDebugCheck(eigenValues[1] >= eigenValues[2]);

    return true;
}

static void EigenSolver3_Tridiagonal(float mat[3][3], float * diag, float * subd)
{
    // Householder reduction T = Q^t M Q
    //   Input:   
    //     mat, symmetric 3x3 matrix M
    //   Output:  
    //     mat, orthogonal matrix Q
    //     diag, diagonal entries of T
    //     subd, subdiagonal entries of T (T is symmetric)
    const float epsilon = 1e-08f;

    float a = mat[0][0];
    float b = mat[0][1];
    float c = mat[0][2];
    float d = mat[1][1];
    float e = mat[1][2];
    float f = mat[2][2];

    diag[0] = a;
    subd[2] = 0.f;
    if (fabsf(c) >= epsilon)
    {
        const float ell = sqrtf(b*b+c*c);
        b /= ell;
        c /= ell;
        const float q = 2*b*e+c*(f-d);
        diag[1] = d+c*q;
        diag[2] = f-c*q;
        subd[0] = ell;
        subd[1] = e-b*q;
        mat[0][0] = 1; mat[0][1] = 0; mat[0][2] = 0;
        mat[1][0] = 0; mat[1][1] = b; mat[1][2] = c;
        mat[2][0] = 0; mat[2][1] = c; mat[2][2] = -b;
    }
    else
    {
        diag[1] = d;
        diag[2] = f;
        subd[0] = b;
        subd[1] = e;
        mat[0][0] = 1; mat[0][1] = 0; mat[0][2] = 0;
        mat[1][0] = 0; mat[1][1] = 1; mat[1][2] = 0;
        mat[2][0] = 0; mat[2][1] = 0; mat[2][2] = 1;
    }
}

static bool EigenSolver3_QLAlgorithm(float mat[3][3], float * diag, float * subd)
{
    // QL iteration with implicit shifting to reduce matrix from tridiagonal
    // to diagonal
    const int maxiter = 32;

    for (int ell = 0; ell < 3; ell++)
    {
        int iter;
        for (iter = 0; iter < maxiter; iter++)
        {
            int m;
            for (m = ell; m <= 1; m++)
            {
                float dd = fabsf(diag[m]) + fabsf(diag[m+1]);
                if ( fabsf(subd[m]) + dd == dd )
                    break;
            }
            if ( m == ell )
                break;

            float g = (diag[ell+1]-diag[ell])/(2*subd[ell]);
            float r = sqrtf(g*g+1);
            if ( g < 0 )
                g = diag[m]-diag[ell]+subd[ell]/(g-r);
            else
                g = diag[m]-diag[ell]+subd[ell]/(g+r);
            float s = 1, c = 1, p = 0;
            for (int i = m-1; i >= ell; i--)
            {
                float f = s*subd[i], b = c*subd[i];
                if ( fabsf(f) >= fabsf(g) )
                {
                    c = g/f;
                    r = sqrtf(c*c+1);
                    subd[i+1] = f*r;
                    c *= (s = 1/r);
                }
                else
                {
                    s = f/g;
                    r = sqrtf(s*s+1);
                    subd[i+1] = g*r;
                    s *= (c = 1/r);
                }
                g = diag[i+1]-p;
                r = (diag[i]-g)*s+2*b*c;
                p = s*r;
                diag[i+1] = g+p;
                g = c*r-b;

                for (int k = 0; k < 3; k++)
                {
                    f = mat[k][i+1];
                    mat[k][i+1] = s*mat[k][i]+c*f;
                    mat[k][i] = c*mat[k][i]-s*f;
                }
            }
            diag[ell] -= p;
            subd[ell] = g;
            subd[m] = 0;
        }

        if ( iter == maxiter )
            // should not get here under normal circumstances
            return false;
    }

    return true;
}
