// This code is in the public domain -- castanyo@yahoo.es

#pragma once
#ifndef NV_MATH_SOLVER_H
#define NV_MATH_SOLVER_H

#include "nvmath.h"

namespace nv
{
    class SparseMatrix;
    class FullVector;


    // Linear solvers.
    bool LeastSquaresSolver(const SparseMatrix & A, const FullVector & b, FullVector & x, float epsilon = 1e-5f);
    bool LeastSquaresSolver(const SparseMatrix & A, const FullVector & b, FullVector & x, const uint * lockedParameters, uint lockedCount, float epsilon = 1e-5f);
    bool SymmetricSolver(const SparseMatrix & A, const FullVector & b, FullVector & x, float epsilon = 1e-5f);
    //void NonSymmetricSolver(const SparseMatrix & A, const FullVector & b, FullVector & x, float epsilon = 1e-5f);

} // nv namespace


#endif // NV_MATH_SOLVER_H
