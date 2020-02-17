#include "OpenNL_psm.h"

/*
 *  Copyright (c) 2004-2010, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     levy@loria.fr
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine,
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX
 *     FRANCE
 *
 */


 /*
  *  This file is a PSM (pluggable software module)
  *   generated from the distribution of Geogram.
  *
  *  See Geogram documentation on:
  *   http://alice.loria.fr/software/geogram/doc/html/index.html
  *
  *  See documentation of the functions bundled in this PSM on:
  *   http://alice.loria.fr/software/geogram/doc/html/nl_8h.html
  */

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#define nl_arg_used(x) (void)x

#ifndef MIN
#define MIN(x,y) (((x) < (y)) ? (x) : (y)) 
#endif

#ifndef MAX
#define MAX(x,y) (((x) > (y)) ? (x) : (y)) 
#endif

#define NL_NEW(T)                (T*)(calloc(1, sizeof(T))) 
#define NL_NEW_ARRAY(T,NB)       (T*)(calloc((size_t)(NB),sizeof(T)))
#define NL_RENEW_ARRAY(T,x,NB)   (T*)(realloc(x,(size_t)(NB)*sizeof(T))) 
#define NL_DELETE(x)             free(x); x = NULL 
#define NL_DELETE_ARRAY(x)       free(x); x = NULL
#define NL_CLEAR_ARRAY(T,x,NB)   memset(x, 0, (size_t)(NB)*sizeof(T)) 

#define NL_NEW_VECTOR(dim) \
    (double*)malloc((size_t)(dim)*sizeof(double))

#define NL_DELETE_VECTOR(ptr) \
    free(ptr)

struct NLMatrixStruct;
typedef struct NLMatrixStruct* NLMatrix;
typedef void(*NLDestroyMatrixFunc)(NLMatrix M);
typedef void(*NLMultMatrixVectorFunc)(NLMatrix M, const double* x, double* y);

#define NL_MATRIX_SPARSE_DYNAMIC 0x1001
#define NL_MATRIX_CRS            0x1002
#define NL_MATRIX_OTHER          0x1006

struct NLMatrixStruct
{
	uint32_t m;
	uint32_t n;
	uint32_t type;
	NLDestroyMatrixFunc destroy_func;
	NLMultMatrixVectorFunc mult_func;
};

/* Dynamic arrays for sparse row/columns */

struct NLCoeff
{
	uint32_t index;
	double value;
};

struct NLRowColumn
{
	uint32_t size;
	uint32_t capacity;
	NLCoeff* coeff;
};

/* Compressed Row Storage */

struct NLCRSMatrix
{
	uint32_t m;
	uint32_t n;
	uint32_t type;
	NLDestroyMatrixFunc destroy_func;
	NLMultMatrixVectorFunc mult_func;
	double* val;
	uint32_t* rowptr;
	uint32_t* colind;
	uint32_t nslices;
	uint32_t* sliceptr;
};

/* SparseMatrix data structure */

struct NLSparseMatrix
{
	uint32_t m;
	uint32_t n;
	uint32_t type;
	NLDestroyMatrixFunc destroy_func;
	NLMultMatrixVectorFunc mult_func;
	uint32_t diag_size;
	uint32_t diag_capacity;
	NLRowColumn* row;
	NLRowColumn* column;
	double*    diag;
	uint32_t row_capacity;
	uint32_t column_capacity;
};

/* NLContext data structure */

typedef void(*NLProgressFunc)(uint32_t cur_iter, uint32_t max_iter, double cur_err, double max_err);

struct NLBufferBinding
{
	void* base_address;
	uint32_t stride;
};

#define NL_BUFFER_ITEM(B,i) *(double*)((void*)((char*)((B).base_address)+((i)*(B).stride)))

struct NLContext
{
	NLBufferBinding* variable_buffer;
	double*        variable_value;
	bool*       variable_is_locked;
	uint32_t*          variable_index;
	uint32_t           n;
	NLMatrix         M;
	NLMatrix         P;
	NLMatrix         B;
	NLRowColumn      af;
	NLRowColumn      al;
	double*        x;
	double*        b;
	uint32_t           nb_variables;
	uint32_t           nb_systems;
	uint32_t           current_row;
	uint32_t           max_iterations;
	bool        max_iterations_defined;
	uint32_t           inner_iterations;
	double         threshold;
	double         omega;
	uint32_t           used_iterations;
	double         error;
	NLProgressFunc   progress_func;

};

void nlDeleteMatrix(NLMatrix M)
{
	if (M == NULL)
		return;
	M->destroy_func(M);
	NL_DELETE(M);
}

void nlMultMatrixVector(NLMatrix M, const double* x, double* y)
{
	M->mult_func(M, x, y);
}

void nlRowColumnConstruct(NLRowColumn* c)
{
	c->size = 0;
	c->capacity = 0;
	c->coeff = NULL;
}

void nlRowColumnDestroy(NLRowColumn* c)
{
	NL_DELETE_ARRAY(c->coeff);
	c->size = 0;
	c->capacity = 0;
}

void nlRowColumnGrow(NLRowColumn* c)
{
	if (c->capacity != 0) {
		c->capacity = 2 * c->capacity;
		c->coeff = NL_RENEW_ARRAY(NLCoeff, c->coeff, c->capacity);
	} else {
		c->capacity = 4;
		c->coeff = NL_NEW_ARRAY(NLCoeff, c->capacity);
	}
}

void nlRowColumnAdd(NLRowColumn* c, uint32_t index, double value)
{
	for (uint32_t i = 0; i < c->size; i++) {
		if (c->coeff[i].index == index) {
			c->coeff[i].value += value;
			return;
		}
	}
	if (c->size == c->capacity)
		nlRowColumnGrow(c);
	c->coeff[c->size].index = index;
	c->coeff[c->size].value = value;
	c->size++;
}

/* Does not check whether the index already exists */
void nlRowColumnAppend(NLRowColumn* c, uint32_t index, double value)
{
	if (c->size == c->capacity)
		nlRowColumnGrow(c);
	c->coeff[c->size].index = index;
	c->coeff[c->size].value = value;
	c->size++;
}

void nlRowColumnZero(NLRowColumn* c)
{
	c->size = 0;
}

void nlRowColumnClear(NLRowColumn* c)
{
	c->size = 0;
	c->capacity = 0;
	NL_DELETE_ARRAY(c->coeff);
}

static int nlCoeffCompare(const void* p1, const void* p2)
{
	return (((NLCoeff*)(p2))->index < ((NLCoeff*)(p1))->index);
}

void nlRowColumnSort(NLRowColumn* c)
{
	qsort(c->coeff, c->size, sizeof(NLCoeff), nlCoeffCompare);
}

/* CRSMatrix data structure */

static void nlCRSMatrixDestroy(NLCRSMatrix* M)
{
	NL_DELETE_ARRAY(M->val);
	NL_DELETE_ARRAY(M->rowptr);
	NL_DELETE_ARRAY(M->colind);
	NL_DELETE_ARRAY(M->sliceptr);
	M->m = 0;
	M->n = 0;
	M->nslices = 0;
}

static void nlCRSMatrixMultSlice(NLCRSMatrix* M, const double* x, double* y, uint32_t Ibegin, uint32_t Iend)
{
	for (uint32_t i = Ibegin; i < Iend; ++i) {
		double sum = 0.0;
		for (uint32_t j = M->rowptr[i]; j < M->rowptr[i + 1]; ++j)
			sum += M->val[j] * x[M->colind[j]];
		y[i] = sum;
	}
}

static void nlCRSMatrixMult(NLCRSMatrix* M, const double* x, double* y)
{
	int nslices = (int)(M->nslices);
	for (int slice = 0; slice < nslices; ++slice)
		nlCRSMatrixMultSlice(M, x, y, M->sliceptr[slice], M->sliceptr[slice + 1]);
}

void nlCRSMatrixConstruct(NLCRSMatrix* M, uint32_t m, uint32_t n, uint32_t nnz, uint32_t nslices)
{
	M->m = m;
	M->n = n;
	M->type = NL_MATRIX_CRS;
	M->destroy_func = (NLDestroyMatrixFunc)nlCRSMatrixDestroy;
	M->mult_func = (NLMultMatrixVectorFunc)nlCRSMatrixMult;
	M->nslices = nslices;
	M->val = NL_NEW_ARRAY(double, nnz);
	M->rowptr = NL_NEW_ARRAY(uint32_t, m + 1);
	M->colind = NL_NEW_ARRAY(uint32_t, nnz);
	M->sliceptr = NL_NEW_ARRAY(uint32_t, nslices + 1);
}

/* SparseMatrix data structure */

static void nlSparseMatrixDestroyRowColumns(NLSparseMatrix* M)
{
	for (uint32_t i = 0; i < M->m; i++)
		nlRowColumnDestroy(&(M->row[i]));
	NL_DELETE_ARRAY(M->row);
}

void nlSparseMatrixDestroy(NLSparseMatrix* M)
{
	assert(M->type == NL_MATRIX_SPARSE_DYNAMIC);
	nlSparseMatrixDestroyRowColumns(M);
	NL_DELETE_ARRAY(M->diag);
}

void nlSparseMatrixAdd(NLSparseMatrix* M, uint32_t i, uint32_t j, double value)
{
	assert(i >= 0 && i <= M->m - 1);
	assert(j >= 0 && j <= M->n - 1);
	if (i == j)
		M->diag[i] += value;
	nlRowColumnAdd(&(M->row[i]), j, value);
}

static void nlSparseMatrixAddSparseMatrix(NLSparseMatrix* M, double mul, const NLSparseMatrix* N)
{
	assert(M->m == N->m);
	assert(M->n == N->n);
	for (uint32_t i = 0; i < N->m; ++i) {
		for (uint32_t jj = 0; jj < N->row[i].size; ++jj) 
			nlSparseMatrixAdd(M, i, N->row[i].coeff[jj].index, mul*N->row[i].coeff[jj].value);
	}
}

static void nlSparseMatrixAddCRSMatrix(NLSparseMatrix* M, double mul, const NLCRSMatrix* N)
{
	assert(M->m == N->m);
	assert(M->n == N->n);
	for (uint32_t i = 0; i < M->m; ++i) {
		for (uint32_t jj = N->rowptr[i]; jj < N->rowptr[i + 1]; ++jj)
			nlSparseMatrixAdd(M, i, N->colind[jj], mul*N->val[jj]);
	}
}

void nlSparseMatrixAddMatrix(NLSparseMatrix* M, double mul, const NLMatrix N)
{
	assert(M->m == N->m);
	assert(M->n == N->n);
	if (N->type == NL_MATRIX_SPARSE_DYNAMIC)
		nlSparseMatrixAddSparseMatrix(M, mul, (const NLSparseMatrix*)N);
	else if (N->type == NL_MATRIX_CRS)
		nlSparseMatrixAddCRSMatrix(M, mul, (const NLCRSMatrix*)N);
}

void nlSparseMatrixZero(NLSparseMatrix* M)
{
	for (uint32_t i = 0; i < M->m; i++)
		nlRowColumnZero(&(M->row[i]));
	NL_CLEAR_ARRAY(double, M->diag, M->diag_size);
}

void nlSparseMatrixClear(NLSparseMatrix* M)
{
	for (uint32_t i = 0; i < M->m; i++)
		nlRowColumnClear(&(M->row[i]));
	NL_CLEAR_ARRAY(double, M->diag, M->diag_size);
}

/* Returns the number of non-zero coefficients */
uint32_t nlSparseMatrixNNZ(NLSparseMatrix* M)
{
	uint32_t nnz = 0;
	for (uint32_t i = 0; i < M->m; i++)
		nnz += M->row[i].size;
	return nnz;
}

void nlSparseMatrixSort(NLSparseMatrix* M)
{
	for (uint32_t i = 0; i < M->m; i++)
		nlRowColumnSort(&(M->row[i]));
}

void nlSparseMatrixMAddRow(NLSparseMatrix* M, uint32_t i1, double s, uint32_t i2)
{
	NLRowColumn* Ri2 = &(M->row[i2]);
	NLCoeff* c = NULL;
	assert(i1 < M->m);
	assert(i2 < M->m);
	for (uint32_t jj = 0; jj < Ri2->size; ++jj) {
		c = &(Ri2->coeff[jj]);
		nlSparseMatrixAdd(M, i1, c->index, s*c->value);
	}
}

void nlSparseMatrixScaleRow(NLSparseMatrix* M, uint32_t i, double s)
{
	NLRowColumn* Ri = &(M->row[i]);
	NLCoeff* c = NULL;
	assert(i < M->m);
	for (uint32_t jj = 0; jj < Ri->size; ++jj) {
		c = &(Ri->coeff[jj]);
		c->value *= s;
	}
	if (i < M->diag_size)
		M->diag[i] *= s;
}

void nlSparseMatrixZeroRow(NLSparseMatrix* M, uint32_t i)
{
	NLRowColumn* Ri = &(M->row[i]);
	assert(i < M->m);
	Ri->size = 0;
	if (i < M->diag_size)
		M->diag[i] = 0.0;
}

/* SparseMatrix x Vector routines, internal helper routines */

static void nlSparseMatrix_mult_rows(NLSparseMatrix* A,	const double* x, double* y)
{
	/*
	 * Note: OpenMP does not like unsigned ints
	 * (causes some floating point exceptions),
	 * therefore I use here signed ints for all
	 * indices.
	 */
	int m = (int)(A->m);
	NLCoeff* c = NULL;
	NLRowColumn* Ri = NULL;
	for (int i = 0; i < m; i++) {
		Ri = &(A->row[i]);
		y[i] = 0;
		for (int ij = 0; ij < (int)(Ri->size); ij++) {
			c = &(Ri->coeff[ij]);
			y[i] += c->value * x[c->index];
		}
	}
}

void nlSparseMatrixMult(NLSparseMatrix* A, const double* x, double* y)
{
	assert(A->type == NL_MATRIX_SPARSE_DYNAMIC);
	nlSparseMatrix_mult_rows(A, x, y);
}

void nlSparseMatrixConstruct(NLSparseMatrix* M, uint32_t m, uint32_t n)
{
	M->m = m;
	M->n = n;
	M->type = NL_MATRIX_SPARSE_DYNAMIC;
	M->destroy_func = (NLDestroyMatrixFunc)nlSparseMatrixDestroy;
	M->mult_func = (NLMultMatrixVectorFunc)nlSparseMatrixMult;
	M->row = NL_NEW_ARRAY(NLRowColumn, m);
	M->row_capacity = m;
	for (uint32_t i = 0; i < n; i++)
		nlRowColumnConstruct(&(M->row[i]));
	M->row_capacity = 0;
	M->column = NULL;
	M->column_capacity = 0;
	M->diag_size = MIN(m, n);
	M->diag_capacity = M->diag_size;
	M->diag = NL_NEW_ARRAY(double, M->diag_size);
}

NLMatrix nlSparseMatrixNew(uint32_t m, uint32_t n)
{
	NLSparseMatrix* result = NL_NEW(NLSparseMatrix);
	nlSparseMatrixConstruct(result, m, n);
	return (NLMatrix)result;
}

static void adjust_diag(NLSparseMatrix* M)
{
	uint32_t new_diag_size = MIN(M->m, M->n);
	if (new_diag_size > M->diag_size) {
		if (new_diag_size > M->diag_capacity) {
			M->diag_capacity *= 2;
			if (M->diag_capacity == 0)
				M->diag_capacity = 16;
			M->diag = NL_RENEW_ARRAY(double, M->diag, M->diag_capacity);
			for (uint32_t i = M->diag_size; i < new_diag_size; ++i)
				M->diag[i] = 0.0;
		}
		M->diag_size = new_diag_size;
	}
}

void nlSparseMatrixAddRow(NLSparseMatrix* M)
{
	++M->m;
	if (M->m > M->row_capacity) {
		M->row_capacity *= 2;
		if (M->row_capacity == 0)
			M->row_capacity = 16;
		M->row = NL_RENEW_ARRAY(NLRowColumn, M->row, M->row_capacity);
	}
	nlRowColumnConstruct(&(M->row[M->m - 1]));
	adjust_diag(M);
}

void nlSparseMatrixAddColumn(NLSparseMatrix* M)
{
	++M->n;
	adjust_diag(M);
}

NLMatrix nlCRSMatrixNewFromSparseMatrix(NLSparseMatrix* M)
{
	uint32_t nnz = nlSparseMatrixNNZ(M);
	uint32_t nslices = 8; /* TODO: get number of cores */
	uint32_t slice, cur_bound, cur_NNZ, cur_row;
	uint32_t k;
	uint32_t slice_size = nnz / nslices;
	NLCRSMatrix* CRS = NL_NEW(NLCRSMatrix);
	nlCRSMatrixConstruct(CRS, M->m, M->n, nnz, nslices);
	nlSparseMatrixSort(M);
	/* Convert matrix to CRS format */
	k = 0;
	for (uint32_t i = 0; i < M->m; ++i) {
		NLRowColumn* Ri = &(M->row[i]);
		CRS->rowptr[i] = k;
		for (uint32_t ij = 0; ij < Ri->size; ij++) {
			NLCoeff* c = &(Ri->coeff[ij]);
			CRS->val[k] = c->value;
			CRS->colind[k] = c->index;
			++k;
		}
	}
	CRS->rowptr[M->m] = k;
	/* Create "slices" to be used by parallel sparse matrix vector product */
	if (CRS->sliceptr != NULL) {
		cur_bound = slice_size;
		cur_NNZ = 0;
		cur_row = 0;
		CRS->sliceptr[0] = 0;
		for (slice = 1; slice < nslices; ++slice) {
			while (cur_NNZ < cur_bound && cur_row < M->m) {
				++cur_row;
				cur_NNZ += CRS->rowptr[cur_row + 1] - CRS->rowptr[cur_row];
			}
			CRS->sliceptr[slice] = cur_row;
			cur_bound += slice_size;
		}
		CRS->sliceptr[nslices] = M->m;
	}
	return (NLMatrix)CRS;
}

void nlMatrixCompress(NLMatrix* M)
{
	NLMatrix CRS = NULL;
	if ((*M)->type != NL_MATRIX_SPARSE_DYNAMIC)
		return;
	CRS = nlCRSMatrixNewFromSparseMatrix((NLSparseMatrix*)*M);
	nlDeleteMatrix(*M);
	*M = CRS;
}

NLContext *nlNewContext()
{
	NLContext* result = NL_NEW(NLContext);
	result->max_iterations = 100;
	result->threshold = 1e-6;
	result->omega = 1.5;
	result->inner_iterations = 5;
	result->progress_func = NULL;
	result->nb_systems = 1;
	return result;
}

void nlDeleteContext(NLContext *context)
{
	nlDeleteMatrix(context->M);
	context->M = NULL;
	nlDeleteMatrix(context->P);
	context->P = NULL;
	nlDeleteMatrix(context->B);
	context->B = NULL;
	nlRowColumnDestroy(&context->af);
	nlRowColumnDestroy(&context->al);
	NL_DELETE_ARRAY(context->variable_value);
	NL_DELETE_ARRAY(context->variable_buffer);
	NL_DELETE_ARRAY(context->variable_is_locked);
	NL_DELETE_ARRAY(context->variable_index);
	NL_DELETE_ARRAY(context->x);
	NL_DELETE_ARRAY(context->b);
	NL_DELETE(context);
}

static double ddot(int n, const double *x, const double *y)
{
	double sum = 0.0;
	for (int i = 0; i < n; i++)
		sum += x[i] * y[i];
	return sum;
}

static void daxpy(int n, double a, const double *x, double *y)
{
	for (int i = 0; i < n; i++)
		y[i] = a * x[i] + y[i];
}

static void dscal(int n, double a, double *x)
{
	for (int i = 0; i < n; i++)
		x[i] *= a;
}

/*
 * The implementation of the solvers is inspired by
 * the lsolver library, by Christian Badura, available from:
 * http://www.mathematik.uni-freiburg.de
 * /IAM/Research/projectskr/lin_solver/
 *
 * About the Conjugate Gradient, details can be found in:
 *  Ashby, Manteuffel, Saylor
 *     A taxononmy for conjugate gradient methods
 *     SIAM J Numer Anal 27, 1542-1568 (1990)
 *
 *  This version is completely abstract, the same code can be used for
 * CPU/GPU, dense matrix / sparse matrix etc...
 *  Abstraction is realized through:
  *   - Abstract matrix interface (NLMatrix), that can implement different
 *     versions of matrix x vector product (CPU/GPU, sparse/dense ...)
 */

static uint32_t nlSolveSystem_PRE_CG(NLContext *context, NLMatrix M, NLMatrix P, double* b, double* x, double eps, uint32_t max_iter, double *sq_bnorm, double *sq_rnorm)
{
	int     N = (int)M->n;
	double* r = NL_NEW_VECTOR(N);
	double* d = NL_NEW_VECTOR(N);
	double* h = NL_NEW_VECTOR(N);
	double *Ad = h;
	uint32_t its = 0;
	double rh, alpha, beta;
	double b_square = ddot(N, b, b);
	double err = eps * eps*b_square;
	double curr_err;
	nlMultMatrixVector(M, x, r);
	daxpy(N, -1., b, r);
	nlMultMatrixVector(P, r, d);
	memcpy(h, d, N * sizeof(double));
	rh = ddot(N, r, h);
	curr_err = ddot(N, r, r);
	while (curr_err > err && its < max_iter) {
		if (context->progress_func != NULL)
			context->progress_func(its, max_iter, curr_err, err);
		nlMultMatrixVector(M, d, Ad);
		alpha = rh / ddot(N, d, Ad);
		daxpy(N, -alpha, d, x);
		daxpy(N, -alpha, Ad, r);
		nlMultMatrixVector(P, r, h);
		beta = 1. / rh;
		rh = ddot(N, r, h);
		beta *= rh;
		dscal(N, beta, d);
		daxpy(N, 1., h, d);
		++its;
		curr_err = ddot(N, r, r);
	}
	NL_DELETE_VECTOR(r);
	NL_DELETE_VECTOR(d);
	NL_DELETE_VECTOR(h);
	*sq_bnorm = b_square;
	*sq_rnorm = curr_err;
	return its;
}

uint32_t nlSolveSystemIterative(NLContext *context, NLMatrix M, NLMatrix P, double* b_in, double* x_in, double eps, uint32_t max_iter, uint32_t inner_iter)
{
	uint32_t N = M->n;
	uint32_t result = 0;
	double rnorm = 0.0;
	double bnorm = 0.0;
	double* b = b_in;
	double* x = x_in;
	assert(M->m == M->n);
	double sq_bnorm, sq_rnorm;
	result = nlSolveSystem_PRE_CG(context, M, P, b, x, eps, max_iter, &sq_bnorm, &sq_rnorm);
	/* Get residual norm and rhs norm */
	bnorm = sqrt(sq_bnorm);
	rnorm = sqrt(sq_rnorm);
	if (bnorm == 0.0)
		context->error = rnorm;
	else
		context->error = rnorm / bnorm;
	context->used_iterations = result;
	return result;
}

static bool nlSolveIterative(NLContext *context)
{
	double* b = context->b;
	double* x = context->x;
	uint32_t n = context->n;
	NLMatrix M = context->M;
	NLMatrix P = context->P;
	for (uint32_t k = 0; k < context->nb_systems; ++k) {
		nlSolveSystemIterative(context, M, P, b, x, context->threshold, context->max_iterations, context->inner_iterations);
		b += n;
		x += n;
	}
	return true;
}

struct NLJacobiPreconditioner
{
	uint32_t m;
	uint32_t n;
	uint32_t type;
	NLDestroyMatrixFunc destroy_func;
	NLMultMatrixVectorFunc mult_func;
	double* diag_inv;
};

static void nlJacobiPreconditionerDestroy(NLJacobiPreconditioner* M)
{
	NL_DELETE_ARRAY(M->diag_inv);
}

static void nlJacobiPreconditionerMult(NLJacobiPreconditioner* M, const double* x, double* y)
{
	for (uint32_t i = 0; i < M->n; ++i)
		y[i] = x[i] * M->diag_inv[i];
}

NLMatrix nlNewJacobiPreconditioner(NLMatrix M_in)
{
	NLSparseMatrix* M = NULL;
	NLJacobiPreconditioner* result = NULL;
	assert(M_in->type == NL_MATRIX_SPARSE_DYNAMIC);
	assert(M_in->m == M_in->n);
	M = (NLSparseMatrix*)M_in;
	result = NL_NEW(NLJacobiPreconditioner);
	result->m = M->m;
	result->n = M->n;
	result->type = NL_MATRIX_OTHER;
	result->destroy_func = (NLDestroyMatrixFunc)nlJacobiPreconditionerDestroy;
	result->mult_func = (NLMultMatrixVectorFunc)nlJacobiPreconditionerMult;
	result->diag_inv = NL_NEW_ARRAY(double, M->n);
	for (uint32_t i = 0; i < M->n; ++i)
		result->diag_inv[i] = (M->diag[i] == 0.0) ? 1.0 : 1.0 / M->diag[i];
	return (NLMatrix)result;
}

void nlSolverParameteri(NLContext *context, uint32_t pname, int param)
{
	if (pname == NL_NB_VARIABLES) {
		assert(param > 0);
		context->nb_variables = (uint32_t)param;
	} else if (pname == NL_MAX_ITERATIONS) {
		assert(param > 0);
		context->max_iterations = (uint32_t)param;
		context->max_iterations_defined = true;
	}
}

void nlSetFunction(NLContext *context, uint32_t pname, NLfunc param)
{
	if (pname == NL_FUNC_PROGRESS)
		context->progress_func = (NLProgressFunc)param;
}

void nlSetVariable(NLContext *context, uint32_t index, double value)
{
	assert(index >= 0 && index <= context->nb_variables - 1);
	NL_BUFFER_ITEM(context->variable_buffer[0], index) = value;
}

double nlGetVariable(NLContext *context, uint32_t index)
{
	assert(index >= 0 && index <= context->nb_variables - 1);
	return NL_BUFFER_ITEM(context->variable_buffer[0], index);
}

void nlLockVariable(NLContext *context, uint32_t index)
{
	assert(index >= 0 && index <= context->nb_variables - 1);
	context->variable_is_locked[index] = true;
}

static void nlVariablesToVector(NLContext *context)
{
	uint32_t n = context->n;
	assert(context->x != NULL);
	for (uint32_t k = 0; k < context->nb_systems; ++k) {
		for (uint32_t i = 0; i < context->nb_variables; ++i) {
			if (!context->variable_is_locked[i]) {
				uint32_t index = context->variable_index[i];
				assert(index < context->n);
				double value = NL_BUFFER_ITEM(context->variable_buffer[k], i);
				context->x[index + k * n] = value;
			}
		}
	}
}

static void nlVectorToVariables(NLContext *context)
{
	uint32_t n = context->n;
	assert(context->x != NULL);
	for (uint32_t k = 0; k < context->nb_systems; ++k) {
		for (uint32_t i = 0; i < context->nb_variables; ++i) {
			if (!context->variable_is_locked[i]) {
				uint32_t index = context->variable_index[i];
				assert(index < context->n);
				double value = context->x[index + k * n];
				NL_BUFFER_ITEM(context->variable_buffer[k], i) = value;
			}
		}
	}
}

void nlCoefficient(NLContext *context, uint32_t index, double value)
{
	assert(index >= 0 && index <= context->nb_variables - 1);
	if (context->variable_is_locked[index]) {
		/*
		 * Note: in al, indices are NLvariable indices,
		 * within [0..nb_variables-1]
		 */
		nlRowColumnAppend(&(context->al), index, value);
	} else {
		/*
		 * Note: in af, indices are system indices,
		 * within [0..n-1]
		 */
		nlRowColumnAppend(&(context->af), context->variable_index[index], value);
	}
}

void nlBegin(NLContext *context, uint32_t prim)
{
	if (prim == NL_SYSTEM) {
		assert(context->nb_variables > 0);
		context->variable_buffer = NL_NEW_ARRAY(NLBufferBinding, context->nb_systems);
		context->variable_value = NL_NEW_ARRAY(double, context->nb_variables * context->nb_systems);
		for (uint32_t k = 0; k < context->nb_systems; ++k) {
			context->variable_buffer[k].base_address =
				context->variable_value +
				k * context->nb_variables;
			context->variable_buffer[k].stride = sizeof(double);
		}
		context->variable_is_locked = NL_NEW_ARRAY(bool, context->nb_variables);
		context->variable_index = NL_NEW_ARRAY(uint32_t, context->nb_variables);
	} else if (prim == NL_MATRIX) {
		if (context->M != NULL)
			return;
		uint32_t n = 0;
		for (uint32_t i = 0; i < context->nb_variables; i++) {
			if (!context->variable_is_locked[i]) {
				context->variable_index[i] = n;
				n++;
			} else
				context->variable_index[i] = (uint32_t)~0;
		}
		context->n = n;
		if (!context->max_iterations_defined)
			context->max_iterations = n * 5;
		context->M = (NLMatrix)(NL_NEW(NLSparseMatrix));
		nlSparseMatrixConstruct((NLSparseMatrix*)(context->M), n, n);
		context->x = NL_NEW_ARRAY(double, n*context->nb_systems);
		context->b = NL_NEW_ARRAY(double, n*context->nb_systems);
		nlVariablesToVector(context);
		nlRowColumnConstruct(&context->af);
		nlRowColumnConstruct(&context->al);
		context->current_row = 0;
	} else if (prim == NL_ROW) {
		nlRowColumnZero(&context->af);
		nlRowColumnZero(&context->al);
	}
}

void nlEnd(NLContext *context, uint32_t prim)
{
	if (prim == NL_MATRIX) {
		nlRowColumnClear(&context->af);
		nlRowColumnClear(&context->al);
	} else if (prim == NL_ROW) {
		NLRowColumn*    af = &context->af;
		NLRowColumn*    al = &context->al;
		NLSparseMatrix* M = (NLSparseMatrix*)context->M;
		double* b = context->b;
		uint32_t nf = af->size;
		uint32_t nl = al->size;
		uint32_t n = context->n;
		uint32_t current_row = context->current_row;
		double S;
		/*
		 * least_squares : we want to solve
		 * A'A x = A'b
		 */
		for (uint32_t i = 0; i < nf; i++) {
			for (uint32_t j = 0; j < nf; j++) {
				nlSparseMatrixAdd(M, af->coeff[i].index, af->coeff[j].index, af->coeff[i].value * af->coeff[j].value);
			}
		}
		for (uint32_t k = 0; k < context->nb_systems; ++k) {
			S = 0.0;
			for (uint32_t jj = 0; jj < nl; ++jj) {
				uint32_t j = al->coeff[jj].index;
				S += al->coeff[jj].value * NL_BUFFER_ITEM(context->variable_buffer[k], j);
			}
			for (uint32_t jj = 0; jj < nf; jj++)
				b[k*n + af->coeff[jj].index] -= af->coeff[jj].value * S;
		}
		context->current_row++;
	}
}

bool nlSolve(NLContext *context)
{
	nlDeleteMatrix(context->P);
	context->P = nlNewJacobiPreconditioner(context->M);
	nlMatrixCompress(&context->M);
	bool result = nlSolveIterative(context);
	nlVectorToVariables(context);
	return result;
}
