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

void nlError(const char* function, const char* message) ;
void nlWarning(const char* function, const char* message) ;

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

extern NLprintfFunc nl_printf;

extern NLfprintfFunc nl_fprintf;

#define NL_NEW_VECTOR(dim) \
    (double*)malloc((size_t)(dim)*sizeof(double))

#define NL_DELETE_VECTOR(ptr) \
    free(ptr)

/* Abstract matrix interface */

struct NLMatrixStruct;
typedef struct NLMatrixStruct* NLMatrix;

typedef void(*NLDestroyMatrixFunc)(NLMatrix M);    

typedef void(*NLMultMatrixVectorFunc)(NLMatrix M, const double* x, double* y);

#define NL_MATRIX_SPARSE_DYNAMIC 0x1001
#define NL_MATRIX_CRS            0x1002
#define NL_MATRIX_OTHER          0x1006
    
struct NLMatrixStruct {
    NLuint m;

    NLuint n;

    NLenum type;

    NLDestroyMatrixFunc destroy_func;

    NLMultMatrixVectorFunc mult_func;
};

/* Dynamic arrays for sparse row/columns */

typedef struct  {
    NLuint index;

    NLdouble value; 
} NLCoeff;

typedef struct {
    NLuint size;
    
    NLuint capacity;

    NLCoeff* coeff;  
} NLRowColumn;

/* Compressed Row Storage */

typedef struct {
    NLuint m;
    
    NLuint n;

    NLenum type;
    
    NLDestroyMatrixFunc destroy_func;

    NLMultMatrixVectorFunc mult_func;
    
    NLdouble* val;    

    NLuint* rowptr;

    NLuint* colind;

    NLuint nslices;

    NLuint* sliceptr;
} NLCRSMatrix;

/* SparseMatrix data structure */

typedef struct {
    NLuint m;
    
    NLuint n;

    NLenum type;
    
    NLDestroyMatrixFunc destroy_func;

    NLMultMatrixVectorFunc mult_func;

    
    NLuint diag_size;

    NLuint diag_capacity;
    
    NLRowColumn* row;

    NLRowColumn* column;

    NLdouble*    diag;

    NLuint row_capacity;

    NLuint column_capacity;
    
} NLSparseMatrix;

/******* extracted from nl_context.h *******/

#ifndef OPENNL_CONTEXT_H
#define OPENNL_CONTEXT_H




/* NLContext data structure */

typedef void(*NLProgressFunc)(
    NLuint cur_iter, NLuint max_iter, double cur_err, double max_err
);

#define NL_STATE_INITIAL                0
#define NL_STATE_SYSTEM                 1
#define NL_STATE_MATRIX                 2
#define NL_STATE_ROW                    3
#define NL_STATE_MATRIX_CONSTRUCTED     4
#define NL_STATE_SYSTEM_CONSTRUCTED     5
#define NL_STATE_SOLVED                 6

typedef struct {
    void* base_address;
    NLuint stride;
} NLBufferBinding;

#define NL_BUFFER_ITEM(B,i) \
    *(double*)((void*)((char*)((B).base_address)+((i)*(B).stride)))


typedef struct {
    NLenum           state;

    NLBufferBinding* variable_buffer;
    
    NLdouble*        variable_value;

    NLboolean*       variable_is_locked;

    NLuint*          variable_index;
    
    NLuint           n;

    NLMatrix         M;

    NLMatrix         P;

    NLMatrix         B;
    
    NLRowColumn      af;

    NLRowColumn      al;

    NLdouble*        x;

    NLdouble*        b;

    NLenum           preconditioner;

    NLboolean        preconditioner_defined;
    
    NLuint           nb_variables;

    NLuint           nb_systems;

    NLuint           current_row;

    NLuint           max_iterations;


    NLboolean        max_iterations_defined;
    
    NLuint           inner_iterations;

    NLdouble         threshold;

    NLboolean        threshold_defined;
    
    NLdouble         omega;

    NLuint           used_iterations;

    NLdouble         error;

    NLProgressFunc   progress_func;
   
} NLContextStruct;

extern NLContextStruct* nlCurrentContext;

void nlCheckState(NLenum state);

void nlTransition(NLenum from_state, NLenum to_state);

#endif

/******* extracted from nl_os.c *******/


#if (defined (WIN32) || defined(_WIN64))
#include <windows.h>
#else
#include <sys/types.h>
#include <sys/times.h> 
#endif

/* Error-reporting functions */

NLprintfFunc nl_printf = printf;
NLfprintfFunc nl_fprintf = fprintf;

void nlError(const char* function, const char* message) {
    nl_fprintf(stderr, "OpenNL error in %s(): %s\n", function, message) ; 
}

void nlWarning(const char* function, const char* message) {
    nl_fprintf(stderr, "OpenNL warning in %s(): %s\n", function, message) ; 
}

void nlPrintfFuncs(NLprintfFunc f1, NLfprintfFunc f2) {
    nl_printf = f1;
    nl_fprintf = f2;
}

void nlDeleteMatrix(NLMatrix M) {
    if(M == NULL) {
        return;
    }
    M->destroy_func(M);
    NL_DELETE(M);
}

void nlMultMatrixVector(
    NLMatrix M, const double* x, double* y
) {
    M->mult_func(M,x,y);
}

void nlRowColumnConstruct(NLRowColumn* c) {
    c->size     = 0;
    c->capacity = 0;
    c->coeff    = NULL;
}

void nlRowColumnDestroy(NLRowColumn* c) {
    NL_DELETE_ARRAY(c->coeff);
    c->size = 0;
    c->capacity = 0;
}

void nlRowColumnGrow(NLRowColumn* c) {
    if(c->capacity != 0) {
        c->capacity = 2 * c->capacity;
        c->coeff = NL_RENEW_ARRAY(NLCoeff, c->coeff, c->capacity);
    } else {
        c->capacity = 4;
        c->coeff = NL_NEW_ARRAY(NLCoeff, c->capacity);
    }
}

void nlRowColumnAdd(NLRowColumn* c, NLuint index, NLdouble value) {
    NLuint i;
    for(i=0; i<c->size; i++) {
        if(c->coeff[i].index == index) {
            c->coeff[i].value += value;
            return;
        }
    }
    if(c->size == c->capacity) {
        nlRowColumnGrow(c);
    }
    c->coeff[c->size].index = index;
    c->coeff[c->size].value = value;
    c->size++;
}

/* Does not check whether the index already exists */
void nlRowColumnAppend(NLRowColumn* c, NLuint index, NLdouble value) {
    if(c->size == c->capacity) {
        nlRowColumnGrow(c);
    }
    c->coeff[c->size].index = index;
    c->coeff[c->size].value = value;
    c->size++;
}

void nlRowColumnZero(NLRowColumn* c) {
    c->size = 0;
}

void nlRowColumnClear(NLRowColumn* c) {
    c->size     = 0;
    c->capacity = 0;
    NL_DELETE_ARRAY(c->coeff);
}

static int nlCoeffCompare(const void* p1, const void* p2) {
    return (((NLCoeff*)(p2))->index < ((NLCoeff*)(p1))->index);
}

void nlRowColumnSort(NLRowColumn* c) {
    qsort(c->coeff, c->size, sizeof(NLCoeff), nlCoeffCompare);
}


/* CRSMatrix data structure */

static void nlCRSMatrixDestroy(NLCRSMatrix* M) {
    NL_DELETE_ARRAY(M->val);
    NL_DELETE_ARRAY(M->rowptr);
    NL_DELETE_ARRAY(M->colind);
    NL_DELETE_ARRAY(M->sliceptr);
    M->m = 0;
    M->n = 0;
    M->nslices = 0;
}

static void nlCRSMatrixMultSlice(
    NLCRSMatrix* M, const double* x, double* y, NLuint Ibegin, NLuint Iend
) {
    NLuint i,j;
    for(i=Ibegin; i<Iend; ++i) {
        double sum=0.0;
        for(j=M->rowptr[i]; j<M->rowptr[i+1]; ++j) {
            sum += M->val[j] * x[M->colind[j]];
        }
        y[i] = sum; 
    }
}

static void nlCRSMatrixMult(
    NLCRSMatrix* M, const double* x, double* y
) {
    int slice;
    int nslices = (int)(M->nslices);
#if defined(_OPENMP)
#pragma omp parallel for private(slice)
#endif
	for(slice=0; slice<nslices; ++slice) {
	    nlCRSMatrixMultSlice(
		M,x,y,M->sliceptr[slice],M->sliceptr[slice+1]
	    );
	}
}

void nlCRSMatrixConstruct(
    NLCRSMatrix* M, NLuint m, NLuint n, NLuint nnz, NLuint nslices
) {
    M->m = m;
    M->n = n;
    M->type = NL_MATRIX_CRS;
    M->destroy_func = (NLDestroyMatrixFunc)nlCRSMatrixDestroy;
	M->mult_func = (NLMultMatrixVectorFunc)nlCRSMatrixMult;
    M->nslices = nslices;
    M->val = NL_NEW_ARRAY(double, nnz);
    M->rowptr = NL_NEW_ARRAY(NLuint, m+1);
    M->colind = NL_NEW_ARRAY(NLuint, nnz);
    M->sliceptr = NL_NEW_ARRAY(NLuint, nslices+1);
}

/* SparseMatrix data structure */


static void nlSparseMatrixDestroyRowColumns(NLSparseMatrix* M) {
    NLuint i;
    for(i=0; i<M->m; i++) {
        nlRowColumnDestroy(&(M->row[i]));
    }
    NL_DELETE_ARRAY(M->row);
}

void nlSparseMatrixDestroy(NLSparseMatrix* M) {
    assert(M->type == NL_MATRIX_SPARSE_DYNAMIC);
    nlSparseMatrixDestroyRowColumns(M);
    NL_DELETE_ARRAY(M->diag);
}

void nlSparseMatrixAdd(NLSparseMatrix* M, NLuint i, NLuint j, NLdouble value) {
    assert(i >= 0 && i <= M->m - 1);
    assert(j >= 0 && j <= M->n - 1);
    if(i == j) {
        M->diag[i] += value;
    }
    nlRowColumnAdd(&(M->row[i]), j, value);
}

static void nlSparseMatrixAddSparseMatrix(
    NLSparseMatrix* M, double mul, const NLSparseMatrix* N    
) {
    NLuint i,jj;
    assert(M->m == N->m);
    assert(M->n == N->n);
	for(i=0; i<N->m; ++i) {
	    for(jj=0; jj<N->row[i].size; ++jj) {
		nlSparseMatrixAdd(
		    M,
		    i, N->row[i].coeff[jj].index,
		    mul*N->row[i].coeff[jj].value
		);
	    }
	}
}

static void nlSparseMatrixAddCRSMatrix(
    NLSparseMatrix* M, double mul, const NLCRSMatrix* N    
) {
    NLuint i,jj;
    assert(M->m == N->m);
    assert(M->n == N->n);
    for(i=0; i<M->m; ++i) {
	for(jj=N->rowptr[i]; jj<N->rowptr[i+1]; ++jj) {
	    nlSparseMatrixAdd(
		M,
		i,
		N->colind[jj],
		mul*N->val[jj]
	    );
	}
    }
}

void nlSparseMatrixAddMatrix(
    NLSparseMatrix* M, double mul, const NLMatrix N
) {
    assert(M->m == N->m);
    assert(M->n == N->n);
    if(N->type == NL_MATRIX_SPARSE_DYNAMIC) {
	nlSparseMatrixAddSparseMatrix(M, mul, (const NLSparseMatrix*)N);
    } else if(N->type == NL_MATRIX_CRS) {
	nlSparseMatrixAddCRSMatrix(M, mul, (const NLCRSMatrix*)N);	
    } else {
	assert(0);
    }
}
    


void nlSparseMatrixZero( NLSparseMatrix* M) {
    NLuint i;
    for(i=0; i<M->m; i++) {
        nlRowColumnZero(&(M->row[i]));
    }
    NL_CLEAR_ARRAY(NLdouble, M->diag, M->diag_size);
}

void nlSparseMatrixClear( NLSparseMatrix* M) {
    NLuint i;
    for(i=0; i<M->m; i++) {
        nlRowColumnClear(&(M->row[i]));
    }
    NL_CLEAR_ARRAY(NLdouble, M->diag, M->diag_size);
}

/* Returns the number of non-zero coefficients */
NLuint nlSparseMatrixNNZ( NLSparseMatrix* M) {
    NLuint nnz = 0;
    NLuint i;
    for(i = 0; i<M->m; i++) {
        nnz += M->row[i].size;
    }
    return nnz;
}

void nlSparseMatrixSort( NLSparseMatrix* M) {
    NLuint i;
    for(i = 0; i<M->m; i++) {
        nlRowColumnSort(&(M->row[i]));                
    }
}

void nlSparseMatrixMAddRow(
    NLSparseMatrix* M, NLuint i1, double s, NLuint i2
) {
    NLuint jj;
    NLRowColumn* Ri2 = &(M->row[i2]);
    NLCoeff* c = NULL;

    assert(i1 < M->m);
    assert(i2 < M->m);
    
    for(jj=0; jj<Ri2->size; ++jj) {
	c = &(Ri2->coeff[jj]);
	nlSparseMatrixAdd(M, i1, c->index, s*c->value);
    }
}

void nlSparseMatrixScaleRow(
    NLSparseMatrix* M, NLuint i, double s
) {
    NLuint jj;
    NLRowColumn* Ri = &(M->row[i]);
    NLCoeff* c = NULL;

    assert(i < M->m);
    
    for(jj=0; jj<Ri->size; ++jj) {
	c = &(Ri->coeff[jj]);
	c->value *= s;
    }
    if(i < M->diag_size) {
	M->diag[i] *= s;
    }
}

void nlSparseMatrixZeroRow(
    NLSparseMatrix* M, NLuint i
) {
    NLRowColumn* Ri = &(M->row[i]);

    assert(i < M->m);
    
    Ri->size = 0;
    if(i < M->diag_size) {
	M->diag[i] = 0.0;
    }
}



/* SparseMatrix x Vector routines, internal helper routines */

static void nlSparseMatrix_mult_rows(
        NLSparseMatrix* A,
        const NLdouble* x,
        NLdouble* y
) {
    /* 
     * Note: OpenMP does not like unsigned ints
     * (causes some floating point exceptions),
     * therefore I use here signed ints for all
     * indices.
     */
    
    int m = (int)(A->m);
    int i,ij;
    NLCoeff* c = NULL;
    NLRowColumn* Ri = NULL;

#if defined(_OPENMP)    
#pragma omp parallel for private(i,ij,c,Ri)
#endif
    
    for(i=0; i<m; i++) {
        Ri = &(A->row[i]);       
        y[i] = 0;
        for(ij=0; ij<(int)(Ri->size); ij++) {
            c = &(Ri->coeff[ij]);
            y[i] += c->value * x[c->index];
        }
    }
}

void nlSparseMatrixMult(
    NLSparseMatrix* A, const NLdouble* x, NLdouble* y
) {
    assert(A->type == NL_MATRIX_SPARSE_DYNAMIC);
    nlSparseMatrix_mult_rows(A, x, y);
}

void nlSparseMatrixConstruct(
    NLSparseMatrix* M, NLuint m, NLuint n
) {
    NLuint i;
    M->m = m;
    M->n = n;
    M->type = NL_MATRIX_SPARSE_DYNAMIC;
    M->destroy_func = (NLDestroyMatrixFunc)nlSparseMatrixDestroy;
    M->mult_func = (NLMultMatrixVectorFunc)nlSparseMatrixMult;
    M->row = NL_NEW_ARRAY(NLRowColumn, m);
	M->row_capacity = m;
    for(i=0; i<n; i++) {
        nlRowColumnConstruct(&(M->row[i]));
    }
	M->row_capacity = 0;
    M->column = NULL;
	M->column_capacity = 0;
    M->diag_size = MIN(m,n);
    M->diag_capacity = M->diag_size;
    M->diag = NL_NEW_ARRAY(NLdouble, M->diag_size);
}

NLMatrix nlSparseMatrixNew(
    NLuint m, NLuint n
) {
    NLSparseMatrix* result = NL_NEW(NLSparseMatrix);
    nlSparseMatrixConstruct(result, m, n);
    return (NLMatrix)result;
}

static void adjust_diag(NLSparseMatrix* M) {
    NLuint new_diag_size = MIN(M->m, M->n);
    NLuint i;
    if(new_diag_size > M->diag_size) {
	if(new_diag_size > M->diag_capacity) {
	    M->diag_capacity *= 2;
	    if(M->diag_capacity == 0) {
		M->diag_capacity = 16;
	    }
	    M->diag = NL_RENEW_ARRAY(double, M->diag, M->diag_capacity);
	    for(i=M->diag_size; i<new_diag_size; ++i) {
		M->diag[i] = 0.0;
	    }
	}
	M->diag_size= new_diag_size;
    }
}

void nlSparseMatrixAddRow( NLSparseMatrix* M) {
    ++M->m;
	if(M->m > M->row_capacity) {
	    M->row_capacity *= 2;
	    if(M->row_capacity == 0) {
		M->row_capacity = 16;
	    }
	    M->row = NL_RENEW_ARRAY(
		NLRowColumn, M->row, M->row_capacity
	    );
	}
	nlRowColumnConstruct(&(M->row[M->m-1]));
    adjust_diag(M);
}

void nlSparseMatrixAddColumn( NLSparseMatrix* M) {
    ++M->n;
    adjust_diag(M);
}




NLMatrix nlCRSMatrixNewFromSparseMatrix(NLSparseMatrix* M) {
    NLuint nnz = nlSparseMatrixNNZ(M);
    NLuint nslices = 8; /* TODO: get number of cores */
    NLuint slice, cur_bound, cur_NNZ, cur_row;
    NLuint i,ij,k; 
    NLuint slice_size = nnz / nslices;
    NLCRSMatrix* CRS = NL_NEW(NLCRSMatrix);

    nlCRSMatrixConstruct(CRS, M->m, M->n, nnz, nslices);
    nlSparseMatrixSort(M);
    /* Convert matrix to CRS format */
    k=0;
    for(i=0; i<M->m; ++i) {
        NLRowColumn* Ri = &(M->row[i]);
        CRS->rowptr[i] = k;
        for(ij=0; ij<Ri->size; ij++) {
            NLCoeff* c = &(Ri->coeff[ij]);
            CRS->val[k] = c->value;
            CRS->colind[k] = c->index;
            ++k;
        }
    }
    CRS->rowptr[M->m] = k;
        
    /* Create "slices" to be used by parallel sparse matrix vector product */
    if(CRS->sliceptr != NULL) {
	cur_bound = slice_size;
	cur_NNZ = 0;
	cur_row = 0;
	CRS->sliceptr[0]=0;
	for(slice=1; slice<nslices; ++slice) {
	    while(cur_NNZ < cur_bound && cur_row < M->m) {
		++cur_row;
		cur_NNZ += CRS->rowptr[cur_row+1] - CRS->rowptr[cur_row];
	    }
	    CRS->sliceptr[slice] = cur_row;
	    cur_bound += slice_size;
	}
	CRS->sliceptr[nslices]=M->m;
    }
    return (NLMatrix)CRS;
}

void nlMatrixCompress(NLMatrix* M) {
    NLMatrix CRS = NULL;
    if((*M)->type != NL_MATRIX_SPARSE_DYNAMIC) {
        return;
    }
    CRS = nlCRSMatrixNewFromSparseMatrix((NLSparseMatrix*)*M);
    nlDeleteMatrix(*M);
    *M = CRS;
}

/******* extracted from nl_context.c *******/


NLContextStruct* nlCurrentContext = NULL;

NLContext nlNewContext() {
    NLContextStruct* result     = NL_NEW(NLContextStruct);
    result->state               = NL_STATE_INITIAL;
    result->max_iterations      = 100;
    result->threshold           = 1e-6;
    result->omega               = 1.5;
    result->inner_iterations    = 5;
    result->progress_func       = NULL;
    result->nb_systems          = 1;
    nlMakeCurrent(result);
    return result;
}

void nlDeleteContext(NLContext context_in) {
    NLContextStruct* context = (NLContextStruct*)(context_in);
    if(nlCurrentContext == context) {
        nlCurrentContext = NULL;
    }

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

void nlMakeCurrent(NLContext context) {
    nlCurrentContext = (NLContextStruct*)(context);
}

NLContext nlGetCurrent() {
    return nlCurrentContext;
}


/* Finite state automaton   */

void nlCheckState(NLenum state) {
    assert(nlCurrentContext->state == state);
}

void nlTransition(NLenum from_state, NLenum to_state) {
    nlCheckState(from_state);
    nlCurrentContext->state = to_state;
}

static double ddot(int n, const double *x, const double *y)
{
	double sum = 0.0;
	for (int i = 0; i < n; i++)
		sum += x[i] * y[i];
	return sum;
}

static void daxpy(int n, double a, const double *x, double *y) {
	for (int i = 0; i < n; i++)
		y[i] = a * x[i] + y[i];
}

static void dscal(int n, double a, double *x) {
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

static NLuint nlSolveSystem_PRE_CG(
    NLMatrix M, NLMatrix P, NLdouble* b, NLdouble* x,
    double eps, NLuint max_iter, double *sq_bnorm, double *sq_rnorm
) {
    NLint     N        = (NLint)M->n;
    NLdouble* r = NL_NEW_VECTOR(N);
    NLdouble* d = NL_NEW_VECTOR(N);
    NLdouble* h = NL_NEW_VECTOR(N);
    NLdouble *Ad = h;
    NLuint its=0;
    NLdouble rh, alpha, beta;
    NLdouble b_square = ddot(N,b,b);
    NLdouble err=eps*eps*b_square;
    NLdouble curr_err;

    nlMultMatrixVector(M,x,r);
    daxpy(N,-1.,b,r);
    nlMultMatrixVector(P,r,d);
	memcpy(h, d, N * sizeof(NLdouble));
    rh=ddot(N,r,h);
    curr_err = ddot(N,r,r);

    while ( curr_err >err && its < max_iter) {
	if(nlCurrentContext != NULL) {
	    if(nlCurrentContext->progress_func != NULL) {
		nlCurrentContext->progress_func(its, max_iter, curr_err, err);
	    }
	}
	nlMultMatrixVector(M,d,Ad);
        alpha=rh/ddot(N,d,Ad);
        daxpy(N,-alpha,d,x);
        daxpy(N,-alpha,Ad,r);
	nlMultMatrixVector(P,r,h);
        beta=1./rh;
	rh=ddot(N,r,h);
	beta*=rh;
        dscal(N,beta,d);
        daxpy(N,1.,h,d);
        ++its;
        curr_err = ddot(N,r,r);
    }
    NL_DELETE_VECTOR(r);
    NL_DELETE_VECTOR(d);
    NL_DELETE_VECTOR(h);
    *sq_bnorm = b_square;
    *sq_rnorm = curr_err;
    return its;
}

/* Main driver routine */

NLuint nlSolveSystemIterative(
    NLMatrix M, NLMatrix P, NLdouble* b_in, NLdouble* x_in,
    double eps, NLuint max_iter, NLuint inner_iter
) {
    NLuint N = M->n;
    NLuint result=0;
    NLdouble rnorm=0.0;
    NLdouble bnorm=0.0; 
    double* b = b_in;
    double* x = x_in;
    assert(M->m == M->n);

	double sq_bnorm, sq_rnorm;
	result = nlSolveSystem_PRE_CG(M,P,b,x,eps,max_iter, &sq_bnorm, &sq_rnorm);

    /* Get residual norm and rhs norm from BLAS context */
    if(nlCurrentContext != NULL) {
	bnorm = sqrt(sq_bnorm);
	rnorm = sqrt(sq_rnorm);
	if(bnorm == 0.0) {
	    nlCurrentContext->error = rnorm;
	} else {
	    nlCurrentContext->error = rnorm/bnorm;
	}
    }
    nlCurrentContext->used_iterations = result;
    return result;
}

static NLboolean nlSolveIterative() {
    NLdouble* b = nlCurrentContext->b;
    NLdouble* x = nlCurrentContext->x;
    NLuint n = nlCurrentContext->n;
    NLuint k;
    NLMatrix M = nlCurrentContext->M;
    NLMatrix P = nlCurrentContext->P;
    
    for(k=0; k<nlCurrentContext->nb_systems; ++k) {
	nlSolveSystemIterative(
	    M,
	    P,
	    b,
	    x,
	    nlCurrentContext->threshold,
	    nlCurrentContext->max_iterations,
	    nlCurrentContext->inner_iterations
	);
	b += n;
	x += n;
    }
    return NL_TRUE;
}

/******* extracted from nl_preconditioners.c *******/




typedef struct {
    NLuint m;

    NLuint n;

    NLenum type;

    NLDestroyMatrixFunc destroy_func;

    NLMultMatrixVectorFunc mult_func;
    
    NLdouble* diag_inv;
    
} NLJacobiPreconditioner;


static void nlJacobiPreconditionerDestroy(NLJacobiPreconditioner* M) {
    NL_DELETE_ARRAY(M->diag_inv);
}

static void nlJacobiPreconditionerMult(
    NLJacobiPreconditioner* M, const double* x, double* y
) {
    NLuint i;
    for(i=0; i<M->n; ++i) {
	y[i] = x[i] * M->diag_inv[i];
    }
}

NLMatrix nlNewJacobiPreconditioner(NLMatrix M_in) {
    NLSparseMatrix* M = NULL;
    NLJacobiPreconditioner* result = NULL;
    NLuint i;
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
    for(i=0; i<M->n; ++i) {
	result->diag_inv[i] = (M->diag[i] == 0.0) ? 1.0 : 1.0/M->diag[i];
    }
    return (NLMatrix)result;
}

static void nlSetupPreconditioner() {
    nlDeleteMatrix(nlCurrentContext->P);
    nlCurrentContext->P = NULL;
    
    switch(nlCurrentContext->preconditioner) {
    case NL_PRECOND_NONE:
        break;
    case NL_PRECOND_JACOBI:
	nlCurrentContext->P = nlNewJacobiPreconditioner(nlCurrentContext->M);
        break;
    default:
        assert(0);
    }
    nlMatrixCompress(&nlCurrentContext->M);
}

/******* extracted from nl_api.c *******/

static NLSparseMatrix* nlGetCurrentSparseMatrix() {
    NLSparseMatrix* result = NULL;
	assert(nlCurrentContext->M != NULL);	    
	assert(nlCurrentContext->M->type == NL_MATRIX_SPARSE_DYNAMIC);
	return (NLSparseMatrix*)(nlCurrentContext->M);
}

/* Get/Set parameters */

void nlSolverParameterd(NLenum pname, NLdouble param) {
    nlCheckState(NL_STATE_INITIAL);
    switch(pname) {
    case NL_THRESHOLD: {
        assert(param >= 0);
        nlCurrentContext->threshold = (NLdouble)param;
        nlCurrentContext->threshold_defined = NL_TRUE;
    } break;
    case NL_OMEGA: {
        nlCurrentContext->omega = (NLdouble)param;
    } break;
    default: {
        nlError("nlSolverParameterd","Invalid parameter");
        assert(0);
    }
    }
}

void nlSolverParameteri(NLenum pname, NLint param) {
    nlCheckState(NL_STATE_INITIAL);
    switch(pname) {
    case NL_NB_VARIABLES: {
        assert(param > 0);
        nlCurrentContext->nb_variables = (NLuint)param;
    } break;
    case NL_NB_SYSTEMS: {
	assert(param > 0);
	nlCurrentContext->nb_systems = (NLuint)param;
    } break;
    case NL_MAX_ITERATIONS: {
        assert(param > 0);
        nlCurrentContext->max_iterations = (NLuint)param;
        nlCurrentContext->max_iterations_defined = NL_TRUE;
    } break;
    case NL_INNER_ITERATIONS: {
        assert(param > 0);
        nlCurrentContext->inner_iterations = (NLuint)param;
    } break;
    case NL_PRECONDITIONER: {
        nlCurrentContext->preconditioner = (NLuint)param;
        nlCurrentContext->preconditioner_defined = NL_TRUE;
    } break;
    default: {
        nlError("nlSolverParameteri","Invalid parameter");
        assert(0);
    }
    }
}

void nlGetDoublev(NLenum pname, NLdouble* params) {
    switch(pname) {
    case NL_THRESHOLD: {
        *params = nlCurrentContext->threshold;
    } break;
    case NL_OMEGA: {
        *params = nlCurrentContext->omega;
    } break;
    case NL_ERROR: {
        *params = nlCurrentContext->error;
    } break;
    default: {
        nlError("nlGetDoublev","Invalid parameter");
        assert(0);
    } 
    }
}

void nlGetIntegerv(NLenum pname, NLint* params) {
    switch(pname) {
    case NL_NB_VARIABLES: {
        *params = (NLint)(nlCurrentContext->nb_variables);
    } break;
    case NL_NB_SYSTEMS: {
	*params = (NLint)(nlCurrentContext->nb_systems);
    } break;
    case NL_MAX_ITERATIONS: {
        *params = (NLint)(nlCurrentContext->max_iterations);
    } break;
    case NL_USED_ITERATIONS: {
        *params = (NLint)(nlCurrentContext->used_iterations);
    } break;
    case NL_PRECONDITIONER: {
        *params = (NLint)(nlCurrentContext->preconditioner);        
    } break;
    default: {
        nlError("nlGetIntegerv","Invalid parameter");
        assert(0);
    } 
    }
}

/* NL functions */

void  nlSetFunction(NLenum pname, NLfunc param) {
    switch(pname) {
    case NL_FUNC_PROGRESS:
        nlCurrentContext->progress_func = (NLProgressFunc)(param);
        break;
    default:
        nlError("nlSetFunction","Invalid parameter");        
        assert(0);
    }
}

/* Get/Set Lock/Unlock variables */

void nlSetVariable(NLuint index, NLdouble value) {
    nlCheckState(NL_STATE_SYSTEM);
    assert(index >= 0 && index <= nlCurrentContext->nb_variables - 1);
    NL_BUFFER_ITEM(nlCurrentContext->variable_buffer[0],index) = value;
}

NLdouble nlGetVariable(NLuint index) {
    assert(nlCurrentContext->state != NL_STATE_INITIAL);
    assert(index >= 0 && index <= nlCurrentContext->nb_variables - 1);
    return NL_BUFFER_ITEM(nlCurrentContext->variable_buffer[0],index);
}

void nlLockVariable(NLuint index) {
    nlCheckState(NL_STATE_SYSTEM);
    assert(index >= 0 && index <= nlCurrentContext->nb_variables - 1);
    nlCurrentContext->variable_is_locked[index] = NL_TRUE;
}

/* System construction */

static void nlVariablesToVector() {
    NLuint n=nlCurrentContext->n;
    NLuint k,i,index;
    NLdouble value;
    
    assert(nlCurrentContext->x != NULL);
    for(k=0; k<nlCurrentContext->nb_systems; ++k) {
	for(i=0; i<nlCurrentContext->nb_variables; ++i) {
	    if(!nlCurrentContext->variable_is_locked[i]) {
		index = nlCurrentContext->variable_index[i];
		assert(index < nlCurrentContext->n);		
		value = NL_BUFFER_ITEM(nlCurrentContext->variable_buffer[k],i);
		nlCurrentContext->x[index+k*n] = value;
	    }
	}
    }
}

static void nlVectorToVariables() {
    NLuint n=nlCurrentContext->n;
    NLuint k,i,index;
    NLdouble value;

    assert(nlCurrentContext->x != NULL);
    for(k=0; k<nlCurrentContext->nb_systems; ++k) {
	for(i=0; i<nlCurrentContext->nb_variables; ++i) {
	    if(!nlCurrentContext->variable_is_locked[i]) {
		index = nlCurrentContext->variable_index[i];
		assert(index < nlCurrentContext->n);
		value = nlCurrentContext->x[index+k*n];
		NL_BUFFER_ITEM(nlCurrentContext->variable_buffer[k],i) = value;
	    }
	}
    }
}


static void nlBeginSystem() {
    NLuint k;
    
    nlTransition(NL_STATE_INITIAL, NL_STATE_SYSTEM);
    assert(nlCurrentContext->nb_variables > 0);

    nlCurrentContext->variable_buffer = NL_NEW_ARRAY(
	NLBufferBinding, nlCurrentContext->nb_systems
    );
    
	nlCurrentContext->variable_value = NL_NEW_ARRAY(
	    NLdouble,
	    nlCurrentContext->nb_variables * nlCurrentContext->nb_systems
	);
	for(k=0; k<nlCurrentContext->nb_systems; ++k) {
	    nlCurrentContext->variable_buffer[k].base_address =
		nlCurrentContext->variable_value +
		k * nlCurrentContext->nb_variables;
	    nlCurrentContext->variable_buffer[k].stride = sizeof(NLdouble);
	}
    
    nlCurrentContext->variable_is_locked = NL_NEW_ARRAY(
	NLboolean, nlCurrentContext->nb_variables
    );
    nlCurrentContext->variable_index = NL_NEW_ARRAY(
	NLuint, nlCurrentContext->nb_variables
    );
}

static void nlEndSystem() {
    nlTransition(NL_STATE_MATRIX_CONSTRUCTED, NL_STATE_SYSTEM_CONSTRUCTED);    
}

static void nlInitializeM() {
    NLuint i;
    NLuint n = 0;

    for(i=0; i<nlCurrentContext->nb_variables; i++) {
        if(!nlCurrentContext->variable_is_locked[i]) {
            nlCurrentContext->variable_index[i] = n;
            n++;
        } else {
            nlCurrentContext->variable_index[i] = (NLuint)~0;
        }
    }

    nlCurrentContext->n = n;
    if(!nlCurrentContext->preconditioner_defined) {
        nlCurrentContext->preconditioner = NL_PRECOND_JACOBI;
    }
    if(!nlCurrentContext->max_iterations_defined) {
        nlCurrentContext->max_iterations = n*5;
    }
    if(!nlCurrentContext->threshold_defined) {
        nlCurrentContext->threshold = 1e-6;
    }

    nlCurrentContext->M = (NLMatrix)(NL_NEW(NLSparseMatrix));
    nlSparseMatrixConstruct(
	     (NLSparseMatrix*)(nlCurrentContext->M), n, n
    );

    nlCurrentContext->x = NL_NEW_ARRAY(
	NLdouble, n*nlCurrentContext->nb_systems
    );
    nlCurrentContext->b = NL_NEW_ARRAY(
	NLdouble, n*nlCurrentContext->nb_systems
    );

    nlVariablesToVector();

    nlRowColumnConstruct(&nlCurrentContext->af);
    nlRowColumnConstruct(&nlCurrentContext->al);

    nlCurrentContext->current_row = 0;
}

static void nlEndMatrix() {
    nlTransition(NL_STATE_MATRIX, NL_STATE_MATRIX_CONSTRUCTED);    
    nlRowColumnClear(&nlCurrentContext->af);
    nlRowColumnClear(&nlCurrentContext->al);
}

static void nlBeginRow() {
    nlTransition(NL_STATE_MATRIX, NL_STATE_ROW);
    nlRowColumnZero(&nlCurrentContext->af);
    nlRowColumnZero(&nlCurrentContext->al);
}

static void nlEndRow() {
    NLRowColumn*    af = &nlCurrentContext->af;
    NLRowColumn*    al = &nlCurrentContext->al;
    NLSparseMatrix* M  = nlGetCurrentSparseMatrix();
    NLdouble* b        = nlCurrentContext->b;
    NLuint nf          = af->size;
    NLuint nl          = al->size;
    NLuint n           = nlCurrentContext->n;
    NLuint current_row = nlCurrentContext->current_row;
    NLuint i,j,jj;
    NLdouble S;
    NLuint k;
    nlTransition(NL_STATE_ROW, NL_STATE_MATRIX);

    /*
     * least_squares : we want to solve
     * A'A x = A'b
     */
    for(i=0; i<nf; i++) {
        for(j=0; j<nf; j++) {
            nlSparseMatrixAdd(
                M, af->coeff[i].index, af->coeff[j].index,
                af->coeff[i].value * af->coeff[j].value
            );
        }
    }
	for(k=0; k<nlCurrentContext->nb_systems; ++k) {
	    S = 0.0;
	    for(jj=0; jj<nl; ++jj) {
		j = al->coeff[jj].index;
		S += al->coeff[jj].value *
		    NL_BUFFER_ITEM(nlCurrentContext->variable_buffer[k],j);
	    }
	    for(jj=0; jj<nf; jj++) {
		b[ k*n+af->coeff[jj].index ] -= af->coeff[jj].value * S;
	    }
	}
    nlCurrentContext->current_row++;
    for(k=0; k<nlCurrentContext->nb_systems; ++k) {
    }
}

void nlCoefficient(NLuint index, NLdouble value) {
    nlCheckState(NL_STATE_ROW);
    assert(index >= 0 && index <= nlCurrentContext->nb_variables - 1);
    if(nlCurrentContext->variable_is_locked[index]) {
	/* 
	 * Note: in al, indices are NLvariable indices, 
	 * within [0..nb_variables-1]
	 */
        nlRowColumnAppend(&(nlCurrentContext->al), index, value);
    } else {
	/*
	 * Note: in af, indices are system indices, 
	 * within [0..n-1]
	 */
        nlRowColumnAppend(
	    &(nlCurrentContext->af),
	    nlCurrentContext->variable_index[index], value
	);
    }
}

void nlBegin(NLenum prim) {
    switch(prim) {
    case NL_SYSTEM: {
        nlBeginSystem();
    } break;
    case NL_MATRIX: {
	nlTransition(NL_STATE_SYSTEM, NL_STATE_MATRIX);
	if(
	    nlCurrentContext->M == NULL
	) {
	    nlInitializeM();
	}
    } break;
    case NL_ROW: {
        nlBeginRow();
    } break;
    default: {
        assert(0);
    }
    }
}

void nlEnd(NLenum prim) {
    switch(prim) {
    case NL_SYSTEM: {
        nlEndSystem();
    } break;
    case NL_MATRIX: {
        nlEndMatrix();
    } break;
    case NL_ROW: {
        nlEndRow();
    } break;
    default: {
        assert(0);
    }
    }
}


/* nlSolve() driver routine */

NLboolean nlSolve() {
    NLboolean result;
    nlCheckState(NL_STATE_SYSTEM_CONSTRUCTED);
	nlSetupPreconditioner();
	result = nlSolveIterative();
    nlVectorToVariables();
    nlTransition(NL_STATE_SYSTEM_CONSTRUCTED, NL_STATE_SOLVED);
    return result;
}
