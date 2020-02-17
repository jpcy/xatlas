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



/******* extracted from nl_private.h *******/

#ifndef OPENNL_PRIVATE_H
#define OPENNL_PRIVATE_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>


#if defined(__APPLE__) && defined(__MACH__)
#define NL_OS_APPLE
#endif

#if defined(__linux__) || defined(__ANDROID__) || defined(NL_OS_APPLE)
#define NL_OS_UNIX
#endif


#if defined(WIN32) || defined(_WIN64)
#define NL_OS_WINDOWS
#endif

#define nl_arg_used(x) (void)x


#if defined(__clang__) || defined(__GNUC__)
#define NL_NORETURN __attribute__((noreturn))
#else
#define NL_NORETURN 
#endif

#if defined(_MSC_VER)
#define NL_NORETURN_DECL __declspec(noreturn) 
#else
#define NL_NORETURN_DECL 
#endif

NL_NORETURN_DECL void nl_assertion_failed(
    const char* cond, const char* file, int line
) NL_NORETURN;

NL_NORETURN_DECL void nl_range_assertion_failed(
    double x, double min_val, double max_val, const char* file, int line
) NL_NORETURN;

NL_NORETURN_DECL void nl_should_not_have_reached(
    const char* file, int line
) NL_NORETURN;

#define nl_assert(x) {                                          \
    if(!(x)) {                                                  \
        nl_assertion_failed(#x,__FILE__, __LINE__) ;            \
    }                                                           \
} 

#define nl_range_assert(x,min_val,max_val) {                    \
    if(((x) < (min_val)) || ((x) > (max_val))) {                \
        nl_range_assertion_failed(x, min_val, max_val,          \
            __FILE__, __LINE__                                  \
        ) ;                                                     \
    }                                                           \
}

#define nl_assert_not_reached {                                 \
    nl_should_not_have_reached(__FILE__, __LINE__) ;            \
}

#ifdef NL_DEBUG
    #define nl_debug_assert(x) nl_assert(x)
    #define nl_debug_range_assert(x,min_val,max_val)            \
                               nl_range_assert(x,min_val,max_val)
#else
    #define nl_debug_assert(x) 
    #define nl_debug_range_assert(x,min_val,max_val) 
#endif

#ifdef NL_PARANOID
    #define nl_parano_assert(x) nl_assert(x)
    #define nl_parano_range_assert(x,min_val,max_val)           \
                               nl_range_assert(x,min_val,max_val)
#else
    #define nl_parano_assert(x) 
    #define nl_parano_range_assert(x,min_val,max_val) 
#endif


void nlError(const char* function, const char* message) ;

void nlWarning(const char* function, const char* message) ;


NLdouble nlCurrentTime(void);

/* classic macros */

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

#define NL_CLEAR(T, x)           memset(x, 0, sizeof(T)) 

#define NL_CLEAR_ARRAY(T,x,NB)   memset(x, 0, (size_t)(NB)*sizeof(T)) 


#define NL_UINT_MAX 0xffffffff

#define NL_USHORT_MAX 0xffff


extern NLprintfFunc nl_printf;

extern NLfprintfFunc nl_fprintf;


#endif

/******* extracted from nl_blas.h *******/



#ifndef OPENNL_BLAS_H
#define OPENNL_BLAS_H

struct NLBlas;

typedef struct NLBlas* NLBlas_t;

typedef enum {
    NoTranspose=0, Transpose=1, ConjugateTranspose=2
} MatrixTranspose ;

typedef enum {
    UnitTriangular=0, NotUnitTriangular=1
} MatrixUnitTriangular ;

typedef enum {
    NL_HOST_MEMORY, NL_DEVICE_MEMORY
} NLmemoryType;

typedef void* (*FUNPTR_malloc)(
    NLBlas_t blas, NLmemoryType type, size_t size
);

typedef void (*FUNPTR_free)(
    NLBlas_t blas, NLmemoryType type, size_t size, void* ptr
);

typedef void (*FUNPTR_memcpy)(
    NLBlas_t blas,
    void* to, NLmemoryType to_type,
    void* from, NLmemoryType from_type,
    size_t size
);

typedef void (*FUNPTR_dscal)(
    NLBlas_t blas, int n, double a, double *x, int incx
);


typedef double (*FUNPTR_ddot)(
    NLBlas_t blas, int n, const double *x, int incx, const double *y, int incy
);

typedef void (*FUNPTR_daxpy)(
    NLBlas_t blas, int n,
    double a, const double *x, int incx, double *y, int incy
);

struct NLBlas {
    FUNPTR_malloc Malloc;
    FUNPTR_free Free;
    FUNPTR_memcpy Memcpy;

    FUNPTR_dscal Dscal;
    FUNPTR_ddot  Ddot;
    FUNPTR_daxpy Daxpy;

    NLboolean has_unified_memory;
    double start_time;
    NLulong flops;
    NLulong used_ram[2];
    NLulong max_used_ram[2];
    
    /* 
     * Used for stats of the linear solver
     * (a bit ugly, should not be here, but
     * more convenient for now...)
     */
    double sq_rnorm; 
    double sq_bnorm;
};

NLboolean nlBlasHasUnifiedMemory(NLBlas_t blas);

void nlBlasResetStats(NLBlas_t blas);

double nlBlasGFlops(NLBlas_t blas);

NLulong nlBlasUsedRam(NLBlas_t blas, NLmemoryType type);

NLulong nlBlasMaxUsedRam(NLBlas_t blas, NLmemoryType type);

NLBlas_t nlHostBlas(void);

#define NL_NEW_VECTOR(blas, memtype, dim) \
    (double*)blas->Malloc(blas,memtype,(size_t)(dim)*sizeof(double))

#define NL_DELETE_VECTOR(blas, memtype, dim, ptr) \
    blas->Free(blas,memtype,(size_t)(dim)*sizeof(double),ptr)



#endif

/******* extracted from nl_matrix.h *******/


#ifndef OPENNL_MATRIX_H
#define OPENNL_MATRIX_H


#ifdef __cplusplus
extern "C" {
#endif


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

NLAPI void NLAPIENTRY nlDeleteMatrix(NLMatrix M);

NLAPI void NLAPIENTRY nlMultMatrixVector(
    NLMatrix M, const double* x, double* y
);
    

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

NLAPI void NLAPIENTRY nlRowColumnConstruct(NLRowColumn* c);

NLAPI void NLAPIENTRY nlRowColumnDestroy(NLRowColumn* c);

NLAPI void NLAPIENTRY nlRowColumnGrow(NLRowColumn* c);

NLAPI void NLAPIENTRY nlRowColumnAdd(
    NLRowColumn* c, NLuint index, NLdouble value
);

NLAPI void NLAPIENTRY nlRowColumnAppend(
    NLRowColumn* c, NLuint index, NLdouble value
);

NLAPI void NLAPIENTRY nlRowColumnZero(NLRowColumn* c);

NLAPI void NLAPIENTRY nlRowColumnClear(NLRowColumn* c);

NLAPI void NLAPIENTRY nlRowColumnSort(NLRowColumn* c);


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

NLAPI void NLAPIENTRY nlCRSMatrixConstruct(
    NLCRSMatrix* M, NLuint m, NLuint n, NLuint nnz, NLuint nslices
);

NLAPI NLuint NLAPIENTRY nlCRSMatrixNNZ(NLCRSMatrix* M);
    

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


NLAPI NLMatrix NLAPIENTRY nlSparseMatrixNew(
    NLuint m, NLuint n
);

NLAPI void NLAPIENTRY nlSparseMatrixConstruct(
    NLSparseMatrix* M, NLuint m, NLuint n
);

NLAPI void NLAPIENTRY nlSparseMatrixDestroy(NLSparseMatrix* M);

NLAPI void NLAPIENTRY nlSparseMatrixMult(
    NLSparseMatrix* A, const NLdouble* x, NLdouble* y
);    
    
NLAPI void NLAPIENTRY nlSparseMatrixAdd(
    NLSparseMatrix* M, NLuint i, NLuint j, NLdouble value
);

NLAPI void NLAPIENTRY nlSparseMatrixAddMatrix(
    NLSparseMatrix* M, double mul, const NLMatrix N
);	
    
NLAPI void NLAPIENTRY nlSparseMatrixZero( NLSparseMatrix* M);

NLAPI void NLAPIENTRY nlSparseMatrixClear( NLSparseMatrix* M);

NLAPI NLuint NLAPIENTRY nlSparseMatrixNNZ( NLSparseMatrix* M);

NLAPI void NLAPIENTRY nlSparseMatrixSort( NLSparseMatrix* M);

NLAPI void NLAPIENTRY nlSparseMatrixAddRow( NLSparseMatrix* M);

NLAPI void NLAPIENTRY nlSparseMatrixAddColumn( NLSparseMatrix* M);

NLAPI void NLAPIENTRY nlSparseMatrixMAddRow(
    NLSparseMatrix* M, NLuint i1, double s, NLuint i2
);

NLAPI void NLAPIENTRY nlSparseMatrixScaleRow(
    NLSparseMatrix* M, NLuint i, double s
);

NLAPI void NLAPIENTRY nlSparseMatrixZeroRow(
    NLSparseMatrix* M, NLuint i
);

    


NLAPI NLMatrix NLAPIENTRY nlCRSMatrixNewFromSparseMatrix(NLSparseMatrix* M);    

NLAPI void NLAPIENTRY nlMatrixCompress(NLMatrix* M);

#ifdef __cplusplus
}
#endif

#endif

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


    NLdouble         start_time;
    
    NLdouble         elapsed_time;

    NLProgressFunc   progress_func;

    NLulong          flops;
   
} NLContextStruct;

extern NLContextStruct* nlCurrentContext;

void nlCheckState(NLenum state);

void nlTransition(NLenum from_state, NLenum to_state);

#endif

/******* extracted from nl_iterative_solvers.h *******/


#ifndef OPENNL_ITERATIVE_SOLVERS_H
#define OPENNL_ITERATIVE_SOLVERS_H



NLAPI NLuint NLAPIENTRY nlSolveSystemIterative(
    NLBlas_t blas,
    NLMatrix M, NLMatrix P, NLdouble* b, NLdouble* x,
    double eps, NLuint max_iter, NLuint inner_iter
);

#endif


/******* extracted from nl_preconditioners.h *******/

#ifndef OPENNL_PRECONDITIONERS_H
#define OPENNL_PRECONDITIONERS_H




/* preconditioners */

NLMatrix nlNewJacobiPreconditioner(NLMatrix M);

#endif

/******* extracted from nl_os.c *******/


#if (defined (WIN32) || defined(_WIN64))
#include <windows.h>
#else
#include <sys/types.h>
#include <sys/times.h> 
#endif

/* Assertions */


void nl_assertion_failed(const char* cond, const char* file, int line) {
    nl_fprintf(
        stderr, 
        "OpenNL assertion failed: %s, file:%s, line:%d\n",
        cond,file,line
    ) ;
    abort() ;
}

void nl_range_assertion_failed(
    double x, double min_val, double max_val, const char* file, int line
) {
    nl_fprintf(
        stderr, 
        "OpenNL range assertion failed: "
	"%f in [ %f ... %f ], file:%s, line:%d\n",
        x, min_val, max_val, file,line
    ) ;
    abort() ;
}

void nl_should_not_have_reached(const char* file, int line) {
    nl_fprintf(
        stderr, 
        "OpenNL should not have reached this point: file:%s, line:%d\n",
        file,line
    ) ;
    abort() ;
}



/* Timing */

#ifdef WIN32
NLdouble nlCurrentTime() {
    return (NLdouble)GetTickCount() / 1000.0 ;
}
#else
double nlCurrentTime() {
    clock_t user_clock ;
    struct tms user_tms ;
    user_clock = times(&user_tms) ;
    return (NLdouble)user_clock / 100.0 ;
}
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





/******* extracted from nl_matrix.c *******/


/*
 Some warnings about const cast in callback for
 qsort() function.
 */

#ifdef __clang__
#pragma GCC diagnostic ignored "-Wcast-qual"
#endif




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

NLuint nlCRSMatrixNNZ(NLCRSMatrix* M) {
    return M->rowptr[M->m];
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
    nlHostBlas()->flops += (NLulong)(2*nlCRSMatrixNNZ(M));
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
    nl_assert(M->type == NL_MATRIX_SPARSE_DYNAMIC);
    nlSparseMatrixDestroyRowColumns(M);
    NL_DELETE_ARRAY(M->diag);
#ifdef NL_PARANOID
    NL_CLEAR(NLSparseMatrix,M);
#endif
}

void nlSparseMatrixAdd(NLSparseMatrix* M, NLuint i, NLuint j, NLdouble value) {
    nl_parano_range_assert(i, 0, M->m - 1);
    nl_parano_range_assert(j, 0, M->n - 1);
    if(i == j) {
        M->diag[i] += value;
    }
    nlRowColumnAdd(&(M->row[i]), j, value);
}

static void nlSparseMatrixAddSparseMatrix(
    NLSparseMatrix* M, double mul, const NLSparseMatrix* N    
) {
    NLuint i,jj;
    nl_assert(M->m == N->m);
    nl_assert(M->n == N->n);
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
    nl_assert(M->m == N->m);
    nl_assert(M->n == N->n);
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
    nl_assert(M->m == N->m);
    nl_assert(M->n == N->n);
    if(N->type == NL_MATRIX_SPARSE_DYNAMIC) {
	nlSparseMatrixAddSparseMatrix(M, mul, (const NLSparseMatrix*)N);
    } else if(N->type == NL_MATRIX_CRS) {
	nlSparseMatrixAddCRSMatrix(M, mul, (const NLCRSMatrix*)N);	
    } else {
	nl_assert_not_reached;
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

    nl_debug_assert(i1 < M->m);
    nl_debug_assert(i2 < M->m);
    
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

    nl_debug_assert(i < M->m);
    
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

    nl_debug_assert(i < M->m);
    
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
    nl_assert(A->type == NL_MATRIX_SPARSE_DYNAMIC);
    nlSparseMatrix_mult_rows(A, x, y);
    nlHostBlas()->flops += (NLulong)(2*nlSparseMatrixNNZ(A));
}

NLMatrix nlSparseMatrixNew(
    NLuint m, NLuint n
) {
    NLSparseMatrix* result = NL_NEW(NLSparseMatrix);
    nlSparseMatrixConstruct(result, m, n);
    return (NLMatrix)result;
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

#ifdef NL_PARANOID
    NL_CLEAR(NLContextStruct, context);
#endif
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
    nl_assert(nlCurrentContext->state == state);
}

void nlTransition(NLenum from_state, NLenum to_state) {
    nlCheckState(from_state);
    nlCurrentContext->state = to_state;
}


/* Preconditioner setup and default solver */

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
        nl_assert_not_reached;
    }
    nlMatrixCompress(&nlCurrentContext->M);
}

static NLboolean nlSolveIterative() {
    NLdouble* b = nlCurrentContext->b;
    NLdouble* x = nlCurrentContext->x;
    NLuint n = nlCurrentContext->n;
    NLuint k;
    NLBlas_t blas = nlHostBlas();
    NLMatrix M = nlCurrentContext->M;
    NLMatrix P = nlCurrentContext->P;
    
    nlCurrentContext->start_time = nlCurrentTime();     
    nlBlasResetStats(blas);
    
    for(k=0; k<nlCurrentContext->nb_systems; ++k) {
	nlSolveSystemIterative(
	    blas,
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

    nlCurrentContext->flops += blas->flops;
    return NL_TRUE;
}

/******* extracted from nl_blas.c *******/


/*
 Many warnings about const double* converted to
 double* when calling BLAS functions that do not
 have the const qualifier in their prototypes.
*/
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wcast-qual"
#pragma GCC diagnostic ignored "-Wcomma"
#endif

#ifndef NL_FORTRAN_WRAP
#define NL_FORTRAN_WRAP(x) x##_
#endif

/* BLAS routines                                                           */
/* copy-pasted from CBLAS (i.e. generated from f2c) */

/*
 * daxpy
 * ddot
     */



typedef NLint     integer ;
typedef NLdouble  doublereal ;
typedef NLboolean logical ;
typedef NLint     ftnlen ;


#ifndef max
#define max(x,y) ((x) > (y) ? (x) : (y))
#endif

/* Subroutine */ static int NL_FORTRAN_WRAP(daxpy)(integer *n, doublereal *da, doublereal *dx, 
        integer *incx, doublereal *dy, integer *incy)
{


    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m, ix, iy, mp1;


/*     constant times a vector plus a vector.   
       uses unrolled loops for increments equal to one.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]


    if (*n <= 0) {
        return 0;
    }
    if (*da == 0.) {
        return 0;
    }
    if (*incx == 1 && *incy == 1) {
        goto L20;
    }

/*        code for unequal increments or equal increments   
            not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
        ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
        iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
        DY(iy) += *da * DX(ix);
        ix += *incx;
        iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1   


          clean-up loop */

L20:
    m = *n % 4;
    if (m == 0) {
        goto L40;
    }
    i__1 = m;
    for (i = 1; i <= m; ++i) {
        DY(i) += *da * DX(i);
/* L30: */
    }
    if (*n < 4) {
        return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= *n; i += 4) {
        DY(i) += *da * DX(i);
        DY(i + 1) += *da * DX(i + 1);
        DY(i + 2) += *da * DX(i + 2);
        DY(i + 3) += *da * DX(i + 3);
/* L50: */
    }
    nl_arg_used(i__1);
    return 0;
} /* daxpy_ */
#undef DY
#undef DX


static doublereal NL_FORTRAN_WRAP(ddot)(integer *n, doublereal *dx, integer *incx, doublereal *dy, 
        integer *incy)
{

    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i, m;
    static doublereal dtemp;
    static integer ix, iy, mp1;


/*     forms the dot product of two vectors.   
       uses unrolled loops for increments equal to one.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]

    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0) {
        return ret_val;
    }
    if (*incx == 1 && *incy == 1) {
        goto L20;
    }

/*        code for unequal increments or equal increments   
            not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
        ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
        iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
        dtemp += DX(ix) * DY(iy);
        ix += *incx;
        iy += *incy;
/* L10: */
    }
    ret_val = dtemp;
    return ret_val;

/*        code for both increments equal to 1   


          clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
        goto L40;
    }
    i__1 = m;
    for (i = 1; i <= m; ++i) {
        dtemp += DX(i) * DY(i);
/* L30: */
    }
    if (*n < 5) {
        goto L60;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= *n; i += 5) {
        dtemp = dtemp + DX(i) * DY(i) + DX(i + 1) * DY(i + 1) + DX(i + 2) * 
                DY(i + 2) + DX(i + 3) * DY(i + 3) + DX(i + 4) * DY(i + 4);
/* L50: */
    }
L60:
    ret_val = dtemp;
    nl_arg_used(i__1);
    return ret_val;
} /* ddot_ */
#undef DY
#undef DX

/* End of BLAS routines */

/* Abstract BLAS interface                                              */


void nlBlasResetStats(NLBlas_t blas) {
    blas->start_time = nlCurrentTime();
    blas->flops = 0;
    blas->used_ram[0] = 0;
    blas->used_ram[1] = 0;
    blas->max_used_ram[0] = 0;
    blas->max_used_ram[1] = 0;
    blas->sq_rnorm = 0.0;
    blas->sq_bnorm = 0.0;
}

double nlBlasGFlops(NLBlas_t blas) {
    double now = nlCurrentTime();
    double elapsed_time = now - blas->start_time;
    return (NLdouble)(blas->flops) / (elapsed_time * 1e9);
}

NLulong nlBlasUsedRam(NLBlas_t blas, NLmemoryType type) {
    return blas->used_ram[type];
}

NLulong nlBlasMaxUsedRam(NLBlas_t blas, NLmemoryType type) {
    return blas->max_used_ram[type];
}

NLboolean nlBlasHasUnifiedMemory(NLBlas_t blas) {
    return blas->has_unified_memory;
}

static void* host_blas_malloc(
    NLBlas_t blas, NLmemoryType type, size_t size
) {
    nl_arg_used(type);
    blas->used_ram[type] += (NLulong)size;
    blas->max_used_ram[type] = MAX(
	blas->max_used_ram[type],blas->used_ram[type]
    );
    return malloc(size);
}

static void host_blas_free(
    NLBlas_t blas, NLmemoryType type, size_t size, void* ptr
) {
    nl_arg_used(type);
    blas->used_ram[type] -= (NLulong)size;
    free(ptr);
}

static void host_blas_memcpy(
    NLBlas_t blas,
    void* to, NLmemoryType to_type,
    void* from, NLmemoryType from_type,
    size_t size
) {
    nl_arg_used(blas);
    nl_arg_used(to_type);
    nl_arg_used(from_type);
    memcpy(to,from,size);
}

static double host_blas_ddot(
    NLBlas_t blas, int n, const double *x, int incx, const double *y, int incy    
) {
    blas->flops += (NLulong)(2*n);
    return NL_FORTRAN_WRAP(ddot)(&n,(double*)x,&incx,(double*)y,&incy);
}

static void host_blas_daxpy(
    NLBlas_t blas, int n, double a, const double *x, int incx, double *y, int incy
) {
    blas->flops += (NLulong)(2*n);
    NL_FORTRAN_WRAP(daxpy)(&n,&a,(double*)x,&incx,y,&incy);
}

static void host_blas_dscal(
    NLBlas_t blas, int n, double a, double *x, int incx
) {
    blas->flops += (NLulong)n;
	for (int i = 0; i < n; i++)
		x[i] *= a;
}

NLBlas_t nlHostBlas() {
    static NLboolean initialized = NL_FALSE;
    static struct NLBlas blas;
    if(!initialized) {
	memset(&blas, 0, sizeof(blas));
	blas.has_unified_memory = NL_TRUE;
	blas.Malloc = host_blas_malloc;
	blas.Free = host_blas_free;
	blas.Memcpy = host_blas_memcpy;
	blas.Ddot = host_blas_ddot;
	blas.Daxpy = host_blas_daxpy;
	blas.Dscal = host_blas_dscal;
	nlBlasResetStats(&blas);
	initialized = NL_TRUE;
    }
    return &blas;
}


/******* extracted from nl_iterative_solvers.c *******/



/* Solvers */

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
 *   - Abstract blas interface (NLBlas_t), that can implement BLAS 
 *     operations on the CPU or on the GPU.
 *   - Abstract matrix interface (NLMatrix), that can implement different
 *     versions of matrix x vector product (CPU/GPU, sparse/dense ...)
 */

static NLuint nlSolveSystem_PRE_CG(
    NLBlas_t blas,
    NLMatrix M, NLMatrix P, NLdouble* b, NLdouble* x,
    double eps, NLuint max_iter
) {
    NLint     N        = (NLint)M->n;
    NLdouble* r = NL_NEW_VECTOR(blas, NL_DEVICE_MEMORY, N);
    NLdouble* d = NL_NEW_VECTOR(blas, NL_DEVICE_MEMORY, N);
    NLdouble* h = NL_NEW_VECTOR(blas, NL_DEVICE_MEMORY, N);
    NLdouble *Ad = h;
    NLuint its=0;
    NLdouble rh, alpha, beta;
    NLdouble b_square = blas->Ddot(blas,N,b,1,b,1);
    NLdouble err=eps*eps*b_square;
    NLdouble curr_err;

    nlMultMatrixVector(M,x,r);
    blas->Daxpy(blas,N,-1.,b,1,r,1);
    nlMultMatrixVector(P,r,d);
	blas->Memcpy(blas, h, NL_HOST_MEMORY, d, NL_HOST_MEMORY, N * sizeof(NLdouble));
    rh=blas->Ddot(blas,N,r,1,h,1);
    curr_err = blas->Ddot(blas,N,r,1,r,1);

    while ( curr_err >err && its < max_iter) {
	if(nlCurrentContext != NULL) {
	    if(nlCurrentContext->progress_func != NULL) {
		nlCurrentContext->progress_func(its, max_iter, curr_err, err);
	    }
	}
	nlMultMatrixVector(M,d,Ad);
        alpha=rh/blas->Ddot(blas,N,d,1,Ad,1);
        blas->Daxpy(blas,N,-alpha,d,1,x,1);
        blas->Daxpy(blas,N,-alpha,Ad,1,r,1);
	nlMultMatrixVector(P,r,h);
        beta=1./rh;
	rh=blas->Ddot(blas,N,r,1,h,1);
	beta*=rh;
        blas->Dscal(blas,N,beta,d,1);
        blas->Daxpy(blas,N,1.,h,1,d,1);
        ++its;
        curr_err = blas->Ddot(blas,N,r,1,r,1);
    }
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, N, r);
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, N, d);
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, N, h);
    blas->sq_bnorm = b_square;
    blas->sq_rnorm = curr_err;
    return its;
}

/* Main driver routine */

NLuint nlSolveSystemIterative(
    NLBlas_t blas,
    NLMatrix M, NLMatrix P, NLdouble* b_in, NLdouble* x_in,
    double eps, NLuint max_iter, NLuint inner_iter
) {
    NLuint N = M->n;
    NLuint result=0;
    NLdouble rnorm=0.0;
    NLdouble bnorm=0.0; 
    double* b = b_in;
    double* x = x_in;
    nl_assert(M->m == M->n);

    if(!nlBlasHasUnifiedMemory(blas)) {
	b = NL_NEW_VECTOR(blas, NL_DEVICE_MEMORY, (int)M->n);
	blas->Memcpy(
	    blas,
	    b, NL_DEVICE_MEMORY,
	    b_in, NL_HOST_MEMORY, (size_t)N*sizeof(double)
	);
	x = NL_NEW_VECTOR(blas, NL_DEVICE_MEMORY, (int)M->n);
	blas->Memcpy(
	    blas,
	    x, NL_DEVICE_MEMORY,
	    x_in, NL_HOST_MEMORY, (size_t)N*sizeof(double)
	);	
    }

	result = nlSolveSystem_PRE_CG(blas,M,P,b,x,eps,max_iter);

    /* Get residual norm and rhs norm from BLAS context */
    if(nlCurrentContext != NULL) {
	bnorm = sqrt(blas->sq_bnorm);
	rnorm = sqrt(blas->sq_rnorm);
	if(bnorm == 0.0) {
	    nlCurrentContext->error = rnorm;
	} else {
	    nlCurrentContext->error = rnorm/bnorm;
	}
    }
    nlCurrentContext->used_iterations = result;

    if(!nlBlasHasUnifiedMemory(blas)) {
	blas->Memcpy(
	    blas,
	    x_in, NL_HOST_MEMORY, x, NL_DEVICE_MEMORY, (size_t)N*sizeof(double)
	);	
	NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, (int)M->n, x);
	NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, (int)M->n, b);
    }
    
    return result;
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
    nlHostBlas()->flops += (NLulong)(M->n);    
}

NLMatrix nlNewJacobiPreconditioner(NLMatrix M_in) {
    NLSparseMatrix* M = NULL;
    NLJacobiPreconditioner* result = NULL;
    NLuint i;
    nl_assert(M_in->type == NL_MATRIX_SPARSE_DYNAMIC);
    nl_assert(M_in->m == M_in->n);
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

/******* extracted from nl_api.c *******/

static NLSparseMatrix* nlGetCurrentSparseMatrix() {
    NLSparseMatrix* result = NULL;
	nl_assert(nlCurrentContext->M != NULL);	    
	nl_assert(nlCurrentContext->M->type == NL_MATRIX_SPARSE_DYNAMIC);
	return (NLSparseMatrix*)(nlCurrentContext->M);
}

/* Get/Set parameters */

void nlSolverParameterd(NLenum pname, NLdouble param) {
    nlCheckState(NL_STATE_INITIAL);
    switch(pname) {
    case NL_THRESHOLD: {
        nl_assert(param >= 0);
        nlCurrentContext->threshold = (NLdouble)param;
        nlCurrentContext->threshold_defined = NL_TRUE;
    } break;
    case NL_OMEGA: {
        nl_range_assert(param,1.0,2.0);
        nlCurrentContext->omega = (NLdouble)param;
    } break;
    default: {
        nlError("nlSolverParameterd","Invalid parameter");
        nl_assert_not_reached;
    }
    }
}

void nlSolverParameteri(NLenum pname, NLint param) {
    nlCheckState(NL_STATE_INITIAL);
    switch(pname) {
    case NL_NB_VARIABLES: {
        nl_assert(param > 0);
        nlCurrentContext->nb_variables = (NLuint)param;
    } break;
    case NL_NB_SYSTEMS: {
	nl_assert(param > 0);
	nlCurrentContext->nb_systems = (NLuint)param;
    } break;
    case NL_MAX_ITERATIONS: {
        nl_assert(param > 0);
        nlCurrentContext->max_iterations = (NLuint)param;
        nlCurrentContext->max_iterations_defined = NL_TRUE;
    } break;
    case NL_INNER_ITERATIONS: {
        nl_assert(param > 0);
        nlCurrentContext->inner_iterations = (NLuint)param;
    } break;
    case NL_PRECONDITIONER: {
        nlCurrentContext->preconditioner = (NLuint)param;
        nlCurrentContext->preconditioner_defined = NL_TRUE;
    } break;
    default: {
        nlError("nlSolverParameteri","Invalid parameter");
        nl_assert_not_reached;
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
    case NL_ELAPSED_TIME: {
        *params = nlCurrentContext->elapsed_time;        
    } break;
    case NL_GFLOPS: {
        if(nlCurrentContext->elapsed_time == 0) {
            *params = 0.0;
        } else {
            *params = (NLdouble)(nlCurrentContext->flops) /
                (nlCurrentContext->elapsed_time * 1e9);
        }
    } break;
    default: {
        nlError("nlGetDoublev","Invalid parameter");
        nl_assert_not_reached;
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
        nl_assert_not_reached;
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
        nl_assert_not_reached;
    }
}

/* Get/Set Lock/Unlock variables */

void nlSetVariable(NLuint index, NLdouble value) {
    nlCheckState(NL_STATE_SYSTEM);
    nl_debug_range_assert(index, 0, nlCurrentContext->nb_variables - 1);
    NL_BUFFER_ITEM(nlCurrentContext->variable_buffer[0],index) = value;
}

NLdouble nlGetVariable(NLuint index) {
    nl_assert(nlCurrentContext->state != NL_STATE_INITIAL);
    nl_debug_range_assert(index, 0, nlCurrentContext->nb_variables - 1);
    return NL_BUFFER_ITEM(nlCurrentContext->variable_buffer[0],index);
}

void nlLockVariable(NLuint index) {
    nlCheckState(NL_STATE_SYSTEM);
    nl_debug_range_assert(index, 0, nlCurrentContext->nb_variables - 1);
    nlCurrentContext->variable_is_locked[index] = NL_TRUE;
}

/* System construction */

static void nlVariablesToVector() {
    NLuint n=nlCurrentContext->n;
    NLuint k,i,index;
    NLdouble value;
    
    nl_assert(nlCurrentContext->x != NULL);
    for(k=0; k<nlCurrentContext->nb_systems; ++k) {
	for(i=0; i<nlCurrentContext->nb_variables; ++i) {
	    if(!nlCurrentContext->variable_is_locked[i]) {
		index = nlCurrentContext->variable_index[i];
		nl_assert(index < nlCurrentContext->n);		
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

    nl_assert(nlCurrentContext->x != NULL);
    for(k=0; k<nlCurrentContext->nb_systems; ++k) {
	for(i=0; i<nlCurrentContext->nb_variables; ++i) {
	    if(!nlCurrentContext->variable_is_locked[i]) {
		index = nlCurrentContext->variable_index[i];
		nl_assert(index < nlCurrentContext->n);
		value = nlCurrentContext->x[index+k*n];
		NL_BUFFER_ITEM(nlCurrentContext->variable_buffer[k],i) = value;
	    }
	}
    }
}


static void nlBeginSystem() {
    NLuint k;
    
    nlTransition(NL_STATE_INITIAL, NL_STATE_SYSTEM);
    nl_assert(nlCurrentContext->nb_variables > 0);

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
    nl_debug_range_assert(index, 0, nlCurrentContext->nb_variables - 1);
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
        nl_assert_not_reached;
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
        nl_assert_not_reached;
    }
    }
}


/* nlSolve() driver routine */

NLboolean nlSolve() {
    NLboolean result;
    nlCheckState(NL_STATE_SYSTEM_CONSTRUCTED);
    nlCurrentContext->start_time = nlCurrentTime();
    nlCurrentContext->elapsed_time = 0.0;
    nlCurrentContext->flops = 0;    
	nlSetupPreconditioner();
	result = nlSolveIterative();
    nlVectorToVariables();
    nlCurrentContext->elapsed_time = nlCurrentTime() - nlCurrentContext->start_time;
    nlTransition(NL_STATE_SYSTEM_CONSTRUCTED, NL_STATE_SOLVED);
    return result;
}
