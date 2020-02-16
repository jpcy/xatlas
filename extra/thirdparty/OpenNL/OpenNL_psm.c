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
    UpperTriangle=0, LowerTriangle=1
} MatrixTriangle ;

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

typedef void (*FUNPTR_dcopy)(
    NLBlas_t blas, int n, const double *x, int incx, double *y, int incy
);

typedef void (*FUNPTR_dscal)(
    NLBlas_t blas, int n, double a, double *x, int incx
);


typedef double (*FUNPTR_ddot)(
    NLBlas_t blas, int n, const double *x, int incx, const double *y, int incy
);

typedef double (*FUNPTR_dnrm2)(NLBlas_t blas, int n, const double *x, int incx);

typedef void (*FUNPTR_daxpy)(
    NLBlas_t blas, int n,
    double a, const double *x, int incx, double *y, int incy
);


typedef void (*FUNPTR_dgemv)( 
    NLBlas_t blas, MatrixTranspose trans, int m, int n, double alpha,
    const double *A, int ldA, const double *x, int incx,
    double beta, double *y, int incy 
);


typedef void (*FUNPTR_dtpsv)(
    NLBlas_t blas, MatrixTriangle uplo, MatrixTranspose trans,
    MatrixUnitTriangular diag, int n, const double *AP,
    double *x, int incx 
);

struct NLBlas {
    FUNPTR_malloc Malloc;
    FUNPTR_free Free;
    FUNPTR_memcpy Memcpy;

    FUNPTR_dcopy Dcopy;
    FUNPTR_dscal Dscal;
    FUNPTR_ddot  Ddot;
    FUNPTR_dnrm2 Dnrm2;
    FUNPTR_daxpy Daxpy;
    FUNPTR_dgemv Dgemv;
    FUNPTR_dtpsv Dtpsv;

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

NLAPI NLboolean NLAPIENTRY nlCRSMatrixLoad(
    NLCRSMatrix* M, const char* filename
);

NLAPI NLboolean NLAPIENTRY nlCRSMatrixSave(
    NLCRSMatrix* M, const char* filename
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

NLAPI NLuint NLAPIENTRY nlMatrixNNZ(NLMatrix M);

NLAPI NLMatrix NLAPIENTRY nlMatrixFactorize(NLMatrix M, NLenum solver);
    


    typedef void(*NLMatrixFunc)(const double* x, double* y);

NLAPI NLMatrix NLAPIENTRY nlMatrixNewFromFunction(
    NLuint m, NLuint n, NLMatrixFunc func
);	     

NLAPI NLMatrix NLAPIENTRY nlMatrixNewFromProduct(
    NLMatrix M, NLboolean product_owns_M,
    NLMatrix N, NLboolean product_owns_N
);


    
#ifdef __cplusplus
}
#endif

#endif

/******* extracted from nl_context.h *******/

#ifndef OPENNL_CONTEXT_H
#define OPENNL_CONTEXT_H




/* NLContext data structure */


typedef NLboolean(*NLSolverFunc)(void);

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


    NLenum           matrix_mode;

    NLMatrix         M;

    NLMatrix         P;

    NLMatrix         B;
    
    NLRowColumn      af;

    NLRowColumn      al;

    NLdouble*        x;

    NLdouble*        b;

    NLenum           solver;

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

    NLSolverFunc     solver_func;

    NLProgressFunc   progress_func;

    NLulong          flops;

    NLenum           eigen_solver;

    NLdouble         eigen_shift;

    NLboolean        eigen_shift_invert;

    NLdouble*        eigen_value;

    NLdouble*        temp_eigen_value;
    
} NLContextStruct;

extern NLContextStruct* nlCurrentContext;

void nlCheckState(NLenum state);

void nlTransition(NLenum from_state, NLenum to_state);

NLboolean nlDefaultSolver(void);

#endif

/******* extracted from nl_iterative_solvers.h *******/


#ifndef OPENNL_ITERATIVE_SOLVERS_H
#define OPENNL_ITERATIVE_SOLVERS_H



NLAPI NLuint NLAPIENTRY nlSolveSystemIterative(
    NLBlas_t blas,
    NLMatrix M, NLMatrix P, NLdouble* b, NLdouble* x,
    NLenum solver,
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


NLboolean nlCRSMatrixSave(NLCRSMatrix* M, const char* filename) {
    NLuint nnz = M->rowptr[M->m];
    FILE* f = fopen(filename, "rb");
    if(f == NULL) {
        nlError("nlCRSMatrixSave", "Could not open file");
        return NL_FALSE;
    }

    fwrite(&M->m, sizeof(NLuint), 1, f);
    fwrite(&M->n, sizeof(NLuint), 1, f);
    fwrite(&nnz, sizeof(NLuint), 1, f);

    fwrite(M->rowptr, sizeof(NLuint), M->m+1, f);
    fwrite(M->colind, sizeof(NLuint), nnz, f);
    fwrite(M->val, sizeof(double), nnz, f);
    
    return NL_TRUE;
}

NLboolean nlCRSMatrixLoad(NLCRSMatrix* M, const char* filename) {
    NLuint nnz = 0;
    FILE* f = fopen(filename, "rb");
    NLboolean truncated = NL_FALSE;
    
    if(f == NULL) {
        nlError("nlCRSMatrixLoad", "Could not open file");
        return NL_FALSE;
    }
    
    truncated = truncated || (
        fread(&M->m, sizeof(NLuint), 1, f) != 1 ||
        fread(&M->n, sizeof(NLuint), 1, f) != 1 ||
        fread(&nnz, sizeof(NLuint), 1, f) != 1
    );

    if(truncated) {
        M->rowptr = NULL;
        M->colind = NULL;
        M->val = NULL;
    } else {
        M->rowptr = NL_NEW_ARRAY(NLuint, M->m+1);
        M->colind = NL_NEW_ARRAY(NLuint, nnz);
        M->val = NL_NEW_ARRAY(double, nnz);
        truncated = truncated || (
            fread(M->rowptr, sizeof(NLuint), M->m+1, f) != M->m+1 ||
            fread(M->colind, sizeof(NLuint), nnz, f) != nnz ||
            fread(M->val, sizeof(double), nnz, f) != nnz
        );
    }

    if(truncated) {
        nlError("nlCRSMatrixSave", "File appears to be truncated");
        NL_DELETE_ARRAY(M->rowptr);
        NL_DELETE_ARRAY(M->colind);
        NL_DELETE_ARRAY(M->val);
        return NL_FALSE;
    } else {
        M->nslices = 1;    
        M->sliceptr = NL_NEW_ARRAY(NLuint, M->nslices+1);
        M->sliceptr[0] = 0;
        M->sliceptr[1] = M->m;
    }

    fclose(f);
    return NL_TRUE;
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

NLuint nlMatrixNNZ(NLMatrix M) {
    if(M->type == NL_MATRIX_SPARSE_DYNAMIC) {
	return nlSparseMatrixNNZ((NLSparseMatrix*)M);
    } else if(M->type == NL_MATRIX_CRS) {
	return nlCRSMatrixNNZ((NLCRSMatrix*)M);	
    }
    return M->m * M->n;
}

NLMatrix nlMatrixFactorize(NLMatrix M, NLenum solver) {
    NLMatrix result = NULL;
    switch(solver) {
	default:
	    nlError("nlMatrixFactorize","unknown solver");
    }
    return result;
}

typedef struct {
    NLuint m;

    NLuint n;

    NLenum type;

    NLDestroyMatrixFunc destroy_func;

    NLMultMatrixVectorFunc mult_func;

    NLMatrixFunc matrix_func;

    NLMatrix M;

    NLboolean owns_M;
    
    NLMatrix N;

    NLboolean owns_N;
    
    NLdouble* work;
} NLMatrixProduct;


static void nlMatrixProductDestroy(NLMatrixProduct* P) {
    NL_DELETE_ARRAY(P->work);
    if(P->owns_M) {
	nlDeleteMatrix(P->M); P->M = NULL;
    }
    if(P->owns_N) {
	nlDeleteMatrix(P->N); P->N = NULL;
    }
}

static void nlMatrixProductMult(
    NLMatrixProduct* P, const NLdouble* x, NLdouble* y
) {
    nlMultMatrixVector(P->N, x, P->work);
    nlMultMatrixVector(P->M, P->work, y);
}

NLMatrix nlMatrixNewFromProduct(
    NLMatrix M, NLboolean owns_M, NLMatrix N, NLboolean owns_N
) {
    NLMatrixProduct* result = NL_NEW(NLMatrixProduct);
    nl_assert(M->n == N->m);
    result->m = M->m;
    result->n = N->n;
    result->type = NL_MATRIX_OTHER;
    result->work = NL_NEW_ARRAY(NLdouble,N->m);
    result->destroy_func = (NLDestroyMatrixFunc)nlMatrixProductDestroy;
    result->mult_func = (NLMultMatrixVectorFunc)nlMatrixProductMult;
    result->M = M;
    result->owns_M = owns_M;
    result->N = N;
    result->owns_N = owns_N;
    return (NLMatrix)result;
}



/******* extracted from nl_context.c *******/


NLContextStruct* nlCurrentContext = NULL;

NLContext nlNewContext() {
    NLContextStruct* result     = NL_NEW(NLContextStruct);
    result->state               = NL_STATE_INITIAL;
    result->solver              = NL_SOLVER_DEFAULT;
    result->max_iterations      = 100;
    result->threshold           = 1e-6;
    result->omega               = 1.5;
    result->inner_iterations    = 5;
    result->solver_func         = nlDefaultSolver;
    result->progress_func       = NULL;
    result->nb_systems          = 1;
    result->matrix_mode         = NL_STIFFNESS_MATRIX;
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

    NL_DELETE_ARRAY(context->eigen_value);
    
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

    if(getenv("NL_LOW_MEM") == NULL) {
        nlMatrixCompress(&nlCurrentContext->M);
    }
}

static NLboolean nlSolveDirect() {
    NLdouble* b = nlCurrentContext->b;
    NLdouble* x = nlCurrentContext->x;
    NLuint n = nlCurrentContext->n;
    NLuint k;
    
    NLMatrix F = nlMatrixFactorize(
	nlCurrentContext->M, nlCurrentContext->solver
    );
    if(F == NULL) {
	return NL_FALSE;
    }
    for(k=0; k<nlCurrentContext->nb_systems; ++k) {
	nlMultMatrixVector(F, b, x);
	b += n;
	x += n;
    }
    nlDeleteMatrix(F);
    return NL_TRUE;
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
	    nlCurrentContext->solver,
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



NLboolean nlDefaultSolver() {
    NLboolean result = NL_TRUE;
    nlSetupPreconditioner();
    switch(nlCurrentContext->solver) {
	case NL_CG:
    {
	    result = nlSolveIterative();
	} break;
	default:
	    nl_assert_not_reached;
    }
    return result;
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

#ifdef NL_USE_ATLAS
int NL_FORTRAN_WRAP(xerbla)(char *srname, int *info) {
    nl_printf(stderr, "** On entry to %6s, parameter number %2d had an illegal value\n",
              srname, *info
    );
    return 0;
} 
#ifndef NL_USE_BLAS
#define NL_USE_BLAS
#endif
#endif

#ifndef NL_USE_BLAS
#define NEEDS_DTPSV
#endif



/* BLAS routines                                                           */
/* copy-pasted from CBLAS (i.e. generated from f2c) */

/*
 * lsame
 * xerbla
 * daxpy
 * ddot
 * dscal
 * dnrm2
 * dcopy
 * dgemv
 * dtpsv
 */



typedef NLint     integer ;
typedef NLdouble  doublereal ;
typedef NLboolean logical ;
typedef NLint     ftnlen ;


#ifndef max
#define max(x,y) ((x) > (y) ? (x) : (y))
#endif

#ifndef NL_USE_BLAS

static int NL_FORTRAN_WRAP(lsame)(const char *ca, const char *cb)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   

    Purpose   
    =======   

    LSAME returns .TRUE. if CA is the same letter as CB regardless of case.   

    Arguments   
    =========   

    CA      (input) CHARACTER*1   
    CB      (input) CHARACTER*1   
            CA and CB specify the single characters to be compared.   

   ===================================================================== 
*/  

    /* System generated locals */
    int ret_val;
    
    /* Local variables */
    int inta, intb, zcode;

    ret_val = *(unsigned char *)ca == *(unsigned char *)cb;
    if (ret_val) {
        return ret_val;
    }

    /* Now test for equivalence if both characters are alphabetic. */

    zcode = 'Z';

    /* Use 'Z' rather than 'A' so that ASCII can be detected on Prime   
       machines, on which ICHAR returns a value with bit 8 set.   
       ICHAR('A') on Prime machines returns 193 which is the same as   
       ICHAR('A') on an EBCDIC machine. */

    inta = *(unsigned char *)ca;
    intb = *(unsigned char *)cb;

    if (zcode == 90 || zcode == 122) {
        /* ASCII is assumed - ZCODE is the ASCII code of either lower or   
          upper case 'Z'. */
        if (inta >= 97 && inta <= 122) inta += -32;
        if (intb >= 97 && intb <= 122) intb += -32;

    } else if (zcode == 233 || zcode == 169) {
        /* EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or   
          upper case 'Z'. */
        if ((inta >= 129 && inta <= 137) || 
            (inta >= 145 && inta <= 153) || 
            (inta >= 162 && inta <= 169)
        )
            inta += 64;
        if (
            (intb >= 129 && intb <= 137) || 
            (intb >= 145 && intb <= 153) || 
            (intb >= 162 && intb <= 169)
        )
            intb += 64;
    } else if (zcode == 218 || zcode == 250) {
        /* ASCII is assumed, on Prime machines - ZCODE is the ASCII code   
          plus 128 of either lower or upper case 'Z'. */
        if (inta >= 225 && inta <= 250) inta += -32;
        if (intb >= 225 && intb <= 250) intb += -32;
    }
    ret_val = inta == intb;
    return ret_val;
    
} /* lsame_ */

/* Subroutine */ static int NL_FORTRAN_WRAP(xerbla)(const char *srname, int *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    XERBLA  is an error handler for the LAPACK routines.   
    It is called by an LAPACK routine if an input parameter has an   
    invalid value.  A message is printed and execution stops.   

    Installers may consider modifying the STOP statement in order to   
    call system-specific exception-handling facilities.   

    Arguments   
    =========   

    SRNAME  (input) CHARACTER*6   
            The name of the routine which called XERBLA.   

    INFO    (input) INT   
            The position of the invalid parameter in the parameter list   

            of the calling routine.   

   ===================================================================== 
*/

    nl_fprintf(stderr, "** On entry to %6s, parameter number %2d had an illegal value\n",
                srname, *info);

/*     End of XERBLA */

    return 0;
} /* xerbla_ */


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

/* Subroutine */ static int NL_FORTRAN_WRAP(dscal)(integer *n, doublereal *da, doublereal *dx, 
    integer *incx)
{


    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i, m, nincx, mp1;


/*     scales a vector by a constant.   
       uses unrolled loops for increment equal to one.   
       jack dongarra, linpack, 3/11/78.   
       modified 3/93 to return if incx .le. 0.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#ifdef DX
#undef DX
#endif
#define DX(I) dx[(I)-1]


    if (*n <= 0 || *incx <= 0) {
        return 0;
    }
    if (*incx == 1) {
        goto L20;
    }

/*        code for increment not equal to 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i = 1; *incx < 0 ? i >= nincx : i <= nincx; i += *incx) {
        DX(i) = *da * DX(i);
/* L10: */
    }
    return 0;

/*        code for increment equal to 1   


          clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
        goto L40;
    }
    i__2 = m;
    for (i = 1; i <= m; ++i) {
        DX(i) = *da * DX(i);
/* L30: */
    }
    if (*n < 5) {
        return 0;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i = mp1; i <= *n; i += 5) {
        DX(i) = *da * DX(i);
        DX(i + 1) = *da * DX(i + 1);
        DX(i + 2) = *da * DX(i + 2);
        DX(i + 3) = *da * DX(i + 3);
        DX(i + 4) = *da * DX(i + 4);
/* L50: */
    }
    nl_arg_used(i__1);
    nl_arg_used(i__2);
    return 0;
} /* dscal_ */
#undef DX

static doublereal NL_FORTRAN_WRAP(dnrm2)(integer *n, doublereal *x, integer *incx)
{


    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1;

    /* Builtin functions */
    /* BL: already declared in the included <math.h>, 
       we do not need it here. */
    /*double sqrt(doublereal); */

    /* Local variables */
    static doublereal norm, scale, absxi;
    static integer ix;
    static doublereal ssq;


/*  DNRM2 returns the euclidean norm of a vector via the function   
    name, so that   

       DNRM2 := sqrt( x'*x )   



    -- This version written on 25-October-1982.   
       Modified on 14-October-1993 to inline the call to DLASSQ.   
       Sven Hammarling, Nag Ltd.   


    
   Parameter adjustments   
       Function Body */
#ifdef X
#undef X
#endif
#define X(I) x[(I)-1]


    if (*n < 1 || *incx < 1) {
        norm = 0.;
    } else if (*n == 1) {
        norm = fabs(X(1));
    } else {
        scale = 0.;
        ssq = 1.;
/*        The following loop is equivalent to this call to the LAPACK 
  
          auxiliary routine:   
          CALL DLASSQ( N, X, INCX, SCALE, SSQ ) */

        i__1 = (*n - 1) * *incx + 1;
        i__2 = *incx;
        for (ix = 1; *incx < 0 ? ix >= (*n-1)**incx+1 : ix <= (*n-1)**incx+1; ix += *incx) {
            if (X(ix) != 0.) {
                absxi = (d__1 = X(ix), fabs(d__1));
                if (scale < absxi) {
/* Computing 2nd power */
                    d__1 = scale / absxi;
                    ssq = ssq * (d__1 * d__1) + 1.;
                    scale = absxi;
                } else {
/* Computing 2nd power */
                    d__1 = absxi / scale;
                    ssq += d__1 * d__1;
                }
            }
/* L10: */
        }
        norm = scale * sqrt(ssq);
    }

    ret_val = norm;

    nl_arg_used(i__1);
    nl_arg_used(i__2);

    return ret_val;

/*     End of DNRM2. */

} /* dnrm2_ */
#undef X

/* Subroutine */ static int NL_FORTRAN_WRAP(dcopy)(integer *n, doublereal *dx, integer *incx, 
        doublereal *dy, integer *incy)
{

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m, ix, iy, mp1;


/*     copies a vector, x, to a vector, y.   
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
        DY(iy) = DX(ix);
        ix += *incx;
        iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1   


          clean-up loop */

L20:
    m = *n % 7;
    if (m == 0) {
        goto L40;
    }
    i__1 = m;
    for (i = 1; i <= m; ++i) {
        DY(i) = DX(i);
/* L30: */
    }
    if (*n < 7) {
        return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= *n; i += 7) {
        DY(i) = DX(i);
        DY(i + 1) = DX(i + 1);
        DY(i + 2) = DX(i + 2);
        DY(i + 3) = DX(i + 3);
        DY(i + 4) = DX(i + 4);
        DY(i + 5) = DX(i + 5);
        DY(i + 6) = DX(i + 6);
/* L50: */
    }
    nl_arg_used(i__1);
    return 0;
} /* dcopy_ */

#undef DX
#undef DY

/* Subroutine */ static int NL_FORTRAN_WRAP(dgemv)(const char *trans, integer *m, integer *n, doublereal *
        alpha, doublereal *a, integer *lda, doublereal *x, integer *incx, 
        doublereal *beta, doublereal *y, integer *incy)
{


    /* System generated locals */
    /* integer a_dim1, a_offset ; */
    integer i__1, i__2; 

    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer lenx, leny, i, j;
/*    extern logical lsame_(char *, char *); */
    static integer ix, iy, jx, jy, kx, ky;
/*    extern int xerbla_(char *, integer *); */


/*  Purpose   
    =======   

    DGEMV  performs one of the matrix-vector operations   

       y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   

    where alpha and beta are scalars, x and y are vectors and A is an   
    m by n matrix.   

    Parameters   
    ==========   

    TRANS  - CHARACTER*1.   
             On entry, TRANS specifies the operation to be performed as   
             follows:   

                TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.   

                TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.   

                TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.   

             Unchanged on exit.   

    M      - INTEGER.   
             On entry, M specifies the number of rows of the matrix A.   
             M must be at least zero.   
             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the number of columns of the matrix A. 
  
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - DOUBLE PRECISION.   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
             Before entry, the leading m by n part of the array A must   
             contain the matrix of coefficients.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, m ).   
             Unchanged on exit.   

    X      - DOUBLE PRECISION array of DIMENSION at least   
             ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'   
             and at least   
             ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.   
             Before entry, the incremented array X must contain the   
             vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    BETA   - DOUBLE PRECISION.   
             On entry, BETA specifies the scalar beta. When BETA is   
             supplied as zero then Y need not be set on input.   
             Unchanged on exit.   

    Y      - DOUBLE PRECISION array of DIMENSION at least   
             ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'   
             and at least   
             ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.   
             Before entry with BETA non-zero, the incremented array Y   
             must contain the vector y. On exit, Y is overwritten by the 
  
             updated vector y.   

    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   


    Level 2 Blas routine.   

    -- Written on 22-October-1986.   
       Jack Dongarra, Argonne National Lab.   
       Jeremy Du Croz, Nag Central Office.   
       Sven Hammarling, Nag Central Office.   
       Richard Hanson, Sandia National Labs.   



       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! NL_FORTRAN_WRAP(lsame)(trans, "N") && ! NL_FORTRAN_WRAP(lsame)(trans, "T") && ! 
            NL_FORTRAN_WRAP(lsame)(trans, "C")) {
        info = 1;
    } else if (*m < 0) {
        info = 2;
    } else if (*n < 0) {
        info = 3;
    } else if (*lda < max(1,*m)) {
        info = 6;
    } else if (*incx == 0) {
        info = 8;
    } else if (*incy == 0) {
        info = 11;
    }
    if (info != 0) {
        NL_FORTRAN_WRAP(xerbla)("DGEMV ", &info);
        return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || (*alpha == 0. && *beta == 1.)) {
        return 0;
    }

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set 
  
       up the start points in  X  and  Y. */

    if (NL_FORTRAN_WRAP(lsame)(trans, "N")) {
        lenx = *n;
        leny = *m;
    } else {
        lenx = *m;
        leny = *n;
    }
    if (*incx > 0) {
        kx = 1;
    } else {
        kx = 1 - (lenx - 1) * *incx;
    }
    if (*incy > 0) {
        ky = 1;
    } else {
        ky = 1 - (leny - 1) * *incy;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A.   

       First form  y := beta*y. */

    if (*beta != 1.) {
        if (*incy == 1) {
            if (*beta == 0.) {
                i__1 = leny;
                for (i = 1; i <= leny; ++i) {
                    Y(i) = 0.;
/* L10: */
                }
            } else {
                i__1 = leny;
                for (i = 1; i <= leny; ++i) {
                    Y(i) = *beta * Y(i);
/* L20: */
                }
            }
        } else {
            iy = ky;
            if (*beta == 0.) {
                i__1 = leny;
                for (i = 1; i <= leny; ++i) {
                    Y(iy) = 0.;
                    iy += *incy;
/* L30: */
                }
            } else {
                i__1 = leny;
                for (i = 1; i <= leny; ++i) {
                    Y(iy) = *beta * Y(iy);
                    iy += *incy;
/* L40: */
                }
            }
        }
    }
    if (*alpha == 0.) {
        return 0;
    }
    if (NL_FORTRAN_WRAP(lsame)(trans, "N")) {

/*        Form  y := alpha*A*x + y. */

        jx = kx;
        if (*incy == 1) {
            i__1 = *n;
            for (j = 1; j <= *n; ++j) {
                if (X(jx) != 0.) {
                    temp = *alpha * X(jx);
                    i__2 = *m;
                    for (i = 1; i <= *m; ++i) {
                        Y(i) += temp * A(i,j);
/* L50: */
                    }
                }
                jx += *incx;
/* L60: */
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= *n; ++j) {
                if (X(jx) != 0.) {
                    temp = *alpha * X(jx);
                    iy = ky;
                    i__2 = *m;
                    for (i = 1; i <= *m; ++i) {
                        Y(iy) += temp * A(i,j);
                        iy += *incy;
/* L70: */
                    }
                }
                jx += *incx;
/* L80: */
            }
        }
    } else {

/*        Form  y := alpha*A'*x + y. */

        jy = ky;
        if (*incx == 1) {
            i__1 = *n;
            for (j = 1; j <= *n; ++j) {
                temp = 0.;
                i__2 = *m;
                for (i = 1; i <= *m; ++i) {
                    temp += A(i,j) * X(i);
/* L90: */
                }
                Y(jy) += *alpha * temp;
                jy += *incy;
/* L100: */
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= *n; ++j) {
                temp = 0.;
                ix = kx;
                i__2 = *m;
                for (i = 1; i <= *m; ++i) {
                    temp += A(i,j) * X(ix);
                    ix += *incx;
/* L110: */
                }
                Y(jy) += *alpha * temp;
                jy += *incy;
/* L120: */
            }
        }
    }

    nl_arg_used(i__1);
    nl_arg_used(i__2);
    return 0;

/*     End of DGEMV . */

} /* dgemv_ */

#undef X
#undef Y
#undef A




#else

extern void NL_FORTRAN_WRAP(daxpy)( 
    int *n, double *alpha, double *x,
    int *incx, double *y, int *incy 
) ;


extern double NL_FORTRAN_WRAP(dnrm2)( int *n, double *x, int *incx ) ;

extern int NL_FORTRAN_WRAP(dcopy)(int* n, double* dx, int* incx, double* dy, int* incy) ;

extern void NL_FORTRAN_WRAP(dscal)(int* n, double* alpha, double *x, int* incx) ;

#ifndef NEEDS_DTPSV
extern void NL_FORTRAN_WRAP(dtpsv)( 
    char *uplo, char *trans, char *diag,
    int *n, double *AP, double *x, int *incx 
) ;
#endif

extern void NL_FORTRAN_WRAP(dgemv)( 
    char *trans, int *m, int *n,
    double *alpha, double *A, int *ldA,
    double *x, int *incx,
    double *beta, double *y, int *incy 
) ;

#endif

#ifdef NEEDS_DTPSV

/* DECK DTPSV */
/* Subroutine */ 
static int NL_FORTRAN_WRAP(dtpsv)(
   const char* uplo, 
   const char* trans, 
   const char* diag, 
   integer* n, 
   doublereal* ap, 
   doublereal* x, 
   integer* incx
) {
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer i__, j, k;
/*    extern logical lsame_(); */
    static integer kk, ix, jx, kx;
/*    extern int xerbla_(); */
    static logical nounit;

/* ***BEGIN PROLOGUE  DTPSV */
/* ***PURPOSE  Solve one of the systems of equations. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1B4 */
/* ***TYPE      DOUBLE PRECISION (STPSV-S, DTPSV-D, CTPSV-C) */
/* ***KEYWORDS  LEVEL 2 BLAS, LINEAR ALGEBRA */
/* ***AUTHOR  Dongarra, J. J., (ANL) */
/*           Du Croz, J., (NAG) */
/*           Hammarling, S., (NAG) */
/*           Hanson, R. J., (SNLA) */
/* ***DESCRIPTION */

/*  DTPSV  solves one of the systems of equations */

/*     A*x = b,   or   A'*x = b, */

/*  where b and x are n element vectors and A is an n by n unit, or */
/*  non-unit, upper or lower triangular matrix, supplied in packed form. */

/*  No test for singularity or near-singularity is included in this */
/*  routine. Such tests must be performed before calling this routine. */

/*  Parameters */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the matrix is an upper or */
/*           lower triangular matrix as follows: */

/*              UPLO = 'U' or 'u'   A is an upper triangular matrix. */

/*              UPLO = 'L' or 'l'   A is a lower triangular matrix. */

/*           Unchanged on exit. */

/*  TRANS  - CHARACTER*1. */
/*           On entry, TRANS specifies the equations to be solved as */
/*           follows: */

/*              TRANS = 'N' or 'n'   A*x = b. */

/*              TRANS = 'T' or 't'   A'*x = b. */

/*              TRANS = 'C' or 'c'   A'*x = b. */

/*           Unchanged on exit. */

/*  DIAG   - CHARACTER*1. */
/*           On entry, DIAG specifies whether or not A is unit */
/*           triangular as follows: */

/*              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */

/*              DIAG = 'N' or 'n'   A is not assumed to be unit */
/*                                  triangular. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  AP     - DOUBLE PRECISION array of DIMENSION at least */
/*           ( ( n*( n + 1))/2). */
/*           Before entry with  UPLO = 'U' or 'u', the array AP must */
/*           contain the upper triangular matrix packed sequentially, */
/*           column by column, so that AP( 1 ) contains a( 1, 1 ), */
/*           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 ) */
/*           respectively, and so on. */
/*           Before entry with UPLO = 'L' or 'l', the array AP must */
/*           contain the lower triangular matrix packed sequentially, */
/*           column by column, so that AP( 1 ) contains a( 1, 1 ), */
/*           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 ) */
/*           respectively, and so on. */
/*           Note that when  DIAG = 'U' or 'u', the diagonal elements of */
/*           A are not referenced, but are assumed to be unity. */
/*           Unchanged on exit. */

/*  X      - DOUBLE PRECISION array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the n */
/*           element right-hand side vector b. On exit, X is overwritten */
/*           with the solution vector x. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/* ***REFERENCES  Dongarra, J. J., Du Croz, J., Hammarling, S., and */
/*                 Hanson, R. J.  An extended set of Fortran basic linear */
/*                 algebra subprograms.  ACM TOMS, Vol. 14, No. 1, */
/*                 pp. 1-17, March 1988. */
/* ***ROUTINES CALLED  LSAME, XERBLA */
/* ***REVISION HISTORY  (YYMMDD) */
/*   861022  DATE WRITTEN */
/*   910605  Modified to meet SLATEC prologue standards.  Only comment */
/*           lines were modified.  (BKS) */
/* ***END PROLOGUE  DTPSV */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/* ***FIRST EXECUTABLE STATEMENT  DTPSV */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --x;
    --ap;

    /* Function Body */
    info = 0;
    if (!NL_FORTRAN_WRAP(lsame)(uplo, "U") && 
        !NL_FORTRAN_WRAP(lsame)(uplo, "L")
    ) {
        info = 1;
    } else if (
        !NL_FORTRAN_WRAP(lsame)(trans, "N") && 
        !NL_FORTRAN_WRAP(lsame)(trans, "T") && 
        !NL_FORTRAN_WRAP(lsame)(trans, "C")
    ) {
        info = 2;
    } else if (
        !NL_FORTRAN_WRAP(lsame)(diag, "U") && 
        !NL_FORTRAN_WRAP(lsame)(diag, "N")
    ) {
        info = 3;
    } else if (*n < 0) {
        info = 4;
    } else if (*incx == 0) {
        info = 7;
    }
    if (info != 0) {
        NL_FORTRAN_WRAP(xerbla)("DTPSV ", &info);
        return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
        return 0;
    }

    nounit = (logical)(NL_FORTRAN_WRAP(lsame)(diag, "N"));

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
        kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
        kx = 1;
    }

/*     Start the operations. In this version the elements of AP are */
/*     accessed sequentially with one pass through AP. */

    if (NL_FORTRAN_WRAP(lsame)(trans, "N")) {

/*        Form  x := inv( A )*x. */

        if (NL_FORTRAN_WRAP(lsame)(uplo, "U")) {
            kk = *n * (*n + 1) / 2;
            if (*incx == 1) {
                for (j = *n; j >= 1; --j) {
                    if (x[j] != 0.) {
                        if (nounit) {
                            x[j] /= ap[kk];
                        }
                        temp = x[j];
                        k = kk - 1;
                        for (i__ = j - 1; i__ >= 1; --i__) {
                            x[i__] -= temp * ap[k];
                            --k;
/* L10: */
                        }
                    }
                    kk -= j;
/* L20: */
                }
            } else {
                jx = kx + (*n - 1) * *incx;
                for (j = *n; j >= 1; --j) {
                    if (x[jx] != 0.) {
                        if (nounit) {
                            x[jx] /= ap[kk];
                        }
                        temp = x[jx];
                        ix = jx;
                        i__1 = kk - j + 1;
                        for (k = kk - 1; k >= i__1; --k) {
                            ix -= *incx;
                            x[ix] -= temp * ap[k];
/* L30: */
                        }
                    }
                    jx -= *incx;
                    kk -= j;
/* L40: */
                }
            }
        } else {
            kk = 1;
            if (*incx == 1) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    if (x[j] != 0.) {
                        if (nounit) {
                            x[j] /= ap[kk];
                        }
                        temp = x[j];
                        k = kk + 1;
                        i__2 = *n;
                        for (i__ = j + 1; i__ <= i__2; ++i__) {
                            x[i__] -= temp * ap[k];
                            ++k;
/* L50: */
                        }
                    }
                    kk += *n - j + 1;
/* L60: */
                }
            } else {
                jx = kx;
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    if (x[jx] != 0.) {
                        if (nounit) {
                            x[jx] /= ap[kk];
                        }
                        temp = x[jx];
                        ix = jx;
                        i__2 = kk + *n - j;
                        for (k = kk + 1; k <= i__2; ++k) {
                            ix += *incx;
                            x[ix] -= temp * ap[k];
/* L70: */
                        }
                    }
                    jx += *incx;
                    kk += *n - j + 1;
/* L80: */
                }
            }
        }
    } else {

/*        Form  x := inv( A' )*x. */

        if (NL_FORTRAN_WRAP(lsame)(uplo, "U")) {
            kk = 1;
            if (*incx == 1) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    temp = x[j];
                    k = kk;
                    i__2 = j - 1;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        temp -= ap[k] * x[i__];
                        ++k;
/* L90: */
                    }
                    if (nounit) {
                        temp /= ap[kk + j - 1];
                    }
                    x[j] = temp;
                    kk += j;
/* L100: */
                }
            } else {
                jx = kx;
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    temp = x[jx];
                    ix = kx;
                    i__2 = kk + j - 2;
                    for (k = kk; k <= i__2; ++k) {
                        temp -= ap[k] * x[ix];
                        ix += *incx;
/* L110: */
                    }
                    if (nounit) {
                        temp /= ap[kk + j - 1];
                    }
                    x[jx] = temp;
                    jx += *incx;
                    kk += j;
/* L120: */
                }
            }
        } else {
            kk = *n * (*n + 1) / 2;
            if (*incx == 1) {
                for (j = *n; j >= 1; --j) {
                    temp = x[j];
                    k = kk;
                    i__1 = j + 1;
                    for (i__ = *n; i__ >= i__1; --i__) {
                        temp -= ap[k] * x[i__];
                        --k;
/* L130: */
                    }
                    if (nounit) {
                        temp /= ap[kk - *n + j];
                    }
                    x[j] = temp;
                    kk -= *n - j + 1;
/* L140: */
                }
            } else {
                kx += (*n - 1) * *incx;
                jx = kx;
                for (j = *n; j >= 1; --j) {
                    temp = x[jx];
                    ix = kx;
                    i__1 = kk - (*n - (j + 1));
                    for (k = kk; k >= i__1; --k) {
                        temp -= ap[k] * x[ix];
                        ix -= *incx;
/* L150: */
                    }
                    if (nounit) {
                        temp /= ap[kk - *n + j];
                    }
                    x[jx] = temp;
                    jx -= *incx;
                    kk -= *n - j + 1;
/* L160: */
                }
            }
        }
    }

    return 0;

/*     End of DTPSV . */

} /* dtpsv_ */

#endif 


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

static void host_blas_dcopy(
    NLBlas_t blas, int n, const double *x, int incx, double *y, int incy    
) {
    nl_arg_used(blas);
    NL_FORTRAN_WRAP(dcopy)(&n,(double*)x,&incx,y,&incy);    
}

static double host_blas_ddot(
    NLBlas_t blas, int n, const double *x, int incx, const double *y, int incy    
) {
    blas->flops += (NLulong)(2*n);
    return NL_FORTRAN_WRAP(ddot)(&n,(double*)x,&incx,(double*)y,&incy);
}

static double host_blas_dnrm2(
    NLBlas_t blas, int n, const double *x, int incx
) {
    blas->flops += (NLulong)(2*n);
    return NL_FORTRAN_WRAP(dnrm2)(&n,(double*)x,&incx);
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
    NL_FORTRAN_WRAP(dscal)(&n,&a,x,&incx);    
}

static void host_blas_dgemv(
    NLBlas_t blas, MatrixTranspose trans, int m, int n, double alpha,
    const double *A, int ldA, const double *x, int incx,
    double beta, double *y, int incy 
) {
    static const char *T[3] = { "N", "T", 0 };
    nl_arg_used(blas);
    NL_FORTRAN_WRAP(dgemv)(
	T[(int)trans],&m,&n,&alpha,(double*)A,&ldA,
	(double*)x,&incx,&beta,y,&incy
    );
    /* TODO: update flops */    
}

static void host_blas_dtpsv(
    NLBlas_t blas, MatrixTriangle uplo, MatrixTranspose trans,
    MatrixUnitTriangular diag, int n, const double *AP,
    double *x, int incx 
) {
    static const char *UL[2] = { "U", "L" };
    static const char *T[3]  = { "N", "T", 0 };
    static const char *D[2]  = { "U", "N" };
    nl_arg_used(blas);    
    NL_FORTRAN_WRAP(dtpsv)(
	UL[(int)uplo],T[(int)trans],D[(int)diag],&n,(double*)AP,x,&incx
    );
    /* TODO: update flops */
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
	blas.Dcopy = host_blas_dcopy;
	blas.Ddot = host_blas_ddot;
	blas.Dnrm2 = host_blas_dnrm2;
	blas.Daxpy = host_blas_daxpy;
	blas.Dscal = host_blas_dscal;
	blas.Dgemv = host_blas_dgemv;
	blas.Dtpsv = host_blas_dtpsv;
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



static NLuint nlSolveSystem_CG(
    NLBlas_t blas,
    NLMatrix M, NLdouble* b, NLdouble* x,
    double eps, NLuint max_iter
) {
    NLint N = (NLint)M->m;

    NLdouble *g = NL_NEW_VECTOR(blas, NL_DEVICE_MEMORY, N);
    NLdouble *r = NL_NEW_VECTOR(blas, NL_DEVICE_MEMORY, N); 
    NLdouble *p = NL_NEW_VECTOR(blas, NL_DEVICE_MEMORY, N);
    NLuint its=0;
    NLdouble t, tau, sig, rho, gam;
    NLdouble b_square=blas->Ddot(blas,N,b,1,b,1);
    NLdouble err=eps*eps*b_square;
    NLdouble curr_err;

    nlMultMatrixVector(M,x,g);
    blas->Daxpy(blas,N,-1.,b,1,g,1);
    blas->Dscal(blas,N,-1.,g,1);
    blas->Dcopy(blas,N,g,1,r,1);
    curr_err = blas->Ddot(blas,N,g,1,g,1);
    while ( curr_err >err && its < max_iter) {
	if(nlCurrentContext != NULL) {
	    if(nlCurrentContext->progress_func != NULL) {
		nlCurrentContext->progress_func(its, max_iter, curr_err, err);
	    }
	}
	nlMultMatrixVector(M,r,p);
        rho=blas->Ddot(blas,N,p,1,p,1);
        sig=blas->Ddot(blas,N,r,1,p,1);
        tau=blas->Ddot(blas,N,g,1,r,1);
        t=tau/sig;
        blas->Daxpy(blas,N,t,r,1,x,1);
        blas->Daxpy(blas,N,-t,p,1,g,1);
        gam=(t*t*rho-tau)/tau;
        blas->Dscal(blas,N,gam,r,1);
        blas->Daxpy(blas,N,1.,g,1,r,1);
        ++its;
        curr_err = blas->Ddot(blas,N,g,1,g,1);
    }
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, N, g);
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, N, r);
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, N, p);
    blas->sq_bnorm = b_square;
    blas->sq_rnorm = curr_err;
    return its;
}

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
    blas->Dcopy(blas,N,d,1,h,1);
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
    NLenum solver,
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

    switch(solver) {
	case NL_CG:
	    if(P == NULL) {
		result = nlSolveSystem_CG(blas,M,b,x,eps,max_iter);
	    } else {
		result = nlSolveSystem_PRE_CG(blas,M,P,b,x,eps,max_iter);
	    }
	    break;
	default:
	    nl_assert_not_reached;
    }


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
    switch(nlCurrentContext->matrix_mode) {
	case NL_STIFFNESS_MATRIX: {
	    nl_assert(nlCurrentContext->M != NULL);	    
	    nl_assert(nlCurrentContext->M->type == NL_MATRIX_SPARSE_DYNAMIC);
	    result = (NLSparseMatrix*)(nlCurrentContext->M);
	} break;
	default:
	    nl_assert_not_reached;
    }
    return result;
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
    case NL_SOLVER: {
        nlCurrentContext->solver = (NLenum)param;
    } break;
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
    case NL_SOLVER: {
        *params = (NLint)(nlCurrentContext->solver);
    } break;
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
    case NL_NNZ: {
        *params = (NLint)(nlMatrixNNZ(nlCurrentContext->M));
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

    /*
     * If the user trusts OpenNL and has left solver as NL_SOLVER_DEFAULT,
     * then we setup reasonable parameters for him.
     */
    if(nlCurrentContext->solver == NL_SOLVER_DEFAULT) {
        nlCurrentContext->solver = NL_CG;
        if(!nlCurrentContext->preconditioner_defined) {
            nlCurrentContext->preconditioner = NL_PRECOND_JACOBI;
        }
        if(!nlCurrentContext->max_iterations_defined) {
            nlCurrentContext->max_iterations = n*5;
        }
        if(!nlCurrentContext->threshold_defined) {
            nlCurrentContext->threshold = 1e-6;
        }
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
	    nlCurrentContext->matrix_mode == NL_STIFFNESS_MATRIX &&
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
    result = nlCurrentContext->solver_func();
    nlVectorToVariables();
    nlCurrentContext->elapsed_time = nlCurrentTime() - nlCurrentContext->start_time;
    nlTransition(NL_STATE_SYSTEM_CONSTRUCTED, NL_STATE_SOLVED);
    return result;
}

/* Eigen solver */

void nlMatrixMode(NLenum matrix) {
    NLuint n = 0;
    nl_assert(
	nlCurrentContext->state == NL_STATE_SYSTEM ||
	nlCurrentContext->state == NL_STATE_MATRIX_CONSTRUCTED
    );
    nlCurrentContext->state = NL_STATE_SYSTEM;
    nlCurrentContext->matrix_mode = matrix;
    nlCurrentContext->current_row = 0;
    switch(matrix) {
	case NL_STIFFNESS_MATRIX: {
	    /* Stiffness matrix is already constructed. */
	} break ;
	default:
	    nl_assert_not_reached;
    }
}


void nlEigenSolverParameterd(
    NLenum pname, NLdouble val
) {
    switch(pname) {
	case NL_EIGEN_SHIFT: {
	    nlCurrentContext->eigen_shift =  val;
	} break;
	case NL_EIGEN_THRESHOLD: {
	    nlSolverParameterd(pname, val);
	} break;
	default:
	    nl_assert_not_reached;
    }
}

void nlEigenSolverParameteri(
    NLenum pname, NLint val
) {
    switch(pname) {
	case NL_EIGEN_SOLVER: {
	    nlCurrentContext->eigen_solver = (NLenum)val;
	} break;
	case NL_NB_VARIABLES:	    
	case NL_NB_EIGENS:
	case NL_EIGEN_MAX_ITERATIONS: {
	    nlSolverParameteri(pname, val);
	} break;
	default:
	    nl_assert_not_reached;
    }
}

void nlEigenSolve() {
    if(nlCurrentContext->eigen_value == NULL) {
	nlCurrentContext->eigen_value = NL_NEW_ARRAY(
	    NLdouble,nlCurrentContext->nb_systems
	);
    }
    
    nlMatrixCompress(&nlCurrentContext->M);
    if(nlCurrentContext->B != NULL) {
	nlMatrixCompress(&nlCurrentContext->B);
    }
    
    switch(nlCurrentContext->eigen_solver) {
	default:
	    nl_assert_not_reached;
    }
}

double nlGetEigenValue(NLuint i) {
    nl_debug_assert(i < nlCurrentContext->nb_variables);
    return nlCurrentContext->eigen_value[i];
}
