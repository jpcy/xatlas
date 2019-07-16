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

typedef void* NLdll;


#define NL_LINK_NOW    1

#define NL_LINK_LAZY   2

#define NL_LINK_GLOBAL 4

#define NL_LINK_QUIET  8

#define NL_LINK_USE_FALLBACK 16

NLdll nlOpenDLL(const char* filename, NLenum flags);

void nlCloseDLL(NLdll handle);

NLfunc nlFindFunction(NLdll handle, const char* funcname);


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
#define NL_MATRIX_SUPERLU_EXT    0x1003    
#define NL_MATRIX_CHOLMOD_EXT    0x1004    
#define NL_MATRIX_FUNCTION       0x1005
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

    NLboolean symmetric_storage;
} NLCRSMatrix;

NLAPI void NLAPIENTRY nlCRSMatrixConstruct(
    NLCRSMatrix* M, NLuint m, NLuint n, NLuint nnz, NLuint nslices
);

NLAPI void NLAPIENTRY nlCRSMatrixConstructSymmetric(
    NLCRSMatrix* M, NLuint n, NLuint nnz
);
    
NLAPI NLboolean NLAPIENTRY nlCRSMatrixLoad(
    NLCRSMatrix* M, const char* filename
);

NLAPI NLboolean NLAPIENTRY nlCRSMatrixSave(
    NLCRSMatrix* M, const char* filename
);

NLAPI NLuint NLAPIENTRY nlCRSMatrixNNZ(NLCRSMatrix* M);
    

/* SparseMatrix data structure */

#define NL_MATRIX_STORE_ROWS          1

#define NL_MATRIX_STORE_COLUMNS       2

#define NL_MATRIX_STORE_SYMMETRIC     4
    
typedef struct {
    NLuint m;
    
    NLuint n;

    NLenum type;
    
    NLDestroyMatrixFunc destroy_func;

    NLMultMatrixVectorFunc mult_func;

    
    NLuint diag_size;

    NLuint diag_capacity;
    
    NLenum storage;

    NLRowColumn* row;

    NLRowColumn* column;

    NLdouble*    diag;

    NLuint row_capacity;

    NLuint column_capacity;
    
} NLSparseMatrix;


NLAPI NLMatrix NLAPIENTRY nlSparseMatrixNew(
    NLuint m, NLuint n, NLenum storage
);

NLAPI void NLAPIENTRY nlSparseMatrixConstruct(
    NLSparseMatrix* M, NLuint m, NLuint n, NLenum storage
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

NLAPI NLMatrix NLAPIENTRY nlCRSMatrixNewFromSparseMatrixSymmetric(
    NLSparseMatrix* M
);    

    
NLAPI void NLAPIENTRY nlMatrixCompress(NLMatrix* M);

NLAPI NLuint NLAPIENTRY nlMatrixNNZ(NLMatrix M);

NLAPI NLMatrix NLAPIENTRY nlMatrixFactorize(NLMatrix M, NLenum solver);
    


    typedef void(*NLMatrixFunc)(const double* x, double* y);

NLAPI NLMatrix NLAPIENTRY nlMatrixNewFromFunction(
    NLuint m, NLuint n, NLMatrixFunc func
);	     

NLAPI NLMatrixFunc NLAPIENTRY nlMatrixGetFunction(NLMatrix M);



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

    NLboolean        user_variable_buffers;
    
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

    NLdouble*        right_hand_side;

    NLdouble         row_scaling;

    NLenum           solver;

    NLenum           preconditioner;

    NLboolean        preconditioner_defined;
    
    NLuint           nb_variables;

    NLuint           nb_systems;

    NLboolean        ij_coefficient_called;
    
    NLuint           current_row;

    NLboolean        least_squares;

    NLboolean        symmetric;

    NLuint           max_iterations;


    NLboolean        max_iterations_defined;
    
    NLuint           inner_iterations;

    NLdouble         threshold;

    NLboolean        threshold_defined;
    
    NLdouble         omega;

    NLboolean        normalize_rows;
    
    NLuint           used_iterations;

    NLdouble         error;


    NLdouble         start_time;
    
    NLdouble         elapsed_time;

    NLSolverFunc     solver_func;

    NLProgressFunc   progress_func;

    NLboolean        verbose;

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

NLMatrix nlNewSSORPreconditioner(NLMatrix M, double omega);

#endif

/******* extracted from nl_superlu.h *******/

#ifndef OPENNL_SUPERLU_H
#define OPENNL_SUPERLU_H



NLAPI NLMatrix NLAPIENTRY nlMatrixFactorize_SUPERLU(
    NLMatrix M, NLenum solver
);

NLboolean nlInitExtension_SUPERLU(void);

NLboolean nlExtensionIsInitialized_SUPERLU(void);

#endif

/******* extracted from nl_cholmod.h *******/

#ifndef OPENNL_CHOLMOD_H
#define OPENNL_CHOLMOD_H



NLAPI NLMatrix NLAPIENTRY nlMatrixFactorize_CHOLMOD(
    NLMatrix M, NLenum solver
);

NLboolean nlInitExtension_CHOLMOD(void);

NLboolean nlExtensionIsInitialized_CHOLMOD(void);

#endif

/******* extracted from nl_arpack.h *******/

#ifndef OPENNL_ARPACK_H
#define OPENNL_ARPACK_H



NLboolean nlInitExtension_ARPACK(void);

NLboolean nlExtensionIsInitialized_ARPACK(void);

void nlEigenSolve_ARPACK(void);


#endif

/******* extracted from nl_mkl.h *******/

#ifndef OPENNL_MKL_H
#define OPENNL_MKL_H



NLboolean nlInitExtension_MKL(void);

NLboolean nlExtensionIsInitialized_MKL(void);

extern NLMultMatrixVectorFunc NLMultMatrixVector_MKL;

#endif

/******* extracted from nl_cuda.h *******/

#ifndef OPENNL_CUDA_EXT_H
#define OPENNL_CUDA_EXT_H



NLboolean nlInitExtension_CUDA(void);

NLboolean nlExtensionIsInitialized_CUDA(void);




NLMatrix nlCUDAMatrixNewFromCRSMatrix(NLMatrix M);

NLMatrix nlCUDAJacobiPreconditionerNewFromCRSMatrix(NLMatrix M);

NLBlas_t nlCUDABlas(void);



#endif

/******* extracted from nl_os.c *******/


#if (defined (WIN32) || defined(_WIN64))
#include <windows.h>
#else
#include <sys/types.h>
#include <sys/times.h> 
#endif

#if defined(GEO_DYNAMIC_LIBS) && defined(NL_OS_UNIX)
#include <dlfcn.h>
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


/* DLLs/shared objects/dylibs */

#if defined(GEO_DYNAMIC_LIBS) 

#  if defined(NL_OS_UNIX)

NLdll nlOpenDLL(const char* name, NLenum flags_in) {
    void* result = NULL;
    int flags = 0;
    if((flags_in & NL_LINK_NOW) != 0) {
	flags |= RTLD_NOW;
    }
    if((flags_in & NL_LINK_LAZY) != 0) {
	flags |= RTLD_LAZY;
    }
    if((flags_in & NL_LINK_GLOBAL) != 0) {
	flags |= RTLD_GLOBAL;
    }
    if((flags_in & NL_LINK_QUIET) == 0) {
	nl_fprintf(stdout,"Trying to load %s\n", name);
    }
    result = dlopen(name, flags);
    if(result == NULL) {
	if((flags_in & NL_LINK_QUIET) == 0) {	
	    nl_fprintf(stderr,"Did not find %s,\n", name);
	    nl_fprintf(stderr,"Retrying with libgeogram_num_3rdparty.so\n");
	}
	if((flags_in & NL_LINK_USE_FALLBACK) != 0) {
	    result=dlopen("libgeogram_num_3rdparty.so", flags);
	    if(result == NULL) {
		if((flags_in & NL_LINK_QUIET) == 0) {		    
		    nlError("nlOpenDLL/dlopen",dlerror());
		}
	    }
        }
    }
    if((flags_in & NL_LINK_QUIET) == 0 && result != NULL) {
	nl_fprintf(stdout,"Loaded %s\n", name);
    }
    
    return result;
}

void nlCloseDLL(void* handle) {
    dlclose(handle);
}

NLfunc nlFindFunction(void* handle, const char* name) {
    /*
     * It is not legal in modern C to cast a void*
     *  pointer into a function pointer, thus requiring this
     *  (quite dirty) function that uses a union.    
     */
    union {
        void* ptr;
        NLfunc fptr;
    } u;
    u.ptr = dlsym(handle, name);
    return u.fptr;
}

#  elif defined(NL_OS_WINDOWS)

NLdll nlOpenDLL(const char* name, NLenum flags) {
    /* Note: NL_LINK_LAZY and NL_LINK_GLOBAL are ignored. */
    void* result = LoadLibrary(name);
    if(result == NULL && ((flags & NL_LINK_USE_FALLBACK) != 0)) {
	if((flags & NL_LINK_QUIET) == 0) {
	    nl_fprintf(stderr,"Did not find %s,\n", name);
	    nl_fprintf(stderr,"Retrying with geogram_num_3rdparty\n");
	}
        result=LoadLibrary("geogram_num_3rdparty.dll");
    }
    return result;
}

void nlCloseDLL(void* handle) {
    FreeLibrary((HMODULE)handle);
}

NLfunc nlFindFunction(void* handle, const char* name) {
    return (NLfunc)GetProcAddress((HMODULE)handle, name);
}

#  endif

#else

NLdll nlOpenDLL(const char* name, NLenum flags) {
    nl_arg_used(name);
    nl_arg_used(flags);
#ifdef NL_OS_UNIX
    nlError("nlOpenDLL","Was not compiled with dynamic linking enabled");
    nlError("nlOpenDLL","(see VORPALINE_BUILD_DYNAMIC in CMakeLists.txt)");        
#else    
    nlError("nlOpenDLL","Not implemented");
#endif    
    return NULL;
}

void nlCloseDLL(void* handle) {
    nl_arg_used(handle);
    nlError("nlCloseDLL","Not implemented");        
}

NLfunc nlFindFunction(void* handle, const char* name) {
    nl_arg_used(handle);
    nl_arg_used(name);
    nlError("nlFindFunction","Not implemented");            
    return NULL;
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
    NLuint i,j,jj;
    NLdouble a;
    
    if(M->symmetric_storage) {
        for(i=0; i<M->m; ++i) {
            y[i] = 0.0;
        }
        for(i=0; i<M->m; ++i) {
            for(jj=M->rowptr[i]; jj<M->rowptr[i+1]; ++jj) {
                a = M->val[jj];
                j = M->colind[jj];
                y[i] += a * x[j];
                if(j != i) {
                    y[j] += a * x[i];
                }
            }
        }
    } else {
    
#if defined(_OPENMP)
#pragma omp parallel for private(slice)
#endif
    
	for(slice=0; slice<nslices; ++slice) {
	    nlCRSMatrixMultSlice(
		M,x,y,M->sliceptr[slice],M->sliceptr[slice+1]
	    );
	}
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
    if(NLMultMatrixVector_MKL != NULL) {
	M->mult_func = (NLMultMatrixVectorFunc)NLMultMatrixVector_MKL;
    } else {
	M->mult_func = (NLMultMatrixVectorFunc)nlCRSMatrixMult;
    }
    M->nslices = nslices;
    M->val = NL_NEW_ARRAY(double, nnz);
    M->rowptr = NL_NEW_ARRAY(NLuint, m+1);
    M->colind = NL_NEW_ARRAY(NLuint, nnz);
    M->sliceptr = NL_NEW_ARRAY(NLuint, nslices+1);
    M->symmetric_storage = NL_FALSE;
}

void nlCRSMatrixConstructSymmetric(
    NLCRSMatrix* M, NLuint n, NLuint nnz
) {
    M->m = n;
    M->n = n;
    M->type = NL_MATRIX_CRS;
    M->destroy_func = (NLDestroyMatrixFunc)nlCRSMatrixDestroy;
    M->mult_func = (NLMultMatrixVectorFunc)nlCRSMatrixMult;
    M->nslices = 0;
    M->val = NL_NEW_ARRAY(double, nnz);
    M->rowptr = NL_NEW_ARRAY(NLuint, n+1);
    M->colind = NL_NEW_ARRAY(NLuint, nnz);
    M->sliceptr = NULL;
    M->symmetric_storage = NL_TRUE;
}


/* SparseMatrix data structure */


static void nlSparseMatrixDestroyRowColumns(NLSparseMatrix* M) {
    NLuint i;
    if(M->storage & NL_MATRIX_STORE_ROWS) {
        for(i=0; i<M->m; i++) {
            nlRowColumnDestroy(&(M->row[i]));
        }
        NL_DELETE_ARRAY(M->row);
    }
    M->storage = (NLenum)((int)(M->storage) & ~NL_MATRIX_STORE_ROWS);
    
    if(M->storage & NL_MATRIX_STORE_COLUMNS) {
        for(i=0; i<M->n; i++) {
            nlRowColumnDestroy(&(M->column[i]));
        }
        NL_DELETE_ARRAY(M->column);
    }
    M->storage = (NLenum)((int)(M->storage) & ~NL_MATRIX_STORE_COLUMNS);    
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
    if((M->storage & NL_MATRIX_STORE_SYMMETRIC) && (j > i)) {
        return;
    }
    if(i == j) {
        M->diag[i] += value;
    }
    if(M->storage & NL_MATRIX_STORE_ROWS) {
        nlRowColumnAdd(&(M->row[i]), j, value);
    }
    if(M->storage & NL_MATRIX_STORE_COLUMNS) {
        nlRowColumnAdd(&(M->column[j]), i, value);
    }
}

static void nlSparseMatrixAddSparseMatrix(
    NLSparseMatrix* M, double mul, const NLSparseMatrix* N    
) {
    NLuint i,j,ii,jj;
    nl_assert(M->m == N->m);
    nl_assert(M->n == N->n);
    if(N->storage & NL_MATRIX_STORE_SYMMETRIC) {
	nl_assert(M->storage & NL_MATRIX_STORE_SYMMETRIC);
    }
    if(N->storage & NL_MATRIX_STORE_ROWS) {
	for(i=0; i<N->m; ++i) {
	    for(jj=0; jj<N->row[i].size; ++jj) {
		nlSparseMatrixAdd(
		    M,
		    i, N->row[i].coeff[jj].index,
		    mul*N->row[i].coeff[jj].value
		);
	    }
	}
    } else {
	nl_assert(N->storage & NL_MATRIX_STORE_COLUMNS);	
	for(j=0; j<N->n; ++j) {
	    for(ii=0; ii<N->column[j].size; ++ii) {
		nlSparseMatrixAdd(
		    M,
		    N->column[j].coeff[ii].index, j,
		    mul*N->column[j].coeff[ii].value
		);
	    }
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
    if(M->storage & NL_MATRIX_STORE_ROWS) {
        for(i=0; i<M->m; i++) {
            nlRowColumnZero(&(M->row[i]));
        }
    }
    if(M->storage & NL_MATRIX_STORE_COLUMNS) {
        for(i=0; i<M->n; i++) {
            nlRowColumnZero(&(M->column[i]));
        }
    }
    NL_CLEAR_ARRAY(NLdouble, M->diag, M->diag_size);
}

void nlSparseMatrixClear( NLSparseMatrix* M) {
    NLuint i;
    if(M->storage & NL_MATRIX_STORE_ROWS) {
        for(i=0; i<M->m; i++) {
            nlRowColumnClear(&(M->row[i]));
        }
    }
    if(M->storage & NL_MATRIX_STORE_COLUMNS) {
        for(i=0; i<M->n; i++) {
            nlRowColumnClear(&(M->column[i]));
        }
    }
    NL_CLEAR_ARRAY(NLdouble, M->diag, M->diag_size);
}

/* Returns the number of non-zero coefficients */
NLuint nlSparseMatrixNNZ( NLSparseMatrix* M) {
    NLuint nnz = 0;
    NLuint i;
    if(M->storage & NL_MATRIX_STORE_ROWS) {
        for(i = 0; i<M->m; i++) {
            nnz += M->row[i].size;
        }
    } else if (M->storage & NL_MATRIX_STORE_COLUMNS) {
        for(i = 0; i<M->n; i++) {
            nnz += M->column[i].size;
        }
    } else {
        nl_assert_not_reached;
    }
    return nnz;
}

void nlSparseMatrixSort( NLSparseMatrix* M) {
    NLuint i;
    if(M->storage & NL_MATRIX_STORE_ROWS) {
        for(i = 0; i<M->m; i++) {
            nlRowColumnSort(&(M->row[i]));                
        }
    } 
    if (M->storage & NL_MATRIX_STORE_COLUMNS) {
        for(i = 0; i<M->n; i++) {
            nlRowColumnSort(&(M->column[i]));
        }
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

    nl_assert(M->storage & NL_MATRIX_STORE_ROWS);
    nl_assert(!(M->storage & NL_MATRIX_STORE_COLUMNS));    
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

static void nlSparseMatrix_mult_rows_symmetric(
    NLSparseMatrix* A,
    const NLdouble* x,
    NLdouble* y
) {
    NLuint m = A->m;
    NLuint i,ij;
    NLCoeff* c = NULL;
    for(i=0; i<m; i++) {
        NLRowColumn* Ri = &(A->row[i]);
        y[i] = 0;
        for(ij=0; ij<Ri->size; ++ij) {
            c = &(Ri->coeff[ij]);
            y[i] += c->value * x[c->index];
            if(i != c->index) {
                y[c->index] += c->value * x[i];
            }
        }
    }
}

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

static void nlSparseMatrix_mult_cols_symmetric(
        NLSparseMatrix* A,
        const NLdouble* x,
        NLdouble* y
) {
    NLuint n = A->n;
    NLuint j,ii;
    NLCoeff* c = NULL;
    for(j=0; j<n; j++) {
        NLRowColumn* Cj = &(A->column[j]);       
        y[j] = 0;
        for(ii=0; ii<Cj->size; ii++) {
            c = &(Cj->coeff[ii]);
            y[c->index] += c->value * x[j];
            if(j != c->index) {
                y[j] += c->value * x[c->index];
            }
        }
    }
}

static void nlSparseMatrix_mult_cols(
        NLSparseMatrix* A,
        const NLdouble* x,
        NLdouble* y
) {
    NLuint n = A->n;
    NLuint j,ii; 
    NLCoeff* c = NULL;
    NL_CLEAR_ARRAY(NLdouble, y, A->m);
    for(j=0; j<n; j++) {
        NLRowColumn* Cj = &(A->column[j]);
        for(ii=0; ii<Cj->size; ii++) {
            c = &(Cj->coeff[ii]);
            y[c->index] += c->value * x[j];
        }
    }
}

void nlSparseMatrixMult(
    NLSparseMatrix* A, const NLdouble* x, NLdouble* y
) {
    nl_assert(A->type == NL_MATRIX_SPARSE_DYNAMIC);
    if(A->storage & NL_MATRIX_STORE_ROWS) {
        if(A->storage & NL_MATRIX_STORE_SYMMETRIC) {
            nlSparseMatrix_mult_rows_symmetric(A, x, y);
        } else {
            nlSparseMatrix_mult_rows(A, x, y);
        }
    } else {
        if(A->storage & NL_MATRIX_STORE_SYMMETRIC) {
            nlSparseMatrix_mult_cols_symmetric(A, x, y);
        } else {
            nlSparseMatrix_mult_cols(A, x, y);
        }
    }
    nlHostBlas()->flops += (NLulong)(2*nlSparseMatrixNNZ(A));
}

NLMatrix nlSparseMatrixNew(
    NLuint m, NLuint n, NLenum storage
) {
    NLSparseMatrix* result = NL_NEW(NLSparseMatrix);
    nlSparseMatrixConstruct(result, m, n, storage);
    return (NLMatrix)result;
}

void nlSparseMatrixConstruct(
    NLSparseMatrix* M, NLuint m, NLuint n, NLenum storage
) {
    NLuint i;
    M->m = m;
    M->n = n;
    M->type = NL_MATRIX_SPARSE_DYNAMIC;
    M->destroy_func = (NLDestroyMatrixFunc)nlSparseMatrixDestroy;
    M->mult_func = (NLMultMatrixVectorFunc)nlSparseMatrixMult;
    M->storage = storage;
    if(storage & NL_MATRIX_STORE_ROWS) {
        M->row = NL_NEW_ARRAY(NLRowColumn, m);
	M->row_capacity = m;
        for(i=0; i<n; i++) {
            nlRowColumnConstruct(&(M->row[i]));
        }
    } else {
        M->row = NULL;
	M->row_capacity = 0;
    }

    if(storage & NL_MATRIX_STORE_COLUMNS) {
        M->column = NL_NEW_ARRAY(NLRowColumn, n);
	M->column_capacity = n;
        for(i=0; i<n; i++) {
            nlRowColumnConstruct(&(M->column[i]));
        }
    } else {
        M->column = NULL;
	M->column_capacity = 0;
    }

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
    if(M->storage & NL_MATRIX_STORE_ROWS) {
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
    }
    adjust_diag(M);
}

void nlSparseMatrixAddColumn( NLSparseMatrix* M) {
    ++M->n;
    if(M->storage & NL_MATRIX_STORE_COLUMNS) {
	if(M->n > M->column_capacity) {
	    M->column_capacity *= 2;
	    if(M->column_capacity == 0) {
		M->column_capacity = 16;
	    }
	    M->column = NL_RENEW_ARRAY(
		NLRowColumn, M->column, M->column_capacity
	    );
	}
	nlRowColumnConstruct(&(M->column[M->n-1]));
    }
    adjust_diag(M);
}




NLMatrix nlCRSMatrixNewFromSparseMatrix(NLSparseMatrix* M) {
    NLuint nnz = nlSparseMatrixNNZ(M);
    NLuint nslices = 8; /* TODO: get number of cores */
    NLuint slice, cur_bound, cur_NNZ, cur_row;
    NLuint i,ij,k; 
    NLuint slice_size = nnz / nslices;
    NLCRSMatrix* CRS = NL_NEW(NLCRSMatrix);

    nl_assert(M->storage & NL_MATRIX_STORE_ROWS);

    if(M->storage & NL_MATRIX_STORE_SYMMETRIC) {
        nl_assert(M->m == M->n);
        nlCRSMatrixConstructSymmetric(CRS, M->n, nnz);        
    } else {
        nlCRSMatrixConstruct(CRS, M->m, M->n, nnz, nslices);
    }
    
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

NLMatrix nlCRSMatrixNewFromSparseMatrixSymmetric(NLSparseMatrix* M) {
    NLuint nnz;
    NLuint i,j,jj,k;
    NLCRSMatrix* CRS = NL_NEW(NLCRSMatrix);
    
    nl_assert(M->storage & NL_MATRIX_STORE_ROWS);
    nl_assert(M->m == M->n);

    nlSparseMatrixSort(M);
    
    if(M->storage & NL_MATRIX_STORE_SYMMETRIC) {
        nnz = nlSparseMatrixNNZ(M);
    } else {
        nnz = 0;
        for(i=0; i<M->n; ++i) {
            NLRowColumn* Ri = &M->row[i];
            for(jj=0; jj<Ri->size; ++jj) {
                j = Ri->coeff[jj].index;
                if(j <= i) {
                    ++nnz;
                }
            }
        }
    }

    nlCRSMatrixConstructSymmetric(CRS, M->n, nnz);        

    k=0;
    for(i=0; i<M->m; ++i) {
        NLRowColumn* Ri = &(M->row[i]);
        CRS->rowptr[i] = k;
        for(jj=0; jj<Ri->size; ++jj) {
            j = Ri->coeff[jj].index;
            if((M->storage & NL_MATRIX_STORE_SYMMETRIC)) {
                nl_debug_assert(j <= i);
            }
            if(j <= i) {
                CRS->val[k] = Ri->coeff[jj].value;
                CRS->colind[k] = j;
                ++k;
            }
        }
    }
    CRS->rowptr[M->m] = k;

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
	case NL_SUPERLU_EXT:
	case NL_PERM_SUPERLU_EXT:      
	case NL_SYMMETRIC_SUPERLU_EXT:
	    result = nlMatrixFactorize_SUPERLU(M,solver);
	    break;
	case NL_CHOLMOD_EXT:
	    result = nlMatrixFactorize_CHOLMOD(M,solver);	    
	    break;
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
} NLFunctionMatrix;

static void nlFunctionMatrixDestroy(NLFunctionMatrix* M) {
    (void)M; /* to avoid 'unused parameter' warning */
    /* 
     * Nothing special to do, 
     * there is no dynamic allocated mem.
     */
}

static void nlFunctionMatrixMult(
    NLFunctionMatrix* M, const NLdouble* x, NLdouble* y
) {
    M->matrix_func(x,y);
}

NLMatrix nlMatrixNewFromFunction(NLuint m, NLuint n, NLMatrixFunc func) {
    NLFunctionMatrix* result = NL_NEW(NLFunctionMatrix);
    result->m = m;
    result->n = n;
    result->type = NL_MATRIX_FUNCTION;
    result->destroy_func = (NLDestroyMatrixFunc)nlFunctionMatrixDestroy;
    result->mult_func = (NLMultMatrixVectorFunc)nlFunctionMatrixMult;
    result->matrix_func = func;
    return (NLMatrix)result;
}

NLMatrixFunc nlMatrixGetFunction(NLMatrix M) {
    if(M == NULL) {
	return NULL;
    }
    if(M->type != NL_MATRIX_FUNCTION) {
	return NULL;
    }
    return ((NLFunctionMatrix*)M)->matrix_func;
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
    result->row_scaling         = 1.0;
    result->inner_iterations    = 5;
    result->solver_func         = nlDefaultSolver;
    result->progress_func       = NULL;
    result->verbose             = NL_FALSE;
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
    NL_DELETE_ARRAY(context->right_hand_side);

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
    /* Check compatibility between solver and preconditioner */
    if(
        nlCurrentContext->solver == NL_BICGSTAB && 
        nlCurrentContext->preconditioner == NL_PRECOND_SSOR
    ) {
        nlWarning(
            "nlSolve", 
            "cannot use SSOR preconditioner with non-symmetric matrix, "
	    "switching to Jacobi"
        );
        nlCurrentContext->preconditioner = NL_PRECOND_JACOBI;        
    }
    if(
        nlCurrentContext->solver == NL_GMRES && 
        nlCurrentContext->preconditioner != NL_PRECOND_NONE
    ) {
        nlWarning("nlSolve", "Preconditioner not implemented yet for GMRES");
        nlCurrentContext->preconditioner = NL_PRECOND_NONE;        
    }
    if(
        nlCurrentContext->solver == NL_SUPERLU_EXT && 
        nlCurrentContext->preconditioner != NL_PRECOND_NONE
    ) {
        nlWarning("nlSolve", "Preconditioner not implemented yet for SUPERLU");
        nlCurrentContext->preconditioner = NL_PRECOND_NONE;        
    }
    if(
        nlCurrentContext->solver == NL_CHOLMOD_EXT && 
        nlCurrentContext->preconditioner != NL_PRECOND_NONE
    ) {
        nlWarning("nlSolve", "Preconditioner not implemented yet for CHOLMOD");
        nlCurrentContext->preconditioner = NL_PRECOND_NONE;        
    }
    if(
        nlCurrentContext->solver == NL_PERM_SUPERLU_EXT && 
        nlCurrentContext->preconditioner != NL_PRECOND_NONE
    ) {
        nlWarning(
	    "nlSolve", "Preconditioner not implemented yet for PERMSUPERLU"
	);
        nlCurrentContext->preconditioner = NL_PRECOND_NONE;        
    }
    if(
        nlCurrentContext->solver == NL_SYMMETRIC_SUPERLU_EXT && 
        nlCurrentContext->preconditioner != NL_PRECOND_NONE
    ) {
        nlWarning(
	    "nlSolve", "Preconditioner not implemented yet for PERMSUPERLU"
	);
        nlCurrentContext->preconditioner = NL_PRECOND_NONE;        
    }

    nlDeleteMatrix(nlCurrentContext->P);
    nlCurrentContext->P = NULL;
    
    switch(nlCurrentContext->preconditioner) {
    case NL_PRECOND_NONE:
        break;
    case NL_PRECOND_JACOBI:
	nlCurrentContext->P = nlNewJacobiPreconditioner(nlCurrentContext->M);
        break;
    case NL_PRECOND_SSOR:
	nlCurrentContext->P = nlNewSSORPreconditioner(
	    nlCurrentContext->M,nlCurrentContext->omega
	);	
        break;
    case NL_PRECOND_USER:
        break;
    default:
        nl_assert_not_reached;
    }

    if(nlCurrentContext->preconditioner != NL_PRECOND_SSOR) {
        if(getenv("NL_LOW_MEM") == NULL) {
            nlMatrixCompress(&nlCurrentContext->M);
        }
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
    NLboolean use_CUDA = NL_FALSE;
    NLdouble* b = nlCurrentContext->b;
    NLdouble* x = nlCurrentContext->x;
    NLuint n = nlCurrentContext->n;
    NLuint k;
    NLBlas_t blas = nlHostBlas();
    NLMatrix M = nlCurrentContext->M;
    NLMatrix P = nlCurrentContext->P;
    
    /*
     * For CUDA: it is implemented for
     *   all iterative solvers except GMRES
     *   Jacobi preconditioner
     */
    if(nlExtensionIsInitialized_CUDA() &&
       (nlCurrentContext->solver != NL_GMRES) && 
       (nlCurrentContext->preconditioner == NL_PRECOND_NONE ||
	nlCurrentContext->preconditioner == NL_PRECOND_JACOBI)
    ) {
	if(nlCurrentContext->verbose) { 
	    nl_printf("Using CUDA\n");
	} 
	use_CUDA = NL_TRUE;
	blas = nlCUDABlas();
	if(nlCurrentContext->preconditioner == NL_PRECOND_JACOBI) {
	    P = nlCUDAJacobiPreconditionerNewFromCRSMatrix(M);
	}
	M = nlCUDAMatrixNewFromCRSMatrix(M);
    }

    /* 
     * We do not count CUDA transfers and CUDA matrix construction
     * when estimating GFlops
     */
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
    
    if(use_CUDA) {
	nlDeleteMatrix(M);
	nlDeleteMatrix(P);
    }
    
    return NL_TRUE;
}



NLboolean nlDefaultSolver() {
    NLboolean result = NL_TRUE;
    nlSetupPreconditioner();
    switch(nlCurrentContext->solver) {
	case NL_CG:
	case NL_BICGSTAB:
	case NL_GMRES: {
	    result = nlSolveIterative();
	} break;

	case NL_SUPERLU_EXT: 
	case NL_PERM_SUPERLU_EXT: 
	case NL_SYMMETRIC_SUPERLU_EXT: 
	case NL_CHOLMOD_EXT: {
	    result = nlSolveDirect();
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

#ifdef NL_USE_SUPERLU
#ifndef NL_USE_BLAS
#define NL_USE_BLAS
/* 
 * The BLAS included in SuperLU does not have DTPSV,
 * we use the DTPSV embedded in OpenNL.
 */
#define NEEDS_DTPSV
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
	    if(nlCurrentContext->verbose && !(its % 100)) {
		nl_printf ( "%d : %.10e -- %.10e\n", its, curr_err, err );
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
	    if( nlCurrentContext->verbose && !(its % 100)) {
		nl_printf ( "%d : %.10e -- %.10e\n", its, curr_err, err );
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

static NLuint nlSolveSystem_BICGSTAB(
    NLBlas_t blas,
    NLMatrix M, NLdouble* b, NLdouble* x,
    double eps, NLuint max_iter
) {
    NLint     N   = (NLint)M->n;
    NLdouble *rT  = NL_NEW_VECTOR(blas, NL_DEVICE_MEMORY, N); 
    NLdouble *d   = NL_NEW_VECTOR(blas, NL_DEVICE_MEMORY, N); 
    NLdouble *h   = NL_NEW_VECTOR(blas, NL_DEVICE_MEMORY, N); 
    NLdouble *u   = NL_NEW_VECTOR(blas, NL_DEVICE_MEMORY, N); 
    NLdouble *Ad  = NL_NEW_VECTOR(blas, NL_DEVICE_MEMORY, N); 
    NLdouble *t   = NL_NEW_VECTOR(blas, NL_DEVICE_MEMORY, N); 
    NLdouble *s   = h;
    NLdouble rTh, rTAd, rTr, alpha, beta, omega, st, tt;
    NLuint its=0;
    NLdouble b_square = blas->Ddot(blas,N,b,1,b,1);
    NLdouble err=eps*eps*b_square;
    NLdouble *r = NL_NEW_VECTOR(blas, NL_DEVICE_MEMORY, N);
    nlMultMatrixVector(M,x,r);
    blas->Daxpy(blas,N,-1.,b,1,r,1);
    blas->Dcopy(blas,N,r,1,d,1);
    blas->Dcopy(blas,N,d,1,h,1);
    blas->Dcopy(blas,N,h,1,rT,1);
    nl_assert( blas->Ddot(blas,N,rT,1,rT,1)>1e-40 );
    rTh=blas->Ddot(blas,N,rT,1,h,1);
    rTr=blas->Ddot(blas,N,r,1,r,1);

    while ( rTr>err && its < max_iter) {
	if(nlCurrentContext != NULL) {
	    if(nlCurrentContext->progress_func != NULL) {
		nlCurrentContext->progress_func(its, max_iter, rTr, err);
	    }
	    if( (nlCurrentContext->verbose) && !(its % 100)) {
		nl_printf ( "%d : %.10e -- %.10e\n", its, rTr, err );
	    }
	}
	nlMultMatrixVector(M,d,Ad);
        rTAd=blas->Ddot(blas,N,rT,1,Ad,1);
        nl_assert( fabs(rTAd)>1e-40 );
        alpha=rTh/rTAd;
        blas->Daxpy(blas,N,-alpha,Ad,1,r,1);
        blas->Dcopy(blas,N,h,1,s,1);
        blas->Daxpy(blas,N,-alpha,Ad,1,s,1);
	nlMultMatrixVector(M,s,t);
        blas->Daxpy(blas,N,1.,t,1,u,1);
        blas->Dscal(blas,N,alpha,u,1);
        st=blas->Ddot(blas,N,s,1,t,1);
        tt=blas->Ddot(blas,N,t,1,t,1);
        if ( fabs(st)<1e-40 || fabs(tt)<1e-40 ) {
            omega = 0.;
        } else {
            omega = st/tt;
        }
        blas->Daxpy(blas,N,-omega,t,1,r,1);
        blas->Daxpy(blas,N,-alpha,d,1,x,1);
        blas->Daxpy(blas,N,-omega,s,1,x,1);
        blas->Dcopy(blas,N,s,1,h,1);
        blas->Daxpy(blas,N,-omega,t,1,h,1);
        beta=(alpha/omega)/rTh;
	rTh=blas->Ddot(blas,N,rT,1,h,1);
	beta*=rTh;
        blas->Dscal(blas,N,beta,d,1);
        blas->Daxpy(blas,N,1.,h,1,d,1);
        blas->Daxpy(blas,N,-beta*omega,Ad,1,d,1);
        rTr=blas->Ddot(blas,N,r,1,r,1);
        ++its;
    }
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, N, r);
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, N, rT);
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, N, d);
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, N, h);
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, N, u);
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, N, Ad);
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, N, t);
    blas->sq_bnorm = b_square;
    blas->sq_rnorm = rTr;
    return its;
}

static NLuint nlSolveSystem_PRE_BICGSTAB(
    NLBlas_t blas,
    NLMatrix M, NLMatrix P, NLdouble* b, NLdouble* x,
    double eps, NLuint max_iter
) {
    NLint     N   = (NLint)M->n;
    NLdouble *rT  = NL_NEW_VECTOR(blas, NL_DEVICE_MEMORY, N);
    NLdouble *d   = NL_NEW_VECTOR(blas, NL_DEVICE_MEMORY, N);
    NLdouble *h   = NL_NEW_VECTOR(blas, NL_DEVICE_MEMORY, N);
    NLdouble *u   = NL_NEW_VECTOR(blas, NL_DEVICE_MEMORY, N);
    NLdouble *Sd  = NL_NEW_VECTOR(blas, NL_DEVICE_MEMORY, N);
    NLdouble *t   = NL_NEW_VECTOR(blas, NL_DEVICE_MEMORY, N);
    NLdouble *aux = NL_NEW_VECTOR(blas, NL_DEVICE_MEMORY, N);
    NLdouble *s   = h;
    NLdouble rTh, rTSd, rTr, alpha, beta, omega, st, tt;
    NLuint its=0;
    NLdouble b_square = blas->Ddot(blas,N,b,1,b,1);
    NLdouble err  = eps*eps*b_square;
    NLdouble *r   = NL_NEW_VECTOR(blas, NL_DEVICE_MEMORY, N);

    nlMultMatrixVector(M,x,r);
    blas->Daxpy(blas,N,-1.,b,1,r,1);
    nlMultMatrixVector(P,r,d);
    blas->Dcopy(blas,N,d,1,h,1);
    blas->Dcopy(blas,N,h,1,rT,1);
    nl_assert( blas->Ddot(blas,N,rT,1,rT,1)>1e-40 );
    rTh=blas->Ddot(blas,N,rT,1,h,1);
    rTr=blas->Ddot(blas,N,r,1,r,1);

    while ( rTr>err && its < max_iter) {
	if(nlCurrentContext != NULL) {	
	    if(nlCurrentContext->progress_func != NULL) {
		nlCurrentContext->progress_func(its, max_iter, rTr, err);
	    }
	    if( (nlCurrentContext->verbose) && !(its % 100)) {
		nl_printf ( "%d : %.10e -- %.10e\n", its, rTr, err );
	    }
	}
	nlMultMatrixVector(M,d,aux);
	nlMultMatrixVector(P,aux,Sd);
        rTSd=blas->Ddot(blas,N,rT,1,Sd,1);
        nl_assert( fabs(rTSd)>1e-40 );
        alpha=rTh/rTSd;
        blas->Daxpy(blas,N,-alpha,aux,1,r,1);
        blas->Dcopy(blas,N,h,1,s,1);
        blas->Daxpy(blas,N,-alpha,Sd,1,s,1);
	nlMultMatrixVector(M,s,aux);
	nlMultMatrixVector(P,aux,t);
        blas->Daxpy(blas,N,1.,t,1,u,1);
        blas->Dscal(blas,N,alpha,u,1);
        st=blas->Ddot(blas,N,s,1,t,1);
        tt=blas->Ddot(blas,N,t,1,t,1);
        if ( fabs(st)<1e-40 || fabs(tt)<1e-40 ) {
            omega = 0.;
        } else {
            omega = st/tt;
        }
        blas->Daxpy(blas,N,-omega,aux,1,r,1);
        blas->Daxpy(blas,N,-alpha,d,1,x,1);
        blas->Daxpy(blas,N,-omega,s,1,x,1);
        blas->Dcopy(blas,N,s,1,h,1);
        blas->Daxpy(blas,N,-omega,t,1,h,1);
        beta=(alpha/omega)/rTh;
	rTh=blas->Ddot(blas,N,rT,1,h,1);
	beta*=rTh;
        blas->Dscal(blas,N,beta,d,1);
        blas->Daxpy(blas,N,1.,h,1,d,1);
        blas->Daxpy(blas,N,-beta*omega,Sd,1,d,1);
        rTr=blas->Ddot(blas,N,r,1,r,1);
        ++its;
    }
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, N, r);
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, N, rT);
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, N, d);
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, N, h);
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, N, u);
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, N, Sd);
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, N, t);
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, N, aux);
    blas->sq_bnorm = b_square;
    blas->sq_rnorm = rTr;
    return its;
}

/* 
 * Note: this one cannot be executed on device (GPU)
 * because it directly manipulates the vectors.
 */
static NLuint nlSolveSystem_GMRES(
    NLBlas_t blas,
    NLMatrix M, NLdouble* b, NLdouble* x,
    double eps, NLuint max_iter, NLuint inner_iter
) {
    NLint    n    = (NLint)M->n;
    NLint    m    = (NLint)inner_iter;
    typedef NLdouble *NLdoubleP;     
    NLdouble *V   = NL_NEW_ARRAY(NLdouble, n*(m+1)   );
    NLdouble *U   = NL_NEW_ARRAY(NLdouble, m*(m+1)/2 );
    NLdouble *r   = NL_NEW_ARRAY(NLdouble, n         );
    NLdouble *y   = NL_NEW_ARRAY(NLdouble, m+1       );
    NLdouble *c   = NL_NEW_ARRAY(NLdouble, m         );
    NLdouble *s   = NL_NEW_ARRAY(NLdouble, m         );
    NLdouble **v  = NL_NEW_ARRAY(NLdoubleP, m+1      );
    NLint i, j, io, uij, u0j; 
    NLint its = -1;
    NLdouble beta, h, rd, dd, nrm2b;

    /* 
     * The way it is written, this routine will not
     * work on the GPU since it directly modifies the
     * vectors.
     */
    nl_assert(nlBlasHasUnifiedMemory(blas));
    
    for ( i=0; i<=m; ++i ){
        v[i]=V+i*n;
    }
    
    nrm2b=blas->Dnrm2(blas,n,b,1);
    io=0;

    do  { /* outer loop */
        ++io;
	nlMultMatrixVector(M,x,r);
        blas->Daxpy(blas,n,-1.,b,1,r,1);
        beta=blas->Dnrm2(blas,n,r,1);
        blas->Dcopy(blas,n,r,1,v[0],1);
        blas->Dscal(blas,n,1./beta,v[0],1);

        y[0]=beta;
        j=0;
        uij=0;
        do { /* inner loop: j=0,...,m-1 */
            u0j=uij;
	    nlMultMatrixVector(M,v[j],v[j+1]);
            blas->Dgemv(
                blas,Transpose,n,j+1,1.,V,n,v[j+1],1,0.,U+u0j,1
            );
            blas->Dgemv(
                blas,NoTranspose,n,j+1,-1.,V,n,U+u0j,1,1.,v[j+1],1
            );
            h=blas->Dnrm2(blas,n,v[j+1],1);
            blas->Dscal(blas,n,1./h,v[j+1],1);
            for (i=0; i<j; ++i ) { /* rotiere neue Spalte */
                double tmp = c[i]*U[uij]-s[i]*U[uij+1];
                U[uij+1]   = s[i]*U[uij]+c[i]*U[uij+1];
                U[uij]     = tmp;
                ++uij;
            }
            { /* berechne neue Rotation */
                rd     = U[uij];
                dd     = sqrt(rd*rd+h*h);
                c[j]   = rd/dd;
                s[j]   = -h/dd;
                U[uij] = dd;
                ++uij;
            }
            { /* rotiere rechte Seite y (vorher: y[j+1]=0) */
                y[j+1] = s[j]*y[j];
                y[j]   = c[j]*y[j];
            }
            ++j;
        } while ( 
            j<m && fabs(y[j])>=eps*nrm2b 
        );
        { /* minimiere bzgl Y */
            blas->Dtpsv(
		blas,
                UpperTriangle,
                NoTranspose,
                NotUnitTriangular,
                j,U,y,1
            );
            /* correct X */
            blas->Dgemv(blas,NoTranspose,n,j,-1.,V,n,y,1,1.,x,1);
        }
    } while ( fabs(y[j])>=eps*nrm2b && (m*(io-1)+j) < (NLint)max_iter);
    
    /* Count the inner iterations */
    its = m*(io-1)+j;
    blas->sq_bnorm = nrm2b*nrm2b;
    blas->sq_rnorm = y[j]*y[j];
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, n, V);
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, n, U);
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, n, r);
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, n, y);
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, n, c);
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, n, s);
    NL_DELETE_VECTOR(blas, NL_DEVICE_MEMORY, n, v);
    return (NLuint)its;
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
	case NL_BICGSTAB:
	    if(P == NULL) {
		result = nlSolveSystem_BICGSTAB(blas,M,b,x,eps,max_iter);
	    } else {
		result = nlSolveSystem_PRE_BICGSTAB(blas,M,P,b,x,eps,max_iter);
	    }
	    break;
	case NL_GMRES:
	    result = nlSolveSystem_GMRES(blas,M,b,x,eps,max_iter,inner_iter);
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
	    if(nlCurrentContext->verbose) {
		nl_printf("in OpenNL : ||Ax-b|| = %e\n",nlCurrentContext->error);
	    }
	} else {
	    nlCurrentContext->error = rnorm/bnorm;
	    if(nlCurrentContext->verbose) {
		nl_printf("in OpenNL : ||Ax-b||/||b|| = %e\n",
		       nlCurrentContext->error
		);
	    }
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




typedef struct {
    NLuint m;

    NLuint n;

    NLenum type;

    NLDestroyMatrixFunc destroy_func;

    NLMultMatrixVectorFunc mult_func;

    NLSparseMatrix* M;

    double omega;
    
    NLdouble* work;
    
} NLSSORPreconditioner;


static void nlSSORPreconditionerDestroy(NLSSORPreconditioner* M) {
    NL_DELETE_ARRAY(M->work);
}



static void nlSparseMatrixMultLowerInverse(
    NLSparseMatrix* A, const NLdouble* x, NLdouble* y, double omega
) {
    NLuint n       = A->n;
    NLdouble* diag = A->diag;
    NLuint i;
    NLuint ij;
    NLCoeff* c = NULL;
    NLdouble S;

    nl_assert(A->storage & NL_MATRIX_STORE_SYMMETRIC);
    nl_assert(A->storage & NL_MATRIX_STORE_ROWS);

    for(i=0; i<n; i++) {
        NLRowColumn*  Ri = &(A->row[i]);       
        S = 0;
        for(ij=0; ij < Ri->size; ij++) {
            c = &(Ri->coeff[ij]);
            nl_parano_assert(c->index <= i); 
            if(c->index != i) {
                S += c->value * y[c->index]; 
            }
        }
        nlHostBlas()->flops += (NLulong)(2*Ri->size);                    
        y[i] = (x[i] - S) * omega / diag[i];
    }
    nlHostBlas()->flops += (NLulong)(n*3);                
}
static void nlSparseMatrixMultUpperInverse(
    NLSparseMatrix* A, const NLdouble* x, NLdouble* y, NLdouble omega
) {
    NLuint n       = A->n;
    NLdouble* diag = A->diag;
    NLint i;
    NLuint ij;
    NLCoeff* c = NULL;
    NLdouble S;

    nl_assert(A->storage & NL_MATRIX_STORE_SYMMETRIC);
    nl_assert(A->storage & NL_MATRIX_STORE_COLUMNS);

    for(i=(NLint)(n-1); i>=0; i--) {
        NLRowColumn*  Ci = &(A->column[i]);       
        S = 0;
        for(ij=0; ij < Ci->size; ij++) {
            c = &(Ci->coeff[ij]);
            nl_parano_assert(c->index >= i); 
            if((NLint)(c->index) != i) {
                S += c->value * y[c->index]; 
            }
        }
        nlHostBlas()->flops += (NLulong)(2*Ci->size);                    
        y[i] = (x[i] - S) * omega / diag[i];
    }
    nlHostBlas()->flops += (NLulong)(n*3);                
}


static void nlSSORPreconditionerMult(
    NLSSORPreconditioner* P, const double* x, double* y
) {
    NLdouble* diag = P->M->diag;
    NLuint i;
    nlSparseMatrixMultLowerInverse(
        P->M, x, P->work, P->omega
    );
    for(i=0; i<P->n; i++) {
        P->work[i] *= (diag[i] / P->omega);
    }
    nlHostBlas()->flops += (NLulong)(P->n);
    nlSparseMatrixMultUpperInverse(
        P->M, P->work, y, P->omega
    );
    nlHostBlas()->Dscal(nlHostBlas(),(NLint)P->n, 2.0 - P->omega, y, 1);
}

NLMatrix nlNewSSORPreconditioner(NLMatrix M_in, double omega) {
    NLSparseMatrix* M = NULL;
    NLSSORPreconditioner* result = NULL;
    nl_assert(M_in->type == NL_MATRIX_SPARSE_DYNAMIC);
    nl_assert(M_in->m == M_in->n);
    M = (NLSparseMatrix*)M_in;
    result = NL_NEW(NLSSORPreconditioner);
    result->m = M->m;
    result->n = M->n;
    result->type = NL_MATRIX_OTHER;
    result->destroy_func = (NLDestroyMatrixFunc)nlSSORPreconditionerDestroy;
    result->mult_func = (NLMultMatrixVectorFunc)nlSSORPreconditionerMult;
    result->M = M;
    result->work = NL_NEW_ARRAY(NLdouble, result->n);
    result->omega = omega;
    return (NLMatrix)result;
}



/******* extracted from nl_superlu.c *******/



#ifdef NL_OS_UNIX
#  ifdef NL_OS_APPLE
#      define SUPERLU_LIB_NAME "libsuperlu_5.dylib"
#  else
#      define SUPERLU_LIB_NAME "libsuperlu.so"
#  endif
#else
#  define SUPERLU_LIB_NAME "libsuperlu.xxx"
#endif









typedef enum {
    SLU_NC,    /* column-wise, no supernode */
    SLU_NCP,   /* column-wise, column-permuted, no supernode 
                  (The consecutive columns of nonzeros, after permutation,
                   may not be stored  contiguously.) */
    SLU_NR,    /* row-wize, no supernode */
    SLU_SC,    /* column-wise, supernode */
    SLU_SCP,   /* supernode, column-wise, permuted */    
    SLU_SR,    /* row-wise, supernode */
    SLU_DN,     /* Fortran style column-wise storage for dense matrix */
    SLU_NR_loc  /* distributed compressed row format  */ 
} Stype_t;

typedef enum {
    SLU_S,     /* single */
    SLU_D,     /* double */
    SLU_C,     /* single complex */
    SLU_Z      /* double complex */
} Dtype_t;


typedef enum {
    SLU_GE,    /* general */
    SLU_TRLU,  /* lower triangular, unit diagonal */
    SLU_TRUU,  /* upper triangular, unit diagonal */
    SLU_TRL,   /* lower triangular */
    SLU_TRU,   /* upper triangular */
    SLU_SYL,   /* symmetric, store lower half */
    SLU_SYU,   /* symmetric, store upper half */
    SLU_HEL,   /* Hermitian, store lower half */
    SLU_HEU    /* Hermitian, store upper half */
} Mtype_t;

typedef int int_t;

typedef struct {
    int_t  nnz;	    /* number of nonzeros in the matrix */
    void *nzval;    /* pointer to array of nonzero values, packed by raw */
    int_t  *colind; /* pointer to array of columns indices of the nonzeros */
    int_t  *rowptr; /* pointer to array of beginning of rows in nzval[] 
		       and colind[]  */
                    /* Note:
		       Zero-based indexing is used;
		       rowptr[] has nrow+1 entries, the last one pointing
		       beyond the last row, so that rowptr[nrow] = nnz. */
} NRformat;

typedef struct {
        Stype_t Stype; /* Storage type: interprets the storage structure 
                          pointed to by *Store. */
        Dtype_t Dtype; /* Data type. */
        Mtype_t Mtype; /* Matrix type: describes the mathematical property of 
                          the matrix. */
        int_t  nrow;   /* number of rows */
        int_t  ncol;   /* number of columns */
        void *Store;   /* pointer to the actual storage of the matrix */
} SuperMatrix;

/* Stype == SLU_DN */
typedef struct {
    int_t lda;    /* leading dimension */
    void *nzval;  /* array of size lda*ncol to represent a dense matrix */
} DNformat;


typedef enum {NO, YES}                                          yes_no_t;
typedef enum {DOFACT, SamePattern, SamePattern_SameRowPerm, FACTORED} fact_t;
typedef enum {NOROWPERM, LargeDiag, MY_PERMR}                   rowperm_t;
typedef enum {NATURAL, MMD_ATA, MMD_AT_PLUS_A, COLAMD,
              METIS_AT_PLUS_A, PARMETIS, ZOLTAN, MY_PERMC}      colperm_t;
typedef enum {NOTRANS, TRANS, CONJ}                             trans_t;
typedef enum {NOEQUIL, ROW, COL, BOTH}                          DiagScale_t;
typedef enum {NOREFINE, SLU_SINGLE=1, SLU_DOUBLE, SLU_EXTRA}    IterRefine_t;
typedef enum {LUSUP, UCOL, LSUB, USUB, LLVL, ULVL}              MemType;
typedef enum {HEAD, TAIL}                                       stack_end_t;
typedef enum {SYSTEM, USER}                                     LU_space_t;
typedef enum {ONE_NORM, TWO_NORM, INF_NORM}                     norm_t;
typedef enum {SILU, SMILU_1, SMILU_2, SMILU_3}                  milu_t;

typedef struct {
    fact_t        Fact;
    yes_no_t      Equil;
    colperm_t     ColPerm;
    trans_t       Trans;
    IterRefine_t  IterRefine;
    double        DiagPivotThresh;
    yes_no_t      SymmetricMode;
    yes_no_t      PivotGrowth;
    yes_no_t      ConditionNumber;
    rowperm_t     RowPerm;
    int           ILU_DropRule;
    double        ILU_DropTol;    /* threshold for dropping */
    double        ILU_FillFactor; /* gamma in the secondary dropping */
    norm_t        ILU_Norm;       /* infinity-norm, 1-norm, or 2-norm */
    double        ILU_FillTol;    /* threshold for zero pivot perturbation */
    milu_t        ILU_MILU;
    double        ILU_MILU_Dim;   /* Dimension of PDE (if available) */
    yes_no_t      ParSymbFact;
    yes_no_t      ReplaceTinyPivot; /* used in SuperLU_DIST */
    yes_no_t      SolveInitialized;
    yes_no_t      RefineInitialized;
    yes_no_t      PrintStat;
    int           nnzL, nnzU;      /* used to store nnzs for now       */
    int           num_lookaheads;  /* num of levels in look-ahead      */
    yes_no_t      lookahead_etree; /* use etree computed from the
                                      serial symbolic factorization */
    yes_no_t      SymPattern;      /* symmetric factorization          */
} superlu_options_t;

typedef void* superlu_options_ptr;

typedef float    flops_t;
typedef unsigned char Logical;

typedef struct {
    int     *panel_histo;    /* histogram of panel size distribution */
    double  *utime;          /* running time at various phases */
    flops_t *ops;            /* operation count at various phases */
    int     TinyPivots;      /* number of tiny pivots */
    int     RefineSteps;     /* number of iterative refinement steps */
    int     expansions;      /* number of memory expansions (SuperLU4) */
} SuperLUStat_t;

/*! \brief Headers for 4 types of dynamatically managed memory */
typedef struct e_node {
    int size;      /* length of the memory that has been used */
    void *mem;     /* pointer to the new malloc'd store */
} ExpHeader;

typedef struct {
    int  size;
    int  used;
    int  top1;  /* grow upward, relative to &array[0] */
    int  top2;  /* grow downward */
    void *array;
} LU_stack_t;

typedef struct {
    int     *xsup;    /* supernode and column mapping */
    int     *supno;   
    int     *lsub;    /* compressed L subscripts */
    int	    *xlsub;
    void    *lusup;   /* L supernodes */
    int     *xlusup;
    void    *ucol;    /* U columns */
    int     *usub;
    int	    *xusub;
    int     nzlmax;   /* current max size of lsub */
    int     nzumax;   /*    "    "    "      ucol */
    int     nzlumax;  /*    "    "    "     lusup */
    int     n;        /* number of columns in the matrix */
    LU_space_t MemModel; /* 0 - system malloc'd; 1 - user provided */
    int     num_expansions;
    ExpHeader *expanders; /* Array of pointers to 4 types of memory */
    LU_stack_t stack;     /* use user supplied memory */
} GlobalLU_t;













typedef void (*FUNPTR_set_default_options)(superlu_options_ptr options);
typedef void (*FUNPTR_ilu_set_default_options)(superlu_options_ptr options);
typedef void (*FUNPTR_StatInit)(SuperLUStat_t *);
typedef void (*FUNPTR_StatFree)(SuperLUStat_t *);

typedef void (*FUNPTR_dCreate_CompCol_Matrix)(
    SuperMatrix *, int, int, int, const double *,
    const int *, const int *, Stype_t, Dtype_t, Mtype_t);

typedef void (*FUNPTR_dCreate_Dense_Matrix)(
    SuperMatrix *, int, int, const double *, int,
    Stype_t, Dtype_t, Mtype_t);

typedef void (*FUNPTR_Destroy_SuperNode_Matrix)(SuperMatrix *);
typedef void (*FUNPTR_Destroy_CompCol_Matrix)(SuperMatrix *);
typedef void (*FUNPTR_Destroy_CompCol_Permuted)(SuperMatrix *);
typedef void (*FUNPTR_Destroy_SuperMatrix_Store)(SuperMatrix *);

typedef void (*FUNPTR_dgssv)(
    superlu_options_ptr, SuperMatrix *, int *, int *, SuperMatrix *,
    SuperMatrix *, SuperMatrix *, SuperLUStat_t *, int *
);

typedef void (*FUNPTR_dgstrs)(
    trans_t, SuperMatrix *, SuperMatrix *, int *, int *,
    SuperMatrix *, SuperLUStat_t*, int *    
);

typedef void (*FUNPTR_get_perm_c)(int, SuperMatrix *, int *);
typedef void (*FUNPTR_sp_preorder)(
   superlu_options_t *, SuperMatrix*, int*, int*, SuperMatrix*
);
typedef int (*FUNPTR_sp_ienv)(int);
typedef int (*FUNPTR_input_error)(const char *, int *);

typedef void (*FUNPTR_dgstrf) (superlu_options_t *options, SuperMatrix *A,
        int relax, int panel_size, int *etree, void *work, int lwork,
        int *perm_c, int *perm_r, SuperMatrix *L, SuperMatrix *U,
    	GlobalLU_t *Glu, /* persistent to facilitate multiple factorizations */
        SuperLUStat_t *stat, int *info
);


typedef struct {
    FUNPTR_set_default_options set_default_options;
    FUNPTR_ilu_set_default_options ilu_set_default_options;    
    FUNPTR_StatInit StatInit;
    FUNPTR_StatFree StatFree;
    FUNPTR_dCreate_CompCol_Matrix dCreate_CompCol_Matrix;
    FUNPTR_dCreate_Dense_Matrix dCreate_Dense_Matrix;
    FUNPTR_Destroy_SuperNode_Matrix Destroy_SuperNode_Matrix;
    FUNPTR_Destroy_CompCol_Matrix Destroy_CompCol_Matrix;
    FUNPTR_Destroy_CompCol_Permuted Destroy_CompCol_Permuted;    
    FUNPTR_Destroy_SuperMatrix_Store Destroy_SuperMatrix_Store;
    FUNPTR_dgssv dgssv;
    FUNPTR_dgstrs dgstrs;
    FUNPTR_get_perm_c get_perm_c;
    FUNPTR_sp_preorder sp_preorder;
    FUNPTR_sp_ienv sp_ienv;
    FUNPTR_dgstrf dgstrf;
    FUNPTR_input_error input_error;
    
    NLdll DLL_handle;
} SuperLUContext;

static SuperLUContext* SuperLU() {
    static SuperLUContext context;
    static NLboolean init = NL_FALSE;
    if(!init) {
        init = NL_TRUE;
        memset(&context, 0, sizeof(context));
    }
    return &context;
}

NLboolean nlExtensionIsInitialized_SUPERLU() {
    return
        SuperLU()->DLL_handle != NULL &&
        SuperLU()->set_default_options != NULL &&
        SuperLU()->ilu_set_default_options != NULL &&   
        SuperLU()->StatInit != NULL &&
        SuperLU()->StatFree != NULL &&
        SuperLU()->dCreate_CompCol_Matrix != NULL &&
        SuperLU()->dCreate_Dense_Matrix != NULL &&
        SuperLU()->Destroy_SuperNode_Matrix != NULL &&
        SuperLU()->Destroy_CompCol_Matrix != NULL &&
        SuperLU()->Destroy_CompCol_Permuted != NULL &&	
        SuperLU()->Destroy_SuperMatrix_Store != NULL &&
        SuperLU()->dgssv != NULL &&
        SuperLU()->dgstrs != NULL &&
	SuperLU()->get_perm_c != NULL &&
	SuperLU()->sp_preorder != NULL &&
	SuperLU()->sp_ienv != NULL &&
	SuperLU()->dgstrf != NULL &&
	SuperLU()->input_error != NULL;
}

static void nlTerminateExtension_SUPERLU(void) {
    if(SuperLU()->DLL_handle != NULL) {
        nlCloseDLL(SuperLU()->DLL_handle);
        SuperLU()->DLL_handle = NULL;
    }
}


#define find_superlu_func(name)                                   \
    if(                                                           \
        (                                                         \
            SuperLU()->name =                                     \
            (FUNPTR_##name)nlFindFunction(SuperLU()->DLL_handle,#name) \
        ) == NULL                                                 \
    ) {                                                           \
        nlError("nlInitExtension_SUPERLU","function not found");  \
        nlError("nlInitExtension_SUPERLU",#name);                 \
        return NL_FALSE;                                          \
    }


NLboolean nlInitExtension_SUPERLU(void) {
    NLenum flags = NL_LINK_NOW | NL_LINK_USE_FALLBACK;
    if(nlCurrentContext == NULL || !nlCurrentContext->verbose) {
	flags |= NL_LINK_QUIET;
    }
    
    if(SuperLU()->DLL_handle != NULL) {
        return nlExtensionIsInitialized_SUPERLU();
    }

    SuperLU()->DLL_handle = nlOpenDLL(SUPERLU_LIB_NAME, flags);
    if(SuperLU()->DLL_handle == NULL) {
        return NL_FALSE;
    }
    
    find_superlu_func(set_default_options);
    find_superlu_func(ilu_set_default_options);    
    find_superlu_func(StatInit);
    find_superlu_func(StatFree);
    find_superlu_func(dCreate_CompCol_Matrix);
    find_superlu_func(dCreate_Dense_Matrix);
    find_superlu_func(Destroy_SuperNode_Matrix);
    find_superlu_func(Destroy_CompCol_Matrix);
    find_superlu_func(Destroy_CompCol_Permuted);        
    find_superlu_func(Destroy_SuperMatrix_Store);
    find_superlu_func(dgssv);
    find_superlu_func(dgstrs);
    find_superlu_func(get_perm_c);    
    find_superlu_func(sp_preorder);
    find_superlu_func(sp_ienv);
    find_superlu_func(dgstrf);
    find_superlu_func(input_error);

    atexit(nlTerminateExtension_SUPERLU);
    return NL_TRUE;
}



typedef struct {
    NLuint m;

    NLuint n;

    NLenum type;

    NLDestroyMatrixFunc destroy_func;

    NLMultMatrixVectorFunc mult_func;

    SuperMatrix L;

    SuperMatrix U;

    int* perm_r;

    int* perm_c;

    trans_t trans;
    
} NLSuperLUFactorizedMatrix;


static void nlSuperLUFactorizedMatrixDestroy(NLSuperLUFactorizedMatrix* M) {
    SuperLU()->Destroy_SuperNode_Matrix(&M->L);
    SuperLU()->Destroy_CompCol_Matrix(&M->U);    
    NL_DELETE_ARRAY(M->perm_r);
    NL_DELETE_ARRAY(M->perm_c);
}

static void nlSuperLUFactorizedMatrixMult(
    NLSuperLUFactorizedMatrix* M, const double* x, double* y
) {
    SuperMatrix B;
    SuperLUStat_t stat;
    int info = 0;
    NLuint i;

    /* Create vector */
    SuperLU()->dCreate_Dense_Matrix(
        &B, (int)(M->n), 1, y, (int)(M->n), 
        SLU_DN, /* Fortran-type column-wise storage */
        SLU_D,  /* doubles */
        SLU_GE  /* general */
    );

    /* copy rhs onto y (superLU matrix-vector product expects it here */
    for(i = 0; i < M->n; i++){
        y[i] = x[i];
    }

    /* Call SuperLU triangular solve */
    SuperLU()->StatInit(&stat) ;

    SuperLU()->dgstrs(
       M->trans, &M->L, &M->U, M->perm_c, M->perm_r, &B, &stat, &info
    );

    SuperLU()->StatFree(&stat) ;
    
    /*  Only the "store" structure needs to be 
     *  deallocated (the array has been allocated
     * by client code).
     */
    SuperLU()->Destroy_SuperMatrix_Store(&B) ;
}

/*
 * Copied from SUPERLU/dgssv.c, removed call to linear solve.
 */
static void dgssv_factorize_only(
      superlu_options_t *options, SuperMatrix *A, int *perm_c, int *perm_r,
      SuperMatrix *L, SuperMatrix *U,
      SuperLUStat_t *stat, int *info, trans_t *trans
) {
    SuperMatrix *AA = NULL;
        /* A in SLU_NC format used by the factorization routine.*/
    SuperMatrix AC; /* Matrix postmultiplied by Pc */
    int      lwork = 0, *etree, i;
    GlobalLU_t Glu; /* Not needed on return. */
    
    /* Set default values for some parameters */
    int      panel_size;     /* panel size */
    int      relax;          /* no of columns in a relaxed snodes */
    int      permc_spec;

    nl_assert(A->Stype == SLU_NR || A->Stype == SLU_NC);
    
    *trans = NOTRANS;

    if ( options->Fact != DOFACT ) *info = -1;
    else if ( A->nrow != A->ncol || A->nrow < 0 ||
	 (A->Stype != SLU_NC && A->Stype != SLU_NR) ||
	 A->Dtype != SLU_D || A->Mtype != SLU_GE )
	*info = -2;
    if ( *info != 0 ) {
	i = -(*info);
	SuperLU()->input_error("SUPERLU/OpenNL dgssv_factorize_only", &i);
	return;
    }

    /* Convert A to SLU_NC format when necessary. */
    if ( A->Stype == SLU_NR ) {
	NRformat *Astore = (NRformat*)A->Store;
	AA = NL_NEW(SuperMatrix);
	SuperLU()->dCreate_CompCol_Matrix(
	    AA, A->ncol, A->nrow, Astore->nnz, 
	    (double*)Astore->nzval, Astore->colind, Astore->rowptr,
	    SLU_NC, A->Dtype, A->Mtype
	);
	*trans = TRANS;
    } else {
        if ( A->Stype == SLU_NC ) AA = A;
    }

    nl_assert(AA != NULL);
    
    /*
     * Get column permutation vector perm_c[], according to permc_spec:
     *   permc_spec = NATURAL:  natural ordering 
     *   permc_spec = MMD_AT_PLUS_A: minimum degree on structure of A'+A
     *   permc_spec = MMD_ATA:  minimum degree on structure of A'*A
     *   permc_spec = COLAMD:   approximate minimum degree column ordering
     *   permc_spec = MY_PERMC: the ordering already supplied in perm_c[]
     */
    permc_spec = options->ColPerm;
    if ( permc_spec != MY_PERMC && options->Fact == DOFACT )
	SuperLU()->get_perm_c(permc_spec, AA, perm_c);
    
    etree = NL_NEW_ARRAY(int,A->ncol);
    SuperLU()->sp_preorder(options, AA, perm_c, etree, &AC);
    panel_size = SuperLU()->sp_ienv(1);
    relax = SuperLU()->sp_ienv(2);
    SuperLU()->dgstrf(options, &AC, relax, panel_size, etree,
            NULL, lwork, perm_c, perm_r, L, U, &Glu, stat, info);

    NL_DELETE_ARRAY(etree);
    SuperLU()->Destroy_CompCol_Permuted(&AC);
    if ( A->Stype == SLU_NR ) {
	SuperLU()->Destroy_SuperMatrix_Store(AA);
	NL_DELETE(AA);
    }
}


NLMatrix nlMatrixFactorize_SUPERLU(
    NLMatrix M, NLenum solver
) {
    NLSuperLUFactorizedMatrix* LU = NULL;
    NLCRSMatrix* CRS = NULL;
    SuperMatrix superM;
    NLuint n = M->n;
    superlu_options_t options;
    SuperLUStat_t     stat;
    NLint info = 0;       /* status code  */
    
    nl_assert(M->m == M->n);

    if(M->type == NL_MATRIX_CRS) {
        CRS = (NLCRSMatrix*)M;
    } else if(M->type == NL_MATRIX_SPARSE_DYNAMIC) {
        CRS = (NLCRSMatrix*)nlCRSMatrixNewFromSparseMatrix((NLSparseMatrix*)M);
    }

    nl_assert(!(CRS->symmetric_storage));
    
    LU = NL_NEW(NLSuperLUFactorizedMatrix);
    LU->m = M->m;
    LU->n = M->n;
    LU->type = NL_MATRIX_OTHER;
    LU->destroy_func = (NLDestroyMatrixFunc)(nlSuperLUFactorizedMatrixDestroy);
    LU->mult_func = (NLMultMatrixVectorFunc)(nlSuperLUFactorizedMatrixMult);
    LU->perm_c = NL_NEW_ARRAY(int, n);
    LU->perm_r = NL_NEW_ARRAY(int, n);    

    SuperLU()->dCreate_CompCol_Matrix(
        &superM, (int)n, (int)n, (int)nlCRSMatrixNNZ(CRS),
        CRS->val, (int*)CRS->colind, (int*)CRS->rowptr, 
        SLU_NR,              /* Row_wise, no supernode */
        SLU_D,               /* doubles                */ 
        CRS->symmetric_storage ? SLU_SYL : SLU_GE
    );

    SuperLU()->set_default_options(&options);
    switch(solver) {
    case NL_SUPERLU_EXT: {
        options.ColPerm = NATURAL;
    } break;
    case NL_PERM_SUPERLU_EXT: {
        options.ColPerm = COLAMD;
    } break;
    case NL_SYMMETRIC_SUPERLU_EXT: {
        options.ColPerm = MMD_AT_PLUS_A;
        options.SymmetricMode = YES;
    } break;
    default: 
        nl_assert_not_reached;
    }
    
    SuperLU()->StatInit(&stat);

    dgssv_factorize_only(
	  &options, &superM, LU->perm_c, LU->perm_r,
	  &LU->L, &LU->U, &stat, &info, &LU->trans
    );

    SuperLU()->StatFree(&stat);
    
    /*
     * Only the "store" structure needs to be deallocated 
     * (the arrays have been allocated by us, they are in CRS).
     */
    SuperLU()->Destroy_SuperMatrix_Store(&superM);
    
    if((NLMatrix)CRS != M) {
        nlDeleteMatrix((NLMatrix)CRS);
    }

    if(info != 0) {
	NL_DELETE(LU);
	LU = NULL;
    }
    return (NLMatrix)LU;
}

/******* extracted from nl_cholmod.c *******/



#ifdef NL_OS_UNIX
#  ifdef NL_OS_APPLE
#      define CHOLMOD_LIB_NAME "libcholmod.dylib"
#  else
#      define CHOLMOD_LIB_NAME "libcholmod.so"
#  endif
#else
#  define CHOLMOD_LIB_NAME "libcholmod.xxx"
#endif


/*            Excerpt from cholmod_core.h                         */


/* A dense matrix in column-oriented form.  It has no itype since it contains
 * no integers.  Entry in row i and column j is located in x [i+j*d].
 */
typedef struct cholmod_dense_struct {
    size_t nrow ;       /* the matrix is nrow-by-ncol */
    size_t ncol ;
    size_t nzmax ;      /* maximum number of entries in the matrix */
    size_t d ;          /* leading dimension (d >= nrow must hold) */
    void *x ;           /* size nzmax or 2*nzmax, if present */
    void *z ;           /* size nzmax, if present */
    int xtype ;         /* pattern, real, complex, or zomplex */
    int dtype ;         /* x and z double or float */
} cholmod_dense ;

/* A sparse matrix stored in compressed-column form. */

typedef struct cholmod_sparse_struct
{
    size_t nrow ;       /* the matrix is nrow-by-ncol */
    size_t ncol ;
    size_t nzmax ;      /* maximum number of entries in the matrix */

    /* pointers to int or SuiteSparse_long: */
    void *p ;           /* p [0..ncol], the column pointers */
    void *i ;           /* i [0..nzmax-1], the row indices */

    /* for unpacked matrices only: */
    void *nz ;          /* nz [0..ncol-1], the # of nonzeros in each col.  In
                         * packed form, the nonzero pattern of column j is in
        * A->i [A->p [j] ... A->p [j+1]-1].  In unpacked form, column j is in
        * A->i [A->p [j] ... A->p [j]+A->nz[j]-1] instead.  In both cases, the
        * numerical values (if present) are in the corresponding locations in
        * the array x (or z if A->xtype is CHOLMOD_ZOMPLEX). */

    /* pointers to double or float: */
    void *x ;           /* size nzmax or 2*nzmax, if present */
    void *z ;           /* size nzmax, if present */

    int stype ;         /* Describes what parts of the matrix are considered:
                         *
        * 0:  matrix is "unsymmetric": use both upper and lower triangular parts
        *     (the matrix may actually be symmetric in pattern and value, but
        *     both parts are explicitly stored and used).  May be square or
        *     rectangular.
        * >0: matrix is square and symmetric, use upper triangular part.
        *     Entries in the lower triangular part are ignored.
        * <0: matrix is square and symmetric, use lower triangular part.
        *     Entries in the upper triangular part are ignored.
        *
        * Note that stype>0 and stype<0 are different for cholmod_sparse and
        * cholmod_triplet.  See the cholmod_triplet data structure for more
        * details.
        */

    int itype ;         /* CHOLMOD_INT:     p, i, and nz are int.
                         * CHOLMOD_INTLONG: p is SuiteSparse_long,
                         *                  i and nz are int.
                         * CHOLMOD_LONG:    p, i, and nz are SuiteSparse_long */

    int xtype ;         /* pattern, real, complex, or zomplex */
    int dtype ;         /* x and z are double or float */
    int sorted ;        /* TRUE if columns are sorted, FALSE otherwise */
    int packed ;        /* TRUE if packed (nz ignored), FALSE if unpacked
                         * (nz is required) */

} cholmod_sparse ;



typedef void* cholmod_common_ptr;
typedef cholmod_dense* cholmod_dense_ptr;
typedef cholmod_sparse* cholmod_sparse_ptr;
typedef void* cholmod_factor_ptr;


typedef enum cholmod_xtype_enum {
    CHOLMOD_PATTERN =0,
    CHOLMOD_REAL    =1,
    CHOLMOD_COMPLEX =2,
    CHOLMOD_ZOMPLEX =3
} cholmod_xtype;


typedef enum cholmod_solve_type_enum {
   CHOLMOD_A    =0,   
   CHOLMOD_LDLt =1,   
   CHOLMOD_LD   =2,   
   CHOLMOD_DLt  =3,   
   CHOLMOD_L    =4,   
   CHOLMOD_Lt   =5,   
   CHOLMOD_D    =6,   
   CHOLMOD_P    =7,   
   CHOLMOD_Pt   =8   
} cholmod_solve_type;
    
typedef int cholmod_stype;

typedef void (*FUNPTR_cholmod_start)(cholmod_common_ptr);

typedef cholmod_sparse_ptr (*FUNPTR_cholmod_allocate_sparse)(
    size_t m, size_t n, size_t nnz, int sorted,
    int packed, int stype, int xtype, cholmod_common_ptr
);

typedef cholmod_dense_ptr (*FUNPTR_cholmod_allocate_dense)(
    size_t m, size_t n, size_t d, int xtype, cholmod_common_ptr
);

typedef cholmod_factor_ptr (*FUNPTR_cholmod_analyze)(
    cholmod_sparse_ptr A, cholmod_common_ptr
);

typedef int (*FUNPTR_cholmod_factorize)(
    cholmod_sparse_ptr A, cholmod_factor_ptr L, cholmod_common_ptr
);

typedef cholmod_dense_ptr (*FUNPTR_cholmod_solve)(
    int solve_type, cholmod_factor_ptr, cholmod_dense_ptr, cholmod_common_ptr
);

typedef void (*FUNPTR_cholmod_free_factor)(
    cholmod_factor_ptr*, cholmod_common_ptr
);

typedef void (*FUNPTR_cholmod_free_dense)(
    cholmod_dense_ptr*, cholmod_common_ptr
);

typedef void (*FUNPTR_cholmod_free_sparse)(
    cholmod_sparse_ptr*, cholmod_common_ptr
);

typedef void (*FUNPTR_cholmod_finish)(cholmod_common_ptr);

typedef struct {
    char cholmod_common[16384];

    FUNPTR_cholmod_start cholmod_start;
    FUNPTR_cholmod_allocate_sparse cholmod_allocate_sparse;
    FUNPTR_cholmod_allocate_dense cholmod_allocate_dense;
    FUNPTR_cholmod_analyze cholmod_analyze;
    FUNPTR_cholmod_factorize cholmod_factorize;
    FUNPTR_cholmod_solve cholmod_solve;
    FUNPTR_cholmod_free_factor cholmod_free_factor;
    FUNPTR_cholmod_free_sparse cholmod_free_sparse;        
    FUNPTR_cholmod_free_dense cholmod_free_dense;
    FUNPTR_cholmod_finish cholmod_finish;
    
    NLdll DLL_handle;
} CHOLMODContext;

static CHOLMODContext* CHOLMOD() {
    static CHOLMODContext context;
    static NLboolean init = NL_FALSE;
    if(!init) {
        init = NL_TRUE;
        memset(&context, 0, sizeof(context));
    }
    return &context;
}

NLboolean nlExtensionIsInitialized_CHOLMOD() {
    return
        CHOLMOD()->DLL_handle != NULL &&
        CHOLMOD()->cholmod_start != NULL &&
        CHOLMOD()->cholmod_allocate_sparse != NULL &&
        CHOLMOD()->cholmod_allocate_dense != NULL &&
        CHOLMOD()->cholmod_analyze != NULL &&
        CHOLMOD()->cholmod_factorize != NULL &&
        CHOLMOD()->cholmod_solve != NULL &&
        CHOLMOD()->cholmod_free_factor != NULL &&
        CHOLMOD()->cholmod_free_sparse != NULL &&
        CHOLMOD()->cholmod_free_dense != NULL &&
        CHOLMOD()->cholmod_finish != NULL ;
}

#define find_cholmod_func(name)                                        \
    if(                                                                \
        (                                                              \
            CHOLMOD()->name =                                          \
            (FUNPTR_##name)nlFindFunction(CHOLMOD()->DLL_handle,#name) \
        ) == NULL                                                      \
    ) {                                                                \
        nlError("nlInitExtension_CHOLMOD","function not found");       \
        return NL_FALSE;                                               \
    }


static void nlTerminateExtension_CHOLMOD(void) {
    if(CHOLMOD()->DLL_handle != NULL) {
        CHOLMOD()->cholmod_finish(&CHOLMOD()->cholmod_common);
        nlCloseDLL(CHOLMOD()->DLL_handle);
        CHOLMOD()->DLL_handle = NULL;
    }
}

NLboolean nlInitExtension_CHOLMOD(void) {
    NLenum flags = NL_LINK_NOW | NL_LINK_USE_FALLBACK;
    if(nlCurrentContext == NULL || !nlCurrentContext->verbose) {
	flags |= NL_LINK_QUIET;
    }
    
    if(CHOLMOD()->DLL_handle != NULL) {
        return nlExtensionIsInitialized_CHOLMOD();
    }

    /*
     *   MKL has a built-in CHOLMOD that conflicts with
     * the CHOLMOD used by OpenNL (to be fixed). For now
     * we simply output a warning message and deactivate the
     * CHOLMOD extension if the MKL extension was initialized
     * before.
     */
    if(NLMultMatrixVector_MKL != NULL) {
	nl_fprintf(
	    stderr,
	    "CHOLMOD extension incompatible with MKL (deactivating)"
	);
	return NL_FALSE;
    }

    
    CHOLMOD()->DLL_handle = nlOpenDLL(CHOLMOD_LIB_NAME,flags);
    if(CHOLMOD()->DLL_handle == NULL) {
        return NL_FALSE;
    }

    find_cholmod_func(cholmod_start);
    find_cholmod_func(cholmod_allocate_sparse);
    find_cholmod_func(cholmod_allocate_dense);
    find_cholmod_func(cholmod_analyze);
    find_cholmod_func(cholmod_factorize);
    find_cholmod_func(cholmod_solve);
    find_cholmod_func(cholmod_free_factor);
    find_cholmod_func(cholmod_free_sparse);
    find_cholmod_func(cholmod_free_dense);
    find_cholmod_func(cholmod_finish);

    CHOLMOD()->cholmod_start(&CHOLMOD()->cholmod_common);

    atexit(nlTerminateExtension_CHOLMOD);
    return NL_TRUE;
}



typedef struct {
    NLuint m;

    NLuint n;

    NLenum type;

    NLDestroyMatrixFunc destroy_func;

    NLMultMatrixVectorFunc mult_func;

    cholmod_factor_ptr L;
    
} NLCholmodFactorizedMatrix;

static void nlCholmodFactorizedMatrixDestroy(NLCholmodFactorizedMatrix* M) {
    CHOLMOD()->cholmod_free_factor(&M->L, &CHOLMOD()->cholmod_common);
}

static void nlCholmodFactorizedMatrixMult(
    NLCholmodFactorizedMatrix* M, const double* x, double* y
) {
    /* 
     * TODO: see whether CHOLDMOD can use user-allocated vectors
     * (and avoid copy)
     */
    cholmod_dense_ptr X=CHOLMOD()->cholmod_allocate_dense(
        M->n, 1, M->n, CHOLMOD_REAL, &CHOLMOD()->cholmod_common
    );
    cholmod_dense_ptr Y=NULL;

    memcpy(X->x, x, M->n*sizeof(double));    
    Y = CHOLMOD()->cholmod_solve(
	CHOLMOD_A, M->L, X, &CHOLMOD()->cholmod_common
    );
    memcpy(y, Y->x, M->n*sizeof(double));    
    
    CHOLMOD()->cholmod_free_dense(&X, &CHOLMOD()->cholmod_common);
    CHOLMOD()->cholmod_free_dense(&Y, &CHOLMOD()->cholmod_common);
}

NLMatrix nlMatrixFactorize_CHOLMOD(
    NLMatrix M, NLenum solver
) {
    NLCholmodFactorizedMatrix* LLt = NULL;
    NLCRSMatrix* CRS = NULL;
    cholmod_sparse_ptr cM= NULL;
    NLuint nnz, cur, i, j, jj;
    int* rowptr = NULL;
    int* colind = NULL;
    double* val = NULL;
    NLuint n = M->n;

    nl_assert(solver == NL_CHOLMOD_EXT);
    nl_assert(M->m == M->n);
    
    if(M->type == NL_MATRIX_CRS) {
        CRS = (NLCRSMatrix*)M;
    } else if(M->type == NL_MATRIX_SPARSE_DYNAMIC) {
	/* 
	 * Note: since we convert once again into symmetric storage,
	 * we could also directly read the NLSparseMatrix there instead
	 * of copying once more...
	 */
        CRS = (NLCRSMatrix*)nlCRSMatrixNewFromSparseMatrix((NLSparseMatrix*)M);
    }

    LLt = NL_NEW(NLCholmodFactorizedMatrix);
    LLt->m = M->m;
    LLt->n = M->n;
    LLt->type = NL_MATRIX_OTHER;
    LLt->destroy_func = (NLDestroyMatrixFunc)(nlCholmodFactorizedMatrixDestroy);
    LLt->mult_func = (NLMultMatrixVectorFunc)(nlCholmodFactorizedMatrixMult);

    /*
     * Compute required nnz, if matrix is not already with symmetric storage,
     * ignore entries in the upper triangular part.
     */
    
    nnz=0;
    for(i=0; i<n; ++i) {
	for(jj=CRS->rowptr[i]; jj<CRS->rowptr[i+1]; ++jj) {
	    j=CRS->colind[jj];
	    if(j <= i) {
		++nnz;
	    }
	}
    }

    /*
     * Copy CRS matrix into CHOLDMOD matrix (and ignore upper trianglar part)
     */
    
    cM = CHOLMOD()->cholmod_allocate_sparse(
        n, n, nnz,    /* Dimensions and number of non-zeros */
        NL_FALSE,     /* Sorted = false */
        NL_TRUE,      /* Packed = true  */
        1,            /* stype (-1 = lower triangular, 1 = upper triangular) */
        CHOLMOD_REAL, /* Entries are real numbers */
        &CHOLMOD()->cholmod_common
    );

    rowptr = (int*)cM->p;
    colind = (int*)cM->i;
    val = (double*)cM->x;
    cur = 0;
    for(i=0; i<n; ++i) {
        rowptr[i] = (int)cur;
	for(jj=CRS->rowptr[i]; jj<CRS->rowptr[i+1]; ++jj) {
            j = CRS->colind[jj];
            if(j <= i) {
		val[cur] = CRS->val[jj];
		colind[cur] = (int)j;
		++cur;
            }
        }
    }
    rowptr[n] = (int)cur;
    nl_assert(cur==nnz);

    LLt->L = CHOLMOD()->cholmod_analyze(cM, &CHOLMOD()->cholmod_common);
    if(!CHOLMOD()->cholmod_factorize(cM, LLt->L, &CHOLMOD()->cholmod_common)) {
        CHOLMOD()->cholmod_free_factor(&LLt->L, &CHOLMOD()->cholmod_common);
	NL_DELETE(LLt);
    }
    
    CHOLMOD()->cholmod_free_sparse(&cM, &CHOLMOD()->cholmod_common);
    
    if((NLMatrix)CRS != M) {
        nlDeleteMatrix((NLMatrix)CRS);
    }

    return (NLMatrix)(LLt);
}


/******* extracted from nl_arpack.c *******/



#ifdef NL_OS_UNIX
#  ifdef NL_OS_APPLE
#      define ARPACK_LIB_NAME "libarpack.dylib"
#  else
#      define ARPACK_LIB_NAME "libarpack.so"
#  endif
#else
#  define ARPACK_LIB_NAME "libarpack.dll"
#endif


typedef int ARint;
typedef int ARlogical;


/* double precision symmetric routines */

typedef void (*FUNPTR_dsaupd)(
    ARint *ido, char *bmat, ARint *n, char *which,
    ARint *nev, double *tol, double *resid,
    ARint *ncv, double *V, ARint *ldv,
    ARint *iparam, ARint *ipntr, double *workd,
    double *workl, ARint *lworkl, ARint *info
);

typedef void (*FUNPTR_dseupd)(
    ARlogical *rvec, char *HowMny, ARlogical *select,
    double *d, double *Z, ARint *ldz,
    double *sigma, char *bmat, ARint *n,
    char *which, ARint *nev, double *tol,
    double *resid, ARint *ncv, double *V,
    ARint *ldv, ARint *iparam, ARint *ipntr,
    double *workd, double *workl,
    ARint *lworkl, ARint *info
);

/* double precision nonsymmetric routines */
    
typedef void (*FUNPTR_dnaupd)(
    ARint *ido, char *bmat, ARint *n, char *which,
    ARint *nev, double *tol, double *resid,
    ARint *ncv, double *V, ARint *ldv,
    ARint *iparam, ARint *ipntr, double *workd,
    double *workl, ARint *lworkl, ARint *info
);

typedef void (*FUNPTR_dneupd)(
    ARlogical *rvec, char *HowMny, ARlogical *select,
    double *dr, double *di, double *Z,
    ARint *ldz, double *sigmar,
    double *sigmai, double *workev,
    char *bmat, ARint *n, char *which,
    ARint *nev, double *tol, double *resid,
    ARint *ncv, double *V, ARint *ldv,
    ARint *iparam, ARint *ipntr,
    double *workd, double *workl,
    ARint *lworkl, ARint *info
);



typedef struct {
    FUNPTR_dsaupd dsaupd;
    FUNPTR_dseupd dseupd;
    FUNPTR_dnaupd dnaupd;
    FUNPTR_dneupd dneupd;
    NLdll DLL_handle;
} ARPACKContext;


static ARPACKContext* ARPACK() {
    static ARPACKContext context;
    static NLboolean init = NL_FALSE;
    if(!init) {
        init = NL_TRUE;
        memset(&context, 0, sizeof(context));
    }
    return &context;
}

NLboolean nlExtensionIsInitialized_ARPACK() {
    return
        ARPACK()->DLL_handle != NULL &&
        ARPACK()->dsaupd != NULL &&
        ARPACK()->dseupd != NULL &&   
        ARPACK()->dnaupd != NULL &&
        ARPACK()->dneupd != NULL;
}

static void nlTerminateExtension_ARPACK(void) {
    if(ARPACK()->DLL_handle != NULL) {
        nlCloseDLL(ARPACK()->DLL_handle);
        ARPACK()->DLL_handle = NULL;
    }
}


static char* u(const char* str) {
    static char buff[1000];
    sprintf(buff, "%s_", str);
    return buff;
}

#define find_arpack_func(name)                                             \
    if(                                                                    \
        (                                                                  \
            ARPACK()->name =                                               \
            (FUNPTR_##name)nlFindFunction(ARPACK()->DLL_handle,u(#name))   \
        ) == NULL                                                          \
    ) {                                                                    \
        nlError("nlInitExtension_ARPACK","function not found");            \
        nlError("nlInitExtension_ARPACK",u(#name));	   		   \
        return NL_FALSE;                                                   \
    }

NLboolean nlInitExtension_ARPACK(void) {
    NLenum flags = NL_LINK_NOW | NL_LINK_USE_FALLBACK;
    if(nlCurrentContext == NULL || !nlCurrentContext->verbose) {
	flags |= NL_LINK_QUIET;
    }
    
    if(ARPACK()->DLL_handle != NULL) {
        return nlExtensionIsInitialized_ARPACK();
    }
    
    ARPACK()->DLL_handle = nlOpenDLL(ARPACK_LIB_NAME, flags);
    if(ARPACK()->DLL_handle == NULL) {
        return NL_FALSE;
    }

    find_arpack_func(dsaupd);
    find_arpack_func(dseupd);
    find_arpack_func(dnaupd);
    find_arpack_func(dneupd);

    atexit(nlTerminateExtension_ARPACK);
    return NL_TRUE;
}



static NLMatrix create_OP(NLboolean symmetric) {
    NLuint n = nlCurrentContext->M->n;
    NLuint i;
    NLMatrix result = NULL;
    
	
    if(nlCurrentContext->eigen_shift != 0.0) {
	/*
	 * A = M
	 */
	NLSparseMatrix* A = NL_NEW(NLSparseMatrix);
	nlSparseMatrixConstruct(A, n, n, NL_MATRIX_STORE_ROWS);
	nlSparseMatrixAddMatrix(A, 1.0, nlCurrentContext->M);
	if(nlCurrentContext->B == NULL) {
	    /*
	     * A = A - shift * Id
	     */
	    for(i=0; i<n; ++i) {
		nlSparseMatrixAdd(A, i, i, -nlCurrentContext->eigen_shift);
	    }
	} else {
	    /*
	     * A = A - shift * B
	     */
	    nlSparseMatrixAddMatrix(
		A, -nlCurrentContext->eigen_shift, nlCurrentContext->B
	    );
	}

	/* 
	 * OP = A^{-1} 
	 */
	if(nlCurrentContext->verbose) {
	    nl_printf("Factorizing matrix...\n");
	}
	result = nlMatrixFactorize(
	    (NLMatrix)A,
	    symmetric ? NL_SYMMETRIC_SUPERLU_EXT : NL_PERM_SUPERLU_EXT
	);
	if(nlCurrentContext->verbose) {
	    if(result == NULL) {
		nl_printf("Could not factorize matrix\n");
	    } else {
		nl_printf("Matrix factorized\n");
	    }
	}
	nlDeleteMatrix((NLMatrix)A);
    } else {
	/* 
	 * OP = M^{-1} 
	 */
	if(nlCurrentContext->verbose) {
	    nl_printf("Factorizing matrix...\n");
	}
	result = nlMatrixFactorize(
	    nlCurrentContext->M,
	    symmetric ? NL_SYMMETRIC_SUPERLU_EXT : NL_PERM_SUPERLU_EXT
	);
	if(nlCurrentContext->verbose) {
	    if(result == NULL) {
		nl_printf("Could not factorize matrix\n");		
	    } else {
		nl_printf("Matrix factorized\n");
	    }
	}
    }

    if(result == NULL) {
	return NULL;
    }
    
    if(nlCurrentContext->B != NULL) {
	/* 
	 * OP = OP * B
	 */	
	result = nlMatrixNewFromProduct(
	    result, NL_TRUE, /* mem. ownership transferred */
	    nlCurrentContext->B, NL_FALSE  /* mem. ownership kept by context */
	);
    }

    return result;
}

static int eigencompare(const void* pi, const void* pj) {
    NLuint i = *(const NLuint*)pi;
    NLuint j = *(const NLuint*)pj;
    double vali = fabs(nlCurrentContext->temp_eigen_value[i]);
    double valj = fabs(nlCurrentContext->temp_eigen_value[j]);
    if(vali == valj) {
	return 0;
    }
    return vali < valj ? -1 : 1;
}

void nlEigenSolve_ARPACK(void) {
    NLboolean symmetric =
	nlCurrentContext->symmetric && (nlCurrentContext->B == NULL); 
    int n = (int)nlCurrentContext->M->n; /* Dimension of the matrix */
    int nev = /* Number of eigenvectors requested */
	(int)nlCurrentContext->nb_systems;
    NLMatrix OP = create_OP(symmetric);
    int ncv = (int)(nev * 2.5); /* Length of Arnoldi factorization */
                 /* Rule of thumb in ARPACK documentation: ncv > 2 * nev */
    int* iparam = NULL;
    int* ipntr  = NULL;
    NLdouble* resid = NULL;
    NLdouble* workev = NULL;
    NLdouble* workd = NULL;
    NLdouble* workl = NULL;
    NLdouble* v = NULL;
    NLdouble* d = NULL;
    ARlogical* select = NULL;
    ARlogical rvec = 1;
    double sigmar = 0.0;
    double sigmai = 0.0;
    int ierr;
    int i,k,kk;
    int ldv = (int)n;
    char* bmat = (char*)"I";   /*Standard problem */
    char* which = (char*)"LM"; /*Largest eigenvalues, but we invert->smallest */
    char* howmny = (char*)"A"; /*which eigens should be computed: all */
    double tol = nlCurrentContext->threshold;
    int ido = 0;  /* reverse communication variable (which operation ?) */
    int info = 1; /* start with initial value of resid */
    int lworkl;   /* size of work array */
    NLboolean converged = NL_FALSE;
    NLdouble value;
    int index;
    int* sorted; /* indirection array for sorting eigenpairs */

    if(OP == NULL) {
	nlError("nlEigenSolve_ARPACK","Could not factorize matrix");
	return;
    }
    
    if(ncv > n) {
	ncv = n;
    }

    if(nev > n) {
	nev = n;
    }

    if(nev + 2 > ncv) {
	nev = ncv  - 2;
    }

    
    if(symmetric) {
	lworkl = ncv * (ncv + 8) ;
    } else {
	lworkl = 3*ncv*ncv + 6*ncv ; 
    }
    iparam = NL_NEW_ARRAY(int, 11);
    ipntr  = NL_NEW_ARRAY(int, 14);

    iparam[1-1] = 1; /* ARPACK chooses the shifts */
    iparam[3-1] = (int)nlCurrentContext->max_iterations;
    iparam[7-1] = 1; /* Normal mode (we do not use 
         shift-invert (3) since we do our own shift-invert */

    workev = NL_NEW_ARRAY(NLdouble, 3*ncv);
    workd = NL_NEW_ARRAY(NLdouble, 3*n);

    resid = NL_NEW_ARRAY(NLdouble, n);
    for(i=0; i<n; ++i) {
	resid[i] = 1.0; /* (double)i / (double)n; */
    }
    v = NL_NEW_ARRAY(NLdouble, ldv*ncv);
    if(symmetric) {
	d = NL_NEW_ARRAY(NLdouble, 2*ncv);
    } else {
	d = NL_NEW_ARRAY(NLdouble, 3*ncv);	
    }
    workl = NL_NEW_ARRAY(NLdouble, lworkl);

    

    if(nlCurrentContext->verbose) {
	if(symmetric) {
	    nl_printf("calling dsaupd()\n");	    
	} else {
	    nl_printf("calling dnaupd()\n");
	}
    }
    while(!converged) {
	/*
	if(nlCurrentContext->verbose) {
	    fprintf(stderr, ".");
	    fflush(stderr);
	}
	*/
	if(symmetric) {
	    ARPACK()->dsaupd(
		&ido, bmat, &n, which, &nev, &tol, resid, &ncv,
		v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info
	    );
	} else {
	    ARPACK()->dnaupd(
		&ido, bmat, &n, which, &nev, &tol, resid, &ncv,
		v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info
	    );
	}
	if(ido == 1) {
	    nlMultMatrixVector(
             OP,
	     workd+ipntr[1-1]-1, /*The "-1"'s are for FORTRAN-to-C conversion */
	     workd+ipntr[2-1]-1  /*to keep the same indices as in ARPACK doc  */
	    );
	} else {
	    converged = NL_TRUE;
	}
    }

    

    if(info < 0) {
	if(symmetric) {
	    nl_fprintf(stderr, "\nError with dsaupd(): %d\n", info);	    
	} else {
	    nl_fprintf(stderr, "\nError with dnaupd(): %d\n", info);
	}
    } else {
	if(nlCurrentContext->verbose) {
	    fprintf(stderr, "\nconverged\n");
	}
	
	select = NL_NEW_ARRAY(ARlogical, ncv);
	for(i=0; i<ncv; ++i) {
	    select[i] = 1;
	}
	
	if(nlCurrentContext->verbose) {
	    if(symmetric) {
		nl_printf("calling dseupd()\n");		
	    } else {
		nl_printf("calling dneupd()\n");
	    }
	}
	
        if(symmetric) {
            ARPACK()->dseupd(
                &rvec, howmny, select, d, v, 
                &ldv, &sigmar, bmat, &n, which, &nev, 
                &tol, resid, &ncv, v, &ldv, 
                iparam, ipntr, workd,
		workl, &lworkl, &ierr 
	    );
        } else {
	    ARPACK()->dneupd(
		&rvec, howmny, select, d, d+ncv,
                v, &ldv, 
                &sigmar, &sigmai, workev, bmat, &n,
		which, &nev, &tol, 
                resid, &ncv, v, &ldv, iparam, 
		ipntr, workd, workl, &lworkl, &ierr 
            ) ;
	}	


	if(nlCurrentContext->verbose) {
	    if(ierr != 0) {		
		if(symmetric) {
		    nl_fprintf(stderr, "Error with dseupd(): %d\n", ierr);		
		} else {
		    nl_fprintf(stderr, "Error with dneupd(): %d\n", ierr);
		}
	    } else {
		if(symmetric) {
		    nl_printf("dseupd() OK, nconv= %d\n", iparam[3-1]);
		} else {
		    nl_printf("dneupd() OK, nconv= %d\n", iparam[3-1]);
		}
	    }
	}
	
	NL_DELETE_ARRAY(select);
    }

    

    for(i=0; i<nev; ++i) {
	d[i] = (fabs(d[i]) < 1e-30) ? 1e30 : 1.0 / d[i] ;
	d[i] += nlCurrentContext->eigen_shift ;
    }            

    
    
    /* Make it visible to the eigen_compare function */
    nlCurrentContext->temp_eigen_value = d;
    sorted = NL_NEW_ARRAY(int, nev);
    for(i=0; i<nev; ++i) {
	sorted[i] = i;
    }
    qsort(sorted, (size_t)nev, sizeof(NLuint), eigencompare);
    nlCurrentContext->temp_eigen_value = NULL;
    
    

    for(k=0; k<nev; ++k) {
	kk = sorted[k];
	nlCurrentContext->eigen_value[k] = d[kk];
	for(i=0; i<(int)nlCurrentContext->nb_variables; ++i) {
	    if(!nlCurrentContext->variable_is_locked[i]) {
		index = (int)nlCurrentContext->variable_index[i];
		nl_assert(index < n);
		value = v[kk*n+index];
		NL_BUFFER_ITEM(
		    nlCurrentContext->variable_buffer[k],(NLuint)i
		) = value;
	    }
	}
    }
    
    

    NL_DELETE_ARRAY(sorted);
    NL_DELETE_ARRAY(workl);
    NL_DELETE_ARRAY(d);
    NL_DELETE_ARRAY(v);
    NL_DELETE_ARRAY(resid);
    NL_DELETE_ARRAY(workd);
    NL_DELETE_ARRAY(workev);
    nlDeleteMatrix(OP);
    NL_DELETE_ARRAY(iparam);
    NL_DELETE_ARRAY(ipntr);
}




/******* extracted from nl_mkl.c *******/



typedef unsigned int MKL_INT;

typedef void (*FUNPTR_mkl_cspblas_dcsrgemv)(
    const char *transa, const MKL_INT *m, const double *a,
    const MKL_INT *ia, const MKL_INT *ja, const double *x, double *y
);

typedef void (*FUNPTR_mkl_cspblas_dcsrsymv)(
    const char *transa, const MKL_INT *m, const double *a,
    const MKL_INT *ia, const MKL_INT *ja, const double *x, double *y
 );

typedef struct {
    NLdll DLL_mkl_intel_lp64;
    NLdll DLL_mkl_intel_thread;
    NLdll DLL_mkl_core;
    NLdll DLL_iomp5;

    FUNPTR_mkl_cspblas_dcsrgemv mkl_cspblas_dcsrgemv;
    FUNPTR_mkl_cspblas_dcsrsymv mkl_cspblas_dcsrsymv;   
} MKLContext;

static MKLContext* MKL() {
    static MKLContext context;
    static NLboolean init = NL_FALSE;
    if(!init) {
        init = NL_TRUE;
        memset(&context, 0, sizeof(context));
    }
    return &context;
}

NLboolean nlExtensionIsInitialized_MKL() {
    if(
	MKL()->DLL_iomp5 == NULL ||
	MKL()->DLL_mkl_core == NULL ||
	MKL()->DLL_mkl_intel_thread == NULL ||
	MKL()->DLL_mkl_intel_lp64 == NULL ||
	MKL()->mkl_cspblas_dcsrgemv == NULL ||
	MKL()->mkl_cspblas_dcsrsymv == NULL 	
    ) {
        return NL_FALSE;
    }
    return NL_TRUE;
}

#define find_mkl_func(name)                                  \
    if(                                                      \
        (                                                    \
            MKL()->name =                                    \
            (FUNPTR_##name)nlFindFunction(                   \
		   MKL()->DLL_mkl_intel_lp64,#name           \
	    )					             \
        ) == NULL                                            \
    ) {                                                      \
        nlError("nlInitExtension_MKL","function not found"); \
        return NL_FALSE;                                     \
    }

static void nlTerminateExtension_MKL(void) {
    if(!nlExtensionIsInitialized_MKL()) {
	return;
    }
    nlCloseDLL(MKL()->DLL_mkl_intel_lp64);
    nlCloseDLL(MKL()->DLL_mkl_intel_thread);
    nlCloseDLL(MKL()->DLL_mkl_core);
    nlCloseDLL(MKL()->DLL_iomp5);
}

NLMultMatrixVectorFunc NLMultMatrixVector_MKL = NULL;

static void NLMultMatrixVector_MKL_impl(NLMatrix M_in, const double* x, double* y) {
    NLCRSMatrix* M = (NLCRSMatrix*)(M_in);
    nl_debug_assert(M_in->type == NL_MATRIX_CRS);
    if(M->symmetric_storage) {
	MKL()->mkl_cspblas_dcsrsymv(
	    "N", /* No transpose */
	    &M->m,
	    M->val,
	    M->rowptr,
	    M->colind,
	    x,
	    y
	);
    } else {
	MKL()->mkl_cspblas_dcsrgemv(
	    "N", /* No transpose */
	    &M->m,
	    M->val,
	    M->rowptr,
	    M->colind,
	    x,
	    y
	);
    }
}


#define INTEL_PREFIX "/opt/intel/"
#define LIB_DIR "lib/intel64/"
#define MKL_PREFIX  INTEL_PREFIX "mkl/" LIB_DIR

NLboolean nlInitExtension_MKL(void) {
    NLenum flags = NL_LINK_LAZY | NL_LINK_GLOBAL;
    if(nlCurrentContext == NULL || !nlCurrentContext->verbose) {
	flags |= NL_LINK_QUIET;
    }
    
    if(MKL()->DLL_mkl_intel_lp64 != NULL) {
        return nlExtensionIsInitialized_MKL();
    }
    
    MKL()->DLL_iomp5 = nlOpenDLL(
	INTEL_PREFIX LIB_DIR "libiomp5.so",
	flags
    );    
    MKL()->DLL_mkl_core = nlOpenDLL(
	MKL_PREFIX "libmkl_core.so",
	flags
    );    
    MKL()->DLL_mkl_intel_thread = nlOpenDLL(
	MKL_PREFIX "libmkl_intel_thread.so",
	flags
    );    
    MKL()->DLL_mkl_intel_lp64 = nlOpenDLL(
	MKL_PREFIX "libmkl_intel_lp64.so",
	flags
    );
    
    if(
	MKL()->DLL_iomp5 == NULL ||
	MKL()->DLL_mkl_core == NULL ||
	MKL()->DLL_mkl_intel_thread == NULL ||
	MKL()->DLL_mkl_intel_lp64 == NULL
    ) {
        return NL_FALSE;
    }

    find_mkl_func(mkl_cspblas_dcsrgemv);
    find_mkl_func(mkl_cspblas_dcsrsymv);

    if(nlExtensionIsInitialized_MKL()) {
	NLMultMatrixVector_MKL = NLMultMatrixVector_MKL_impl;
    }
    
    atexit(nlTerminateExtension_MKL);
    return NL_TRUE;
}




/******* extracted from nl_cuda.c *******/




/*      CUDA structures and functions                     */
/* Repeated here so that one can compile OpenNL without   */
/* requiring CUDA to be installed in the system.          */


struct cudaDeviceProp {
    char name[256];
    size_t totalGlobalMem;
    size_t sharedMemPerBlock;
    int regsPerBlock;
    int warpSize;
    size_t memPitch;
    int maxThreadsPerBlock;
    int maxThreadsDim[3];
    int maxGridSize[3];
    int clockRate;
    size_t totalConstMem;
    int major;
    int minor;
    size_t textureAlignment;
    size_t texturePitchAlignment;
    int deviceOverlap;
    int multiProcessorCount;
    int kernelExecTimeoutEnabled;
    int integrated;
    int canMapHostMemory;
    int computeMode;
    int maxTexture1D;
    int maxTexture1DMipmap;
    int maxTexture1DLinear;
    int maxTexture2D[2];
    int maxTexture2DMipmap[2];
    int maxTexture2DLinear[3];
    int maxTexture2DGather[2];
    int maxTexture3D[3];
    int maxTexture3DAlt[3];
    int maxTextureCubemap;
    int maxTexture1DLayered[2];
    int maxTexture2DLayered[3];
    int maxTextureCubemapLayered[2];
    int maxSurface1D;
    int maxSurface2D[2];
    int maxSurface3D[3];
    int maxSurface1DLayered[2];
    int maxSurface2DLayered[3];
    int maxSurfaceCubemap;
    int maxSurfaceCubemapLayered[2];
    size_t surfaceAlignment;
    int concurrentKernels;
    int ECCEnabled;
    int pciBusID;
    int pciDeviceID;
    int pciDomainID;
    int tccDriver;
    int asyncEngineCount;
    int unifiedAddressing;
    int memoryClockRate;
    int memoryBusWidth;
    int l2CacheSize;
    int maxThreadsPerMultiProcessor;
    int streamPrioritiesSupported;
    int globalL1CacheSupported;
    int localL1CacheSupported;
    size_t sharedMemPerMultiprocessor;
    int regsPerMultiprocessor;
    int managedMemSupported;
    int isMultiGpuBoard;
    int multiGpuBoardGroupID;
    int singleToDoublePrecisionPerfRatio;
    int pageableMemoryAccess;
    int concurrentManagedAccess;
    char padding[1024]; /* More room for future evolutions */
};

enum cudaComputeMode {
    cudaComputeModeDefault          = 0, 
    cudaComputeModeExclusive        = 1, 
    cudaComputeModeProhibited       = 2, 
    cudaComputeModeExclusiveProcess = 3   
};

enum cudaMemcpyKind {
    cudaMemcpyHostToHost          =   0, 
    cudaMemcpyHostToDevice        =   1, 
    cudaMemcpyDeviceToHost        =   2, 
    cudaMemcpyDeviceToDevice      =   3, 
    cudaMemcpyDefault             =   4  
};

typedef int cudaError_t;

typedef cudaError_t (*FUNPTR_cudaGetDeviceCount)(int* device_count);
typedef cudaError_t (*FUNPTR_cudaGetDeviceProperties)(
    struct cudaDeviceProp *props, int device
);
typedef cudaError_t (*FUNPTR_cudaDeviceReset)(void);
typedef cudaError_t (*FUNPTR_cudaMalloc)(void **devPtr, size_t size);
typedef cudaError_t (*FUNPTR_cudaFree)(void* devPtr);
typedef cudaError_t (*FUNPTR_cudaMemcpy)(
    void *dst, const void *src, size_t count, enum cudaMemcpyKind kind
);

#define find_cuda_func(name)                                  \
    if(                                                       \
        (                                                     \
            CUDA()->name =                                    \
            (FUNPTR_##name)nlFindFunction(                    \
		   CUDA()->DLL_cudart,#name                   \
	    )					              \
        ) == NULL                                             \
    ) {                                                       \
        nlError("nlInitExtension_CUDA: function not found", #name); \
        return NL_FALSE;                                      \
    }


/*      CUBLAS structures and functions                   */


struct cublasContext;
typedef struct cublasContext *cublasHandle_t;
typedef int cublasStatus_t;

typedef enum {
    CUBLAS_SIDE_LEFT =0, 
    CUBLAS_SIDE_RIGHT=1
} cublasSideMode_t; 

typedef enum {
    CUBLAS_FILL_MODE_LOWER=0, 
    CUBLAS_FILL_MODE_UPPER=1
} cublasFillMode_t;

typedef enum {
    CUBLAS_OP_N=0,  
    CUBLAS_OP_T=1,  
    CUBLAS_OP_C=2  
} cublasOperation_t;

typedef enum {
    CUBLAS_DIAG_NON_UNIT=0, 
    CUBLAS_DIAG_UNIT=1
} cublasDiagType_t; 

typedef cublasStatus_t (*FUNPTR_cublasCreate)(cublasHandle_t* handle);
typedef cublasStatus_t (*FUNPTR_cublasDestroy)(cublasHandle_t handle);

typedef cublasStatus_t (*FUNPTR_cublasGetVersion)(
    cublasHandle_t handle, int* version
);

typedef cublasStatus_t (*FUNPTR_cublasDdot)(
    cublasHandle_t handle, int n,
    const double *x, int incx,
    const double *y, int incy,
    double *result
);

typedef cublasStatus_t (*FUNPTR_cublasDcopy)(
    cublasHandle_t handle, int n,
    const double *x, int incx,
    const double *y, int incy
);

typedef cublasStatus_t (*FUNPTR_cublasDaxpy)(
    cublasHandle_t handle, int n,
    const double* alpha, 
    const double *x, int incx,
    const double *y, int incy
);

typedef cublasStatus_t (*FUNPTR_cublasDscal)(
    cublasHandle_t handle, int n,
    const double* alpha, 
    const double *x, int incx
);

typedef cublasStatus_t (*FUNPTR_cublasDnrm2)(
    cublasHandle_t handle, int n,
    const double *x, int incx,
    double* result
);
                                
typedef cublasStatus_t (*FUNPTR_cublasDdgmm)(
    cublasHandle_t handle, cublasSideMode_t mode,
    int m, int n,
    const double* A, int lda,
    const double* x, int incx,
    double* C, int ldc
);

typedef cublasStatus_t (*FUNPTR_cublasDgemv)(
    cublasHandle_t handle, 
    cublasOperation_t trans, 
    int m,
    int n,
    const double *alpha,
    const double *A,
    int lda,
    const double *x,
    int incx,
    const double *beta, 
    double *y, 
    int incy    
);

typedef cublasStatus_t (*FUNPTR_cublasDtpsv)(
    cublasHandle_t handle, cublasFillMode_t uplo,
    cublasOperation_t trans, cublasDiagType_t diag,
    int n, const double *AP,
    double* x, int incx
);


#define find_cublas_func(name)	    		              \
    if(                                                       \
        (                                                     \
            CUDA()->name =                                    \
            (FUNPTR_##name)nlFindFunction(                    \
		   CUDA()->DLL_cublas,#name "_v2"             \
	    )					              \
        ) == NULL                                             \
    ) {                                                       \
        nlError("nlInitExtension_CUDA: function not found", #name); \
        return NL_FALSE;                                      \
    }

#define find_cublas_func_v1(name)                             \
    if(                                                       \
        (                                                     \
            CUDA()->name =                                    \
            (FUNPTR_##name)nlFindFunction(                    \
		   CUDA()->DLL_cublas,#name                   \
	    )					              \
        ) == NULL                                             \
    ) {                                                       \
        nlError("nlInitExtension_CUDA: function not found", #name); \
        return NL_FALSE;                                      \
    }




/*      CUSPARSE structures and functions                 */


struct cusparseContext;
typedef struct cusparseContext *cusparseHandle_t;
typedef int cusparseStatus_t;
struct cusparseMatDescr;
typedef struct cusparseMatDescr *cusparseMatDescr_t;

typedef enum {
    CUSPARSE_MATRIX_TYPE_GENERAL = 0, 
    CUSPARSE_MATRIX_TYPE_SYMMETRIC = 1,     
    CUSPARSE_MATRIX_TYPE_HERMITIAN = 2, 
    CUSPARSE_MATRIX_TYPE_TRIANGULAR = 3 
} cusparseMatrixType_t;

typedef enum {
    CUSPARSE_INDEX_BASE_ZERO = 0, 
    CUSPARSE_INDEX_BASE_ONE = 1
} cusparseIndexBase_t;

typedef enum {
    CUSPARSE_OPERATION_NON_TRANSPOSE = 0,  
    CUSPARSE_OPERATION_TRANSPOSE = 1,  
    CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE = 2  
} cusparseOperation_t;

struct cusparseHybMat;
typedef struct cusparseHybMat *cusparseHybMat_t;

typedef enum {
    CUSPARSE_HYB_PARTITION_AUTO = 0,  
    CUSPARSE_HYB_PARTITION_USER = 1,  
    CUSPARSE_HYB_PARTITION_MAX = 2    
} cusparseHybPartition_t;

typedef cusparseStatus_t (*FUNPTR_cusparseCreate)(cusparseHandle_t* handle);
typedef cusparseStatus_t (*FUNPTR_cusparseDestroy)(cusparseHandle_t handle);
typedef cusparseStatus_t (*FUNPTR_cusparseGetVersion)(
    cusparseHandle_t handle, int* version
);
typedef cusparseStatus_t (*FUNPTR_cusparseCreateMatDescr)(
    cusparseMatDescr_t* descr
);
typedef cusparseStatus_t (*FUNPTR_cusparseDestroyMatDescr)(
    cusparseMatDescr_t descr
);
typedef cusparseStatus_t (*FUNPTR_cusparseSetMatType)(
    cusparseMatDescr_t descr, cusparseMatrixType_t mtype
);
typedef cusparseStatus_t (*FUNPTR_cusparseSetMatIndexBase)(
    cusparseMatDescr_t descr, cusparseIndexBase_t ibase
);
typedef cusparseStatus_t (*FUNPTR_cusparseDcsrmv)(
     cusparseHandle_t handle, cusparseOperation_t transA, 
     int m, int n, int nnz,
     const double *alpha, const cusparseMatDescr_t descrA, 
     const double *csrSortedValA, const int *csrSortedRowPtrA, 
     const int *csrSortedColIndA, const double *x, 
     const double *beta, double *y
);
typedef cusparseStatus_t (*FUNPTR_cusparseCreateHybMat)(
    cusparseHybMat_t *hybA
);
typedef cusparseStatus_t (*FUNPTR_cusparseDestroyHybMat)(
    cusparseHybMat_t hybA
);
typedef cusparseStatus_t (*FUNPTR_cusparseDcsr2hyb)(
    cusparseHandle_t handle,
    int m,
    int n,
    const cusparseMatDescr_t descrA,
    const double *csrSortedValA,
    const int *csrSortedRowPtrA,
    const int *csrSortedColIndA,
    cusparseHybMat_t hybA,
    int userEllWidth,
    cusparseHybPartition_t partitionType
);
typedef cusparseStatus_t (*FUNPTR_cusparseDhybmv)(
    cusparseHandle_t handle,
    cusparseOperation_t transA,
    const double *alpha,
    const cusparseMatDescr_t descrA,
    const cusparseHybMat_t hybA,
    const double *x,
    const double *beta,
    double *y
);

#define find_cusparse_func(name)                              \
    if(                                                       \
        (                                                     \
            CUDA()->name =                                    \
            (FUNPTR_##name)nlFindFunction(                    \
		   CUDA()->DLL_cusparse,#name                 \
	    )					              \
        ) == NULL                                             \
    ) {                                                       \
        nlError("nlInitExtension_CUDA : function not found", #name);	\
        return NL_FALSE;                                      \
    }




typedef struct {
    NLdll DLL_cudart;
    FUNPTR_cudaGetDeviceCount cudaGetDeviceCount;
    FUNPTR_cudaGetDeviceProperties cudaGetDeviceProperties;
    FUNPTR_cudaDeviceReset cudaDeviceReset;
    FUNPTR_cudaMalloc cudaMalloc;
    FUNPTR_cudaFree cudaFree;
    FUNPTR_cudaMemcpy cudaMemcpy;

    NLdll DLL_cublas;
    cublasHandle_t HNDL_cublas;
    FUNPTR_cublasCreate cublasCreate;
    FUNPTR_cublasDestroy cublasDestroy;
    FUNPTR_cublasGetVersion cublasGetVersion;
    FUNPTR_cublasDdot cublasDdot;
    FUNPTR_cublasDcopy cublasDcopy;
    FUNPTR_cublasDaxpy cublasDaxpy;
    FUNPTR_cublasDscal cublasDscal;    
    FUNPTR_cublasDnrm2 cublasDnrm2;
    FUNPTR_cublasDdgmm cublasDdgmm;
    FUNPTR_cublasDgemv cublasDgemv;
    FUNPTR_cublasDtpsv cublasDtpsv;
    
    NLdll DLL_cusparse;
    cusparseHandle_t HNDL_cusparse;
    FUNPTR_cusparseCreate cusparseCreate;
    FUNPTR_cusparseDestroy cusparseDestroy;
    FUNPTR_cusparseGetVersion cusparseGetVersion;
    FUNPTR_cusparseCreateMatDescr cusparseCreateMatDescr;
    FUNPTR_cusparseDestroyMatDescr cusparseDestroyMatDescr;    
    FUNPTR_cusparseSetMatType cusparseSetMatType;
    FUNPTR_cusparseSetMatIndexBase cusparseSetMatIndexBase;    
    FUNPTR_cusparseDcsrmv cusparseDcsrmv;
    FUNPTR_cusparseCreateHybMat cusparseCreateHybMat;
    FUNPTR_cusparseDestroyHybMat cusparseDestroyHybMat;
    FUNPTR_cusparseDcsr2hyb cusparseDcsr2hyb;
    FUNPTR_cusparseDhybmv cusparseDhybmv;
    
    int devID;
} CUDAContext;

static CUDAContext* CUDA() {
    static CUDAContext context;
    static NLboolean init = NL_FALSE;
    if(!init) {
        init = NL_TRUE;
        memset(&context, 0, sizeof(context));
    }
    return &context;
}

NLboolean nlExtensionIsInitialized_CUDA() {
    if(
	CUDA()->DLL_cudart == NULL ||
	CUDA()->cudaGetDeviceCount == NULL ||
	CUDA()->cudaGetDeviceProperties == NULL ||
	CUDA()->cudaDeviceReset == NULL ||
	CUDA()->cudaMalloc == NULL ||
	CUDA()->cudaFree == NULL ||
	CUDA()->cudaMemcpy == NULL ||
	
	CUDA()->DLL_cublas == NULL ||
	CUDA()->HNDL_cublas == NULL ||
	CUDA()->cublasCreate == NULL ||
	CUDA()->cublasDestroy == NULL ||
	CUDA()->cublasGetVersion == NULL ||
	CUDA()->cublasDdot == NULL ||
	CUDA()->cublasDcopy == NULL ||
	CUDA()->cublasDaxpy == NULL ||
	CUDA()->cublasDscal == NULL ||
	CUDA()->cublasDnrm2 == NULL ||
	CUDA()->cublasDdgmm == NULL ||
	
	CUDA()->DLL_cusparse == NULL ||
	CUDA()->HNDL_cusparse == NULL ||
	CUDA()->cusparseCreate == NULL ||
	CUDA()->cusparseDestroy == NULL ||
	CUDA()->cusparseGetVersion == NULL ||
	CUDA()->cusparseCreateMatDescr == NULL ||
	CUDA()->cusparseDestroyMatDescr == NULL ||	
	CUDA()->cusparseSetMatType == NULL ||
	CUDA()->cusparseSetMatIndexBase == NULL ||
	CUDA()->cusparseDcsrmv == NULL ||
	CUDA()->cusparseCreateHybMat == NULL ||
	CUDA()->cusparseDestroyHybMat == NULL ||
	CUDA()->cusparseDcsr2hyb == NULL ||
	CUDA()->cusparseDhybmv == NULL
    ) {
        return NL_FALSE;
    }
    return NL_TRUE;
}

static void nlTerminateExtension_CUDA(void) {
    if(!nlExtensionIsInitialized_CUDA()) {
	return;
    }

    CUDA()->cusparseDestroy(CUDA()->HNDL_cusparse);    
    nlCloseDLL(CUDA()->DLL_cusparse);
    
    CUDA()->cublasDestroy(CUDA()->HNDL_cublas);
    nlCloseDLL(CUDA()->DLL_cublas);

    CUDA()->cudaDeviceReset();    
    nlCloseDLL(CUDA()->DLL_cudart);
}



static int ConvertSMVer2Cores(int major, int minor) {
    /* Defines for GPU Architecture types (using the SM version 
       to determine the # of cores per SM */
    typedef struct {
        int SM; /* 0xMm (hexadecimal notation), 
                    M = SM Major version, 
                    and m = SM minor version */
        int Cores;
    } sSMtoCores;

    sSMtoCores nGpuArchCoresPerSM[] = {
        { 0x10,  8 }, /* Tesla Generation   (SM 1.0) G80 class    */
        { 0x11,  8 }, /* Tesla Generation   (SM 1.1) G8x class    */
        { 0x12,  8 }, /* Tesla Generation   (SM 1.2) G9x class    */
        { 0x13,  8 }, /* Tesla Generation   (SM 1.3) GT200 class  */
        { 0x20, 32 }, /* Fermi Generation   (SM 2.0) GF100 class  */
        { 0x21, 48 }, /* Fermi Generation   (SM 2.1) GF10x class  */
        { 0x30, 192}, /* Kepler Generation  (SM 3.0) GK10x class  */
        { 0x35, 192}, /* Kepler Generation  (SM 3.5) GK11x class  */
	{ 0x50, 128}, /* Maxwell Generation (SM 5.0) GM10x class 
                             (yes, #cores smaller than with 3.x)  */
	{ 0x52, 128}, /* Maxwell Generation (SM 5.2) GM20x class  */
	{ 0x60, 64 }, /* Pascal Generation  (SM 6.0) GP100,GP102  
              (yes, 64, but GP100 has superfast double precision) */
	{ 0x61, 128}, /* Pascal Generation  (SM 6.1) GP104 class  
                               (but FP64 runs as 1/32 FP32 speed) */ 	
        {   -1, -1 }
    };
    int index = 0;
    while (nGpuArchCoresPerSM[index].SM != -1) {
        if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor)) {
            return nGpuArchCoresPerSM[index].Cores;
        }
        index++;
    }
    /* If we don't find the values, we default use the 
       previous one to run properly */
    nl_printf(
      "MapSMtoCores for SM %d.%d is undefined.  Default to use %d Cores/SM\n",
      major, minor, nGpuArchCoresPerSM[8].Cores
    );
    return nGpuArchCoresPerSM[8].Cores;
}

static int getBestDeviceID() {
    int current_device     = 0, sm_per_multiproc  = 0;
    int max_compute_perf   = 0, max_perf_device   = 0;
    int device_count       = 0, best_SM_arch      = 0;
    int compute_perf       = 0;
    struct cudaDeviceProp deviceProp;
    CUDA()->cudaGetDeviceCount(&device_count);
    /* Find the best major SM Architecture GPU device */
    while (current_device < device_count) {
        CUDA()->cudaGetDeviceProperties(&deviceProp, current_device);
        /* If this GPU is not running on Compute Mode prohibited, 
           then we can add it to the list */
        if (deviceProp.computeMode != cudaComputeModeProhibited) {
            if (deviceProp.major > 0 && deviceProp.major < 9999) {
                best_SM_arch = MAX(best_SM_arch, deviceProp.major);
            }
        }
        current_device++;
    }
    /* Find the best CUDA capable GPU device */
    current_device = 0;
    while (current_device < device_count) {
        CUDA()->cudaGetDeviceProperties(&deviceProp, current_device);
        /* If this GPU is not running on Compute Mode prohibited, 
           then we can add it to the list */
        if (deviceProp.computeMode != cudaComputeModeProhibited) {
            if (deviceProp.major == 9999 && deviceProp.minor == 9999) {
                sm_per_multiproc = 1;
            } else {
                sm_per_multiproc = ConvertSMVer2Cores(
		    deviceProp.major, deviceProp.minor
		);
            }
            compute_perf  =
		deviceProp.multiProcessorCount *
		sm_per_multiproc * deviceProp.clockRate;
            if (compute_perf  > max_compute_perf) {
                /* If we find GPU with SM major > 2, search only these */
                if (best_SM_arch > 2) {
                    /* If our device==dest_SM_arch, choose this, or else pass */
                    if (deviceProp.major == best_SM_arch) {
                        max_compute_perf  = compute_perf;
                        max_perf_device   = current_device;
                    }
                } else {
                    max_compute_perf  = compute_perf;
                    max_perf_device   = current_device;
                }
            }
        }
        ++current_device;
    }

 
    return max_perf_device;
}




#ifdef NL_OS_UNIX
#  define LIBPREFIX "lib"
#  ifdef NL_OS_APPLE
#      define LIBEXTENSION ".dylib"
#  else
#      define LIBEXTENSION ".so"
#  endif
#else
#  define LIBPREFIX 
#  define LIBEXTENSION ".dll"
#endif


NLboolean nlInitExtension_CUDA(void) {
    struct cudaDeviceProp deviceProp;
    int cublas_version;
    int cusparse_version;
    NLenum flags = NL_LINK_LAZY | NL_LINK_GLOBAL;
    if(nlCurrentContext == NULL || !nlCurrentContext->verbose) {
	flags |= NL_LINK_QUIET;
    }
    
    if(nlExtensionIsInitialized_CUDA()) {
	return NL_TRUE;
    }

    CUDA()->DLL_cudart = nlOpenDLL(
	LIBPREFIX "cudart" LIBEXTENSION, flags
    );

    find_cuda_func(cudaGetDeviceCount);
    find_cuda_func(cudaGetDeviceProperties);
    find_cuda_func(cudaDeviceReset);        
    find_cuda_func(cudaMalloc);
    find_cuda_func(cudaFree);
    find_cuda_func(cudaMemcpy);
    
    CUDA()->devID = getBestDeviceID();

    if(CUDA()->cudaGetDeviceProperties(&deviceProp, CUDA()->devID)) {
	nl_fprintf(stderr,"OpenNL CUDA: could not find a CUDA device\n");
	return NL_FALSE;
    }
    
    nl_printf("OpenNL CUDA: Device ID = %d\n", CUDA()->devID);
    nl_printf("OpenNL CUDA: Device name=%s\n", deviceProp.name);
    nl_printf(
	"OpenNL CUDA: Device has %d Multi-Processors, "
	"%d cores per Multi-Processor, SM %d.%d compute capabilities\n",
	deviceProp.multiProcessorCount,
	ConvertSMVer2Cores(deviceProp.major, deviceProp.minor),
	deviceProp.major, deviceProp.minor
    );

    nl_printf(
	"OpenNL CUDA: %d kB shared mem. per block, %d per MP\n",
	(int)(deviceProp.sharedMemPerBlock / 1024),
	(int)(deviceProp.sharedMemPerMultiprocessor / 1024)
    );
    
    nl_printf(
	"OpenNL CUDA: %d regs. per block, %d per MP\n",
	deviceProp.regsPerBlock,
	deviceProp.regsPerMultiprocessor	
    );

    nl_printf(
	"OpenNL CUDA: warpsize=%d\n",
	deviceProp.warpSize
    );
    
    if ((deviceProp.major * 0x10 + deviceProp.minor) < 0x11) {
        nl_fprintf(stderr, "OpenNL CUDA requires a minimum CUDA compute 1.1 capability\n");
        CUDA()->cudaDeviceReset();
	return NL_FALSE;
    }
    
    CUDA()->DLL_cublas = nlOpenDLL(
	LIBPREFIX "cublas" LIBEXTENSION, flags
    );

    find_cublas_func(cublasCreate);
    find_cublas_func(cublasDestroy);
    find_cublas_func(cublasGetVersion);    
    find_cublas_func(cublasDdot);
    find_cublas_func(cublasDaxpy);
    find_cublas_func(cublasDcopy);
    find_cublas_func(cublasDscal);
    find_cublas_func(cublasDnrm2);
    find_cublas_func(cublasDgemv);        
    find_cublas_func(cublasDtpsv);    
    find_cublas_func_v1(cublasDdgmm);

    
    if(CUDA()->cublasCreate(&CUDA()->HNDL_cublas)) {
	return NL_FALSE;
    }

    if(CUDA()->cublasGetVersion(CUDA()->HNDL_cublas, &cublas_version)) {
	return NL_FALSE;
    }
    nl_printf("OpenNL CUDA: cublas version = %d\n", cublas_version);
    
    CUDA()->DLL_cusparse = nlOpenDLL(
	LIBPREFIX "cusparse" LIBEXTENSION, flags
    );
    find_cusparse_func(cusparseCreate);
    find_cusparse_func(cusparseDestroy);
    find_cusparse_func(cusparseGetVersion);
    find_cusparse_func(cusparseCreateMatDescr);
    find_cusparse_func(cusparseDestroyMatDescr);    
    find_cusparse_func(cusparseSetMatType);
    find_cusparse_func(cusparseSetMatIndexBase);
    find_cusparse_func(cusparseDcsrmv);
    find_cusparse_func(cusparseCreateHybMat);
    find_cusparse_func(cusparseDestroyHybMat);
    find_cusparse_func(cusparseDcsr2hyb);
    find_cusparse_func(cusparseDhybmv);                    
    
    if(CUDA()->cusparseCreate(&CUDA()->HNDL_cusparse)) {
	return NL_FALSE;
    }
    if(CUDA()->cusparseGetVersion(CUDA()->HNDL_cusparse, &cusparse_version)) {
	return NL_FALSE;
    }
    nl_printf("OpenNL CUDA: cusparse version = %d\n", cusparse_version);
    
    if(!nlExtensionIsInitialized_CUDA()) {
	return NL_FALSE;
    }

    atexit(nlTerminateExtension_CUDA);
    return NL_TRUE;
    
}

static void nlCUDACheckImpl(int status, int line) {
    if(status != 0) {
	nl_fprintf(stderr,"nl_cuda.c:%d fatal error %d\n",line, status);
	CUDA()->cudaDeviceReset();    	
	exit(-1);
    }
}

#define nlCUDACheck(status) nlCUDACheckImpl(status, __LINE__)



typedef struct {
    NLuint m;
    NLuint n;
    NLenum type;
    NLDestroyMatrixFunc destroy_func;
    NLMultMatrixVectorFunc mult_func;
    cusparseMatDescr_t descr;
    NLuint nnz;
    int* colind;
    int* rowptr;
    double* val;
    cusparseHybMat_t hyb;
} NLCUDASparseMatrix;


static void nlCRSMatrixCUDADestroyCRS(NLCUDASparseMatrix* Mcuda) {
    if(Mcuda->colind != NULL) {
	nlCUDACheck(CUDA()->cudaFree(Mcuda->colind));
	Mcuda->colind = NULL;
    }
    if(Mcuda->rowptr != NULL) {
	nlCUDACheck(CUDA()->cudaFree(Mcuda->rowptr));
	Mcuda->rowptr = NULL;
    }
    if(Mcuda->val != NULL) {
	nlCUDACheck(CUDA()->cudaFree(Mcuda->val));
	Mcuda->val = NULL;
    }
}

static void nlCRSMatrixCUDADestroy(NLCUDASparseMatrix* Mcuda) {
    if(Mcuda->hyb != NULL) {
	nlCUDACheck(CUDA()->cusparseDestroyHybMat(Mcuda->hyb));
    }
    nlCRSMatrixCUDADestroyCRS(Mcuda);
    nlCUDACheck(CUDA()->cusparseDestroyMatDescr(Mcuda->descr));
    memset(Mcuda, 0, sizeof(*Mcuda));
}

static void nlCRSMatrixCUDAMult(
    NLCUDASparseMatrix* Mcuda, const double* x, double* y
) {
    const double one = 1;
    const double zero = 0;
    if(Mcuda->hyb != NULL) {
	nlCUDACheck(
	    CUDA()->cusparseDhybmv(
		CUDA()->HNDL_cusparse,
		CUSPARSE_OPERATION_NON_TRANSPOSE,
		&one,
		Mcuda->descr,
		Mcuda->hyb,
		x,
		&zero,
		y
	    )
	);
    } else {
	nlCUDACheck(
	    CUDA()->cusparseDcsrmv(
		CUDA()->HNDL_cusparse,
		CUSPARSE_OPERATION_NON_TRANSPOSE,
		(int)Mcuda->m,
		(int)Mcuda->n,
		(int)Mcuda->nnz,
		&one,
		Mcuda->descr,
		Mcuda->val,
		Mcuda->rowptr,
		Mcuda->colind,
		x,
		&zero,
		y
		)
	    );
    }
    nlCUDABlas()->flops += (NLulong)(2*Mcuda->nnz);
}

NLMatrix nlCUDAMatrixNewFromCRSMatrix(NLMatrix M_in) {
    NLCUDASparseMatrix* Mcuda = NL_NEW(NLCUDASparseMatrix);
    NLCRSMatrix* M = (NLCRSMatrix*)(M_in);
    size_t colind_sz, rowptr_sz, val_sz;
    nl_assert(M_in->type == NL_MATRIX_CRS);
    nlCUDACheck(CUDA()->cusparseCreateMatDescr(&Mcuda->descr));
    if(M->symmetric_storage) {
	nlCUDACheck(CUDA()->cusparseSetMatType(
			Mcuda->descr, CUSPARSE_MATRIX_TYPE_SYMMETRIC)
	);
    } else {
	nlCUDACheck(CUDA()->cusparseSetMatType(
			Mcuda->descr, CUSPARSE_MATRIX_TYPE_GENERAL)
	);	
    }
    nlCUDACheck(CUDA()->cusparseSetMatIndexBase(
		    Mcuda->descr, CUSPARSE_INDEX_BASE_ZERO)
    );	
    Mcuda->m = M->m;
    Mcuda->n = M->n;
    Mcuda->nnz = nlCRSMatrixNNZ(M);

    colind_sz = (size_t)Mcuda->nnz*sizeof(int);
    rowptr_sz = (size_t)(Mcuda->m+1)*sizeof(int);
    val_sz    = (size_t)Mcuda->nnz*sizeof(double);

    nlCUDACheck(CUDA()->cudaMalloc((void**)&Mcuda->colind,colind_sz));
    nlCUDACheck(CUDA()->cudaMalloc((void**)&Mcuda->rowptr,rowptr_sz));
    nlCUDACheck(CUDA()->cudaMalloc((void**)&Mcuda->val,val_sz));
    nlCUDACheck(CUDA()->cudaMemcpy(
      Mcuda->colind, M->colind, colind_sz, cudaMemcpyHostToDevice)
    );
    nlCUDACheck(CUDA()->cudaMemcpy(
      Mcuda->rowptr, M->rowptr, rowptr_sz, cudaMemcpyHostToDevice)
    );
    nlCUDACheck(CUDA()->cudaMemcpy(
      Mcuda->val, M->val, val_sz, cudaMemcpyHostToDevice)
    );
    Mcuda->hyb=NULL;
    if(!M->symmetric_storage) {
	nlCUDACheck(CUDA()->cusparseCreateHybMat(&Mcuda->hyb));
	nlCUDACheck(CUDA()->cusparseDcsr2hyb(
			CUDA()->HNDL_cusparse,
			(int)M->m,
			(int)M->n,
			Mcuda->descr,
			Mcuda->val,
			Mcuda->rowptr,
			Mcuda->colind,
			Mcuda->hyb,
			0,
			CUSPARSE_HYB_PARTITION_AUTO 
	));
	/* We no longer need the CRS part */
	nlCRSMatrixCUDADestroyCRS(Mcuda);	
    }
    Mcuda->type=NL_MATRIX_OTHER;
    Mcuda->destroy_func=(NLDestroyMatrixFunc)nlCRSMatrixCUDADestroy;
    Mcuda->mult_func=(NLMultMatrixVectorFunc)nlCRSMatrixCUDAMult;
    return (NLMatrix)Mcuda;
}



typedef struct {
    NLuint m;
    NLuint n;
    NLenum type;
    NLDestroyMatrixFunc destroy_func;
    NLMultMatrixVectorFunc mult_func;
    double* val;
} NLDiagonalMatrixCUDA;

static void nlDiagonalMatrixCUDADestroy(NLDiagonalMatrixCUDA* Mcuda) {
    nlCUDACheck(CUDA()->cudaFree(Mcuda->val));
    memset(Mcuda, 0, sizeof(*Mcuda));
}

static void nlDiagonalMatrixCUDAMult(
    NLDiagonalMatrixCUDA* Mcuda, const double* x, double* y
) {
    int N = (int)Mcuda->n;
    /*
     * vector x vector component-wise product implemented
     * using diagonal matrix x matrix function.
     */
    nlCUDACheck(CUDA()->cublasDdgmm(
	CUDA()->HNDL_cublas, CUBLAS_SIDE_LEFT,
	N, 1,
	x, N,
	Mcuda->val, 1,
	y, N
    ));
    nlCUDABlas()->flops += (NLulong)N;
}

static NLMatrix nlDiagonalMatrixCUDANew(const double* diag, NLuint n) {
    NLDiagonalMatrixCUDA* Mcuda = NL_NEW(NLDiagonalMatrixCUDA);
    Mcuda->m = n;
    Mcuda->n = n;
    Mcuda->type = NL_MATRIX_OTHER;
    nlCUDACheck(CUDA()->cudaMalloc(
       (void**)&Mcuda->val, n*sizeof(double))
    );
    nlCUDACheck(CUDA()->cudaMemcpy(
       Mcuda->val, diag, n*sizeof(double), cudaMemcpyHostToDevice)
    );
    Mcuda->destroy_func=(NLDestroyMatrixFunc)nlDiagonalMatrixCUDADestroy;
    Mcuda->mult_func=(NLMultMatrixVectorFunc)nlDiagonalMatrixCUDAMult;
    return (NLMatrix)Mcuda;
}

NLMatrix nlCUDAJacobiPreconditionerNewFromCRSMatrix(NLMatrix M_in) {
    NLuint N = M_in->n;
    NLuint i,jj;
    double* diag = NULL;
    NLMatrix result = NULL;
    NLCRSMatrix* M = (NLCRSMatrix*)(M_in);
    nl_assert(M_in->type == NL_MATRIX_CRS);
    diag = NL_NEW_ARRAY(double,N);
    for(i=0; i<N; ++i) {
	for(jj=M->rowptr[i]; jj<M->rowptr[i+1]; ++jj) {
	    if(M->colind[jj] == i) {
		diag[i] = M->val[jj];
	    }
	}
    }
    for(i=0; i<N; ++i) {
	diag[i] = ((diag[i] == 0.0) ? 1.0 : 1.0 / diag[i]);
    }
    result = nlDiagonalMatrixCUDANew(diag, N);
    NL_DELETE_ARRAY(diag);
    return result;
}



static void* cuda_blas_malloc(
    NLBlas_t blas, NLmemoryType type, size_t size
) {
    void* result = NULL;
    blas->used_ram[type] += (NLulong)size;
    blas->max_used_ram[type] = MAX(
	blas->max_used_ram[type],blas->used_ram[type]
    );
    if(type == NL_HOST_MEMORY) {
	result = malloc(size);
    } else {
	nlCUDACheck(CUDA()->cudaMalloc(&result,size));	
    }
    return result;
}

static void cuda_blas_free(
    NLBlas_t blas, NLmemoryType type, size_t size, void* ptr
) {
    blas->used_ram[type] -= (NLulong)size;
    if(type == NL_HOST_MEMORY) {
	free(ptr);
    } else {
	nlCUDACheck(CUDA()->cudaFree(ptr));	
    }
}

static void cuda_blas_memcpy(
    NLBlas_t blas,
    void* to, NLmemoryType to_type,
    void* from, NLmemoryType from_type,
    size_t size
) {
    enum cudaMemcpyKind kind = cudaMemcpyDefault;
    nl_arg_used(blas);
    if(from_type == NL_HOST_MEMORY) {
	if(to_type == NL_HOST_MEMORY) {
	    kind = cudaMemcpyHostToHost;
	} else {
	    kind = cudaMemcpyHostToDevice;	    
	}
    } else {
	if(to_type == NL_HOST_MEMORY) {
	    kind = cudaMemcpyDeviceToHost;
	} else {
	    kind = cudaMemcpyDeviceToDevice;	    
	}
    }
    nlCUDACheck(CUDA()->cudaMemcpy(to, from, size, kind));
}

static void cuda_blas_dcopy(
    NLBlas_t blas, int n, const double *x, int incx, double *y, int incy    
) {
    nl_arg_used(blas);
    CUDA()->cublasDcopy(CUDA()->HNDL_cublas,n,x,incx,y,incy);
}

static double cuda_blas_ddot(
    NLBlas_t blas, int n, const double *x, int incx, const double *y, int incy
) {
    double result = 0.0;
    blas->flops += (NLulong)(2*n);
    CUDA()->cublasDdot(CUDA()->HNDL_cublas,n,x,incx,y,incy,&result);
    return result;
}

static double cuda_blas_dnrm2(
    NLBlas_t blas, int n, const double *x, int incx
) {
    double result = 0.0;
    blas->flops += (NLulong)(2*n);
    CUDA()->cublasDnrm2(CUDA()->HNDL_cublas,n,x,incx,&result);    
    return result;
}

static void cuda_blas_daxpy(
    NLBlas_t blas, int n,
    double a, const double *x, int incx, double *y, int incy
) {
    blas->flops += (NLulong)(2*n);
    CUDA()->cublasDaxpy(CUDA()->HNDL_cublas,n,&a,x,incx,y,incy);        
}

static void cuda_blas_dscal(
    NLBlas_t blas, int n, double a, double *x, int incx    
) {
    blas->flops += (NLulong)n;
    CUDA()->cublasDscal(CUDA()->HNDL_cublas,n,&a,x,incx);            
}


static void cuda_blas_dgemv(
    NLBlas_t blas, MatrixTranspose trans, int m, int n, double alpha,
    const double *A, int ldA, const double *x, int incx,
    double beta, double *y, int incy 
) {
    nl_arg_used(blas);
    /* TODO: update FLOPS */
    CUDA()->cublasDgemv(
	CUDA()->HNDL_cublas, (cublasOperation_t)trans,
	m, n, &alpha, A, ldA, x, incx, &beta, y, incy
    );
}

static void cuda_blas_dtpsv(
    NLBlas_t blas, MatrixTriangle uplo, MatrixTranspose trans,
    MatrixUnitTriangular diag, int n, const double *AP,
    double *x, int incx 
) {
    nl_arg_used(blas);
    /* TODO: update FLOPS */
    CUDA()->cublasDtpsv(
	CUDA()->HNDL_cublas,
	(cublasFillMode_t)uplo,
	(cublasOperation_t)trans,
	(cublasDiagType_t)diag, n,
	AP, x, incx	
    );
}


NLBlas_t nlCUDABlas() {
    static NLboolean initialized = NL_FALSE;
    static struct NLBlas blas;
    if(!initialized) {
	memset(&blas, 0, sizeof(blas));
	blas.has_unified_memory = NL_FALSE;
	blas.Malloc = cuda_blas_malloc;
	blas.Free = cuda_blas_free;
	blas.Memcpy = cuda_blas_memcpy;
	blas.Dcopy = cuda_blas_dcopy;
	blas.Ddot = cuda_blas_ddot;
	blas.Dnrm2 = cuda_blas_dnrm2;
	blas.Daxpy = cuda_blas_daxpy;
	blas.Dscal = cuda_blas_dscal;
	blas.Dgemv = cuda_blas_dgemv;
	blas.Dtpsv = cuda_blas_dtpsv;
	nlBlasResetStats(&blas);
	initialized = NL_TRUE;
    }
    return &blas;
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
	case NL_MASS_MATRIX: {
	    nl_assert(nlCurrentContext->B != NULL);
	    nl_assert(nlCurrentContext->B->type == NL_MATRIX_SPARSE_DYNAMIC);
	    result = (NLSparseMatrix*)(nlCurrentContext->B);
	} break;
	default:
	    nl_assert_not_reached;
    }
    return result;
}




NLboolean nlInitExtension(const char* extension) {
    if(!strcmp(extension, "SUPERLU")) {
        return nlInitExtension_SUPERLU();
    } else if(!strcmp(extension, "CHOLMOD")) {
        return nlInitExtension_CHOLMOD();
    } else if(!strcmp(extension, "ARPACK")) {
	/* 
	 * SUPERLU is needed by OpenNL's ARPACK driver
	 * (factorizes the matrix for the shift-invert spectral
	 *  transform).
	 */
	return nlInitExtension_SUPERLU() && nlInitExtension_ARPACK();
    } else if(!strcmp(extension, "MKL")) {
	return nlInitExtension_MKL();
    } else if(!strcmp(extension, "CUDA")) {
	return nlInitExtension_CUDA();
    }
    return NL_FALSE;
}

NLboolean nlExtensionIsInitialized(const char* extension) {
    if(!strcmp(extension, "SUPERLU")) {
        return nlExtensionIsInitialized_SUPERLU();
    } else if(!strcmp(extension, "CHOLMOD")) {
        return nlExtensionIsInitialized_CHOLMOD();
    } else if(!strcmp(extension, "ARPACK")) {
	/* 
	 * SUPERLU is needed by OpenNL's ARPACK driver
	 * (factorizes the matrix for the shift-invert spectral
	 *  transform).
	 */
	return nlExtensionIsInitialized_SUPERLU() && nlExtensionIsInitialized_ARPACK();
    } else if(!strcmp(extension, "MKL")) {
	return nlExtensionIsInitialized_MKL();
    } else if(!strcmp(extension, "CUDA")) {
	return nlExtensionIsInitialized_CUDA();
    }
    return NL_FALSE;
}

void nlInitialize(int argc, char** argv) {
    int i=0;
    char* ptr=NULL;
    char extension[255];
    /* Find all the arguments with the form:
     * nl:<extension>=true|false
     * and try to activate the corresponding extensions.
     */
    for(i=1; i<argc; ++i) {
	ptr = strstr(argv[i],"=true");
	if(!strncmp(argv[i], "nl:", 3) &&
	   (strlen(argv[i]) > 3) &&
	   (ptr != NULL)) {
	    strncpy(extension, argv[i]+3, (size_t)(ptr-argv[i]-3));
	    extension[(size_t)(ptr-argv[i]-3)] = '\0';
	    if(nlInitExtension(extension)) {
		nl_fprintf(stdout,"OpenNL %s: initialized\n", extension);
	    } else {
		nl_fprintf(stderr,"OpenNL %s: could not initialize\n", extension);		
	    }
	}
    }
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
    case NL_LEAST_SQUARES: {
        nlCurrentContext->least_squares = (NLboolean)param;
    } break;
    case NL_MAX_ITERATIONS: {
        nl_assert(param > 0);
        nlCurrentContext->max_iterations = (NLuint)param;
        nlCurrentContext->max_iterations_defined = NL_TRUE;
    } break;
    case NL_SYMMETRIC: {
        nlCurrentContext->symmetric = (NLboolean)param;        
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

void nlGetBooleanv(NLenum pname, NLboolean* params) {
    switch(pname) {
    case NL_LEAST_SQUARES: {
        *params = nlCurrentContext->least_squares;
    } break;
    case NL_SYMMETRIC: {
        *params = nlCurrentContext->symmetric;
    } break;
    default: {
        nlError("nlGetBooleanv","Invalid parameter");
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
    case NL_LEAST_SQUARES: {
        *params = (NLint)(nlCurrentContext->least_squares);
    } break;
    case NL_MAX_ITERATIONS: {
        *params = (NLint)(nlCurrentContext->max_iterations);
    } break;
    case NL_SYMMETRIC: {
        *params = (NLint)(nlCurrentContext->symmetric);
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


/* Enable / Disable */

void nlEnable(NLenum pname) {
    switch(pname) {
	case NL_NORMALIZE_ROWS: {
	    nl_assert(nlCurrentContext->state != NL_STATE_ROW);
	    nlCurrentContext->normalize_rows = NL_TRUE;
	} break;
	case NL_VERBOSE: {
	    nlCurrentContext->verbose = NL_TRUE;
	} break;
	case NL_VARIABLES_BUFFER: {
	    nlCurrentContext->user_variable_buffers = NL_TRUE;
	} break;
    default: {
        nlError("nlEnable","Invalid parameter");        
        nl_assert_not_reached;
    }
    }
}

void nlDisable(NLenum pname) {
    switch(pname) {
	case NL_NORMALIZE_ROWS: {
	    nl_assert(nlCurrentContext->state != NL_STATE_ROW);
	    nlCurrentContext->normalize_rows = NL_FALSE;
	} break;
	case NL_VERBOSE: {
	    nlCurrentContext->verbose = NL_FALSE;
	} break;
	case NL_VARIABLES_BUFFER: {
	    nlCurrentContext->user_variable_buffers = NL_FALSE;
	} break;
	default: {
	    nlError("nlDisable","Invalid parameter");                
	    nl_assert_not_reached;
	}
    }
}

NLboolean nlIsEnabled(NLenum pname) {
    NLboolean result = NL_FALSE;
    switch(pname) {
	case NL_NORMALIZE_ROWS: {
	    result = nlCurrentContext->normalize_rows;
	} break;
	case NL_VERBOSE: {
	    result = nlCurrentContext->verbose;
	} break;
	case NL_VARIABLES_BUFFER: {
	    result = nlCurrentContext->user_variable_buffers;
	} break;
	default: {
	    nlError("nlIsEnables","Invalid parameter");
	    nl_assert_not_reached;
	}
    }
    return result;
}


/* NL functions */

void  nlSetFunction(NLenum pname, NLfunc param) {
    switch(pname) {
    case NL_FUNC_SOLVER:
        nlCurrentContext->solver_func = (NLSolverFunc)(param);
        nlCurrentContext->solver = NL_SOLVER_USER;	
        break;
    case NL_FUNC_MATRIX:
	nlDeleteMatrix(nlCurrentContext->M);
	nlCurrentContext->M = nlMatrixNewFromFunction(
	    nlCurrentContext->n, nlCurrentContext->n,
	    (NLMatrixFunc)param
	);
        break;
    case NL_FUNC_PRECONDITIONER:
	nlDeleteMatrix(nlCurrentContext->P);
	nlCurrentContext->P = nlMatrixNewFromFunction(
	    nlCurrentContext->n, nlCurrentContext->n,
	    (NLMatrixFunc)param
	);
        nlCurrentContext->preconditioner = NL_PRECOND_USER;
        break;
    case NL_FUNC_PROGRESS:
        nlCurrentContext->progress_func = (NLProgressFunc)(param);
        break;
    default:
        nlError("nlSetFunction","Invalid parameter");        
        nl_assert_not_reached;
    }
}

void nlGetFunction(NLenum pname, NLfunc* param) {
    switch(pname) {
    case NL_FUNC_SOLVER:
        *param = (NLfunc)(nlCurrentContext->solver_func);
        break;
    case NL_FUNC_MATRIX:
        *param = (NLfunc)(nlMatrixGetFunction(nlCurrentContext->M));
        break;
    case NL_FUNC_PRECONDITIONER:
        *param = (NLfunc)(nlMatrixGetFunction(nlCurrentContext->P));
        break;
    default:
        nlError("nlGetFunction","Invalid parameter");                
        nl_assert_not_reached;
    }
}


/* Get/Set Lock/Unlock variables */

void nlSetVariable(NLuint index, NLdouble value) {
    nlCheckState(NL_STATE_SYSTEM);
    nl_debug_range_assert(index, 0, nlCurrentContext->nb_variables - 1);
    NL_BUFFER_ITEM(nlCurrentContext->variable_buffer[0],index) = value;
}

void nlMultiSetVariable(NLuint index, NLuint system, NLdouble value) {
    nlCheckState(NL_STATE_SYSTEM);
    nl_debug_range_assert(index, 0, nlCurrentContext->nb_variables-1);
    nl_debug_range_assert(system, 0, nlCurrentContext->nb_systems-1);    
    NL_BUFFER_ITEM(nlCurrentContext->variable_buffer[system],index) = value;
}

NLdouble nlGetVariable(NLuint index) {
    nl_assert(nlCurrentContext->state != NL_STATE_INITIAL);
    nl_debug_range_assert(index, 0, nlCurrentContext->nb_variables - 1);
    return NL_BUFFER_ITEM(nlCurrentContext->variable_buffer[0],index);
}

NLdouble nlMultiGetVariable(NLuint index, NLuint system) {
    nl_assert(nlCurrentContext->state != NL_STATE_INITIAL);
    nl_debug_range_assert(index, 0, nlCurrentContext->nb_variables-1);
    nl_debug_range_assert(system, 0, nlCurrentContext->nb_systems-1);
    return NL_BUFFER_ITEM(nlCurrentContext->variable_buffer[system],index);    
}


void nlLockVariable(NLuint index) {
    nlCheckState(NL_STATE_SYSTEM);
    nl_debug_range_assert(index, 0, nlCurrentContext->nb_variables - 1);
    nlCurrentContext->variable_is_locked[index] = NL_TRUE;
}

void nlUnlockVariable(NLuint index) {
    nlCheckState(NL_STATE_SYSTEM);
    nl_debug_range_assert(index, 0, nlCurrentContext->nb_variables - 1);
    nlCurrentContext->variable_is_locked[index] = NL_FALSE;
}

NLboolean nlVariableIsLocked(NLuint index) {
    nl_assert(nlCurrentContext->state != NL_STATE_INITIAL);
    nl_debug_range_assert(index, 0, nlCurrentContext->nb_variables - 1);
    return nlCurrentContext->variable_is_locked[index];
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
    
    if(nlCurrentContext->user_variable_buffers) {
	nlCurrentContext->variable_value = NULL;
    } else {
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
    NLenum storage = NL_MATRIX_STORE_ROWS;


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
        if(nlCurrentContext->least_squares || nlCurrentContext->symmetric) {
            nlCurrentContext->solver = NL_CG;
            if(!nlCurrentContext->preconditioner_defined) {
                nlCurrentContext->preconditioner = NL_PRECOND_JACOBI;
            }
        } else {
            nlCurrentContext->solver = NL_BICGSTAB;
        }
        if(!nlCurrentContext->max_iterations_defined) {
            nlCurrentContext->max_iterations = n*5;
        }
        if(!nlCurrentContext->threshold_defined) {
            nlCurrentContext->threshold = 1e-6;
        }
    }

    
    /* SSOR preconditioner requires rows and columns */
    if(nlCurrentContext->preconditioner == NL_PRECOND_SSOR) {
        storage = (storage | NL_MATRIX_STORE_COLUMNS);
    }

    /* a least squares problem results in a symmetric matrix */
    if(nlCurrentContext->least_squares) {
        nlCurrentContext->symmetric = NL_TRUE;
    }

    if(
	nlCurrentContext->symmetric &&
        nlCurrentContext->preconditioner == NL_PRECOND_SSOR 
    ) {
	/* 
	 * For now, only used with SSOR preconditioner, because
	 * for other modes it is either unsupported (SUPERLU) or
	 * causes performance loss (non-parallel sparse SpMV)
	 */
        storage = (storage | NL_MATRIX_STORE_SYMMETRIC);
    }

    nlCurrentContext->M = (NLMatrix)(NL_NEW(NLSparseMatrix));
    nlSparseMatrixConstruct(
	     (NLSparseMatrix*)(nlCurrentContext->M), n, n, storage
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

    nlCurrentContext->right_hand_side = NL_NEW_ARRAY(
	double, nlCurrentContext->nb_systems
    );
    nlCurrentContext->current_row = 0;
}

static void nlEndMatrix() {
    nlTransition(NL_STATE_MATRIX, NL_STATE_MATRIX_CONSTRUCTED);    

    nlRowColumnClear(&nlCurrentContext->af);
    nlRowColumnClear(&nlCurrentContext->al);
    
    if(!nlCurrentContext->least_squares) {
        nl_assert(
            nlCurrentContext->ij_coefficient_called || (
                nlCurrentContext->current_row == 
                nlCurrentContext->n
            )
        );
    }
}

static void nlBeginRow() {
    nlTransition(NL_STATE_MATRIX, NL_STATE_ROW);
    nlRowColumnZero(&nlCurrentContext->af);
    nlRowColumnZero(&nlCurrentContext->al);
}

static void nlScaleRow(NLdouble s) {
    NLRowColumn*    af = &nlCurrentContext->af;
    NLRowColumn*    al = &nlCurrentContext->al;
    NLuint nf            = af->size;
    NLuint nl            = al->size;
    NLuint i,k;
    for(i=0; i<nf; i++) {
        af->coeff[i].value *= s;
    }
    for(i=0; i<nl; i++) {
        al->coeff[i].value *= s;
    }
    for(k=0; k<nlCurrentContext->nb_systems; ++k) {
	nlCurrentContext->right_hand_side[k] *= s;
    }
}

static void nlNormalizeRow(NLdouble weight) {
    NLRowColumn*    af = &nlCurrentContext->af;
    NLRowColumn*    al = &nlCurrentContext->al;
    NLuint nf            = af->size;
    NLuint nl            = al->size;
    NLuint i;
    NLdouble norm = 0.0;
    for(i=0; i<nf; i++) {
        norm += af->coeff[i].value * af->coeff[i].value;
    }
    for(i=0; i<nl; i++) {
        norm += al->coeff[i].value * al->coeff[i].value;
    }
    norm = sqrt(norm);
    nlScaleRow(weight / norm);
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

    if(nlCurrentContext->normalize_rows) {
        nlNormalizeRow(nlCurrentContext->row_scaling);
    } else if(nlCurrentContext->row_scaling != 1.0) {
        nlScaleRow(nlCurrentContext->row_scaling);
    }
    /*
     * if least_squares : we want to solve
     * A'A x = A'b
     */

    if(nlCurrentContext->least_squares) {
        for(i=0; i<nf; i++) {
            for(j=0; j<nf; j++) {
                nlSparseMatrixAdd(
                    M, af->coeff[i].index, af->coeff[j].index,
                    af->coeff[i].value * af->coeff[j].value
                );
            }
        }
	for(k=0; k<nlCurrentContext->nb_systems; ++k) {
	    S = -nlCurrentContext->right_hand_side[k];
	    for(jj=0; jj<nl; ++jj) {
		j = al->coeff[jj].index;
		S += al->coeff[jj].value *
		    NL_BUFFER_ITEM(nlCurrentContext->variable_buffer[k],j);
	    }
	    for(jj=0; jj<nf; jj++) {
		b[ k*n+af->coeff[jj].index ] -= af->coeff[jj].value * S;
	    }
	}
    } else {
        for(jj=0; jj<nf; ++jj) {
            nlSparseMatrixAdd(
                M, current_row, af->coeff[jj].index, af->coeff[jj].value
            );
        }
	for(k=0; k<nlCurrentContext->nb_systems; ++k) {
	    b[k*n+current_row] = nlCurrentContext->right_hand_side[k];
	    for(jj=0; jj<nl; ++jj) {
		j = al->coeff[jj].index;
		b[k*n+current_row] -= al->coeff[jj].value *
		    NL_BUFFER_ITEM(nlCurrentContext->variable_buffer[k],j);
	    }
	}
    }
    nlCurrentContext->current_row++;
    for(k=0; k<nlCurrentContext->nb_systems; ++k) {
	nlCurrentContext->right_hand_side[k] = 0.0;
    }
    nlCurrentContext->row_scaling = 1.0;
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

void nlAddIJCoefficient(NLuint i, NLuint j, NLdouble value) {
    NLSparseMatrix* M  = nlGetCurrentSparseMatrix();    
    nlCheckState(NL_STATE_MATRIX);
    nl_debug_range_assert(i, 0, nlCurrentContext->nb_variables - 1);
    nl_debug_range_assert(j, 0, nlCurrentContext->nb_variables - 1);
#ifdef NL_DEBUG
    for(NLuint i=0; i<nlCurrentContext->nb_variables; ++i) {
        nl_debug_assert(!nlCurrentContext->variable_is_locked[i]);
    }
#endif    
    nlSparseMatrixAdd(M, i, j, value);
    nlCurrentContext->ij_coefficient_called = NL_TRUE;
}

void nlAddIRightHandSide(NLuint i, NLdouble value) {
    nlCheckState(NL_STATE_MATRIX);
    nl_debug_range_assert(i, 0, nlCurrentContext->nb_variables - 1);
#ifdef NL_DEBUG
    for(NLuint i=0; i<nlCurrentContext->nb_variables; ++i) {
        nl_debug_assert(!nlCurrentContext->variable_is_locked[i]);
    }
#endif
    nlCurrentContext->b[i] += value;
    nlCurrentContext->ij_coefficient_called = NL_TRUE;
}

void nlMultiAddIRightHandSide(NLuint i, NLuint k, NLdouble value) {
    NLuint n = nlCurrentContext->n;
    nlCheckState(NL_STATE_MATRIX);
    nl_debug_range_assert(i, 0, nlCurrentContext->nb_variables - 1);
    nl_debug_range_assert(k, 0, nlCurrentContext->nb_systems - 1);
#ifdef NL_DEBUG
    for(NLuint i=0; i<nlCurrentContext->nb_variables; ++i) {
        nl_debug_assert(!nlCurrentContext->variable_is_locked[i]);
    }
#endif
    nlCurrentContext->b[i + k*n] += value;
    nlCurrentContext->ij_coefficient_called = NL_TRUE;
}

void nlRightHandSide(NLdouble value) {
    nlCurrentContext->right_hand_side[0] = value;
}

void nlMultiRightHandSide(NLuint k, NLdouble value) {
    nl_debug_range_assert(k, 0, nlCurrentContext->nb_systems - 1);
    nlCurrentContext->right_hand_side[k] = value;
}

void nlRowScaling(NLdouble value) {
    nlCheckState(NL_STATE_MATRIX);
    nlCurrentContext->row_scaling = value;
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

void nlUpdateRightHandSide(NLdouble* values) {
    /*
     * If we are in the solved state, get back to the
     * constructed state.
     */
    nl_assert(nlCurrentContext->nb_systems == 1);
    if(nlCurrentContext->state == NL_STATE_SOLVED) {
        nlTransition(NL_STATE_SOLVED, NL_STATE_SYSTEM_CONSTRUCTED);
    }
    nlCheckState(NL_STATE_SYSTEM_CONSTRUCTED);
    memcpy(nlCurrentContext->x, values, nlCurrentContext->n * sizeof(double));
}


/* Buffers management */

void nlBindBuffer(
    NLenum buffer, NLuint k, void* addr, NLuint stride
) {
    nlCheckState(NL_STATE_SYSTEM);    
    nl_assert(nlIsEnabled(buffer));
    nl_assert(buffer == NL_VARIABLES_BUFFER);
    nl_assert(k<nlCurrentContext->nb_systems);
    if(stride == 0) {
	stride = sizeof(NLdouble);
    }
    nlCurrentContext->variable_buffer[k].base_address = addr;
    nlCurrentContext->variable_buffer[k].stride = stride;
}


/* Eigen solver */

void nlMatrixMode(NLenum matrix) {
    NLuint n = 0;
    NLuint i;
    nl_assert(
	nlCurrentContext->state == NL_STATE_SYSTEM ||
	nlCurrentContext->state == NL_STATE_MATRIX_CONSTRUCTED
    );
    nlCurrentContext->state = NL_STATE_SYSTEM;
    nlCurrentContext->matrix_mode = matrix;
    nlCurrentContext->current_row = 0;
    nlCurrentContext->ij_coefficient_called = NL_FALSE;
    switch(matrix) {
	case NL_STIFFNESS_MATRIX: {
	    /* Stiffness matrix is already constructed. */
	} break ;
	case NL_MASS_MATRIX: {
	    if(nlCurrentContext->B == NULL) {
		for(i=0; i<nlCurrentContext->nb_variables; ++i) {
		    if(!nlCurrentContext->variable_is_locked[i]) {
			++n;
		    }
		}
		nlCurrentContext->B = (NLMatrix)(NL_NEW(NLSparseMatrix));
		nlSparseMatrixConstruct(
		    (NLSparseMatrix*)(nlCurrentContext->B),
		    n, n, NL_MATRIX_STORE_ROWS
		);
	    }
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
	case NL_SYMMETRIC:
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
	case NL_ARPACK_EXT:
	    nlEigenSolve_ARPACK();
	    break;
	default:
	    nl_assert_not_reached;
    }
}

double nlGetEigenValue(NLuint i) {
    nl_debug_assert(i < nlCurrentContext->nb_variables);
    return nlCurrentContext->eigen_value[i];
}
