#define GEOGRAM_PSM
#ifndef GEO_STATIC_LIBS
#define GEO_DYNAMIC_LIBS
#endif
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

/******* extracted from nl.h *******/

#ifndef OPENNL_H
#define OPENNL_H

#include <stdio.h>

typedef unsigned int    NLenum;

typedef unsigned char   NLboolean;

typedef unsigned int    NLbitfield;

typedef void            NLvoid;

typedef signed char     NLbyte;

typedef short           NLshort;

typedef int             NLint; 

typedef unsigned char   NLubyte;

typedef unsigned short  NLushort;

typedef unsigned int    NLuint;  

typedef unsigned long   NLulong; 

typedef int             NLsizei;

typedef float           NLfloat;

typedef double          NLdouble;

typedef void(*NLfunc)(void);

struct NLContext; 

#define NL_FALSE   0x0
#define NL_TRUE    0x1

#define NL_NB_VARIABLES     0x101

#define NL_MAX_ITERATIONS   0x103

    NLContext *nlNewContext(void);

    void nlDeleteContext(NLContext *context);

    void nlSolverParameteri(NLContext *context, NLenum pname, NLint param);

#define NL_FUNC_PROGRESS       0x603

    void nlSetFunction(NLContext *context, NLenum pname, NLfunc param);

    void nlSetVariable(NLContext *context, NLuint i, NLdouble value);

    NLdouble nlGetVariable(NLContext *context, NLuint i);

    void nlLockVariable(NLContext *context, NLuint index);

#define NL_SYSTEM  0x0

#define NL_MATRIX  0x1

#define NL_ROW     0x2

    void nlBegin(NLContext *context, NLenum primitive);

    void nlEnd(NLContext *context, NLenum primitive);

    void nlCoefficient(NLContext *context, NLuint i, NLdouble value);

    NLboolean nlSolve(NLContext *context);

#endif
