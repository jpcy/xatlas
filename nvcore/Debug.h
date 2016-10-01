// This code is in the public domain -- Ignacio Castaño <castano@gmail.com>

#pragma once
#ifndef NV_CORE_DEBUG_H
#define NV_CORE_DEBUG_H

#include "nvcore.h"
#include <assert.h>
#include <stdarg.h> // va_list

#define nvCheck(exp)     if (!(exp)) { nvDebugPrint("%s %s %s\n", #exp, __FILE__, __LINE__); }
#define nvDebugCheck(exp) assert(exp)
#define nvDebug(...)    nvDebugPrint(__VA_ARGS__)
void nvDebugPrint( const char *msg, ... ) __attribute__((format (printf, 1, 2)));

#endif // NV_CORE_DEBUG_H
