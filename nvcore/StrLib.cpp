// This code is in the public domain -- Ignacio Castaño <castano@gmail.com>

#include "StrLib.h"

#include "Memory.h"
#include "Utils.h" // swap

#include <math.h>   // log
#include <stdio.h>  // vsnprintf
#include <string.h> // strlen, strcmp, etc.

#if NV_CC_MSVC
#include <stdarg.h> // vsnprintf
#endif

using namespace nv;

namespace 
{
    static char * strAlloc(uint size)
    {
        return malloc<char>(size);
    }

    static char * strReAlloc(char * str, uint size)
    {
        return realloc<char>(str, size);
    }

    static void strFree(const char * str)
    {
        return free<char>(str);
    }

    /*static char * strDup( const char * str )
    {
        nvDebugCheck( str != NULL );
        uint len = uint(strlen( str ) + 1);
        char * dup = strAlloc( len );
        memcpy( dup, str, len );
        return dup;
    }*/

    // helper function for integer to string conversion.
    static char * i2a( uint i, char *a, uint r )
    {
        if( i / r > 0 ) {
            a = i2a( i / r, a, r );
        }
        *a = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"[i % r];
        return a + 1;
    }

    // Locale independent functions.
    static inline char toUpper( char c ) {
        return (c<'a' || c>'z') ? (c) : (c+'A'-'a');
    }
    static inline char toLower( char c ) {
        return (c<'A' || c>'Z') ? (c) : (c+'a'-'A');
    }
    static inline bool isAlpha( char c ) {
        return (c>='a' && c<='z') || (c>='A' && c<='Z');
    }
    static inline bool isDigit( char c ) {
        return c>='0' && c<='9';
    }
    static inline bool isAlnum( char c ) {
        return (c>='a' && c<='z') || (c>='A' && c<='Z') || (c>='0' && c<='9');
    }

}

uint nv::strLen(const char * str)
{
    nvDebugCheck(str != NULL);
    return U32(strlen(str));
}

int nv::strDiff(const char * s1, const char * s2)
{
    nvDebugCheck(s1 != NULL);
    nvDebugCheck(s2 != NULL);
    return strcmp(s1, s2);
}

int nv::strCaseDiff(const char * s1, const char * s2)
{
    nvDebugCheck(s1 != NULL);
    nvDebugCheck(s1 != NULL);
#if NV_CC_MSVC
    return _stricmp(s1, s2);
#else
    return strcasecmp(s1, s2);
#endif
}

bool nv::strEqual(const char * s1, const char * s2)
{
    if (s1 == s2) return true;
    if (s1 == NULL || s2 == NULL) return false;
    return strcmp(s1, s2) == 0;
}

bool nv::strCaseEqual(const char * s1, const char * s2)
{
    if (s1 == s2) return true;
    if (s1 == NULL || s2 == NULL) return false;
    return strCaseDiff(s1, s2) == 0;
}

bool nv::strBeginsWith(const char * str, const char * prefix)
{
    //return strstr(str, prefix) == dst;
    return strncmp(str, prefix, strlen(prefix)) == 0;
}

bool nv::strEndsWith(const char * str, const char * suffix)
{
    uint ml = strLen(str);
    uint sl = strLen(suffix);
    if (ml < sl) return false;
    return strncmp(str + ml - sl, suffix, sl) == 0;
}

// @@ Add asserts to detect overlap between dst and src?
void nv::strCpy(char * dst, uint size, const char * src)
{
    nvDebugCheck(dst != NULL);
    nvDebugCheck(src != NULL);
#if NV_CC_MSVC && _MSC_VER >= 1400
    strcpy_s(dst, size, src);
#else
    NV_UNUSED(size);
    strcpy(dst, src);
#endif
}

void nv::strCpy(char * dst, uint size, const char * src, uint len)
{
    nvDebugCheck(dst != NULL);
    nvDebugCheck(src != NULL);
#if NV_CC_MSVC && _MSC_VER >= 1400
    strncpy_s(dst, size, src, len);
#else
    int n = min(len+1, size);
    strncpy(dst, src, n);
    dst[n-1] = '\0';
#endif
}

void nv::strCat(char * dst, uint size, const char * src)
{
    nvDebugCheck(dst != NULL);
    nvDebugCheck(src != NULL);
#if NV_CC_MSVC && _MSC_VER >= 1400
    strcat_s(dst, size, src);
#else
    NV_UNUSED(size);
    strcat(dst, src);
#endif
}

NVCORE_API const char * nv::strSkipWhiteSpace(const char * str)
{
    nvDebugCheck(str != NULL);
    while (*str == ' ') str++;
    return str;
}

NVCORE_API char * nv::strSkipWhiteSpace(char * str)
{
    nvDebugCheck(str != NULL);
    while (*str == ' ') str++;
    return str;
}


/** Pattern matching routine. I don't remember where did I get this. */
bool nv::strMatch(const char * str, const char * pat)
{
    nvDebugCheck(str != NULL);
    nvDebugCheck(pat != NULL);

    char c2;

    while (true) {
        if (*pat==0) {
            if (*str==0) return true;
            else         return false;
        }
        if ((*str==0) && (*pat!='*')) return false;
        if (*pat=='*') {
            pat++;
            if (*pat==0) return true;
            while (true) {
                if (strMatch(str, pat)) return true;
                if (*str==0) return false;
                str++;
            }
        }
        if (*pat=='?') goto match;
        if (*pat=='[') {
            pat++;
            while (true) {
                if ((*pat==']') || (*pat==0)) return false;
                if (*pat==*str) break;
                if (pat[1] == '-') {
                    c2 = pat[2];
                    if (c2==0) return false;
                    if ((*pat<=*str) && (c2>=*str)) break;
                    if ((*pat>=*str) && (c2<=*str)) break;
                    pat+=2;
                }
                pat++;
            }
            while (*pat!=']') {
                if (*pat==0) {
                    pat--;
                    break;
                }
                pat++;
            }
            goto match;
        }

        if (*pat == NV_PATH_SEPARATOR) {
            pat++;
            if (*pat==0) return false;
        }
        if (*pat!=*str) return false;

match:
        pat++;
        str++;
    }
}

bool nv::isNumber(const char * str) {
    while(*str != '\0') {
        if (!isDigit(*str)) return false;
        str++;
    }
    return true;
}

/// Clone this string
String String::clone() const
{
    String str(data);
    return str;
}

void String::setString(const char * str)
{
    if (str == NULL) {
        data = NULL;
    }
    else {
        allocString( str );
        addRef();
    }
}

void String::setString(const char * str, uint length)
{
    nvDebugCheck(str != NULL);

    allocString(str, length);
    addRef();
}

// Add reference count.
void String::addRef()
{
    if (data != NULL)
    {
        setRefCount(getRefCount() + 1);
    }
}

// Decrease reference count.
void String::release()
{
    if (data != NULL)
    {
        const uint16 count = getRefCount();
        setRefCount(count - 1);
        if (count - 1 == 0) {
            free(data - 2);
            data = NULL;
        }
    }
}

void String::allocString(const char * str, uint len)
{
    const char * ptr = malloc<char>(2 + len + 1);

    setData( ptr );
    setRefCount( 0 );

    // Copy string.
    strCpy(const_cast<char *>(data), len+1, str, len);

    // Add terminating character.
    const_cast<char *>(data)[len] = '\0';
}

void nv::swap(String & a, String & b) {
    swap(a.data, b.data);
}
