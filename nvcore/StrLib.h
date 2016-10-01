// This code is in the public domain -- Ignacio Castaño <castano@gmail.com>

#pragma once
#ifndef NV_CORE_STRING_H
#define NV_CORE_STRING_H

#include "Debug.h"
#include "Hash.h" // hash

//#include <string.h> // strlen, etc.

#if NV_OS_WIN32
#define NV_PATH_SEPARATOR '\\'
#else
#define NV_PATH_SEPARATOR '/'
#endif

namespace nv
{

    NVCORE_API uint strHash(const char * str, uint h) NV_PURE;

    /// String hash based on Bernstein's hash.
    inline uint strHash(const char * data, uint h = 5381)
    {
        uint i = 0;
        while(data[i] != 0) {
            h = (33 * h) ^ uint(data[i]);
            i++;
        }
        return h;
    }

    template <> struct Hash<const char *> {
        uint operator()(const char * str) const { return strHash(str); }
    };

    NVCORE_API uint strLen(const char * str) NV_PURE;                       // Asserts on NULL strings.

    NVCORE_API int strDiff(const char * s1, const char * s2) NV_PURE;       // Asserts on NULL strings.
    NVCORE_API int strCaseDiff(const char * s1, const char * s2) NV_PURE;   // Asserts on NULL strings.
    NVCORE_API bool strEqual(const char * s1, const char * s2) NV_PURE;     // Accepts NULL strings.
    NVCORE_API bool strCaseEqual(const char * s1, const char * s2) NV_PURE; // Accepts NULL strings.

    template <> struct Equal<const char *> {
        bool operator()(const char * a, const char * b) const { return strEqual(a, b); }
    };

    NVCORE_API bool strBeginsWith(const char * dst, const char * prefix) NV_PURE;
    NVCORE_API bool strEndsWith(const char * dst, const char * suffix) NV_PURE;


    NVCORE_API void strCpy(char * dst, uint size, const char * src);
    NVCORE_API void strCpy(char * dst, uint size, const char * src, uint len);
    NVCORE_API void strCat(char * dst, uint size, const char * src);

    NVCORE_API const char * strSkipWhiteSpace(const char * str);
    NVCORE_API char * strSkipWhiteSpace(char * str);

    NVCORE_API bool strMatch(const char * str, const char * pat) NV_PURE;

    NVCORE_API bool isNumber(const char * str) NV_PURE;

    template <int count> void strCpySafe(char (&buffer)[count], const char *src) {
        strCpy(buffer, count, src);
    }

    template <int count> void strCatSafe(char (&buffer)[count], const char * src) {
        strCat(buffer, count, src);
    }
	
    /// String class.
    class NVCORE_CLASS String
    {
    public:

        /// Constructs a null string. @sa isNull()
        String()
        {
            data = NULL;
        }

        /// Constructs a shared copy of str.
        String(const String & str)
        {
            data = str.data;
            if (data != NULL) addRef();
        }

        /// Constructs a shared string from a standard string.
        String(const char * str)
        {
            setString(str);
        }

        /// Constructs a shared string from a standard string.
        String(const char * str, int length)
        {
            setString(str, length);
        }

        /// Dtor.
        ~String()
        {
            release();
        }

        String clone() const;

        /// Release the current string and allocate a new one.
        const String & operator=( const char * str )
        {
            release();
            setString( str );
            return *this;
        }

        /// Implement value semantics.
        String & operator=( const String & str )
        {
            if (str.data != data)
            {
                release();
                data = str.data;
                addRef();
            }
            return *this;
        }

        /// Equal operator.
        bool operator==( const String & str ) const
        {
            return strMatch(str.data, data);
        }

        /// Equal operator.
        bool operator==( const char * str ) const
        {
            return strMatch(str, data);
        }

        /// Not equal operator.
        bool operator!=( const String & str ) const
        {
            return !strMatch(str.data, data);
        }

        /// Not equal operator.
        bool operator!=( const char * str ) const
        {
            return !strMatch(str, data);
        }

        /// Returns true if this string is the null string.
        bool isNull() const { return data == NULL; }

        /// Return the exact length.
        uint length() const { nvDebugCheck(data != NULL); return strLen(data); }

        /// Return the hash of the string.
        uint hash() const { nvDebugCheck(data != NULL); return strHash(data); }

        /// const char * cast operator.
        operator const char * () const { return data; }

        /// Get string pointer.
        const char * str() const { return data; }


    private:

        // Add reference count.
        void addRef();

        // Decrease reference count.
        void release();

        uint16 getRefCount() const
        {
            nvDebugCheck(data != NULL);
            return *reinterpret_cast<const uint16 *>(data - 2);
        }

        void setRefCount(uint16 count) {
            nvDebugCheck(data != NULL);
            nvCheck(count < 0xFFFF);
            *reinterpret_cast<uint16 *>(const_cast<char *>(data - 2)) = uint16(count);
        }

        void setData(const char * str) {
            data = str + 2;
        }

        void allocString(const char * str)
        {
            allocString(str, strLen(str));
        }

        void allocString(const char * str, uint length);

        void setString(const char * str);
        void setString(const char * str, uint length);

        // Swap strings.
        friend void swap(String & a, String & b);

    private:

        const char * data;

    };

    template <> struct Hash<String> {
        uint operator()(const String & str) const { return str.hash(); }
    };

} // nv namespace

#endif // NV_CORE_STRING_H
