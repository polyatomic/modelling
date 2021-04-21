#ifndef _ARCHIVE_H_
#define _ARCHIVE_H_

#include "PersistentStream.h"

typedef void* PtrType;

#ifdef DECLSPEC
   #undef DECLSPEC
#endif

#ifdef _WIN32
   #ifdef MATH_EXPORTS
      #define DECLSPEC __declspec(dllexport)
   #else
      #define DECLSPEC __declspec(dllimport)
   #endif
#else
   #define DECLSPEC
#endif

DECLSPEC void archive(int&, PersistentStream&);
DECLSPEC void archive(bool&, PersistentStream&);
DECLSPEC void archive(double&, PersistentStream&);
DECLSPEC void archive(string&, PersistentStream&);
DECLSPEC void archive(unsigned int&, PersistentStream&);
DECLSPEC void archive(unsigned long&, PersistentStream&);
DECLSPEC void archive(long&, PersistentStream&);
DECLSPEC void archive(PtrType&, PersistentStream&);

#endif
