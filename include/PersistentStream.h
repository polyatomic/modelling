#ifndef _PERSISTENTSTREAM_H_
#define _PERSISTENTSTREAM_H_

#include <strstream>

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

using namespace std;

class DECLSPEC PersistentStream
{
public:
   PersistentStream(bool);
   ~PersistentStream();
   bool IsSaving() const;
   void write(const char *pch, int nCount);
   void read(char *pch, int nCount);
   void clear();
   int size() const;
   const char *GetBuf();

   static const bool save;
   static const bool load;

private:
   bool m_mode;
   bool m_frozen;
   strstream *m_stream;
};

#endif
