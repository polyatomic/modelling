#ifndef _MCITERATOR_H_
#define _MCITERATOR_H_

#pragma warning(disable: 4231)

#include <string>

#ifdef DECLSPEC
   #undef DECLSPEC
#endif

#ifdef EXPIMP_TEMPLATE
   #undef EXPIMP_TEMPLATE
#endif

#ifdef _WIN32
   #ifdef MATH_EXPORTS
      #define DECLSPEC __declspec(dllexport)
      #define EXPIMP_TEMPLATE
   #else
      #define DECLSPEC __declspec(dllimport)
      #define EXPIMP_TEMPLATE extern
   #endif
#else
   #define DECLSPEC
   #ifdef MATH_EXPORTS
      #define EXPIMP_TEMPLATE
   #else
      #define EXPIMP_TEMPLATE extern
   #endif
#endif

template<class T>
class MCIterator
{
public:
   MCIterator(unsigned long s = 0):
   m_size(s),
   m_counter(1),
   m_inside(true)
   {
      if (s == 0)
      {
         m_inside = false;
         m_type = m_curpos = 0;
         return;
      }
      m_type = new T[s];
      m_curpos = m_type;
   }

   MCIterator(const MCIterator &rhs)
   {
      m_size = rhs.m_size;
      m_counter = rhs.m_counter;
      m_inside = rhs.m_inside;
      if (m_size == 0)
      {
         m_type = m_curpos = 0;
         return;
      }
      m_type = new T[m_size];
      for (unsigned long i=0; i < m_size; i++)
         m_type[i] = rhs.m_type[i];
      m_curpos = m_type + (m_counter-1);
   }

   MCIterator& operator=(const MCIterator &rhs)
   {
      if (this == &rhs)
         return *this;
      delete [] m_type;
      m_size = rhs.m_size;
      m_counter = rhs.m_counter;
      m_inside = rhs.m_inside;
      if (m_size == 0)
      {
         m_type = m_curpos = 0;
         return *this;
      }
      m_type = new T[m_size];
      for (unsigned long i=0; i < m_size; i++)
         m_type[i] = rhs.m_type[i];
      m_curpos = m_type + (m_counter-1);
      return *this;
   }

   ~MCIterator()
   {
      delete [] m_type;
   }

   const MCIterator& operator++()
   {
      if (m_counter >= m_size)
         m_inside = false;
      else
      {
         ++m_curpos;
         ++m_counter;
      }
      return *this;
   }

   operator bool() const
   {
      return m_inside;
   }

   T& operator *()
   {
      return *m_curpos;
   }

   const T& operator *() const
   {
      return *m_curpos;
   }

   T& operator [](unsigned long offSet)
   {
      return m_curpos[offSet];
   }

   const T& operator [](unsigned long offSet) const
   {
      return m_curpos[offSet];
   }

   MCIterator& ToFirst()
   {
      if (m_size == 0)
         return *this;
      m_counter = 1;
      m_curpos = m_type;
      m_inside = true;
      return *this;
   }

   unsigned long Size() const
   {
      return m_size;
   }

   void Resize(unsigned long s)
   {
      unsigned long min, i;
      m_counter = 1;
      m_inside = true;
      if (s == 0)
      {
         delete [] m_type;
         m_size = s;
         m_inside = false;
         m_type = m_curpos = 0;
         return;
      }
      T *type = new T[s];
      if (s < m_size) min = s;
      else min = m_size;
      m_size = s;
      for (i=0; i < min; i++)
      {
         type[i] = m_type[i];
      }
      delete [] m_type;
      m_type = type;
      m_curpos = m_type;
   }

private:
   unsigned long m_counter;
   T *m_curpos;
   bool m_inside;

protected:
   T *m_type;
   unsigned long m_size;
};

EXPIMP_TEMPLATE template class DECLSPEC MCIterator<std::string>;
EXPIMP_TEMPLATE template class DECLSPEC MCIterator<void *>;

#endif
