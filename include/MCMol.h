#ifndef _MCMOL_H_
#define _MCMOL_H_

#include "MCAtom.h"
#include "Persistent.h"

#ifdef DECLSPEC
   #undef DECLSPEC
#endif

#ifdef _WIN32
   #ifdef CHEMISTRY_EXPORTS
      #define DECLSPEC __declspec(dllexport)
   #else
      #define DECLSPEC __declspec(dllimport)
   #endif
#else
   #define DECLSPEC
#endif

class DECLSPEC MCMol : public Persistent
{
public:
   MCMol();
   ~MCMol();
   void Init(int, int);
   void SetCoords(const double *);
   void SetCoords(int, const double *);
   void SetAtom(int, const string&);
   void SetComment(const string&);
   string GetComment() const;
   void CreateBond(int, int, int, int);
   void SetFormalCharge(int, int);
   int GetTotalCharge() const;
   unsigned int NAtoms() const;
   unsigned int FindNBonds();
   unsigned int NBonds() const;
   void GetCoordinates(double *) const;
   void GetCoordinates(int, double *) const;
   const char *GetAtomicSymbol(int) const;
   unsigned int GetAtomicNumber(int) const;
   MCAtom *GetAtom(int) const;
   int GetAtomPos(MCAtom *) const;
   MCMol& operator=(const MCMol&);
   virtual void archive(PersistentStream&);

private:
   void Clear();

   int m_size;
   int m_nBonds;
   double *m_AtomCoordinates;
   MCAtom *m_Atoms;
   string m_comment;
};

#endif
