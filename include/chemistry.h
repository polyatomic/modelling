#ifndef _CHEMISTRY_H_
#define _CHEMISTRY_H_

#include "MCIterator.h"
#include "MCMol.h"

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

DECLSPEC MCIterator<MCMol> MCReadXYZFile(const char *fileName);
DECLSPEC void MCReadXYZFile(MCIterator<MCMol>& molecules, istream& InputStream);
DECLSPEC void MCWriteMolBinary(MCMol& mol, ostream& ofile);
DECLSPEC void MCGetNextMolBinary(MCMol& molecule, istream& InputStream);
DECLSPEC void MCWriteXYZFile(const MCMol& mol, ostream& ofile);
DECLSPEC void MCWriteGamFile(const MCMol& mol, ostream& ofile, string& kw);

#endif
