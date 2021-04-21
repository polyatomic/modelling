#pragma warning(disable: 4786 4251)

#include <fstream>
#include <iomanip>

#include "chemistry.h"

using namespace std;

void MCWriteMolBinary(MCMol& mol, ostream& ofile)
{
   int length;
   PersistentStream ps(PersistentStream::save);
   mol.archive(ps);
   length = ps.size();
   ofile.write((char *)&length, sizeof(int));
   ofile.write(ps.GetBuf(), ps.size());
}

void MCWriteXYZFile(const MCMol& mol, ostream& ofile)
{
   int i;
   double *coords;
   ofile << setw(6) << mol.NAtoms() << endl;
   ofile << endl;
   coords = new double[3*mol.NAtoms()];
   mol.GetCoordinates(coords);
   ofile.precision(9);
   ofile.setf(ios::scientific, ios::floatfield);
   for (i=0; i < mol.NAtoms(); i++)
   {
      ofile << mol.GetAtomicSymbol(i) << setw(22) << coords[3*i] << setw(18) << coords[3*i+1] << setw(18) << coords[3*i+2] << endl;
   }
   delete [] coords;
}

void MCWriteGamFile(const MCMol& mol, ostream& ofile, string& kw)
{
   int i;
   double *coords;
   coords = new double[3*mol.NAtoms()];
   mol.GetCoordinates(coords);
   ofile << fixed;
   ofile << kw;
   for (i=0; i < mol.NAtoms(); i++)
   {
      ofile << setprecision(1) << setw(2) << setiosflags(ios::left) << mol.GetAtomicSymbol(i) << resetiosflags(ios::left)
            << setw(6) << (double)mol.GetAtomicNumber(i) << setprecision(9) << setw(18) << coords[i*3]
            << setw(18) << coords[i*3+1] << setw(18) << coords[i*3+2] << endl;
   }
   ofile << " $end" << endl;
   delete [] coords;
}
