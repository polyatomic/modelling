#pragma warning(disable: 4786 4251)

#include <string>
#include <fstream>
#include <sstream>

#include "chemistry.h"

MCIterator<MCMol> MCReadXYZFile(const char *fileName)
{
   MCIterator<MCMol> ret(0);
   ifstream ifile(fileName);
   if (!ifile)
   {
      goto end;
   }
   MCReadXYZFile(ret, ifile);
end:
   ifile.close();
   return ret;
}

void MCReadXYZFile(MCIterator<MCMol>& molecules, istream& InputStream)
{
   int i, j, na, nmols;
   string line, ele;
   double coords[3];
   istringstream isstream;
   getline(InputStream, line, '\n');
   na = atoi(line.c_str());
   for (nmols=0;;)
   {
      if (InputStream.eof()) break;
      getline(InputStream, line, '\n');
      nmols++;
   }
   nmols /= (na+2);
   InputStream.clear();
   InputStream.seekg(0, ios::beg);
   molecules.Resize(nmols);
   for (j=0; j < nmols; j++)
   {
      getline(InputStream, line, '\n');
      getline(InputStream, line, '\n');
      molecules[j].Init(na, 0);
      for (i=0; i < na; i++)
      {
         getline(InputStream, line, '\n');
         isstream.str(line);
         isstream >> ele >> coords[0] >> coords[1] >> coords[2];
         isstream.clear();
         molecules[j].SetCoords(i, coords);
         molecules[j].SetAtom(i, ele);
      }
   }
}

void MCGetNextMolBinary(MCMol& molecule, istream& InputStream)
{
   int length;
   char *buf;
   PersistentStream ps(PersistentStream::load);
   InputStream.read((char *)&length, sizeof(int));
   buf = new char[length];
   InputStream.read(buf, length);
   ps.write((const char *)buf, length);
   delete [] buf;
   molecule.archive(ps);
   ps.clear();
}
