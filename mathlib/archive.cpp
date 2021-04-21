#include "PersistentStream.h"
#include "archive.h"

void archive(int& i, PersistentStream& ps)
{
   if (ps.IsSaving())
      ps.write(reinterpret_cast<char *>(&i), sizeof(i));
   else
      ps.read(reinterpret_cast<char *>(&i), sizeof(i));
}

void archive(double& d, PersistentStream& ps)
{
   if (ps.IsSaving())
      ps.write(reinterpret_cast<char *>(&d), sizeof(d));
   else
      ps.read(reinterpret_cast<char *>(&d), sizeof(d));
}

void archive(string& s, PersistentStream& ps)
{
   int i;
   if (ps.IsSaving())
   {
      i = s.size();
      ps.write(reinterpret_cast<char *>(&i), sizeof(i));
      ps.write(s.c_str(), i);
   }
   else
   {
      ps.read(reinterpret_cast<char *>(&i), sizeof(i));
      char *p = new char[i];
      ps.read(p, i);
      s.assign(p, i);
      delete [] p;
   }
}

void archive(unsigned int& i, PersistentStream& ps)
{
   if (ps.IsSaving())
      ps.write(reinterpret_cast<char *>(&i), sizeof(i));
   else
      ps.read(reinterpret_cast<char *>(&i), sizeof(i));
}

void archive(bool& b, PersistentStream& ps)
{
   if (ps.IsSaving())
      ps.write(reinterpret_cast<char *>(&b), sizeof(b));
   else
      ps.read(reinterpret_cast<char *>(&b), sizeof(b));
}

void archive(unsigned long& i, PersistentStream& ps)
{
   if (ps.IsSaving())
      ps.write(reinterpret_cast<char *>(&i), sizeof(i));
   else
      ps.read(reinterpret_cast<char *>(&i), sizeof(i));
}

void archive(long& i, PersistentStream& ps)
{
   if (ps.IsSaving())
      ps.write(reinterpret_cast<char *>(&i), sizeof(i));
   else
      ps.read(reinterpret_cast<char *>(&i), sizeof(i));
}

void archive(PtrType& p, PersistentStream& ps)
{
   if (ps.IsSaving())
      ps.write(reinterpret_cast<char *>(&p), sizeof(p));
   else
      ps.read(reinterpret_cast<char *>(&p), sizeof(p));
}
