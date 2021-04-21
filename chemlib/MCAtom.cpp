#pragma warning(disable: 4786 4251)

#include "archive.h"
#include "MCAtom.h"

MCAtom::MCAtom():
m_elemno(0),
m_FormalCharge(0)
{
}

MCAtom::~MCAtom()
{
}

void MCAtom::SetElement(const string& element)
{
   map<string, unsigned int>::iterator it;
   it = m_elemnos.find(element);
   if (it != m_elemnos.end())
   {
      m_elemno = it->second;
   }
}

unsigned int MCAtom::GetAtomicNumber() const
{
   return m_elemno;
}

int MCAtom::GetCoordinationNumber() const
{
   return m_sigmaBonds.size();
}

int MCAtom::NumberOfNeighbors(int n) const
{
   int i, ret=0;
   for (i=0; i < m_sigmaBonds.size(); i++)
   {
      if (m_sigmaBonds[i]->GetAtomicNumber() == n) ret++;
   }
   return ret;
}

MCAtom *MCAtom::GetNeighbor(int i) const
{
   if (i >= 0 && i < m_sigmaBonds.size())
      return m_sigmaBonds[i];
   return NULL;
}

int MCAtom::GetBondType(int i) const
{
   if (i >= 0 && i < m_sigmaBonds.size())
      return m_bondTypes[i];
   return 0;
}

int MCAtom::GetBondStereo(int i) const
{
   if (i >= 0 && i < m_sigmaBonds.size())
      return m_bondStereos[i];
   return 0;
}

void MCAtom::SetBondStereo(int i, int st)
{
   if (i >= 0 && i < m_sigmaBonds.size())
      m_bondStereos[i] = st;
}

int MCAtom::GetFormalCharge() const
{
   return m_FormalCharge;
}

bool MCAtom::Contains(MCAtom *a) const
{
   for (int i=0; i < m_sigmaBonds.size(); i++)
   {
      MCAtom *atm = m_sigmaBonds[i];
      if (a == atm) return true;
   }
   return false;
}

void MCAtom::archive(PersistentStream& ps)
{
   int cn, i, tmp;
   ::archive(m_elemno, ps);
   ::archive(m_FormalCharge, ps);
   if (ps.IsSaving())
   {
      cn = m_sigmaBonds.size();
      ::archive(cn, ps);
   } else
   {
      ::archive(cn, ps);
   }
   for (i=0; i < cn; i++)
   {
      if (ps.IsSaving())
      {
         tmp = m_bondTypes[i];
         ::archive(tmp, ps);
         tmp = m_bondStereos[i];
         ::archive(tmp, ps);
      } else
      {
         ::archive(tmp, ps);
         m_bondTypes.push_back(tmp);
         ::archive(tmp, ps);
         m_bondStereos.push_back(tmp);
      }
   }
}

const char *MCAtom::GetAtomicSymbol(unsigned int atomicNumber)
{
   map<unsigned int, string>::iterator it;
   it = m_elemsymbols.find(atomicNumber);
   if (it != m_elemsymbols.end())
   {
      return (it->second).c_str();
   }
   return m_dummyName.c_str();
}

unsigned int MCAtom::GetAtomicNumber(const char *symbol)
{
   unsigned int elemno = 0;
   map<string, unsigned int>::iterator it;
   it = m_elemnos.find(symbol);
   if (it != m_elemnos.end())
   {
      elemno = it->second;
   }
   return elemno;
}

double MCAtom::GetAveragedWeight(unsigned int atomicNumber)
{
   map<unsigned int, double>::iterator it;
   it = m_averagedweights.find(atomicNumber);
   if (it != m_averagedweights.end())
   {
      return it->second;
   }
   return 0.0;
}

map<string, unsigned int> MCAtom::m_elemnos = GetElemNos();
map<unsigned int, string> MCAtom::m_elemsymbols = GetElemSymbols();
map<unsigned int, double> MCAtom::m_averagedweights = GetAveragedWeights();
const string MCAtom::m_dummyName = "";
