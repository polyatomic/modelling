#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <set>

#include "calc.h"
#include "util.h"

#include "Descriptors.h"

using namespace std;

Descriptors::Descriptors():
m_na(0),
m_na2(0),
m_na3(0),
m_na4(0),
m_nt2(0),
m_nt3(0),
m_nt4(0),
m_ndesc(0),
m_ngp2sum(0),
m_ngp3sum(0),
m_nblocks(0),
m_t2(0),
m_t3(0),
m_t4(0),
m_t2b(0),
m_t3b(0),
m_t4b(0),
m_st3(0),
m_st4(0),
m_ngp2(0),
m_ngp3(0),
m_ngp4(0),
m_svp3i(0),
m_svp4i(0),
m_gamma(0),
m_svp2(0),
m_svp3(0),
m_svp4(0),
m_dist(0)
{
   int i;
   for (i=0; i < 20; i++) m_blocksizes[i] = 0;
}

Descriptors::~Descriptors()
{
   delete [] m_ngp2;
   delete [] m_t2b;
   delete [] m_dist;
   delete [] m_t3b;
   delete [] m_t4b;
   ReleaseGrid();
}

bool Descriptors::Init(int na, int nt2, int *t2, int nt3, int *st3, int *t3, int nt4, int *st4, int *t4, int *ngp, double *gamma, int nthreads)
{
   int i;
   bool res = true;
   ReleaseGrid();
   if (!t2)
   {
      cerr << "\"Descriptors\" object initialization: nothing to do" << endl;
      return false;
   }
   m_t2 = t2;
   cout << "Number of atoms: " << na << endl;
   m_na = na;
   m_na2 = na*(na-1)/2;
   cout << "Number of atom pairs: " << m_na2 << endl;
   cout << "Number of atom type pairs: " << nt2 << endl;
   m_nt2 = nt2;
   delete [] m_ngp2;
   m_ngp2 = new int[nt2];
   for (i=0; i < nt2; i++) m_ngp2[i] = ngp[i];
   m_gamma = gamma;
   delete [] m_t2b;
   m_t2b = new int[nt2];
   delete [] m_dist;
   m_dist = new double[nthreads*m_na2];
   m_nt3 = nt3;
   if (!nt3) return res;
   m_na3 = na*(na-1)*(na-2)/6;
   cout << "Number of atom triplets: " << m_na3 << endl;
   cout << "Number of atom type triplets: " << nt3 << endl;
   m_t3 = t3;
   m_st3 = st3;
   delete [] m_t3b;
   m_t3b = new int[nt3];
   m_nt4 = nt4;
   if (!nt4) return res;
   m_na4 = na*(na-1)*(na-2)*(na-3)/24;
   cout << "Number of atom quadruplets: " << m_na4 << endl;
   cout << "Number of atom type quadruplets: " << nt4 << endl;
   m_t4 = t4;
   m_st4 = st4;
   delete [] m_t4b;
   m_t4b = new int[nt4];
   return res;
}

void Descriptors::GetDistanceRanges(double *r, double *dmin, double *dmax)
{
   int i, j, k, tp;
   double tmp, dist;
   for (j=0,i=0; j < m_na-1; j++)
   {
      for (k=j+1; k < m_na; k++)
      {
         tp = m_t2[i];
         tmp = r[3*j] - r[3*k];
         dist = tmp*tmp;
         tmp = r[3*j+1] - r[3*k+1];
         dist += tmp*tmp;
         tmp = r[3*j+2] - r[3*k+2];
         dist += tmp*tmp;
         dist = sqrt(dist);
         if (dmin[tp] > dist) dmin[tp] = dist;
         if (dmax[tp] < dist) dmax[tp] = dist;
         i++;
      }
   }
}

void Descriptors::FindClosestMultiplets(double *r, set<int> *pss, set<Triplet, TripletCompare> *tss, set<Sixtuplet, SixtupletCompare> *sss)
{
   int i, j, k, i1, i2, i3, i4, tp, itmp;
   int *dloc;
   double d, tmp;
   double *x1, *x2, *dist;
   double d2[6];
   int idx[6];
   int idxsrc[6];
   Triplet t;
   Sixtuplet s;
   if (!m_nt3) return;
   dist = m_dist;
   for (i=0,k=0; i < m_na-1; i++)
   {
      x1 = r + 3*i;
      for (j=i+1; j < m_na; j++)
      {
         tp = m_t2[k];
         x2 = r + 3*j;
         tmp = x1[0] - x2[0];
         d = tmp*tmp;
         tmp = x1[1] - x2[1];
         d += tmp*tmp;
         tmp = x1[2] - x2[2];
         d += tmp*tmp;
         d2[0] = dist[k] = 1.0/sqrt(d);
         GetIdx(d2, idxsrc, tp, 2);
         idx[0] = idxsrc[0];
         pss[tp].insert(idx[0]);
         k++;
      }
   }
   dloc = m_t3 + m_na3;
   for (i1=0,k=0; i1 < m_na-2; i1++)
   {
      for (i2=i1+1; i2 < m_na-1; i2++)
      {
         for (i3=i2+1; i3 < m_na; i3++)
         {
            tp = m_t3[k];
            d2[0] = dist[dloc[3*k]];
            d2[1] = dist[dloc[3*k+1]];
            d2[2] = dist[dloc[3*k+2]];
            GetIdx(d2, idxsrc, tp, 3);
            idx[0] = idxsrc[0];
            idx[1] = idxsrc[1];
            idx[2] = idxsrc[2];
            if (tp == 0 || tp == 3) sort(idx, idx+3);
            else sort(idx, idx+2);
            t.v[0] = idx[0];
            t.v[1] = idx[1];
            t.v[2] = idx[2];
            tss[tp].insert(t);
            k++;
         }
      }
   }
   if (!m_nt4) return;
   dloc = m_t4 + m_na4;
   for (i1=0,k=0; i1 < m_na3-3; i1++)
   {
      for (i2=i1+1; i2 < m_na-2; i2++)
      {
         for (i3=i2+1; i3 < m_na-1; i3++)
         {
            for (i4=i3+1; i4 < m_na; i4++)
            {
               tp = m_t4[k];
               d2[0] = dist[dloc[6*k]];
               d2[1] = dist[dloc[6*k+1]];
               d2[2] = dist[dloc[6*k+2]];
               d2[3] = dist[dloc[6*k+3]];
               d2[4] = dist[dloc[6*k+4]];
               d2[5] = dist[dloc[6*k+5]];
               GetIdx(d2, idxsrc, tp, 4);
               idx[0] = idxsrc[0];
               idx[1] = idxsrc[1];
               idx[2] = idxsrc[2];
               idx[3] = idxsrc[3];
               idx[4] = idxsrc[4];
               idx[5] = idxsrc[5];
               if (tp == 0 || tp == 4)
               {
                  if (idx[0] > idx[5])
                  {
                     itmp = idx[0];
                     idx[0] = idx[5];
                     idx[5] = itmp;
                     itmp = idx[2];
                     idx[2] = idx[3];
                     idx[3] = itmp;
                  }
                  if (idx[1] > idx[4])
                  {
                     itmp = idx[1];
                     idx[1] = idx[4];
                     idx[4] = itmp;
                     itmp = idx[2];
                     idx[2] = idx[3];
                     idx[3] = itmp;
                  }
                  if (idx[0] > idx[1])
                  {
                     itmp = idx[0];
                     idx[0] = idx[1];
                     idx[1] = itmp;
                     itmp = idx[4];
                     idx[4] = idx[5];
                     idx[5] = itmp;
                  }
                  if (idx[0] > idx[2])
                  {
                     itmp = idx[0];
                     idx[0] = idx[2];
                     idx[2] = itmp;
                     itmp = idx[3];
                     idx[3] = idx[5];
                     idx[5] = itmp;
                  }
                  if (idx[1] > idx[2])
                  {
                     itmp = idx[1];
                     idx[1] = idx[2];
                     idx[2] = itmp;
                     itmp = idx[3];
                     idx[3] = idx[4];
                     idx[4] = itmp;
                  }
                  if (idx[0] == idx[1])
                  {
                     if (idx[1] == idx[2])
                     {
                        if (idx[3] > idx[4])
                        {
                           itmp = idx[3];
                           idx[3] = idx[4];
                           idx[4] = itmp;
                        }
                        if (idx[3] > idx[5])
                        {
                           itmp = idx[3];
                           idx[3] = idx[5];
                           idx[5] = itmp;
                        }
                        if (idx[4] > idx[5])
                        {
                           itmp = idx[4];
                           idx[4] = idx[5];
                           idx[5] = itmp;
                        }
                     } else if (idx[4] > idx[5])
                     {
                        itmp = idx[4];
                        idx[4] = idx[5];
                        idx[5] = itmp;
                     }
                  }
               } else if (tp == 1 || tp == 2)
               {
                  if (idx[0] > idx[1])
                  {
                     itmp = idx[0];
                     idx[0] = idx[1];
                     idx[1] = itmp;
                     itmp = idx[4];
                     idx[4] = idx[5];
                     idx[5] = itmp;
                  }
                  if (idx[0] > idx[2])
                  {
                     itmp = idx[0];
                     idx[0] = idx[2];
                     idx[2] = itmp;
                     itmp = idx[3];
                     idx[3] = idx[5];
                     idx[5] = itmp;
                  }
                  if (idx[1] > idx[2])
                  {
                     itmp = idx[1];
                     idx[1] = idx[2];
                     idx[2] = itmp;
                     itmp = idx[3];
                     idx[3] = idx[4];
                     idx[4] = itmp;
                  }
                  if (idx[0] == idx[1])
                  {
                     if (idx[1] == idx[2])
                     {
                        if (idx[3] > idx[4])
                        {
                           itmp = idx[3];
                           idx[3] = idx[4];
                           idx[4] = itmp;
                        }
                        if (idx[3] > idx[5])
                        {
                           itmp = idx[3];
                           idx[3] = idx[5];
                           idx[5] = itmp;
                        }
                        if (idx[4] > idx[5])
                        {
                           itmp = idx[4];
                           idx[4] = idx[5];
                           idx[5] = itmp;
                        }
                     } else if (idx[4] > idx[5])
                     {
                        itmp = idx[4];
                        idx[4] = idx[5];
                        idx[5] = itmp;
                     }
                  }
               } else
               {
                  itmp = (idx[1] < idx[2]) ? idx[1] : idx[2];
                  if (idx[3] < itmp || idx[4] < itmp)
                  {
                     itmp = idx[1];
                     idx[1] = idx[3];
                     idx[3] = itmp;
                     itmp = idx[2];
                     idx[2] = idx[4];
                     idx[4] = itmp;
                  }
                  if (idx[1] > idx[2])
                  {
                     itmp = idx[1];
                     idx[1] = idx[2];
                     idx[2] = itmp;
                     itmp = idx[3];
                     idx[3] = idx[4];
                     idx[4] = itmp;
                  } else if (idx[1] == idx[2])
                  {
                     if (idx[3] > idx[4])
                     {
                        itmp = idx[3];
                        idx[3] = idx[4];
                        idx[4] = itmp;
                     }
                  }
               }
               s.v[0] = idx[0];
               s.v[1] = idx[1];
               s.v[2] = idx[2];
               s.v[3] = idx[3];
               s.v[4] = idx[4];
               s.v[5] = idx[5];
               sss[tp].insert(s);
               k++;
            }
         }
      }
   }
}

void Descriptors::SetGrid(double *dmin, double *dmax)
{
   int i, j, nb, dif;
   double b, e, stp;
   ReleaseGrid();
   m_svp2 = new double*[m_nt2];
   for (i=0; i < m_nt2; i++) m_svp2[i] = 0;
   for (i=0; i < m_nt2; i++)
   {
      b = 1.0/dmax[i];
      e = 1.0/dmin[i];
      stp = (e - b)/(m_ngp2[i] - 1);
      m_svp2[i] = new double[m_ngp2[i]];
      for (j=0; j < m_ngp2[i]; j++)
      {
         m_svp2[i][j] = b + j*stp;
      }
   }
   m_nblocks = 0;
   for (i=0,m_ndesc=0; i < m_nt2; i++)
   {
      m_t2b[i] = m_ndesc;
      m_blocksizes[m_nblocks++] = m_ngp2[i];
      m_ndesc += m_ngp2[i];
   }
   m_ngp2sum = m_ndesc;
   if (!m_nt3) goto end;
   m_svp3 = new double*[m_nt3];
   for (i=0; i < m_nt3; i++) m_svp3[i] = 0;
   m_svp3i = new int*[m_nt3];
   for (i=0; i < m_nt3; i++) m_svp3i[i] = 0;
   m_ngp3 = new int[m_nt3];
   if (!ReadGridPoints(0, 0, 0, 0, "g30.txt")) AssignGridPoints(0, 0, 0, 0);
   if (!ReadGridPoints(1, 1, 2, 1, "g31.txt")) AssignGridPoints(1, 1, 2, 1);
   if (!ReadGridPoints(1, 1, 0, 2, "g32.txt")) AssignGridPoints(1, 1, 0, 2);
   if (!ReadGridPoints(2, 2, 2, 3, "g33.txt")) AssignGridPoints(2, 2, 2, 3);
   for (i=0,j=0; i < m_nt3; i++)
   {
      m_t3b[i] = j;
      m_blocksizes[m_nblocks++] = m_ngp3[i];
      j += m_ngp3[i];
   }
   m_ndesc += j;
   m_ngp3sum = m_ndesc;
   if (!m_nt4) goto end;
   m_svp4 = new double*[m_nt4];
   for (i=0; i < m_nt4; i++) m_svp4[i] = 0;
   m_svp4i = new int*[m_nt4];
   for (i=0; i < m_nt4; i++) m_svp4i[i] = 0;
   m_ngp4 = new int[m_nt4];
   if (!ReadGridPoints(0, 0, 0, 0, 0, 0, 0, "g40.txt") ||
       !ReadGridPoints(1, 1, 1, 2, 2, 2, 1, "g41.txt") ||
       !ReadGridPoints(1, 1, 1, 0, 0, 0, 2, "g42.txt") ||
       !ReadGridPoints(0, 1, 1, 1, 1, 2, 3, "g43.txt") ||
       !ReadGridPoints(2, 2, 2, 2, 2, 2, 4, "g44.txt"))
   {
      m_nt4 = 0;
      delete [] m_ngp4;
      delete [] m_svp4;
      m_ngp4 = 0;
      m_svp4 = 0;
      return;
   }
   for (i=0,j=0; i < m_nt4; i++)
   {
      m_t4b[i] = j;
      m_blocksizes[m_nblocks++] = m_ngp4[i];
      j += m_ngp4[i];
   }
   m_ndesc += j;
end:
   ReadGridPoints(0, "g20.txt");
   ReadGridPoints(1, "g21.txt");
   ReadGridPoints(2, "g22.txt");
   nb = 0;
   for (i=0,j=0; i < m_nt2; i++)
   {
      m_t2b[i] = j;
      m_blocksizes[nb++] = m_ngp2[i];
      j += m_ngp2[i];
   }
   dif = m_ngp2sum - j;
   m_ndesc -= dif;
   m_ngp2sum -= dif;
   m_ngp3sum -= dif;
}

void Descriptors::ReleaseGrid()
{
   int i;
   if (!m_svp2) return;
   for (i=0; i < m_nt2; i++)
   {
      delete [] m_svp2[i];
   }
   delete [] m_svp2;
   m_svp2 = 0;
   if (!m_svp3) return;
   for (i=0; i < m_nt3; i++)
   {
      delete [] m_svp3[i];
      delete [] m_svp3i[i];
   }
   delete [] m_svp3;
   delete [] m_svp3i;
   delete [] m_ngp3;
   m_svp3 = 0;
   m_svp3i = 0;
   m_ngp3 = 0;
   if (!m_svp4) return;
   for (i=0; i < m_nt4; i++)
   {
      delete [] m_svp4[i];
      delete [] m_svp4i[i];
   }
   delete [] m_svp4;
   delete [] m_svp4i;
   delete [] m_ngp4;
   m_svp4 = 0;
   m_svp4i = 0;
   m_ngp4 = 0;
}

void Descriptors::GetIdx(double *d, int *idxs, int tp, int bo)
{
   int i, j;
   if (bo == 2)
   {
      idxs[0] = m_ngp2[tp] - 2;
      for (i=0; i < m_ngp2[tp]; i++)
      {
         if (d[0] < m_svp2[tp][i])
         {
            idxs[0] = i-1;
            break;
         }
      }
      if (idxs[0] < 0) idxs[0] = 0;
      if (fabs(m_svp2[tp][idxs[0]]-d[0]) > fabs(m_svp2[tp][idxs[0]+1]-d[0])) idxs[0]++;
   } else if (bo == 3)
   {
      switch (tp)
      {
         case 0:
            for (j=0; j < 3; j++)
            {
               idxs[j] = m_ngp2[0] - 2;
               for (i=0; i < m_ngp2[0]; i++)
               {
                  if (d[j] < m_svp2[0][i])
                  {
                     idxs[j] = i-1;
                     break;
                  }
               }
               if (idxs[j] < 0) idxs[j] = 0;
               if (fabs(m_svp2[0][idxs[j]]-d[j]) > fabs(m_svp2[0][idxs[j]+1]-d[j])) idxs[j]++;
//               for (i=0,mindif=INF; i < m_ngp2[0]; i++)
//               {
//                  dif = fabs(m_svp2[0][i]-d[j]);
//                  if (dif < mindif)
//                  {
//                     mindif = dif;
//                     idxs[j] = i;
//                  }
//               }
            }
            break;
         case 1:
            for (j=0; j < 2; j++)
            {
               idxs[j] = m_ngp2[1] - 2;
               for (i=0; i < m_ngp2[1]; i++)
               {
                  if (d[j] < m_svp2[1][i])
                  {
                     idxs[j] = i-1;
                     break;
                  }
               }
               if (idxs[j] < 0) idxs[j] = 0;
               if (fabs(m_svp2[1][idxs[j]]-d[j]) > fabs(m_svp2[1][idxs[j]+1]-d[j])) idxs[j]++;
            }
            idxs[2] = m_ngp2[2] - 2;
            for (i=0; i < m_ngp2[2]; i++)
            {
               if (d[2] < m_svp2[2][i])
               {
                  idxs[2] = i-1;
                  break;
               }
            }
            if (idxs[2] < 0) idxs[2] = 0;
            if (fabs(m_svp2[2][idxs[2]]-d[2]) > fabs(m_svp2[2][idxs[2]+1]-d[2])) idxs[2]++;
            break;
         case 2:
            for (j=0; j < 2; j++)
            {
               idxs[j] = m_ngp2[1] - 2;
               for (i=0; i < m_ngp2[1]; i++)
               {
                  if (d[j] < m_svp2[1][i])
                  {
                     idxs[j] = i-1;
                     break;
                  }
               }
               if (idxs[j] < 0) idxs[j] = 0;
               if (fabs(m_svp2[1][idxs[j]]-d[j]) > fabs(m_svp2[1][idxs[j]+1]-d[j])) idxs[j]++;
            }
            idxs[2] = m_ngp2[0] - 2;
            for (i=0; i < m_ngp2[0]; i++)
            {
               if (d[2] < m_svp2[0][i])
               {
                  idxs[2] = i-1;
                  break;
               }
            }
            if (idxs[2] < 0) idxs[2] = 0;
            if (fabs(m_svp2[0][idxs[2]]-d[2]) > fabs(m_svp2[0][idxs[2]+1]-d[2])) idxs[2]++;
            break;
         case 3:
            for (j=0; j < 3; j++)
            {
               idxs[j] = m_ngp2[2] - 2;
               for (i=0; i < m_ngp2[2]; i++)
               {
                  if (d[j] < m_svp2[2][i])
                  {
                     idxs[j] = i-1;
                     break;
                  }
               }
               if (idxs[j] < 0) idxs[j] = 0;
               if (fabs(m_svp2[2][idxs[j]]-d[j]) > fabs(m_svp2[2][idxs[j]+1]-d[j])) idxs[j]++;
            }
            break;
      }
   } else if (bo == 4)
   {
      switch (tp)
      {
         case 0:
            for (j=0; j < 6; j++)
            {
               idxs[j] = m_ngp2[0] - 2;
               for (i=0; i < m_ngp2[0]; i++)
               {
                  if (d[j] < m_svp2[0][i])
                  {
                     idxs[j] = i-1;
                     break;
                  }
               }
               if (idxs[j] < 0) idxs[j] = 0;
               if (fabs(m_svp2[0][idxs[j]]-d[j]) > fabs(m_svp2[0][idxs[j]+1]-d[j])) idxs[j]++;
            }
            break;
         case 1:
            for (j=0; j < 3; j++)
            {
               idxs[j] = m_ngp2[1] - 2;
               for (i=0; i < m_ngp2[1]; i++)
               {
                  if (d[j] < m_svp2[1][i])
                  {
                     idxs[j] = i-1;
                     break;
                  }
               }
               if (idxs[j] < 0) idxs[j] = 0;
               if (fabs(m_svp2[1][idxs[j]]-d[j]) > fabs(m_svp2[1][idxs[j]+1]-d[j])) idxs[j]++;
            }
            for (j=3; j < 6; j++)
            {
               idxs[j] = m_ngp2[2] - 2;
               for (i=0; i < m_ngp2[2]; i++)
               {
                  if (d[j] < m_svp2[2][i])
                  {
                     idxs[j] = i-1;
                     break;
                  }
               }
               if (idxs[j] < 0) idxs[j] = 0;
               if (fabs(m_svp2[2][idxs[j]]-d[j]) > fabs(m_svp2[2][idxs[j]+1]-d[j])) idxs[j]++;
            }
            break;
         case 2:
            for (j=0; j < 3; j++)
            {
               idxs[j] = m_ngp2[1] - 2;
               for (i=0; i < m_ngp2[1]; i++)
               {
                  if (d[j] < m_svp2[1][i])
                  {
                     idxs[j] = i-1;
                     break;
                  }
               }
               if (idxs[j] < 0) idxs[j] = 0;
               if (fabs(m_svp2[1][idxs[j]]-d[j]) > fabs(m_svp2[1][idxs[j]+1]-d[j])) idxs[j]++;
            }
            for (j=3; j < 6; j++)
            {
               idxs[j] = m_ngp2[0] - 2;
               for (i=0; i < m_ngp2[0]; i++)
               {
                  if (d[j] < m_svp2[0][i])
                  {
                     idxs[j] = i-1;
                     break;
                  }
               }
               if (idxs[j] < 0) idxs[j] = 0;
               if (fabs(m_svp2[0][idxs[j]]-d[j]) > fabs(m_svp2[0][idxs[j]+1]-d[j])) idxs[j]++;
            }
            break;
         case 3:
            idxs[0] = m_ngp2[0] - 2;
            for (i=0; i < m_ngp2[0]; i++)
            {
               if (d[0] < m_svp2[0][i])
               {
                  idxs[0] = i-1;
                  break;
               }
            }
            if (idxs[0] < 0) idxs[0] = 0;
            if (fabs(m_svp2[0][idxs[0]]-d[0]) > fabs(m_svp2[0][idxs[0]+1]-d[0])) idxs[0]++;
            for (j=1; j < 5; j++)
            {
               idxs[j] = m_ngp2[1] - 2;
               for (i=0; i < m_ngp2[1]; i++)
               {
                  if (d[j] < m_svp2[1][i])
                  {
                     idxs[j] = i-1;
                     break;
                  }
               }
               if (idxs[j] < 0) idxs[j] = 0;
               if (fabs(m_svp2[1][idxs[j]]-d[j]) > fabs(m_svp2[1][idxs[j]+1]-d[j])) idxs[j]++;
            }
            idxs[5] = m_ngp2[2] - 2;
            for (i=0; i < m_ngp2[2]; i++)
            {
               if (d[5] < m_svp2[2][i])
               {
                  idxs[5] = i-1;
                  break;
               }
            }
            if (idxs[5] < 0) idxs[5] = 0;
            if (fabs(m_svp2[2][idxs[5]]-d[5]) > fabs(m_svp2[2][idxs[5]+1]-d[5])) idxs[5]++;
            break;
         case 4:
            for (j=0; j < 6; j++)
            {
               idxs[j] = m_ngp2[2] - 2;
               for (i=0; i < m_ngp2[2]; i++)
               {
                  if (d[j] < m_svp2[2][i])
                  {
                     idxs[j] = i-1;
                     break;
                  }
               }
               if (idxs[j] < 0) idxs[j] = 0;
               if (fabs(m_svp2[2][idxs[j]]-d[j]) > fabs(m_svp2[2][idxs[j]+1]-d[j])) idxs[j]++;
            }
            break;
      }
   }
}

void Descriptors::Calculate(double *r, bool *bid, double *x, int wsi)
{
   int i, j, k, l, tp, tpb, i1, i2, i3, i4;
   int *dloc;
   double d, tmp;
   double *x1, *x2, *dist, *xp;
   bool *bidp;
   double d1[6];
   double d2[6];
   double gammas[6];
   double (*symmetrized_kernels3[2])(double *, double *, double *, double (*)(double, double, double));
   symmetrized_kernels3[0] = symmetrized_kernel3_0;
   symmetrized_kernels3[1] = symmetrized_kernel3_1;
   double (*symmetrized_kernels4[3])(double *, double *, double *, double (*)(double, double, double));
   symmetrized_kernels4[0] = symmetrized_kernel4_0;
   symmetrized_kernels4[1] = symmetrized_kernel4_1;
   symmetrized_kernels4[2] = symmetrized_kernel4_2;
   dist = m_dist + wsi*m_na2;
   for (i=0; i < m_ndesc; i++) x[i] = 0.0;
   for (i=0,k=0; i < m_na-1; i++)
   {
      x1 = r + 3*i;
      for (j=i+1; j < m_na; j++)
      {
         x2 = r + 3*j;
         tmp = x1[0] - x2[0];
         d = tmp*tmp;
         tmp = x1[1] - x2[1];
         d += tmp*tmp;
         tmp = x1[2] - x2[2];
         d += tmp*tmp;
         dist[k] = 1.0/sqrt(d);
         k++;
      }
   }
   for (i=0; i < m_na2; i++)
   {
      tp = m_t2[i];
      tpb = m_t2b[tp];
      for (l=0; l < m_ngp2[tp]; l++)
      {
         if (bid[tpb+l]) x[tpb+l] += kernel(m_svp2[tp][l], dist[i], m_gamma[tp]);
      }
   }
   if (!m_nt3) goto end;
   xp = x + m_ngp2sum;
   bidp = bid + m_ngp2sum;
   dloc = m_t3 + m_na3;
   for (i1=0,k=0; i1 < m_na-2; i1++)
   {
      for (i2=i1+1; i2 < m_na-1; i2++)
      {
         for (i3=i2+1; i3 < m_na; i3++)
         {
            tp = m_t3[k];
            tpb = m_t3b[tp];
            GetGammas(gammas, tp, 3);
            d2[0] = dist[dloc[3*k]];
            d2[1] = dist[dloc[3*k+1]];
            d2[2] = dist[dloc[3*k+2]];
            for (l=0; l < m_ngp3[tp]; l++)
            {
               if (bidp[tpb+l])
               {
                  d1[0] = m_svp3[tp][3*l];
                  d1[1] = m_svp3[tp][3*l+1];
                  d1[2] = m_svp3[tp][3*l+2];
                  xp[tpb+l] += symmetrized_kernels3[m_st3[tp]](d1, d2, gammas, kernel);
               }
            }
            k++;
         }
      }
   }
   if (!m_nt4) goto end;
   xp = x + m_ngp3sum;
   bidp = bid + m_ngp3sum;
   dloc = m_t4 + m_na4;
   for (i1=0,k=0; i1 < m_na-3; i1++)
   {
      for (i2=i1+1; i2 < m_na-2; i2++)
      {
         for (i3=i2+1; i3 < m_na-1; i3++)
         {
            for (i4=i3+1; i4 < m_na; i4++)
            {
               tp = m_t4[k];
               tpb = m_t4b[tp];
               GetGammas(gammas, tp, 4);
               d2[0] = dist[dloc[6*k]];
               d2[1] = dist[dloc[6*k+1]];
               d2[2] = dist[dloc[6*k+2]];
               d2[3] = dist[dloc[6*k+3]];
               d2[4] = dist[dloc[6*k+4]];
               d2[5] = dist[dloc[6*k+5]];
               for (l=0; l < m_ngp4[tp]; l++)
               {
                  if (bidp[tpb+l])
                  {
                     d1[0] = m_svp4[tp][6*l];
                     d1[1] = m_svp4[tp][6*l+1];
                     d1[2] = m_svp4[tp][6*l+2];
                     d1[3] = m_svp4[tp][6*l+3];
                     d1[4] = m_svp4[tp][6*l+4];
                     d1[5] = m_svp4[tp][6*l+5];
                     xp[tpb+l] += symmetrized_kernels4[m_st4[tp]](d1, d2, gammas, kernel);
                  }
               }
               k++;
            }
         }
      }
   }
end:
   for (i=0,j=0; i < m_ndesc; i++)
   {
      if (bid[i])
      {
         if (i == j)
            j++;
         else
            x[j++] = x[i];
      }
   }
}

int Descriptors::GetNDescriptors()
{
   return m_ndesc;
}

int Descriptors::GetDescriptorBlockSizes(int *blocks)
{
   int i;
   for (i=0; i < m_nblocks; i++) blocks[i] = m_blocksizes[i];
   return m_nblocks;
}

void Descriptors::Get3BGridPoints(set<Triplet, TripletCompare> *tss)
{
   int i, j;
   Triplet t;
   for (i=0; i < m_nt3; i++)
   {
      for (j=0; j < m_ngp3[i]; j++)
      {
         t.v[0] = m_svp3i[i][3*j];
         t.v[1] = m_svp3i[i][3*j+1];
         t.v[2] = m_svp3i[i][3*j+2];
         tss[i].insert(t);
      }
   }
}

void Descriptors::Get4BGridPoints(set<Sixtuplet, SixtupletCompare> *sss)
{
   int i, j;
   Sixtuplet s;
   for (i=0; i < m_nt4; i++)
   {
      for (j=0; j < m_ngp4[i]; j++)
      {
         s.v[0] = m_svp4i[i][6*j];
         s.v[1] = m_svp4i[i][6*j+1];
         s.v[2] = m_svp4i[i][6*j+2];
         s.v[3] = m_svp4i[i][6*j+3];
         s.v[4] = m_svp4i[i][6*j+4];
         s.v[5] = m_svp4i[i][6*j+5];
         sss[i].insert(s);
      }
   }
}

void Descriptors::SetKernel(int kernel_type)
{
   switch (kernel_type)
   {
      case 0:
         kernel = gaussian_kernel;
         cout << "Using squared exponential kernel" << endl;
         break;
      case 1:
         kernel = laplacian_kernel;
         cout << "Using laplacian kernel" << endl;
         break;
      case 2:
         kernel = matern32_kernel;
         cout << "Using Matern 3/2 kernel" << endl;
         break;
      case 3:
         kernel = matern52_kernel;
         cout << "Using Matern 5/2 kernel" << endl;
         break;
   }
}

void Descriptors::GetGammas(double *gammas, int tp, int bo)
{
   int nt22;
   if (bo == 3)
   {
      switch (tp)
      {
         case 0:
            gammas[0] = m_gamma[m_nt2];
            gammas[1] = m_gamma[m_nt2];
            gammas[2] = m_gamma[m_nt2];
            break;
         case 1:
            gammas[0] = m_gamma[m_nt2+1];
            gammas[1] = m_gamma[m_nt2+1];
            gammas[2] = m_gamma[m_nt2+2];
            break;
         case 2:
            gammas[0] = m_gamma[m_nt2+3];
            gammas[1] = m_gamma[m_nt2+3];
            gammas[2] = m_gamma[m_nt2+4];
            break;
         case 3:
            gammas[0] = m_gamma[m_nt2+5];
            gammas[1] = m_gamma[m_nt2+5];
            gammas[2] = m_gamma[m_nt2+5];
            break;
      }
   } else if (bo == 4)
   {
      nt22 = m_nt2+6;
      switch (tp)
      {
         case 0:
            gammas[0] = m_gamma[nt22];
            gammas[1] = m_gamma[nt22];
            gammas[2] = m_gamma[nt22];
            gammas[3] = m_gamma[nt22];
            gammas[4] = m_gamma[nt22];
            gammas[5] = m_gamma[nt22]; 
            break;
         case 1:
            gammas[0] = m_gamma[nt22+1];
            gammas[1] = m_gamma[nt22+1];
            gammas[2] = m_gamma[nt22+1];
            gammas[3] = m_gamma[nt22+2];
            gammas[4] = m_gamma[nt22+2];
            gammas[5] = m_gamma[nt22+2]; 
            break;
         case 2:
            gammas[0] = m_gamma[nt22+3];
            gammas[1] = m_gamma[nt22+3];
            gammas[2] = m_gamma[nt22+3];
            gammas[3] = m_gamma[nt22+4];
            gammas[4] = m_gamma[nt22+4];
            gammas[5] = m_gamma[nt22+4]; 
            break;
         case 3:
            gammas[0] = m_gamma[nt22+5];
            gammas[1] = m_gamma[nt22+6];
            gammas[2] = m_gamma[nt22+6];
            gammas[3] = m_gamma[nt22+6];
            gammas[4] = m_gamma[nt22+6];
            gammas[5] = m_gamma[nt22+7]; 
            break;
         case 4:
            gammas[0] = m_gamma[nt22+8];
            gammas[1] = m_gamma[nt22+8];
            gammas[2] = m_gamma[nt22+8];
            gammas[3] = m_gamma[nt22+8];
            gammas[4] = m_gamma[nt22+8];
            gammas[5] = m_gamma[nt22+8]; 
            break;
      }
   }
}

void Descriptors::AssignGridPoints(int isrc1, int isrc2, int isrc3, int idst)
{
   int i1, i2, i3, count, inxt;
   double d1, d2, d3, d1inv, d2inv, d3inv, hs;
   for (i1=0,count=0; i1 < m_ngp2[isrc1]; i1++)
   {
      d1 = 1.0/m_svp2[isrc1][i1];
      if (isrc1 == isrc2) inxt = i1; else inxt = 0;
      for (i2=inxt; i2 < m_ngp2[isrc2]; i2++)
      {
         d2 = 1.0/m_svp2[isrc2][i2];
         if (isrc2 == isrc3) inxt = i2; else inxt = 0;
         for (i3=inxt; i3 < m_ngp2[isrc3]; i3++)
         {
            d3 = 1.0/m_svp2[isrc3][i3];
            hs = 0.5*(d1 + d2 + d3);
            if ((hs-d1) >= 0.0 && (hs-d2) >= 0.0 && (hs-d3) >= 0.0) count++;
         }
      }
   }
   m_ngp3[idst] = count;
   m_svp3[idst] = new double[3*count];
   m_svp3i[idst] = new int[3*count];
   for (i1=0,count=0; i1 < m_ngp2[isrc1]; i1++)
   {
      d1inv = m_svp2[isrc1][i1];
      d1 = 1.0/d1inv;
      if (isrc1 == isrc2) inxt = i1; else inxt = 0;
      for (i2=inxt; i2 < m_ngp2[isrc2]; i2++)
      {
         d2inv = m_svp2[isrc2][i2];
         d2 = 1.0/d2inv;
         if (isrc2 == isrc3) inxt = i2; else inxt = 0;
         for (i3=inxt; i3 < m_ngp2[isrc3]; i3++)
         {
            d3inv = m_svp2[isrc3][i3];
            d3 = 1.0/d3inv;
            hs = 0.5*(d1 + d2 + d3);
            if ((hs-d1) >= 0.0 && (hs-d2) >= 0.0 && (hs-d3) >= 0.0)
            {
               m_svp3[idst][3*count] = d1inv; m_svp3[idst][3*count+1] = d2inv; m_svp3[idst][3*count+2] = d3inv;
               m_svp3i[idst][3*count] = i1; m_svp3i[idst][3*count+1] = i2; m_svp3i[idst][3*count+2] = i3;
//cout << i1 << " " << i2 << " " << i3 << endl;
               count++;
            }
         }
      }
   }
   cout << "Generated " << count << " grid points for triplets of type " << idst << endl;
}

bool Descriptors::ReadGridPoints(int idst, const char *fn)
{
   int i, ngp, gp;
   ifstream ifile;
   istringstream isstream;
   string line;
   ifile.open(fn);
   if (ifile)
   {
      getline(ifile, line, '\n');
      isstream.str(line);
      isstream >> ngp;
      isstream.clear();
      m_ngp2[idst] = ngp;
      for (i=0; i < ngp; i++)
      {
         getline(ifile, line, '\n');
         isstream.str(line);
         isstream >> gp;
         m_svp2[idst][i] = m_svp2[idst][gp];
         isstream.clear();
      }
      cout << "Read " << ngp << " grid points for pairs of type " << idst << endl;
   } else
   {
      cerr << "Unable to open " << fn << endl;
      return false;
   }
   ifile.close();
   return true;
}

bool Descriptors::ReadGridPoints(int isrc1, int isrc2, int isrc3, int idst, const char *fn)
{
   int i, ngp, gp1, gp2, gp3;
   ifstream ifile;
   istringstream isstream;
   string line;
   ifile.open(fn);
   if (ifile)
   {
      getline(ifile, line, '\n');
      isstream.str(line);
      isstream >> ngp;
      isstream.clear();
      m_ngp3[idst] = ngp;
      if (ngp)
      {
         m_svp3[idst] = new double[3*ngp];
         m_svp3i[idst] = new int[3*ngp];
      }
      for (i=0; i < ngp; i++)
      {
         getline(ifile, line, '\n');
         isstream.str(line);
         isstream >> gp1 >> gp2 >> gp3;
         m_svp3[idst][3*i] = m_svp2[isrc1][gp1];
         m_svp3[idst][3*i+1] = m_svp2[isrc2][gp2];
         m_svp3[idst][3*i+2] = m_svp2[isrc3][gp3];
         m_svp3i[idst][3*i] = gp1;
         m_svp3i[idst][3*i+1] = gp2;
         m_svp3i[idst][3*i+2] = gp3;
         isstream.clear();
      }
      cout << "Read " << ngp << " grid points for triplets of type " << idst << endl;
   } else
   {
      cerr << "Unable to open " << fn << endl;
      return false;
   }
   ifile.close();
   return true;
}

bool Descriptors::ReadGridPoints(int isrc1, int isrc2, int isrc3, int isrc4, int isrc5, int isrc6, int idst, const char *fn)
{
   int i, ngp, gp1, gp2, gp3, gp4, gp5, gp6;
   ifstream ifile;
   istringstream isstream;
   string line;
   ifile.open(fn);
   if (ifile)
   {
      getline(ifile, line, '\n');
      isstream.str(line);
      isstream >> ngp;
      isstream.clear();
      m_ngp4[idst] = ngp;
      if (ngp)
      {
         m_svp4[idst] = new double[6*ngp];
         m_svp4i[idst] = new int[6*ngp];
      }
      for (i=0; i < ngp; i++)
      {
         getline(ifile, line, '\n');
         isstream.str(line);
         isstream >> gp1 >> gp2 >> gp3 >> gp4 >> gp5 >> gp6;
         m_svp4[idst][6*i] = m_svp2[isrc1][gp1];
         m_svp4[idst][6*i+1] = m_svp2[isrc2][gp2];
         m_svp4[idst][6*i+2] = m_svp2[isrc3][gp3];
         m_svp4[idst][6*i+3] = m_svp2[isrc4][gp4];
         m_svp4[idst][6*i+4] = m_svp2[isrc5][gp5];
         m_svp4[idst][6*i+5] = m_svp2[isrc6][gp6];
         m_svp4i[idst][6*i] = gp1;
         m_svp4i[idst][6*i+1] = gp2;
         m_svp4i[idst][6*i+2] = gp3;
         m_svp4i[idst][6*i+3] = gp4;
         m_svp4i[idst][6*i+4] = gp5;
         m_svp4i[idst][6*i+5] = gp6;
         isstream.clear();
      }
      cout << "Read " << ngp << " grid points for quadruplets of type " << idst << endl;
   } else
   {
      cerr << "Unable to open " << fn << endl;
      return false;
   }
   ifile.close();
   return true;
}

double (*Descriptors::kernel)(double, double, double) = 0;
