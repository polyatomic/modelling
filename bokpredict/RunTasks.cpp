#pragma warning(disable: 4251)

#include <iostream>
#include <iomanip>
#include <set>
#include <algorithm>
#include <thread>

#include "chemistry.h"
#include "util.h"
#include "calc.h"

#include "RunTasks.h"
#include "Descriptors.h"
#include "tasks.h"

extern PType g_params;

void analyze_triplet(Triplet& t, int na, int i, int j, int k, int pids[])
{
   int itmp;
   if (t.v[0] == t.v[1])
   {
      if (t.v[1] == t.v[2])
      {
         pids[0] = IDX(na,i,j);
         pids[1] = IDX(na,i,k);
         pids[2] = IDX(na,j,k);
      } else
      {
         itmp = t.v[0];
         t.v[0] = t.v[2];
         t.v[2] = itmp;
         pids[0] = IDX(na,i,k);
         pids[1] = IDX(na,j,k);
         pids[2] = IDX(na,i,j);
      }
   } else
   {
      if (t.v[1] == t.v[2])
      {
         pids[0] = IDX(na,i,j);
         pids[1] = IDX(na,i,k);
         pids[2] = IDX(na,j,k);
      } else
      {
         itmp = t.v[0];
         t.v[0] = t.v[1];
         t.v[1] = itmp;
         pids[0] = IDX(na,i,j);
         pids[1] = IDX(na,j,k);
         pids[2] = IDX(na,i,k);
      }
   }
}

void analyze_quadruplet(Quadruplet& q, int na, int i, int j, int k, int l, int pids[])
{
   int itmp;
   if (q.v[0] == q.v[1])
   {
      if (q.v[1] == q.v[2])
      {
         if (q.v[2] == q.v[3])
         {
            pids[0] = IDX(na,i,j);
            pids[1] = IDX(na,i,k);
            pids[2] = IDX(na,i,l);
            pids[3] = IDX(na,j,k);
            pids[4] = IDX(na,j,l);
            pids[5] = IDX(na,k,l);
         } else
         {
            itmp = q.v[0];
            q.v[0] = q.v[3];
            q.v[3] = itmp;
            pids[0] = IDX(na,i,l);
            pids[1] = IDX(na,j,l);
            pids[2] = IDX(na,k,l);
            pids[3] = IDX(na,i,j);
            pids[4] = IDX(na,i,k);
            pids[5] = IDX(na,j,k);
         }
      } else
      {
         if (q.v[2] == q.v[3])
         {
            if (q.v[1] > q.v[2])
            {
               itmp = q.v[0];
               q.v[0] = q.v[2];
               q.v[2] = itmp;
               itmp = q.v[1];
               q.v[1] = q.v[3];
               q.v[3] = itmp;
               pids[0] = IDX(na,k,l);
               pids[1] = IDX(na,i,k);
               pids[2] = IDX(na,j,k);
               pids[3] = IDX(na,i,l);
               pids[4] = IDX(na,j,l);
               pids[5] = IDX(na,i,j);
            } else
            {
               pids[0] = IDX(na,i,j);
               pids[1] = IDX(na,i,k);
               pids[2] = IDX(na,i,l);
               pids[3] = IDX(na,j,k);
               pids[4] = IDX(na,j,l);
               pids[5] = IDX(na,k,l);
            }
         } else
         {
            itmp = q.v[0];
            q.v[0] = q.v[2];
            q.v[2] = itmp;
            pids[0] = IDX(na,j,k);
            pids[1] = IDX(na,i,k);
            pids[2] = IDX(na,k,l);
            pids[3] = IDX(na,i,j);
            pids[4] = IDX(na,j,l);
            pids[5] = IDX(na,i,l);
         }
      }
   } else
   {
      if (q.v[1] == q.v[2])
      {
         if (q.v[2] == q.v[3])
         {
            pids[0] = IDX(na,i,j);
            pids[1] = IDX(na,i,k);
            pids[2] = IDX(na,i,l);
            pids[3] = IDX(na,j,k);
            pids[4] = IDX(na,j,l);
            pids[5] = IDX(na,k,l);
         } else
         {
            if (q.v[0] < q.v[1])
            {
               itmp = q.v[1];
               q.v[1] = q.v[3];
               q.v[3] = itmp;
               pids[0] = IDX(na,i,l);
               pids[1] = IDX(na,i,j);
               pids[2] = IDX(na,i,k);
               pids[3] = IDX(na,j,l);
               pids[4] = IDX(na,k,l);
               pids[5] = IDX(na,j,k);
            } else
            {
               itmp = q.v[0];
               q.v[0] = q.v[2];
               q.v[2] = itmp;
               pids[0] = IDX(na,j,k);
               pids[1] = IDX(na,i,k);
               pids[2] = IDX(na,k,l);
               pids[3] = IDX(na,i,j);
               pids[4] = IDX(na,j,l);
               pids[5] = IDX(na,i,l);
            }
         }
      } else
      {
         if (q.v[2] == q.v[3])
         {
            itmp = q.v[1];
            q.v[1] = q.v[0];
            q.v[0] = itmp;
            pids[0] = IDX(na,i,j);
            pids[1] = IDX(na,j,k);
            pids[2] = IDX(na,j,l);
            pids[3] = IDX(na,i,k);
            pids[4] = IDX(na,i,l);
            pids[5] = IDX(na,k,l);
         } else
         {
            if (q.v[1] < q.v[2])
            {
               itmp = q.v[0];
               q.v[0] = q.v[3];
               q.v[3] = itmp;
               pids[0] = IDX(na,j,l);
               pids[1] = IDX(na,k,l);
               pids[2] = IDX(na,i,l);
               pids[3] = IDX(na,j,k);
               pids[4] = IDX(na,i,j);
               pids[5] = IDX(na,i,k);
            } else
            {
               itmp = q.v[1];
               q.v[1] = q.v[2];
               q.v[2] = itmp;
               pids[0] = IDX(na,i,k);
               pids[1] = IDX(na,i,j);
               pids[2] = IDX(na,i,l);
               pids[3] = IDX(na,j,k);
               pids[4] = IDX(na,k,l);
               pids[5] = IDX(na,j,l);
            }
         }
      }
   }
}

int GetSubType(Triplet& t)
{
   if (t.v[0] == t.v[1])
   {
      if (t.v[1] == t.v[2])
      {
         return 0;
      }
   }
   return 1;
}

int GetSubType(Quadruplet& q)
{
   int i, comp, count;
   comp = q.v[0];
   for (i=1,count=1; i < 4; i++)
   {
      if (q.v[i] == comp) count++;
   }
   if (count == 2) return 2;
   if (count == 4) return 0;
   return 1;
}

void run_DIIDFunctor2(int ivars[], bool *bb, double *dvars[], DIIDFunctor *dcalc)
{
   int j, k;
   double *rcoords;
   int istart = ivars[0];
   int iend = ivars[1];
   int olength = ivars[2];
   int nidi = ivars[3];
   int wsi = ivars[4];
   int ii = ivars[5];
   int r = ivars[6];
   double *rp = dvars[0];
   double *wp = dvars[1];
   double *B = dvars[2];
   double *Ypred = dvars[3];
   for (j=istart; j < iend; j++)
   {
      rcoords = rp + j*olength;
      (*dcalc)(rcoords, bb, wp, wsi);
      for (k=0; k < nidi; k++) Ypred[ii] += B[k]*wp[k];
      Ypred[ii] += B[nidi];
      ii++;
   }
}

void RunTasks()
{
   int i, j, k, l, m, o, mbo, nstoich, nstr, na, n, sz1, sz2, nstrtotal, nt1, nt2, nt3, nt4, lt, tp, ndesc, na3, na4, ii, jj, kk, nb, lc, npar, jstart, jend, tc;
   int *nstrs, *ns, *nas, *itypes, *ngp, *stypes, *idinit, *tsizes;
   int pids[10];
   int bsizes[20];
   char stlabels[20][10];
   int **ltypes1, **ltypes2, **ltypes3, **ltypes4, **stypes3, **stypes4;
   double dval, sum, tmp, maesum;
   double *Y, *yptr, *v, *vptr, *dmin, *dmax, *rcoords, *gamma, *B, *Ypred, *wrow;
   bool *bb;
   vector<int> t1;
   vector<int>::iterator t1pos;
   vector<Pair> t2;
   vector<Triplet> t3;
   vector<Quadruplet> t4;
   vector<Pair>::iterator t2pos;
   vector<Triplet>::iterator t3pos;
   vector<Quadruplet>::iterator t4pos;
   MCIterator<MCMol> mols, molstmp, molsTest;
   set<Triplet, TripletCompare> tss[4];
   set<Triplet, TripletCompare> tss_ref[4];
   set<Triplet, TripletCompare>::iterator it, it2;
   set<Sixtuplet, SixtupletCompare> sss[5];
   set<Sixtuplet, SixtupletCompare> sss_ref[5];
   set<Sixtuplet, SixtupletCompare>::iterator its, its2;
   set<int> pss[3];
   set<int>::iterator itp;
   Pair p;
   Triplet t;
   Quadruplet q;
   int prec = numeric_limits<double>::max_digits10;
   streamsize oldprec = cout.precision();
   B = 0;
   bb = 0;
   idinit = 0;
   wrow = 0;
   mbo = g_params.mbo;
   if (mbo < 1)
   {
      cerr << "Undefined maximum body order" << endl;
      return;
   } else if (mbo > 5)
   {
      cerr << "Body order exceeding 5 not supported" << endl;
      return;
   }
   string s = g_params.molFile;
   vector<string> tkns2;
   vector<string> tkns = MCTokenize(s, "\r\n");
   nstoich = tkns.size();
   if (nstoich > 0)
      cout << "Reading " << nstoich << " stoichiometrie(s)" << endl;
   else
   {
      cerr << "Configurations undefined" << endl;
      return;
   }
   nstrs = new int[nstoich];
   ns = new int[nstoich];
   nas = new int[nstoich];
   for (i=0; i < nstoich; i++)
   {
      s = tkns[i];
      molstmp = MCReadXYZFile(s.c_str());
      if (!molstmp)
      {
         cerr << "Error opening " << s << endl;
         delete [] nstrs;
         delete [] ns;
         delete [] nas;
         return;
      }
      nstr = molstmp.Size();
      nstrs[i] = nstr;
      cout << "Number of configurations = " << nstr << endl;
      na = molstmp[0].NAtoms();
      cout << "Number of atoms in configuration = " << na << endl << endl;
      n = 3*na;
      ns[i] = n;
      nas[i] = na;
      sz1 = mols.Size();
      sz2 = mols.Size() + nstr;
      mols.Resize(sz2);
      for (k=0,j=sz1; j < sz2; j++,k++)
         mols[j] = molstmp[k];
   }
   molstmp.Resize(0);
   nstrtotal = mols.Size();
   s = g_params.eFile;
   tkns = MCTokenize(s, "\r\n");
   if (tkns.size() > 0 && tkns.size() != nstoich)
   {
      cerr << "Input error" << endl;
      delete [] nstrs;
      delete [] ns;
      delete [] nas;
      return;
   }
   Y = 0;
   if (tkns.size() > 0)
   {
      Y = new double[nstrtotal];
      for (i=0,sz1=0; i < nstoich; i++)
      {
         s = tkns[i];
         if (s.size() > 0)
         {
            yptr = Y + sz1;
            if (!File2Array(s, yptr))
            {
               cerr << "Error reading " << s << endl;
               delete [] nstrs;
               delete [] ns;
               delete [] nas;
               delete [] Y;
               return;
            }
         } else
         {
            cerr << "Input error" << endl;
            delete [] nstrs;
            delete [] ns;
            delete [] nas;
            delete [] Y;
            return;
         }
         sz1 += nstrs[i];
      }
   }
   tsizes = new int[g_params.nthreads];
   int ivars[MAX_THREADS][6];
   double *dvars[MAX_THREADS][4];
   thread threads[MAX_THREADS];
   tkns = MCTokenize(g_params.t1, "\r\n");
   nt1 = tkns.size();
   cout << "Number of atom types = " << nt1 << endl;
   cout << "Atom types:" << endl;
   for (i=0; i < nt1; i++)
   {
      j = atoi(tkns[i].c_str());
      t1.push_back(j);
      cout << MCAtom::GetAtomicSymbol(j) << endl;
   }
   cout << endl;
   tkns = MCTokenize(g_params.t2, "\r\n");
   nt2 = tkns.size();
   lc = 0;
   cout << "Number of atom type pairs = " << nt2 << endl;
   cout << "Atom type pairs, numbers of grid points for distance and shape factors:" << endl;
   ngp = new int[nt2];
   switch (mbo)
   {
      case 2:
         npar = nt2;
         break;
      case 3:
         npar = nt2 + 6;
         break;
      case 4:
         npar = nt2 + 15;
         break;
   }
   gamma = new double[npar];
   for (i=0; i < nt2; i++)
   {
      ngp[i] = 0;
   }
   for (i=0; i < npar; i++)
   {
      gamma[i] = 1.0;
   }
   for (i=0; i < nt2; i++)
   {
      tkns2 = MCTokenize(tkns[i], " ");
      j = atoi(tkns2[0].c_str());
      p.v[0] = j;
      j = atoi(tkns2[1].c_str());
      p.v[1] = j;
      t2.push_back(p);
      if (tkns2.size() > 2)
      {
         j = atoi(tkns2[2].c_str());
         ngp[i] = j;
      }
      if (tkns2.size() > 3)
      {
         dval = atof(tkns2[3].c_str());
         gamma[i] = dval;
      }
      strcpy(stlabels[lc], MCAtom::GetAtomicSymbol(p.v[0]));
      strcat(stlabels[lc], MCAtom::GetAtomicSymbol(p.v[1]));
      lc++;
      cout << MCAtom::GetAtomicSymbol(p.v[0]) << " " << MCAtom::GetAtomicSymbol(p.v[1]) << " " << ngp[i] << " " << gamma[i];
      cout << endl;
   }
   k = nt2;
   cout << endl;
   tkns = MCTokenize(g_params.t3, "\r\n");
   nt3 = tkns.size();
   cout << "Number of atom type triplets = " << nt3 << endl;
   cout << "Atom type triplets and shape factors:" << endl;
   for (i=0; i < nt3; i++)
   {
      tkns2 = MCTokenize(tkns[i], " ");
      j = atoi(tkns2[0].c_str());
      t.v[0] = j;
      j = atoi(tkns2[1].c_str());
      t.v[1] = j;
      j = atoi(tkns2[2].c_str());
      t.v[2] = j;
      t3.push_back(t);
      if (tkns2.size() > 3)
      {
         dval = atof(tkns2[3].c_str());
         gamma[k++] = dval;
         if (tkns2.size() > 4)
         {
            dval = atof(tkns2[4].c_str());
            gamma[k++] = dval;
         }
      }
      strcpy(stlabels[lc], MCAtom::GetAtomicSymbol(t.v[0]));
      strcat(stlabels[lc], MCAtom::GetAtomicSymbol(t.v[1]));
      strcat(stlabels[lc], MCAtom::GetAtomicSymbol(t.v[2]));
      lc++;
      cout << MCAtom::GetAtomicSymbol(t.v[0]) << " " << MCAtom::GetAtomicSymbol(t.v[1]) << " " << MCAtom::GetAtomicSymbol(t.v[2]);
      if (tkns2.size() == 4)
      {
         cout << " " << gamma[k-1];
      } else if (tkns2.size() == 5)
      {
         cout << " " << gamma[k-2] << " " << gamma[k-1];
      }
      cout << endl;
   }
   cout << endl;
   tkns = MCTokenize(g_params.t4, "\r\n");
   nt4 = tkns.size();
   cout << "Number of atom type quadruplets = " << nt4 << endl;
   cout << "Atom type quadruplets and shape factors:" << endl;
   for (i=0; i < nt4; i++)
   {
      tkns2 = MCTokenize(tkns[i], " ");
      j = atoi(tkns2[0].c_str());
      q.v[0] = j;
      j = atoi(tkns2[1].c_str());
      q.v[1] = j;
      j = atoi(tkns2[2].c_str());
      q.v[2] = j;
      j = atoi(tkns2[3].c_str());
      q.v[3] = j;
      t4.push_back(q);
      if (tkns2.size() > 4)
      {
         dval = atof(tkns2[4].c_str());
         gamma[k++] = dval;
         if (tkns2.size() > 5)
         {
            dval = atof(tkns2[5].c_str());
            gamma[k++] = dval;
            if (tkns2.size() > 6)
            {
               dval = atof(tkns2[6].c_str());
               gamma[k++] = dval;
            }
         }
      }
      strcpy(stlabels[lc], MCAtom::GetAtomicSymbol(q.v[0]));
      strcat(stlabels[lc], MCAtom::GetAtomicSymbol(q.v[1]));
      strcat(stlabels[lc], MCAtom::GetAtomicSymbol(q.v[2]));
      strcat(stlabels[lc], MCAtom::GetAtomicSymbol(q.v[3]));
      lc++;
      cout << MCAtom::GetAtomicSymbol(q.v[0]) << " " << MCAtom::GetAtomicSymbol(q.v[1]) << " " << MCAtom::GetAtomicSymbol(q.v[2]) << " " << MCAtom::GetAtomicSymbol(q.v[3]);
      if (tkns2.size() == 5)
      {
         cout << " " << gamma[k-1];
      } else if (tkns2.size() == 6)
      {
         cout << " " << gamma[k-2] << " " << gamma[k-1];
      } else if (tkns2.size() == 7)
      {
         cout << " " << gamma[k-3] << " " << gamma[k-2] << " " << gamma[k-1];
      }
      cout << endl;
   }
   cout << endl;
   for (i=0,j=0; i < nstoich; i++) j += nstrs[i]*ns[i];
   v = new double[j];
   for (i=0,sz1=0,sz2=0; i < nstoich; i++)
   {
      vptr = v + sz1;
      for (j=0; j < nstrs[i]; j++) mols[sz2+j].GetCoordinates(vptr+j*ns[i]);
      if (i > 0) mols[i] = mols[sz2];
      sz1 += nstrs[i]*ns[i];
      sz2 += nstrs[i];
   }
   mols.Resize(nstoich);
   Descriptors *descs = new Descriptors[nstoich];
   ltypes1 = new int*[nstoich];
   ltypes2 = new int*[nstoich];
   for (i=0; i < nstoich; i++) ltypes2[i] = 0;
   ltypes3 = new int*[nstoich];
   for (i=0; i < nstoich; i++) ltypes3[i] = 0;
   ltypes4 = new int*[nstoich];
   for (i=0; i < nstoich; i++) ltypes4[i] = 0;
   stypes3 = new int*[nstoich];
   for (i=0; i < nstoich; i++) stypes3[i] = 0;
   stypes4 = new int*[nstoich];
   for (i=0; i < nstoich; i++) stypes4[i] = 0;
   dmin = new double[nt2];
   dmax = new double[nt2];
   tkns = MCTokenize(g_params.dist, "\r\n");
   if (tkns.size() == nt2)
   {
      for (i=0; i < nt2; i++)
      {
         tkns2 = MCTokenize(tkns[i], " ");
         dmin[i] = atof(tkns2[0].c_str());
         dmax[i] = atof(tkns2[1].c_str());
      }
   }
   TDIIDFunctor<Descriptors> *dcalcs = new TDIIDFunctor<Descriptors>[nstoich];
   DIIDFunctor **dcalcsp = new DIIDFunctor*[nstoich];
   for (lt=0,sz1=0; lt < nstoich; lt++)
   {
      cout << "Stoichiometry #" << lt+1 << endl;
      na = nas[lt];
      itypes = ltypes1[lt] = new int[na];
      for (i=0; i < na; i++)
      {
         k = mols[lt].GetAtomicNumber(i);
         t1pos = lower_bound(t1.begin(), t1.end(), k);
         tp = t1pos - t1.begin();
         itypes[i] = tp;
      }
      if (mbo > 1)
      {
         itypes = ltypes2[lt] = new int[na*(na-1)/2];
         for (i=0,m=0; i < na-1; i++)
         {
            k = mols[lt].GetAtomicNumber(i);
            for (j=i+1; j < na; j++)
            {
               p.v[0] = k;
               p.v[1] = mols[lt].GetAtomicNumber(j);
               if (p.v[0] > p.v[1]) { l = p.v[0]; p.v[0] = p.v[1]; p.v[1] = l; }
               t2pos = lower_bound(t2.begin(), t2.end(), p, PairCompare());
               tp = t2pos - t2.begin();
               itypes[m++] = tp;
            }
         }
      }
      if (mbo > 2)
      {
         if (nt3)
         {
            stypes = stypes3[lt] = new int[nt3];
            for (i=0; i < nt3; i++) stypes[i] = GetSubType(t3[i]);
         }
         na3 = na*(na-1)*(na-2)/6;
         if (na3) itypes = ltypes3[lt] = new int[4*na3];
         for (i=0,o=0; i < na-2; i++)
         {
            l = mols[lt].GetAtomicNumber(i);
            for (j=i+1; j < na-1; j++)
            {
               m = mols[lt].GetAtomicNumber(j);
               for (k=j+1; k < na; k++)
               {
                  t.v[0] = l;
                  t.v[1] = m;
                  t.v[2] = mols[lt].GetAtomicNumber(k);
                  analyze_triplet(t, na, i, j, k, pids);
                  t3pos = find(t3.begin(), t3.end(), t);
                  tp = t3pos - t3.begin();
                  itypes[o] = tp;
                  itypes[na3+3*o] = pids[0]; itypes[na3+3*o+1] = pids[1]; itypes[na3+3*o+2] = pids[2];
                  o++;
               }
            }
         }
      }
      if (mbo > 3)
      {
         if (nt4)
         {
            stypes = stypes4[lt] = new int[nt4];
            for (i=0; i < nt4; i++) stypes[i] = GetSubType(t4[i]);
         }
         na4 = na*(na-1)*(na-2)*(na-3)/24;
         if (na4) itypes = ltypes4[lt] = new int[7*na4];
         for (i=0,o=0; i < na-3; i++)
         {
            ii = mols[lt].GetAtomicNumber(i);
            for (j=i+1; j < na-2; j++)
            {
               jj = mols[lt].GetAtomicNumber(j);
               for (k=j+1; k < na-1; k++)
               {
                  kk = mols[lt].GetAtomicNumber(k);
                  for (l=k+1; l < na; l++)
                  {
                     q.v[0] = ii;
                     q.v[1] = jj;
                     q.v[2] = kk;
                     q.v[3] = mols[lt].GetAtomicNumber(l);
                     analyze_quadruplet(q, na, i, j, k, l, pids);
                     t4pos = find(t4.begin(), t4.end(), q);
                     tp = t4pos - t4.begin();
                     itypes[o] = tp;
                     itypes[na4+6*o] = pids[0]; itypes[na4+6*o+1] = pids[1]; itypes[na4+6*o+2] = pids[2]; itypes[na4+6*o+3] = pids[3]; itypes[na4+6*o+4] = pids[4]; itypes[na4+6*o+5] = pids[5];
                     o++;
                  }
               }
            }
         }
      }
      if (!descs[lt].Init(na, nt2, ltypes2[lt], nt3, stypes3[lt], ltypes3[lt], nt4, stypes4[lt], ltypes4[lt], ngp, gamma, g_params.nthreads)) goto end;
      if (g_params.dist.size() == 0)
      {
         for (i=0; i < nt2; i++)
         {
            dmin[i] = 1000.0;
            dmax[i] = 0.0;
         }
         vptr = v + sz1;
         for (i=0; i < nstrs[lt]; i++)
         {
            rcoords = vptr + i*ns[lt];
            descs[lt].GetDistanceRanges(rcoords, dmin, dmax);
         }
      }
      cout << "Minimum and maximum distances and their differences:" << endl;
      for (i=0; i < nt2; i++)
         cout << MCAtom::GetAtomicSymbol(t2[i].v[0]) << " " << MCAtom::GetAtomicSymbol(t2[i].v[1]) << " " << setprecision(prec) << dmin[i] << " " << dmax[i] << " " << setprecision(oldprec) << dmax[i]-dmin[i] << endl;
      descs[lt].SetGrid(dmin, dmax);
      sz1 += nstrs[lt]*ns[lt];
      cout << endl;
   }
   ndesc = descs[0].GetNDescriptors();
   cout << "Number of descriptors: " << ndesc << endl << endl;
   nb = descs[0].GetDescriptorBlockSizes(bsizes);
   cout << "Descriptor block sizes for different subsystems:" << endl;
   for (i=0; i < nb; i++)
   {
      cout << stlabels[i] << " " << bsizes[i] << endl;
   }
   cout << endl;
   for (i=0; i < nstoich; i++)
   {
      dcalcs[i].init(&descs[i], &Descriptors::Calculate);
      dcalcsp[i] = &dcalcs[i];
   }
   Descriptors::SetKernel(g_params.kernel_type);
   bb = new bool[ndesc];
   for (i=0; i < ndesc; i++) bb[i] = true;
   wrow = new double[ndesc*g_params.nthreads];
   B = new double[ndesc+1];
   s = g_params.linearCoeffs;
   if (s.size() > 0)
   {
      cout << "Reading descriptor coeffs" << endl;
      File2Array(s, B);
   } else
   {
      cerr << "Descriptor coeffs undefined" << endl;
      goto end;
   }
   if (g_params.nthreads > 1) cout << "Using " << g_params.nthreads << " threads" << endl;
   Ypred = new double[nstrtotal];
   for (j=0; j < nstrtotal; j++) Ypred[j] = 0.0;
   for (i=0,sz1=0,ii=0; i < nstoich; i++)
   {
      vptr = v + sz1;
      divide(tsizes, nstrs[i], g_params.nthreads);
      jstart = 0;
      jend = tsizes[0];
      for (tc=0; tc < g_params.nthreads; tc++)
      {
         ivars[tc][0] = jstart;
         ivars[tc][1] = jend;
         ivars[tc][2] = ns[i];
         ivars[tc][3] = ndesc;
         ivars[tc][4] = tc;
         ivars[tc][5] = ii;
         dvars[tc][0] = vptr;
         dvars[tc][1] = wrow + tc*ndesc;
         dvars[tc][2] = B;
         dvars[tc][3] = Ypred;
         threads[tc] = thread(run_DIIDFunctor2, ivars[tc], bb, dvars[tc], dcalcsp[i]);
         ii += tsizes[tc];
         jstart = jend;
         if (tc < g_params.nthreads-1) jend = jstart + tsizes[tc+1];
      }
      for (tc=0; tc < g_params.nthreads; tc++)
         threads[tc].join();
      sz1 += nstrs[i]*ns[i];
   }
   cout << endl << "Observed\tPredicted\tError" << endl;
   for (j=0,sum=0.0,maesum=0.0; j < nstrtotal; j++)
   {
      tmp = Ypred[j] - Y[j];
      cout << Y[j] << "\t" << Ypred[j] << "\t" << setprecision(prec) << tmp << setprecision(oldprec) << endl;
      sum += tmp*tmp;
      maesum += fabs(tmp);
   }
   cout << "RMSE(test) = " << sqrt(sum/nstrtotal) << endl;
   cout << "MAE(test) = " << maesum/nstrtotal << endl;
   delete [] Ypred;
end:
   delete [] nstrs;
   delete [] ns;
   delete [] nas;
   delete [] Y;
   delete [] tsizes;
   delete [] ngp;
   delete [] gamma;
   delete [] v;
   delete [] descs;
   for (i=0; i < nstoich; i++) delete [] ltypes1[i];
   delete [] ltypes1;
   for (i=0; i < nstoich; i++) delete [] ltypes2[i];
   delete [] ltypes2;
   for (i=0; i < nstoich; i++) delete [] ltypes3[i];
   delete [] ltypes3;
   for (i=0; i < nstoich; i++) delete [] ltypes4[i];
   delete [] ltypes4;
   for (i=0; i < nstoich; i++) delete [] stypes3[i];
   delete [] stypes3;
   for (i=0; i < nstoich; i++) delete [] stypes4[i];
   delete [] stypes4;
   delete [] dmin;
   delete [] dmax;
   delete [] dcalcs;
   delete [] dcalcsp;
   delete [] idinit;
   delete [] bb;
   delete [] wrow;
   delete [] B;
}
