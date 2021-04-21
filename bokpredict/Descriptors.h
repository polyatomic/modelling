#ifndef _DESCRIPTORS_H_
#define _DESCRIPTORS_H_

class Descriptors
{
public:
   Descriptors();
   ~Descriptors();
   bool Init(int na, int nt2, int *t2, int nt3, int *st3, int *t3, int nt4, int *st4, int *t4, int *ngp, double *gamma, int nthreads);
   void GetDistanceRanges(double *r, double *dmin, double *dmax);
   void FindClosestMultiplets(double *r, set<int> *pss, set<Triplet, TripletCompare> *tss, set<Sixtuplet, SixtupletCompare> *sss);
   void SetGrid(double *dmin, double *dmax);
   void Calculate(double *r, bool *bid, double *x, int wsi);
   int GetNDescriptors();
   int GetDescriptorBlockSizes(int *blocks);
   void Get3BGridPoints(set<Triplet, TripletCompare> *tss);
   void Get4BGridPoints(set<Sixtuplet, SixtupletCompare> *sss);

   static void SetKernel(int kernel_type);

private:
   void AssignGridPoints(int isrc1, int isrc2, int isrc3, int idst);
   bool ReadGridPoints(int idst, const char *fn);
   bool ReadGridPoints(int isrc1, int isrc2, int isrc3, int idst, const char *fn);
   bool ReadGridPoints(int isrc1, int isrc2, int isrc3, int isrc4, int isrc5, int isrc6, int idst, const char *fn);
   void GetGammas(double *gammas, int tp, int bo);
   void ReleaseGrid();
   void GetIdx(double *d, int *idxs, int tp, int bo);

   int m_na;
   int m_na2;
   int m_na3;
   int m_na4;
   int m_nt2;
   int m_nt3;
   int m_nt4;
   int m_ndesc;
   int m_ngp2sum;
   int m_ngp3sum;
   int m_nblocks;
   int *m_t2;
   int *m_t3;
   int *m_t4;
   int *m_t2b;
   int *m_t3b;
   int *m_t4b;
   int *m_st3;
   int *m_st4;
   int *m_ngp2;
   int *m_ngp3;
   int *m_ngp4;
   int **m_svp3i;
   int **m_svp4i;
   int m_blocksizes[20];
   double *m_gamma;
   double **m_svp2;
   double **m_svp3;
   double **m_svp4;
   double *m_dist;

   static double (*kernel)(double, double, double);
};

#endif
