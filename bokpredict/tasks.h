#ifndef _TASKS_H_
#define _TASKS_H_

class DIIDFunctor
{
public:
   virtual void operator()(double *darg1, bool *barg, double *darg2, int iarg) = 0;
   virtual void Call(double *darg1, bool *barg, double *darg2, int iarg) = 0;
};

template<class FClass>
class TDIIDFunctor : public DIIDFunctor
{
public:
   TDIIDFunctor()
   {
      pt2Object = 0;
      fpt = 0;
   }

   TDIIDFunctor(FClass *_pt2Object, void(FClass::*_fpt)(double *darg1, bool *barg, double *darg2, int iarg))
   {
      pt2Object = _pt2Object;
      fpt = _fpt;
   }

   virtual void init(FClass *_pt2Object, void(FClass::*_fpt)(double *darg1, bool *barg, double *darg2, int iarg))
   {
      pt2Object = _pt2Object;
      fpt = _fpt;
   }

   virtual void operator()(double *darg1, bool *barg, double *darg2, int iarg)
   {
      (*pt2Object.*fpt)(darg1, barg, darg2, iarg);
   }

   virtual void Call(double *darg1, bool *barg, double *darg2, int iarg)
   {
      (*pt2Object.*fpt)(darg1, barg, darg2, iarg);
   }

private:
   void (FClass::*fpt)(double *darg1, bool *barg, double *darg2, int iarg);
   FClass *pt2Object;
};

#define MAX_THREADS 16

void divide(int *tsizes, int total, int nthreads);

#endif
