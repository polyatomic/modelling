#pragma warning(disable: 4251)

#include "tasks.h"

void divide(int *tsizes, int total, int nthreads)
{
   int i, rem, ts;
   ts = total/nthreads;
   rem = total%nthreads;
   for (i=0; i < rem; i++) tsizes[i] = ts+1;
   for (i=rem; i < nthreads; i++) tsizes[i] = ts;
}
