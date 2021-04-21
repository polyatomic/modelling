#ifndef _RUNTASKS_H_
#define _RUNTASKS_H_

#include <string>

using namespace std;

struct PType
{
   string molFile;
   string eFile;
   string t1;
   string t2;
   string t3;
   string t4;
   string t5;
   string dist;
   string linearCoeffs;
   int mbo;
   int kernel_type;
   int nthreads;
};

void RunTasks();

#endif
