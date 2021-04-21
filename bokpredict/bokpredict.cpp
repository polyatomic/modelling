#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#ifdef _WIN32
   #include <windows.h>
#endif

#include "RunTasks.h"

PType g_params;

bool ReadParams(const char *fn)
{
   int startpos, len, ival;
   bool bt, res;
   string line, par, val;
   res = true;
   g_params.mbo = 0;
   g_params.kernel_type = 0;
   g_params.nthreads = 1;
   ifstream ifile(fn);
   if (!ifile)
   {
      res = false;
      goto end;
   }
   for (;;)
   {
      if (ifile.eof()) break;
      getline(ifile, line, '\n');
      len = line.size();
      if (len == 0) continue;
      if (line[0] == '#') continue;
      startpos = line.find("=");
      if (startpos != string::npos)
      {
         par = line.substr(0, startpos);
         while (par.size() > 0 && par[par.size()-1] == ' ') par.erase(par.size()-1, 1);
         startpos++;
         val = line.substr(startpos, len-startpos);
         while (val.size() > 0 && val[0] == ' ') val.erase(0, 1);
         bt = false;
         while (val.size() > 0 && val[val.size()-1] == '\\')
         {
            val.erase(val.size()-1, 1);
            bt = true;
         }
         if (bt) continue;
      } else if (line[len-1] == '\\')
      {
         while (line.size() > 0 && line[line.size()-1] == '\\') line.erase(line.size()-1, 1);
         val += '\n';
         val += line;
         continue;
      } else if (par.size() > 0)
      {
         val += '\n';
         val += line;
      }
      if (par == "GEOMETRIES")
      {
         g_params.molFile = val;
      } else if (par == "ENERGIES")
      {
         g_params.eFile = val;
      } else if (par == "MAX_BODY_ORDER")
      {
         ival = atoi(val.c_str());
         g_params.mbo = ival;
      } else if (par == "ATOM_TYPES")
      {
         g_params.t1 = val;
      } else if (par == "ATOM_TYPE_PAIRS")
      {
         g_params.t2 = val;
      } else if (par == "ATOM_TYPE_TRIPLETS")
      {
         g_params.t3 = val;
      } else if (par == "ATOM_TYPE_QUADRUPLETS")
      {
         g_params.t4 = val;
      } else if (par == "ATOM_TYPE_QUINTUPLETS")
      {
         g_params.t5 = val;
      } else if (par == "KERNEL_TYPE")
      {
         ival = atoi(val.c_str());
         g_params.kernel_type = ival;
      } else if (par == "DISTANCES")
      {
         g_params.dist = val;
      } else if (par == "N_THREADS")
      {
         ival = atoi(val.c_str());
         g_params.nthreads = ival;
      } else if (par == "LINEAR_COEFFS")
      {
         g_params.linearCoeffs = val;
      }
      par = "";
   }
end:
   ifile.close();
   return res;
}

int main(int argc, char *argv[])
{
   int i;
   string pfile;
   double sUser=0;
   if (argc < 3)
   {
      cerr << "bokfit version 1.0.0.1\n";
      cerr << "Usage: " << argv[0] << " -p <parameter file>\n";
      return 1;
   }
#ifdef _WIN32
   HANDLE hp;
   FILETIME ftimec, ftimee, fsys, fuser;
   __int64 i64User, i64UserPrev;
   hp = GetCurrentProcess();
   i64UserPrev = 0;
#else
   clock_t c_start = clock();
#endif
   chrono::steady_clock::time_point t_start = chrono::high_resolution_clock::now();
   for (i=1; i < 3; i++)
   {
      if (string(argv[i]) == "-p")
      {
         if (i+1 == argc)
         {
            cerr << "Errors in command line" << endl;
            return 1;
         }
         pfile = argv[i+1];
      }
   }
   const char *pfname = pfile.c_str();
   if (!ReadParams(pfname))
   {
      cerr << "Error reading parameter file\n";
      return 1;
   }
   RunTasks();
   chrono::steady_clock::time_point t_end = chrono::high_resolution_clock::now();
#ifdef _WIN32
   GetProcessTimes(hp, &ftimec, &ftimee, &fsys, &fuser);
   i64User = *((__int64 *)&fuser);
   sUser = ((double)(i64User-i64UserPrev))/10000000;
#else
   clock_t c_end = clock();
   sUser = (c_end - c_start) / CLOCKS_PER_SEC;
#endif
   cout << setiosflags(ios::fixed) << setprecision(1);
   cout << "CPU time: " << sUser << " seconds" << endl;
   cout << "Wall clock time: " << chrono::duration<double>(t_end - t_start).count() << " seconds" << endl;
   return 0;
}
