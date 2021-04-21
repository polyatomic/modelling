#include <string>
#include <iostream>
#include <fstream>

#include "util.h"

bool File2Array(string& fn, double arr[])
{
   int i, len;
   double fval;
   string line;
   ifstream efile(fn.c_str());
   i = 0;
   if (efile)
   {
      for (;;)
      {
         if (efile.eof()) break;
         getline(efile, line, '\n');
         len = line.size();
         if (len == 0) continue;
         fval = atof(line.c_str());
         arr[i++] = fval;
      }
      efile.close();
   } else
   {
      cerr << "Error opening " << fn << endl;
      return false;
   }
   return true;
}

vector<string> MCTokenize(const string& str, const string& delimiters)
{
   vector<string> tokens;
   string::size_type lastPos = str.find_first_not_of(delimiters, 0);
   string::size_type pos = str.find_first_of(delimiters, lastPos);
   while (string::npos != pos || string::npos != lastPos)
   {
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      lastPos = str.find_first_not_of(delimiters, pos);
      pos = str.find_first_of(delimiters, lastPos);
   }
   return tokens;
}
