#ifndef _PERSISTENT_H_
#define _PERSISTENT_H_

#include "PersistentStream.h"

class Persistent
{
public:
   virtual void archive(PersistentStream&) = 0;
};

#endif
