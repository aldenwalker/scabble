#ifndef __scallop__
#define __scallop__

#include <vector>

struct arc{
	int first;
	int last;
	
	int first_word;
	int last_word;
};

struct polygon {
    int size;
	  std::vector<int> arc;
};


extern "C" {
  //#include "sssgmp.h"
}
#include <gmp.h>
//#include "gmp/gmp-exec/include/gmp.h"

extern "C" {
  #include "matrix.h"
}

#endif
