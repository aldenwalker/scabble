#ifndef __lp__
#define __lp__


#include <vector>

#include "scabble.h"
#include "rational.h"


using namespace std;

enum scallop_lp_solver {GLPK_DOUBLE, GLPK_EXACT, QSOPT_EXACT, EXLP};

void do_linear_program(   vector<string> &w,
                          vector<double> &weight,
                          vector<arc> &arc_list, 
                          vector<polygon> &polygon_list, 
                          double &scl, 
                          vector<double> &solutionVector,
                          scallop_lp_solver solver, 
                          int VERBOSE);

void linear_program_from_ratmat( vector<polygon>& polygon_list,
                                 vector<rational>& solutionVector,
                                 rational& scl,
                                 RatMat* constraints,
                                 vector<int>& equalityType,
                                 scallop_lp_solver solver, int VERBOSE );

void init_lp();

#endif
