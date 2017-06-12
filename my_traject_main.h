#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <list>
#include <cassert>
#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
#include "compute_A.h"
#include "compute_Constraint.h"
// choose exact integral type
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz ET;
#else
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
#endif

 //program and solution types
typedef CGAL::Quadratic_program_from_iterators
<float**,                                             // for A
 float*,                                              // for b
 //CGAL::Const_oneset_iterator<CGAL::Comparison_result>,// for r
 CGAL::Comparison_result*,                            // for r
 bool*,                                               // for fl
 float*,                                              // for l
 bool*,                                               // for fu
 float*,                                              // for u
 float**,                                             // for D
 float*>                                              // for c
Program;

//typedef CGAL::Nonnegative_linear_program_from_iterators
//<float**,                                             // for A
 //float*,                                              // for b
 //CGAL::Const_oneset_iterator<CGAL::Comparison_result>,// for r
 //float*>                                              // for c
//Program;

typedef CGAL::Quadratic_program_solution<ET> Solution;

int main(int argc, char **argv);
