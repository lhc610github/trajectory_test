#ifndef COMPUTE_CONSTRAINT_H_
#define COMPUTE_CONSTRAINT_H_
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <list>
class compute_Constraint
{
    public:
        compute_Constraint();
        ~compute_Constraint();
    private:
        int order;
        int m;
        int n;
        float k_r;
        float k_psi;
        int n_intermediate;
        float corridor_width;
};
#endif
