#ifndef COMPUTE_A_H_
#define COMPUTE_A_H_
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <list>
class compute_A
{
    public:
        compute_A(int temp_order, int temp_m);
        ~compute_A();
        void run_compute(float mu_r, float mu_psi, int k_r, int k_psi, float *t);
        void print_A();
        void print_vector(int *sth, int length);
        void print_matrix(float **sth, int row, int col);
    private:
        int order;
        int m;
        int A_col_row;
        float ** c_A_res;
        void init(int temp_order, int temp_m);
};
#endif
