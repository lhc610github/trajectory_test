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
        compute_Constraint(int temp_order, int temp_n, int temp_m, int temp_kr, int temp_kpsi);
        ~compute_Constraint();
        bool compute_waypoint_C(float ** temp_keyframe,float *temp_t);
        bool compute_pos_D_C(float temp_eps);
        void print_Matrix(float **sth,int row,int col);
        void print_Vector(float *sth,int length);
        void combine_Cb();
        float ** C;
        float * b;

    private:
        int order;
        int m;
        int n;
        int k_r;
        int k_psi;
        float **keyframe; //new in compute_waypoint_C
        float *t; //new in compute_waypoint_C
        int n_intermediate;
        float corridor_width;
        float ** C1; //new in compute_waypoint_C
        float * b1; //new in compute_waypoint_C
        int C1_size[2]; //[row,col]
        int b1_size[2]; //[row,col]
        float ** C2; //new in compute_waypoint_C
        float * b2; //new in compute_waypoint_C
        int C2_size[2]; //[row,col]
        int b2_size[2]; //[row,col]
        int C_size[2]; //[row,col]
        int b_size[2]; //[row,col]
};
#endif
