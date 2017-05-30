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

    private:
        int order;
        int m;
        int n;
        float k_r;
        float k_psi;
        float **keyframe; //new in compute_waypoint_C
        float *t; //new in compute_waypoint_C
        int n_intermediate;
        float corridor_width;
        float ** C1; //new in compute_waypoint_C
        float * b1; //new in compute_waypoint_C
        float ** C2; //new in compute_waypoint_C
        float * b2; //new in compute_waypoint_C
};
#endif
