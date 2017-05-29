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
        compute_Constraint(int temp_order, int temp_n, int temp_m);
        ~compute_Constraint();
        bool compute_waypoint_C(float ** temp_keyframe,float *temp_t);

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
};
#endif
