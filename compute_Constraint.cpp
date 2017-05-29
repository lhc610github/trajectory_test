#include "compute_Constraint.h"
compute_Constraint::
compute_Constraint()
{
}

compute_Constraint::
~compute_Constraint()
{
}
bool
compute_Constraint::
compute_waypoint_C(float ** keyframe,int temp_m)
{
    int i,j,k;
    m = temp_m;//
    float ** C1;
    C1 = new float *[2*m*n];
    for (i=0 ;i < 2*m*n ;i++)
    {
        C1[i] = new float[n*(order+1)*m];
        memset(&(*C1[i]),0,sizeof(float)*(n*(order+1)*m));
    }
    float * b1;
    b1 = new float[2*m*n];
    float *waypoint;
    waypoint new float *[3]

    for (i=0 ;i < m ;i++)
    {
            waypoint = keyframe[];
    }
}
