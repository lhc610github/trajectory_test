#include "my_traject_main.h"
int
main(int argc, char **argv)
{
   int order = 6;
   int m = 5;
   int n = 4;
   float mu_r = 1;
   float mu_psi = 1;
   int k_r = 4;
   int k_psi = 2;
   float t_index[6] = {0,2,4,6,8,10};
   float ** keyframe;
   keyframe = new float *[6];
   for(int i=0;i<6;i++)
   {
       keyframe[i] = new float[4];
       memset(&keyframe[i][0],0,sizeof(float)*4);
   }
   keyframe[3][0] = 1;
   keyframe[3][1] = 2;
   keyframe[3][2] = 3;
   keyframe[3][3] = 4;
   compute_A C_A(order,m);
   C_A.run_compute(mu_r, mu_psi,  k_r, k_psi, t_index);
   //C_A.print_A();
   compute_Constraint C_C(order,n,m,k_r,k_psi);
   C_C.compute_waypoint_C(keyframe,t_index);
   printf("done\n");
   for (int i=0;i <6;i++)
   {
       delete[] keyframe[i];
   }
   delete[] keyframe;
   return 0;
}
