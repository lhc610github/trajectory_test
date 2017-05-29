#include "my_traject_main.h"
int
main(int argc, char **argv)
{
   int order = 6;
   int m = 5;
   float mu_r = 1;
   float mu_psi = 1;
   int k_r = 4;
   int k_psi = 2;
   float t_index[6] = {0,2,4,6,8,10};
   compute_A C_A(order,m);
   C_A.run_compute(mu_r, mu_psi,  k_r, k_psi, t_index);
   C_A.print_A();
   return 0;
}
