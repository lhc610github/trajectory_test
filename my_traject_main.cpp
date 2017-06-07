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
   float temp_eps = 0.01;
   compute_A C_A(order,m);
   C_A.run_compute(mu_r, mu_psi,  k_r, k_psi, t_index);
   //C_A.print_A();
   compute_Constraint C_C(order,n,m,k_r-1,k_psi);//don't know the reason why k_r-1
   C_C.compute_waypoint_C(keyframe,t_index);
   C_C.compute_pos_D_C(temp_eps);
   C_C.combine_Cb();
   CGAL::Const_oneset_iterator<CGAL::Comparison_result>
         r(   CGAL::SMALLER);
   bool fl[C_C.C_size[0]];
   float l[C_C.C_size[0]];
   for (int i=0; i<C_C.C_size[0]; i++)
   {
       fl[i] = false;
       l[i] = 0.0;
   }
   bool fu[C_C.C_size[0]];
   float u[C_C.C_size[0]];
   float c[C_C.C_size[0]];
   for (int i=0; i<C_C.C_size[0]; i++)
   {
       fu[i] = false;
       u[i] = 0.0;
       c[i] = 0.0;
   }
   printf("initiall done\n");
   printf("sizeof A is [%d] X [%d]\n",C_C.C_size[0],C_C.C_size[1]);
   printf("sizeof b is [%d] \n",C_C.b_size[1]);
   printf("sizeof D is [%d] X [%d]\n",C_A.A_col_row,C_A.A_col_row);
   //Program qp (C_C.C_size[0],C_C.C_size[1],C_C.C,C_C.b,r,fl,l,fu,u,C_A.c_A_res,c);
   Program qp (C_C.C_size[0],C_C.C_size[1],C_C.C,C_C.b,r,fl,NULL,fu,NULL,C_A.c_A_res,c);
   CGAL::print_linear_program(std::cout, qp, "test_traject");
   Solution s = CGAL::solve_quadratic_program(qp, ET());
   std::cout << s;
   Solution::Index_iterator it = s.basic_variable_indices_begin();
   Solution::Index_iterator end = s.basic_variable_indices_end();
   for (; it != end; ++it)
       std::cout << *it << " ";
   std::cout << std::endl;
   for (int i=0;i <6;i++)
   {
       delete[] keyframe[i];
   }
   delete[] keyframe;
   return 0;
}
