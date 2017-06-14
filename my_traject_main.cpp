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
   keyframe[0][0] = 0;
   keyframe[0][1] = 0;
   keyframe[0][2] = 0;
   keyframe[0][3] = 0;
   keyframe[1][0] = 2;
   keyframe[1][1] = 3;
   keyframe[1][2] = 5;
   keyframe[1][3] = 1;
   keyframe[2][0] = 1;
   keyframe[2][1] = 2;
   keyframe[2][2] = 3;
   keyframe[2][3] = 2;
   keyframe[3][0] = 1;
   keyframe[3][1] = 4;
   keyframe[3][2] = 2;
   keyframe[3][3] = 3;
   keyframe[4][0] = 7;
   keyframe[4][1] = 3;
   keyframe[4][2] = 2;
   keyframe[4][3] = 4;
   keyframe[5][0] = keyframe[0][0];
   keyframe[5][1] = keyframe[0][1];
   keyframe[5][2] = keyframe[0][2];
   keyframe[5][3] = keyframe[0][3];
   float temp_eps = FLT_MIN;
   printf("eps is %f \n",temp_eps);
   compute_A C_A(order,m);
   C_A.run_compute(mu_r, mu_psi,  k_r, k_psi, t_index);
   //C_A.print_A();
   compute_Constraint C_C(order,n,m,k_r-1,k_psi);//don't know the reason why k_r-1
   C_C.compute_waypoint_C(keyframe,t_index);
   C_C.compute_pos_D_C(temp_eps);
   C_C.combine_Cb();
   //CGAL::Const_oneset_iterator<CGAL::Comparison_result>
   CGAL::Comparison_result r[C_C.C_size[0]]; //r(   CGAL::EQUAL); //r(   CGAL::SMALLER);
   CGAL::Quadratic_program_pricing_strategy strategy[] = {
       CGAL::QP_CHOOSE_DEFAULT,
       CGAL::QP_DANTZIG,
       CGAL::QP_PARTIAL_DANTZIG,
       CGAL::QP_BLAND,
       CGAL::QP_FILTERED_DANTZIG,
       CGAL::QP_PARTIAL_FILTERED_DANTZIG
   };
   bool fl[C_C.C_size[0]];
   float l[C_C.C_size[0]];
   for (int i=0; i<C_C.C_size[0]; i++)
   {
       fl[i] = false;
       l[i] = -20.1;
       if (i<40)
            r[i] = CGAL::EQUAL;
       else
            r[i] = CGAL::SMALLER;
   }
   bool fu[C_C.C_size[0]];
   float u[C_C.C_size[0]];
   float c[C_C.C_size[0]];
   float c0;
   for (int i=0; i<C_C.C_size[0]; i++)
   {
       fu[i] = false;
       u[i] = 20.8;
       c[i] = 0.0;
   }
    c0 = 1.0;
   printf("initiall done\n");
   printf("sizeof A is [%d] X [%d]\n",C_C.C_size[0],C_C.C_size[1]);
   printf("sizeof b is [%d] \n",C_C.b_size[1]);
   printf("sizeof D is [%d] X [%d]\n",C_A.A_col_row,C_A.A_col_row);
   //C_A.print_2D();
   //C_C.print_C();
   //C_C.print_b();

//solve Quadratic Program problem
   Program qp (C_C.C_size[0],C_C.C_size[1],C_C.C,C_C.b,r,fl,l,fu,u,C_A.D2_res,c,c0);
   //Program qp (C_C.C_size[0],C_C.C_size[1],C_C.C,C_C.b,r,fl,NULL,fu,NULL,C_A.D2_res,c,c0);

   CGAL::Quadratic_program_options options;
   options.set_verbosity(1);
   options.set_pricing_strategy(strategy[5]);
   options.set_auto_validation(true);
   //return 1;
   Solution s = CGAL::solve_quadratic_program(qp, ET(),options);
   assert (s.is_valid());
   std::cout << s;
   //Solution::Index_iterator it = s.basic_variable_indices_begin();
   //Solution::Index_iterator end = s.basic_variable_indices_end();
   //for (; it != end; ++it)
       //std::cout << *it << " ";
   //std::cout << std::endl;
   //Solution::Optimality_certificate_iterator result_it = s.optimality_certificate_begin();
   //Solution::Optimality_certificate_iterator result_end = s.optimality_certificate_end();
   //for (; result_it != result_end; ++result_it)
       //std::cout << *result_it << " ";
   //std::cout << std::endl;
   std::cout << " result: " << std::endl;
   Solution::Variable_value_iterator result_it = s.variable_values_begin();
   Solution::Variable_value_iterator result_end = s.variable_values_end();
   for (; result_it != result_end; ++result_it)
       std::cout << *result_it << " ";
   std::cout << std::endl;
   printf("done\n");
   for (int i=0;i <6;i++)
   {
       delete[] keyframe[i];
   }
   delete[] keyframe;
   return 0;
}
