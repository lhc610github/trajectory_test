//#include <CGAL/basic.h>
//#include <CGAL/QP_models.h>
//#include <CGAL/QP_functions.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <list>
class my_traject
{
    public:
        my_traject(int temp_order, int temp_m);
        ~my_traject();
        void computeA(float mu_r, float mu_psi, int k_r, int k_psi, float *t);
        void print_A();
    private:
        int order;
        int m;
        int A_col_row;
        float ** c_A_res;
        void init(int temp_order, int temp_m);
};

my_traject::
my_traject(int temp_order,int temp_m)
{
    init(temp_order,temp_m);
    c_A_res = new float *[A_col_row];
    for (int i=0 ; i < A_col_row ; i++)
    {
        c_A_res[i] = new float[A_col_row];
    }
    std::cout<<A_col_row<<std::endl;
}

my_traject::
~my_traject()
{
    delete c_A_res;
}


void
my_traject::init(int temp_order, int temp_m)
{
    order = temp_order;
    m = temp_m;
    A_col_row = 4*temp_m*(temp_order+1);
}

void
my_traject::computeA(float mu_r, float mu_psi, int k_r, int k_psi, float *t)
{


    int l_polyn_r = order - k_r + 1;
    int polynomial_r[l_polyn_r]; 
    int i,j,k;
    for (i = 0; i < l_polyn_r; i++ )
    {
        for (j = 0; j < k_r; j++ )
        {
            polynomial_r[i] *= order - i - j;
        }
    }



    int l_polyn_psi = order - k_psi + 1;
    int polynomial_psi[l_polyn_psi]; 
    for (i = 0; i < l_polyn_psi; i++ )
    {
        for (j = 0; j < k_psi; j++ )
        {
            polynomial_psi[i] *= order - i - j;
        }
    }

    float A[m][4][order+1][order+1];
    memset(A,0,sizeof(A));
    //float A_res[A_col_row][A_col_row];
    //printf("size A : %d \n",(int)sizeof(A_res));
    for (i = 0; i < A_col_row ; i++)
    memset(&(*c_A_res[i]),0,sizeof(float)*A_col_row);

    int order_t_r,order_t_psi;
    int count_state;
    int base_count_i;
    int base_count_state;
    float mu_x;
    
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < (order+1); j++)
        {
            for (k = j; k < (order+1); k++)
            {
                // Position
                if(j < l_polyn_r && k < l_polyn_r)
                {
                    order_t_r = ((order-k_r-j)+(order-k_r-k));
                    if (j == k)
                    {
                        A[i][0][j][k] = pow(polynomial_r[j],2)/(order_t_r + 1) \
                                    *(pow(t[i+1],(order_t_r + 1)) - pow(t[i],(order_t_r + 1)));
                        A[i][1][j][k] = pow(polynomial_r[j],2)/(order_t_r + 1) \
                                    *(pow(t[i+1],(order_t_r + 1)) - pow(t[i],(order_t_r + 1)));
                        A[i][2][j][k] = pow(polynomial_r[j],2)/(order_t_r + 1) \
                                    *(pow(t[i+1],(order_t_r + 1)) - pow(t[i],(order_t_r + 1)));
                    }
                    else
                    {
                        A[i][0][j][k] = 2*polynomial_r[j]*polynomial_r[k]/(order_t_r + 1) \
                                    *(pow(t[i+1],(order_t_r + 1)) - pow(t[i],(order_t_r + 1)));
                        A[i][1][j][k] = 2*polynomial_r[j]*polynomial_r[k]/(order_t_r + 1) \
                                    *(pow(t[i+1],(order_t_r + 1)) - pow(t[i],(order_t_r + 1)));
                        A[i][2][j][k] = 2*polynomial_r[j]*polynomial_r[k]/(order_t_r + 1) \
                                    *(pow(t[i+1],(order_t_r + 1)) - pow(t[i],(order_t_r + 1)));
                    }

                }
                else
                {
                          A[i][0][j][k] = 0;
                          A[i][1][j][k] = 0;
                          A[i][2][j][k] = 0;
                }


                // Yaw

                if(j < l_polyn_psi && k < l_polyn_psi)
                {
                    order_t_psi = ((order - k_psi - j)+(order - k_psi - k));
                    if (j == k)
                    {
                        A[i][3][j][k] = pow(polynomial_psi[j],2)/(order_t_psi + 1) \
                                    *(pow(t[i+1],(order_t_psi + 1)) - pow(t[i],(order_t_psi + 1)));
                    }
                    else
                    {
                        A[i][3][j][k] = 2*polynomial_psi[j]*polynomial_psi[k]/(order_t_psi + 1) \
                                    *(pow(t[i+1],(order_t_psi + 1)) - pow(t[i],(order_t_psi + 1)));
                    }

                }
                else
                {
                        A[i][3][j][k] = 0;
                }
            }
        }
        // TODO:  Make A symetric

        base_count_i = i*4*(order+1);
        for (count_state = 0; count_state < 4 ;count_state ++)
        {
            if(count_state == 3)
            {
                mu_x = mu_psi;
            }
            else
            {
                mu_x = mu_r;
            }
            base_count_state = count_state*(order+1);
            for ( j = 0 ; j < (order+1) ; j++ )
            {
                for (k = 0 ; k <= j ; k++)
                {
                    if ( k == j )
                    {
                        c_A_res[base_count_i + base_count_state + j][base_count_i + base_count_state + k] = A[i][count_state][j][k] * mu_x;
                    }
                    else
                    {
                        c_A_res[base_count_i + base_count_state + j][base_count_i + base_count_state + k] = (A[i][count_state][j][k]+A[i][count_state][k][j])/2*mu_x;
                        c_A_res[base_count_i + base_count_state + k][base_count_i + base_count_state + j] = c_A_res[base_count_i + base_count_state + j][base_count_i + base_count_state + k];
                    }
                }
            }
        }
    }

}

void
my_traject::print_A()
{
    for(int i=0;i < A_col_row;i++)
    {
        for(int j=0;j < A_col_row;j++)
            printf("%.1f ",c_A_res[i][j]);
        std::cout<<std::endl;
    }
}


int 
main()
{
   int order = 6;
   int m = 5;
   float mu_r = 1;
   float mu_psi = 1;
   int k_r = 4;
   int k_psi = 2;
   float t_index[6] = {0,2,4,6,8,10};
   my_traject C_A(order,m);
   C_A.computeA(mu_r, mu_psi,  k_r, k_psi, t_index);
   C_A.print_A();
   return 0;
}
