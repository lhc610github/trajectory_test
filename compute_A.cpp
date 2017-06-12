#include "compute_A.h"
compute_A::
compute_A(int temp_order,int temp_m)
{
    init(temp_order,temp_m);
    c_A_res = new float *[A_col_row];
    for (int i=0 ; i < A_col_row ; i++)
    {
        c_A_res[i] = new float [A_col_row];
        memset(c_A_res[i],0,sizeof(float)*A_col_row);
    }
    D2_res = new float *[A_col_row];
    for (int i=0 ; i < A_col_row ; i++)
    {
        D2_res[i] = new float [A_col_row];//[i+1];
        memset(D2_res[i],0,sizeof(float)*(A_col_row));//(i+1));
    }
}

compute_A::
~compute_A()
{
    //for (int i=0 ; i < A_col_row ; i++)
    //{
       //delete[] D2_res[i];
    //}
    delete[] D2_res;
    //for (int i=0 ; i < A_col_row ; i++)
    //{
        //printf("%d  ",i);
        //delete[] c_A_res[i];
    //}
    delete[] c_A_res;
}


void
compute_A::init(int temp_order, int temp_m)
{
    order = temp_order;
    m = temp_m;
    A_col_row = 4*temp_m*(temp_order+1);
}

void
compute_A::run_compute(float mu_r, float mu_psi, int k_r, int k_psi, float *t)
{


    int l_polyn_r = order - k_r + 1;
    int * polynomial_r;
    polynomial_r = new int[l_polyn_r];
    int i,j,k;
    for (i = 0; i < l_polyn_r; i++ )
    {
        for (j = 0; j < k_r; j++ )
        {
            if (j == 0)
                polynomial_r[i] = (order - i - j);
            else
                polynomial_r[i] *= (order - i - j);
        }
    }
    //print_vector(polynomial_r,l_polyn_r);



    int l_polyn_psi = order - k_psi + 1;
    int *polynomial_psi;
    polynomial_psi = new int[l_polyn_psi]; 
    for (i = 0; i < l_polyn_psi; i++ )
    {
        for (j = 0; j < k_psi; j++ )
        {
            if (j == 0)
                polynomial_psi[i] = order - i - j;
            else
                polynomial_psi[i] *= order - i - j;
        }
    }
    //print_vector(polynomial_psi,l_polyn_psi);

    float ****A;
    A = new float ***[m];
    for (i = 0; i<m ;i++)
    {
        A[i] = new float **[4];
        for (j = 0; j < 4; j++)
        {
            A[i][j] = new float *[order+1];
            for (k=0; k<(order+1); k++)
            {
                A[i][j][k] = new float[order+1];
                memset(A[i][j][k],0,sizeof(float)*(order+1));
            }
        }
    }

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
            //printf("step: %d  state: %d \n",i+1,count_state+1);
            //print_matrix(A[i][count_state],(order+1),(order+1));
            for ( j = 0 ; j < (order+1) ; j++ )
            {
                for (k = 0 ; k <= j ; k++)
                {
                    if ( k == j )
                    {
                        c_A_res[base_count_i + base_count_state + j][base_count_i + base_count_state + k] = (float)A[i][count_state][j][k] * mu_x;
                    }
                    else
                    {
                        c_A_res[base_count_i + base_count_state + j][base_count_i + base_count_state + k] = ((float)A[i][count_state][j][k]+(float)A[i][count_state][k][j])/2 * mu_x;
                        c_A_res[base_count_i + base_count_state + k][base_count_i + base_count_state + j] = (float)c_A_res[base_count_i + base_count_state + j][base_count_i + base_count_state + k];
                    }
                }
            }
        }
    }
    get_2D();

    for (i = 0; i<m ;i++)
    {
        for (j = 0; j < 4; j++)
        {
            for (k=0; k<(order+1); k++)
            {
               delete[] A[i][j][k];
            }
            delete[] A[i][j];
        }
        delete[] A[i];
    }
    delete[] A;
    delete[] polynomial_r;
    delete[] polynomial_psi;
}
void
compute_A::print_2D()
{
    printf("2D :\n");
    for (int i=0; i< A_col_row; i++)
    {
        for (int j=0; j<i+1; j++)
        {
           printf("%.2f ",D2_res[i][j]);
        }
        printf("\n");
    }
}
void
compute_A::print_A()
{
    print_matrix(c_A_res,A_col_row,A_col_row);
}

void
compute_A::print_vector(int *sth, int length)
{
    printf("vector: ");
    for (int i=0; i<length ;i++)
        printf("%d ",sth[i]);
    printf("\n");
}

void 
compute_A::print_matrix(float **sth, int row, int col)
{
    printf("matrix: \n");
    for (int i=0; i<row ;i++)
    {
        printf("| ");
        for (int j=0; j<col ;j++)
        {
            printf("%.2f ",sth[i][j]);
        }
        printf("|\n");
    }
    printf("\n");
}

void
compute_A::get_2D()
{
    float max_2D = -1;
    float min_2D = 1;
    for (int i=0; i< A_col_row; i++)
    {
        for (int j=0; j< A_col_row; j++)//i+1; j++)
        {
            D2_res[i][j] = 2*c_A_res[i][j];
            max_2D = (D2_res[i][j]>max_2D)?D2_res[i][j]:max_2D;
            min_2D = (D2_res[i][j]<min_2D)?D2_res[i][j]:min_2D;
        }
    }
    printf("get 2D \n in 2D max: %.2f   min: %.2f \n",max_2D,min_2D);
}

