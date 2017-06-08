#include "compute_Constraint.h"
compute_Constraint::
compute_Constraint(int temp_order, int temp_n, int temp_m, int temp_kr, int temp_kpsi)
{
    order = temp_order;
    n = temp_n;
    m = temp_m;
    k_r = temp_kr;
    k_psi = temp_kpsi;

    keyframe = new float *[m+1];
    for (int i=0 ;i < m+1 ;i++)
    {
        keyframe[i] = new float[n];
    }
    t = new float[m+1];
    C1 = new float *[2*m*n];
    for (int i=0 ;i < 2*m*n ;i++)
    {
        C1[i] = new float[n*(order+1)*m];
        memset(C1[i],0,sizeof(float)*(n*(order+1)*m));
    }
    C1_size[0] = 2*m*n;
    C1_size[1] = n*(order+1)*m;
    b1 = new float[2*m*n];
    b1_size[0] = 1;
    b1_size[1] = 2*m*n;

    int num = 2*m*(n-1)*k_r;
    C2 = new float *[num];
    for (int i=0 ;i < num ;i++)
    {
        C2[i] = new float[n*(order+1)*m];
        memset(C2[i],0,sizeof(float)*(n*(order+1)*m));
    }
    C2_size[0] = num;
    C2_size[1] = n*(order+1)*m;
    b2 = new float[num];
    b2_size[0] = 1;
    b2_size[1] = num;
    // C's row is num of variables
    // C's col is num of equations 
    C_size[1] = C1_size[0] + C2_size[0];
    C_size[0] = C1_size[1];//same as C2
    C = new float *[C_size[0]];
    for (int i=0 ;i < C_size[0] ;i++)
    {
        C[i] = new float[C_size[1]];
        memset(C[i],0,sizeof(float)*(C_size[1]));
    }
    b_size[0] = 1;
    b_size[1] = b1_size[1] + b2_size[1];
    b = new float[b_size[1]];
    memset(b,0,sizeof(float)*(b_size[1]));
}

compute_Constraint::
~compute_Constraint()
{
    for (int i=0 ;i < m+1 ;i++)
    {
        delete[] keyframe[i];
    }
    delete[] keyframe;
    delete[] t;
    for (int i=0 ;i < 2*m*n ;i++)
    {
        delete[] C1[i];
    }
    delete[] C1;
    delete[] b1;
    //for (int i=0 ;i < C2_size[0] ;i++)
    //{
        //delete[] C2[i];
    //}
    delete[] C2;
    delete[] b2;
    //for (int i=0 ;i < C_size[0];i++)
    //{
        //delete[] C[i];
    //}
    delete[] C;
    delete[] b;
}
bool
compute_Constraint::
compute_waypoint_C(float ** temp_keyframe,float * temp_t) //keyframe m*4 matrix
{
    int i,j,k;
    for (i=0 ;i < m+1 ;i++)
    {
        t[i] = temp_t[i];
        for (j=0 ;j < n ;j++)
        {
            keyframe[i][j] = temp_keyframe[i][j];
        }
    }
    float *waypoint;
    waypoint = new float [4];
    float *values;
    values = new float [order+1];
    for (i=0 ;i < m ;i++)
    {
        waypoint = keyframe[i];
        //printf("waypoint ");
        //print_Vector(waypoint,4);
        if(i == 0) // Initial and  Final Position
        {   // Initial
            //memset(values,0,sizeof(float)*(order+1));
            for (j=0 ;j<(order+1) ;j++)
            {
                // computation of polynomials
                values[j] = pow(t[i],(order-j));
            }
            for (k=0 ;k<n ;k++ )
            {
                for (j=0 ;j<(order+1) ;j++)
                {
                    C1[k][k*(order+1) + j] = values[j];
                }
            }
            for (j=0 ;j<n ;j++ )
            {
                b1[j]=waypoint[j];
            }

            // Final
            for (j=0 ;j<(order+1) ;j++)
            {
                // computation of polynomials
                values[j] = pow(t[m],(order-j));
            }
            for (k=0 ;k<n ;k++ )
            {
                for (j=0 ;j<(order+1) ;j++)
                {
                    C1[k+n][(m-1)*(order+1)*n + k*(order+1) + j] = values[j];
                }
            }
            for (j=0 ;j<n ;j++)
            {
                b1[n+j] = waypoint[j];
            }
        }
        else
        {
            //Elsewhere
            for (j=0 ;j<(order+1) ;j++)
            {
                // computation of polynomials
                values[j] = pow(t[i],(order-j));
            }

            for (k=0 ;k<n ;k++ )
            {
                for (j=0 ;j<(order+1) ;j++)
                {
                    C1[k + 2*n*i][(i-1)*(order+1)*n + k*(order+1) + j] = values[j];
                }
            }
            for (j=0 ;j<n ;j++)
            {
                b1[2*n*i + j] = waypoint[j];
            }

            for (k=0 ;k<n ;k++ )
            {
                for (j=0 ;j<(order+1) ;j++)
                {
                    C1[k + 2*n*i+n][(i)*(order+1)*n + k*(order+1) + j] = values[j];
                }
            }
            for (j=0 ;j<n ;j++)
            {
                b1[2*n*i + n+j] = waypoint[j];
            }
        }
    }
    //printf("C1 ");
    //print_Matrix(C1,2*m*n,m*(order+1)*n);
    //printf("b1 ");
    //print_Vector(b1,2*m*n);
    delete[] waypoint;
    delete[] values;
    return true;
}

bool
compute_Constraint::
compute_pos_D_C(float temp_eps)
{
    int num = 2*m*(n-1)*k_r;
    for (int i=0 ;i < num ;i++)
    {
        b2[i] = temp_eps;
    }
    //int i,j,k;
    // constraintData.m
        float *** C_r; // constraintData_r
        C_r = new float **[m];
        for (int i = 0;i < m;i++)
        {
            C_r[i] = new float *[k_r];
            for (int j=0;j<k_r;j++)
            {
                C_r[i][j] = new float [3];
                memset(C_r[i][j],0,sizeof(float)*3);
            }
        }
        //velocity
        if(k_r >= 1)
        {
            for (int i = 0;i < 3;i ++)
            {
                C_r[0][0][i] = 0; //At starting position
                C_r[0][1][i] = 0; //At starting position
                C_r[0][2][i] = 0; //At starting position
                C_r[1][2][i] = 0; //At starting position
                C_r[2][2][i] = 0; //At starting position
                C_r[3][2][i] = 0; //At starting position
            }
            for (int i = 1;i < m;i ++)
            {
                C_r[i][0][0] = temp_eps; // x,y velocities
                C_r[i][0][1] = temp_eps;
                C_r[i][0][2] = temp_eps; // z velocities
            }
        }
        //acceleration
        if(k_r >= 2)
        {
                C_r[0][1][2] = 0; //At starting position
            for (int i = 1;i < m;i ++)
            {
                C_r[i][1][0] = temp_eps; // x,y velocities
                C_r[i][1][1] = temp_eps;
                C_r[i][1][2] = temp_eps; // z velocities
            }
        }
        //jerk
        //snap
        float ** C_psi; // constraintData_psi
        C_psi = new float *[m];
        for (int i = 0;i < m;i++)
        {
            C_psi[i] = new float [k_psi];
            memset(C_psi[i],0,sizeof(float)*k_psi);
        }
        //velocity
        if(k_psi >= 1)
        {
            C_psi[0][0] = 0;
        }
        //acceleration
        //jerk
        //snap
    float values[order+1];
    bool continuity[n-1];
    for (int i = 0;i < m;i++)
    {
        for (int k = 0;k < k_r;k++)
        {
            if(i == 0)
            {
                // Initial
                for (int j = 0;j < order+1;j++)
                {
                    if((j+k+1)>order)
                    {
                        values[j] = 0;
                    }
                    else
                    {
                        int tempcoeffs = (order-j);
                        for (int time_coeffs = 1; time_coeffs<=k; time_coeffs++)
                        {
                            tempcoeffs *= (order-j-time_coeffs);
                        }
                        values[j] = tempcoeffs*pow(t[i],(order-j-k-1));
                    }
                }
                //printf("Values ");
                //print_Vector(values,order+1);
                
                for (int l_con = 0;l_con < (n-1);l_con++)
                {
                    continuity[l_con] = (C_r[i][k][l_con] == temp_eps)?true:false;
                    if(continuity[l_con])
                    {
                        for (int j = 0;j < order+1;j++)
                        {
                            C2[l_con+k*(n-1)][l_con*(order+1)+j] = values[j];
                            C2[l_con+k*(n-1)][(m-1)*(order+1)*n + l_con*(order+1)+j] = -values[j];
                        }
                        b2[l_con+k*(n-1)] = 0;
                    }
                    else
                    {
                        for (int j = 0;j < order+1;j++)
                        {
                            C2[l_con+k*(n-1)][l_con*(order+1)+j] = values[j];
                        }
                        b2[l_con+k*(n-1)] = C_r[i][k][l_con];

                    }
                }
                // Final
                for (int j = 0;j < order+1;j++)
                {
                    if((j+k+1)>order)
                    {
                        values[j] = 0;
                    }
                    else
                    {
                        int tempcoeffs = (order-j);
                        for (int time_coeffs = 1; time_coeffs<=k; time_coeffs++)
                        {
                            tempcoeffs *= (order-j-time_coeffs);
                        }
                        values[j] = tempcoeffs*pow(t[i],(order-j-k-1));
                    }
                }
                for (int l_con = 0;l_con < (n-1);l_con++)
                {
                    if(!continuity[l_con])
                    {
                        for (int j = 0;j < order+1;j++)
                        {
                            C2[l_con+k*(n-1)+(n-1)*k_r][(m-1)*(order+1)*n + l_con*(order+1)+j] = values[j];
                        }
                        b2[l_con+k*(n-1)+(n-1)*k_r] = C_r[i][k][l_con];
                    }
                }
            }
            else
            { // Elsewhere
                for (int j = 0;j < order+1;j++)
                {
                    if((j+k+1)>order)
                    {
                        values[j] = 0;
                    }
                    else
                    {
                        int tempcoeffs = (order-j);
                        for (int time_coeffs = 1; time_coeffs<=k; time_coeffs++)
                        {
                            tempcoeffs *= (order-j-time_coeffs);
                        }
                        values[j] = tempcoeffs*pow(t[i],(order-j-k-1));
                    }
                }

                for (int l_con = 0;l_con < (n-1);l_con ++)
                {
                    continuity[l_con] = (C_r[i][k][l_con]==temp_eps)? true:false;
                    if (continuity[l_con])
                    {
                        for (int j = 0;j < order+1;j++)
                        {
                            C2[l_con + k*(n-1) + 2*i*(n-1)*k_r][(i-1)*(order+1)*n + l_con*(order+1)+j] = values[j];
                            C2[l_con + k*(n-1) + 2*i*(n-1)*k_r][i*(order+1)*n + l_con*(order+1)+j] = -values[j];
                        }
                        b2[l_con + k*(n-1) + 2*i*(n-1)*k_r] = 0;
                    }
                    else
                    {
                        for (int j = 0;j < order+1;j++)
                        {
                            C2[l_con + k*(n-1) + 2*i*(n-1)*k_r][(i-1)*(order+1)*n + l_con*(order+1)+j] = values[j];
                            C2[l_con + k*(n-1) + 2*i*(n-1)*k_r + (n-1)*k_r][i*(order+1)*n + l_con*(order+1)+j] = values[j];
                        }
                        b2[l_con + k*(n-1) + 2*i*(n-1)*k_r] = C_r[i][k][l_con];
                        b2[l_con + k*(n-1) + 2*i*(n-1)*k_r + (n-1)*k_r] = C_r[i][k][l_con];
                    }
                    
                }
            }
        }
    }
    //printf("C2 ");
    //print_Matrix(C2,num,n*(order+1)*m);
    //printf("b2 ");
    //print_Vector(b2,num);

        for (int i = 0;i < m;i++)
        {
            for (int j=0;j < k_r;j++)
            {
                delete[] C_r[i][j];
            }
            delete[] C_r[i];
        }
        delete[] C_r;

        for (int i = 0;i < m;i++)
        {
            delete[] C_psi[i];
        }
        delete[] C_psi;
    return true;
}

void
compute_Constraint::
print_Matrix(float **sth,int row,int col)
{
    printf("Matrix:\n");
    for(int i=0 ;i < row; i++)
    {
        printf("| ");
        for(int j=0 ;j < col; j++)
        {
            printf("%.2f ",sth[i][j]);
        }
        printf("|\n");
    }
}

void
compute_Constraint::
print_Vector(float *sth,int length)
{
    printf("Vector:\n");
    printf("| ");
    for(int i=0 ;i < length; i++)
    {
        printf("%.2f ",sth[i]);
    }
    printf("|\n");
}

void
compute_Constraint::
print_C()
{
    for(int i=0; i < C_size[0]/2; i++)
    {
        for(int j=0; j < C_size[1]; j++)
            printf("%.2f ",C[i][j]);
        printf("\n");
    }
}

void
compute_Constraint::
print_b()
{
    for(int j=0; j < b_size[1]; j++)
        printf("%.2f ",b[j]);
    printf("\n");
}

void
compute_Constraint::
combine_Cb()
{
// C
    // C1
float max_C = -1;
float min_C = 1;
float max_b = -1;
float min_b = 1;
    printf("start combine\n");
    for(int i=0; i < C1_size[0]; i++)
    {
        for(int j=0; j < C1_size[1]; j++)
        {
        C[j][i] = C1[i][j];
        max_C = (C[j][i] > max_C)? C[j][i]:max_C;
        min_C = (C[j][i] < min_C)? C[j][i]:min_C;
        }
    }
    // C2
    for(int i=0; i < C2_size[0]; i++)
    {
        for(int j=0; j < C2_size[1]; j++)
        {
        int row_C = i+C1_size[0];
        C[j][row_C] = C2[i][j];
        max_C = (C[j][row_C] > max_C)? C[j][row_C]:max_C;
        min_C = (C[j][row_C] < min_C)? C[j][row_C]:min_C;
        }
    }
    printf("C done\n");
// b
    // b1
    for(int i=0; i < b1_size[1]; i++)
    {
        b[i] = b1[i];
        max_b = (b[i] > max_b)? b[i]:max_b;
        min_b = (b[i] < min_b)? b[i]:min_b;
    }
    // b2
    for(int i=0; i < b2_size[1]; i++)
    {
        int row_b = i+b1_size[1];
        b[row_b] = b2[i];
        max_b = (b[row_b] > max_b)? b[row_b]:max_b;
        min_b = (b[row_b] < min_b)? b[row_b]:min_b;
    }
    printf("C[%d][%d]  combine done \n in C   max: %.2f  min: %.2f \n",C_size[0],C_size[1],max_C,min_C);
    printf("b[%d][%d]  combine done \n in b   max: %.2f  min: %.2f \n",b_size[0],b_size[1],max_b,min_b);
    
}
