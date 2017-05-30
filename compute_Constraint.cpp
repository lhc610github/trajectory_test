#include "compute_Constraint.h"
compute_Constraint::
compute_Constraint(int temp_order, int temp_n, int temp_m)
{
    order = temp_order;
    n = temp_n;
    m = temp_m;
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
        memset(&(*C1[i]),0,sizeof(float)*(n*(order+1)*m));
    }
    b1 = new float[2*m*n];
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
    delete[] waypoint;
    delete[] values;
    return true;
}
