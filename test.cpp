#include <iostream>
#include <stdio.h>
#include <omp.h>
#include <time.h>

using namespace std;

#define SIZE 10000
#define N 100
#define S 25

int n = N;
int s = S;

double a[SIZE],b[SIZE],c[SIZE];

// Initializing the matrices

void mat_init(double *a, double *b, int n)
{
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
            a[i*n + j] = 1;

    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
            b[i*n + j] = 2;

}


void mat_multi(double *a, double *b, double *c, int n)
{
    //double start_t = omp_get_wtime();
    clock_t start=clock();

    int i,j,k;

    #pragma omp num_threads(5) for private(i,j,k)
    for( i=0; i<n; i++)
        for( j=0; j<n; j++)
            for( k=0; k<n; k++)
                c[i*n+j]+=a[i*n+k]*b[k*n+j];

    start = clock() - start;

    double ms = ((double)(start)*1000)/CLOCKS_PER_SEC;
    //double stop_t = omp_get_wtime();
    cout<<"Naive multiplication requires "<<ms<<"ms"<<endl;
}

void mat_print(double *a, int n)
{
    cout<<endl<<endl<<endl<<"************************************************************"<<endl;

    for (int i = 0; i < n; ++i)
    {
        cout<<endl;
        for (int j = 0; j < n; ++j)
        {
            /* code */
            cout<<a[i*n+j]<<" ";
        }
    }

    cout<<endl<<endl<<endl<<"************************************************************"<<endl;
}

void mat_empty(double *a, int n)
{
    for (int i = 0; i < n; ++i)
    {
        /* code */
        for (int j = 0; j < n; ++j)
        {
            /* code */
            c[i*n+j]=0;
        }
    }
}

void tiled_mat_multiply(double *a, double *b, double *c, int n)
{
    int i,j,k,i1,j1,k1,tid;

    clock_t start = clock();

    double start_t,stop_t;

    omp_set_nested(1);
    #pragma omp parallel shared(a,b,c) private(i1,j1,k1,i,j,k,tid) num_threads(omp_get_num_procs())
    {

        /*
        tid = omp_get_thread_num();

        if(tid == 0)
        {
            cout<<"Master thread encountered "<<endl<<endl;
            start_t = omp_get_wtime();
        }
        */

        #pragma omp for
        for ( i1 = 0; i1 < n; i1+=s)
            for ( j1 = 0; j1 < n; j1+=s)
                for ( k1 = 0; k1 < n; k1+=s)
                    for( i=i1; i <i1+s && i<n; i++)
                        for ( j=j1; j< j1+s && j<n; ++j)
                            for( k=k1; k< k1+s && k<n; ++k)
                                c[i*n+j]+=a[i*n+k]*b[k*n+j];
    }


    /*if(tid==0)
    {
        stop_t = omp_get_wtime();
    }*/

    start = clock() - start;
    double ms = ((double)(start)*1000)/CLOCKS_PER_SEC;

    cout<<"Tiled matrix multiplication requires "<<ms<<"ms"<<endl;

}

int main()
{
    mat_init(a,b,n);
    mat_multi(a,b,c,n);
    //mat_print(c,n);
    mat_empty(c,n);
    tiled_mat_multiply(a,b,c,n);
    //mat_print(c,n);
    return 0;
}
