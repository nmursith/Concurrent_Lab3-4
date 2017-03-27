#include <iostream>
#include <omp.h>
#include "paralleloptimized.h"

double doParallelImprovedMultiplication(int n){
 double start = omp_get_wtime();

 //double matrix2[n][n];
 //double results[n][n];
 int row, col,k;
 int c;
 #pragma omp parallel for collapse(2) private(row,col,k) shared(matrix1,matrix2,resultMatrix) //schedule(dynamic,50)
  for( row =0; row<n; row++){
   for( k=0; k<n; k++){
     for( col =0,c=0; col <n; col++){
         resultMatrix[row][col] += matrix1[row][k] * matrix2[k][col];

       }
       //cout<<resultMatrix[row][col] <<endl;
       //resultMatrix[row][col] = sum;
       //cout << resultMatrix[row][col] <<"  " ;
     }
   //cout <<endl;
 }


  double end = omp_get_wtime();
 //cout << endl << "Time for Parallel Multiplication " << n<< "x"<<n <<": " << (end-start)*1000 <<" milliseconds" << endl;
 return (end-start) * 1000;
}
