#include <iostream>
#include <omp.h>
#include <ctime>
#include <cstdlib>
#include "paralleloptimized.h"
#include <fstream>
#include <cmath>

#define SAMPLES 30
using namespace std;
double ** matrix1;//= populate(n,n);
double ** matrix2;//= populate(n,n);
double ** resultMatrix;
double ** transposedMatrix;
double timeResults[3][500] = {0};



double **  populate(int row, int col);
double doSerialMultiplication(int n);
double doParallelMultiplication(int n);
//double doParallelImprovedMultiplication(int n);
void transpose(int n);
void calculateTime(double  elapsed [][SAMPLES], int step);
void doPopulateall(int n);
void freeMemory();
int main()

{
	//#pragma omp parallel num_threads(4)
// 	#pragma omp parallel
// 	#pragma omp_set_num_threads(4)
// //     cout << "Hello, world.\n";
// //
// // //		#ifdef _OPENMP
// // cout << " Number of processors available:" << omp_get_num_procs() << " MAX Number of threads " << omp_get_max_threads() << endl;
// // cout.flush();
//#endif
		double  timeElapsed [3][SAMPLES]= {0};

		ofstream myfile;
    myfile.open ("speedup108.csv");
	//populate(matrix2,n,n);
	for(int step=100; step<=1000; step+=100){
		for(int samples =0; samples <SAMPLES ; samples++){
						doPopulateall(step);
						double time1 = doSerialMultiplication(step);
					//std::cout << '\n';
						double time2 = doParallelMultiplication(step);
						//transpose(step);
					//std::cout << '\n';
						double time3 = doParallelImprovedMultiplication(step);

						// cout << "Results for n = "<<step << endl;
						// cout << "	Serial 			= " << time1<<endl;
						// cout << "	Parallel 		= " <<time2<<endl;
						//
						// cout << "	Parallel Improved	= " <<time3<<endl;
						// cout << "	Speedup 		= "<< time1/time2<<endl;
						// cout << "	Speedup Improved	= "<< time1/time3<<endl;
						// myfile << time1/time3 <<endl;
						//cout << endl;

						timeElapsed[0][samples] = time1;
						timeElapsed[1][samples] = time1/time2;
						timeElapsed[2][samples] = time1/time3;
						//cout << "	Improved = " << results[2][r][k]<<endl<<endl;

						freeMemory();
		}
		calculateTime(timeElapsed, step);


	}
	myfile.close();
}

void calculateTime(double elapsed[][SAMPLES], int step){
//  double elapsed[3][samples];
  double average[3] = {0,0,0};
  double variance[3] = {0,0,0};

  int i,CASE;

	#pragma omp parallel for schedule(dynamic,50)  private(i,CASE) shared(average, variance)
	for(int CASE = 0; CASE<3; CASE++){
  	for(i=0; i<SAMPLES; i++){
    		average[CASE]+=elapsed[CASE][i];
  	}
		  average[CASE] = average[CASE]/SAMPLES;

	}




#pragma omp parallel for  private(i,CASE) shared(average, variance)
  for( CASE=1; CASE<3; CASE++){

  	for(i=0; i<SAMPLES; i++){
      variance[CASE] += pow((elapsed[CASE][i] - average[CASE]), 2);
    	}
    	variance[CASE] = variance[CASE]/SAMPLES;
        //printf("%f\t", sqrt(variance[j]));
      variance[CASE] = sqrt(variance[CASE]);

      printf("sample needed CASE %d Step %d  average %lf ---> %f\n", CASE, step, average[CASE], pow(((100 * variance[CASE] * 1.96)/(5* average[CASE])),2));

  }
  printf("\n" );
}

void doPopulateall(int n){
	matrix1= populate(n,n);
	matrix2= populate(n,n);
	resultMatrix =  (double **)malloc(n * sizeof(double *)); //new double*[numRow];
	transposedMatrix =  (double **)malloc(n * sizeof(double *)); //new double*[numRow];

	for (int i = 0; i < n; ++i)
	{
		resultMatrix[i] =   new double[n];//(double *)malloc(n * sizeof(double *));// new double[numCol];
		//resultMatrix[i] = {0};
		transposedMatrix[i] =  new double[n];//(double *)malloc(n * sizeof(double *));// new double[numCol];
		//transposedMatrix[i]={0};
	}

}

void transpose(int n) {
    int i,j;
	#pragma omp parallel for collapse(2) schedule(dynamic,50)  private(i,j) shared(matrix1,transposedMatrix) num_threads(omp_get_num_procs())
    for(i=0; i<n; i++) {
        for(j=0; j<n; j++) {
            transposedMatrix[i][j] = matrix2[j][i];

        }
    }
}

void freeMemory(){
	free(matrix1);
	free(matrix2);
	free(resultMatrix);
	free(transposedMatrix);

}
//  double doParallelImprovedMultiplication(int n){
//
//
// 	double start = omp_get_wtime();
//
// 	//double matrix2[n][n];
// 	//double results[n][n];
// 	int row, col,k;
//
// 	#pragma omp parallel for collapse(2) private(row,col,k) shared(matrix1,matrix2,resultMatrix) //schedule(dynamic,50)
// 	//collapse(2)
// 	for( row =0; row<n; row++){
//
// 		for( k=0; k<n; k++){
// 			for( col =0; col <n; col++){
// 				//matrix1[row][col]  = rand();
// 				//double sum =0;
// 				//#pragma omp parallel for reduction(+:sum)
//
//
// 				//	cout << "** " <<matrix1[row][k] << "**"<<  transposedMatrix[col][k]<<endl;
// 					resultMatrix[row][col] += matrix1[row][k] * matrix2[k][col];
//
// 				}
// 				//cout<<resultMatrix[row][col] <<endl;
// 				//resultMatrix[row][col] = sum;
// 				//cout << resultMatrix[row][col] <<"  " ;
// 			}
// 		//cout <<endl;
// 	}
//
// /*
// 	int TILE = 16;
// 	for ( int i=0; i<n; i+=TILE )
// 	        for ( int j=0; j<n; j+=TILE )
// 	            for ( int k=0; k<n; k+=TILE )
// 	            // Regular multiply inside the tiles
// 	                for ( int y=i; y<i+TILE; y++ ){
// 	                    for ( int x=j; x<j+TILE; x++ ){
// 														double sum=0;
// 		                        for ( int z=k; z<k+TILE; z++ ){
// 		                             sum += matrix1[y][z] * transposedMatrix[x][z]; //A(y,z)*B(z,x);
// 															 }
// 															 cout<<sum <<endl;
// 															 resultMatrix[y][x] = sum;
// 														}
//
// 												}
//
//
// */
//
//
// 	 double end = omp_get_wtime();
// 	//cout << endl << "Time for Parallel Multiplication " << n<< "x"<<n <<": " << (end-start)*1000 <<" milliseconds" << endl;
// 	return (end-start) * 1000;
//  }



 double doParallelMultiplication(int n){


	double start = omp_get_wtime();

	//double matrix2[n][n];
	//double results[n][n];
int row, col, k;
	#pragma omp parallel for  collapse(2) private(row,col,k) shared(matrix1,matrix2,resultMatrix) //schedule(dynamic,50)
	//collapse(2)
	for(row  =0; row<n; row++){
		for( col =0; col <n; col++){
			//matrix1[row][col]  = rand();
			//double sum =0;
			resultMatrix[row][col] = 0;
			for( k=0; k<n; k++){
				resultMatrix[row][col] += matrix1[row][k] * matrix2[k][col];

			}
			//cout << sum << "   ";
			//cout<<resultMatrix[row][col] <<endl;
			//resultMatrix[row][col] = sum;
			//cout << resultMatrix[row][col] <<"  " ;
		}
		//cout <<endl;
	}
	//cout <<endl;


	 double end = omp_get_wtime();
	//cout << endl << "Time for Parallel Multiplication " << n<< "x"<<n <<": " << (end-start)*1000 <<" milliseconds" << endl;
	return (end-start) * 1000;
 }

 double doSerialMultiplication(int n){
	double start = omp_get_wtime();
	 int row, col, k;
	//double matrix2[n][n];
	//double results[n][n];


	for( row =0; row<n; row++){
		for( col =0; col <n; col++){
			//matrix1[row][col]  = rand();
			double sum =0;
			resultMatrix[row][col]=0;
			for( k=0; k<n; k++){
				resultMatrix[row][col] += matrix1[row][k] * matrix2[k][col];

			}
			//cout << sum << "   ";
			//cout<<resultMatrix[row][col] <<endl;
			//resultMatrix[row][col] = sum;
			//cout << resultMatrix[row][col] <<"  " ;
		}
		//cout <<endl;
	}
	//cout <<endl;


	 double end = omp_get_wtime();
	//cout << endl << "Time for Serial Multiplication " << n<< "x"<<n <<": " << (end-start)*1000 <<" milliseconds" << endl;
	return (end-start) * 1000;
 }


 double ** populate(int numRow, int numCol){
	srand( (unsigned)time( NULL ) );
	///double matrix[numRow][numCol]; //[numRow][numCol];

	 double ** matrix =  (double **)malloc(numRow * sizeof(double *)); //new double*[numRow];

	#pragma omp parallel for
	for (int i = 0; i < numCol; ++i)
	{
		matrix[i] = (double *)malloc(numCol * sizeof(double *));// new double[numCol];
	}

	#pragma omp parallel for collapse(2)
	for(int row =0; row<numRow; row++){

		for(int col =0; col <numCol; col++){
			matrix[row][col] = rand()%10000;

		}

	}

	return matrix;

}
