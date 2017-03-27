// paralleloptimized
#ifndef PARALLELOPTIMIZED_H
#define PARALLELOPTIMIZED_H

extern double ** matrix1;//= populate(n,n);
extern double ** matrix2;//= populate(n,n);
extern double ** resultMatrix;
//extern double ** transposedMatrix;
//extern double timeResults[3][500] = {0};

double doParallelImprovedMultiplication(int n);

#endif
