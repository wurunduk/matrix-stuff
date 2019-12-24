#ifndef _MATRIX_H_
#define _MATRIX_H_

#include "matrix_exception.h" 

typedef struct{
	double value;
	int count;
} deigen_value;

MatrixException FillMatrix(double* m, const int w, const int h, const int type);
MatrixException ReadMatrix(double* matrix, const int w, const int h, const char* file_name);

double Length(const double* matrix, const int w, const int h);
void PrintClean(const double* matrix, const int w, const int h);
void PrintMatrix(const double* matrix, const int size, int print_size);
void PrintEigenVector(const deigen_value* vector, const int size, int print_size);

#endif
