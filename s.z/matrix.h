#ifndef _MATRIX_H_
#define _MATRIX_H_

#include "matrix_exception.h" 

MatrixException FillMatrix(double* m, const int w, const int h);
MatrixException ReadMatrix(double* matrix, const int w, const int h, const char* file_name);

double Length(const double* matrix, const int w, const int h);
void PrintClean(const double* matrix, const int w, const int h);
void PrintMatrix(const double* matrix, const int size, int print_size);
void PrintVector(const double* vector, const int size, int print_size);

#endif
