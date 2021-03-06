#ifndef _MATRIX_H_
#define _MATRIX_H_

#include "matrix_exception.h"
#include <memory>

namespace Matrix {
	void FillMatrix(double* matrix, const int w, const int h);
	void GetAnswerVector(double* vector, const int size);
    void PrintClean(const double* matrix, const int w, const int h);
    void Print(const double* matrix, const int size, int print_size = 10);
	void Print(const double* matrix, const int size, const int* indexes, int print_size = 10);
	void PrintVector(const double* vector, const int size, int print_size = 10);
	void PrintVector(const double* vector, const int size, const int* indexes, int print_size = 10);
	double Length(const double* matrix, const int w, const int h);
	double LengthVector(const double* vector, const int size);

    void NullMatrix(double* matrix, const int size);
	void EMatrix(double* matrix, const int size);
    MatrixException ReadMatrix(double* matrix, const int w, const int h, const char* file_name);
    MatrixException InitMatrix(double* matrix, const int w, const int h, const char* file_name);

	void GetRHSVector(const double* matrix, double* RHSVector, const int size);

    void GetBlock(const double* A, double* block, const int x, const int y, const int x1, const int y1, const int matrix_size);
    void PutBlock(double* A, const double* block, const int x, const int y, const int x1, const int y1, const int matrix_size);

	MatrixException CreateVector(double** vector, const int size);
    MatrixException CreateMatrix(double** matrix, const int w, const int h);
    MatrixException CreateMatrix(double** matrix, const int w, const int h, const char* file_name);

	int GetInverseMatrix(double* matrix, double* matrixReversed, int m, double norm, int* transposition_m);

	double* MultiplyMatrices(double* matrix1, double* matrix2, double* res, const int h, const int wh, const int w);
	double* MultiplyMatricesNMK(double* matrix1, double* matrix2, double* res, const int h, const int wh, const int w);
	void MultiplyMatricesNN_4(double* a, double* b, double* c, const int m);
	void MultiplyMatricesNN(double* a, double* b, double* c, const int n);
    double* SubstractMatrices(double* matrix, const double* matrix2, const int w, const int h);

    double GetError(double* vector, const int size);

	int SolveBlock(double* matrix, double* rhs, double* answer, const int size, const int block_size, double* array[]);
}

#endif
