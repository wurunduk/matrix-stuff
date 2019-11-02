#ifndef _MATRIX_H_
#define _MATRIX_H_

#include "matrix_exception.h"
#include <memory>

namespace Matrix{
	void FillMatrix(double* matrix, const int w, const int h);	
	void GetAnswerVector(double* vector, const int size);
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

	int GetInverseMatrix(double* matrix, double* inverse, int size);

	double* MultiplyMatrixByVector(const double* matrix, const double* vector, double* answer, const int size);	
	double* MultiplyMatrices(double* matrix, const double* matrix2, const int q, const int n, const int m);
    double* SubstractMatrices(double* matrix, const double* matrix2, const int n, const int m);
    
	double* SubstractVectors(double* v1, const double* v2, const int size);
    
    double GetError(double* vector, const int size);

	void SolveBlock(double* matrix, double* rhs, double* answer, const int size, const int block_size, double* array[]);

	//Matrix operator-(const Matrix& m);

    //Matrix operator*(const Matrix& m);
    //Matrix operator*(const double& k);
    //Matrix& operator*=(const double& k);
}

#endif
