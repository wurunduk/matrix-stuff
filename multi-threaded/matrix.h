#ifndef _MATRIX_H_
#define _MATRIX_H_

#include "matrix_exception.h"
#include <memory>

typedef struct{
	double* inverse_block;
	double* block;
    double* block_temp;
    double* block_temp_im;
    double* block_temp_sub;
    double* block_me;
    double* block_me_temp;
    double* block_me_temp_im;
    double* block_me_temp_sub;
    double* vector_block;
    double* vector_block_temp;
    double* vector_block_temp_im;
	double* block_ee;
	double* block_ee_temp;
	double* vector_e;
	double* vector_e_temp;
} address_array;

typedef struct{
    double* matrix;
    double* rhs;
	double* answer;
    char* file_name;
	int size;
	int block_size;

	int return_value;

	int thread_count;
	int thread_id;
    
    double work_time;
    
    pthread_barrier_t* barrier;
} arg;

namespace Matrix{
	void InitializeTempAddresses(addresses_array* adr, int block_size, int end);
	void DeleteTempAddresses(addresses_array* adr);
    
    void AttachMatrices(arg* in);

    void PrintClean(const double* matrix, const int w, const int h);
    void Print(const double* matrix, const int size, int print_size = 10);
	void Print(const double* matrix, const int size, const int* indexes, int print_size = 10);
	void PrintVector(const double* vector, const int size, int print_size = 10);
	void PrintVector(const double* vector, const int size, const int* indexes, int print_size = 10);
	double Length(const double* matrix, const int w, const int h);
	double LengthVector(const double* vector, const int size);
    
    void GetAnswerVector(double* vector, const int size);
    void GetRHSVector(const double* matrix, double* RHSVector, const int size);
    
    void NullMatrix(double* matrix, const int size);
	void EMatrix(double* matrix, const int size);
    MatrixException FillMatrix(double* matrix, const int w, const int h);	
    MatrixException ReadMatrix(double* matrix, const int w, const int h, const char* file_name);
    
    void GetBlock(const double* A, double* block, const int x, const int y, const int x1, const int y1, const int matrix_size);
    void PutBlock(double* A, const double* block, const int x, const int y, const int x1, const int y1, const int matrix_size);

	int GetInverseMatrix(double* matrix, double* matrixReversed, int m, double norm, int* transposition_m);

	double* MultiplyMatrices(double* matrix1, double* matrix2, double* res, const int h, const int wh, const int w);
	double* MultiplyMatricesNMK(double* matrix1, double* matrix2, double* res, const int h, const int wh, const int w);
	void MultiplyMatricesNN_4(double* a, double* b, double* c, const int m);
	void MultiplyMatricesNN(double* a, double* b, double* c, const int n);
    double* SubstractMatrices(double* matrix, const double* matrix2, const int w, const int h);
    
    double GetError(double* vector, const int size);

	void* SolveBlock(void* in);
}

#endif
