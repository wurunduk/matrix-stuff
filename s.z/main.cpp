#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "matrix.h"
#include "alg.h"
#include "matrix_exception.h"

void set_fpu_exception_mask (void);

static void ReportError(MatrixException ex){
	switch (ex){
		case NO_ERROR:
		return;
		case UNKNOWN_ERROR:
		printf("unknown error happened while working with matrix, error code %d\n", static_cast<int>(ex));
		return;
		case NOT_ENOUGH_MEMORY:
		printf("Could not allocate enough memory for the matrix.\n");
		return;
		case CAN_NOT_OPEN:
		printf("Could not open the matrix file.\n");
		return;
		case FILE_CORRUPT:
		printf("Matrix file did not have enough numbers to fill the given array.\n");
		return;
        case NON_SYMMETRICAL:
		printf("Matrix is not symmetrical.\n");
		return;
	}
}

void PrintUsage(const char* executable_name){
    printf("Usage: %s n eps for auto generated matrix of size n with precision eps\n %s n eps 'file' to read matrix of size n from file precision eps\n", executable_name, executable_name);
}

int main(int argc, char* argv[])
{
	set_fpu_exception_mask();
	char* file_name = nullptr;
	int matrix_size = 0;
	double* a = nullptr;
	deigen_value* values = nullptr;
	int found_values_amount = 0;
	double eps;
	int e = 0;			//MatrixException

	if(!(argc==4 || argc==3) || !(matrix_size = atoi(argv[1]))){
        PrintUsage(argv[0]);        
        return 1;
    }
    
    eps = atof(argv[2]); 

	if(argc == 4) file_name = argv[3];

	printf("Loading file %s with matrix size %d, precision %e\n", file_name, matrix_size, eps);

	a = new double[matrix_size*matrix_size];
	values = new deigen_value[matrix_size];
	
	if(file_name){
		e |= ReadMatrix(a, matrix_size, matrix_size, file_name);
	}else{
		int type = 0;
		printf("Input matrix type 0 - |i-j|, 1 - Hilbert, 2 - diagonal ones, 3 - diagonal twoes\n");
		while(scanf("%d", &type) != 1){}
		e |= FillMatrix(a, matrix_size, matrix_size, type);
	}
	
	ReportError((MatrixException)e);
	if(e != NO_ERROR){
		delete[] a;
		delete[] values;
		return 1;
	}

	PrintMatrix(a, matrix_size, 10);
	
	double trace = .0;
    double inv1 = .0;
	double inv2 = .0;

	for (int i = 0; i < matrix_size; ++i)
	{
		trace += a[i * matrix_size + i];
		for (int j = 0; j < matrix_size; ++j)
			inv1 += a[i * matrix_size + j] * a[j * matrix_size + i];
	}

	auto t = clock();
	int total_iterations = FindValues(a, values, matrix_size, &found_values_amount, eps);
	t = clock() - t;

	for (int i = 0; i < matrix_size;)
	{
		trace -= values[i].value*values[i].count;
		inv2 += values[i].value*values[i].value*values[i].count;
		i += values[i].count;
	}

	printf("Found %d eigen values:\n", found_values_amount);
    PrintEigenVector(values, matrix_size, 12);

	printf("Elapsed time: %.2f, Iterations: %d, Trace difference: %e, inv2: %e\n", (double)t / CLOCKS_PER_SEC, total_iterations, fabs(trace), fabs(sqrt(inv2) - sqrt(inv1)));

	delete[] a;
	delete[] values;

	return 0;
}
