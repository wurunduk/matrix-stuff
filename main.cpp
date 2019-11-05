#include "matrix.h"
#include "matrix_exception.h"
#include <ctime>
#include <iostream>

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
	}
}

void PrintUsage(const char* executable_name){
    printf("Usage: %s n m for auto generated matrix of size n with block size m\n %s n m 'file' to read matrix of size n from file with block size m\n", executable_name, executable_name);
}

int main(int argc, char* argv[]){
    
    char* file_name = nullptr;
    int matrix_size = 0;
    int block_size = 0;
    int e = 0;
    double *A, *x, *b, *Ax;
    double *block_m, *block_me, *block_ee, *vec_m, *vec_e;
    double* temps[12];

    if(!(argc==4 || argc==3) || !(matrix_size = atoi(argv[1])) || !(block_size = atoi(argv[2])) ){
        PrintUsage(argv[0]);        
        return 1;
    }

    if(argc == 3) file_name = argv[3];

    printf("loading file %s with matrix size %d, block size %d\n", file_name, matrix_size, block_size);

	auto t = clock();
	
    //initialize used variables
    e |= Matrix::CreateMatrix(&A, matrix_size, matrix_size, file_name);
    e |= Matrix::CreateVector(&x, matrix_size);
    e |= Matrix::CreateVector(&b, matrix_size);

    e |= Matrix::CreateVector(&Ax, matrix_size);

    e |= Matrix::CreateMatrix(&block_m, block_size, block_size);
    e |= Matrix::CreateMatrix(&block_me, block_size, matrix_size%block_size);
    e |= Matrix::CreateMatrix(&block_ee, matrix_size%block_size, matrix_size%block_size);
    e |= Matrix::CreateVector(&vec_m, block_size);
    e |= Matrix::CreateVector(&vec_e, matrix_size%block_size);

    ReportError((MatrixException)e);
    if(e != NO_ERROR)
    {
        if(A) delete[] A;
        if(x) delete[] x;
        if(b) delete[] b;
        if(Ax) delete[] Ax;
		if(block_m) delete[] block_m;
		if(block_me) delete[] block_me;
		if(block_ee) delete[] block_ee;
		if(vec_m) delete[] vec_m;
		if(vec_e) delete[] vec_e;
        return 2;
    }
	Matrix::GetRHSVector(A, b, matrix_size);

	//Matrix::Solve(A, b, x, matrix_size);
    Matrix::SolveBlock(A, b, x, matrix_size, block_size, temps);
    
    //at this point matrix A and vector b are wrong, but we can reinitialize them;
    Matrix::InitMatrix(A, matrix_size, matrix_size, file_name);
    Matrix::GetRHSVector(A, b, matrix_size);
    
	Matrix::MultiplyMatrixByVector(A, x, Ax, matrix_size);
    
	//printf("Solution found: \n");
	//Matrix::PrintVector(x, matrix_size);
	//printf("Resulting Ax= vector: \n");
	//Matrix::PrintVector(Ax, matrix_size);	
	//printf("Intendent b vector: \n");
	//Matrix::PrintVector(b, matrix_size);

	double residual = Matrix::LengthVector(Matrix::SubstractVectors(Ax, b, matrix_size), matrix_size);
	//double error = Matrix::LengthVector(Matrix::SubstractVectors(x, x1, matrix_size), matrix_size);
    double error = Matrix::GetError(x, matrix_size);

	printf("Residual=%e\nerror=%e\n", residual, error);
	t = clock()-t;
    printf("Elapsed time: %f\n", static_cast<float>(t)/CLOCKS_PER_SEC);
    
	if(A) delete[] A;
	if(x) delete[] x;
	if(b) delete[] b;
	if(Ax) delete[] Ax;
	if(block_m) delete[] block_m;
	if(block_me) delete[] block_me;
	if(block_ee) delete[] block_ee;
	if(vec_m) delete[] vec_m;
	if(vec_e) delete[] vec_e;
    
    return 0;
}
