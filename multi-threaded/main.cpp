#include "matrix.h"
#include "matrix_exception.h"
#include <ctime>
#include <iostream>
#include <pthread.h>

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
	}
}

void PrintUsage(const char* executable_name){
    printf("Usage: %s n m p for auto generated matrix of size n with block size m and p threads\n %s n m p 'file' to read matrix of size n from file with block size m and p threads\n", executable_name, executable_name);
}

int main(int argc, char* argv[]){
    set_fpu_exception_mask();
    char* file_name = nullptr;
    int matrix_size = 0;
    int block_size = 0;
	int thread_count = 0;
    int e = 0;				//exceptions count to track errors
    double *A, *x, *b, *Ax;	
	arg* args = nullptr;
	pthread_barrier_t barrier;

    if(!(argc==5 || argc==4) || !(matrix_size = atoi(argv[1])) || !(block_size = atoi(argv[2])) || !(thread_count = atoi(argv[3]))){
        PrintUsage(argv[0]);
        return 1;
    }

    if(argc == 5) file_name = argv[4];

    printf("loading file %s with matrix size %d, block size %d, thread count %d\n", file_name, matrix_size, block_size, thread_count);

	pthread_barrier_init(&barrier, 0, thread_count);

	//size of end block of the matrix
    int end = matrix_size - (matrix_size/block_size)*block_size;
	
    //initialize used variables
    e |= Matrix::CreateMatrix(&A, matrix_size, matrix_size, file_name);
    e |= Matrix::CreateVector(&x, matrix_size);
    e |= Matrix::CreateVector(&b, matrix_size);

    e |= Matrix::CreateVector(&Ax, matrix_size);

    args = new arg[thread_count];

	for(int i = 0; i < thread_count; i++){
		//InitializeTempAddresses(&(args[i]->adr), &e);
		args[i].matrix = A;
		args[i].rhs = b;
		args[i].x = x;
		args[i].size = matrix_size;
		args[i].block_size = block_size;

		args[i].return_value = 0;

		args[i].thread_count = thread_count;
		args[i].thread_id = id;
		args[i].work_time = 0;
		args[i].matrix = A;

		args[i].barrier = &barrier;
	}
    
    ReportError((MatrixException)e);
    if(e != NO_ERROR)
    {
        delete[] A;
        delete[] x;
        delete[] b;
        delete[] Ax;
        
        for(int i = 0; i < thread_count; i++)
            DeleteTempAddresses(&(args[i]->adr));

		delete[] args;
        return 2;
    }

	Matrix::GetRHSVector(A, b, matrix_size);

    auto t = clock();
    
    int res = Matrix::SolveBlock(A, b, x, matrix_size, block_size, temps);
    
    t = clock()-t;
    
    if(res == 0){
        //at this point matrix A and vector b are wrong, but we can reinitialize them;
        Matrix::InitMatrix(A, matrix_size, matrix_size, file_name);
        Matrix::GetRHSVector(A, b, matrix_size);
        
        Matrix::MultiplyMatrices(A, x, Ax, matrix_size, matrix_size, 1);
        
        printf("Solution found: \n");
        Matrix::PrintVector(x, matrix_size);
        printf("Resulting Ax= vector: \n");
        Matrix::PrintVector(Ax, matrix_size);	
        printf("Intendent b vector: \n");
        Matrix::PrintVector(b, matrix_size);
        
        double residual = Matrix::LengthVector(Matrix::SubstractMatrices(Ax, b, matrix_size, 1), matrix_size);
        double error = Matrix::GetError(x, matrix_size);

        printf("Residual=%e Error=%e n=%d m=%d k=%d elapsed=%.2f\n", residual, error, matrix_size, block_size, thread_count, static_cast<float>(t)/CLOCKS_PER_SEC);
    }
    else{
        printf("Matrix could not be inverted.\n");
    }
    
	delete[] A;
	delete[] x;
	delete[] b;
	delete[] Ax;
    
    for(int i = 0; i < thread_count; i++)
    	DeleteTempAddresses(&(args[i]->adr));

	delete[] args;

	pthread_barrier_destroy(&barrier);
    
    return 0;
}
