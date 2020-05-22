#include "matrix.h"
#include "matrix_exception.h"
#include "mpi.h"
#include <ctime>
#include <iostream>
#include <sys/sysinfo.h>

#include <cmath>
#include <string.h>

void set_fpu_exception_mask (void);

static void ReportError(MatrixException ex){
	switch (ex){
		case NO_ERROR:
		return;
		case UNKNOWN_ERROR:
		printf("Unknown error happened while working with matrix, error code %d\n", static_cast<int>(ex));
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
		case UNINVERTABLE:
		printf("Matrix could not be inverted.\n");
		return;
	}
}

void PrintUsage(const char* executable_name){
	printf("Usage: mpirun -n 4 %s n m for auto generated matrix of size n with block size m\n mpirun -n 4 %s n m p 'file' to read matrix of size n from file with block size m\n", executable_name, executable_name);
}

int main(int argc, char* argv[]){
	MPI_Init (&argc, &argv);
	set_fpu_exception_mask();
	
	char* file_name = nullptr;
	int matrix_size = 0;
	int block_size = 0;
	int process_count = 0;
	int process_id = 0;
	int local_error_id = 0;		//store error id of this process if any occured
	int global_error_id = 0;

	double *A, *x, *b, *Ax;	

	MPI_Comm_size(MPI_COMM_WORLD, &process_count);
	MPI_Comm_rank(MPI_COMM_WORLD, &process_id);

	if(!(argc == 4 || argc == 3) || !(matrix_size = atoi(argv[1])) || !(block_size = atoi(argv[2]))){
		if(process_id == 0) PrintUsage(argv[0]);
		MPI_Finalize();
		return 1;
	}

	if(argc == 4) file_name = argv[3];

	if(process_id == 0)
		printf("Loading file %s with matrix size %d, block size %d\n", file_name, matrix_size, block_size);


	A = new double[matrix_size*matrix_size];
	x = new double[matrix_size*3];
	b = x + matrix_size;
	Ax = b + matrix_size;

	arg args;

	args.matrix = A;
	args.rhs = b;
	args.answer = x;
	args.size = matrix_size;
	args.block_size = block_size;
	args.file_name = file_name;

	args.minimal_index = 0;
	args.minimal_norm = 0;

	args.return_value = 0;

	args.process_count = process_count;
	args.process_id = process_id;
	args.work_time = .0;
	args.cpu_time = .0;
	args.matrix = A;
	
	Matrix::SolveBlock(&args);

	//get any errors
	/*for(int i = 0; i < thread_count; i++){
		e |= args[i].return_value;
		if(args[i].return_value != 0){
			printf("Thread %d error:\n\t", i);
			ReportError((MatrixException)args[i].return_value);
		}
	}*/
	
	if(e == 0){
		//at this point matrix A and vector b are wrong, but we can reinitialize them;
		if(file_name) 
			Matrix::ReadMatrix(A, matrix_size, matrix_size, file_name);
		else
			Matrix::FillMatrix(A, matrix_size, matrix_size);
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

		printf("Residual=%e Error=%e Time=%.2lf n=%d m=%d k=%d\n", residual, error, args[0].work_time, matrix_size, block_size, thread_count);

		for(int i = 0; i < thread_count; i++){
			printf ("Cpu time of thread %d: %.2lf, Time of thread %d: %.2lf\n", i, args[i].cpu_time, i, args[i].work_time);
		}
	}
	
	delete[] A;
	delete[] x;

	delete[] args;
	
	MPI_Finalize ();
	return 0;
}
