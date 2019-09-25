#include "matrix.h"
#include "matrix_exception.h"
#include <ctime>
#include <iostream>

void PrintUsage(const char* executable_name){
    printf("Usage: %s n for auto generated matrix of size n\n %s n file to read matrix of size n from file\n", executable_name, executable_name);
}

int main(int argc, char* argv[]){
    
    char* file_name = nullptr;
    int matrix_size = 0;
	double *A, *x, *b, *A1, *x1, *b1, *Ax;

    if(!(argc==3 || argc==2) || !(matrix_size = atoi(argv[1]))){
        PrintUsage(argv[0]);        
        return 1;
    }

    if(argc == 3) file_name = argv[2];

    printf("loading file %s with matrix size %d\n", file_name, matrix_size);

	auto t = clock();
	
    MatrixException r1 = Matrix::CreateMatrix(&A, matrix_size, file_name);
	MatrixException r2 = Matrix::CreateMatrix(&A1, matrix_size, file_name);

	MatrixException r3 = Matrix::CreateVector(&b, matrix_size);
	MatrixException r4 = Matrix::CreateVector(&b1, matrix_size);

	MatrixException r5 = Matrix::CreateVector(&x, matrix_size);
	MatrixException r6 = Matrix::CreateVector(&x1, matrix_size);

	MatrixException r7 = Matrix::CreateVector(&Ax, matrix_size);

	ReportError(r1);
	ReportError(r2);
	ReportError(r3);
	ReportError(r4);
	ReportError(r5);
	ReportError(r6);
	ReportError(r7);
	if((r1|r2|r3|r4|r5|r6|r7) != NO_ERROR) return 2;
	
	Matrix::GetRHSVector(A, b, matrix_size);
	Matrix::GetRHSVector(A1, b1, matrix_size);

	Matrix::GetAnswerVector(x1, matrix_size);

	Matrix::Solve(A, b, x, matrix_size);

	Matrix::MultiplyMatrixByVector(A1, x, Ax, matrix_size);

	printf("Solution found: \n");
	Matrix::PrintVector(x, matrix_size);
	printf("Resulting Ax= vector: \n");
	Matrix::PrintVector(Ax, matrix_size);	
	printf("Intendent b vector: \n");
	Matrix::PrintVector(b1, matrix_size);

	double residual = Matrix::LengthVector(Matrix::SubstractVectors(Ax, b, matrix_size), matrix_size);
	double error = Matrix::LengthVector(Matrix::SubstractVectors(x, x1, matrix_size), matrix_size);

	printf("residual=%e\nerror=%e\n", residual, error);
	t = clock()-t;
    printf("Elapsed time: %f\n", static_cast<float>(t)/CLOCKS_PER_SEC);
    
    return 0;
}
