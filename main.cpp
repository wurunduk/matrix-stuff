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
    Matrix A, A1;    

    if(!(argc==3 || argc==2) || !(matrix_size = atoi(argv[1]))){
        PrintUsage(argv[0]);        
        return 1;
    }

    if(argc == 3) file_name = argv[2];


    auto t = clock();

    printf("loading file %s with matrix size %d\n", file_name, matrix_size);

    MatrixException result = A.CreateMatrix(matrix_size, matrix_size, file_name);
	MatrixException result2 = A1.CreateMatrix(matrix_size, matrix_size, file_name);

	ReportError(result);
	if(result != NO_ERROR) return 2;

	A.Print(12);
	auto b = A.GetRHSVector();
	auto x = A.Solve(&b);
	//got the solution, let's check how close it is now
	b = A1.GetRHSVector();
	auto b1 = A1*x;
	x.Print();
	printf("vector we want:\n");
	b.Print();
	printf("vector Ax=\n");
	b1.Print();

    t = clock()-t;
    printf("Elapsed time: %f\n", static_cast<float>(t)/CLOCKS_PER_SEC);
    return 0;
}
