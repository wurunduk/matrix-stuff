#include "matrix.h"
#include "matrix_exception.h"
#include <ctime>
#include <iostream>

void set_fpu_exception_mask (void);

static void ReportError(MatrixException ex) {
    switch (ex) {
        case NO_ERROR:
        return;
        case UNKNOWN_ERROR:
        printf("Unknown error happened while working with matrix, error code %d.\n", static_cast<int>(ex));
        return;
        case NOT_ENOUGH_MEMORY:
        printf("Could not allocate enough memory for the matrix.\n");
        return;
        case CAN_NOT_OPEN:
        printf("Could not open matrix file.\n");
        return;
        case FILE_CORRUPT:
        printf("Matrix file did not have enough numbers to fill the given array.\n");
        return;
    }
}

void PrintUsage(const char* executable_name) {
    printf("Usage: %s n m for auto generated matrix of size n with block size m.\n %s n m 'file' to read matrix of size n from file with block size m.\n", executable_name, executable_name);
}

int main(int argc, char* argv[]) {
    set_fpu_exception_mask();

    char* file_name = nullptr;
    int matrix_size = 0;
    int block_size = 0;
    int e = 0;
    double *A, *x, *b, *Ax;
    const int temp_matrix_count = 17;
    double* temps[temp_matrix_count];

    if (!(argc==4 || argc==3) || !(matrix_size = atoi(argv[1])) || !(block_size = atoi(argv[2]))) {
        PrintUsage(argv[0]);
        return 1;
    }

    if (argc == 4) file_name = argv[3];

    printf("Loading file %s with matrix size %d, block size %d.\n", file_name, matrix_size, block_size);

    int end = matrix_size - (matrix_size/block_size)*block_size;

    // Initialize used variables
    e |= Matrix::CreateMatrix(&A, matrix_size, matrix_size, file_name);
    e |= Matrix::CreateVector(&x, matrix_size);
    e |= Matrix::CreateVector(&b, matrix_size);

    e |= Matrix::CreateVector(&Ax, matrix_size);

    e |= Matrix::CreateMatrix(&temps[0], block_size, block_size);//temp m*m blocks
    e |= Matrix::CreateMatrix(&temps[1], block_size, block_size);
    e |= Matrix::CreateMatrix(&temps[2], block_size, block_size);
    e |= Matrix::CreateMatrix(&temps[3], block_size, block_size);
    e |= Matrix::CreateMatrix(&temps[4], block_size, end);//temp m*e blocks
    e |= Matrix::CreateMatrix(&temps[5], block_size, end);
    e |= Matrix::CreateMatrix(&temps[6], block_size, end);
    e |= Matrix::CreateMatrix(&temps[7], block_size, end);
    e |= Matrix::CreateMatrix(&temps[8], block_size, block_size);//inversed block
    e |= Matrix::CreateVector(&temps[9], block_size);//temp b vector blocks
    e |= Matrix::CreateVector(&temps[10], block_size);
    e |= Matrix::CreateVector(&temps[11], block_size);
    e |= Matrix::CreateVector(&temps[12], block_size);
    e |= Matrix::CreateMatrix(&temps[13], end, end);//temp e*e blocks
    e |= Matrix::CreateMatrix(&temps[14], end, end);
    e |= Matrix::CreateVector(&temps[15], end);//temp vector end blocks
    e |= Matrix::CreateVector(&temps[16], end);

    ReportError((MatrixException)e);
    if (e != NO_ERROR)
    {
        if (A) delete[] A;
        if (x) delete[] x;
        if (b) delete[] b;
        if (Ax) delete[] Ax;

        for (int i = 0; i < temp_matrix_count; i++)
            if (temps[i]) delete[] temps[i];
        return 2;
    }
    Matrix::GetRHSVector(A, b, matrix_size);

    auto t = clock();

    int res = Matrix::SolveBlock(A, b, x, matrix_size, block_size, temps);

    t = clock()-t;

    if (res == 0) {
        // At this point matrix A and vector b are wrong, reinitialize them
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

        printf("Residual=%e\nerror=%e\n", residual, error);
        printf("Elapsed time: %f\n", static_cast<float>(t)/CLOCKS_PER_SEC);
    } else {
        printf("Matrix could not be inverted.\n");
    }

    if (A) delete[] A;
    if (x) delete[] x;
    if (b) delete[] b;
    if (Ax) delete[] Ax;

    for (int i = 0; i < temp_matrix_count; i++)
        if (temps[i]) delete[] temps[i];

    return 0;
}
