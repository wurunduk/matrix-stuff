#include "matrix.h"
#include "matrix_exception.h"

#include <cmath>

void Matrix::FillMatrix(double* m, const int size){
	for(int y = 0; y < size; y++)
		for(int x = 0; x < size; x++)
			m[y*size + x] = 1./(x+y+1);
}

void Matrix::GetAnswerVector(double* vector, const int size){
	for(int i = 0; i < size; i++)
		vector[i] = ((i+1)&1);
}

void Matrix::Print(const double* matrix, const int size, int print_size){
	int n = size;
    if(n > print_size) n = print_size;

    for(int y = 0; y < n-1; y++){
        for(int x = 0; x < n-1; x++)
            printf("%lf ", matrix[x + y*size]);
		printf("... %lf ", matrix[size-1 + y*size]);
        printf("\n");
    }

	printf("...\n");

	for(int x = 0; x < n-1; x++)
            printf("%lf ", matrix[x + (size-1)*size]);
	
	printf("... %lf ", matrix[size*size-1]);
    printf("\n");

}

void Matrix::Print(const double* matrix, const int size, const int* indexes, int print_size){
	int n = size;
    if(n > print_size) n = print_size;

    for(int y = 0; y < n-1; y++){
        for(int x = 0; x < n-1; x++)
            printf("%lf ", matrix[x + indexes[y]*size]);
		printf("... %lf ", matrix[size-1 + indexes[y]*size]);
        printf("\n");
    }

	printf("...\n");

	for(int x = 0; x < n-1; x++)
            printf("%lf ", matrix[x + indexes[(size-1)]*size]);
	
	printf("... %lf ", matrix[indexes[size-1]*size + (size-1)]);
    printf("\n");

}

void Matrix::PrintVector(const double* vector, const int size, int print_size){
	int n = size;
    if(n > print_size) n = print_size;

    for(int y = 0; y < n-1; y++){
        printf("%lf", vector[y]);
        printf("\n");
    }

	printf("...\n");
    printf("%lf\n", vector[size-1]);
}

void Matrix::PrintVector(const double* vector, const int size, const int* indexes, int print_size){
	int n = size;
    if(n > print_size) n = print_size;

    for(int y = 0; y < n-1; y++){
        printf("%lf", vector[indexes[y]]);
        printf("\n");
    }

	printf("...\n");
    printf("%lf\n", vector[indexes[size-1]]);
}

double Matrix::Length(const double* matrix, const int size){
	double max, sum;
	for(int y = 0; y < size; y++){
		sum = 0;
		for(int x = 0; x < size; x++){
			sum += fabs(matrix[x + y*size]);
		}
		if(y==0) max = sum;
		if(sum > max) max = sum;
	}

	return max;
}

double Matrix::LengthVector(const double* vector, const int size){
	double max=fabs(vector[0]);
	for(int x = 0; x < size; x++){
		if(fabs(vector[x]) > max) max = fabs(vector[x]);
	}

	return max;
}

void Matrix::GetRHSVector(const double* matrix, double* RHSVector, const int size){
	for(int y = 0; y < size; y++){
		double sum = 0;
		for(int x = 0; x < size; x+=2){
			sum += matrix[x + y*size];
		}
		RHSVector[y] = sum;
	}
}

void Matrix::Solve(double* matrix, double* rhs, double* answer, const int size){
	int offset = 0;
	//create a permutation, so we dont take time to actually move elements in the matrix
	auto indexes = new int[size];
	for(int i = 0; i < size; i++)
		indexes[i] = i;

	//when we complete a step of Gaussian algorithm, we should apply it again to the matrix of size m-1. 
	//For that we will just think of the next element on the diagonal as the first one.
	while(offset != size){
		//we want to find the element with lowest inverse length
		//with numbers it will be just 1/k, so we can search for the largest element in the column
        //Print(matrix, size, indexes, 15);
        //printf("\n");
		int minimal_length_index = offset;
		double max_length = matrix[offset + indexes[offset]*size];
		for(int y = offset; y < size; y++){
			double k = matrix[offset+indexes[y]*size];
			if(k>max_length){
				max_length = k;
				minimal_length_index = y;
			} 
		}
		
		//now we can swap the top and the lowest line
		int tempIndex = indexes[offset];
		indexes[offset] = indexes[minimal_length_index];
		indexes[minimal_length_index] = tempIndex;

		//now we want to normalize our line using the first element
		//this will make the first element of current top line = 1
		double k = 1.0/max_length;
		for(int x = offset; x < size; x++){
			matrix[x+indexes[offset]*size] *= k;
		}
		rhs[indexes[offset]] *= k;

		//substract the top line from all the lines below it
		for(int y = offset+1; y < size; y++){
			rhs[indexes[y]] -= rhs[indexes[offset]]*matrix[offset+indexes[y]*size];
			for(int x = size-1; x >= offset; x--){
				matrix[x+indexes[y]*size] -= matrix[x+indexes[offset]*size]*matrix[offset+indexes[y]*size];
			}
		}

		//repeat for the sub matrix
		offset += 1;
	}
	

	//at this point we have an upper diagonal matrix and we can get the answer
	//known as reverse step of Gauss algorithm
	for(int y = size-1; y >= 0; y--){
		double sum = 0.0;
		//technically this is an x coordinate
		for(int i = size-1; i > y; i--){
			sum += answer[i]*matrix[indexes[y]*size + i];
		}

		answer[y] = rhs[indexes[y]] - sum;
	}
	
	delete[] indexes;
}

void Matrix::SolveBlock(double* matrix, double* rhs, double* answer, const int size, const int block_size){
    int test = block_size;
    test += 1;
	int offset = 0;
    int index_size = (size+2)/block_size;
    
    double* block = nullptr;
    
	//create a permutation, so we dont take time to actually move elements in the matrix
	auto indexes = new int[index_size];
	for(int i = 0; i < index_size; i++)
		indexes[i] = i;

	//when we complete a step of Gaussian algorithm, we should apply it again to the matrix of size m-1. 
	//For that we will just think of the next element on the diagonal as the first one.
	while(offset != index_size){
		//we want to find the element with lowest inverse length
		//with numbers it will be just 1/k, so we can search for the largest element in the column
        //Print(matrix, size, indexes, 15);
        //printf("\n");
		int minimal_length_index = offset;
        GetBlock(matrix, block, 0, 0, size, block_size);
		double max_length = Length(block, block_size);
		for(int y = offset; y < index_size; y++){
            GetBlock(matrix, block, offset, indexes[y], size, block_size);
            //this will crash on last smaller blocks
			double k = Length(block, block_size);
			if(k>max_length){
				max_length = k;
				minimal_length_index = y;
			} 
		}
		
		//now we can swap the top and the lowest line(size+2)/block_size
		int tempIndex = indexes[offset];
		indexes[offset] = indexes[minimal_length_index];
		indexes[minimal_length_index] = tempIndex;

		//now we want to normalize our line using the first element
		//this will make the first element of current top line = 1
		double k = 1.0/max_length;
		for(int x = offset; x < size; x++){
			matrix[x+indexes[offset]*size] *= k;
		}
		rhs[indexes[offset]] *= k;

		//substract the top line from all the lines below it
		for(int y = offset+1; y < size; y++){
			rhs[indexes[y]] -= rhs[indexes[offset]]*matrix[offset+indexes[y]*size];
			for(int x = size-1; x >= offset; x--){
				matrix[x+indexes[y]*size] -= matrix[x+indexes[offset]*size]*matrix[offset+indexes[y]*size];
			}
		}

		//repeat for the sub matrix
		offset += 1;
	}
	
	

	//at this point we have an upper diagonal matrix and we can get the answer
	//known as reverse step of Gauss algorithm
	for(int y = size-1; y >= 0; y--){
		double sum = 0.0;
		//technically this is an x coordinate
		for(int i = size-1; i > y; i--){
			sum += answer[i]*matrix[indexes[y]*size + i];
		}

		answer[y] = rhs[indexes[y]] - sum;
	}
	
	delete[] indexes;
}

void Matrix::GetBlockSize(int* q, int* m, const int x, const int y, const int matrix_size, const int block_size){
    *q = x*block_size;
    if(*q + block_size > matrix_size) *q = matrix_size%block_size;
    
    *m = y*block_size;
    if(*m + block_size > matrix_size) *m = matrix_size%block_size;
}

void Matrix::GetBlock(const double* A, double* block, const int x, const int y, const int matrix_size, const int block_size){
    int q, m;
    
    GetBlockSize(&q, &m, x, y, matrix_size, block_size);
    
    block = new double[q*m];
    
    for(int i = 0; i < m; i++){
        for(int j = 0; j < q-4; j+=4){
            block[j + i*q] = A[x*block_size + j + y*block_size*matrix_size + i*matrix_size];
            block[j + i*q + 1] = A[x*block_size + j + 1 + y*block_size*matrix_size + i*matrix_size];
            block[j + i*q + 2] = A[x*block_size + j + 2 + y*block_size*matrix_size + i*matrix_size];
            block[j + i*q + 3] = A[x*block_size + j + 3 + y*block_size*matrix_size + i*matrix_size];
        }
    }
}

void Matrix::PutBlock(double* A, const double* block, const int x, const int y, const int matrix_size, const int block_size){
    int q, m;
    
    GetBlockSize(&q, &m, x, y, matrix_size, block_size);
    
    for(int i = 0; i < m; i++)
        for(int j = 0; j < q-4; j+= 4){
            A[x*block_size + j + y*block_size*matrix_size + i*matrix_size] = block[j + i*q];
            A[x*block_size + j + 1 + y*block_size*matrix_size + i*matrix_size] = block[j + 1 + i*q];
            A[x*block_size + j + 2 + y*block_size*matrix_size + i*matrix_size] = block[j + 2 + i*q];
            A[x*block_size + j + 3 + y*block_size*matrix_size + i*matrix_size] = block[j + 3 + i*q];
        }
}

void Matrix::NullMatrix(double* matrix, const int size){
    for(int x = 0; x < size; x++)
        matrix[x] = 0;
}

MatrixException Matrix::ReadMatrix(double* matrix, const int size, const char* file_name){
    FILE* f = fopen(file_name, "r");
    if(!f) return CAN_NOT_OPEN;
    
    for(int i = 0; i < size*size; i++){
        if(fscanf(f, "%lf", &matrix[i]) != 1){
            fclose(f);
            return FILE_CORRUPT;
        }
    }
    
    return NO_ERROR;
}

MatrixException Matrix::InitMatrix(double* matrix, const int size, const char* file_name){
    if(file_name == nullptr) {
        FillMatrix(matrix, size);
        return NO_ERROR;
    }
    else return ReadMatrix(matrix, size, file_name);
}

MatrixException Matrix::CreateVector(double** vector, const int size){
    *vector = new double[size];
    if(!*vector) return MatrixException::NOT_ENOUGH_MEMORY;
    
    NullMatrix(*vector, size);

    return NO_ERROR;
}

MatrixException Matrix::CreateMatrix(double** matrix, const int size){
    *matrix = new double[size*size];
    if(!*matrix) return MatrixException::NOT_ENOUGH_MEMORY;
    
    NullMatrix(*matrix, size*size);
    
    return NO_ERROR;
}

MatrixException Matrix::CreateMatrix(double** matrix, const int size, const char* file_name){
    //if no name was supplied, we want to initialize and generate the matrix
    MatrixException res = CreateMatrix(matrix, size);
	if(res == NO_ERROR){
		return InitMatrix(*matrix, size, file_name);
	}
	else return res;
}

double* Matrix::MultiplyMatrixByVector(const double* matrix, const double* vector, double* answer, const int size){    
    for(int y = 0; y < size; y++){
		double sum = 0.0;
        for(int x = 0; x < size; x++){
            sum += matrix[x + y*size]*vector[x];
        }
		answer[y] = sum;
    }
	return answer;
}

double* Matrix::SubstractVectors(double* v1, const double* v2, const int size){
	for(int x = 0; x < size; x++)
		v1[x] -= v2[x];

	return v1;
}


double Matrix::GetError(double* vector, const int size){
    double max=fabs(vector[0]-1);
	for(int x = 0; x < size; x++){
		if(fabs(vector[x] - ((x+1)&1)) > max) max = fabs(vector[x] - ((x+1)&1));
	}

	return max;
}





























