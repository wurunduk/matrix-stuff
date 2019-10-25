#include "matrix.h"
#include "matrix_exception.h"

#include <cmath>

void Matrix::FillMatrix(double* m, const int w, const int h){
	for(int y = 0; y < h; y++)
		for(int x = 0; x < w; x++)
			m[y*w + x] = fabs(y-x) + 1;
}

void Matrix::GetAnswerVector(double* vector, const int size){
	int i = 0;
	for(; i < size-3; i+=4){
		vector[i] = 1;
		vector[i+1] = 0;
		vector[i+2] = 1;
		vector[i+3] = 0;
	}
	for(; i < size; i++)
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
		int x = 0;
		for(; x < size-3; x+=4){
			sum += fabs(matrix[x + y*size]);
			sum += fabs(matrix[x + 1 + y*size]);
			sum += fabs(matrix[x + 2 + y*size]);
			sum += fabs(matrix[x + 3 + y*size]);
		}
		for(; x < size; x++)
			sum += fabs(matrix[x + y*size]);
		if(y==0) max = sum;
		if(sum > max) max = sum;
	}

	return max;
}

double Matrix::LengthVector(const double* vector, const int size){
	double max = 0;
	int x = 0;
	//count four elements at once
	for(; x < size-3; x+=4){
		max += fabs(vector[x]);
		max += fabs(vector[x+1]);
		max += fabs(vector[x+2]);
		max += fabs(vector[x+3]);
	}
	//add the last elements
	for(; x < size; x++) max += fabs(vector[x]);

	return max;
}

void Matrix::GetRHSVector(const double* matrix, double* RHSVector, const int size){
	for(int y = 0; y < size; y++){
		RHSVector[y] = 0;
		int x = 0;
		for(; x < size-7; x+=8){
			RHSVector[y] += matrix[x + y*size];
			RHSVector[y] += matrix[x + 2 + y*size];
			RHSVector[y] += matrix[x + 4 + y*size];
			RHSVector[y] += matrix[x + 6 + y*size];
		}
		for(; x < size; x+=2)
			RHSVector[y] += matrix[x + y*size];
	}
}

void Matrix::GetBlock(const double* A, double* block, const int x, const int y, const int x1, const int y1, const int matrix_size){
	int h = y1-y;
	int w = x1-x;
    for(int i = 0; i < h; i++){
		int j = 0;
        for(; j < w-4; j+=4){
            block[j + i*w] = A[(y+i)*matrix_size + x + j];
            block[j + i*w + 1] = A[(y+i)*matrix_size + x + j + 1];
            block[j + i*w + 2] = A[(y+i)*matrix_size + x + j + 2];
            block[j + i*w + 3] = A[(y+i)*matrix_size + x + j + 3];
        }
		for(; j < w; j++)
			block[j + i*w] = A[(y+i)*matrix_size + x + j];
    }
}

void Matrix::PutBlock(double* A, const double* block, const int x, const int y, const int x1, const int y1, const int matrix_size){
    int h = y1-y;
	int w = x1-x;
    for(int i = 0; i < h; i++){
		int j = 0;
        for(; j < w-4; j+= 4){
            A[(y+i)*matrix_size + x + j] = block[j + i*w];
            A[(y+i)*matrix_size + x + j + 1] = block[j + 1 + i*w];
            A[(y+i)*matrix_size + x + j + 2] = block[j + 2 + i*w];
            A[(y+i)*matrix_size + x + j + 3] = block[j + 3 + i*w];
        }
		for(; j < w; j++)
			A[(y+i)*matrix_size + x + j] = block[j + i*w];
	}
}

void Matrix::NullMatrix(double* matrix, const int size){
	int x = 0;
    for(; x < size-3; x+=4){
        matrix[x] = 0;
		matrix[x+1] = 0;
		matrix[x+2] = 0;
		matrix[x+3] = 0;
	}
	for(; x < size; x++)
		matrix[x] = 0;
}

void Matrix::EMatrix(double* matrix, const int size){
	NullMatrix(matrix, size);
	int x = 0;
	for(; x < size-7; x+= 8){
		matrix[x*size + x] = 1;
		matrix[x*size + x+1] = 1;
		matrix[x*size + x+2] = 1;
		matrix[x*size + x+3] = 1;
		matrix[x*size + x+4] = 1;
		matrix[x*size + x+5] = 1;
		matrix[x*size + x+6] = 1;
		matrix[x*size + x+7] = 1;
	}
	for(; x < size; x++)
		matrix[x*size + x] = 1;
}

MatrixException Matrix::ReadMatrix(double* matrix, const int w, const int h, const char* file_name){
    FILE* f = fopen(file_name, "r");
    if(!f) return CAN_NOT_OPEN;
    
    for(int i = 0; i < w*h; i++){
        if(fscanf(f, "%lf", &matrix[i]) != 1){
            fclose(f);
            return FILE_CORRUPT;
        }
    }
    
    return NO_ERROR;
}

MatrixException Matrix::InitMatrix(double* matrix, const int w, const int h, const char* file_name){
    if(file_name == nullptr) {
        FillMatrix(matrix, w, h);
        return NO_ERROR;
    }
    else return ReadMatrix(matrix, w, h, file_name);
}

MatrixException Matrix::CreateVector(double** vector, const int size){
    *vector = new double[size];
    if(!*vector) return MatrixException::NOT_ENOUGH_MEMORY;
    
    NullMatrix(*vector, size);

    return NO_ERROR;
}

MatrixException Matrix::CreateMatrix(double** matrix, const int w, const int h){
    *matrix = new double[w*h];
    if(!*matrix) return MatrixException::NOT_ENOUGH_MEMORY;
    
    NullMatrix(*matrix, w*h);
    
    return NO_ERROR;
}

MatrixException Matrix::CreateMatrix(double** matrix, const int w, const int h, const char* file_name){
    //if no name was supplied, we want to initialize and generate the matrix
    MatrixException res = CreateMatrix(matrix, w, h);
	if(res == NO_ERROR){
		return InitMatrix(*matrix, w, h, file_name);
	}
	else return res;
}

int Matrix::GetInverseMatrix(double* matrix, double* inverse, int size){
	double temp_matrix[16];
	double t = 0;
	const double EPS = 1e-16;
	for(int x = 0; x < size; x++){
		int y = x;
		for(; y < size; y++){
			if(fabs(matrix[y*size + x]) > EPS){
				int n = 0;
				for(; n < size-3; n += 4){
					//swap elements
					temp_matrix[0] = matrix[x*size + n];
					matrix[x*size+n] = matrix[y*size+n];
					matrix[y*size+n] = temp_matrix[0];

					temp_matrix[1] = matrix[x*size + n+1];
					matrix[x*size+n+1] = matrix[y*size+n+1];
					matrix[y*size+n+1] = temp_matrix[1];

					temp_matrix[2] = matrix[x*size + n+2];
					matrix[x*size+n+2] = matrix[y*size+n+2];
					matrix[y*size+n+2] = temp_matrix[2];

					temp_matrix[3] = matrix[x*size + n+3];
					matrix[x*size+n+3] = matrix[y*size+n+3];
					matrix[y*size+n+3] = temp_matrix[3];
					
					temp_matrix[0] = inverse[x*size + n];
					inverse[x*size+n] = inverse[y*size+n];
					inverse[y*size+n] = temp_matrix[0];

					temp_matrix[1] = inverse[x*size + n+1];
					inverse[x*size+n+1] = inverse[y*size+n+1];
					inverse[y*size+n+1] = temp_matrix[1];

					temp_matrix[2] = inverse[x*size + n+2];
					inverse[x*size+n+2] = inverse[y*size+n+2];
					inverse[y*size+n+2] = temp_matrix[2];

					temp_matrix[3] = inverse[x*size + n+3];
					inverse[x*size+n+3] = inverse[y*size+n+3];
					inverse[y*size+n+3] = temp_matrix[3];
				}
				for(;n<size;n++){
					temp_matrix[0] = matrix[x*size + n];
					matrix[x*size+n] = matrix[y*size+n];
					matrix[y*size+n] = temp_matrix[0];

					temp_matrix[0] = inverse[x*size + n];
					inverse[x*size+n] = inverse[y*size+n];
					inverse[y*size+n] = temp_matrix[0];
				}
				break;
			}
		}
		if(y==size) return 0;

		t = matrix[x*size+x];
		int n = 0;
		for(; n < size-3; n += 4){
			matrix[x*size + n] /= t;
			matrix[x*size + n+1] /= t;
			matrix[x*size + n+2] /= t;
			matrix[x*size + n+3] /= t;

			inverse[x*size + n] /= t;
			inverse[x*size + n+1] /= t;
			inverse[x*size + n+2] /= t;
			inverse[x*size + n+3] /= t;
		}
		for(; n < size; n++){
			matrix[x*size + n] /= t;
			inverse[x*size + n] /= t;
		}

		for(n = 0; n < size; n++){
			t = matrix[n*size + x];
			int k = 0;
			for(; k < size-3; k+=4){
				if(n != x){
					inverse[n*size + k] -= inverse[x*size + k] * t;
					inverse[n*size + k+1] -= inverse[x*size + k+1] * t;
					inverse[n*size + k+2] -= inverse[x*size + k+2] * t;
					inverse[n*size + k+3] -= inverse[x*size + k+3] * t;

					matrix[n*size + k] -= matrix[x*size + k] * t;
					matrix[n*size + k+1] -= matrix[x*size + k+1] * t;
					matrix[n*size + k+2] -= matrix[x*size + k+2] * t;
					matrix[n*size + k+3] -= matrix[x*size + k+3] * t;
				}
			}
			for(; k < size; k++){
				if(n != x){
					inverse[n*size + k] -= inverse[x*size + k] * t;
					matrix[n*size + k] -= matrix[x*size + k] * t;
				}
			}
		}
	}
	return 1;
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

void Matrix::SolveBlock(double* matrix, double* rhs, double* answer, const int size, const int block_size){
    int step = size/block_size;
	int end = size - step*block_size;

	int offset = 0;
    
    double* block = nullptr;
    
	//create a permutation, so we dont take time to actually move elements in the matrix
	auto indexes = new int[size];
	for(int i = 0; i < size; i++)
		indexes[i] = i;

	//when we complete a step of Gaussian algorithm, we should apply it again to the matrix of size m-1. 
	//For that we will just think of the next element on the diagonal as the first one.
	while(offset != step){
		//we want to find the element with lowest inverse length
		int minimal_length_index = offset;
        GetBlock(matrix, block, offset*block_size, offset*block_size, size, block_size);
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

		
		double minimal_norm = 10e300;
		double minimal_norm_index = offset;
		for(int y = offset; y < step; y++){
			GetBlock(matrix, block, offset*block_size, offset*block_size, size, block_size);
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


























