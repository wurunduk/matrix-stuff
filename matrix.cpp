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
	double sum=0;
	for(int x = 0; x < size; x++){
		sum += fabs(vector[x]);
	}

	return sum;
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
/*
void Matrix::Solve(const Matrix* rhs, Matrix* answer){
	int offset = 0;
	//create a permutation, so we dont take time to actually move elements in the matrix
	int* indexes = reinterpret_cast<int*>(malloc(height*sizeof(int)));
	for(int i = 0; i < height; i++)
		indexes[i] = i;

	//when we complete a step of Gaussian algorithm, we should apply it again to the matrix of size m-1. 
	//For that we will just think of the next element on the diagonal as the first one.
	while(offset != width){
		//we want to find the element with lowest inverse length
		//with numbers it will be just 1/k, so we can search for the largest element in the column
		int minimal_length_index = offset;
		double max_length = matrix[offset + indexes[offset]*width];
		for(int y = offset; y < height; y++){
			double k = matrix[offset+indexes[y]*width];
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
		for(int x = offset; x < width; x++){
			matrix[x+indexes[offset]*width] *= k;
		}
		rhs->matrix[indexes[offset]] *= k;

		//substract the top line from all the lines below it
		for(int y = offset+1; y < height; y++){
			rhs->matrix[indexes[y]] -= rhs->matrix[indexes[offset]]*matrix[offset+indexes[y]*width];
			for(int x = width-1; x >= offset; x--){
				matrix[x+indexes[y]*width] -= matrix[x+indexes[offset]*width]*matrix[offset+indexes[y]*width];
			}
		}

		//repeat for the sub matrix
		offset += 1;
	}

	//at this point we have an upper diagonal matrix and we can get the answer
	for(int y = height-1; y >= 0; y--){
		double sum = 0.0;
		//technically this is an x coordinate
		for(int i = width-1; i > y; i--){
			sum += answer->matrix[i]*matrix[indexes[y]*width + i];
		}

		answer->matrix[y] = rhs->matrix[indexes[y]] - sum;
	}
	
	free(indexes);
}*/

MatrixException Matrix::CreateVector(double** vector, const int size){
    *vector = new double[size];
    if(!*vector) return MatrixException::NOT_ENOUGH_MEMORY;

    return NO_ERROR;
}

MatrixException Matrix::CreateMatrix(double** matrix, const int size){
    *matrix = new double[size*size];
    if(!*matrix) return MatrixException::NOT_ENOUGH_MEMORY;

    return NO_ERROR;
}

MatrixException Matrix::CreateMatrix(double** matrix, const int size, const char* file_name){
    //if no name was supplied, we want to initialize and generate the matrix
    MatrixException res = CreateMatrix(matrix, size);
	if(res == NO_ERROR){
		if(file_name == nullptr) {
			FillMatrix(*matrix, size);
			return NO_ERROR;
		}
	}
	else return res;

    FILE* f = fopen(file_name, "r");
    if(!f) return CAN_NOT_OPEN;
    
    for(int i = 0; i < size*size; i++){
        if(fscanf(f, "%lf", &(*matrix)[i]) != 1){
            fclose(f);
            return FILE_CORRUPT;
        }
    }

    return NO_ERROR;
}

double* Matrix::MultiplyMatrixByVector(const double* matrix, const double* vector, double* answer, const int size){    
    for(int y = 0; y < size; y++){
		double sum = 0;
        for(int x = 0; x < size; x++){
            sum += matrix[y*size + x]*vector[x];
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

/*Matrix Matrix::operator-(const Matrix& m){
	Matrix result;
	result.CreateMatrix(m.width, m.height);

	if(m.height != this->height || m.width != this->width){
		printf("Can not substract matrices of different sizes %dx%d - %dx%d.\n", this->width, this->height, m.width, m.height);
		return result;
	}

	for(int i = 0; i < m.height*m.width; i++){
			result.matrix[i] = this->matrix[i] - m.matrix[i];
	}

	return result;
}

Matrix Matrix::operator*(const Matrix& m){

    Matrix result;
    result.CreateMatrix(m.width, this->height);

    if(m.height != this->width){
        printf("Can not multiply matrix of width %d by matrix of height %d", this->width, m.height);
        return result;
    }
    
    for(int y = 0; y < this->height; y++){
        for(int x = 0; x < m.width; x++){
            double sum = 0;
            //i - x coordinate in this matrix, y coordinate in matrix m
            for(int i = 0; i < this->width; i++)
                sum += this->matrix[y*this->width + i]*m.matrix[x + i*m.width];
            result.matrix[x + y*result.width] = sum;
        }
    }
    
    return result;
}

Matrix Matrix::operator*(const double& k){
    Matrix result;
    result.CreateMatrix(this->width, this->height);
    for(int i = 0; i < this->width*this->height; i++)
        result.matrix[i] = this->matrix[i]*k;
    return result;
}

Matrix& Matrix::operator*=(const double& k){
    for(int i = 0; i < this->width*this->height; i++)
        matrix[i]*=k;
    return *this;
}*/
































