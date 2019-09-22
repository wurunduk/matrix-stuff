#include "matrix.h"
#include "matrix_exception.h"

Matrix::Matrix(){
    width = 0;
    height = 0;
}

void Matrix::FillMatrix(){
	for(int y = 0; y < height; y++)
		for(int x = 0; x < width; x++)
			matrix[x + y*width] = 1.0/(x+y+1.0);
}

void Matrix::Print(int size_) const{
	int sizeX = size_;
	int sizeY = size_;
    if(sizeX > width) sizeX = width;
    if(sizeY > height) sizeY = height;

    for(int y = 0; y < sizeY; y++){
        for(int x = 0; x < sizeX; x++)
            printf("%lf ", matrix[x + y*width]);
        printf("\n");
    }
}

double Matrix::Length() const{
	double max, sum;
	for(int y = 0; y < width; y++){
		sum = 0;
		for(int x = 0; x < height; x++){
			sum += matrix[x + y*width];
		}
		if(y==0) max = sum;
		if(sum > max) max = sum;
	}

	return max;
}


Matrix Matrix::GetRHSVector(){
	Matrix result;
	result.CreateMatrix(1, height);

	for(int y = 0; y < height; y++){
		double sum = 0;
		for(int x = 0; x < width; x+=2){
			sum += matrix[x + y*width];
		}
		result.matrix[y] = sum;
	}
	return result;
}

Matrix Matrix::Solve(const Matrix* rhs){
	int offset = 0;
	//create a permutation, so we dont take time to actually move elements in the matrix
	auto indexes = new int[height];
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
		for(int x = offset; x < width; x++){
			matrix[x+indexes[offset]*width] /= max_length;
			rhs->matrix[indexes[offset]] /= max_length;
		}

		//substract the top line from all the lines below it
		for(int y = offset+1; y < height; y++){
			rhs->matrix[indexes[y]] -= rhs->matrix[indexes[offset]]*matrix[offset+indexes[y]*width];
			for(int x = width-1; x >= offset; x--){
				matrix[x+indexes[y]*width] -= matrix[x+indexes[offset]*width]*matrix[offset+indexes[y]*width];
			}
		}	

		printf("Step %d\n", offset);
		Print(12);

		//repeat for the sub matrix
		offset += 1;
	}

	Matrix answer;
	answer.CreateMatrix(1, height);

	//at this point we have an upper diagonal matrix and we can get the answer
	for(int y = height-1; y >= 0; y--){
		double sum = 0.0;
		//technically this is an x coordinate
		for(int i = width-1; i > indexes[y]; i--)
			sum += answer.matrix[i]*matrix[indexes[y]*width + i];

		answer.matrix[y] = rhs->matrix[indexes[y]] - sum;
	}
	
	delete[] indexes;

	return answer;
}

MatrixException Matrix::CreateMatrix(int width_, int height_){
    width = width_;
    height = height_;
    double* temp = new double[height_*width_];
    if(!temp) return MatrixException::NOT_ENOUGH_MEMORY;
    
    matrix = std::unique_ptr<double[]>(temp);
    return NO_ERROR;
}

MatrixException Matrix::CreateMatrix(int width_, int height_, const char* file_name){
    //if no name was supplied, we want to initialize and generate the matrix
    MatrixException res = CreateMatrix(width_, height_);
    if(file_name == nullptr){
		if(res == NO_ERROR) FillMatrix();
        return res;
    }
	//initialize matrix if we are reading from file
    else if((res = CreateMatrix(width_, height_)) && res != NO_ERROR) return res;

	if(res == NO_ERROR){
		if(file_name == nullptr) {
			FillMatrix();
			return NO_ERROR;
		}
	}
	else return res;

    FILE* f = fopen(file_name, "r");
    if(!f) return CAN_NOT_OPEN;
    
    for(int i = 0; i < width_*height_; i++){
        if(fscanf(f, "%lf", &matrix[i]) != 1){
            fclose(f);
            return FILE_CORRUPT;
        }
    }

    return NO_ERROR;
}

Matrix Matrix::operator-(const Matrix& m){
	Matrix result;
	result.CreateMatrix(m.width, m.height);

	if(m.height != this->height || m.width != this->width){
		printf("Can not substract matrices of different sizes %dx%d - %dx%d.\n", this->width, this->height, m.width, m.height);
		return result;
	}

	for(int i = 0; i < m.height*m.width; i++)
			result.matrix[i] = this->matrix[i] - m.matrix[i];

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
            for(int i = 0; i < this->width/*we know m.height=this.width at this moment*/; i++)
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
}
































