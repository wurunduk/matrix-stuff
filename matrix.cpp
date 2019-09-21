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
    if(size_ > width) size_ = width;
    if(size_ > height) size_ = height;

    for(int y = 0; y < height; y++){
        for(int x = 0; x < width; x++)
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

Matrix Matrix::Solve(const Matrix& rhs){

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
































