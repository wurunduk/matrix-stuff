#ifndef _MATRIX_H_
#define _MATRIX_H_

#include "matrix_exception.h"
#include <memory>

class Matrix{
    public:
    std::unique_ptr<double[]> matrix;
    int width;
    int height;

    Matrix();

	void FillMatrix();	
    void Print(int size_ = 10) const;
	double Length() const;

	Matrix GetRHSVector();

	Matrix Solve(const Matrix* rhs);

    MatrixException CreateMatrix(int width_, int height_);
    MatrixException CreateMatrix(int width_, int height_, const char* file_name);

	Matrix operator-(const Matrix& m);

    Matrix operator*(const Matrix& m);
    Matrix operator*(const double& k);
    Matrix& operator*=(const double& k);
};

#endif
