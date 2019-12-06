#include "matrix.h"
#include "matrix_exception.h"

#include <cmath>
#include <string.h>

void Matrix::FillMatrix(double* m, const int w, const int h){
	for(int y = 0; y < h; y++)
		for(int x = 0; x < w; x++)
			m[y*w + x] = fabs(x-y);
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

void Matrix::PrintClean(const double* matrix, const int w, const int h){
    for(int y = 0; y < h; y++){
        for(int x = 0; x < w; x++)
            printf("%lf ", matrix[x + y*w]);
        printf("\n");
    }
    printf("\n");
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

double Matrix::Length(const double* matrix, const int w, const int h){
	double max, sum;
	for(int y = 0; y < h; y++){
		sum = 0;
		int x = 0;
		for(; x < w-3; x+=4){
			sum += fabs(matrix[x + y*w]);
			sum += fabs(matrix[x + 1 + y*w]);
			sum += fabs(matrix[x + 2 + y*w]);
			sum += fabs(matrix[x + 3 + y*w]);
		}
		for(; x < w; x++)
			sum += fabs(matrix[x + y*w]);
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
	NullMatrix(matrix, size*size);
	int x = 0;
	for(; x < size-7; x+= 8){
		matrix[x*size + x] = 1;
		matrix[(x+1)*size + x+1] = 1;
		matrix[(x+2)*size + x+2] = 1;
		matrix[(x+3)*size + x+3] = 1;
		matrix[(x+4)*size + x+4] = 1;
		matrix[(x+5)*size + x+5] = 1;
		matrix[(x+6)*size + x+6] = 1;
		matrix[(x+7)*size + x+7] = 1;
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

int Matrix::GetInverseMatrix(double* matrix, double* matrixReversed, int m, double norm, int* transposition_m){
	double temp = 0;
  	double temp_m[4];
	const double EPS = 1e-16;

  	for (int i = 0; i < m; i++)
    	transposition_m[i] = i;

  	for (int k = 0; k < m; k++){
      	double max_norm = -10e300;
      	int i_temp = k, j_temp = k;
      	for (int i = k; i < m; i++){
          	for (int j = k; j < m; j++){
              	temp = matrix[i*m + j];
              	if (max_norm < fabs (temp)){
                  	max_norm = fabs (temp);
                  	i_temp = i;
                  	j_temp = j;
                }
            }
        }
      	if (fabs (max_norm) < EPS * norm)
         	return 0;

      	if (i_temp != k){
          	int j = 0;
          	for (; j < m - 3; j += 4){
		        temp_m[0] = matrix[k*m+j];
		        matrix[k*m+j] = matrix[i_temp*m+j];
		        matrix[i_temp*m+j] = temp_m[0];

				temp_m[1] = matrix[k*m+j + 1];
		        matrix[k*m+j+1] = matrix[i_temp*m+j+1];
		        matrix[i_temp*m+j+1] = temp_m[1];

				temp_m[2] = matrix[k*m+j + 2];
		        matrix[k*m+j+2] = matrix[i_temp*m+j+2];
		        matrix[i_temp*m+j+2] = temp_m[2];

				temp_m[3] = matrix[k*m+j + 3];
		        matrix[k*m+j+3] = matrix[i_temp*m+j+3];
		        matrix[i_temp*m+j+3] = temp_m[3];

				temp_m[0] = matrixReversed[k*m+j];
		        matrixReversed[k*m+j] = matrixReversed[i_temp*m+j];
		        matrixReversed[i_temp*m+j] = temp_m[0];

		        temp_m[1] = matrixReversed[k*m+j + 1];
		        matrixReversed[k*m+j+1] = matrixReversed[i_temp*m+j+1];
		        matrixReversed[i_temp*m+j+1] = temp_m[1];

				temp_m[2] = matrixReversed[k*m+j + 2];
		        matrixReversed[k*m+j+2] = matrixReversed[i_temp*m+j+2];
		        matrixReversed[i_temp*m+j+2] = temp_m[2];

				temp_m[3] = matrixReversed[k*m+j + 3];
		        matrixReversed[k*m+j+3] = matrixReversed[i_temp*m+j+3];
		        matrixReversed[i_temp*m+j+3] = temp_m[3];;
            }
          	for (; j < m; j++){
		        temp = matrix[k*m+j];
		        matrix[k*m+j] = matrix[i_temp*m+j];
		        matrix[i_temp*m+j] = temp;
		        temp = matrixReversed[k*m+j];
		        matrixReversed[k*m+j] = matrixReversed[i_temp*m+j];
		        matrixReversed[i_temp*m+j] = temp;
            }

        }

		if (j_temp != k){
          	for (int i = 0; i < m; i++){
		          temp = matrix[i*m+k];
		          matrix[i*m+k] = matrix[i*m+j_temp];
		          matrix[i*m+j_temp] = temp;
            }
          	transposition_m[k] = j_temp;
        }

		temp = matrix[k*m+k];
		int l = 0;
		for (; l < m - 3; l += 4){
			matrix[k*m+l] /= temp;
			matrix[k*m+l+1] /= temp;
			matrix[k*m+l+2] /= temp;
			matrix[k*m+l+3] /= temp;

			matrixReversed[k*m+l] /= temp;
			matrixReversed[k*m+l+1] /= temp;
			matrixReversed[k*m+l+2] /= temp;
			matrixReversed[k*m+l+3] /= temp;
        }
      	for (; l < m; l++){
          	matrix[k*m+l] /= temp;
          	matrixReversed[k*m+l] /= temp;
        }
      	for (l = 0; l < m; l++){
			temp = matrix[l*m + k];
		    int p = 0;
		    for (; p < m - 3; p += 4){
            	if (l != k){
					matrixReversed[l*m + p] -= matrixReversed[k*m + p]*temp;
					matrixReversed[l*m + p+1] -= matrixReversed[k*m + p+1]*temp;
					matrixReversed[l*m + p+2] -= matrixReversed[k*m + p+2]*temp;
					matrixReversed[l*m + p+3] -= matrixReversed[k*m + p+3]*temp;

					matrix[l*m + p] -= matrix[k*m + p]*temp;
					matrix[l*m + p+1] -= matrix[k*m + p+1]*temp;
					matrix[l*m + p+2] -= matrix[k*m + p+2]*temp;
					matrix[l*m + p+3] -= matrix[k*m + p+3]*temp;
                }
            }
          	for (; p < m; p++){
            	if (l != k){
		       		matrixReversed[l*m + p] -= matrixReversed[k*m + p]*temp;
		            matrix[l*m + p] -= matrix[k*m + p]*temp;
                }
            }
        }
    }

  	for (int i = m - 1; i >= 0; i--){
      	if (i != transposition_m[i]){
          	for (int j = 0; j < m; j++){
              	temp = matrixReversed[i * m + j];
             	matrixReversed[i * m + j] = matrixReversed[transposition_m[i] * m + j];
             	matrixReversed[transposition_m[i] * m +j] = temp;
            }
        }
    }

  return 1;
}

double* Matrix::MultiplyMatrices(double* matrix1, double* matrix2, double* res, const int h, const int wh, const int w){
	if(h == wh && wh == w){
		if(h % 4 == 0){
			MultiplyMatricesNN_4(matrix1, matrix2, res, h);
			return res;
		}
		else{
			MultiplyMatricesNN(matrix1, matrix2, res, h);
			return res;
		}
	}
	else
		MultiplyMatricesNMK(matrix1, matrix2, res, h, wh, w);

	return res;
}

double* Matrix::MultiplyMatricesNMK(double* matrix1, double* matrix2, double* res, const int h, const int wh, const int w){
	int q = 0;
	double sum = 0;
	for(int y = 0; y < h; y++){
		for(int x = 0; x < w; x++){
			sum = 0;
			for(q = 0; q < wh-3; q+=4){
				sum += matrix1[y*wh + q]*matrix2[q*w + x];
				sum += matrix1[y*wh + q+1]*matrix2[(q+1)*w + x];
				sum += matrix1[y*wh + q+2]*matrix2[(q+2)*w + x];
				sum += matrix1[y*wh + q+3]*matrix2[(q+3)*w + x];
			}
			for(;q < wh; q++){
				sum += matrix1[y*wh + q]*matrix2[q*w + x];
			}
			res[y*w + x] = sum;
		}
	}
	return res;
}

void Matrix::MultiplyMatricesNN_4(double* a, double* b, double* c, const int m){
	const double * pa, * pb;
  double * pc;
  double s00, s01, s02, s03;
  double s10, s11, s12, s13;
  double s20, s21, s22, s23;
  double s30, s31, s32, s33;
  double a0, a1, a2, a3;
  double b0, b1, b2, b3;
  int i, j, k;
  bzero (c, m * m * sizeof (double));
  // matrix_0 (c, m);
  for (i = 0; i < m; i += 4)
    {
      for (j = 0; j < m; j += 4)
        {
          s00 = 0;
          s01 = 0;
          s02 = 0;
          s03 = 0;

          s10 = 0;
          s11 = 0;
          s12 = 0;
          s13 = 0;

          s20 = 0;
          s21 = 0;
          s22 = 0;
          s23 = 0;

          s30 = 0;
          s31 = 0;
          s32 = 0;
          s33 = 0;

          pa = a + i * m;
          pb = b + j;

          for (k = 0; k < m; k++)
            {
              a0 = pa[0];
              a1 = pa[m];
              a2 = pa[2 * m];
              a3 = pa[3 * m];

              b0 = pb[0];
              b1 = pb[1];
              b2 = pb[2];
              b3 = pb[3];

              s00 += a0 * b0;
              s01 += a0 * b1;
              s02 += a0 * b2;
              s03 += a0 * b3;

              s10 += a1 * b0;
              s11 += a1 * b1;
              s12 += a1 * b2;
              s13 += a1 * b3;

              s20 += a2 * b0;
              s21 += a2 * b1;
              s22 += a2 * b2;
              s23 += a2 * b3;

              s30 += a3 * b0;
              s31 += a3 * b1;
              s32 += a3 * b2;
              s33 += a3 * b3;

              pa += 1;
              pb += m;
            }
          pc = c + i * m + j;

          pc[0] += s00;
          pc[1] += s01;
          pc[2] += s02;
          pc[3] += s03;

          pc[m] += s10;
          pc[m + 1] += s11;
          pc[m + 2] += s12;
          pc[m + 3] += s13;

          pc[2 * m] += s20;
          pc[2 * m + 1] += s21;
          pc[2 * m + 2] += s22;
          pc[2 * m + 3] += s23;

          pc[3 * m] += s30;
          pc[3 * m + 1] += s31;
          pc[3 * m + 2] += s32;
          pc[3 * m + 3] += s33;
        }
    }
}

void Matrix::MultiplyMatricesNN(double* a, double* b, double* c, const int n){
	int bm, bi, nbm, nbi, N = 16;
  int l, nl;
  int i, j, m;
  double *pa, *pb, *pc, s;
  for (bm = 0; bm < n; bm += N)
    {
      nbm = (bm + N <= n ? bm + N : n);
      for (bi = 0; bi < n; bi += N)
        {
          nbi = (bi + N <= n ? bi + N : n);
          for (m = bm, pc = c + bm; m < nbm; m++, pc++)
            for (i = bi; i < nbi; i++)
              pc[i * n] = 0.;
          for (l = 0; l < n; l += N)
            {
              nl = (l + N <= n ? l + N : n);
              for (m = bm, pc = c+bm; m < nbm; m++, pc++)
                for (i = bi, pb = b + m; i < nbi; i++)
                  {
                    pa = a + l + i * n;
                    for (s = 0., j = l; j < nl; j++)
                      s += *(pa++) * pb[j * n];
                    pc[i * n] += s;
                  }
            }
        }
    }
}

double* Matrix::SubstractMatrices(double* matrix, const double* matrix2, const int w, const int h){
	for(int y = 0; y < h; y++){
		int x = 0;
		for(; x < w-3; x+=4){
			matrix[y*w + x] -= matrix2[y*w+x];
			matrix[y*w + x+1] -= matrix2[y*w+x+1];
			matrix[y*w + x+2] -= matrix2[y*w+x+2];
			matrix[y*w + x+3] -= matrix2[y*w+x+3];
		}
		for(;x<w;x++)
			matrix[y*w + x] -= matrix2[y*w+x];
	}

	return matrix;
}


double Matrix::GetError(double* vector, const int size){
    double max=fabs(vector[0]-1);
	for(int x = 0; x < size; x++){
		if(fabs(vector[x] - ((x+1)&1)) > max) max = fabs(vector[x] - ((x+1)&1));
	}

	return max;
}

void Matrix::SolveBlock(double* matrix, double* rhs, double* answer, const int size, const int block_size, double* temps[17]){
    int step = size/block_size;
	int end = size - step*block_size;

	int offset = 0;
    
    double* block = temps[0];
    double* block_temp = temps[1];
    double* block_temp_im = temps[2];
    double* block_temp_sub = temps[3];
	block_temp_sub += 1;
    block_temp_sub -= 1;
    double* block_me = temps[4];
    double* block_me_temp = temps[5];
    double* block_me_temp_im = temps[6];
    double* block_me_temp_sub = temps[7];
	block_me_temp_sub += 1;
    block_me_temp_sub -= 1;
	double* inverse_block = temps[8];
    double* vector_block = temps[9];
    double* vector_block_temp = temps[10];
    double* vector_block_temp_im = temps[11];
    double* vector_block_temp_sub = temps[12];
	vector_block_temp_sub += 1;
    vector_block_temp_sub -= 1;
	double* block_ee = temps[13];
	double* block_ee_temp = temps[14];
	double* vector_e = temps[15];
	double* vector_e_temp = temps[16];
	
    
	//create a permutation, so we dont take time to actually move elements in the matrix
	auto indexes = new int[step];
	auto indexes_m = new int[size];
	for(int i = 0; i < step; i++)
		indexes[i] = i;

	double n = Length(matrix, size, size);
    
    //PrintClean(matrix, size, size);
    //PrintClean(rhs, 1, size);

	//when we complete a step of Gaussian algorithm, we should apply it again to the matrix of size m-1. 
	//For that we will just think of the next element on the diagonal as the first one.
	while(offset != step){
		double minimal_norm = 10e300;
		int minimal_norm_index = offset;
        int found_inversable = 0;
		for(int y = offset; y < step; y++){
			GetBlock(matrix, block, offset*block_size, indexes[y]*block_size, 
									offset*block_size + block_size, indexes[y]*block_size + block_size, size);
			EMatrix(inverse_block, block_size);
			if(GetInverseMatrix(block, inverse_block, block_size, n, indexes_m)){
                found_inversable = 1;
                double k = Length(inverse_block, block_size, block_size);
                if(minimal_norm > k){
                    minimal_norm = k;
                    minimal_norm_index = y;
                }
            }
		}
		if(found_inversable == 0){
            printf("No matrices can be inverted on column %d\n", offset);
            return;
        }
		
		//swap the line to the top
		int temp_index = indexes[offset];
		indexes[offset] = indexes[minimal_norm_index];
		indexes[minimal_norm_index] = temp_index;
        
		//get inverse block
        GetBlock(matrix, block, offset*block_size, indexes[offset]*block_size, 
								offset*block_size + block_size, indexes[offset]*block_size + block_size, size);
        EMatrix(inverse_block, block_size);
        //we know this matrix exists, no need to check it
        GetInverseMatrix(block, inverse_block, block_size, n, indexes_m);


        //firstly normalize the rhs vector
        GetBlock(rhs, vector_block, 0, indexes[offset]*block_size, 
									1, indexes[offset]*block_size + block_size, 1);
		MultiplyMatrices(inverse_block, vector_block, vector_block_temp, block_size, block_size, 1);
        PutBlock(rhs, vector_block_temp, 0, indexes[offset]*block_size, 
										 1, indexes[offset]*block_size + block_size, 1);		 

		//if unfull end block exists norm it too
		if(end > 0){
			GetBlock(matrix, block_me, step*block_size, indexes[offset]*block_size, 
									   step*block_size + end, indexes[offset]*block_size + block_size, size);
			MultiplyMatrices(inverse_block, block_me, block_me_temp, block_size, block_size, end);
            PutBlock(matrix, block_me_temp, step*block_size, indexes[offset]*block_size, 
											step*block_size + end, indexes[offset]*block_size + block_size, size);
		}
		
        //normalize the first row of the matrix
        for(int x = offset+1; x < step; x++){
            //normalize current block
            GetBlock(matrix, block, x*block_size, indexes[offset]*block_size, 
									x*block_size + block_size, indexes[offset]*block_size + block_size, size);
			MultiplyMatrices(inverse_block, block, block_temp, block_size, block_size, block_size);
            PutBlock(matrix, block_temp, x*block_size, indexes[offset]*block_size, 
										 x*block_size + block_size, indexes[offset]*block_size + block_size, size);
        }
        
		//everything normalized, block_me_temp has ready me normalized block of the top row
		//vector_block_temp has the same for rhs vector

        //substract top line of blocks from the all the bottom ones of block_size
        for(int y = offset+1; y < step; y++){
            //first element of the current row   
            GetBlock(matrix, block, offset*block_size, indexes[y]*block_size, 
									offset*block_size + block_size, indexes[y]*block_size + block_size, size);


			MultiplyMatrices(block, vector_block_temp, vector_block_temp_im, block_size, block_size, 1);
			//rhs block
            GetBlock(rhs, vector_block, 0, indexes[y]*block_size, 
										1, indexes[y]*block_size + block_size, 1);
			SubstractMatrices(vector_block, vector_block_temp_im, 1, block_size);
            PutBlock(rhs, vector_block, 0, indexes[y]*block_size, 
										1, indexes[y]*block_size + block_size, 1);
            
            //end block if it exists
            if(end > 0){
				MultiplyMatrices(block, block_me_temp, block_me_temp_im, block_size, block_size, end);
                GetBlock(matrix, block_me, step*block_size, indexes[y]*block_size, 
										   step*block_size + end, indexes[y]*block_size + block_size, size);
				SubstractMatrices(block_me, block_me_temp_im, block_size, end);
                PutBlock(matrix, block_me, step*block_size, indexes[y]*block_size, 
										   step*block_size + end, indexes[y]*block_size + block_size, size);
            }
            
            //substract all other blocks
            for(int x = offset+1; x < step; x++){
                GetBlock(matrix, block_temp, x*block_size, indexes[offset]*block_size, 
											 x*block_size + block_size, indexes[offset]*block_size + block_size, size);
				MultiplyMatrices(block, block_temp, block_temp_im, block_size, block_size, block_size);
                GetBlock(matrix, block_temp, x*block_size, indexes[y]*block_size, 
											 x*block_size + block_size, indexes[y]*block_size + block_size, size);
				SubstractMatrices(block_temp, block_temp_im, block_size, block_size);
                PutBlock(matrix, block_temp, x*block_size, indexes[y]*block_size, 
											 x*block_size + block_size, indexes[y]*block_size + block_size, size);
            }
        }

		//if the last row exists, substract from them too
		if(end > 0){
			//first element of the last end row   
            GetBlock(matrix, block_me, offset*block_size, step*block_size, 
									   offset*block_size + block_size, step*block_size + end, size);
            
            //rhs vector end block
            GetBlock(rhs, vector_e, 0, step*block_size, 
									1, step*block_size + end, 1);
			MultiplyMatrices(block_me, vector_block_temp, vector_e_temp, end, block_size, 1);
			SubstractMatrices(vector_e, vector_e_temp, 1, end);
            PutBlock(rhs, vector_e, 0, step*block_size, 
									1, step*block_size + end, 1);

            //ee block
            GetBlock(matrix, block_ee, step*block_size, step*block_size, 
									   step*block_size + end, step*block_size + end, size);
			MultiplyMatrices(block_me, block_me_temp, block_ee_temp, end, block_size, end);
			SubstractMatrices(block_ee, block_ee_temp, end, end);
            PutBlock(matrix, block_ee, step*block_size, step*block_size, 
									   step*block_size + end, step*block_size + end, size);
            
            //substract all blocks
            for(int x = offset+1; x < step ; x++){
                GetBlock(matrix, block_temp, x*block_size, indexes[offset]*block_size, 
											 x*block_size + block_size, indexes[offset]*block_size + block_size, size);
				MultiplyMatrices(block_me, block_temp, block_me_temp_im, end, block_size, block_size);
                GetBlock(matrix, block_me_temp_sub, x*block_size, step*block_size, 
											 		x*block_size + block_size, step*block_size + end, size);
				SubstractMatrices(block_me_temp_sub, block_me_temp_im, block_size, end);
                PutBlock(matrix, block_me_temp_sub, x*block_size, step*block_size, 
											 		x*block_size + block_size, step*block_size + end, size);
            }
		}
		
		offset += 1;
	}

	if(end > 0){
		GetBlock(matrix, block_ee, step*block_size, step*block_size, 
									   step*block_size + end, step*block_size + end, size);
		EMatrix(block_ee_temp, end);
		if(!GetInverseMatrix(block_ee, block_ee_temp, end, n, indexes_m)){
        	printf("End part of the matrix could not be inverted\n");
            return;
        }

        GetBlock(rhs, vector_e, 0, step*block_size, 
                                1, step*block_size + end, 1);
		MultiplyMatrices(block_ee_temp, vector_e, vector_e_temp, end, end, 1);

		PutBlock(rhs, vector_e_temp, 0, step*block_size, 
                                     1, step*block_size + end, 1);

        for(int y = 0; y < step; y++){
            GetBlock(matrix, block_me, step*block_size, indexes[y]*block_size, 
									step*block_size + end, indexes[y]*block_size + block_size, size);
            MultiplyMatrices(block_me, vector_e_temp, vector_block, block_size, end, 1);
            GetBlock(rhs, vector_block_temp, 0, indexes[y]*block_size,
                                            1, indexes[y]*block_size + block_size, 1);
            
            SubstractMatrices(vector_block_temp, vector_block, 1, block_size);
            PutBlock(rhs, vector_block_temp, 0, indexes[y]*block_size,
                                            1, indexes[y]*block_size + block_size, 1);
        }
    } 

	
	//backwards step
	for (int x = step-1; x > 0; x--){
     	GetBlock(rhs, vector_block, 0, indexes[x] * block_size, 
								    1, indexes[x] * block_size + block_size, 1);
      	for (int y = 0; y < x; y++){
          	GetBlock (matrix, block_temp, x * block_size, indexes[y] * block_size, 
										x * block_size + block_size, indexes[y] * block_size + block_size, size);
          	//multiplication block_temp by vector_block put into vector_block_temp
			MultiplyMatrices(block_temp, vector_block, vector_block_temp, block_size, block_size, 1);
		
          	GetBlock (rhs, vector_block, 0, indexes[y] * block_size, 
										 1, indexes[y] * block_size + block_size, 1);
          	//substract vector_block and vector_block_temp put into vector_block_temp_sub
			SubstractMatrices(vector_block, vector_block_temp, 1, block_size);
          	PutBlock (rhs, vector_block, 0, indexes[y] * block_size, 
										 1, indexes[y] * block_size + block_size, 1);
        }
    }
	
	//swap the answer vector back
	for(int y = 0; y < size; y++){
        if(y < step*block_size)
            answer[y] = rhs[indexes[y/block_size]*block_size + y%block_size];
        else
            answer[y] = rhs[y];
    }
	
	
	delete[] indexes;
	delete[] indexes_m;
}


























