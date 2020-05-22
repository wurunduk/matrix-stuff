#include "matrix.h"
#include "matrix_exception.h"
#include <sys/sysinfo.h>

#include <cmath>
#include <string.h>

void Matrix::InitializeTempAddresses(address_array* adr, int size, int block_size, int end){

	//memset addresses
	adr->inverse_block = new double[block_size*block_size*5];
	adr->block = adr->inverse_block + block_size*block_size;
	adr->block_temp = adr->block + block_size*block_size;
	adr->block_temp_im = adr->block_temp + block_size*block_size;
	adr->block_temp_sub = adr->block_temp_im + block_size*block_size;
/*if(id == 0){
		memset(in->matrix + step*m*n, 0, end*n*sizeof(double));
		memset(in->rhs + n - end, 0, end*sizeof(double));
		memset(in->answer + n - end, 0, end*sizeof(double));
	}

	for(int i = 0; (i + id)*m < n; i += in->thread_count){
		if((i + id)*m + m < n){
			memset(in->matrix + (i + id)*m*n, 0, m*n*sizeof(double));
			memset(in->rhs + (i+id)*m, 0, m*sizeof(double));
			memset(in->answer + (i+id)*m, 0, m*sizeof(double));
		}
	}*/


	adr->temp_row = new double[block_size*size];
	
	adr->block_me = new double[block_size*end*4];
	adr->block_me_temp = adr->block_me + block_size*end;
	adr->block_me_temp_im = adr->block_me_temp + block_size*end;
	adr->block_me_temp_sub = adr->block_me_temp_im + block_size*end;
	
	adr->vector_block = new double[block_size*3];
	adr->vector_block_temp = adr->vector_block + block_size;
	adr->vector_block_temp_im = adr->vector_block_temp + block_size;
	
	adr->block_ee = new double[end*end*2];
	adr->block_ee_temp = adr->block_ee + end*end;
	
	adr->vector_e = new double[end*2];
	adr->vector_e_temp = adr->vector_e + end;
}

void Matrix::DeleteTempAddresses(address_array* adr){
	delete[] adr->inverse_block;
	delete[] adr->temp_row;
	delete[] adr->block_me;
	delete[] adr->vector_block;
	delete[] adr->block_ee;
	delete[] adr->vector_e;
}

MatrixException Matrix::AttachMatrices(arg* in){
	int n = in->size;
	
	auto e = in->file_name ? ReadMatrix(in->matrix, n, n, in->file_name) : FillMatrix(in->matrix, n, n);
	if(e == NO_ERROR) GetRHSVector(in->matrix, in->rhs, n);
	return e;
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
	double max = 0.0, sum;
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
	if(size > 0) max = fabs(vector[x]);
	//add the last elements
	for(; x < size; x++) 
		if(fabs(vector[x]) > max) max = fabs(vector[x]);

	return max;
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

MatrixException Matrix::FillMatrix(double* m, const int w, const int h){
	for(int y = 0; y < h; y++)
		for(int x = 0; x < w; x++)
			m[y*w + x] = fabs(x-y);
				
	return NO_ERROR;
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

int Matrix::GetInverseMatrix(double* matrix, double* matrixReversed, int m, double norm, int* transposition_m){
	double temp = 0;
	double temp_m[4];
	const double EPS = 1e-16;

	for (int i = 0; i < m; i++)
		transposition_m[i] = i;

		for (int k = 0; k < m; k++){
				double max_norm = 0.0;
				int i_temp = k, j_temp = k;
				for (int i = k; i < m; i++){
						for (int j = k; j < m; j++){
								temp = matrix[i*m + j];
								if(i == k && j == k) max_norm = fabs (temp);
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

void* Matrix::SolveBlock(void* in){
	cpu_set_t cpu;

	CPU_ZERO (&cpu);
		CPU_SET (get_nprocs () - ((arg*)in)->thread_id - 1, &cpu);
		pthread_setaffinity_np (pthread_self (), sizeof (cpu), &cpu);

	arg* args = (arg*)in;
	int size = args->size;
	int block_size = args->block_size;
		int step = size/block_size;
	int end = size - step*block_size;

	int id = args->thread_id;
	int p = args->thread_count;

	double* rhs = args->rhs;
	double* answer = args->answer;
	double* matrix = args->matrix;

	int index_start;

	address_array a;
	InitializeTempAddresses(&a, size, block_size, end);
		
	//attaches matrice rows to processor caches and reades the matrix from file\equation
		int e = AttachMatrices(args);

	args->return_value = e;

	//no barrier in attach matricces, wait for first thread to read the matrix
		pthread_barrier_wait(args->barrier);
		
	for(int i = 0; i < p; i++){
		if((args + i - id)->return_value != NO_ERROR){
			//other thread had an error, terminate
			DeleteTempAddresses(&a);
						return nullptr;
		}
	}

	pthread_barrier_wait(args->barrier);

	//create a permutation, so we dont take time to actually move elements in the matrix
	auto indexes = new int[step];
	auto indexes_m = new int[size];
	for(int i = 0; i < step; i++)
		indexes[i] = i;

	args->work_time = get_full_time();
	args->cpu_time = get_cpu_time();

	double n = Length(matrix, size, size);

	for(int offset = 0; offset < step; offset++){
		double minimal_norm = 0.0;
		int minimal_norm_index = offset;
				int found_inversable = 0;

		index_start = id < offset % p ? id + p - offset % p : id - offset % p;

		for(int y = offset + index_start; y < step; y+= p){
			GetBlock(matrix, a.block, offset*block_size, y*block_size, 
									offset*block_size + block_size, y*block_size + block_size, size);
			EMatrix(a.inverse_block, block_size);
			if(GetInverseMatrix(a.block, a.inverse_block, block_size, n, indexes_m)){
								found_inversable = 1;
								double k = Length(a.inverse_block, block_size, block_size);
								if(y == offset) minimal_norm = k;
								if(minimal_norm > k){
										minimal_norm = k;
										minimal_norm_index = y;
								}
						}
		}
		if(found_inversable == 0){
			args->return_value = 1;
				}

		args->minimal_norm = minimal_norm;
		args->minimal_index = minimal_norm_index;

		pthread_barrier_wait(args->barrier);

		found_inversable = 0;
		for(int i = 0; i < p; i++){
			if((args + i - id)->return_value != 1){
				if(found_inversable == 0){
					minimal_norm = (args + i - id)->minimal_norm;
					minimal_norm_index = (args + i - id)->minimal_index;
					found_inversable = 1;
				}

				if((args + i - id)->minimal_norm < minimal_norm){
					minimal_norm = (args + i - id)->minimal_norm;
					minimal_norm_index = (args + i - id)->minimal_index;
				}
			}
		}

		if(found_inversable == 0){
			args->return_value = UNINVERTABLE;
			delete[] indexes;
						delete[] indexes_m;
			DeleteTempAddresses(&a);
						return nullptr;
		}

		pthread_barrier_wait(args->barrier);

		args->return_value = 0;
		
		//save indexes to rearrange the answer in the end
		int temp_index = indexes[offset];
		indexes[offset] = indexes[minimal_norm_index];
		indexes[minimal_norm_index] = temp_index;

		if(id == offset%p){
			//physically move the rows 
			memcpy (a.temp_row, matrix + offset * block_size * size, block_size * size * sizeof(double));
			memcpy (matrix + offset * block_size * size, matrix + minimal_norm_index * block_size * size, block_size * size * sizeof(double));
			memcpy (matrix + minimal_norm_index * block_size * size, a.temp_row, block_size * size * sizeof(double));

			memcpy (a.vector_block, rhs + offset * block_size, block_size * sizeof(double));
			memcpy (rhs + offset * block_size, rhs + minimal_norm_index * block_size, block_size * sizeof(double));
			memcpy (rhs + minimal_norm_index * block_size, a.vector_block, block_size * sizeof(double));

			//get inverse block
				GetBlock(matrix, a.block, offset*block_size, offset*block_size, 
									offset*block_size + block_size, offset*block_size + block_size, size);
				EMatrix(a.inverse_block, block_size);
				//we know this matrix exists, no need to check it
				GetInverseMatrix(a.block, a.inverse_block, block_size, n, indexes_m);
				//firstly normalize the rhs vector
				GetBlock(rhs, a.vector_block, 0, offset*block_size, 
										1, offset*block_size + block_size, 1);
			MultiplyMatrices(a.inverse_block, a.vector_block, a.vector_block_temp, block_size, block_size, 1);
				PutBlock(rhs, a.vector_block_temp, 0, offset*block_size, 
											 1, offset*block_size + block_size, 1);

			for(int x = offset+1; x < step; x++){
						GetBlock(matrix, a.block_temp_sub, x*block_size, offset*block_size, 
																		x*block_size + block_size, offset*block_size + block_size, size);
						MultiplyMatrices(a.inverse_block, a.block_temp_sub, a.block_temp, block_size, block_size, block_size);
						PutBlock(matrix, a.block_temp, x*block_size, offset*block_size, 
																		x*block_size + block_size, offset*block_size + block_size, size);
						}

			//if unfull end block exists norm it too
			if(end > 0){
				GetBlock(matrix, a.block_me, step*block_size, offset*block_size, 
											 step*block_size + end, offset*block_size + block_size, size);
				MultiplyMatrices(a.inverse_block, a.block_me, a.block_me_temp, block_size, block_size, end);
						PutBlock(matrix, a.block_me_temp, step*block_size, offset*block_size, 
												step*block_size + end, offset*block_size + block_size, size);
			}
		}

		pthread_barrier_wait(args->barrier);

		//copy the first row, rhs and end part
		memcpy (a.temp_row, matrix + offset * block_size * size, block_size * size * sizeof(double));
		GetBlock(rhs, a.vector_block_temp, 0, offset*block_size, 
										1, offset*block_size + block_size, 1);

		if(end > 0)
			GetBlock(matrix, a.block_me_temp, step*block_size, offset*block_size, 
											 step*block_size + end, offset*block_size + block_size, size);
		
		//everything normalized, block_me_temp has ready me normalized block of the top row
		//vector_block_temp has the same for rhs vector

		index_start = id < (offset+1) % p ? id + p - (offset+1) % p : id - (offset+1) % p;

				//substract top line of blocks from the all the bottom ones of block_size
				for(int y = offset+1+index_start; y < step; y+= p){
						//first element of the current row   
						GetBlock(matrix, a.block, offset*block_size, y*block_size, 
									offset*block_size + block_size, y*block_size + block_size, size);

			MultiplyMatrices(a.block, a.vector_block_temp, a.vector_block_temp_im, block_size, block_size, 1);
			//rhs block
						GetBlock(rhs, a.vector_block, 0, y*block_size, 
										1, y*block_size + block_size, 1);
			SubstractMatrices(a.vector_block, a.vector_block_temp_im, 1, block_size);
						PutBlock(rhs, a.vector_block, 0, y*block_size, 
										1, y*block_size + block_size, 1);
						
						//end block if it exists
						if(end > 0){
				MultiplyMatrices(a.block, a.block_me_temp, a.block_me_temp_im, block_size, block_size, end);
								GetBlock(matrix, a.block_me, step*block_size, y*block_size, 
											 step*block_size + end, y*block_size + block_size, size);
				SubstractMatrices(a.block_me, a.block_me_temp_im, block_size, end);
								PutBlock(matrix, a.block_me, step*block_size, y*block_size, 
											 step*block_size + end, y*block_size + block_size, size);
						}

						//substract all other blocks
						for(int x = offset+1; x < step; x++){
								GetBlock(a.temp_row, a.block_temp, x*block_size, 0, 
																						x*block_size + block_size, 0 + block_size, size);
				MultiplyMatrices(a.block, a.block_temp, a.block_temp_im, block_size, block_size, block_size);
								GetBlock(matrix, a.block_temp, x*block_size, y*block_size, 
											 x*block_size + block_size, y*block_size + block_size, size);
				SubstractMatrices(a.block_temp, a.block_temp_im, block_size, block_size);
								PutBlock(matrix, a.block_temp, x*block_size, y*block_size, 
											 x*block_size + block_size, y*block_size + block_size, size);
						}
				}

		//if the last row exists, substract from them too
		if(id == 0 && end > 0){
			//first element of the last end row   
						GetBlock(matrix, a.block_me, offset*block_size, step*block_size, 
										 offset*block_size + block_size, step*block_size + end, size);
						
						//rhs vector end block
						GetBlock(rhs, a.vector_e, 0, step*block_size, 
									1, step*block_size + end, 1);
			MultiplyMatrices(a.block_me, a.vector_block_temp, a.vector_e_temp, end, block_size, 1);
			SubstractMatrices(a.vector_e, a.vector_e_temp, 1, end);
						PutBlock(rhs, a.vector_e, 0, step*block_size, 
									1, step*block_size + end, 1);

						//ee block
						GetBlock(matrix, a.block_ee, step*block_size, step*block_size, 
										 step*block_size + end, step*block_size + end, size);
			MultiplyMatrices(a.block_me, a.block_me_temp, a.block_ee_temp, end, block_size, end);
			SubstractMatrices(a.block_ee, a.block_ee_temp, end, end);
						PutBlock(matrix, a.block_ee, step*block_size, step*block_size, 
										 step*block_size + end, step*block_size + end, size);
						
						//substract all blocks
						
						for(int x = offset+1; x < step ; x++){
								GetBlock(a.temp_row, a.block_temp, x*block_size, 0, 
											 x*block_size + block_size, block_size, size);
				MultiplyMatrices(a.block_me, a.block_temp, a.block_me_temp_im, end, block_size, block_size);
								GetBlock(matrix, a.block_me_temp_sub, x*block_size, step*block_size, 
													x*block_size + block_size, step*block_size + end, size);
				SubstractMatrices(a.block_me_temp_sub, a.block_me_temp_im, block_size, end);
								PutBlock(matrix, a.block_me_temp_sub, x*block_size, step*block_size, 
													x*block_size + block_size, step*block_size + end, size);
						}
		}

		pthread_barrier_wait(args->barrier);
	}


	if(end > 0 && id == 0){
		GetBlock(matrix, a.block_ee, step*block_size, step*block_size, 
										 step*block_size + end, step*block_size + end, size);
		EMatrix(a.block_ee_temp, end);
		if(!GetInverseMatrix(a.block_ee, a.block_ee_temp, end, n, indexes_m)){
					printf("End part of the matrix could not be inverted\n");
						delete[] indexes;
						delete[] indexes_m;
			DeleteTempAddresses(&a);
			args->return_value = UNINVERTABLE;
						return nullptr;
				}

				GetBlock(rhs, a.vector_e, 0, step*block_size, 
																1, step*block_size + end, 1);
		MultiplyMatrices(a.block_ee_temp, a.vector_e, a.vector_e_temp, end, end, 1);
		PutBlock(rhs, a.vector_e_temp, 0, step*block_size, 
																		 1, step*block_size + end, 1);

				for(int y = 0; y < step; y++){
						GetBlock(matrix, a.block_me, step*block_size, y*block_size, 
										 step*block_size + end, y*block_size + block_size, size);
						MultiplyMatrices(a.block_me, a.vector_e_temp, a.vector_block, block_size, end, 1);
						GetBlock(rhs, a.vector_block_temp, 0, y*block_size,
																						1, y*block_size + block_size, 1);
						SubstractMatrices(a.vector_block_temp, a.vector_block, 1, block_size);
						PutBlock(rhs, a.vector_block_temp, 0, y*block_size,
																						1, y*block_size + block_size, 1);
				}
		} 

	pthread_barrier_wait(args->barrier);

	if((args - id)->return_value == UNINVERTABLE){
		delete[] indexes;
		delete[] indexes_m;
		DeleteTempAddresses(&a);
		args->return_value = 0;
		return nullptr;
	}

	//backwards step
	for (int x = step-1; x > 0; x--){
			GetBlock(rhs, a.vector_block, 0, x * block_size, 
										1, x * block_size + block_size, 1);
				for (int y = id; y < x; y+=p){
						GetBlock (matrix, a.block_temp, x * block_size, y * block_size, 
										x * block_size + block_size, y * block_size + block_size, size);
			MultiplyMatrices(a.block_temp, a.vector_block, a.vector_block_temp, block_size, block_size, 1);
						GetBlock (rhs, a.vector_block_temp_im, 0, y * block_size, 
										 1, y * block_size + block_size, 1);
			SubstractMatrices(a.vector_block_temp_im, a.vector_block_temp, 1, block_size);
						PutBlock (rhs, a.vector_block_temp_im, 0, y * block_size, 
										 1, y * block_size + block_size, 1);
				}
		pthread_barrier_wait(args->barrier);
		}

	if(id == 0){
		//fill the answer
		for(int y = 0; y < size; y++){
				answer[y] = rhs[y];
		}
	}	

	args->work_time = get_full_time() - args->work_time;
	args->cpu_time = get_cpu_time() - args->cpu_time;	


	delete[] indexes;
	delete[] indexes_m;

	DeleteTempAddresses(&a);
		args->return_value = 0;
		return nullptr;
}


























