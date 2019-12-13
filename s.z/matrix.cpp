#include <stdio.h>
#include <math.h>

#include "matrix.h"

MatrixException FillMatrix(double* m, const int w, const int h){
	for(int y = 0; y < h; y++)
		for(int x = 0; x < w; x++)
			m[y*w + x] = fabs(x-y);

	return NO_ERROR;
}

MatrixException ReadMatrix(double* matrix, const int w, const int h, const char* file_name){
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

double Length(const double* matrix, const int w, const int h){
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

void PrintClean(const double* matrix, const int w, const int h){
    for(int y = 0; y < h; y++){
        for(int x = 0; x < w; x++)
            printf("%lf ", matrix[x + y*w]);
        printf("\n");
    }
    printf("\n");
}

void PrintMatrix(const double* matrix, const int size, int print_size = 10){
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

void PrintVector(const double* vector, const int size, int print_size = 10){
	int n = size;
    if(n > print_size) n = print_size;

    for(int y = 0; y < n-1; y++){
        printf("%lf", vector[y]);
        printf("\n");
    }

	printf("...\n");
    printf("%lf\n", vector[size-1]);
}
