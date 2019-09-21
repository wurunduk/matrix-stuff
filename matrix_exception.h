#ifndef _MATRIX_EXCEPTION_H_
#define _MATRIX_EXCEPTION_H_

#include <iostream>

enum MatrixException{
    NO_ERROR = 0,
    UNKNOWN_ERROR,
    NOT_ENOUGH_MEMORY,
    CAN_NOT_OPEN,
    FILE_CORRUPT
};

static void ReportError(MatrixException ex){
	switch (ex){
		case NO_ERROR:
		return;
		case UNKNOWN_ERROR:
		printf("unknown error happened while working with matrix, error code %d\n", static_cast<int>(ex));
		return;
		case NOT_ENOUGH_MEMORY:
		printf("Could not allocate enough memory for the matrix.\n");
		return;
		case CAN_NOT_OPEN:
		printf("Could not open the matrix file.\n");
		return;
		case FILE_CORRUPT:
		printf("Matrix file did not have enough numbers to fill the given array.\n");
		return;
	}
}

#endif
