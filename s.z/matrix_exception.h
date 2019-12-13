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
#endif
