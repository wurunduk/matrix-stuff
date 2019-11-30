#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <iostream>
#include <pthread.h>
#include <cmath>

typedef struct{
    double* a;
    int n1, n2;
    int m;
    int p;
    
    double left_num;
    double right_num;
    
    double work_time = 0;
    
    pthread_barrier_t* barrier;
} arg;

double get_full_time(){
    struct timeval buffer;
    gettimeofday(&buffer, 0);
    return (double)buffer.tv_sec + (double)buffer.tv_usec/1000000.;
}

double get_time(){
    struct rusage buffer;
    getrusage(RUSAGE_THREAD, &buffer);
    return (double)buffer.ru_utime.tv_sec + (double)buffer.ru_utime.tv_usec/1000000.;
}

void* thread_function(void* in);
void fill_function(double* matrix, int n1, int n2);
int read_file(double* matrix, int n1, int n2, char* file_name);
void print_error(int id);



void* thread_function(void* in){
    arg* args = (arg*)in;
    
    args->work_time = get_time();
    
    int p = args->p;
    int m = args->m;
    double* a = args->a;
    
    int length = args->n1*args->n2/p;
    double cur_num = 0;

    int end_index = length*(m+1);
    //last thread needs to finish the end part aswell
    if(args->m == p - 1){
        end_index += args->n1*args->n2 % p;
    }
    
    args->left_num = a[length*m];
    args->right_num = a[length*m];
    double temp = 0;
    if(m == 0) temp = a[0];
    else temp = a[length*m - 1];
    
    //printf("thread %d works from %d to %d\n", m, length*m, end_index);
    
    for(int i = length*m; i < end_index - 1; i++){
        //no neighbours for end elements
        if(i == 0 || i == args->n1*args->n2 - 1) continue;
        //save the last element of thsi threads array so we can count it later using other thread's result
        if(i == end_index - 2){
            args->right_num = a[i];
        }
        cur_num = a[i];
        a[i] = (temp + a[i+1])/2.0;
        temp = cur_num;
    }
    
    pthread_barrier_wait(args->barrier);
    
    if(length == 1){
        if(m != 0 && m != p-1) a[m] = ((args-1)->right_num + (args+1)->left_num)/2.0;
    }
    else if(m != p-1){
        a[end_index-1] = (args->right_num + (args+1)->left_num)/2.0;
    }
        
    args->work_time = get_time() - args->work_time;
    
    return 0;
}


void print_error(int id){
    switch(id){
        case -1:
            printf("File could not be opened\n");
            break;
        case -2:
            printf("File read error\n");
            break;
        default:
            printf("Unknown error");
    }
}

void print_matrix(double* a, int width, int height, int print_size = 10){
    int pw = width > print_size ? print_size : width;
    int ph = height > print_size ? print_size : height;
    for(int y = 0; y < ph; y++){
        for(int x = 0; x < pw; x++)
            printf("%lf ", a[x+y*width]);
    printf("\n");
    }
}

void fill_function(double* a, int n1, int n2){
    for(int y = 0; y < n2; y++)
        for(int x = 0; x < n1; x++)
            a[x + n1*y] = 1.0/(x+y+1);
}

int read_file(double* matrix, int n1, int n2, char* file_name){
    FILE* f = fopen(file_name, "r");
    if(!f) return -1;
    
    for(int x = 0; x < n1*n2; x++){
        if(fscanf(f, "%lf ", &matrix[x]) != 1){
            return -2;
        }
    }
    
    return 1;
}

int main(int argc, char* argv[]) {
    char* file_name = nullptr;
    int n1=0, n2=0;
    int p = 0;

    if(!(argc==5 || argc==4) || !(p = atoi(argv[1])) || !(n1 = atoi(argv[2])) || !(n2 = atoi(argv[3]))){
        printf("Wrong input format\nUse ./a.out p n1 n2 filename or ./a.out p n1 n2\n");     
        return 1;
    }

    printf("Threads: %d, size %dx%d", p, n1, n2);
    if(file_name) printf(", file %s.", argv[4]);
    printf("\n");
    
    if(p > n1*n2){
        printf("Tried to use more threads than makes sense, using %d threads.\n", n1*n2);
        p = n1*n2;
    }
    
    double* a = new double[n1*n2];
    
    if(argc == 5){
        file_name = argv[4];
        int res = read_file(a, n1, n2, file_name);
        if(res < 0){
            print_error(res);
            delete[] a;
            return 2;
        }
    }
    else{
        fill_function(a, n1, n2);
    }
    
    print_matrix(a, n1, n2);
    
    pthread_barrier_t barrier;
    
    pthread_t* threads = new pthread_t[p];
    arg* args = new arg[p];
    
    pthread_barrier_init(&barrier, 0, p);
    
    double t = get_full_time();
    /*
    double* a;
    int n1, n2;
    int m;
    int p;
    
    double left_num;
    double right_num;
    
    double work_time = 0;
    
    pthread_barrier_t* barrier;
     */
    
    for(int i = 0; i < p; i++){
        args[i].a = a;
        args[i].n1 = n1;
        args[i].n2 = n2;
        args[i].m = i;
        args[i].p = p;
        args[i].left_num = 0;
        args[i].right_num = 0;
        args[i].work_time = 0;
        args[i].barrier = &barrier;
        if(pthread_create(&threads[i], 0, &thread_function, args+i))
        {
            delete[] a;
            delete[] threads;
            delete[] args;
            printf("Could not create a thread %d.\n", i);
            return -1;
        }
    }
    
    for(int i = 0; i < p; i++){
        pthread_join(threads[i], 0);
        printf("Total time taken by thread %d: %lf\n", i, args[i].work_time);
    }
    
    printf("Total time taken: %lf\n", get_full_time() - t);
    
    print_matrix(a, n1, n2);
    
    pthread_barrier_destroy(&barrier);
    
    delete[] a;
    delete[] threads;
    delete[] args;
    
    return 0;
}

























