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
    
    double* lefts;
    double* rights;
    
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
    
    int n1 = args->n1;
    int n2 = args->n2;
    int p = args->p;    //total threads amount
    int m = args->m;    //current thread id
    double* a = args->a;
    
    double* lefts = args->lefts;
    double* rights = args->rights;
    
    double* temps;
    
    int length = n1/p;
    double cur_num = 0;

    int end_length = length;
    //last thread needs to finish the end part aswell
    if(m == p - 1){
        end_length += n1 % p;
        temps = new double[length + n1%p];
    }
    else{
        temps = new double[length];
    }
    
    for(y = 0; y < n2; y++){
        //FIX HERE SOMETHING PLEASE
        //if thread has a single element, make sure left and right elemt is this element
        lefts[y] = a[y*n1 + m*length];
        rights[y] = a[y*n1 + m*length];
        
        double temps[0];
        if(m == 0) temps[0] = a[y*n1];
        else temps[0] = a[y*n1 + length*m - 1];
        
        //printf("thread %d works from %d to %d\n", m, length*m, end_index);
        
        for(int x = 0; x < end_length - 1; x++){
            //save the last element of this threads array so we can count it later using other thread's result
            if(x == end_length - 2){
                rights[y] = a[y*n1 + m*length + x];
            }
            
            cur_num = a[y*n1 + m*length + x];
            
            double c = 0;
            double sum = 0;
            
            if(y*n1 + m*length + x - 1 > 0 ){
                c+=1;
                sum += temps[0];
            }
            
            if((y-1)*n1 + m*length + x > 0){
                c+=1;
                sum += temps[x];
            }
            
            if(y*n1 + m*length + x + 1 < (y+1)*n1){
                c+= 1;
                sum+= a[y*n1 + m*length + x + 1];
            }
            
            if((y+1)*n1 < n2*n1){
                c+= 1;
                sum += a[(y+1)*n1 + m*length + x];
            }
            
            a[y*n1 + m*length + x] = sum/c;
            temps[x] = cur_num;
        }
    }
    
    
    pthread_barrier_wait(args->barrier);
    
    
        for(int y = 0; y < n2; y++){
            double c = 0;
            double sum = 0;
            
            if(y*n1 + m*length + end_length-1 > 0 ){
                c+=1;
                sum += temps[0];
            }
            
            if((y-1)*n1 + m*length + x > 0){
                c+=1;
                sum += temps[x];
            }
            
            if(y*n1 + m*length + x + 1 < (y+1)*n1){
                c+= 1;
                sum+= a[y*n1 + m*length + x + 1];
            }
            
            if((y+1)*n1 < n2*n1){
                c+= 1;
                sum += a[(y+1)*n1 + m*length + x];
            }
            
            a[y*n1 + m*length + end_length-1] = sum/c;
            temps[x] = cur_num;
        }
    
    
    if(length == 1){
        //if(m != 0 && m != p-1) a[m] = ((args-1)->right_num + (args+1)->left_num)/2.0;
    }
    else if(m != p-1){
        for(int y = 0; y < n2; y++){
            a[y*n1 + m*length + end_length-1] = 
            double c = 0;
            double sum = 0;
            
            if(y*n1 + m*length + end_length-2 > 0 ){
                c+=1;
                sum += temps[0];
            }
            
            if((y-1)*n1 + m*length + x > 0){
                c+=1;
                sum += temps[x];
            }
            
            if(y*n1 + m*length + x + 1 < (y+1)*n1){
                c+= 1;
                sum+= a[y*n1 + m*length + x + 1];
            }
            
            if((y+1)*n1 < n2*n1){
                c+= 1;
                sum += a[(y+1)*n1 + m*length + x];
            }
            
            a[y*n1 + m*length + end_length-1] = sum/c;
            temps[x] = cur_num;
        }
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
    
    double** lefts_array = new double*[p];
    double** rights_array = new double*[p];
    
    for(int i = 0; i < p; i ++){
        lefts_array[i] = new double[n2];
        rights_array[i] = new double[n2];
    }
    
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
    
    double* lefts;
    double* rights;
    
    double work_time = 0;
    
    pthread_barrier_t* barrier;
     */
    
    for(int i = 0; i < p; i++){
        args[i].a = a;
        args[i].n1 = n1;
        args[i].n2 = n2;
        args[i].m = i;
        args[i].p = p;
        args[i].lefts = lefts_array[i];
        args[i].rights = rights_array[i];
        args[i].work_time = 0;
        args[i].barrier = &barrier;
        if(pthread_create(&threads[i], 0, &thread_function, args+i))
        {
            delete[] a;
            delete[] threads;
            delete[] args;
            
            for(int i = 0; i < p; i ++){
                delete[] lefts_array;
                delete[] rights_array;
            }
            
            delete[] lefts_array;
            delete[] rights_array;
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
    
    for(int i = 0; i < p; i ++){
        delete[] lefts_array;
        delete[] rights_array;
    }
    
    delete[] lefts_array;
    delete[] rights_array;
    
    return 0;
}

























