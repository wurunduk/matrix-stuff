#include <stdio.h>
#include <pthread.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <cmath>
#include <iostream>

typedef struct{
    double* array;
    int size;
    
    int thread_id;
    int total_threads;

    int found_pairs;
    double sum_found;
    
    double time;

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
void print_error(int id);double current_num = 0;
void print_array(double* a, int width, int print_size = 10);

void* thread_function(void* in){
    arg* args = (arg*)in;
    
    args->time = get_time();
    
    double* a = args->array;
    int n = args->size;
    
    int id = args->thread_id;
    int p = args->total_threads;
    
    int length = n/p;
    
    int count = 0;
    double sum = .0;
    
    int end_length = length;
    if(id==p-1){
        end_length += n%p;
    }
    
    double left1 = .0, left2 = .0;
    double right1 = .0, right2 = .0;
    
    //search out division for numbers
    for(int x = 0; x < end_length; x++){
        if(x + length*id - 2 >= 0 && x + length*id + 2 < n){
            if(x==0) {
                left1 = a[length*id + x - 2];
                left2 = a[length*id + x - 1];
            }
            
            if(id == 0){
                left1 = a[0];
                left2 = a[1];
            }
            
            if(x==end_length - 1) {
                right1 = a[length*id + x + 1];
                right2 = a[length*id + x + 2];
            }
            
            if(a[length*id + x-2]*a[length*id + x+2] >= 0 && a[length*id + x] < sqrt( a[length*id + x-2]*a[length*id + x+2] )){
                sum += a[length*id + x];
                count += 1;
            }
        }
    }
    
    //store in the args array
    args->sum_found = sum;
    args->found_pairs = count;
    
    pthread_barrier_wait(args->barrier);
    
    //when everyone stored their results, add them to our result
    for(int i = 0; i < p; i++){
        if(i == id) continue;
        sum += (args+i-id)->sum_found;
        count += (args+i-id)->found_pairs;
    }
    
    if(count == 0){
        args->time = get_time() - args->time;
        return 0;
    }
    
    double new_element = sum/count;
    
    for(int x = 0; x < end_length; x++){
        if(x + length*id  - 2 >= 0 && x + length*id  + 2 < n){
            
            double r = a[length*id + x + 2];
            
            if(x == end_length - 2) r = right1;
            if(x == end_length - 1) r = right2;
            
            if(left1*r >= 0 && a[length*id + x] < sqrt(left1*r)){
                left1 = left2;
                left2 = a[length*id + x];
                a[length*id + x] = new_element;
            }
        }
    }
    
    args->time = get_time() - args->time;
    
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
            printf("Unknown error\n");
    }
}

void print_array(double* a, int width, int print_size){
    int pw = width > print_size ? print_size : width;
    
    for(int x = 0; x < pw; x++)
        printf("%lf ", a[x]);
    printf("\n");
}

void fill_function(double* a, int w){
    for(int x = 0; x < w; x++)
        a[x] = 1.0/(x + 1);
}

int read_file(double* matrix, int w, char* file_name){
    FILE* f = fopen(file_name, "r");
    if(!f) return -1;
    
    for(int x = 0; x < w; x++){
        if(fscanf(f, "%lf ", &matrix[x]) != 1){
            return -2;
        }
    }
    
    return 1;
}

int main(int argc, char* argv[]) {
    char* file_name = nullptr;
    int thread_count = 0;
    int array_size = 0;
    
    pthread_barrier_t barrier;
    
    if(!(argc==3 || argc==4) || !(thread_count = atoi(argv[1])) || !(array_size = atoi(argv[2]))){
        printf("Wrong input format\nUse ./a.out p n or ./a.out p n file_name\n");     
        return 1;
    }
    
    if(array_size < 5) printf("Nothing to do.\n");
    
    if(argc==4) file_name = argv[3];
    
    double* a = new double[array_size];
    pthread_t* threads = new pthread_t[thread_count];
    arg* args = new arg[thread_count];
    
    if(file_name){
        int res = read_file(a, array_size, file_name);
        if(res < 0){
            print_error(res);
            delete[] a;
            delete[] threads;
            delete[] args;
        }
    }
    else fill_function(a, array_size);
    
    print_array(a, array_size);
    
    pthread_barrier_init(&barrier, 0, thread_count);
    
    double t = get_full_time();
    
    for(int i = 0; i < thread_count; i++){
        args[i].array = a;
        
        args[i].size = array_size;
        
        args[i].thread_id = i;
        args[i].total_threads = thread_count;

        args[i].found_pairs = 0;
        args[i].sum_found = .0;
        
        args[i].time = .0;
        
        args[i].barrier = &barrier;
        
        if(pthread_create(&threads[i], 0, &thread_function, args+i))
        {
            delete[] threads;
            delete[] args;
            delete[] a;
            printf("Could not create a thread %d\n", i);
            return -1;
        }
    }

    
    for(int i = 0; i < thread_count; i++){
        pthread_join(threads[i], 0);
        printf("Total time taken by thread %d: %lf\n", i, args[i].time);
    }
    
    printf("Total time taken: %lf\n", get_full_time() - t);
    
    print_array(a, array_size);
    
    pthread_barrier_destroy(&barrier);
    
    delete[] threads;
    delete[] args;
    delete[] a;
    
    return 0;
}

























