#include <stdio.h>
#include <pthread.h>

typedef struct{
    char* file_name;
    int thread_id;
    int total_threads;
    int return_value;
    int* answer;
    
    int* max_nums;
    int* flags;
    pthread_mutex_t* mutex;
    pthread_barrier_t* barrier;
} arg;


void* thread_function(void* in);

void* thread_function(void* in){
    arg* args = (arg*)in;
    
    FILE* f = fopen(args->file_name, "r");
    if(!f) {
        args->return_value=-1;
        printf("Thread %d could not open the file!\n", args->thread_id);
        pthread_barrier_wait(args->barrier);
        pthread_barrier_wait(args->barrier);
        return 0;
    }
    
    double current_num = 0;
    double max_num = 0;
    int k = fscanf(f, "%lf ", &current_num);
    if(k == 1){
        args->flags[args->thread_id] = 1;
    }else{
        //file was empty
        args->flags[args->thread_id] = 0;
        if(!feof(f))
        {
            args->return_value = -2;
            printf("Thread %d could not read a number!\n", args->thread_id);
            fclose(f);
            pthread_barrier_wait(args->barrier);
            pthread_barrier_wait(args->barrier);
            return 0;
        }
    }
    
    max_num = current_num;
    
    while(fscanf(f, "%lf ", &current_num) == 1){
        if(current_num > max_num) max_num = current_num;
    }
    
    if(!feof(f)){
        args->return_value = -3;
        printf("Thread %d did not reach end of the file!\n", args->thread_id);
        fclose(f);
        pthread_barrier_wait(args->barrier);
        pthread_barrier_wait(args->barrier);
        return 0;
    }
    
    args->max_nums[args->thread_id] = max_num;
    
    //wait until other threads get the biggest number on their array
    pthread_barrier_wait(args->barrier);
    
    //got the biggest number in the whole array
    for(int i = 0; i < args->total_threads; i++){
        if(args->max_nums[i] > max_num && args->flags[i]) max_num = args->max_nums[i]; 
    }
    
    max_num /= 2;
    
    int n = 0;
    
    rewind(f);
    
    while(fscanf(f, "%lf ", &current_num) == 1){
        if(current_num > max_num) n += 1;
    }
    
    if(!feof(f)){
        args->return_value = -3;
        printf("Thread %d did not reach end of the file!\n", args->thread_id);
        fclose(f);
        pthread_barrier_wait(args->barrier);
        return 0;
    }
    
    pthread_mutex_lock(args->mutex);
    *args->answer += n;
    pthread_mutex_unlock(args->mutex);
    
    pthread_barrier_wait(args->barrier);
    
    args->return_value = *(args->answer);
    
    fclose(f);
    return 0;
}

int main(int argc, char* argv[]) {
    int files_num = argc-1;
    int answer = 0;
    
    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_barrier_t barrier;
    
    pthread_t* threads = new pthread_t[files_num];
    arg* args = new arg[files_num];
    int* max_nums = new int[files_num];
    int* flags = new int[files_num];
    
    pthread_mutex_init(&mutex, 0);
    pthread_barrier_init(&barrier, 0, files_num);
    
    for(int i = 0; i < files_num; i++){
        args[i].thread_id = i;
        args[i].total_threads = files_num;
        args[i].file_name = argv[i+1];
        args[i].return_value = 0;
        args[i].answer = &answer;
        
        args[i].max_nums = max_nums;
        args[i].flags = flags;
        args[i].mutex = &mutex;
        args[i].barrier = &barrier;
        if(pthread_create(&threads[i], 0, &thread_function, args+i))
        {
            delete[] threads;
            delete[] args;
            delete[] max_nums;
            delete[] flags;
            printf("Could not create a thread %d\n", i);
            return -1;
        }
    }
    
    int error = 0;
    for(int i = 0; i < files_num; i++){
        pthread_join(threads[i], 0);
        printf("Return value of thread %d is %d\n", i, args[i].return_value);
        if(args[i].return_value < 0) error = args[i].return_value;
    }
    if(!error) printf("Answer %d\n", answer);
    else printf("Error %d occured\n", error);
    
    pthread_mutex_destroy(&mutex);
    pthread_barrier_destroy(&barrier);
    
    delete[] threads;
    delete[] args;
    delete[] max_nums;
    delete[] flags;
    
    return 0;
}

























