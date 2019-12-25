#include <stdio.h>
#include <pthread.h>

typedef struct{
    char* file_name;
    int thread_id;
    int total_threads;
    int max_count;
    int local_max;
    int found_max;
    int* answer;
    int return_value;
    pthread_barrier_t* barrier;
} arg;


void* thread_function(void* in);

void* thread_function(void* in){
    arg* args = (arg*)in;
    
    FILE* f = fopen(args->file_name, "r");
    if(!f) {
        args->return_value=-1;
        pthread_barrier_wait(args->barrier);
        return 0;
    }
    
    int p = args->total_threads;
    int id = args->thread_id;
    
    double cur = 0;
    double max;
    double max_count = 0;
    
    //read first element of the array
    int k = fscanf(f, "%lf ", &cur);
    if(k == 1){
        max = cur;
        max_count = 1;
        args->found_max = 1;
        
        while(fscanf(f, "%lf ", &cur) == 1){
            if(cur > max){
                max = cur;
                max_count = 1;
            }
        }
        if(!feof(f)){
            args->return_value = -2;
            printf("Thread %d could not parse a number!\n", args->thread_id);
            fclose(f);
            pthread_barrier_wait(args->barrier);
            return 0;
        }
    }
    else{
        args->max_count = 0;
        args->found_max = 0;
        //save error value of the current threads, 
        //fill the berriers so other threads could finish calculations
        //args->return_value = -2;
        printf("Thread %d had no numbers!\n", args->thread_id);
        fclose(f);
        pthread_barrier_wait(args->barrier);
        return 0;
    }
    
    args->max_count = max_count;
    args->local_max = max;
    
    //wait for the sync before continuing
    pthread_barrier_wait(args->barrier);
    
    int new_max = 0;
    int found_max_total = 0;
    
    if(id == 0){
        for(int i = 0; i < p; i++){
            if((args+i)->found_max == 1){
                if((args+i)->local_max == new_max){
                    found_max_total += (args+i)->max_count;
                }
                else if((args+i)->local_max > new_max){
                    new_max = (args+i)->local_max;
                    found_max_total = (args+i)->max_count;
                }
            }
        }
        
        *(args->answer) = found_max_total;
    }
    
    fclose(f);
    return 0;
}

int main(int argc, char* argv[]) {
    int files_num = argc-1;
    
    if(files_num <= 0) {
        printf("No files specified\n");
        return 1;
    }
    
    //used by threads for syncing answers and longest arrays
    int answer = 0;
    
    pthread_barrier_t barrier;
    
    pthread_t* threads = new pthread_t[files_num];
    arg* args = new arg[files_num];
    
    pthread_barrier_init(&barrier, 0, files_num);
    /*
     * typedef struct{
    char* file_name;
    int thread_id;
    int total_threads;
    int max_count;
    int local_max;
    int* answer;
    pthread_barrier_t* barrier;
    } arg;*/
    for(int i = 0; i < files_num; i++){
        args[i].thread_id = i;
        args[i].total_threads = files_num;
        args[i].file_name = argv[i+1];
        args[i].answer = &answer;
        args[i].local_max = 0;
        args[i].found_max = 0;
        args[i].max_count = 0;
        args[i].return_value = 0;    
        args[i].barrier = &barrier;
        if(pthread_create(&threads[i], 0, &thread_function, args+i))
        {
            delete[] threads;
            delete[] args;
            printf("Could not create a thread %d\n", i);
            return -1;
        }
    }
    
    int error = 0;
    int found_files = 0;
    for(int i = 0; i < files_num; i++){
        pthread_join(threads[i], 0);
        found_files += args[i].found_max;
        printf("Return value of thread %d is %d\n", i, args[i].return_value);
        if(args[i].return_value < 0) error = args[i].return_value;
    }
    if(!error && found_files > 0) printf("Answer %d\n", answer);
    else if(found_files == 0) printf("No numbers were found in any files\n");
    else printf("Error %d occured\n", error);
    
    pthread_barrier_destroy(&barrier);
    
    delete[] threads;
    delete[] args;
    
    return 0;
}

























