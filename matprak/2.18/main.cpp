#include <stdio.h>
#include <pthread.h>

typedef struct{
    char* file_name;
    int thread_id;
    int total_threads;
    int return_value;
    int* answer;
    int* lowest_num;
    int* longest_array;
    pthread_mutex_t* mutex;
    pthread_barrier_t* barrier;
} arg;


void* thread_function(void* in);

void* thread_function(void* in){
    arg* args = (arg*)in;
    
    FILE* f = fopen(args->file_name, "r");
    if(!f) {
        args->return_value=-1;
        pthread_barrier_wait(args->barrier);
        pthread_barrier_wait(args->barrier);
        return 0;
    }
    
    double last = 0, cur = 0;
    //elements in the longest array
    int longest_length = 1;
    double lowest_on_longest;
    
    //elements in the current array
    int length = 1;
    double lowest_element;
    
    //total amount of elements in the file lower then given number
    int n = 0;
    
    //read first two elements of the array
    int k = fscanf(f, "%lf %lf ", &last, &cur);
    if(k == 2){
        //check if the first two elements create a descending array
        lowest_element = last;
        if(cur < lowest_element){
            length += 1;
            lowest_on_longest = lowest_element = cur;
        }
        else{
            length = 1;
        }
        lowest_on_longest = lowest_element;
        
        last = cur;
        
        while(fscanf(f, "%lf ", &cur) == 1){
            //printf("%lf ", cur);
            if(cur < lowest_element){
                lowest_element = cur;
                length += 1;
            }
            else{
                //the decreasing stoped, if its longer then the previous, save it
                if(longest_length < length){
                    longest_length = length;
                    lowest_on_longest = lowest_element;
                }
                else{
                    length = 1;
                    lowest_element = cur;
                }
            }
            last = cur;
        }
        if(!feof(f)){
            args->return_value = -2;
            printf("Thread %d could not parse a number!\n", args->thread_id);
            fclose(f);
            pthread_barrier_wait(args->barrier);
            pthread_barrier_wait(args->barrier);
            return 0;
        }
        
        //check if the array which ended on the file end is the lowest
        if(longest_length < length){
            longest_length = length;
            lowest_on_longest = lowest_element;
        }
        else{
            length = 1;
            lowest_element = cur;
        }
    }
    else if(k==1){
        longest_length = 0;
        lowest_element = last;
    }
    else{
        //save error value of the current threads, 
        //fill the berriers so other threads could finish calculations
        args->return_value = -2;
        printf("Thread %d had no numbers!\n", args->thread_id);
        fclose(f);
        pthread_barrier_wait(args->barrier);
        pthread_barrier_wait(args->barrier);
        return 0;
    }
    
    //let every thread save their lowest element and largest array length
    pthread_mutex_lock(args->mutex);
    if(*(args->longest_array) < longest_length){
        *(args->longest_array) = longest_length;
        *(args->lowest_num) = lowest_on_longest;
    }
    pthread_mutex_unlock(args->mutex);
    
    //wait for the sync before continuing
    pthread_barrier_wait(args->barrier);
    
    rewind(f);
    
    double min = *(args->lowest_num);
    
    //iterate through the array and find the elements lower than the given one
    while(fscanf(f, "%lf ", &cur) == 1){
        if(cur < min) n++;
    }
    if(!feof(f)){
        args->return_value = -2;
        printf("Thread %d could not parse a number!\n", args->thread_id);
        fclose(f);
        pthread_barrier_wait(args->barrier);
        return 0;
    }
    
    pthread_mutex_lock(args->mutex);
    *(args->answer) += n;
    pthread_mutex_unlock(args->mutex);
    
    pthread_barrier_wait(args->barrier);
    
    args->return_value = *(args->answer);
    
    fclose(f);
    return 0;
}

int main(int argc, char* argv[]) {
    int files_num = argc-1;

    //used by threads for syncing answers and longest arrays
    int answer = 0;
    int longest_array = 0;
    int lowest_num = 0;
    
    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_barrier_t barrier;
    
    pthread_t* threads = new pthread_t[files_num];
    arg* args = new arg[files_num];
    
    pthread_mutex_init(&mutex, 0);
    pthread_barrier_init(&barrier, 0, files_num);
    
    for(int i = 0; i < files_num; i++){
        args[i].thread_id = i;
        args[i].total_threads = files_num;
        args[i].file_name = argv[i+1];
        args[i].return_value = 0;
        args[i].answer = &answer;
        args[i].lowest_num = &lowest_num;
        args[i].longest_array = &longest_array;
        args[i].mutex = &mutex;
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
    for(int i = 0; i < files_num; i++){
        pthread_join(threads[i], 0);
        printf("Return value of thread %d is %d\n", i, args[i].return_value);
        if(args[i].return_value < 0) error = args[i].return_value;
    }
    if(!error) printf("Answer %d\n", answer);
    else printf("Error %d occured\n", error);
    
    pthread_mutex_destroy(&mutex);
    
    delete[] threads;
    delete[] args;
    
    return 0;
}

























