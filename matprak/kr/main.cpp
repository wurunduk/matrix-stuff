#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <iostream>
#include <pthread.h>
#include <cmath>

typedef struct{
    int n;
    int m;
    int p;
    
    int cur_found_pairs = 0;
    int sum_found_pairs = 0;
    int index;
    int length;
    int* total_found_pairs;
    int* sum_of_found_pairs;
    
    double work_time = 0;
    
    int first_prime;
    
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

void* thread_function(void* in){
    arg* args = (arg*)in;
    
    args->work_time = get_time();
    
    while(*(args->total_found_pairs) != args->n){
        int first_prime = 0;
    
        int last_prime = 0;
    
        int is_prime = 1;
        
        args->cur_found_pairs = 0;
        args->sum_found_pairs = 0;
        
        for(int prime = args->index*args->length; prime < (args->index + 1)*args->length; prime++){
            is_prime = 1;
            for(int i = 2; i < sqrt(prime); i++){
                if(prime%i == 0){
                    //this one is devidable
                    is_prime = 0;
                    break;
                }
            }
            if(is_prime){
                if(first_prime == 0){
                    first_prime = prime;
                }
                //we know that thedifference is 6 between this two
                if(prime - last_prime == 6){
                    args->cur_found_pairs += 1;
                    args->sum_found_pairs += prime + last_prime;
                }
                last_prime = prime;
            }
        }
        
        args->first_prime = first_prime;
        
        pthread_barrier_wait(args->barrier);
            
            if(args->m < ((args->p) - 1)){
                //if(args->index < 100)
                //        printf("yay %d %d\n", (args+1)->first_prixme, last_prime);
                if((args+1)->first_prime - last_prime == 6){
                    
                    args->cur_found_pairs += 1;
                    args->sum_found_pairs += (args+1)->first_prime + last_prime ;
                }
            }
            
            //if(args->index < 100)
                //printf("Thread %d found %d pairs on interval [%d, %d]\n", args->m, args->cur_found_pairs, args->index*args->length, (args->index + 1)*args->length);
            
            pthread_barrier_wait(args->barrier);
            //do the counting in the main thread
            if(args->m == 0){
                for(int i = 0; i < args->p; i++){
                    int current_pairs = *(args->total_found_pairs);
                    //printf("pairs found on previous steps: %d\n", current_pairs);
                    //int current_answer = *(args->sum_of_found_pairs);
                    //adding the result of this thread will send over the limit, whioch means
                    //the right pair amount is reached in this thread
                    if(current_pairs + (args+i)->cur_found_pairs > args->n){
                        printf("thread %d went over the limit with %d found pairs, pairs found %d, sum %d\n", i, (args+i)->cur_found_pairs, current_pairs, *(args->sum_of_found_pairs));
                        //recalculate answer of this thread until total_found_pairs = n exactly
                        int prime = ((args->index)+i)*args->length;
                        last_prime = 0;
                        (args+i)->cur_found_pairs = 0;
                        (args+i)->sum_found_pairs = 0;
                        while(current_pairs + (args+i)->cur_found_pairs != args->n){
                            is_prime = 1;
                            for(int i = 2; i < sqrt(prime); i++){
                                if(prime%i == 0){
                                    //this one is devidable
                                    is_prime = 0;
                                    break;
                                }
                            }
                            if(is_prime){
                                //we know that thedifference is 6 between this two
                                if(last_prime == 0){
                                    last_prime = prime;
                                }
                                else if(prime - last_prime == 6){
                                    (args+i)->cur_found_pairs += 1;
                                    (args+i)->sum_found_pairs += prime + last_prime;
                                }
                                last_prime = prime;
                            }
                            prime++;
                        }
                        *(args->total_found_pairs) += (args+i)->cur_found_pairs;
                        *(args->sum_of_found_pairs) += (args+i)->sum_found_pairs;
                        break;
                    }
                    else if(current_pairs + (args+i)->cur_found_pairs < args->n){
                        *(args->total_found_pairs) += (args+i)->cur_found_pairs;
                        *(args->sum_of_found_pairs) += (args+i)->sum_found_pairs;
                    }
                    else {
                        printf("thread %d got exactly the limit with %d found pairs, pairs found %d, sum %d\n", i, (args+i)->cur_found_pairs, current_pairs, *(args->sum_of_found_pairs));
                        *(args->total_found_pairs) += (args+i)->cur_found_pairs;
                        *(args->sum_of_found_pairs) += (args+i)->sum_found_pairs;
                        printf("thread %d got exactly the limit with %d found pairs, pairs found %d, sum %d\n", i, (args+i)->cur_found_pairs, (args+i)->sum_found_pairs, *(args->sum_of_found_pairs));
                        break;
                    }
                }
            }
            
        pthread_barrier_wait(args->barrier);
        args->index += args->p;
    }
    
    args->work_time = get_time() - args->work_time;
    
    return 0;
}

int main(int argc, char* argv[]) {
    int p = 1;
    int n = 1;
    
    //double time = get_time();
    
    if(argc!=3){
        printf("Error reading parametres\n");
        return 1;
    }
    
    p = atoi(argv[1]);
    n = atoi(argv[2]);
    
    //used by threads for syncing answers and longest arrays
    int answer = 0;
    int pairs_found = 0;
    
    double t = get_full_time();
    
    int length = 1000;
    if(length < 10) length = 10;
    printf("Length %d %d \n ", length, n);
    
    pthread_barrier_t barrier;
    
    pthread_t* threads = new pthread_t[p];
    arg* args = new arg[p];
    
    pthread_barrier_init(&barrier, 0, p);
    /*
     * int n;
    int m;
    int p;
    
    int cur_found_pairs = 0;
    int index;
    int* total_found_pairs;
    
    int first_prime_distance;
    int first_prime;
    
    pthread_mutex_t* mutex;
    pthread_barrier_t* barrier;
     */
    
    for(int i = 0; i < p; i++){
        args[i].n = n;
        args[i].m = i;
        args[i].p = p;
        args[i].cur_found_pairs = 0;
        args[i].sum_found_pairs = 0;
        args[i].index = i;
        args[i].length = length;
        args[i].total_found_pairs = &pairs_found;
        args[i].sum_of_found_pairs = &answer;
        args[i].first_prime = 0;
        args[i].work_time = 0;
        args[i].barrier = &barrier;
        if(pthread_create(&threads[i], 0, &thread_function, args+i))
        {
            delete[] threads;
            delete[] args;
            printf("Could not create a thread %d\n", i);
            return -1;
        }
    }
    
    for(int i = 0; i < p; i++){
        pthread_join(threads[i], 0);
        printf("Total time taken by thread %d: %lf\n", i, args[i].work_time);
    }
    
    printf("Total time taken: %lf\n", get_full_time() - t);
    
    printf("Found %d pairs summed up to %d\n", pairs_found, answer);
    
    pthread_barrier_destroy(&barrier);
    
    delete[] threads;
    delete[] args;
    
    return 0;
}

























