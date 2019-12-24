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
    
    int largest_found_prime = 0;
    
    int cur_found_pairs = 0;
    
    int index;
    int length;
    int* largest_answer;
    int* total_found_pairs;
    
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

int is_prime(int num){
    int is_prime = 1;
    double c = sqrt(num) + 1;
    if(num%2 == 0) return 0;
    for(int i = 3; i < c; i+=2){
        if(num%i == 0){
            //this one is devidable
            is_prime = 0;
            break;
        }
    }
    return is_prime;
}

void* thread_function(void* in);

void* thread_function(void* in){
    arg* args = (arg*)in;
    
    int m = args->m;
    int n = args->n;
    int p = args->p;
    
    int length = args->length;
    
    args->work_time = get_time();
    
    while(*(args->total_found_pairs) < n){        
        args->cur_found_pairs = 0;
        
        int start = args->index*length + ((args->index*length+1)&1);
        
        int next_prime_1 = is_prime(start);
        int next_prime_2 = is_prime(start+2);
        int next_prime_3 = is_prime(start+4);
        
        for(int prime = start; prime < (args->index + 1)*length; prime+=2){
            int next_is_prime = is_prime(prime+6);
            if(next_prime_1){
                if(next_is_prime){
                    args->largest_found_prime = prime+6;
                    args->cur_found_pairs += 1;
                    
                }
            }
            next_prime_1 = next_prime_2;
            next_prime_2 = next_prime_3;
            next_prime_3 = next_is_prime;
        }
        
        pthread_barrier_wait(args->barrier);
        //do the counting in the main thread
        if(m == 0){
            //printf("step %llu\n", *(args->total_found_pairs));
            for(int i = 0; i < p; i++){
                int current_pairs = *(args->total_found_pairs);
                //adding the result of this thread will send over the limit, which means
                //the right pair amount is reached in this thread
                if(current_pairs + (args+i)->cur_found_pairs > n){
                    printf("thread %d went over the limit with %d found pairs, pairs found %d\n", i, (args+i)->cur_found_pairs, current_pairs);
                    start = (args+i)->index*length + (((args+i)->index*length+1)&1);
                    next_prime_1 = is_prime(start);
                    next_prime_2 = is_prime(start+2);
                    next_prime_3 = is_prime(start+4);
                    //recalculate answer of this thread until total_found_pairs = n exactly
                    for(int prime = start; prime < ((args+i)->index + 1)*length; prime+=2){
                        int next_is_prime = is_prime(prime+6); 
                        if(next_prime_1){
                            if(next_is_prime){
                                *(args->largest_answer) = prime+6;
                                *(args->total_found_pairs) += 1;
                                current_pairs += 1;
                                
                                
                                if(current_pairs == n) break;
                            }
                        }
                        next_prime_1 = next_prime_2;
                        next_prime_2 = next_prime_3;
                        next_prime_3 = next_is_prime;
                    }
                    break;
                }
                else if(current_pairs + (args+i)->cur_found_pairs < n){
                    *(args->total_found_pairs) += (args+i)->cur_found_pairs;
                }
                else {
                    *(args->total_found_pairs) += (args+i)->cur_found_pairs;
                    *(args->largest_answer) = (args+i)->largest_found_prime;
                    break;
                }
            }
        }
            
        //printf("%d didthe meme %d %d\n", m, *(args->total_found_pairs), args->cur_found_pairs);
            
        pthread_barrier_wait(args->barrier);
        args->index += p;
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
    
    if(p < 0){
        printf("You need more than 0 threads\n");
        return 1;
    }
    
    n = atoi(argv[2]);
    
    //used by threads for syncing answers and longest arrays
    int answer = 0;
    int total_found_pairs = 0;
    
    int length = 10000;
    
    pthread_barrier_t barrier;
    
    pthread_t* threads = new pthread_t[p];
    arg* args = new arg[p];
    
    pthread_barrier_init(&barrier, 0, p);
    
    double t = get_full_time();
    
    for(int i = 0; i < p; i++){
        args[i].n = n;
        args[i].m = i;
        args[i].p = p;
        args[i].cur_found_pairs = 0;
        args[i].index = i;
        args[i].length = length;
        args[i].largest_answer = &answer;
        args[i].total_found_pairs = &total_found_pairs;
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
    
    printf("%d-th pair number: %d\n", n, answer);
    
    pthread_barrier_destroy(&barrier);
    
    delete[] threads;
    delete[] args;
    
    return 0;
}

























