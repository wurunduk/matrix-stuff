#include <stdio.h>
#include <pthread.h>
#include <cmath>

typedef struct{
    char* file_name;
    int thread_id;
    int total_threads;
    int max_count;
    int local_max;
    int did_find_element;
    int* answer;
    int return_value;
} arg;

const double epsilon = 1e-16;

void* thread_function(void* in);

void* thread_function(void* in){
    arg* args = (arg*)in;
    
    FILE* f = fopen(args->file_name, "r");
    if(!f) {
        args->return_value=-1;
        return 0;
    }
    
    double cur = 0;
    double max;
    double max_count = 0;
    
    //read first element of the array
    int k = fscanf(f, "%lf ", &cur);
    if(k == 1){
        max = cur;
        max_count = 1;
        args->did_find_element = 1;
        
        while(fscanf(f, "%lf ", &cur) == 1){
            if(cur > max){
                max = cur;
                max_count = 1;
            }
			else if(fabs(cur - max) < epsilon){
				max_count +=1;
			}
        }
        if(!feof(f)){
            args->return_value = -2;
            printf("Thread %d could not parse a number!\n", args->thread_id);
            fclose(f);
            return 0;
        }
    }
    else{
        args->max_count = 0;
        args->did_find_element = 0;
        //save error value of the current threads, 
        //fill the berriers so other threads could finish calculations
        //args->return_value = -2;
        printf("Thread %d had no numbers!\n", args->thread_id);
        fclose(f);
        return 0;
    }
    
    args->max_count = max_count;
    args->local_max = max;
    fclose(f);
    return 0;
}

void ReportError(int e){
    switch(e){
        case -1:
            printf("File could not be opened.\n");
        break;
        case -2:
            printf("Could not read a number.\n");
        break;
        default:
            printf("I should not be here.\n");
        break;
    }
}

int main(int argc, char* argv[]) {
    int files_num = argc-1;
    
    if(files_num <= 0) {
        printf("No files specified\n");
        return 1;
    }
    
    int answer = 0;
    
    pthread_t* threads = new pthread_t[files_num];
    arg* args = new arg[files_num];
    
    for(int i = 0; i < files_num; i++){
        args[i].thread_id = i;
        args[i].total_threads = files_num;
        args[i].file_name = argv[i+1];
        args[i].answer = &answer;
        args[i].local_max = 0;
        args[i].did_find_element = 0;
        args[i].max_count = 0;
        args[i].return_value = 0;    
        if(pthread_create(&threads[i], 0, &thread_function, args+i))
        {
            delete[] threads;
            delete[] args;
            printf("Could not create a thread %d\n", i);
            return -1;
        }
    }
    
    int error = 0;
    int found_files_amount = 0;
    for(int i = 0; i < files_num; i++){
        pthread_join(threads[i], 0);
        found_files_amount += args[i].did_find_element;
		if(args[i].did_find_element) answer += args[i].max_count;
        printf("Return value of thread %d is %d, nums found %d\n", i, args[i].return_value, args[i].max_count);
        if(args[i].return_value < 0) error = args[i].return_value;
    }

    if(!error && found_files_amount > 0) printf("Answer %d\n", answer);
    else if(found_files_amount == 0) printf("No numbers were found in any file\n");
    else ReportError(error);
    
    delete[] threads;
    delete[] args;
    
    return 0;
}

























