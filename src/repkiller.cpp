#define __STDC_FORMAT_MACROS

#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <ctype.h>
#include <pthread.h>
#include <float.h>
#include "structs.h"
#include "commonFunctions.h"
#include "comparisonFunctions.h"

#define LOGT(s) fprintf(stdout, "[INFO] %s T = %e\n", s, (double)(end-begin)/CLOCKS_PER_SEC);

#define LOGI(s) fprintf(stdout, "[INFO] %s\n", s);

uint64_t total_bytes_in_use = 0;
int DEBUG_ACTIVE = 0;
int HARD_DEBUG_ACTIVE = 0;

void print_all();
void init_args(int argc, char ** av, FILE ** multifrags, char ** out_file_base_path,
    char * path_frags);
int main(int ac, char **av) {
    // Clocks to measure time
    clock_t begin, end;
    // The sequences manager to store ids, lengths, etc
    sequence_manager seq_manager;

    char multifrags_path[512]; multifrags_path[0] = '\0';
    char * out_file_base_path = new char[MAX_LINE];
    char * out_file_name = new char[MAX_LINE];

    // Open frags file, lengths file and output files %%%%%%%%%%%%%%%%%%%%%%%%%%%
    FILE * frags_file = NULL, * lengths_file = NULL;
    init_args(ac, av, &frags_file, &out_file_base_path, multifrags_path);

    // Concat .lengths to path of multifrags %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    char path_lengths[READLINE];
    path_lengths[0]='\0';
    strcpy(path_lengths, multifrags_path);
    strcat(path_lengths, ".lengths");
    lengths_file = fopen64(path_lengths, "rb");
    if(lengths_file == NULL) terror("Could not open input lengths file");

    // Load fragments into array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LOGI("Loading fragments into memory...");
    begin = clock();
    FragmentsDatabase frag_db(frags_file, lengths_file, seq_manager);
    end = clock();
    LOGT("Fragments loaded into memory succesfully!");
    fclose(frags_file);
    fclose(lengths_file);

    // Generate fragment groups %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LOGI("Generating fragment groups...");
    begin = clock();
    FGList * efrags_groups = new FGList();
    generate_fragment_groups(frag_db, efrags_groups, 0.8, 8.0);
    end = clock();
    LOGT("Fragment groups generated succesfully!");

    // Save results in file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LOGI("Saving results into secondary memory...");
    begin = clock();
    save_all_frag_pairs(out_file_base_path, seq_manager, *efrags_groups);
    end = clock();
    LOGT("Fragments saved into csv file.");

    // Final cleanup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    delete [] out_file_base_path;
    delete [] out_file_name;
}

void print_all(){
    fprintf(stdout, "USAGE:\n");
    fprintf(stdout, "           repkiller -multifrags [query] -out [results]\n");
    fprintf(stdout, "OPTIONAL:\n");
    fprintf(stdout, "           -hash_table_divisor [Integer:   1<=X] (default 100)\n");
    fprintf(stdout, "           --help      Shows the help for program usage\n");
}

void init_args(int argc, char ** av, FILE ** multifrags, char ** out_file_base_path,
    char * path_frags){

    int pNum = 0;
    while(pNum < argc){
        if(strcmp(av[pNum], "--debug") == 0) DEBUG_ACTIVE = 1;
        if(strcmp(av[pNum], "--hdebug") == 0) HARD_DEBUG_ACTIVE = 1;
        if(strcmp(av[pNum], "--help") == 0){
            print_all();
            exit(1);
        }
        if(strcmp(av[pNum], "-multifrags") == 0){
            *multifrags = fopen64(av[pNum+1], "rb");
            strncpy(path_frags, av[pNum+1], strlen(av[pNum+1]));
            path_frags[strlen(av[pNum+1])] = '\0';
            if(multifrags==NULL) terror("Could not open multifrags file");
        }
        if(strcmp(av[pNum], "-out") == 0){
            if(av[pNum+1]==NULL) terror("ERR¿?");
            strcpy(*out_file_base_path, av[pNum+1]);
        }
        pNum++;
    }

    if(*multifrags==NULL ||  *out_file_base_path==NULL || path_frags[0] == '\0'){
        print_all();
        terror("A frags file and an output file must be specified");
    }
}
