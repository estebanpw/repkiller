#define __STDC_FORMAT_MACROS

#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <ctype.h>
#include <pthread.h>
#include <float.h>
#include <set>
#include "structs.h"
#include "commonFunctions.h"
#include "comparisonFunctions.h"
#include "evolutionaryEventsFunctions.h"

#define LOGT(s) fprintf(stdout, "[INFO] %s T = %e\n", s, (double)(end-begin)/CLOCKS_PER_SEC);

using namespace std;

uint64_t total_bytes_in_use = 0;
int DEBUG_ACTIVE = 0;
int HARD_DEBUG_ACTIVE = 0;

void print_all();
void init_args(int argc, char ** av, FILE ** multifrags, char ** out_file_base_path,
    char * path_frags, uint64_t * ht_size);
int main(int ac, char **av) {
    //Iterator
    uint64_t i;
    //Number of frags loaded, number of sequences loaded
    uint64_t total_frags, n_files;
    //Last ID of synteny lists
    uint64_t last_s_id;
    //Array of fragments to hold them all, pointer to free when changing pointer
    struct FragFile * loaded_frags = NULL;
    //Clocks to measure time
    clock_t begin, end;
    //The sequences manager to store ids, lengths, etc
    sequence_manager * seq_manager = new sequence_manager();
    //The number of iterations to trimm
    char multifrags_path[512]; multifrags_path[0] = '\0';
    //Initial hash table size (divisor of the longest sequence)
    uint64_t ht_size = 100; //Default
    //Default behaviour is the trimmed frags do not already exist
    char * out_file_base_path = new char[MAX_LINE];
    //
    char * out_file_name = new char[MAX_LINE];

    //Open frags file, lengths file and output files %%%%%%%%%%%%%%%%%%%%%%%%%%%
    FILE * frags_file = NULL, * lengths_file = NULL;
    init_args(ac, av, &frags_file, &out_file_base_path, multifrags_path, &ht_size);

    //Concat .lengths to path of multifrags %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    char path_lengths[READLINE];
    path_lengths[0]='\0';
    strcpy(path_lengths, multifrags_path);
    strcat(path_lengths, ".lengths");
    lengths_file = fopen64(path_lengths, "rb");
    if(lengths_file == NULL) terror("Could not open input lengths file");

    //Load lengths and substract accumulated length %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    n_files =  seq_manager->load_sequences_descriptors(lengths_file);
    end = clock();
    LOGT("Loaded sequence descriptors.");
    fclose(lengths_file);

    //Load fragments into array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    load_fragments_local(frags_file, &total_frags, &loaded_frags);
    end = clock();
    print_memory_usage();
    LOGT("Loaded fragments into memory.");
    fclose(frags_file);

    //Frags to blocks conversion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    memory_pool * mp = new memory_pool(POOL_SIZE);
    uint64_t max_len_sequence = seq_manager->get_maximum_length();
    uint64_t coord_aux;
    // ram debug over here
    hash_table * ht = new hash_table(mp, max_len_sequence/ht_size, seq_manager, max_len_sequence);
    for(i=0;i<total_frags;i++){
        //Switch coordinates of reversed fragments. This can only be done at the end of trimming and not meanwhile!
        if(loaded_frags[i].strand == 'r'){ coord_aux = loaded_frags[i].yStart; loaded_frags[i].yStart = loaded_frags[i].yEnd; loaded_frags[i].yEnd = coord_aux;}
        ht->insert_block(&loaded_frags[i]);
    }
    //ht->print_hash_table(2);

    compute_order_of_blocks(ht, n_files);
    end = clock();
    //fprintf(stdout, "[INFO] Insertion of fragments into hash table completed. Load factor = %e. T = %e\n", ht->get_load_factor(), (double)(end-begin)/CLOCKS_PER_SEC);
    print_memory_usage();
    LOGT("Insertion of fragments into hash table completed.");

    //Generate synteny blocks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    auto synteny_block_list = compute_synteny_list(ht, n_files, mp, &last_s_id);
    //traverse_synteny_list_and_write(synteny_block_list, n_files, "init");
    //traverse_synteny_list(synteny_block_list);
    end = clock();
    print_memory_usage();
    LOGT("Generated synteny blocks.");

    // Redo synteny blocks into frags group list %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    auto frags_groups = redo_synteny_system(synteny_block_list);
    end = clock();
    //print_frags_grops_list(ptr_fgl);
    print_memory_usage();
    LOGT("Frags groups list created.");

    // Extend frags groups using parameter 'e' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    begin = clock();
    FGList efrags_groups;
    extend_groups(*frags_groups, &efrags_groups);
    end = clock();
    LOGT("Extensions calculated.");

    begin = clock();
    save_all_frag_pairs(out_file_base_path, seq_manager, efrags_groups);
    end = clock();
    LOGT("Fragments saved into csv file.");

    std::free(loaded_frags);

    delete [] out_file_base_path;
    delete [] out_file_name;

    delete ht;
    delete mp;
    delete seq_manager;

    return 0;
}

void print_all(){
    fprintf(stdout, "USAGE:\n");
    fprintf(stdout, "           repkiller -multifrags [query] -out [results]\n");
    fprintf(stdout, "OPTIONAL:\n");
    fprintf(stdout, "           -hash_table_divisor [Integer:   1<=X] (default 100)\n");
    fprintf(stdout, "           --help      Shows the help for program usage\n");
}

void init_args(int argc, char ** av, FILE ** multifrags, char ** out_file_base_path,
    char * path_frags, uint64_t * ht_size){

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
        if(strcmp(av[pNum], "-hash_table_divisor") == 0){
            *ht_size = (uint64_t) atoi(av[pNum+1]);
            if(*ht_size < 1) terror("The hash table divisor must be one at least");
        }
        pNum++;
    }

    if(*multifrags==NULL ||  *out_file_base_path==NULL || path_frags[0] == '\0'){
        print_all();
        terror("A frags file and an output file must be specified");
    }
}
