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
#include "evolutionaryEventsFunctions.h"

#define LOG(s) fprintf(stdout, "[INFO] %s T = %e\n", s, (double)(end-begin)/CLOCKS_PER_SEC);

uint64_t total_bytes_in_use = 0;
int DEBUG_ACTIVE = 0;
int HARD_DEBUG_ACTIVE = 0;

void print_all();
void init_args(int argc, char ** av, FILE ** multifrags, char ** out_file_base_path,
    uint64_t * min_len_trimming, uint64_t * min_trim_itera, char * path_frags, uint64_t * ht_size,
    FILE ** trim_frags_file, bool * trim_frags_file_write);
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
    //Minimum length to fragment to accept a trimming
    uint64_t min_len = 50; //Default
    //The sequences manager to store ids, lengths, etc
    sequence_manager * seq_manager = new sequence_manager();
    //The number of iterations to trimm
    uint64_t N_ITERA = 0; //Default
    //Path to the multifrags file
    char multifrags_path[512]; multifrags_path[0] = '\0';
    //Initial hash table size (divisor of the longest sequence)
    uint64_t ht_size = 100; //Default
    //Default behaviour is the trimmed frags do not already exist
    bool trim_frags_file_write = true;
    //
    char * out_file_base_path = new char[MAX_LINE];
    //
    char * out_file_name = new char[MAX_LINE];

    //Open frags file, lengths file and output files %%%%%%%%%%%%%%%%%%%%%%%%%%%
    FILE * frags_file = NULL, * lengths_file = NULL, * trim_frags_file = NULL;
    init_args(ac, av, &frags_file, &out_file_base_path, &min_len, &N_ITERA, multifrags_path, &ht_size, &trim_frags_file, &trim_frags_file_write);

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
    LOG("Loaded sequence descriptors.");
    fclose(lengths_file);

    //Load fragments into array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    load_fragments_local(frags_file, &total_frags, &loaded_frags);
    end = clock();
    print_memory_usage();
    LOG("Loaded fragments into memory.");
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
    LOG("Insertion of fragments into hash table completed.");

    //Generate synteny blocks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    Synteny_list * synteny_block_list = compute_synteny_list(ht, n_files, mp, &last_s_id);
    //traverse_synteny_list_and_write(synteny_block_list, n_files, "init");
    //traverse_synteny_list(synteny_block_list);
    end = clock();
    print_memory_usage();
    LOG("Generated synteny blocks.");

    // Redo synteny blocks into frags group list %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    Frags_Groups_List * ptr_fgl = redo_synteny_system(synteny_block_list);
    end = clock();
    //print_frags_grops_list(ptr_fgl);
    print_memory_usage();
    LOG("Frags groups list created.");

    begin = clock();
    save_all_frag_pairs(out_file_base_path, seq_manager, ptr_fgl);
    end = clock();
    LOG("Fragments saved into csv file.");

    // DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //ht->print_hash_table(2);



    //print_maptable_portion(map_table, 0, 194, 50, 0);
    //print_maptable_portion(map_table, 0, 231, 50, 1);

    /*
    print_maptable_portion(map_table, 833568  , 833983, 60, 1);
    print_maptable_portion(map_table, 834344  , 834759, 60, 2);
    print_maptable_portion(map_table, 834361  , 834776, 60, 3);
    */
    /*
    getchar();
    */
    /*
    Bucket * b = ht->get_value_at(208419);
    if(b == NULL) terror("Could not find requested");

    Synteny_list * sbl = find_synteny_block_from_block(sbl, &b->b);
    if(sbl != NULL) printSyntenyListNode(sbl);
    getchar();
    */

    std::free(loaded_frags);

    //fclose(out_file);

    if(trim_frags_file != NULL){
        fclose(trim_frags_file);
    }

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
    fprintf(stdout, "           -min_len_trimming   [Integer:   0<=X] (default 50)\n");
    fprintf(stdout, "           -min_trim_itera     [Integer:   0<=X] (default 500)\n");
    fprintf(stdout, "           -hash_table_divisor [Integer:   1<=X] (default 100)\n");
    fprintf(stdout, "           -reuse_trim_frags   [Path to destination output file]\n");
    fprintf(stdout, "                               Notice that this will save the file if\n");
    fprintf(stdout, "                               it does not exist and load it if it does\n");
    fprintf(stdout, "           --help      Shows the help for program usage\n");
}

void init_args(int argc, char ** av, FILE ** multifrags, char ** out_file_base_path,
    uint64_t * min_len_trimming, uint64_t * min_trim_itera, char * path_frags, uint64_t * ht_size,
    FILE ** trim_frags_file, bool * trim_frags_file_write){

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
        if(strcmp(av[pNum], "-reuse_trim_frags") == 0){
            if(exists_file(av[pNum+1]) == 0){
                *trim_frags_file = fopen64(av[pNum+1], "wb");
                *trim_frags_file_write = true;
                if(trim_frags_file == NULL) terror("Could not open trimmed frags file");
            }else{
                *trim_frags_file = fopen64(av[pNum+1], "rb");
                *trim_frags_file_write = false;
                if(trim_frags_file == NULL) terror("Could not open existing trimmed frags file");
            }
        }
        if(strcmp(av[pNum], "-out") == 0){
            if(av[pNum+1]==NULL) terror("ERR¿?");
            strcpy(*out_file_base_path, av[pNum+1]);
        }
        if(strcmp(av[pNum], "-min_len_trimming") == 0){
            *min_len_trimming = (uint64_t) atoi(av[pNum+1]);
            if(*min_len_trimming < 0) terror("Minimum trimming length must be zero or more");
        }
        if(strcmp(av[pNum], "-min_trim_itera") == 0){
            *min_trim_itera = (uint64_t) atoi(av[pNum+1]);
            if(*min_trim_itera < 0) terror("Minimum number of trimming iterations must be zero or more");
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
