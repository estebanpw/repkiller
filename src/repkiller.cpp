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

uint64_t total_bytes_in_use = 0;
int DEBUG_ACTIVE = 0;
int HARD_DEBUG_ACTIVE = 0;

void print_all();
void init_args(int argc, char ** av, FILE ** multifrags, FILE ** out_file,
    uint64_t * min_len_trimming, uint64_t * min_trim_itera, char * path_frags, uint64_t * ht_size, 
    char * path_files, FILE ** trim_frags_file, bool * trim_frags_file_write);
void repetitions_detector(Synteny_list * synteny_block_list, sequence_manager * seq_manager, FILE * out_file);

int main(int ac, char **av) {

    
    //Iterator
    uint64_t i;
    //Number of frags loaded, number of sequences loaded
    uint64_t total_frags, n_files;
    //Last ID of synteny lists
    uint64_t last_s_id;
    //Array of fragments to hold them all, pointer to free when changing pointer
    struct FragFile * loaded_frags = NULL, * aux_pointer = NULL;
    //Clocks to measure time
    clock_t begin, end;
    //Minimum length to fragment to accept a trimming
    uint64_t min_len = 50; //Default
    //The sequences manager to store ids, lengths, etc
    sequence_manager * seq_manager = new sequence_manager();
    //The number of iterations to trimm
    uint64_t N_ITERA = 500; //Default
    //Path to the multifrags file
    char multifrags_path[512];
    multifrags_path[0] = '\0';
    //Path to the fastas
    char fastas_path[512];
    fastas_path[0] = '\0';
    //Initial hash table size (divisor of the longest sequence)
    uint64_t ht_size = 100; //Default
    //Default behaviour is the trimmed frags do not already exist
    bool trim_frags_file_write = true; 
    


    //Open frags file, lengths file and output files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FILE * frags_file = NULL, * lengths_file = NULL, * out_file = NULL, * trim_frags_file = NULL;
    init_args(ac, av, &frags_file, &out_file, &min_len, &N_ITERA, multifrags_path, &ht_size, fastas_path, &trim_frags_file, &trim_frags_file_write);

    //Concat .lengths to path of multifrags %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    char path_lengths[READLINE];
    path_lengths[0]='\0';
    strcpy(path_lengths, multifrags_path);
    strcat(path_lengths, ".lengths");
    lengths_file = fopen64(path_lengths, "rb");
    if(lengths_file == NULL) terror("Could not open input lengths file");

    //Load lengths and substract accumulated length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    n_files =  seq_manager->load_sequences_descriptors(lengths_file);
    end = clock();
    fprintf(stdout, "[INFO] Loaded sequence descriptors. T = %e\n", (double)(end-begin)/CLOCKS_PER_SEC);
    fclose(lengths_file);

    //Load fragments into array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    load_fragments_local(frags_file, &total_frags, &loaded_frags);
    end = clock();
    print_memory_usage();
    fprintf(stdout, "[INFO] Loaded fragments into memory. T = %e\n", (double)(end-begin)/CLOCKS_PER_SEC);
    fclose(frags_file); 


    //Initial mapping of fragments to table of genomes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    unsigned char ** map_table = (unsigned char **) std::calloc(n_files, sizeof(unsigned char *));
    total_bytes_in_use += n_files*sizeof(unsigned char *);
    //Allocate map table
    for(i=0; i<n_files; i++){ 
        map_table[i] = (unsigned char *) std::calloc(seq_manager->get_sequence_by_label(i)->len, sizeof(unsigned char)); 
        total_bytes_in_use += seq_manager->get_sequence_by_label(i)->len * sizeof(unsigned char);
        if(map_table[i] == NULL) terror("Could not allocate map table"); 
    }
    map_frags_to_genomes(map_table, loaded_frags, total_frags, seq_manager);
    //Compute initial coverage
    get_coverage_from_genome_grid(map_table, seq_manager, n_files, min_len);
    seq_manager->print_sequences_data();
    end = clock();
    print_memory_usage();
    fprintf(stdout, "[INFO] Initial mapping completed. T = %e\n", (double)(end-begin)/CLOCKS_PER_SEC);

    //Trimming %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    uint64_t ratio_itera = N_ITERA/5;
    //Check if frags had already been computed and exist 
    if(trim_frags_file != NULL && trim_frags_file_write == false){
        // Load fragments trimmed
        total_bytes_in_use -= total_frags * sizeofFragment();
        if(0 == fread(&total_frags, sizeof(uint64_t), 1, trim_frags_file)) terror("Incorrect number of fragments to load");
        loaded_frags = (struct FragFile *) std::realloc(loaded_frags, total_frags*sizeofFragment());
        total_bytes_in_use += total_frags * sizeofFragment();
        if(loaded_frags == NULL) terror("Could not load existing trimmed frags");
        if(0 == fread(loaded_frags, sizeofFragment(), total_frags, trim_frags_file)) terror("Specified fragments are of size zero");
    }else{
        // Compute fragments trim
        total_bytes_in_use -= total_frags * sizeofFragment();
        for(i=0;i<N_ITERA;i++){
            if(i % ratio_itera == 0) fprintf(stdout, "[INFO] Iteration %"PRIu64"\n", i);
            aux_pointer = trim_fragments_and_map(map_table, loaded_frags, &total_frags, min_len);
            //fprintf(stdout, "FRAGS: %"PRIu64"\n", total_frags);
            free(loaded_frags); //A new list is being allocated in the function
            loaded_frags = aux_pointer;
        }
        total_bytes_in_use += total_frags * sizeofFragment();
    }
    
    
    //Compute final coverage
    get_coverage_from_genome_grid(map_table, seq_manager, n_files, min_len);
    //Write frags if it was set
    if(trim_frags_file != NULL && trim_frags_file_write == true){
        //Write number of frags and frags
        fwrite(&total_frags, sizeof(uint64_t), 1, trim_frags_file);
        fwrite(loaded_frags, sizeofFragment(), total_frags, trim_frags_file);
    }
    //seq_manager->print_sequences_data();
    end = clock();
    print_memory_usage();
    fprintf(stdout, "[INFO] Trimming of fragments completed after %"PRIu64" iteration(s).\n       Number of final fragments: %"PRIu64". T = %e\n", N_ITERA, total_frags, (double)(end-begin)/CLOCKS_PER_SEC);

    
    //Frags to blocks conversion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    fprintf(stdout, "[INFO] Insertion of fragments into hash table completed. T = %e\n", (double)(end-begin)/CLOCKS_PER_SEC);
     
    //Generate synteny blocks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    Synteny_list * synteny_block_list = compute_synteny_list(ht, n_files, mp, &last_s_id);
    //traverse_synteny_list_and_write(synteny_block_list, n_files, "init");
    traverse_synteny_list(synteny_block_list);
    end = clock();
    print_memory_usage();
    fprintf(stdout, "[INFO] Generated synteny blocks. T = %e\n", (double)(end-begin)/CLOCKS_PER_SEC);

    //Load sequences %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    seq_manager->read_dna_sequences(fastas_path);
    end = clock();
    print_memory_usage();
    fprintf(stdout, "[INFO] Loaded DNA sequences. T = %e\n", (double)(end-begin)/CLOCKS_PER_SEC);

    //Detect repetitions
    repetitions_detector(synteny_block_list, seq_manager, out_file);
    printf("End of repetitions_detector\n");

    
    
    
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
    
    
    for(i=0; i<n_files; i++){
        std::free(map_table[i]);
    }
    std::free(map_table);
    std::free(loaded_frags);

    fclose(out_file);
    
    if(trim_frags_file != NULL){
        fclose(trim_frags_file);
    }
    

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

void init_args(int argc, char ** av, FILE ** multifrags, FILE ** out_file,
    uint64_t * min_len_trimming, uint64_t * min_trim_itera, char * path_frags, uint64_t * ht_size,
    char * path_files, FILE ** trim_frags_file, bool * trim_frags_file_write){
    
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
            *out_file = fopen64(av[pNum+1], "wt");
            if(out_file==NULL) terror("Could not open output file");
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
        if (strcmp(av[pNum], "-pathfiles") == 0) {
            strncpy(path_files, av[pNum+1], strlen(av[pNum+1]));
            path_files[strlen(av[pNum+1])] = '\0';
        }
        pNum++;
    }
    
    if(*multifrags==NULL || *out_file==NULL || path_frags[0] == '\0'){
        print_all();
        terror("A frags file and an output file must be specified");
    }
}

void repetitions_detector(Synteny_list * synteny_block_list, sequence_manager * seq_manager, FILE * out_file) {
    Synteny_list * pointer_sbl = synteny_block_list;
    Synteny_block * pointer_sb;
    bool advance = true;
    bool found = false;
    uint64_t index = 0;
    uint64_t id;
    uint64_t i;
    uint64_t synteny_block_size;
    uint64_t repetition_index = 0;

    fprintf(out_file, "Repetition;Sequence;Starting position;Ending position;Repetition\n");

    while(pointer_sbl != NULL){
        advance = true;
        pointer_sb = pointer_sbl->sb;
        synteny_block_size = get_synteny_block_size(pointer_sb);

        uint64_t id_repetitions_matrix[synteny_block_size][2];
        for (i = 0; i < synteny_block_size; i++) {
            id_repetitions_matrix[i][1] = 0;
        }

        while(pointer_sb != NULL && advance){
            id = pointer_sb->b->genome->id;
            index = 0;
            while (!found && id_repetitions_matrix[index][1] != 0 && index < synteny_block_size) {
                found = id_repetitions_matrix[index][0] == id;
                index++;
            }
            id_repetitions_matrix[index-1][0] = id;
            id_repetitions_matrix[index-1][1]++;
            if (id_repetitions_matrix[index-1][1] > 1) {
                advance = false;
                Synteny_block * aux_pointer = pointer_sbl->sb;
                while (aux_pointer != NULL) {
                    fprintf(out_file, "%"PRIu64";%"PRIu64";%"PRIu64";%"PRIu64";", repetition_index, aux_pointer->b->genome->id, aux_pointer->b->start, aux_pointer->b->end);
                    seq_manager->print_sequence_region(out_file, aux_pointer->b->genome->id, aux_pointer->b->start, aux_pointer->b->end);
                    aux_pointer = aux_pointer->next;
                }
                repetition_index++;
            }
            found = false;
            pointer_sb = pointer_sb->next;
        }
        pointer_sbl = pointer_sbl->next;
    }
}

