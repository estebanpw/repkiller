#pragma once

#include <pthread.h>
#include <stdio.h>
#include <cstdlib>
#include <inttypes.h>
#include <ctype.h>
#include <cfloat>
#include <string.h>
#include <math.h>

#include "structs.h"
#include "structs.h"
#include "comparisonFunctions.h"

using namespace std;

/**
 * Print the error message 's' and exit(-1)
 */
void terror(const char *s);

/**
 * Function to read char by char buffered from a FILE
 */
char buffered_fgetc(char *buffer, uint64_t *pos, uint64_t *read, FILE *f);

/*
	Checks if a file exists (1) if so (0) otherwise
*/
int exists_file(const char * file_name);


void generate_fragment_groups(const FragmentsDatabase & frags_db, FGList * efrags_groups, double len_pos_ratio, double threshold);

/*
	Prints memory requested from heap
*/
void print_memory_usage();

/*
  Prints the standard csv header
*/
void write_header(FILE * f, uint64_t sx_len, uint64_t sy_len);

/*
 TODO desc
*/
void save_frags_from_group(FILE * out_file, FragsGroup & fg, heuristic_sorted_list * hsl, uint64_t gid);

/*
 TODO desc
*/
void save_frag_pair(FILE * out_file, uint64_t seq1_label, uint64_t seq2_label, const sequence_manager & seq_mngr, FGList * fgl);

/*
 TODO desc
*/
void save_all_frag_pairs(char * out_file_base_path, const sequence_manager & seq_manager, FGList & fgl);

void print_load(double percentage);
