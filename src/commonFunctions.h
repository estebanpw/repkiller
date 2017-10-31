#pragma once

#include <memory>
#include <algorithm>
#include <fstream>

#include "structs.h"
#include "FragmentsDatabase.h"
#include "SequenceOcupationList.h"

using namespace std;

void print_help();

void init_args(const vector<string> & args, ifstream & multifrags, ifstream & lengths_file, ifstream & inf_file, string & out_file_base_path,
               string & path_frags, double & len_pos_ratio, double & pos_ratio);
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

void printFragment(const FragFile & f);

size_t generate_fragment_groups(const FragmentsDatabase & frags_db, FGList & efrags_groups, const sequence_manager & seq_manager, double len_pos_ratio, double pos_ratio);

/*
   Prints the standard csv header
 */
void write_header(ofstream & f, uint64_t sx_len, uint64_t sy_len);

/*
   TODO desc
 */
void save_frags_from_group(ofstream & out_file, FragsGroup & fg, uint64_t gid);

/*
   TODO desc
 */
void save_frag_pair(ofstream & out_file, uint64_t seq1_label, uint64_t seq2_label, const sequence_manager & seq_mngr, const FGList & fgl);

/*
   TODO desc
 */
void save_all_frag_pairs(const string & out_file_base_path, const sequence_manager & seq_manager, const FGList & fgl);

void print_load(double percentage);

void sort_groups(FGList & fgl, const size_t * diag_func);

void generate_diagonal_func(const FragmentsDatabase & fdb, size_t * diag_func);
