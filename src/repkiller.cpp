#define __STDC_FORMAT_MACROS

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cstdlib>
#include <string>
#include <ctype.h>
#include <float.h>
#include <thread>

#include "FragmentsDatabase.h"
#include "structs.h"
#include "commonFunctions.h"

#define LOGT(s) fprintf(stdout, "[INFO] %s T = %e\n", s, (double)(end-begin)/CLOCKS_PER_SEC);
#define LOGI(s) fprintf(stdout, "[INFO] %s\n", s);

using namespace std;

int main(int argc, char * argv []) {
  string out_file_base_path, multifrags_path;
  double len_pos_ratio, pos_ratio;
  vector<string> args(argc);
  args.assign(argv, argv + argc);

  // Open frags file, lengths file and output files %%%%%%%%%%%%%%%%%%%%%%%%%%%
  FILE * frags_file = nullptr, * lengths_file = nullptr;
  try {
    init_args(args, frags_file, lengths_file, out_file_base_path, multifrags_path, len_pos_ratio, pos_ratio);
  } catch (const invalid_argument & e) {
    cerr << e.what() << endl;
    print_all();
    throw e;
  }

  cout << "--- Running REPKILLER v0.8 by Carles Bordas ---\n     Bitlab - Arquitectura de Computadores\n           Universidad de Málaga 2017\n\n" << flush;

  // Clocks to measure time
  clock_t begin, end;
  double elapsed;
  // The sequences manager to store ids, lengths, etc
  sequence_manager seq_manager;

  // Load fragments into array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cout << "Loading fragments into memory..." << endl;
  begin = clock();
  FragmentsDatabase frag_db(frags_file, lengths_file, seq_manager);
  end = clock();
  elapsed = (double)(end - begin) / CLOCKS_PER_SEC;
  cout << "Fragments loaded into memory succesfully!" << endl;
  cout << "\t# Total fragments loaded: " << frag_db.getTotalFrags() << endl;
  cout << "\t# Elapsed time: " << elapsed << "s" << endl;
  fclose(frags_file);
  fclose(lengths_file);

  // Generate fragment groups %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cout << "Generating fragment groups... (this might take a while)" << endl;
  begin = clock();
  FGList efrags_groups;
  (void) generate_fragment_groups(frag_db, efrags_groups, seq_manager, len_pos_ratio, pos_ratio);
  end = clock();
  cout << "Groups created succesfully!" << endl;
  cout << "\t# Number of groups: " << efrags_groups.size() << endl;
  cout << "\t# Elapsed time: " << ((double)(end - begin) / CLOCKS_PER_SEC) << " s" << endl;

  // Generating diagonal function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cout << "Generating diagonal function for heuristic calculations" << endl;
  begin = clock();
  size_t * diag_func = new size_t[frag_db.getA()];
  generate_diagonal_func(frag_db, diag_func);
  end = clock();
  cout << "Diagonal function created succesfully!" << endl;
  cout << "\t# Elapsed time: " << ((double)(end - begin) / CLOCKS_PER_SEC) << " s" << endl;

  // Sort groups by heuristic value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cout << "Sorting fragments by heuristic value" << endl;
  begin = clock();
  sort_groups(efrags_groups, diag_func);
  end = clock();
  cout << "Fragments groups sorted succesfully!" << endl;
  cout << "\t# Elapsed time: " << ((double)(end - begin) / CLOCKS_PER_SEC) << " s" << endl;

  // Save results in file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cout << "Saving results into secondary memory..." << endl;
  begin = clock();
  try {
    save_all_frag_pairs(out_file_base_path, seq_manager, efrags_groups);
  } catch (const runtime_error & e) {
    cout << "Coult not access specified output path " << out_file_base_path << ", saving results into ./repkillerresults.csv" << endl;
    save_all_frag_pairs("./repkillerresults.csv", seq_manager, efrags_groups);
  }
  end = clock();
  cout << "Fragments saved into csv file." << endl;
  cout << "\t# Elapsed time: " << ((double)(end - begin) / CLOCKS_PER_SEC) << " s" << endl;

  delete[] diag_func;
}
