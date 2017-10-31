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

using namespace std;

int main(int argc, char * argv []) {
  string out_file_base_path, multifrags_path;
  double len_pos_ratio, pos_ratio;

  // Open frags file, lengths file and output files %%%%%%%%%%%%%%%%%%%%%%%%%%%
  ifstream frags_file, lengths_file, inf_file;
  try {
    vector<string> args(argc);
    args.assign(argv, argv + argc);
    init_args(args, frags_file, lengths_file, inf_file, out_file_base_path, multifrags_path, len_pos_ratio, pos_ratio);
  } catch (const invalid_argument & e) {
    cerr << e.what() << endl;
    print_help();
    exit(1);
  }

  cout << "--- Running REPKILLER v0.8.b by Carles Bordas ---\n"
          "     Bitlab - Arquitectura de Computadores\n"
          "           Universidad de Málaga 2017\n"
          "\n" << flush;

  // Clocks to measure time
  clock_t begin, end;
  double elapsed;
  // The sequences manager to store ids, lengths, etc
  sequence_manager seq_manager;

  // Load fragments into array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cout << "Loading fragments into memory...\n";
  cout << flush;
  begin = clock();
  FragmentsDatabase frag_db(frags_file, lengths_file, seq_manager);
  end = clock();
  elapsed = (double)(end - begin) / CLOCKS_PER_SEC;
  cout << "Fragments loaded into memory succesfully!\n";
  cout << "\t# Total fragments loaded: " << frag_db.getTotalFrags() << "\n";
  cout << "\t# Elapsed time: " << elapsed << "s\n";
  cout << flush;
  frags_file.close();
  lengths_file.close();
  inf_file.close();

  // Generate fragment groups %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cout << "Generating fragment groups... (this might take a while)\n";
  cout << flush;
  begin = clock();
  FGList efrags_groups;
  (void) generate_fragment_groups(frag_db, efrags_groups, seq_manager, len_pos_ratio, pos_ratio);
  end = clock();
  cout << "Groups created succesfully!\n";
  cout << "\t# Number of groups: " << efrags_groups.size() << "\n";
  cout << "\t# Elapsed time: " << ((double)(end - begin) / CLOCKS_PER_SEC) << " s\n";
  cout << flush;

  // Generating diagonal function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cout << "Generating diagonal function for heuristic calculations\n" << flush;
  begin = clock();
  unique_ptr<size_t[]> diag_func(new size_t[frag_db.getA()]);
  generate_diagonal_func(frag_db, diag_func.get());
  end = clock();
  cout << "Diagonal function created succesfully!\n";
  cout << "\t# Elapsed time: " << ((double)(end - begin) / CLOCKS_PER_SEC) << " s\n" << flush;

  // Sort groups by heuristic value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cout << "Sorting fragments by heuristic value\n" << flush;
  begin = clock();
  sort_groups(efrags_groups, diag_func.get());
  end = clock();
  cout << "Fragments groups sorted succesfully!\n";
  cout << "\t# Elapsed time: " << ((double)(end - begin) / CLOCKS_PER_SEC) << " s\n" << flush;

  // Save results in file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cout << "Saving results into secondary memory...\n";
  cout << flush;
  begin = clock();
  try {
    save_all_frag_pairs(out_file_base_path, seq_manager, efrags_groups);
  } catch (const runtime_error & e) {
    cout << "Couldn't access " << out_file_base_path << ", saving results into ./repkiller_results.csv\n" << flush;
    save_all_frag_pairs("./repkiller_results.csv", seq_manager, efrags_groups);
  }
  end = clock();
  cout << "Fragments saved into csv file.\n";
  cout << "\t# Elapsed time: " << ((double)(end - begin) / CLOCKS_PER_SEC) << " s\n" << flush;

  for (auto fg : efrags_groups) delete fg;
}
