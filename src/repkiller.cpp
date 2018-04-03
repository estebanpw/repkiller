#define __STDC_FORMAT_MACROS

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <thread>
#include <utility>
#include <queue>

#include "FragmentsDatabase.h"
#include "structs.h"
#include "commonFunctions.h"
#include "SaverQueue.h"

#define MAX_THREADS 3
#if MAX_THREADS <= 0
  #error "MAX THREADS MUST BE GREATER THAN ZERO"
#endif

using namespace std;

void execWithParams(const FragmentsDatabase & frag_db, pair<double, double> param,
  const sequence_manager & seq_manager, const string & out_file_base_path, SaverQueue & sq);

int main(int argc, char * argv []) {
  string out_file_base_path, multifrags_path;
  queue<pair<double, double>> params;

  // Open frags file, lengths file and output files %%%%%%%%%%%%%%%%%%%%%%%%%%%
  ifstream frags_file;
  try {
    vector<string> args(argc);
    args.assign(argv, argv + argc);
    init_args(args, frags_file, out_file_base_path, multifrags_path, params);
  } catch (const invalid_argument & e) {
    cerr << e.what() << endl;
    print_help();
    exit(1);
  }

  // Printing the header of the program
  cout << "--- Running REPKILLER v0.9.b ---\n"
          "     Bitlab - Arquitectura de Computadores\n"
          "           Universidad de Málaga 2018\n"
          "\n" << flush;

  // The sequences manager to store ids, lengths, etc
  sequence_manager seq_manager;

  // Load fragments into array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  FragmentsDatabase frag_db(frags_file, seq_manager);

  // Read header
  frags_file.close();

  SaverQueue sq(seq_manager);
  sq.start();

  std::vector<std::thread> threads;
  std::mutex mtx;
  for (size_t i = 0; i < MAX_THREADS; ++i) {
    threads.emplace_back([&] {
      while (true) {
        std::unique_lock<std::mutex> lck(mtx);
        if (params.empty()) return;
        auto param = params.front(); params.pop();
        lck.unlock();
        execWithParams(frag_db, param, seq_manager, out_file_base_path, sq);
      }
    });
  }

  for (auto & thread : threads) thread.join();
  sq.stop();

  cout << "Repkiller finished with no errors\n";
}

void execWithParams(const FragmentsDatabase & frag_db, pair<double, double> param,
  const sequence_manager & seq_manager, const string & out_file_base_path, SaverQueue & sq) {
  // Generate fragment groups %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  FGList * efrags_groups = new FGList;
  (void) generate_fragment_groups(frag_db, *efrags_groups, seq_manager, param.first, param.second);
  // Generating diagonal function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  size_t * diag_func = new size_t[frag_db.getA()];
  generate_diagonal_func(frag_db, diag_func);
  // Sort groups by heuristic value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  sort_groups(*efrags_groups, diag_func);
  delete[] diag_func;
  // Save results in file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  auto out_file_path = out_file_base_path;
  sq.addRequest(out_file_path, efrags_groups);
}
