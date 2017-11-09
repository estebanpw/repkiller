#define __STDC_FORMAT_MACROS

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <thread>
#include <utility>
#include <queue>
#include <limits>

#include "FragmentsDatabase.h"
#include "structs.h"
#include "commonFunctions.h"
#include "SaverQueue.h"

#define MAX_THREADS 1
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
  ifstream frags_file, lengths_file, inf_file, fasta_file_0, fasta_file_1;
  try {
    vector<string> args(argc);
    args.assign(argv, argv + argc);
    init_args(args, frags_file, lengths_file, inf_file, fasta_file_0, fasta_file_1, out_file_base_path, multifrags_path, params);
  } catch (const invalid_argument & e) {
    cerr << e.what() << endl;
    print_help();
    exit(1);
  }

  // Printing the header of the program
  cout << "--- Running REPKILLER v0.9.a by Carles Bordas ---\n"
          "     Bitlab - Arquitectura de Computadores\n"
          "           Universidad de Málaga 2017\n"
          "\n" << flush;

  // The sequences manager to store ids, lengths, etc
  sequence_manager seq_manager;
  (void) seq_manager.load_sequences_descriptors(lengths_file);
  lengths_file.close();
  seq_manager.loadSequence(0, fasta_file_0);
  seq_manager.loadSequence(1, fasta_file_1);
  fasta_file_0.close();
  fasta_file_1.close();

  // Load fragments into array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  FragmentsDatabase frag_db(frags_file, lengths_file, seq_manager, 100);
  frags_file.close();
  inf_file.close();

  cout << "FRAGMENTS AND SEQUENCES LOADED" << endl;
  // Generate subfragments from  overlapping ones %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  auto last = frag_db.getTotalFrags();
  generateSubfragments(frag_db, 100, seq_manager);
  cout << "WENT FROM " << last << " TO " << frag_db.getTotalFrags() << endl;

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
  cout << "GENERATING GROUPS" << endl;
  (void) generate_fragment_groups(frag_db, *efrags_groups, seq_manager, param.first, param.second);

  // Remove 1-element groups %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  size_t psize = efrags_groups->size();
  efrags_groups->erase(std::remove_if(efrags_groups->begin(), efrags_groups->end(), [](const FragsGroup * fgp){ return fgp->size() <= 1; }), efrags_groups->end());
  for (auto fg : *efrags_groups) {
    fg->erase(
      std::unique(
        fg->begin(),
        fg->end(),
        [](auto f1, auto f2){ return f1->xStart == f2->xStart and f1->yStart == f2->yStart and f1->length == f2->length; }
      ),
      fg->end()
    );
  }

  printf("Went from %lu groups to %lu groups\n", psize, efrags_groups->size());

  size_t i = 0;
  for (auto fg : *efrags_groups) {
    std::cout << "GROUP " << i++ << std::endl;
    for (auto fp : *fg) {
      const FragFile & f = *fp;
      printFragment(f);
      auto pair = seq_manager.getSequencePair(f);
      std::cout << "\t" << pair.first << "\n" << std::flush;
      std::cout << "\t" << pair.second << "\n" << std::flush;
      std::cout << std::endl;
    }
    getchar();
  }
  //divideGroups(*efrags_groups, 0.8);

  //exit(0);

  // Generating diagonal function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //size_t * diag_func = new size_t[frag_db.getA()];
  //generate_diagonal_func(frag_db, diag_func);

  // Sort groups by heuristic value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //sort_groups(*efrags_groups, diag_func);
  //delete[] diag_func;

  // Save results in file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  auto out_file_path = out_file_base_path + "-" + std::to_string(param.first) + "-" + std::to_string(param.second) + ".csv";
  sq.addRequest(out_file_path, efrags_groups);
}
