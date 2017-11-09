#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <utility>
#include <algorithm>

#include "structs.h"

#define DIV 60000

using namespace std;

class FragmentsDatabase {
private:
  uint64_t frags_count = 0;
public:
  vector<FragFile> ** loaded_frags;
  size_t nc, nr;
  FragmentsDatabase(ifstream & frags_file, ifstream & lengths_file, sequence_manager & seq_manager, uint64_t threshold);
  ~FragmentsDatabase();
  inline void add(FragFile f) { loaded_frags[f.xStart / DIV][f.yStart / DIV].push_back(f); frags_count++; }
  uint64_t getTotalFrags() const { return frags_count; }

};
