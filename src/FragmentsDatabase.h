#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <utility>
#include <algorithm>

#include "structs.h"

using namespace std;

class FragmentsDatabase {
private:
  unique_ptr<vector<FragFile>[]> loaded_frags;
  size_t vsize;
  uint64_t frags_count = 0;
public:
  FragmentsDatabase(ifstream & frags_file, ifstream & lengths_file, sequence_manager & seq_manager);
  auto getA() const{
    return vsize;
  }
  auto begin() const {
    return loaded_frags.get();
  };
  auto end() const {
    return loaded_frags.get() + vsize - 1;
  };
  uint64_t getTotalFrags() const {
    return frags_count;
  }
};
