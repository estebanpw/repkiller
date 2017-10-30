#pragma once

#include <iostream>
#include <fstream>
#include <arpa/inet.h>
#include <vector>
#include <utility>

#include "structs.h"

using namespace std;

class FragmentsDatabase {
private:
vector<FragFile> * loaded_frags;
size_t vsize;
uint64_t frags_count = 0;
public:
FragmentsDatabase(ifstream & frags_file, ifstream & lengths_file, sequence_manager & seq_manager);
auto getA() {
  return vsize;
}
auto begin() const {
  return loaded_frags;
};
auto end() const {
  return loaded_frags + vsize - 1;
};
uint64_t getTotalFrags() const {
  return frags_count;
}
~FragmentsDatabase();
};
