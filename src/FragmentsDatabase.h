#pragma once

#include <iostream>
#include <arpa/inet.h>
#include <vector>

#include "structs.h"

using namespace std;

class FragmentsDatabase {
private:
vector<FragFile*> * loaded_frags;
size_t vsize;
uint64_t num_frags = 0;
static bool readFragment(struct FragFile *frag, FILE *f);
public:
FragmentsDatabase(FILE * frags_file, FILE * lengths_file, sequence_manager & seq_manager);
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
        return num_frags;
}
~FragmentsDatabase();
};
