#pragma once

#include <iostream>
#include <arpa/inet.h>

#include "structs.h"

using namespace std;

class FragmentsDatabase {
private:
FragFile * loaded_frags;
uint64_t num_frags = 0;
static bool readFragment(struct FragFile *frag, FILE *f);
public:
FragmentsDatabase(FILE * frags_file, FILE * lengths_file, sequence_manager & seq_manager);
FragFile * getFragAt(size_t index) const {
        return &loaded_frags[index];
}
auto cbegin() const {
        return loaded_frags;
};
auto cend() const {
        return loaded_frags + num_frags;
};
uint64_t getTotalFrags() const {
        return num_frags;
}
~FragmentsDatabase();
};
