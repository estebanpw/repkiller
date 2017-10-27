#pragma once

#include <forward_list>
#include <vector>
#include <inttypes.h>
#include <limits>
#include <cmath>

#include "structs.h"

#define DIVISOR 100

class SequenceOcupationList {
  private:
    struct Ocupation {
      uint64_t center;
      uint64_t length;
      FragsGroup * group;
      Ocupation(uint64_t center, uint64_t length, FragsGroup * group) : center(center), length(length), group(group) {};
    };

    double len_ratio, pos_ratio;
    size_t max_index;
    forward_list<Ocupation> ** ocupations;
    const forward_list<Ocupation> * get_suitable_indices(uint64_t center) const;
    double deviation(Ocupation oc, uint64_t center, uint64_t length) const;

  public:
    SequenceOcupationList(double len_pos_ratio, double pos_ratio, uint64_t max_length);
    FragsGroup * get_associated_group(uint64_t center, uint64_t length) const;
    void insert(uint64_t center, uint64_t length, FragsGroup * group);
};
