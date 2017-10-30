#include "SequenceOcupationList.h"

SequenceOcupationList::SequenceOcupationList(double len_ratio, double pos_ratio, uint64_t seq_size)
  : len_ratio(len_ratio), pos_ratio(pos_ratio), max_index(seq_size / DIVISOR),
  ocupations(new forward_list<SequenceOcupationList::Ocupation>*[max_index + 1]) {
  for (size_t i = 0; i < max_index + 1; i++)
    ocupations[i] = new forward_list<Ocupation>();
}

SequenceOcupationList::~SequenceOcupationList() {
  for (size_t i = 0; i < max_index + 1; i++)
    delete ocupations[i];
  delete[] ocupations;
}

const forward_list<SequenceOcupationList::Ocupation> * SequenceOcupationList::get_suitable_indices(uint64_t center) const {
  return ocupations[center / DIVISOR];
}

double SequenceOcupationList::deviation(SequenceOcupationList::Ocupation oc, uint64_t center, uint64_t length) const {
  uint64_t dif_len = length > oc.length ? length - oc.length : oc.length - length;
  // Linear transformation of dif len so 0 is TOO OUT and 1 is EXACT MATCH
  auto sim_len = -fabs(dif_len / (length * len_ratio)) + 1.0;
  if (sim_len < 0) return 0.0;

  uint64_t dif_cen = center > oc.center ? center - oc.center : oc.center - center;
  // Linear transformation of dif pos so 0 is TOO DEVIATED and 1 is EXACT MATCH
  auto sim_pos = -fabs(dif_cen / (length * pos_ratio)) + 1.0;
  if (sim_pos < 0) return 0.0;
  return sim_len * 0.4 + sim_pos * 0.6;
}

FragsGroup * SequenceOcupationList::get_associated_group(uint64_t center, uint64_t length) const {
  FragsGroup * fg = nullptr;
  double d = 0;
  {
    auto sind = get_suitable_indices(center);
    for (auto oc : *sind) {
      double curr_dev = deviation(oc, center, length);
      if (curr_dev > d) {
        d = curr_dev;
        fg = oc.group;
      }
    }
  }

  if (center > 0) {
    auto sind = get_suitable_indices(center - 1);
    for (auto oc : *sind) {
      double curr_dev = deviation(oc, center, length);
      if (curr_dev > d) {
        d = curr_dev;
        fg = oc.group;
      }
    }
  }

  if (center < max_index) {
    auto sind = get_suitable_indices(center + 1);
    for (auto oc : *sind) {
      double curr_dev = deviation(oc, center, length);
      if (curr_dev > d) {
        d = curr_dev;
        fg = oc.group;
      }
    }
  }

  if (center > 1) {
    auto sind = get_suitable_indices(center - 2);
    for (auto oc : *sind) {
      double curr_dev = deviation(oc, center, length);
      if (curr_dev > d) {
        d = curr_dev;
        fg = oc.group;
      }
    }
  }

  if (center < max_index - 1) {
    auto sind = get_suitable_indices(center + 2);
    for (auto oc : *sind) {
      double curr_dev = deviation(oc, center, length);
      if (curr_dev > d) {
        d = curr_dev;
        fg = oc.group;
      }
    }
  }
  return fg;
}

void SequenceOcupationList::insert(uint64_t center, uint64_t length, FragsGroup * group) {
  SequenceOcupationList::Ocupation noc(center, length, group);
  ocupations[center / DIVISOR]->push_front(noc);
}
