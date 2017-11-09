#include "structs.h"

#define __STDC_FORMAT_MACROS


uint64_t sequence_manager::load_sequences_descriptors(ifstream & lengths_file) {
  //Calculate number of sequences according to size of lengths file
  lengths_file.seekg(0, ios::end);
  size_t n_seqs = lengths_file.tellg() / sizeof(uint64_t);
  lengths_file.seekg(0, ios::beg);

  //Load sequence data into sequences descriptors
  uint64_t acum = 0;
  for(size_t i = 0; i < n_seqs; i++) {
    sequences[i].id = i;
    sequences[i].acum = acum;
    lengths_file.read((char *) &sequences[i].len, sizeof(uint64_t));
    sequences[i].len++;
    acum += sequences[i].len;
    //fprintf(stdout, "[INFO] Sequence %"PRIu64" has length %"PRIu64"\n", i, st[i].len);
  }

  return sequences.size();
}

void sequence_manager::loadSequence(size_t index, std::ifstream & ifile) {
  Sequence & seq = sequences[index];
  seq.nucleotides.clear();
  uint64_t seq_len = seq.len;
  seq.nucleotides.reserve(seq_len);
  ifile.ignore(std::numeric_limits<streamsize>::max(), (int) '\n');
  while (seq.nucleotides.size() < seq_len and ifile.good()) {
    char c = ifile.get();
    if (c == 'A' or c == 'C' or c == 'G' or c == 'T' or c == 'N') {
      seq.nucleotides += c;
    }
  }
}

std::string sequence_manager::getBasesFromFrag(size_t seq_id, const FragFile & f) const {
  if (seq_id == 0)
    return getBasesAt(seq_id, f.xStart, f.xEnd);
  else  {
    if (f.strand == 'f') return getBasesAt(seq_id, f.yStart, f.yEnd);
    else {
      std::string r = getBasesAt(seq_id, f.yEnd, f.yStart);
      std::reverse(r.begin(), r.end());
      return r;
    }
  }

}

std::pair<std::string, std::string> sequence_manager::getSequencePair(const FragFile & f) const {
  return make_pair(getBasesFromFrag(0, f), getBasesFromFrag(1, f));
}

std::string sequence_manager::getBasesAt(size_t seq_id, uint64_t start, uint64_t end) const {
  return sequences[seq_id].nucleotides.substr(start, end - start);
}

uint64_t sequence_manager::get_maximum_length() const {
        uint64_t m_len = 0;
        for(size_t i = 0; i < sequences.size(); i++) {
                if(m_len < this->get_sequence_by_label(i).len) {
                        m_len = this->get_sequence_by_label(i).len;
                }
        }
        return m_len;
}


void sequence_manager::print_sequences_data() const {
        fprintf(stdout, "[INFO] Sequences data:\n");
        for(size_t i = 0; i < sequences.size(); i++) {
                fprintf(stdout, "\t(%" PRIu64 ")\tL:%" PRIu64 "\tC:%" PRIu32 "\tF:%" PRIu64 "\n", i, this->sequences[i].len, this->sequences[i].coverage, this->sequences[i].n_frags);
        }
}

const Sequence & sequence_manager::get_sequence_by_label(uint64_t label) const {
        return sequences.at(label);
}
