#include "structs.h"

void sequence_manager::read_header(string input_header) {
        //Allocate heap for sequences struct to hold lengths and ids
       header = input_header;
}

const void sequence_manager::write_header(ofstream & out_file) const {
        // Write headers to output file
        out_file << header;
}

uint64_t sequence_manager::get_maximum_length() const {
        uint64_t i;
        uint64_t m_len = 0;
        for(i=0; i<sequences.size(); i++) {
                if(m_len < this->get_sequence_by_label(i).len) {
                        m_len = this->get_sequence_by_label(i).len;
                }
        }
        return m_len;
}

const Sequence & sequence_manager::get_sequence_by_label(uint64_t label) const {
        return this->sequences[label];
}
