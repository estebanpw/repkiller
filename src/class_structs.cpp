#include "structs.h"

#define __STDC_FORMAT_MACROS

sequence_manager::sequence_manager(){
        this->n_sequences = 0;
}

uint64_t sequence_manager::load_sequences_descriptors(FILE * lengths_file){


        //Calculate number of sequences according to size of lengths file
        fseeko(lengths_file, 0L, SEEK_END);
        this->n_sequences = ftello(lengths_file)/sizeof(uint64_t);
        fseeko(lengths_file, 0L, SEEK_SET);

        //Allocate heap for sequences struct to hold lengths and ids
        this->sequences = (Sequence *) malloc(n_sequences*sizeof(Sequence));

        if(this->sequences == NULL) throw runtime_error("Could not allocate memory for sequences");

        //Load sequence data into sequences descriptors
        uint64_t i=0, acum = 0;
        while(i<this->n_sequences) {
                this->sequences[i].id = i;
                this->sequences[i].acum = acum;
                if(1 != fread(&this->sequences[i].len, sizeof(uint64_t), 1, lengths_file)) throw runtime_error("Wrong number of sequences or sequence file corrupted");
                this->sequences[i].len = this->sequences[i].len + 1;
                acum += this->sequences[i].len;
                //fprintf(stdout, "[INFO] Sequence %"PRIu64" has length %"PRIu64"\n", i, st[i].len);
                i++;
        }

        return this->n_sequences;
}

uint64_t sequence_manager::get_maximum_length() const {
        uint64_t i;
        uint64_t m_len = 0;
        for(i=0; i<this->n_sequences; i++) {
                if(m_len < this->get_sequence_by_label(i)->len) {
                        m_len = this->get_sequence_by_label(i)->len;
                }
        }
        return m_len;
}


void sequence_manager::print_sequences_data() const {
        uint64_t i;
        fprintf(stdout, "[INFO] Sequences data:\n");
        for(i=0; i<this->n_sequences; i++) {
                fprintf(stdout, "\t(%" PRIu64 ")\tL:%" PRIu64 "\tC:%" PRIu32 "\tF:%" PRIu64 "\n", i, this->sequences[i].len, this->sequences[i].coverage, this->sequences[i].n_frags);
        }
}

Sequence * sequence_manager::get_sequence_by_label(uint64_t label) const {
        if(label < this->n_sequences) return &this->sequences[label]; else return NULL;
}

sequence_manager::~sequence_manager(){
        free(this->sequences);
}

SequenceOcupationList::SequenceOcupationList(double len_pos_ratio, double threshold, uint64_t seq_size)
        : len_pos_ratio(len_pos_ratio),
        threshold(threshold),
        max_index(seq_size / 10),
        ocupations(new forward_list<Ocupation>*[max_index + 1]) {
        for (size_t i = 0; i < max_index + 1; i++)
                ocupations[i] = new forward_list<Ocupation>();
}

const forward_list<Ocupation> * SequenceOcupationList::get_suitable_indices(uint64_t center) const {
        return ocupations[center / 10];
}

double SequenceOcupationList::deviation(Ocupation oc, uint64_t center, uint64_t length) const {
        uint64_t dif_len = length > oc.length ? length - oc.length : oc.length - length;
        uint64_t dif_cen = center > oc.center ? center - oc.center : oc.center - center;
        return dif_len * len_pos_ratio + dif_cen * (1 - len_pos_ratio);
}

FragsGroup * SequenceOcupationList::get_associated_group(uint64_t center, uint64_t length) const {
        FragsGroup * fg = nullptr;
        double d = threshold;
        auto sind = get_suitable_indices(center);
        for (auto oc : *sind) {
                double curr_dev = deviation(oc, center, length);
                if (curr_dev < d) {
                        d = curr_dev;
                        fg = oc.group;
                }
        }
        return fg;
}

void SequenceOcupationList::insert(uint64_t center, uint64_t length, FragsGroup * group) {
        // Obtain index range
        size_t min_idx = center, max_idx = center;
        Ocupation noc(center, length, group);
        double d;
        // Brute force min index calc
        while ((d = deviation(noc, min_idx, length)) < threshold and min_idx > 0) min_idx--;
        // Brute force max index calc
        while ((d = deviation(noc, max_idx, length)) < threshold and max_idx < max_index - 1) max_idx++;
        for (size_t i = min_idx / 10; i <= max_idx / 10; i++) ocupations[i]->push_front(noc);
}
