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
