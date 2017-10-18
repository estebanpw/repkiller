#define __STDC_FORMAT_MACROS

#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <ctype.h>
#include "structs.h"
#include "commonFunctions.h"
#include "comparisonFunctions.h"

FragmentsDatabase::FragmentsDatabase(FILE * frags_file, FILE * lengths_file, sequence_manager & seq_manager) {
        uint64_t unused_len, total_frags;

        (void) seq_manager.load_sequences_descriptors(lengths_file);

        //Compute number of fragments in file
        fseeko(frags_file, 0L, SEEK_END);
        total_frags = ftello(frags_file) - 2*(sizeof(uint64_t)); //Remove the headers
        fseeko(frags_file, 0L, SEEK_SET);

        //Divide by size of frag to get the number of fragments
        //Plus one because it might have padding, thus rounding up to bottom and missing 1 struct
        total_frags = 1 + total_frags/sizeof(FragFile);
        //Allocate memory for all frags
        loaded_frags = (struct FragFile *) std::malloc(total_frags * sizeof(FragFile));
        if(loaded_frags == NULL) terror("Could not allocate heap for all fragments");

        //Skip headers
        readSequenceLength(&unused_len, frags_file);
        readSequenceLength(&unused_len, frags_file);

        //To keep track of current frag
        num_frags = 0;

        while(!feof(frags_file)) {
                if (!readFragment(&loaded_frags[num_frags], frags_file)) break;
                if (loaded_frags[num_frags].seqX != 0) continue;
                //printFragment(loaded_frags[num_frags]);
                ++num_frags;
                if(num_frags > total_frags) terror("Something went wrong. More fragments than expected");
        }
        printf("Loaded %" PRIu64 " fragments.\n", num_frags);
}

FragmentsDatabase::~FragmentsDatabase() {
        std::free(loaded_frags);
}

heuristic_sorted_list::heuristic_sorted_list()
{
}

void heuristic_sorted_list::insert_new(uint64_t id, uint64_t value, hs_bucket * prev)
{
        hs_bucket * nb;

        nb = (hs_bucket *)malloc(sizeof(hs_bucket));
        nb->id = id;
        nb->value = value;
        if (prev == NULL) {
                nb->next = this->head;
                head = nb;
        } else {
                nb->next = prev->next;
                prev->next = nb;
        }
}

void heuristic_sorted_list::print() {
        hs_bucket * hsb = this->head;

        while (hsb != NULL) {
                printf("ID: %" PRIu64 ": VAL: %" PRIu64 " -> ", hsb->id, hsb->value);
                hsb = hsb->next;
        }
        printf("\n");
}

void heuristic_sorted_list::insert(uint64_t id, uint64_t value)
{
        hs_bucket * phsb;
        hs_bucket * hsb;

        phsb = NULL;
        hsb = this->head;
        while (hsb != NULL) {
                if (hsb->value >= value) {
                        insert_new(id, value, phsb);
                        return;
                }else{
                        phsb = hsb;
                        hsb = hsb->next;
                }
        }
        insert_new(id, value, hsb);
}

bool heuristic_sorted_list::contains(uint64_t id) {
        hs_bucket * hsb;

        hsb = this->head;
        while (hsb != NULL) {
                if (hsb->id == id) return true;
                else hsb = hsb->next;
        }
        return false;
}

void heuristic_sorted_list::clear() {
        hs_bucket * aux;
        while (this->head != NULL) {
                aux = this->head->next;
                free(aux);
                this->head = aux;
        }
}

heuristic_sorted_list::~heuristic_sorted_list()
{
        clear();
}

sequence_manager::sequence_manager(){
        this->n_sequences = 0;
}

uint64_t sequence_manager::load_sequences_descriptors(FILE * lengths_file){


        //Calculate number of sequences according to size of lengths file
        fseeko(lengths_file, 0L, SEEK_END);
        this->n_sequences = ftello(lengths_file)/sizeof(uint64_t);
        fseeko(lengths_file, 0L, SEEK_SET);

        //Allocate heap for sequences struct to hold lengths and ids
        this->sequences = (Sequence *) std::malloc(n_sequences*sizeof(Sequence));

        if(this->sequences == NULL) terror("Could not allocate memory for sequence descriptors");

        //Load sequence data into sequences descriptors
        uint64_t i=0, acum = 0;
        while(i<this->n_sequences) {
                this->sequences[i].id = i;
                this->sequences[i].acum = acum;
                if(1 != fread(&this->sequences[i].len, sizeof(uint64_t), 1, lengths_file)) terror("Wrong number of sequences or sequence file corrupted");
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
        uint64_t i;
        for(i=0; i<this->n_sequences; i++) {
                if(this->sequences[i].seq != NULL) std::free(this->sequences[i].seq);
        }
        std::free(this->sequences);
}

SequenceOcupationList::SequenceOcupationList(double len_pos_ratio, double threshold, uint64_t seq_size)
        : len_pos_ratio(len_pos_ratio),
        threshold(threshold),
        max_index(seq_size / 100),
        ocupations(new std::forward_list<Ocupation>*[max_index + 1]) {
        for (size_t i = 0; i < max_index + 1; i++)
                ocupations[i] = new std::forward_list<Ocupation>();
}

const std::forward_list<Ocupation> * SequenceOcupationList::get_suitable_indices(uint64_t center) const {
        return ocupations[center / 100];
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

        // Brute force min index calc
        double d = deviation(noc, min_idx, length) * 10;
        while (d < threshold and min_idx > 0) {
                min_idx--;
                d = deviation(noc, min_idx, length) * 10;
        }

        // Brute force max index calc
        d = deviation(noc, max_idx, length) * 10;
        while (d < threshold and max_idx < max_index - 1) {
                max_idx++;
                d = deviation(noc, max_idx, length) * 10;
        }

        for (size_t i = min_idx / 100; i <= max_idx / 100; i++) ocupations[i]->push_front(noc);
}
