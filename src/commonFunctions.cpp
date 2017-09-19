#define __STDC_FORMAT_MACROS

#include <stdio.h>
#include <cstdlib>
#include <sys/time.h>
#include <inttypes.h>
#include <ctype.h>
#include <cfloat>
#include <string.h>
#include <math.h>
#include "structs.h"
#include "comparisonFunctions.h"
#include "alignment_functions.h"


void terror(const char *s) {
    printf("ERR**** %s ****\n", s);
    exit(-1);
}

char buffered_fgetc(char *buffer, uint64_t *pos, uint64_t *read, FILE *f) {
    if (*pos >= READBUF) {
        *pos = 0;
        memset(buffer, 0, READBUF);
        *read = fread(buffer, 1, READBUF, f);
    }
    *pos = *pos + 1;
    return buffer[*pos-1];
}

int exists_file(const char * file_name){
    FILE * file;
    if ((file = fopen64(file_name, "r")) != NULL){
        fclose(file);
        return 1;
    }
    return 0;
}

void load_fragments_local(FILE * fragsfile, uint64_t * n_frags, struct FragFile ** loaded_frags){

    struct FragFile temp_frag;
    uint64_t unused_len, total_frags;

    //Compute number of fragments in file
    fseeko(fragsfile, 0L, SEEK_END);
    total_frags = ftello(fragsfile) - 2*(sizeof(uint64_t)); //Remove the headers
    fseeko(fragsfile, 0L, SEEK_SET);

    //Divide by size of frag to get the number of fragments
    //Plus one because it might have padding, thus rounding up to bottom and missing 1 struct
    total_frags = 1 + total_frags/sizeofFragment();

    //Allocate memory for all frags
    struct FragFile * temp_frags_array = (struct FragFile *) std::malloc(total_frags * sizeofFragment());
    if(temp_frags_array == NULL) terror("Could not allocate heap for all fragments");
    total_bytes_in_use += total_frags * sizeofFragment();

    fprintf(stdout, "[INFO] There are %"PRIu64" fragments to be loaded, requiring %"PRIu64" Megabyte(s) of RAM\n", total_frags, (total_frags*sizeof(struct FragFile))/(1024*1024));

    //Skip headers
    readSequenceLength(&unused_len, fragsfile);
    readSequenceLength(&unused_len, fragsfile);

    //To keep track of current frag
    *n_frags = 0;


    while(!feof(fragsfile)){
        readFragment(&temp_frag, fragsfile);

        //printFragment(&temp_frag); getchar();
        //Transform coordinates to local
        //Actually, it is not required anymore.

        /*if(prevx != temp_frag.seqX || prevy != temp_frag.seqY){
            printf("Frag: (%"PRIu64", %"PRIu64") coord: (%"PRIu64", %"PRIu64") diag:: %"PRId64"\n", temp_frag.seqX, temp_frag.seqY, temp_frag.xStart, temp_frag.yStart, temp_frag.diag);
            getchar();
        }
        prevx = temp_frag.seqX; prevy = temp_frag.seqY;*/


        //printf("SeqX SeqY: (%"PRIu64", %"PRIu64")\n", temp_frag.seqX, temp_frag.seqY);
        //printf("Frags. (%"PRIu64", %"PRIu64")\n", sequences[temp_frag.seqX].acum, sequences[temp_frag.seqY].acum);

        temp_frag.xStart = temp_frag.xStart;
        temp_frag.xEnd = temp_frag.xEnd;


        temp_frag.yStart = temp_frag.yStart;
        temp_frag.yEnd = temp_frag.yEnd;


        //temp_frag.xStart = temp_frag.xStart - sequences[temp_frag.seqX].acum;
        //temp_frag.xEnd = temp_frag.xEnd - sequences[temp_frag.seqX].acum;

        //temp_frag.yStart = temp_frag.yStart - sequences[temp_frag.seqY].acum;
        //temp_frag.yEnd = temp_frag.yEnd - sequences[temp_frag.seqY].acum;


        //Copy temp fragment into array
        memcpy(&temp_frags_array[*n_frags], &temp_frag, sizeofFragment());

        *n_frags = *n_frags + 1;

        if(*n_frags > total_frags){ terror("Something went wrong. More fragments than expected");}
    }

    //Copy pointer of array heap
    *loaded_frags = temp_frags_array;
}

void get_coverage_from_genome_grid(unsigned char ** maptable, sequence_manager * seq_manager, uint64_t n_seqs, uint64_t min_len_without_breaks){
    uint64_t i, j;

    uint64_t sum_bases;
    uint64_t current;

    for(i=0;i<n_seqs;i++){
        current = 0;
        sum_bases = 0;
        for(j=0;j<seq_manager->get_sequence_by_label(i)->len;j++){
            if(maptable[i][j] != COVERFRAG){
                if(current >= min_len_without_breaks){
                    sum_bases += current;
                }
                current = 0;
            }else{
                current++;
            }
        }
        seq_manager->get_sequence_by_label(i)->coverage = (100*sum_bases)/seq_manager->get_sequence_by_label(i)->len;
    }
}


void write_maptable_to_disk(unsigned char ** maptable, sequence_manager * seq_manager, const char * out_file_path){
    uint64_t i, j;
    FILE * out_file;

    char output_file[512];

    for(i=0; i<seq_manager->get_number_of_sequences(); i++){

        output_file[0]='\0';
        sprintf(output_file, "%s_%d", out_file_path, (int) i);
        out_file = fopen64(output_file, "wt");
        if(out_file == NULL) terror("Could not open output debug file");

        for(j=0; j<seq_manager->get_sequence_by_label(i)->len; j++){
            fprintf(out_file, "%u", maptable[i][j]);
            if(j != 0 && j % 50 == 0) fprintf(out_file, "\n");
        }
        fprintf(out_file, "\n");

        fclose(out_file);
    }

}

void print_maptable_portion(unsigned char ** maptable, uint64_t from, uint64_t to, uint64_t rate, uint64_t seq){
    uint64_t i, acu = 0;
    static uint64_t CALL = 0;
    fprintf(stdout, "===========(%"PRIu64" -> %"PRIu64")@%"PRIu64"\n", from, to, seq);
    for(i=from;i<to;i++){
        fprintf(stdout, "%u", maptable[seq][i]);
        if(acu != 0 && acu % rate == 0) fprintf(stdout, "\n");
        acu++;
    }
    fprintf(stdout, "\n===========CALL: %"PRIu64"\n", CALL);
    CALL++;
}

void traverse_synteny_list(Synteny_list * sbl){
    Synteny_list * ptr_sbl = sbl;
    Synteny_block * ptr_sb;
    //bool forward = true;
    while(ptr_sbl != NULL){
        fprintf(stdout, "SBL:%"PRIu64"\n", ptr_sbl->id);
        ptr_sb = ptr_sbl->sb;
        while(ptr_sb != NULL){
            fprintf(stdout, "\t");printBlock(ptr_sb->b);
            ptr_sb = ptr_sb->next;
        }
        /*
        if(ptr_sbl->next == NULL){
            forward = !forward;
            printf("And changing direction!!!!!!!\n");
        }

        if(forward){
            ptr_sbl = ptr_sbl->next;
        }else{
            ptr_sbl = ptr_sbl->prev;
        }
        */
        ptr_sbl = ptr_sbl->next;
        //getchar();
    }
}

void traverse_synteny_list_and_write(Synteny_list * sbl, uint64_t n_sequences, char * tag){
    FILE * writer;
    for(uint64_t i=0;i<n_sequences;i++){
        char name[MAX_LINE];
        name[0] = '\0';
        sprintf(name, "%s_%" PRIu64"_%s_%s", "blocks_", i, tag, ".txt");
        writer = fopen64(name, "wt");
        if(writer == NULL) throw "Could not open synteny list write output";
        fprintf(writer, "bID\tini\tfin\tlong\tgen\tsyn\n");
        Synteny_list * ptr_sbl = sbl;
        Synteny_block * ptr_sb;
        //bool forward = true;
        while(ptr_sbl != NULL){
            //fprintf(stdout, "SBL:\n");
            ptr_sb = ptr_sbl->sb;
            while(ptr_sb != NULL){
                //fprintf(stdout, "\t");printBlock(ptr_sb->b);
                if(i==ptr_sb->b->genome->id) printBlockJoseMode(ptr_sb->b, writer);
                ptr_sb = ptr_sb->next;
            }

            ptr_sbl = ptr_sbl->next;
        }
        fclose(writer);
    }
}



Synteny_list * find_synteny_block_from_block(Synteny_list * sbl, Block * b){
    Synteny_list * ptr_sbl = sbl;
    Synteny_block * ptr_sb;
    while(ptr_sbl != NULL){
        ptr_sb = ptr_sbl->sb;
        while(ptr_sb != NULL){
            if(isBlockEqualTo(b, ptr_sb->b) == 1){
                return ptr_sbl;
            }
            ptr_sb = ptr_sb->next;
        }
        ptr_sbl = ptr_sbl->next;
    }
    return NULL;
}

void find_fragments_from_maptable(unsigned char ** maptable, uint64_t start, uint64_t end, uint64_t seq, struct FragFile * frags, uint64_t n_frags){
    uint64_t i = 0;
    struct FragFile target;
    target.xStart = start;
    target.xEnd = end;
    target.seqX = seq;

    while(i<n_frags){
        if(isFragmentEqualTo(&target, &frags[i]) == 1){
            printFragment(&frags[i]);
            //getchar();
        }
        i++;
    }
}

int compare_ranges(Annotation * a, Annotation * b){
    if(a->start <= b->end && b->start <= a->end) return 0; //Overlap
    if(a->start < b->start) return -1;
    return 1;
}

int compare_two_annotations(Annotation * a, Annotation * b){
    if(a->start == b->start) return 0;
    if(a->start < b->start) return -1;
    return 1;
}

Annotation * binary_search_annotations(uint64_t start, uint64_t end, Annotation * anot, uint64_t n_annots){

   uint64_t first = 0;
   uint64_t last = n_annots - 1;
   uint64_t middle = (first+last)/2;
   int compare;
   Annotation aux;
   aux.start = start;
   aux.end = end;

   while (first <= last) {

       compare = compare_ranges(&aux, &anot[middle]);
       if (compare == 0) return &anot[middle];
       if (compare<0) last = middle - 1;
       else first = middle + 1;
       middle = (first + last)/2;
   }
   return NULL; // Not found
}

void quick_sort_annotations(Annotation * array, uint64_t x, uint64_t y) {

    Annotation pivot, aux;
    uint64_t x1, y1;

    memcpy(&pivot, &array[(x+y)/2], sizeofAnnotation());
    x1 = x;
    y1 = y;

    do{

        while (compare_two_annotations(&pivot, &array[x1]) > 0) x1++;
        while (compare_two_annotations(&pivot, &array[y1]) < 0) y1--;
        if (x1 < y1) {
            memcpy(&aux, &array[x1], sizeofAnnotation());
            memcpy(&array[x1], &array[y1], sizeofAnnotation());
            memcpy(&array[y1], &aux, sizeofAnnotation());
            x1++;
            y1--;
        }
        else if (x1 == y1) x1++;
    } while (x1 <= y1);

    if (x < y1) quick_sort_annotations(array, x, y1);
    if (x1 < y) quick_sort_annotations(array, x1, y);
}

uint64_t quick_pow4(uint32_t n){
    static uint64_t pow4[33]={1L, 4L, 16L, 64L, 256L, 1024L, 4096L, 16384L, 65536L,
    262144L, 1048576L, 4194304L, 16777216L, 67108864L, 268435456L, 1073741824L, 4294967296L,
    17179869184L, 68719476736L, 274877906944L, 1099511627776L, 4398046511104L, 17592186044416L,
    70368744177664L, 281474976710656L, 1125899906842624L, 4503599627370496L, 18014398509481984L,
    72057594037927936L, 288230376151711744L, 1152921504606846976L, 4611686018427387904L};
    return pow4[n];
}

uint64_t quick_pow4byLetter(uint32_t n, const char c){
    static uint64_t pow4_G[33]={2*1L, 2*4L, 2*16L, 2*64L, 2*256L, 2*1024L, 2*4096L, 2*16384L, 2*65536L,
    (uint64_t)2*262144L, (uint64_t)2*1048576L,(uint64_t)2*4194304L, (uint64_t)2*16777216L, (uint64_t)2*67108864L, (uint64_t)2*268435456L, (uint64_t)2*1073741824L, (uint64_t)2*4294967296L,
    (uint64_t)2*17179869184L, (uint64_t)2*68719476736L, (uint64_t)2*274877906944L, (uint64_t)2*1099511627776L, (uint64_t)2*4398046511104L, (uint64_t)2*17592186044416L,
    (uint64_t)2*70368744177664L, (uint64_t)2*281474976710656L, (uint64_t)2*1125899906842624L, (uint64_t)2*4503599627370496L, (uint64_t)2*18014398509481984L,
    (uint64_t)2*72057594037927936L, (uint64_t) 2*288230376151711744L, (uint64_t) 2*1152921504606846976L, (uint64_t) 2*4611686018427387904L};

    static uint64_t pow4_T[33]={3*1L, 3*4L, 3*16L, 3*64L, 3*256L, 3*1024L, 3*4096L, 3*16384L, 3*65536L,
    (uint64_t)3*262144L, (uint64_t) 3*1048576L, (uint64_t)3*4194304L, (uint64_t)3*16777216L, (uint64_t)3*67108864L, (uint64_t)3*268435456L, (uint64_t)3*1073741824L, (uint64_t)3*4294967296L,
    (uint64_t)3*17179869184L, (uint64_t)3*68719476736L, (uint64_t)3*274877906944L, (uint64_t)3*1099511627776L, (uint64_t)3*4398046511104L, (uint64_t)3*17592186044416L,
    (uint64_t)3*70368744177664L, (uint64_t)3*281474976710656L, (uint64_t)3*1125899906842624L, (uint64_t)3*4503599627370496L, (uint64_t)3*18014398509481984L,
    (uint64_t)3*72057594037927936L, (uint64_t) 3*288230376151711744L, (uint64_t) 3*1152921504606846976L, (uint64_t) 3*4611686018427387904L};

    if(c == 'A') return 0;
    if(c == 'C') return quick_pow4(n);
    if(c == 'G') return pow4_G[n];
    if(c == 'T') return pow4_T[n];
    return 0;
}

uint64_t hashOfWord(const char * word, uint32_t k){

    uint64_t value = 0, jIdx;
    for(jIdx=0;jIdx<k;jIdx++){
        value += quick_pow4byLetter(k-(jIdx+1), word[jIdx]);
    }
    return value;

}

inline int64_t compare_letters(char a, char b){
    if(a != 'N') return (a == b) ? POINT : -POINT;
    return -POINT;
}

int overlapped_words(uint64_t xstart, uint64_t xend, uint64_t ystart, uint64_t yend){
    if(xstart <= yend && ystart <= xend) return 0; //Overlap
    if(xstart < ystart) return -1;
    return 1;
}

void alignment_from_hit(sequence_manager * seq_man, Word * a, Word * b, Quickfrag * qf, uint64_t kmer_size){

    //@Important: a->pos and b->pos should be ending of the hit + 1

    /*
    printf("From hits\n");
    seq_man->print_sequence_region(a->genome->id, a->pos - kmer_size, a->pos);
    seq_man->print_sequence_region(b->genome->id, b->pos - kmer_size, b->pos);
    */

    int64_t curr_pos_a = (int64_t) a->pos;
    int64_t curr_pos_b = (int64_t) b->pos;
    int64_t final_end_a = (int64_t) a->pos - 1, final_start_a = final_end_a - kmer_size + 1, final_start_b = curr_pos_b - kmer_size;
    int64_t score_right = (int64_t) kmer_size * POINT;
    int64_t score_left = score_right;
    int64_t high_left = score_left, high_right = score_right;
    int64_t start_block_a = (int64_t) a->b->start;
    int64_t end_block_a = (int64_t) a->b->end;
    int64_t start_block_b = (int64_t) b->b->start;
    int64_t end_block_b = (int64_t) b->b->end;
    qf->t_len = kmer_size;
    uint64_t idents = kmer_size;
    uint64_t final_idents = 0;

    int keep_going = 1;

    //Forward search
    while(keep_going == 1){


        if(score_right > 0 && curr_pos_a < end_block_a && curr_pos_b < end_block_b){
            if(curr_pos_a  > end_block_a ||  curr_pos_b > end_block_b) break;
            if(compare_letters(a->b->genome->seq[curr_pos_a], b->b->genome->seq[curr_pos_b]) == POINT){ score_right+=POINT; idents++; }else{ score_right-=POINT;}
            if(high_right <= score_right){
                final_end_a = curr_pos_a;
                high_right = score_right;
                final_idents = idents;
            }
            curr_pos_a++;
            curr_pos_b++;
        }else{
            keep_going = 0;
        }
    }

    keep_going = 1;
    curr_pos_a = a->pos - kmer_size - 1;
    curr_pos_b = b->pos - kmer_size - 1;

    score_left = high_right;

    //Backward search
    while(keep_going == 1){

        if(score_left > 0 && curr_pos_a >= start_block_a && curr_pos_b >= start_block_b){
            if(curr_pos_a < start_block_a || curr_pos_b < start_block_b ) break;
            if(compare_letters(a->b->genome->seq[curr_pos_a], b->b->genome->seq[curr_pos_b]) == POINT){ score_left+=POINT; idents++; }else{ score_left-=POINT;}
            if(high_left <= score_left){
                final_start_a = curr_pos_a;
                final_start_b = curr_pos_b;
                high_left = score_left;
                final_idents = idents;
            }
            curr_pos_a--;
            curr_pos_b--;
        }else{
            keep_going = 0;
        }
    }



    qf->t_len = final_end_a - final_start_a + 1;
    qf->sim = (final_idents / (long double) qf->t_len)*100;
    qf->x_start = final_start_a;
    qf->y_start = final_start_b;
    qf->diag = (int64_t)qf->x_start - (int64_t)qf->y_start;
    qf->x = a->b->genome;
    qf->y = b->b->genome;

    /*
    printf("Aligned\n");
    seq_man->print_sequence_region(a->genome->id, final_start_a, final_end_a);
    seq_man->print_sequence_region(b->genome->id, final_start_b, final_start_b+qf->t_len);
    getchar();
    */

}


inline char complement(char c){

    switch(c){
    case ('A'): return 'T';
    break;
    case ('C'): return 'G';
    break;
    case ('G'): return 'C';
    break;
    case ('T'): return 'A';
    break;
    case ('N'): return 'N';
    break;
    case ('-'): return '-'; // Gaps
    break;
    case ('\0'): return '\0';

    }

    printf("T B: (%c) ____\n", c);
    terror("Unrecognized nucleotide in hit");
    return 0;
}


void alignment_from_hit_reverse(sequence_manager * seq_man, Word * a, Word * b, Quickfrag * qf, uint64_t kmer_size){

    //@Important: a->pos and b->pos should be ending of the hit + 1
    //@Important: a will be reversed
    //@Important: Everything is the same except how the access is done

    /*
    printf("From hits\n");
    seq_man->print_sequence_region(a->genome->id, a->pos - kmer_size, a->pos);
    seq_man->print_sequence_region(b->genome->id, b->pos - kmer_size, b->pos);
    */

    int64_t curr_pos_a = (int64_t) a->pos - kmer_size - 1;
    int64_t curr_pos_b = (int64_t) b->pos;
    int64_t final_end_a = (int64_t) curr_pos_a, final_start_a = a->pos - 1, final_start_b = curr_pos_b - kmer_size;
    int64_t score_right = (int64_t) kmer_size * POINT;
    int64_t score_left = score_right;
    int64_t high_left = score_left, high_right = score_right;
    int64_t start_block_a = (int64_t) a->b->start;
    int64_t end_block_a = (int64_t) a->b->end;
    int64_t start_block_b = (int64_t) b->b->start;
    int64_t end_block_b = (int64_t) b->b->end;
    qf->t_len = kmer_size;
    uint64_t idents = kmer_size;
    uint64_t final_idents = 0;

    int keep_going = 1;


    //printf("From: %"PRId64" to %"PRId64"\n", start_block_a, end_block_a);
    //printf("currspos. %"PRId64"\n", curr_pos_a);
    //getchar();

    //Forward search
    while(keep_going == 1){


        if(score_right > 0 && curr_pos_a >= start_block_a && curr_pos_b < end_block_b){
            if(curr_pos_a  < start_block_a ||  curr_pos_b > end_block_b) break;
            if(compare_letters(complement(a->b->genome->seq[curr_pos_a]), b->b->genome->seq[curr_pos_b]) == POINT){ score_right+=POINT; idents++; }else{ score_right-=POINT;}
            if(high_right <= score_right){
                final_end_a = curr_pos_a;
                high_right = score_right;
                final_idents = idents;
            }
            curr_pos_a--;
            curr_pos_b++;
        }else{
            keep_going = 0;
        }
    }

    keep_going = 1;
    curr_pos_a = a->pos;
    curr_pos_b = b->pos - kmer_size - 1;

    score_left = high_right;

    //printf("From: %"PRId64" to %"PRId64"\n", start_block_a, end_block_a);
    //printf("currspos. %"PRId64"\n", curr_pos_a);
    //printf("Aktuelle: %"PRId64" and %"PRId64"\n", final_start_a, final_end_a);
    //getchar();

    //Backward search
    while(keep_going == 1){

        if(score_left > 0 && curr_pos_a < end_block_a && curr_pos_b >= start_block_b){
            if(curr_pos_a >= end_block_a || curr_pos_b < start_block_b ) break;
            if(compare_letters(complement(a->b->genome->seq[curr_pos_a]), b->b->genome->seq[curr_pos_b]) == POINT){ score_left+=POINT; idents++; }else{ score_left-=POINT;}
            if(high_left <= score_left){
                final_start_a = curr_pos_a;
                final_start_b = curr_pos_b;
                high_left = score_left;
                final_idents = idents;
            }
            curr_pos_a++;
            curr_pos_b--;
        }else{
            keep_going = 0;
        }
    }

    //printf("Ending\n");
    //printf("Aktuelle: %"PRId64" and %"PRId64"\n", final_start_a, final_end_a);
    //getchar();

    qf->t_len = final_start_a - final_end_a + 1;
    qf->sim = (final_idents / (long double) qf->t_len)*100;
    qf->x_start = final_end_a;
    qf->y_start = final_start_b;
    qf->diag = (int64_t)qf->x_start - (int64_t)qf->y_start;
    qf->x = a->b->genome;
    qf->y = b->b->genome;
    /*
    printf("Aligned\n");
    seq_man->print_sequence_region(a->genome->id, final_start_a, final_end_a);
    seq_man->print_sequence_region(b->genome->id, final_start_b, final_start_b+qf->t_len);
    getchar();
    */

}

void inplace_dna_switch(char *a, char *b, uint64_t l1, uint64_t l2, uint64_t t){
    uint64_t i;
    char aux;
    for(i=0;i<t;i++){
        aux = a[i+l1];
        a[i+l1] = b[i+l2];
        b[i+l2] = aux;
    }
}

void inplace_reverse_and_complement(char *d, uint64_t l){
    uint64_t i;
    char c;
    for(i=0;i<l/2;i++){
        c = d[i];
        d[i] = d[l-i-1];
        d[l-i-1] = c;
        //printf("%c%c", d[i], d[l-i-1]);
    }

    //printf("\n");
    for(i=0;i<l;i++){
        //printf("eeee: %d", (int)d[i]);
        d[i] = complement(d[i]);


    }
}


inline void strrev(char *p, char *d, uint32_t k){
    uint32_t i;
    char c;
    for(i=0;i<k;i++){
    	c = p[k-i-1];
    	switch(c){
    	case ('A'): c=('T');
    	break;
    	case ('C'): c=('G');
    	break;
    	case ('G'): c=('C');
    	break;
    	case ('T'): c=('A');
    	break;
    	}
        d[i] = c;
    }

}

void align_region(sequence_manager * seq_man, Synteny_block * unlinked_sb, uint32_t kmer_size, dictionary_hash * dhw, Quickfrag * result){
    //Restore dictionary
    dhw->clear();

    //Kmer reading
    char curr_kmer[kmer_size], rev_kmer[kmer_size];
    uint64_t kmer_index = 0;
    char c;
    uint64_t precomputed_sizeofQuickfrag = sizeofQuickfrag();
    bool first = true;

    //To get hits
    Wordbucket * hit;
    Wordbucket ** hit_list;

    //Current pos in sequence
    uint64_t advanced_steps = 0;
    uint64_t i;

    //Results of alignment
    Quickfrag aligned_qf;
    int64_t curr_diag;
    Word align_word;

    //Current unlinked block
    Synteny_block * unlinked = unlinked_sb;

    //For all sequences in the unlinked synteny block
    while(unlinked != NULL){
        while(advanced_steps < unlinked->b->end){
            //Get nucleotide
            c = unlinked->b->genome->seq[advanced_steps];
            advanced_steps++;

            if(c == 'A' || c == 'C' || c == 'T' || c == 'G'){
                curr_kmer[kmer_index] = c;
                kmer_index++;
            }else{
                kmer_index = 0;
            }

            //Check if we have a kmer big enough
            if(kmer_index == kmer_size){

                //Insert word in dictionary
                hit_list = dhw->put_and_hit(curr_kmer, 'f', advanced_steps, unlinked->b);

                for(i=0;i<dhw->get_candidates_number();i++){
                    //For every hit candidate
                    hit = hit_list[i];

                    if(hit != NULL){

                        //We have a hit we should try an alignment
                        //only if the hit is not overlapping
                        //or it is overlapping but on different diagonal
                        if(first){

                            //There is no frag yet, so try first one
                            align_word.pos = advanced_steps;
                            align_word.strand = 'f';
                            align_word.b = unlinked->b;
                            if(hit->w.strand == 'f'){
                                alignment_from_hit(seq_man, &align_word, &hit->w, result, kmer_size);
                            }else{
                                alignment_from_hit_reverse(seq_man, &hit->w, &align_word, result, kmer_size);
                            }
                            first = false;
                        }else{

                            //There is already a frag, check for overlapping and diagonal
                            if(overlapped_words(result->x_start, result->x_start+result->t_len, advanced_steps-kmer_size-1, advanced_steps-1) != 0){
                                //If it is not overlapped
                                curr_diag = (int64_t) advanced_steps - (int64_t) hit->w.pos;
                                if(curr_diag != result->diag){
                                    //We can try new alignment
                                    align_word.pos = advanced_steps;
                                    align_word.strand = 'f';
                                    align_word.b->genome = unlinked->b->genome;
                                    if(hit->w.strand == 'f'){
                                        alignment_from_hit(seq_man, &align_word, &hit->w, &aligned_qf, kmer_size);
                                    }else{
                                        alignment_from_hit_reverse(seq_man, &hit->w, &align_word, &aligned_qf, kmer_size);
                                    }

                                    //Only copy if new alignment is better

                                    if(aligned_qf.sim > result->sim){
                                        memcpy(result, &aligned_qf, precomputed_sizeofQuickfrag);
                                    }
                                }
                            }
                        }

                    }
                }

                //Insert reversed word in dictionary
                strrev(curr_kmer, rev_kmer, kmer_size);

                hit_list = dhw->put_and_hit(rev_kmer, 'r', advanced_steps, unlinked->b);

                for(i=0;i<dhw->get_candidates_number();i++){
                    //For every hit candidate
                    hit = hit_list[i];

                    if(hit->w.strand == 'r') continue; //Skip since a 'r' - 'r' word is the same as 'f' - 'f'

                    if(hit != NULL){

                        //We have a hit we should try an alignment
                        //only if the hit is not overlapping
                        //or it is overlapping but on different diagonal
                        if(first){

                            //There is no frag yet, so try first one
                            align_word.pos = advanced_steps;
                            align_word.strand = 'r';
                            align_word.b = unlinked->b;
                            alignment_from_hit_reverse(seq_man, &align_word, &hit->w, result, kmer_size);
                            //printQuickfrag(&aligned_qf);
                            memcpy(result, &aligned_qf, precomputed_sizeofQuickfrag);
                            first = false;
                        }else{
                            //printf("on the else track\n");
                            //getchar();
                            //There is already a frag, check for overlapping and diagonal
                            if(overlapped_words(result->x_start, result->x_start+result->t_len, advanced_steps-kmer_size-1, advanced_steps-1) != 0){
                                //If it is not overlapped
                                curr_diag = (int64_t) advanced_steps - (int64_t) hit->w.pos;
                                if(curr_diag != result->diag){
                                    //We can try new alignment
                                    align_word.pos = advanced_steps;
                                    align_word.strand = 'r';
                                    align_word.b->genome = unlinked->b->genome;
                                    alignment_from_hit_reverse(seq_man, &align_word, &hit->w, &aligned_qf, kmer_size);

                                    if(aligned_qf.sim > result->sim){
                                        memcpy(result, &aligned_qf, precomputed_sizeofQuickfrag);
                                    }
                                }
                            }
                        }

                    }
                }
                //Displace
                memmove(&curr_kmer[0], &curr_kmer[1], kmer_size-1);
                memmove(&rev_kmer[0], &rev_kmer[1], kmer_size-1);
                kmer_index--;
            }

        }
        //Advance block
        unlinked = unlinked->next;
        kmer_index = 0;
    }





}

void read_words_from_synteny_block_and_align(sequence_manager * seq_man, Synteny_list * sbl, uint32_t kmer_size, dictionary_hash * dhw, Quickfrag ** qfmat, unsigned char ** qfmat_state){

    //@Important: A single hit between two sequences should be enough to align two blocks

    //Erase what we had for previous alignments
    dhw->clear();
    uint64_t i,j;
    for(i=0;i<seq_man->get_number_of_sequences();i++){
        for(j=0;j<seq_man->get_number_of_sequences();j++) qfmat_state[i][j] = 0;
    }
    Synteny_list * sbl_ptr = sbl;
    Synteny_block * sb_ptr;
    //To keep track of where we are
    uint64_t advanced_steps;

    //Kmer reading
    char curr_kmer[kmer_size], rev_kmer[kmer_size];
    uint64_t kmer_index = 0;
    char c;

    //To make things clearer
    Quickfrag * qf;
    Quickfrag aligned_qf;
    int64_t curr_diag;
    Word align_word;

    //And to speed things up
    uint64_t precomputed_sizeofQuickfrag = sizeofQuickfrag();
    Wordbucket * hit;
    Wordbucket ** hit_list;

    sb_ptr = sbl_ptr->sb;

    //For all blocks in the synteny
    while(sb_ptr != NULL){

        //Get next word of current block
        advanced_steps = sb_ptr->b->start;
        while(advanced_steps < sb_ptr->b->end){

            //Get nucleotide
            c = sb_ptr->b->genome->seq[advanced_steps];
            advanced_steps++;

            if(c == 'A' || c == 'C' || c == 'T' || c == 'G'){
                curr_kmer[kmer_index] = c;
                kmer_index++;
            }else{
                kmer_index = 0;
            }

            //Check if we have a kmer big enough
            if(kmer_index == kmer_size){

                //printf("Putting %s at %"PRIu64"\n", curr_kmer, advanced_steps);

                //Insert word in dictionary
                //printf("%s\n", curr_kmer);
                hit_list = dhw->put_and_hit(curr_kmer, 'f', advanced_steps, sb_ptr->b);

                for(i=0;i<dhw->get_candidates_number();i++){
                    //For every hit candidate
                    hit = hit_list[i];

                    if(hit != NULL){
                        //printf("Got hit and the state is %u. The synteny level is %"PRIu64"\n", (int)qfmat_state[sb_ptr->b->genome->id][hit->w.b->genome->id], sbl->synteny_level);
                        //printf("saux_id: %"PRIu64" and hit genome id: %"PRIu64"\n", sb_ptr->b->genome->id, hit->w.b->genome->id);
                        //getchar();
                        //In the sake of clarity
                        qf = &qfmat[sb_ptr->b->genome->id][hit->w.b->genome->id];

                        //We have a hit we should try an alignment
                        //only if the hit is not overlapping
                        //or it is overlapping but on different diagonal
                        if(qfmat_state[sb_ptr->b->genome->id][hit->w.b->genome->id] == 0){
                            //printf("BINGO!!!!!!!!!!!!!\n");

                            //There is no frag yet, so try first one
                            align_word.pos = advanced_steps;
                            align_word.strand = 'f';
                            align_word.b = sb_ptr->b;
                            if(hit->w.strand == 'f'){
                                alignment_from_hit(seq_man, &align_word, &hit->w, &aligned_qf, kmer_size);
                            }else{
                                alignment_from_hit_reverse(seq_man, &hit->w, &align_word, &aligned_qf, kmer_size);
                            }

                            //printQuickfrag(&aligned_qf);
                            qfmat_state[sb_ptr->b->genome->id][hit->w.b->genome->id] = 1;
                            qfmat_state[hit->w.b->genome->id][sb_ptr->b->genome->id] = 1;
                            memcpy(qf, &aligned_qf, precomputed_sizeofQuickfrag);
                            memcpy(&qfmat[hit->w.b->genome->id][sb_ptr->b->genome->id], &aligned_qf, precomputed_sizeofQuickfrag);

                        }else{
                            //printf("on the else track\n");
                            //getchar();
                            //There is already a frag, check for overlapping and diagonal
                            if(overlapped_words(qf->x_start, qf->x_start+qf->t_len, advanced_steps-kmer_size-1, advanced_steps-1) != 0){
                                //If it is not overlapped
                                curr_diag = (int64_t) advanced_steps - (int64_t) hit->w.pos;
                                if(curr_diag != qf->diag){
                                    //We can try new alignment
                                    align_word.pos = advanced_steps;
                                    align_word.strand = 'f';
                                    align_word.b->genome = sb_ptr->b->genome;
                                    if(hit->w.strand == 'f'){
                                        alignment_from_hit(seq_man, &align_word, &hit->w, &aligned_qf, kmer_size);
                                    }else{
                                        alignment_from_hit_reverse(seq_man, &hit->w, &align_word, &aligned_qf, kmer_size);
                                    }
                                    //printQuickfrag(&aligned_qf);
                                    //printf("HAPPENS AS WELL\n");
                                    //getchar();
                                    //getchar();

                                    //Only copy if new alignment is better

                                    if(aligned_qf.sim > qf->sim){
                                        memcpy(qf, &aligned_qf, precomputed_sizeofQuickfrag);
                                        memcpy(&qfmat[hit->w.b->genome->id][sb_ptr->b->genome->id], &aligned_qf, precomputed_sizeofQuickfrag);
                                    }
                                }
                            }
                        }

                    }
                }



                //Insert reversed word in dictionary
                strrev(curr_kmer, rev_kmer, kmer_size);

                hit_list = dhw->put_and_hit(rev_kmer, 'r', advanced_steps, sb_ptr->b);

                for(i=0;i<dhw->get_candidates_number();i++){
                    //For every hit candidate
                    hit = hit_list[i];

                    if(hit->w.strand == 'r') continue; //Skip since a 'r' - 'r' word is the same as 'f' - 'f'

                    if(hit != NULL){
                        //printf("Got hit and the state is %u. The synteny level is %"PRIu64"\n", (int)qfmat_state[sb_ptr->b->genome->id][hit->w.b->genome->id], sbl->synteny_level);
                        //printf("saux_id: %"PRIu64" and hit genome id: %"PRIu64"\n", sb_ptr->b->genome->id, hit->w.b->genome->id);
                        //getchar();
                        //In the sake of clarity
                        qf = &qfmat[sb_ptr->b->genome->id][hit->w.b->genome->id];

                        //We have a hit we should try an alignment
                        //only if the hit is not overlapping
                        //or it is overlapping but on different diagonal
                        if(qfmat_state[sb_ptr->b->genome->id][hit->w.b->genome->id] == 0){
                            //printf("BINGO!!!!!!!!!!!!!\n");

                            //There is no frag yet, so try first one
                            align_word.pos = advanced_steps;
                            align_word.strand = 'r';
                            align_word.b = sb_ptr->b;
                            alignment_from_hit_reverse(seq_man, &align_word, &hit->w, &aligned_qf, kmer_size);
                            //printQuickfrag(&aligned_qf);
                            qfmat_state[sb_ptr->b->genome->id][hit->w.b->genome->id] = 1;
                            qfmat_state[hit->w.b->genome->id][sb_ptr->b->genome->id] = 1;
                            memcpy(qf, &aligned_qf, precomputed_sizeofQuickfrag);
                            memcpy(&qfmat[hit->w.b->genome->id][sb_ptr->b->genome->id], &aligned_qf, precomputed_sizeofQuickfrag);

                        }else{
                            //printf("on the else track\n");
                            //getchar();
                            //There is already a frag, check for overlapping and diagonal
                            if(overlapped_words(qf->x_start, qf->x_start+qf->t_len, advanced_steps-kmer_size-1, advanced_steps-1) != 0){
                                //If it is not overlapped
                                curr_diag = (int64_t) advanced_steps - (int64_t) hit->w.pos;
                                if(curr_diag != qf->diag){
                                    //We can try new alignment
                                    align_word.pos = advanced_steps;
                                    align_word.strand = 'r';
                                    align_word.b->genome = sb_ptr->b->genome;
                                    alignment_from_hit_reverse(seq_man, &align_word, &hit->w, &aligned_qf, kmer_size);
                                    //printQuickfrag(&aligned_qf);
                                    //printf("HAPPENS AS WELL\n");
                                    //getchar();
                                    //getchar();

                                    //Only copy if new alignment is better

                                    if(aligned_qf.sim > qf->sim){
                                        memcpy(qf, &aligned_qf, precomputed_sizeofQuickfrag);
                                        memcpy(&qfmat[hit->w.b->genome->id][sb_ptr->b->genome->id], &aligned_qf, precomputed_sizeofQuickfrag);
                                    }
                                }
                            }
                        }

                    }
                }




                //Displace
                memmove(&curr_kmer[0], &curr_kmer[1], kmer_size-1);
                memmove(&rev_kmer[0], &rev_kmer[1], kmer_size-1);
                kmer_index--;
            }

        }

        //Advance block
        sb_ptr = sb_ptr->next;
        kmer_index = 0;
    }
    fprintf(stdout, "Synteny level: %"PRIu64"\n", sbl->synteny_level);
    printSyntenyBlock(sbl->sb);
    for(i=0;i<seq_man->get_number_of_sequences();i++){
        for(j=0;j<seq_man->get_number_of_sequences();j++){
            if(qfmat_state[i][j] == 1) qfmat[i][j].sim = 100.0 - qfmat[i][j].sim;
        }
    }
    printQuickFragMatrix(qfmat, qfmat_state, seq_man->get_number_of_sequences());
    //getchar();

}

/*
struct alignment_arguments{
    char * seq_A;
    uint64_t start_A;
    uint64_t end_A;
    char * seq_B;
    uint64_t start_B;
    uint64_t end_B;
    int64_t iGap;
    int64_t eGap;
    struct cell * mc;
    struct cell * f0;
    struct cell * f1;
    struct cell * alignment;
    struct cell * alignment_reverse;
    char * seq_for_reverse;
};
*/

void * compute_NW_on_pthreads(void * a){

    alignment_arguments * aa = (alignment_arguments *) a;

    *aa->alignment = NWscore2rows(aa->seq_A, 0, aa->end_A, aa->seq_B, 0, aa->end_B, aa->iGap, aa->eGap, aa->mc, aa->f0, aa->f1);
    memcpy(aa->seq_for_reverse, aa->seq_B, aa->end_B);
    inplace_reverse_and_complement(aa->seq_for_reverse, aa->end_B - 2);
    *aa->alignment_reverse = NWscore2rows(aa->seq_A, 0, aa->end_A, aa->seq_for_reverse, 0, aa->end_B, aa->iGap, aa->eGap, aa->mc, aa->f0, aa->f1);

    return NULL;
}


void fill_quickfrag_matrix_NW(sequence_manager * seq_man, char ** seq_for_reverse, Synteny_list * sbl, Quickfrag ** qfmat, unsigned char ** qfmat_state, int iGap, int eGap, struct cell ** mc, struct cell ** f0, struct cell ** f1, pthread_t * threads){
    Synteny_block * sb_ptr = sbl->sb;
    Synteny_block * sb_ptr_intern;
    uint64_t n_sequences = seq_man->get_number_of_sequences();
    struct cell alignment[n_sequences];
    struct cell alignment_reverse[n_sequences];
    Quickfrag qf[n_sequences];
    unsigned char active_threads[n_sequences];
    alignment_arguments aa[n_sequences];

    uint64_t i,j;
    for(i=0;i<seq_man->get_number_of_sequences();i++){
        active_threads[i] = 0;
        for(j=0;j<seq_man->get_number_of_sequences();j++){

            qfmat_state[i][j] = 0;
        }
    }

    uint64_t index = 0;
    while(sb_ptr != NULL){
        sb_ptr_intern = sb_ptr->next;
        while(sb_ptr_intern != NULL){


            #ifdef VERBOSE
            printf("aligning :"); printBlock(sb_ptr->b); printBlock(sb_ptr_intern->b);
            #endif
            /*
            for(uint64_t z=0;z<(sb_ptr->b->end - sb_ptr->b->start); z++){
                printf("%c", sb_ptr->b->genome->seq[sb_ptr->b->start+z]);
                getchar();
            }
            */


            aa[index].seq_A = &sb_ptr->b->genome->seq[sb_ptr->b->start];
            aa[index].s_x = sb_ptr->b->genome;
            aa[index].start_A = 0;
            aa[index].end_A = sb_ptr->b->end - sb_ptr->b->start;
            aa[index].seq_B = &sb_ptr_intern->b->genome->seq[sb_ptr_intern->b->start];
            aa[index].s_y = sb_ptr_intern->b->genome;
            aa[index].start_B = 0;
            aa[index].end_B = sb_ptr_intern->b->end - sb_ptr_intern->b->start;
            aa[index].iGap = (int64_t) iGap;
            aa[index].eGap = (int64_t) eGap;
            printf("I am index: %"PRIu64"\n", index);
            aa[index].mc = mc[index];
            aa[index].f0 = f0[index];
            aa[index].f1 = f1[index];
            aa[index].alignment = &alignment[index];
            aa[index].alignment_reverse = &alignment_reverse[index];
            aa[index].seq_for_reverse = seq_for_reverse[index];

            int error;
            if( 0 != (error = pthread_create(&threads[index], NULL, compute_NW_on_pthreads, (void *) &aa[index] ))){
                fprintf(stdout, "Thread %"PRIu64" returned %d:", index, error); terror("Could not launch");
            }
            /*
            compute_NW_on_pthreads(&sb_ptr->b->genome->seq[sb_ptr->b->start], 0, sb_ptr->b->end - sb_ptr->b->start, &sb_ptr_intern->b->genome->seq[sb_ptr_intern->b->start], 0, sb_ptr_intern->b->end - sb_ptr_intern->b->start, (int64_t) iGap, (int64_t) eGap, mc[index], f0[index], f1[index], &alignment[index], &alignment_reverse[index], seq_for_reverse[index]);
            */

            /*
            alignment = NWscore2rows(&sb_ptr->b->genome->seq[sb_ptr->b->start], 0, sb_ptr->b->end - sb_ptr->b->start, &sb_ptr_intern->b->genome->seq[sb_ptr_intern->b->start], 0, sb_ptr_intern->b->end - sb_ptr_intern->b->start, (int64_t) iGap, (int64_t) eGap, mc, f0, f1);
            // problem here with \0 (solved)
            memcpy(&seq_for_reverse[sb_ptr_intern->b->start], &sb_ptr_intern->b->genome->seq[sb_ptr_intern->b->start], sb_ptr_intern->b->end - sb_ptr_intern->b->start);


            inplace_reverse_and_complement(&seq_for_reverse[sb_ptr_intern->b->start], sb_ptr_intern->b->end - sb_ptr_intern->b->start -2);
            alignment_reverse = NWscore2rows(&sb_ptr->b->genome->seq[sb_ptr->b->start], 0, sb_ptr->b->end - sb_ptr->b->start, &seq_for_reverse[sb_ptr_intern->b->start], 0, sb_ptr_intern->b->end - sb_ptr_intern->b->start, iGap, eGap, mc, f0, f1);
            */
            active_threads[index] = 1;
            index++; // The index controls the pthread that is launched
            sb_ptr_intern = sb_ptr_intern->next;
        }
        sb_ptr = sb_ptr->next;
        // Make wait/join here
        for(i=0;i<n_sequences;i++){
            if(active_threads[i] == 1) pthread_join(threads[i], NULL);
        }
        for(i=0;i<index;i++){
            // Paste data into quickfrag matrix
            if(alignment_reverse[i].ident > alignment[i].ident) alignment[i] = alignment_reverse[i];

            qf[i].x_start = alignment[i].xs;
            qf[i].y_start = alignment[i].ys;
            qf[i].t_len = MAX(alignment[i].xe - alignment[i].xs, alignment[i].ye - alignment[i].ys);
            qf[i].diag = 0; // Not needed
            //qf[i].sim = 100 - (long double) (MIN(100, (uint64_t)alignment[i].ident * 100)) / (long double) qf[i].t_len;
            qf[i].sim = 100 - (long double) MIN(100, (  (alignment[i].ident * 100)/(long double) qf[i].t_len  ));
            qf[i].x = aa[i].s_x;
            qf[i].y = aa[i].s_y;

            #ifdef VERBOSE
            printf("Best is:\n");
            printCell(&alignment[i]); printf("Has len: %"PRIu64" and sim %Le\n", MAX(alignment[i].xe - alignment[i].xs, alignment[i].ye - alignment[i].ys), qf[i].sim);
            //getchar();
            #endif

            memcpy(&qfmat[qf[i].x->id][qf[i].y->id], &qf[i], sizeofQuickfrag());
            memcpy(&qfmat[qf[i].y->id][qf[i].x->id], &qf[i], sizeofQuickfrag());
            qfmat_state[qf[i].x->id][qf[i].y->id] = 1;
            qfmat_state[qf[i].y->id][qf[i].x->id] = 1;
        }
        index = 0; // Restart index to align more
    }


    /*
    // Restart pthread index
    index = 0;
    sb_ptr = sbl->sb;
    while(sb_ptr != NULL){
        sb_ptr_intern = sb_ptr->next;
        while(sb_ptr_intern != NULL){



            if(alignment_reverse[index].ident > alignment[index].ident) alignment[index] = alignment_reverse[index];
            qf[index].x_start = alignment[index].xs;
            qf[index].y_start = alignment[index].ys;
            qf[index].t_len = alignment[index].xe - alignment[index].xs;
            qf[index].diag = 0; // Not needed
            qf[index].sim = 100 - (long double)(MIN(100, (uint64_t)alignment[index].ident * 100)) / (long double) qf[index].t_len;
            qf[index].x = sb_ptr->b->genome;
            qf[index].y = sb_ptr_intern->b->genome;

            #ifdef VERBOSE
            printf("Best is:\n");
            printCell(&alignment[index]);
            getchar();
            #endif

            memcpy(&qfmat[sb_ptr->b->genome->id][sb_ptr_intern->b->genome->id], &qf[index], sizeofQuickfrag());
            memcpy(&qfmat[sb_ptr_intern->b->genome->id][sb_ptr->b->genome->id], &qf[index], sizeofQuickfrag());
            qfmat_state[sb_ptr->b->genome->id][sb_ptr_intern->b->genome->id] = 1;
            qfmat_state[sb_ptr_intern->b->genome->id][sb_ptr->b->genome->id] = 1;
            sb_ptr_intern = sb_ptr_intern->next;

            // Controls pthread
            index++;
        }
        sb_ptr = sb_ptr->next;
    }
    */


    #ifdef VERBOSE
    printQuickFragMatrix(qfmat, qfmat_state, seq_man->get_number_of_sequences());
    #endif
}

Slist * UPGMA_joining_clustering(Quickfrag ** M, double ** submat, unsigned char ** qfmat_state, uint64_t N, memory_pool * mp){


    //printQuickFragMatrix(M, qfmat_state, 2);

    uint64_t i, j, i_p, j_p;
    double f_min = -1;
    uint64_t i_min = 0, j_min = 1;
    uint64_t nodes_in_dendro = 1;

    Slist ** dendrogram = (Slist **) mp->request_bytes(N*sizeof(Slist*));

    uint64_t * dendro_lengths = (uint64_t *) mp->request_bytes(N*sizeof(uint64_t));

    unsigned char * skip_i = (unsigned char *) mp->request_bytes(N*sizeof(unsigned char));
    for(i=0;i<N;i++){
        dendrogram[i] = NULL;
        skip_i[i] = 0;
        dendro_lengths[i] = 1;
    }


    //Generate submatrix, since some rows or columns will be null
    i_p = 0; j_p = 0;
    for(i=0;i<N;i++){
        j_p = 0;
        for(j=0;j<N;j++){
            if(qfmat_state[i][j] == 1){
                // Copy value
                if(i_p == j_p) j_p++;
                submat[i_p][j_p] = M[i][j].sim;
                if(dendrogram[i_p] == NULL){
                    dendrogram[i_p] = (Slist *) mp->request_bytes(sizeofSlist());
                    if(M[i][j].x->id == i) dendrogram[i_p]->s = M[i][j].x;
                    if(M[i][j].y->id == i) dendrogram[i_p]->s = M[i][j].y;
                    dendrogram[i_p]->next = NULL;

                }
                j_p++;
            }
        }

        if(j_p > 0) i_p++;
    }

    /*
    for(i=0;i<i_p;i++){
        printf("ID: %"PRIu64"\n", dendrogram[i]->s->id);
    }
    */

    j_p = i_p;
    while(nodes_in_dendro < i_p){
        #ifdef VERBOSE
        printUnstatedDoubleMatrix(submat, i_p, skip_i);
        #endif

        //Find smallest
        for(i=0;i<i_p;i++){
            //Jump those that have been removed
            if(skip_i[i] == 1) continue;

            for(j=0;j<i;j++){
                //Jump removed & Find minimum
                if(skip_i[j] == 0 && (f_min == -1 || f_min >= submat[i][j])){
                    f_min = submat[i][j];
                    i_min = i;
                    j_min = j;
                }


            }
        }

        //Join minimums
        //printf("Min on %"PRIu64", %"PRIu64"\n", i_min, j_min);


        if( dendrogram[j_min]->next == NULL || dendrogram[i_min]->next == NULL){
            //Add one makes a tuple
            Slist * aux = dendrogram[i_min];
            while(aux->next != NULL){
                aux = aux->next;
            }

            aux->next = dendrogram[j_min];
            //last_d->next = dendrogram[j_min];
            //last_d = dendrogram[j_min];
            dendro_lengths[i_min]++;
        }else{
            Slist * null_seq = (Slist *) mp->request_bytes(sizeofSlist());
            null_seq->s = NULL;
            null_seq->next = dendrogram[j_min];

            Slist * aux = dendrogram[i_min];
            while(aux->next != NULL){
                aux = aux->next;
            }
            aux->next = null_seq;
            dendro_lengths[i_min] += dendro_lengths[j_min];

        }




        //Remove existing rows and cols
        skip_i[j_min] = 1;

        //Update new distances
        for(i=0;i<i_p;i++){
            if(skip_i[i] == 0){
                submat[i_min][i] = ((submat[i_min][i] + submat[j_min][i]) * dendro_lengths[i_min]) /2;
                submat[i][i_min] = submat[i_min][i];
            }
        }
        /*
        for(i=0;i<i_p;i++){
            if(skip_i[i] == 0){
                submat[j_min][i] = ((submat[i_min][i] + submat[j_min][i]) * dendro_lengths[i_min]) /2;
            }
        }
        */

        nodes_in_dendro++;


        f_min = -1; //Restart minimum search

    }

    //Print clusters
    #ifdef VERBOSE
    printf("Phylogenetic Clustering: \n");printDendrogramList(dendrogram[i_min]);
    printf("\n");
    #endif


    return dendrogram[i_min];

}


uint64_t generate_multiple_alignment(arguments_multiple_alignment * arg_mul_al){
    /*
    struct arguments_multiple_alignment{
    sequence_manager * seq_man;
    char ** seq_for_reverse;
    Synteny_list * sbl;
    Quickfrag ** qfmat;
    unsigned char ** qfmat_state;
    int iGap;
    int eGap;
    struct cell ** mc;
    struct cell ** f0;
    struct cell ** f1;
    pthread_t * threads;
    double ** submat;
    memory_pool * mp;
    // For full NW
    char * aux_dummy_sequence; // The one that is not used in the backtrackings
    uint64_t * sequence_ids;
    char ** recon_X;
    char ** recon_Y;
    char ** recon_Z;
    char ** seq_X;
    char ** seq_Y;
    char ** seq_Z;
    int64_t * cell_path_y;
    struct positioned_cell * mc_f;
    struct cell_f ** table_f;
    char * writing_buffer_alignment;
    long double window;

    */
    // First generate the guidance tree to tell how to pair sequences
    //fill_quickfrag_matrix_NW(sequence_manager * seq_man, char ** seq_for_reverse, Synteny_list * sbl, Quickfrag ** qfmat, unsigned char ** qfmat_state, int iGap, int eGap, struct cell ** mc, struct cell ** f0, struct cell ** f1, pthread_t * threads){
    fill_quickfrag_matrix_NW(arg_mul_al->seq_man, arg_mul_al->seq_for_reverse, arg_mul_al->sbl, arg_mul_al->qfmat, arg_mul_al->qfmat_state, arg_mul_al->iGap, arg_mul_al->eGap, arg_mul_al->mc, arg_mul_al->f0, arg_mul_al->f1, arg_mul_al->threads);
    Slist * tree = UPGMA_joining_clustering(arg_mul_al->qfmat, arg_mul_al->submat, arg_mul_al->qfmat_state, arg_mul_al->seq_man->get_number_of_sequences(), arg_mul_al->mp);

    uint64_t current_sequences_in_X = 0;
    uint64_t current_sequences_in_Y = 0;
    //uint64_t current_sequences_in_Z = 0;
    uint64_t final_len;
    memset(arg_mul_al->sequence_ids, 0, arg_mul_al->n_sequences * sizeof(uint64_t));

    // Link IDs to synteny blocks
    uint64_t starting_pos[arg_mul_al->n_sequences];
    uint64_t ending_pos[arg_mul_al->n_sequences];
    Synteny_block * sb_ptr = arg_mul_al->sbl->sb;
    while(sb_ptr != NULL){
        starting_pos[sb_ptr->b->genome->id] = sb_ptr->b->start;
        ending_pos[sb_ptr->b->genome->id] = sb_ptr->b->end;
        sb_ptr = sb_ptr->next;
    }

    // Traverse tree and make pairwise alignments
    BasicAlignment ba;
    Slist * d = tree;
    while(d != NULL){
        if(d->s == NULL) printf(" ] "); else printf(" %"PRIu64" ", d->s->id);
        //getchar();

        if(d->s != NULL){

            if(current_sequences_in_X == 0){
                memcpy(&arg_mul_al->seq_X[current_sequences_in_X][0], &d->s->seq[starting_pos[d->s->id]], ending_pos[d->s->id] - starting_pos[d->s->id]);
                arg_mul_al->sequence_ids[current_sequences_in_X] = d->s->id;
                current_sequences_in_X++;
            }else{
                memcpy(&arg_mul_al->seq_Y[current_sequences_in_Y][0], &d->s->seq[starting_pos[d->s->id]], ending_pos[d->s->id] - starting_pos[d->s->id]);
                current_sequences_in_Y++;
            }

            if(current_sequences_in_Y == 1){

                //printf("I am going to gap-bounded align: \nX[-1]%s\nY[0]%s\n", arg_mul_al->seq_X[current_sequences_in_X-1], arg_mul_al->seq_Y[0]);


                build_unanchored_alignment(arg_mul_al->cell_path_y, arg_mul_al->seq_X, arg_mul_al->seq_Y, current_sequences_in_X, current_sequences_in_Y);
                uint64_t xlen = get_max_length_of_sequences(arg_mul_al->seq_X, current_sequences_in_X);
                uint64_t ylen = get_max_length_of_sequences(arg_mul_al->seq_Y, current_sequences_in_Y);

                final_len = build_multiple_alignment(arg_mul_al->recon_X, arg_mul_al->recon_Y, arg_mul_al->seq_X, arg_mul_al->seq_Y, current_sequences_in_X, current_sequences_in_Y, arg_mul_al->table_f, arg_mul_al->mc_f, arg_mul_al->writing_buffer_alignment, xlen, ylen, arg_mul_al->cell_path_y, &arg_mul_al->window, arg_mul_al->iGap, arg_mul_al->eGap, &ba, arg_mul_al->aux_dummy_sequence);

                // After aligning group X with group Y, move Y to X


                memcpy(&arg_mul_al->seq_X[current_sequences_in_X][0], &arg_mul_al->seq_Y[0][0], final_len);
                arg_mul_al->sequence_ids[current_sequences_in_X] = d->s->id;
                current_sequences_in_X++;
                //printf("CE IMPRESSIVE -> %s\n", arg_mul_al->seq_Y[0]);




                current_sequences_in_Y = 0;
            }



        }else{

        }

        d = d->next;
    }
    // Finished
    printf("Done with multiple alignment\n");
    /*
    for(uint64_t i=0; i<current_sequences_in_X; i++){
        printf("X: %"PRIu64" -> %s\n", arg_mul_al->sequence_ids[i], arg_mul_al->seq_X[i]);
    }
    */
    //getchar();

    return final_len;

}


void find_event_location(Slist * dendrogram, Event e, void * data, Event_handling * genomes_affected, Slist * last_dendro_pointer){
    Slist * d = dendrogram;
    Slist * prev, * next;

    if(last_dendro_pointer != NULL) d = last_dendro_pointer;

    while(d != NULL){
        //if(d->s == NULL) printf(" ] "); else printf(" %"PRIu64" ", d->s->id);
        prev = d;
        d = d->next;
        if(d != NULL) next = d->next;

        if(prev != NULL && d != NULL && next != NULL){
            if(prev->s != NULL && d->s != NULL && next->s != NULL){

                last_dendro_pointer = next;

                if(e == inversion){
                    //check for strands

                    strand_matrix * sm_B = (strand_matrix *) data;
                    printf("%u %u %u %u\n", sm_B->get_strands(next->s->id, d->s->id), sm_B->get_strands(d->s->id, prev->s->id), sm_B->get_strands(next->s->id, prev->s->id), sm_B->get_strands(prev->s->id, d->s->id));
                    sm_B->print_strand_matrix();
                    if(sm_B->get_strands(next->s->id, d->s->id) == sm_B->get_strands(d->s->id, prev->s->id)){
                        //A reversion happened in 'prev'
                        genomes_affected->genomes_affected[prev->s->id] = true;
                    }
                    if(sm_B->get_strands(next->s->id, prev->s->id) == sm_B->get_strands(prev->s->id, d->s->id)){
                        //A reversion happened in 'd'
                        genomes_affected->genomes_affected[d->s->id] = true;
                    }
                    genomes_affected->type_of_event = inversion;

                }
                else if(e == transposition){
                    strand_matrix * sm_B = (strand_matrix *) data;
                    sm_B->print_strand_matrix();
                    if(sm_B->get_strands(next->s->id, d->s->id) == sm_B->get_strands(d->s->id, prev->s->id)){
                        //A reversion happened in 'prev'
                        genomes_affected->genomes_affected[prev->s->id] = true;
                    }
                    if(sm_B->get_strands(next->s->id, prev->s->id) == sm_B->get_strands(prev->s->id, d->s->id)){
                        //A reversion happened in 'd'
                        genomes_affected->genomes_affected[d->s->id] = true;
                    }
                    genomes_affected->type_of_event = transposition;

                }
                else if(e == indel){

                    // TODO this should be remade


                    Indel_handling * indel_hd = (Indel_handling *) data;



                    int64_t d_next_prev = (int64_t) labs((int64_t) indel_hd->bp_lengths[next->s->id] - (int64_t) indel_hd->bp_lengths[prev->s->id]);
                    //int64_t d_next_d = (int64_t) labs((int64_t) indel_hd->bp_lengths[next->s->id] - (int64_t) indel_hd->bp_lengths[d->s->id]);
                    #define RELATIVELY_SIMILAR 0.15
                    #define VERY_SIMILAR 0.9


                    //
                    printf("printing values of importance:\n");
                    printf("%"PRId64" (d_next_prev) %Le (C2-1) %Le (C1-1) %Le (C2-2) %Le (C1-2)\n", d_next_prev, (long double) indel_hd->bp_lengths[d->s->id] / (long double) indel_hd->bp_lengths[next->s->id], (long double) indel_hd->bp_lengths[prev->s->id] / (long double) indel_hd->bp_lengths[next->s->id], (long double) indel_hd->bp_lengths[prev->s->id] / (long double) indel_hd->bp_lengths[next->s->id],  (long double) indel_hd->bp_lengths[d->s->id] / (long double) indel_hd->bp_lengths[next->s->id]);

                    if(d_next_prev >= 0){
                        // possible deletion
                        if((long double) indel_hd->bp_lengths[d->s->id] / (long double) indel_hd->bp_lengths[next->s->id] > VERY_SIMILAR){
                            if(! ((long double) indel_hd->bp_lengths[prev->s->id] / (long double) indel_hd->bp_lengths[next->s->id] > VERY_SIMILAR)){
                                if(indel_hd->bp_lengths[prev->s->id] > indel_hd->bp_lengths[next->s->id]){
                                    // There is insertion
                                    printf("Def an insertion in %"PRIu64"\n", prev->s->id);
                                    genomes_affected->type_of_event = insertion;
                                    genomes_affected->genomes_affected[prev->s->id] = true;
                                }else{
                                    // There is deletion
                                    genomes_affected->type_of_event = deletion;
                                    printf("Def a deletion in %"PRIu64"\n", prev->s->id);
                                    genomes_affected->genomes_affected[prev->s->id] = true;
                                }

                                return;
                            }
                        }
                    }else{
                        if((long double) indel_hd->bp_lengths[prev->s->id] / (long double) indel_hd->bp_lengths[next->s->id] > VERY_SIMILAR){
                            if(! ((long double) indel_hd->bp_lengths[d->s->id] / (long double) indel_hd->bp_lengths[next->s->id] > VERY_SIMILAR)){
                                if(indel_hd->bp_lengths[prev->s->id] > indel_hd->bp_lengths[next->s->id]){
                                    // There is insertion
                                    genomes_affected->type_of_event = insertion;
                                    printf("Def an insertion in %"PRIu64"\n", d->s->id);
                                    genomes_affected->genomes_affected[d->s->id] = true;
                                }else{
                                    // There is deletion
                                    genomes_affected->type_of_event = deletion;
                                    printf("Def a deletion in %"PRIu64"\n", d->s->id);
                                    genomes_affected->genomes_affected[d->s->id] = true;
                                }

                                return;
                            }
                        }
                    }
                }


                /*
                switch(e){
                    case inversion:
                    {
                        //check for strands

                        strand_matrix * sm_B = (strand_matrix *) data;
                        printf("%u %u %u %u\n", sm_B->get_strands(next->s->id, d->s->id), sm_B->get_strands(d->s->id, prev->s->id), sm_B->get_strands(next->s->id, prev->s->id), sm_B->get_strands(prev->s->id, d->s->id));
                        sm_B->print_strand_matrix();
                        if(sm_B->get_strands(next->s->id, d->s->id) == sm_B->get_strands(d->s->id, prev->s->id)){
                            //A reversion happened in 'prev'
                            genomes_affected[prev->s->id] = true;
                        }
                        if(sm_B->get_strands(next->s->id, prev->s->id) == sm_B->get_strands(prev->s->id, d->s->id)){
                            //A reversion happened in 'd'
                            genomes_affected[d->s->id] = true;
                        }

                    }
                    break;
                    case transposition:
                    {
                        strand_matrix * h = (strand_matrix *) data;
                        h->print_strand_matrix();
                        if((h->get_strands(next->s->id, next->s->id) == h->get_strands(d->s->id, d->s->id)) == (h->get_strands(d->s->id, d->s->id) == h->get_strands(prev->s->id, prev->s->id))){
                            //A transposition happened in 'prev'
                            genomes_affected[prev->s->id] = true;
                        }
                        if((h->get_strands(next->s->id, next->s->id) == h->get_strands(prev->s->id, prev->s->id)) == (h->get_strands(prev->s->id, prev->s->id) == h->get_strands(d->s->id, d->s->id))){
                            //A transposition  happened in 'd'
                            genomes_affected[d->s->id] = true;
                        }

                    }

                    case duplication: {}
                    break;
                    case insertion: {}
                    break;
                    case deletion: {}
                    break;

                }
                */
            }else{
                //Collapse
                // TODO

            }
        }
    }
}

int compare_distances_indel(const void * a, const void * b){

    if( *(uint64_t *) a > *(uint64_t *) b) return 1;
    if(*(uint64_t *) a == *(uint64_t *) b) return 0;
    return -1;
}

long double median_from_vector(uint64_t * v, uint64_t l){
    if(l % 2 != 0){
        //its not pair
        return (long double) v[l/2];
    }else{
        //Pair
        return ((long double) v[(l-1)/2] + (long double) v[(l-1)/2 + 1])/2;
    }
}

void print_memory_usage(){
    //fprintf(stdout, "[INFO] Current RAM usage: %"PRIu64" (MB)\n", total_bytes_in_use/(1024*1024));
    fprintf(stdout, "[INFO] Current RAM usage: %"PRIu64" bytes, which are %"PRIu64" Megabytes. \n", total_bytes_in_use, total_bytes_in_use/(1024*1024));
}

void write_header(FILE * f, uint64_t sx_len, uint64_t sy_len){
    fprintf(f, "CSV file\n");
    fprintf(f, "[Jul.15 -- < bitlab - Departamento de Arquitectura de Computadores >\n");
    fprintf(f, "SeqX filename	: DATA1.dat\n");
    fprintf(f, "SeqY filename	: DATA2.dat\n");
    fprintf(f, "SeqX name	: S1\n");
    fprintf(f, "SeqY name	: S2\n");
    fprintf(f, "SeqX length	: %"PRIu64"\n", sx_len);
    fprintf(f, "SeqY length	: %"PRIu64"\n", sy_len);
    fprintf(f, "Min.fragment.length	: 0\n");
    fprintf(f, "Min.Identity	: 0.0\n");
    fprintf(f, "Total hits	: 0\n");
    fprintf(f, "Total hits (used)	: 0\n");
    fprintf(f, "Total fragments	: 0\n");
    fprintf(f, "Total CSBs:	: 0\n");
    fprintf(f, "Frag/CSB,xStart,yStart,xEnd,yEnd,strand,block,length,score,ident,similarity,identity,geneX,geneY\n");
    fprintf(f, "========================================================\n");
}


void save_frags_from_block(FILE * out_file, Block * b, heuristic_sorted_list * hsl) {
  Frags_list * fl = b->f_list;
  FragFile f;
  uint64_t i, t;
  int v;

  i = 0;
  while (fl != NULL) {
    f = *fl->f;
    hsl->insert(i++, abs(f.xStart - f.yStart));
    fl = fl->next;
  }
  t = hsl->get_first();
  v = i == 1 ? 0 : 1;
  i = 0;
  fl = b->f_list;
  while (fl != NULL){
    f = *fl->f;
    if (v != 0)
      v = i == t ? 1 : 2;
    fprintf(out_file, "Frag,%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",%c,", f.xStart, f.yStart, f.xEnd, f.yEnd, f.strand);
    fprintf(out_file, "%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",%.2f,%.2f,0,%d\n", b->id, f.length, f.score, f.ident, f.similarity, ((float)f.ident * 100 / (float)f.length), v);
    fl = fl->next;
    i++;
  }
  hsl->clear();
}

/*int contains_repetitions(Synteny_block * ptr_sb, uint64_t seq1_label, uint64_t seq2_label) {
  uint8_t repetitions;
  Block * ptr_b;

  repetitions = 0b00;

  while (ptr_sb != NULL) {
    ptr_b = ptr_sb->b;
    if (ptr_b->genome->id == seq1_label) {
      if (repetitions & 0b01) return 1;
      else repetitions |= 0b01;
    } else if (ptr_b->genome->id == seq2_label) {
      if (repetitions & 0b10) return 1;
      else repetitions |= 0b10;
    }
    ptr_sb = ptr_sb->next;
  }

  return 0;
}*/

void save_frag_pair(FILE * out_file, uint64_t seq1_label, uint64_t seq2_label, sequence_manager * seq_mngr, Synteny_list * sbl) {
  Sequence * seq1, * seq2;
  Synteny_list * ptr_sbl = sbl;
  Synteny_block * ptr_sb;
  Block * ptr_b;
  heuristic_sorted_list hsl;
  uint64_t gen_id;
  //int repetitions;

  seq1 = seq_mngr->get_sequence_by_label(seq1_label);
  seq2 = seq_mngr->get_sequence_by_label(seq2_label);

  write_header(out_file, seq1->len, seq2->len);
  while (ptr_sbl != NULL){
    ptr_sb = ptr_sbl->sb;
    // Check if SB contains repetitions
    //repetitions = contains_repetitions(ptr_sb, seq1_label, seq2_label);

    // Write actual values
    while (ptr_sb != NULL){
        ptr_b = ptr_sb->b;
        gen_id = ptr_b->genome->id;
        if (gen_id == seq1_label){
            save_frags_from_block(out_file, ptr_b, &hsl);
        }
        ptr_sb = ptr_sb->next;
    }

    ptr_sbl = ptr_sbl->next;
  }
}

void save_all_frag_pairs(char * out_file_base_path, sequence_manager * seq_manager, Synteny_list * sbl){
  // TODO better descriptions
  // Iterators
  uint64_t i, j;
  // Path to svc
  char out_file_name[512];
  // Number of sequences involved
  uint64_t n_seq;
  //
  FILE * out_file;

  n_seq = seq_manager->get_number_of_sequences();

  // For each pair of sequences
  for(i=0;i<n_seq;i++) for(j=i+1;j<n_seq;j++){
    sprintf(out_file_name, "%s-%"PRIu64"-%"PRIu64".csv", out_file_base_path, i, j);
    out_file = fopen64(out_file_name, "wt");
    if (out_file == NULL) continue;
    save_frag_pair(out_file, i, j, seq_manager, sbl);
    fclose(out_file);
  }

}
