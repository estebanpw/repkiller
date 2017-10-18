#define __STDC_FORMAT_MACROS

#include <stdio.h>
#include <cstdlib>
#include <string>
#include <ctype.h>
#include <pthread.h>
#include <float.h>
#include "structs.h"
#include "commonFunctions.h"
#include "comparisonFunctions.h"

#define LOGT(s) fprintf(stdout, "[INFO] %s T = %e\n", s, (double)(end-begin)/CLOCKS_PER_SEC);
#define LOGI(s) fprintf(stdout, "[INFO] %s\n", s);

int main(int argc, char * argv []) {
        // Clocks to measure time
        clock_t begin, end;
        // The sequences manager to store ids, lengths, etc
        sequence_manager seq_manager;
        std::string out_file_base_path, multifrags_path;
        std::vector<std::string> args(argc);
        args.assign(argv, argv + argc);

        // Open frags file, lengths file and output files %%%%%%%%%%%%%%%%%%%%%%%%%%%
        FILE * frags_file = nullptr, * lengths_file = nullptr;
        init_args(args, frags_file, out_file_base_path, multifrags_path);

        // Concat .lengths to path of multifrags %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        std::string path_lengths = multifrags_path + ".lengths";
        lengths_file = fopen64(path_lengths.c_str(), "rb");
        if(lengths_file == nullptr) terror("Could not open input lengths file");

        // Load fragments into array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        LOGI("Loading fragments into memory...");
        begin = clock();
        FragmentsDatabase frag_db(frags_file, lengths_file, seq_manager);
        end = clock();
        LOGT("Fragments loaded into memory succesfully!");
        fclose(frags_file);
        fclose(lengths_file);

        // Generate fragment groups %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        LOGI("Generating fragment groups...");
        begin = clock();
        FGList * efrags_groups = new FGList();
        generate_fragment_groups(frag_db, efrags_groups, seq_manager, 0.853, 10.0);
        end = clock();
        LOGT("Fragment groups generated succesfully!");

        // Save results in file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        LOGI("Saving results into secondary memory...");
        begin = clock();
        save_all_frag_pairs(out_file_base_path, seq_manager, *efrags_groups);
        end = clock();
        LOGT("Fragments saved into csv file.");

}
