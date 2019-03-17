#ifndef __GEMINO_SEQUENCER_H
#define __GEMINO_SEQUENCER_H



#include <RcppArmadillo.h>
#include <vector>  // vector class
#include <string>  // string class

#include "gemino_types.hpp"  // uint32

using namespace Rcpp;



namespace sequencer {

// Goes from character (coerced to integer) to index from 0:3 (4 for non-nucleotide)
std::vector<uint8> nt_map = {
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
    4,2,4,1,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,4,0,4,4,4,4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4};

std::vector<std::string> mm_nucleos = {"CAG", "TAG", "TCG", "TCA", "NNN"};

}




/*
 Way to store info for the sequence identifier line in FASTQ file.
 It's organized as such:
 @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>
    <read>:<is filtered>:<control number>:<sample number>
 (where the newline and tab are ' ')
 (from https://help.basespace.illumina.com/articles/descriptive/fastq-files/)
 Example: "@SIM:1:FCX:1:15:6329:1045 1:N:0:2"
 */
struct SequenceIdentifierInfo {

    // (Numerical) Read number. 1 can be single read or read 2 of paired-end:
    uint32 read;

    SequenceIdentifierInfo() {}
    SequenceIdentifierInfo(
        const std::string& instrument,
        const uint32& run_number,
        const std::string& flowcell_ID,
        const uint32& lane,
        const uint32& tile,
        const uint32& x_pos,
        const uint32& y_pos,
        const uint32& read_,
        const std::string& is_filtered,
        const uint32& control_number,
        const uint32& sample_number
    )
        : read(read_),
          before_read(),
          after_read() {

        // Before read:
        before_read = "@";
        before_read += instrument + ':' + std::to_string(run_number) + ':' + flowcell_ID + ':';
        before_read += std::to_string(lane) + ':' + std::to_string(tile) + ':';
        before_read += std::to_string(x_pos) + ':' + std::to_string(y_pos) + ' ';
        // After read
        after_read = ':' + is_filtered + ':';
        after_read += std::to_string(control_number) + ':';
        after_read += std::to_string(sample_number);

    };

    std::string get_line() {
        std::string out = before_read + std::to_string(read) + after_read;
        return out;
    }

private:

    std::string before_read;
    std::string after_read;

};




#endif
