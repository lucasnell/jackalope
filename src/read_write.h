#ifndef __JACKAL_READ_WRITE_H
#define __JACKAL_READ_WRITE_H

#include <RcppArmadillo.h>
#include <vector>               // vector class
#include <string>               // string class
#include <ctime>

#include <fstream>
#include <zlib.h>

// #ifdef _OPENMP
// #include <omp.h>  // omp
// #endif


#include "jackal_types.h"  // integer types
#include "seq_classes_ref.h"  // Ref* classes
#include "seq_classes_var.h"  // Var* classes


//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(cpp11)]]



using namespace Rcpp;



// Size of the block of memory to use for reading non-indexed fasta files and fai files.
#define LENGTH 0x1000 // hexadecimel for 4096.
// Maximum uint32 value:
#define MAX_INT 4294967295




// For parsing segregating sites from ms-style output files
namespace parse_ms {
    const std::string site = "segsites:";
    const std::string pos = "positions:";
}



inline void expand_path(std::string& file_name) {
    // Obtain environment containing function
    Environment base("package:base");
    // Make function callable from C++
    Function pe_r = base["path.expand"];
    // Call the function and receive its list output
    SEXP fai_file_exp = pe_r(file_name);
    file_name = as<std::string>(fai_file_exp);
    return;
}




// Info for one sequence on one variant
class OneVarSeqVCF {

public:

    // Whether to include this variant in the current line of VCF file:
    bool include = false;
    // Starting/ending indices on mutations vector:
    std::pair<uint32, uint32> ind = std::make_pair(0, 0);
    // Starting/ending position on ref. sequence:
    std::pair<uint32, uint32> pos;

    OneVarSeqVCF(const VarSequence& var_seq_)
        : var_seq(&var_seq_) {
        reset_pos();
    };


    /*
    Determine whether this variant should be included in a VCF line for given
    sequence starting and ending positions.
    If this variant has a deletion at the input position, this method updates that
    and the boolean for whether the line is still expanding (changes it to true).
    */
    void check(const uint32& pos_start,
               uint32& pos_end,
               bool& still_growing) {

        if (pos_end >= pos.first) {

            include = true;
            const Mutation* mut(&(var_seq->mutations[ind.second]));

            while (ind.second < var_seq->mutations.size() &&
                   get_first_pos(var_seq->mutations[ind.second]) < pos_end) {

                ind.second++;

            }

            if (ind.second >= var_seq->mutations.size() ||
                (ind.second < var_seq->mutations.size() &&
                get_first_pos(var_seq->mutations[ind.second]) > pos_end)) {

                ind.second--;

            }

            mut = &(var_seq->mutations[ind.second]);
            set_second_pos(*mut);

            if (pos.second > pos_end) {
                pos_end = pos.second;
                still_growing = true;
            }

        }

        return;
    }


    /*
    This "dumps" the necessary haploid information for the VCF's `ALT` string,
    then iterates to the next mutation information
    */
    void dump(std::string& alt_str,
              const uint32& pos_start,
              const uint32& pos_end) {

        if (include) {

            // Fill with reference sequence:
            alt_str.reserve((pos_end - pos_start + 1) * 2); // `* 2` in case of insertions
            for (uint32 i = pos_start; i <= pos_end; i++) {
                alt_str.push_back(var_seq->ref_seq->nucleos[i]);
            }

            // Add mutations from back:
            const Mutation* mut;
            uint32 pos;
            uint32 n_muts = ind.second - ind.first + 1;
            for (uint32 i = 0; i < n_muts; i++) {
                mut = &(var_seq->mutations[ind.second - i]);
                pos = mut->old_pos - pos_start;
                if (pos >= alt_str.size()) {
                    stop(std::string("\nPosition ") + std::to_string(pos) +
                        std::string(" on alt. string is too high for total ") +
                        std::string("alt. string length of ") +
                        std::to_string(alt_str.size()));
                }
                if (mut->size_modifier == 0) { // substitution
                    alt_str[pos] = mut->nucleos[0];
                } else if (mut->size_modifier > 0) { // insertion
                    // Copy so we can remove last nucleotide before inserting:
                    std::string nts = mut->nucleos;
                    alt_str[pos] = nts.back();
                    nts.pop_back();
                    alt_str.insert(pos, nts);  // inserts before `pos`
                } else {  // deletion
                    alt_str.erase(pos, static_cast<size_t>(std::abs(mut->size_modifier)));
                }
            }

            // Now iterate:
            ind.second++;
            ind.first = ind.second;
            reset_pos();
            include = false;

        } else alt_str.resize(0);

        return;
    }


    // Reset to new variant sequence
    void new_var(const VarSequence& var_seq_) {
        var_seq = &var_seq_;
        include = false;
        ind = std::make_pair(0, 0);
        reset_pos();
        return;
    }

private:

    const VarSequence* var_seq;

    // Set positions when indices are the same (when initializing, iterating, new variant)
    void reset_pos() {
        if (ind.first >= var_seq->mutations.size()) {
            pos = std::make_pair(MAX_INT, MAX_INT); // max uint32 values
        } else {
            const Mutation& mut(var_seq->mutations[ind.first]);
            set_first_pos(mut);
            set_second_pos(mut);
        }
        return;
    }

    /*
    Gets first ref. sequence position for a mutation, compensating for the fact
    that deletions have to be treated differently
    */
    inline void set_first_pos(const Mutation& mut) {
        pos.first = mut.old_pos;
        if (mut.size_modifier < 0 && mut.old_pos > 0) pos.first--;
        return;
    }
    // Same as above, but returns the integer rather than setting it
    inline uint32 get_first_pos(const Mutation& mut) {
        uint32 pos_first = mut.old_pos;
        if (mut.size_modifier < 0 && mut.old_pos > 0) pos_first--;
        return pos_first;
    }
    /*
    Gets last ref. sequence position for a mutation, compensating for the fact
    that deletions have to be treated differently
    */
    inline void set_second_pos(const Mutation& mut) {
        pos.second = mut.old_pos;
        if (mut.size_modifier < 0) {
            if (mut.old_pos > 0) {
                pos.second -= (1 + mut.size_modifier);
            } else {
                pos.second -= mut.size_modifier;
            }
        }
        return;
    }
    // Same as above, but returns the integer rather than setting it
    inline uint32 get_second_pos(const Mutation& mut) {
        uint32 pos_second = mut.old_pos;
        if (mut.size_modifier < 0) {
            if (mut.old_pos > 0) {
                pos_second -= (1 + mut.size_modifier);
            } else {
                pos_second -= mut.size_modifier;
            }
        }
        return pos_second;
    }

};

// Map mutations among all variants for one sequence
class OneSeqMutMapper {

public:

    const VarSet* var_set;
    uint32 seq_ind;
    std::vector<uint32> near_mut_ind; // index for nearest mutation for each variant
    /* starting and ending positions for nearest mutation for each variant: */
    std::vector<uint32> near_mut_start;
    std::vector<uint32> near_mut_end;
    std::vector<sint32> near_mut_smod; // size mod. for nearest mutation for each variant
    uint32 nearest_mut;  // position for overall nearest mutation

    OneSeqMutMapper(const VarSet& var_set_,
                    const uint32& seq_ind_)
        : var_set(&var_set_),
          seq_ind(seq_ind_),
          near_mut_ind(var_set_.size(), 0),
          near_mut_start(var_set_.size()),
          near_mut_end(var_set_.size()),
          near_mut_smod(var_set_.size(), 0),
          nearest_mut() {

        construct();

        // So they don't move around in memory
        to_add_inds.reserve(var_set_.size());
        to_add_smods.reserve(var_set_.size());

    };


    /*
    Set the reference sequence and info for mutations to add to a new line in the
    VCF file.
    It returns a boolean for whether it has more line(s) to add.
    */
    bool iterate(std::string& ref_str,
                 std::string& alt_str) {

        // ref_str;
        // Last end position for any that overlapping mutations:
        uint32 last_end = nearest_mut;
        // index for last mutation to include for each variant
        std::vector<uint32> last_mut_ind = near_mut_ind;
        // boolean for whether we're still merging mutations
        bool still_merging = true;
        // Vector of booleans for whether to include each variant:
        std::vector<bool> include(var_set->size(), false);

        while (still_merging) {
            still_merging = false;
            for (uint32 i = 0; i < var_set->size(); i++) {
                // const VarSequence& vs((*var_set)[i][seq_ind]);
                if (near_mut_start[i] <= last_end) {
                    if (near_mut_end[i] > last_end) {
                        last_end = near_mut_end[i];
                        still_merging = true;
                    }
                    ;
                    // if (near_mut_smod[i] < 0) {
                    //     if (near_mut_start[i] > 0) near_mut_start[i]--;
                    //     near_mut_end[i] = near_mut_start[i] - near_mut_smod[i];
                    // }
                }
            }
        }

        // for (uint32 i = 0; i < var_set->size(); i++) {
        //     const VarSequence& vs((*var_set)[i][seq_ind]);
        //     if (near_mut_pos[i] == nearest_mut) {
        //         const Mutation& mut(vs.mutations[near_mut_ind[i]]);
        //         near_mut_pos[i] = mut.old_pos;
        //         near_mut_smod[i] = mut.size_modifier;
        //     } else near_mut_pos[i] = MAX_INT; // max uint32 value
        // }


        // // Now go back through and calculate for next iteration
        // for (uint32 i = 0; i < var_set->size(); i++) {
        //     const VarSequence& vs((*var_set)[i][seq_ind]);
        //     if (near_mut_pos[i] == nearest_mut) {
        //         const Mutation& mut(vs.mutations[near_mut_ind[i]]);
        //         near_mut_pos[i] = mut.old_pos;
        //         // Deletions are treated differently, unless they are at pos 0:
        //         if (mut.size_modifier < 0 && mut.old_pos > 0) {
        //             near_mut_pos[i]--;
        //         }
        //     } else near_mut_pos[i] = MAX_INT; // max uint32 value
        // }
        // nearest_mut = *std::min_element(near_mut_pos.begin(), near_mut_pos.end());
        if (nearest_mut == MAX_INT) return false;
        return true;
    }


    // Change the sequence this object refers to
    void new_seq(const uint& seq_ind_) {
        seq_ind = seq_ind_;
        near_mut_ind = std::vector<uint32>(var_set->size(), 0);
        construct();
        return;
    }


private:

    std::vector<uint32> to_add_inds;
    std::vector<sint32> to_add_smods;  // size modifiers for mutations to add

    void construct() {

        for (uint32 i = 0; i < var_set->size(); i++) {
            const VarSequence& vs((*var_set)[i][seq_ind]);
            if (!vs.mutations.empty()) {
                const Mutation& mut(vs.mutations.front());
                near_mut_start[i] = mut.old_pos;
                near_mut_smod[i] = mut.size_modifier;
                // Deletions have to be treated differently
                if (near_mut_smod[i] < 0) {
                    if (near_mut_start[i] > 0) near_mut_start[i]--;
                    near_mut_end[i] = near_mut_start[i] - near_mut_smod[i];
                } else near_mut_end[i] = near_mut_start[i];
            } else {
                near_mut_start[i] = MAX_INT; // max uint32 value
                near_mut_end[i] = MAX_INT;
            }
        }

        nearest_mut = *std::min_element(near_mut_start.begin(), near_mut_start.end());

        return;
    }
};



//' Write `variants` to VCF file, uncompressed.
//'
//' Starting out doing for only one sequence.
//'
//' @noRd
//'
//[[Rcpp::export]]
void write_vcf_uncomp(std::string file_name,
                      SEXP var_set_ptr) {

    XPtr<VarSet> var_set(var_set_ptr);

    expand_path(file_name);

    // std::ofstream out_file(file_name);
    //
    // if (!out_file.is_open()) {
    //     Rcout << "Unable to open file " << file_name << std::endl;
    // }
    uint32 n_vars = var_set->size();
    uint32 n_seqs = var_set->reference->size();

    Rcout << "##fileformat=VCFv4.3\n";
    Rcout << "##fileDate="; // << date_str << '\n';
    std::time_t t = std::time(nullptr);
    std::tm tm = *std::localtime(&t);
    Rcout << std::put_time(&tm, "%Y%m%d") << '\n';
    Rcout << "##source=jackal\n";
    for (uint32 i = 0; i < n_seqs; i++) {
        const RefSequence& rs(var_set->reference->operator[](i));
        Rcout << "##contig=<ID=" << rs.name << ',';
        Rcout << "length=" << rs.size() << ">\n";
    }
    Rcout << "##phasing=full\n";
    Rcout << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    Rcout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (uint32 i = 0; i < n_vars; i++) {
        Rcout << '\t' << var_set->operator[](i).name;
    }
    Rcout << '\n';

    // for (uint32 i = 0; i < ref.size(); i++) {
    //     out_file << '>';
    //     out_file << ref[i].name;
    //     out_file << '\n';
    //
    //     const std::string& seq_str(ref[i].nucleos);
    //
    //     for (uint32 pos = 0; pos < seq_str.size(); pos++) {
    //         out_file << seq_str[pos];
    //         if ((pos % text_width) == (text_width - 1)) out_file << '\n';
    //     }
    //     out_file << '\n';
    // }
    //
    // out_file.close();

    return;

}




#endif
