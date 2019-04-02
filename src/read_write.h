#ifndef __JACKAL_READ_WRITE_H
#define __JACKAL_READ_WRITE_H

#include <RcppArmadillo.h>
#include <vector>               // vector class
#include <string>               // string class
#include <iomanip>
#include <ctime>

#include <fstream>
#include <zlib.h>

// #ifdef _OPENMP
// #include <omp.h>  // omp
// #endif


#include "jackalope_types.h"  // integer types
#include "seq_classes_ref.h"  // Ref* classes
#include "seq_classes_var.h"  // Var* classes


//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(cpp11)]]



using namespace Rcpp;



// Size of the block of memory to use for reading non-indexed fasta files and fai files.
#define LENGTH 0x1000 // hexadecimel for 4096.
// Maximum uint32 value:
#define MAX_INT 4294967295UL




// For parsing segregating sites from ms-style output files
namespace parse_ms {
const std::string site = "segsites:";
const std::string pos = "positions:";
}


// Formatted date as "20190331"
inline std::string vcf_date() {
    // Obtain environment containing function
    Environment base("package:base");
    // Make function callable from C++
    Function sys_date = base["Sys.Date"];
    Function format = base["format"];
    // Call the function and change its output using format
    SEXP date = sys_date();
    SEXP fmt_date = format(date, "%Y%m%d");
    std::string vcf_date_str = as<std::string>(fmt_date);
    return vcf_date_str;
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

    /*
     Integer indicating whether to include this variant in the current line of VCF file,
     and, if so, which alternate string to use.
     0 indicates not to include it at all.
     After running `check` but before `dump`, it's set ambiguously to 1 if it should
     be included. The `dump` method assigns the actual index to the alternate string
     to use.
     */
    uint32 gt_index = 0;
    // Starting/ending indices on mutations vector:
    std::pair<uint32, uint32> ind = std::make_pair(0, 0);
    // Starting/ending position on ref. sequence:
    std::pair<uint32, uint32> pos;

    OneVarSeqVCF() {};


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

            gt_index = 1;
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

            /*
             Checking for a deletion right after the current mutation:
             (the second part of this statement is added because contiguous deletions
             are prevented)
             */
            if (ind.second < var_seq->mutations.size() &&
                var_seq->mutations[ind.second].size_modifier >= 0) {
                const Mutation& next_mut(var_seq->mutations[ind.second + 1]);
                if (next_mut.size_modifier < 0 &&
                    next_mut.old_pos == (var_seq->mutations[ind.second].old_pos + 1)) {
                    ind.second++;
                }
            }

            mut = &(var_seq->mutations[ind.second]);
            set_second_pos(*mut);

            if (pos.second > pos_end) {
                pos_end = pos.second;
                still_growing = true;
            }

        } else {

            gt_index = 0;

        }

        return;
    }


    /*
     This "dumps" the necessary haploid information for the VCF's `ALT` string,
     then iterates to the next mutation information
     */
    void dump(std::vector<std::string>& unq_alts,
              uint32& gt_tmp,
              const uint32& pos_start,
              const uint32& pos_end,
              const std::string& ref_str) {

        if (gt_index > 0) {

            /*
             First create alternate string:
             */
            // Fill with reference sequence:
            std::string alt_str = ref_str;

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

            /*
             Double-check that two mutations didn't combine to turn it back into the
             reference string:
             */
            if (alt_str != ref_str) {
                /*
                 Now see if that string exists already, and assign `gt_index` accordingly
                 */
                auto iter = std::find(unq_alts.begin(), unq_alts.end(), alt_str);
                // If it doesn't already exist, we add it:
                if (iter == unq_alts.end()) {
                    gt_index = unq_alts.size();
                    unq_alts.push_back(alt_str);
                } else {
                    gt_index = iter - unq_alts.begin();
                }
                gt_index++;  // <-- because alt. indices start at 1
                gt_tmp = gt_index;  // this stores `gt_index` before resetting below
            } else {
                gt_tmp = 0;
            }

            // Now iterate:
            ind.second++;
            ind.first = ind.second;
            reset_pos();
            gt_index = 0;

        } else {

            gt_tmp = 0;

        }

        return;
    }


    // Reset to new variant sequence
    void set_var(const VarSequence& var_seq_) {
        var_seq = &var_seq_;
        gt_index = 0;
        ind = std::make_pair(0, 0);
        reset_pos();
        return;
    }

    void compare_pos(uint32& pos_start,
                     uint32& pos_end) const {

        // If this is the new nearest mutation, override both positions:
        if (pos.first < pos_start) {
            pos_start = pos.first;
            pos_end = pos.second;
        }
        /*
         If this one ties with a previous mutation
         and the ending position of this one is further along than the original,
         then we override that.
         */
        if (pos.first == pos_start && pos.second > pos_end) {
            pos_end = pos.second;
        }

        return;

    }



private:

    const VarSequence* var_seq;

    // Set positions when indices are the same (when initializing, iterating, new variant)
    void reset_pos() {
        if (ind.first >= var_seq->mutations.size()) {
            pos = std::make_pair(MAX_INT, MAX_INT); // max uint32 values
        } else {
            const Mutation* mut(&(var_seq->mutations[ind.first]));
            set_first_pos(*mut);
            /*
             Checking for a deletion right after the current mutation:
             (the second part of this statement is added because contiguous deletions
             are prevented elsewhere)
             */
            if (ind.second < var_seq->mutations.size() && mut->size_modifier >= 0) {
                const Mutation& next_mut(var_seq->mutations[ind.second + 1]);
                if (next_mut.size_modifier < 0 &&
                    next_mut.old_pos == (mut->old_pos + 1)) {
                    ind.second++;
                    mut = &(var_seq->mutations[ind.second]);
                }
            }
            set_second_pos(*mut);
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
class WriterVCF {

public:

    const VarSet* var_set;
    uint32 seq_ind;
    const std::string* ref_nts;

    std::vector<OneVarSeqVCF> var_infos;
    // Starting/ending positions on reference sequence for overall nearest mutation:
    std::pair<uint32,uint32> mut_pos = std::make_pair(MAX_INT, MAX_INT);
    // Strings for all unique alt. strings among  variants. Grouping is not relevant here.
    std::vector<std::string> unq_alts;
    // Indices for how to group each variant's genotype information:
    arma::umat sample_groups;
    // Names for each sample:
    std::vector<std::string> sample_names;

    WriterVCF(const VarSet& var_set_,
              const uint32& seq_ind_,
              const IntegerMatrix& sample_groups_)
        : var_set(&var_set_),
          seq_ind(seq_ind_),
          ref_nts(),
          var_infos(var_set_.size()),
          unq_alts(),
          sample_groups(as<arma::umat>(sample_groups_) - 1),
          gt_indexes(var_set_.size()) {

        // Now checking sequence index:
        if (seq_ind >= var_set->reference->size()) {
            str_stop({"\nWhen specifying a sequence index for VCF output, ",
                     "you must provide an integer < the number of sequences."});
        }

        unq_alts.reserve(var_set_.size());
        construct();
        make_names();

    };


    /*
     Set the strings for the sequence position (`POS`), reference sequence (`REF`),
     alternative alleles (`ALT`), and genotype information (`GT` format field)
     to add to a new line in the VCF file.
     */
    void iterate(std::string& pos_str,
                 std::string& ref_str,
                 std::string& alt_str,
                 std::vector<std::string>& gt_strs) {

        // Reset all strings
        if (ref_str.size() > 0) ref_str.clear();
        if (alt_str.size() > 0) alt_str.clear();
        for (std::string& gt : gt_strs) if (gt.size() > 0) gt.clear();

        /*
         Boolean for whether we're still merging mutations.
         Only deletions can change this from false to true.
         */
        bool still_growing = true;
        /*
         Now going through sequences until it's no longer merging, updating the starting
         and ending positions each time:
         */
        while (still_growing) {
            still_growing = false;
            for (uint32 i = 0; i < var_infos.size(); i++) {
                var_infos[i].check(mut_pos.first, mut_pos.second, still_growing);
            }
        }

        // Create reference sequence:
        ref_str.reserve(mut_pos.second - mut_pos.first + 1);
        if (mut_pos.second >= ref_nts->size()) {
            stop(std::string("\nPosition ") + std::to_string(mut_pos.second) +
                std::string(" on ref. string is too high for total ") +
                std::string("ref. string length of ") +
                std::to_string(ref_nts->size()));
        }
        for (uint32 i = mut_pos.first; i <= mut_pos.second; i++) {
            ref_str.push_back(ref_nts->at(i));
        }

        /*
         Go back through and collect information for each variant that's
         getting included:
         */
        pos_str = std::to_string(mut_pos.first + 1);  //bc it's 1-based indexing
        unq_alts.clear();
        for (uint32 i = 0; i < var_infos.size(); i++) {
            var_infos[i].dump(unq_alts, gt_indexes[i], mut_pos.first, mut_pos.second,
                              ref_str);
        }
        if (unq_alts.empty()) stop("unq_alts.empty()");
        // Fill alt. string:
        alt_str += unq_alts[0];
        for (uint32 i = 1; i < unq_alts.size(); i++) alt_str += ',' + unq_alts[i];


        /*
         Now fill genotype (`GT`) info, using `sample_groups` to group them
         */
        if (gt_strs.size() != sample_groups.n_rows) {
            str_stop({"\nInput vector for GT field info isn't the same size ",
                     "as the number of rows in the `sample_matrix` argument."});
        }
        uint32 gt_i;
        for (uint32 i = 0; i < sample_groups.n_rows; i++) {
            std::string& gt(gt_strs[i]);
            gt_i = gt_indexes[sample_groups(i,0)];
            gt = std::to_string(gt_i);
            for (uint32 j = 1; j < sample_groups.n_cols; j++) {
                gt_i = gt_indexes[sample_groups(i,j)];
                gt += '|';
                gt += std::to_string(gt_i);
            }
        }


        // Check for the new nearest mutation position:
        mut_pos = std::make_pair(MAX_INT, MAX_INT);
        for (uint32 i = 0; i < var_infos.size(); i++) {
            var_infos[i].compare_pos(mut_pos.first, mut_pos.second);
        }

        return;
    }


    // Change the sequence this object refers to
    void new_seq(const uint& seq_ind_) {
        seq_ind = seq_ind_;
        construct();
        return;
    }

    // Fill a string with the header info
    void fill_header(std::string& chunk) {
        chunk = "##fileformat=VCFv4.3\n";
        chunk += "##fileDate=";
        chunk += vcf_date();
        chunk += '\n';
        chunk += "##source=jackalope\n";
        for (uint32 i = 0; i < var_set->reference->size(); i++) {
            const RefSequence& rs(var_set->reference->operator[](i));
            chunk += "##contig=<ID=" + rs.name + ',';
            chunk += "length=" + std::to_string(rs.size()) + ">\n";
        }
        chunk += "##phasing=full\n";
        chunk += "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number ";
        chunk +=    "of Samples With Data\">\n";
        chunk += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
        chunk += "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype";
        chunk +=    "Quality\">\n";
        chunk += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        for (uint32 i = 0; i < sample_names.size(); i++) {
            chunk += '\t' + sample_names[i];
        }
        chunk += '\n';
        return;
    }


private:

    std::vector<uint32> gt_indexes;  // temporarily stores gt info

    void construct() {

        ref_nts = &(var_set->reference->sequences[seq_ind].nucleos);

        /*
         Set pointer for the focal sequence in each variant
         and set positions in `mut_pos` field
         */
        for (uint32 i = 0; i < var_infos.size(); i++) {
            var_infos[i].set_var((*var_set)[i][seq_ind]);
            var_infos[i].compare_pos(mut_pos.first, mut_pos.second);
        }

        return;
    }

    // Creating vector of names for each sample:
    void make_names() {
        uint32 n_samples = sample_groups.n_rows;
        sample_names = std::vector<std::string>(n_samples, "");
        for (uint32 i = 0; i < n_samples; i++) {
            std::string& sn(sample_names[i]);
            sn = var_set->operator[](sample_groups(i,0)).name;
            for (uint32 j = 1; j < sample_groups.n_cols; j++) {
                sn += "__";
                sn += var_set->operator[](sample_groups(i,j)).name;
            }
        }
        return;
    }
};


// Overloaded for writing to stdout, uncompressed file, or gzipped file
inline void chunk_to_output(const std::string& null_str,
                            const std::string& chunk) {
    Rcout << chunk;
    return;
}
inline void chunk_to_output(std::ofstream& out_file,
                            const std::string& chunk) {
    out_file << chunk;
    return;
}
inline void chunk_to_output(gzFile& out_file,
                            const std::string& chunk) {
    gzwrite(out_file, chunk.c_str(), chunk.size());
    return;
}






//' Template doing most of the work for writing to a VCF file.
//'
//' `T` should be `std::string`, `std::ofstream`, or `gzFile`, for the three
//' specializations of the `chunk_to_output` function above.
//'
//' @noRd
//'
template <typename T>
void write_vcf_(XPtr<VarSet> var_set,
                T& out_file,
                WriterVCF writer) {

    // Very high quality that will essentially round to Pr(correct) = 1
    // (only needed as string):
    std::string max_qual = "441453";

    uint32 n_seqs = var_set->reference->size();
    uint32 n_samples = writer.sample_groups.n_rows;

    // String of text to append to, then to insert into output:
    std::string chunk;


    /*
     Header
     */
    writer.fill_header(chunk);
    chunk_to_output(out_file, chunk);

    /*
     Data lines
     */

    std::string pos_str = "";
    std::string ref_str = "";
    std::string alt_str = "";
    std::vector<std::string> gt_strs(n_samples, "");

    for (uint32 seq = 0; seq < n_seqs; seq++) {
        writer.new_seq(seq);
        while (writer.mut_pos.first < MAX_INT) {
            // Set information for this line:
            writer.iterate(pos_str, ref_str, alt_str, gt_strs);
            // CHROM
            chunk = var_set->reference->operator[](writer.seq_ind).name;
            // POS
            chunk += '\t' + pos_str;
            // ID
            chunk += "\t.";
            // REF
            chunk += '\t' + ref_str;
            // ALT
            chunk += '\t' + alt_str;
            // QUAL (setting to super high value)
            chunk += '\t' + max_qual;
            // FILTER
            chunk += "\tPASS";
            // INFO
            chunk += "\tNS=" + std::to_string(n_samples);
            // FORMAT
            chunk += "\tGT:GQ";
            // Sample info (setting GQ to super high value)
            for (uint32 i = 0; i < n_samples; i++) {
                chunk += '\t' + gt_strs[i];
                chunk += ':' + max_qual;
            }
            chunk += '\n';
            chunk_to_output(out_file, chunk);
        }
    }


    return;

}







#endif
