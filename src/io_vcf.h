#ifndef __JACKAL_VCF_IO_H
#define __JACKAL_VCF_IO_H

#include <RcppArmadillo.h>
#include <vector>               // vector class
#include <string>               // string class

#include <fstream>
#include "zlib.h"

#include "htslib/bgzf.h"  // BGZF

#include "jackalope_types.h"  // integer types
#include "seq_classes_ref.h"  // Ref* classes
#include "seq_classes_var.h"  // Var* classes
#include "io.h"  // File* classes



using namespace Rcpp;


// Maximum uint64 value:
#define MAX_INT 18446744073709551615ULL



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
    uint64 gt_index = 0;
    // Starting/ending indices on mutations vector:
    std::pair<uint64, uint64> ind = std::make_pair(0, 0);
    // Starting/ending position on ref. sequence:
    std::pair<uint64, uint64> pos;

    OneVarSeqVCF() : var_seq(nullptr) {};


    /*
     Determine whether this variant should be included in a VCF line for given
     sequence starting and ending positions.
     If this variant has a deletion at the input position, this method updates that
     and the boolean for whether the line is still expanding (changes it to true).
     */
    void check(const uint64& pos_start,
               uint64& pos_end,
               bool& still_growing);


    /*
     This "dumps" the necessary haploid information for the VCF's `ALT` string,
     then iterates to the next mutation information
     */
    void dump(std::vector<std::string>& unq_alts,
              uint64& gt_tmp,
              const uint64& pos_start,
              const uint64& pos_end,
              const std::string& ref_str);


    // Reset to new variant sequence
    void set_var(const VarSequence& var_seq_) {
        var_seq = &var_seq_;
        gt_index = 0;
        ind = std::make_pair(0, 0);
        reset_pos();
        return;
    }

    void compare_pos(uint64& pos_start,
                     uint64& pos_end) const {

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
            pos = std::make_pair(MAX_INT, MAX_INT); // max uint64 values
        } else {
            const Mutation* mut(&(var_seq->mutations[ind.first]));
            set_first_pos(*mut);
            /*
             Checking for a deletion right after the current mutation:
             (the second part of this statement is added because contiguous deletions
             are prevented elsewhere)
             */
            if (ind.second < (var_seq->mutations.size() - 1) &&
                mut->size_modifier >= 0) {
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
    inline uint64 get_first_pos(const Mutation& mut) {
        uint64 pos_first = mut.old_pos;
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
    inline uint64 get_second_pos(const Mutation& mut) {
        uint64 pos_second = mut.old_pos;
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
    uint64 seq_ind;
    const std::string* ref_nts;

    std::vector<OneVarSeqVCF> var_infos;
    // Starting/ending positions on reference sequence for overall nearest mutation:
    std::pair<uint64,uint64> mut_pos = std::make_pair(MAX_INT, MAX_INT);
    // Strings for all unique alt. strings among  variants. Grouping is not relevant here.
    std::vector<std::string> unq_alts;
    // Indices for how to group each variant's genotype information:
    arma::umat sample_groups;
    // Names for each sample:
    std::vector<std::string> sample_names;

    WriterVCF(const VarSet& var_set_,
              const uint64& seq_ind_,
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
     Returns false if you shouldn't write to file for this iteration (if all mutations
     by chance have been cancelled out result in the reference sequence).
     Returns true otherwise.
     */
    bool iterate(std::string& pos_str,
                 std::string& ref_str,
                 std::string& alt_str,
                 std::vector<std::string>& gt_strs);


    // Change the sequence this object refers to
    void new_seq(const uint64& seq_ind_) {
        seq_ind = seq_ind_;
        construct();
        return;
    }

    // Fill a string with the header info
    void fill_header(std::string& pool) {
        pool = "##fileformat=VCFv4.3\n";
        pool += "##fileDate=";
        pool += vcf_date();
        pool += '\n';
        pool += "##source=jackalope\n";
        for (uint64 i = 0; i < var_set->reference->size(); i++) {
            const RefSequence& rs(var_set->reference->operator[](i));
            pool += "##contig=<ID=" + rs.name + ',';
            pool += "length=" + std::to_string(rs.size()) + ">\n";
        }
        pool += "##phasing=full\n";
        pool += "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number ";
        pool +=    "of Samples With Data\">\n";
        pool += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
        pool += "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype";
        pool +=    "Quality\">\n";
        pool += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        for (uint64 i = 0; i < sample_names.size(); i++) {
            pool += '\t' + sample_names[i];
        }
        pool += '\n';
        return;
    }


private:

    std::vector<uint64> gt_indexes;  // temporarily stores gt info

    void construct() {

        ref_nts = &(var_set->reference->sequences[seq_ind].nucleos);

        /*
         Set pointer for the focal sequence in each variant
         and set positions in `mut_pos` field
         */
        for (uint64 i = 0; i < var_infos.size(); i++) {
            var_infos[i].set_var((*var_set)[i][seq_ind]);
            var_infos[i].compare_pos(mut_pos.first, mut_pos.second);
        }

        return;
    }

    // Creating vector of names for each sample:
    void make_names() {
        uint64 n_samples = sample_groups.n_rows;
        sample_names = std::vector<std::string>(n_samples, "");
        for (uint64 i = 0; i < n_samples; i++) {
            std::string& sn(sample_names[i]);
            sn = var_set->operator[](sample_groups(i,0)).name;
            for (uint64 j = 1; j < sample_groups.n_cols; j++) {
                sn += "__";
                sn += var_set->operator[](sample_groups(i,j)).name;
            }
        }
        return;
    }
};





//' Template doing most of the work for writing to a VCF file.
//'
//' `T` should be `FileUncomp` or `FileBGZF` from `io.h`
//'
//' @noRd
//'
template <typename T>
inline void write_vcf_(XPtr<VarSet> var_set,
                       const std::string& file_name,
                       const int& compress,
                       WriterVCF writer) {

    T out_file(file_name, compress);

    // Very high quality that will essentially round to Pr(correct) = 1
    // (only needed as string):
    std::string max_qual = "441453";

    uint64 n_seqs = var_set->reference->size();
    uint64 n_samples = writer.sample_groups.n_rows;

    // String of text to append to, then to insert into output:
    std::string pool;


    /*
     Header
     */
    writer.fill_header(pool);
    out_file.write(pool);

    /*
     Data lines
     */

    std::string pos_str = "";
    std::string ref_str = "";
    std::string alt_str = "";
    std::vector<std::string> gt_strs(n_samples, "");

    for (uint64 seq = 0; seq < n_seqs; seq++) {
        writer.new_seq(seq);
        while (writer.mut_pos.first < MAX_INT) {
            Rcpp::checkUserInterrupt();
            /*
             Set information for this line, unless by chance multiple mutations
             cause it to revert back to the reference. This would result in
             `writer.iterate` to return false. It should occur very rarely.
             */
            if (writer.iterate(pos_str, ref_str, alt_str, gt_strs)) {
                // CHROM
                pool = var_set->reference->operator[](writer.seq_ind).name;
                // POS
                pool += '\t' + pos_str;
                // ID
                pool += "\t.";
                // REF
                pool += '\t' + ref_str;
                // ALT
                pool += '\t' + alt_str;
                // QUAL (setting to super high value)
                pool += '\t' + max_qual;
                // FILTER
                pool += "\tPASS";
                // INFO
                pool += "\tNS=" + std::to_string(n_samples);
                // FORMAT
                pool += "\tGT:GQ";
                // Sample info (setting GQ to super high value)
                for (uint64 i = 0; i < n_samples; i++) {
                    pool += '\t' + gt_strs[i];
                    pool += ':' + max_qual;
                }
                pool += '\n';
                out_file.write(pool);
            }
        }
    }

    out_file.close();

    return;

}







#endif
