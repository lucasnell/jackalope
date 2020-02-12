#ifndef __JACKALOPE_VCF_IO_H
#define __JACKALOPE_VCF_IO_H

#include "jackalope_config.h" // controls debugging and diagnostics output

#include <RcppArmadillo.h>
#include <vector>               // vector class
#include <string>               // string class

#include <fstream>
#include "zlib.h"

#include "htslib/bgzf.h"  // BGZF

#include "jackalope_types.h"  // integer types
#include "ref_classes.h"  // Ref* classes
#include "hap_classes.h"  // Hap* classes
#include "io.h"  // File* classes



using namespace Rcpp;


// Maximum uint64 value:
#define MAX_INT 18446744073709551615ULL



/*
 This function produces a vector of indices that map the VCF indices onto the
 original chromosome names from the `ref_genome` object.
 This is in case the chromosomes are in a different order in the VCF file.
 */
inline std::vector<uint64> match_chrom_names(const std::vector<std::string>& from_ref,
                                             const std::vector<std::string>& from_vcf,
                                             const bool& print_names) {

    std::vector<uint64> order_(from_ref.size());

    for (uint64 i = 0; i < order_.size(); i++) {
        auto iter = std::find(from_vcf.begin(), from_vcf.end(),
                              from_ref[i]);
        if (iter == from_vcf.end()) {
            std::vector<std::string> err_msg;
            if (print_names) {
                for (const std::string& s : from_vcf) err_msg.push_back(s + '\n');
            }
            err_msg.push_back("\nChromosome name(s) in VCF file don't match those in ");
            err_msg.push_back("the `ref_genome` object. It's probably easiest ");
            err_msg.push_back("to manually change the `ref_genome` object ");
            err_msg.push_back("(using `$set_names()` method) to have the same names ");
            err_msg.push_back("as the VCF file.");
            str_stop(err_msg);
        }
        order_[i] = iter - from_vcf.begin();
    }

    return order_;
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





// Info for one chromosome on one haplotype
class OneHapChromVCF {

public:

    /*
     Integer indicating whether to include this haplotype in the current line of VCF file,
     and, if so, which alternate string to use.
     0 indicates not to include it at all.
     After running `check` but before `dump`, it's set ambiguously to 1 if it should
     be included. The `dump` method assigns the actual index to the alternate string
     to use.
     */
    uint64 gt_index = 0;
    // Starting/ending indices on mutations vector:
    std::pair<uint64, uint64> mut_ind = std::make_pair(0, 0);
    // Starting/ending position on ref. chromosome:
    std::pair<uint64, uint64> ref_pos;

    OneHapChromVCF() : hap_chrom(nullptr) {};


    /*
     Determine whether this haplotype should be included in a VCF line for given
     chromosome starting and ending positions.
     If this haplotype has a deletion at the input position, this method updates that
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


    // Reset to new haplotype chromosome
    void set_hap(const HapChrom& hap_chrom_) {
        hap_chrom = &hap_chrom_;
        gt_index = 0;
        mut_ind = std::make_pair(0, 0);
        reset_pos();
        return;
    }

    void compare_pos(uint64& pos_start,
                     uint64& pos_end) const {

        // If this is the new nearest mutation, override both positions:
        if (ref_pos.first < pos_start) {
            pos_start = ref_pos.first;
            pos_end = ref_pos.second;
        }
        /*
         If this one ties with a previous mutation
         and the ending position of this one is further along than the original,
         then we override that.
         */
        if (ref_pos.first == pos_start && ref_pos.second > pos_end) {
            pos_end = ref_pos.second;
        }

        return;

    }




private:

    const HapChrom* hap_chrom;

    // Set positions when indices are the same (when initializing, iterating, new haplotype)
    void reset_pos() {

        if (mut_ind.first >= hap_chrom->mutations.size()) {
            ref_pos = std::make_pair(MAX_INT, MAX_INT); // max uint64 values
        } else {
            uint64 index = mut_ind.first;
            set_first_pos(mut_ind.first);
            /*
             Checking for a deletion right after the current mutation:
             (the second part of this statement is added because contiguous deletions
             are prevented elsewhere)
             */
            if (mut_ind.second < (hap_chrom->mutations.size()-1) &&
                hap_chrom->size_modifier(mut_ind.first) >= 0) {
                if (hap_chrom->size_modifier(mut_ind.second + 1) < 0 &&
                    hap_chrom->mutations.old_pos[mut_ind.second + 1] ==
                    (hap_chrom->mutations.old_pos[mut_ind.first] + 1)) {
                    mut_ind.second++;
                    index = mut_ind.second;
                }
            }
            set_second_pos(index);
        }
        return;
    }

    /*
     Gets first ref. chromosome position for a mutation, compensating for the fact
     that deletions have to be treated differently
     */
    inline void set_first_pos(const uint64& index) {
        ref_pos.first = hap_chrom->mutations.old_pos[index];
        if (hap_chrom->size_modifier(index) < 0 &&
            hap_chrom->mutations.old_pos[index] > 0) ref_pos.first--;
        return;
    }
    // Same as above, but returns the integer rather than setting it
    inline uint64 get_first_pos(const uint64& index) {
        uint64 pos_first = hap_chrom->mutations.old_pos[index];
        if (hap_chrom->size_modifier(index) < 0 &&
            hap_chrom->mutations.old_pos[index] > 0) pos_first--;
        return pos_first;
    }
    /*
     Gets last ref. chromosome position for a mutation, compensating for the fact
     that deletions have to be treated differently
     */
    inline void set_second_pos(const uint64& index) {
        ref_pos.second = hap_chrom->mutations.old_pos[index];
        if (hap_chrom->size_modifier(index) < 0) {
            if (hap_chrom->mutations.old_pos[index] > 0) {
                ref_pos.second -= (1 + hap_chrom->size_modifier(index));
            } else {
                ref_pos.second -= hap_chrom->size_modifier(index);
            }
        }
        return;
    }
    // Same as above, but returns the integer rather than setting it
    inline uint64 get_second_pos(const uint64& index) {
        uint64 pos_second = hap_chrom->mutations.old_pos[index];
        if (hap_chrom->size_modifier(index) < 0) {
            if (hap_chrom->mutations.old_pos[index] > 0) {
                pos_second -= (1 + hap_chrom->size_modifier(index));
            } else {
                pos_second -= hap_chrom->size_modifier(index);
            }
        }
        return pos_second;
    }

};




// Map mutations among all haplotypes for one chromosome
class WriterVCF {

public:

    const HapSet* hap_set;
    uint64 chrom_ind;
    const std::string* ref_nts;

    std::vector<OneHapChromVCF> hap_infos;
    // Starting/ending positions on reference chromosome for overall nearest mutation:
    std::pair<uint64,uint64> mut_pos = std::make_pair(MAX_INT, MAX_INT);
    // Strings for all unique alt. strings among  haplotypes. Grouping is not relevant here.
    std::vector<std::string> unq_alts;
    // Indices for how to group each haplotype's genotype information:
    arma::umat sample_groups;
    // Names for each sample:
    std::vector<std::string> sample_names;

    WriterVCF(const HapSet& hap_set_,
              const uint64& chrom_ind_,
              const IntegerMatrix& sample_groups_)
        : hap_set(&hap_set_),
          chrom_ind(chrom_ind_),
          ref_nts(),
          hap_infos(hap_set_.size()),
          unq_alts(),
          sample_groups(as<arma::umat>(sample_groups_) - 1),
          gt_indexes(hap_set_.size()) {

        // Now checking chromosome index:
        if (chrom_ind >= hap_set->reference->size()) {
            str_stop({"\nWhen specifying a chromosome index for VCF output, ",
                     "you must provide an integer < the number of chromosomes."});
        }

        unq_alts.reserve(hap_set_.size());
        construct();
        make_names();

    };


    /*
     Set the strings for the chromosome position (`POS`), reference chromosome (`REF`),
     alternative alleles (`ALT`), and genotype information (`GT` format field)
     to add to a new line in the VCF file.
     Returns false if you shouldn't write to file for this iteration (if all mutations
     by chance have been cancelled out result in the reference chromosome).
     Returns true otherwise.
     */
    bool iterate(std::string& pos_str,
                 std::string& ref_str,
                 std::string& alt_str,
                 std::vector<std::string>& gt_strs);


    // Change the chromosome this object refers to
    void new_chrom(const uint64& chrom_ind_) {
        chrom_ind = chrom_ind_;
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
        for (uint64 i = 0; i < hap_set->reference->size(); i++) {
            const RefChrom& rs(hap_set->reference->operator[](i));
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

        ref_nts = &(hap_set->reference->chromosomes[chrom_ind].nucleos);

        /*
         Set pointer for the focal chromosome in each haplotype
         and set positions in `mut_pos` field
         */
        for (uint64 i = 0; i < hap_infos.size(); i++) {
            hap_infos[i].set_hap((*hap_set)[i][chrom_ind]);
            hap_infos[i].compare_pos(mut_pos.first, mut_pos.second);
        }

        return;
    }

    // Creating vector of names for each sample:
    void make_names() {
        uint64 n_samples = sample_groups.n_rows;
        sample_names = std::vector<std::string>(n_samples, "");
        for (uint64 i = 0; i < n_samples; i++) {
            std::string& sn(sample_names[i]);
            sn = hap_set->operator[](sample_groups(i,0)).name;
            for (uint64 j = 1; j < sample_groups.n_cols; j++) {
                sn += "__";
                sn += hap_set->operator[](sample_groups(i,j)).name;
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
inline void write_vcf_(XPtr<HapSet> hap_set,
                       const std::string& file_name,
                       const int& compress,
                       WriterVCF writer) {

    T out_file(file_name, compress);

    // Very high quality that will essentially round to Pr(correct) = 1
    // (only needed as string):
    std::string max_qual = "441453";

    uint64 n_chroms = hap_set->reference->size();
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

    for (uint64 chrom = 0; chrom < n_chroms; chrom++) {
        writer.new_chrom(chrom);
        while (writer.mut_pos.first < MAX_INT) {
            Rcpp::checkUserInterrupt();
            /*
             Set information for this line, unless by chance multiple mutations
             cause it to revert back to the reference. This would result in
             `writer.iterate` to return false. It should occur very rarely.
             */
            if (writer.iterate(pos_str, ref_str, alt_str, gt_strs)) {
                // CHROM
                pool = hap_set->reference->operator[](writer.chrom_ind).name;
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
