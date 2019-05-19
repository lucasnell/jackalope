/*
 Functions to read/write to/from VCF files
*/


#include <RcppArmadillo.h>

#include <fstream>
#include <string>
#include <vector>
#include "zlib.h"
#ifdef _OPENMP
#include <omp.h>  // omp
#endif

#include <progress.hpp>  // for the progress bar


#include "jackalope_types.h"  // integer types
#include "seq_classes_ref.h"  // Ref* classes
#include "seq_classes_var.h"  // Var* classes
#include "str_manip.h"  // filter_nucleos
#include "util.h"  // str_stop, thread_check
#include "io.h"
#include "io_vcf.h"

using namespace Rcpp;





/*
 Determine whether this variant should be included in a VCF line for given
 sequence starting and ending positions.
 If this variant has a deletion at the input position, this method updates that
 and the boolean for whether the line is still expanding (changes it to true).
 */
void OneVarSeqVCF::check(const uint64& pos_start,
                         uint64& pos_end,
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
void OneVarSeqVCF::dump(std::vector<std::string>& unq_alts,
                        uint64& gt_tmp,
                        const uint64& pos_start,
                        const uint64& pos_end,
                        const std::string& ref_str) {

    if (gt_index > 0) {

        /*
         First create alternate string:
         */
        // Fill with reference sequence:
        std::string alt_str = ref_str;

        // Add mutations from back:
        const Mutation* mut;
        uint64 pos;
        uint64 n_muts = ind.second - ind.first + 1;
        for (uint64 i = 0; i < n_muts; i++) {
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










/*
 Set the strings for the sequence position (`POS`), reference sequence (`REF`),
 alternative alleles (`ALT`), and genotype information (`GT` format field)
 to add to a new line in the VCF file.
 */
bool WriterVCF::iterate(std::string& pos_str,
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
        for (uint64 i = 0; i < var_infos.size(); i++) {
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
    for (uint64 i = mut_pos.first; i <= mut_pos.second; i++) {
        ref_str.push_back(ref_nts->at(i));
    }

    /*
     Go back through and collect information for each variant that's
     getting included:
     */
    pos_str = std::to_string(mut_pos.first + 1);  //bc it's 1-based indexing
    unq_alts.clear();
    for (uint64 i = 0; i < var_infos.size(); i++) {
        var_infos[i].dump(unq_alts, gt_indexes[i], mut_pos.first, mut_pos.second,
                          ref_str);
    }
    /*
     `do_write` will be false if overlapping mutations result in the reference
     sequence again.
     It being false should be a very rare occurrence.
     */
    bool do_write = !unq_alts.empty();

    if (do_write) {

        // Fill alt. string:
        alt_str += unq_alts[0];
        for (uint64 i = 1; i < unq_alts.size(); i++) alt_str += ',' + unq_alts[i];


        /*
         Now fill genotype (`GT`) info, using `sample_groups` to group them
         */
        if (gt_strs.size() != sample_groups.n_rows) {
            str_stop({"\nInput vector for GT field info isn't the same size ",
                     "as the number of rows in the `sample_matrix` argument."});
        }
        uint64 gt_i;
        for (uint64 i = 0; i < sample_groups.n_rows; i++) {
            std::string& gt(gt_strs[i]);
            gt_i = gt_indexes[sample_groups(i,0)];
            gt = std::to_string(gt_i);
            for (uint64 j = 1; j < sample_groups.n_cols; j++) {
                gt_i = gt_indexes[sample_groups(i,j)];
                gt += '|';
                gt += std::to_string(gt_i);
            }
        }

    }

    /*
     Regardless of whether or not to write these mutations, we need to
     check for the new nearest mutation position.
     Otherwise, we'll be stuck in an infinite loop.
     */
    mut_pos = std::make_pair(MAX_INT, MAX_INT);
    for (uint64 i = 0; i < var_infos.size(); i++) {
        var_infos[i].compare_pos(mut_pos.first, mut_pos.second);
    }

    return do_write;
}







/*
 ==================================================================
                READ
 ==================================================================
 */

//' Read VCF from a vcfR object.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
SEXP read_vcfr(SEXP reference_ptr,
               const std::vector<std::string>& var_names,
               const std::vector<std::vector<std::string>>& haps_list,
               const std::vector<uint64>& seq_inds,
               const std::vector<uint64>& pos,
               const std::vector<std::string>& ref_seq) {

    XPtr<RefGenome> reference(reference_ptr);
    uint64 n_muts = haps_list.size();
    uint64 n_vars = var_names.size();
    uint64 n_seqs = reference->size();

    Mutation new_mut;

    XPtr<VarSet> var_set(new VarSet(*reference, var_names));

    for (uint64 mut_i = 0; mut_i < n_muts; mut_i++) {

        const std::string& ref(ref_seq[mut_i]);
        const std::vector<std::string>& haps(haps_list[mut_i]);
        const uint64& seq_i(seq_inds[mut_i]);

        for (uint64 var_i = 0; var_i < n_vars; var_i++) {

            const std::string& alt(haps[var_i]);

            // If it's blank or if it's the same as the reference, move on:
            if (alt.size() == 0 || alt == ref) continue;

            // Else, mutate accordingly:
            VarSequence& var_seq((*var_set)[var_i][seq_i]);

            if (alt.size() == ref.size()) {
                /*
                 ------------
                 substitution(s)
                 ------------
                 */
                for (uint64 i = 0; i < ref.size(); i++) {
                    if (alt[i] != ref[i]) {
                        new_mut = Mutation(pos[mut_i] + i, pos[mut_i] + i, alt[i]);
                        var_seq.mutations.push_back(new_mut);
                    }
                }
            } else if (alt.size() > ref.size()) {
                /*
                 ------------
                 insertion
                 ------------
                 */
                // Copy the string so it can be manipulated
                std::string alt_copy = alt;

                /*
                 For all sequences but the last in the REF string, just make
                 them substitutions if they differ from ALT.
                 */
                uint64 i = 0;
                for (; i < (ref.size()-1); i++) {
                    if (alt[i] != ref[i]) {
                        new_mut = Mutation(pos[mut_i] + i, pos[mut_i] + i, alt_copy[i]);
                        var_seq.mutations.push_back(new_mut);
                    }
                }
                // Erase all the nucleotides that have already been added (if any):
                if (ref.size() > 1) alt_copy.erase(0, ref.size() - 1U);
                /*
                 Make the last one an insertion proper
                 */
                new_mut = Mutation(pos[mut_i] + i, pos[mut_i] + i, alt_copy);
                var_seq.mutations.push_back(new_mut);

            } else {
                /*
                 ------------
                 deletion
                 ------------
                 */
                /*
                 For all sequences in the ALT string, just make them substitutions
                 if they differ from REF.
                 (Note that this goes to the end of ALT, not REF, as it does for
                 insertions.)
                 */
                uint64 i = 0;
                for (; i < alt.size(); i++) {
                    if (alt[i] != ref[i]) {
                        new_mut = Mutation(pos[mut_i] + i, pos[mut_i] + i, alt[i]);
                        var_seq.mutations.push_back(new_mut);
                    }
                }

                // size modifier:
                sint64 sm = static_cast<sint64>(alt.size()) -
                    static_cast<sint64>(ref.size());

                new_mut = Mutation(pos[mut_i] + i, pos[mut_i] + i, sm);
                var_seq.mutations.push_back(new_mut);

            }

        }

    }


    /*
     Go back and re-calculate positions and variant sequence sizes
     */
    for (uint64 seq_i = 0; seq_i < n_seqs; seq_i++) {
        for (uint64 var_i = 0; var_i < n_vars; var_i++) {
            VarSequence& var_seq((*var_set)[var_i][seq_i]);
            var_seq.calc_positions();
        }
    }

    return var_set;
}









/*
 ==================================================================
 WRITE
 ==================================================================
 */



//' Write `variants` to VCF file.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
void write_vcf_cpp(std::string out_prefix,
                   const int& compress,
                   SEXP var_set_ptr,
                   const IntegerMatrix& sample_matrix,
                   const bool& show_progress) {

    XPtr<VarSet> var_set(var_set_ptr);

    expand_path(out_prefix);

    if (any(sample_matrix < 1).is_true()) {
        str_stop({"\nIn the input matrix specifying which samples each ",
                 "variant belongs to, there are values < 1."});
    }
    if (any(sample_matrix > var_set->size()).is_true()) {
        str_stop({"\nIn the input matrix specifying which samples each ",
                 "variant belongs to, there are values > the number of variants."});
    }
    if (any(is_na(sample_matrix)).is_true()) {
        str_stop({"\nIn the input matrix specifying which samples each ",
                 "variant belongs to, there are missing values."});
    }


    // Start the `WriterVCF` object
    WriterVCF writer(*var_set, 0, sample_matrix);

    std::string file_name = out_prefix + ".vcf";

    if (compress > 0) {

        // Use wrapper of `BGZF` to write to compressed VCF file
        write_vcf_<FileBGZF>(var_set, file_name, compress, writer);

    } else {

        // Use wrapper of `std::ofstream` to write to uncompressed VCF file
        write_vcf_<FileUncomp>(var_set, file_name, compress, writer);

    }


    return;

}



