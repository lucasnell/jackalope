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
               const std::vector<uint32>& seq_inds,
               const std::vector<uint32>& pos,
               const std::vector<std::string>& ref_seq) {

    XPtr<RefGenome> reference(reference_ptr);
    uint32 n_muts = haps_list.size();
    uint32 n_vars = var_names.size();
    uint32 n_seqs = reference->size();

    Mutation new_mut;

    XPtr<VarSet> var_set(new VarSet(*reference, var_names));

    for (uint32 mut_i = 0; mut_i < n_muts; mut_i++) {

        const std::string& ref(ref_seq[mut_i]);
        const std::vector<std::string>& haps(haps_list[mut_i]);
        const uint32& seq_i(seq_inds[mut_i]);

        for (uint32 var_i = 0; var_i < n_vars; var_i++) {

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
                for (uint32 i = 0; i < ref.size(); i++) {
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
                uint32 i = 0;
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
                uint32 i = 0;
                for (; i < alt.size(); i++) {
                    if (alt[i] != ref[i]) {
                        new_mut = Mutation(pos[mut_i] + i, pos[mut_i] + i, alt[i]);
                        var_seq.mutations.push_back(new_mut);
                    }
                }

                // size modifier:
                sint32 sm = static_cast<sint32>(alt.size()) -
                    static_cast<sint32>(ref.size());

                new_mut = Mutation(pos[mut_i] + i, pos[mut_i] + i, sm);
                var_seq.mutations.push_back(new_mut);

            }

        }

    }


    /*
     Go back and re-calculate positions and variant sequence sizes
     */
    for (uint32 seq_i = 0; seq_i < n_seqs; seq_i++) {
        for (uint32 var_i = 0; var_i < n_vars; var_i++) {
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



