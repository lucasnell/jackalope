/*
 Functions to read/write to/from VCF files
*/

#include "jackalope_config.h" // controls debugging and diagnostics output

#include <RcppArmadillo.h>

#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include "zlib.h"
#ifdef _OPENMP
#include <omp.h>  // omp
#endif

#include <progress.hpp>  // for the progress bar


#include "jackalope_types.h"  // integer types
#include "ref_classes.h"  // Ref* classes
#include "hap_classes.h"  // Hap* classes
#include "str_manip.h"  // filter_nucleos, cpp_str_split_delim_str, count_substr
#include "util.h"  // str_stop, thread_check
#include "io.h"
#include "io_vcf.h"
#include "htslib/vcf.h"
#include "htslib/vcfutils.h"





#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1


using namespace Rcpp;





/*
 Determine whether this haplotype should be included in a VCF line for given
 chromosome starting and ending positions.
 If this haplotype has a deletion at the input position, this method updates that
 and the boolean for whether the line is still expanding (changes it to true).
 */
void OneHapChromVCF::check(const uint64& pos_start,
                         uint64& pos_end,
                         bool& still_growing) {

    if (pos_end >= ref_pos.first) {

        gt_index = 1;

        while (mut_ind.second < hap_chrom->mutations.size() &&
               get_first_pos(mut_ind.second) < pos_end) {

            mut_ind.second++;

        }

        if (mut_ind.second >= hap_chrom->mutations.size() ||
            (mut_ind.second < hap_chrom->mutations.size() &&
            get_first_pos(mut_ind.second) > pos_end)) {

            mut_ind.second--;

        }

        /*
         Checking for a deletion right after the current mutation:
         (the second part of this statement is added because contiguous deletions
         are prevented)
         */
        if (mut_ind.second < (hap_chrom->mutations.size() - 1) &&
            hap_chrom->size_modifier(mut_ind.second) >= 0) {
            if (hap_chrom->size_modifier(mut_ind.second + 1) < 0 &&
                hap_chrom->mutations.old_pos[mut_ind.second + 1] ==
                (hap_chrom->mutations.old_pos[mut_ind.second] + 1)) {
                mut_ind.second++;
            }
        }

        set_second_pos(mut_ind.second);

        if (ref_pos.second > pos_end) {
            pos_end = ref_pos.second;
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
void OneHapChromVCF::dump(std::vector<std::string>& unq_alts,
                          uint64& gt_tmp,
                          const uint64& pos_start,
                          const uint64& pos_end,
                          const std::string& ref_str) {

    const AllMutations& mutations(hap_chrom->mutations);

    if (gt_index > 0) {

        /*
         First create alternate string:
         */
        // Fill with reference chromosome:
        std::string alt_str = ref_str;

        // Add mutations from back:
        uint64 pos;
        uint64 n_muts = mut_ind.second - mut_ind.first + 1;
        for (uint64 i = 0; i < n_muts; i++) {
            uint64 index = mut_ind.second - i;
            pos = mutations.old_pos[index] - pos_start;
            if (pos >= alt_str.size()) {
                stop(std::string("\nPosition ") + std::to_string(pos) +
                    std::string(" on alt. string is too high for total ") +
                    std::string("alt. string length of ") +
                    std::to_string(alt_str.size()));
            }
            if (hap_chrom->size_modifier(index) == 0) { // substitution
                alt_str[pos] = mutations.nucleos[index][0];
            } else if (hap_chrom->size_modifier(index) > 0) { // insertion
                // Copy so we can remove last nucleotide before inserting:
                std::string nts(mutations.nucleos[index]);
                alt_str[pos] = nts.back();
                nts.pop_back();
                alt_str.insert(pos, nts);  // inserts before `pos`
            } else {  // deletion
                alt_str.erase(pos, static_cast<size_t>(
                        std::abs(hap_chrom->size_modifier(index))));
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
        mut_ind.second++;
        mut_ind.first = mut_ind.second;
        reset_pos();
        gt_index = 0;

    } else {

        gt_tmp = 0;

    }

    return;
}










/*
 Set the strings for the chromosome position (`POS`), reference chromosome (`REF`),
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
     Now going through chromosomes until it's no longer merging, updating the starting
     and ending positions each time:
     */
    while (still_growing) {
        still_growing = false;
        for (uint64 i = 0; i < hap_infos.size(); i++) {
            hap_infos[i].check(mut_pos.first, mut_pos.second, still_growing);
        }
    }


    // Create reference chromosome:
    if (mut_pos.second >= ref_nts->size()) {
        str_stop({"\nPosition ", std::to_string(mut_pos.second),
            " on ref. string is too high for total ",
            "ref. string length of ",
            std::to_string(ref_nts->size()), ". ",
            "For debugging, mut_pos.first = ", std::to_string(mut_pos.first)});
    }
    ref_str.reserve(mut_pos.second - mut_pos.first + 1);

    for (uint64 i = mut_pos.first; i <= mut_pos.second; i++) {
        ref_str.push_back(ref_nts->at(i));
    }

    /*
     Go back through and collect information for each haplotype that's
     getting included:
     */
    pos_str = std::to_string(mut_pos.first + 1);  //bc it's 1-based indexing
    unq_alts.clear();
    for (uint64 i = 0; i < hap_infos.size(); i++) {
        hap_infos[i].dump(unq_alts, gt_indexes[i], mut_pos.first, mut_pos.second,
                          ref_str);
    }


    /*
     `do_write` will be false if overlapping mutations result in the reference
     chromosome again.
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
    for (uint64 i = 0; i < hap_infos.size(); i++) {
        hap_infos[i].compare_pos(mut_pos.first, mut_pos.second);
    }

    return do_write;
}







/*
 ==================================================================
                READ
 ==================================================================
 */


/*
 Make haplotype names from vector of sample names and ploidy info.
 */
void make_hap_names(std::vector<std::string>& hap_names,
                    const std::vector<std::string>& samp_names,
                    const int& ploidy) {

    if (ploidy == 1) {

        hap_names = samp_names;

    } else {

        hap_names.reserve(samp_names.size() * ploidy);

        /*
         Check for whether they're output from jackalope.
         If so, then split by "__". Otherwise add "_1", "_2", etc.
         */
        bool from_jlp = count_substr(samp_names[0], "__") ==
            static_cast<uint32>(ploidy - 1);
        for (uint32 i = 1; i < samp_names.size(); i++) {
            const std::string& s(samp_names[i]);
            if (s.size() == 0) stop("Can't have zero-sized sample names in VCF files.");
            if (!from_jlp) continue;
            if (s.size() < 3) {
                from_jlp = false;
                break;
            }
            uint32 n_dunders = count_substr(s, "__");
            from_jlp = n_dunders == static_cast<uint32>(ploidy - 1);
            // Can't have "__" on either end:
            if (from_jlp) {
                from_jlp = !(s.front() == '_' && s[1] == '_') &&
                    !(s[s.size() - 2] == '_' && s.back() == '_');
            }
        }

        if (from_jlp) {

            for (const std::string& samp : samp_names) {
                std::vector<std::string> sub_samps = cpp_str_split_delim_str(samp, "__");
                for (const std::string& s : sub_samps) hap_names.push_back(s);
            }

        } else {

            for (const std::string& samp : samp_names) {
                for (int j = 0; j < ploidy; j++) {
                    hap_names.push_back(samp + '_' + std::to_string(j + 1));
                }
            }

        }

    }

    return;
}



/*

 Fill vectors of info from VCF file.

 Used info from http://wresch.github.io/2014/11/18/process-vcf-file-with-htslib.html
 and
 <https://github.com/samtools/htslib/blob/dd6f0b72c92591252bb77818663629cc1a129949/
 htslib/vcf.h#L835>

 */


int fill_vcf_info(const std::string& fn,
                  std::vector<std::string>& chrom_names,
                  std::vector<std::string>& hap_names,
                  std::vector<std::vector<std::string>>& alts_list,
                  std::vector<uint64>& chrom_inds,
                  std::vector<uint64>& positions,
                  std::vector<std::string>& ref_chrom) {

    int n_chroms = 0;

    // genotype data for each call
    int ngt_arr = 0;
    int ngt     = 0;
    int *gt     = NULL;

    // open VCF/BCF file
    htsFile * inf = bcf_open(fn.c_str(), "r");
    if (inf == NULL) {
        return EXIT_FAILURE;
    }

    // read header
    bcf_hdr_t *hdr = bcf_hdr_read(inf);
    int n_samps = bcf_hdr_nsamples(hdr);

    // Read sample names, to be used later for `hap_names`
    std::vector<std::string> samp_names;
    samp_names.reserve(n_samps);
    for (int k = 0; k < n_samps; k++) {
        samp_names.push_back(std::string(hdr->samples[k]));
    }

    // report names of all the chromosomes in the VCF file
    const char **c_names = NULL;
    c_names = bcf_hdr_seqnames(hdr, &n_chroms);
    if (c_names == NULL) {
        bcf_close(inf);
        bcf_hdr_destroy(hdr);
        return EXIT_FAILURE;
    }
    chrom_names.reserve(n_chroms);
    for (uint32 i = 0; i < static_cast<uint32>(n_chroms); i++) {
        chrom_names.push_back(std::string(c_names[i]));
    }

    // struc for storing each record
    bcf1_t *rec = bcf_init();
    if (rec == NULL) {
        free(c_names);
        bcf_close(inf);
        bcf_hdr_destroy(hdr);
        return EXIT_FAILURE;
    }

    uint64 chrom, pos;
    std::string ref;
    int ploidy = -1;

    while (bcf_read(inf, hdr, rec) == 0) {

        ngt = bcf_get_genotypes(hdr, rec, &gt, &ngt_arr);
        if (ngt <= 0) continue; // GT not present

        // This needs to come before making `ref` so that it fills `rec->d`
        bcf_unpack(rec, BCF_UN_ALL);

        chrom = static_cast<uint64>(rec->rid);
        pos = static_cast<uint64>(rec->pos);
        ref = std::string(rec->d.allele[0]);

        chrom_inds.push_back(chrom);
        positions.push_back(pos);
        ref_chrom.push_back(ref);

        if (ploidy != -1 && ploidy != static_cast<int>(ngt / n_samps)) {
            stop("All ploidy must be the same in VCF files.");
        }
        ploidy = ngt / n_samps;

        alts_list.push_back(std::vector<std::string>(0));
        std::vector<std::string>& alts(alts_list.back());
        alts.reserve(ngt);

        for (int i = 0; i < n_samps; i++) {
            int32_t *ptr = gt + i*ploidy;
            for (int j = 0; j < ploidy; j++) {
                // if true, the sample has smaller ploidy
                if (ptr[j] == bcf_int32_vector_end) {
                    stop("All samples must have the same ploidy");
                }

                // missing allele
                if (bcf_gt_is_missing(ptr[j])) {
                    alts.push_back("");
                    continue;
                }

                // the VCF 0-based allele index
                int allele_index = bcf_gt_allele(ptr[j]);

                alts.push_back(std::string(rec->d.allele[allele_index]));
            }
        }

    }


    // Memory management
    free(gt);
    free(c_names);
    bcf_hdr_destroy(hdr);
    bcf_close(inf);
    bcf_destroy(rec);

    // Create list of haplotype names:
    make_hap_names(hap_names, samp_names, ploidy);


    return EXIT_SUCCESS;
}




/*
 Add mutations to a HapSet object based on VCF-file info vectors.
 */

void add_vcf_mutations(HapSet& hap_set,
                       const std::vector<std::vector<std::string>>& alts_list,
                       const std::vector<uint64>& chrom_inds,
                       const std::vector<uint64>& positions,
                       const std::vector<std::string>& ref_chrom,
                       const std::vector<uint64>& ind_map) {

    uint64 n_muts = alts_list.size();
    uint64 n_haps = hap_set.size();

    arma::Mat<sint64> size_mods(n_haps, hap_set.reference->size(), arma::fill::zeros);

    sint64 size_mod_i; // used temporarily for each deletion and insertion

    uint64 new_pos;

    for (uint64 mut_i = 0; mut_i < n_muts; mut_i++) {

        const std::string& ref(ref_chrom[mut_i]);
        const std::vector<std::string>& alts(alts_list[mut_i]);
        const uint64& chrom_i(ind_map[chrom_inds[mut_i]]);

        for (uint64 hap_i = 0; hap_i < n_haps; hap_i++) {

            const std::string& alt(alts[hap_i]);

            // If it's blank or if it's the same as the reference, move on:
            if (alt.size() == 0 || alt == ref) continue;

            // Else, mutate accordingly:
            HapChrom& hap_chrom(hap_set[hap_i][chrom_i]);
            AllMutations& mutations(hap_chrom.mutations);
            sint64& size_mod(size_mods(hap_i, chrom_i));

            // Make sure that positions are never before any existing mutations
            if (!mutations.empty() && mutations.old_pos.back() >= positions[mut_i]) {
                str_stop({"\nFor VCF files, \"Positions are sorted numerically, in ",
                         "increasing order, within each reference sequence CHROM.\" ",
                         "(VCFv4.3 specification). ",
                         "In jackalope, multiple records with the same POS are also ",
                         "not permitted"});
            }

            if (alt.size() == ref.size()) {
                /*
                 ------------
                 substitution(s)
                 ------------
                 */
                for (uint64 i = 0; i < ref.size(); i++) {
                    if (alt[i] != ref[i]) {
                        new_pos = positions[mut_i] + i + size_mod;
                        mutations.push_back(positions[mut_i] + i, new_pos, alt[i]);
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
                 For all chromosomes but the last in the REF string, just make
                 them substitutions if they differ from ALT.
                 */
                uint64 i = 0;
                for (; i < (ref.size()-1); i++) {
                    if (alt[i] != ref[i]) {
                        new_pos = positions[mut_i] + i + size_mod;
                        mutations.push_back(positions[mut_i] + i, new_pos, alt_copy[i]);
                    }
                }
                // Erase all the nucleotides that have already been added (if any):
                if (ref.size() > 1) alt_copy.erase(0, ref.size() - 1U);
                /*
                 Make the last one an insertion proper
                 */
                size_mod_i = alt_copy.size() - 1;
                new_pos = positions[mut_i] + i + size_mod;
                mutations.push_back(positions[mut_i] + i, new_pos,
                                    alt_copy.c_str());
                size_mod += size_mod_i;
                hap_chrom.chrom_size += size_mod_i;

            } else {
                /*
                 ------------
                 deletion
                 ------------
                 */
                /*
                 For all chromosomes in the ALT string, just make them substitutions
                 if they differ from REF.
                 (Note that this goes to the end of ALT, not REF, as it does for
                 insertions.)
                 */
                uint64 i = 0;
                for (; i < alt.size(); i++) {
                    if (alt[i] != ref[i]) {
                        new_pos = positions[mut_i] + i + size_mod;
                        mutations.push_back(positions[mut_i] + i, new_pos, alt[i]);
                    }
                }

                size_mod_i = static_cast<sint64>(alt.size()) -
                    static_cast<sint64>(ref.size());

                new_pos = positions[mut_i] + i + size_mod;
                mutations.push_back(positions[mut_i] + i, new_pos, nullptr);
                size_mod += size_mod_i;
                hap_chrom.chrom_size += size_mod_i;

            }

        }

    }


    return;
}





//[[Rcpp::export]]
SEXP read_vcf_cpp(SEXP reference_ptr,
                  const std::string& fn,
                  const bool& print_names) {

    /*
     ------------
     Fill vectors of info from VCF file
     ------------
     */

    std::vector<std::string> chrom_names;
    std::vector<std::string> hap_names;
    std::vector<std::vector<std::string>> alts_list;
    std::vector<uint64> chrom_inds;
    std::vector<uint64> positions;
    std::vector<std::string> ref_chrom;

    /*
     Count # lines in file. This is an over-estimation bc it includes the header,
     but I'm okay with this.
     */
    uint64 n_lines;
    if (true) {  // nested scope to close inFile after it's done
        std::ifstream inFile(fn);
        n_lines = std::count(std::istreambuf_iterator<char>(inFile),
                             std::istreambuf_iterator<char>(), '\n');
    }

    // Reserve memory based on # lines:
    alts_list.reserve(n_lines);
    chrom_inds.reserve(n_lines);
    positions.reserve(n_lines);
    ref_chrom.reserve(n_lines);

    int status = fill_vcf_info(fn, chrom_names, hap_names, alts_list, chrom_inds,
                               positions, ref_chrom);

    if (status != EXIT_SUCCESS) {
        stop(std::string("Error reading file ") + fn);
    }


    /*
     ------------
     Now add VCF info to a new HapSet object
     ------------
     */
    XPtr<RefGenome> reference(reference_ptr);

    // Verify that names in the VCF file match those in the reference genome
    if (chrom_names.size() != reference->size()) {
        str_stop({"\nThe number of chromosomes in the VCF file doesn't match ",
                 "that for the `ref_genome` object."});
    }
    std::vector<std::string> ref_names;
    ref_names.reserve(reference->size());
    for (uint32 i = 0; i < reference->size(); i++) {
        ref_names.push_back(reference->chromosomes[i].name);
    }
    // Vector to map indices for position on `chrom_names` to position on
    // `reference->chromosomes`:
    std::vector<uint64> ind_map = match_chrom_names(ref_names, chrom_names, print_names);

    // Finally create HapSet
    XPtr<HapSet> hap_set(new HapSet(*reference, hap_names));
    // ...and add mutations:
    add_vcf_mutations(*hap_set, alts_list, chrom_inds, positions, ref_chrom, ind_map);

    return hap_set;

}













/*
 ==================================================================
 WRITE
 ==================================================================
 */



//' Write `haplotypes` to VCF file.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
void write_vcf_cpp(std::string out_prefix,
                   const int& compress,
                   SEXP hap_set_ptr,
                   const IntegerMatrix& sample_matrix,
                   const bool& show_progress) {

    XPtr<HapSet> hap_set(hap_set_ptr);

    expand_path(out_prefix);

    if (any(sample_matrix < 1).is_true()) {
        str_stop({"\nIn the input matrix specifying which samples each ",
                 "haplotype belongs to, there are values < 1."});
    }
    if (any(sample_matrix > hap_set->size()).is_true()) {
        str_stop({"\nIn the input matrix specifying which samples each ",
                 "haplotype belongs to, there are values > the number of haplotypes."});
    }
    if (any(is_na(sample_matrix)).is_true()) {
        str_stop({"\nIn the input matrix specifying which samples each ",
                 "haplotype belongs to, there are missing values."});
    }


    // Start the `WriterVCF` object
    WriterVCF writer(*hap_set, 0, sample_matrix);

    std::string file_name = out_prefix + ".vcf";

    if (compress > 0) {

        // Use wrapper of `BGZF` to write to compressed VCF file
        write_vcf_<FileBGZF>(hap_set, file_name, compress, writer);

    } else {

        // Use wrapper of `std::ofstream` to write to uncompressed VCF file
        write_vcf_<FileUncomp>(hap_set, file_name, compress, writer);

    }


    return;

}



