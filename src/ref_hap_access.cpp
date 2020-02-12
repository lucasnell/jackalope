
/*
 ********************************************************

 Retrieving and printing info from Ref* and Hap* classes, from R environment.
 Many of these functions are used for testing.

 ********************************************************
 */


#include "jackalope_config.h" // controls debugging and diagnostics output

#include <RcppArmadillo.h>
#include <vector>  // vector class
#include <string>  // string class
#include <algorithm>  // lower_bound, sort
#include <deque>  // deque
#include <cmath>  // pow, log, exp
#include <pcg/pcg_random.hpp> // pcg prng


#include "jackalope_types.h"  // integer types
#include "ref_classes.h"  // Ref* classes
#include "hap_classes.h"  // Hap* classes
#include "pcg.h"  // pcg seeding
#include "phylogenomics.h"  // match_ and template functions



using namespace Rcpp;



/*
 ========================================================================================
 ========================================================================================

 Printing

 ========================================================================================
 ========================================================================================
 */



//' Function to print info on a `RefGenome`.
//'
//' Access `RefGenome` class's print method from R.
//'
//' @noRd
//'
//[[Rcpp::export]]
void print_ref_genome(SEXP ref_genome_ptr) {
    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    ref_genome->print();
    return;
}



//' Function to print info on a HapSet.
//'
//' Access `HapSet` class's print method from R.
//'
//' @noRd
//'
//[[Rcpp::export]]
void print_hap_set(SEXP hap_set_ptr) {
    XPtr<HapSet> hap_set(hap_set_ptr);
    hap_set->print();
    return;
}




/*
 ========================================================================================
 ========================================================================================

 Making new

 ========================================================================================
 ========================================================================================
 */



//' Make a RefGenome object from a set of chromosomes.
//'
//' Used for testing.
//'
//' @noRd
//'
//[[Rcpp::export]]
SEXP make_ref_genome(const std::vector<std::string>& chroms) {

    // Make pointer:
    XPtr<RefGenome> ref_genome(new RefGenome(), true);

    // Reference RefGenome fields:
    uint64 n_chroms = chroms.size();
    std::deque<RefChrom>& chromosomes(ref_genome->chromosomes);
    uint64& total_size(ref_genome->total_size);

    // Add to fields:
    chromosomes = std::deque<RefChrom>(n_chroms, RefChrom());
    for (uint64 i = 0; i < n_chroms; i++) {
        chromosomes[i].nucleos = chroms[i];
        chromosomes[i].name = "chrom" + std::to_string(i);
        total_size += chroms[i].size();
    }

    return ref_genome;
}


//' Make a HapSet object from a RefGenome pointer and # haplotypes.
//'
//' Used for testing.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
SEXP make_hap_set(SEXP ref_genome_ptr, const uint64& n_haps) {
    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    // FYI, `n_haps` can be zero:
    XPtr<HapSet> hap_set(new HapSet(*ref_genome, n_haps), true);
    return hap_set;
}




/*
 ========================================================================================
 ========================================================================================

 Viewing numbers of chromosomes/haplotypes

 ========================================================================================
 ========================================================================================
 */


//[[Rcpp::export]]
IntegerVector view_ref_genome_nchroms(SEXP ref_genome_ptr) {
    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    IntegerVector out(1);
    out[0] = ref_genome->size();
    return out;
}


// Number of chromosomes
//[[Rcpp::export]]
IntegerVector view_hap_set_nchroms(SEXP hap_set_ptr) {
    XPtr<HapSet> hap_set(hap_set_ptr);
    IntegerVector out(1);
    out[0] = hap_set->reference->size();
    return out;
}

// Number of haplotypes
//[[Rcpp::export]]
IntegerVector view_hap_set_nhaps(SEXP hap_set_ptr) {
    XPtr<HapSet> hap_set(hap_set_ptr);
    IntegerVector out(1);
    out[0] = hap_set->size();
    return out;
}


/*
 ========================================================================================
 ========================================================================================

 Viewing chromosome sizes

 ========================================================================================
 ========================================================================================
 */


//[[Rcpp::export]]
IntegerVector view_ref_genome_chrom_sizes(SEXP ref_genome_ptr) {
    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    std::vector<uint64> tmp = ref_genome->chrom_sizes();
    IntegerVector out(tmp.size());
    for (uint64 i = 0; i < tmp.size(); i++) out[i] = tmp[i];
    return out;
}


//' See all chromosome sizes in a HapGenome object within a HapSet.
//'
//' @noRd
//'
//[[Rcpp::export]]
IntegerVector view_hap_genome_chrom_sizes(SEXP hap_set_ptr,
                                        const uint64& hap_ind) {

    XPtr<HapSet> hap_set(hap_set_ptr);
    const HapGenome& hap_genome((*hap_set)[hap_ind]);

    IntegerVector out(hap_genome.size());
    for (uint64 i = 0; i < hap_genome.size(); i++) {
        const HapChrom& hap_chrom(hap_genome.chromosomes[i]);
        out[i] = hap_chrom.chrom_size;
    }
    return out;
}







/*
 ========================================================================================
 ========================================================================================

 Viewing one chromosome from a genome

 ========================================================================================
 ========================================================================================
 */



//[[Rcpp::export]]
std::string view_ref_genome_chrom(SEXP ref_genome_ptr, const uint64& chrom_ind) {
    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    std::string out = (*ref_genome)[chrom_ind].nucleos;
    return out;
}


//' Function to piece together the strings for one chromosome in a HapGenome.
//'
//' @noRd
//'
//[[Rcpp::export]]
std::string view_hap_genome_chrom(SEXP hap_set_ptr,
                               const uint64& hap_ind,
                               const uint64& chrom_ind) {

    XPtr<HapSet> hap_set(hap_set_ptr);
    const HapChrom& hap_chrom((*hap_set)[hap_ind][chrom_ind]);
    std::string out = hap_chrom.get_chrom_full();
    return out;
}




/*
 ========================================================================================
 ========================================================================================

 Viewing an entire genome

 ========================================================================================
 ========================================================================================
 */


//[[Rcpp::export]]
std::vector<std::string> view_ref_genome(SEXP ref_genome_ptr) {
    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    std::vector<std::string> out(ref_genome->size(), "");
    for (uint64 i = 0; i < ref_genome->size(); i++) {
        const RefChrom& ref_chrom((*ref_genome)[i]);
        out[i] = ref_chrom.nucleos;
    }
    return out;
}

//' Function to piece together the strings for all chromosomes in a HapGenome.
//'
//' @noRd
//'
//[[Rcpp::export]]
std::vector<std::string> view_hap_genome(SEXP hap_set_ptr,
                                        const uint64& hap_ind) {

    XPtr<HapSet> hap_set(hap_set_ptr);
    const HapGenome& hap_genome((*hap_set)[hap_ind]);

    std::vector<std::string> out(hap_genome.size(), "");
    for (uint64 i = 0; i < hap_genome.size(); i++) {
        const HapChrom& hap_chrom(hap_genome[i]);
        out[i] = hap_chrom.get_chrom_full();
    }
    return out;
}





/*
 ========================================================================================
 ========================================================================================

 Viewing names

 ========================================================================================
 ========================================================================================
 */


//[[Rcpp::export]]
std::vector<std::string> view_ref_genome_chrom_names(SEXP ref_genome_ptr) {
    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    std::vector<std::string> out;
    out.reserve(ref_genome->size());
    for (const RefChrom& chrom : (*ref_genome).chromosomes) out.push_back(chrom.name);
    return out;
}

//' See all haplotype-genome names in a HapSet object.
//'
//' @noRd
//'
//[[Rcpp::export]]
std::vector<std::string> view_hap_set_hap_names(SEXP hap_set_ptr) {
    XPtr<HapSet> hap_set(hap_set_ptr);
    std::vector<std::string> out;
    out.reserve(hap_set->size());
    for (const HapGenome& vg : hap_set->haplotypes) out.push_back(vg.name);
    return out;
}


/*
 ========================================================================================
 ========================================================================================

 Viewing GC and other nucleotide content

 ========================================================================================
 ========================================================================================
 */


//' See GC content in a RefGenome object.
//'
//' @noRd
//'
//[[Rcpp::export]]
double view_ref_genome_gc_content(SEXP ref_genome_ptr,
                                  const uint64& chrom_ind,
                                  const uint64& start,
                                  const uint64& end) {

    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    const std::string& chrom = (*ref_genome)[chrom_ind].nucleos;
    double gc = gc_prop(chrom, start, end);
    return gc;
}

//' See GC content in a HapSet object.
//'
//' @noRd
//'
//[[Rcpp::export]]
double view_hap_set_gc_content(SEXP hap_set_ptr,
                               const uint64& chrom_ind,
                               const uint64& hap_ind,
                               const uint64& start,
                               const uint64& end) {
    XPtr<HapSet> hap_set(hap_set_ptr);
    const HapChrom& hap_chrom((*hap_set)[hap_ind][chrom_ind]);
    std::string chrom;
    uint64 mut_i = 0;
    hap_chrom.set_chrom_chunk(chrom, start, end - start + 1, mut_i);
    double gc = gc_prop(chrom);
    return gc;
}

//' See any nucleotide's content in a RefGenome object.
//'
//' @noRd
//'
//[[Rcpp::export]]
double view_ref_genome_nt_content(SEXP ref_genome_ptr,
                                  const char& nt,
                                  const uint64& chrom_ind,
                                  const uint64& start,
                                  const uint64& end) {

    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    const std::string& chrom = (*ref_genome)[chrom_ind].nucleos;
    double ntp = nt_prop(chrom, nt, start, end);
    return ntp;
}

//' See any nucleotide's content in a HapSet object.
//'
//' @noRd
//'
//[[Rcpp::export]]
double view_hap_set_nt_content(SEXP hap_set_ptr,
                               const char& nt,
                               const uint64& chrom_ind,
                               const uint64& hap_ind,
                               const uint64& start,
                               const uint64& end) {
    XPtr<HapSet> hap_set(hap_set_ptr);
    const HapChrom& hap_chrom((*hap_set)[hap_ind][chrom_ind]);
    std::string chrom;
    uint64 mut_i = 0;
    hap_chrom.set_chrom_chunk(chrom, start, end - start + 1, mut_i);
    double ntp = nt_prop(chrom, nt);
    return ntp;
}




/*
 ========================================================================================
 ========================================================================================

 Setting names

 ========================================================================================
 ========================================================================================
 */


//[[Rcpp::export]]
void set_ref_genome_chrom_names(
        SEXP ref_genome_ptr,
        const std::vector<uint64>& chrom_inds,
        const std::vector<std::string>& names) {
    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    if (names.size() != chrom_inds.size()) stop("names and chrom_inds aren't the same size");
    if (*std::max_element(chrom_inds.begin(), chrom_inds.end()) >= ref_genome->size()) {
        stop("at least one value in chrom_inds is too large");
    }
    for (uint64 i = 0; i < chrom_inds.size(); i++) {
        (*ref_genome)[chrom_inds[i]].name = names[i];
    }
    return;
}



//[[Rcpp::export]]
void clean_ref_genome_chrom_names(SEXP ref_genome_ptr) {

    XPtr<RefGenome> ref_genome(ref_genome_ptr);

    std::string char_map;
    char_map.reserve(256);
    for (uint64 i = 0; i < 256; i++) char_map += static_cast<char>(i);
    std::string bad_chars = " :;=%,\\|/\"\'";
    for (char& c : bad_chars) char_map[c] = '_';

    for (uint64 i = 0; i < ref_genome->size(); i++) {
        for (char& c : (*ref_genome)[i].name) {
            c = char_map[c];
        }
    }
    return;
}


//[[Rcpp::export]]
void set_hap_set_hap_names(
        SEXP hap_set_ptr,
        const std::vector<uint64>& hap_inds,
        const std::vector<std::string>& names) {
    XPtr<HapSet> hap_set(hap_set_ptr);
    if (names.size() != hap_inds.size()) stop("names and hap_inds aren't the same size");
    if (*std::max_element(hap_inds.begin(), hap_inds.end()) >= hap_set->size()) {
        stop("at least one value in hap_inds is too large");
    }
    for (uint64 i = 0; i < hap_inds.size(); i++) {
        (*hap_set)[hap_inds[i]].name = names[i];
    }
    return;
}



/*
 ========================================================================================
 ========================================================================================

 Remove elements

 ========================================================================================
 ========================================================================================
 */


//[[Rcpp::export]]
void remove_ref_genome_chroms(
        SEXP ref_genome_ptr,
        std::vector<uint64> chrom_inds) {

    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    std::deque<RefChrom>& chromosomes(ref_genome->chromosomes);

    // Checking for duplicates:
    std::sort(chrom_inds.begin(), chrom_inds.end());
    if (adjacent_find(chrom_inds.begin(), chrom_inds.end()) != chrom_inds.end()) {
        stop("duplicates detected in chrom_inds");
    }

    // Number of deleted nucleotides:
    uint64 n_del = 0;

    for (uint64 i = 1; i <= chrom_inds.size(); i++) {
        // Going backward so I don't have to update later ones each time:
        uint64 j = chrom_inds[(chrom_inds.size() - i)];
        n_del += chromosomes[j].size();
        chromosomes.erase(chromosomes.begin() + j);
    }
    clear_memory<std::deque<RefChrom>>(chromosomes);

    // Update total size:
    ref_genome->total_size -= n_del;

    return;
}



//[[Rcpp::export]]
void remove_hap_set_haps(
        SEXP hap_set_ptr,
        std::vector<uint64> hap_inds) {

    XPtr<HapSet> hap_set(hap_set_ptr);
    std::vector<HapGenome>& haplotypes(hap_set->haplotypes);

    // Checking for duplicates:
    std::sort(hap_inds.begin(), hap_inds.end());
    if (adjacent_find(hap_inds.begin(), hap_inds.end()) != hap_inds.end()) {
        stop("duplicates detected in hap_inds");
    }

    for (uint64 i = 1; i <= hap_inds.size(); i++) {
        // Going backward so I don't have to update later ones each time:
        uint64 j = hap_inds[(hap_inds.size() - i)];
        haplotypes.erase(haplotypes.begin() + j);
    }
    clear_memory<std::vector<HapGenome>>(haplotypes);
    return;
}





/*
 ========================================================================================
 ========================================================================================

 Add elements

 ========================================================================================
 ========================================================================================
 */


//[[Rcpp::export]]
void add_ref_genome_chroms(
        SEXP ref_genome_ptr,
        const std::vector<std::string>& new_chroms,
        const std::vector<std::string>& new_names) {

    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    std::deque<RefChrom>& chromosomes(ref_genome->chromosomes);

    if (new_chroms.size() != new_names.size()) {
        stop("In `add_ref_genome_chroms`, `new_chroms` must be the same size as `new_names`");
    }

    for (uint64 i = 0; i < new_chroms.size(); i++) {
        chromosomes.push_back(RefChrom(new_names[i], new_chroms[i]));
        // Update total size:
        ref_genome->total_size += new_chroms[i].size();
    }


    return;
}


// Add blank, named haplotypes
//[[Rcpp::export]]
void add_hap_set_haps(
        SEXP hap_set_ptr,
        const std::vector<std::string>& new_names) {

    XPtr<HapSet> hap_set(hap_set_ptr);
    std::vector<HapGenome>& haplotypes(hap_set->haplotypes);
    const RefGenome& ref(*(hap_set->reference));

    for (uint64 i = 0; i < new_names.size(); i++) {
        haplotypes.push_back(HapGenome(new_names[i], ref));
    }

    return;
}
// Duplicate existing haplotype(s):
//[[Rcpp::export]]
void dup_hap_set_haps(
        SEXP hap_set_ptr,
        const std::vector<uint64>& hap_inds,
        const std::vector<std::string>& new_names) {

    XPtr<HapSet> hap_set(hap_set_ptr);
    std::vector<HapGenome>& haplotypes(hap_set->haplotypes);
    const RefGenome& ref(*(hap_set->reference));

    if (hap_inds.size() != new_names.size()) {
        stop("In `dup_hap_set_haps`, `hap_inds` must be the same size as `new_names`");
    }
    if (*std::max_element(hap_inds.begin(), hap_inds.end()) >= haplotypes.size()) {
        stop("In `dup_hap_set_haps`, one or more `hap_inds` is too large");
    }

    for (uint64 i = 0; i < new_names.size(); i++) {
        // Add blank haplotype:
        haplotypes.push_back(HapGenome(new_names[i], ref));
        // Add mutation information:
        HapGenome& new_vg(haplotypes.back());
        const HapGenome& old_vg(haplotypes[hap_inds[i]]);
        for (uint64 j = 0; j < new_vg.chromosomes.size(); j++) {
            HapChrom& new_vs(new_vg.chromosomes[j]);
            const HapChrom& old_vs(old_vg.chromosomes[j]);
            new_vs.mutations = old_vs.mutations;
            new_vs.chrom_size = old_vs.chrom_size;
        }
    }

    return;
}




/*
 ========================================================================================
 ========================================================================================

 For molecular evolution testing

 ========================================================================================
 ========================================================================================
 */




//' Turns a HapGenome's mutations into a list of data frames.
//'
//' Internal function for testing.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
DataFrame view_mutations(SEXP hap_set_ptr, const uint64& hap_ind) {

    XPtr<HapSet> hap_set(hap_set_ptr);
    const HapGenome& hap_genome((*hap_set)[hap_ind]);

    uint64 n_muts = 0;
    for (const HapChrom& vs : hap_genome.chromosomes) n_muts += vs.mutations.size();

    std::vector<sint64> size_mod;
    size_mod.reserve(n_muts);
    std::vector<uint64> old_pos;
    old_pos.reserve(n_muts);
    std::vector<uint64> new_pos;
    new_pos.reserve(n_muts);
    std::vector<std::string> nucleos;
    nucleos.reserve(n_muts);
    std::vector<uint64> haps(n_muts, hap_ind);
    std::vector<uint64> chroms;
    chroms.reserve(n_muts);

    for (uint64 i = 0; i < hap_genome.size(); i++) {
        const HapChrom& hap_chrom(hap_genome.chromosomes[i]);
        uint64 n_muts_i = hap_chrom.mutations.size();
        for (uint64 j = 0; j < n_muts_i; ++j) {
            size_mod.push_back(hap_chrom.size_modifier(j));
            old_pos.push_back(hap_chrom.mutations.old_pos[j]);
            new_pos.push_back(hap_chrom.mutations.new_pos[j]);
            nucleos.push_back("");
            if (hap_chrom.mutations.nucleos[j] != nullptr) {
                nucleos.back() = std::string(hap_chrom.mutations.nucleos[j]);
            }
            chroms.push_back(i);
        }
    }

    DataFrame out = DataFrame::create(
        _["hap"] = haps,
        _["chrom"] = chroms,
        _["size_mod"] = size_mod,
        _["old_pos"] = old_pos,
        _["new_pos"] = new_pos,
        _["nucleos"] = nucleos);

    return out;
}


//' Turns a HapGenome's mutations into a list of data frames.
//'
//' Internal function for testing.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
List examine_mutations(SEXP hap_set_ptr, const uint64& hap_ind, const uint64& chrom_ind) {

    XPtr<HapSet> hap_set_xptr(hap_set_ptr);
    const HapGenome& hap_genome((*hap_set_xptr)[hap_ind]);
    const HapChrom& hap_chrom(hap_genome[chrom_ind]);
    const AllMutations& muts(hap_chrom.mutations);

    std::string bases = "TCAG";
    std::vector<uint64> base_inds(85);
    uint64 j = 0;
    for (const char& c : bases) {
        base_inds[static_cast<uint64>(c)] = j;
        j++;
    }

    uint64 n_muts = hap_chrom.mutations.size();
    arma::mat sub_mat(4, 4, arma::fill::zeros);
    uint64 max_ins = 0;
    uint64 max_del = 0;
    for (uint64 i = 0; i < n_muts; i++) {
        sint64 mi = hap_chrom.size_modifier(i);
        if (mi == 0) continue;
        if (mi > 0) {
            if (mi > static_cast<sint64>(max_ins)) max_ins = mi;
        } else {
            uint64 mid = static_cast<uint64>(std::abs(mi));
            if (mid > max_del) max_del = mid;
        }
    }
    arma::mat ins_mat(4, max_ins, arma::fill::zeros);
    arma::mat del_mat(4, max_del, arma::fill::zeros);
    std::vector<uint64> pos_vec(n_muts);

    for (uint64 mut_i = 0; mut_i < n_muts; mut_i++) {

        char c = (*(hap_chrom.ref_chrom))[muts.old_pos[mut_i]];
        uint64 i = base_inds[static_cast<uint64>(c)];
        sint64 smod = hap_chrom.size_modifier(mut_i);
        if (smod == 0) {
            uint64 j = base_inds[static_cast<uint64>(muts.nucleos[mut_i][0])];
            sub_mat(i, j)++;
        } else if (smod > 0) {
            uint64 j = static_cast<uint64>(smod - 1);
            ins_mat(i, j)++;
        } else {
            uint64 j = static_cast<uint64>(std::abs(smod + 1));
            del_mat(i, j)++;
        }

        pos_vec[mut_i] = hap_chrom.mutations.old_pos[mut_i];
    }

    List out = List::create(
        _["sub"] = wrap(sub_mat),
        _["ins"] = wrap(ins_mat),
        _["del"] = wrap(del_mat),
        _["pos"] = pos_vec);

    return out;
}




//' Add mutations manually from R.
//'
//' This section applies to the next 3 functions.
//'
//' Note that all indices are in 0-based C++ indexing. This means that the first
//' item is indexed by `0`, and so forth.
//'
//' @param hap_set_ptr External pointer to a C++ `HapSet` object
//' @param hap_ind Integer index to the desired haplotype. Uses 0-based indexing!
//' @param chrom_ind Integer index to the desired chromosome. Uses 0-based indexing!
//' @param new_pos_ Integer index to the desired subsitution location.
//'     Uses 0-based indexing!
//'
//' @noRd
NULL_ENTRY;

//' @describeIn add_mutations Add a substitution.
//'
//' @inheritParams add_mutations
//' @param nucleo_ Character to substitute for existing one.
//'
//' @noRd
//'
//[[Rcpp::export]]
void add_substitution(SEXP hap_set_ptr, const uint64& hap_ind,
                      const uint64& chrom_ind,
                      const char& nucleo_,
                      const uint64& new_pos_) {
    XPtr<HapSet> hap_set(hap_set_ptr);
    HapGenome& hap_genome((*hap_set)[hap_ind]);
    HapChrom& hap_chrom(hap_genome[chrom_ind]);
    hap_chrom.add_substitution(nucleo_, new_pos_);
    return;
}
//' @describeIn add_mutations Add an insertion.
//'
//' @inheritParams add_mutations
//' @param nucleos_ Nucleotides to insert at the desired location.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
void add_insertion(SEXP hap_set_ptr, const uint64& hap_ind,
                   const uint64& chrom_ind,
                   const std::string& nucleos_,
                   const uint64& new_pos_) {
    XPtr<HapSet> hap_set(hap_set_ptr);
    HapGenome& hap_genome((*hap_set)[hap_ind]);
    HapChrom& hap_chrom(hap_genome[chrom_ind]);
    hap_chrom.add_insertion(nucleos_, new_pos_);
    return;
}
//' @describeIn add_mutations Add a deletion.
//'
//' @inheritParams add_mutations
//' @param size_ Size of deletion.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
void add_deletion(SEXP hap_set_ptr,
                  const uint64& hap_ind,
                  const uint64& chrom_ind,
                  const uint64& size_,
                  const uint64& new_pos_) {
    XPtr<HapSet> hap_set(hap_set_ptr);
    HapGenome& hap_genome((*hap_set)[hap_ind]);
    HapChrom& hap_chrom(hap_genome[chrom_ind]);
    hap_chrom.add_deletion(size_, new_pos_);
    return;
}


