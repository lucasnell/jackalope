#include <RcppArmadillo.h>
#include <string> // string class
#include <algorithm> // sort
#include <vector> // vector class

#include "gemino_types.h" // integer types, VariantSet, SequenceSet
#include "str_manip.h" // cpp_to_upper
#include "util.h" // gc_content

using namespace Rcpp;






// ======================================================================================
// ======================================================================================

//      SequenceSet: Info for sets of DNA sequences

// ======================================================================================
// ======================================================================================



//'
//' Create a SequenceSet from vectors for sequences and sequence names
//'
//' @param sequences Character vector of sequences.
//' @param seq_names Character vector of sequence names.
//'
//' @return A pointer to a \code{SequenceSet} class in C++ (class \code{Rcpp::XPtr}).
//'
//'
// [[Rcpp::export]]
XPtr<SequenceSet> SequenceSet_characters(std::vector<std::string>& sequences,
                                         std::vector<std::string>& seq_names,
                                         bool merged = true) {

    XPtr<SequenceSet> ss(new SequenceSet, true);
    ss->merged = merged;
    ss->total_size = 0;
    ss->seq_sizes = std::vector<uint>(0);
    ss->seq_names = seq_names;
    ss->sequences = sequences;

    uint size_i;

    for (int i = 0; i < seq_names.size(); i++) {
        size_i = sequences[i].size();
        ss->total_size += size_i;
        ss->seq_sizes.push_back(size_i);
    }

    return ss;
}


// [[Rcpp::export]]
uint64 total_size_SequenceSet(XPtr<SequenceSet> ss) {
    return ss->total_size;
}

// [[Rcpp::export]]
std::vector<uint> seq_sizes_SequenceSet(XPtr<SequenceSet> ss) {
    return ss->seq_sizes;
}

// [[Rcpp::export]]
char one_nucleo_SequenceSet(XPtr<SequenceSet> ss, int scaff, int pos) {
    return ss->sequences[scaff][pos];
}

// [[Rcpp::export]]
std::string one_scaff_SequenceSet(XPtr<SequenceSet> ss, int scaff) {
    return ss->sequences[scaff];
}

// [[Rcpp::export]]
bool is_merged_SequenceSet(XPtr<SequenceSet> ss) {
    return ss->merged;
}

// [[Rcpp::export]]
std::vector<std::string> seq_names_SequenceSet(XPtr<SequenceSet> ss) {
    return ss->seq_names;
}

// [[Rcpp::export]]
std::vector<double> all_scaff_gc_SequenceSet(XPtr<SequenceSet> ss) {
    std::vector<double> gc_prop_vector;
    for (std::string& s : ss->sequences) {
        double gc = gc_prop(s);
        gc_prop_vector.push_back(gc);
    }
    return gc_prop_vector;
}

// [[Rcpp::export]]
double one_scaff_gc_SequenceSet(XPtr<SequenceSet> ss, const std::string& seq_name) {
    auto iter = std::find(ss->seq_names.begin(), ss->seq_names.end(), seq_name);
    if (iter == ss->seq_names.end()) return NA_REAL;
    uint ind = iter - ss->seq_names.begin();
    std::string& s = ss->sequences[ind];
    double gc = gc_prop(s);
    return gc;
}

// [[Rcpp::export]]
double one_range_gc_SequenceSet(XPtr<SequenceSet> ss, const std::string& seq_name,
                                const uint& start, const uint& stop) {
    auto iter = std::find(ss->seq_names.begin(), ss->seq_names.end(), seq_name);
    if (iter == ss->seq_names.end()) return NA_REAL;
    uint ind = iter - ss->seq_names.begin();
    std::string& s = ss->sequences[ind];
    double gc = gc_prop(s, start, stop);
    return gc;
}

// [[Rcpp::export]]
std::vector<double> mult_ranges_gc_SequenceSet(XPtr<SequenceSet> ss,
                                          const std::vector<std::string>& seq_names,
                                          const std::vector<uint>& starts,
                                          const std::vector<uint>& stops) {
    uint n_ranges = seq_names.size();
    std::vector<double> gc_prop_vector(n_ranges);
    if (starts.size() != n_ranges || stops.size() != n_ranges) {
        stop("seq_names, starts, and stops must all be of same length in "
                 "mult_ranges_gc_SequenceSet");
    }
    for (uint i = 0; i < n_ranges; i++) {
        double gc = one_range_gc_SequenceSet(ss, seq_names[i], starts[i], stops[i]);
        gc_prop_vector[i] = gc;
    }
    return gc_prop_vector;
}







//'
//' Create a new SequenceSet from an existing one's pointer.
//'
//' @param input Input sequence set object.
//'
//' @return A pointer to a \code{SequenceSet} class in C++ (class \code{Rcpp::XPtr}).
//'
//'
// [[Rcpp::export]]
XPtr<SequenceSet> SequenceSet_copy(XPtr<SequenceSet> input) {

    XPtr<SequenceSet> ss(new SequenceSet, true);
    ss->merged = input->merged;
    ss->total_size = input->total_size;
    ss->seq_sizes = input->seq_sizes;
    ss->seq_names = input->seq_names;
    ss->sequences = input->sequences;

    return ss;
}




// [[Rcpp::export]]
void SummarizeSequenceSet(XPtr<SequenceSet> ss, int console_width) {

    // 32 characters is the narrowest I'll allow
    // (I'd start getting negative lengths and other problems otherwise)
    console_width = std::max(console_width, 32);

    int num_seqs = ss->sequences.size();
    std::vector<int> inds;
    if (num_seqs <= 10) {
        for (int i = 0; i < num_seqs; i++) inds.push_back(i);
    } else {
        for (int i = 0; i < 5; i++) inds.push_back(i);
        inds.push_back(-1);
        for (int i = (num_seqs - 5 + 1); i < num_seqs; i++) inds.push_back(i);
    }

    Rcout.imbue(std::locale(""));
    Rcout << "< Set of " << num_seqs << " sequences >" << std::endl;
    Rcout << "# Total size: " << ss->total_size << " bp" << std::endl;
    // Rcout << "# Sequences:" << std::endl;

    int ind_i, name_width = 10, length_width = 9;
    // Console width minus name width AND length width AND spaces between
    int seq_print_len = console_width - name_width - length_width - 2;
    // Number of chars print before and after elipses for a long string
    int before_elips = std::ceil((seq_print_len - 3) / 2);
    int after_elips = seq_print_len - 3 - before_elips;

    Rprintf("%-*s %s%-*s %*s\n", name_width, "  name",
            std::string(before_elips - 4, ' ').c_str(),
            seq_print_len - (before_elips - 4), "sequence",
            length_width, "length");

    for (int i = 0; i < inds.size(); ++i) {
        ind_i = inds[i];
        if (ind_i == -1) {
            Rprintf("%-10s %-*s %9s\n", "...", seq_print_len, "...", "...");
            continue;
        }
        // Print name
        const std::string& name_i = ss->seq_names[ind_i];
        Rprintf("%-10.10s ", name_i.c_str());
        // Print sequence
        const std::string& seq_i = ss->sequences[ind_i];
        if (seq_i.size() > seq_print_len){
            for (int j = 0; j < before_elips; j++) Rcout << seq_i[j];
            Rcout << "...";
            for (int j = (seq_i.size() - after_elips); j < seq_i.size(); j++) {
                Rcout << seq_i[j];
            }
        } else {
            Rprintf("%-*s", seq_print_len, seq_i.c_str());
        }
        // Print width
        if (ss->seq_sizes[ind_i] > 999999999) {
            Rprintf(" %9.2E", ss->seq_sizes[ind_i]);
        } else {
            Rprintf(" %9i", ss->seq_sizes[ind_i]);
        }
        Rcout << std::endl;
    }
}



















// ======================================================================================
// ======================================================================================

//      VariantSet: Info for sets of genomic variants

// ======================================================================================
// ======================================================================================


//'
//' Create a VariantSet from input vectors.
//'
//' Note: The estimate of total segregating sites will be underestimated because
//' it counts all insertions as one site, even if they're > 1 bp long.
//'
//' @param nucleos A list of character vectors, each string representing all the
//'     nucleotides present at segregating sites for a particular variant on a
//'     particular scaffold.
//' @param sites A list of integer vectors, each vector representing all the
//'     segregating sites among all variants for one scaffold.
//' @param scaffold_lengths A list of numeric or integer vectors, each vector
//'     representing all the scaffold lengths for a particular variant.
//'
//'
//' @return A pointer to a \code{VariantSet} class in C++ (class \code{Rcpp::XPtr}).
//'
//'
// [[Rcpp::export]]
XPtr<VariantSet> VariantSet_vectors(
        const std::vector< std::vector<std::string> >& nucleos,
        const std::vector< std::vector< std::vector<uint> > >& sites,
        const std::vector< std::vector<uint> >& scaffold_lengths) {

    uint n_vars = sites.size();
    uint n_scaffs = sites[0].size();
    for (uint i = 0; i < n_vars; i++) {
        if (sites[i].size() != n_scaffs) {
            stop("All variants should have the same number of scaffolds");
        }
    }

    XPtr<VariantSet> vs(new VariantSet(n_vars, n_scaffs), true);

    for (uint i = 0; i < n_scaffs; i++) {
        std::vector<uint> sites_all;
        for (uint j = 0; j < n_vars; j++) {
            vs->variant_info[j].sites[i] = sites[j][i];
            for (uint k = 0; k < nucleos[j][i].size(); k++) {
                vs->variant_info[j].nucleos[i].push_back(nucleos[j][i][k]);
            }
            vs->variant_info[j].scaffold_lengths[i] = scaffold_lengths[j][i];
            sites_all.insert(sites_all.end(), sites[j][i].begin(),
                                 sites[j][i].end());
            std::sort(sites_all.begin(), sites_all.end());
            sites_all.erase(unique(sites_all.begin(), sites_all.end()),
                            sites_all.end());
        }
        vs->total_segr_sites += sites_all.size();
    }

    return vs;
}

// ----------
// Accessing VariantSet fields
// ----------

// [[Rcpp::export]]
int n_variants_VS(const XPtr<VariantSet>& vs) {
    return vs->n_variants;
}
// [[Rcpp::export]]
uint total_segr_sites_VS(const XPtr<VariantSet>& vs) {
    return vs->total_segr_sites;
}

// One variant's nucleos for one scaffold
// [[Rcpp::export]]
std::vector<char> nucleos_VS(const XPtr<VariantSet>& vs,
                             uint variant_index,
                             uint scaff_index) {

    const std::vector<char>& vni = vs->variant_info[variant_index].nucleos[scaff_index];
    std::vector<char> out = vni;
    // std::string out = "";
    // for (uint i = 0; i < vni.size(); i++) {
    //     if (vni[i] != '\0') out += vni[i];
    // }

    return out;
}

// One variant's sites for one scaffold
// [[Rcpp::export]]
std::vector<uint> sites_VS(const XPtr<VariantSet>& vs,
                           uint variant_index,
                           uint scaff_index) {

    std::vector<uint> out = vs->variant_info[variant_index].sites[scaff_index];

    return out;
}

// One variant's scaffold length for one scaffold
// [[Rcpp::export]]
uint scaffold_length_VS(XPtr<VariantSet> vs, uint variant_index, uint scaff_index) {

    uint out = vs->variant_info[variant_index].scaffold_lengths[scaff_index];

    return out;
}



// ----------

//'
//' Create a new VariantSet from an existing one's pointer.
//'
//' @param input Input variant set object.
//'
//' @return A pointer to a \code{VariantSet} class in C++ (class \code{Rcpp::XPtr}).
//'
//'
// [[Rcpp::export]]
XPtr<VariantSet> VariantSet_copy(XPtr<VariantSet> input) {

    XPtr<VariantSet> vs(new VariantSet, true);
    vs->n_variants = input->n_variants;
    vs->total_segr_sites = input->total_segr_sites;
    vs->variant_info = input->variant_info;

    return vs;
}




// ----------
// Accessing variant scaffolds from a VariantSet object
// ----------


/*
 Retrieve one scaffold from one variant
 */

// [[Rcpp::export]]
std::string variants_retr_scaff(uint scaff_num,
                                uint variant_num,
                                const XPtr<SequenceSet>& ss,
                                const XPtr<VariantSet>& vs) {

    // Going from R to C++ indexing
    uint variant_index = variant_num - 1;
    uint scaff_index = scaff_num - 1;

    const std::vector<char>& nucleos = vs->variant_info[variant_index].nucleos[scaff_index];
    const std::vector<uint>& sites = vs->variant_info[variant_index].sites[scaff_index];
    const std::string& ref = ss->sequences[scaff_index];

    std::string scaffold_out = "";
    uint sites_index = 0;
    for (uint i = 0; i < ref.size(); i++) {
        if (sites[sites_index] == i) {
            while (sites[sites_index] == i) {
                if (nucleos[sites_index] != '\0') scaffold_out += nucleos[sites_index];
                sites_index++;
            }
        } else {
            scaffold_out += ref[i];
        }
    }

    return scaffold_out;
}



// Same as above, but only for C++ code
// The reason for this one is that it can be run in parallel

std::string cpp_retr_scaff(uint scaff_index,
                           uint variant_index,
                           const XPtr<SequenceSet>& ss,
                           const XPtr<VariantSet>& vs) {

    const std::vector<char>& nucleos = vs->variant_info[variant_index].nucleos[scaff_index];
    const std::vector<uint>& sites = vs->variant_info[variant_index].sites[scaff_index];
    const std::string& ref = ss->sequences[scaff_index];

    std::string scaffold_out = "";
    uint sites_index = 0;
    for (uint i = 0; i < ref.size(); i++) {
        if (sites[sites_index] == i) {
            while (sites[sites_index] == i) {
                if (nucleos[sites_index] != '\0') scaffold_out += nucleos[sites_index];
                sites_index++;
            }
        } else {
            scaffold_out += ref[i];
        }
    }

    return scaffold_out;
}


/*
 Retrieve a portion of a variant's scaffold
 */

// [[Rcpp::export]]
std::string variants_retr_seq(const size_t& start_pos,
                              const size_t& length_out,
                              uint scaff_num,
                              uint variant_num,
                              const XPtr<SequenceSet>& ss,
                              const XPtr<VariantSet>& vs) {

    std::string scaff = variants_retr_scaff(scaff_num, variant_num, ss, vs);

    std::string out_string = scaff.substr(start_pos - 1, length_out);

    return out_string;
}




/*
 Retrieve all scaffolds from one variant
 */

// [[Rcpp::export]]
std::vector<std::string> variants_retr_var(const uint& variant_num,
                                           const XPtr<SequenceSet>& ss,
                                           const XPtr<VariantSet>& vs) {

    uint n_scaffs = ss->sequences.size();

    std::vector<std::string> variant_out(n_scaffs);

    for (uint i = 0; i < n_scaffs; i++) {
        variant_out[i] = variants_retr_scaff(i + 1, variant_num, ss, vs);
    }

    return variant_out;
}




// [[Rcpp::export]]
void SummarizeVariantSet(XPtr<VariantSet> vs) {

    Rcout.imbue(std::locale(""));
    Rcout << "# Variants: " << vs->n_variants << std::endl;
    Rcout << "# Segregating sites: " << vs->total_segr_sites << std::endl;
}

