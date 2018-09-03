
/*
 ********************************************************

 Print methods for RefGenome and VarSet objects

 ********************************************************
 */


#include <RcppArmadillo.h>

#include "sequence_classes.h" // Ref* and Var* classes

using namespace Rcpp;




/*
 Calling `base::options("width")$width`
 */
int get_width() {
    // Obtain environment containing function
    Environment base("package:base");
    // Make function callable from C++
    Function opt_r = base["options"];
    // Call the function and receive its list output
    List width_list = opt_r("width");
    int console_width = width_list["width"];
    return console_width;
}






/*
 ==========================================
 ==========================================

 Printing reference genome info

 ==========================================
 ==========================================
 */

void RefGenome::print() const {

    int console_width = get_width();

    // 32 characters is the narrowest I'll allow
    // (I'd start getting negative lengths and other problems otherwise)
    console_width = std::max(console_width, 32);

    int num_seqs = size();
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
    Rcout << "# Total size: " << total_size << " bp" << std::endl;
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

    for (int i = 0; i < inds.size(); i++) {
        ind_i = inds[i];
        if (ind_i == -1) {
            Rprintf("%-10s %-*s %9s\n", "...", seq_print_len, "...", "...");
            continue;
        }
        const RefSequence& rs(sequences[ind_i]);
        const std::string& name_i(rs.name);
        const std::string& seq_i(rs.nucleos);
        // Print name
        Rprintf("%-10.10s ", name_i.c_str());
        // Print sequence
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
        if (rs.size() > 999999999) {
            Rprintf(" %9.2E", rs.size());
        } else {
            Rprintf(" %9i", rs.size());
        }
        Rcout << std::endl;
    }
}




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






/*
 ==========================================
 ==========================================

 For printing info on a set of variants

 ==========================================
 ==========================================
 */


// For printing info on a set of variants
void VarSet::print() const noexcept {

    uint32 total_muts = 0;
    for (const VarGenome& vg : variants) {
        for (const VarSequence& vs : vg.var_genome) {
            total_muts += vs.mutations.size();
        }
    }

    int console_width = get_width();

    int n_spaces = static_cast<int>(
        std::ceil(static_cast<double>(console_width - 21) / 2)
    );

    for (int i = 0; i < n_spaces; i++) Rcout << ' ';
    Rcout << "<< Variants object >>" << std::endl;

    Rcout.imbue(std::locale(""));
    Rcout << "# Variants: " << size() << std::endl;
    Rcout << "# Mutations: " << total_muts << std::endl;
    Rcout << std::endl;

    n_spaces = static_cast<int>(
        std::ceil(static_cast<double>(console_width - 28) / 2)
    );

    for (int i = 0; i < n_spaces; i++) Rcout << ' ';
    Rcout << "<< Reference genome info: >>" << std::endl;
    reference.print();
}





//' Function to print info on a VarSet.
//'
//' Access `VarSet` class's print method from R.
//'
//' @noRd
//'
//[[Rcpp::export]]
void print_var_set(SEXP var_set_ptr) {
    XPtr<VarSet> var_set(var_set_ptr);
    var_set->print();
    return;
}
