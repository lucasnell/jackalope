#include <RcppArmadillo.h>

#include <string>
#include <vector>
#include <sstream>  // ostringstream
#include <iomanip> // setting precision on ostringstream


#include "jackalope_types.h"  // integer types
#include "util.h"  // str_stop
#include "io.h"  // File* classes

using namespace Rcpp;




template <typename T>
void write_bed__(const std::string& file_name,
                 const std::vector<arma::mat>& gamma_mats,
                 const std::vector<std::string>& seq_names,
                 const int& compress) {

    // Open file:
    T out_file(file_name, compress);

    // ostringstream to make sure doubles are printed with many digits
    std::ostringstream dbl_stream;
    // Set fixed-point notation and precision
    dbl_stream << std::fixed << std::setprecision(12);

    // Write to file
    std::string line;
    line.reserve(500);  // <- should be more than enough
    for (uint64 i = 0; i < gamma_mats.size(); i++) {
        if (gamma_mats[i].n_rows == 0) continue;
        const arma::mat& gm(gamma_mats[i]);
        // First line's a bit different:
        line = seq_names[i] + '\t';
        line += "0\t";
        line += std::to_string(static_cast<int>(gm(0,0))) + '\t';
        line += seq_names[i] + "_0" + '\t';
        dbl_stream << gm(0,1);
        line += dbl_stream.str() + '\n';
        dbl_stream.str(std::string()); // clear stream
        out_file.write(line);
        // The rest of the lines:
        for (uint64 j = 1; j < gm.n_rows; j++) {
            line = seq_names[i] + '\t';
            line += std::to_string(static_cast<int>(gm(j-1,0))) + '\t';
            line += std::to_string(static_cast<int>(gm(j,0))) + '\t';
            line += seq_names[i] + '_' + std::to_string(j) + '\t';
            dbl_stream << gm(j,1);
            line += dbl_stream.str() + '\n';
            dbl_stream.str(std::string()); // clear stream
            out_file.write(line);
        }
    }

    // Close file:
    out_file.close();

    return;

}


//' Write Gamma matrix info to a tab-delimited BED file.
//'
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
void write_bed(std::string out_prefix,
               const std::vector<arma::mat>& gamma_mats,
               const std::vector<std::string>& seq_names,
               const int& compress,
               const std::string& comp_method) {

    if (seq_names.size() != gamma_mats.size()) {
        str_stop({"\nIn internal function `write_bed`, the list of gamma matrices ",
                 "is not the same length as the vector of sequence names. ",
                 "Since the check for this is already done in `site_var`, ",
                 "I believe that something very fishy is going on..."});
    }

    expand_path(out_prefix);

    std::string file_name = out_prefix + ".bed";

    if (compress > 0) {

        if (comp_method == "gzip") {
            write_bed__<FileGZ>(file_name, gamma_mats, seq_names, compress);
        } else if (comp_method == "bgzip") {
            write_bed__<FileBGZF>(file_name, gamma_mats, seq_names, compress);
        } else stop("\nUnrecognized compression method.");

    } else {

        write_bed__<FileUncomp>(file_name, gamma_mats, seq_names, compress);

    }

    return;
}





