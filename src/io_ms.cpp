/*
 Read info from ms-style output files
 */

#include "jackalope_config.h" // controls debugging and diagnostics output

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
#include "ref_classes.h"  // Ref* classes
#include "hap_classes.h"  // Hap* classes
#include "str_manip.h"  // filter_nucleos
#include "util.h"  // str_stop, thread_check
#include "io.h"

using namespace Rcpp;




// For parsing segregating sites from ms-style output files
namespace parse_ms {
    const std::string site = "segsites:";
    const std::string pos = "positions:";
}

/*
 Parse from gene trees in ms-style output
 */

void ms_parse_tree_line(std::string& line,
                        std::vector<std::vector<std::string>>& newick_strings) {

    if (line[0] == '/' && line[1] == '/') {
        newick_strings.push_back(std::vector<std::string>(0));
        return;
    }
    if (line[0] == '[' || line[0] == '(') {
        if (newick_strings.empty()) {
            str_stop({"\nIn the input ms-style output file containing gene trees, ",
                     "the first gene tree is not preceded with a line containing \"//\"."});
        }
        newick_strings.back().push_back(line);
    }
    return;
}

//' Read a ms output file with newick gene trees and return the gene tree strings.
//'
//' @param ms_file File name of the ms output file.
//'
//' @return A vector of strings for each set of gene trees.
//'
//' @noRd
//'
//[[Rcpp::export]]
std::vector<std::vector<std::string>> read_ms_trees_(std::string ms_file) {

    std::vector<std::vector<std::string>> newick_strings;

    expand_path(ms_file);

    gzFile file;
    file = gzopen(ms_file.c_str(), "rb");
    if (! file) {
        std::string e = "gzopen of " + ms_file + " failed: " + strerror(errno) + ".\n";
        Rcpp::stop(e);
    }

    // Scroll through buffers
    std::string lastline = "";
    char *buffer = new char[LENGTH];

    while (1) {
        Rcpp::checkUserInterrupt();
        int err;
        int bytes_read;
        bytes_read = gzread(file, buffer, LENGTH - 1);
        buffer[bytes_read] = '\0';

        // Recast buffer as a std::string:
        std::string mystring(buffer);
        mystring = lastline + mystring;

        // std::vector of strings for parsed buffer:
        std::vector<std::string> svec = cpp_str_split_newline(mystring);

        // Scroll through lines derived from the buffer.
        for (uint64 i = 0; i < svec.size() - 1; i++){
            ms_parse_tree_line(svec[i], newick_strings);
        }
        // Manage the last line.
        lastline = svec.back();

        // Check for end of file (EOF) or errors.
        if (bytes_read < LENGTH - 1) {
            if ( gzeof(file) ) {
                ms_parse_tree_line(lastline, newick_strings);
                break;
            } else {
                std::string error_string = gzerror(file, &err);
                if (err) {
                    std::string e = "Error: " + error_string + ".\n";
                    stop(e);
                }
            }
        }
    }
    delete[] buffer;
    gzclose (file);

    return newick_strings;
}




/*
 Parse from segregating sites in ms-style output
 */

// For organizing info from each chromosome
struct MS_SitesInfo {
    uint64 n_sites;
    std::vector<double> positions;
    std::vector<std::vector<bool>> segr_bools;

    arma::mat to_mat(const uint64& chrom_i) {

        // Now create and fill the output matrix
        arma::mat M(n_sites, segr_bools.size() + 1);
        if (positions.size() != n_sites) {
            str_stop({"\nIn creation of segregation-sites info ",
                     "for chromosome number ", std::to_string(chrom_i + 1),
                     ", the listed positions for each site (line starting with ",
                     "'positions:') does not have a length that's the same "
                     "as the # sites as given by the line starting with 'segsites:'."});
        }
        M.col(0) = arma::conv_to<arma::vec>::from(positions);
        for (uint64 i = 0; i < segr_bools.size(); i++) {
            if (segr_bools[i].size() != n_sites) {
                str_stop({"\nIn creation of segregation-sites info ",
                         "for chromosome number ", std::to_string(chrom_i + 1),
                         ", the listed number of sites (line starting with ",
                         "'segsites:') does not agree with the number of "
                         "items in the ", std::to_string(i + 1), "th line ",
                         "of segregating sites info (ones filled with 0s and 1s)."});
            }
            for (uint64 j = 0; j < segr_bools[i].size(); j++) {
                M(j, i+1) = static_cast<double>(segr_bools[i][j]);
            }
        }

        return M;
    }
};

// For parsing a single line from the file
void ms_parse_sites_line(std::string& line,
                         std::vector<MS_SitesInfo>& sites_infos) {

    if (line[0] == '0' || line[0] == '1') {
        if (sites_infos.empty()) return; // sometimes it has a header that starts with 1/0
        trimws(line);
        std::vector<std::vector<bool>>& segr_bools(sites_infos.back().segr_bools);
        segr_bools.push_back(std::vector<bool>());
        segr_bools.back().reserve(line.size());
        for (char& c : line) segr_bools.back().push_back(c == '1');
    } else if (line[0] == '/' || line[0] == '/') {
        sites_infos.push_back(MS_SitesInfo());
    } else if (line.compare(0, parse_ms::site.size(), parse_ms::site) == 0) {
        line.erase(0, parse_ms::site.size());
        trimws(line);
        sites_infos.back().n_sites = std::stoi(line);
    } else if (line.compare(0, parse_ms::pos.size(), parse_ms::pos) == 0) {
        line.erase(0, parse_ms::pos.size());
        trimws(line);
        std::vector<std::string> pos_str = cpp_str_split_delim(line, ' ');
        if (sites_infos.empty()) {
            str_stop({"\nIn parsing of segregation-sites info from a file, ",
                     "a line starting with '//' should always appear before ",
                     "one starting with 'positions:'."});
        }
        std::vector<double>& pos(sites_infos.back().positions);
        if (!pos.empty()) stop("multiple positions for the same locus.");
        pos.reserve(pos_str.size());
        for (std::string& ps : pos_str) pos.push_back(std::stod(ps));
    }
    return;
}


//' Read a ms output file with segregating sites and return the matrices of site info.
//'
//' @param ms_file File name of the ms output file.
//'
//' @return A vector of strings for each set of gene trees.
//'
//' @noRd
//'
//[[Rcpp::export]]
arma::field<arma::mat> coal_file_sites(std::string ms_file) {

    std::vector<MS_SitesInfo> sites_infos;

    expand_path(ms_file);

    gzFile file;
    file = gzopen(ms_file.c_str(), "rb");
    if (! file) {
        std::string e = "gzopen of " + ms_file + " failed: " + strerror(errno) + ".\n";
        Rcpp::stop(e);
    }

    // Scroll through buffers
    std::string lastline = "";

    char *buffer = new char[LENGTH];

    while (1) {
        Rcpp::checkUserInterrupt();
        int err;
        int bytes_read;
        bytes_read = gzread(file, buffer, LENGTH - 1);
        buffer[bytes_read] = '\0';

        // Recast buffer as a std::string:
        std::string mystring(buffer);
        mystring = lastline + mystring;

        // std::vector of strings for parsed buffer:
        std::vector<std::string> svec = cpp_str_split_newline(mystring);

        // Scroll through lines derived from the buffer.
        for (uint64 i = 0; i < svec.size() - 1; i++){
            ms_parse_sites_line(svec[i], sites_infos);
        }
        // Manage the last line.
        lastline = svec.back();

        // Check for end of file (EOF) or errors.
        if (bytes_read < LENGTH - 1) {
            if ( gzeof(file) ) {
                ms_parse_sites_line(lastline, sites_infos);
                break;
            } else {
                std::string error_string = gzerror(file, &err);
                if (err) {
                    std::string e = "Error: " + error_string + ".\n";
                    stop(e);
                }
            }
        }
    }
    delete[] buffer;
    gzclose (file);

    arma::field<arma::mat> sites_mats(sites_infos.size());
    for (uint64 i = 0; i < sites_infos.size(); i++) {
        sites_mats[i] = sites_infos[i].to_mat(i);
    }

    return sites_mats;
}


