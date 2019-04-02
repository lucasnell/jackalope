#include <RcppArmadillo.h>

#include <fstream>
#include <string>
#include <vector>
#include <zlib.h>

#include "jackal_types.h"  // integer types
#include "seq_classes_ref.h"  // Ref* classes
#include "seq_classes_var.h"  // Var* classes
#include "str_manip.h"  // filter_nucleos
#include "util.h"  // str_stop
#include "read_write.h"

using namespace Rcpp;





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

    while (1) {
        Rcpp::checkUserInterrupt();
        int err;
        int bytes_read;
        char buffer[LENGTH];
        bytes_read = gzread(file, buffer, LENGTH - 1);
        buffer[bytes_read] = '\0';

        // Recast buffer as a std::string:
        std::string mystring(reinterpret_cast<char*>(buffer));
        mystring = lastline + mystring;

        char split = '\n';
        // std::vector of strings for parsed buffer:
        std::vector<std::string> svec = cpp_str_split_delim(mystring, split);

        // Scroll through lines derived from the buffer.
        for (uint32 i = 0; i < svec.size() - 1; i++){
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
    gzclose (file);

    return newick_strings;
}




/*
 Parse from segregating sites in ms-style output
 */

// For organizing info from each sequence
struct MS_SitesInfo {
    uint32 n_sites;
    std::vector<double> positions;
    std::vector<std::vector<bool>> segr_bools;

    arma::mat to_mat(const uint32& seq) {

        // Now create and fill the output matrix
        arma::mat M(n_sites, segr_bools.size() + 1);
        if (positions.size() != n_sites) {
            str_stop({"\nIn creation of segregation-sites info ",
                     "for sequence number ", std::to_string(seq + 1),
                     ", the listed positions for each site (line starting with ",
                     "'positions:') does not have a length that's the same "
                     "as the # sites as given by the line starting with 'segsites:'."});
        }
        M.col(0) = arma::conv_to<arma::vec>::from(positions);
        for (uint32 i = 0; i < segr_bools.size(); i++) {
            if (segr_bools[i].size() != n_sites) {
                str_stop({"\nIn creation of segregation-sites info ",
                         "for sequence number ", std::to_string(seq + 1),
                         ", the listed number of sites (line starting with ",
                         "'segsites:') does not agree with the number of "
                         "items in the ", std::to_string(i + 1), "th line ",
                         "of segregating sites info (ones filled with 0s and 1s)."});
            }
            for (uint32 j = 0; j < segr_bools[i].size(); j++) {
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

    while (1) {
        Rcpp::checkUserInterrupt();
        int err;
        int bytes_read;
        char buffer[LENGTH];
        bytes_read = gzread(file, buffer, LENGTH - 1);
        buffer[bytes_read] = '\0';

        // Recast buffer as a std::string:
        std::string mystring(reinterpret_cast<char*>(buffer));
        mystring = lastline + mystring;

        char split = '\n';
        // std::vector of strings for parsed buffer:
        std::vector<std::string> svec = cpp_str_split_delim(mystring, split);

        // Scroll through lines derived from the buffer.
        for (uint32 i = 0; i < svec.size() - 1; i++){
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
    gzclose (file);

    arma::field<arma::mat> sites_mats(sites_infos.size());
    for (uint32 i = 0; i < sites_infos.size(); i++) {
        sites_mats[i] = sites_infos[i].to_mat(i);
    }

    return sites_mats;
}








/*
 ==================================================================
 ==================================================================

                READ VCF

 ==================================================================
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

    XPtr<VarSet> var_set(new VarSet(*reference, var_names));

    for (uint32 mut_i = 0; mut_i < n_muts; mut_i++) {

        const std::string& ref(ref_seq[mut_i]);
        const std::vector<std::string>& haps(haps_list[mut_i]);
        const uint32& seq_i(seq_inds[mut_i]);

        for (uint32 var_i = 0; var_i < n_vars; var_i++) {

            // If it's blank or if it's the same as the reference, move on:
            if (haps[var_i].size() == 0 || haps[var_i] == ref) continue;

            VarSequence& var_seq((*var_set)[var_i][seq_i]);

            // Else, mutate accordingly:
            Mutation new_mut;
            if (ref.size() <= haps[var_i].size()) {
                /*
                ------------
                substitution and/or insertion
                ------------
                */
                new_mut = Mutation(pos[mut_i], pos[mut_i], haps[var_i]);
            } else {
                /*
                ------------
                deletion
                ------------
                */
                // size modifier:
                sint32 sm = static_cast<sint32>(haps[var_i].size()) -
                    static_cast<sint32>(ref.size());
                /* Position (have to do this bc VCFs include non-deleted nucleotide when
                calculating deltion positions) */
                uint32 np = pos[mut_i] + 1;
                new_mut = Mutation(np, np, sm);
            }
            var_seq.mutations.push_back(new_mut);

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




// ==================================================================
// ==================================================================

//                          READ FASTA - NON-INDEXED

// ==================================================================
// ==================================================================


// Parse one line of input from a file and add to output

void parse_fasta_line(const std::string& line, const bool& cut_names,
                      RefGenome& ref) {

    if (line.find(">") != std::string::npos) {
        std::string name_i = "";
        if (cut_names) {
            std::string::size_type spc = line.find(' ', 2);
            if (spc == std::string::npos) spc = line.size();
            name_i = line.substr(1, spc);
            // Remove any spaces if they exist (they would occur at the beginning)
            name_i.erase(std::remove_if(name_i.begin(), name_i.end(), ::isspace),
                         name_i.end());
        } else {
            name_i = line.substr(1, line.size());
        }
        RefSequence seq(name_i, "");
        ref.sequences.push_back(seq);
    } else {
        ref.sequences.back().nucleos += line;
        ref.total_size += line.size();
    }
    return;
}




/*
 C++ function to add to a RefGenome object from a non-indexed fasta file.
 Does most of the work for `read_fasta_noind` below.
 */
void append_ref_noind(RefGenome& ref,
                      std::string fasta_file,
                      const bool& cut_names,
                      const bool& remove_soft_mask) {

    expand_path(fasta_file);

    gzFile file;
    file = gzopen(fasta_file.c_str(), "rb");
    if (! file) {
        std::string e = "gzopen of " + fasta_file + " failed: " + strerror(errno) + ".\n";
        Rcpp::stop(e);
    }

    // Scroll through buffers
    std::string lastline = "";

    while (1) {
        Rcpp::checkUserInterrupt();
        int err;
        int bytes_read;
        char buffer[LENGTH];
        bytes_read = gzread(file, buffer, LENGTH - 1);
        buffer[bytes_read] = '\0';

        // Recast buffer as a std::string:
        std::string mystring(reinterpret_cast<char*>(buffer));
        mystring = lastline + mystring;

        char split = '\n'; // Must be single quotes!
        // std::vector of strings for parsed buffer:
        std::vector<std::string> svec = cpp_str_split_delim(mystring, split);

        // Scroll through lines derived from the buffer.
        for (uint32 i = 0; i < svec.size() - 1; i++){
            parse_fasta_line(svec[i], cut_names, ref);
        }
        // Manage the last line.
        lastline = svec.back();

        // Check for end of file (EOF) or errors.
        if (bytes_read < LENGTH - 1) {
            if ( gzeof(file) ) {
                parse_fasta_line(lastline, cut_names, ref);
                break;
            } else {
                std::string error_string = gzerror (file, & err);
                if (err) {
                    std::string e = "Error: " + error_string + ".\n";
                    stop(e);
                }
            }
        }
    }
    gzclose (file);

    // Remove weird characters and remove soft masking if desired:
    for (uint32 i = 0; i < ref.size(); i++) {
        filter_nucleos(ref.sequences[i].nucleos, remove_soft_mask);
    }

    return;

}



//' Read a non-indexed fasta file to a \code{RefGenome} object.
//'
//' @param file_names File names of the fasta file(s).
//' @param cut_names Boolean for whether to cut sequence names at the first space.
//'     Defaults to \code{TRUE}.
//' @param remove_soft_mask Boolean for whether to remove soft-masking by making
//'    sequences all uppercase. Defaults to \code{TRUE}.
//'
//' @return Nothing.
//'
//' @noRd
//'
//[[Rcpp::export]]
SEXP read_fasta_noind(const std::vector<std::string>& fasta_files,
                      const bool& cut_names,
                      const bool& remove_soft_mask) {

    XPtr<RefGenome> ref_xptr(new RefGenome(), true);
    RefGenome& ref(*ref_xptr);

    for (const std::string& fasta : fasta_files) {
        append_ref_noind(ref, fasta, cut_names, remove_soft_mask);
    }

    return ref_xptr;
}







// ==================================================================
// ==================================================================

//                          READ FASTA - INDEXED

// ==================================================================
// ==================================================================


// Parse one line of input from a fasta index file and add to output

void parse_line_fai(const std::string& line,
                    std::vector<uint64>& offsets,
                    std::vector<std::string>& names,
                    std::vector<uint64>& lengths,
                    std::vector<uint32>& line_lens) {

    char split = '\t';

    if (line != "") {
        std::vector<std::string> split_line = cpp_str_split_delim(line, split);
        names.push_back(split_line[0]);
        lengths.push_back(std::stoull(split_line[1]));
        offsets.push_back(std::stoull(split_line[2]));
        line_lens.push_back(std::stoul(split_line[3]));
    }
    return;
}


// Get info from a fasta index file

void read_fai(const std::string& fai_file,
              std::vector<uint64>& offsets,
              std::vector<std::string>& names,
              std::vector<uint64>& lengths,
              std::vector<uint32>& line_lens) {


    gzFile file;
    file = gzopen(fai_file.c_str(), "rb");
    if (! file) {
        std::string e = "gzopen of " + fai_file + " failed: " + strerror(errno) +
            ".\n";
        Rcpp::stop(e);
    }

    // Scroll through buffers
    std::string lastline = "";

    while (1) {
        Rcpp::checkUserInterrupt();
        int err;
        int bytes_read;
        char buffer[LENGTH];
        bytes_read = gzread(file, buffer, LENGTH - 1);
        buffer[bytes_read] = '\0';

        // Recast buffer as a std::string:
        std::string mystring(reinterpret_cast<char*>(buffer));
        mystring = lastline + mystring;

        char split = '\n'; // Must be single quotes!
        // std::vector of strings for parsed buffer:
        std::vector<std::string> svec = cpp_str_split_delim(mystring, split);

        // Scroll through lines derived from the buffer.
        for (uint32 i = 0; i < svec.size() - 1; i++){
            parse_line_fai(svec[i], offsets, names, lengths, line_lens);
        }
        // Manage the last line.
        lastline = svec.back();

        // Check for end of file (EOF) or errors.
        if (bytes_read < LENGTH - 1) {
            if ( gzeof(file) ) {
                parse_line_fai(lastline, offsets, names, lengths, line_lens);
                break;
            } else {
                std::string error_string = gzerror (file, & err);
                if (err) {
                    std::string e = "Error: " + error_string + ".\n";
                    stop(e);
                }
            }
        }
    }
    gzclose (file);

    return;
}






/*
 C++ function to add to a RefGenome object from an indexed fasta file.
 Does most of the work for `read_fasta_ind` below.
 */
void append_ref_ind(RefGenome& ref,
                    std::string fasta_file,
                    std::string fai_file,
                    const bool& remove_soft_mask) {

    std::vector<uint64> offsets;
    std::vector<std::string> names;
    std::vector<uint64> lengths;
    std::vector<uint32> line_lens;

    expand_path(fasta_file);
    expand_path(fai_file);

    // Fill info from index file:
    read_fai(fai_file, offsets, names, lengths, line_lens);

    if (offsets.size() != names.size() || names.size() != lengths.size() ||
        lengths.size() != line_lens.size()) {
        stop("Wrong sizes.");
    }

    // Process fasta file:
    gzFile file;
    file = gzopen(fasta_file.c_str(), "rb");
    if (! file) {
        std::string e = "gzopen of " + fasta_file + " failed: " + strerror(errno) + ".\n";
        Rcpp::stop(e);
    }

    const uint32 n_seqs0 = ref.size(); // starting # sequences
    uint32 n_new_seqs = offsets.size();
    uint64 LIMIT = 4194304;
    ref.sequences.resize(n_seqs0 + n_new_seqs, RefSequence());

    for (uint32 i = 0; i < n_new_seqs; i++) {

        Rcpp::checkUserInterrupt();

        RefSequence& rs(ref.sequences[i+n_seqs0]);
        rs.name = names[i];

        // Length of the whole sequence including newlines
        uint64 len = lengths[i] + lengths[i] / line_lens[i] + 1;

        sint64 bytes_read;

        for (uint64 j = 0; j < len; j += (LIMIT-1)) {
            gzseek(file, offsets[i] + j, SEEK_SET);
            uint32 partial_len = LIMIT;
            if (len - j < LIMIT) partial_len = len - j;
            char buffer[partial_len];
            bytes_read = gzread(file, buffer, partial_len - 1);
            buffer[bytes_read] = '\0';

            // Recast buffer as a std::string:
            std::string seq_str(static_cast<char*>(buffer));

            // Remove newlines
            seq_str.erase(remove(seq_str.begin(), seq_str.end(), '\n'),
                          seq_str.end());

            // Filter out weird characters and remove soft masking if requested
            filter_nucleos(seq_str, remove_soft_mask);

            rs.nucleos += seq_str;
            ref.total_size += seq_str.size();

            // Check for errors.
            if (bytes_read < partial_len) {
                if ( gzeof(file) ) {
                    warning("fai file lengths appear incorrect; re-index or "
                                "check output manually for accuracy");
                    break;
                } else {
                    int err;
                    std::string error_string = gzerror(file, &err);
                    if (err) {
                        std::string e = "Error: " + error_string + ".\n";
                        stop(e);
                    }
                }
            }
        }

    }

    gzclose(file);

    return;
}

//' Read an indexed fasta file to a \code{RefGenome} object.
//'
//' @param file_name File name of the fasta file.
//' @param remove_soft_mask Boolean for whether to remove soft-masking by making
//'    sequences all uppercase. Defaults to \code{TRUE}.
//' @param offsets Vector of sequence offsets from the fasta index file.
//' @param names Vector of sequence names from the fasta index file.
//' @param lengths Vector of sequence lengths from the fasta index file.
//' @param line_lens Vector of sequence line lengths from the fasta index file.
//'
//' @return Nothing.
//'
//' @noRd
//'
//'
//[[Rcpp::export]]
SEXP read_fasta_ind(const std::vector<std::string>& fasta_files,
                    const std::vector<std::string>& fai_files,
                    const bool& remove_soft_mask) {

    XPtr<RefGenome> ref_xptr(new RefGenome(), true);
    RefGenome& ref(*ref_xptr);

    if (fasta_files.size() != fai_files.size()) {
        str_stop({"\nThe vector of fasta index files must be the same length as ",
                 "the vector of fasta files."});
    }

    for (uint32 i = 0; i < fasta_files.size(); i++) {
        append_ref_ind(ref, fasta_files[i], fai_files[i], remove_soft_mask);
    }

    return ref_xptr;

}





// ==================================================================
// ==================================================================

//                          WRITE

// ==================================================================
// ==================================================================






//' Write \code{RefGenome} to an uncompressed fasta file.
//'
//' @param file_name File name of output fasta file.
//' @param ref_genome_ptr An external pointer to a \code{RefGenome} C++ object.
//' @param text_width The number of characters per line in the output fasta file.
//'
//' @return Nothing.
//'
//' @noRd
//'
//'
//[[Rcpp::export]]
void write_fasta_fa(std::string file_name,
                    SEXP ref_genome_ptr,
                    const uint32& text_width){

    XPtr<RefGenome> ref_xptr(ref_genome_ptr);
    RefGenome& ref(*ref_xptr);

    expand_path(file_name);

    std::ofstream out_file(file_name);

    if (out_file.is_open()) {

        for (uint32 i = 0; i < ref.size(); i++) {
            out_file << '>';
            out_file << ref[i].name;
            out_file << '\n';

            const std::string& seq_str(ref[i].nucleos);

            for (uint32 pos = 0; pos < seq_str.size(); pos++) {
                out_file << seq_str[pos];
                if ((pos % text_width) == (text_width - 1)) out_file << '\n';
            }
            out_file << '\n';
        }
        out_file.close();

    } else {
        Rcout << "Unable to open file " << file_name << std::endl;
    }

    return;
}




//' Write \code{RefGenome} to a compressed fasta file.
//'
//' @inheritParams write_fasta_fa
//'
//' @return Nothing.
//'
//' @noRd
//'
//[[Rcpp::export]]
void write_fasta_gz(std::string file_name,
                    SEXP ref_genome_ptr,
                    const uint32& text_width){

    XPtr<RefGenome> ref_xptr(ref_genome_ptr);
    RefGenome& ref(*ref_xptr);

    expand_path(file_name);

    // Initialize filehandle.
    gzFile fi;

    // Initialize file.
    // Note that gzfile does not tolerate initializing an empty file.
    // Use ofstream instead.
    if (!std::ifstream(file_name)){
        std::ofstream myfile;
        myfile.open(file_name, std::ios::out | std::ios::binary);
        myfile.close();
    }


    fi = gzopen(file_name.c_str(), "wb");
    if (!fi) {
        std::string e = "gzopen of " + file_name + " failed: " + strerror (errno) + ".\n";
        Rcpp::stop(e);
    }

    for (uint32 i = 0; i < ref.size(); i++) {
        std::string name = '>' + ref[i].name + '\n';
        gzwrite(fi, name.c_str(), name.size());

        const std::string& seq_str(ref[i].nucleos);
        uint32 num_rows = seq_str.length() / text_width;

        for (uint32 i = 0; i < num_rows; i++) {
            std::string one_line = seq_str.substr(i * text_width, text_width);
            one_line += '\n';
            gzwrite(fi, one_line.c_str(), one_line.size());
        }

        // If there are leftover characters, create a shorter item at the end.
        if (seq_str.length() % text_width != 0) {
            std::string one_line = seq_str.substr(text_width * num_rows);
            one_line += '\n';
            gzwrite(fi, one_line.c_str(), one_line.size());
        }
    }
    gzclose(fi);

    return;
}





//' Write `variants` to VCF file.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
void write_vcf_cpp(std::string out_prefix,
                   const bool& compress,
                   SEXP var_set_ptr,
                   const IntegerMatrix& sample_matrix,
                   const bool& testing) {

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


    /*
     The `testing` argument overrides everything else.
     It results in all output being streamed to Rcout.
     */
    if (testing) {

        std::string out_file = "";

        write_vcf_<std::string>(var_set, out_file, writer);

        return;
    }

    if (compress) {
        std::string file_name = out_prefix + ".vcf.gz";
        /*
         Initialize file.
         Note that gzfile does not tolerate initializing an empty file.
         Use ofstream instead.
         */
        if (!std::ifstream(file_name)){
            std::ofstream tmp_file;
            tmp_file.open(file_name, std::ios::out | std::ios::binary);
            tmp_file.close();
        }

        gzFile out_file = gzopen(file_name.c_str(), "wb");
        if (!out_file) {
            std::string e = "gzopen of " + file_name + " failed: " +
                strerror(errno) + ".\n";
            Rcpp::stop(e);
        }

        write_vcf_<gzFile>(var_set, out_file, writer);

        gzclose(out_file);

    } else {

        std::string file_name = out_prefix + ".vcf";
        std::ofstream out_file(file_name);

        if (!out_file.is_open()) {
            Rcout << "Unable to open file " << file_name << std::endl;
        }

        write_vcf_<std::ofstream>(var_set, out_file, writer);

        out_file.close();
    }


    return;

}

