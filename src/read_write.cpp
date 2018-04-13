#include <RcppArmadillo.h>

#include <fstream>
#include <string>
#include <vector>
#include <zlib.h>

#include "sequence_classes.h"  // RefGenome and RefSequence classes
#include "str_manip.h"
#include "read_write.h"

using namespace Rcpp;




/*
 Calling `base::path.expand(file_name)` and changing input string
*/
void expand_path(std::string& file_name) {
    // Obtain environment containing function
    Environment base("package:base");
    // Make function callable from C++
    Function pe_r = base["path.expand"];
    // Call the function and receive its list output
    SEXP fai_file_exp = pe_r(file_name);
    file_name = as<std::string>(fai_file_exp);
    return;
}


// Size of the block of memory to use for reading non-indexed fasta files and fai files.
#define LENGTH 0x1000 // hexadecimel for 4096.



// ==================================================================
// ==================================================================

//                          READ - NON-INDEXED

// ==================================================================
// ==================================================================


// Parse one line of input from a file and add to output

void parse_line(const std::string& line, const bool& cut_names,
                RefGenome& ref) {

    if (line.find(">") != std::string::npos) {
        std::string name_i = "";
        if (cut_names) {
            std::string::size_type spc = line.find(' ', 2);
            if (spc == std::string::npos) spc = line.size();
            name_i = line.substr(1, spc);
            // Remove any spaces if they exist (they would occur at the beginning)
            name_i.erase(remove_if(name_i.begin(), name_i.end(), ::isspace),
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
 C++ function to fill an empty RefGenome object from a non-indexed fasta file.
 Does most of the work for `read_fasta_noind` below.
 */
void fill_ref_noind(RefGenome& ref,
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
        for (uint i = 0; i < svec.size() - 1; i++){
            parse_line(svec[i], cut_names, ref);
        }
        // Manage the last line.
        lastline = svec.back();

        // Check for end of file (EOF) or errors.
        if (bytes_read < LENGTH - 1) {
            if ( gzeof(file) ) {
                parse_line(lastline, cut_names, ref);
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

    if (remove_soft_mask) {
        for (int i = 0; i < ref.size(); i++) {
            cpp_to_upper(ref.sequences[i].nucleos);
        }
    }

    return;

}



//' Read a non-indexed fasta file to a \code{RefGenome} object.
//'
//' @param file_name File name of the fasta file.
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
SEXP read_fasta_noind(const std::string& fasta_file,
                      const bool& cut_names,
                      const bool& remove_soft_mask) {

    XPtr<RefGenome> ref_xptr(new RefGenome(), true);
    RefGenome& ref(*ref_xptr);

    fill_ref_noind(ref, fasta_file, cut_names, remove_soft_mask);

    return ref_xptr;
}






// ==================================================================
// ==================================================================

//                          READ - INDEXED

// ==================================================================
// ==================================================================


// Parse one line of input from a fasta index file and add to output

void parse_line_fai(const std::string& line,
                    std::vector<uint64>& offsets,
                    std::vector<std::string>& names,
                    std::vector<uint64>& lengths,
                    std::vector<uint>& line_lens) {

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
              std::vector<uint>& line_lens) {


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
        for (uint i = 0; i < svec.size() - 1; i++){
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
 C++ function to fill an empty RefGenome object from an indexed fasta file.
 Does most of the work for `read_fasta_ind` below.
 */
void fill_ref_ind(RefGenome& ref,
                  std::string fasta_file,
                  std::string fai_file,
                  const bool& remove_soft_mask) {

    std::vector<uint64> offsets;
    std::vector<std::string> names;
    std::vector<uint64> lengths;
    std::vector<uint> line_lens;

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

    uint n_seqs = offsets.size();
    uint64 LIMIT = 4194304;
    ref.sequences = std::deque<RefSequence>(n_seqs, RefSequence());

    for (uint i = 0; i < n_seqs; i++) {

        Rcpp::checkUserInterrupt();

        RefSequence& rs(ref.sequences[i]);
        rs.name = names[i];

        // Length of the whole sequence including newlines
        uint64 len = lengths[i] + lengths[i] / line_lens[i] + 1;

        uint64 bytes_read;

        for (uint64 j = 0; j < len; j += (LIMIT-1)) {
            gzseek(file, offsets[i] + j, SEEK_SET);
            uint partial_len = LIMIT;
            if (len - j < LIMIT) partial_len = len - j;
            char buffer[partial_len];
            bytes_read = gzread(file, buffer, partial_len - 1);
            buffer[bytes_read] = '\0';

            // Recast buffer as a std::string:
            std::string seq_str(static_cast<char*>(buffer));

            // Remove newlines
            seq_str.erase(remove(seq_str.begin(), seq_str.end(), '\n'),
                          seq_str.end());

            if (remove_soft_mask) cpp_to_upper(seq_str);

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
    gzclose (file);
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
SEXP read_fasta_ind(const std::string& fasta_file,
                    const std::string& fai_file,
                    const bool& remove_soft_mask) {

    XPtr<RefGenome> ref_xptr(new RefGenome(), true);
    RefGenome& ref(*ref_xptr);

    fill_ref_ind(ref, fasta_file, fai_file, remove_soft_mask);

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
//' @param ref_ An external pointer to a \code{RefGenome} C++ object.
//' @param text_width The number of characters per line in the output fasta file.
//'
//' @return Nothing.
//'
//' @noRd
//'
//'
//[[Rcpp::export]]
void write_fasta_fa(std::string file_name,
                    SEXP ref_,
                    const uint& text_width){

    XPtr<RefGenome> ref_xptr(ref_);
    RefGenome& ref(*ref_xptr);

    std::ofstream out_file(file_name);

    if (out_file.is_open()) {

        for (uint i = 0; i < ref.size(); i++) {
            out_file << '>';
            out_file << ref[i].name;
            out_file << '\n';

            const std::string& seq_str(ref[i].nucleos);

            for (uint pos = 0; pos < seq_str.size(); pos++) {
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
void write_fasta_gz(const std::string& file_name,
                    SEXP ref_,
                    const uint& text_width){

    XPtr<RefGenome> ref_xptr(ref_);
    RefGenome& ref(*ref_xptr);

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

    for (uint i = 0; i < ref.size(); i++) {
        std::string name = '>' + ref[i].name + '\n';
        gzwrite(fi, name.c_str(), name.size());

        const std::string& seq_str(ref[i].nucleos);
        uint num_rows = seq_str.length() / text_width;

        for (uint i = 0; i < num_rows; i++) {
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



