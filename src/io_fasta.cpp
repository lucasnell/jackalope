
/*
 Functions to read and write to/from FASTA files
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
#include "io.h"   // expand_path, File* classes, `LENGTH`

using namespace Rcpp;



/*
 ==================================================================
 ==================================================================

 READ FASTA - NON-INDEXED

 ==================================================================
 ==================================================================
 */

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
        RefChrom chrom(name_i, "");
        ref.chromosomes.push_back(chrom);
    } else {
        ref.chromosomes.back().nucleos += line;
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
    delete[] buffer;
    gzclose (file);

    // Remove weird characters and remove soft masking if desired:
    for (uint64 i = 0; i < ref.size(); i++) {
        filter_nucleos(ref.chromosomes[i].nucleos, remove_soft_mask);
    }

    return;

}



//' Read a non-indexed fasta file to a \code{RefGenome} object.
//'
//' @param file_names File names of the fasta file(s).
//' @param cut_names Boolean for whether to cut chromosome names at the first space.
//'     Defaults to \code{TRUE}.
//' @param remove_soft_mask Boolean for whether to remove soft-masking by making
//'    chromosomes all uppercase. Defaults to \code{TRUE}.
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
                    std::vector<uint64>& line_lens) {

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
              std::vector<uint64>& line_lens) {


    gzFile file;
    file = gzopen(fai_file.c_str(), "rb");
    if (! file) {
        std::string e = "gzopen of " + fai_file + " failed: " + strerror(errno) +
            ".\n";
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

    delete[] buffer;
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
    std::vector<uint64> line_lens;

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

    const uint64 n_chroms0 = ref.size(); // starting # chromosomes
    uint64 n_new_chroms = offsets.size();
    uint64 LIMIT = 4194304;
    ref.chromosomes.resize(n_chroms0 + n_new_chroms, RefChrom());

    for (uint64 i = 0; i < n_new_chroms; i++) {

        Rcpp::checkUserInterrupt();

        RefChrom& rs(ref.chromosomes[i+n_chroms0]);
        rs.name = names[i];

        // Length of the whole chromosome including newlines
        uint64 len = lengths[i] + lengths[i] / line_lens[i] + 1;

        sint64 bytes_read;

        for (uint64 j = 0; j < len; j += (LIMIT-1)) {
            gzseek(file, offsets[i] + j, SEEK_SET);
            uint64 partial_len = LIMIT;
            if (len - j < LIMIT) partial_len = len - j;
            char *buffer = new char[partial_len];
            bytes_read = gzread(file, buffer, partial_len - 1);
            buffer[bytes_read] = '\0';

            // Recast buffer as a std::string:
            std::string chrom_str(buffer);

            // Remove newlines
            chrom_str.erase(remove(chrom_str.begin(), chrom_str.end(), '\n'),
                          chrom_str.end());

            // Filter out weird characters and remove soft masking if requested
            filter_nucleos(chrom_str, remove_soft_mask);

            rs.nucleos += chrom_str;
            ref.total_size += chrom_str.size();

            // Check for errors.
            if (bytes_read < static_cast<sint64>(partial_len)) {
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

            delete[] buffer;


        }

    }

    gzclose(file);

    return;
}

//' Read an indexed fasta file to a \code{RefGenome} object.
//'
//' @param file_name File name of the fasta file.
//' @param remove_soft_mask Boolean for whether to remove soft-masking by making
//'    chromosomes all uppercase. Defaults to \code{TRUE}.
//' @param offsets Vector of chromosome offsets from the fasta index file.
//' @param names Vector of chromosome names from the fasta index file.
//' @param lengths Vector of chromosome lengths from the fasta index file.
//' @param line_lens Vector of chromosome line lengths from the fasta index file.
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

    for (uint64 i = 0; i < fasta_files.size(); i++) {
        append_ref_ind(ref, fasta_files[i], fai_files[i], remove_soft_mask);
    }

    return ref_xptr;

}





// ==================================================================
// ==================================================================

//                          WRITE

// ==================================================================
// ==================================================================




/*
 Template that does most of the work to write from RefGenome to FASTA files of
 varying formats (gzip, bgzip, uncompressed).
 `T` should be `FileUncomp` or `FileBGZF` from `io.h`
 */
template <typename T>
inline void write_ref_fasta__(const std::string& file_name,
                              const int& compress,
                              const RefGenome& ref,
                              const uint64& text_width,
                              const bool& show_progress) {

    T file(file_name, compress);

    Progress prog_bar(ref.total_size, show_progress);

    std::string one_line;
    one_line.reserve(text_width + 2);

    for (uint64 i = 0; i < ref.size(); i++) {

        if (prog_bar.check_abort()) break;

        std::string name = '>' + ref[i].name + '\n';
        file.write(name);

        const std::string& chrom_str(ref[i].nucleos);
        uint64 num_rows = chrom_str.length() / text_width;
        uint64 n_chars = 0;

        for (uint64 i = 0; i < num_rows; i++) {
            // Check every 10,000 characters for user interrupt:
            if (n_chars > 10000) {
                if (prog_bar.check_abort()) break;
                n_chars = 0;
            }
            one_line = chrom_str.substr(i * text_width, text_width);
            one_line += '\n';
            file.write(one_line);
            n_chars += text_width;
        }

        if (prog_bar.is_aborted() || prog_bar.check_abort()) break;

        // If there are leftover characters, create a shorter item at the end.
        if (chrom_str.length() % text_width != 0) {
            one_line = chrom_str.substr(text_width * num_rows);
            one_line += '\n';
            file.write(one_line);
        }

        prog_bar.increment(chrom_str.size());

    }

    file.close();

    return;
}

//' Write \code{RefGenome} to an uncompressed fasta file.
//'
//' @param out_prefix Prefix to file name of output fasta file.
//' @param ref_genome_ptr An external pointer to a \code{RefGenome} C++ object.
//' @param text_width The number of characters per line in the output fasta file.
//' @param compress Boolean for whether to compress output.
//'
//' @return Nothing.
//'
//' @noRd
//'
//'
//[[Rcpp::export]]
void write_ref_fasta(const std::string& out_prefix,
                     SEXP ref_genome_ptr,
                     const uint64& text_width,
                     const int& compress,
                     const std::string& comp_method,
                     const bool& show_progress) {

    XPtr<RefGenome> ref_xptr(ref_genome_ptr);
    RefGenome& ref(*ref_xptr);

    std::string file_name = out_prefix + ".fa";

    expand_path(file_name);

    if (compress > 0) {

        if (comp_method == "gzip") {
            write_ref_fasta__<FileGZ>(file_name, compress, ref, text_width,
                                      show_progress);
        } else if (comp_method == "bgzip") {
            write_ref_fasta__<FileBGZF>(file_name, compress, ref, text_width,
                                        show_progress);
        } else stop("\nUnrecognized compression method.");

    } else {
        write_ref_fasta__<FileUncomp>(file_name, compress, ref, text_width,
                                      show_progress);
    }

    return;
}






template <typename T>
void write_haps_fasta__(const std::string& out_prefix,
                        const HapSet& hap_set,
                        const uint64& text_width,
                        const int& compress,
                        const uint64& n_threads,
                        const bool& show_progress) {

    Progress prog_bar(hap_set.reference->size() * hap_set.size(), show_progress);

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads) if (n_threads > 1)
{
#endif
    std::string line;
    line.reserve(text_width + 1);
    std::string name;
    name.reserve(text_width + 1);

    // Parallelize the Loop
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (uint64 v = 0; v < hap_set.size(); v++) {

        if (prog_bar.is_aborted() || prog_bar.check_abort()) continue;

        std::string file_name = out_prefix + "__" + hap_set[v].name + ".fa";
        T out_file(file_name, compress);

        for (uint64 s = 0; s < hap_set.reference->size(); s++) {

            if (prog_bar.is_aborted() || prog_bar.check_abort()) break;

            name = '>';
            name += (*hap_set.reference)[s].name;
            name += '\n';
            out_file.write(name);

            const HapChrom& hap_chrom(hap_set[v][s]);
            uint64 mut_i = 0;
            uint64 line_start = 0;
            uint64 n_chars = 0;

            while (line_start < hap_chrom.chrom_size) {
                // Check every 10,000 characters for user interrupt:
                if (n_chars > 10000) {
                    if (prog_bar.check_abort()) break;
                    n_chars = 0;
                }
                hap_chrom.set_chrom_chunk(line, line_start,
                                      text_width, mut_i);
                line += '\n';
                out_file.write(line);
                line_start += text_width;
                n_chars += text_width;
            }

            prog_bar.increment(hap_set.reference->operator[](s).size());

        }

        out_file.close();

    }

#ifdef _OPENMP
}
#endif

}



//' Write \code{HapSet} to an uncompressed fasta file.
//'
//' @param out_prefix Prefix to file name of output fasta file.
//' @param hap_set_ptr An external pointer to a \code{HapSet} C++ object.
//' @param text_width The number of characters per line in the output fasta file.
//' @param compress Boolean for whether to compress output.
//'
//' @return Nothing.
//'
//' @noRd
//'
//'
//[[Rcpp::export]]
void write_haps_fasta(std::string out_prefix,
                      SEXP hap_set_ptr,
                      const uint64& text_width,
                      const int& compress,
                      const std::string& comp_method,
                      uint64 n_threads,
                      const bool& show_progress) {

    XPtr<HapSet> haps_xptr(hap_set_ptr);
    HapSet& hap_set(*haps_xptr);

    // Check that # threads isn't too high and change to 1 if not using OpenMP
    thread_check(n_threads);

    expand_path(out_prefix);

    if (compress > 0) {

        if (comp_method == "gzip") {
            write_haps_fasta__<FileGZ>(out_prefix, hap_set, text_width, compress,
                                       n_threads, show_progress);
        } else if (comp_method == "bgzip") {
            write_haps_fasta__<FileBGZF>(out_prefix, hap_set, text_width, compress,
                                         n_threads, show_progress);
        } else stop("\nUnrecognized compression method.");

    } else {

        write_haps_fasta__<FileUncomp>(out_prefix, hap_set, text_width, compress,
                                       n_threads, show_progress);

    }


    return;
}






