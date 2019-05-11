#ifndef __JACKAL_SEQ_CLASSES_REF_H
#define __JACKAL_SEQ_CLASSES_REF_H


/*
 ********************************************************

 Classes to store reference-genome sequence info

 ********************************************************
 */



#include <RcppArmadillo.h>
#include <vector>  // vector class
#include <string>  // string class
#include <deque>  // deque class

#include "jackalope_types.h"  // integer types
#include "util.h"  // clear_memory, get_width

using namespace Rcpp;




/*
 =========================================
 One reference-genome sequence (e.g., chromosome, scaffold)
 =========================================
 */
struct RefSequence {

    // Member variables
    std::string name;
    std::string nucleos;

    // Constructors
    RefSequence() : name(""), nucleos("") {};
    RefSequence(const std::string& name_, const std::string& nucleos_)
        : name(name_), nucleos(nucleos_) {};
    RefSequence(const std::string& nucleos_)
        : name(""), nucleos(nucleos_) {};

    // Overloaded operator so nucleotides can be easily extracted
    char operator[](const uint64& idx) const {
        if (idx >= nucleos.size()) {
            stop("Trying to extract nucleotide that doesn't exist");
        }
        return nucleos[idx];
    }
    char& operator[](const uint64& idx) {
        if (idx >= nucleos.size()) {
            stop("Trying to extract nucleotide that doesn't exist");
        }
        return nucleos[idx];
    }
    // To resize this sequence
    void reserve(const uint64& n) {
        nucleos.reserve(n);
        return;
    }
    // To resize this sequence
    void resize(const uint64& n, const char& x) {
        // nucleos.resize(n, x); // the below way should be faster based on testing
        nucleos = std::string(n, x);
        return;
    }
    // To add character to sequence
    void push_back(const char& nt) {
        nucleos.push_back(nt);
        return;
    }
    // To add character to sequence
    void operator+=(const char& nt) {
        nucleos += nt;
        return;
    }
    // To return the size of this sequence
    uint64 size() const noexcept {
        return nucleos.size();
    }
    // For sorting from largest to smallest sequence
    bool operator > (const RefSequence& other) const noexcept {
        return size() > other.size();
    }

    /*
     ------------------
     For filling a read at a given starting position from a sequence of a
     given starting position and size.
     Used only for sequencer.
     ------------------
     */
    void fill_read(std::string& read,
                   const uint64& read_start,
                   const uint64& seq_start,
                   uint64 n_to_add) const {
        // Making sure end doesn't go beyond the sequence bounds
        if ((seq_start + n_to_add - 1) >= nucleos.size()) {
            n_to_add = nucleos.size() - seq_start;
        }
        // Make sure the read is long enough (this fxn should never shorten it):
        if (read.size() < n_to_add + read_start) read.resize(n_to_add + read_start, 'N');
        for (uint64 i = 0; i < n_to_add; i++) {
            read[(read_start + i)] = this->nucleos[(seq_start + i)];
        }
        return;
    }

};




/*
 =========================================
 One full reference genome
 =========================================
 */

struct RefGenome {

    // Member variables
    uint64 total_size = 0;
    std::deque<RefSequence> sequences;
    bool merged = false;
    // For storing original names if merged:
    std::deque<std::string> old_names;
    // Only added for compatibility with templates in sequencing code:
    std::string name = "REF";

    // Constructors
    RefGenome()
        : sequences(std::deque<RefSequence>(0)) {};
    RefGenome(const RefGenome& ref_)
        : total_size(ref_.total_size), sequences(ref_.sequences),
          merged(ref_.merged), old_names(ref_.old_names) {};
    RefGenome(const uint64& N)
        : sequences(std::deque<RefSequence>(N, RefSequence())) {};
    RefGenome(const std::deque<std::string>& seqs) {
        uint64 n_seqs = seqs.size();
        sequences = std::deque<RefSequence>(n_seqs, RefSequence());
        for (uint64 i = 0; i < n_seqs; i++) {
            sequences[i].nucleos = seqs[i];
            sequences[i].name = "seq" + std::to_string(i);
            total_size += seqs[i].size();
        }
    }
    // Overloaded operator so sequences can be easily extracted
    // It returns a reference so no copying is done and so changes can be made
    RefSequence& operator[](const uint64& idx) {
        if (idx >= sequences.size()) {
            stop("Trying to extract sequence that doesn't exist");
        }
        return sequences[idx];
    }
    const RefSequence& operator[](const uint64& idx) const {
        if (idx >= sequences.size()) {
            stop("Trying to extract sequence that doesn't exist");
        }
        return sequences[idx];
    }
    // To return the number of sequences
    uint64 size() const noexcept {
        return sequences.size();
    }
    // To return the sequence sizes
    std::vector<uint64> seq_sizes() const {
        std::vector<uint64> out(size());
        for (uint64 i = 0; i < out.size(); i++) out[i] = sequences[i].size();
        return out;
    }
    // For printing reference genome info
    void print() const {

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

        Rcout << "< Set of " << big_int_format<int>(num_seqs) << " sequences >";
        Rcout << std::endl;
        Rcout << "# Total size: " << big_int_format<uint64>(total_size) << " bp";
        Rcout << std::endl;

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

        for (int i = 0; i < static_cast<int>(inds.size()); i++) {
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
            int seq_i_size = static_cast<int>(seq_i.size());
            if (seq_i_size > seq_print_len){
                for (int j = 0; j < before_elips; j++) Rcout << seq_i[j];
                Rcout << "...";
                for (int j = (seq_i_size - after_elips); j < seq_i_size; j++) {
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

};




#endif
