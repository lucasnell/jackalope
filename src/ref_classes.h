#ifndef __JACKALOPE_REF_CLASSES_H
#define __JACKALOPE_REF_CLASSES_H


/*
 ********************************************************

 Classes to store reference-genome chromosome info

 ********************************************************
 */



#include "jackalope_config.h" // controls debugging and diagnostics output

#include <RcppArmadillo.h>
#include <vector>  // vector class
#include <string>  // string class
#include <deque>  // deque class

#include "jackalope_types.h"  // integer types
#include "util.h"  // clear_memory, get_width

using namespace Rcpp;




/*
 =========================================
 One reference-genome chromosome (e.g., chromosome, scaffold)
 =========================================
 */
struct RefChrom {

    // Member variables
    std::string name;
    std::string nucleos;

    // Constructors
    RefChrom() : name(""), nucleos("") {};
    RefChrom(const std::string& name_, const std::string& nucleos_)
        : name(name_), nucleos(nucleos_) {};
    RefChrom(const std::string& nucleos_)
        : name(""), nucleos(nucleos_) {};

    // Overloaded operator so nucleotides can be easily extracted
    char operator[](const uint64& idx) const {
#ifdef __JACKALOPE_DEBUG
        if (idx >= nucleos.size()) {
            stop("Trying to extract nucleotide that doesn't exist");
        }
#endif
        return nucleos[idx];
    }
    char& operator[](const uint64& idx) {
#ifdef __JACKALOPE_DEBUG
        if (idx >= nucleos.size()) {
            stop("Trying to extract nucleotide that doesn't exist");
        }
#endif
        return nucleos[idx];
    }
    // To resize this chromosome
    void reserve(const uint64& n) {
        nucleos.reserve(n);
        return;
    }
    // To resize this chromosome
    void resize(const uint64& n, const char& x) {
        // nucleos.resize(n, x); // the below way should be faster based on testing
        nucleos = std::string(n, x);
        return;
    }
    // To add character to chromosome
    void push_back(const char& nt) {
        nucleos.push_back(nt);
        return;
    }
    // To add character to chromosome
    void operator+=(const char& nt) {
        nucleos += nt;
        return;
    }
    // To return the size of this chromosome
    uint64 size() const noexcept {
        return nucleos.size();
    }
    // For sorting from largest to smallest chromosome
    bool operator > (const RefChrom& other) const noexcept {
        return size() > other.size();
    }

    /*
     ------------------
     For filling a read at a given starting position from a chromosome of a
     given starting position and size.
     Used only for sequencer.
     ------------------
     */
    void fill_read(std::string& read,
                   const uint64& read_start,
                   const uint64& chrom_start,
                   uint64 n_to_add) const {
        // Making sure end doesn't go beyond the chromosome bounds
        if ((chrom_start + n_to_add - 1) >= nucleos.size()) {
            n_to_add = nucleos.size() - chrom_start;
        }
        // Make sure the read is long enough (this fxn should never shorten it):
        if (read.size() < n_to_add + read_start) read.resize(n_to_add + read_start, 'N');
        for (uint64 i = 0; i < n_to_add; i++) {
            read[(read_start + i)] = nucleos[(chrom_start + i)];
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
    std::deque<RefChrom> chromosomes;
    bool merged = false;
    // For storing original names if merged:
    std::deque<std::string> old_names = std::deque<std::string>(0);
    // Only added for compatibility with templates in sequencing code:
    std::string name = "REF";

    // Constructors
    RefGenome()
        : chromosomes(std::deque<RefChrom>(0)) {};
    RefGenome(const RefGenome& ref_)
        : total_size(ref_.total_size), chromosomes(ref_.chromosomes),
          merged(ref_.merged), old_names(ref_.old_names) {};
    RefGenome(const uint64& N)
        : chromosomes(std::deque<RefChrom>(N, RefChrom())) {};
    RefGenome(const std::deque<std::string>& chroms) {
        uint64 n_chroms = chroms.size();
        chromosomes = std::deque<RefChrom>(n_chroms, RefChrom());
        for (uint64 i = 0; i < n_chroms; i++) {
            chromosomes[i].nucleos = chroms[i];
            chromosomes[i].name = "chrom" + std::to_string(i);
            total_size += chroms[i].size();
        }
    }
    // Overloaded operator so chromosomes can be easily extracted
    // It returns a reference so no copying is done and so changes can be made
    RefChrom& operator[](const uint64& idx) {
#ifdef __JACKALOPE_DEBUG
        if (idx >= chromosomes.size()) {
            stop("Trying to extract chromosome that doesn't exist");
        }
#endif
        return chromosomes[idx];
    }
    const RefChrom& operator[](const uint64& idx) const {
#ifdef __JACKALOPE_DEBUG
        if (idx >= chromosomes.size()) {
            stop("Trying to extract chromosome that doesn't exist");
        }
#endif
        return chromosomes[idx];
    }
    // To return the number of chromosomes
    uint64 size() const noexcept {
        return chromosomes.size();
    }
    // To return the chromosome sizes
    std::vector<uint64> chrom_sizes() const {
        std::vector<uint64> out(size());
        for (uint64 i = 0; i < out.size(); i++) out[i] = chromosomes[i].size();
        return out;
    }
    // For printing reference genome info
    void print() const {

        int console_width = get_width();

        // 32 characters is the narrowest I'll allow
        // (I'd start getting negative lengths and other problems otherwise)
        console_width = std::max(console_width, 32);

        int num_chroms = size();
        std::vector<int> inds;
        if (num_chroms <= 10) {
            for (int i = 0; i < num_chroms; i++) inds.push_back(i);
        } else {
            for (int i = 0; i < 5; i++) inds.push_back(i);
            inds.push_back(-1);
            for (int i = (num_chroms - 5 + 1); i < num_chroms; i++) inds.push_back(i);
        }

        Rcout << "< Set of " << big_int_format<int>(num_chroms) << " chromosomes >";
        Rcout << std::endl;
        Rcout << "# Total size: " << big_int_format<uint64>(total_size) << " bp";
        Rcout << std::endl;

        int ind_i, name_width = 10, length_width = 9;
        // Console width minus name width AND length width AND spaces between
        int chrom_print_len = console_width - name_width - length_width - 2;
        // Number of chars print before and after elipses for a long string
        int before_elips = std::ceil((chrom_print_len - 3) / 2);
        int after_elips = chrom_print_len - 3 - before_elips;

        Rprintf("%-*s %s%-*s %*s\n", name_width, "  name",
                std::string(before_elips - 4, ' ').c_str(),
                chrom_print_len - (before_elips - 4), "chromosome",
                length_width, "length");

        for (int i = 0; i < static_cast<int>(inds.size()); i++) {
            ind_i = inds[i];
            if (ind_i == -1) {
                Rprintf("%-10s %-*s %9s\n", "...", chrom_print_len, "...", "...");
                continue;
            }
            const RefChrom& rs(chromosomes[ind_i]);
            const std::string& name_i(rs.name);
            const std::string& chrom_i(rs.nucleos);
            // Print name
            Rprintf("%-10.10s ", name_i.c_str());
            // Print chromosome
            int chrom_i_size = static_cast<int>(chrom_i.size());
            if (chrom_i_size > chrom_print_len){
                for (int j = 0; j < before_elips; j++) Rcout << chrom_i[j];
                Rcout << "...";
                for (int j = (chrom_i_size - after_elips); j < chrom_i_size; j++) {
                    Rcout << chrom_i[j];
                }
            } else {
                Rprintf("%-*s", chrom_print_len, chrom_i.c_str());
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
