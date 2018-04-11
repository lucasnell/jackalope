#ifndef _VARIANT_CLASS_H
#define _VARIANT_CLASS_H


/*
 * THIS FILE SEEKS TO REPLACE THE OLD CLASSES
 */



#include <RcppArmadillo.h>
#include <vector>  // vector class
#include <string>  // string class
#include <algorithm>  // lower_bound, sort
#include <deque>  // deque

#include "gemino_types.h"  // integer types

using namespace Rcpp;



/*
 ========================================================================
 ========================================================================

 One reference-genome sequence (e.g., chromosome, scaffold)

 ========================================================================
 ========================================================================
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
    char operator[](const uint& idx) const {
        if (idx >= nucleos.size()) {
            stop("Trying to extract nucleotide that doesn't exist");
        }
        return nucleos[idx];
    }
    // To return the size of this sequence
    uint size() const noexcept {
        return nucleos.size();
    }
    // For sorting from largest to smallest sequence
    bool operator > (const RefSequence& other) const noexcept {
        return size() > other.size();
    }
};





/*
 ========================================================================
 ========================================================================

 One reference genome

 ========================================================================
 ========================================================================
 */

struct RefGenome {

    // Member variables
    uint64 total_size = 0;
    std::deque<RefSequence> sequences;
    bool merged = false;
    // For storing original names if merged:
    std::deque<std::string> old_names;

    // Constructors
    RefGenome()
        : sequences(std::deque<RefSequence>(0)) {};
    RefGenome(const RefGenome& ref_)
        : total_size(ref_.total_size), sequences(ref_.sequences),
          merged(ref_.merged), old_names(ref_.old_names) {};
    RefGenome(const uint& N)
        : sequences(std::deque<RefSequence>(N, RefSequence())) {};
    RefGenome(const std::deque<std::string>& seqs) {
        uint n_seqs = seqs.size();
        sequences = std::deque<RefSequence>(n_seqs, RefSequence());
        for (uint i = 0; i < n_seqs; i++) {
            sequences[i].nucleos = seqs[i];
            sequences[i].name = "seq" + std::to_string(i);
            total_size += seqs[i].size();
        }
    }
    // Overloaded operator so sequences can be easily extracted
    // It returns a const reference so no copying is done and no changes can be made
    const RefSequence& operator[](const uint& idx) const {
        if (idx >= sequences.size()) {
            stop("Trying to extract sequence that doesn't exist");
        }
        return sequences[idx];
    }
    // To return the number of sequences
    uint size() const noexcept {
        return sequences.size();
    }
    // For printing reference genome info
    void print() const;
};






/*
 ========================================================================
 ========================================================================

 One mutation (substitution, insertion, or deletion)

 ========================================================================
 ========================================================================
 */

struct Mutation {

    // How this mutation changes the overall sequence size:
    sint size_modifier;
    // Position on the old (i.e., reference) sequence:
    uint old_pos;
    // Position on the new, variant sequence:
    uint new_pos;
    // Nucleotides associated with this mutation:
    std::string nucleos;

    // Constructors
    Mutation() {};
    Mutation(uint old_pos_, uint new_pos_, std::string nucleos_)
        : size_modifier(nucleos_.size() - 1), old_pos(old_pos_),
          new_pos(new_pos_), nucleos(nucleos_) {};
    // For deletions:
    Mutation(uint old_pos_, uint new_pos_, sint size_modifier_)
        : size_modifier(size_modifier_), old_pos(old_pos_),
          new_pos(new_pos_), nucleos("") {};

    /*
     Is this Mutation "less than" another Mutation object?
     This is used to `std::sort` Mutation objects as well as to compare them like you
     would an integer.
     I'm using the old_pos field bc it should never change based on new mutations.
     */
    bool operator<(const Mutation& other) const {
        return old_pos < other.old_pos;
    }
    // Same for greater than and equal to
    bool operator>(const Mutation& other) const {
        return old_pos > other.old_pos;
    }
    bool operator==(const Mutation& other) const {
        return old_pos == other.old_pos;
    }
    // For easily outputting mutation sequence
    const char& operator[](const uint& idx) const {
        return nucleos[idx];
    }
};



/*
 ========================================================================
 ========================================================================

 One sequence from one variant haploid genome

 ========================================================================
 ========================================================================
 */

class VarSequence {
public:

    const RefSequence& ref_seq;
    mutable std::deque<Mutation> mutations;
    uint seq_size;

    // Constructor
    VarSequence(const RefSequence& ref)
        : ref_seq(ref), mutations(std::deque<Mutation>()),
          seq_size(ref.size()) {};

    /*
     Since all other classes have a size() method, I'm including this here:
     */
    uint size() const noexcept {
        return seq_size;
    }


    /*
     ------------------
     Re-calculate new positions (and total sequence size)
     ------------------
     */
    void calc_positions(bool sort_first = true);
    void calc_positions(uint mut_i);
    void calc_positions(uint mut_i, const sint& modifier);



    /*
     ------------------
     Retrieve a nucleotide (char type) from the variant sequence
     based on the position in the new, variant sequence
     ------------------
     */
    char get_nt(const uint& new_pos) const;

    /*
     ------------------
     Retrieve all nucleotides (i.e., the full sequence; std::string type) from
     the variant sequence
     ------------------
     */
    std::string get_seq_full() const;


    /*
     ------------------
     Retrieve the first part of a sequence from the variant sequence.
     ------------------
     */
    std::string get_seq_start(uint out_length) const;

    /*
     ------------------
     Set a string object to a chunk of a sequence from the variant sequence.
     ------------------
     */
    void set_seq_chunk(std::string& chunk_str,
                       const uint& start,
                       const uint& chunk_size,
                       uint& mut_i) const;

    /*
     ------------------
     Adding mutations somewhere in the deque
     ------------------
     */
    void add_deletion(const uint& size_, const uint& new_pos_);
    void add_insertion(const std::string& nucleos_, const uint& new_pos_);
    void add_substitution(const char& nucleo, const uint& new_pos_);



private:

    /*
     -------------------
     Internal function to "blowup" mutation(s) due to a deletion.
     By "blowup", I mean it removes substitutions and insertions if they're covered
     entirely by the deletion, and it merges any deletions that are contiguous.
     -------------------
     */
    void deletion_blowup_(
            std::deque<Mutation>::iterator& iter,
            uint& deletion_start, uint& deletion_end, sint& size_mod);




    /*
     -------------------
     Inner function to merge an insertion and deletion.
     -------------------
     */
    void merge_del_ins_(std::deque<Mutation>::iterator& insertion,
                        uint& deletion_start,
                        uint& deletion_end,
                        sint& size_mod,
                        sint& iterd);




    /*
     -------------------
     Inner function to remove Mutation and keep iterator from being invalidated.
     -------------------
     */
    void remove_mutation_(std::deque<Mutation>::iterator& mutation);
    void remove_mutation_(std::deque<Mutation>::iterator& mutation1,
                          std::deque<Mutation>::iterator& mutation2);


    /*
     ------------------
     Internal function for finding character of either mutation or reference
     given an index (in the "new", variant sequence) and an iterator pointing to
     a single Mutation object.
     ------------------
     */
    char get_char_(const uint& new_pos, const uint& mut) const;

    /*
     ------------------
     Inner function to return an iterator to the Mutation object nearest to
     (without being past) an input position on the "new", variant sequence.
     ------------------
     */
    uint get_mut_(const uint& new_pos) const;

};



/*
 ========================================================================
 ========================================================================

 One variant haploid genome

 ========================================================================
 ========================================================================
 */

class VarGenome {
public:

    // Fields
    std::string name;
    std::deque<VarSequence> var_genome;

    // Constructors
    VarGenome(const RefGenome& ref) {
        name = "";
        for (uint i = 0; i < ref.size(); i++) {
            VarSequence vs(ref[i]);
            var_genome.push_back(vs);
        }
    };
    VarGenome(const std::string& name_, const RefGenome& ref) {
        name = name_;
        for (uint i = 0; i < ref.size(); i++) {
            VarSequence vs(ref[i]);
            var_genome.push_back(vs);
        }
    };

    // For easily outputting a reference to a VarSequence
    VarSequence& operator[](const uint& idx) {
        VarSequence& vs(var_genome[idx]);
        return vs;
    }
    // To return the number of sequences
    uint size() const noexcept {
        return var_genome.size();
    }

private:

};






/*
 ========================================================================
 ========================================================================

 Multiple variant haploid genomes (based on the same reference)

 ========================================================================
 ========================================================================
 */


class VarSet {
public:
    std::deque<VarGenome> variants;
    RefGenome reference;

    /*
     Constructors:
     */
    VarSet() {};
    VarSet(const RefGenome& ref, const uint& n_vars)
        : variants(std::deque<VarGenome>(n_vars, VarGenome(ref))),
          reference(ref) {
        for (uint i = 0; i < n_vars; i++) variants[i].name = "var" + std::to_string(i);
    };
    VarSet(const std::string& fasta_file, const uint& n_vars,
           const bool& cut_names = true, const bool& remove_soft_mask = true);
    VarSet(const std::string& fasta_file, const std::string& fai_file,
           const uint& n_vars,
           const bool& remove_soft_mask = true);
    VarSet(const std::deque<std::string>& seqs, const uint& n_vars);

    // For easily outputting a reference to a VarGenome
    VarGenome& operator[](const uint& idx) {
        if (idx >= variants.size()) {
            stop("trying to access a VarGenome that doesn't exist");
        }
        VarGenome& vg(variants[idx]);
        return vg;
    }
    // To return the number of variants
    uint size() const noexcept {
        return variants.size();
    }

    // For printing output
    void print() const noexcept;

    /*
     Fill VarGenome objects after the reference has been filled
     */
    void fill_vars(const uint& n_vars) {
        VarGenome vg(reference);
        for (uint i = 0; i < n_vars; i++) variants.push_back(vg);
        return;
    }
    // Overloaded for if you want to provide names
    void fill_vars(const std::vector<std::string>& names) {
        for (uint i = 0; i < names.size(); i++) {
            VarGenome vg(names[i], reference);
            variants.push_back(vg);
        }
        return;
    }
};

#endif
