#ifndef __JACKAL_SEQ_CLASSES_VAR_H
#define __JACKAL_SEQ_CLASSES_VAR_H


/*
 ********************************************************

 Classes to store variant sequence info.

 ********************************************************
 */



#include <RcppArmadillo.h>
#include <vector>  // vector class
#include <string>  // string class
#include <deque>  // deque class

#include "jackalope_types.h"  // integer types
#include "seq_classes_ref.h"  // Ref* classes
#include "util.h"  // clear_memory

using namespace Rcpp;




/*
 ========================================================================================
 ========================================================================================

 One mutation (substitution, insertion, or deletion)

 ========================================================================================
 ========================================================================================
 */

struct Mutation {

    // How this mutation changes the overall sequence size:
    sint64 size_modifier;
    // Position on the old (i.e., reference) sequence:
    uint64 old_pos;
    // Position on the new, variant sequence:
    uint64 new_pos;
    // Nucleotides associated with this mutation:
    std::string nucleos;

    // Constructors
    Mutation() {};
    Mutation(uint64 old_pos_, uint64 new_pos_, std::string nucleos_)
        : size_modifier(nucleos_.size() - 1), old_pos(old_pos_),
          new_pos(new_pos_), nucleos(nucleos_) {};
    Mutation(uint64 old_pos_, uint64 new_pos_, char nucleo_)
        : size_modifier(0), old_pos(old_pos_),
          new_pos(new_pos_), nucleos(1, nucleo_) {};
    Mutation(const Mutation& other)
        : size_modifier(other.size_modifier), old_pos(other.old_pos),
          new_pos(other.new_pos), nucleos(other.nucleos) {};
    // For deletions:
    Mutation(uint64 old_pos_, uint64 new_pos_, sint64 size_modifier_)
        : size_modifier(size_modifier_), old_pos(old_pos_),
          new_pos(new_pos_), nucleos("") {};

    /*
     These operators compare Mutation objects.
     If there is any overlap between the two mutations, then neither > nor <
     will return true.
     They are used to determine how and whether to merge VarSequence objects.
    */
    bool operator<(const Mutation& other) const {
        // Check for a deletion spanning the distance between the two:
        if (size_modifier < 0) {
            uint64 end_pos = old_pos + static_cast<uint64>(std::abs(size_modifier)) - 1;
            return end_pos < other.old_pos;
        }
        return old_pos < other.old_pos;
    }
    bool operator>(const Mutation& other) const {
        // Check for a deletion spanning the distance between the two:
        if (other.size_modifier < 0) {
            uint64 other_end_pos = other.old_pos +
                static_cast<uint64>(std::abs(other.size_modifier)) - 1;
            return old_pos > other_end_pos;
        }
        return old_pos > other.old_pos;
    }


    // For easily outputting mutation sequence
    const char& operator[](const uint64& idx) const {
        return nucleos[idx];
    }
};








/*
 ========================================================================================
 ========================================================================================

 Variant genomes

 ========================================================================================
 ========================================================================================
 */

// (This class will later need access to private members of VarSequence.)
class LocationSampler;


/*
 =========================================
 One sequence from one variant haploid genome
 =========================================
 */

class VarSequence {

    friend class LocationSampler;

public:

    const RefSequence* ref_seq;  // pointer to const RefSequence
    std::deque<Mutation> mutations;
    uint64 seq_size;
    std::string name;

    // Constructors
    VarSequence() {};
    VarSequence(const RefSequence& ref)
        : ref_seq(&ref),
          mutations(),
          seq_size(ref.size()),
          name(ref.name) {};

    /*
     Since all other classes have a size() method, I'm including this here:
     */
    uint64 size() const noexcept {
        return seq_size;
    }

    // Clear mutation info and restore RAM
    void clear() {
        mutations.clear();
        clear_memory<std::deque<Mutation>>(mutations);
        seq_size = ref_seq->size();
        return;
    }

    // Replace existing mutation information in this VarSequence with another
    void replace(const VarSequence& other) {
        ref_seq = other.ref_seq;
        mutations = other.mutations;
        seq_size = other.seq_size;
        name = other.name;
        return;
    }

    // Add existing mutation information in another `VarSequence` to this one
    VarSequence& operator+=(const VarSequence& other);


    /*
     ------------------
     Re-calculate new positions (and total sequence size)
     ------------------
     */
    void calc_positions();
    void calc_positions(uint64 mut_i);
    void calc_positions(uint64 mut_i, const sint64& modifier);



    /*
     ------------------
     Retrieve a nucleotide (char type) from the variant sequence
     based on the position in the new, variant sequence
     ------------------
     */
    char get_nt(const uint64& new_pos) const;

    /*
     ------------------
     Retrieve all nucleotides (i.e., the full sequence; std::string type) from
     the variant sequence
     ------------------
     */
    std::string get_seq_full() const;



    /*
     ------------------
     Set a string object to a chunk of a sequence from the variant sequence.
     ------------------
     */
    void set_seq_chunk(std::string& chunk_str,
                       const uint64& start,
                       const uint64& chunk_size,
                       uint64& mut_i) const;

    /*
     ------------------
     Adding mutations somewhere in the deque
     ------------------
     */
    void add_deletion(const uint64& size_, const uint64& new_pos_);
    void add_insertion(const std::string& nucleos_, const uint64& new_pos_);
    void add_substitution(const char& nucleo, const uint64& new_pos_);


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
                   uint64 n_to_add) const;



private:

    /*
     -------------------
     Internal function to "blowup" mutation(s) due to a deletion.
     By "blowup", I mean it removes substitutions and insertions if they're covered
     entirely by the deletion, and it merges any deletions that are contiguous.
     -------------------
     */
    void deletion_blowup_(uint64& mut_i, uint64& deletion_start, uint64& deletion_end,
                          sint64& size_mod);




    /*
     -------------------
     Inner function to merge an insertion and deletion.
     -------------------
     */
    void merge_del_ins_(uint64& insert_i,
                        uint64& deletion_start,
                        uint64& deletion_end,
                        sint64& size_mod);




    /*
     -------------------
     Inner function to remove Mutation and keep iterator from being invalidated.
     -------------------
     */
    void remove_mutation_(uint64& mut_i);
    void remove_mutation_(uint64& mut_i1, uint64& mut_i2);


    /*
     ------------------
     Internal function for finding character of either mutation or reference
     given an index (in the "new", variant sequence) and an iterator pointing to
     a single Mutation object.
     ------------------
     */
    char get_char_(const uint64& new_pos, const uint64& mut) const;

    /*
     ------------------
     Inner function to return an iterator to the Mutation object nearest to
     (without being past) an input position on the "new", variant sequence.
     ------------------
     */
    uint64 get_mut_(const uint64& new_pos) const;

};



/*
 =========================================
 One variant haploid genome
 =========================================
 */

class VarGenome {
public:

    // Fields
    std::string name;
    std::deque<VarSequence> var_genome;

    // Constructors
    VarGenome() {};
    VarGenome(const RefGenome& ref) {
        name = "";
        for (uint64 i = 0; i < ref.size(); i++) {
            VarSequence var_seq(ref[i]);
            var_genome.push_back(var_seq);
        }
    };
    VarGenome(const std::string& name_, const RefGenome& ref) {
        name = name_;
        for (uint64 i = 0; i < ref.size(); i++) {
            VarSequence var_seq(ref[i]);
            var_genome.push_back(var_seq);
        }
    };

    // For easily outputting a reference to a VarSequence
    VarSequence& operator[](const uint64& idx) {
        VarSequence& var_seq(var_genome[idx]);
        return var_seq;
    }
    // const version
    const VarSequence& operator[](const uint64& idx) const {
        const VarSequence& var_seq(var_genome[idx]);
        return var_seq;
    }
    // To return the number of sequences
    uint64 size() const noexcept {
        return var_genome.size();
    }
    // To return all the sequence sizes
    std::vector<uint64> seq_sizes() const noexcept {
        std::vector<uint64> out(size());
        for (uint64 i = 0; i < out.size(); i++) out[i] = var_genome[i].size();
        return out;
    }

private:

};





/*
 =========================================
 Multiple variant haploid genomes (based on the same reference)
 =========================================
 */

class VarSet {
public:
    std::deque<VarGenome> variants;
    const RefGenome* reference;  // pointer to const RefGenome

    /*
     Constructors:
     */
    VarSet(const RefGenome& ref) : variants(), reference(&ref) {};
    VarSet(const RefGenome& ref, const uint64& n_vars)
        : variants(n_vars, VarGenome(ref)),
          reference(&ref) {
        for (uint64 i = 0; i < n_vars; i++) variants[i].name = "var" + std::to_string(i);
    };
    // If you already have the names:
    VarSet(const RefGenome& ref, const std::vector<std::string>& names_)
        : variants(names_.size(), VarGenome(ref)),
          reference(&ref) {
        for (uint64 i = 0; i < names_.size(); i++) variants[i].name = names_[i];
    };

    // For easily outputting a reference to a VarGenome
    VarGenome& operator[](const uint64& idx) {
        if (idx >= variants.size()) {
            stop("trying to access a VarGenome that doesn't exist");
        }
        VarGenome& vg(variants[idx]);
        return vg;
    }
    // const version of above
    const VarGenome& operator[](const uint64& idx) const {
        if (idx >= variants.size()) {
            stop("trying to access a VarGenome that doesn't exist");
        }
        const VarGenome& vg(variants[idx]);
        return vg;
    }
    // To return the number of variants
    uint64 size() const noexcept {
        return variants.size();
    }
    // To return the minimum size of a given sequence
    uint64 min_size(const uint64& i) const {
        uint64 ms = variants[0][i].size();
        for (const VarGenome vg : variants) {
            if (vg[i].size() < ms) ms = vg[i].size();
        }
        return ms;
    }

    /*
     Fill VarGenome objects after the reference has been filled
     */
    void fill_vars(const uint64& n_vars) {
        VarGenome vg(*reference);
        for (uint64 i = 0; i < n_vars; i++) variants.push_back(vg);
        return;
    }
    // Overloaded for if you want to provide names
    void fill_vars(const std::vector<std::string>& names) {
        for (uint64 i = 0; i < names.size(); i++) {
            VarGenome vg(names[i], *reference);
            variants.push_back(vg);
        }
        return;
    }

    // For printing output
    void print() const noexcept;

};

#endif
