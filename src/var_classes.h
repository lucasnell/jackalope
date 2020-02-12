#ifndef __JACKALOPE_VAR_CLASSES_H
#define __JACKALOPE_VAR_CLASSES_H


/*
 ********************************************************

 Classes to store haplotype chromosome info.

 ********************************************************
 */


#include "jackalope_config.h" // controls debugging and diagnostics output

#include <RcppArmadillo.h>
#include <vector>  // vector class
#include <string>  // string class
#include <cstring> // for std::strcpy
#include <deque>  // deque class

#include "jackalope_types.h"  // integer types
#include "ref_classes.h"  // Ref* classes
#include "util.h"  // clear_memory

using namespace Rcpp;




/*
 ========================================================================================
 ========================================================================================

 One mutation (substitution, insertion, or deletion)

 ========================================================================================
 ========================================================================================
 */

// struct Mutation {
//
//     // How this mutation changes the overall chromosome size:
//     sint64 size_modifier;
//     // Position on the old (i.e., reference) chromosome:
//     uint64 old_pos;
//     // Position on the new, haplotype chromosome:
//     uint64 new_pos;
//     // Nucleotides associated with this mutation:
//     std::string nucleos;
//
//     // Constructors
//     Mutation(uint64 old_pos_, uint64 new_pos_, std::string nucleos_)
//         : size_modifier(nucleos_.size() - 1), old_pos(old_pos_),
//           new_pos(new_pos_), nucleos(nucleos_) {};
//     Mutation(uint64 old_pos_, uint64 new_pos_, char nucleo_)
//         : size_modifier(0), old_pos(old_pos_),
//           new_pos(new_pos_), nucleos(1, nucleo_) {};
//     Mutation(const Mutation& other)
//         : size_modifier(other.size_modifier), old_pos(other.old_pos),
//           new_pos(other.new_pos), nucleos(other.nucleos) {};
//     // For deletions:
//     Mutation(uint64 old_pos_, uint64 new_pos_, sint64 size_modifier_)
//         : size_modifier(size_modifier_), old_pos(old_pos_),
//           new_pos(new_pos_), nucleos("") {};
//
//     /*
//      These operators compare Mutation objects.
//      If there is any overlap between the two mutations, then neither > nor <
//      will return true.
//      They are used to determine how and whether to merge HapChrom objects.
//     */
//     bool operator<(const Mutation& other) const {
//         // Check for a deletion spanning the distance between the two:
//         if (size_modifier < 0) {
//             uint64 end_pos = old_pos + static_cast<uint64>(std::abs(size_modifier)) - 1;
//             return end_pos < other.old_pos;
//         }
//         return old_pos < other.old_pos;
//     }
//     bool operator>(const Mutation& other) const {
//         // Check for a deletion spanning the distance between the two:
//         if (other.size_modifier < 0) {
//             uint64 other_end_pos = other.old_pos +
//                 static_cast<uint64>(std::abs(other.size_modifier)) - 1;
//             return old_pos > other_end_pos;
//         }
//         return old_pos > other.old_pos;
//     }
//
//
//     // For easily outputting mutation chromosome
//     const char& operator[](const uint64& idx) const {
//         return nucleos[idx];
//     }
// };



struct AllMutations {

    std::deque<uint64> old_pos;
    std::deque<uint64> new_pos;
    std::deque<char*> nucleos;

    AllMutations() : old_pos(), new_pos(), nucleos() {}

    AllMutations(const AllMutations& other)
        : old_pos(other.old_pos),
          new_pos(other.new_pos),
          nucleos(other.nucleos.size(), nullptr) {
        for (uint64 i = 0; i < nucleos.size(); i++) {
            fill_one_nucleos__(other.nucleos[i], &nucleos[i]);
        }
    }
    AllMutations& operator=(const AllMutations& other) {

        old_pos = other.old_pos;
        new_pos = other.new_pos;
        for (uint64 i = 0; i < nucleos.size(); i++) delete [] nucleos[i];
        nucleos = std::deque<char*>(other.nucleos.size(), nullptr);
        for (uint64 i = 0; i < nucleos.size(); i++) {
            fill_one_nucleos__(other.nucleos[i], &nucleos[i]);
        }
        return *this;
    }

    ~AllMutations() {
        for (uint64 i = 0; i < nucleos.size(); i++) {
            delete [] nucleos[i];
        }
    }


    inline size_t size() const noexcept {
        return old_pos.size();
    }

    inline bool empty() const noexcept {
        return old_pos.empty();
    }

    inline void clear() {

        if (old_pos.size() == 0) return;

        old_pos.clear();
        new_pos.clear();
        nucleos.clear();

        return;
    }
    // Add to front
    inline void push_front(const uint64& op,
                           const uint64& np,
                           const char* nts) {
        old_pos.push_front(op);
        new_pos.push_front(np);
        nucleos.push_front(nullptr);
        fill_one_nucleos__(nts, &nucleos.front());
        return;
    }
    inline void push_front(const uint64& op,
                           const uint64& np,
                           const char& nt) {
        old_pos.push_front(op);
        new_pos.push_front(np);
        nucleos.push_front(nullptr);
        fill_one_nucleos__(nt, &nucleos.front());
        return;
    }
    // Add to back
    inline void push_back(const uint64& op,
                          const uint64& np,
                          const char* nts) {
        old_pos.push_back(op);
        new_pos.push_back(np);
        nucleos.push_back(nullptr);
        fill_one_nucleos__(nts, &nucleos.back());
        return;
    }
    inline void push_back(const uint64& op,
                          const uint64& np,
                          const char& nt) {
        old_pos.push_back(op);
        new_pos.push_back(np);
        nucleos.push_back(nullptr);
        fill_one_nucleos__(nt, &nucleos.back());
        return;
    }
    // Add to middle
    inline void insert(const uint64& ind,
                       const uint64& op,
                       const uint64& np,
                       const char* nts) {

        old_pos.insert(old_pos.begin() + ind, op);
        new_pos.insert(new_pos.begin() + ind, np);
        nucleos.insert(nucleos.begin() + ind, nullptr);
        fill_one_nucleos__(nts, &nucleos[ind]);
        return;
    }
    inline void insert(const uint64& ind,
                       const uint64& op,
                       const uint64& np,
                       const char& nt) {

        old_pos.insert(old_pos.begin() + ind, op);
        new_pos.insert(new_pos.begin() + ind, np);
        nucleos.insert(nucleos.begin() + ind, nullptr);
        fill_one_nucleos__(nt, &nucleos[ind]);
        return;
    }
    // Remove from position
    inline void erase(const uint64& ind) {

        old_pos.erase(old_pos.begin() + ind);
        new_pos.erase(new_pos.begin() + ind);
        delete [] nucleos[ind];
        nucleos.erase(nucleos.begin() + ind);
        return;
    }

    // Remove between positions
    inline void erase(const uint64& ind1, const uint64& ind2) {

        old_pos.erase(old_pos.begin() + ind1, old_pos.begin() + ind2);
        new_pos.erase(new_pos.begin() + ind1, new_pos.begin() + ind2);
        for (uint64 ind = ind1; ind < ind2; ind++) delete [] nucleos[ind];
        nucleos.erase(nucleos.begin() + ind1, nucleos.begin() + ind2);
        return;
    }


private:


    // Fill one `char *` in `nucleos` with values from a `const char *` object
    inline void fill_one_nucleos__(const char* nts, char** one_nucleos) {
        if (nts != nullptr) {
            size_t nts_size = std::strlen(nts);
            *one_nucleos = new char[nts_size + 1];
            for (size_t i = 0; i < nts_size; i++) {
                (*one_nucleos)[i] = nts[i];
            }
            (*one_nucleos)[nts_size] = '\0';
        }
        return;
    }
    // Same but for a reference to one character
    inline void fill_one_nucleos__(const char& nt, char** one_nucleos) {
        *one_nucleos = new char[2];
        (*one_nucleos)[0] = nt;
        (*one_nucleos)[1] = '\0';
        return;
    }

};








/*
 ========================================================================================
 ========================================================================================

 Single haplotypes

 ========================================================================================
 ========================================================================================
 */

// (These classes will later need access to private members of HapChrom.)
class SubMutator;


/*
 =========================================
 One chromosome from one haplotype haploid genome
 =========================================
 */

class HapChrom {

    friend SubMutator;

public:

    const RefChrom* ref_chrom;  // pointer to const RefChrom
    AllMutations mutations;
    uint64 chrom_size;
    std::string name;

    // Constructors
    HapChrom() : ref_chrom(nullptr) {};
    HapChrom(const RefChrom& ref)
        : ref_chrom(&ref),
          mutations(),
          chrom_size(ref.size()),
          name(ref.name) {};

    /*
     Since all other classes have a size() method, I'm including this here:
     */
    uint64 size() const noexcept {
        return chrom_size;
    }

    // Size modifier for a mutation
    sint64 size_modifier(const uint64& ind) const {

        sint64 size_mod;

        if (ind < (mutations.new_pos.size() - 1)) {
            size_mod = mutations.new_pos[ind+1] - mutations.old_pos[ind+1];
#ifdef __JACKALOPE_DEBUG
        } else if (ind == (mutations.size() - 1)) {
            size_mod = chrom_size - ref_chrom->size();
        } else stop("ind too large inside size_modifier function.");
#else
        } else {
            size_mod = chrom_size - ref_chrom->size();
        }
#endif

        size_mod += static_cast<sint64>(mutations.old_pos[ind] - mutations.new_pos[ind]);

        return size_mod;
}

    // Add existing mutation information in another `HapChrom` to this one,
    // adding to the back of `mutations`, with a starting mutation index
    sint64 add_to_back(const HapChrom& other, const uint64& mut_i);


    /*
     ------------------
     Re-calculate new positions (and total chromosome size)
     ------------------
     It re-calculates positions for all mutations after AND INCLUDING the given index
     */
    inline void calc_positions(uint64 mut_i, const sint64& modifier) {
        // Updating individual Mutation objects
        for (; mut_i < mutations.size(); ++mut_i) {
            mutations.new_pos[mut_i] += modifier;
        }
        // Updating full chromosome size
        chrom_size += modifier;

        return;
    }


    /*
     ------------------
     Retrieve all nucleotides (i.e., the full chromosome; std::string type) from
     the haplotype chromosome
     ------------------
     */
    std::string get_chrom_full() const;



    /*
     ------------------
     Set a string object to a chunk of a chromosome from the haplotype chromosome.
     ------------------
     */
    void set_chrom_chunk(std::string& chunk_str,
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
     For filling a read at a given starting position from a chromosome of a
     given starting position and size.
     Used only for sequencer.
     ------------------
     */
    void fill_read(std::string& read,
                   const uint64& read_start,
                   const uint64& chrom_start,
                   uint64 n_to_add) const;



private:


    /*
     -------------------
     Inner function to get old position for deletion.
     -------------------
     */
    uint64 deletion_old_pos_(const uint64& deletion_start,
                             const uint64& deletion_end,
                             const uint64& mut_i) const;



    /*
     -------------------
     Inner function to adjust a single mutation for a deletion.
     -------------------
     */
    void deletion_one_mut_(const uint64& mut_i,
                           const uint64& deletion_start,
                           const uint64& deletion_end,
                           const sint64& full_size_mod,
                           sint64& new_size_mod,
                           std::vector<uint64>& rm_inds);



    /*
     ------------------
     Internal function for finding character of either mutation or reference
     given an index (in the "new", haplotype chromosome) and an index for a
     single Mutation object.
     This only works if you've already narrowed it down to the Mutation object
     that is directly previous to the index position.
     ------------------
     */
    inline char get_char_(const uint64& new_pos,
                          const uint64& mut_i) const {
        char out;
        uint64 ind = new_pos - mutations.new_pos[mut_i];
        if (static_cast<sint64>(ind) > size_modifier(mut_i)) {
            ind += (mutations.old_pos[mut_i] - size_modifier(mut_i));
            out = (*ref_chrom)[ind];
        } else {
            if (mutations.nucleos[mut_i] == nullptr) {
                std::string err_msg = "mutations.nucleos[mut_i] == nullptr at ";
                err_msg += std::to_string(mut_i);
                stop(err_msg.c_str());
            }
            out = mutations.nucleos[mut_i][ind];
        }
        return out;
    }




    /*
     ------------------
     Inner function to return an iterator to the Mutation object nearest to
     (without being past) an input position on the "new", haplotype chromosome.
     ------------------
     */
    uint64 get_mut_(const uint64& new_pos) const;

};



/*
 =========================================
 One haplotype haploid genome
 =========================================
 */

class HapGenome {
public:

    // Fields
    std::string name;
    std::vector<HapChrom> chromosomes;

    // Constructors
    HapGenome() {};
    HapGenome(const RefGenome& ref) {
        name = "";
        for (uint64 i = 0; i < ref.size(); i++) {
            HapChrom hap_chrom(ref[i]);
            chromosomes.push_back(hap_chrom);
        }
    };
    HapGenome(const std::string& name_, const RefGenome& ref) {
        name = name_;
        for (uint64 i = 0; i < ref.size(); i++) {
            HapChrom hap_chrom(ref[i]);
            chromosomes.push_back(hap_chrom);
        }
    };

    // For easily outputting a reference to a HapChrom
    HapChrom& operator[](const uint64& idx) {
        HapChrom& hap_chrom(chromosomes[idx]);
        return hap_chrom;
    }
    // const version
    const HapChrom& operator[](const uint64& idx) const {
        const HapChrom& hap_chrom(chromosomes[idx]);
        return hap_chrom;
    }
    // To return the number of chromosomes
    uint64 size() const noexcept {
        return chromosomes.size();
    }
    // To return all the chromosome sizes
    std::vector<uint64> chrom_sizes() const noexcept {
        std::vector<uint64> out(size());
        for (uint64 i = 0; i < out.size(); i++) out[i] = chromosomes[i].size();
        return out;
    }

private:

};





/*
 =========================================
 Multiple haplotype haploid genomes (based on the same reference)
 =========================================
 */

class HapSet {
public:
    std::vector<HapGenome> haplotypes;
    const RefGenome* reference;  // pointer to const RefGenome

    /*
     Constructors:
     */
    HapSet(const RefGenome& ref) : haplotypes(), reference(&ref) {};
    HapSet(const RefGenome& ref, const uint64& n_haps)
        : haplotypes(n_haps, HapGenome(ref)),
          reference(&ref) {
        for (uint64 i = 0; i < n_haps; i++) haplotypes[i].name = "hap" + std::to_string(i);
    };
    // If you already have the names:
    HapSet(const RefGenome& ref, const std::vector<std::string>& names_)
        : haplotypes(names_.size(), HapGenome(ref)),
          reference(&ref) {
        for (uint64 i = 0; i < names_.size(); i++) haplotypes[i].name = names_[i];
    };

    // For easily outputting a reference to a HapGenome
    HapGenome& operator[](const uint64& idx) {
#ifdef __JACKALOPE_DEBUG
        if (idx >= haplotypes.size()) {
            stop("trying to access a HapGenome that doesn't exist");
        }
#endif
        HapGenome& vg(haplotypes[idx]);
        return vg;
    }
    // const version of above
    const HapGenome& operator[](const uint64& idx) const {
#ifdef __JACKALOPE_DEBUG
        if (idx >= haplotypes.size()) {
            stop("trying to access a HapGenome that doesn't exist");
        }
#endif
        const HapGenome& vg(haplotypes[idx]);
        return vg;
    }
    // To return the number of haplotypes
    uint64 size() const noexcept {
        return haplotypes.size();
    }
    // To return the minimum size of a given chromosome
    uint64 min_size(const uint64& i) const {
        uint64 ms = haplotypes[0][i].size();
        for (const HapGenome vg : haplotypes) {
            if (vg[i].size() < ms) ms = vg[i].size();
        }
        return ms;
    }

    /*
     Fill HapGenome objects after the reference has been filled
     */
    void fill_haps(const uint64& n_haps) {
        HapGenome vg(*reference);
        for (uint64 i = 0; i < n_haps; i++) haplotypes.push_back(vg);
        return;
    }
    // Overloaded for if you want to provide names
    void fill_haps(const std::vector<std::string>& names) {
        for (uint64 i = 0; i < names.size(); i++) {
            HapGenome vg(names[i], *reference);
            haplotypes.push_back(vg);
        }
        return;
    }

    // For printing output
    void print() const noexcept;

};

#endif
