#ifndef __GEMINO_DIGEST_H
#define __GEMINO_DIGEST_H

#include <RcppArmadillo.h>
#include <string> // find, string::npos, string::size_type, ...
#include <vector> // vector class
#include <deque> // deque class

#include "sequence_classes.h" // classes RefGenome, VarSet, VarSequence



/*
 Organizes info needed for digestion
 */
struct DigestInfo {
    std::vector<std::string> bind_sites;
    std::vector<uint> len5s;
    uint max_size;

    // Constructors
    DigestInfo() : bind_sites(), len5s(), max_size() {};

    DigestInfo(const std::vector<std::string>& bind_sites_,
               const std::vector<uint>& len5s_)
        : bind_sites(bind_sites_), len5s(len5s_), max_size(bind_sites_[0].size()) {

        for (uint i = 1; i < bind_sites_.size(); i++) {
            if (bind_sites_[i].size() > max_size) max_size = bind_sites_[i].size();
        }

    };

    // Get the # bp to skip before looking for another match
    uint to_skip(const uint& j) const noexcept {
        return bind_sites[j].size();
    }
};





/*
 When searching for multiple strings simultaneously, this is a way to return the
 first position where one of them match and the string match's index
 */
struct MultiOut {
    // fields
    std::string::size_type new_pos;
    uint j;
    // Constructor
    MultiOut() {
        new_pos = std::string::npos;
        j = 1000000;
    }
    // reset to default values
    void reset() {
        new_pos = std::string::npos;
        j = 1000000;
    }
};



// Info on the digestion of a sequence's chunk

class RefSeqChunk {

public:

    const uint seq;    // index of the sequence this refers to
    const uint start;  // index of the starting position in the sequence of this chunk
    const uint end;    // index of the ending position in the sequence of this chunk
    std::deque<uint> cut_sites;  // cut sites found
    /*
     Earliest possible start for the next chunk in this sequence. This is used
     when merging the "seam" between this chunk and the one after it.
     */
    uint next_start;

    RefSeqChunk(const std::string& ref_,
                const uint& seq_,
                const uint& start_,
                const uint& chunk_size)
        : seq(seq_), start(start_),
          end(std::min(static_cast<uint>(start_ + chunk_size - 1),
                       static_cast<uint>(ref_.size() - 1))),
          cut_sites() ,
          ref(ref_) {};

    // To return the number of cut sites
    uint size() const noexcept {
        return cut_sites.size();
    }
    // For combining with an output deque
    void combine(std::deque<uint>& out_dq) {
        for (uint j = 0; j < cut_sites.size(); j++) {
            out_dq.push_back(cut_sites[j]);
        }
    }

    // Digest this chunk
    void digest(const DigestInfo& dinfo);
    // Make sure no positions overlap in the "seam" between two chunks.
    void merge_seam(const RefSeqChunk& prev, const DigestInfo& dinfo);

private:
    const std::string& ref;  // reference genome string this refers to
};




// Info on the digestion of a variant's sequence
class VarSeqDigest {

public:

    const VarSequence& var_info;  // info on this variant's sequence
    const uint var;    // index of the variant this refers to
    const uint seq;    // index of the sequence this refers to
    const uint chunk_size;    // size of chunks to digest sequence by
    std::deque<uint>& cut_sites;  // cut sites found

    VarSeqDigest(VarSet& var_info_,
                 const uint& var_,
                 const uint& seq_,
                 const uint& chunk_size_,
                 std::deque<uint>& cut_sites_)
        : var_info(var_info_[var_][seq_]), var(var_), seq(seq_),
          chunk_size(chunk_size_), cut_sites(cut_sites_), start(0) {};

    // Digest this sequence
    void digest(const DigestInfo& dinfo);

private:
    uint start;  // starting positions of chunks used as it's digested
};




#endif
