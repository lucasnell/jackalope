//
// This file digests sequences from variant set and sequence set objects
//


#include <RcppArmadillo.h>
#include <string> // find, string::npos, string::size_type, ...
#include <vector> // vector class
#include <deque> // deque class
#ifdef _OPENMP
#include <omp.h>  // omp
#endif


#include "gemino_types.h" // integer types
#include "sequence_classes.h" // classes RefGenome and VarSet
#include "digest.h" // DigestInfo, MultiOut,

using namespace Rcpp;





/*
 Search for multiple substrings simulateously, change the `search_obj` object so that it
 stores...
 (1) the first position in `s` that matches one of the strings in `bind_sites` and
 (2) the index in `bind_sites` for the site that matches to a position in `s`
 */
void multi_search(
        MultiOut& search_obj,
        const DigestInfo& dinfo,
        const std::string& s,
        const uint& last_pos) {

    search_obj.reset();

    std::string::size_type new_pos;

    for (uint j = 0; j < dinfo.bind_sites.size(); j++) {
        new_pos = s.find(dinfo.bind_sites[j], last_pos);
        if (new_pos < search_obj.new_pos) {
            search_obj.new_pos = new_pos;
            search_obj.j = j;
        }
    }

    return;
}


// ------------------
// Digest one string
// ------------------

/*
 (`pos_mod` is to modify positions because `s` is a chunk not from the start of the
  sequence.)
 */
void digest_seq(std::deque<uint>& pos_vec,
                uint& start_pos,
                const std::string& s,
                const DigestInfo& dinfo,
                const uint& pos_mod = 0,
                const uint& end = 0) {

    uint cut_site_pos;
    MultiOut search_obj;

    multi_search(search_obj, dinfo, s, start_pos);

    if (end == 0) {
        while (search_obj.new_pos != std::string::npos) {
            cut_site_pos = search_obj.new_pos + dinfo.len5s[search_obj.j] + pos_mod;
            pos_vec.push_back(cut_site_pos);
            start_pos = search_obj.new_pos + dinfo.to_skip(search_obj.j);
            multi_search(search_obj, dinfo, s, start_pos);
        }
    } else {
        while (search_obj.new_pos != std::string::npos && search_obj.new_pos <= end) {
            cut_site_pos = search_obj.new_pos + dinfo.len5s[search_obj.j] + pos_mod;
            pos_vec.push_back(cut_site_pos);
            start_pos = search_obj.new_pos + dinfo.to_skip(search_obj.j);
            multi_search(search_obj, dinfo, s, start_pos);
        }
    }

    return;
}




/*
 =============================================================================
 =============================================================================
 =============================================================================
 =============================================================================

 Variant digestion

 =============================================================================
 =============================================================================
 =============================================================================
 =============================================================================
 */




// Digest one variant genome's sequence

void VarSeqDigest::digest(const DigestInfo& dinfo) {

    // Reset these if it's been digested before
    start = 0;
    if (!cut_sites.empty()) cut_sites.clear();

    if (chunk_size >= var_info.seq_size) {
        std::string seq = var_info.get_seq_full();
        digest_seq(cut_sites, start, seq, dinfo);
        return;
    }

    // index to the nearest mutation for a given chunk's starting position:
    uint mut_i = 0;

    // Initialize chunk sequence string
    std::string chunk_str(chunk_size, 'x');

    /*
     I'm iterating by `chunk_size - (max_size - 1)` to make sure I catch any
    matches that overlap the "seam" between chunks.
    */
    uint it = chunk_size - (dinfo.max_size - 1);

    /*
     The last chunk start position in this sequence (without allowing only overlap
     at the end)
     */
    uint last_chunk = var_info.seq_size - dinfo.max_size;

    // Digesting each chunk:
    for (uint chunk_start = 0; chunk_start <= last_chunk; chunk_start += it) {
        // Filling new chunk string:
        var_info.set_seq_chunk(chunk_str, chunk_start, chunk_size, mut_i);
        // Compensating for matches at the end of the previous chunk:
        if (start > it) {
            start -= it;
        } else start = 0;
        // Digest chunk:
        digest_seq(cut_sites, start, chunk_str, dinfo, chunk_start);
    }

    return;
}




// ------------------
// Digest all sequences, all variants
// ------------------



//' Internal C++ function to digest all sequences for all variants in a variant set.
//'
//'
//'
//' @param var_ An external pointer to a C++ \code{VarSet} object
//'     representing variants from the reference genome.
//' @param bind_sites Vector of enzyme full recognition site(s).
//' @param len5s A vector of the numbers of characters of the prime5 sites for each
//'     recognition site.
//' @param chunk_size The size of chunks to divide sequences into when digesting.
//' @param n_cores The number of cores to use for processing.
//'
//' @return A list of lists, each sub-list containing multiple vectors representing
//'     the locations of cut sites for a given variant on a given sequence.
//'     Indexing the output list would be done as such:
//'     \code{output_list[[variant_index]][[sequence_index]][position_index]}.
//'
//' @noRd
//'
// [[Rcpp::export]]
std::vector< std::vector< std::deque<uint> > > digest_var(
        SEXP var_,
        const std::vector<std::string>& bind_sites,
        const std::vector<uint>& len5s,
        const uint& chunk_size,
        const uint& n_cores) {

    XPtr<VarSet> var_xptr(var_);
    VarSet& var(*var_xptr);

    const uint n_vars = var.size();
    const uint n_seqs = var.reference.size();

    const DigestInfo dinfo(bind_sites, len5s);

    if (chunk_size < (2 * dinfo.max_size - 1)) {
        stop("chunk size needs to be at least 2 times the maximum binding "
                 "site length minus 1 (i.e., `chunk_size >= 2 * "
                 "max(<binding site sizes>) - 1`)");
    }

    std::vector< std::vector< std::deque<uint> > > out_vec(
            n_vars, std::vector< std::deque<uint> >(
                    n_seqs, std::deque<uint>(0)));

    #ifdef _OPENMP
    #pragma omp parallel for default(shared) num_threads(n_cores) schedule(dynamic) collapse(2)
    #endif
    for (uint v = 0; v < n_vars; v++) {
        for (uint s = 0; s < n_seqs; s++) {
            VarSeqDigest vsd(var, v, s, chunk_size, out_vec[v][s]);
            vsd.digest(dinfo);
        }
    }

    return out_vec;
}











/*
 =============================================================================
 =============================================================================
 =============================================================================
 =============================================================================

 Reference genome digestion

 =============================================================================
 =============================================================================
 =============================================================================
 =============================================================================
 */




// Digest chunk of a sequence

void RefSeqChunk::digest(const DigestInfo& dinfo) {

    // Reset this if it's been digested already
    if (!cut_sites.empty()) cut_sites.clear();

    /*
     `next_start` is in the RefSeqChunk class and keeps track of earliest possible
     start for the next chunk in this sequence (`start` is const, `next_start` isn't):
     */
    next_start = start;

    digest_seq(cut_sites, next_start, ref, dinfo, 0, end);

    return;
}


/*
 Make sure no positions overlap in the "seam" between two chunks.
 `prev` should come before the focal `RefSeqChunk`
 */
void RefSeqChunk::merge_seam(const RefSeqChunk& prev, const DigestInfo& dinfo) {

    // If they're not on the same sequence, don't do anything
    if (prev.seq != seq) return;
    // If there were no cut sites on this chunk or the previous one, don't do anything
    if (cut_sites.empty() || prev.cut_sites.empty()) return;

    std::deque<uint>& pos_vec2(cut_sites);

    // Earliest possible starting position in `pos_vec2`
    uint start_pos = prev.next_start;

    if (pos_vec2.front() >= start_pos) return;

    // First remove position(ref) from `pos_vec2` that are < `start_pos`
    while (pos_vec2.front() < start_pos) {
        pos_vec2.pop_front();
        if (pos_vec2.empty()) break;
    }

    // Now search through sequence again until positions are "merged"
    MultiOut search_obj;
    multi_search(search_obj, dinfo, ref, start_pos);
    while (search_obj.new_pos != std::string::npos) {
        uint cut_site_pos = search_obj.new_pos + dinfo.len5s[search_obj.j];
        while (pos_vec2.front() < cut_site_pos) pos_vec2.pop_front();
        if (pos_vec2.front() == cut_site_pos) {  // <---------- how I'm defining "merged"
            break;
        } else {
            pos_vec2.push_front(cut_site_pos);
        }
        start_pos = search_obj.new_pos + dinfo.to_skip(search_obj.j);
        multi_search(search_obj, dinfo, ref, start_pos);
    }

    /*
     If we've changed the whole deque of cut sites, then adjust the `next_start` field
     (this should virtually never happen with adequately large chunk sizes)
     */
    if (start_pos > next_start) next_start = start_pos;

    return;
}






//' Internal C++ function to digest all sequences in a reference genome.
//'
//' This function shouldn't be exported after testing is finished.
//'
//'
//' @param ref_ An external pointer to a C++ \code{RefGenome} object
//'     representing the reference genome.
//' @param bind_sites Vector of enzyme full recognition site(s).
//' @param len5s A vector of the numbers of characters of the prime5 sites for each
//'     recognition site.
//' @param n_cores The number of cores to use for processing. This value is ignored
//'     if the input reference genome is merged and \code{chunk_size == 0}.
//'     Defaults to \code{1}.
//' @param chunk_size Size of chunks to break sequences into for processing.
//'     This value is ignored if it's set to zero.
//'     Ideally this is set to a value that results in a number of chunks divisible by
//'     the number of cores you're using, and is most useful when `n_cores` is greater
//'     than the number of scaffolds.
//'     Breaking into increasingly small chunks results in increasing overhead, so
//'     beware of making this argument very small.
//'     Reference genome sequences are not copied during this function, so using
//'     this argument for a reference genome does NOT decrease memory usage
//'     appreciably.
//'     Defaults to \code{0}.
//'
//' @return A list of vectors, each vector representing the locations of cut sites
//'     on a given sequence.
//'     Indexing the output list would be done as such:
//'     \code{output_list[[sequence_index]][position_index]}.
//'
//' @noRd
//'
//[[Rcpp::export]]
std::vector< std::deque<uint> > digest_ref(
        SEXP ref_,
        const std::vector<std::string>& bind_sites,
        const std::vector<uint>& len5s,
        const uint& n_cores = 1,
        const uint& chunk_size = 0) {

    const XPtr<RefGenome> ref(ref_);
    const RefGenome& reference(*ref);

    std::vector< std::deque<uint> > out_vec(reference.size());
    DigestInfo dinfo(bind_sites, len5s);

    if (chunk_size > 0) {
        if (chunk_size < (2 * dinfo.max_size - 1)) {
            stop("chunk size needs to be at least 2 times the maximum binding "
                     "site length minus 1 (i.e., `chunk_size >= 2 * "
                     "max(<binding site sizes>) - 1`)");
        }
        /*
         I'm iterating by `chunk_size - (max_size - 1)` to make sure I catch any
        matches that overlap the "seam" between chunks.
        */
        uint iter_by = chunk_size - (dinfo.max_size - 1);
        // Set up temporary deque containing a deque of RefSeqChunk objects:
        std::deque<RefSeqChunk> chunk_dq;
        for (uint i = 0; i < reference.size(); i++) {
            const std::string& seq(reference[i].nucleos);

            if (chunk_size >= seq.size()) {
                RefSeqChunk tmp_sc(seq, i, 0, chunk_size);
                chunk_dq.push_back(tmp_sc);
            } else {
                /*
                 The last chunk start position (without allowing only overlap
                 at the end):
                 */
                uint last_chunk = seq.size() - dinfo.max_size;
                for (uint s = 0; s < last_chunk; s += iter_by) {
                    RefSeqChunk tmp_sc(seq, i, s, chunk_size);
                    chunk_dq.push_back(tmp_sc);
                }
            }
        }

        #ifdef _OPENMP
        #pragma omp parallel for default(shared) num_threads(n_cores) schedule(static)
        #endif
        for (uint i = 0; i < chunk_dq.size(); i++) {
            chunk_dq[i].digest(dinfo);
        }

        /*
         Now check for points over the seams and make sure they aren't too close to
         each other.
         */
        #ifdef _OPENMP
        #pragma omp parallel for default(shared) num_threads(n_cores) schedule(static)
        #endif
        for (uint i = 1; i < chunk_dq.size(); i++) {
            chunk_dq[i].merge_seam(chunk_dq[(i-1)], dinfo);
        }

        /*
         Lastly combine all RefSeqChunk objects for a given sequence together into
         one deque.
         */
        for (uint i = 0; i < reference.size(); i++) {
            std::deque<uint>& out_dq(out_vec[i]);
            while (chunk_dq.front().seq == i) {
                chunk_dq.front().combine(out_dq);
                chunk_dq.pop_front();
                if (chunk_dq.empty()) break;
            }
        }
    } else if (reference.merged) {
        uint i = 0, j = 0; // j is only used for the call to digest_seq
        const std::string& seq(reference[i].nucleos);
        digest_seq(out_vec[i], j, seq, dinfo);
    } else {
        #ifdef _OPENMP
        #pragma omp parallel for default(shared) num_threads(n_cores) schedule(dynamic)
        #endif
        for (uint i = 0; i < reference.size(); i++) {
            const std::string& seq_s(reference[i].nucleos);
            uint j = 0; // j is only used for the call to digest_seq
            digest_seq(out_vec[i], j, seq_s, dinfo);
        }
    }

    return out_vec;
}

