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


#include "jackal_types.h" // integer types
#include "seq_classes_ref.h"  // Ref* classes
#include "seq_classes_var.h"  // Var* classes
#include "digest.h" // DigestInfo, MultiOut,
#include "str_manip.h" // rev_comp

using namespace Rcpp;



/*
 Remove duplicates from a vector of strings in place.
 */
void unique_(std::vector<std::string>& x) {
    std::sort(x.begin(), x.end());
    std::vector<std::string>::iterator iter = std::unique(x.begin(), x.end());
    x.resize(iter - x.begin());
    return;
}




//' Calculate how many bases come before a cleavage site.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
std::vector<uint32> get_precleavage_lens(const std::vector<std::string>& seqs) {

    std::vector<uint32> out(seqs.size());

    for (uint32 i = 0; i < seqs.size(); i++) {
        uint32 cleav = seqs[i].find('/');
        out[i] = cleav;
    }

    return out;
}


/*
 Expand sequences with non-specific nucleobases.
 */
void expand_sites(const std::vector<std::string>& sites,
                  std::vector<std::string>& seqs_out,
                  const bool& add_rev_comp = true) {

    uint32 n_combs = 1;
    for (uint32 i = 0; i < sites.size(); i++) n_combs *= sites[i].size();
    seqs_out.reserve(n_combs);

    for (uint32 i = 0; i < sites[0].size(); i++) {
        seqs_out.push_back(std::string(1, sites[0][i]));
    }

    for (uint32 i = 1; i < sites.size(); i++) {
        const std::string& site_i(sites[i]);
        uint32 n = seqs_out.size();
        for (uint32 j = 1; j < site_i.size(); j++) {
            for (uint32 k = 0; k < n; k++) {
                seqs_out.push_back(seqs_out[k] + site_i[j]);
            }
        }
        for (uint32 k = 0; k < n; k++) {
            seqs_out[k] += site_i[0];
        }
    }
    // Add all reverse complements, whether or not they already exist in the output vector
    if (add_rev_comp) {
        seqs_out.reserve(n_combs * 2);
        for (uint32 i = 0; i < n_combs; i++) {
            std::string rc = seqs_out[i];
            rev_comp(rc);
            seqs_out.push_back(rc);
        }
    }

    return;
}




//' Expand sequences for reverse complements and for non-specific nucleobases.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
std::vector<std::string> expand_seqs(const std::vector<std::string>& seqs) {

    std::vector<std::string> out;
    std::string bases = "TCAG";

    for (uint32 i = 0; i < seqs.size(); i++) {
        const std::string& seq(seqs[i]);
        // See if the sequence is degenerate (i.e., ambiguous) or not:
        bool degenerate = false;
        for (const char& c : seq) {
            if (bases.find(c) == std::string::npos) degenerate = true;
        }
        if (!degenerate) {
            out.push_back(seq);
            std::string rc = seq;
            rev_comp(rc);
            if (rc != seq) out.push_back(rc);
        } else {
            std::vector<std::string> sites(seq.size());
            for (uint32 j = 0; j < seq.size(); j++) {
                sites[j] = digest::deg_bases[seq[j]];
            }
            // Create all combos of these sites and add to `out`
            expand_sites(sites, out);
        }
    }
    // Remove duplicates:
    unique_(out);

    return out;
}






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
        const uint32& last_pos) {

    search_obj.reset();

    std::string::size_type new_pos;

    for (uint32 j = 0; j < dinfo.bind_sites.size(); j++) {
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
void digest_seq(std::deque<uint32>& pos_vec,
                uint32& start_pos,
                const std::string& s,
                const DigestInfo& dinfo,
                const uint32& pos_mod = 0,
                const uint32& end = 0) {

    uint32 cut_site_pos;
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
    uint32 mut_i = 0;

    // Initialize chunk sequence string
    std::string chunk_str(chunk_size, 'x');

    /*
     I'm iterating by `chunk_size - (max_size - 1)` to make sure I catch any
    matches that overlap the "seam" between chunks.
    */
    uint32 it = chunk_size - (dinfo.max_size - 1);

    /*
     The last chunk start position in this sequence (without allowing only overlap
     at the end)
     */
    uint32 last_chunk = var_info.seq_size - dinfo.max_size;

    // Digesting each chunk:
    for (uint32 chunk_start = 0; chunk_start <= last_chunk; chunk_start += it) {
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
//' @param var_set_ptr An external pointer to a C++ \code{VarSet} object
//'     representing variants from the reference genome.
//' @inheritParams bind_sites digest_ref
//' @inheritParams len5s digest_ref
//' @param chunk_size The size of chunks to divide sequences into when digesting.
//' @param n_cores The number of cores to use for processing. Defaults to \code{1}.
//'
//' @return A list of lists, each sub-list containing multiple vectors representing
//'     the locations of cut sites for a given variant on a given sequence.
//'     Indexing the output list would be done as such:
//'     \code{output_list[[variant_index]][[sequence_index]][position_index]}.
//'
//' @noRd
//'
// [[Rcpp::export]]
SEXP digest_var_set(
        SEXP var_set_ptr,
        const std::vector<std::string>& bind_sites,
        const std::vector<uint32>& len5s,
        const uint32& chunk_size,
        const uint32& n_cores = 1) {

    XPtr<VarSet> var_set(var_set_ptr);

    const uint32 n_vars = var_set->size();
    const uint32 n_seqs = var_set->reference->size();

    XPtr<GenomeSetDigest> out_digests(new GenomeSetDigest(n_vars, n_seqs));

    const DigestInfo dinfo(bind_sites, len5s);

    if (chunk_size < (2 * dinfo.max_size - 1)) {
        stop("chunk size needs to be at least 2 times the maximum binding "
                 "site length minus 1 (i.e., `chunk_size >= 2 * "
                 "max(<binding site sizes>) - 1`)");
    }


    #ifdef _OPENMP
    #pragma omp parallel for default(shared) num_threads(n_cores) schedule(dynamic) collapse(2)
    #endif
    for (uint32 v = 0; v < n_vars; v++) {
        for (uint32 s = 0; s < n_seqs; s++) {
            VarSeqDigest vsd(*var_set, v, s, chunk_size, (*out_digests)[v][s]);
            vsd.digest(dinfo);
        }
    }

    return out_digests;
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

    std::deque<uint32>& pos_vec2(cut_sites);

    // Earliest possible starting position in `pos_vec2`
    uint32 start_pos = prev.next_start;

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
        uint32 cut_site_pos = search_obj.new_pos + dinfo.len5s[search_obj.j];
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
//'
//'
//' @param ref_genome_ptr An external pointer to a C++ \code{RefGenome} object
//'     representing the reference genome.
//' @param bind_sites Vector of enzyme full recognition site(s).
//' @param len5s A vector of the numbers of characters of the prime5 sites for each
//'     recognition site.
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
//' @param n_cores The number of cores to use for processing. This value is ignored
//'     if the input reference genome is merged and \code{chunk_size == 0}.
//'     Defaults to \code{1}.
//'
//' @return A list of vectors, each vector representing the locations of cut sites
//'     on a given sequence.
//'     Indexing the output list would be done as such:
//'     \code{output_list[[sequence_index]][position_index]}.
//'
//' @noRd
//'
//[[Rcpp::export]]
SEXP digest_ref(
        SEXP ref_genome_ptr,
        const std::vector<std::string>& bind_sites,
        const std::vector<uint32>& len5s,
        const uint32& chunk_size = 0,
        const uint32& n_cores = 1) {

    const XPtr<RefGenome> ref_genome(ref_genome_ptr);

    XPtr<GenomeDigest> out_digest_xptr(new GenomeDigest(ref_genome->size()));
    GenomeDigest& out_digest(*out_digest_xptr);

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
        uint32 iter_by = chunk_size - (dinfo.max_size - 1);
        // Set up temporary deque containing a deque of RefSeqChunk objects:
        std::deque<RefSeqChunk> chunk_dq;
        for (uint32 i = 0; i < ref_genome->size(); i++) {
            const std::string& seq((*ref_genome)[i].nucleos);

            if (chunk_size >= seq.size()) {
                RefSeqChunk tmp_sc(seq, i, 0, chunk_size);
                chunk_dq.push_back(tmp_sc);
            } else {
                /*
                 The last chunk start position (without allowing only overlap
                 at the end):
                 */
                uint32 last_chunk = seq.size() - dinfo.max_size;
                for (uint32 s = 0; s < last_chunk; s += iter_by) {
                    RefSeqChunk tmp_sc(seq, i, s, chunk_size);
                    chunk_dq.push_back(tmp_sc);
                }
            }
        }

        #ifdef _OPENMP
        #pragma omp parallel for default(shared) num_threads(n_cores) schedule(static)
        #endif
        for (uint32 i = 0; i < chunk_dq.size(); i++) {
            chunk_dq[i].digest(dinfo);
        }

        /*
         Now check for points over the seams and make sure they aren't too close to
         each other.
         */
        #ifdef _OPENMP
        #pragma omp parallel for default(shared) num_threads(n_cores) schedule(static)
        #endif
        for (uint32 i = 1; i < chunk_dq.size(); i++) {
            chunk_dq[i].merge_seam(chunk_dq[(i-1)], dinfo);
        }

        /*
         Lastly combine all RefSeqChunk objects for a given sequence together into
         one deque.
         */
        for (uint32 i = 0; i < ref_genome->size(); i++) {
            std::deque<uint32>& out_dq(out_digest[i]);
            while (chunk_dq.front().seq == i) {
                chunk_dq.front().combine(out_dq);
                chunk_dq.pop_front();
                if (chunk_dq.empty()) break;
            }
        }
    } else if (ref_genome->merged) {
        uint32 i = 0, j = 0; // j is only used for the call to digest_seq
        const std::string& seq((*ref_genome)[i].nucleos);
        digest_seq(out_digest[i], j, seq, dinfo);
    } else {
        #ifdef _OPENMP
        #pragma omp parallel for default(shared) num_threads(n_cores) schedule(dynamic)
        #endif
        for (uint32 i = 0; i < ref_genome->size(); i++) {
            const std::string& seq_s((*ref_genome)[i].nucleos);
            uint32 j = 0; // j is only used for the call to digest_seq
            digest_seq(out_digest[i], j, seq_s, dinfo);
        }
    }

    return out_digest_xptr;
}

