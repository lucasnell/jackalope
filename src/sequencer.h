#ifndef __JACKAL_SEQUENCER_H
#define __JACKAL_SEQUENCER_H



#include <RcppArmadillo.h>
#include <vector>  // vector class
#include <string>  // string class
#include <pcg/pcg_random.hpp> // pcg prng
#include <fstream> // for writing FASTQ files
#include "zlib.h"  // for writing to compressed FASTQ
#include <progress.hpp>  // for the progress bar




#include "jackalope_types.h"  // uint32
#include "pcg.h"  // ruinf_01
#include "util.h"  // str_stop
#include "read_write.h"  // File* types


using namespace Rcpp;



namespace sequencer {

// Goes from character (coerced to integer) to index from 0:3 (4 for non-nucleotide)
const std::vector<uint8> nt_map = {
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
    4,2,4,1,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,4,0,4,4,4,4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4};

const std::vector<std::string> mm_nucleos = {"CAG", "TAG", "TCG", "TCA", "NNN"};

}





/*

Info for making reads and writing them, for one thread.

`T` can be `[Illumina/PacBio]Reference` or `[Illumina/PacBio]Variants`.
 `F` should be `FileUncomp`, `FileGZ`, or `FileBGZF`.
*/

template <typename T, typename F>
class ReadWriterOneThread {

public:

    T read_filler;
    const uint32 n_reads;           // # reads to create
    const uint32 read_pool_size;    // reads per pool
    uint32 reads_made;              // Number of reads already made
    uint32 reads_in_pool;           // Number of reads in current pool
    bool do_write;                  // Whether to write to file
    const double prob_dup;          // probability of duplicate
    const uint32 n_read_ends;       // (1 for SE Illumina or PacBio, 2 for PE Illumina)
    std::vector<std::vector<char>> fastq_pools;

    ReadWriterOneThread(const T& read_filler_base,
                        const uint32& n_reads_,
                        const uint32& read_pool_size_,
                        const double& prob_dup_,
                        const uint32& n_read_ends_)
        : read_filler(read_filler_base),
          n_reads(n_reads_),
          read_pool_size(read_pool_size_),
          reads_made(0),
          reads_in_pool(0),
          do_write(false),
          prob_dup(prob_dup_),
          n_read_ends(n_read_ends_),
          fastq_pools(n_read_ends_) {};


    // Write contents in `fastq_pools` to UNcompressed file(s).
    void write(std::vector<F>& files) {
        for (uint32 i = 0; i < fastq_pools.size(); i++) {
            files[i].write(fastq_pools[i]);
            fastq_pools[i].clear();
        }
        reads_in_pool = 0;
        do_write = false;
        return;
    }
    /* Overloaded for one file: */
    void write(F& file) {
        file.write(fastq_pools[0]);
        fastq_pools[0].clear();
        reads_in_pool = 0;
        do_write = false;
        return;
    }


    /*
     Add new read(s) to `fastq_pools`, and update bool for whether you should
     write to file
     */
    void create_reads(pcg64& eng) {
        read_filler.template one_read<std::vector<char>>(fastq_pools, eng);
        reads_made += n_read_ends;
        reads_in_pool += n_read_ends;
        double dup = runif_01(eng);
        while (dup < prob_dup && reads_made < n_reads &&
               reads_in_pool < read_pool_size) {
            read_filler.template re_read<std::vector<char>>(fastq_pools, eng);
            reads_made += n_read_ends;
            reads_in_pool += n_read_ends;
            dup = runif_01(eng);
        }
        do_write = reads_in_pool >= read_pool_size || reads_made >= n_reads;
        return;
    }



};








//' Split number of reads by number of threads.
//'
//' @noRd
//'
inline std::vector<uint32> split_n_reads(const uint32& n_reads,
                                         const uint32& n_threads) {
    std::vector<uint32> out(n_threads, n_reads / n_threads);
    uint32 sum_reads = std::accumulate(out.begin(), out.end(), 0U);
    uint32 i = 0;
    while (sum_reads < n_reads) {
        out[i]++;
        i++;
        sum_reads++;
    }
    return out;
}






/*
 ======================================================================================
 ======================================================================================

 Input/output

 ======================================================================================
 ======================================================================================
 */





/*
 For one file type and read filler type, make Illumina reads and write them to file(s).
 Does most of the work of `write_reads_cpp_` below.
 This should only be called inside that function.

 `T` should be `[Illumina|PacBio]Reference` or `[Illumina|PacBio]Variants`.
 `F` should be `FileUncomp`, `FileGZ`, or `FileBGZF`.

 */
template <typename T, typename F>
inline void write_reads_one_filetype_(const T& read_filler_base,
                                      const std::string& out_prefix,
                                      const uint32& n_reads,
                                      const double& prob_dup,
                                      const uint32& read_pool_size,
                                      const uint32& n_read_ends,
                                      const uint32& n_threads,
                                      const bool& show_progress,
                                      const int& compress) {

    const std::vector<uint32> reads_per_thread = split_n_reads(n_reads, n_threads);

    // Generate seeds for random number generators (1 RNG per thread)
    const std::vector<std::vector<uint64>> seeds = mt_seeds(n_threads);

    // Progress bar
    Progress prog_bar(n_reads, show_progress);

    // Create and open files:
    std::vector<F> files;
    for (uint32 i = 0; i < n_read_ends; i++) {
        std::string file_name = out_prefix + "_R" + std::to_string(i+1)+ ".fq";
        files.push_back(F(file_name, compress));
    }


#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads) shared(files) if (n_threads > 1)
{
#endif

    std::vector<uint64> active_seeds;

    // Write the active seed per thread or just write one of the seeds.
#ifdef _OPENMP
    uint32 active_thread = omp_get_thread_num();
#else
    uint32 active_thread = 0;
#endif
    active_seeds = seeds[active_thread];

    pcg64 eng = seeded_pcg(active_seeds);

    uint32 reads_this_thread = reads_per_thread[active_thread];

    ReadWriterOneThread<T,F> writer(read_filler_base, reads_this_thread,
                                    read_pool_size, prob_dup, n_read_ends);

    uint32 reads_written;
    uint32 n_chars = 0;

    while (writer.reads_made < reads_this_thread) {

        /*
         Every 10,000 characters written, check that the user hasn't
         interrupted the process.
         (Doing it this way makes the check approximately the same between
          illumina and pacbio.)
         */
        if (n_chars > 10000) {
            if (prog_bar.check_abort()) break;
            n_chars = 0;
        }

        writer.create_reads(eng);

        if (writer.do_write) {
            // Save info for progress bar:
            reads_written = writer.reads_in_pool;
            // Write to files:
#ifdef _OPENMP
#pragma omp critical
{
#endif
            writer.write(files);
#ifdef _OPENMP
}
#endif
            // Increment progress bar
            prog_bar.increment(reads_written);
            // Increment # chars written
            n_chars += (writer.fastq_pools[0].size() * writer.fastq_pools.size());
        }
    }

#ifdef _OPENMP
}
#endif

    // Close files
    for (uint32 i = 0; i < files.size(); i++) {
        files[i].close();
    }

    return;
};




/*
 For one file type, make Illumina reads and write them to file(s).
 Does most of the work of `write_reads_cpp_` below.

 `T` should be `[Illumina|PacBio]Reference` or `[Illumina|PacBio]Variants`.
 `F` should be `FileUncomp`, `FileGZ`, or `FileBGZF`.

 */

template <typename T>
inline void write_reads_cpp_(const T& read_filler_base,
                             std::string out_prefix,
                             const uint32& n_reads,
                             const double& prob_dup,
                             const uint32& read_pool_size,
                             const uint32& n_read_ends,
                             uint32 n_threads,
                             const bool& show_progress,
                             const int& compress,
                             const std::string& comp_method) {

    expand_path(out_prefix);

    // Compressed output run serially
    if (compress > 0 && n_threads == 1) {

        if (comp_method == "gzip") {
            write_reads_one_filetype_<T, FileGZ>(
                    read_filler_base, out_prefix, n_reads, prob_dup,
                    read_pool_size, n_read_ends, n_threads, show_progress, compress);
        } else if (comp_method == "bgzip") {
            write_reads_one_filetype_<T, FileBGZF>(
                    read_filler_base, out_prefix, n_reads, prob_dup,
                    read_pool_size, n_read_ends, n_threads, show_progress, compress);
        } else stop("\nUnrecognized compression method.");

    /*
     Compressed output run in parallel.
     The only way I've found to make this actually have a speed advantage over
     running serially is to make it first write to uncompressed output in parallel,
     then do the compression using `BGZF` also in parallel.
     */
    } else if (compress > 0 && n_threads > 1) {

        if (show_progress) {
            Rcout << "Progress for writing uncompressed reads..." << std::endl;
        }
        // First do it uncompressed:
        write_reads_one_filetype_<T, FileUncomp>(
                read_filler_base, out_prefix, n_reads, prob_dup,
                read_pool_size, n_read_ends, n_threads, show_progress, compress);
        // Then compress using multiple threads:
        if (show_progress) Rcout << std::endl << "Now compressing reads..." << std::endl;
        for (uint32 i = 0; i < n_read_ends; i++) {
            std::string file_name = out_prefix + "_R" + std::to_string(i+1)+ ".fq";
            bgzip_file(file_name, static_cast<int>(n_threads),
                       static_cast<int>(compress));
            if (show_progress) Rcout << "finished compressing " << file_name << std::endl;
        }
    // Uncompressed output run in serial or parallel
    } else {
        write_reads_one_filetype_<T, FileUncomp>(
                read_filler_base, out_prefix, n_reads, prob_dup,
                read_pool_size, n_read_ends, n_threads, show_progress, compress);
    }


}








#endif
