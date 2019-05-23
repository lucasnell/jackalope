#ifndef __JACKAL_SEQUENCER_H
#define __JACKAL_SEQUENCER_H



#include <RcppArmadillo.h>
#include <vector>  // vector class
#include <string>  // string class
#include <pcg/pcg_random.hpp> // pcg prng
#include <fstream> // for writing FASTQ files
#include "zlib.h"  // for writing to compressed FASTQ
#include <progress.hpp>  // for the progress bar


// for bgzip_file
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "htslib/bgzf.h"

#include "jackalope_types.h"  // uint64
#include "pcg.h"  // ruinf_01
#include "util.h"  // str_stop, thread_check, split_int
#include "io.h"  // File* types
#include "alias_sampler.h"  // Alias sampler


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



//' bgzip a file, potentially using multiple threads.
//'
//' @noRd
//'
inline int bgzip_file(const std::string& file_name,
                      const int& n_threads,
                      const int& compress) {

    // Make writer object:
    FileBGZF bgzf(file_name, n_threads, compress);
    // Others used below:
    void *buffer;
    int c;
    struct stat sbuf;
    int f_src = fileno(stdin);

    const int WINDOW_SIZE = 64 * 1024;


    // Checking source file:
    if (stat(file_name.c_str(), &sbuf) < 0) {
        str_stop({"\nIn bgzip step, file ", file_name,
                 " had non-zero status: ", strerror(errno), "."});
    }
    if ((f_src = open(file_name.c_str(), O_RDONLY)) < 0) {
        str_stop({"\nIn bgzip step, file ", file_name, " could not be opened."});
    }

    // Create buffer of info to pass between:
    buffer = malloc(WINDOW_SIZE);
#ifdef _WIN32
    _setmode(f_src, O_BINARY);
#endif
    // Writing info from one to another
    while ((c = read(f_src, buffer, WINDOW_SIZE)) > 0) {
        bgzf.write(buffer, c);
    }

    bgzf.close();
    unlink(file_name.c_str());
    free(buffer);
    close(f_src);
    return 0;

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
    const uint64 n_reads;           // # reads to create
    const uint64 read_pool_size;    // reads per pool
    uint64 reads_made;              // Number of reads already made
    uint64 reads_in_pool;           // Number of reads in current pool
    bool do_write;                  // Whether to write to file
    const double prob_dup;          // probability of duplicate
    const uint64 n_read_ends;       // (1 for SE Illumina or PacBio, 2 for PE Illumina)
    std::vector<std::vector<char>> fastq_pools;

    ReadWriterOneThread(const T& read_filler_base,
                        const uint64& n_reads_,
                        const uint64& read_pool_size_,
                        const double& prob_dup_,
                        const uint64& n_read_ends_)
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
        for (uint64 i = 0; i < fastq_pools.size(); i++) {
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


    inline uint64 pool_size() {
        if (fastq_pools.size() == 0) return 0ULL;
        return fastq_pools[0].size() * fastq_pools.size();
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
                                      const uint64& n_reads,
                                      const double& prob_dup,
                                      const uint64& read_pool_size,
                                      const uint64& n_read_ends,
                                      const uint64& n_threads,
                                      const int& compress,
                                      Progress& prog_bar) {

    const std::vector<uint64> reads_per_thread = split_int(n_reads, n_threads);

    // Generate seeds for random number generators (1 RNG per thread)
    const std::vector<std::vector<uint64>> seeds = mt_seeds(n_threads);

    // Create and open files:
    std::vector<F> files(n_read_ends);
    for (uint64 i = 0; i < n_read_ends; i++) {
        std::string file_name = out_prefix + "_R" + std::to_string(i+1)+ ".fq";
        files[i].set(file_name, compress);
    }


#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads) shared(files) if (n_threads > 1)
{
#endif

    std::vector<uint64> active_seeds;

    // Write the active seed per thread or just write one of the seeds.
#ifdef _OPENMP
    uint64 active_thread = omp_get_thread_num();
#else
    uint64 active_thread = 0;
#endif
    active_seeds = seeds[active_thread];

    pcg64 eng = seeded_pcg(active_seeds);

    uint64 reads_this_thread = reads_per_thread[active_thread];

    ReadWriterOneThread<T,F> writer(read_filler_base, reads_this_thread,
                                    read_pool_size, prob_dup, n_read_ends);

    uint64 reads_written;
    uint64 old_reads = 0;
    uint64 new_reads = 0;

    while (writer.reads_made < reads_this_thread) {

        old_reads = writer.pool_size();

        writer.create_reads(eng);

        new_reads += (writer.pool_size() - old_reads);

        /*
         Every 10,000 characters created, check that the user hasn't
         interrupted the process.
         (Doing it this way makes the check approximately the same between
          illumina and pacbio.)
         */
        if (new_reads > 10000) {
            if (prog_bar.check_abort()) break;
            new_reads = 0;
        }

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

        }
    }

#ifdef _OPENMP
}
#endif

    // Close files
    for (uint64 i = 0; i < files.size(); i++) {
        files[i].close();
    }

    return;
};




/*
 For one sequencing-filler type, make reads and write them to file(s).

 `T` should be `[Illumina|PacBio]Reference` or `[Illumina|PacBio]Variants`.

 */

template <typename T>
inline void write_reads_cpp_(const T& read_filler_base,
                             std::string out_prefix,
                             const uint64& n_reads,
                             const double& prob_dup,
                             const uint64& read_pool_size,
                             const uint64& n_read_ends,
                             uint64 n_threads,
                             const int& compress,
                             const std::string& comp_method,
                             Progress& prog_bar) {

    expand_path(out_prefix);

    // To make sure the following if statements are still accurate if OpenMP not used,
    // check that # threads isn't too high and change to 1 if not using OpenMP:
    thread_check(n_threads);

    // Compressed output run serially
    if (compress > 0 && n_threads == 1) {

        if (comp_method == "gzip") {
            write_reads_one_filetype_<T, FileGZ>(
                    read_filler_base, out_prefix, n_reads, prob_dup,
                    read_pool_size, n_read_ends, n_threads, compress, prog_bar);
        } else if (comp_method == "bgzip") {
            write_reads_one_filetype_<T, FileBGZF>(
                    read_filler_base, out_prefix, n_reads, prob_dup,
                    read_pool_size, n_read_ends, n_threads, compress, prog_bar);
        } else stop("\nUnrecognized compression method.");

    /*
     Compressed output run in parallel.
     The only way I've found to make this actually have a speed advantage over
     running serially is to make it first write to uncompressed output in parallel,
     then do the compression using `BGZF` also in parallel.
     */
    } else if (compress > 0 && n_threads > 1) {

        // First do it uncompressed:
        write_reads_one_filetype_<T, FileUncomp>(
                read_filler_base, out_prefix, n_reads, prob_dup,
                read_pool_size, n_read_ends, n_threads, compress, prog_bar);
        for (uint64 i = 0; i < n_read_ends; i++) {
            std::string file_name = out_prefix + "_R" + std::to_string(i+1)+ ".fq";
            bgzip_file(file_name, static_cast<int>(n_threads),
                       static_cast<int>(compress));
            // For progress bar, I assume compression takes half as long:
            prog_bar.increment(n_reads / (n_read_ends * 2));
        }
    // Uncompressed output run in serial or parallel
    } else {
        write_reads_one_filetype_<T, FileUncomp>(
                read_filler_base, out_prefix, n_reads, prob_dup,
                read_pool_size, n_read_ends, n_threads, compress, prog_bar);
    }


}



/*
 Samples for # reads per variant.
 It uses binomial distribution for a `n_reads` that decreases each iteration and a
 `variant_probs` vector whose values go up.
 This is ~600x faster than doing separate samples (via AliasSampler) for every read.
 */
inline std::vector<uint64> read_per_var(uint64 n_reads,
                                        std::vector<double> variant_probs) {

    pcg64 eng = seeded_pcg();

    double sum_probs = std::accumulate(variant_probs.begin(), variant_probs.end(), 0.0);
    for (double& p : variant_probs) p /= sum_probs;

    std::vector<uint64> out;
    out.reserve(variant_probs.size());

    std::binomial_distribution<uint64> distr(n_reads, 0.5);

    for (uint64 i = 0; i < variant_probs.size(); i++) {

        // Update distribution:
        distr.param(std::binomial_distribution<uint64>::param_type(
                n_reads, variant_probs[i]));

        // Sample from binomial distribution:
        out.push_back(distr(eng));

        // Decrease `n_reads` based on # sampled:
        n_reads -= out.back();
        // If `n_reads` is zero, we can resize `out` and stop now:
        if (n_reads == 0) {
            out.resize(variant_probs.size(), 0ULL);
            break;
        }

        // Increase probabilities:
        sum_probs = 1 - variant_probs[i];
        for (uint64 j = (i+1); j < variant_probs.size(); j++) {
            variant_probs[j] /= sum_probs;
        }
    }

    return out;

}




/*
 Same as above, but for when you want a separate file per variant
 */

template <typename T>
inline void write_reads_cpp_sep_files_(const VarSet& var_set,
                                       const std::vector<double>& variant_probs,
                                       T& read_filler_base,
                                       const std::string& out_prefix,
                                       const uint64& n_reads,
                                       const double& prob_dup,
                                       const uint64& read_pool_size,
                                       const uint64& n_read_ends,
                                       const uint64& n_threads,
                                       const int& compress,
                                       const std::string& comp_method,
                                       Progress& prog_bar) {

    // Sample for reads per file:
    std::vector<uint64> reads_per_file = read_per_var(n_reads / n_read_ends,
                                                      variant_probs);
    if (n_read_ends > 1) for (uint64& rpf : reads_per_file) rpf *= n_read_ends;

    // Now do each variant as separate file:
    std::vector<double> var_probs_(var_set.size(), 0.0);

    for (uint64 i = 0; i < var_set.size(); i++) {

        if (prog_bar.check_abort()) break;

        var_probs_[i] = 1;
        read_filler_base.variant_sampler = AliasSampler(var_probs_);

        std::string out_prefix_ = out_prefix + '_' + var_set[i].name;

        write_reads_cpp_<T>(
            read_filler_base, out_prefix_, reads_per_file[i], prob_dup,
            read_pool_size, n_read_ends, n_threads, compress, comp_method,
            prog_bar);

        var_probs_[i] = 0;
    }

    return;
}





#endif
