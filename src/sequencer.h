#ifndef __GEMINO_SEQUENCER_H
#define __GEMINO_SEQUENCER_H



#include <RcppArmadillo.h>
#include <vector>  // vector class
#include <string>  // string class
#include <pcg/pcg_random.hpp> // pcg prng
#include <fstream> // for writing FASTQ files
#include <zlib.h>  // for writing to compressed FASTQ
#include <progress.hpp>  // for the progress bar

#include "gemino_types.h"  // uint32
#include "pcg.h"  // ruinf_01

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

Info for making reads and writing them, for one core.

`T` can be `[Illumina/PacBio]Reference` or `[Illumina/PacBio]Variants`.
*/

template <typename T>
class ReadWriterOneCore {

public:

    T read_filler;
    const uint32 n_reads;           // # reads to create
    const uint32 read_chunk_size;   // reads per chunk
    uint32 reads_made;              // Number of reads already made
    uint32 reads_in_chunk;          // Number of reads in current chunk
    bool do_write;                  // Whether to write to file
    const double prob_dup;          // probability of duplicate
    const uint32 n_read_ends;       // (1 for SE Illumina or PacBio, 2 for PE Illumina)
    std::vector<std::string> fastq_chunks;

    ReadWriterOneCore(const T& read_filler_base,
                      const uint32& n_reads_,
                      const uint32& read_chunk_size_,
                      const double& prob_dup_,
                      const uint32& n_read_ends_)
        : read_filler(read_filler_base),
          n_reads(n_reads_),
          read_chunk_size(read_chunk_size_),
          reads_made(0),
          reads_in_chunk(0),
          do_write(false),
          prob_dup(prob_dup_),
          n_read_ends(n_read_ends_),
          fastq_chunks(n_read_ends_, "") {};


    // Write contents in `fastq_chunks` to UNcompressed file(s).
    void write(std::vector<std::ofstream>& files) {
        for (uint32 i = 0; i < fastq_chunks.size(); i++) {
            files[i] << fastq_chunks[i];
            fastq_chunks[i].clear();
        }
        reads_in_chunk = 0;
        do_write = false;
        return;
    }
    // Write contents in `fastq_chunks` to compressed file(s).
    void write(std::vector<gzFile>& files) {
        for (uint32 i = 0; i < fastq_chunks.size(); i++) {
            gzwrite(files[i], fastq_chunks[i].c_str(), fastq_chunks[i].size());
            fastq_chunks[i].clear();
        }
        reads_in_chunk = 0;
        do_write = false;
        return;
    }
    /* Overloaded for one file: */
    void write(std::ofstream& file) {
        file << fastq_chunks[0];
        fastq_chunks[0].clear();
        reads_in_chunk = 0;
        do_write = false;
        return;
    }
    void write(gzFile& file) {
        gzwrite(file, fastq_chunks[0].c_str(), fastq_chunks[0].size());
        fastq_chunks[0].clear();
        reads_in_chunk = 0;
        do_write = false;
        return;
    }

    /*
     Add new read(s) to `fastq_chunks`, and update bool for whether you should
     write to file
     */
    void create_reads(pcg64& eng) {
        read_filler.one_read(fastq_chunks, eng);
        reads_made += n_read_ends;
        reads_in_chunk += n_read_ends;
        double dup = runif_01(eng);
        while (dup < prob_dup && reads_made < n_reads &&
               reads_in_chunk < read_chunk_size) {
            read_filler.re_read(fastq_chunks, eng);
            reads_made += n_read_ends;
            reads_in_chunk += n_read_ends;
            dup = runif_01(eng);
        }
        do_write = reads_in_chunk >= read_chunk_size || reads_made >= n_reads;
        return;
    }



};








//' Split number of reads by number of cores.
//'
//' @noRd
//'
inline std::vector<uint32> split_n_reads(const uint32& n_reads,
                                         const uint32& n_cores) {
    std::vector<uint32> out(n_cores, n_reads / n_cores);
    uint32 sum_reads = std::accumulate(out.begin(), out.end(), 0U);
    uint i = 0;
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
 Create and open files:
 */

// Uncompressed version:
inline void open_fastq(std::vector<std::ofstream>& files,
                             const std::string& out_prefix) {

    for (uint32 i = 0; i < files.size(); i++) {

        std::string file_name = out_prefix + "_R" + std::to_string(i+1)+ ".fq";

        files[i].open(file_name);

        if (!files[i].is_open()) {
            std::string e = "Unable to open file " + file_name + ".\n";
            Rcpp::stop(e);
        }

    }

    return;

}
// Compressed version:
inline void open_fastq(std::vector<gzFile>& files,
                             const std::string& out_prefix) {

    for (uint32 i = 0; i < files.size(); i++) {

        std::string file_name = out_prefix + "_R" + std::to_string(i+1)+ ".fq.gz";

        /*
         Initialize file.
         Note that gzfile does not tolerate initializing an empty file.
         Use ofstream instead.
         */
        if (!std::ifstream(file_name)){
            std::ofstream myfile;
            myfile.open(file_name, std::ios::out | std::ios::binary);
            myfile.close();
        }

        files[i] = gzopen(file_name.c_str(), "wb");
        if (!files[i]) {
            std::string e = "gzopen of " + file_name + " failed: " +
                strerror(errno) + ".\n";
            Rcpp::stop(e);
        }

    }

    return;
}
/*
 Close FASTQ files
 */
// Uncompressed version:
inline void close_fastq(std::vector<std::ofstream>& files) {
    for (uint32 i = 0; i < files.size(); i++) files[i].close();
    return;
}
// Compressed version:
inline void close_fastq(std::vector<gzFile>& files) {
    for (uint32 i = 0; i < files.size(); i++) gzclose(files[i]);
    return;
}





/*
 Make Illumina reads and write them to file(s).

 `T` should be `[Illumina|PacBio]Reference` or `[Illumina|PacBio]Variants`.
 `U` should be `gzFile` or `std::ofstream`.

 */
template <typename T, typename U>
void write_reads_cpp_(const T& read_filler_base,
                      const std::string& out_prefix,
                      const uint32& n_reads,
                      const double& prob_dup,
                      const uint32& read_chunk_size,
                      const uint32& n_read_ends,
                      uint32 n_cores,
                      const bool& show_progress) {


    // To make sure reads_per_core is still accurate if OpenMP not used:
#ifndef _OPENMP
    n_cores = 1;
#endif

    const std::vector<uint32> reads_per_core = split_n_reads(n_reads, n_cores);

    // Generate seeds for random number generators (1 RNG per core)
    const std::vector<std::vector<uint64>> seeds = mc_seeds(n_cores);

    /*
    Create and open files:
    */
    std::vector<U> files(n_read_ends);
    open_fastq(files, out_prefix);

    // Progress bar
    Progress prog_bar(n_reads, show_progress);


#ifdef _OPENMP
#pragma omp parallel num_threads(n_cores) if (n_cores > 1)
{
#endif

    std::vector<uint64> active_seeds;

    // Write the active seed per core or just write one of the seeds.
#ifdef _OPENMP
    uint32 active_thread = omp_get_thread_num();
#else
    uint32 active_thread = 0;
#endif
    active_seeds = seeds[active_thread];

    pcg64 eng = seeded_pcg(active_seeds);

    uint32 reads_this_core = reads_per_core[active_thread];

    ReadWriterOneCore<T> writer(read_filler_base, reads_this_core,
                                read_chunk_size, prob_dup, n_read_ends);

    uint32 reads_written;

    while (writer.reads_made < reads_this_core) {

        // Every 10 reads, check that the user hasn't interrupted the process
        if (writer.reads_made % 10 == 0 && prog_bar.is_aborted()) break;

        writer.create_reads(eng);

        if (writer.do_write) {
#ifdef _OPENMP
#pragma omp critical
#endif
            reads_written = writer.reads_in_chunk;
            writer.write(files);
            prog_bar.increment(reads_written);
        }
    }

#ifdef _OPENMP
}
#endif


    // Close files
    close_fastq(files);

    return;
};







#endif
