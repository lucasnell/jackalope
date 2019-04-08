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
*/

template <typename T>
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
    void write(std::vector<std::ofstream>& files) {
        for (uint32 i = 0; i < fastq_pools.size(); i++) {
            files[i].write(fastq_pools[i].data(), fastq_pools[i].size());
            fastq_pools[i].clear();
        }
        reads_in_pool = 0;
        do_write = false;
        return;
    }
    // Write contents in `fastq_pools` to compressed file(s).
    void write(std::vector<gzFile>& files) {
        for (uint32 i = 0; i < fastq_pools.size(); i++) {
            gzwrite(files[i], fastq_pools[i].data(), fastq_pools[i].size());
            fastq_pools[i].clear();
        }
        reads_in_pool = 0;
        do_write = false;
        return;
    }
    // // Write contents in `fastq_pools` to compressed file(s), thread-safe method.
    // void write(std::vector<gzlog *> files) {
    //     for (uint32 i = 0; i < fastq_pools.size(); i++) {
    //         // char *cstr = new char[fastq_pools[i].size() + 1];
    //         // strcpy(cstr, fastq_pools[i].c_str());
    //         // gzlog_write(files[i], cstr, fastq_pools[i].size());
    //         // delete [] cstr;
    //         gzlog_write(files[i], reinterpret_cast<char*>(fastq_pools[i].data()),
    //                     fastq_pools[i].size());
    //         fastq_pools[i].clear();
    //     }
    //     reads_in_pool = 0;
    //     do_write = false;
    //     return;
    // }

    /* Overloaded for one file: */
    void write(std::ofstream& file) {
        file.write(fastq_pools[0].data(), fastq_pools[0].size());
        fastq_pools[0].clear();
        reads_in_pool = 0;
        do_write = false;
        return;
    }
    void write(gzFile& file) {
        gzwrite(file, fastq_pools[0].data(), fastq_pools[0].size());
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
 Test opening files (meant to be done NOT in parallel):
 */
template<typename W>
inline void test_open_fastq(const uint32& n_read_ends,
                            const std::string& out_prefix) {
    return;
}

template<>
inline void test_open_fastq<std::ofstream>(const uint32& n_read_ends,
                                           const std::string& out_prefix) {

    for (uint32 i = 0; i < n_read_ends; i++) {

        std::string file_name = out_prefix + "_R" + std::to_string(i+1)+ ".fq";

        std::ofstream out_file(file_name);

        if (!out_file.is_open()) {
            std::string e = "Unable to open file " + file_name + ".\n";
            Rcpp::stop(e);
        }

    }

    return;

}

template<>
inline void test_open_fastq<gzFile>(const uint32& n_read_ends,
                                    const std::string& out_prefix) {

    for (uint32 i = 0; i < n_read_ends; i++) {

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

        gzFile file = gzopen(file_name.c_str(), "wb");
        if (!file) {
            str_stop({"gzopen of ", file_name, " failed: ", strerror(errno), ".\n"});
        }

    }

    return;
}




/*
 Create and open files:
 */

// Uncompressed version:
inline void open_fastq(std::vector<std::ofstream>& files,
                             const std::string& out_prefix) {

    for (uint32 i = 0; i < files.size(); i++) {

        std::string file_name = out_prefix + "_R" + std::to_string(i+1)+ ".fq";

        files[i].open(file_name, std::ofstream::out | std::ofstream::app);

    }

    return;

}
// Compressed version:
inline void open_fastq(std::vector<gzFile>& files,
                             const std::string& out_prefix) {

    for (uint32 i = 0; i < files.size(); i++) {

        std::string file_name = out_prefix + "_R" + std::to_string(i+1)+ ".fq.gz";

        files[i] = gzopen(file_name.c_str(), "ab");

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
 `W` should be `gzFile` or `std::ofstream`.

 */
template <typename T, typename W>
void write_reads_cpp_(const T& read_filler_base,
                      const std::string& out_prefix,
                      const uint32& n_reads,
                      const double& prob_dup,
                      const uint32& read_pool_size,
                      const uint32& n_read_ends,
                      uint32 n_threads,
                      const bool& show_progress) {


    // To make sure reads_per_thread is still accurate if OpenMP not used:
#ifndef _OPENMP
    n_threads = 1;
#endif

    const std::vector<uint32> reads_per_thread = split_n_reads(n_reads, n_threads);

    // Generate seeds for random number generators (1 RNG per thread)
    const std::vector<std::vector<uint64>> seeds = mt_seeds(n_threads);

    // Test that files can be opened before going into parallel chunk:
    test_open_fastq<W>(n_read_ends, out_prefix);

    // Progress bar
    Progress prog_bar(n_reads, show_progress);

    // Create and open files:
    std::vector<W> files(n_read_ends);
    open_fastq(files, out_prefix);


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

    ReadWriterOneThread<T> writer(read_filler_base, reads_this_thread,
                                  read_pool_size, prob_dup, n_read_ends);

    uint32 reads_written;

    while (writer.reads_made < reads_this_thread) {

        // Every 10 reads, check that the user hasn't interrupted the process
        if (writer.reads_made % 10 == 0 && prog_bar.is_aborted()) break;

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
        }
    }

#ifdef _OPENMP
}
#endif

    // Close files
    close_fastq(files);

    return;
};








// /*
// Make Illumina reads and write them to file(s).
//
// `T` should be `[Illumina|PacBio]Reference` or `[Illumina|PacBio]Variants`.
//
// */
// template <typename T>
// void write_reads_cpp_gz_(const T& read_filler_base,
//                          const std::string& out_prefix,
//                          const uint32& n_reads,
//                          const double& prob_dup,
//                          const uint32& read_pool_size,
//                          const uint32& n_read_ends,
//                          uint32 n_threads,
//                          const bool& show_progress) {
//
//
//     // To make sure reads_per_thread is still accurate if OpenMP not used:
// #ifndef _OPENMP
//     stop("\nwrite_reads_cpp_gz_ should never be run when OpenMP is not enabled");
// #endif
//
//     if (n_threads <= 1) stop("\nwrite_reads_cpp_gz_ should never be run with 1 thread");
//
//     const std::vector<uint32> reads_per_thread = split_n_reads(n_reads, n_threads);
//
//     // Allow for pool size to get 10x as large as specified in case of conflicts
//     // regarding access to output file
//     const uint32 max_read_pool_size = 10 * read_pool_size;
//
//     // Generate seeds for random number generators (1 RNG per thread)
//     const std::vector<std::vector<uint64>> seeds = mt_seeds(n_threads);
//
//     // Test that files can be opened before going into parallel chunk:
//     test_open_fastq<gzFile>(n_read_ends, out_prefix);
//
//     // Now create gzlog* object(s)
//     std::vector<gzlog *> files(n_read_ends);
//     for (uint32 i = 0; i < n_read_ends; i++) {
//         // (gzlog_open adds the ".gz")
//         std::string path_str = out_prefix + "_R" + std::to_string(i+1)+ ".fq";
//         char *path = new char[path_str.length() + 1];
//         // (gzlog class will delete pointer after you close the file.)
//         strcpy(path, path_str.c_str());
//         files[i] = gzlog_open(path);
//         if (files[i] == NULL) {
//             str_stop({"\nThere was an error creating the file ", path_str + ".gz"});
//         }
//     }
//
//     // Progress bar
//     Progress prog_bar(n_reads, show_progress);
//
//     // Whether writing is currently happening:
//     bool writing = false;
//
//
// #ifdef _OPENMP
// #pragma omp parallel num_threads(n_threads) shared(files, writing) if (n_threads > 1)
// {
// #endif
//
// #ifdef _OPENMP
//     uint32 active_thread = omp_get_thread_num();
// #else
//     uint32 active_thread = 0;
// #endif
//
//     std::vector<uint64> active_seeds;
//
//     // Write the active seed per thread or just write one of the seeds.
//     active_seeds = seeds[active_thread];
//
//     pcg64 eng = seeded_pcg(active_seeds);
//
//     uint32 reads_this_thread = reads_per_thread[active_thread];
//
//     ReadWriterOneThread<T> writer(read_filler_base, reads_this_thread,
//                                   read_pool_size, prob_dup, n_read_ends);
//
//     uint32 reads_written;
//
//     while (writer.reads_made < reads_this_thread) {
//
//         // Every 10 reads, check that the user hasn't interrupted the process
//         if (writer.reads_made % 10 == 0 && prog_bar.is_aborted()) break;
//
//         writer.create_reads(eng);
//
//         if (writer.reads_in_pool > max_read_pool_size ||
//             writer.reads_made >= reads_this_thread ||
//             (writer.do_write && !writing)) {
// #ifdef _OPENMP
// #pragma omp critical
// {
// #endif
//             writing = true;
//             reads_written = writer.reads_in_pool;
//             writer.write(files);
//             prog_bar.increment(reads_written);
//             writing = false;
// #ifdef _OPENMP
// }
// #endif
//         }
//     }
//
// #ifdef _OPENMP
// }
// #endif
//
//
//
//     // Finish compressing, close files, and print warnings if bad things happened
//     std::vector<int> out_codes_compr(n_read_ends);
//     std::vector<int> out_codes_close(n_read_ends);
//     for (uint32 i = 0; i < n_read_ends; i++) {
//         out_codes_compr[i] = gzlog_compress(files[i]);
//         out_codes_close[i] = gzlog_close(files[i]);
//     }
//     for (uint32 i = 0; i < n_read_ends; i++) {
//         std::string fn = out_prefix + "_R" + std::to_string(i+1)+ ".fq.gz";
//
//         if (out_codes_compr[i] == -1) {
//             str_warn({"\nWhen compressing the file ", fn,", there was a file I/O error. ",
//                      "This can happen when the device has run out of space or when ",
//                      "there are permissions / ownership issues."});
//         }
//         if (out_codes_compr[i] == -2) {
//             str_warn({"\nWhen compressing the file ", fn,", there was a memory ",
//                      "allocation failure. You may want to check the ",
//                      "file for nonsense."});
//         }
//         if (out_codes_compr[i] == 3 || out_codes_close[i] == -3) {
//             str_warn({"\nWhen closing the file ", fn,", the C++ object ",
//                      "used to create it was invalid. You may want to check the ",
//                      "file for nonsense."});
//         }
//     }
//
//
//     return;
// };






#endif
