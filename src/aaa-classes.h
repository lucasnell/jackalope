# ifndef __GEMINO_DEFUNT_CLASSES_H
# define __GEMINO_DEFUNT_CLASSES_H

#include <RcppArmadillo.h>
#include <fstream> // fstream
#include <string> // string class
#include <algorithm> // sort
#include <vector> // vector class

// #include <boost/iostreams/filtering_stream.hpp>
// #include <boost/iostreams/filter/gzip.hpp>

#include "gemino_types.h" // integer types, VariantSet, SequenceSet
#include "str_manip.h" // cpp_to_upper


using namespace Rcpp;





XPtr<SequenceSet> SequenceSet_characters(std::vector<std::string>& seq_names,
                                         std::vector<std::string>& sequences,
                                         bool merged);

uint64 total_size_SequenceSet(XPtr<SequenceSet> ss);

std::vector<uint> seq_sizes_SequenceSet(XPtr<SequenceSet> ss);

char one_nucleo_SequenceSet(XPtr<SequenceSet> ss, int scaff, int pos);

std::string one_scaff_SequenceSet(XPtr<SequenceSet> ss, int scaff);

XPtr<SequenceSet> SequenceSet_fasta(std::string file_name, bool cut_names,
                                    bool remove_soft_mask);

// XPtr<SequenceSet> SequenceSet_fastagz(std::string file_name, bool cut_names,
//                                       bool remove_soft_mask);

void SummarizeSequenceSet(XPtr<SequenceSet> ss, int console_width);






/*
 ---------------------
 ---------------------

 For variants objects

 ---------------------
 ---------------------
 */





XPtr<VariantSet> VariantSet_vectors(const std::vector< std::vector<std::string> >& nucleos,
                                    const std::vector< std::vector< std::vector<uint> > >& sites,
                                    const std::vector< std::vector<uint> >& scaffold_lengths);

int n_variants_VS(const XPtr<VariantSet>& vs);

uint total_segr_sites_VS(const XPtr<VariantSet>& vs);

std::vector<char> nucleos_VS(const XPtr<VariantSet>& vs, uint variant_index, uint scaff_index);

std::vector<uint> sites_VS(const XPtr<VariantSet>& vs, uint variant_index, uint scaff_index);

uint scaffold_length_VS(XPtr<VariantSet> vs, uint variant_index, uint scaff_index);


std::string variants_retr_scaff(uint scaff_num, uint variant_num,
                           const XPtr<SequenceSet>& ss,
                           const XPtr<VariantSet>& vs);

std::string cpp_retr_scaff(uint scaff_index, uint variant_index,
                      const XPtr<SequenceSet>& ss,
                      const XPtr<VariantSet>& vs);

std::string variants_retr_seq(const size_t& start_pos, const size_t& length_out,
                         uint scaff_num, uint variant_num,
                         const XPtr<SequenceSet>& ss,
                         const XPtr<VariantSet>& vs);

std::vector<std::string> variants_retr_var(const uint& variant_num,
                                 const XPtr<SequenceSet>& ss,
                                 const XPtr<VariantSet>& vs);


# endif
