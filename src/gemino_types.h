# ifndef __GEMINO_TYPES_H
# define __GEMINO_TYPES_H


#include <RcppArmadillo.h>
#include <cstdint>
#include <vector>
#include <string>
#include <deque>


using namespace Rcpp;

typedef uint_fast8_t uint8;
typedef uint_fast32_t uint;
typedef int_fast32_t sint;
typedef uint_fast64_t uint64;
typedef int_fast64_t sint64;



struct SequenceSet {

    // Member variables
    uint64 total_size;
    std::vector<uint> seq_sizes;
    std::vector<std::string> seq_names;
    std::vector<std::string> sequences;
    bool merged;

    // Constructors
    SequenceSet() {
        total_size = 0;
        seq_sizes = std::vector<uint>(0);
        seq_names = std::vector<std::string>(0);
        sequences = std::vector<std::string>(0);
        merged = false;
    }
};

// Struct for one variant
struct OneVariant {

    // Member variables
    std::vector< std::vector<char> > nucleos;
    std::vector< std::vector<uint> > sites;
    std::vector< uint > scaffold_lengths;

    // Constructors
    OneVariant()
    {
        nucleos = std::vector< std::vector<char> >(0, std::vector<char>(0));
        sites = std::vector< std::vector<uint> >(0, std::vector<uint>(0));
        scaffold_lengths = std::vector< uint >(0);
    }
    OneVariant(uint n_scaffs)
    {
        nucleos = std::vector< std::vector<char> >(n_scaffs, std::vector<char>(0));
        sites = std::vector< std::vector<uint> >(n_scaffs, std::vector<uint>(0));
        scaffold_lengths = std::vector< uint >(n_scaffs, 0);
    }
    OneVariant(std::vector<uint> scaff_lens)
    {
        nucleos = std::vector< std::vector<char> >(scaff_lens.size(), std::vector<char>(0));
        sites = std::vector< std::vector<uint> >(scaff_lens.size(), std::vector<uint>(0));
        scaffold_lengths = scaff_lens;
    }
};



// Struct for multiple variants
struct VariantSet {

    // Member variables
    int n_variants;
    uint total_segr_sites;
    std::vector< OneVariant > variant_info;

    // Constructors
    VariantSet()
    {
        n_variants = 0;
        total_segr_sites = 0;
        variant_info = std::vector<OneVariant>(0, OneVariant());
    }
    VariantSet(int n_variants_)
    {
        n_variants = n_variants_;
        total_segr_sites = 0;
        variant_info = std::vector<OneVariant>(n_variants_, OneVariant());
    }
    VariantSet(int n_variants_, uint n_scaffs)
    {
        n_variants = n_variants_;
        total_segr_sites = 0;
        variant_info = std::vector<OneVariant>(n_variants_, OneVariant(n_scaffs));
    }
    VariantSet(int n_variants_, uint n_scaffs, uint total_segr_sites_)
    {
        n_variants = n_variants_;
        total_segr_sites = total_segr_sites_;
        variant_info = std::vector<OneVariant>(n_variants_, OneVariant(n_scaffs));
    }
    VariantSet(int n_variants_, uint total_segr_sites_, std::vector<uint> scaff_lens)
    {
        n_variants = n_variants_;
        total_segr_sites = total_segr_sites_;
        variant_info = std::vector<OneVariant>(n_variants_, OneVariant(scaff_lens));
    }
};






#endif
