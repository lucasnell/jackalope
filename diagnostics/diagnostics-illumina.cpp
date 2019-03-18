#include <Rcpp.h>
#include <string>
#include <vector>

//[[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

#define NUCLEO_BASES "TCAG"


//' Groups vectors of qualities by position.
//'
//'
//[[Rcpp::export]]
std::vector<std::vector<int>> quals_by_pos(const std::vector<std::string>& input) {

    size_t n_reads = input.size();
    size_t n_pos = input[0].size();
    std::vector<std::vector<int>> output(n_pos);

    for (size_t i = 0; i < n_pos; i++) {
        output[i].reserve(n_reads);
        for (size_t j = 0; j < n_reads; j++) {
            output[i].push_back(static_cast<int>(input[j][i]) - static_cast<int>('!'));
        }
    }

    return output;
}


//' Indicates which nucleotide mostly comprises each read.
//' Useful for distinguishing reverse complements.
//'
//[[Rcpp::export]]
std::vector<int> which_nt(const std::vector<std::string>& reads) {

    size_t n_reads = reads.size();
    size_t n_pos = reads[0].size();
    std::vector<int> output(n_reads);
    std::vector<int> nt_map(256, 0);
    std::string bases = NUCLEO_BASES;
    for (int i = 0; i < bases.size(); i++) nt_map[bases[i]] = i;

    std::vector<int> nt_counts(4);

    for (size_t i = 0; i < n_reads; i++) {
        for (int& c : nt_counts) c = 0;
        for (size_t j = 0; j < n_pos; j++) {
            nt_counts[nt_map[reads[i][j]]]++;
        }
        int max_nt_ind = std::max_element(nt_counts.begin(), nt_counts.end()) -
            nt_counts.begin();
        output[i] = max_nt_ind;
    }

    return output;
}


//[[Rcpp::export]]
std::vector<int> mm_by_qual(const std::vector<std::string>& reads,
                            const std::vector<std::string>& quals,
                            const int& max_qual) {

    size_t n_reads = reads.size();
    size_t n_pos = reads[0].size();
    std::vector<int> output(max_qual, 0);
    std::vector<int> nt_map(256, 0);
    std::string bases = NUCLEO_BASES;
    for (int i = 0; i < bases.size(); i++) nt_map[bases[i]] = i;

    std::vector<int> nt_counts(4);
    std::string max_nts(n_reads,'A');

    for (size_t i = 0; i < n_reads; i++) {
        for (int& c : nt_counts) c = 0;
        for (size_t j = 0; j < n_pos; j++) {
            nt_counts[nt_map[reads[i][j]]]++;
        }
        int max_nt_ind = std::max_element(nt_counts.begin(), nt_counts.end()) -
            nt_counts.begin();
        max_nts[i] = bases[max_nt_ind];
    }

    for (size_t i = 0; i < n_reads; i++) {
        for (size_t j = 0; j < n_pos; j++) {
            if (reads[i][j] != max_nts[i]) {
                int qual = static_cast<int>(quals[i][j]) - static_cast<int>('!');
                output[qual]++;
            }
        }
    }

    return output;
}






//' Count indels when read should only have one two-character repeating motif.
//' This function is for one read only.
//'
void count_indels_one_(const std::string& read,
                       const char& motif1,
                       const char& motif2,
                       int& insertions,
                       int& deletions) {

    if (read.size() < 10) return;

    std::vector<int> output(2, 0);
    auto iter1 = read.begin();
    auto iter2 = read.begin() + 1;

    // Compensate for if read started halfway through motif:
    if (*(iter1+1) == motif1 && *(iter2+1) == motif2) {
        iter1++;
        iter2++;
    }

    while (iter2 < (read.end() - 1)) {

        if (*iter1 != motif1 || *iter2 != motif2) {

            // Insertion before first position
            if (*(iter1+1) == motif1 && *(iter2+1) == motif2) insertions++;
            // Insertion before second position
            if (*iter1 == motif1 && *(iter2+1) == motif2) insertions++;
            // Deletion of first position
            if (*iter1 == motif2 && *iter2 == motif2) deletions++;
            // Deletion of second position
            if (*iter1 == motif2 && *iter2 == motif1) deletions++;
            /*
             Now reset to get to a point where the iterator reaches the first motif
             and second one.
             You'll miss some indels doing this, but if you keep the indel rates pretty
             low, it should be a rare event and good enough for rough diagnostics.
             */
            while (iter2 < read.end() && (*iter1 != motif1 || *iter2 != motif2)) {
                iter1++;
                iter2++;
            }
        } else {
            if (iter2 > (read.end() - 2)) break;
            iter1 += 2;
            iter2 += 2;
        }

    }

    return;
}

//' Count indels when reads should only have one two-character repeating motif.
//' This function is for multiple reads.
//'
//' It outputs first the vector of insertions, then deletions
//'
//'
//[[Rcpp::export]]
std::vector<std::vector<int>> count_indels(const std::vector<std::string>& reads,
                                           const char& motif1,
                                           const char& motif2) {

    std::vector<std::vector<int>> output(2, std::vector<int>(reads.size(), 0));

    for (size_t i = 0; i < reads.size(); i++) {
        count_indels_one_(reads[i], motif1, motif2,
                          output[0][i], output[1][i]);
    }

    return output;

}
