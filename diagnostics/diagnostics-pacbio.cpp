#include <Rcpp.h>
#include <string>
#include <vector>

//[[Rcpp::plugins(cpp11)]]
using namespace Rcpp;



//' Mean quality for each read.
//'
//'
//[[Rcpp::export]]
std::vector<double> mean_quals(const std::vector<std::string>& quals) {

    size_t n_reads = quals.size();
    std::vector<double> output;
    output.reserve(n_reads);
    int qual_start = static_cast<int>('!');

    for (const std::string& q : quals) {
        double sum = 0;
        double n = static_cast<double>(q.size());
        for (const char& c : q) {
            sum += static_cast<double>(static_cast<int>(c) - qual_start);
        }
        output.push_back(sum / n);
    }

    return output;
}
