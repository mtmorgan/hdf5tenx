#ifndef MARGIN_H
#define MARGIN_H

#include <Rcpp.h>
#include <vector>

class margin {

public:

    std::vector<int> n;
    std::vector<double> sum, sumsq;

    margin(hsize_t n) : n(n), sum(n), sumsq(n) {};

    inline void update(int i, double d) {
        n[i] += 1;
        sum[i] += d;
        sumsq[i] += d * d;
    }

    Rcpp::List as_list() {
        return Rcpp::List::create(
            Rcpp::IntegerVector(n.begin(), n.end()),
            Rcpp::NumericVector(sum.begin(), sum.end()),
            Rcpp::NumericVector(sumsq.begin(), sumsq.end())
            );
    }

};

#endif
