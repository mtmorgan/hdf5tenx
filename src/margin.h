#ifndef MARGIN_H
#define MARGIN_H

#include <Rcpp.h>
#include <vector>

class margin {

public:

    int offset;
    std::vector<int> n;
    std::vector<double> sum, sumsq;

  margin( int offset, hsize_t n) : offset(offset), n(n), sum(n), sumsq(n) {};

    inline void update(int i, double d) {
        n[i] += 1;
        sum[i] += d;
        sumsq[i] += d * d;
    }

    Rcpp::List as_list() {
        Rcpp::IntegerVector r_offset(1);
        r_offset[0] = offset;
        return Rcpp::List::create(
            r_offset,
            Rcpp::IntegerVector(n.begin(), n.end()),
            Rcpp::NumericVector(sum.begin(), sum.end()),
            Rcpp::NumericVector(sumsq.begin(), sumsq.end())
            );
    }

};

#endif
