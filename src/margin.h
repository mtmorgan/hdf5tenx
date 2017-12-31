#ifndef MARGIN_H
#define MARGIN_H

#include <Rcpp.h>
#include <vector>

class margin {

public:

    std::vector<int> n;
    std::vector<double> sum, sumsq;

  margin( const int n ) : n(n), sum(n), sumsq(n) {};

    inline void update( const int i, const double d ) {
        if (d == 0)
            return;
        n[i] += 1;
        sum[i] += d;
        sumsq[i] += d * d;
    }

    Rcpp::List as_list() const {
        return Rcpp::List::create(
            Rcpp::_["n"] = Rcpp::IntegerVector(n.begin(), n.end()),
            Rcpp::_["sum"] = Rcpp::NumericVector(sum.begin(), sum.end()),
            Rcpp::_["sumsq"] = Rcpp::NumericVector(sumsq.begin(), sumsq.end())
            );
    }

};

#endif
