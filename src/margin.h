#ifndef MARGIN_H
#define MARGIN_H

#include <Rcpp.h>
#include <vector>

class margin {

public:

    std::vector<int> n;
    std::vector<double> sum, sumsq;

    margin(int n) : n(n), sum(n), sumsq(n) {};

    void reset() {
        for (int i = 0; i < n.size(); ++i) {
            n[i] = sum[i] = sumsq[i] = 0;
        }
    }

    void update(int i, double d) {
        n[i] += 1;
        sum[i] += d;
        sumsq[i] += d * d;
    }

    void join( const margin &rhs ) {
        const int length = n.size();
        for (int i = 0; i < length; ++i) {
            n[i] += rhs.n[i];
            sum[i] += rhs.sum[i];
            sumsq[i] += rhs.sumsq[i];
        }
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
