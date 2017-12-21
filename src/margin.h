#ifndef MARGIN_H
#define MARGIN_H

#include <Rcpp.h>
#include <vector>

class margin {

private:

    std::vector<int> n;
    std::vector<double> sum, sumsq;

public:

    margin(int n) : n(n), sum(n), sumsq(n) {};

    void update(int i, double d) {
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
