#include "c++/H5Cpp.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List hello_world() {

    CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
    NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
    List z            = List::create( x, y ) ;

    return z ;
}
