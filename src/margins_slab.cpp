#include <string>
#include <Rcpp.h>
#include "c++/H5Cpp.h"

#include "margin.h"
#include "slab.h"

std::vector<int> tenx_dim(H5::H5File h5, const H5std_string group)
{
    const H5std_string
        barcodes( group + "/barcodes"),
        genes( group + "/genes");
    std::vector<hsize_t> h5dim(2);

    h5.openDataSet(genes).getSpace().getSimpleExtentDims( &h5dim[0], NULL);
    h5.openDataSet(barcodes).getSpace().getSimpleExtentDims( &h5dim[1], NULL);

    std::vector<int> result(2);
    std::copy(h5dim.begin(), h5dim.end(), result.begin());

    return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector indptr( const std::string fname, const std::string group )
{
    H5::H5File h5( fname, H5F_ACC_RDONLY );
    H5::DataSet dataset = h5.openDataSet( group + "/indptr" );
    hsize_t n;
    dataset.getSpace().getSimpleExtentDims(&n, NULL);
    std::vector<int64_t> vec(n);
    dataset.read(&vec[0], H5::PredType::NATIVE_LONG);
    Rcpp::NumericVector result(n);
    std::copy(vec.begin(), vec.end(), result.begin());
    return result;
}

// [[Rcpp::export]]
Rcpp::List margins_slab(
    const std::string fname, const std::string group,
    const std::vector<double> indptr, const int offset, const int count
    )
{
    const H5std_string
        indices_name( group + "/indices" ),
        data_name( group + "/data" );

    H5::H5File h5( fname, H5F_ACC_RDONLY );

    // indices + data
    std::vector<int64_t> indices_v( count ), data_v( count );
    slab( h5, indices_name ).read( indices_v, count, indptr[offset] );
    slab( h5, data_name ).read( data_v, count, indptr[offset] );

    // summarize
    std::vector<int> dim = tenx_dim( h5, group );
    margin gene( dim[0] ), cell( dim[1] );
    int l = offset, i, j = l - 1;
    for (int k = 0; k < count; ++k) {
        i = indices_v[k];
        if (indptr[offset] + k == indptr[l]) {
            ++j;
            ++l;
        }
        gene.update(i, data_v[k]);
        cell.update(j, data_v[k]);
    }

    return Rcpp::List::create(gene.as_list(), cell.as_list());
}
