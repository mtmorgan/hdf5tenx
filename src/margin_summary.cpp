#include <string>
#include <Rcpp.h>
#include "c++/H5Cpp.h"

using namespace Rcpp;
using namespace H5;

std::vector<int> tenx_dim(H5File h5, const H5std_string group)
{
    const H5std_string
        barcodes( group + "/barcodes"),
        genes( group + "/genes");
    std::vector<hsize_t> h5dim(2);

    h5.openDataSet(genes).getSpace().getSimpleExtentDims( &h5dim[0], NULL);
    h5.openDataSet(barcodes).getSpace().getSimpleExtentDims( &h5dim[1], NULL);

    std::vector<int> result(2);
    copy(h5dim.begin(), h5dim.end(), result.begin());

    return result;
}

std::vector<int64_t> tenx_indptr(H5File h5, H5std_string indptr)
{
    hsize_t n;
    DataSet dataset = h5.openDataSet(indptr);
    dataset.getSpace().getSimpleExtentDims(&n, NULL);
    std::vector<int64_t> result(n);
    dataset.read(&result[0], PredType::NATIVE_LONG);
    return result;
}

// [[Rcpp::export]]
List tenx_margins(
    CharacterVector r_fname, CharacterVector r_group, IntegerVector r_bufsize
    )
{
    const H5std_string
        h5_name( Rcpp::as<std::string>(r_fname) ),
        group_name( Rcpp::as<std::string>(r_group) ),
        indptr_name( group_name + "/indptr" ),
        indices_name( group_name + "/indices" ),
        data_name( group_name + "/data" );

    H5File h5( h5_name, H5F_ACC_RDONLY );
    std::vector<int> dim = tenx_dim( h5, group_name );
    std::vector<int64_t> indptr = tenx_indptr( h5, indptr_name );

    // indices + data
    DataSet
        indices = h5.openDataSet( indices_name ),
        data = h5.openDataSet( data_name );
    DataSpace
        indices_dataspace = indices.getSpace(),
        data_dataspace = data.getSpace();
    hsize_t
        rank = indices_dataspace.getSimpleExtentNdims(),
        offset[] = { 0 },
        count[] = { Rcpp::as<int>(r_bufsize) };
    hsize_t indices_n;
    indices.getSpace().getSimpleExtentDims( &indices_n, NULL );

    std::vector<int64_t> indices_v( Rcpp::as<int>(r_bufsize) );
    std::vector<int64_t> data_v( Rcpp::as<int>(r_bufsize) );

    hsize_t i, j = -1, k = 0, l = 0;
    std::vector<int> n_i( dim[0] ), n_j( dim[1] );
    std::vector<double>
        sum_i( dim[0] ), sum_j( dim[1] ), sumsq_i( dim[0] ), sumsq_j( dim[1] );
    std::adjacent_difference(indptr.begin() + 1, indptr.end(), n_j.begin());

    while ( count[0] > 0 ) {
        DataSpace memspace( rank, count );

        indices_dataspace.selectHyperslab( H5S_SELECT_SET, count, offset );
        indices.read(
            &indices_v[0], PredType::NATIVE_LONG, memspace, indices_dataspace
            );

        data_dataspace.selectHyperslab( H5S_SELECT_SET, count, offset );
        data.read(
            &data_v[0], PredType::NATIVE_LONG, memspace, data_dataspace
            );

        for (k = 0; k < count[0]; ++k) {
            i = indices_v[k];
            if (offset[0] + k == indptr[l]) {
                ++j;
                ++l;
            }

            const double d = data_v[k], d2 = d * d;
            n_i[i] += 1;
            sum_i[i] += d; sum_j[j] += d;
            sumsq_i[i] += d2; sumsq_j[j] += d2;
        }

        Rcpp::Rcout << j << " / " << dim[1] << std::endl;
        offset[0] += count[0];
        count[0] = std::min( count[0], indices_n - offset[0] );
    }

    List result = List::create(
        List::create(
            IntegerVector(n_i.begin(), n_i.end()),
            NumericVector(sum_i.begin(), sum_i.end()),
            NumericVector(sumsq_i.begin(), sumsq_i.end())
            ),
        List::create(
            IntegerVector(n_j.begin(), n_j.end()),
            NumericVector(sum_j.begin(), sum_j.end()),
            NumericVector(sumsq_j.begin(), sumsq_j.end())
            )
        );

    return result;
}
