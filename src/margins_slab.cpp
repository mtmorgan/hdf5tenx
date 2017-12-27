#include <string>
#include <Rcpp.h>
#include "c++/H5Cpp.h"

#include "margin.h"
#include "slab.h"

using namespace Rcpp;
using namespace H5;

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

std::vector<int64_t> tenx_indptr(H5::H5File h5, const H5std_string indptr)
{
    hsize_t n;
    H5::DataSet dataset = h5.openDataSet(indptr);
    dataset.getSpace().getSimpleExtentDims(&n, NULL);
    std::vector<int64_t> result(n);
    dataset.read(&result[0], H5::PredType::NATIVE_LONG);
    return result;
}

// [[Rcpp::export]]
List margins_slab(
    const std::string fname, const std::string group,
    const int offset, const int count
    )
{
    const H5std_string
        h5_name( fname ),
        group_name( group ),
        indptr_name( group_name + "/indptr" ),
        indices_name( group_name + "/indices" ),
        data_name( group_name + "/data" );

    H5File h5( h5_name, H5F_ACC_RDONLY );
    std::vector<int> dim = tenx_dim( h5, group_name );
    std::vector<int64_t> indptr = tenx_indptr( h5, indptr_name );

    // indices + data
    slab indices( h5, indices_name ), data( h5, data_name );
    std::vector<int64_t> indices_v( count ), data_v( count );
    indices.read( indices_v, count, indptr[offset] );
    data.read( data_v, count, indptr[offset] );

    // summarize
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

    return List::create(gene.as_list(), cell.as_list());
}
