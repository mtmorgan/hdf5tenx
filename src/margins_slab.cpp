#include <string>
#include <Rcpp.h>
#include "c++/H5Cpp.h"

#include "margin.h"

std::vector<hsize_t> margins_dim(H5::H5File h5, const H5std_string group)
{
    const H5std_string
        barcodes( group + "/barcodes"),
        genes( group + "/genes");
    std::vector<hsize_t> dim(2);

    h5.openDataSet(genes).getSpace().getSimpleExtentDims( &dim[0], NULL);
    h5.openDataSet(barcodes).getSpace().getSimpleExtentDims( &dim[1], NULL);

    return dim;
}

void slab_read(
    std::vector<int64_t>& data, H5::H5File& h5,
    const H5std_string group, const hsize_t c_count, const hsize_t c_start
    )
{
    H5::DataSet dataset = h5.openDataSet( group );
    H5::DataSpace dataspace = dataset.getSpace();
    const hsize_t count[] = { c_count }, start[] = { c_start };
    H5::DataSpace memspace( 1, count );
    dataspace.selectHyperslab( H5S_SELECT_SET, count, start );
    dataset.read( &data[0], H5::PredType::NATIVE_LONG, memspace, dataspace );
}

// [[Rcpp::export]]
Rcpp::NumericVector indptr( const std::string fname, const std::string group )
{
    H5::H5File h5( fname, H5F_ACC_RDONLY );
    H5::DataSet dataset = h5.openDataSet( group + "/indptr" );
    hsize_t n;

    dataset.getSpace().getSimpleExtentDims( &n, NULL );
    std::vector<int64_t> vec(n);
    Rcpp::NumericVector result(n);

    dataset.read( &vec[0], H5::PredType::NATIVE_LONG );
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
    slab_read( indices_v, h5, indices_name, count, indptr[offset] );
    slab_read( data_v, h5, data_name, count, indptr[offset] );

    // summarize
    std::vector<hsize_t> dim = margins_dim( h5, group );
    margin gene( dim[0] ), cell( dim[1] );
    int i, j = offset;
    for (int k = 0; k < count; ++k) {
        i = indices_v[k];
        if (indptr[offset] + k == indptr[j + 1])
            ++j;
        gene.update(i, data_v[k]);
        cell.update(j, data_v[k]);
    }

    return Rcpp::List::create(gene.as_list(), cell.as_list());
}
