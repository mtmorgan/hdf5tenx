#include <string>
#include <Rcpp.h>
#include "c++/H5Cpp.h"

#include "margin.h"

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

std::vector<hsize_t> margins_dim_internal(
    H5::H5File h5, const H5std_string group
    )
{
    const H5std_string
        barcodes( group + "/barcodes"),
        genes( group + "/genes");
    std::vector<hsize_t> dim(2);

    h5.openDataSet(genes).getSpace().getSimpleExtentDims( &dim[0], NULL);
    h5.openDataSet(barcodes).getSpace().getSimpleExtentDims( &dim[1], NULL);

    return dim;
}

// [[Rcpp::export]]
Rcpp::IntegerVector margins_dim( std::string fname, std::string group )
{
    H5::H5File h5( fname, H5F_ACC_RDONLY );
    std::vector<hsize_t> dim = margins_dim_internal( h5, group );

    Rcpp::IntegerVector result( dim.size() );
    std::copy( dim.begin(), dim.end(), result.begin() );
    return result;
}

void slab_read(
    std::vector<int64_t>& data, H5::H5File& h5, const H5std_string group,
    const hsize_t rank, const hsize_t *count, const hsize_t *start
    )
{
    H5::DataSet dataset = h5.openDataSet( group );
    H5::DataSpace dataspace = dataset.getSpace();
    H5::DataSpace memspace( rank, count );

    dataspace.selectHyperslab( H5S_SELECT_SET, count, start );
    dataset.read( &data[0], H5::PredType::NATIVE_LONG, memspace, dataspace );
}

// [[Rcpp::export]]
Rcpp::List margins_slab(
    const std::string fname, const std::string group,
    const std::vector<double> indptr, const int begin, const int end
    )
{
    const H5std_string
        indices_name( group + "/indices" ),
        data_name( group + "/data" );
    const hsize_t
        rank = 1,
        count[] = { indptr[end] - indptr[begin] },
        start[] = { indptr[begin] } ;

    H5::H5File h5( fname, H5F_ACC_RDONLY );

    // indices + data
    std::vector<int64_t> indices( count[0] ), data( count[0] );
    slab_read( indices, h5, indices_name, rank, count, start );
    slab_read( data, h5, data_name, rank, count,  start );

    // summarize
    const std::vector<hsize_t> dim = margins_dim_internal( h5, group );
    margin gene( 0, dim[0] ), cell( begin, end - begin );
    int i, j = 0;
    for (int k = 0; k < count[0]; ++k) {
        i = indices[k];
        if (indptr[begin] + k == indptr[begin + j + 1])
            ++j;
        gene.update(i, data[k]);
        cell.update(j, data[k]);
    }

    return Rcpp::List::create(gene.as_list(), cell.as_list());
}
