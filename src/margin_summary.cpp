#include <string>
#include <Rcpp.h>
#include <progress.hpp>
#include "c++/H5Cpp.h"

#include "margin.h"

using namespace Rcpp;
using namespace H5;

class slab {
private:
    DataSet dataset;
    DataSpace dataspace;

public:
    slab(H5File& h5, H5std_string name) {
        dataset = h5.openDataSet( name );
        dataspace = dataset.getSpace();
    }

    void read(std::vector<int64_t>& data, hsize_t c_count, hsize_t c_offset) {
        hsize_t count[] = { c_count }, start[] = { c_offset };
        DataSpace memspace( 1, count );
        dataspace.selectHyperslab( H5S_SELECT_SET, count, start );
        dataset.read(&data[0], PredType::NATIVE_LONG, memspace, dataspace);
    }
};

std::vector<int> tenx_dim(H5File h5, const H5std_string group)
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
    indices_dataspace.getSimpleExtentDims( &indices_n, NULL );

    std::vector<int64_t> indices_v( Rcpp::as<int>(r_bufsize) );
    std::vector<int64_t> data_v( Rcpp::as<int>(r_bufsize) );

    int i, j = -1, l = 0;
    margin gene( dim[0] ), cell( dim[1] );

    Progress progress( dim[1] );
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

        for (int k = 0; k < count[0]; ++k) {
            i = indices_v[k];
            if (offset[0] + k == indptr[l]) {
                ++j;
                ++l;
            }

            gene.update(i, data_v[i]);
            cell.update(j, data_v[j]);
        }

        if (progress.check_abort()) {
            return List();
        }
        progress.update( j );

        offset[0] += count[0];
        count[0] = std::min( count[0], indices_n - offset[0] );
    }

    return List::create( gene.as_list(), cell.as_list() );
}

// [[Rcpp::export]]
List tenx_margins_slab(
    CharacterVector r_fname, CharacterVector r_group,
    NumericVector r_offset, NumericVector r_count
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
    const int
        c_count = Rcpp::as<int>(r_count),
        c_offset = Rcpp::as<int>(r_offset) - 1;

    // indices + data
    slab indices( h5, indices_name ), data( h5, data_name );
    std::vector<int64_t> indices_v( c_count ), data_v( c_count );
    indices.read( indices_v, c_count, indptr[c_offset] );
    data.read( data_v, c_count, indptr[c_offset] );

    // summarize
    margin gene( dim[0] ), cell( dim[1] );
    int l = c_offset, i, j = l - 1;
    for (int k = 0; k < c_count; ++k) {
        i = indices_v[k];
        if (indptr[c_offset] + k == indptr[l]) {
            ++j;
            ++l;
        }
        gene.update(i, data_v[k]);
        cell.update(j, data_v[k]);
    }

    return List::create(gene.as_list(), cell.as_list());
}
