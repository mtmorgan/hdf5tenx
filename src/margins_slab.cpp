#include <string>
#include <Rcpp.h>
#include "c++/H5Cpp.h"

#include "margin.h"

std::vector<hsize_t> get_dim( H5::H5File h5, const H5std_string group )
{
    H5::DataSet dataset = h5.openDataSet( group );
    H5::DataSpace dataspace = dataset.getSpace();
    const int rank = dataspace.getSimpleExtentNdims();
    std::vector<hsize_t> dim( rank );
    dataspace.getSimpleExtentDims( &dim[0], NULL );
    return dim;
}

// [[Rcpp::export]]
Rcpp::IntegerVector margins_dim( std::string fname, std::string group )
{
    H5::H5File h5( fname, H5F_ACC_RDONLY );
    std::vector<hsize_t> dim = get_dim( h5, group );

    Rcpp::IntegerVector result( dim.size() );
    std::reverse_copy( dim.begin(), dim.end(), result.begin() );

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

Rcpp::List as_result(const int begin, const margin gene, const margin cell)
{
    Rcpp::Environment env = Rcpp::new_env();
    char id[20];
    std::sprintf(id, "%d", begin);
    env[ id ] = cell.as_list();
    return Rcpp::List::create(
        Rcpp::_["row"] = gene.as_list(),
        Rcpp::_["column"] = env
        );
}

// [[Rcpp::export]]
Rcpp::List margins_rle_slab(
    const std::string fname, const std::string group,
    const int bufsize, const int begin, const int end
    )
{
    H5::H5File h5( fname, H5F_ACC_RDONLY );
    const H5std_string
        indptr_path( group + "/indptr" ), genes_path( group + "/genes" ),
        indices_path( group + "/indices"), data_path( group + "/data" );
    const hsize_t
        nrow = get_dim( h5, genes_path )[0], ncol = end - begin;
    hsize_t rank = 1, count, start;
    std::vector<int64_t>
        indptr( ncol + 1 ), indices( bufsize ), data( bufsize );
    margin gene( nrow ), cell( ncol );

    count = ncol + 1; start = begin;
    slab_read( indptr, h5, indptr_path, rank, &count, &start );

    count = bufsize; start = indptr[0];
    const hsize_t nelt = indptr[indptr.size() - 1];
    int i, j = -1, l = 0;
    while ( nelt - start > 0 ) {
        count = std::min( count, nelt - start);
        slab_read( indices, h5, indices_path, rank, &count, &start );
        slab_read( data, h5, data_path, rank, &count,  &start );
        for (int k = 0; k < count; ++k) {
            i = indices[k];
            if (indptr[0] + l == indptr[j + 1])
                ++j;
            ++l;
            gene.update(i, data[k]);
            cell.update(j, data[k]);
        }

        start += count;
    }

    return as_result( begin, gene, cell );
}

// [[Rcpp::export]]
Rcpp::List margins_dense_slab(
    const std::string fname, const std::string group,
    const int bufsize, const int begin, const int end
    )
{
    H5::H5File h5( fname, H5F_ACC_RDONLY );
    // transposed!
    const hsize_t nrow = get_dim( h5, group )[1], ncol = end - begin;
    hsize_t rank = 2, count[] = { bufsize, nrow }, start[] = { begin, 0 };
    std::vector<int64_t> data( bufsize * nrow );
    margin gene( nrow ), cell( ncol );

    int l = 0;
    while ( end - start[0] > 0 ) {
        count[0] = std::min( count[0], end - start[0] );
        slab_read( data, h5, group, rank, count, start );
        for (int  j = 0; j < count[0]; ++j) {
            for (int i = 0; i < nrow; ++i) {
                const double d = data[j * nrow + i];
                gene.update(i, d);
                cell.update(l + j, d);
            }
        }
        start[0] += count[0];
        l += count[0];
    }
    
    return as_result( begin, gene, cell );
}
