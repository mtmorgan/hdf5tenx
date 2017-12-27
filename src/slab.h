#ifndef SLAB_H
#define SLAB_H

#include <vector>
#include "c++/H5Cpp.h"

class slab {
private:
    H5::DataSet dataset;
    H5::DataSpace dataspace;

public:
    slab(H5::H5File& h5, H5std_string name) {
        dataset = h5.openDataSet( name );
        dataspace = dataset.getSpace();
    }

    void read(
        std::vector<int64_t>& data, hsize_t c_count, hsize_t c_start
    ) {
        hsize_t count[] = { c_count }, start[] = { c_start };
        H5::DataSpace memspace( 1, count );
        dataspace.selectHyperslab( H5S_SELECT_SET, count, start );
        dataset.read(&data[0], H5::PredType::NATIVE_LONG, memspace, dataspace);
    }

};

#endif
