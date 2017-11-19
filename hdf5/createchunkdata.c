#include <stdio.h>
#include <hdf5.h>

int main(void) {
    hid_t   file_id, dset_id, space_id, dcpl_id;
    herr_t status;
    hsize_t chunk_dims[2] = {4, 4};
    hsize_t dset_dims[2] = {12, 12};
    int     buffer[12][12],i,j;

    /* Create the file */
    file_id = H5Fcreate("file.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* Create a dataset creation property list and set it to use chunking
    *      */
    dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl_id, 2, chunk_dims);

    /* Create the dataspace and the chunked dataset */
    space_id = H5Screate_simple(2, dset_dims, NULL);
    dset_id = H5Dcreate(file_id, "dataset", H5T_NATIVE_INT, space_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);

    /* Write to the dataset */
    for(i = 0; i < 12; i++)
        for(j = 0; j < 12; j++)
        {
            buffer[i][j]=i+j;
        }
    status = H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);

    /* Close */
    H5Dclose(dset_id);
    H5Sclose(space_id);
    H5Pclose(dcpl_id);
    H5Fclose(file_id);
    return 0;
}
