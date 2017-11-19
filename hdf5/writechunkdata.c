#include <stdio.h>
#include <hdf5.h>
#include <stdlib.h>
#include <sys/time.h>

double Time()
{
    struct timeval  t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec * 1.e-6;
}

int main(void) {
    hid_t   file_id, dset_id, fspace_id,fspace_id1,  mspace_id, dcpl_id;
    herr_t status;
    hsize_t chunk_dims[2] = {10, 1};
    hsize_t dset_dims[2] = {10, 10};
    hsize_t mem_dims[1] = {5};
    hsize_t start1[2] = {3, 2}, start[2] = {0, 0};
    hsize_t count1[2] = {5, 1}, count[2] = {1, 5};
    int     buffer[5], i;
    double t1,t2,t3;


    /* Create the file */
    file_id = H5Fcreate("file.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* Create a dataset creation property list and set it to use chunking
    *      * with a chunk size of 10x1 */
    dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl_id, 2, chunk_dims);

    /* Create the dataspace and the chunked dataset */
    fspace_id = H5Screate_simple(2, dset_dims, NULL);
   // fspace_id1 = H5Screate_simple(2, dset_dims, NULL);
    dset_id = H5Dcreate(file_id, "dataset", H5T_NATIVE_INT, fspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);

    /* Select the elements from 3, 2 to 7, 2 */
    H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, start, NULL, count, NULL);
    //H5Sselect_hyperslab(fspace_id1, H5S_SELECT_SET, start1, NULL, count1, NULL);

    /* Create the memory dataspace */
    mspace_id = H5Screate_simple(1, mem_dims, NULL);

    /*for(i = 0; i < 5 ; i++)
    {
        printf("%d\n", buffer[i]);
    }*/
    /* Write to the dataset */
    t1 = Time();
    status = H5Dwrite(dset_id, H5T_NATIVE_INT, mspace_id, fspace_id, H5P_DEFAULT, buffer);
    //t2 = Time();
    //status = H5Dwrite(dset_id, H5T_NATIVE_INT, mspace_id, fspace_id1, H5P_DEFAULT, buffer);
    t3 = Time();
    printf("write time = %f \n", t3 - t1);

    /* Close */
    H5Dclose(dset_id);
    H5Sclose(fspace_id);
    H5Sclose(mspace_id);
    H5Pclose(dcpl_id);
    H5Fclose(file_id);
    return 0;
}

