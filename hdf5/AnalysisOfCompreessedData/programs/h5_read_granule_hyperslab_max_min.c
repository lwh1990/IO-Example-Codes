#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h>

#define DATASETNAME "/All_Data/CrIS-SDR_All/ES_ImaginaryLW"

#define difference(a, b, result)\
        do {                            \
        (result).tv_sec = (a).tv_sec - (b).tv_sec; \
        (result).tv_usec = (a).tv_usec - (b).tv_usec;  \
        if ((result).tv_usec < 0) {  \
                --(result).tv_sec;   \
                (result).tv_usec += 1000000; \
                } \
        } while (0)


int
main (int argc, char *argv[]) {

    hid_t       file;
    hid_t       fapl;
    hid_t       dset;
    hid_t       space, dataspace;
    char       *inputfile;
    char       *dsetname;
    int         ndims;
    int         i, j, k, l, m, status;
    hsize_t     dims[4];
    hsize_t     count[4];
    hssize_t    offset[4];
    float       rdata[4][30][9][163];
//    float       rdata[4][30][9][717];
    float       max = 0.0, min = 0.0;

    struct timeval start, finish;
    struct timeval startw, finishw;
    struct timeval startw_ch, finishw_ch;
    struct timeval diff, diffw, diffw_ch;
    double total, totalw, totalw_ch, min_totalw_ch;

/* We can get all these dimensions from the dataset in the file, but for now I'm setting them here. */
    dims[0] = 60;
    dims[1] = 30;
    dims[2] = 9;
    dims[3] = 163;
//    dims[3] = 717;
    count[0] = 4;
    count[1] = 30;
    count[2] = 9;
    count[3] = 163;
//    count[3] = 717;
    offset[0] = 0;
    offset[1] = 0;
    offset[2] = 0;
    offset[3] = 0;

    /*
     * Open file and dataset using the default properties.
     */

    dsetname = argv[1];
    
    inputfile = argv[2];

    gettimeofday(&startw, NULL);

    fapl = H5Pcreate (H5P_FILE_ACCESS);
 
    file = H5Fopen (inputfile, H5F_ACC_RDONLY, fapl);

    if (file < 0) {
        fprintf(stderr, "Error:  Failed to open file %s.\n", inputfile);     
    } else {
        printf("Opened input file %s.\n", inputfile);
    }
 
    dset = H5Dopen (file, dsetname,  H5P_DEFAULT);

    if (dset < 0) {
        fprintf(stderr, "Error:  Failed to open dataset %s.\n", dsetname);
    } else {
        printf("Opened dataset %s.\n", dsetname);
    }

    dataspace = H5Screate_simple(4, count, NULL);

    space = H5Dget_space (dset);

    for (i=0; i<(dims[0]/count[0]); i++) {
        offset[0] = i* count[0];

        status = H5Sselect_hyperslab (space, H5S_SELECT_SET,
                                      offset, NULL, count, NULL);

        status = H5Dread (dset, H5T_IEEE_F32LE, dataspace, space, H5P_DEFAULT,
                          (void *) rdata);

        for (j = 0; j<4; ++j) {
        for (k=0; k<dims[1]; ++k) {
            for (l=0; l<dims[2]; ++l) {
                for (m=0; m<dims[3]; ++m) {
                 
                    if (rdata[j][k][l][m]) {
                        if (max < rdata[j][k][l][m]) 
                            max = rdata[j][k][l][m];
                        if (min > rdata[j][k][l][m])
                            min = rdata[j][k][l][m];
                    }
                }
            }
            
        }
        }
        if (i%5 == 0) {
            printf("Hyperslab %d value (0,0,0,0):  %f\n", i, rdata[0][0][0][0]);
        }

    
    } // i-loop

    gettimeofday(&finishw, NULL);
 
    printf("\nMax and min for %s:\n", dsetname);
    printf("    Max:  %f\n", max);
    printf("    Min:  %f\n\n", min);

    status = H5Sclose (dataspace);
    status = H5Sclose (space);
    H5Dclose(dset);
    H5Pclose(fapl);
    H5Fclose(file);

    difference(finishw, startw, diffw);
    totalw = (double)diffw.tv_sec+(double)diffw.tv_usec/(double)1000000.0;
    printf("\nRead time (by granule hyperslab):  %f sec\n", totalw);

}
