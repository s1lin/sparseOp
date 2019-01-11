//
// Created by shilei on 1/8/19.
//

#ifndef SPARSEOP2_SPARSEMATRIX_H
#define SPARSEOP2_SPARSEMATRIX_H

#include <cstdio>
#include <cstdlib>
#include <omp.h>

#include "../../helper/mmio.h"


namespace DataStructure {

    template<class T>
    class SparseMatrix {

        //return code from mmio
        int ret_code;

        //type of the matrix
        MM_typecode matcode;

        //.mtx file name
        const char *mtx_file;

        /*
         * M: num of rows
         * N: num of columns;
         * nz: num of non zeros;
         */
        int M, N, nz;

        /*
         * Lp: the column pointer of L
         * Ij: the index of column pointer
         * Li: the row index of L
         */
        //  and
        int i, *Lp, Ij, *Li;

        //value of L
        T *Lx;


    public:

        //default constructor
        SparseMatrix() {}

        //constructor
        SparseMatrix(const char *mtx) {
            this->mtx_file = mtx;
        }

        //read routine and storage L into CCS format
        void read() {
            //.mtx file
            FILE *f;

            //Check the file:
            if ((f = fopen(mtx_file, "r")) == nullptr)
                exit(1);

            //Go through banner:
            if (mm_read_banner(f, &matcode) != 0) {
                printf("Could not process Matrix Market banner.\n");
                exit(1);
            }

            //Go through matrix type:
            if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
                mm_is_sparse(matcode)) {
                printf("Not support Matrix type: [%s]\n", mm_typecode_to_str(matcode));
                exit(1);
            }

            //read the dimensions and number of nonzeros of the matrix
            if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) != 0)
                exit(1);

            //start the compressed column storage
            //allocate memory
            Lp = (int *) malloc((M + 1) * sizeof(int));
            Li = (int *) malloc(nz * sizeof(int));
            Lx = (T *) malloc(nz * sizeof(T));

            //init variables
            Ij = 1;
            Lp[0] = 0;
            int prev_j = 0;

            //scan the .mtx line by line from 0 to the num of non zeros.
            for (i = 0; i < nz; i++) {

                int curr_j, curr_i;

                /*
                 * Since the matrix for this time is already lower triangular matrix.
                 * there is no need to compare the row index and the column index because
                 * ri >= ci
                 */
                fscanf(f, "%d %d %lg\n", &curr_i, &curr_j, &Lx[i]);

                //the array starts with 0
                Li[i] = curr_i - 1;

                //Check whether the new column starts
                if (prev_j != curr_j && curr_j != 1) {
                    Lp[Ij] = i;
                    Ij++;
                    prev_j = curr_j;
                }
            }

            //set the last column pointer to the size of the nzs
            Lp[Ij] = nz;
            Ij++;

            //close file
            if (f != stdin)
                fclose(f);
        }

        //get the number of nonzeros
        int getNz() {
            return this->nz;
        }

        //get the column pointer
        int *getLp() {
            return this->Lp;
        }

        //get the row index
        int *getLi() {
            return this->Li;
        }

        //get the real values
        T *getLx() {
            return this->Lx;
        }

        //get the size of L
        int getSize() {
            return this->M;
        }

        //print matrix in ccs format
        void print() {

            fprintf(stdout, "Lx:");
            for (i = 0; i < nz; i++) {
                fprintf(stdout, " %2g ", Lx[i]);
            }
            fprintf(stdout, "\n Li:");
            for (i = 0; i < nz; i++) {
                fprintf(stdout, " %d ", Li[i]);
            }
            fprintf(stdout, "\n Lp:");
            for (i = 0; i < Ij; i++) {
                fprintf(stdout, " %d ", Lp[i]);
            }
            fprintf(stdout, "\n");
            fprintf(stdout, "\n");

        }

        //deconstructor
        ~SparseMatrix() {
        }

    };
}


#endif //SPARSEOP2_SPARSEMATRIX_H
