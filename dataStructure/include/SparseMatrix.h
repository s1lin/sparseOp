//
// Created by shilei on 1/8/19.
//

#ifndef SPARSEOP2_SPARSEMATRIX_H
#define SPARSEOP2_SPARSEMATRIX_H

#include <cstdio>
#include <cstdlib>

#include "../../helper/mmio.h"


namespace DataStructure {

    template<class T>
    class SparseMatrix {

        int ret_code;
        MM_typecode matcode;
        FILE *f;
        int M, N, nz;
        //Ij:index of current Lp: Column start Index and Li:Row Indices
        int i, Ij, *Lp, *Li;
        T *Lx;
        const char *mtx_file;

    public:
        SparseMatrix() {

        }

        SparseMatrix(const char *mtx) {
            this->mtx_file = mtx;
        }

        void read() {
            if ((f = fopen(mtx_file, "r")) == nullptr)
                exit(1);


            if (mm_read_banner(f, &matcode) != 0) {
                printf("Could not process Matrix Market banner.\n");
                exit(1);
            }


            if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
                mm_is_sparse(matcode)) {
                printf("Not support Matrix type: [%s]\n", mm_typecode_to_str(matcode));
                exit(1);
            }

            /* dimension*/

            if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) != 0)
                exit(1);


            /* compressed column storage */

            Lp = (int *) malloc((M + 1) * sizeof(int));
            Li = (int *) malloc(nz * sizeof(int));
            Lx = (T *) malloc(nz * sizeof(T));
            Ij = 1;
            Lp[0] = 0;
            int prev_j = 0;
            for (i = 0; i < nz; i++) {
                int curr_j, curr_i;
                T curr_val;
                fscanf(f, "%d %d %lg\n", &curr_i, &curr_j, &Lx[i]);
                //todo:Only read the lower part for now
                Li[i] = curr_i - 1;
                if (prev_j != curr_j && curr_j != 1) {
                    Lp[Ij] = i;
                    Ij++;
                    prev_j = curr_j;
                }
            }
            Lp[Ij] = nz;
            Ij++;

            if (f != stdin) {
                fclose(f);
            }
        }

        int getNz() {
            return this->nz;
        }

        int *getLp() {
            return this->Lp;
        }

        int *getLi() {
            return this->Li;
        }

        T *getLx() {
            return this->Lx;
        }

        int getSize() {
            return this->M;
        }

        void print() {
            /* print matrix */
            for (i = 0; i < nz; i++) {
                fprintf(stdout, " Lx: %2g ", Lx[i]);
            }
            fprintf(stdout, "\n");
            for (i = 0; i < nz; i++) {
                fprintf(stdout, " Li: %d ", Li[i]);
            }
            fprintf(stdout, "\n");
            for (i = 0; i < Ij; i++) {
                fprintf(stdout, " Lp: %d ", Lp[i]);
            }
            fprintf(stdout, "\n");
            fprintf(stdout, "\n");

        }

        ~SparseMatrix() {}

    };
}


#endif //SPARSEOP2_SPARSEMATRIX_H
