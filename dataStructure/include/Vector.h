//
// Created by shilei on 1/8/19.
//

#ifndef SPARSEOP2_VECTOR_H
#define SPARSEOP2_VECTOR_H

#include <cstdio>
#include <set>
#include <iterator>
#include "Constants.h"

extern "C" {
#include "../../helper/mmio.h"
};

using namespace std;

namespace DataStructure {

    template<unsigned int VT, class T>
    class Vector {

        int ret_code, i, M, nz;

        MM_typecode matcode;

        FILE *f;

        T *Lx;

        const char *mtx_file;

        std::set<int> nzB;

    public:
        Vector() {

        }

        Vector(const char *mtx) {
            this->mtx_file = mtx;
        }

        void read() {
            switch (VT) {
                case 0:
                    readDense(this->mtx_file);
                    break;
                case 1:
                    readSparse(this->mtx_file);
                    break;
                default:
                    printf("No matching Vector Type.\n");
                    break;
            }
        }

        int getSize() {
            return this->M;
        }

        T *getLx() {
            return this->Lx;
        }

        std::set<int> getNzB(){
            return this->nzB;
        }


        void print() {
            /* print matrix */
            for (i = 0; i < M; i++) {
                if (Lx[i] != 0.0)
                    fprintf(stdout, "%2g ", Lx[i]);
            }
            fprintf(stdout, "\n");
            fprintf(stdout, "\n");
        }

        ~Vector() {}

    private:

        void readDense(const char *mtx) {
            mtx_file = mtx;
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
            int N;
            if ((ret_code = mm_read_mtx_array_size(f, &M, &N)) != 0)
                exit(1);

            /* compressed column storage */

            Lx = (T *) malloc(M * sizeof(T));

            for (i = 0; i < M; i++) {
                T curr_val;
                fscanf(f, "%lg\n", &curr_val);
                Lx[i] = curr_val;
            }

            if (f != stdin) {
                fclose(f);
            }
        }

        void readSparse(const char *mtx) {
            mtx_file = mtx;
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
            int N;
            if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) != 0)
                exit(1);


            /* compressed column storage */

            Lx = (T *) malloc(M * sizeof(T));

            for (i = 0; i < nz; i++) {
                int curr_j, curr_i;
                T curr_val;
                fscanf(f, "%d %d %lg\n", &curr_i, &curr_j, &curr_val);
                Lx[curr_i - 1] = curr_val;
                nzB.insert(curr_i - 1);//B = {i|bi != 0} the index of nonzero elements of b.
            }

            if (f != stdin) {
                fclose(f);
            }

        }
    };
}


#endif //SPARSEOP2_VECTOR_H
