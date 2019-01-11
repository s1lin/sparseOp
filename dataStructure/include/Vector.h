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

    template<class T>
    class Vector {

        //return code from mmio
        int ret_code;

        //type of the matrix
        MM_typecode matcode;

        //.mtx file name
        const char *mtx_file;

        /*
         * M: num of rows
         * nz: num of non zeros;
         */
        int M, nz;

        //value of RHS
        T *Lx;

        //store .mtx file
        FILE *f;

        //type of the Vector
        int VT;//dense or sparse

        //the indices of nonzeros of RHS
        std::set<int> nzB;

    public:

        //default constructor
        Vector() {}

        //constructor
        Vector(const char *mtx) {
            this->mtx_file = mtx;
        }

        //read routine and depends on the type of the vector and invokes different routines.
        void read() {

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
            int N;
            if ((ret_code = mm_read_mtx_array_size(f, &M, &N)) != 0)
                exit(1);

            //allocate memory
            Lx = (T *) malloc(M * sizeof(T));

            //chose different reader based on the vector type
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

        //get the size of the RHS
        int getSize() {
            return this->M;
        }

        //get the value of RHS
        T *getLx() {
            return this->Lx;
        }

        //get the indices of the nonzeros in RHS
        std::set<int> getNzB() {
            return this->nzB;
        }

        //set the type of the vector
        void setVT(VectorType vt) {
            this->VT = vt;
        }

        //print vector
        void print() {

            for (int i = 0; i < M; i++) {

                fprintf(stdout, "%2g ", Lx[i]);
            }
            fprintf(stdout, "\n");
            fprintf(stdout, "\n");
        }

        ~Vector() {}

    private:

        void readDense(const char *mtx) {

            //scan the .mtx line by line from 0 to the size of the RHS
            for (int i = 0; i < M; i++)
                fscanf(f, "%lg\n", &Lx[i]);

            //close file
            if (f != stdin)
                fclose(f);
        }

        void readSparse(const char *mtx) {

            //initialize variables
            T curr_val;

            for (int i = 0; i < nz; i++) {

                int curr_i, curr_j;

                //scan line by line
                fscanf(f, "%d %d %lg\n", &curr_i, &curr_j, &curr_val);

                //store
                Lx[curr_i - 1] = curr_val;

                //also put in the current index of non zero element
                nzB.insert(curr_i - 1);//B = {i|bi != 0} the index of nonzero elements of b.
            }

            //close file
            if (f != stdin)
                fclose(f);
        }

    };


}


#endif //SPARSEOP2_VECTOR_H
