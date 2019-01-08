//
// Created by shilei on 1/8/19.
//

#ifndef SPARSEOP2_SPARSEMATRIX_H
#define SPARSEOP2_SPARSEMATRIX_H

#include <cstdio>

extern "C" {
#include "../../helper/mmio.h"
};

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

        explicit SparseMatrix(const char *);

        void read();

        int *getLp();

        int *getLi();

        T *getLx();

        long getSize();

        void print();

        ~SparseMatrix() {}

    };
}


#endif //SPARSEOP2_SPARSEMATRIX_H
