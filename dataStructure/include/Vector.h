//
// Created by shilei on 1/8/19.
//

#ifndef SPARSEOP2_VECTOR_H
#define SPARSEOP2_VECTOR_H

#include <cstdio>
#include "Constants.h"

extern "C" {
#include "../../helper/mmio.h"
};

using namespace std;

namespace DataStructure {
    
    template<unsigned int VT, class T>
    class Vector {
        int ret_code;
        MM_typecode matcode;
        FILE *f;
        int i, M, nz;
        T *Lx;
        const char *mtx_file;

    public:

        explicit Vector(const char *);

        void read();

        int getSize();

        T *getLx();

        void print();

        ~Vector(){}

    private:

        void readDense(const char *);

        void readSparse(const char *);

    };

}

#endif //SPARSEOP2_VECTOR_H
