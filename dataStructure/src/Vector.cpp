//
// Created by shilei on 1/8/19.
//

#include <cstdlib>
#include <Vector.h>

using namespace DataStructure;

template<unsigned int VT, class T>
Vector<VT, T>::Vector(const char *mtx) {
    this->mtx_file = mtx;
}

template<unsigned int VT, class T>
void Vector<VT, T>::read() {

    switch (VT) {
        case 0:
            Vector::readDense(this->mtx_file);
            break;
        case 1:
            Vector::readSparse(this->mtx_file);
            break;
        default:
            printf("No matching Vector Type.\n");
            break;
    }
}

template<unsigned int VT, class T>
void Vector<VT, T>::readDense(const char *mtx) {
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
    }

    if (f != stdin) {
        fclose(f);
    }
}

template<unsigned int VT, class T>
void Vector<VT, T>::readSparse(const char *mtx) {
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
    }

    if (f != stdin) {
        fclose(f);
    }
}

template<unsigned int VT, class T>
int Vector<VT, T>::getSize() {
    return this->M;
}

template<unsigned int VT, class T>
T *Vector<VT, T>::getLx() {
    return this->Lx;
}

template<unsigned int VT, class T>
void Vector<VT, T>::print() {

    /* print matrix */
    for (i = 0; i < M; i++) {
        if (Lx[i] != 0.0)
            fprintf(stdout, "%2g ", Lx[i]);
    }
    fprintf(stdout, "\n");
    fprintf(stdout, "\n");

}
