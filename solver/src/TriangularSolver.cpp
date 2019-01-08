//
// Created by shilei on 1/8/19.
//

#include <SparseMatrix.h>
#include <Vector.h>
#include <TriangularSolver.h>
#include <cstdio>

using namespace DataStructure;

template<unsigned int VT, class T>
TriangularSolve<VT, T>::TriangularSolve(SparseMatrix <T> A, Vector <VT, T> x) {
    this->A = A;
    this->x = x;
}


template<unsigned int VT, class T>
void TriangularSolve<VT, T>::setA(SparseMatrix <T> A) {
    this->A = A;
}

template<unsigned int VT, class T>
void TriangularSolve<VT, T>::setx(Vector <VT, T> x) {
    this->x = x;
}

/*
* Lower triangular solver Lx=b
* L is stored in the compressed column storage format
* Inputs are:
* n : the matrix dimension
* Lp : the column pointer of L
* Li : the row index of L
* Lx : the values of L
* In/Out:
* x : the right hand-side b at start and the solution x at the end.
*/
template<unsigned int VT, class T>
int TriangularSolve<VT, T>::lsolve() {


    /* check inputs */
//    if (!A.getLp() || !A.getLi() || x.getLx() != nullptr)
//        exit(1);

    T *Lx = A.getLx();
    T *Lxx = x.getLx();

    int *Lp = A.getLp();
    int *Li = A.getLi();

    for (long j = 0; j < A.getSize(); j++) {

        Lxx[j] /= Lx[Lp[j]];

        for (int p = Lp[j] + 1; p < Lp[j + 1]; p++) {
            Lxx[Li[p]] -= Lx[p] * Lxx[j];
        }

    }

    fprintf(stdout, "Solution:");
    for (int i = 0; i < x.getSize(); i++) {
        fprintf(stdout, " %f ", Lxx[i]);
    }

}
