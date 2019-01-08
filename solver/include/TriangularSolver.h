//
// Created by shilei on 1/8/19.
//

#ifndef SPARSEOP2_TRIANGULARSOLVER_H
#define SPARSEOP2_TRIANGULARSOLVER_H
using namespace DataStructure;

template<unsigned int VectorType, class T>
class TriangularSolve {

    SparseMatrix<T> A;

    Vector<VectorType, T> x;

public:

    TriangularSolve(SparseMatrix<T> A, Vector<VectorType, T> x);

    void setA(SparseMatrix<T> A);

    void setx(Vector<VectorType, T> x);

    int lsolve();

    int verify();
};


#endif //SPARSEOP2_TRIANGULARSOLVER_H
