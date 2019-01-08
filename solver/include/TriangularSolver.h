//
// Created by shilei on 1/8/19.
//

#ifndef SPARSEOP2_TRIANGULARSOLVER_H
#define SPARSEOP2_TRIANGULARSOLVER_H

#include <SparseMatrix.h>
#include <Vector.h>
#include <Eigen/SparseCore>
#include <Eigen/Sparse>
#include <iostream>

using namespace DataStructure;

template<unsigned int VectorType, class T>
class TriangularSolve {

    SparseMatrix<T> A;

    Vector<VectorType, T> x;

public:

    TriangularSolve(SparseMatrix<T> A, Vector<VectorType, T> x) {
        this->A = A;
        this->x = x;
    }

    void setA(SparseMatrix<T> A) {
        this->A = A;
    }

    void setx(Vector<VectorType, T> x) {
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
    int lsolve() {

        /* check inputs */
//    if (!A.getLp() || !A.getLi() || x.getLx() != nullptr)
//        exit(1);

        T *Lx = A.getLx();
        T *Lxx = x.getLx();

        int *Lp = A.getLp();
        int *Li = A.getLi();

        int nz = A.getNz();

        for (int j = 0; j < A.getSize(); j++) {

            Lxx[j] /= Lx[Lp[j]];

            for (int p = Lp[j] + 1; p < Lp[j + 1]; p++) {
                Lxx[Li[p]] -= Lx[p] * Lxx[j];
            }

        }

//        fprintf(stdout, "Solution:");
//        for (int i = 0; i < x.getSize(); i++) {
//            fprintf(stdout, " %f ", Lxx[i]);
//        }

    }


    int verify() {

        //Reinitialize X
        T *Lxx = x.getLx();
        x.read();

        int M = A.getSize();
        int nz = A.getNz();

        T *Lx = A.getLx();
        T *Lxv = x.getLx();

        int *Lp = A.getLp();
        int *Li = A.getLi();

        Eigen::SparseMatrix<double, Eigen::RowMajor, long int> A(M, M);
        Eigen::VectorXd b(M, 1), xV(M, 1);

        typedef  Eigen::Triplet<double> Triplet;
        std::vector<Triplet> triplets;

        triplets.reserve(nz);

        for (int j = 0; j < M; j++) {
            for (int p = Lp[j]; p < Lp[j + 1]; p++) {
//                printf("\n(%d %d %f)", Li[p], j, Lx[p]);
//                printf("\n(%d)", Li[p]);
                triplets.push_back(Triplet(Li[p],j,Lx[p]));
//                A.coeffRef(Li[p], j) = Lx[p];
            }
        }

        A.setFromTriplets(triplets.begin(), triplets.end());

        for (int i = 0; i < M; i++) {
            b(i) = Lxv[i];
        }

//        cout<< A<< endl;
        xV = A.triangularView<Eigen::Lower>().solve(b);

        printf("Verification:");
        for (int i = 0; i < M; i++) {
            if (abs(xV[i] - Lxx[i]) > 1e-10) {
                printf("\n(%d %f)", i, Lxx[i]);
            }
        }

    }
};

#endif //SPARSEOP2_TRIANGULARSOLVER_H
