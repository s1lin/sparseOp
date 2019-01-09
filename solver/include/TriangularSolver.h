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
#include <sys/time.h>
#include <omp.h>
#include <set>
#include <list>
#include <iterator>
#include <algorithm>

using namespace DataStructure;

template<unsigned int VectorType, class T>
class TriangularSolve {

    SparseMatrix<T> A;

    Vector<VectorType, T> x;

    struct ReachSet {
        std::set<int> rs;
        std::list<int> os;
    } reachSet;

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

//    int cholesky() {
//
//        T *Lx = A.getLx();
//        int *Lp = A.getLp();
//        int *Li = A.getLi();
//
//        int nz = A.getNz();
//        int M = A.getSize();
//
//        for (int k = 0; k < M; k++) {
//            if (k == Li[k]) {
//                Lx[Lp[k]] = sqrt(Lx[Lp[k]]); //A(k,k) = Sqrt(A(k,k))
//
//                for (int i = Lp[k] + 1; i < Lp[k+1]; i++) {
//                    Lx[Lp[i]] /= Lx[Lp[k]]; //A(i,k) /= A(k,k)
//                }
//            }
//
//            for (int j = Lp[k] + 1; j < Lp[k+1]; j++) {
//                for (int i = Lp[k]; i < M; i++) {
//                    Lx[] -= Lx[] * Lx[]; //A(i,j) -= A(i,k) * A(j,k)
//                }
//            }
//
//        }
//    }

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

        if (VectorType == VectorType::sparse) {
            analysis();
            lsolve_sparse();
        }

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

    }

    int lsolve_sparse() {

        T *Lx = A.getLx();
        T *Lxx = x.getLx();

        int *Lp = A.getLp();
        int *Li = A.getLi();

        int nz = A.getNz();

        for (int j : reachSet.os) {

            Lxx[j] /= Lx[Lp[j]];

            for (int p = Lp[j] + 1; p < Lp[j + 1]; p++) {
                Lxx[Li[p]] -= Lx[p] * Lxx[j];
            }

        }

    }

    int lsolve_parallel(int num_threads) {

        if (VectorType == VectorType::sparse) {
            analysis();
        }
        T *Lx = A.getLx();
        T *Lxx = x.getLx();

        int *Lp = A.getLp();
        int *Li = A.getLi();

        int nz = A.getNz();
        int j, p;

        omp_set_num_threads(num_threads);


        for (j = 0; j < A.getSize(); j++) {
            Lxx[j] /= Lx[Lp[j]];

#pragma omp parallel shared(Lx, Lxx, Lp, Li) private(p)
            cout << omp_get_thread_num() << " ";
#pragma omp for

            for (p = Lp[j] + 1; p < Lp[j + 1]; p++) {
                Lxx[Li[p]] -= Lx[p] * Lxx[j];
            }

        }
    }


    int analysis() {

        int e = 1; //level number

        int M = A.getSize();//Matrix Size
        int nz = A.getNz();//Matrix Size
        T *Lx = A.getLx();
        int *Lp = A.getLp();
        int *Li = A.getLi();

        std::set<int> nzB = x.getNzB();

        for (int i : nzB) {

            if (this->reachSet.rs.insert(i).second) {
                this->reachSet.os.push_back(i);
            }

            int index = i;
            int j = Lp[index];

            while (j < Lp[index + 1]) {

                if (Li[j] == index) {
                    if ((Lp[index + 1] - Lp[index] == 1)) {
                        break;
                    }
                    j++;
                    continue;
                }
                cout << "Li[" << j << "]:" << Li[j] << "     ";
                if (!this->reachSet.rs.insert(Li[j]).second) {
                    j++;
                    continue;
                }
                this->reachSet.os.push_back(Li[j]);

                for (int p = j + 1; p < nz; p++) {
                    if (Li[p] == Li[j] && !this->reachSet.rs.insert(Li[p]).second) {
                        cout << "Inner: Li[" << p << "]:" << Li[p] << " \n";
                        index = Li[p];
                        j = Lp[index];
                        break;
                    }
                }

            }
        }

        cout << endl << "ReachSet: Size=" << reachSet.rs.size() << "\n value: ";

        for (int i:reachSet.os) {
            cout << i << " ";
        }
        cout << endl;
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

        typedef Eigen::Triplet<double> Triplet;
        std::vector<Triplet> triplets;

        triplets.reserve(nz);

        for (int j = 0; j < M; j++) {
            for (int p = Lp[j]; p < Lp[j + 1]; p++) {
                triplets.push_back(Triplet(Li[p], j, Lx[p]));
            }
        }

        A.setFromTriplets(triplets.begin(), triplets.end());

//        cout << A;

        for (int i = 0; i < M; i++) {
            b(i) = Lxv[i];
        }

        struct timeval tim;

        cout << "Eigen-------" << endl;
        gettimeofday(&tim, NULL);
        double t1 = tim.tv_sec + (tim.tv_usec / 1e+6);

        xV = A.triangularView<Eigen::Lower>().solve(b);

        gettimeofday(&tim, NULL);
        double t2 = tim.tv_sec + (tim.tv_usec / 1e+6);
        cout << "Eigen Solve Used:" << t2 - t1 << "s." << endl;

        printf("Verification:");
        int index = 0;
        for (index = 0; index < M; index++) {
            if (abs(xV[index] - Lxx[index]) > 1e-10) {
                printf("\n(%d %f)", index, abs(xV[index] - Lxx[index]));
            }
        }

        if (index == M) {
            printf("Clear!\n");
        }
    }
};

#endif //SPARSEOP2_TRIANGULARSOLVER_H
