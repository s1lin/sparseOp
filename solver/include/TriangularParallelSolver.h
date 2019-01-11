//
// Created by shilei on 1/8/19.
//

#ifndef SPARSEOP2_TRIANGULARPARALLELSOLVER_H
#define SPARSEOP2_TRIANGULARPARALLELSOLVER_H

#include <SparseMatrix.h>
#include <Vector.h>
#include <Graph.h>

#include <Eigen/SparseCore>
#include <Eigen/Sparse>

#include <iostream>
#include <sys/time.h>
#include <omp.h>
#include <set>
#include <list>
#include <iterator>
#include <stack>
#include <algorithm>

using namespace std;
using namespace DataStructure;

template<unsigned int VectorType, class T>
class TriangularParallelSolver {

    SparseMatrix<T> A;

    Vector<VectorType, T> x;

    stack<int> reachSet;

    multiset<int> usedNodes;

public:

    TriangularParallelSolver(SparseMatrix<T> A, Vector<VectorType, T> x) {
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

        T *Lx = A.getLx();
        T *Lxx = x.getLx();

        int *Lp = A.getLp();
        int *Li = A.getLi();

        for (int j = 0; j < A.getSize(); j++) {

            Lxx[j] /= Lx[Lp[j]];
            int p;

#pragma omp parallel for shared(Lx, Lxx, Li, Lp) private(p)
            for (p = Lp[j] + 1; p < Lp[j + 1]; p++) {
                Lxx[Li[p]] -= Lx[p] * Lxx[j];
            }

        }
    }

    int lsolve_sparse() {

        T *Lx = A.getLx();
        T *Lxx = x.getLx();

        int *Lp = A.getLp();
        int *Li = A.getLi();

        while (!reachSet.empty()) {
            int j = reachSet.top();

            Lxx[j] /= Lx[Lp[j]];
            int p;
//#pragma omp parallel for shared(Lx, Lxx) private(p)
            for (p = Lp[j] + 1; p < Lp[j + 1]; p++) {
                Lxx[Li[p]] -= Lx[p] * Lxx[j];
            }

            reachSet.pop();
        }
    }

    int analysis() {

        int *Lp = A.getLp();
        int *Li = A.getLi();

        set<int> nzB = x.getNzB();

        Graph graph(A.getSize(), nzB);

        #pragma omp parallel for
        for (int j = 0; j < A.getSize(); j++) {
            for (int p = Lp[j] + 1; p < Lp[j + 1]; p++) {
                graph.addEdge(j, Li[p]);
            }
        }

        this->reachSet = graph.dfs();
    }

    int verify() {

        T *Lxx = x.getLx();

        //Reinitialize X
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
        vector<Triplet> triplets;

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

#endif //SPARSEOP2_TRIANGULARPARALLELSOLVER_H
