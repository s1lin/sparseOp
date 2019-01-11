//
// Created by shilei on 1/8/19.
//

#ifndef SPARSEOP2_TRIANGULARSOLVER_H
#define SPARSEOP2_TRIANGULARSOLVER_H

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
class TriangularSolve {

    //Sparse Matrix
    SparseMatrix<T> A;

    //Right hand side
    Vector<VectorType, T> x;

    //Reach set generates from Matrix
    stack<int> reachSet;

public:

    //Constructor
    TriangularSolve(SparseMatrix<T> A, Vector<VectorType, T> x) {
        this->A = A;
        this->x = x;
    }

    //set Matrix
    void setA(SparseMatrix<T> A) {
        this->A = A;
    }

    //set RHS
    void setx(Vector<VectorType, T> x) {
        this->x = x;
    }

    /*
    * Lower triangular solver Lx=b
    * L is stored in the compressed column storage format
    * In/Out:
    * x : the right hand-side b at start and the solution x at the end.
    */
    void lsolve() {

        //Value of :
        T *Lx = A.getLx();

        //value of Right hand side
        T *Lxx = x.getLx();

        //Lp : the column pointer of L
        int *Lp = A.getLp();

        //Li : the row index of L
        int *Li = A.getLi();

        //backward substitute
        for (int j = 0; j < A.getSize(); j++) {

            Lxx[j] /= Lx[Lp[j]];

            for (int p = Lp[j] + 1; p < Lp[j + 1]; p++) {
                Lxx[Li[p]] -= Lx[p] * Lxx[j];
            }
        }
    }

    /*
    * Lower triangular solver Lx=b
    * A is stored in the compressed column storage format
    * x : the right hand-side b at start and the solution x at the end.
    */
    void lsolve_sparse() {

        //Value of :
        T *Lx = A.getLx();

        //value of Right hand side
        T *Lxx = x.getLx();

        //Lp : the column pointer of L
        int *Lp = A.getLp();

        //Li : the row index of L
        int *Li = A.getLi();

        //backward substitute
        while (!reachSet.empty()) {


            int j = reachSet.top();

            Lxx[j] /= Lx[Lp[j]];

            for (int p = Lp[j] + 1; p < Lp[j + 1]; p++) {
                Lxx[Li[p]] -= Lx[p] * Lxx[j];
            }

            //pop the first element.
            reachSet.pop();
        }
    }

    /*
     * Generate reach set;
     * A is stored in the compressed column storage format;
     * g is the Directed Acyclic Graph for adjacency;
     * reachset generates from depth-first search on g.
    */
    int analysis() {

        //Lp : the column pointer of L
        int *Lp = A.getLp();

        //Li : the row index of L
        int *Li = A.getLi();

        //nzB: the nonzero elements of x.
        set<int> nzB = x.getNzB();

        //create Graph objects
        Graph g(A.getSize(), nzB);

        //adding edges
        for (int j = 0; j < A.getSize(); j++)
            for (int p = Lp[j] + 1; p < Lp[j + 1]; p++)
                g.addEdge(j, Li[p]);

        //generates reachSet
        this->reachSet = g.dfs();
    }

    /*
     * Verify Solution
     * A is stored in the compressed column storage format;
     * g is the Directed Acyclic Graph for adjacency;
     * reachset generates from depth-first search on g.
    */
    int verify() {

        //The solution from lsolve()
        T *Lxx = x.getLx();

        //Reinitialize X
        x.read();

        //The size of Matrix
        int M = A.getSize();

        //The nonzero elements of A
        int nz = A.getNz();

        //Value of :
        T *Lx = A.getLx();

        //value of Right hand side
        T *Lxv = x.getLx();

        //Lp : the column pointer of L
        int *Lp = A.getLp();

        //Li : the row index of L
        int *Li = A.getLi();

        //create verify matrix Av
        Eigen::SparseMatrix<double, Eigen::RowMajor, long int> Av(M, M);

        //create right-hand side vector b
        Eigen::VectorXd b(M, 1);

        //create left-hand side xV
        Eigen::VectorXd xV(M, 1);

        /*
         * To initialize Av from Eigen::Triplet which were pushed from the nonzeros A
         * triplets: vectors to hold nonzeros of A
         */
        typedef Eigen::Triplet<double> Triplet;
        vector<Triplet> triplets;

        triplets.reserve(nz);

        //push all nonzeros to triplets
        for (int j = 0; j < M; j++) {
            for (int p = Lp[j]; p < Lp[j + 1]; p++) {
                triplets.push_back(Triplet(Li[p], j, Lx[p]));
            }
        }

        //set Av from triplets
        Av.setFromTriplets(triplets.begin(), triplets.end());

        //set b from Lxv
        for (int i = 0; i < M; i++) {
            b(i) = Lxv[i];
        }


        //Solve By Eigen triangular solver
        struct timeval tim;

        cout << "Eigen-------" << endl;
        gettimeofday(&tim, NULL);
        double t1 = tim.tv_sec + (tim.tv_usec / 1e+6);

        xV = Av.triangularView<Eigen::Lower>().solve(b);

        gettimeofday(&tim, NULL);
        double t2 = tim.tv_sec + (tim.tv_usec / 1e+6);
        cout << "Eigen Solve Used:" << t2 - t1 << "s." << endl;

        /*
         * To verify the result, I did xV[i] - Lxx[i] to see whether the value is larger than 1e-10.
         * xV: The solution generates from Eigen triangular solver
         * Lxx: The solution generates from lsolve()
         */
        printf("Verification:");
        for (int index = 0; index < M; index++) {
            if (abs(xV[index] - Lxx[index]) > 1e-10) {
                printf("\n(%d %f)", index, abs(xV[index] - Lxx[index]));
            }
        }

        printf("Clear!\n");

    }
};

#endif //SPARSEOP2_TRIANGULARSOLVER_H
