//
// Created by shilei on 1/8/19.
//

#include "SparseMatrix.h"
#include "Vector.h"
#include "TriangularSolver.h"

using namespace DataStructure;
using namespace std;

int main() {

    struct timeval tim;

    SparseMatrix<double> A("/home/shilei/CLionProjects/sparseOp/matrix/af_0_k101.mtx");
//    Vector<VectorType::sparse, double> b("/home/shilei/CLionProjects/sparseOp/matrix/b_sparse_af_0_k101.mtx");
//    Vector<VectorType::sparse, double> b2("/home/shilei/CLionProjects/sparseOp/matrix/b_sparse_af_0_k101.mtx");
    Vector<VectorType::sparse, double> b3("/home/shilei/CLionProjects/sparseOp/matrix/b_sparse_af_0_k101.mtx");
//    Vector<VectorType::dense, double> b("/home/shilei/CLionProjects/sparseOp/matrix/b_dense_af_0_k101.mtx");
////
////    SparseMatrix<double> A("/home/shilei/CLionProjects/sparseOp/matrix/testA.mtx");
////    SparseMatrix<double> A("/home/shilei/CLionProjects/sparseOp/matrix/b1_ss.mtx");
////    Vector<VectorType::sparse, double> b("/home/shilei/CLionProjects/sparseOp/matrix/b1_ss_b_sparse.mtx");
////    Vector<VectorType::sparse, double> b2("/home/shilei/CLionProjects/sparseOp/matrix/b1_ss_b_sparse.mtx");
////    Vector<VectorType::sparse, double> b("/home/shilei/CLionProjects/sparseOp/matrix/testA_b_sparse.mtx");
////    Vector<VectorType::sparse, double> b2("/home/shilei/CLionProjects/sparseOp/matrix/testA_b_sparse.mtx");
//
//    cout << "-------Begin initialize A-------" << endl;
//
    gettimeofday(&tim, NULL);
    double t1 = tim.tv_sec + (tim.tv_usec / 1e+6);

    A.read();

    gettimeofday(&tim, NULL);
    double t2 = tim.tv_sec + (tim.tv_usec / 1e+6);
//    cout << "Initializing A Used:" << t2 - t1 << "s." << endl;
//    cout << "Finish initializing A" << "\n\n";
//
////    A.print();
//
//    cout << "-------Begin initialize b-------" << endl;
//    gettimeofday(&tim, NULL);
//    t1 = tim.tv_sec + (tim.tv_usec / 1e+6);
//
//    b.read();
//    b2.read();
    b3.read();
//
//    gettimeofday(&tim, NULL);
//    t2 = tim.tv_sec + (tim.tv_usec / 1e+6);
//    cout << "Initializing B Used:" << t2 - t1 << "s." << endl;
//    cout << "Finish initializing b" << "\n\n";
//
//    TriangularSolve<VectorType::sparse, double> triangularSolve1(A, b);
////    TriangularSolve<VectorType::dense, double> triangularSolve1(A, b);
//
//    cout << "-------Begin Serial solve-------" << endl;
//    gettimeofday(&tim, NULL);
//    t1 = tim.tv_sec + (tim.tv_usec / 1e+6);
//
//    triangularSolve1.lsolve();
//
//    gettimeofday(&tim, NULL);
//    t2 = tim.tv_sec + (tim.tv_usec / 1e+6);
//
//    double origin_t = t2- t1;
//
//    cout << "Solve Used:" << origin_t << "s." << endl;
//    cout << "Finish Analysis." << "\n\n";
//
////    b.print();
//
//    TriangularSolve<VectorType::sparse, double> triangularSolve2(A, b2);
//    cout << "-------Begin optimized solve-------" << endl;
//
//    gettimeofday(&tim, NULL);
//    t1 = tim.tv_sec + (tim.tv_usec / 1e+6);
//
//    triangularSolve2.analysis();
//
//    gettimeofday(&tim, NULL);
//    t2 = tim.tv_sec + (tim.tv_usec / 1e+6);
//
//    cout << "Analysis Used:" << t2-t1 << "s." << endl;
//    cout << "Finish Analysis." << "\n\n";
//
//    t1 = tim.tv_sec + (tim.tv_usec / 1e+6);
//
//    triangularSolve2.lsolve_sparse();
//
//    gettimeofday(&tim, NULL);
//    t2 = tim.tv_sec + (tim.tv_usec / 1e+6);
//
//    double ops_t = t2- t1;
//
//    cout << "Solve Used:" << ops_t << "s." << endl;
//    cout << "improve:" << origin_t - ops_t << "s." << endl;
//    cout << "Finish Solving." << "\n\n";


    TriangularSolve<VectorType::sparse, double> triangularSolve3(A, b3);
    cout << "-------Begin Parallel solve-------" << endl;

    gettimeofday(&tim, NULL);
    t1 = tim.tv_sec + (tim.tv_usec / 1e+6);

    triangularSolve3.analysis();

    gettimeofday(&tim, NULL);
    t2 = tim.tv_sec + (tim.tv_usec / 1e+6);

    cout << "Analysis Used:" << t2-t1 << "s." << endl;
    cout << "Finish Solving." << "\n\n";

    t1 = tim.tv_sec + (tim.tv_usec / 1e+6);

    triangularSolve3.lsolve_parallel(2);

    gettimeofday(&tim, NULL);
    t2 = tim.tv_sec + (tim.tv_usec / 1e+6);

    double p_t = t2- t1;

    cout << "Solve Used:" << p_t << "s." << endl;
//    cout << "improve:" << origin_t - p_t << "s." << endl;
    cout << "Finish Solving." << "\n\n";

//    b.print();

    cout << "-------Begin to verify-------" << endl;
//    triangularSolve1.verify();
//    triangularSolve2.verify();
    triangularSolve3.verify();
    cout << "Finish Verifing." << endl;

    return 0;
}