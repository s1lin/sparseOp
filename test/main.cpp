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
    Vector<VectorType::sparse, double> b("/home/shilei/CLionProjects/sparseOp/matrix/b_sparse_af_0_k101.mtx");

//    SparseMatrix<double> A("/home/shilei/CLionProjects/sparseOp/matrix/b1_ss.mtx");
//    Vector<VectorType::dense, double> b("/home/shilei/CLionProjects/sparseOp/matrix/b1_ss_b.mtx");

    cout << "-------Begin initialize A-------" << endl;

    gettimeofday(&tim, NULL);
    double t1 = tim.tv_sec + (tim.tv_usec / 1e+6);

    A.read();

    gettimeofday(&tim, NULL);
    double t2 = tim.tv_sec + (tim.tv_usec / 1e+6);
    cout << "Initializing A Used:" << t2 - t1 << "s." << endl;
    cout << "Finish initializing A" << "\n\n";


    cout << "-------Begin initializing b-------" << endl;
    gettimeofday(&tim, NULL);
    t1 = tim.tv_sec + (tim.tv_usec / 1e+6);

    b.read();

    gettimeofday(&tim, NULL);
    t2 = tim.tv_sec + (tim.tv_usec / 1e+6);
    cout << "Initializing B Used:" << t2 - t1 << "s." << endl;
    cout << "Finish initializing b" << "\n\n";

//    A.print();
//    b.print();

    TriangularSolve<VectorType::sparse, double> triangularSolve1(A, b);


    cout << "-------Begin to solve-------" << endl;
    gettimeofday(&tim, NULL);
    t1 = tim.tv_sec + (tim.tv_usec / 1e+6);

    triangularSolve1.lsolve();

    gettimeofday(&tim, NULL);
    t2 = tim.tv_sec + (tim.tv_usec / 1e+6);
    cout << "Solve Used:" << t2 - t1 << "s." << endl;
    cout << "Finish Solving." << "\n\n";


    cout << "-------Begin to verify-------" << endl;
    triangularSolve1.verify();
    cout << "Finish Verifing." << endl;

    return 0;
}