//
// Created by shilei on 1/8/19.
//

#include "SparseMatrix.h"
#include "Vector.h"
#include "TriangularSolver.h"

using namespace DataStructure;
using namespace std;

int main() {

    SparseMatrix<double> A("/home/shilei/CLionProjects/sparseOp/matrix/af_0_k101.mtx");
    Vector<VectorType::dense, double> b("/home/shilei/CLionProjects/sparseOp/matrix/b_dense_af_0_k101.mtx");

//    SparseMatrix<double> A("/home/shilei/CLionProjects/sparseOp/matrix/b1_ss.mtx");
//    Vector<VectorType::dense, double> b("/home/shilei/CLionProjects/sparseOp/matrix/b1_ss_b.mtx");

    cout << "Begin initialize A" << endl;
    A.read();
    cout << "Finish initializing A" << endl;
    cout << "Begin initializing b" << endl;
    b.read();
    cout << "Finish initializing b" << endl;

//    A.print();
//    b.print();

    TriangularSolve<VectorType::dense, double> triangularSolve1(A, b);

    cout << "start to solve by the triangular solver" << endl;
    triangularSolve1.lsolve();
    cout << "finish solving." << endl;

    cout << "start to verify solution by Eigen" << endl;
    triangularSolve1.verify();
    cout << "finish verify." << endl;

    return 0;
}