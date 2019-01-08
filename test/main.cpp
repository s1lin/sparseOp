//
// Created by shilei on 1/8/19.
//

#include "SparseMatrix.h"
#include "Vector.h"
#include "TriangularSolver.h"

#include <Eigen/Dense>
using namespace DataStructure;


int main() {

//    SparseMatrix<double> A("/home/shilei/CLionProjects/sparseOp/matrix/af_0_k101.mtx");
//    Vector<VectorType::dense, double> b("/home/shilei/CLionProjects/sparseOp/matrix/b_dense_af_0_k101.mtx");


    SparseMatrix<double> A("/home/shilei/CLionProjects/sparseOp/matrix/b1_ss.mtx");
    Vector<VectorType::dense, double> b("/home/shilei/CLionProjects/sparseOp/matrix/b1_ss_b.mtx");

    A.read();
    b.read();

    A.print();
    b.print();

    TriangularSolve<VectorType::dense, double> triangularSolve1(A, b);

    triangularSolve1.lsolve();
    triangularSolve1.verify();

    return 0;
}