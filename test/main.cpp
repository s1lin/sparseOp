//
// Created by shilei on 1/8/19.
//
//#include "structure/SparseMatrix.h"
//#include "structure/Vector.h"

#include "SparseMatrix.h"
#include "Vector.h"
#include "TriangularSolver.h"

//#include <Eigen/Dense>
using namespace DataStructure;


int main() {
    SparseMatrix<double> A("/home/shilei/CLionProjects/SparseTriangularSolve/matrix/b1_ss.mtx");
    Vector<VectorType::sparse, double> b("/home/shilei/CLionProjects/SparseTriangularSolve/matrix/b1_ss_b.mtx");

    A.read();
    b.read();

    A.print();
    b.print();

    TriangularSolve<VectorType::sparse, double> triangularSolve1(A, b);

    triangularSolve1.setA(A);
    triangularSolve1.setx(b);
    triangularSolve1.lsolve();

//    Eigen::MatrixXd m(7, 7);
//    m(4, 0) = -.03599942;
//    m(0, 0) = 1;
//    m(5, 0) = -.0176371;
//    m(6, 0) = -.007721779;
//    m(1, 1) = -2;
//    m(2, 2) = -2;
//    m(3, 3) = -2;
//    m(4, 4) = 1;
//    m(5, 5) = 1;
//    m(6, 6) = 1;
//
//    Eigen::VectorXd bt(7), xT(7);
//    bt(0) = -.0001;
//    bt(1) = .1167;
//    bt(2) = -.2333;
//    bt(3) = .1167;
//    bt(4) = -.4993128;
//    bt(5) = .3435885;
//    bt(6) = .7467878;
////
////    std::cout << bt << std::endl;
////
//    xT = m.triangularView<Eigen::Lower>().solve(bt);
////
//    printf("\n(%g, %g, %g, %g, %g, %g, %g)", xT(0), xT(1), xT(2), xT(3), xT(4), xT(5), xT(6));

    return 0;
}