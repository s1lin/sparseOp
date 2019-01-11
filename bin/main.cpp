//
// Created by shilei on 1/8/19.
//

#include "SparseMatrix.h"
#include "Vector.h"
#include "TriangularSolver.h"
#include <string.h>

using namespace DataStructure;
using namespace std;


int main(int argc, char *argv[]) {

    //init time struct
    struct timeval tim;

    //init A and x from the input
    SparseMatrix<double> A(argv[1]);

    Vector<double> b1(argv[3]);
    Vector<double> b2(argv[3]);

    //set the type of vector
    if (strcmp(argv[3], "sparse") != 0) {
        b1.setVT(VectorType::sparse);
        b2.setVT(VectorType::sparse);
    }
    else {
        b1.setVT(VectorType::dense);
        b2.setVT(VectorType::dense);
    }


    cout << "-------Begin initialize A-------" << endl;
    gettimeofday(&tim, NULL);
    double t1 = tim.tv_sec + (tim.tv_usec / 1e+6);

    A.read();

    gettimeofday(&tim, NULL);
    double t2 = tim.tv_sec + (tim.tv_usec / 1e+6);
    cout << "Initializing A Used:" << t2 - t1 << "s." << endl;
    cout << "Finish initializing A" << "\n\n";


    for (int i = 0; i < 5; i++) {
        cout << "*******Num of Execution " << i+1 << "***********\n";
        cout << "-------Begin initialize b-------" << endl;
        gettimeofday(&tim, NULL);
        t1 = tim.tv_sec + (tim.tv_usec / 1e+6);

        b1.read();
        b2.read();

        gettimeofday(&tim, NULL);
        t2 = tim.tv_sec + (tim.tv_usec / 1e+6);
        cout << "Initializing B Used:" << t2 - t1 << "s." << endl;
        cout << "Finish initializing b" << "\n\n";

        TriangularSolve<double> triangularSolve1(A, b1);

        cout << "-------Begin Serial solve-------" << endl;
        gettimeofday(&tim, NULL);
        t1 = tim.tv_sec + (tim.tv_usec / 1e+6);

        triangularSolve1.lsolve();

        gettimeofday(&tim, NULL);
        t2 = tim.tv_sec + (tim.tv_usec / 1e+6);

        double origin_t = t2 - t1;

        cout << "Solve Used:" << origin_t << endl;
        cout << "Finish Solving." << "\n\n";


        TriangularSolve<double> triangularSolve2(A, b2);
        cout << "-------Begin optimized solve-------" << endl;

        gettimeofday(&tim, NULL);
        t1 = tim.tv_sec + (tim.tv_usec / 1e+6);

        triangularSolve2.analysis();

        gettimeofday(&tim, NULL);
        t2 = tim.tv_sec + (tim.tv_usec / 1e+6);

        cout << "Analysis Used:" << t2 - t1 << "s." << endl;
        cout << "Finish Analysis." << "\n\n";

        gettimeofday(&tim, NULL);
        t1 = tim.tv_sec + (tim.tv_usec / 1e+6);

        triangularSolve2.lsolve_sparse();

        gettimeofday(&tim, NULL);
        t2 = tim.tv_sec + (tim.tv_usec / 1e+6);

        double ops_t = t2 - t1;

        cout << "Solve Used:" << ops_t << endl;
        cout << "improve:" << origin_t - ops_t << "s." << endl;

        cout << "-------Begin to verify-------" << endl;
        triangularSolve1.verify();
        triangularSolve2.verify();
        cout << "Finish Verifing.\n\n";
    }

    return 0;
}

