#include "sdlp.h"
#include "BaseGeoOpUtils.h"
#include "Mesh.h"

using namespace sdlp;

double* sdlpMain(double extremeDirection[3], HalfSpace** halfSpaceSet, int numOfHalfSpaces) {

    int d = 3;
    int m = numOfHalfSpaces;
    Eigen::VectorXd x(d);    // decision variables
    Eigen::VectorXd c(d);    // objective coefficients
    Eigen::MatrixXd A(m, d); // constraint matrix
    Eigen::VectorXd b(m);    // constraint bound

    //c << 0.0, 0.0, 1.0;   // default
    c << extremeDirection[0], extremeDirection[1], extremeDirection[2];

    for (int i = 0; i < numOfHalfSpaces; i++) {
        HalfSpace* hp = halfSpaceSet[i];
        A.row(i) << hp->ABCD[0], hp->ABCD[1], hp->ABCD[2];
        b(i) = -hp->ABCD[3];
    }

    double minobj = linprog(c, A, b, x);

    if (minobj == numeric_limits<double>::infinity()) {
        //cout << "INFEASIBLE" << endl;
        return NULL;
    }
    if (minobj == -numeric_limits<double>::infinity()) {
        //cout << "UNBOUNDED" << endl;
        return NULL;
    }

    double* kernel_point = new double[d];
    for (int i = 0; i < d; i++)
        kernel_point[i] = x[i];
    return kernel_point;

}

double* sdlpMain(Mesh* hostMeshptr, double extremeDirection[3]) {

    vector<HalfSpace> halfSpaceSet;
    computeHalfSpacesFromTriangles(hostMeshptr->getAllTris(), hostMeshptr->getAllVerts(), halfSpaceSet);
    HalfSpace** halfSpaces = new HalfSpace * [halfSpaceSet.size()];
    for (int i = 0; i < halfSpaceSet.size(); i++)
        halfSpaces[i] = &halfSpaceSet[i];

    double* kernel_point = sdlpMain(extremeDirection, halfSpaces, halfSpaceSet.size());
    delete[] halfSpaces;

    return kernel_point;
}
