#pragma once

#include <cmath>

/*  * ****************************************************** *
    * Taken from:                                            *
    * http://paulbourke.net/geometry/polygonise/             *
    *                                                        *
    * ****************************************************** *
*/

struct XYZ {
    double x;
    double y;
    double z;
    double scalarValue;

    XYZ();
    XYZ(double* xyz_s);
    double* coords();
    bool operator==(XYZ p);
};

struct TRIANGLE {
    XYZ p[3];
};

struct GRIDCELL {
    XYZ p[8];
    double val[8];
};

XYZ VertexInterp(double isolevel, XYZ p1, XYZ p2, double valp1, double valp2);
int Polygonise(GRIDCELL grid, double isolevel, TRIANGLE* triangles);

