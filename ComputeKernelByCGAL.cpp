// @taken_from https://doc.cgal.org/latest/Convex_hull_3/index.html

#include <CGAL/Convex_hull_3/dual/halfspace_intersection_with_constructions_3.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>
#include <list>

#include "BaseGeoOpUtils.h"
#include "ComputeKernelByCGAL.h"
#include "CGALUtils.h"

Mesh* computeKernelByCGAL(Mesh& hostMesh, double* kernelPoint) {

    vector<double*> halfSpaceCoeffs = computeHalfSpaceCoeffsFromTriangles(hostMesh.getAllTris(), hostMesh.getAllVerts());

    std::list<CGALPlane> planes;
    for (int i = 0; i < halfSpaceCoeffs.size(); i++) {
        double* coeffs = halfSpaceCoeffs[i];
        typename K::Plane_3 plane(coeffs[0], coeffs[1], coeffs[2], coeffs[3]);
        planes.push_back(plane);
        delete[] coeffs;
    }
    halfSpaceCoeffs.clear();
  
    CGALMesh chull;
    //CGAL::halfspace_intersection_with_constructions_3(planes.begin(), planes.end(), chull);
    //CGAL::halfspace_intersection_3(planes.begin(), planes.end(), chull);    // if no point inside the intersection is provided, one will be automatically found using linear programming
    CGAL::halfspace_intersection_3(planes.begin(), planes.end(), chull, CGALPoint(kernelPoint[0], kernelPoint[1], kernelPoint[2]));   
    //std::cout << "The convex hull contains " << num_vertices(chull) << " vertices" << std::endl;

    Mesh* kernel = convertCGALMeshToMesh(chull);
    return kernel;
}