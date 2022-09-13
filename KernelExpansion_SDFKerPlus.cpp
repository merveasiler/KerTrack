// @author Merve Asiler

#include "KernelExpansion_SDFKerPlus.h"
#include "BaseGeoOpUtils.h"
#include "MarchingCubes2.h"

KernelExpansion_SDFKerPlus::KernelExpansion_SDFKerPlus(double* extremeDirection, const Mesh& hostMesh, Grid grid) :
	KernelExpansion(extremeDirection, hostMesh, grid) {};

void KernelExpansion_SDFKerPlus::expandKernel() {

	// is kernel empty?
	if (this->initialPoint == nullptr) {
		delete kernel;
		kernel = NULL;
		return;	// NOT STAR-SHAPED!
	}

	// find the cell which the initial point belongs to and append it into the queue
	int* initialCell = grid.findHomeCell(initialPoint);
	possibleKernelCells.push(initialCell);	// this is a queue consisting of cells whose corners will be checked for being inside/outside kernel
	handledCells.push_back(initialCell);

	double p[3] = { initialPoint[0], initialPoint[1], initialPoint[2] };
	//cout << "Is this really a kernel point?: " << findClosestValueSatisfiedByPoint(p, halfSpaceSet) << endl;

	// check each cell in the queue through flood fill style
	while (!possibleKernelCells.empty()) {
		int* cell_indices = possibleKernelCells.front();
		possibleKernelCells.pop();

		// compute the 0th corner coordinates of the cell
		double initialCorner[3], finalCorner[3];
		for (int j = 0; j < 3; j++) {
			initialCorner[j] = grid.minGridCoords[j] + (cell_indices[j] * grid.cellSize[j]);
			finalCorner[j] = grid.cellSize[j] + initialCorner[j];
		}

		// ---> handle the corners in the ordering taken by first change in x, then in z, then in y 
		double coords[8][4] = { {initialCorner[0], initialCorner[1], initialCorner[2], 1.0},
								{finalCorner[0], initialCorner[1], initialCorner[2], 1.0},
								{finalCorner[0], initialCorner[1], finalCorner[2], 1.0},
								{initialCorner[0], initialCorner[1], finalCorner[2], 1.0},
								{initialCorner[0], finalCorner[1], initialCorner[2], 1.0},
								{finalCorner[0], finalCorner[1], initialCorner[2], 1.0},
								{finalCorner[0], finalCorner[1], finalCorner[2], 1.0},
								{initialCorner[0], finalCorner[1], finalCorner[2], 1.0}
		};
		double* scalarVectors[8];

		// for each corner compute the scalar value to be used in marching cubes by putting it into the system of inequalities
		for (int i = 0; i < 8; i++)
			scalarVectors[i] = findValueVectorSatisfiedByPoint(coords[i], halfSpaceSet);

		computeScalarsForMarchingCube(scalarVectors, coords);

		// compute the triangles inside the corresponding marching cube cell by using scalar values of cell corners obtained just above
		computeMarchingCubeCell(coords, scalarVectors);

		// Decide on neighbors for appending into queue
		int i = cell_indices[0], j = cell_indices[1], k = cell_indices[2];
		decideOnNeighbor(i - 1, j, k, coords[0][3], coords[3][3], coords[4][3], coords[7][3]);	// Neighbor at x-1
		decideOnNeighbor(i + 1, j, k, coords[1][3], coords[2][3], coords[5][3], coords[6][3]);	// Neighbor at x+1
		decideOnNeighbor(i, j - 1, k, coords[0][3], coords[1][3], coords[2][3], coords[3][3]);	// Neighbor at y-1
		decideOnNeighbor(i, j + 1, k, coords[4][3], coords[5][3], coords[6][3], coords[7][3]);	// Neighbor at y+1
		decideOnNeighbor(i, j, k - 1, coords[0][3], coords[1][3], coords[4][3], coords[5][3]);	// Neighbor at z-1
		decideOnNeighbor(i, j, k + 1, coords[2][3], coords[3][3], coords[6][3], coords[7][3]);	// Neighbor at z+1

	}

	for (int i = 0; i < kernel->getNumOfTris(); i++) {
		if (kernel->getTriangle(i)->corners[0] == kernel->getTriangle(i)->corners[1] ||
			kernel->getTriangle(i)->corners[1] == kernel->getTriangle(i)->corners[2] ||
			kernel->getTriangle(i)->corners[2] == kernel->getTriangle(i)->corners[0])
			cout << "errorrrrr: triangle no: " << i << endl;
	}

	int num_of_non_kernel_points = 0;
	for (int i = 0; i < kernel->getNumOfVerts(); i++) {
		if (findClosestValueSatisfiedByPoint(kernel->getVertex(i)->coords, halfSpaceSet) > 3 * EPSILON ||
			findClosestValueSatisfiedByPoint(kernel->getVertex(i)->coords, halfSpaceSet) > 3 * EPSILON ||
			findClosestValueSatisfiedByPoint(kernel->getVertex(i)->coords, halfSpaceSet) > 3 * EPSILON) {
			cout << "non-kernel point" << endl;
			num_of_non_kernel_points++;
		}
	}
	cout << "num of non-kernel points is: " << num_of_non_kernel_points << endl;

}

void KernelExpansion_SDFKerPlus::computeScalarsForMarchingCube(double* scalarVectors[8], double coords[8][4]) {

	for (int i = 0; i < 8; i++) {
		int j = 0;
		for (; j < halfSpaceSet.size(); j++) {
			if (scalarVectors[i][j] > 3 * EPSILON) {
				coords[i][3] = 1.0;	// store anything positive
				break;
			}
		}
		if (j == halfSpaceSet.size())
			coords[i][3] = -1.0;	// store anything negative

	}

}

void KernelExpansion_SDFKerPlus::computeMarchingCubeCell(double coords[8][4], double* scalarVectors[8]) {

	GRIDCELL2 gridCell;
	for (int c_id = 0; c_id < 8; c_id++) {
		XYZ2 xyz(coords[c_id], scalarVectors[c_id]);
		gridCell.p[c_id] = xyz;
		gridCell.val[c_id] = xyz.scalarValue;
		gridCell.vect[c_id] = xyz.scalarVector;
	}
	gridCell.vectSize = halfSpaceSet.size();

	TRIANGLE2 triangles[5];
	int number = Polygonise2(gridCell, 0, triangles);
	if (number > 0) {
		for (int t = 0; t < number; t++) { // triangles
			int corners[3];
			if (triangles[t].p[0] == triangles[t].p[1] || triangles[t].p[0] == triangles[t].p[2] || triangles[t].p[1] == triangles[t].p[2]) {
				//cout << "invalid triangle caught!" << endl;
				continue;	// invalid triangle
			}

			for (int v = 0; v < 3; v++) { // triangle vertices
				XYZ2 xyz = triangles[t].p[v];
				bool previously_saved_vertex = false;
				for (int pv = 0; pv < kernel->getNumOfVerts(); pv++) {
					Vertex* prevVertex = kernel->getVertex(pv);
					if (abs(xyz.x - (double)prevVertex->coords[0]) < 2 * EPSILON && abs(xyz.y - (double)prevVertex->coords[1]) < 2 * EPSILON && abs(xyz.z - (double)prevVertex->coords[2]) < 2 * EPSILON) {
						corners[v] = pv;
						previously_saved_vertex = true;
						break;
					}
				}
				if (!previously_saved_vertex) {
					corners[v] = kernel->getNumOfVerts();
					kernel->addVertex(xyz.x, xyz.y, xyz.z);
				}
			}

			kernel->addTriangle(corners[0], corners[1], corners[2]);
		}
	}
}