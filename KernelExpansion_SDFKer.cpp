// @author Merve Asiler

#include "KernelExpansion_SDFKer.h"
#include "BaseGeoOpUtils.h"
#include "MarchingCubes.h"

KernelExpansion_SDFKer::KernelExpansion_SDFKer(const Mesh& hostMesh, double cellSizeRatio) :
	KernelExpansion(hostMesh, cellSizeRatio) {};

void KernelExpansion_SDFKer::expandKernel() {

	// is kernel empty?
	if (this->extremeCorners[0] == nullptr && this->extremeCorners[1] == nullptr) {
		delete kernel;
		kernel = NULL;
		return;	// NOT STAR-SHAPED!
	}

	// reserve a 2D array of size <grid.numOfCells[1] x grid.numOfCells[2]> for back face and front face of a cell separately
	double** backFaceArray = new double* [grid.numOfCells[1]], ** frontFaceArray = new double* [grid.numOfCells[1]];
	for (int j = 0; j < grid.numOfCells[1]; j++) {
		backFaceArray[j] = new double[grid.numOfCells[2]];
		frontFaceArray[j] = new double[grid.numOfCells[2]];
	}

	// initialize back face of the initial cell
	// ... which means; when i = 0: 
	for (int j = 0; j < grid.numOfCells[1]; j++) {
		for (int k = 0; k < grid.numOfCells[2]; k++)
			backFaceArray[j][k] = evaluateCorner(0, j, k);
	}

	// now; when i > 0: 
	for (int i = 1; i < grid.numOfCells[0]; i++) {

		// when j = 0:
		for (int k = 0; k < grid.numOfCells[2]; k++)
			frontFaceArray[0][k] = evaluateCorner(i, 0, k);
		// when j > 0:
		for (int j = 1; j < grid.numOfCells[1]; j++) {

			// when k = 0:
			frontFaceArray[j][0] = evaluateCorner(i, j, 0);
			// when k > 0:
			for (int k = 1; k < grid.numOfCells[2]; k++) {
				frontFaceArray[j][k] = evaluateCorner(i, j, k);
				evaluateCell(backFaceArray, frontFaceArray, i, j, k);
			}
		}

		// transfer from front face to back face, initialize front face
		for (int j = 0; j < grid.numOfCells[1]; j++) {
			for (int k = 0; k < grid.numOfCells[2]; k++)
				backFaceArray[j][k] = frontFaceArray[j][k];
		}

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

	// Clean-up
	for (int j = 0; j < grid.numOfCells[1]; j++) {
		delete[] frontFaceArray[j];
		delete[] backFaceArray[j];
	}

	delete[] frontFaceArray;
	delete[] backFaceArray;
}

double KernelExpansion_SDFKer::evaluateCorner(int i, int j, int k) {

	// compute the corner coordinates
	int cell_indices[3] = { i, j, k };
	double cornerCoords[3];
	for (int c = 0; c < 3; c++)
		cornerCoords[c] = grid.minGridCoords[c] + (cell_indices[c] * grid.cellSize[c]);

	// for the corner compute the scalar value to be used by putting it into the system of inequalities
	double scalarValue = findClosestValueSatisfiedByPoint(cornerCoords, halfSpaceSet);

	return scalarValue;
}

void KernelExpansion_SDFKer::evaluateCell(double** backFaceArray, double** frontFaceArray, int i, int j, int k) {

	// compute the 0th corner coordinates of the cell
	int cell_indices[3] = { i - 1, j - 1, k - 1 };
	double initialCorner[3], finalCorner[3];
	for (int d = 0; d < 3; d++) {
		initialCorner[d] = grid.minGridCoords[d] + (cell_indices[d] * grid.cellSize[d]);
		finalCorner[d] = grid.cellSize[d] + initialCorner[d];
	}

	// ---> handle the corners in the ordering taken by first change in x, then in z, then in y 
	double coords[8][4] = { {initialCorner[0], initialCorner[1], initialCorner[2], backFaceArray[j - 1][k - 1]},
							{finalCorner[0], initialCorner[1], initialCorner[2], frontFaceArray[j - 1][k - 1]},
							{finalCorner[0], initialCorner[1], finalCorner[2], frontFaceArray[j - 1][k]},
							{initialCorner[0], initialCorner[1], finalCorner[2], backFaceArray[j - 1][k]},
							{initialCorner[0], finalCorner[1], initialCorner[2], backFaceArray[j][k - 1]},
							{finalCorner[0], finalCorner[1], initialCorner[2], frontFaceArray[j][k - 1]},
							{finalCorner[0], finalCorner[1], finalCorner[2], frontFaceArray[j][k]},
							{initialCorner[0], finalCorner[1], finalCorner[2], backFaceArray[j][k]}
	};
	// compute the triangles inside the corresponding marching cube cell by using scalar values of cell corners obtained just above
	computeMarchingCubeCell(coords);

}

void KernelExpansion_SDFKer::computeMarchingCubeCell(double coords[8][4]) {

	GRIDCELL gridCell;
	for (int c_id = 0; c_id < 8; c_id++) {
		XYZ xyz(coords[c_id]);
		gridCell.p[c_id] = xyz;
		gridCell.val[c_id] = xyz.scalarValue;
	}

	TRIANGLE triangles[5];
	int number = Polygonise(gridCell, 0, triangles);
	if (number > 0) {
		for (int t = 0; t < number; t++) { // triangles
			int corners[3];
			if (triangles[t].p[0] == triangles[t].p[1] || triangles[t].p[0] == triangles[t].p[2] || triangles[t].p[1] == triangles[t].p[2]) {
				//cout << "invalid triangle caught!" << endl;
				continue;	// invalid triangle
			}

			for (int v = 0; v < 3; v++) { // triangle vertices
				XYZ xyz = triangles[t].p[v];
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

