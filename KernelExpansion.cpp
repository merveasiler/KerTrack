// @author Merve Asiler

#include "KernelExpansion.h"
#include "BaseGeoOpUtils.h"
#include "sdlp.h"

KernelExpansion::KernelExpansion(const Mesh& hostMesh, double* extremeDirection) {

	this->hostMeshptr = &hostMesh;

	computeHalfSpacesFromTriangles(hostMeshptr->getAllTris(), hostMeshptr->getAllVerts(), this->halfSpaceSet);

	if (extremeDirection) {
		// FIND POINT BY USING SINGLE EXTREME DIRECTION
		this->extremeCorners[0] = nullptr;	this->extremeCorners[1] = nullptr;
		this->initialPoint = sdlpMain(extremeDirection, halfSpaceSet);	// compute initial kernel point
	}
	else {
		// FIND POINT BY USING ALL OF THE 6 EXTREME DIRECTIONS, TAKE THEIR AVERAGE
		extremeCorners[0] = new double[3];	extremeCorners[1] = new double[3];
		double extremeDirections[6][3] = { {-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1} };
		for (int i = 0; i < 6; i++) {
			this->initialPoint = sdlpMain(extremeDirections[i], halfSpaceSet);	// compute initial kernel point at the given extreme direction
			if (!this->initialPoint)
				return;
			this->extremeCorners[i % 2][int(i / 2)] = this->initialPoint[int(i / 2)];
			delete[] this->initialPoint;
		}
		this->initialPoint = nullptr;
	}

}

KernelExpansion::KernelExpansion(const Mesh& hostMesh) {

	this->hostMeshptr = &hostMesh;
	this->initialPoint = nullptr;
	this->extremeCorners[0] = nullptr;	this->extremeCorners[1] = nullptr;
	//computeHalfSpacesFromTriangles(hostMesh.getAllTris(), hostMesh.getAllVerts(), this->halfSpaceSet);

}

KernelExpansion::~KernelExpansion() {
	if (initialPoint) {
		delete[] initialPoint;
		initialPoint = nullptr;
	}

	if (extremeCorners[0]) {
		delete[] extremeCorners[0];
		extremeCorners[0] = nullptr;
		delete[] extremeCorners[1];
		extremeCorners[1] = nullptr;
	}
	
	halfSpaceSet.clear();

	for (int i = 0; i < handledCells.size(); i++) {
		delete[] handledCells[i];
		handledCells[i] = nullptr;
	}
	handledCells.clear();
}

Mesh& KernelExpansion::getKernel() {
	return kernel;
}

double* KernelExpansion::getInitialKernelPoint() {
	return initialPoint;
}

vector<HalfSpace>& KernelExpansion::getHalfSpaceSet() {
	return halfSpaceSet;
}

Grid KernelExpansion::getGrid() {
	return grid;
}

double** KernelExpansion::computeCellCoords(double initialCorner[3], double finalCorner[3]) {

	double** coords = new double* [8];
	// ---> handle the corners in the ordering taken by first change in x, then in z, then in y 
	coords[0] = new double[3]{ initialCorner[0], initialCorner[1], initialCorner[2] };
	coords[1] = new double[3]{ finalCorner[0], initialCorner[1], initialCorner[2] };
	coords[2] = new double[3]{ finalCorner[0], initialCorner[1], finalCorner[2] };
	coords[3] = new double[3]{ initialCorner[0], initialCorner[1], finalCorner[2] };
	coords[4] = new double[3]{ initialCorner[0], finalCorner[1], initialCorner[2] };
	coords[5] = new double[3]{ finalCorner[0], finalCorner[1], initialCorner[2] };
	coords[6] = new double[3]{ finalCorner[0], finalCorner[1], finalCorner[2] };
	coords[7] = new double[3]{ initialCorner[0], finalCorner[1], finalCorner[2] };

	return coords;
}



void KernelExpansion::decideOnNeighbor(int neigh_i, int neigh_j, int neigh_k, double neighVert1_scalar, double neighVert2_scalar, double neighVert3_scalar, double neighVert4_scalar) {

	if (isCellValidForQueue(neigh_i, neigh_j, neigh_k)) {

		// if any vertex of the given face is a kernel point, then we should check this neighbor cell sharing the face
		if (neighVert1_scalar <= 0 ||
			neighVert2_scalar <= 0 ||
			neighVert3_scalar <= 0 ||
			neighVert4_scalar <= 0) {	// scalarValue is located.

			int* neigh_ijk = new int[3]{ neigh_i, neigh_j, neigh_k };
			possibleKernelCells.push(neigh_ijk);
			handledCells.push_back(neigh_ijk);
		}
	}
}

bool KernelExpansion::isCellValidForQueue(int  neighbor_i, int neighbor_j, int neighbor_k) {

	int neighbor_ijk[3] = { neighbor_i, neighbor_j, neighbor_k };

	for (int i = 0; i < 3; i++)
		if (neighbor_ijk[i] >= grid.numOfCells[i] || neighbor_ijk[i] < 0)	// there is no such a neighbor because of invalid index
			return false;

	// did we previously check this cell?
	bool is_found = false;
	for (int i = 0; i < handledCells.size(); i++) {
		is_found = true;
		for (int j = 0; j < 3; j++) {
			if (handledCells[i][j] != neighbor_ijk[j]) {
				is_found = false;
				break;
			}
		}
		if (is_found)
			return false;
	}

	return true;
}

bool KernelExpansion::isInKernel(double* scalarVector) {
	for (int i = 0; i < halfSpaceSet.size(); i++)
		if (scalarVector[i] > EPSILON)
			return false;	// out
	return true;		// in
}

void KernelExpansion::checkKernelForNonKernelVertices() {

	int num_of_non_kernel_points = 0;
	for (int i = 0; i < kernel.getNumOfVerts(); i++) {
		if (findClosestValueSatisfiedByPoint(kernel.getVertex(i).coords, halfSpaceSet) > 3 * EPSILON ||
			findClosestValueSatisfiedByPoint(kernel.getVertex(i).coords, halfSpaceSet) > 3 * EPSILON ||
			findClosestValueSatisfiedByPoint(kernel.getVertex(i).coords, halfSpaceSet) > 3 * EPSILON) {
			cout << "non-kernel point" << endl;
			num_of_non_kernel_points++;
		}
	}
	cout << "num of non-kernel points is: " << num_of_non_kernel_points << endl;

}

void KernelExpansion::checkKernelForIrregularTriangles() {

	int num_of_irregular_triangles = 0;
	for (int i = 0; i < kernel.getNumOfTris(); i++) {
		if (kernel.getTriangle(i).corners[0] == kernel.getTriangle(i).corners[1] ||
			kernel.getTriangle(i).corners[1] == kernel.getTriangle(i).corners[2] ||
			kernel.getTriangle(i).corners[2] == kernel.getTriangle(i).corners[0]) {
			cout << "errorrrrr: triangle no: " << i << endl;
			num_of_irregular_triangles++;
		}
	}
	cout << "num of irregular triangles is: " << num_of_irregular_triangles << endl;

}



Grid::Grid(const Mesh& hostMesh, double* leastCoordinates, double* mostCoordinates, int* gridDimension, double cellSizeRatio) {

	if (leastCoordinates == NULL)
		this->defineGridByMeshCoordinates(hostMesh);								// USE MESH BOUNDARY BOX TO DEFINE GRID
	else
		this->defineGridByExternalCoordinates(leastCoordinates, mostCoordinates);	// USE THE GIVEN BOUNDARY BOX TO DEFINE GRID

	if (gridDimension)
		this->partitionGridByNumOfCells(gridDimension);
	else
		this->partitionGridByDiagonalRatio(cellSizeRatio);
}

void Grid::defineGridByMeshCoordinates(const Mesh& hostMesh) {

	for (int j = 0; j < 3; j++) {
		this->minGridCoords[j] = numeric_limits<double>::infinity();
		this->maxGridCoords[j] = -numeric_limits<double>::infinity();
	}

	// find the "most" and "least" coordinates of the mesh
	for (int i = 0; i < hostMesh.getNumOfVerts(); i++) {
		Vertex vertex = hostMesh.getVertex(i);
		for (int j = 0; j < 3; j++) {
			if (this->minGridCoords[j] > vertex.coords[j])
				this->minGridCoords[j] = vertex.coords[j];
			if (this->maxGridCoords[j] < vertex.coords[j])
				this->maxGridCoords[j] = vertex.coords[j];
		}
	}

}

void Grid::defineGridByExternalCoordinates(double* leastCoordinates, double* mostCoordinates) {

	for (int j = 0; j < 3; j++) {
		this->minGridCoords[j] = leastCoordinates[j];
		this->maxGridCoords[j] = mostCoordinates[j];
	}

}

void Grid::partitionGridByCellSize(double cellSize) {

	for (int j = 0; j < 3; j++)
		this->cellSize[j] = cellSize;

	// compute number of cells of grid at each direction
	for (int j = 0; j < 3; j++) {
		this->numOfCells[j] = (this->maxGridCoords[j] - this->minGridCoords[j]) / this->cellSize[j];
		//this->maxGridCoords[j] = this->minGridCoords[j] + (this->cellSize[j] * ++this->numOfCells[j]);
	}

	/*
	// CHECK IF DEFINED POINT IS A VALID KERNEL POINT
	cout << "MAX Coords: " << maxGridCoords[0] << " " << maxGridCoords[1] << " " << maxGridCoords[2] << endl;
	cout << "MIN Coords: " << minGridCoords[0] << " " << minGridCoords[1] << " " << minGridCoords[2] << endl;
	cout << "DIFFERENCE: ";
	for (int j = 0; j < 3; j++)
		cout << (maxGridCoords[j] - minGridCoords[j]) << " ";
	cout << endl;
	*/

}

void Grid::partitionGridByDiagonalRatio(double ratio) {

	// compute each edge length of the grid
	double edgeLength[3];
	for (int j = 0; j < 3; j++)
		edgeLength[j] = this->maxGridCoords[j] - this->minGridCoords[j];

	// compute diagonal's length & cellSize
	double diagonal = sqrt((edgeLength[0] * edgeLength[0]) + (edgeLength[1] * edgeLength[1]) + (edgeLength[2] * edgeLength[2]));
	double cellSize = diagonal * ratio;

	partitionGridByCellSize(cellSize);

}

void Grid::partitionGridByNumOfCells(int gridDimension[3]) {

	for (int d = 0; d < 3; d++) {
		this->numOfCells[d] = gridDimension[d];
		this->cellSize[d] = (this->maxGridCoords[d] - this->minGridCoords[d]) / this->numOfCells[d];
	}

}

int Grid::getNumOfCells() {

	return numOfCells[0] * numOfCells[1] * numOfCells[2];
}

double Grid::getCellSize() {

	return cellSize[0];
}

double* Grid::getMinCoords() {

	return minGridCoords;
}

double* Grid::getMaxCoords() {

	return maxGridCoords;
}

int* Grid::findHomeCell(double* point) {

	int* cell_indices = new int[3];

	for (int d = 0; d < 3; d++) {	// d: dimension
		double coordValue = minGridCoords[d];
		for (int i = 0; i < numOfCells[d] + 1; i++) {
			if (point[d] < coordValue) {
				cell_indices[d] = i - 1;
				break;
			}
			coordValue += cellSize[d];
		}
	}

	return cell_indices;
}

double* Grid::findCellCoords(int i, int j, int k) {

	int cell_indices[3] = { i, j, k };
	double* cornerCoords = new double[3];
	for (int c = 0; c < 3; c++)
		cornerCoords[c] = minGridCoords[c] + (cell_indices[c] * cellSize[c]);

	return cornerCoords;
}

