// @author Merve Asiler

#pragma once

#include "Mesh.h"
#include "BasicGeometricElements.h"
#include <queue>

enum class NeighborPolarity {
	ONE_SIDED,
	DISTINCT_SIDED_INCLUSION,
	DISTINCT_SIDED_EXCLUSION,
	NO_POLARITY
};

struct Grid {
	int numOfCells[3];	// number of cells for each of the three dimension
	double cellSize[3];	// length of the edges for each of the three dimension, normally equal because of cube
	double minGridCoords[3], maxGridCoords[3];

	Grid() {};
	Grid(const Mesh& hostMesh, double* leastCoordinates, double* mostCoordinates, int* gridDimension, double cellSizeRatio);
	void defineGridByMeshCoordinates(const Mesh& hostMesh);
	void defineGridByExternalCoordinates(double* leastCoordinates, double* mostCoordinates);
	void partitionGridByCellSize(double cellSize);
	void partitionGridByDiagonalRatio(double ratio);
	void partitionGridByNumOfCells(int gridDimension[3]);
	int getNumOfCells();
	double getCellSize();
	double* getMinCoords();
	double* getMaxCoords();
	int* findHomeCell(double* point);
	double* findCellCoords(int i, int j, int k);
	
};

class KernelExpansion {

protected:
	double* initialPoint = nullptr;
	double* extremeCorners[2] = { nullptr, nullptr };
	const Mesh* hostMeshptr;
	Mesh kernel;

	Grid grid;

	queue<int*> possibleKernelCells;		// the cells which may include some pieces of the kernel
	vector<int*> handledCells;				// the cells which were previously checked or in queue to be checked for being inside kernel 
	vector<HalfSpace> halfSpaceSet;

	double** computeCellCoords(double initialCorner[3], double finalCorner[3]);

	bool isCellValidForQueue(int  neighbor_i, int neighbor_j, int neighbor_k);
	void decideOnNeighbor(int neigh_i, int neigh_j, int neigh_k, double neighVert1_scalar, double neighVert2_scalar, double neighVert3_scalar, double neighVert4_scalar);
	bool isInKernel(double* scalarVector);

public:
	KernelExpansion(const Mesh& hostMesh, double* extremeDirection);
	KernelExpansion(const Mesh& hostMesh);
	~KernelExpansion();
	Mesh& getKernel();
	double* getInitialKernelPoint();
	vector<HalfSpace>& getHalfSpaceSet();
	Grid getGrid();
	virtual void expandKernel() = 0;
	void checkKernelForNonKernelVertices();
	void checkKernelForIrregularTriangles();
};
