// @author Merve Asiler

#pragma once

#include "KernelExpansion.h"

class KernelExpansion_MDFKerPlus : public KernelExpansion {

	vector<int*> allTriedPlaneCombinations;

	void sendRay(double source[3], double directions[3]);
	void findPossibleKernelVertex(double* scalars, double source_coords[3]);
	void updateClosest3Planes(double* values, int currentPlaneIds[3], int previousPlaneIds[3]);
	double* computePossibleKernelVertex(double* scalarsVector, int closestPlaneIds[3]);
	double* computePossibleKernelVertex(int closestPlaneIds[3]);
	int addIfNoFrontierObstacle(double candidateKernelVertex[3], int candidateVertexProducerIds[3], double source_coords[3]);
	bool isAlreadyKernelVertex(double coords[3]);
	int findIfFrontierObstacle(double* candidateKernelVertex, int candidateVertexProducerIds[3], double* sourcePoint);
	void findByFrontierPlanes(int closestPlaneIds[3], int obstaclePlaneId, double source_coords[3]);
	void replaceOnePlaneBy(queue<int*>& possiblePlaneCombinations, int planeIds[3], int replacerId);

	void decideOnNeighborCells(int* cell_indices, double* corner_in_outs, double** scalarsVector);
	void addNeighborCellIntoQueue(int neigh_i, int neigh_j, int neigh_k,
		double neighVert1_scalar, double neighVert2_scalar, double neighVert3_scalar, double neighVert4_scalar,
		NeighborPolarity np1, NeighborPolarity np2, NeighborPolarity np3, NeighborPolarity np4);

public:
	KernelExpansion_MDFKerPlus(const Mesh& hostMesh, double cellSizeRatio);
	~KernelExpansion_MDFKerPlus();
	void expandKernel();
	void findMissingKernelVertices(Mesh* incompleteKernel, double cellSize);
};