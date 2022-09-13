// @author Merve Asiler

#pragma once

#include "KernelExpansion.h"

class KernelExpansion_SDFKer : public KernelExpansion {

	double evaluateCorner(int i, int j, int k);
	void evaluateCell(double** backFaceArray, double** frontFaceArray, int i, int j, int k);
	void computeMarchingCubeCell(double coords[8][4]);

public:
	KernelExpansion_SDFKer(const Mesh& hostMesh, double cellSizeRatio);
	void expandKernel();
};