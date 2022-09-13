// @author Merve Asiler

#pragma once

#include "KernelExpansion.h"

class KernelExpansion_SDFKerPlus : public KernelExpansion {

	void computeScalarsForMarchingCube(double* scalarVectors[8], double coords[8][4]);
	void computeMarchingCubeCell(double coords[8][4], double* scalarVectors[8]);

public:
	KernelExpansion_SDFKerPlus(double* extremeDirection, const Mesh& hostMesh, Grid grid);
	void expandKernel();
};