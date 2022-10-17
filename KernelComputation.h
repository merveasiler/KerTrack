// @author Merve Asiler

#pragma once

#include <string>
using namespace std;

class Mesh;

void doExperimentForPaper(string meshName);

void ComputeKernel(string meshName, string algoType);

/*
Mesh* ComputeKernelBySDFKer(Mesh* mesh, double cellSizeRatio);

Mesh* ComputeKernelBySDFKerPlus(Mesh* mesh, double cellSizeRatio, double* extremeDirection);

Mesh* ComputeKernelByMDFKer(Mesh* mesh, int gridDimension[3]);

Mesh* ComputeKernelByMDFKerPlus(Mesh* mesh, double cellSizeRatio);
*/

Mesh* ComputeKernelByKerTrack(Mesh& mesh);

Mesh* ComputeKernelByCGAL(Mesh& mesh, double* extremeDirection);

void FindKernelPoint_SDLP(string meshName);

void SphericalParametrize(string meshName);

void ShapeMorphByKernel(string sourceMeshName, string targetMeshName);

void ShapeMorphByLerp(string sourceMeshName, string targetMeshName);



