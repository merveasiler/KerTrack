// @author Merve Asiler

#pragma once
#pragma comment(lib, "boost_filesystem-vc140-mt.lib")

#include <string>
using namespace std;

class Mesh;

void doExperimentForPaper(string meshName);

void ComputeKernel(string meshName, string algoType);

void ComputeBatchKernel(string inputFolderName, string outputFolderName, string algoType);

Mesh ComputeKernelByKerTrack(Mesh& mesh);

Mesh ComputeKernelByCGAL(Mesh& mesh, double* extremeDirection);

void FindKernelPoint_SDLP(string meshName);

void SphericalParametrize(string meshName);

void ShapeMorphByKernel(string sourceMeshName, string targetMeshName);

void ShapeMorphByLerp(string sourceMeshName, string targetMeshName);



