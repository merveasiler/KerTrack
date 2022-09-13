// @author Merve Asiler

#include "KernelComputation.h"
#include "KernelExpansion_SDFKer.h"
#include "KernelExpansion_SDFKerPlus.h"
#include "KernelExpansion_MDFKer.h"
#include "KernelExpansion_MDFKerPlus.h"
#include "KernelExpansion_KerTrack.h"
#include "SomehowKernel.h"
#include "ComputeKernelByCGAL.h"
#include "sdlp.h"
#include "SceneManager.h"
#include "CommonUtils.h"
#include "CGALUtils.h"

#include <ctime>

/* ***************************** GLOBAL VARIABLES ***************************** */
/* */ vector<double> produceColorSource(Mesh* ground_truth, Mesh* exp_mesh);
/* */ string fileNames[9] = {"321.off", "325.off", "326.off", "350.off", "cube.off", "cube_broken.off", "liver.obj", "Rock_2.obj", "Rock_6.obj"};
/* */ double extremeDirections[6][3] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
/* */ double cellSizeRatios[9][4] = {{0.07, 0.05, 0.035, 0.01}, {0.07, 0.05, 0.035, 0.01}, {0.07, 0.05, 0.035, 0.01}, {0.01, 0.005, 0.0045, 0.0025},
/* */								{0.07, 0.05, 0.035, 0.01}, {0.07, 0.05, 0.035, 0.01}, {0.07, 0.035, 0.025, 0.01}, {0.07, 0.035, 0.025, 0.01}, {0.07, 0.035, 0.025, 0.01}};
/* */ double executionCount = 10;
/* **************************************************************************** */

void doExperimentForPaper(string meshName) {

	/************************************************ READ MESH ************************************************/
	Mesh* mesh = new Mesh();
	if (meshName.substr(meshName.length() - 3, 3) == "off")
		mesh->loadOff(meshName.c_str());
	else
		mesh->loadObj(meshName.c_str());

	/*************************************** INGREDIENTS / PREPARATIONS  ***************************************/
		// ... for kernel computation
	double extremeDirection[3] = { 0, 0, 1 };
	vector<double> cellSizeRatios{ 0.075}; // the ideal scale is the 2nd one 
	vector<Mesh*> kernels;
	kernels.resize(cellSizeRatios.size() * 2 + 1);	// number of experiments: 1 + 4 + 4
		// ... for positions of the shapes on the scene
	vector<double*> positions;
	Grid boundingBox(mesh);	// compute the bounding box of the shape
	double* minborder = boundingBox.getMinCoords();
	double* maxborder = boundingBox.getMaxCoords();
	double objectWidth[3], breakWidth[3];
	for (int j = 0; j < 3; j++) {
		objectWidth[j] = maxborder[j] - minborder[j];
		breakWidth[j] = objectWidth[j] / 2.0;
	}
	double sceneSizeUnit = objectWidth[0] > objectWidth[1] ? objectWidth[0] : objectWidth[1];
	double totalSceneSize = (sceneSizeUnit + sceneSizeUnit / 10.0) * (kernels.size() / 2);
		// ... for shape materials
	vector<tuple<tuple<Mesh*, MaterialSetting*>, tuple<Mesh*, MaterialSetting*>>> outputs;
	MaterialSetting* kernelMatSetting = new MaterialSetting(0, 0, 1, 0);
	MaterialSetting* meshMatSetting = new MaterialSetting(1, 1, 1, 0.5);

	/************************************************ ALGORITHMS ************************************************/
		// ... CGAL
	Mesh* groundTruth = kernels[0] = ComputeKernelByCGAL(mesh, extremeDirection);
	positions.push_back(new double[3]{objectWidth[0] + breakWidth[0], 0, 0 });	// on the corner scene
	outputs.push_back(make_tuple(make_tuple(kernels[0], kernelMatSetting), make_tuple(mesh, meshMatSetting)));
	cout << endl;

		// ... KERTRACK
	for (int j = 1, i = 0; i < kernels.size()/2; i++, j++) {
		kernels[j] = ComputeKernelByKerTrack(mesh);
		double* position = new double[3]{ 0, 0, 0 };	// on the lower row
		position[0] = i * (objectWidth[0] + breakWidth[0]);
		positions.push_back(position);
		outputs.push_back(make_tuple(make_tuple(kernels[j], kernelMatSetting), make_tuple(mesh, meshMatSetting)));
		cout << endl;
	}

		// ... SDF-KER
	for (int j = 0, i = 1 + kernels.size() / 2; i < kernels.size(); i++, j++) {
		kernels[i] = ComputeKernelBySDFKer(mesh, cellSizeRatios[j]);
		double* position = new double[3]{ 0, objectWidth[1] + breakWidth[1], 0 };	// on the upper row
		position[0] = j * (objectWidth[0] + breakWidth[0]);
		positions.push_back(position);
		outputs.push_back(make_tuple(make_tuple(kernels[i], kernelMatSetting), make_tuple(mesh, meshMatSetting)));
		cout << endl;
	}

	if (groundTruth) {

		// POLYHEDRON-KERNEL (ITALIAN JOB)
		Mesh* polyhedronKernel = new Mesh();
		polyhedronKernel->loadOff((meshName.substr(meshName.length() - 3, 3) + "_kernel.off").c_str());
		kernels.push_back(polyhedronKernel);
		positions.push_back(new double[3]{ -1, -1, -1 });	// not to be used
		cout << endl;

		/*********************************** COMPARE HAUSDORFF DISTANCES & VOLUMES ***********************************/
		double groundTruthVolume = groundTruth->computeVolume();
		for (int j = 0, mod = (kernels.size() - 1) / 2, i = 1; i < kernels.size(); i++, j++) {
			double volume = kernels[i]->computeVolume();
			if (volume < EPSILON)
				continue;
			cout << "Volume of the kernel of " << cellSizeRatios[j % mod] << ": " << volume << " out of " << groundTruthVolume << "." << endl;
			cout << "Volume difference between the one <" << cellSizeRatios[j % mod] << "> and ground truth (CGAL): " << groundTruthVolume - volume << endl;

			double* hd = computeHausdorffDistance(groundTruth, kernels[i]);
			cout << "Hausdorff distance between the one <" << cellSizeRatios[j % mod] << "> and ground truth (CGAL): <" << hd[0] << ", " << hd[1] << ", " << hd[2] << ">." << endl << endl;
			delete[] hd;
		}

		/*************************************************** DRAW ***************************************************/
		vector<tuple<tuple<Mesh*, MaterialSetting*>, tuple<Mesh*, MaterialSetting*>>> outputs_to_draw;
		vector<double*> positions_to_draw;
		for (int i = 0; i < cellSizeRatios.size() + 1; i++) {
			outputs_to_draw.push_back(outputs[i]);
			positions_to_draw.push_back(positions[i]);
		}
		drawMultipleScenes(outputs_to_draw, positions_to_draw, totalSceneSize);
	}

	/************************************************* CLEAN-UP *************************************************/
	for (int i = 0; i < kernels.size(); i++) {
		delete kernels[i];
		delete[] positions[i];
	}
	delete kernelMatSetting;
	delete meshMatSetting;
	delete mesh;

}

void ComputeKernel(string meshName, string algoType) {

	// Read mesh
	Mesh* mesh = new Mesh();
	if (meshName.substr(meshName.length() - 3, 3) == "off")
		mesh->loadOff(meshName.c_str());
	else
		mesh->loadObj(meshName.c_str());

	// Ingredients:
	double extremeDirection[3] = { 0, -1, 0 };
	double cellSizeRatio = 0.07;
	int gridDimension[3] = {2, 2, 2};
	Mesh* kernel = nullptr;

	// Compute kernel
	if (algoType == "sdfker")
		kernel = ComputeKernelBySDFKer(mesh, cellSizeRatio);
	else if (algoType == "sdfker_plus")
		kernel = ComputeKernelBySDFKerPlus(mesh, cellSizeRatio, extremeDirection);
	else if (algoType == "mdfker")
		kernel = ComputeKernelByMDFKer(mesh, gridDimension);
	else if (algoType == "mdfker_plus")
		kernel = ComputeKernelByMDFKerPlus(mesh, cellSizeRatio);
	else if (algoType == "kertrack")
		kernel = ComputeKernelByKerTrack(mesh);
	else if (algoType == "kernel_by_cgal")
		kernel = ComputeKernelByCGAL(mesh, extremeDirection);
	else;

	// Draw
	if (kernel) {
		MaterialSetting* kernelMatSetting = new MaterialSetting(0, 0, 1, 0);
		MaterialSetting* meshMatSetting = new MaterialSetting(1, 1, 1, 0.5);
		vector<tuple<Mesh*, MaterialSetting*>> mesh_mat_set = { make_tuple(kernel, kernelMatSetting), make_tuple(mesh, meshMatSetting) };
		drawMultipleMeshToScene(mesh_mat_set);
		delete kernel;
	}
	delete mesh;

}

Mesh* ComputeKernelBySDFKer(Mesh* mesh, double cellSizeRatio) {

	// Execute by <executionCount>-many times
	double totalTime = 0;
	for (int i = 0; i < executionCount; i++) {
		clock_t begin = clock();
		KernelExpansion* kernelExpansion = new KernelExpansion_SDFKer(*mesh, cellSizeRatio);
		kernelExpansion->expandKernel();
		Mesh* kernel = kernelExpansion->getKernel();
		clock_t end = clock();
		totalTime += double(end - begin) / CLOCKS_PER_SEC;

		if (kernel)
			delete kernel;
		delete kernelExpansion;
	}

	// Execute to see the results
	KernelExpansion* kernelExpansion = new KernelExpansion_SDFKer(*mesh, cellSizeRatio);
	kernelExpansion->expandKernel();
	Mesh* kernel = kernelExpansion->getKernel();
	if (kernel) {
		double* kernelPoint = kernelExpansion->getInitialKernelPoint();
		cout << "Used kernel point: " << kernelPoint[0] << " " << kernelPoint[1] << " " << kernelPoint[2] << endl;
		cout << "Cell size is " << cellSizeRatio << " of the grid's diagonal: " << kernelExpansion->getGrid().getCellSize() << endl;
		cout << "Number of cells in grid: " << kernelExpansion->getGrid().getNumOfCells() << endl;
		cout << "Mesh: [faces: " << mesh->getNumOfTris() << "], [edges: " << mesh->getNumOfEdges() << "], [vertices: " << mesh->getNumOfVerts() << "]" << endl;
		delete kernel;
	}
	else {
		kernel = NULL;
		cout << "Kernel is empty!" << endl;
		cout << "Mesh: [faces: " << mesh->getNumOfTris() << "]\n";
	}

	// Print the average time
	double elapsed_secs = totalTime / executionCount;
	cout << "Kernel computation has been completed in " << elapsed_secs << " second(s) by SDF-KER-PLUS with cell size ratio = " << cellSizeRatio << "." << endl;
	delete kernelExpansion;
	return kernel;

}

Mesh* ComputeKernelBySDFKerPlus(Mesh* mesh, double cellSizeRatio, double* extremeDirection) {

	// Execute by <executionCount>-many times
	double totalTime = 0;
	for (int i = 0; i < executionCount; i++) {
		clock_t begin = clock();
		Grid grid(mesh);
		grid.constructGridByDiagonalRatio(cellSizeRatio);
		KernelExpansion* kernelExpansion = new KernelExpansion_SDFKerPlus(extremeDirection, *mesh, grid);
		kernelExpansion->expandKernel();
		Mesh* kernel = kernelExpansion->getKernel();
		clock_t end = clock();
		totalTime += double(end - begin) / CLOCKS_PER_SEC;

		if (kernel)
			delete kernel;
		delete kernelExpansion;
	}

	// Execute to see the results
	Grid grid(mesh);
	grid.constructGridByDiagonalRatio(cellSizeRatio);
	KernelExpansion* kernelExpansion = new KernelExpansion_SDFKerPlus(extremeDirection, *mesh, grid);
	kernelExpansion->expandKernel();
	Mesh* kernel = kernelExpansion->getKernel();
	if (kernel) {
		double* kernelPoint = kernelExpansion->getInitialKernelPoint();
		cout << "Used kernel point: " << kernelPoint[0] << " " << kernelPoint[1] << " " << kernelPoint[2] << endl;
		cout << "Cell size is " << cellSizeRatio << " of the grid's diagonal: " << grid.getCellSize() << endl;
		cout << "Number of cells in grid: " << grid.getNumOfCells() << endl;
		cout << "Mesh: [faces: " << mesh->getNumOfTris() << "], [edges: " << mesh->getNumOfEdges() << "], [vertices: " << mesh->getNumOfVerts() << "]" << endl;
		delete kernel;
	}
	else {
		kernel = NULL;
		cout << "Kernel is empty!" << endl;
		cout << "Mesh: [faces: " << mesh->getNumOfTris() << "]\n";
	}

	// Print the average time
	double elapsed_secs = totalTime / executionCount;
	cout << "Kernel computation has been completed in " << elapsed_secs << " second(s)." << endl;
	delete kernelExpansion;
	return kernel;

}

Mesh* ComputeKernelByMDFKer(Mesh* mesh, int gridDimension[3]) {

	// Execute by <executionCount>-many times
	double totalTime = 0;
	for (int i = 0; i < executionCount; i++) {
		clock_t begin = clock();
		Grid grid(mesh);
		KernelExpansion* kernelExpansion = new KernelExpansion_MDFKer(*mesh, gridDimension);
		kernelExpansion->expandKernel();
		Mesh* kernel = NULL, *incompleteKernel = kernelExpansion->getKernel();
		if (incompleteKernel)
			kernel = computeConvexHull(incompleteKernel->getAllVerts());
		clock_t end = clock();
		totalTime += double(end - begin) / CLOCKS_PER_SEC;

		if (incompleteKernel) {
			delete incompleteKernel;
			delete kernel;
		}
		delete kernelExpansion;
	}

	// Execute to see the results
	Grid grid(mesh);
	KernelExpansion* kernelExpansion = new KernelExpansion_MDFKer(*mesh, gridDimension);
	kernelExpansion->expandKernel();
	Mesh* kernel, *incompleteKernel = kernelExpansion->getKernel();
	if (incompleteKernel) {
		kernel = computeConvexHull(incompleteKernel->getAllVerts());
		double* kernelPoint = kernelExpansion->getInitialKernelPoint();
		cout << "Used kernel point: " << kernelPoint[0] << " " << kernelPoint[1] << " " << kernelPoint[2] << endl;
		cout << "Cell dimension is: " << gridDimension[0] << " x " << gridDimension[1] << " x " << gridDimension[2] << ".\n";
		cout << "Number of found kernel vertices: " << incompleteKernel->getNumOfVerts() << "." << endl;
		cout << "Mesh: [faces: " << mesh->getNumOfTris() << "], [edges: " << mesh->getNumOfEdges() << "], [vertices: " << mesh->getNumOfVerts() << "]" << endl;
		cout << "Kernel: [faces: " << kernel->getNumOfTris() << "], [edges: " << kernel->getNumOfEdges() << "], [vertices: " << kernel->getNumOfVerts() << "]" << endl;
		delete incompleteKernel;
	}
	else {
		kernel = NULL;
		cout << "Kernel is empty!" << endl;
		cout << "Mesh: [faces: " << mesh->getNumOfTris() << "]\n";
	}
	
	// Print the average time
	double elapsed_secs = totalTime / executionCount;
	cout << "Kernel computation has been completed in " << elapsed_secs << " second(s) by MDF-KER." << endl;
	delete kernelExpansion;
	return kernel;
}

Mesh* ComputeKernelByMDFKerPlus(Mesh* mesh, double cellSizeRatio) {

	// Execute by <executionCount>-many times
	double totalTime = 0;
	for (int i = 0; i < executionCount; i++) {
		clock_t begin = clock();
		KernelExpansion* kernelExpansion = new KernelExpansion_MDFKerPlus(*mesh, cellSizeRatio);
		kernelExpansion->expandKernel();
		Mesh* kernel = NULL, *incompleteKernel = kernelExpansion->getKernel();
		if (incompleteKernel)
			kernel = computeConvexHull(incompleteKernel->getAllVerts());
		clock_t end = clock();
		totalTime += double(end - begin) / CLOCKS_PER_SEC;

		if (incompleteKernel) {
			delete incompleteKernel;
			delete kernel;
		}
		delete kernelExpansion;
	}

	// Execute to see the results
	KernelExpansion* kernelExpansion = new KernelExpansion_MDFKerPlus(*mesh, cellSizeRatio);
	kernelExpansion->expandKernel();
	Mesh* kernel, *incompleteKernel = kernelExpansion->getKernel();
	if (incompleteKernel) {
		kernel = computeConvexHull(incompleteKernel->getAllVerts());
		double* kernelPoint = kernelExpansion->getInitialKernelPoint();
		cout << "Used kernel point: " << kernelPoint[0] << " " << kernelPoint[1] << " " << kernelPoint[2] << endl;
		cout << "Cell size is " << cellSizeRatio << " of the grid's diagonal: " << kernelExpansion->getGrid().getCellSize() << endl;
		cout << "Number of cells in grid: " << kernelExpansion->getGrid().getNumOfCells() << endl;
		cout << "Number of found kernel vertices: " << incompleteKernel->getNumOfVerts() << "." << endl;
		cout << "Mesh: [faces: " << mesh->getNumOfTris() << "], [edges: " << mesh->getNumOfEdges() << "], [vertices: " << mesh->getNumOfVerts() << "]" << endl;
		cout << "Kernel: [faces: " << kernel->getNumOfTris() << "], [edges: " << kernel->getNumOfEdges() << "], [vertices: " << kernel->getNumOfVerts() << "]" << endl;
		delete incompleteKernel;
	}
	else {
		kernel = NULL;
		cout << "Kernel is empty!" << endl;
		cout << "Mesh: [faces: " << mesh->getNumOfTris() << "]\n";
	}

	// Print the average time
	double elapsed_secs = totalTime / executionCount;
	cout << "Kernel computation has been completed in " << elapsed_secs << " second(s) by MDF-KER-PLUS." << endl;
	delete kernelExpansion;
	return kernel;

}

Mesh* ComputeKernelByKerTrack(Mesh* mesh) {

	// Execute by <executionCount>-many times
	double totalTime = 0;
	for (int i = 0; i < executionCount; i++) {
		clock_t begin = clock();
		KernelExpansion_KerTrack kertrack(*mesh);
		kertrack.expandKernel();
		Mesh* kernel = NULL, *incompleteKernel = kertrack.getKernel();
		if (incompleteKernel)
			kernel = computeConvexHull(incompleteKernel->getAllVerts());
		clock_t end = clock();
		totalTime += double(end - begin) / CLOCKS_PER_SEC;
		
		if (incompleteKernel) {
			delete kernel;
			delete incompleteKernel;
		}
	}
	
	// Execute to see the results
	KernelExpansion* kernelExpansion = new KernelExpansion_KerTrack(*mesh);
	kernelExpansion->expandKernel();
	Mesh* kernel, *incompleteKernel = kernelExpansion->getKernel();
	if (incompleteKernel) {
		kernel = computeConvexHull(incompleteKernel->getAllVerts());
		cout << "Number of found kernel vertices: " << incompleteKernel->getNumOfVerts() << "." << endl;
		cout << "Mesh: [faces: " << mesh->getNumOfTris() << "], [edges: " << mesh->getNumOfEdges() << "], [vertices: " << mesh->getNumOfVerts() << "]" << endl;
		cout << "Kernel: [faces: " << kernel->getNumOfTris() << "], [edges: " << kernel->getNumOfEdges() << "], [vertices: " << kernel->getNumOfVerts() << "]" << endl;
		delete incompleteKernel;
	}
	else {
		kernel = NULL;
		cout << "Kernel is empty!" << endl;
		cout << "Mesh: [faces: " << mesh->getNumOfTris() << "]\n";
	}		

	// Print the average time
	double elapsed_secs = totalTime / executionCount;
	cout << "Kernel computation has been completed in " << elapsed_secs << " second(s) by KerTrack." << endl;
	delete kernelExpansion;
	return kernel;
}

Mesh* ComputeKernelByCGAL(Mesh* mesh, double* extremeDirection) {

	// Execute by <executionCount>-many times
	double totalTime = 0;
	for (int i = 0; i < executionCount; i++) {
		clock_t begin = clock();
		double* kernelPoint = sdlpMain(mesh, extremeDirection);
		Mesh* kernel = NULL;
		if (kernelPoint)
			kernel = computeKernelByCGAL(mesh, kernelPoint);
		clock_t end = clock();
		totalTime += double(end - begin) / CLOCKS_PER_SEC;

		if (kernelPoint) {
			delete[] kernelPoint;
			delete kernel;
		}
	}

	// Execute to see the results
	double* kernelPoint = sdlpMain(mesh, extremeDirection);
	Mesh* kernel;
	if (kernelPoint) {
		kernel = computeKernelByCGAL(mesh, kernelPoint);
		cout << "Mesh: [faces: " << mesh->getNumOfTris() << "], [edges: " << mesh->getNumOfEdges() << "], [vertices: " << mesh->getNumOfVerts() << "]" << endl;
		cout << "Kernel: [faces: " << kernel->getNumOfTris() << "], [edges: " << kernel->getNumOfEdges() << "], [vertices: " << kernel->getNumOfVerts() << "]" << endl;
		delete[] kernelPoint;
	}
	else {
		kernel = NULL;
		cout << "Kernel is empty!" << endl;
		cout << "Mesh: [faces: " << mesh->getNumOfTris() << "]\n";
	}

	// Print the average time
	double elapsed_secs = totalTime / executionCount;
	cout << "Kernel computation has been completed in " << elapsed_secs << " second(s) by CGAL." << endl;
	return kernel;
}

void CompareFloodFillAndCGAL(string meshName) {
	Mesh* mesh = new Mesh();
	if (meshName.substr(meshName.length() - 3, 3) == "off")
		mesh->loadOff(meshName.c_str());
	else
		mesh->loadObj(meshName.c_str());

	double* extremeDirection = extremeDirections[4];
	double cellSizeRatio = 0.035;// cellSizeRatios[6][1];
	Grid grid(mesh);
	grid.constructGridByDiagonalRatio(cellSizeRatio);
	//KernelExpansion* kernelExpansion = new KernelExpansion_MDFKer(extremeDirection, *mesh, grid);	// ALGO: MINIMUM_DISTANCED_PLANE_TRIPLE
	KernelExpansion* kernelExpansion = new KernelExpansion_SDFKer(*mesh, cellSizeRatio);

	// Compute kernel by flood fill (ff)
	clock_t begin_ff = clock();
	kernelExpansion->expandKernel();
	clock_t end_ff = clock();
	double elapsed_secs_ff = double(end_ff - begin_ff) / CLOCKS_PER_SEC;
	cout << "Flood fill kernel computation has been completed in " << elapsed_secs_ff << " second(s)." << endl;
	Mesh* kernel_ff = kernelExpansion->getKernel();
	//Mesh* incompleteKernel = kernelExpansion->getKernel();												// ALGO: MINIMUM_DISTANCED_PLANE_TRIPLE
	//Mesh* kernel_ff = computeConvexHull(incompleteKernel->getAllVerts());									// ALGO: MINIMUM_DISTANCED_PLANE_TRIPLE


	// Compute kernel by cgal
	clock_t begin_cgal = clock();
	double* kernelPoint = kernelExpansion->getInitialKernelPoint();
	Mesh* kernel_cgal = computeKernelByCGAL(mesh, kernelPoint);
	clock_t end_cgal = clock();
	double elapsed_secs_cgal = double(end_cgal - begin_cgal) / CLOCKS_PER_SEC;
	cout << "CGAL kernel computation has been completed in " << elapsed_secs_cgal << " second(s)." << endl << endl;
	
	//cout << "num of found kernel vertices: " << incompleteKernel->getNumOfVerts() << endl;	// ALGO: MINIMUM_DISTANCED_PLANE_TRIPLE
	//delete incompleteKernel;																	// ALGO: MINIMUM_DISTANCED_PLANE_TRIPLE
	cout << "num of convex hull kernel vertices: " << kernel_ff->getNumOfVerts() << endl;

	// Produce color map
	vector<double> distances = produceColorSource(kernel_cgal, kernel_ff);
	drawDoubleMeshToScene(kernel_cgal, mesh, distances);

	delete kernel_ff;
	delete kernel_cgal;
	delete kernelExpansion;
	delete mesh;
}

void FindKernelPoint_SDLP(string meshName) {

	Mesh* mesh = new Mesh();
	if (meshName.substr(meshName.length() - 3, 3) == "off")
		mesh->loadOff(meshName.c_str());
	else
		mesh->loadObj(meshName.c_str());

	double extremeDirection[3] = { 0, 0, 1 };

	double* kernel_point = sdlpMain(mesh, extremeDirection);
	if (kernel_point != NULL) {
		cout << "Final kernel point: " << kernel_point[0] << " " << kernel_point[1] << " " << kernel_point[2] << endl;
		delete[] kernel_point;
	}
	delete mesh;
}

void SphericalParametrize(string meshName) {

	Mesh* mesh = new Mesh();
	if (meshName.substr(meshName.length() - 3, 3) == "off")
		mesh->loadOff(meshName.c_str());
	else
		mesh->loadObj(meshName.c_str());

	double extremeDirection[3] = { 0, 0, 1 };

	double* kernelPoint = sdlpMain(mesh, extremeDirection);
	double* center = kernelPoint;
	float radius = 1.0;

	// project mesh vertices onto the sphere
	for (int i = 0; i < mesh->getNumOfVerts(); i++) {
		Vertex* v = mesh->getVertex(i);
		double rayDirection[3];
		for (int j = 0; j < 3; j++)
			rayDirection[j] = v->coords[j] - center[j];
		normalize(rayDirection);
		// new coordinates:
		for (int j = 0; j < 3; j++)
			v->coords[j] = center[j] + rayDirection[j] * radius;
	}

	drawSphereOnMesh(mesh, center, radius);

	delete[] center;
	delete mesh;
}

vector<double> produceColorSource(Mesh* ground_truth, Mesh* exp_mesh) {

	vector<double> distances;
	double maxDistance = -numeric_limits<double>::infinity();

	for (int i = 0; i < ground_truth->getNumOfVerts(); i++) {

		Vertex* v1 = ground_truth->getVertex(i);		
		double minDistance = numeric_limits<double>::infinity();

		for (int j = 0; j < exp_mesh->getNumOfVerts(); j++) {
	
			Vertex* v2 = exp_mesh->getVertex(j);
			double* diff = diffVects(v1->coords, v2->coords);
			double distance = computeLength(diff);
			delete[] diff;

			if (distance < minDistance)
				minDistance = distance;
		}

		if (minDistance < EPSILON)
			minDistance = 0;
		distances.push_back(minDistance);
		if (minDistance > maxDistance)
			maxDistance = minDistance;
	}

	if (!(maxDistance < EPSILON)) {
		// normalize
		for (int i = 0; i < distances.size(); i++)
			distances[i] /= maxDistance;
	}

	return distances;
}

