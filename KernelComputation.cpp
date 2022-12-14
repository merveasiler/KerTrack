// @author Merve Asiler

#include "KernelComputation.h"
#include "SphericalParametrization.h"
#include "ShapeMorphing.h"
#include "KernelExpansion_KerTrack.h"
#include "ComputeKernelByCGAL.h"
#include "sdlp.h"
#include "SceneManager.h"
#include "CommonUtils.h"
#include "CGALUtils.h"

#include <boost/filesystem.hpp>
#include <ctime>
#include <fstream>

using namespace boost::filesystem;

/* ***************************** GLOBAL VARIABLES ***************************** */
/* */ vector<double> produceColorSource(Mesh* ground_truth, Mesh* exp_mesh);
/* */ void computeCenterOfKernel(Mesh& mesh, double center[3]);
/* */ string fileNames[9] = {"321.off", "325.off", "326.off", "350.off", "cube.off", "cube_broken.off", "liver.obj", "Rock_2.obj", "Rock_6.obj"};
/* */ double extremeDirections[6][3] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
/* */ double cellSizeRatios[9][4] = {{0.07, 0.05, 0.035, 0.01}, {0.07, 0.05, 0.035, 0.01}, {0.07, 0.05, 0.035, 0.01}, {0.01, 0.005, 0.0045, 0.0025},
/* */								{0.07, 0.05, 0.035, 0.01}, {0.07, 0.05, 0.035, 0.01}, {0.07, 0.035, 0.025, 0.01}, {0.07, 0.035, 0.025, 0.01}, {0.07, 0.035, 0.025, 0.01}};
/* */ double executionCount = 10;
/* **************************************************************************** */

void doExperimentForPaper(string meshName) {

	/************************************************ READ MESH ************************************************/
	Mesh mesh;
	if (meshName.substr(meshName.length() - 3, 3) == "off")
		mesh.loadOff(meshName.c_str());
	else
		mesh.loadObj(meshName.c_str());

	/*************************************** INGREDIENTS / PREPARATIONS  ***************************************/
		// ... for kernel computation
	double extremeDirection[3] = { 0, 0, 1 };
	vector<double> cellSizeRatios{ 0.075}; // the ideal scale is the 2nd one 
	vector<Mesh> kernels;
	kernels.resize(cellSizeRatios.size() * 2 + 1);	// number of experiments: 1 + 4 + 4
		// ... for positions of the shapes on the scene
	vector<double*> positions;
	Grid boundingBox(mesh, NULL, NULL, NULL, 1);	// compute the bounding box of the shape
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
	Mesh groundTruth = kernels[0] = ComputeKernelByCGAL(mesh, extremeDirection);
	positions.push_back(new double[3]{objectWidth[0] + breakWidth[0], 0, 0 });	// on the corner scene
	outputs.push_back(make_tuple(make_tuple(&kernels[0], kernelMatSetting), make_tuple(&mesh, meshMatSetting)));
	cout << endl;
	
		// ... KERTRACK
	for (int j = 1, i = 0; i < kernels.size()/2; i++, j++) {
		kernels[j] = ComputeKernelByKerTrack(mesh);
		double* position = new double[3]{ 0, 0, 0 };	// on the lower row
		position[0] = i * (objectWidth[0] + breakWidth[0]);
		positions.push_back(position);
		outputs.push_back(make_tuple(make_tuple(&kernels[j], kernelMatSetting), make_tuple(&mesh, meshMatSetting)));
		cout << endl;
	}

	if (groundTruth.getNumOfVerts() > 0) {

		// POLYHEDRON-KERNEL (ITALIAN JOB)
		Mesh polyhedronKernel;
		polyhedronKernel.loadOff((meshName.substr(meshName.length() - 3, 3) + "_kernel.off").c_str());
		kernels[2] = polyhedronKernel;
		positions.push_back(new double[3]{ -1, -1, -1 });	// not to be used
		cout << endl;

		/*********************************** COMPARE HAUSDORFF DISTANCES & VOLUMES ***********************************/
		double groundTruthVolume = groundTruth.computeVolume();
		for (int j = 0, mod = (kernels.size() - 1) / 2, i = 1; i < kernels.size(); i++, j++) {
			double volume = 0;
			if (kernels[i].getNumOfVerts() > 0)
				volume = kernels[i].computeVolume();
			if (volume < EPSILON)
				continue;
			cout << "Volume of the kernel of " << cellSizeRatios[j % mod] << ": " << volume << " out of " << groundTruthVolume << "." << endl;
			double volumeDiff = groundTruthVolume - volume;
			if (volumeDiff < 1e9*EPSILON)
				volumeDiff = 0;
			cout << "Volume difference between the one <" << cellSizeRatios[j % mod] << "> and ground truth (CGAL): " << volumeDiff << endl;

			double* hd = computeHausdorffDistance(groundTruth, kernels[i]);
			for (int k=0; k<3; k++)
				if (hd[k] < 1e10*EPSILON)
					hd[k] = 0;
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
	for (int i = 0; i < kernels.size(); i++) 
		delete[] positions[i];

	delete kernelMatSetting;
	delete meshMatSetting;

}

void ComputeKernel(string meshName, string algoType) {

	// Read mesh
	Mesh mesh;
	if (meshName.substr(meshName.length() - 3, 3) == "off")
		mesh.loadOff(meshName.c_str());
	else
		mesh.loadObj(meshName.c_str());

	// Ingredients:
	double extremeDirection[3] = { 0, -1, 0 };
	double cellSizeRatio = 0.07;
	int gridDimension[3] = {2, 2, 2};
	Mesh kernel;

	// Compute kernel
	
	/*
	if (algoType == "sdfker")
		kernel = ComputeKernelBySDFKer(mesh, cellSizeRatio);
	else if (algoType == "sdfker_plus")
		kernel = ComputeKernelBySDFKerPlus(mesh, cellSizeRatio, extremeDirection);
	else if (algoType == "mdfker")
		kernel = ComputeKernelByMDFKer(mesh, gridDimension);
	else if (algoType == "mdfker_plus")
		kernel = ComputeKernelByMDFKerPlus(mesh, cellSizeRatio);
	*/
	if (algoType == "kertrack")
		kernel = ComputeKernelByKerTrack(mesh);
	else if (algoType == "kernel_by_cgal")
		kernel = ComputeKernelByCGAL(mesh, extremeDirection);
	else;

	// Draw
	if (kernel.getNumOfVerts() > 0) {
		MaterialSetting* kernelMatSetting = new MaterialSetting(0, 0, 1, 0);
		MaterialSetting* meshMatSetting = new MaterialSetting(1, 1, 1, 0.5);
		vector<tuple<Mesh*, MaterialSetting*>> mesh_mat_set = { make_tuple(&kernel, kernelMatSetting), make_tuple(&mesh, meshMatSetting) };
		drawMultipleMeshToScene(mesh_mat_set);
	}

}

void ComputeBatchKernel(string inputFolderName, string outputFolderName, string algoType) {

	// Read folder, fecth mesh names
	vector<string> meshNames;
	path p(inputFolderName);
	for (auto i = directory_iterator(p); i != directory_iterator(); i++)
	{
		if (!is_directory(i->path())) //we eliminate directories
			meshNames.push_back(i->path().filename().string());
	}

	double avgTime_star = 0, avgTime_nonstar = 0;
	int numOfStarShapes = 0, numOfNonStarShapes = 0;
	string extension;
	std::ofstream outputFile;
	if (algoType == "batch_kernel_kertrack") {
		outputFile.open(outputFolderName + "/KernelResults_KerTrack.txt");
		extension = "_kernel_by_kertrack.off";
	}
	else if (algoType == "batch_kernel_cgal") {
		outputFile.open(outputFolderName + "/KernelResults_CGAL.txt");
		extension = "_kernel_by_cgal.off";
	}

	// Compute batch kernel
	for (int i = 0; i < meshNames.size(); i++) {
		// Read mesh
		string meshName = meshNames[i];
		Mesh mesh;
		if (meshName.substr(meshName.length() - 3, 3) == "off")
			mesh.loadOff((inputFolderName + "/" + meshName).c_str());
		else
			mesh.loadObj((inputFolderName + "/" + meshName).c_str());

		if (mesh.isManifold() == false) {
			cout << "NOT MANIFOLD: " << i << ": " << meshName << " !!!!!!!!!!!!!!!!!!" << endl;
			continue;
		}
		cout << i << ": " << meshName << endl;

		// kernel computation
		Mesh kernel;
		clock_t	begin = clock();
		if (algoType == "batch_kernel_kertrack") {
			KernelExpansion_KerTrack kertrack(mesh);
			kertrack.expandKernel();
			Mesh & incompleteKernel = kertrack.getKernel();
			if (incompleteKernel.getNumOfVerts() > 0)
				kernel = computeConvexHull(incompleteKernel.getAllVerts());
		}
		else if (algoType == "batch_kernel_cgal") {
			kernel = computeKernelByCGAL(mesh, NULL);
		}
		else;
		clock_t	end = clock();

		outputFile << meshName << endl;
		double totalTime = double(end - begin) / CLOCKS_PER_SEC;

		// output notes
		if (kernel.getNumOfVerts() > 0) {
			outputFile << "\t" << "Number of found kernel vertices: " << kernel.getNumOfVerts() << "." << endl;
			outputFile << "\t" << "Mesh: [faces: " << mesh.getNumOfTris() << "], [edges: " << mesh.getNumOfEdges() << "], [vertices: " << mesh.getNumOfVerts() << "]" << endl;
			outputFile << "\t" << "Kernel: [faces: " << kernel.getNumOfTris() << "], [edges: " << kernel.getNumOfEdges() << "], [vertices: " << kernel.getNumOfVerts() << "]" << endl;
			//kernel->writeOff(outputFolderName + "/" + meshName.substr(0, meshName.length() - 4) + extension);
			avgTime_star += totalTime;
			numOfStarShapes++;
		}
		else {
			outputFile << "\t" << "Kernel is empty!" << endl;
			outputFile << "\t" << "Mesh: [faces: " << mesh.getNumOfTris() << "]\n";
			avgTime_nonstar += totalTime;
			numOfNonStarShapes++;
		}

		outputFile << "\t" << "Kernel computation has been completed in " << totalTime << " second(s) by KerTrack." << endl;
	}

	outputFile << endl << endl;
	outputFile << "AVERAGE Kernel computation has been completed in " << avgTime_star / numOfStarShapes << " second(s) for " << numOfStarShapes << " star-shapes." << endl;
	outputFile << "AVERAGE Kernel computation has been completed in " << avgTime_nonstar / numOfNonStarShapes << " second(s) for " << numOfNonStarShapes << " nonstar-shapes." << endl;
	outputFile.close();

}

Mesh ComputeKernelByKerTrack(Mesh& mesh) {

	// Execute by <executionCount>-many times
	double totalTime = 0;
	for (int i = 0; i < executionCount; i++) {
		KernelExpansion_KerTrack kertrack(mesh);
		clock_t begin = clock();
		kertrack.expandKernel(); 
		clock_t end = clock(); 
		Mesh* kernel = NULL, &incompleteKernel = kertrack.getKernel();
		totalTime += double(end - begin) / CLOCKS_PER_SEC;
	}
	
	// Execute to see the results
	KernelExpansion_KerTrack kertrack(mesh);
	kertrack.expandKernel();
	Mesh kernel, &incompleteKernel = kertrack.getKernel();
	if (incompleteKernel.getNumOfVerts() > 0) {
		kernel = computeConvexHull(incompleteKernel.getAllVerts());
		cout << "Number of found kernel vertices: " << incompleteKernel.getNumOfVerts() << "." << endl;
		cout << "Mesh: [faces: " << mesh.getNumOfTris() << "], [edges: " << mesh.getNumOfEdges() << "], [vertices: " << mesh.getNumOfVerts() << "]" << endl;
		cout << "Kernel: [faces: " << kernel.getNumOfTris() << "], [edges: " << kernel.getNumOfEdges() << "], [vertices: " << kernel.getNumOfVerts() << "]" << endl;
	}
	else {
		cout << "Kernel is empty!" << endl;
		cout << "Mesh: [faces: " << mesh.getNumOfTris() << "]\n";
	}		

	// Print the average time
	double elapsed_secs = totalTime / executionCount;
	cout << "Kernel computation has been completed in " << elapsed_secs << " second(s) by KerTrack." << endl;

	return kernel;
}

Mesh ComputeKernelByCGAL(Mesh& mesh, double* extremeDirection) {

	// Execute by <executionCount>-many times
	double totalTime = 0;
	for (int i = 0; i < executionCount; i++) {
		clock_t begin = clock();
		Mesh kernel = computeKernelByCGAL(mesh, NULL);
		clock_t end = clock();
		totalTime += double(end - begin) / CLOCKS_PER_SEC;
	}

	// Execute to see the results
	Mesh kernel;
	kernel = computeKernelByCGAL(mesh, NULL);
	if (kernel.getNumOfVerts() > 0) {
		cout << "Mesh: [faces: " << mesh.getNumOfTris() << "], [edges: " << mesh.getNumOfEdges() << "], [vertices: " << mesh.getNumOfVerts() << "]" << endl;
		cout << "Kernel: [faces: " << kernel.getNumOfTris() << "], [edges: " << kernel.getNumOfEdges() << "], [vertices: " << kernel.getNumOfVerts() << "]" << endl;
	}
	else {
		cout << "Kernel is empty!" << endl;
		cout << "Mesh: [faces: " << mesh.getNumOfTris() << "]\n";
	}

	// Print the average time
	double elapsed_secs = totalTime / executionCount;
	cout << "Kernel computation has been completed in " << elapsed_secs << " second(s) by CGAL." << endl;
	return kernel;
}

void FindKernelPoint_SDLP(string meshName) {

	Mesh* mesh = new Mesh();
	if (meshName.substr(meshName.length() - 3, 3) == "off")
		mesh->loadOff(meshName.c_str());
	else
		mesh->loadObj(meshName.c_str());

	double extremeDirection[3] = { 0, 0, 1 };

	double* kernel_point = sdlpMain(*mesh, extremeDirection);
	if (kernel_point != NULL) {
		cout << "Final kernel point: " << kernel_point[0] << " " << kernel_point[1] << " " << kernel_point[2] << endl;
		delete[] kernel_point;
	}
	delete mesh;
}

void SphericalParametrize(string meshName) {

	Mesh mesh;
	if (meshName.substr(meshName.length() - 3, 3) == "off")
		mesh.loadOff(meshName.c_str());
	else
		mesh.loadObj(meshName.c_str());

	int resolution = 10;
	double radius[1] = { 1.0 };
	
	double center[3];
	computeCenterOfKernel(mesh, center);

	/*
	double extremeDirection[3] = { 0, 0, 1 };
	double* kernelPoint = sdlpMain(mesh, extremeDirection);
	double* center = kernelPoint;
	*/

	Mesh sphericalMesh;
	parametrizeByKernel(mesh, sphericalMesh, center, radius, resolution);
	drawMeshOnSphere(&sphericalMesh, &mesh, center, radius[0]);

}

void ShapeMorphByKernel(string sourceMeshName, string targetMeshName) {

	Mesh* sourceMesh = new Mesh();
	Mesh* targetMesh = new Mesh();

	if (sourceMeshName.substr(sourceMeshName.length() - 3, 3) == "off")
		sourceMesh->loadOff(sourceMeshName.c_str());
	else
		sourceMesh->loadObj(sourceMeshName.c_str());
	
	if (targetMeshName.substr(targetMeshName.length() - 3, 3) == "off")
		targetMesh->loadOff(targetMeshName.c_str());
	else
		targetMesh->loadObj(targetMeshName.c_str());

	int numOfInterMeshes = 5;
	vector<Mesh*> interMeshes;
	for (int i = 0; i < numOfInterMeshes; i++)
		interMeshes.push_back(new Mesh());

	double centerSource[3] = { 0, 0, 0 };
	double centerTarget[3] = { 0, 0, 0 };
	//double centerSource[3], centerTarget[3];
	//computeCenterOfKernel(sourceMesh, centerSource);
	//computeCenterOfKernel(targetMesh, centerTarget);
	
	morphByKernel(sourceMesh, targetMesh, interMeshes, centerSource, centerTarget);

	vector<tuple<Mesh*, MaterialSetting*>> outputs;
	MaterialSetting* meshMatSetting = new MaterialSetting(1, 1, 1, 0);
	outputs.push_back(make_tuple(sourceMesh, meshMatSetting));
	for (int i = 0; i < interMeshes.size(); i++)
		outputs.push_back(make_tuple(interMeshes[i], meshMatSetting));

	drawMultipleMeshToScene(outputs);

	delete sourceMesh;
	delete targetMesh;	
}

void ShapeMorphByLerp(string sourceMeshName, string targetMeshName) {

	Mesh* sourceMesh = new Mesh();
	Mesh* targetMesh = new Mesh();

	if (sourceMeshName.substr(sourceMeshName.length() - 3, 3) == "off")
		sourceMesh->loadOff(sourceMeshName.c_str());
	else
		sourceMesh->loadObj(sourceMeshName.c_str());

	if (targetMeshName.substr(targetMeshName.length() - 3, 3) == "off")
		targetMesh->loadOff(targetMeshName.c_str());
	else
		targetMesh->loadObj(targetMeshName.c_str());

	int numOfInterMeshes = 5;
	vector<Mesh*> interMeshes;
	for (int i = 0; i < numOfInterMeshes; i++)
		interMeshes.push_back(new Mesh());

	morphByLerp(sourceMesh, targetMesh, interMeshes);

	vector<tuple<Mesh*, MaterialSetting*>> outputs;
	MaterialSetting* meshMatSetting = new MaterialSetting(1, 1, 1, 0);
	outputs.push_back(make_tuple(sourceMesh, meshMatSetting));
	for (int i = 0; i < interMeshes.size(); i++)
		outputs.push_back(make_tuple(interMeshes[i], meshMatSetting));

	drawMultipleMeshToScene(outputs);

	delete sourceMesh;
	delete targetMesh;
}

void computeCenterOfKernel(Mesh& mesh, double center[3]) {

	KernelExpansion_KerTrack kt(mesh);
	kt.expandKernel();
	Mesh& kernel = kt.getKernel();

	for (int k = 0; k < 3; k++)
		center[k] = 0;
	for (int i = 0; i < kernel.getNumOfVerts(); i++) {
		for (int k = 0; k < 3; k++)
			center[k] += kernel.getVertex(i).coords[k];
	}
	for (int k = 0; k < 3; k++)
		center[k] /= kernel.getNumOfVerts();

}

vector<double> produceColorSource(Mesh* ground_truth, Mesh* exp_mesh) {

	vector<double> distances;
	double maxDistance = -numeric_limits<double>::infinity();

	for (int i = 0; i < ground_truth->getNumOfVerts(); i++) {

		Vertex v1 = ground_truth->getVertex(i);		
		double minDistance = numeric_limits<double>::infinity();

		for (int j = 0; j < exp_mesh->getNumOfVerts(); j++) {
	
			Vertex v2 = exp_mesh->getVertex(j);
			double* diff = diffVects(v1.coords, v2.coords);
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

