// @author Merve Asiler

#include "KernelExpansion_MDFKer.h"
#include "BaseGeoOpUtils.h"

KernelExpansion_MDFKer::KernelExpansion_MDFKer(const Mesh& hostMesh, int gridDimension[3]) :
	KernelExpansion(hostMesh, gridDimension) {};

KernelExpansion_MDFKer::~KernelExpansion_MDFKer() {
	for (int i = 0; i < allTriedPlaneCombinations.size(); i++)
		delete[] allTriedPlaneCombinations[i];
	allTriedPlaneCombinations.clear();
}

void KernelExpansion_MDFKer::expandKernel() {

	// is kernel empty?
	if (this->extremeCorners[0] == nullptr && this->extremeCorners[1] == nullptr) {
		delete kernel;
		kernel = NULL;
		return;	// NOT STAR-SHAPED!
	}

	bool** backFace = new bool* [grid.numOfCells[1]], ** frontFace = new bool* [grid.numOfCells[1]];
	for (int j = 0; j < grid.numOfCells[1]; j++) {
		backFace[j] = new bool[grid.numOfCells[2]];
		frontFace[j] = new bool[grid.numOfCells[2]];
		for (int k = 0; k < grid.numOfCells[2]; k++) {
			backFace[j][k] = true;
			frontFace[j][k] = true;
		}
	}

	for (int i = 0; i < grid.numOfCells[0]; i++) {
		for (int j = 0; j < grid.numOfCells[1]; j++) {

			bool loop_broken = false;
			vector<int> preExcluders;

			for (int k = 0; k < grid.numOfCells[2]; k++) {
				vector<int> excluders = evaluateCorner(i, j, k);
				if (excluders.size() == 0) {
					frontFace[j][k] = true;
					//if (k!=0)
					//	evaluateCorner(i, j, k-0.5);

					if (i == 0 || j == 0 || k == 0);
					else {
						double prevIndices[3] = { i - 1.0, j - 1.0, k - 1.0 };
						double nextIndices[3] = { i, j, k };
						bool prevVals[3] = { backFace[j][k], frontFace[j - 1][k], frontFace[j][k - 1] };
						expandByZoom(prevIndices, nextIndices, prevVals, 1);
					}
					loop_broken = true;
					break;
				}
				else {
					frontFace[j][k] = false;
					if (excluders.size() == 1 && preExcluders.size() == 1)
						if (dotProduct(halfSpaceSet[excluders[0]]->ABCD, halfSpaceSet[preExcluders[0]]->ABCD) <= 0) {
							double prevIndices[3] = { i - 1.0, j - 1.0, k - 1.0 };
							double nextIndices[3] = { i, j, k };
							bool prevVals[3] = { backFace[j][k], frontFace[j - 1][k], frontFace[j][k - 1] };
							expandByZoom(prevIndices, nextIndices, prevVals, 1);
						}
				}
				preExcluders.clear();
				preExcluders = excluders;
			}

			if (loop_broken) {
				vector<int> preExcluders;
				for (int k = grid.numOfCells[2] - 1; k >= 0; k--) {
					vector<int> excluders = evaluateCorner(i, j, k);
					if (excluders.size() == 0) {
						frontFace[j][k] = true;
						//if (k != grid.numOfCells[2] - 1)
						//	evaluateCorner(i, j, k+0.5);

						if (i == 0 || j == 0 || k == grid.numOfCells[2] - 1);
						else {
							double prevIndices[3] = { i - 1.0, j - 1.0, k + 1.0 };
							double nextIndices[3] = { i, j, k };
							bool prevVals[3] = { backFace[j][k], frontFace[j - 1][k], frontFace[j][k + 1] };
							expandByZoom(prevIndices, nextIndices, prevVals, 1);
						}
						break;
					}
					else {
						frontFace[j][k] = false;
						if (excluders.size() == 1 && preExcluders.size() == 1)
							if (dotProduct(halfSpaceSet[excluders[0]]->ABCD, halfSpaceSet[preExcluders[0]]->ABCD) <= 0) {
								double prevIndices[3] = { i - 1.0, j - 1.0, k - 1.0 };
								double nextIndices[3] = { i, j, k };
								bool prevVals[3] = { backFace[j][k], frontFace[j - 1][k], frontFace[j][k - 1] };
								expandByZoom(prevIndices, nextIndices, prevVals, 1);
							}
					}
					preExcluders.clear();
					preExcluders = excluders;
				}
			}
		}

		for (int j = 0; j < grid.numOfCells[1]; j++) {
			for (int k = 0; k < grid.numOfCells[2]; k++) {
				backFace[j][k] = frontFace[j][k];
				frontFace[j][k] = true;
			}
		}

	}

	for (int j = 0; j < grid.numOfCells[1]; j++) {
		delete[] backFace[j];
		delete[] frontFace[j];
	}
	delete[] backFace;
	delete[] frontFace;

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

}

void KernelExpansion_MDFKer::expandByZoom(double prevIndices[3], double nextIndices[3], bool prevVals[3], int depth) {

	double new_prevIndices[3], new_nextIndices[3];
	for (int d = 0; d < 3; d++) {
		new_prevIndices[d] = prevIndices[d];
		new_nextIndices[d] = nextIndices[d];
	}
	int prev_num_of_kernel_verts = kernel->getNumOfVerts();

	for (int d = 0; d < 3; d++) {
		if (prevVals[d])
			new_prevIndices[d] = nextIndices[d]; // there won't be zooming at that axis
		else {
			double midIndices[3] = { nextIndices[0], nextIndices[1], nextIndices[2] };
			midIndices[d] = (prevIndices[d] + nextIndices[d]) / 2.0;

			vector<int> excluders = evaluateCorner(midIndices[0], midIndices[1], midIndices[2]);
			if (excluders.size() == 0)
				new_nextIndices[d] = midIndices[d];
			else
				new_prevIndices[d] = midIndices[d];
		}
	}

	//if (kernel->getNumOfVerts() > prev_num_of_kernel_verts)
	if (depth < 3)
		expandByZoom(new_prevIndices, new_nextIndices, prevVals, depth + 1);


}

vector<int> KernelExpansion_MDFKer::evaluateCorner(double i, double j, double k) {

	// compute the corner coordinates
	double cell_indices[3] = { i, j, k };
	double cornerCoords[3];
	for (int c = 0; c < 3; c++)
		cornerCoords[c] = grid.minGridCoords[c] + (cell_indices[c] * grid.cellSize[c]);

	// for the corner compute the scalar value vector to be used by putting it into the system of inequalities
	double* scalarsVector = findValueVectorSatisfiedByPoint(cornerCoords, halfSpaceSet);
	vector<int> excluders = kernelInOutState(scalarsVector);

	if (excluders.size() == 0)
		findPossibleKernelVertex(scalarsVector, cornerCoords);

	delete[] scalarsVector;
	return excluders;
}

vector<int> KernelExpansion_MDFKer::kernelInOutState(double* scalarVector) {

	vector<int> excluders;

	for (int i = 0; i < halfSpaceSet.size(); i++) {
		if (scalarVector[i] > 0) {
			for (int j = 0; j < excluders.size(); j++) {
				double cosAngle = dotProduct(halfSpaceSet[excluders[j]]->ABCD, halfSpaceSet[i]->ABCD);
				if (cosAngle < 0) {
					excluders.push_back(i);
					return excluders;
				}
			}
			if (excluders.size() == 0)
				excluders.push_back(i);
		}
	}

	return excluders;

}

void KernelExpansion_MDFKer::findPossibleKernelVertex(double* scalarsVector, double source_coords[3]) {

	int closestPlaneIds[3] = { -1, -1, -1 };
	updateClosest3Planes(scalarsVector, closestPlaneIds, closestPlaneIds);							// find the planes satisfying the minimum 3 distance value by selecting from scalarsVector
	double* candidateKernelVertex = computePossibleKernelVertex(scalarsVector, closestPlaneIds);	// however those 3 planes must be able to create a vertex

	if (candidateKernelVertex) {
		int obstaclePlaneId = addIfNoFrontierObstacle(candidateKernelVertex, closestPlaneIds, source_coords);
		if (obstaclePlaneId < 0)	// no obstacle
			return;
		findByFrontierPlanes(closestPlaneIds, obstaclePlaneId, source_coords);
	}
	else {
		cout << closestPlaneIds[0] << " " << closestPlaneIds[1] << " " << closestPlaneIds[2] << endl;
		cout << "EXPECTED TO BE IMPOSSIBLE!\n";
	}

}

/*
	It takes an array representing some distance values to a determined point and finds the minimum 3 of those values.
	Optionally, 1 or 2 of the minimum values could be given and only the remaining(s) may be asked to be found.
	CurrentPlaneIds are filled with the already given minimum(s), and for the unknown one(s) it includes -1.
	If it includes -1, the previous value of the corresponding minimum is lost (if it existed before). Therefore;
	PreviousPlaneIds is required to hold the last plane indices selected as the minimum distanced ones.
*/
void KernelExpansion_MDFKer::updateClosest3Planes(double* values, int currentPlaneIds[3], int previousPlaneIds[3]) {

	double previousDistances[3];
	double lowerLimit = -1.0;	// this will be an indicator which value has been used at most as the minimum, previously.
	int maxId = -1;				// this will be an indicator to understand which plane was used among the ones equal distanced.
	for (int i = 0; i < 3; i++) {
		if (previousPlaneIds[i] < 0)
			previousDistances[i] = -1.0;
		else
			previousDistances[i] = abs(values[previousPlaneIds[i]]);

		if (previousDistances[i] >= lowerLimit) {
			lowerLimit = previousDistances[i];
			if (currentPlaneIds[i] == -1 && previousPlaneIds[i] > maxId)
				maxId = previousPlaneIds[i];
		}
	}

	double minDistances[3];
	for (int i = 0; i < 3; i++) {
		if (currentPlaneIds[i] < 0)
			minDistances[i] = numeric_limits<double>::infinity();
		else
			minDistances[i] = abs(values[currentPlaneIds[i]]);
	}

	for (int i = 0; i < halfSpaceSet.size(); i++) {
		if (i == previousPlaneIds[0] || i == previousPlaneIds[1] || i == previousPlaneIds[2])
			continue;
		if (abs(values[i]) < lowerLimit)
			continue;
		if (abs(values[i]) == lowerLimit && i <= maxId)
			continue;

		if (abs(values[i]) < minDistances[0] && abs(values[i]) >= lowerLimit) {
			minDistances[2] = minDistances[1];
			currentPlaneIds[2] = currentPlaneIds[1];
			minDistances[1] = minDistances[0];
			currentPlaneIds[1] = currentPlaneIds[0];
			minDistances[0] = abs(values[i]);
			currentPlaneIds[0] = i;
		}
		else if (abs(values[i]) < minDistances[1] && abs(values[i]) >= lowerLimit) {
			minDistances[2] = minDistances[1];
			currentPlaneIds[2] = currentPlaneIds[1];
			minDistances[1] = abs(values[i]);
			currentPlaneIds[1] = i;
		}
		else if (abs(values[i]) < minDistances[2] && abs(values[i]) >= lowerLimit) {
			minDistances[2] = abs(values[i]);
			currentPlaneIds[2] = i;
		}
	}
}

double* KernelExpansion_MDFKer::computePossibleKernelVertex(double* scalarsVector, int closestPlaneIds[3]) {

	Plane* planes[3];
	int previousPlaneIds[3];

	while (true) {
		for (int i = 0; i < 3; i++) {
			previousPlaneIds[i] = closestPlaneIds[i];	// save the current ones as the previous since they may change in the following steps.
			planes[i] = halfSpaceSet[closestPlaneIds[i]];
		}

		// First check any 2 of the 3 planes are parallel (or the same) or not
		if (abs(abs(dotProduct(planes[0]->ABCD, planes[2]->ABCD)) - 1) < 3 * EPSILON) {
			closestPlaneIds[2] = -1;
		}
		else if (abs(abs(dotProduct(planes[0]->ABCD, planes[1]->ABCD)) - 1) < 3 * EPSILON) {
			closestPlaneIds[1] = closestPlaneIds[2];
			closestPlaneIds[2] = -1;
		}
		else if (abs(abs(dotProduct(planes[1]->ABCD, planes[2]->ABCD)) - 1) < 3 * EPSILON) {
			//closestPlaneIds[1] = closestPlaneIds[2];
			closestPlaneIds[2] = -1;
		}
		else {
			// Second, check if the 3 planes coincide on a line, or a point
			double* point = find3PlaneIntersection(planes);
			if (point)
				return point;
			else
				closestPlaneIds[2] = -1;
		}

		updateClosest3Planes(scalarsVector, closestPlaneIds, previousPlaneIds);
	}

}

double* KernelExpansion_MDFKer::computePossibleKernelVertex(int closestPlaneIds[3]) {

	Plane* closestPlanes[3] = { halfSpaceSet[closestPlaneIds[0]], halfSpaceSet[closestPlaneIds[1]], halfSpaceSet[closestPlaneIds[2]] };
	double* candidateKernelVertex = find3PlaneIntersection(closestPlanes);
	return candidateKernelVertex;
}

int KernelExpansion_MDFKer::addIfNoFrontierObstacle(double candidateKernelVertex[3], int candidateVertexProducerIds[3], double source_coords[3]) {

	int frontierObstacle = -1;
	if (isAlreadyKernelVertex(candidateKernelVertex) == false) {
		int obstaclePlaneId = findIfFrontierObstacle(candidateKernelVertex, candidateVertexProducerIds, source_coords);
		if (obstaclePlaneId < 0)
			kernel->addVertex(candidateKernelVertex[0], candidateKernelVertex[1], candidateKernelVertex[2]);
		else
			frontierObstacle = obstaclePlaneId;
	}

	delete[] candidateKernelVertex;
	return frontierObstacle;
}

bool KernelExpansion_MDFKer::isAlreadyKernelVertex(double coords[3]) {

	for (int i = 0; i < kernel->getNumOfVerts(); i++) {
		Vertex* v = kernel->getVertex(i);
		if (abs(v->coords[0] - coords[0]) < EPSILON && abs(v->coords[1] - coords[1]) < EPSILON && abs(v->coords[2] - coords[2]) < EPSILON)
			return true;
	}
	return false;
}

int KernelExpansion_MDFKer::findIfFrontierObstacle(double* candidateKernelVertex, int candidateVertexProducerIds[3], double* sourcePoint) {

	int frontierHalfSpaceId = -1;

	double* direction = new double[3];
	for (int i = 0; i < 3; i++)
		direction[i] = candidateKernelVertex[i] - sourcePoint[i];
	double candidateDistance = computeLength(direction);
	normalize(direction);

	Line* line = new Line(sourcePoint, direction);
	delete[] direction;

	double* scalarVector = findValueVectorSatisfiedByPoint(candidateKernelVertex, halfSpaceSet);

	for (int i = 0; i < halfSpaceSet.size(); i++) {
		if (scalarVector[i] > 0) {
			double distance = findLinePlaneIntersection(line, halfSpaceSet[i]);
			if (distance >= 0 && distance < candidateDistance) {
				frontierHalfSpaceId = i;
				candidateDistance = distance;
			}
		}
	}

	delete line;
	delete[] scalarVector;
	return frontierHalfSpaceId;

}

void KernelExpansion_MDFKer::findByFrontierPlanes(int closestPlaneIds[3], int obstaclePlaneId, double source_coords[3]) {

	queue<int*> possiblePlaneCombinations;
	replaceOnePlaneBy(possiblePlaneCombinations, closestPlaneIds, obstaclePlaneId);

	while (!possiblePlaneCombinations.empty()) {
		int* planeIds = possiblePlaneCombinations.front();
		double* candidateKernelVertex = computePossibleKernelVertex(planeIds);
		if (candidateKernelVertex) {
			int obstacle_id = addIfNoFrontierObstacle(candidateKernelVertex, planeIds, source_coords);
			if (obstacle_id >= 0)
				replaceOnePlaneBy(possiblePlaneCombinations, planeIds, obstacle_id);
		}
		possiblePlaneCombinations.pop();
	}

	for (int i = 0; i < allTriedPlaneCombinations.size(); i++)
		delete[] allTriedPlaneCombinations[i];
	allTriedPlaneCombinations.clear();
}

void KernelExpansion_MDFKer::replaceOnePlaneBy(queue<int*>& possiblePlaneCombinations, int planeIds[3], int replacerId) {

	int* candidates[3];
	candidates[0] = new int[3]{ planeIds[0], planeIds[1], replacerId };
	candidates[1] = new int[3]{ planeIds[0], planeIds[2], replacerId };
	candidates[2] = new int[3]{ planeIds[1], planeIds[2], replacerId };

	for (int c = 0; c < 3; c++) {
		int* cand = candidates[c];
		bool is_found = false;
		for (int i = 0; i < allTriedPlaneCombinations.size(); i++) {
			int* comb = allTriedPlaneCombinations[i];
			for (int j = 0; j < 3; j++) {
				if (comb[j] == cand[0] && comb[(j + 1) % 3] == cand[1] && comb[(j + 2) % 3] == cand[2]) {
					is_found = true;
					break;
				}
				if (comb[j] == cand[0] && comb[(j + 1) % 3] == cand[2] && comb[(j + 2) % 3] == cand[1]) {
					is_found = true;
					break;
				}
			}

			if (is_found)
				break;
		}

		if (!is_found) {
			allTriedPlaneCombinations.push_back(cand);
			possiblePlaneCombinations.push(cand);
		}
	}
}

void KernelExpansion_MDFKer::decideOnNeighborCells(int* cell_indices, double* corner_in_outs, double** scalarsVector) {

	// define if the neighbor cells have different in/out feature for different half-spaces
	NeighborPolarity* edgePolarityTypes = determineNeighborPolarity(scalarsVector);

	// Decide on neighbors for appending into queue
	int i = cell_indices[0], j = cell_indices[1], k = cell_indices[2];
	addNeighborCellIntoQueue(i - 1, j, k, corner_in_outs[0], corner_in_outs[3], corner_in_outs[4], corner_in_outs[7], edgePolarityTypes[1], edgePolarityTypes[2], edgePolarityTypes[7], edgePolarityTypes[9]);	// Neighbor at x-1
	addNeighborCellIntoQueue(i + 1, j, k, corner_in_outs[1], corner_in_outs[2], corner_in_outs[5], corner_in_outs[6], edgePolarityTypes[3], edgePolarityTypes[4], edgePolarityTypes[6], edgePolarityTypes[10]);	// Neighbor at x+1
	addNeighborCellIntoQueue(i, j - 1, k, corner_in_outs[0], corner_in_outs[1], corner_in_outs[2], corner_in_outs[3], edgePolarityTypes[0], edgePolarityTypes[1], edgePolarityTypes[3], edgePolarityTypes[5]);	// Neighbor at y-1
	addNeighborCellIntoQueue(i, j + 1, k, corner_in_outs[4], corner_in_outs[5], corner_in_outs[6], corner_in_outs[7], edgePolarityTypes[8], edgePolarityTypes[9], edgePolarityTypes[10], edgePolarityTypes[11]);// Neighbor at y+1
	addNeighborCellIntoQueue(i, j, k - 1, corner_in_outs[0], corner_in_outs[1], corner_in_outs[4], corner_in_outs[5], edgePolarityTypes[0], edgePolarityTypes[2], edgePolarityTypes[4], edgePolarityTypes[8]);	// Neighbor at z-1
	addNeighborCellIntoQueue(i, j, k + 1, corner_in_outs[2], corner_in_outs[3], corner_in_outs[6], corner_in_outs[7], edgePolarityTypes[5], edgePolarityTypes[6], edgePolarityTypes[7], edgePolarityTypes[11]);	// Neighbor at z+1

}

void KernelExpansion_MDFKer::addNeighborCellIntoQueue(int neigh_i, int neigh_j, int neigh_k,
	double neighVert1_scalar, double neighVert2_scalar, double neighVert3_scalar, double neighVert4_scalar,
	NeighborPolarity np1, NeighborPolarity np2, NeighborPolarity np3, NeighborPolarity np4) {

	if (isCellValidForQueue(neigh_i, neigh_j, neigh_k)) {

		bool add_neighbor_into_queue;

		if (neighVert1_scalar <= 0 || neighVert2_scalar <= 0 || neighVert3_scalar <= 0 || neighVert4_scalar <= 0)
			add_neighbor_into_queue = true;/*
		else if (np1 == NeighborPolarity::DISTINCT_SIDED_INCLUSION || np2 == NeighborPolarity::DISTINCT_SIDED_INCLUSION ||
				 np3 == NeighborPolarity::DISTINCT_SIDED_INCLUSION || np4 == NeighborPolarity::DISTINCT_SIDED_INCLUSION)
			add_neighbor_into_queue = true;*/
		else
			add_neighbor_into_queue = false;

		// it satisfies the above conditions, add it into the queue
		if (add_neighbor_into_queue) {
			int* neighbor_ijk = new int[3]{ neigh_i, neigh_j, neigh_k };
			possibleKernelCells.push(neighbor_ijk);
			handledCells.push_back(neighbor_ijk);
		}
	}
}



void KernelExpansion_MDFKer::findMissingKernelVertices(Mesh* incompleteKernel, double cellSize) {

	for (int t = 0; t < incompleteKernel->getNumOfTris(); t++) {
		Triangle* triangle = incompleteKernel->getTriangle(t);
		triangle->computeNormal(incompleteKernel->getVertex(triangle->corners[0]), incompleteKernel->getVertex(triangle->corners[1]), incompleteKernel->getVertex(triangle->corners[2]));
		bool is_matched = false;
		for (int hs = 0; hs < halfSpaceSet.size(); hs++) {
			if (abs(halfSpaceSet[hs]->ABCD[0] - triangle->normal[0]) < EPSILON && abs(halfSpaceSet[hs]->ABCD[1] - triangle->normal[1]) < EPSILON && abs(halfSpaceSet[hs]->ABCD[2] - triangle->normal[2]) < EPSILON) {
				is_matched = true;
				break;
			}
		}
		if (is_matched)
			continue;

		cout << "THERE HAS BEEN MISSED ONE!" << endl;

		double center[3];
		for (int i = 0; i < 3; i++) {
			center[i] = 0;
			for (int j = 0; j < 3; j++)
				center[i] += incompleteKernel->getVertex(triangle->corners[j])->coords[i];
			center[i] /= 3.0;
			center[i] += triangle->normal[i] * cellSize;
		}

		double* center_ScalarsVector = findValueVectorSatisfiedByPoint(center, halfSpaceSet);
		findPossibleKernelVertex(center_ScalarsVector, center);
	}

}
















