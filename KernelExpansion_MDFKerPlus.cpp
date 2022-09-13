// @author Merve Asiler

#include "KernelExpansion_MDFKerPlus.h"
#include "BaseGeoOpUtils.h"

KernelExpansion_MDFKerPlus::KernelExpansion_MDFKerPlus(const Mesh& hostMesh, double cellSizeRatio) :
	KernelExpansion(hostMesh, cellSizeRatio) {};

KernelExpansion_MDFKerPlus::~KernelExpansion_MDFKerPlus() {
	for (int i = 0; i < allTriedPlaneCombinations.size(); i++)
		delete[] allTriedPlaneCombinations[i];
	allTriedPlaneCombinations.clear();
}

void KernelExpansion_MDFKerPlus::expandKernel() {

	// is kernel empty?
	if (this->extremeCorners[0] == nullptr && this->extremeCorners[1] == nullptr) {
		delete kernel;
		kernel = NULL;
		return;	// NOT STAR-SHAPED!
	}

	this->initialPoint = new double[3];
	for (int i = 0; i < 3; i++)
		this->initialPoint[i] = (grid.getMinCoords()[i] + grid.getMaxCoords()[i]) / 2.0;

	double loop_limits[2] = { grid.getMaxCoords()[0], grid.getMaxCoords()[1] };
	double loop_test_elements[2] = { this->initialPoint[0], this->initialPoint[1] };
	double loop_step_sizes[2] = { grid.cellSize[0], grid.cellSize[1] };
	double* prevScalarsVector = NULL;
	double prevCoords[3];

	for (int i = 0; i < 2; i++) {
		for (double x = this->initialPoint[0]; loop_test_elements[0] <= loop_limits[0]; loop_test_elements[0] += abs(loop_step_sizes[0]), x += loop_step_sizes[0]) {

			for (int j = 0; j < 2; j++) {
				for (double y = this->initialPoint[1]; loop_test_elements[1] <= loop_limits[1]; loop_test_elements[1] += abs(loop_step_sizes[1]), y += loop_step_sizes[1]) {
					double z = this->initialPoint[2];
					double coords[3] = { x, y, z };
					double direction[3] = { 0, 0, 1.0 };
					sendRay(coords, direction);
				}

				loop_test_elements[1] = -initialPoint[1];
				loop_limits[1] = -grid.getMinCoords()[1];
				loop_step_sizes[1] = -grid.cellSize[1];
			}

		}
		loop_test_elements[1] = initialPoint[1];
		loop_limits[1] = grid.getMaxCoords()[1];
		loop_step_sizes[1] = grid.cellSize[1];
		loop_test_elements[0] = -initialPoint[0];
		loop_limits[0] = -grid.getMinCoords()[0];
		loop_step_sizes[0] = -grid.cellSize[0];
	}

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

void KernelExpansion_MDFKerPlus::sendRay(double source[3], double direction[3]) {

	// compute the distances to each halfspace at the given direction
	// detect the farthest excluder halfspace
	double theFarthestExcluderDistance = -numeric_limits<double>::infinity();
	int theFarthestExcluderId = -1;
	double* distances = new double[halfSpaceSet.size()];
	double POS_EPSILON = 3 * EPSILON;
	double NEG_EPSILON = -POS_EPSILON;
	Line* ray = new Line(source, direction);
	for (int i = 0; i < halfSpaceSet.size(); i++) {
		distances[i] = findLinePlaneIntersection(ray, halfSpaceSet[i]);
		if (distances[i] > POS_EPSILON && dotProduct(direction, halfSpaceSet[i]->ABCD) < NEG_EPSILON) {
			if (distances[i] > theFarthestExcluderDistance) {
				theFarthestExcluderDistance = distances[i];
				theFarthestExcluderId = i;
			}
		}
		else if (distances[i] < NEG_EPSILON && dotProduct(direction, halfSpaceSet[i]->ABCD) > POS_EPSILON) {
			if (abs(distances[i]) > theFarthestExcluderDistance) {
				theFarthestExcluderDistance = -distances[i];
				theFarthestExcluderId = i;
			}
		}
	}

	delete ray;
	double point1[3], point2[3];

	// project onto kernel border
	if (theFarthestExcluderId >= 0) {
		for (int i = 0; i < halfSpaceSet.size(); i++)
			distances[i] -= distances[theFarthestExcluderId];
		for (int i = 0; i < 3; i++)
			point1[i] = source[i] + direction[i] * distances[theFarthestExcluderId];

		double theClosestDistance;
		if (distances[theFarthestExcluderId] > 0) {
			theClosestDistance = numeric_limits<double>::infinity();
			for (int i = 0; i < halfSpaceSet.size(); i++)
				if (distances[i] > POS_EPSILON && distances[i] < theClosestDistance)
					theClosestDistance = distances[i];
		}
		else {
			theClosestDistance = -numeric_limits<double>::infinity();
			for (int i = 0; i < halfSpaceSet.size(); i++)
				if (distances[i] < NEG_EPSILON && distances[i] > theClosestDistance)
					theClosestDistance = distances[i];
		}
		for (int i = 0; i < 3; i++)
			point2[i] = point1[i] + direction[i] * theClosestDistance;
	}

	else {
		double theClosestForwardDistance = numeric_limits<double>::infinity();
		double theClosestBackwardDistance = -numeric_limits<double>::infinity();
		int forwardHalfSpaceId = -1, backwardHalfSpaceId = -1, onHalfSpaceId = -1;
		for (int i = 0; i < halfSpaceSet.size(); i++) {
			if (distances[i] > POS_EPSILON && distances[i] < theClosestForwardDistance)
				theClosestForwardDistance = distances[i];
			else if (distances[i] < NEG_EPSILON && distances[i] > theClosestBackwardDistance)
				theClosestBackwardDistance = distances[i];
			else if (abs(distances[i]) <= POS_EPSILON)
				;	//?
		}
		for (int i = 0; i < 3; i++) {
			point1[i] = source[i] + direction[i] * theClosestForwardDistance;
			point2[i] = source[i] + direction[i] * theClosestBackwardDistance;
		}
	}

	delete[] distances;

	double* scalarsVector = findValueVectorSatisfiedByPoint(point1, halfSpaceSet);
	findPossibleKernelVertex(scalarsVector, point1);

	if (abs(point1[0] - point2[0]) < EPSILON && abs(point1[1] - point2[1]) < EPSILON && abs(point1[2] - point2[2]) < EPSILON)
		return;

	scalarsVector = findValueVectorSatisfiedByPoint(point2, halfSpaceSet);
	findPossibleKernelVertex(scalarsVector, point2);

}

void KernelExpansion_MDFKerPlus::findPossibleKernelVertex(double* scalarsVector, double source_coords[3]) {

	int closestPlaneIds[3] = { -1, -1, -1 };
	this->updateClosest3Planes(scalarsVector, closestPlaneIds, closestPlaneIds);							// find the planes satisfying the minimum 3 distance value by selecting from scalarsVector
	double* candidateKernelVertex = computePossibleKernelVertex(scalarsVector, closestPlaneIds);	// however those 3 planes must be able to create a vertex
	delete[] scalarsVector;

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
void KernelExpansion_MDFKerPlus::updateClosest3Planes(double* values, int currentPlaneIds[3], int previousPlaneIds[3]) {

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

double* KernelExpansion_MDFKerPlus::computePossibleKernelVertex(double* scalarsVector, int closestPlaneIds[3]) {

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

double* KernelExpansion_MDFKerPlus::computePossibleKernelVertex(int closestPlaneIds[3]) {

	Plane* closestPlanes[3] = { halfSpaceSet[closestPlaneIds[0]], halfSpaceSet[closestPlaneIds[1]], halfSpaceSet[closestPlaneIds[2]] };
	double* candidateKernelVertex = find3PlaneIntersection(closestPlanes);
	return candidateKernelVertex;
}

int KernelExpansion_MDFKerPlus::addIfNoFrontierObstacle(double candidateKernelVertex[3], int candidateVertexProducerIds[3], double source_coords[3]) {

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

bool KernelExpansion_MDFKerPlus::isAlreadyKernelVertex(double coords[3]) {

	for (int i = 0; i < kernel->getNumOfVerts(); i++) {
		Vertex* v = kernel->getVertex(i);
		if (abs(v->coords[0] - coords[0]) < EPSILON && abs(v->coords[1] - coords[1]) < EPSILON && abs(v->coords[2] - coords[2]) < EPSILON)
			return true;
	}
	return false;
}

int KernelExpansion_MDFKerPlus::findIfFrontierObstacle(double* candidateKernelVertex, int candidateVertexProducerIds[3], double* sourcePoint) {

	int frontierHalfSpaceId = -1;

	double* direction = new double[3];
	for (int i = 0; i < 3; i++)
		direction[i] = candidateKernelVertex[i] - sourcePoint[i];
	double candidateDistance = computeLength(direction);
	normalize(direction);

	Line* line = new Line(sourcePoint, direction);
	delete[] direction;

	for (int i = 0; i < halfSpaceSet.size(); i++) {
		if (i == candidateVertexProducerIds[0] || i == candidateVertexProducerIds[1] || i == candidateVertexProducerIds[2])
			continue;
		double distance = findLinePlaneIntersection(line, halfSpaceSet[i]);
		if (distance >= 0 && (distance - candidateDistance) < EPSILON) {
			frontierHalfSpaceId = i;
			candidateDistance = distance;
		}
	}

	delete line;
	return frontierHalfSpaceId;
}

void KernelExpansion_MDFKerPlus::findByFrontierPlanes(int closestPlaneIds[3], int obstaclePlaneId, double source_coords[3]) {

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
}

void KernelExpansion_MDFKerPlus::replaceOnePlaneBy(queue<int*>& possiblePlaneCombinations, int planeIds[3], int replacerId) {

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

void KernelExpansion_MDFKerPlus::decideOnNeighborCells(int* cell_indices, double* corner_in_outs, double** scalarsVector) {

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

void KernelExpansion_MDFKerPlus::addNeighborCellIntoQueue(int neigh_i, int neigh_j, int neigh_k,
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



void KernelExpansion_MDFKerPlus::findMissingKernelVertices(Mesh* incompleteKernel, double cellSize) {

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
















