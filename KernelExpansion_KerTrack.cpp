// @author Merve Asiler

#include "KernelExpansion_KerTrack.h"
#include "BaseGeoOpUtils.h"

KernelExpansion_KerTrack::KernelExpansion_KerTrack(const Mesh& hostMesh) :
	KernelExpansion(hostMesh) {};

KernelExpansion_KerTrack::~KernelExpansion_KerTrack() {
}

void KernelExpansion_KerTrack::expandKernel() {

	//cout << "Number of bounding planes: " << halfSpaceSet.size() << " out of " << hostMeshptr->getNumOfTris() << " triangles." << endl;

	for (int i = 0; i < halfSpaceSet.size(); i++) {
		edgePartners.push_back(queue<EdgePartnerTriple*>());
		isKernelFace.push_back(ProcessColor::RED);
		edgePartnerIds.push_back(vector<int>());
	}

	clock_t begin = clock();
	double point[3];
	findInitialPoint_3(point);
	clock_t end = clock();
	cout << "Point finding time: " << double(end - begin) / CLOCKS_PER_SEC << endl;
	if (point[0] == numeric_limits<double>::infinity()) {
		kernel = NULL;
		return;	// NOT STAR-SHAPED!
	}

	double* newpoint1, *newpoint2;
	int theClosestId1, theClosestId2;
	tie(theClosestId1, newpoint1) = findTheClosestHalfSpace(point);
	tie(theClosestId2, newpoint2) = findTheClosestHalfSpace(newpoint1, theClosestId1);
	delete[] newpoint1;
	delete[] newpoint2;
	int number_of_kernel_faces = 0;
	
	for (bool any_white_left = false; ; any_white_left = false) {
		for (int i = 0; i < halfSpaceSet.size(); i++) {
			if (isKernelFace[i] == ProcessColor::WHITE) {
				while(!edgePartners[i].empty()) {
					// initializations
					EdgePartnerTriple* ept = edgePartners[i].front();
					int partnerId = ept->partnerId;
					double* edgeDirection = ept->edgeDirection;
					double* startPoint = ept->startPoint;
					int startPointId = ept->startPointId;
					edgePartners[i].pop();

					// find the corner on the given edge
					vector<int> theClosestIds;
					double* kernelVertex;
					int previousNumOfVerts = kernel->getNumOfVerts();
					tie(theClosestIds, kernelVertex) = findTheClosestHalfSpace(startPointId, edgeDirection, i, partnerId);	

					if (kernel->getNumOfVerts() > previousNumOfVerts || isTripleSame(kernelVertex, startPoint))
						orderTheFaces(i, partnerId, theClosestIds, kernelVertex, edgeDirection);
					
					//if (startPointId >= 0)
					//	kernel->addEdge(startPointId, kernel->getNumOfVerts() - 1);

					delete[] kernelVertex;
					delete ept;
				}
				any_white_left = true;
				isKernelFace[i] = ProcessColor::GREEN;
				number_of_kernel_faces++;
				break;
			}
		}
		if (!any_white_left)
			break;
	}

	/*
	int num_of_non_kernel_points = 0;
	for (int i = 0; i < kernel->getNumOfVerts(); i++) {
		if (findClosestValueSatisfiedByPoint(kernel->getVertex(i)->coords, halfSpaceSet) > 3 * EPSILON ||
			findClosestValueSatisfiedByPoint(kernel->getVertex(i)->coords, halfSpaceSet) > 3 * EPSILON ||
			findClosestValueSatisfiedByPoint(kernel->getVertex(i)->coords, halfSpaceSet) > 3 * EPSILON) {
			//cout << "non-kernel point" << endl;
			num_of_non_kernel_points++;
		}
	}
	cout << "num of non-kernel points is: " << num_of_non_kernel_points << endl;
	*/

}

tuple<int, double*> KernelExpansion_KerTrack::findTheClosestHalfSpace(double point[3]) {
	
	double* scalarsArray = findValueVectorSatisfiedByPoint(point, halfSpaceSet);
	vector<double> scalarsVector;
	filterRepetitions(scalarsArray, scalarsVector);

	double theClosestDistance = numeric_limits<double>::infinity();
	int theClosestId = -1;
	for (int i = 0; i < halfSpaceSet.size(); i++) {
		if (abs(scalarsVector[i]) < theClosestDistance) {
			theClosestId = i;
			theClosestDistance = abs(scalarsVector[i]);
		}
	}

	scalarsVector.clear();
	double* newpoint = new double[3];
	for (int i = 0; i < 3; i++)
		newpoint[i] = point[i] + halfSpaceSet[theClosestId].ABCD[i] * theClosestDistance;
	return make_tuple(theClosestId, newpoint);
}

tuple<int, double*> KernelExpansion_KerTrack::findTheClosestHalfSpace(double point[3], int id) {

	double theClosestDistance = numeric_limits<double>::infinity();
	int theClosestId = -1;
	double* theClosestDirection = NULL, *intersectionLineDirection = NULL;
	for (int s = 0, i = 0, bound = id; s < 2; s++) {
		for (; i < bound; i++) {
			
			double* lineDirection = crossProduct(halfSpaceSet[i].ABCD, halfSpaceSet[id].ABCD);
			normalize(lineDirection);
			double* rayDirection = crossProduct(halfSpaceSet[id].ABCD, lineDirection);
			normalize(rayDirection);
			Line ray(point, rayDirection);
			double t = findLinePlaneIntersection(ray, halfSpaceSet[i]);
			if (t < theClosestDistance) {
				theClosestId = i;
				theClosestDistance = t;
				if (theClosestDirection) {
					delete[] theClosestDirection;
					delete[] intersectionLineDirection;
				}
				theClosestDirection = rayDirection;
				intersectionLineDirection = lineDirection;
			}
			else {
				delete[] rayDirection;
				delete[] lineDirection;
			}
		}
		i = id + 1;
		bound = halfSpaceSet.size();
	}

	double* newpoint = new double[3];
	for (int i = 0; i < 3; i++)
		newpoint[i] = point[i] + theClosestDirection[i] * theClosestDistance;
	delete[] theClosestDirection;

	edgePartners[id].push(new EdgePartnerTriple(theClosestId, intersectionLineDirection, newpoint, 0));
	edgePartnerIds[id].push_back(theClosestId);
	edgePartnerIds[theClosestId].push_back(id);
	isKernelFace[id] = ProcessColor::WHITE;
	isKernelFace[theClosestId] = ProcessColor::WHITE;
	
	vector<int> parentIdsForNewVertex;
	parentIdsForNewVertex.push_back(id);
	parentIdsForNewVertex.push_back(theClosestId);
	vertexParentIds.push_back(parentIdsForNewVertex);
	kernel->addVertex(newpoint[0], newpoint[1], newpoint[2]);

	delete[] intersectionLineDirection;
	return make_tuple(theClosestId, newpoint);
}

tuple<vector<int>, double*> KernelExpansion_KerTrack::findTheClosestHalfSpace(int vertexId, double* lineDirection, int lineParent1Id, int lineParent2Id) {

	double theClosestDistance = numeric_limits<double>::infinity();
	vector<int> theClosestIds;
	Line line(kernel->getVertex(vertexId)->coords, lineDirection);
	
	for (int i = 0; i < halfSpaceSet.size(); i++) {
		bool previous_parent = false;
		for (int j = 0; j < vertexParentIds[vertexId].size(); j++)
			if (i == vertexParentIds[vertexId][j]) {
				previous_parent = true;
				break;
			}
		if (previous_parent)
			continue;

		double t = findLinePlaneIntersection(line, halfSpaceSet[i]);
		if (t >= 0 && t <= theClosestDistance) {
			if (t < theClosestDistance)
				theClosestIds.clear();
			theClosestIds.push_back(i);
			theClosestDistance = t;
		}
	}

	double* newpoint = new double[3];
	for (int i = 0; i < 3; i++)
		newpoint[i] = kernel->getVertex(vertexId)->coords[i] + lineDirection[i] * theClosestDistance;

	if (theClosestDistance > 0) {
		bool previousVertex = false;
		for (int v = 0; v < kernel->getNumOfVerts(); v++)
			if (isTripleSame(kernel->getVertex(v)->coords, newpoint)) {
				previousVertex = true;
				vertexId = v;
				break;
			}
		if (!previousVertex) {
			kernel->addVertex(newpoint[0], newpoint[1], newpoint[2]);
			vector<int> parentIdsForNewVertex;
			parentIdsForNewVertex.push_back(lineParent1Id);
			parentIdsForNewVertex.push_back(lineParent2Id);
			vertexParentIds.push_back(parentIdsForNewVertex);
			vertexId = kernel->getNumOfVerts() - 1;
		}
	}

	for (int i = 0; i < theClosestIds.size(); i++)
		vertexParentIds[vertexId].push_back(theClosestIds[i]);

	return make_tuple(theClosestIds, newpoint);

}

void KernelExpansion_KerTrack::orderTheFaces(int base_id, int partner_id, vector<int> next_partner_ids, double* startPoint, double* currentEdgeDirection) {

	vector<double*> edgeDirections;
	vector<double> cosAngles;
	bool should_revert = false;

	double* currentEdgeDirection2 = crossProduct(halfSpaceSet[base_id].ABCD, halfSpaceSet[partner_id].ABCD);
	if (dotProduct(currentEdgeDirection, currentEdgeDirection2) < 0) {
		delete[] currentEdgeDirection2;
		currentEdgeDirection2 = multVect(currentEdgeDirection, -1.0);
		currentEdgeDirection = currentEdgeDirection2;
		should_revert = true;
	}
	else
		delete[] currentEdgeDirection2;

	// ordering
	for (int i = 0; i < next_partner_ids.size(); i++) {
		double* nextEdgeDirection = crossProduct(halfSpaceSet[base_id].ABCD, halfSpaceSet[next_partner_ids[i]].ABCD);
		normalize(nextEdgeDirection);
		edgeDirections.push_back(nextEdgeDirection);
		double cosAngle = dotProduct(currentEdgeDirection, nextEdgeDirection);
		cosAngles.push_back(cosAngle);

		for (int j = 0; j < cosAngles.size(); j++) {
			if (cosAngle < cosAngles[j]) {
				for (int k = cosAngles.size() - 1; k >= j; k--) {
					cosAngles[k] = cosAngles[k - 1];
					edgeDirections[k] = edgeDirections[k - 1];
				}
				cosAngles[j] = cosAngle;
				edgeDirections[j] = nextEdgeDirection;
				break;
			}
		}
	}
	
	if (should_revert) {
		for (int i = 0; i < next_partner_ids.size(); i++) {
			double* reversedEdgeDirection = multVect(edgeDirections[i], -1.0);
			delete edgeDirections[i];
			edgeDirections[i] = reversedEdgeDirection;
		}

	}

	if (shouldSaveEdge(base_id, next_partner_ids[0]))
		edgePartners[base_id].push(new EdgePartnerTriple(next_partner_ids[0], edgeDirections[0], startPoint, kernel->getNumOfVerts() - 1));
	
	next_partner_ids.push_back(partner_id);
	for (int i = 0, j = 1; j < next_partner_ids.size(); j++) {
		if (next_partner_ids[j] == -1)
			continue;
		if (shouldSaveEdge(next_partner_ids[i], next_partner_ids[j]))
			saveFoundEdgeToProcess(next_partner_ids[i], next_partner_ids[j], startPoint, halfSpaceSet[base_id].ABCD);
		i = j;
	}

	edgeDirections.clear();
	cosAngles.clear();
}

bool KernelExpansion_KerTrack::shouldSaveEdge(int hs1_id, int hs2_id) {

	bool walkedEdge = false;
	for (int j = 0; j < edgePartnerIds[hs1_id].size(); j++)
		if (edgePartnerIds[hs1_id][j] == hs2_id) {
			walkedEdge = true;
			break;
		}
	if (!walkedEdge) {
		edgePartnerIds[hs1_id].push_back(hs2_id);
		edgePartnerIds[hs2_id].push_back(hs1_id);
		if (isKernelFace[hs1_id] != ProcessColor::GREEN)
			isKernelFace[hs1_id] = ProcessColor::WHITE;
		if (isKernelFace[hs2_id] != ProcessColor::GREEN)
			isKernelFace[hs2_id] = ProcessColor::WHITE;
		return true;
	}
	return false;

}

void KernelExpansion_KerTrack::saveFoundEdgeToProcess(int base_id, int partner_id, double* startPoint, double* directioner) {

	double* edgeDirection = crossProduct(halfSpaceSet[base_id].ABCD, halfSpaceSet[partner_id].ABCD);
	normalize(edgeDirection);
	double* nextEdgeDirection = new double[3]{edgeDirection[0], edgeDirection[1], edgeDirection[2] };
	if (dotProduct(edgeDirection, directioner) > 0)
		nextEdgeDirection = multVect(edgeDirection, -1);

	edgePartners[base_id].push(new EdgePartnerTriple(partner_id, nextEdgeDirection, startPoint, kernel->getNumOfVerts() - 1));
	delete[] nextEdgeDirection;
	delete[] edgeDirection;

}

double* KernelExpansion_KerTrack::findInitialPoint_1() {
		
	if (this->initialPoint) {
		// find the center of the kernel's bonding box
		double* point = new double[3]{ 0, 0, 0 };
		for (int i = 0; i < 6; i++)
			point[i / 2] += this->initialPoint[i * 3 + i / 2];
		for (int i = 0; i < 3; i++)
			point[i] /= 2.0;
		return point;
	}
	return NULL;
	
}

double* KernelExpansion_KerTrack::findInitialPoint_2() {

	Vertex* v = hostMeshptr->getVertex(0);
	double* point = new double[3]{ v->coords[0], v->coords[1], v->coords[2] };
	vector<int> satisfied_face_ids;

	for (int i = 0; i < halfSpaceSet.size(); i++) {
		double distance = halfSpaceSet[i].ABCD[3];
		for (int k = 0; k < 3; k++)
			distance += halfSpaceSet[i].ABCD[k] * point[k];
		if (distance <= 0)
			satisfied_face_ids.push_back(i);
		else {
			for (int j = 0; j < satisfied_face_ids.size(); j++) {
				double in_distance = halfSpaceSet[satisfied_face_ids[j]].ABCD[3];
				for (int k = 0; k < 3; k++)
					in_distance += halfSpaceSet[satisfied_face_ids[j]].ABCD[k] * point[k];
				double in_point[3];
				for (int k = 0; k < 3; k++)
					in_point[k] = halfSpaceSet[satisfied_face_ids[j]].ABCD[k] * in_distance + point[k];
				for (int c = 0; c < satisfied_face_ids.size(); c++) {
					double check_distance = halfSpaceSet[satisfied_face_ids[c]].ABCD[3];
					for (int k = 0; k < 3; k++)
						check_distance += halfSpaceSet[satisfied_face_ids[c]].ABCD[k] * in_point[k];
					if (check_distance > 0) {
						for (int k = 0; k < 3; k++)
							in_point[k] = -halfSpaceSet[satisfied_face_ids[c]].ABCD[k] * check_distance + in_point[k];
					}
				}
				double* in_direction = diffVects(in_point, point);
				in_distance = computeLength(in_direction);
				normalize(in_direction);
				distance = -halfSpaceSet[i].ABCD[3];
				double payda = 0;
				for (int k = 0; k < 3; k++) {
					distance -= halfSpaceSet[i].ABCD[k] * point[k];
					payda += halfSpaceSet[i].ABCD[k] * in_direction[k];
				}
				distance /= payda;
				if (distance < in_distance)
					for (int k = 0; k < 3; k++)
						point[k] += in_direction[k] * distance;
				delete[] in_direction;
			}
		}
	}

	return point;
}

void KernelExpansion_KerTrack::findInitialPoint_3(double* point) {

	Vertex* v = hostMeshptr->getVertex(0);
	for (int k = 0; k < 3; k++)
		point[k] = v->coords[k];

	for (int i = 0; i < halfSpaceSet.size(); i++) {
		// calculate the distance of current point to i'th plane
		double distance = halfSpaceSet[i].ABCD[3];
		for (int k = 0; k < 3; k++)
			distance += halfSpaceSet[i].ABCD[k] * point[k];
		// if distance > 0, then it means point does not satisfy the i'th halfspace
		if (distance > 0) {
			// project the point onto i'th plane, now it is the new point
			for (int k = 0; k < 3; k++)
				point[k] = -halfSpaceSet[i].ABCD[k] * distance + point[k];

			// we made the previous point satisfy all the halfpaces upto i'th, check them for the new point
			for (int j = 0; j < i; j++) {
				// calculate the distance of the new point to j'th plane
				double distance = halfSpaceSet[j].ABCD[3];
				for (int k = 0; k < 3; k++)
					distance += halfSpaceSet[j].ABCD[k] * point[k];
				// if distance > 0, then it means that new point does not satisfy the j'th halfspace
				// but the previous point was satisfying j'th half-space
				// then we understand that j'th halfspace passes through between the previous point and the new point (current point)
				if (distance > EPSILON) {
					double lineDir[3], goingDir[3];
					crossProduct(halfSpaceSet[j].ABCD, halfSpaceSet[i].ABCD, lineDir);
					crossProduct(lineDir, halfSpaceSet[i].ABCD, goingDir);
					Line goingLine(point, goingDir);
					double t = findLinePlaneIntersection(goingLine, halfSpaceSet[j]);

					double zeroVector[3] = { 0, 0, 0 };
					if (isTripleSame(lineDir, zeroVector) || isTripleSame(goingDir, zeroVector) || t == numeric_limits<double>::infinity()) {
						for (int k = 0; k < 3; k++)
							point[k] = numeric_limits<double>::infinity();
						return;
					}

					for (int k = 0; k < 3; k++)
						point[k] = goingDir[k] * t + point[k];
				
					for (int m = 0; m < j; m++) {
						double distance = halfSpaceSet[m].ABCD[3];
						for (int k = 0; k < 3; k++)
							distance += halfSpaceSet[m].ABCD[k] * point[k];
						if (distance > 0) {
							for (int k = 0; k < 3; k++)
								point[k] = numeric_limits<double>::infinity();
							return;
						}
					}
				}
			}
		}
	}
}

double* KernelExpansion_KerTrack::findInitialPoint_4() {
	
	Vertex* v = hostMeshptr->getVertex(0);
	double* point = new double[3]{ v->coords[0], v->coords[1], v->coords[2] };
	for (int i = 0; i < halfSpaceSet.size(); i++) {
		// calculate the distance of current point to i'th plane
		double distance = halfSpaceSet[i].ABCD[3];
		for (int k = 0; k < 3; k++)
			distance += halfSpaceSet[i].ABCD[k] * point[k];
		// if distance > 0, then it means point does not satisfy the i'th halfspace
		if (distance > 0) {
			// project the point onto i'th plane, now it is the new point
			for (int k = 0; k < 3; k++)
				point[k] = -halfSpaceSet[i].ABCD[k] * distance + point[k];

			// we made the previous point satisfy all the halfpaces upto i'th, check them for the new point
			double therightmostDist = -numeric_limits<double>::infinity();
			double upperLimit = numeric_limits<double>::infinity();
			for (int j = 0; j < i; j++) {
				// calculate the distance of the new point to j'th plane
				double distance = halfSpaceSet[j].ABCD[3];
				for (int k = 0; k < 3; k++)
					distance += halfSpaceSet[j].ABCD[k] * point[k];
				// if distance > 0, then it means that new point does not satisfy the j'th halfspace
				// but the previous point was satisfying j'th half-space
				// then we understand that j'th halfspace passes through between the previous point and the new point (current point)
				if (distance > EPSILON) {
					double cosAngle = dotProduct(halfSpaceSet[j].ABCD, halfSpaceSet[i].ABCD);
					double tanAngle = tan(acos(cosAngle));
					Line upLine(point, halfSpaceSet[i].ABCD);
					double upDist = findLinePlaneIntersection(upLine, halfSpaceSet[j]);
					double rightDist;
					if (upDist == numeric_limits<double>::infinity()) {
						Line rightLine(point, halfSpaceSet[j].ABCD);
						rightDist = findLinePlaneIntersection(rightLine, halfSpaceSet[j]);
						if (rightDist < upperLimit)
							upperLimit = rightDist;
						else
							continue;
					}
					else
						rightDist = upDist / -tanAngle;

					if (rightDist > therightmostDist)
						therightmostDist = rightDist;
				}
			}
			for (int k = 0; k < 3; k++)
				point[k] += halfSpaceSet[i].ABCD[k] * therightmostDist;
		}
	}

	for (int i = 0; i < halfSpaceSet.size(); i++) {
		double distance = halfSpaceSet[i].ABCD[3];
		for (int k = 0; k < 3; k++)
			distance += halfSpaceSet[i].ABCD[k] * point[k];
		if (distance > 0) {
			cout << "NOT STAR-SHAPED!!!" << endl;
			return NULL;
		}
	}
	return point;
	
}

void KernelExpansion_KerTrack::filterRepetitions(double* distances,  vector<double>& scalarsVector) {

	double _EPSILON = EPSILON * 3;
	vector<HalfSpace> temp_halfSpaceSet;

	for (int i = 0; i < halfSpaceSet.size(); i++) {
		bool isTheSame = false;
		for (int j = 0; j < temp_halfSpaceSet.size(); j++) {
			if (abs(distances[i] - scalarsVector[j]) < _EPSILON) {
				if (isTripleSame(halfSpaceSet[i].ABCD, temp_halfSpaceSet[j].ABCD)) {
					isTheSame = true;
					break;
				}
			}
		}
		if (!isTheSame) {
			temp_halfSpaceSet.push_back(halfSpaceSet[i]);
			scalarsVector.push_back(distances[i]);
		}
	}

	halfSpaceSet.clear();
	delete[] distances;
	distances = NULL;
	for (int j = 0; j < temp_halfSpaceSet.size(); j++)
		halfSpaceSet.push_back(temp_halfSpaceSet[j]);
	temp_halfSpaceSet.clear();

}
