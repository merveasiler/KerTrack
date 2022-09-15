// @author Merve Asiler

#include "KernelExpansion_KerTrack.h"
#include "BaseGeoOpUtils.h"

KernelExpansion_KerTrack::KernelExpansion_KerTrack(const Mesh& hostMesh) :
	KernelExpansion(hostMesh) {};

KernelExpansion_KerTrack::~KernelExpansion_KerTrack() {
}

void KernelExpansion_KerTrack::expandKernel() {

	//cout << "Number of bounding planes: " << halfSpaceSet.size() << " out of " << hostMeshptr->getNumOfTris() << " triangles." << endl;

	double point[3];
	findInitialPoint_3(point);
	if (point[0] == numeric_limits<double>::infinity())
		return;	// NOT STAR-SHAPED!

	vector<double> scalarsVector = initialize(point);
	int theClosestId1 = findTheClosestHalfSpace(point, scalarsVector);
	int theClosestId2 = findTheClosestHalfSpace(point, theClosestId1);
	int number_of_kernel_faces = 0;
	
	for (bool any_white_left = false; ; any_white_left = false) {
		for (int i = 0; i < halfSpaceSet.size(); i++) {
			if (isKernelFace[i] == ProcessColor::WHITE) {
				while(!edgePartners[i].empty()) {
					// initializations
					EdgePartnerTriple& ept = edgePartners[i].front();
					int partnerId = ept.partnerId;
					double edgeDirection[3] = { ept.edgeDirection[0], ept.edgeDirection[1], ept.edgeDirection[2] };
					double startPoint[3] = { ept.startPoint[0], ept.startPoint[1], ept.startPoint[2] };
					int startPointId = ept.startPointId;
					edgePartners[i].pop();

					// find the corner on the given edge
					int previousNumOfVerts = kernel.getNumOfVerts();
					double kernelVertex[3];
					vector<int> theClosestIds = findTheClosestHalfSpace(startPointId, edgeDirection, i, partnerId, kernelVertex);	

					if (kernel.getNumOfVerts() > previousNumOfVerts || isTripleSame(kernelVertex, startPoint))
						orderTheFaces(i, partnerId, theClosestIds, kernelVertex, edgeDirection);

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

	
	int num_of_non_kernel_points = 0;
	for (int i = 0; i < kernel.getNumOfVerts(); i++) {
		if (findClosestValueSatisfiedByPoint(kernel.getVertex(i).coords, halfSpaceSet) > 3 * EPSILON ||
			findClosestValueSatisfiedByPoint(kernel.getVertex(i).coords, halfSpaceSet) > 3 * EPSILON ||
			findClosestValueSatisfiedByPoint(kernel.getVertex(i).coords, halfSpaceSet) > 3 * EPSILON) {
			num_of_non_kernel_points++;
		}
	}
	cout << "num of non-kernel points is: " << num_of_non_kernel_points << endl;
	

}

vector<double> KernelExpansion_KerTrack::initialize(double* point) {

	double* scalarsArray = findValueVectorSatisfiedByPoint(point, halfSpaceSet);
	vector<double> scalarsVector;
	filterRepetitions(scalarsArray, scalarsVector);

	for (int i = 0; i < halfSpaceSet.size(); i++) {
		edgePartners.push_back(queue<EdgePartnerTriple>());
		isKernelFace.push_back(ProcessColor::RED);
		edgePartnerIds.push_back(vector<int>());
	}
	
	return scalarsVector;
}

int KernelExpansion_KerTrack::findTheClosestHalfSpace(double* point, vector<double>& scalarsVector) {

	double theClosestDistance = numeric_limits<double>::infinity();
	int theClosestId = -1;
	for (int i = 0; i < halfSpaceSet.size(); i++) {
		if (abs(scalarsVector[i]) < theClosestDistance) {
			theClosestId = i;
			theClosestDistance = abs(scalarsVector[i]);
		}
	}

	scalarsVector.clear();
	for (int i = 0; i < 3; i++)
		point[i] = point[i] + halfSpaceSet[theClosestId].ABCD[i] * theClosestDistance;
	return theClosestId;
}

int KernelExpansion_KerTrack::findTheClosestHalfSpace(double* point, int id) {

	double theClosestDistance = numeric_limits<double>::infinity();
	int theClosestId = -1;
	double theClosestDirection[3], intersectionLineDirection[3];;
	for (int s = 0, i = 0, bound = id; s < 2; s++) {
		for (; i < bound; i++) {
			
			double lineDirection[3], rayDirection[3];
			crossProduct(halfSpaceSet[i].ABCD, halfSpaceSet[id].ABCD, lineDirection);
			normalize(lineDirection);
			crossProduct(halfSpaceSet[id].ABCD, lineDirection, rayDirection);
			normalize(rayDirection);
			Line ray(point, rayDirection);
			double t = findLinePlaneIntersection(ray, halfSpaceSet[i]);
			if (t < theClosestDistance) {
				theClosestId = i;
				theClosestDistance = t;
				for (int k = 0; k < 3; k++) {
					theClosestDirection[k] = rayDirection[k];
					intersectionLineDirection[k] = lineDirection[k];
				}
			}
		}
		i = id + 1;
		bound = halfSpaceSet.size();
	}

	for (int i = 0; i < 3; i++)
		point[i] = point[i] + theClosestDirection[i] * theClosestDistance;

	edgePartners[id].push(EdgePartnerTriple(theClosestId, intersectionLineDirection, point, 0));
	edgePartnerIds[id].push_back(theClosestId);
	edgePartnerIds[theClosestId].push_back(id);
	isKernelFace[id] = ProcessColor::WHITE;
	isKernelFace[theClosestId] = ProcessColor::WHITE;
	
	vector<int> parentIdsForNewVertex;
	parentIdsForNewVertex.push_back(id);
	parentIdsForNewVertex.push_back(theClosestId);
	vertexParentIds.push_back(parentIdsForNewVertex);
	kernel.addVertex(point[0], point[1], point[2]);

	return theClosestId;
}

vector<int> KernelExpansion_KerTrack::findTheClosestHalfSpace(int vertexId, double* lineDirection, int lineParent1Id, int lineParent2Id, double* newpoint) {

	double theClosestDistance = numeric_limits<double>::infinity();
	vector<int> theClosestIds;
	Line line(kernel.getVertex(vertexId).coords, lineDirection);
	
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

	for (int i = 0; i < 3; i++)
		newpoint[i] = kernel.getVertex(vertexId).coords[i] + lineDirection[i] * theClosestDistance;
	
	int startId = vertexId;
	if (theClosestDistance > 0) {
		bool previousVertex = false;
		for (int v = 0; v < kernel.getNumOfVerts(); v++)
			if (isTripleSame(kernel.getVertex(v).coords, newpoint)) {
				previousVertex = true;
				vertexId = v;
				break;
			}
		if (!previousVertex) {
			kernel.addVertex(newpoint[0], newpoint[1], newpoint[2]);
			vector<int> parentIdsForNewVertex;
			parentIdsForNewVertex.push_back(lineParent1Id);
			parentIdsForNewVertex.push_back(lineParent2Id);
			vertexParentIds.push_back(parentIdsForNewVertex);
			vertexId = kernel.getNumOfVerts() - 1;
		}
	}

	for (int i = 0; i < theClosestIds.size(); i++)
		vertexParentIds[vertexId].push_back(theClosestIds[i]);

	kernel.addEdge(startId, vertexId);

	return theClosestIds;

}

void KernelExpansion_KerTrack::orderTheFaces(int base_id, int partner_id, vector<int> next_partner_ids, double* startPoint, double* currentEdgeDirection) {

	vector<int> ordered_next_partner_ids;
	vector<double> cosAngles;
	bool should_revert = false;

	double currentEdgeDirection2[3];
	crossProduct(halfSpaceSet[base_id].ABCD, halfSpaceSet[partner_id].ABCD, currentEdgeDirection2);
	if (dotProduct(currentEdgeDirection, currentEdgeDirection2) < 0) {
		multVect(currentEdgeDirection, -1.0, currentEdgeDirection2);
		for (int k=0; k<3; k++)
			currentEdgeDirection[k] = currentEdgeDirection2[k];
		should_revert = true;
	}
	
	// ordering
	for (int i = 0; i < next_partner_ids.size(); i++) {
		
		double nextEdgeDirection[3];
		crossProduct(halfSpaceSet[base_id].ABCD, halfSpaceSet[next_partner_ids[i]].ABCD, nextEdgeDirection);
		normalize(nextEdgeDirection);
		double cosAngle = dotProduct(currentEdgeDirection, nextEdgeDirection);
		cosAngles.push_back(cosAngle);
		ordered_next_partner_ids.push_back(next_partner_ids[i]);

		for (int j = 0; j < cosAngles.size(); j++) {
			if (cosAngle < cosAngles[j]) {
				for (int k = cosAngles.size() - 1; k >= j; k--) {
					cosAngles[k] = cosAngles[k - 1];
					ordered_next_partner_ids[k] = ordered_next_partner_ids[k - 1];
				}
				cosAngles[j] = cosAngle;
				ordered_next_partner_ids[j] = next_partner_ids[i];
				break;
			}
		}
	}

	ordered_next_partner_ids.insert(ordered_next_partner_ids.begin(), base_id);
	ordered_next_partner_ids.push_back(partner_id);

	for (int i = 0, j = 1; j < ordered_next_partner_ids.size(); j++) {
		if (isWalkedEdge(ordered_next_partner_ids[i], ordered_next_partner_ids[j])) {
			double edgeDirection[3];
			if (i == 0) 
				findEdgeDirection(ordered_next_partner_ids[i], ordered_next_partner_ids[j], should_revert, NULL, edgeDirection);	
			else 
				findEdgeDirection(ordered_next_partner_ids[i], ordered_next_partner_ids[j], should_revert, halfSpaceSet[base_id].ABCD, edgeDirection);
			
			if (isValidEdge(startPoint, edgeDirection)) {
				edgePartners[ordered_next_partner_ids[i]].push(EdgePartnerTriple(ordered_next_partner_ids[j], edgeDirection, startPoint, kernel.getNumOfVerts() - 1));
				i = j;
			}
		}
		else
			i = j;
	}

	ordered_next_partner_ids.clear();
	cosAngles.clear();
}

bool KernelExpansion_KerTrack::isWalkedEdge(int hs1_id, int hs2_id) {

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

void KernelExpansion_KerTrack::findEdgeDirection(int hs1_id, int hs2_id, bool should_revert, double* directioner, double* edgeDirection) {

	crossProduct(halfSpaceSet[hs1_id].ABCD, halfSpaceSet[hs2_id].ABCD, edgeDirection);
	normalize(edgeDirection);
	if (directioner) {
		if (dotProduct(edgeDirection, directioner) > 0)
			multVect(edgeDirection, -1, edgeDirection);
	}
	else {
		if (should_revert)
			multVect(edgeDirection, -1.0, edgeDirection);
	}

}

bool KernelExpansion_KerTrack::isValidEdge(double* startPoint, double* edgeDirection) {

	double testPoint[3];
	for (int k = 0; k < 3; k++)
		testPoint[k] = startPoint[k] + EPSILON * edgeDirection[k];
	double t = findClosestValueSatisfiedByPoint(testPoint, halfSpaceSet);
	if (t > EPSILON)
		return false;
	return true;
	
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

	Vertex v = hostMeshptr->getVertex(0);
	double* point = new double[3]{ v.coords[0], v.coords[1], v.coords[2] };
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

	Vertex v = hostMeshptr->getVertex(0);
	for (int k = 0; k < 3; k++)
		point[k] = v.coords[k];

	for (int i = 0; i < hostMeshptr->getNumOfTris(); i++) {

		// construct the halfspace
		Triangle tri = hostMeshptr->getTriangle(i);
		HalfSpace halfSpace(hostMeshptr->getVertex(tri.corners[0]).coords, tri.normal, false);
		halfSpaceSet.push_back(halfSpace);

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
	
	Vertex v = hostMeshptr->getVertex(0);
	double* point = new double[3]{ v.coords[0], v.coords[1], v.coords[2] };
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
