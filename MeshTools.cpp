// @author Merve Asiler

#include "Mesh.h"
#include "BasicGeometricElements.h"
#include "BaseGeoOpUtils.h"

void Mesh::computeTrisAngles() {

	for (unsigned int i = 0; i < tris.size(); i++) {

		Triangle* triangle = tris[i];

		Edge* a = edges[triangle->edgeList[0]];
		Edge* b = edges[triangle->edgeList[1]];
		Edge* c = edges[triangle->edgeList[2]];

		// law of cosines
		double cos_a = (pow((b->length), 2) + pow((c->length), 2) - pow((a->length), 2)) / (2 * b->length * c->length);
		double cos_b = (pow((a->length), 2) + pow((c->length), 2) - pow((b->length), 2)) / (2 * a->length * c->length);
		double cos_c = (pow((a->length), 2) + pow((b->length), 2) - pow((c->length), 2)) / (2 * a->length * b->length);

		triangle->angleList.push_back(acos(cos_a));
		triangle->angleList.push_back(acos(cos_b));
		triangle->angleList.push_back(acos(cos_c));

	}

}

void Mesh::tetrahedralizeSurface() {

	// find each triangle's surface tetrahedron 
	for (int i = 0; i < tris.size(); i++) {
		
		//compute the peak point of the tetrahedron (the fourth point)
		Triangle* triangle = tris[i];
		double* peakCoords = new double[3];

		// get the pointers for the triangle corners
		Vertex* corners[3];
		for (int j = 0; j < 3; j++)
			corners[j] = verts[triangle->corners[j]];

		// find the center of the triangle
		triangle->computeCenter(corners[0], corners[1], corners[2]);

		// find the vector in the direction of the normal which will point to the peak point
		double edgeVector1[3], edgeVector2[3];
		for (int j = 0; j < 3; j++) {
			edgeVector1[j] = corners[1]->coords[j] - corners[0]->coords[j];
			edgeVector2[j] = corners[2]->coords[j] - corners[1]->coords[j];
		}
		double* dirVect = crossProduct(edgeVector1, edgeVector2);
		double scaleFactor = sqrt(computeLength(dirVect));
		for (int j = 0; j < 3; j++)
			dirVect[j] /= scaleFactor;

		// find the peak point
		for (int j = 0; j < 3; j++)
			peakCoords[j] = triangle->center[j] + dirVect[j];
		delete[] dirVect;

		Vertex* peak = new Vertex(i, peakCoords);
		Tetrahedron* tetrahedron = new Tetrahedron(i, peak);
		tets.push_back(tetrahedron);

	}

}

// Writes the tetrahedralized surface version of a mesh into a file 
void Mesh::writeSurfaceTetrahedralizedVersion(const char* name) {

	int nVerts = verts.size();
	int nTris = tris.size();

	ofstream offFile(name);
	if (offFile.is_open()) {

		// Initialize wrtiting, write "OFF", number of vertices and triangles to .off file
		offFile << "OFF" << endl;
		offFile << nVerts + nTris << " " << nTris * 4 << " " << "0" << endl;

		// Write the vertices into .off file
		for (int i = 0; i < nVerts; i++) {
			double* coords = verts[i]->coords;
			offFile << coords[0] << " " << coords[1] << " " << coords[2] << " " << endl;	// x, y, z coordinates of the i^{th} vertex
		}

		for (int i = 0; i < nTris; i++) {
			double* coords = tets[i]->peak->coords;
			offFile << coords[0] << " " << coords[1] << " " << coords[2] << " " << endl;
		}
		// Write the triangles into .off file
		for (int i = 0; i < nTris; i++) {
			int* corners = tris[i]->corners;
			offFile << "3" << " " << corners[0] << " " << corners[1] << " " << corners[2] << endl;
			offFile << "3" << " " << corners[0] << " " << corners[1] << " " << nVerts + i << endl;
			offFile << "3" << " " << corners[1] << " " << corners[2] << " " << nVerts + i << endl;
			offFile << "3" << " " << corners[2] << " " << corners[0] << " " << nVerts + i << endl;
		}

		offFile.close();
	}
}

double Mesh::computeVolume() {

	double volume = 0;
	for (int i = 0; i < tris.size(); i++) {

		Triangle* triangle = tris[i];
		Vertex* corners[3];

		for (int j = 0; j < 3; j++)
			corners[j] = verts[triangle->corners[j]];
		double* a_cross_b = crossProduct(corners[0]->coords, corners[1]->coords);

		double a_cross_b_dot_c = dotProduct(a_cross_b, corners[2]->coords);
		delete[] a_cross_b;

		volume += a_cross_b_dot_c / 6.0;
	}

	return abs(volume);
}

vector<double> Mesh::computeCurvaturePerVertex() {

	vector<double> curvatures;
	vector< vector <double> > aroundAngles;	// for each vertex it holds the angles around it
	aroundAngles.resize(verts.size());

	for (int i = 0; i < tris.size(); i++) {
		Triangle* t = tris[i];
		for (int c = 0; c < 3; c++) {
			Edge* e1 = edges[t->edgeList[c]];
			Edge* e2 = edges[t->edgeList[(c+1)%3]];
			double* vect1 = diffVects(verts[e1->endVerts[0]]->coords, verts[e1->endVerts[1]]->coords);
			double* vect2 = diffVects(verts[e2->endVerts[1]]->coords, verts[e2->endVerts[0]]->coords);
			normalize(vect1);
			normalize(vect2);
			if (e1->endVerts[1] == e2->endVerts[0])
				aroundAngles[e1->endVerts[1]].push_back( acos(dotProduct(vect1, vect2)) );
			else if (e1->endVerts[1] == e2->endVerts[1])
				aroundAngles[e1->endVerts[1]].push_back( PI - acos(dotProduct(vect1, vect2)) );
			else if (e1->endVerts[0] == e2->endVerts[1])
				aroundAngles[e1->endVerts[0]].push_back( acos(dotProduct(vect1, vect2)) );
			else // (e1->endVerts[0] == e2->endVerts[0])
				aroundAngles[e1->endVerts[0]].push_back( PI - acos(dotProduct(vect1, vect2)) );
		}
	}

	for (int i = 0; i < verts.size(); i++) {
		double angleSum = 0;
		for (int j = 0; j < aroundAngles[i].size(); j++)
			angleSum += aroundAngles[i][j];
		curvatures.push_back(2*PI - angleSum);
	}

	return curvatures;
}

int findDistinctEnd(int* c1, int* c2) {

	int i, j;
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++)
			if (c1[i] == c2[j])
				break;
		if (j == 3)
			return i;
	}
	return -1;
}

bool isEdgeConvex(Triangle* t1, Triangle* t2, vector<Vertex*> vertexSet) {

	// check whether they are on the same plane
	bool theSame = true;
	for (int i=0; i<3; i++)
		if (abs(t1->normal[i] - t2->normal[i]) > EPSILON) {
			theSame = false;
			break;
		}
	if (theSame)
		return true;
	/*
	theSame = true;
	for (int i = 0; i < 3; i++)
		if (abs(t1->normal[i] + t2->normal[i]) > EPSILON) {
			theSame = false;
			break;
		}
	if (theSame)
		return true;
		*/
	// find the distinct vertex on t1
	int p_id = findDistinctEnd(t1->corners, t2->corners);
	Vertex* p = vertexSet[t1->corners[p_id]];
	Line* t1_line;
	Plane* t2_plane;

	// check for perpendicular planes:
	double cosAngle = findCosAngleBetween(t1->normal, t2->normal);	// then the cos of angle between the planes is -cosAngle.
	if (abs(cosAngle) < EPSILON) {	// if they are perpendicular
		double plane_normal[3];
		for (int i = 0; i < 3; i++)
			plane_normal[i] = (t1->normal[i] + t2->normal[i]) / 2;
		normalize(plane_normal);
		Vertex* q = vertexSet[t1->corners[(p_id+1)%3]];
		t2_plane = new Plane(q->coords, plane_normal);
	}
	else							// if they are not perpendicular
		t2_plane = new Plane(vertexSet[t2->corners[0]]->coords, t2->normal);

	// find the intersection
	t1_line = new Line(p->coords, t1->normal);
	double parameter = findLinePlaneIntersection(t1_line, t2_plane);
	delete t1_line;
	delete t2_plane;

	if (abs(cosAngle) < EPSILON && parameter < 0)
		return false;
	if (abs(cosAngle) < EPSILON && parameter > 0)
		return true;
	if (cosAngle < 0 && parameter > 0) 
		return false;
	if (cosAngle > 0 && parameter < 0)
		return false;
	return true;
}

void Mesh::makeConvex() {

	bool need_to_check = true;
	int counter = 0;
	while (need_to_check) {
		if (counter == 0)
			break;
		counter++;
		need_to_check = false;
		for (int i = 0; i < edges.size(); i++) {
			Triangle* t1 = tris[edges[i]->triList[0]];
			Triangle* t2 = tris[edges[i]->triList[1]];

			t1->computeNormal(verts[t1->corners[0]], verts[t1->corners[1]], verts[t1->corners[2]]);
			t2->computeNormal(verts[t2->corners[0]], verts[t2->corners[1]], verts[t2->corners[2]]);

			if (isEdgeConvex(t1, t2, verts))
				continue;

			// concave
			need_to_check = true;	// after all the edges are handled, we will re-check in one more loop

			double n[3];
			for (int j = 0; j < 3; j++)
				n[j] = ((double)t1->normal[j] + t2->normal[j]) / 2.0;
			normalize(n);
			int ind = findDistinctEnd(t1->corners, t2->corners);
			Vertex* p = verts[t1->corners[ind]];
			Plane* plane = new Plane(p->coords, n);

			Vertex* v1 = verts[t1->corners[(ind + 1) % 3]];
			double* new_v1_coords = projectVertexOnPlane(plane, v1->coords);

			Vertex* v2 = verts[t1->corners[(ind + 2) % 3]];
			double* new_v2_coords = projectVertexOnPlane(plane, v2->coords);

			for (int j = 0; j < 3; j++) {
				v1->coords[j] = new_v1_coords[j];
				v2->coords[j] = new_v2_coords[j];
			}

			delete[] new_v1_coords;
			delete[] new_v2_coords;
			delete plane;
			
		}
	}
	
	for (int i = 0; i < tris.size(); i++) {
		delete[] tris[i]->normal;
		tris[i]->normal = NULL;
	}
}

void Mesh::removeAbnormalTris() {
	int initial_num_of_edges = edges.size();
	for (int i = 0, removed = 0; i < initial_num_of_edges; i++) {
		Edge* edge = edges[i - removed];
		if (edge->triList.size() == 1) {
			cout << "1-neighbor-triangle" << endl;
			edges.erase(edges.begin() + (i - removed));
			removed++;
		}
		else if (edge->triList.size() > 2)
			cout << "more than 2-neighbor-triangle" << endl;

	}
}
