// @author Merve Asiler

#include "BasicMeshElements.h"

void Vertex::setNormal(double x, double y, double z) {
	normal[0] = x;
	normal[1] = y;
	normal[2] = z;
};

Vertex::~Vertex() {
	delete[] coords;
	delete[] normal;
	vertList.clear();
	triList.clear();
	edgeList.clear();
};


Edge::~Edge() {
	delete[] endVerts;
	triList.clear();
};

void Edge::computeLength(Vertex* v1, Vertex* v2) {
	double* coords1 = v1->coords;
	double* coords2 = v2->coords;
	length = sqrt(pow(coords1[0] - coords2[0], 2) + pow(coords1[1] - coords2[1], 2) + pow(coords1[2] - coords2[2], 2));
};


Triangle::Triangle(const Triangle& t) {
	this->idx = t.idx;
	this->corners = t.corners;		// shallow copy
	for (int i = 0; i < t.edgeList.size(); i++)
		this->edgeList.push_back(t.edgeList[i]);
	for (int i = 0; i < t.triList.size(); i++)
		this->triList.push_back(t.triList[i]);
	if (t.normal != NULL)
		this->normal = t.normal;	// shallow copy
}

Triangle::~Triangle() {
	delete[] corners;
	if (center != NULL)
		delete[] center;
	if (normal != NULL)
		delete[] normal;
	if (areaVect != NULL)
		delete[] areaVect;
	edgeList.clear();
	triList.clear();
	angleList.clear();
};

void Triangle::computeCenter(Vertex* v1, Vertex* v2, Vertex* v3) {
	center = new double[3];
	for (int i = 0; i < 3; i++)
		center[i] = (v1->coords[i] + v2->coords[i] + v3->coords[i]) / 3;
}

void Triangle::computeNormal(Vertex* v1, Vertex* v2, Vertex* v3) {

	double edgeVect1[3], edgeVect2[3];
	for (int i = 0; i < 3; i++) {
		edgeVect1[i] = v2->coords[i] - v1->coords[i];
		edgeVect2[i] = v3->coords[i] - v2->coords[i];
	}

	normal = crossProduct(edgeVect1, edgeVect2);
	normalize(normal);
}

void Triangle::computeNormal(Vertex* v1, Vertex* v2, Vertex* v3, double* normal) {

	double edgeVect1[3], edgeVect2[3];
	for (int i = 0; i < 3; i++) {
		edgeVect1[i] = v2->coords[i] - v1->coords[i];
		edgeVect2[i] = v3->coords[i] - v2->coords[i];
	}

	crossProduct(edgeVect1, edgeVect2, normal);
	normalize(normal);
}

void Triangle::computeAreaVector(Vertex* v1, Vertex* v2, Vertex* v3) {

	double vect1[3], vect2[3];
	for (int i = 0; i < 3; i++)
		vect1[i] = v3->coords[i] - v2->coords[i];
	for (int i = 0; i < 3; i++)
		vect2[i] = v1->coords[i] - v2->coords[i];

	// normalize
	double length1 = sqrt(pow(vect1[0], 2) + pow(vect1[1], 2) + pow(vect1[2], 2));
	double length2 = sqrt(pow(vect2[0], 2) + pow(vect2[1], 2) + pow(vect2[2], 2));
	for (int i = 0; i < 3; i++) {
		vect1[i] = vect1[i] / length1;
		vect2[i] = vect2[i] / length2;
	}

	// cross product
	areaVect = new double[3];
	areaVect[0] = vect1[1] * vect2[2] - vect1[2] * vect2[1];
	areaVect[1] = vect1[2] * vect2[0] - vect1[0] * vect2[2];
	areaVect[2] = vect1[0] * vect2[1] - vect1[1] * vect2[0];

	// normalize
	double length = sqrt(pow(areaVect[0], 2) + pow(areaVect[1], 2) + pow(areaVect[2], 2));
	for (int i = 0; i < 3; i++)
		areaVect[i] = areaVect[i] / length;
}

void Triangle::setAngSkewness(double angSkewness) {
	angularSkewness = angSkewness;
}

double Triangle::getAngSkewness() {
	return angularSkewness;
}

void Triangle::setSquish(double squish) {
	this->squish = squish;
}

double Triangle::getSquish() {
	return squish;
}

TriangleWithVerts::~TriangleWithVerts() {
	for (int i = 0; i < 3; i++)
		vertices[i] = NULL;
	this->corners = NULL;
	this->normal = NULL;
	this->edgeList.clear();
	this->triList.clear();
};

double* TriangleWithVerts::computeBarycentricCoords(double point[3]) {	// assumes the triangle's normal is ready

	double* operand1 = diffVects(vertices[1]->coords, vertices[0]->coords);
	double* operand2 = diffVects(vertices[2]->coords, vertices[0]->coords);
	double* product = crossProduct(operand1, operand2);
	double denominator = dotProduct(product, normal);
	delete[] operand1;
	delete[] operand2;
	delete[] product;
	
	double* baryCoords = new double[3];
	for (int k = 0; k < 3; k++) {
		double* operand1 = diffVects(vertices[(2 + k) % 3]->coords, vertices[(1 + k) % 3]->coords);
		double* operand2 = diffVects(point, vertices[(1 + k) % 3]->coords);
		double* product = crossProduct(operand1, operand2);
		baryCoords[(2 + k) % 3] = dotProduct(product, normal) / denominator;
		delete[] operand1;
		delete[] operand2;
		delete[] product;
	}

	// check if this is a valid barycentric (the intersection point is on the triangle or not)
	if (abs((baryCoords[0] + baryCoords[1] + baryCoords[2]) - 1.0) < EPSILON) {
		for (int k = 0; k < 3; k++)
			if (baryCoords[k] < 0 || baryCoords[k] > 1.0) {
				delete[] baryCoords;
				return NULL;
			}

	}
	else {
		delete[] baryCoords;
		return NULL;
	}

	return baryCoords;
}
