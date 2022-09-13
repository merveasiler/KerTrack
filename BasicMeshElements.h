// @author Merve Asiler

#pragma once

#include "BaseMathOpUtils.h"

struct Vertex
{
	int idx;
	double* coords = NULL;
	double* normal = NULL;

	// adjacencies
	vector< int > vertList;
	vector< int > triList;
	vector< int > edgeList;

	Vertex(int i, double* c) : idx(i), coords(c) { normal = new double[3]; };
	void setNormal(double x, double y, double z);
	~Vertex();
};

struct Edge
{
	int idx;
	int* endVerts = NULL;
	double length;

	// adjacency
	vector< int > triList;

	Edge(int i, int* c) : idx(i), endVerts(c) {};
	~Edge();
	void computeLength(Vertex* v1, Vertex* v2);
};

struct Triangle
{
	int idx;
	int* corners = NULL;
	double* center = NULL;	// center of gravity
	double* normal = NULL;
	double* areaVect = NULL; // unit area vector
	vector< int > edgeList;
	vector< int > triList; // neighbor tris
	vector< double > angleList;
	double angularSkewness;
	double squish;

	Triangle() {};
	Triangle(int i, int* c) : idx(i), corners(c) {};
	Triangle(const Triangle& t);	// shallow copy
	~Triangle();

	void computeCenter(Vertex* v1, Vertex* v2, Vertex* v3);

	void computeNormal(Vertex* v1, Vertex* v2, Vertex* v3);

	void computeNormal(Vertex* v1, Vertex* v2, Vertex* v3, double* normal);

	void computeAreaVector(Vertex* v1, Vertex* v2, Vertex* v3);

	void setAngSkewness(double angSkewness);

	double getAngSkewness();

	void setSquish(double squish);

	double getSquish();

};

struct TriangleWithVerts : public Triangle {
	Vertex* vertices[3];
	TriangleWithVerts(const Triangle& t, Vertex* v1, Vertex* v2, Vertex* v3) : Triangle(t) { 
		vertices[0] = v1; 
		vertices[1] = v2;
		vertices[2] = v3;
	}
	~TriangleWithVerts();
	double* computeBarycentricCoords(double point[3]);
};

struct Tetrahedron
{
	int idx;			// also the id of base triangle of the tetrahedron
	Vertex* peak;		// peak point of the tetrahedron
	double transformationDiff;

	Tetrahedron(int i, Vertex* p) : idx(i), peak(p) {};
	~Tetrahedron() {
		delete peak;
	};

	void setTransformationDiff(double transformationDiff) {
		this->transformationDiff = transformationDiff;
	}

	double getTransformationDiff() {
		return transformationDiff;
	}

};