// @author Merve Asiler

#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
#include "BasicMeshElements.h"

class Mesh
{

private:

	vector< Vertex* > verts;
	vector< Triangle* > tris;
	vector< Edge* > edges;
	vector< Tetrahedron* > tets;

public:
	Mesh() {};
	~Mesh();
	void loadObj(const char* name);
	void loadOff(const char* name);

	// methods to construct mesh
	void addTriangle(int v1, int v2, int v3);
	void addEdge(int v1, int v2);
	void addVertex(double x, double y, double z);
	int makeVertsNeighbor(int v1i, int v2i);

	// methods to reach mesh elements
	int getNumOfVerts() const { return verts.size(); };
	int getNumOfTris() const { return tris.size(); };
	int getNumOfEdges() const { return edges.size(); };
	int getNumOfTets() { return tets.size(); };
	Vertex* getVertex(int i) const { return verts[i]; };
	Triangle* getTriangle(int i) const { return tris[i]; };
	Edge* getEdge(int i) const { return edges[i]; };
	Tetrahedron* getTetrahedron(int i) { return tets[i]; };
	const vector<Vertex*>& getAllVerts() const { return verts; };
	const vector<Triangle*>& getAllTris() const { return tris; };
	const vector<Edge*>& getAllEdges() const { return edges; };
	const vector<Tetrahedron*>& getAllTetrahedra() const { return tets; };

	// methods to compute mesh features
	void computeTrisAngles();
	double computeVolume();
	void tetrahedralizeSurface();
	void writeSurfaceTetrahedralizedVersion(const char* name);
	vector<double> computeCurvaturePerVertex();
	void makeConvex();
	void removeAbnormalTris();
	void loadOffFromMatrices();
};

#pragma once