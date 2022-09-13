// @author Merve Asiler

#pragma once

#include "KernelExpansion.h"

enum class ProcessColor {
	RED = 2,
	WHITE = 4,
	GREEN = 6
};

struct EdgePartnerTriple {

	int partnerId;
	double* edgeDirection;
	double* startPoint;
	int startPointId;

	EdgePartnerTriple(int partnerId, double* edgeDirection, double* startPoint, int startPointId) {
		this->partnerId = partnerId;
		this->edgeDirection = new double[3];
		this->startPoint = new double[3];
		for (int i = 0; i < 3; i++) {
			this->edgeDirection[i] = edgeDirection[i];
			this->startPoint[i] = startPoint[i];
		}
		this->startPointId = startPointId;
	};

	~EdgePartnerTriple() {
		delete[] edgeDirection;
		delete[] startPoint;
	}

};

class KernelExpansion_KerTrack : public KernelExpansion {

	vector<queue<EdgePartnerTriple*>> edgePartners;	// "id" of the other plane to construct and edge line with this & "direction" to walk on the line
	vector<vector<int>> edgePartnerIds;
	vector<vector<int>> vertexParentIds;
	vector<ProcessColor> isKernelFace;
	double _EPSILON = 3 * EPSILON;

	tuple<int, double*> findTheClosestHalfSpace(double point[3]);
	tuple<int, double*> findTheClosestHalfSpace(double point[3], int id);
	tuple<vector<int>, double*> findTheClosestHalfSpace(int vertexId, double* lineDirection, int lineParent1Id, int lineParent2Id);
	void orderTheFaces(int base_id, int partner_id, vector<int> next_partner_ids, double* startPoint, double* currentEdgeDirection);
	bool shouldSaveEdge(int hs1_id, int hs2_id);
	void saveFoundEdgeToProcess(int base_id, int next_partner_id, double* startPoint, double* edgeDirectioner);
	double* findInitialPoint_1();
	double* findInitialPoint_2();
	double* findInitialPoint_3();
	double* findInitialPoint_4();

public:
	KernelExpansion_KerTrack(const Mesh& hostMesh);
	~KernelExpansion_KerTrack();
	void expandKernel();
};