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
	double edgeDirection[3];
	double startPoint[3];
	int startPointId;

	EdgePartnerTriple(int partnerId, double* edgeDirection, double* startPoint, int startPointId) {
		this->partnerId = partnerId;
		for (int i = 0; i < 3; i++) {
			this->edgeDirection[i] = edgeDirection[i];
			this->startPoint[i] = startPoint[i];
		}
		this->startPointId = startPointId;
	};

	~EdgePartnerTriple() {
	}

};

class KernelExpansion_KerTrack : public KernelExpansion {

	vector<queue<EdgePartnerTriple>> edgePartners;	// "id" of the other plane to construct and edge line with this & "direction" to walk on the line
	vector<vector<int>> edgePartnerIds;
	vector<vector<int>> vertexParentIds;
	vector<ProcessColor> isKernelFace;
	double _EPSILON = 3 * EPSILON;

	vector<double> initialize(double* point);
	int findTheClosestHalfSpace(double* point, vector<double>& scalarsVector);
	int findTheClosestHalfSpace(double* point, int id);
	vector<int> findTheClosestHalfSpace(int vertexId, double* lineDirection, int lineParent1Id, int lineParent2Id, double* newpoint);
	void orderTheFaces(int base_id, int partner_id, vector<int> next_partner_ids, double* startPoint, double* currentEdgeDirection);
	bool isWalkedEdge(int hs1_id, int hs2_id);
	void findEdgeDirection(int hs1_id, int hs2_id, bool should_revert, double* directioner, double* edgeDirection);
	bool isValidEdge(double* startPoint, double* edgeDirection);
	double* findInitialPoint_1();
	double* findInitialPoint_2();
	void findInitialPoint_3(double* point);
	double* findInitialPoint_4();
	void filterRepetitions(double* distances, vector<double>& scalarsVector);

public:
	KernelExpansion_KerTrack(const Mesh& hostMesh);
	~KernelExpansion_KerTrack();
	void expandKernel();
};