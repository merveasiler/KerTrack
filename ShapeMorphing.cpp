// @author Merve Asiler

#include "ShapeMorphing.h"

void morphByKernel(Mesh* sourceMesh, Mesh* targetMesh, vector<Mesh*>& interMeshes, double centerSource[3], double centerTarget[3]) {

	double time_step = 1.0 / (interMeshes.size() + 1);
	double refDir[3] = { 0, 0, 1 };

	// define source mesh vertices w.r.t. center of the kernel
	vector <double*> directions_distances_source;
	for (int i = 0; i < sourceMesh->getNumOfVerts(); i++) {
		Vertex v = sourceMesh->getVertex(i);
		double* rayDirection = new double[4];
		for (int j = 0; j < 3; j++)
			rayDirection[j] = v.coords[j] - centerSource[j];
		rayDirection[3] = computeLength(rayDirection);	// only for the first 3 elements
		normalize(rayDirection);						// only for the first 3 elements
		directions_distances_source.push_back(rayDirection);
	}

	// define target mesh vertices w.r.t. center of the kernel
	vector <double*> directions_distances_target;
	for (int i = 0; i < targetMesh->getNumOfVerts(); i++) {
		Vertex v = targetMesh->getVertex(i);
		double* rayDirection = new double[4];
		for (int j = 0; j < 3; j++)
			rayDirection[j] = v.coords[j] - centerTarget[j];
		rayDirection[3] = computeLength(rayDirection);	// only for the first 3 elements
		normalize(rayDirection);						// only for the first 3 elements
		directions_distances_target.push_back(rayDirection);
	}

	// interpolate
	vector<double> angleDiffs, distanceDiffs;
	for (int i = 0; i < sourceMesh->getNumOfVerts(); i++) {
		// angle
		double* sourceDir = directions_distances_source[i];
		double* targetDir = directions_distances_target[i];


		double cosAngle = dotProduct(sourceDir, targetDir);
		double angle = acos(cosAngle), totalAngleDiff;
		double* axis1 = crossProduct(refDir, sourceDir);
		double* axis2 = crossProduct(refDir, targetDir);
		normalize(axis1);
		normalize(axis2);
		double s = dotProduct(axis1, axis2);
		if (s < 0)
			totalAngleDiff = angle;
		else 
			totalAngleDiff = 2 * PI - angle;

		double angleDiff = time_step * totalAngleDiff;
		angleDiffs.push_back(angleDiff);

		double totalDistanceDiff = directions_distances_target[i][3] - directions_distances_source[i][3];
		double distanceDiff = time_step * totalDistanceDiff;
		distanceDiffs.push_back(distanceDiff);
	}

	double translation[3], centerDiff[3];
	for (int i = 0; i < 3; i++) {
		centerDiff[i] = centerTarget[i] - centerSource[i];
		centerDiff[i] *= time_step;
	}

	// construct intermediate meshes
	for (int m = 0; m < interMeshes.size(); m++) {
		Mesh* mesh = interMeshes[m];
		double main_time_step = m + 1;
		for (int v = 0; v < sourceMesh->getNumOfVerts(); v++) {

			double* sourceDir = directions_distances_source[v];
			double* targetDir = directions_distances_target[v];

			double* axis = crossProduct(sourceDir, targetDir);
			double s = computeLength(axis);
			double angle = main_time_step * angleDiffs[v];
			double Vx[3][3] = { {0, -axis[2], axis[1]},
								{axis[2], 0, -axis[0]},
								{-axis[1], axis[0], 0} };

			double toplam[3][3] = { {1, -axis[2], axis[1]},
								{axis[2], 1, -axis[0]},
								{-axis[1], axis[0], 1} };

			double Vx2[3][3] = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };

			for (int k = 0; k < 3; k++)
				for (int i = 0; i < 3; i++)
					for (int j = 0; j < 3; j++)
						Vx2[k][j] += Vx[k][i] * Vx[i][j];

			double cosAngle = cos(angle);
			double coeff = (1.0 - cosAngle) / (s * s);
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					Vx2[i][j] *= coeff;

			double R[3][3];
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					R[i][j] = toplam[i][j] + Vx2[i][j];

			double dir[3] = { 0, 0, 0 };
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					dir[i] += R[i][j] * sourceDir[j];

			normalize(dir);
			double distance = sourceDir[3] + main_time_step * distanceDiffs[v];

			double coords[3];
			for (int i = 0; i < 3; i++)
				coords[i] = centerDiff[i] * main_time_step + centerSource[i] + dir[i] * distance;

			mesh->addVertex(coords[0], coords[1], coords[2]);
		}

		for (int t = 0; t < sourceMesh->getNumOfVerts(); t++) {
			Triangle triangle = sourceMesh->getTriangle(t);
			mesh->addTriangle(triangle.corners[0], triangle.corners[1], triangle.corners[2]);
		}

	}

}



/*
void morphByKernel(Mesh* sourceMesh, Mesh* targetMesh, vector<Mesh*>& interMeshes, double centerSource[3], double centerTarget[3]) {

	double time_step = 1.0 / (interMeshes.size()+1);
	double refDir[3] = { 0, 0, 1 };
	
	// define source mesh vertices w.r.t. center of the kernel
	vector <double*> directions_distances_source;
	for (int i = 0; i < sourceMesh->getNumOfVerts(); i++) {
		Vertex v = sourceMesh->getVertex(i);
		double* rayDirection = new double[4];
		for (int j = 0; j < 3; j++)
			rayDirection[j] = v.coords[j] - centerSource[j];
		rayDirection[3] = computeLength(rayDirection);	// only for the first 3 elements
		normalize(rayDirection);						// only for the first 3 elements
		directions_distances_source.push_back(rayDirection);
	}

	// define target mesh vertices w.r.t. center of the kernel
	vector <double*> directions_distances_target;
	for (int i = 0; i < targetMesh->getNumOfVerts(); i++) {
		Vertex v = targetMesh->getVertex(i);
		double* rayDirection = new double[4];
		for (int j = 0; j < 3; j++)
			rayDirection[j] = v.coords[j] - centerTarget[j];
		rayDirection[3] = computeLength(rayDirection);	// only for the first 3 elements
		normalize(rayDirection);						// only for the first 3 elements
		directions_distances_target.push_back(rayDirection);
	}

	// interpolate
	vector<double> angleDiffs, distanceDiffs;
	for (int i = 0; i < sourceMesh->getNumOfVerts(); i++) {
		// angle
		double* sourceDir = directions_distances_source[i];
		double* targetDir = directions_distances_target[i];
		
		double cosAngle = dotProduct(sourceDir, targetDir);
		double angle = acos(cosAngle), totalAngleDiff;
		double* cross = crossProduct(sourceDir, targetDir);
		double sinAngle = computeLength(cross);
		if (sinAngle > 0)
			totalAngleDiff = angle;
		else if (cosAngle > 0)
			totalAngleDiff = 2 * PI - angle;
		else
			totalAngleDiff = (PI - angle) + PI;
		
		double angleDiff = time_step * totalAngleDiff;
		angleDiffs.push_back(angleDiff);

		double totalDistanceDiff = directions_distances_target[i][3] - directions_distances_source[i][3];
		double distanceDiff = time_step * totalDistanceDiff;
		distanceDiffs.push_back(distanceDiff);
	}

	double translation[3], centerDiff[3];
	for (int i = 0; i < 3; i++) {
		centerDiff[i] = centerTarget[i] - centerSource[i];
		centerDiff[i] *= time_step;
	}

	// construct intermediate meshes
	for (int m = 0; m < interMeshes.size(); m++) {
		Mesh* mesh = interMeshes[m];
		double main_time_step = m + 1;
		for (int v = 0; v < sourceMesh->getNumOfVerts(); v++) {
		
			double* sourceDir = directions_distances_source[v];
			double* targetDir = directions_distances_target[v];

			double* axis = crossProduct(sourceDir, targetDir);
			double s = computeLength(axis);
			double angle = main_time_step * angleDiffs[v];
			double Vx[3][3] = { {0, -axis[2], axis[1]},
								{axis[2], 0, -axis[0]},
								{-axis[1], axis[0], 0} };

			double toplam[3][3] = { {1, -axis[2], axis[1]},
								{axis[2], 1, -axis[0]},
								{-axis[1], axis[0], 1} };

			double Vx2[3][3] = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };

			for (int k = 0; k < 3; k++)
				for (int i = 0; i < 3; i++)
					for (int j = 0; j < 3; j++)
						Vx2[k][j] += Vx[k][i] * Vx[i][j];

			double cosAngle = cos(angle);
			double coeff = (1.0 - cosAngle) / (s * s);
			for (int i = 0; i < 3; i++) 
				for (int j = 0; j < 3; j++)
					Vx2[i][j] *= coeff;

			double R[3][3];
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					R[i][j] = toplam[i][j] + Vx2[i][j];

			double dir[3] = { 0, 0, 0 };
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					dir[i] += R[i][j] * sourceDir[j];

			normalize(dir);
			double distance = sourceDir[3] + main_time_step * distanceDiffs[v];

			double coords[3];
			for (int i = 0; i < 3; i++)
				coords[i] = centerDiff[i] * main_time_step + centerSource[i] + dir[i] * distance;

			mesh->addVertex(coords[0], coords[1], coords[2]);
		}

		for (int t = 0; t < sourceMesh->getNumOfVerts(); t++) {
			Triangle triangle = sourceMesh->getTriangle(t);
			mesh->addTriangle(triangle.corners[0], triangle.corners[1], triangle.corners[2]);
		}

	}

}
*/

