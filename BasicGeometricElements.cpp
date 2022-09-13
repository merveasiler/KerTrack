// @author Merve Asiler

#include "BasicGeometricElements.h"

#pragma region LINE

/*
	@param A point on the line:
	@param The direction vector of the line
*/
Line::Line(double* point, double* directionVector) {
	this->point = new double[3];
	this->directionVector = new double[3];

	for (unsigned int i = 0; i < 3; i++) {
		this->point[i] = point[i];
		this->directionVector[i] = directionVector[i];
	}
}

Line::~Line() {
	delete[] point;
	delete[] directionVector;
	point = nullptr;
	directionVector = nullptr;
}

#pragma endregion

#pragma region HALFPLANE

HalfPlane::HalfPlane(double* point, double* directionVector, double* normalVector, bool doesCoverNormalSide) : Line(point, directionVector) {

	this->normalOnPlane = new double[3];
	for (int i = 0; i < 3; i++)
		this->normalOnPlane[i] = normalVector[i];
	this->doesCoverNormalSide = doesCoverNormalSide;
}

HalfPlane::~HalfPlane() {
	delete[] normalOnPlane;
	normalOnPlane = nullptr;
}

#pragma endregion

#pragma region PLANE

/*
	@param A point on the plane
	@param The normal vector of the plane
*/
Plane::Plane(double* point, double* normalVector) {

	this->ABCD = new double[4];
	this->point = new double[3];

	for (unsigned int i = 0; i < 3; i++) {
		this->ABCD[i] = normalVector[i];
		this->point[i] = point[i];
	}

	this->ABCD[3] = -dotProduct(normalVector, point);
}

Plane::~Plane() {
	delete[] ABCD;
	delete[] point;
	ABCD = nullptr;
	point = nullptr;
}

void Plane::setId(int idx) {
	
	this->idx = idx;
}

#pragma endregion

#pragma region HALFSPACE

HalfSpace::HalfSpace(double* point, double* normalVector, bool isLargerThanZero) : Plane(point, normalVector) {

	this->isLargerThanZero = isLargerThanZero;
}

HalfSpace::~HalfSpace() {

}

#pragma endregion