#pragma once
#ifndef LCOOR_H
#define LCOOR_H

#include "Point2D.h"
#include <math.h>
#include "Data.h"

class LCoordinate2D
{
public:
	double L1, L2, L3;
	Point2D *pointI, *pointJ, *pointK;
	double  *a, *b, *c;

	double **matrK;
	double *vecF;

	double detA = 0;

	LCoordinate2D(Point2D *pointI, Point2D *pointJ, Point2D *pointK) {
		this->pointI = new Point2D(*pointI);
		this->pointJ = new Point2D(*pointJ);
		this->pointK = new Point2D(*pointK);

		this->matrK = new double*[3];
		for (int i = 0; i < 3; i++)
		{
			this->matrK[i] = new double[3];
		}

		this->vecF = new double[3];

		a = new double[3];
		b = new double[3];
		c = new double[3];
	}

	void evalParams() {
		
		detA = 0.5 * ((pointJ->x * pointK->y - pointJ->y * pointK->x)
			- (pointI->x * pointK->y - pointI->y * pointK->x)
			+ (pointI->x * pointJ->y - pointI->y * pointJ->x));

		defineKoeff();
		
		defineMatrK();
		defineVecF();
	}

	double** getMatrK()
	{
		return matrK;
	}

	double* getVecF()
	{
		return vecF;
	}

	LCoordinate2D();
	~LCoordinate2D();

private:
	void defineKoeff();
	void defineMatrK();
	void defineVecF();
	double ** multiply(double ** matrix, double alpha);
	double integral(double(*func_k)(double , double));
	double k(int i, int j);


};

LCoordinate2D::LCoordinate2D()
{
	this->pointI = nullptr;
	this->pointJ = nullptr;
	this->pointK = nullptr;

}
LCoordinate2D::~LCoordinate2D()
{
}

void LCoordinate2D::defineKoeff() {
	a[0] = pointJ->x*pointK->y - pointK->x*pointJ->y;
	b[0] = pointJ->y - pointK->y;
	c[0] = pointK->x - pointJ->x;

	a[1] = pointK->x*pointI->y - pointI->x*pointK->y;
	b[1] = pointK->y - pointI->y;
	c[1] = pointI->x - pointK->y;

	a[2] = pointI->x*pointJ->y - pointJ->x*pointI->y;
	b[2] = pointI->y - pointJ->y;
	c[2] = pointJ->x - pointI->x;
/*
	L1 = 1 / (2 * detA)*(a[0] + b[0] * pointI->x + c[0] * pointI->y);
	L2 = 1 / (2 * detA)*(a[1] + b[1] * pointJ->x + c[1] * pointJ->y);
	L3 = 1 / (2 * detA)*(a[2] + b[2] * pointK->x + c[2] * pointK->y);
*/
	L1 = 1 / (2 * detA)*(a[0] + b[0] * (pointK->getX() + pointJ->getX()) / 2 + c[0] * (pointK->getY() + pointJ->getY()) / 2);
	L2 = 1 / (2 * detA)*(a[1] + b[1] * (pointK->getX() + pointI->getX()) / 2 + c[1] * (pointK->getY() + pointI->getY()) / 2);
	L3 = 1 / (2 * detA)*(a[2] + b[2] * (pointI->getX() + pointJ->getX()) / 2 + c[2] * (pointI->getY() + pointJ->getY()) / 2);
}

void LCoordinate2D::defineMatrK()
{
	for (int i = 0; i < COUNT_CONUS; i++)
 		for (int j = 0; j < COUNT_CONUS; j++)
			matrK[i][j] = k(i, j);
	
}

void LCoordinate2D::defineVecF()
{
	double integral_Q = (-1)*integral(func_Q);
	vecF[0] = L1*integral_Q;
	vecF[1] = L2*integral_Q;
	vecF[2] = L3*integral_Q;
}

//for matrix
double** LCoordinate2D::multiply(double** matrix, double alpha) {
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			matrix[i][j] = matrix[i][j] * alpha;
		}
	}
	return matrix;
}

double LCoordinate2D::integral(double (*func_k)(double, double ))
{
 	double a = func_k((pointI->x + pointK->x)/2, (pointI->y + pointK->y)/2);
	double b = func_k((pointI->x - pointJ->x)/2, (pointI->y - pointJ->y)/2);
	double c = func_k((pointJ->x - pointK->x)/2, (pointJ->y - pointK->y)/2);

	return (detA / 3)*(a + b + c);
}


double LCoordinate2D::k(int i, int j)
{
	return (b[i] * b[j] /(4*detA) * integral(func_kxx) + c[i]*c[j]/(4*detA)*integral(func_kyy));
}

#endif // !LCOOR_H