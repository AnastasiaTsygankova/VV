#pragma once
#ifndef POINT_H
#define POINT_H
class Point2D
{
public:
	double x, y;
	int num;

	Point2D();
	Point2D(double x, double y);
	Point2D(const Point2D &point);
	Point2D(double x, double y, int num);
	Point2D(const Point2D &point, int num);

	double getX() {
		return this->x;
	}
	double getY() {
		return this->y;
	}

	void setX(double x) {
		this->x = x;
	}
	void setY(double y) {
		this->y = y;
	}

	~Point2D();
};

Point2D::Point2D()
{
	this->x = 0;
	this->y = 0;
}
Point2D::Point2D(double x, double y) {
	this->x = x;
	this->y = y;
}
Point2D::Point2D(const Point2D &point) {
	this->x = point.x;
	this->y = point.y;
	this->num = point.num;
}
Point2D::Point2D(const Point2D &point, int num) {
	this->x = point.x;
	this->y = point.y;
	this->num = num;
}
Point2D::Point2D(double x, double y, int num)
{
	this->x = x;
	this->y = y;
	this->num = num;
}
Point2D::~Point2D()
{

}
#endif