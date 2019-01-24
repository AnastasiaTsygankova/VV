#pragma once

double					// »сходные значени¤ мировых координат границ окна L, R, B, T
L = -1,
R = 4,
B = -1,
T = 4;

#define COUNT_CONUS  3

double x_start = 0.0, y_start = 0.0, x_end = 10.0, y_end = 10.0;
double hx = 1, hy = 1;

bool needToDraw = true;

double rPoint = 0.2;
int minColorVal = 50, maxColorVal = 255;

double func_Q(double x, double y) {
	return x;
}

double func_phi(double x, double y) {
	return 1;
}
double
func_kxx(double x, double y)
{
	return 1;
}

double
func_kyy(double x, double y)
{
	return 1;
}