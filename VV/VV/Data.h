#pragma once

double					// �������� ������� ������� ��������� ������ ���� L, R, B, T
L = -1,
R = 4,
B = -1,
T = 4;

#define COUNT_CONUS  3

double x_start = 0.0, y_start = 0.0, x_end = 4.0, y_end = 4.0;
double hx = 2.0, hy = 2.0;

bool needToDraw = true;

double rPoint = 0.2;
int minColorVal = 50, maxColorVal = 255;

double func_Q(double x, double y) {
	return 0;
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