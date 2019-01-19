#pragma once
#ifndef RECT_H
#define RECT_H

#include <fstream>
#include <list>
#include <algorithm>
#include "Point2D.h"
#include "LCoordinate2D.h"
#include "SolveSystemOfLinearEquation.h"



class triangle2D
{	
public:
	Point2D *points;

	bool thread[COUNT_CONUS];

	double **matr_k;
	double *f;

	triangle2D(double ix, double iy, double jx, double jy, double kx, double ky, bool i, bool j, bool k, int ni, int nj, int nk)
	{
		points = (Point2D*)std::malloc(3 * sizeof(Point2D));//new Point2D *[3];


		this->points[0] = Point2D(ix, iy, ni);
		this->points[1] = Point2D(jx, jy, nj);
		this->points[2] = Point2D(kx, ky, nk);

		thread[0] = i;
		thread[1] = j;
		thread[2] = k;
	}

	Point2D getPoint(int num)
	{
		return points[num];
	}

};


class rectangle
{
public:
	double x0, y0, x1, y1, hx, hy;
	int count_x, count_y;

	double **k;
	double *f;

	double *U;

	int count_points;

	std::list<triangle2D> *list_of_triangles;


	rectangle(double x0, double x1, double y0, double y1, double hx, double hy)
	{
		this->x0 = x0;
		this->x1 = x1;
		this->y0 = y0;
		this->y1 = y1;
		this->hx = hx;
		this->hy = hy;

		list_of_triangles = new list <triangle2D>;

		detect_greed();

		k = new double*[count_points];
		for (int i = 0; i < count_points; i++)
			k[i] = new double[count_points];

		f = new double[count_points];

		for (int i = 0; i < count_points; i++)
		{
			f[i] = 0;
			for (int j = 0; j < count_points; j++)
				k[i][j] = 0;
		}

	}


	void
	detect_greed()
	{
		double x = x0, y = y1;
	
		count_x = count_points_(x0, x1, hx);
		count_y = count_points_(y0, y1, hy);

		for(int i_x = 0; i_x < (count_x - 1); i_x++)
		{
			for(int i_y = 0; i_y < (count_y - 1);)
			{
				list_of_triangles->push_back
				(*new triangle2D(x, y, x, max(y - hy, y0), min(x + hx, x1), y,
					true, false, false, 
					i_x * count_y + i_y, i_x * count_y + (i_y + 1), (i_x + 1) * count_y + i_y));

				y = max(y - hy, y0);
				i_y++;
				list_of_triangles->push_back
				(*new triangle2D(x, y, min(x + hx, x1), y, min(x + hx, x1), y + hy,
					true, false, true,
					i_x * count_y + i_y, (i_x + 1) * count_y + i_y, (i_x + 1) * count_y + (i_y - 1)));
			}
			x += hx;
			y = y1;
		}

		this->count_points = count_x * count_y;
	}
	
	int	
	count_points_(double z0, double z1, double hz)
	{
		int cnt = (int)((z1 - z0) / hz);
		cnt += ((z0 + cnt * hz) == z1) ? 1 : 2;
		return cnt;
	}
	
	void
	formation_triangles()
	{
		fstream f_;

		f_.open("out_triangles.txt", std::ios_base::out, std::ios_base::trunc);
		
		f_ << "count triangles = " << list_of_triangles->size() << endl << endl;

		for (std::list<triangle2D>::iterator it = list_of_triangles->begin(); it != list_of_triangles->end(); it++)
		{
			/*printf("triangle \npoint %d %f %f \npoint %d  %f %f \npoint %d  %f %f \n",
				(*it).points[0].num, (*it).points[0].getX(), (*it).points[0].getY(),
				(*it).points[1].num, (*it).points[1].getX(), (*it).points[1].getY(),
				(*it).points[2].num, (*it).points[2].getX(), (*it).points[2].getY());
				*/
		
			LCoordinate2D *triangle =
 				new LCoordinate2D(&(*it).points[0], &(*it).points[1], &(*it).points[2]);

			triangle->evalParams();

			(*it).matr_k = triangle->getMatrK();
			(*it).f = triangle->vecF;
			
			for (int i = 0; i < 3; i++)
				f_ << "point " << (*it).points[i].num << " x " << (*it).points[i].getX() << " y " << (*it).points[i].getY() << endl;
			f_ << endl;
							
		}
		f_.close();	

	}

	double *
	calculate_temp()
	{
		formation_triangles();

		for (std::list<triangle2D>::iterator tr = list_of_triangles->begin(); tr != list_of_triangles->end(); tr++)
			compare_triangle(*tr);

		for (int i = 0; i < count_points; i++, printf("  =  %2.f  \n", f[i]))
			for (int j = 0; j < count_points; j++)
				printf("%2.f ", k[i][j]);
		
		border_conditions();
		
		for (int i = 0; i < count_points; i++, printf("  =  %2.f  \n", f[i]))
			for (int j = 0; j < count_points; j++)
				printf("%2.f ", k[i][j]);

		U = gauss(k, f, count_points);

		return U;
	}


	void
	compare_triangle(triangle2D tr)
	{
		for (int i = 0; i < count_points; i++)
			for (int c1 = 0; c1 < 3; c1++)
				if (i == tr.points[c1].num)
					for (int j = 0; j < count_points; j++)
						for (int c2 = 0; c2 < 3; c2++)
							if (j == tr.points[c2].num)
							{
								k[i][j] += tr.matr_k[c1][c2];
								f[i] = tr.f[c1];//??
							}
	}

	void
	sub_from_matrix(double T, int stb)
	{
		double _k = k[stb][stb];
		printf("\n k\n");

		for (int i = 0; i < count_points; i++)
		{
			k[stb][i] /= _k;
			printf("%2.f ", k[stb][i]);
		}
		f[stb] /= _k;

		printf("\n f\n");

		for (int i = 0; i < count_points; k[i][stb] = 0, i++)
		{
			f[i] -= k[i][stb] * T;
			printf("%2.f ", f[i]);

		}

		k[stb][stb] = 1;
	}

	void
	border_conditions()
	{
		for (int i = 0; i < count_y; i++)
		{
			sub_from_matrix(func_phi(x0, std::max(y0, y1 - hy * i)), i);
			sub_from_matrix(func_phi(x1, std::max(y0, y1 - hy * i)), i + ((count_y) * (count_x - 1)));

		}


		for (int j = 1; j < count_x - 1; j++)
		{
			sub_from_matrix(func_phi(std::min(x1, x0 + j*hx), y0), j * (count_y));
			sub_from_matrix(func_phi(x1, std::max(y0, y1 - hy * j)), j * (count_y) + (count_y -1));

		}

	}

};

#endif // !RECT_H
