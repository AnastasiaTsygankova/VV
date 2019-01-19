#include "stdafx.h"
#include <iostream>
#include "Data.h"
#include "LCoordinate2D.h"
#include "Point2D.h"
#include "SolveSystemOfLinearEquation.h"
#include "rectangle.h"
#include "fstream"
#include "math.h"

int main()
{
	double  denom = 0.00001; //желаемая точность
	rectangle *rec = new rectangle(x_start, x_end, y_start, y_end, hx, hy);
	double *U = rec->calculate_temp();

	fstream fst;

	fst.open("out_U.txt", std::ios_base::out, std::ios_base::trunc);
	
	//for (int i = 0; i < rec->count_x; fst << "\n", i++)
	//	for (int j = 0; j < rec->count_y; j++)
	//		//printf("%d: %.2f  ", i + j*(rec->count_y), U[i + j*(rec->count_y)]);
	//		fst << /*i + j*(rec->count_y) << ": " <<*/ (U[i + j*(rec->count_y)] - remainder(U[i + j*(rec->count_y)], denom))<< "  ";
	//
	
	
	for (int i = 0; i < rec->count_points; i++)
		fst << (U[i] - remainder(U[i], denom)) << endl;

	
	fst.close();

	/*for (std::list<triangle2D>::iterator it = rec->list_of_triangles.begin(); it != rec->list_of_triangles.end(); it++)
	{
		printf("triangle \ni %f %f \nj %f %f \nk %f %f \n",
			(*it).points[0].getX(), (*it).points[0].getY(),
			(*it).points[1].getX(), (*it).points[1].getY(),
			(*it).points[2].getX(), (*it).points[2].getY());
	}
*/
	system("pause");
}