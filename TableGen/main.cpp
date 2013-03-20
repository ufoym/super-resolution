#include <iostream>
#include <cstring>
#include <vector>
#include <cmath>
#include <ctime>
#include "../image.h"
#include "table.h"




int main()
{
	
	TableGenerator gen;
	gen.exec();

	gen.savePatchTableBinary("patch_table.dat");
	gen.savePolygonTableBinary("polygon_table.dat");


	int start = clock();

	gen.loadPatchTableBinary("patch_table.dat");
	gen.loadPolygonTableBinary("polygon_table.dat");

	int elapsed = (clock() - start) / 1000;
	std::cout << "elapsed time: " << elapsed << "s\n";

	gen.visualize();
	

	
	return 0;
}









//void generateTable()
//{
//	for (float d1 = -28 * 0.3; d1 <= 28 * 0.3; d1 += 0.3) {
//		for (float theta1 = 0.0; theta1 <= 175.0; theta1 += 5.0) {
//			for (float d2 = -28 * 0.3; d2 <= 28 * 0.3; d2 += 0.3) {
//				for (float theta2 = 0.0; theta2 <= 175.0; theta2 += 5.0) {
//
//					// case 1: the two lines coincide.
//					if (d1 == d2 && theta1 == theta2) {
//
//						continue;
//					}
//
//					// case 2: the two lines intersect.
//
//					bool verti1 = (std::abs(theta1 - 90) < 0.1);
//					bool verti2 = (std::abs(theta2 - 90) < 0.1);
//
//					double x1, y1, x2, y2, k1, k2, x, y,
//						t1 = theta1 / PI,
//						t2 = theta2 / PI;
//					std::vector<float> markers;
//
//					if (verti1 || verti2) {
//						if (verti1) {
//							std::swap(d1, d2);
//							std::swap(theta1, theta2);
//						}
//						x1 = d1 * std::sin(t1);
//						y1 = -d1 * std::cos(t1);
//						k1 = std::tan(t1);
//						x2 = d2;
//						x = x2;
//						y = y1 + k1 * (x - x1);
//					}
//					else {
//
//					}
//
//				}
//			}
//		}
//	}
//}
//
//bool isAccepted(float d1, float theta1, float d2, float theta2) {
//	bool vertical1 = (std::abs(theta1 - 90) < 0.1);
//	bool vertical2 = (std::abs(theta2 - 90) < 0.1);
//
//	if (vertical1 && vertical2) {
//		return std::abs(d1 - d2) < 0.1;
//	}
//
//	if (vertical1 || vertical2) {
//		if (vertical1) {
//			std::swap(d1, d2);
//			std::swap(theta1, theta2);
//		}
//		double x2 = d2;
//		double t1 = theta1 / PI;
//		double x1 = d1 * std::sin(t1);
//		double y1 = -d1 * std::cos(t1);
//		double k1 = std::tan(t1);
//		double x12 = x2;
//		double y12 = y1 + k1 * (x12 - x1);
//
//		return 0. <= x12 && x12 <= 6 
//			&& 0. <= y12 && y12 <= 6;
//	}
//
//	// now non of the two lines are vertical.
//	double t1 = theta1 / PI,
//		   t2 = theta2 / PI,
//		   x1 = d1 * std::sin(t1),
//		   x2 = d2 * std::sin(t2),
//		   y1 = -d1 * std::cos(t1),
//		   y2 = -d2 * std::cos(t2),
//		   k1 = std::tan(t1),
//		   k2 = std::tan(t2);
//
//	//if (d1 == d2 && theta1 == theta2)
//	//	return true;
//
//	double x12 = ((y1 - y2) + (k2*x2 - k1*x1)) / (k2 - k1);
//	double y12 = y1 + k1 * (x12 - x1);
//
//	return 0. <= x12 && x12 <= 6 
//		&& 0. <= y12 && y12 <= 6;
//}
//
//int old_main()
//{
//	int num_corners = 0;
//	for (float d1 = -28 * 0.3; d1 <= 28 * 0.3; d1 += 0.3) {
//		for (float theta1 = 0.0; theta1 <= 175.0; theta1 += 5.0) {
//			for (float d2 = -28 * 0.3; d2 <= 28 * 0.3; d2 += 0.3) {
//				for (float theta2 = 0.0; theta2 <= 175.0; theta2 += 5.0) {
//					if (isAccepted(d1, theta1, d2, theta2)) {
//						num_corners++;
//					}
//				}
//			}
//		}
//	}
//
//	std::cout << num_corners << std::endl;
//	return 0;
//}