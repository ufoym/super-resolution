#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <iostream>


struct BSpline
{
	int k;
	float s;
	std::vector<float> t;
	std::vector<float> c;
	std::vector<float> u;
	int id_shape;
	int id_node;
	std::vector<float> data;
	std::vector<float> samples;
};

class BSplineFitter
{
private:
	#define NUM_CURVE_TYPES		3
	#define NUM_NEIGHBORS		3
	#define SAMPLE_INTERVAL		1
	std::vector<BSpline> splines;

public:
	BSplineFitter(std::vector<float>&			nodes,
				std::vector<std::vector<int>>&	indices,
				std::vector<bool>&				junction_map);

	void saveToSVG(std::string& filename);

	void getSplines(std::vector<BSpline>& splines);

	void clear();

private:
	void _fitBSpline(long k, float s, std::vector<float>& points, BSpline& spline);

	void _sampleBSpline(BSpline& spline, std::vector<float>& samples);

	void _saveBSpline();

};