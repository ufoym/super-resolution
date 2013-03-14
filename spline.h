#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <iostream>


struct BSpline
{
	int k;
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
#define NUM_NEIGHBORS		10
#define SAMPLE_INTERVAL		1
private:
	std::vector<BSpline> splines;

public:
	BSplineFitter(std::vector<float>&			nodes,
				std::vector<std::vector<int>>&	indices,
				std::vector<bool>&				junction_map);

	void saveToSVG(std::string& filename);

	void getSplines(std::vector<BSpline>& splines);

	void clear();

private:
	void _fitBSpline(std::vector<float>& points, BSpline& spline);

	void _sampleBSpline(BSpline& spline, std::vector<float>& samples);

	void _saveBSpline();

};