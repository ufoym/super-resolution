#pragma once
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include "image.h"
#include "helper.h"

//-------------------------------------------------------------------------
// helper: constant variable
//-------------------------------------------------------------------------

// suppose white background
const rgbf vBG = rgbf(255, 255, 255); 

// partition 3x3 window to 2 halves in four ways: - | / \
// for each way:
//     estimate min1/min2 and max1/max2 of either half window
//     set min=min1 max=max1 or min=min2 max=max2
//     compare local min/max and set global min/max
const int PAR_POS9[4][8] = { 
	{5,2,7,1,4,6,3,8}, 
	{5,1,6,2,3,7,4,8}, 
	{5,1,2,6,7,8,4,3}, 
	{6,1,3,5,8,7,2,4}
};

const int NEIB_POS9[9][2] = { 
	{ 0, 0}, {-1,-1}, { 0,-1}, 
	{ 1,-1}, { 1, 0}, { 1, 1}, 
	{ 0, 1}, {-1, 1}, {-1, 0}
};

const int NEIB_POS[] = {
	 0, 0,	-1, 0,	 0,-1,	 0, 1,	 1, 0,	-1,-1,	-1, 1,	 1,-1,	
	 1, 1,	-2, 0,	 0,-2,	 0, 2,	 2, 0,	-2,-1,	-2, 1,	-1,-2,	
	-1, 2,	 1,-2,	 1, 2,	 2,-1,	 2, 1,	-2,-2,	-2, 2,	 2,-2,	
	 2, 2,	-3, 0,	 0,-3,	 0, 3,	 3, 0,	-3,-1,	-3, 1,	-1,-3,	
	-1, 3,	 1,-3,	 1, 3,	 3,-1,	 3, 1,	-3,-2,	-3, 2,	-2,-3,	
	-2, 3,	 2,-3,	 2, 3,	 3,-2,	 3, 2,	-4, 0,	 0,-4,	 0, 4,	
	 4, 0,	-4,-1,	-4, 1,	-1,-4,	-1, 4,	 1,-4,	 1, 4,	 4,-1,	
	 4, 1,	-3,-3,	-3, 3,	 3,-3,	 3, 3,	-4,-2,	-4, 2,	-2,-4,	
	-2, 4,	 2,-4,	 2, 4,	 4,-2,	 4, 2,	-5, 0,	-4,-3,	-4, 3,	
	-3,-4,	-3, 4,	 0,-5,	 0, 5,	 3,-4,	 3, 4,	 4,-3,	 4, 3,	
	 5, 0,	-5,-1,	-5, 1,	-1,-5,	-1, 5,	 1,-5,	 1, 5,	 5,-1,	
	 5, 1,	-5,-2,	-5, 2,	-2,-5,	-2, 5,	 2,-5,	 2, 5,	 5,-2,	
	 5, 2,	-4,-4,	-4, 4,	 4,-4,	 4, 4,	-5,-3,	-5, 3,	-3,-5,	
	-3, 5,	 3,-5,	 3, 5,	 5,-3,	 5, 3,	-6, 0,	 0,-6,	 0, 6,	
	 6, 0,	-6,-1,	-6, 1,	-1,-6,	-1, 6,	 1,-6,	 1, 6,	 6,-1,	
	 6, 1,	-6,-2,	-6, 2,	-2,-6,	-2, 6,	 2,-6,	 2, 6,	 6,-2
};

//-------------------------------------------------------------------------
// helper: type conversion
//-------------------------------------------------------------------------

inline rgb rgbf2rgb(const rgbf& x)
{rgb v((int)x.r, (int)x.g, (int)x.b); return v;}
inline rgbf rgb2rgbf(const rgb& x) 
{rgbf v(x.r, x.g, x.b); return v;}

//-------------------------------------------------------------------------
// disjoint-set forests using union-by-rank and path compression (sort of)
//-------------------------------------------------------------------------

typedef struct 
{
	int rank;
	int p;
	int size;
	rgbf sum;
} uni_elt;

//-------------------------------------------------------------------------

class universe 
{
public:

	universe(image<rgb> *im);

	~universe();

	int find(int x);

	void join(int x, int y);

	int size(int x) const;
	int num_sets() const;
	rgb color(int x);

private:
	uni_elt *elts;
	int num;
};


//-------------------------------------------------------------------------
// helper of segmentation core
//-------------------------------------------------------------------------

struct edge 
{
	float w;
	int a, b;
	edge() {w=0;}
};

static inline bool operator<(const edge &a, const edge &b) 
{
	return a.w < b.w;
}

// dissimilarity measure between pixels
static inline float diff(image<rgb> *im, int x1, int y1, int x2, int y2) 
{
	return (rgb2rgbf(im->get(x1,y1)) - rgb2rgbf(im->get(x2,y2))).length();
}

static inline float diff(const rgbf & e1, const rgbf &e2)
{
	//float rmean = (e1.r + e2.r) / 2;
	//rgbf tmp = e1 - e2;
	//return sqrt((int((512+rmean)*tmp.r*tmp.r)>>8) + 4*tmp.g*tmp.g \
	//	+ (int((767-rmean)*tmp.b*tmp.b)>>8));
	return (e1 - e2).length();
}

//-------------------------------------------------------------------------

struct bd_info 
{
	int a;
	int b;
	int first_a;
	int first_b;
	float cost;
	int edge_num;
	bool operator<(const bd_info & v) {	return cost < v.cost; }
};

//-------------------------------------------------------------------------
// segmentation core
//-------------------------------------------------------------------------

class Segment
{
	// --------------------------------------------------------------------
private:
	std::map<int,rgb> color_map;
	image<int> *seg_map;

	// --------------------------------------------------------------------
public:
	Segment(image<rgb> *im);

	~Segment();;
	int count() const;
	std::map<int,rgb> get_color_map();  // added by chizhang for the ease of aliased rendering, 2013-3-12.

	image<rgb> *vis();

	image<rgb> *vis_pseudo() const;

	// --------------------------------------------------------------------

	void trace( std::vector<float>&				nodes,
				std::vector<std::vector<int>>&	indices,
				std::vector<bool>&				junction_map ) const;

	// --------------------------------------------------------------------

private:

	// ------------------------------------anti-aliasing artifact removal

	float _project_dist(const rgbf color, const rgbf vCenter, 
		const rgbf vDir, float &t) const;

	void _compute_line_model_array(image<rgb> *im, 
		const int x, const int y,
		rgbf *aNeib, const int nNeib, rgbf& vDir) const;

	float _linear_combinability(const rgbf &c1, const rgbf &c2, 
		const rgbf &c, float g_fVarTol, float &t) const;


	void _find_local_min_max(
		rgbf *aNeib, const rgbf vDir, const float g_fVarTol, 
		const int par, const int start, const int end,
		float &tmin, float &tmax, int &idmin, int &idmax) const;

	bool _select_extrema_with_spatial_constraint(
		rgbf *aNeib, const rgbf vDir, float g_fVarTol, 
		int &idmin, int &idmax) const;

	float _blending_alpha(image<rgb> *im, const int x, const int y,
		const int nNeib, rgbf &p_min, rgbf &p_max, float &prob,
		float g_fVarTol = 0.1f) const;

	image<rgb> *remove_blending(image<rgb> *im) const;

	// ---------------------------------------modified pff's segmentation

	universe *_segment_graph(image<rgb> *im, 
		int num_edges, edge *edges) const;

	universe *segment(image<rgb> *im) const;

	// -----------------------------------------------------merge segment

	float _merge_thres(const int size_a, const int size_b) const;

	void _color_merge( int num, universe * u, edge * edges ) const;

	void _boundary_merge( int num, universe * u, edge * edges, 
		image<rgb> * im ) const;

	// -----------------------------------------------------diagonal fuse

	rgbf _get_px(image<rgb> *im, const int x, const int y,
		const rgbf &replaceColor = rgbf(255, 255, 255)) const;

	double _kappa_t_px(image<rgb> *im,const int x,const int y) const;
	double _kappa_m_px(image<rgb> *im,const int x,const int y) const;
	void diagonal_fuse(image<rgb> *im, const double tau_k = 1.0);
};