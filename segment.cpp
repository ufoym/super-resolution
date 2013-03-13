#include "segment.h"

void Segment::diagonal_fuse( image<rgb> *im, const double tau_k /*= 1.0*/ )
{
	int w = seg_map->width(), h = seg_map->height();
	for (int y = 0; y < h-1; y++) {
		int *row_seg_map = seg_map->access[y]; 
		for (int x = 0; x < w-1; x++) {
			if(! ( row_seg_map[x] != seg_map->get(x, y+1)
				&& row_seg_map[x] != seg_map->get(x+1, y)
				&& seg_map->get(x+1, y+1) != seg_map->get(x, y+1)
				&& seg_map->get(x+1, y+1) != seg_map->get(x+1, y)))
				continue;
			//Compute psi_f(i,j,P) according to Eq.4.44
			double psi_f = 0; {                 
				psi_f += _kappa_m_px(im, x,  y);    //kappa^m_{i,j}                 
				psi_f += _kappa_m_px(im, x+1,y+1);  //kappa^m_{i+1,j+1}                 
				psi_f -= _kappa_t_px(im, x,  y+1);  //kappa^t_{i+1,j}                   
				psi_f -= _kappa_t_px(im, x+1,y);    //kappa^t_{i,j+1}
			}
			if(abs(psi_f) >= tau_k) {
				//Compute (delta_i, delta_j) according to Eq.4.46
				int delta_y=0, delta_x=0; {
					if(psi_f > tau_k) {
						delta_y = 1;
						delta_x = -1;
					}
					else if(psi_f < -tau_k) {
						delta_y = -1;
						delta_x = -1;
					}
				}

				//Compute omega_0, omega_1 according to Eq.4.45
				bool omega_0 = \
					(row_seg_map[x] == seg_map->get((x+delta_x+w)%w, y)
					&&row_seg_map[x] == seg_map->get(x,(y+delta_y+h)%h)
					&&row_seg_map[x] != seg_map->get((x+delta_x+w)%w, 
					(y+delta_y+h)%h));
				bool omega_1 = \
					(row_seg_map[x] == seg_map->get((x-delta_x+w)%w, y)
					&&row_seg_map[x] == seg_map->get(x,(y-delta_y+h)%h)
					&&row_seg_map[x] != seg_map->get((x-delta_x+w)%w, 
					(y-delta_y+h)%h));

				//Compute (y_0~3, x_0~3) according to Eq.4.47
				int y_0=0, y_1=0, y_2=0, y_3=0, 
					x_0=0, x_1=0, x_2=0, x_3=0; {
						if(psi_f > tau_k) {//ij_2 <-> ij_3 to keep symmetry 
							//(different from thesis, but similar to VM)
							y_0 = y;    x_0 = x;
							y_1 = y+1;  x_1 = x+1;
							y_2 = y;    x_2 = x+1;  
							y_3 = y+1;  x_3 = x;
						}
						else if(psi_f < -tau_k) {
							y_0 = y+1;  x_0 = x;
							y_1 = y;    x_1 = x+1;
							y_2 = y;    x_2 = x;
							y_3 = y+1;  x_3 = x+1;
						}
				}
				//Compute psi_y according to Eq.4.48
				double psi_y = 0; {
					rgbf tmp;
					tmp = _get_px(im,x_2,y_2) - _get_px(im,x_0,y_0);
					psi_y += tmp.dot(tmp);
					tmp = _get_px(im,x_2,y_2) - _get_px(im,x_1,y_1);
					psi_y += tmp.dot(tmp);
					tmp = _get_px(im,x_3,y_3) - _get_px(im,x_0,y_0);
					psi_y -= tmp.dot(tmp);
					tmp = _get_px(im,x_3,y_3) - _get_px(im,x_1,y_1);
					psi_y -= tmp.dot(tmp);
				}

				if(psi_y < 0) {
					if(omega_0 == false)
						seg_map->set(x_2,y_2, seg_map->get(x_0,y_0));
					else if(omega_1 == false)
						seg_map->set(x_3,y_3, seg_map->get(x_0,y_0));
				}
				else {
					if(omega_1 == false)
						seg_map->set(x_3,y_3, seg_map->get(x_0,y_0));
					else if(omega_0 == false)
						seg_map->set(x_2,y_2, seg_map->get(x_0,y_0));
				}
			} //end if
		} //end for_j
	} //end for_i
}

double Segment::_kappa_m_px( image<rgb> *im,const int x,const int y ) const
{
	rgbf tmp = _get_px(im, x-1, y+1) \
		- _get_px(im, x, y)*2 + _get_px(im, x+1, y-1);
	return tmp.dot(tmp);
}

double Segment::_kappa_t_px( image<rgb> *im,const int x,const int y ) const
{
	rgbf tmp = _get_px(im, x-1, y-1) \
		- _get_px(im, x, y)*2 + _get_px(im, x+1, y+1);
	return tmp.dot(tmp);
}

rgbf Segment::_get_px( image<rgb> *im, const int x, const int y, const rgbf &replaceColor /*= rgbf(255, 255, 255)*/ ) const
{
	return im->check(x, y) ? rgb2rgbf(im->get(x, y)) : replaceColor;
}

void Segment::_boundary_merge( int num, universe * u, edge * edges, image<rgb> * im ) const
{
	int width = im->width();
	std::vector<bd_info> co_boundary_info;
	std::map<std::pair<int,int>, int> co_boundary_map;
	for (int i = 0; i < num; i++) {
		int a = u->find(edges[i].a);
		int b = u->find(edges[i].b);
		if (a != b) {
			float color_dist = diff(
				rgb2rgbf(im->get(edges[i].a%width, edges[i].a/width)),
				rgb2rgbf(im->get(edges[i].b%width, edges[i].b/width))
				);
			std::pair<int,int> p(a,b);
			std::map<std::pair<int,int>, int>::const_iterator found \
				= co_boundary_map.find(p);
			if (found == co_boundary_map.end()) {
				bd_info bd = {a,b,edges[i].a,edges[i].b,color_dist,1};
				co_boundary_map[p] = (int)co_boundary_info.size();
				co_boundary_info.push_back(bd);
			}
			else {
				int idx = found->second;
				co_boundary_info[idx].cost += color_dist;
				co_boundary_info[idx].edge_num += 1;
			}
		}
	}
	sort(co_boundary_info.begin(), co_boundary_info.end());
	std::vector<bd_info>::const_iterator info_iter \
		= co_boundary_info.begin();
	for (; info_iter != co_boundary_info.end(); info_iter++) {
		int a = u->find(info_iter->first_a);
		int b = u->find(info_iter->first_b);
		if (a!=b) {
			int size_a = u->size(a);
			int size_b = u->size(b);
			if (info_iter->cost/info_iter->edge_num \
				< 0.5*_merge_thres(size_a,size_b))
				u->join(a, b);
		}
	}
}

void Segment::_color_merge( int num, universe * u, edge * edges ) const
{
	for (int i = 0; i < num; i++) {
		int a = u->find(edges[i].a);
		int b = u->find(edges[i].b);
		if (a != b) {
			int size_a = u->size(a);
			int size_b = u->size(b);
			float color_dist = diff(rgb2rgbf(u->color(a)), 
				rgb2rgbf(u->color(b)));
			if ( color_dist < _merge_thres(size_a, size_b) )
				u->join(a, b);
		}
	}
}

float Segment::_merge_thres( const int size_a, const int size_b ) const
{
	int min_size = std::min(size_a,size_b);
	return std::max(200.f/min_size, std::max(60.f-0.1f*min_size,15.f));
}

universe * Segment::segment( image<rgb> *im ) const
{
	int width = im->width();
	int height = im->height();

	// build graph
	edge *edges = new edge[width*height*4];
	int num = 0;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			for(int i = 3; i <= 6; ++i) {
				int ox = NEIB_POS9[i][0], oy = NEIB_POS9[i][1];
				if(!im->check(x+ox, y+oy))
					continue;
				edges[num].a = y * width + x;
				edges[num].b = (y+oy) * width + (x+ox);
				edges[num].w = diff(im, x, y, x+ox, y+oy);
				num++;
			}
		}
	}
	// segment
	universe *u = _segment_graph(im, num, edges);

	// post process: merge components
	_color_merge(num, u, edges);			// similar segment color
	_boundary_merge(num, u, edges, im);		// low co-boundary cost

	delete [] edges;
	return u;
}

universe * Segment::_segment_graph( image<rgb> *im, int num_edges, edge *edges ) const
{
	// sort edges by weight
	std::sort(edges, edges + num_edges);

	// make a disjoint-set forest
	universe *u = new universe(im);

	// for each edge, in non-decreasing weight order...
	for (int i = 0; i < num_edges; i++) {
		edge *pedge = &edges[i];

		// components connected by this edge
		int a = u->find(pedge->a);
		int b = u->find(pedge->b);
		if (a != b && pedge->w == 0) 
			u->join(a, b);
	}

	// free up
	return u;
}

image<rgb> * Segment::remove_blending( image<rgb> *im ) const
{
	int width = im->width();
	int height = im->height();
	image<rgb> *output = new image<rgb>(width, height);
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			rgbf p_min, p_max;
			float prob = 0;
			float alpha = _blending_alpha(im, x,y, 9,p_min,p_max,prob);
			output->set(x, y,
				(alpha < 0 || alpha > 1 || prob < 0.8) ? im->get(x, y) 
				: ( (alpha < 0.5) ? rgbf2rgb(p_min) : rgbf2rgb(p_max) )
				);
		}
	}  
	return output;
}

float Segment::_blending_alpha( image<rgb> *im, const int x, const int y, const int nNeib, rgbf &p_min, rgbf &p_max, float &prob, float g_fVarTol /*= 0.1f*/ ) const
{
	rgbf vDir;
	rgbf *vNeib = new rgbf[nNeib];
	_compute_line_model_array(im, x, y, vNeib, nNeib, vDir);

	rgbf res;
	int idmin, idmax;
	float alpha = -1;
	if(_select_extrema_with_spatial_constraint(vNeib, vDir, 
		g_fVarTol, idmin, idmax)) {
			prob = _linear_combinability(vNeib[idmin],vNeib[idmax], 
				vNeib[0], g_fVarTol, alpha);
			int ox_min = NEIB_POS[idmin*2], oy_min = NEIB_POS[idmin*2+1];
			int ox_max = NEIB_POS[idmax*2], oy_max = NEIB_POS[idmax*2+1];
			p_min = im->check(x+ox_min, y+oy_min) ? 
				rgb2rgbf(im->get(x+ox_min, y+oy_min)) : vBG;
			p_max = im->check(x+ox_max, y+oy_max) ? 
				rgb2rgbf(im->get(x+ox_max, y+oy_max)) : vBG;
	}
	delete [] vNeib;
	return alpha;
}

bool Segment::_select_extrema_with_spatial_constraint( rgbf *aNeib, const rgbf vDir, float g_fVarTol, int &idmin, int &idmax ) const
{
	float max_min = -1.f;
	for (int i = 0; i < 4; ++i) {

		int idmin1 = 0, idmax1 = 0; 
		float tmin1 = 1.f, tmax1 = -1.f;
		_find_local_min_max(aNeib, vDir, g_fVarTol, i, 
			0,3, tmin1, tmax1, idmin1, idmax1);

		int idmin2 = 0, idmax2 = 0; 
		float tmin2 = 1.f, tmax2 = -1.f;
		_find_local_min_max(aNeib, vDir, g_fVarTol, i, 
			5,8, tmin2, tmax2, idmin2, idmax2);

		if ((tmax1 < 0 || tmin2 > 0) && (tmax2 < 0 || tmin1 > 0))
			continue;

		int idmin_local = 0, idmax_local = 0; 
		float max_min_local = -1.f;
		float max1_min2 = tmax1 - tmin2, max2_min1 = tmax2 - tmin1;
		// prefer valid one with large max_min
		if (max1_min2 > max2_min1 || (tmax2 < 0 || tmin1 > 0)) {
			idmin_local = idmin2;
			idmax_local = idmax1;
			max_min_local = max1_min2;
		}
		else {
			idmin_local = idmin1;
			idmax_local = idmax2;
			max_min_local = max2_min1;
		}
		if (max_min_local > max_min) {
			max_min = max_min_local;
			idmin = idmin_local;
			idmax = idmax_local;
		}
	}
	return (max_min >= 0);
}

void Segment::_find_local_min_max( rgbf *aNeib, const rgbf vDir, const float g_fVarTol, const int par, const int start, const int end, float &tmin, float &tmax, int &idmin, int &idmax ) const
{
	for(int j = start; j < end; ++j) {
		float t;
		int cur = PAR_POS9[par][j];
		float dist2 = _project_dist(aNeib[cur], aNeib[0], vDir, t);
		if(dist2 > 585225*g_fVarTol)
			continue;
		if( t < tmin) {
			tmin = t;
			idmin = cur;
		}
		if(t > tmax) {
			tmax = t;
			idmax = cur;
		}
	}
}

float Segment::_linear_combinability( const rgbf &c1, const rgbf &c2, const rgbf &c, float g_fVarTol, float &t ) const
{
	rgbf c2c1 = c2 - c1;
	float c2c1len = c2c1.length();
	rgbf c2c1norm = c2c1 / c2c1len;
	t = (c-c1).dot(c2c1norm) / c2c1len;
	rgbf distVec = c - (c2c1 * t) - c1;
	float dist2 = distVec.dot(distVec);

	if(t < 0) {
		rgbf vToEndPt1 = c-c1;
		dist2 = std::max(dist2, vToEndPt1.dot(vToEndPt1));
	}
	if(t > 1) {
		rgbf vToEndPt2 = c-c2;
		dist2 = std::max(dist2, vToEndPt2.dot(vToEndPt2));
	}
	t = saturate(t);

	float linComb = exp(-dist2/65025.0f * 0.5f / g_fVarTol) ;
	return linComb;
}

void Segment::_compute_line_model_array( image<rgb> *im, const int x, const int y, rgbf *aNeib, const int nNeib, rgbf& vDir ) const
{
	rgbf *vNeibMean = new rgbf[nNeib];
	rgbf vCenter;
	for(int i = 0; i < nNeib; ++i) {
		int ox = NEIB_POS[i*2], oy = NEIB_POS[i*2+1];
		aNeib[i] = im->check(x+ox, y+oy) ? 
			rgb2rgbf(im->get(x+ox, y+oy)) : vBG;
		vCenter += aNeib[i];
	}
	vCenter /= float(nNeib);
	for(int i = 0; i < nNeib; ++i)
		vNeibMean[i] = aNeib[i] - vCenter;

	//EM iterations
	const int NITER = 3;
	rgbf p = vCenter.normalize(0.001f);
	for(int iter = 0; iter < NITER; ++iter) {
		rgbf t = rgbf(0.0001f,0.0001f,0.0001f);
		for(int i = 0; i < nNeib; ++i)
			t += vNeibMean[i] * vNeibMean[i].dot(p);
		p = t.normalize();
	}
	vDir = p;
	delete [] vNeibMean;
}

float Segment::_project_dist( const rgbf color, const rgbf vCenter, const rgbf vDir, float &t ) const
{
	rgbf v = color - vCenter;
	t = vDir.dot(v);
	rgbf vDist = color - (vCenter + vDir * t);
	return vDist.dot(vDist);
}

void Segment::trace( std::vector<float>& nodes, std::vector<std::vector<int>>& indices, std::vector<bool>& junction_map ) const
{
	int width = seg_map->width();
	int height = seg_map->height();
	//Compute F according to Eq. 5.2
	image<int> *F = new image<int>(width+1, height+1);
	for (int y = 0; y < height+1; y++) {
		int *row_px_F = F->access[y]; 
		int *row_px_seg = seg_map->access[y]; 
		for (int x = 0; x < width+1; x++) {
			if(y==0 || y==height || x==0 || x==width) {
				row_px_F[x] = 1;
			}
			else {
				int seg_val = row_px_seg[x];
				for(int delta_y=-1; delta_y<=0; delta_y++) {
					for(int delta_x=-1; delta_x<=0; delta_x++) {
						if(seg_map->get(x+delta_x,y+delta_y)!=seg_val){
							row_px_F[x] = 1;
							break;
						}
					}
				}
			}			
		}
	}

	int c = 0;
	image<int> *E = new image<int>(width+1, height+1);
	//scan all nodes to index the nodes
	for (int y = 0; y < height+1; y++) 
	{
		int *row_px_E = E->access[y]; 
		int *row_px_F = F->access[y]; 
		for (int x = 0; x < width+1; x++)
		{
			if(row_px_F[x] == 1) 
			{
				row_px_E[x] = c;
				nodes.push_back(float(x));
				nodes.push_back(float(y));
				c++;
			}
		}
	}

	static int offset_M[4][6] = {		
		-1,	0,	-1,	-1,	-1, 0,			
		0,	1,	-1,	0,	0,	0,			
		1,	0,	0,	0,	0,	-1,			
		0,	-1,	0,	-1,	-1,	-1 };
		indices.resize(count());
		for (int y = 0; y < height; y++) 
		{
			int *row_px_seg = seg_map->access[y]; 
			for (int x = 0; x < width; x++)
			{
				int k = row_px_seg[x];
				if(indices[k].size()==0) 
				{
					//backup F
					image<int> *temp_F = F->copy();
					//trace boundary_k
					int temp_y = y;
					int temp_x = x;
					indices[k].push_back(E->get(temp_x, temp_y));

					int cur_idx = 2;
					temp_F->set(temp_x, temp_y, 0);

					while (true) {
						temp_y += offset_M[cur_idx][0];
						temp_x += offset_M[cur_idx][1];
						indices[k].push_back(E->get(temp_x, temp_y));

						//compute M <- M_k,i,j(S,F) according to Eq. 5.4
						std::vector<int> possible_off_idx;
						for (int index = 0; index < 4; index++) {
							int val_y = temp_y + offset_M[index][0];
							int val_x = temp_x + offset_M[index][1];
							if (val_y<0 || val_y>height 
								|| val_x<0 || val_x>width)
								continue;
							int f_val = temp_F->get(val_x, val_y);

							val_y = temp_y + offset_M[index][2];
							val_x = temp_x + offset_M[index][3];

							int k_val_1 = -9999;	
							//assign a special value when out of bound

							if (val_y>=0 && val_y<height 
								&& val_x>=0 && val_x<width)
								k_val_1 = seg_map->get(val_x, val_y);

							val_y = temp_y + offset_M[index][4];
							val_x = temp_x + offset_M[index][5];

							int k_val_2 = -9999;	
							//assign a special value when out of bound

							if (val_y>=0 && val_y<height 
								&& val_x>=0 && val_x<width)
								k_val_2 = seg_map->get(val_x, val_y);

							if(f_val==1 && k==k_val_1 && k!=k_val_2)
								possible_off_idx.push_back(index);
						}

						if(possible_off_idx.size()==1) {
							temp_F->set(temp_x, temp_y, 0);
							cur_idx = possible_off_idx[0];
						}
						else if (possible_off_idx.size()>1)	{
							int min_val = 999;
							int min_idx = -1;
							for(int l_idx=0; 
								l_idx<(int)possible_off_idx.size(); 
								l_idx++) {
									int temp_val = (
										possible_off_idx.at(l_idx)-cur_idx+4
										)%4;
									if(temp_val < min_val) {
										min_val = temp_val;
										min_idx = l_idx;
									}
							}
							cur_idx = possible_off_idx.at(min_idx);
						}
						else break;
					}
					delete temp_F;
				}
			}
		}
		delete F;
		delete E;


		// Generate junction map for nodes
		const int patch_location[] = {-1,-1, -1,0, 0,-1, 0,0};
		const unsigned nodes_num = nodes.size() / 2;
		junction_map.resize(nodes_num, false);
		for (unsigned i = 0; i < nodes_num; i++) {
			std::set<int> segs;
			for (unsigned j = 0; j < 8; j+=2) {
				int x = (int)nodes[i*2] + patch_location[j];
				int y = (int)nodes[i*2+1] + patch_location[j+1];
				segs.insert(
					seg_map->check(x, y) ? seg_map->get(x, y) : -1
					);
			}
			if (segs.size() > 2)
				junction_map[i] = true;
		}
}

image<rgb> * Segment::vis_pseudo() const
{
	int width = seg_map->width();
	int height = seg_map->height();
	image<rgb> *output = new image<rgb>(width, height);
	rgb *colors = new rgb[width*height];
	for (int i = 0; i < width*height; i++) {
		colors[i].r = (unsigned char)rand();
		colors[i].g = (unsigned char)rand();
		colors[i].b = (unsigned char)rand();
	}
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
			output->set(x, y, colors[seg_map->get(x, y)]);
	delete [] colors;  
	return output;
}

image<rgb> * Segment::vis()
{
	int width = seg_map->width();
	int height = seg_map->height();
	image<rgb> *output = new image<rgb>(width, height);
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
			output->set(x, y, color_map[seg_map->get(x, y)]);
	return output;
}

std::map<int,rgb> Segment::get_color_map()
{
	return color_map;
}

int Segment::count() const
{
	return (int)color_map.size();
}

Segment::~Segment()
{
	delete seg_map;
}

Segment::Segment( image<rgb> *im )
{
	image<rgb> *blend_free_img = remove_blending(im);
	universe *seg = segment(blend_free_img); 

	// construct segmentation map
	int width = im->width();
	int height = im->height();
	seg_map = new image<int>(width, height);
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
			seg_map->set(x, y, seg->find(y * width + x));

	diagonal_fuse(im);

	// normalize segmentation map and build color map
	std::map<int,int> comp_to_index;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			int comp = seg_map->get(x, y);
			int index = comp_to_index[comp];				
			if (index == 0) {	// first meet this comp
				index = comp_to_index.size();
				// int index = comp_to_index[comp]; has actually
				// inserted this comp, so there's no need to minus
				// 1 here

				comp_to_index[comp] = index;
				color_map[index-1] = seg->color(comp);
			}	// else:		index = 1~count
			seg_map->set(x, y, index-1);
		}
	}  

	delete blend_free_img;
	delete seg;
}


rgb universe::color( int x )
{
	return rgbf2rgb(elts[x].sum / float(elts[x].size));
}

int universe::num_sets() const
{
	return num;
}

int universe::size( int x ) const
{
	return elts[x].size;
}

void universe::join( int x, int y )
{
	if (elts[x].rank > elts[y].rank) {
		elts[y].p = x;
		elts[x].size += elts[y].size;
		elts[x].sum += elts[y].sum;
	} else {
		elts[x].p = y;
		elts[y].size += elts[x].size;
		elts[y].sum += elts[x].sum;
		if (elts[x].rank == elts[y].rank)
			elts[y].rank++;
	}
	num--;
}

int universe::find( int x )
{
	int y = x;
	while (y != elts[y].p)
		y = elts[y].p;
	elts[x].p = y;
	return y;
}

universe::~universe()
{
	delete [] elts;
}

universe::universe( image<rgb> *im )
{
	int width = im->width();
	int height = im->height();
	num = width * height;
	elts = new uni_elt[num];
	for (int y = 0, i = 0; y < height; y++) {
		for (int x = 0; x < width; x++, i++) {
			elts[i].rank = 0;
			elts[i].size = 1;
			elts[i].p = i;
			elts[i].sum = rgb2rgbf(im->get(x, y));
		}
	}
}
