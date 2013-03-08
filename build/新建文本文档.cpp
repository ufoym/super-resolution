chizhangsysu
	杨铭 16:19:24
https://github.com/ufoym/super-resolution
你看看能否访问

	另外，
http://rogerdudler.github.com/git-guide/index.zh.html

这是git的一个非常清爽的教程，你可以花几分钟看一下，命令记不记住没关系，关键是理解git的工作流程，有什么疑问可以问我~
	杨铭 16:27:31
http://windows.github.com/
杨铭 16:27:37
	install this
	暂时无名 19:40:32

	nodes里面的元素是按x1,y1,x2,y2,...,xn,yn,...这样储存的吗？
	杨铭 19:49:16
	yes
	杨铭 20:47:28


	//---------------------------------------------------------------------
	// visualization

	void _visualize_poly(std::ofstream &file) const
{
	std::vector<Loop>::const_iterator loop_iter;
	std::vector<unsigned>::const_iterator node_idx_iter;
	for (loop_iter = loops.begin();
		loop_iter != loops.end(); loop_iter++) {
			file << "<polygon points=\"";
			for (node_idx_iter = loop_iter->node_indice.begin();
				node_idx_iter != loop_iter->node_indice.end();
				node_idx_iter++)
				file<< nodes[*node_idx_iter].x << ","
				<< nodes[*node_idx_iter].y << " ";
			file<< "\" style=\"fill-opacity:0; "
				<< "stroke-opacity:0.8; "
				<< "stroke:#00ff00; stroke-width:0.2\"/>"
				<< std::endl;
	}
}

void _visualize_node(std::ofstream &file) const
{
	for (unsigned node_idx = 0; node_idx<nodes.size(); ++node_idx)
		file<< "<circle cx=\"" << nodes[node_idx].x << "\" "
		<< "cy=\"" << nodes[node_idx].y << "\" " << "r=\""
		<< nodes[node_idx].section_indice.size() / 10.f
		<< "\" "
		<< "fill=\"red\" />\n"

		<< "<text x=\"" << nodes[node_idx].x << "\" "
		<< "y=\"" << nodes[node_idx].y << "\" "
		<< "font-size=\"0.5\" fill=\"red\">"
		<< node_idx << "</text>"
		<< std::endl;
}

void _visualize_bezier(std::ofstream &file) const
{
	std::vector<std::vector<float>> bezier_contour \
		= _bezier_representation();
	std::vector<std::vector<float>>::iterator bezier_contour_iter;
	std::vector<float>::iterator node_iter;
	for (bezier_contour_iter = bezier_contour.begin();
		bezier_contour_iter != bezier_contour.end();
		bezier_contour_iter++) {
			file << "<path d=\"";
			unsigned num_seg = bezier_contour_iter->size() / 8;
			node_iter = bezier_contour_iter->begin();
			for (unsigned i = 0; i < num_seg; ++i) {
				if (i == 0) {
					file<< "M"
						<< bezier_contour_iter->at(0) << ","
						<< bezier_contour_iter->at(1) << " ";
				}
				file<< "C"
					<< bezier_contour_iter->at(i*8+2) << ","
					<< bezier_contour_iter->at(i*8+3) << " "
					<< bezier_contour_iter->at(i*8+4) << ","
					<< bezier_contour_iter->at(i*8+5) << " "
					<< bezier_contour_iter->at(i*8+6) << ","
					<< bezier_contour_iter->at(i*8+7) << " ";
			}

			file<< "\" style=\"fill-opacity:0; "
				<< "stroke-opacity:0.8; "
				<< "stroke:#0000ff; stroke-width:0.2\"/>"
				<< std::endl;
	}
}
