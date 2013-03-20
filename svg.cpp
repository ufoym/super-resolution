#include "svg.h"

void SVGWriter::close()
{
	file << "</svg>" << std::endl;
	file.close();
}

void SVGWriter::writeText( float x, float y, std::string text, std::string color )
{
	file << "<text x=\"" << x << "\" y=\"" << y 
		<< "\" font-size=\"0.5\" fill=\"" << color 
		<< "\">" << text << "</text>";
}

void SVGWriter::writePolyline( std::vector<float>& points, std::string color )
{
	file << "<polyline points=\"";
	for (int i = 0; i < points.size(); i += 2) {
		file << points[i] << ","
			<< points[i + 1] << " ";
	}
	file << "\" style=\"fill-opacity:0; "
		<< "stroke-opacity:0.8; "
		<< "stroke:" << color << "; stroke-width:0.1\"/>"
		<< std::endl;
}

void SVGWriter::writePolygon( std::vector<float>& points, std::string color )
{
	file << "<polygon points=\"";
	for (int i = 0; i < points.size(); i += 2) {
		file << points[i] << ","
			<< points[i + 1] << " ";
	}
	file << "\" style=\"fill-opacity:0; "
		<< "stroke-opacity:0.8; "
		<< "stroke:" << color << "; stroke-width:0.1\"/>"
		<< std::endl;
}


void SVGWriter::writePolygonWithText( std::vector<float>& points, std::string color )
{
	file << "<polygon points=\"";
	for (int i = 0; i < points.size(); i += 2) {
		file << points[i] << ","
			 << points[i + 1] << " ";
	}
	file << "\" style=\"fill-opacity:0; "
		<< "stroke-opacity:0.8; "
		<< "stroke:" << color << "; stroke-width:0.1\"/>"
		<< std::endl;
	for (int i = 0; i < points.size(); i += 2) {
		file << "<text x=\"" << points[i] << "\" y=\"" << points[i + 1]
			 << "\" font-size=\"0.5\" fill=\"" << color 
			 << "\">" << i / 2 << "</text>\n";
	}
}


void SVGWriter::writeDots( std::vector<float>& points, std::string color )
{
	for (int i = 0 ; i < points.size(); i += 2) {
		file << "<circle cx=\"" << points[i] 
		<< "\" cy=\"" << points[i + 1] 
		<< "\" r=\"0.1\" fill=\"" << color << "\" />"
			<< std::endl;
	}
}

void SVGWriter::writeDot( float x, float y, std::string color )
{
	file << "<circle cx=\"" << x
		<< "\" cy=\"" << y
		<< "\" r=\"0.1\" fill=\"" << color << "\" />"
		<< std::endl;
}

SVGWriter::SVGWriter( std::string filename, std::string bg_img_path )
{
	file.open(filename, std::ios::out | std::ios::binary);
	file << "<?xml version=\"1.0\" standalone=\"no\"?>\n"
		"<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n"
		"\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"
		"<svg width=\"30\" height=\"30\" version=\"1.1\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns=\"http://www.w3.org/2000/svg\">\n";
	if (bg_img_path != "") {
		file << "<image y=\"0\" x=\"0\" xlink:href=\"" << bg_img_path << "\" />\n";
	}
}
