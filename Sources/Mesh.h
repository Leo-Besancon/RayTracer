#pragma once

#include "Box.h"

class Mesh : public Box
{
	Vector _center;
	double _scale;

public:
	void constr(Material material, std::string filename);
	Mesh(Material material, std::string filename, std::string texture_filename = "PROCEDURAL");
	~Mesh();

	std::vector<Vector*> get_v() const { return _v; };
	std::vector<std::pair<double, double>*> get_vt() const { return _vt; };
	std::vector<Vector*> get_vn() const { return _vn; };
	std::vector<Face*> get_f() const { return _f; };

	Vector get_center() { return _center; };
	double get_scale() { return _scale; };

	void normalize();
	void translate(Vector vec);
	void scale(double a);

	void loadTexturefromBMP(std::string filename);
};