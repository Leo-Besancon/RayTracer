#pragma once

#include "Ray.h"
#include "IObject.h"


struct Face {
	int v1, v2, v3, vt1, vt2, vt3, vn1, vn2, vn3;
};

class Box : public IObject
{

protected :
	std::vector<Box*> _boxes;
	
	Vector _min, _max;
	int _recursion_level;

	std::vector<Vector*> _v;
	std::vector<std::pair<double, double>*> _vt;
	std::vector<Vector*> _vn;
	std::vector<Face*> _f;

	unsigned char* _texture;
	int  _texture_width, _texture_height;

	bool faceInsideBox(Face* f) const;
	bool HitBoundingBox(Ray r, double t) const;
	bool planeBoxOverlap(double normal[3], double vert[3], double maxbox[3]) const;

public:
	Box(Material material, Vector min, Vector max, const std::vector<Vector*>& v, const std::vector<std::pair<double, double>*>& vt, const std::vector<Vector*>& vn, std::vector<Face*> f, int recursion_level, unsigned char * texture = nullptr, int texture_width = 0, int texture_height = 0);
	~Box();

	virtual Vector get_center() const { return (_min + _max) * 0.5; };

	void clearBoxes();

	size_t nb_faces() const { return _f.size(); };

	void computeBoxes(int n_box);

	void computeMinMax();

	double intersection_plan(Ray ray, Vector P, Vector n) ;

	std::tuple<double, Vector*, Vector> intersection_face(Ray ray, Face* f) ;

	std::tuple<double, Vector*, Vector> intersection(Ray ray, double t) ;

	Vector get_texture_color(double u, double v) const;

protected:
	Box() {};
};

