#pragma once

#include "IObject.h"

class Sphere : public IObject
{
	double _radius;
	Vector _center;

public:
	Sphere(Material material, Vector color, double radius);
	~Sphere() {};
	
	double get_rayon() const { return _radius; };
	void set_rayon(double radius) { _radius = radius; };

	Vector get_center() const { return _center; };
	void set_center(Vector center) { _center = center; };

	double get_surfaceArea() const { return (4.0 * M_PI * _radius * _radius); };

	std::tuple<double, Vector*, Vector> intersection(Ray ray, double t = std::numeric_limits<double>::max()) ;

};
