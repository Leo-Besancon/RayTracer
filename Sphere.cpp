#include "stdafx.h"
#include "Sphere.h"
#include "Scene.h"


Sphere::Sphere(Material material, Vector center, double radius) : _center(center), _radius(radius)
{
	_material = material;
}

std::tuple<double, Vector*, Vector> Sphere::intersection(Ray ray, double t) 
{
	t = -1;

	Vector* n = nullptr;

	Vector CO = ray.get_center() - get_center();
	double a = ray.get_dir().get_norm2();
	double b = 2 * ray.get_dir().scalar(CO);
	double c = CO.get_norm2() - _radius * _radius;

	double delta = b * b - 4 * a*c;
	
	if (delta >= 0)
	{
		double t1 = (- b - sqrt(delta)) / (2 * a);
		double t2 = (- b + sqrt(delta)) / (2 * a);
		if (t2 >= 0)
		{
			if (t1 < 0)
			{
				t = t2;
			}
			else
			{
				t = t1;
			}
			n = new Vector();
			*n = (ray.getPoint(t) - get_center()).normalize();
		}
	}

	return std::tuple<double, Vector*, Vector>(t,n,_material._color);
}
