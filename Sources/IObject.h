#pragma once

#include "IAnimatable.h"

#include "Ray.h"
#include "Material.h"

class Scene;

class IObject : public IAnimatable
{
public:
	virtual ~IObject() {};
	virtual std::tuple<double, Vector*, Vector> intersection(Ray r, double t = std::numeric_limits<double>::max()) = 0;
	virtual double get_surfaceArea() const { return 1.0; };
	virtual Vector get_center() const = 0;

	Material get_material() const { return _material; };
	void set_material(Material material) { _material = material; };

protected:
	IObject() {};
	Material _material;
};
