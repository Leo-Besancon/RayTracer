#pragma once

#include "IAnimatable.h"

class Light : public IAnimatable
{
	Vector _center;
	Vector _intensity;

public:
	Light(Vector center, Vector intensity) : _center(center), _intensity(intensity) {};
	~Light() { };

	Vector get_center() const { return _center; };
	void set_center(Vector c) { _center = c; };

	Vector get_intensity() const { return _intensity; };
	void set_intensity(Vector i) { _intensity = i; };

};
