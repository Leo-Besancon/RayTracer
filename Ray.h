#pragma once

#include "IAnimatable.h"

class Ray : public IAnimatable
{
	Vector _center;
	Vector _dir;

public:
	Ray() {};
	Ray(Vector center, Vector dir) : _center(center), _dir(dir) {};
	~Ray() {};

	Vector get_center() const { return _center; };
	void set_center(Vector c) { _center = c; };

	Vector get_dir() const { return _dir; };
	void set_dir(Vector dir) { _dir = dir; };

	double get_norm() const { return _dir.get_norm(); };
	Vector getPoint(double t) const { return (_center + _dir * t); };
	
	void normalize() { _dir = _dir.normalize(); };

	static Ray getNonAARay(int i, int j, int H, int W, double depth, Vector pos, Vector up = Vector(0, 1, 0), Vector direction = Vector(0, 0, -1));
	static Ray getAARay(int i, int j, int H, int W, double depth, Vector pos, Vector up = Vector(0, 1, 0), Vector direction = Vector(0, 0, -1));
	static Ray getAADoFRay(int i, int j, int H, int W, double depth, Vector pos, double focal, Vector up = Vector(0, 1, 0), Vector direction = Vector(0, 0, -1));

	static Ray getRandRay(Vector center, Vector dir);
	static Ray getRandRayAngleUniform(Vector center, double surface, Vector dir);

    static Ray getRandRayPhong(Vector center, Vector dir, double phong_exponent);

	void translate(Vector vec) { _center = _center + vec; };
	void rotate_x(double theta, Vector rotation_center);
	void rotate_y(double theta, Vector rotation_center);
	void rotate_z(double theta, Vector rotation_center);

	void apply_animations(std::vector<Animation> animations, double time);
	void reverse_animations(std::vector<Animation> animations, double time);

};
