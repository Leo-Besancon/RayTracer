#pragma once

#include "IAnimatable.h"

class Camera : public IAnimatable
{
	Vector _center;
	Vector _direction;
	Vector _up;
	double _focal;
	double _fov_degrees;
	int _H;
	int _W;

public:
	Camera::Camera(Vector center, Vector direction, Vector up, double fov_degrees, double focal, int H, int W) : _center(center), _direction(direction), _up(up), _focal(focal), _fov_degrees(fov_degrees), _H(H), _W(W) { }
	~Camera() {};

	Vector get_center() const { return _center; };
	void set_center(Vector center) { _center = center; };
	Vector get_direction() const { return _direction; };
	void set_direction(Vector direction) { _direction = direction; };
	Vector get_up() const { return _up; };
	void set_up(Vector up) { _up = up; };
	double get_focal() const { return _focal; };
	void set_focal(double focal) { _focal = focal; };
	double get_fov_degrees() const { return _fov_degrees; };
	void set_fov_degrees(double fov_degrees) { _fov_degrees = fov_degrees; };
	double get_fov_rad() const { return (_fov_degrees * M_PI / 180.0); };
	double get_depth() const { return ( _H / (2 * tan(_fov_degrees * M_PI / 180.0 / 2))); };
	int get_H() const { return _H; };
	void set_H(int H) { _H = H; };
	int get_W() const { return _W; };
	void set_W(int W) { _W = W; };

};