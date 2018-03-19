#pragma once

struct Animation
{
	double _start_time;
	double _end_time;

	Vector _translation;
	double _scale;
	double _rotation_x;
	Vector _rotation_center_x;
	double _rotation_y;
	Vector _rotation_center_y;
	double _rotation_z;
	Vector _rotation_center_z;

	Animation(double start_time = 0., double end_time = 0., Vector translation = Vector(0., 0., 0.), double scale = 1., double rotation_x = 0., Vector rotation_center_x = Vector(0., 0., 0.), double rotation_y = 0., Vector rotation_center_y = Vector(0., 0., 0.), double rotation_z = 0., Vector rotation_center_z = Vector(0., 0., 0.)) :
		_start_time(start_time), _end_time(end_time), _translation(translation), _scale(scale), _rotation_x(rotation_x), _rotation_center_x(rotation_center_x), _rotation_y(rotation_y), _rotation_center_y(rotation_center_y), _rotation_z(rotation_z), _rotation_center_z(rotation_center_z) { };
};

class IAnimatable
{
public:
	virtual ~IAnimatable() { };

	void addAnimation(Animation a) { _animations.push_back(a); };
	std::vector<Animation> get_animations() const { return _animations; };

	virtual Vector get_center() const = 0;

	virtual Vector get_center(double time) const final;

protected:
	std::vector<Animation> _animations;

	IAnimatable() { };
};
