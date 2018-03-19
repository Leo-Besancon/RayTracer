#include "stdafx.h"
#include "Ray.h"


Ray Ray::getNonAARay(int i, int j, int H, int W, double depth, Vector pos, Vector up, Vector direction) {

	Vector right = direction.cross(up);

	Vector dir = right * (j - W / 2) + up * (i - H / 2) + direction * depth;
	Ray r = Ray(pos, dir);
	r.normalize();
	return r;
};

Ray Ray::getAARay(int i, int j, int H, int W, double depth, Vector pos, Vector up, Vector direction) {

	Vector right = direction.cross(up);

	double x = distrib(engine);
	double y = distrib(engine);
	double R = sqrt(-2 * log(x));
	double u = R * cos(2 * M_PI *y)*0.5;
	double v = R * sin(2 * M_PI *y)*0.5;

	Vector dir = right * (j - W / 2 + u - 0.5) + up * (i - H / 2 + v - 0.5) + direction * depth;
	Ray r = Ray(pos, dir);
	r.normalize();
	return r;
};

Ray Ray::getAADoFRay(int i, int j, int H, int W, double depth, Vector pos, double focal, Vector up, Vector direction) {

	Ray r1 = getAARay(i, j, H, W, depth, pos, up, direction);

	Vector right = direction.cross(up);


	Vector dir = r1.get_dir();
	double px = (distrib(engine) - 0.5) * 5.;
	double py = (distrib(engine) - 0.5) * 5.;

	Vector pos2 = pos + (right * px) + (up * py);

	Vector dir2 = (pos + ( dir * focal)) - pos2;

	Ray r = Ray(pos2, dir2);
	r.normalize();
	return r;
};

Ray Ray::getRandRay(Vector center, Vector n)
{
	double rand1 = distrib(engine);
	double rand2 = distrib(engine);

	double sqrt1 = sqrt(1 - rand2);

	double x_local = cos(2 * M_PI * rand1) * sqrt1;
	double y_local = sin(2 * M_PI * rand1) * sqrt1;
	double z_local = sqrt(rand2);

	Vector nx = Vector(n.get_z() * (-1), n.get_y(), n.get_x());
	Vector ny = n.cross(nx);

	Vector dir = nx * x_local + ny * y_local + n * z_local;

	Ray ray = Ray(center, dir);

	return ray;
}

Ray Ray::getRandRayAngleUniform(Vector center, double surface, Vector dir)
{
	double rayon = sqrt(surface / (4.0 * M_PI));

	double rand1 = distrib(engine);
	double rand2 = distrib(engine);

	double sqrt1 = sqrt(rand2);

	double x_local = cos(2 * M_PI * rand1) * sqrt1;
	double y_local = sin(2 * M_PI * rand1) * sqrt1;
	double z_local = sqrt(1 - rand2);

	Vector nx = Vector(dir.get_z() * (-1), dir.get_y(), dir.get_x());
	Vector ny = dir.cross(nx);

	Vector dir2 = (nx * x_local + ny * y_local + dir * z_local).normalize();
	Vector center2 = center + dir2 * rayon;

	Ray ray = Ray(center2, dir2);

	return ray;
};

Ray Ray::getRandRayPhong(Vector center, Vector dir, double phong_exponent)
{
	double rand1 = distrib(engine);
	double rand2 = distrib(engine);

	double phong_term = pow(rand2, 1. / (phong_exponent + 1.));

	double sqrt1 = sqrt(1-(phong_term * phong_term));

	double x_local = cos(2 * M_PI * rand1) * sqrt1;
	double y_local = sin(2 * M_PI * rand1) * sqrt1;
	double z_local = sqrt(1 - phong_term);

	Vector nx = Vector(dir.get_z() * (-1), dir.get_y(), dir.get_x());
	Vector ny = dir.cross(nx);

	Vector dir2 = (nx * x_local + ny * y_local + dir * z_local).normalize();

	Ray ray = Ray(center, dir2);

	return ray;
};

void Ray::rotate_x(double theta, Vector rotation_center)
{
	_dir = _dir.rotate_x(theta);

	translate(rotation_center * (-1.));
	_center = _center.rotate_x(theta);
	translate(rotation_center);
}

void Ray::rotate_y(double theta, Vector rotation_center)
{
	_dir = _dir.rotate_y(theta);

	translate(rotation_center * (-1.));
	_center = _center.rotate_y(theta);
	translate(rotation_center);
}
void Ray::rotate_z(double theta, Vector rotation_center)
{
	_dir = _dir.rotate_z(theta);

	translate(rotation_center * (-1.));
	_center = _center.rotate_z(theta);
	translate(rotation_center);
}

void Ray::apply_animations(std::vector<Animation> animations, double time)
{
	for (int i = 0; i < animations.size(); i++)
	{
		Animation a = animations[i];
		if (a._start_time >= time || a._start_time > a._end_time)
		{
			continue;
		}
		else
		{
			double progress;

			if (a._end_time > time)
				progress = (time - a._start_time) / (a._end_time - a._start_time);
			else
				progress = 1.;

			translate(a._translation * progress);
			rotate_x(a._rotation_x * progress, a._rotation_center_x);
			rotate_y(a._rotation_y * progress, a._rotation_center_y);
			rotate_z(a._rotation_z * progress, a._rotation_center_z);
		}
	}
}

void Ray::reverse_animations(std::vector<Animation> animations, double time)
{
	for (int i = 0; i < animations.size(); i++)
	{
		animations[i]._rotation_x = animations[i]._rotation_x * (-1);
		animations[i]._rotation_y = animations[i]._rotation_y * (-1);
		animations[i]._rotation_z = animations[i]._rotation_z * (-1);
		animations[i]._translation = animations[i]._translation * (-1);
	}
	apply_animations(animations, time);
}
