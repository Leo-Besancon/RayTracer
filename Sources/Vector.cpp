#include "stdafx.h"

Vector Vector::cross(Vector a) const
{
	double u1 = get_x();
	double u2 = get_y();
	double u3 = get_z();
	double v1 = a.get_x();
	double v2 = a.get_y();
	double v3 = a.get_z();

	return Vector(u2*v3-u3*v2,u3*v1-u1*v3,u1*v2-u2*v1);
}

Vector Vector::normalize()
{
	double norm = get_norm();
	Vector vect = Vector(_x / norm, _y / norm, _z / norm);

	return vect;
}

std::ostream& operator << (std::ostream& os, const Vector& A)
{
	os << "( " << A.get_x() << ", " << A.get_y() << ", " << A.get_z() << " )";
	return os;
}

void Vector::operator+=(const Vector& a)
{
	_x += a.get_x();
	_y += a.get_y();
	_z += a.get_z();
}

Vector Vector::rotate_x(double theta)
{
	double theta_rad = theta * M_PI / 180.;
	double x = _x;
	double y = cos(theta_rad) * _y - sin(theta_rad) * _z;
	double z = sin(theta_rad) * _y + cos(theta_rad) * _z;

	return Vector(x, y, z);
}

Vector Vector::rotate_y(double theta)
{
	double theta_rad = theta * M_PI / 180.;

	double x = cos(theta_rad) * _x + sin(theta_rad) * _z;
	double y = _y;
	double z = - sin(theta_rad) * _x + cos(theta_rad) * _z;

	return Vector(x, y, z);
}

Vector Vector::rotate_z(double theta)
{
	double theta_rad = theta * M_PI / 180.;

	double x = cos(theta_rad) * _x - sin(theta_rad) * _y;
	double y = sin(theta_rad) * _x + cos(theta_rad) * _y;
	double z = _z;

	return Vector(x, y, z);
}
