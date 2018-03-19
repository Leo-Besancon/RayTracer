#pragma once

class Vector
{
	double _x;
	double _y;
	double _z;

public:
	Vector() { };
	Vector(double x, double y, double z) : _x(x), _y(y), _z(z) { };
	~Vector() { };

	double get_x() const { return _x; };
	double get_y() const { return _y; };
	double get_z() const { return _z; };

	void set_x(double x) { _x = x; };
	void set_y(double y) { _y = y; };
	void set_z(double z) { _z = z; };

	double get_norm() const { return sqrt(_x*_x + _y * _y + _z * _z); };
	double get_norm2() const { return (_x*_x + _y * _y + _z * _z); };

	double min() const { return std::min(std::min(_x, _y), _z); }
	double max() const { return std::max(std::max(_x, _y), _z); }
	double moy() const { return ((_x + _y + _z) / 3.); }


	Vector min(double a) const { return Vector(std::min(a, _x), std::min(a, _y), std::min(a, _z)); }
	Vector max(double a) const { return Vector(std::max(a, _x), std::max(a, _y), std::max(a, _z)); }

	double scalar(Vector a) const { return (_x * a.get_x() + _y * a.get_y() + _z * a.get_z()); };

	Vector cross(Vector a) const;

	Vector operator+(const Vector& a) const { return  Vector(_x + a.get_x(), _y + a.get_y(), _z + a.get_z()); };
	void operator+=(const Vector& a);
	Vector operator-(const Vector& a) const { return Vector(_x - a.get_x(), _y - a.get_y(), _z - a.get_z()); };
	Vector operator*(double a) const { return Vector(_x*a, _y*a, _z*a); };
	Vector operator/(double a) const { return Vector(_x/a, _y/a, _z/a); };
	Vector operator*(const Vector& a) const { return Vector(_x * a.get_x(), _y * a.get_y(), _z * a.get_z());};

	friend std::ostream& operator <<(std::ostream& sortie, const Vector& A);

	Vector normalize();
	Vector rotate_x(double theta);
	Vector rotate_y(double theta);
	Vector rotate_z(double theta);
};
