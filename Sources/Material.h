#pragma once

struct Material
{
	Vector _color;
	bool _mirror;
	Vector _specular_color;
	bool _transparent;
	double _n_object;
	bool _emissive;
	double _emissivity;
	bool _phong;
	double _phong_exponent;

	Material() :
		_color(Vector(0.,0.,0.)), _mirror(false), _specular_color(Vector(0., 0., 0.)), _transparent(false), _n_object(1.0), _emissive(false), _emissivity(0.), _phong(false), _phong_exponent(1.) {};
	Material(Vector color, bool mirror, Vector specular_color, bool transparent, double n_object, bool emissive, double emissivity, bool phong, double phong_exponent) :
		_color(color), _mirror(mirror), _specular_color(specular_color), _transparent(transparent), _n_object(n_object), _emissive(emissive), _emissivity(emissivity), _phong(phong), _phong_exponent(phong_exponent) {};

	static Material Material::createMirror(Vector specular_color);
	static Material Material::createTransparent(Vector specular_color, double n_object);
	static Material Material::createEmissive(Vector color, double emissivity);
	static Material Material::createDiffuse(Vector color);
	static Material Material::createPhong(Vector color, Vector specular_color, double phong_exponent);
};