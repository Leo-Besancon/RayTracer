#include "stdafx.h"
#include "Material.h"


Material Material::createMirror(Vector specular_color)
{
	Material m = Material();
	m._mirror = true;
	m._specular_color = specular_color;
	return m;
}

Material Material::createTransparent(Vector specular_color, double n_object)
{
	Material m = Material();
	m._transparent = true;
	m._specular_color = specular_color;
	m._n_object = n_object;
	return m;
}

Material Material::createEmissive(Vector color, double emissivity)
{
	Material m = Material();
	m._color = color;
	m._emissive = true;
	m._emissivity = emissivity;
	return m;
}

Material Material::createDiffuse(Vector color)
{
	Material m = Material();
	m._color = color;
	return m;
}

Material Material::createPhong(Vector color, Vector specular_color, double phong_exponent)
{
	Material m = Material();
	m._color = color;
	m._specular_color = specular_color;
	m._phong = true;
	m._phong_exponent = phong_exponent;
	return m;
}
