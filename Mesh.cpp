#include "stdafx.h"

#include "Mesh.h"
#include "Scene.h"

#include <fstream>

std::string replace_slashes_in_string(std::string str) {
	std::string ret = "";
	for (int i = 0; i < str.length(); ++i) {
		if (str[i] == '/')
		{
			if (i < str.length() - 1 && str[i + 1] == '/')
			{
				ret.append(" 0 ");
			}
			else
			{
				ret.append(" ");
			}
		}
		else
		{
			ret.append(1, str[i]);
		}
	}
	return ret;
}

Mesh::Mesh(Material material, std::string filename, std::string texture_filename)
{
#ifdef LOGS
	std::cout << "Trying to parse : " << filename << ".obj" << std::endl;
#endif

	_material = material;

	std::ifstream mesh_file(filename + ".obj", std::ios::in);
	if (mesh_file)
	{
		for (std::string line; getline(mesh_file, line); )
		{
			if (line.length() < 3)
			{
				continue;
			}

			double a, b, c;
			std::istringstream toParse(line.substr(2));

			if (line[0] == 'v' && line[1] == ' ')
			{
				toParse >> a >> b >> c;
				Vector* v = new Vector(a, b, c);
				_v.push_back(v);

#ifdef LOGS
				std::cout << "Parsed v : " << a << " " << b << " " << c << std::endl;
#endif
			}
			else if (line[0] == 'v' && line[1] == 't')
			{
				static int cpt = 0;
				toParse >> a >> b;
				std::pair<double, double>* vt = new std::pair<double, double>(a, b);
				_vt.push_back(vt);

#ifdef LOGS
				std::cout << "Parsed vt : " << a << " " << b << std::endl;
#endif
			}
			else if (line[0] == 'v' && line[1] == 'n')
			{
				toParse >> a >> b >> c;
				Vector* vn = new Vector(a, b, c);
				_vn.push_back(vn);

#ifdef LOGS
				std::cout << "Parsed vn : " << a << " " << b << " " << c << std::endl;
#endif
			}
			else if (line[0] == 'f' && line[1] == ' ')
			{
				Face* f = new Face();
				std::istringstream toParse2(replace_slashes_in_string(line.substr(2)));

				toParse2 >> f->v1 >> f->vt1 >> f->vn1 >> f->v2 >> f->vt2 >> f->vn2 >> f->v3 >> f->vt3 >> f->vn3;

				_f.push_back(f);

#ifdef LOGS
				std::cout << "Parsed f : " << f.v1 << "/" << f.vt1 << "/" << f.vn1 << " " << f.v2 << "/" << f.vt2 << "/" << f.vn2 << " " << f.v3 << "/" << f.vt3 << "/" << f.vn3 << std::endl;
#endif
			}
			else {
#ifdef LOGS
				std::cout << "Could not parse line ! " << std::endl;
#endif
			}
		}
	}
	else
	{
		std::cout << "Unable to load the mesh file : " + filename + ".obj !" << std::endl;
	}

	computeMinMax();

	_center = (_min + _max) * 0.5;
	_scale = ((_max - _min) * 0.5).max();

	_recursion_level = 0;

	loadTexturefromBMP(texture_filename);
}

Mesh::~Mesh()
{
	free(_texture);
}

void Mesh::normalize()
{
	translate(_center * (-1.));
	scale(1. / _scale);

	computeMinMax();
}

void Mesh::translate(Vector vec)
{
	_center = _center + vec;

	for (int i = 0; i < _v.size(); i++)
	{
		*_v[i] = *_v[i] + vec;
	}
	computeMinMax();
}

void Mesh::scale(double a)
{
	_scale = _scale * a;

	Vector center = _center;

	translate(center * (-1.));
	for (int i = 0; i < _v.size(); i++)
	{
		*_v[i] = *_v[i] * a;
	}
	translate(center);

	computeMinMax();
}

void Mesh::loadTexturefromBMP(std::string filename)
{
	FILE *fd;
	fd = fopen(filename.c_str(), "rb");
	if (fd == NULL)
	{
		printf("Error: fopen failed \n");
		return;
	}

	unsigned char header[54];

	// Read header
	fread(header, sizeof(unsigned char), 54, fd);

	// Capture dimensions
	_texture_width = *(int*)&header[18];
	_texture_height = *(int*)&header[22];

	int padding = 0;

	// Calculate padding
	while ((_texture_width * 3 + padding) % 4 != 0)
	{
		padding++;
	}

	// Compute new width, which includes padding
	int widthnew = _texture_width * 3 + padding;

	// Allocate memory to store image data (non-padded)
	_texture = (unsigned char *)malloc(_texture_width *  _texture_height * 3 * sizeof(unsigned char));
	if (_texture == NULL)
	{
		printf("Error: Malloc failed\n");
		return;
	}

	// Allocate temporary memory to read widthnew size of data
	unsigned char* data = (unsigned char *)malloc(widthnew * sizeof(unsigned int));

	// Read row by row of data and remove padded data.
	for (int i = 0; i < _texture_height; i++)
	{
		// Read widthnew length of data
		fread(data, sizeof(unsigned char), widthnew, fd);

		// Retain width length of data, and swizzle RB component.
		// BMP stores in BGR format, my usecase needs RGB format
		for (int j = 0; j < _texture_width * 3; j += 3)
		{
			int index = (i *  _texture_width * 3) + (j);
			_texture[index + 0] = data[j + 2];
			_texture[index + 1] = data[j + 1];
			_texture[index + 2] = data[j + 0];
		}
	}

	free(data);
	fclose(fd);
}