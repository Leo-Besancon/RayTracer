#include "stdafx.h"
#include "Box.h"

#include "macros_snippets.h"

Box::Box(Material material, Vector min, Vector max, const std::vector<Vector*> & v, const std::vector<std::pair<double, double>*> & vt, const std::vector<Vector*> & vn, std::vector<Face*> f, int recursion_level, unsigned char* texture, int  texture_width, int texture_height, bool texture_procedural) :
	_min(min), _max(max), _recursion_level(recursion_level), _v(v), _vt(vt), _vn(vn), _texture(texture), _texture_width(texture_width), _texture_height(texture_height), _texture_procedural(texture_procedural)
{
	_material = material;

	for (int i = 0; i < f.size(); i++)
	{
		if (faceInsideBox(f[i]))
		{
			_f.push_back(f[i]);
		}
	}
}

Box::~Box()
{
	clearBoxes();
}

void Box::clearBoxes()
{
	for (int i = 0; i < _boxes.size(); i++)
	{
		if (_boxes[i] != nullptr)
			delete _boxes[i];
	}
	_boxes.clear();
}

void Box::computeBoxes(int n_box)
{
	Vector min, max;

	clearBoxes();

	if (nb_faces() >= 20)
	{
		for (int i = 0; i < n_box; i++)
		{
			for (int j = 0; j < n_box; j++)
			{
				for (int k = 0; k < n_box; k++)
				{
					min.set_x(_min.get_x() + (i * 1.0 / n_box) * (_max.get_x() - _min.get_x()));
					min.set_y(_min.get_y() + (j * 1.0 / n_box) * (_max.get_y() - _min.get_y()));
					min.set_z(_min.get_z() + (k * 1.0 / n_box) * (_max.get_z() - _min.get_z()));
					max.set_x(_min.get_x() + ((i + 1.0) * 1.0 / n_box) * (_max.get_x() - _min.get_x()));
					max.set_y(_min.get_y() + ((j + 1.0)  * 1.0 / n_box) * (_max.get_y() - _min.get_y()));
					max.set_z(_min.get_z() + ((k + 1.0) * 1.0 / n_box) * (_max.get_z() - _min.get_z()));

					Box* box = new Box(_material, min, max, _v, _vt, _vn, _f, _recursion_level + 1, _texture, _texture_width, _texture_height, _texture_procedural);

					if (box->nb_faces() == 0)
					{
						delete box;
					}
					else
					{
						_boxes.push_back(box);
					}

				}
			}
		}

		if (_boxes.size() <= 2)
		{
			clearBoxes();
		}

		for (int i = 0; i < _boxes.size(); i++)
		{
			_boxes[i]->computeBoxes(n_box);
		}

	}
}

void Box::computeMinMax()
{
	double cur_min_x = std::numeric_limits<double>::max();
	double cur_min_y = std::numeric_limits<double>::max();
	double cur_min_z = std::numeric_limits<double>::max();
	double cur_max_x = std::numeric_limits<double>::min();
	double cur_max_y = std::numeric_limits<double>::min();
	double cur_max_z = std::numeric_limits<double>::min();

	for (int i = 0; i < _f.size(); i++)
	{
		double candidate_min_x = std::min(std::min(_v[_f[i]->v1 - 1]->get_x(), _v[_f[i]->v2 - 1]->get_x()), _v[_f[i]->v3 - 1]->get_x());
		double candidate_min_y = std::min(std::min(_v[_f[i]->v1 - 1]->get_y(), _v[_f[i]->v2 - 1]->get_y()), _v[_f[i]->v3 - 1]->get_y());
		double candidate_min_z = std::min(std::min(_v[_f[i]->v1 - 1]->get_z(), _v[_f[i]->v2 - 1]->get_z()), _v[_f[i]->v3 - 1]->get_z());
		double candidate_max_x = std::max(std::max(_v[_f[i]->v1 - 1]->get_x(), _v[_f[i]->v2 - 1]->get_x()), _v[_f[i]->v3 - 1]->get_x());
		double candidate_max_y = std::max(std::max(_v[_f[i]->v1 - 1]->get_y(), _v[_f[i]->v2 - 1]->get_y()), _v[_f[i]->v3 - 1]->get_y());
		double candidate_max_z = std::max(std::max(_v[_f[i]->v1 - 1]->get_z(), _v[_f[i]->v2 - 1]->get_z()), _v[_f[i]->v3 - 1]->get_z());

		if (cur_min_x > candidate_min_x)
			cur_min_x = candidate_min_x;
		if (cur_min_y > candidate_min_y)
			cur_min_y = candidate_min_y;
		if (cur_min_z > candidate_min_z)
			cur_min_z = candidate_min_z;
		if (cur_max_x < candidate_max_x)
			cur_max_x = candidate_max_x;
		if (cur_max_y < candidate_max_y)
			cur_max_y = candidate_max_y;
		if (cur_max_z < candidate_max_z)
			cur_max_z = candidate_max_z;
	}

	_min = Vector(cur_min_x, cur_min_y, cur_min_z);
	_max = Vector(cur_max_x, cur_max_y, cur_max_z);

	for (int i = 0; i < _boxes.size(); i++)
	{
		_boxes[i]->computeMinMax();
	}
}

double Box::intersection_plan(Ray ray, Vector P, Vector n) 
{
	double t = -1;
	double t2;
	double scalar = ray.get_dir().scalar(n);

	if (scalar != 0)
	{
		Vector OP = ray.get_center() - P;
		t2 = -OP.scalar(n) / scalar;
		if (t2 > 0)
			t = t2;
	}

	return t;
}

std::tuple<double, Vector*, Vector> Box::intersection_face(Ray ray, Face* f) 
{
	double t = -1;
	Vector* n_bary = nullptr;
	Vector color_bary = _material._color;

	// Compute the face normal

	Vector v12 = *_v[f->v2 - 1] - *_v[f->v1 - 1];
	Vector v13 = *_v[f->v3 - 1] - *_v[f->v1 - 1];

	Vector n = v12.cross(v13);

	if (n.scalar(*_vn[f->vn1 - 1]) < 0)
		n = n * (-1);

	double t_plan = intersection_plan(ray, *_v[f->v1 - 1], n);

	if (t_plan != -1)
	{

		Vector P = ray.getPoint(t_plan);

		Vector v0 = *_v[f->v3 - 1] - *_v[f->v1 - 1];
		Vector v1 = *_v[f->v2 - 1] - *_v[f->v1 - 1];
		Vector v2 = P - *_v[f->v1 - 1];

		// Compute dot products
		double dot00 = v0.scalar(v0);
		double dot01 = v0.scalar(v1);
		double dot02 = v0.scalar(v2);
		double dot11 = v1.scalar(v1);
		double dot12 = v1.scalar(v2);

		// Compute barycentric coordinates
		double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
		double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
		double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

		// Check if point is in triangle
		if ((u >= 0) && (v >= 0) && (u + v < 1))
		{
			t = t_plan;
			n_bary = new Vector();
			*n_bary = *_vn[f->vn1 - 1] * (1 - u - v) + *_vn[f->vn2 - 1] * v + *_vn[f->vn3 - 1] * u;

			if (_texture_procedural == true)
			{
				color_bary = get_texture_procedural_color(*n_bary);
			}
			else
			{
				if (_texture != nullptr)
				{
					if (f->vt1 - 1 >= 0 && f->vt2 - 1 >= 0 && f->vt3 - 1 >= 0)
					{
						std::pair<double, double> uv_color1 = *_vt[f->vt1 - 1];
						std::pair<double, double> uv_color2 = *_vt[f->vt2 - 1];
						std::pair<double, double> uv_color3 = *_vt[f->vt3 - 1];

						double u_color1 = uv_color1.first;
						double u_color2 = uv_color2.first;
						double u_color3 = uv_color3.first;

						double v_color1 = uv_color1.second;
						double v_color2 = uv_color2.second;
						double v_color3 = uv_color3.second;

						double u_bary = u_color1 * (1 - u - v) + u_color2 * v + u_color3 * u;;
						double v_bary = v_color1 * (1 - u - v) + v_color2 * v + v_color3 * u;;

						color_bary = get_texture_color(u_bary, v_bary);
					}
					else
					{
						color_bary = _material._color;
					}
			}

			}
		}
	}
	return std::tuple<double, Vector*, Vector>(t, n_bary, color_bary);
}

std::tuple<double, Vector*, Vector> Box::intersection(Ray ray, double t) 
{
	double t2;

	Vector color = _material._color;
	Vector color2;

	Vector* n = nullptr;
	Vector* n2 = nullptr;

	// Returns false if the discovered t is smaller than the distance to the box
	if (HitBoundingBox(ray, t) == false)
	{
		t = -1;
		return std::tuple<double, Vector*, Vector>(t, n, color);
	}

	// If it crosses the box, check with children boxes

	for (int i = 0; i < _boxes.size(); i++)
	{
		if (_boxes[i] != nullptr)
		{
			std::tuple<double, Vector*, Vector> tuple = _boxes[i]->intersection(ray,t);
			t2 = std::get<0>(tuple);
			n2 = std::get<1>(tuple);
			color2 = std::get<2>(tuple);

			if (t2 != -1 && t2 < t)
			{
				t = t2;
				n = n2;
				color = color2;
			}
			else if (n2 != nullptr)
			{
				delete n2;
			}
		}
	}

	// If no children boxes, Compute with all the faces

	if (_boxes.size() == 0)
	{
		for (int i = 0; i < _f.size(); i++)
		{
			std::tuple<double, Vector*, Vector> tuple = intersection_face(ray, _f[i]);
			t2 = std::get<0>(tuple);
			n2 = std::get<1>(tuple);
			color2 = std::get<2>(tuple);

			if (t2 != -1 && t2 < t)
			{
				t = t2;
				n = n2;
				color = color2;
			}
			else if (n2 != nullptr)
			{
				delete n2;
			}
		}
	}

	if (t == std::numeric_limits<double>::max())
		t = -1;

	return std::tuple<double, Vector*, Vector>(t, n, color);
}

// Quick check with the 6 planes, cf snippet : https://github.com/erich666/GraphicsGems/blob/master/gems/RayBox.c
bool Box::HitBoundingBox(Ray r, double t) const
{
	double minB[NUMDIM], maxB[NUMDIM];		/*box */
	double origin[NUMDIM], dir[NUMDIM];		/*ray */
	double coord[NUMDIM];

	minB[0] = _min.get_x();
	minB[1] = _min.get_y();
	minB[2] = _min.get_z();
	maxB[0] = _max.get_x();
	maxB[1] = _max.get_y();
	maxB[2] = _max.get_z();
	origin[0] = r.get_center().get_x();
	origin[1] = r.get_center().get_y();
	origin[2] = r.get_center().get_z();
	dir[0] = r.get_dir().get_x();
	dir[1] = r.get_dir().get_y();
	dir[2] = r.get_dir().get_z();

	bool inside = true;
	char quadrant[NUMDIM];
	register int i;
	int whichPlane;
	double maxT[NUMDIM];
	double candidatePlane[NUMDIM];

	/* Find candidate planes; this loop can be avoided if
	rays cast all from the eye(assume perpsective view) */
	for (i = 0; i < NUMDIM; i++)
	{
		if (origin[i] < minB[i]) {
			quadrant[i] = LEFT;
			candidatePlane[i] = minB[i];
			inside = false;
		}
		else if (origin[i] > maxB[i]) {
			quadrant[i] = RIGHT;
			candidatePlane[i] = maxB[i];
			inside = false;
		}
		else {
			quadrant[i] = MIDDLE;
		}
	}

	/* Ray origin inside bounding box */
	if (inside) {
		return true;
	}

	/* Calculate T distances to candidate planes */
	for (i = 0; i < NUMDIM; i++)
	{
		if (quadrant[i] != MIDDLE && dir[i] != 0.)
			maxT[i] = (candidatePlane[i] - origin[i]) / dir[i];
		else
			maxT[i] = -1.;
	}

	/* Get largest of the maxT's for final choice of intersection */
	whichPlane = 0;

	for (i = 1; i < NUMDIM; i++)
	{
		if (maxT[whichPlane] < maxT[i])
			whichPlane = i;
	}

	/* Check final candidate actually inside box */
	if (maxT[whichPlane] < 0.)
		return false;

	for (i = 0; i < NUMDIM; i++)
	{
		if (whichPlane != i) 
		{
			coord[i] = origin[i] + maxT[whichPlane] * dir[i];
			if (coord[i] < minB[i] || coord[i] > maxB[i])
				return false;
		}
		else {
			coord[i] = candidatePlane[i];
		}
	}

	/* Check if the intersection is behind the discovered point*/
	if ((Vector(coord[0], coord[1], coord[2]) + (r.get_center() * (-1))).get_norm2() > t*t)
		return false;

	return true;				/* ray hits box */
}

bool Box::planeBoxOverlap(double normal[3], double vert[3], double maxbox[3]) const	// -NJMP-
{
	int q;
	double vmin[3], vmax[3], v;
	for (q = X; q <= Z; q++)
	{
		v = vert[q];					// -NJMP-
		if (normal[q]>0.0f)
		{
			vmin[q] = -maxbox[q] - v;	// -NJMP-
			vmax[q] = maxbox[q] - v;	// -NJMP-
		}
		else
		{
			vmin[q] = maxbox[q] - v;	// -NJMP-
			vmax[q] = -maxbox[q] - v;	// -NJMP-
		}
	}
	if (DOT(normal, vmin) > 0.0) return false;	// -NJMP-
	if (DOT(normal, vmax) >= 0.0) return true;	// -NJMP-

	return false;
}

// Compute faces inside Box; Cf snippet :  http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/tribox3.txt
bool Box::faceInsideBox(Face* f) const
{
	double boxcenter[3];
	double boxhalfsize[3];
	double triverts[3][3];

	Vector center = (_min + _max) * 0.5;
	Vector halfsize = (_max - _min) * 0.5;

	boxcenter[0] = center.get_x();
	boxcenter[1] = center.get_y();
	boxcenter[2] = center.get_z();

	boxhalfsize[0] = halfsize.get_x();
	boxhalfsize[1] = halfsize.get_y();
	boxhalfsize[2] = halfsize.get_z();

	triverts[0][0] = _v[f->v1 - 1]->get_x();
	triverts[0][1] = _v[f->v1 - 1]->get_y();
	triverts[0][2] = _v[f->v1 - 1]->get_z();
	triverts[1][0] = _v[f->v2 - 1]->get_x();
	triverts[1][1] = _v[f->v2 - 1]->get_y();
	triverts[1][2] = _v[f->v2 - 1]->get_z();
	triverts[2][0] = _v[f->v3 - 1]->get_x();
	triverts[2][1] = _v[f->v3 - 1]->get_y();
	triverts[2][2] = _v[f->v3 - 1]->get_z();


	/*    use separating axis theorem to test overlap between triangle and box */
	/*    need to test for overlap in these directions: */
	/*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
	/*       we do not even need to test these) */
	/*    2) normal of the triangle */
	/*    3) crossproduct(edge from tri, {x,y,z}-directin) */
	/*       this gives 3x3=9 more tests */

	double v0[3], v1[3], v2[3];
	//   float axis[3];

	double min, max, p0, p1, p2, rad, fex, fey, fez;		// -NJMP- "d" local variable removed
	double normal[3], e0[3], e1[3], e2[3];

	/* This is the fastest branch on Sun */
	/* move everything so that the boxcenter is in (0,0,0) */

	SUB(v0, triverts[0], boxcenter);
	SUB(v1, triverts[1], boxcenter);
	SUB(v2, triverts[2], boxcenter);

	/* compute triangle edges */

	SUB(e0, v1, v0);      /* tri edge 0 */
	SUB(e1, v2, v1);      /* tri edge 1 */
	SUB(e2, v0, v2);      /* tri edge 2 */
						  /* Bullet 3:  */
						  /*  test the 9 tests first (this was faster) */
	fex = fabs(e0[X]);
	fey = fabs(e0[Y]);
	fez = fabs(e0[Z]);
	AXISTEST_X01(e0[Z], e0[Y], fez, fey);
	AXISTEST_Y02(e0[Z], e0[X], fez, fex);
	AXISTEST_Z12(e0[Y], e0[X], fey, fex);

	fex = fabs(e1[X]);
	fey = fabs(e1[Y]);
	fez = fabs(e1[Z]);

	AXISTEST_X01(e1[Z], e1[Y], fez, fey);
	AXISTEST_Y02(e1[Z], e1[X], fez, fex);
	AXISTEST_Z0(e1[Y], e1[X], fey, fex);

	fex = fabs(e2[X]);
	fey = fabs(e2[Y]);
	fez = fabs(e2[Z]);

	AXISTEST_X2(e2[Z], e2[Y], fez, fey);
	AXISTEST_Y1(e2[Z], e2[X], fez, fex);
	AXISTEST_Z12(e2[Y], e2[X], fey, fex);

	/* Bullet 1: */
	/*  first test overlap in the {x,y,z}-directions */
	/*  find min, max of the triangle each direction, and test for overlap in */
	/*  that direction -- this is equivalent to testing a minimal AABB around */
	/*  the triangle against the AABB */

	/* test in X-direction */

	FINDMINMAX(v0[X], v1[X], v2[X], min, max);

	if (min>boxhalfsize[X] || max<-boxhalfsize[X]) return false;

	/* test in Y-direction */
	FINDMINMAX(v0[Y], v1[Y], v2[Y], min, max);
	if (min>boxhalfsize[Y] || max<-boxhalfsize[Y]) return false;

	/* test in Z-direction */
	FINDMINMAX(v0[Z], v1[Z], v2[Z], min, max);
	if (min>boxhalfsize[Z] || max<-boxhalfsize[Z]) return false;

	/* Bullet 2: */
	/*  test if the box intersects the plane of the triangle */
	/*  compute plane equation of triangle: normal*x+d=0 */

	CROSS(normal, e0, e1);

	// -NJMP- (line removed here)

	if (!planeBoxOverlap(normal, v0, boxhalfsize)) return false;	// -NJMP-

	return true;   /* box and triangle overlaps */
}

Vector Box::get_texture_color(double u, double v) const
{
	Vector color;

	int i = (int)round((v *_texture_height));
	int j = (int)round((u *_texture_width)) * 3;
	int index = (i * _texture_width * 3) + j;

	color.set_x(_texture[index] / 255.);
	color.set_y(_texture[index + 1] / 255.);
	color.set_z(_texture[index + 2] / 255.);

	return color;
}

Vector Box::get_texture_procedural_color(Vector n) const
{
	Vector color;

	color.set_x(abs(n.get_x() - 2 * n.get_y()));
	color.set_y(abs(n.get_x() - 2 * cos(n.get_y())));
	color.set_z(abs(n.get_z() - 2 * sin(n.get_y())));

	color = color.normalize();

	return color;
}