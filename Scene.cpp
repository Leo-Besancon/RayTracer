#include "stdafx.h"
#include "Scene.h"


Scene::~Scene()
{
	// Delete Lights
	for (int i = 0; i < _lights.size(); i++)
	{
		if (_lights[i] != nullptr)
			delete _lights[i];
	}
	_lights.clear();

	// Delete Objects
	for (int i = 0; i < _objects.size(); i++)
	{
		if (_objects[i] != nullptr)
			delete _objects[i];
	}
	_objects.clear();
}

std::tuple<Vector*, Vector*, IObject*, Vector> Scene::computeIntersection(Ray r, double time) const
{
	IObject* current_Object = nullptr;
	Vector* intersection = nullptr;
	Vector* current_n = nullptr;
	Vector current_color;

	double current_t = std::numeric_limits<double>::max();

	for (int k = 0; k < _objects.size(); k++)
	{
		IObject* o = _objects[k];

		r.reverse_animations(o->get_animations(), time);
		std::tuple<double, Vector*, Vector> tuple = o->intersection(r);
		r.apply_animations(o->get_animations(), time);

		double t = std::get<0>(tuple);
		Vector* n = std::get<1>(tuple);
		Vector color = std::get<2>(tuple);

		if (t != -1 && t < current_t)
		{
			current_t = t;
			current_Object = o;
			current_color = color;

			if (current_n != nullptr)
			{
				delete current_n;
			}
			current_n = n;
		}
		else
		{
			if (n != nullptr)
			{
				delete n;
			}
		}
	}

	if (current_Object != nullptr)
	{
		r.reverse_animations(current_Object->get_animations(), time);
		intersection = new Vector(r.getPoint(current_t));
		r.apply_animations(current_Object->get_animations(), time);

		Ray r2 = Ray(*intersection, *current_n);
		r2.apply_animations(current_Object->get_animations(), time);
		*intersection = r2.get_center();
		*current_n = r2.get_dir();
	}
	return std::tuple<Vector*, Vector*, IObject*, Vector>(intersection, current_n, current_Object, current_color);
}

double Scene::computeFresnel(Ray r, double n_air, double n_object, Vector* n) const
{
	double k0, R, T;
	Vector i;

	k0 = pow((n_air - n_object) / (n_air + n_object), 2.);

	if (r.get_dir().scalar(*n) < 0)
	{
		i = r.get_dir() * (-1);
	}
	else
	{

		Vector n2 = *n * (-1);

		i = r.get_dir() * (n_object / n_air) - n2 * (n_object / n_air * r.get_dir().scalar(n2) + sqrt(1 - n_object * n_object / (n_air * n_air) * (1 - r.get_dir().scalar(n2) * r.get_dir().scalar(n2))));
	}

	i.normalize();
	R = k0 + (1.0 - k0) * pow((1 - i.scalar(*n)), 5.);
	T = 1. - R;
	return T;
}

bool Scene::computeShadows(Vector nudged, Vector center, double time) const
{
	bool b = false;

#ifdef SHADOWS
	Vector dir = center - nudged;

	Ray r = Ray(nudged, dir);
	r.normalize();

	std::tuple<Vector*, Vector*, IObject*, Vector> tuple = computeIntersection(r, time);
	Vector* intersection2 = std::get<0>(tuple);
	Vector* n2 = std::get<1>(tuple);

	if (intersection2 != nullptr)
	{
		Vector dir2 = *intersection2 - nudged;

		if (dir2.get_norm2() + _epsilon < dir.get_norm2())
		{
			b = true;
		}
	}
	delete n2;
	delete intersection2;
#endif

	return b;
}

Vector Scene::computeIntensityMirror(Ray r, Vector* intersection, IObject* o, Vector* n, Vector color, int nombre_rebonds, double time, bool showEmissiveSurfaces, bool fromTransparent) const
{
	Vector intensity = Vector(0., 0., 0.);

	if (o->get_material()._mirror || fromTransparent)
	{
		Vector nudged2 = *intersection + *n * _epsilon;
		Vector dir2 = r.get_dir() - (*n * 2 * r.get_dir().scalar(*n));

		Ray r2 = Ray(nudged2, dir2);

		std::tuple<Vector*, Vector*, IObject*, Vector> tuple2 = computeIntersection(r2, time);
		Vector* intersection2 = std::get<0>(tuple2);
		Vector* n2 = std::get<1>(tuple2);
		IObject* o2 = std::get<2>(tuple2);
		Vector color2 = std::get<3>(tuple2);

		if (intersection2 != nullptr)
		{
#ifdef RUSSIAN_ROULETTE
			double rand = distrib(engine);
			if (rand < 0.9)
				intensity = computeIntensity(r2, intersection2, o2, n2, color2, nombre_rebonds - 1, time) * o->get_material()._specular_color / 0.9;
#else
			intensity = computeIntensity(r2, intersection2, o2, n2, color2, nombre_rebonds - 1, time) * o->get_material()._specular_color;
#endif
		}
		delete n2;
		delete intersection2;
	}
	return intensity;
}

Vector Scene::computeIntensityTransparent(Ray r, Vector* intersection, IObject* o, Vector* n, Vector color, int nombre_rebonds, double time, bool showEmissiveSurfaces) const
{
	Vector intensity = Vector(0, 0, 0);

	if (o->get_material()._transparent)
	{
		double n_object = o->get_material()._n_object;
		double T = 1.0;
		double rand = distrib(engine);

#ifdef FRESNEL
			T = computeFresnel(r, _n_air, n_object, n);
#endif

		if (rand < T)
		{
			double n_1 = _n_air;
			double n_2 = n_object;
			Vector* normal = n;

			if (r.get_dir().scalar(*n) >= 0)
			{
				n_1 = n_object;
				n_2 = _n_air;
				*normal = *normal * (-1.);
			}

			if (1 - n_1 * n_1 / (n_2 * n_2) * (1 - r.get_dir().scalar(*n) * r.get_dir().scalar(*n)) >= 0)
			{
				Vector* intersection2 = nullptr;
				Vector* n2 = nullptr;
				IObject* o2 = nullptr;
				Vector color2;
				Ray r2;

				Vector dir2 = r.get_dir() * (n_1 / n_2) - *normal * (n_1 / n_2 * r.get_dir().scalar(*normal) + sqrt(1 - n_1 * n_1 / (n_2 * n_2) * (1 - r.get_dir().scalar(*normal) * r.get_dir().scalar(*normal))));
				Vector nudged2 = *intersection - *normal * _epsilon;

				r2 = Ray(nudged2, dir2);

				std::tuple<Vector*, Vector*, IObject*, Vector> tuple3 = computeIntersection(r2, time);
				intersection2 = std::get<0>(tuple3);
				n2 = std::get<1>(tuple3);
				o2 = std::get<2>(tuple3);
				color2 = std::get<3>(tuple3);

				if (intersection2 != nullptr)
				{
#ifdef RUSSIAN_ROULETTE
					double rand2 = distrib(engine);
					if (rand2 < 0.9)
						intensity = computeIntensity(r2, intersection2, o2, n2, color2, nombre_rebonds - 1, time, showEmissiveSurfaces) / 0.9;
#else
					intensity = computeIntensity(r2, intersection2, o2, n2, color2, nombre_rebonds - 1, time, showEmissiveSurfaces);
#endif
				}
				delete n2;
				delete intersection2;
			}
			else
			{
#ifdef RUSSIAN_ROULETTE
				double rand2 = distrib(engine);
				double proba = 0.9;
				if (rand2 < proba)
					intensity = computeIntensityMirror(r, intersection, o, n, color, nombre_rebonds - 1, time, showEmissiveSurfaces, true) / 0.9;
#else
				intensity = computeIntensityMirror(r, intersection, o, n, color, nombre_rebonds - 1, time, showEmissiveSurfaces, true);
#endif
			}
		}
		else
		{
#ifdef RUSSIAN_ROULETTE
			double rand = distrib(engine);
			double proba = 0.9;
			if (rand < proba)
				intensity = computeIntensityMirror(r, intersection, o, n, color, nombre_rebonds - 1, time, showEmissiveSurfaces, true) / 0.9;
#else
			intensity = computeIntensityMirror(r, intersection, o, n, color, nombre_rebonds - 1, time, showEmissiveSurfaces, true);
#endif
		}
	}
	return intensity;
}

Vector Scene::computeIntensityIndirect(Ray r, Vector * intersection, IObject * o, Vector * n, Vector color, int nombre_rebonds, double time, bool showEmissiveSurfaces) const
{
	Vector intensity = Vector(0., 0., 0.);

#ifdef INDIRECT_LIGHT
	Vector* intersection2 = nullptr;
	Vector* n2 = nullptr;
	IObject* o2 = nullptr;
	Vector color2;

	Vector nudged2 = *intersection + *n * _epsilon;

	Ray r2;

	double p = o->get_material()._phong ? 0.5 : 1.;

	Vector reflected = r.get_dir() - (*n * 2 * r.get_dir().scalar(*n));
	double exponent = o->get_material()._phong_exponent;

	double rand = distrib(engine);

	if (o->get_material()._phong == true && rand > p)
	{
		r2 = Ray::getRandRayPhong(*intersection, reflected, exponent);
		if (r2.get_dir().scalar(*n) <= 0.)
			return Vector(0., 0., 0.);
		if (r2.get_dir().scalar(reflected) <= 0.)
			return Vector(0., 0., 0.);
	}
	else
	{
		r2 = Ray::getRandRay(nudged2, *n);
	}

	std::tuple<Vector*, Vector*, IObject*, Vector> tuple2 = computeIntersection(r2, time);
	intersection2 = std::get<0>(tuple2);
	n2 = std::get<1>(tuple2);
	o2 = std::get<2>(tuple2);
	color2 = std::get<3>(tuple2);

	Vector indirect_intensity = Vector(0, 0, 0);

	if (intersection2 != nullptr)
		{
#ifdef RUSSIAN_ROULETTE
		double rand2 = distrib(engine);
		double proba = (color.get_x() + color.get_y() + color.get_z()) / 3.5;
		if (rand2 < proba)
			indirect_intensity = computeIntensity(r2, intersection2, o2, n2, color2, nombre_rebonds - 1, time, false) / proba;
#else
		indirect_intensity = computeIntensity(r2, intersection2, o2, n2, color2, nombre_rebonds - 1, time, false);
#endif
		
		double proba_phong = (exponent + 1.) * std::pow((r2.get_dir().scalar(reflected)), exponent);
		double proba_diffuse = n->scalar(r2.get_dir());

		double proba = p * proba_diffuse + (1. - p) * proba_phong;

		if (o->get_material()._phong == true && rand > p)
			intensity += indirect_intensity * n->scalar(r2.get_dir()) * (exponent + 2.) * std::pow((r2.get_dir().scalar(reflected)), exponent) * o->get_material()._specular_color / (2. * M_PI) / proba;
		else
			intensity += indirect_intensity * n->scalar(r2.get_dir()) * color / (2. * M_PI) / proba;

	}

	delete n2;
	delete intersection2;
#endif

	return intensity;
}

Vector Scene::computeIntensityEmissive(Ray r, Vector * intersection, IObject * o, Vector * n, Vector color, int nombre_rebonds, double time, bool showEmissiveSurfaces) const
{
	Vector intensity = Vector(0., 0., 0.);

#ifdef EMISSIVE_LIGHT
	if (o->get_material()._emissive && showEmissiveSurfaces)
	{
		intensity += o->get_material()._color * o->get_material()._emissivity;
	}
#endif

	return intensity;
}

Vector Scene::computeIntensityDirect(Ray r, Vector * intersection, IObject * o, Vector * n, Vector color, int nombre_rebonds, double time, bool showEmissiveSurfaces) const
{
	Vector intensity = Vector(0., 0., 0.);

#ifdef DIRECT_LIGHT
	Ray r2;

	// We aim one of the emissive object (with chances proportional to total light intensity of the object)

	double rand = distrib(engine);
	double sum = 0;

	std::vector<double> proba;

	for (int i = 0; i < _lights_objects.size(); i++)
	{
		proba.push_back(_lights_objects[i]->get_material()._emissivity / _lights_objects[i]->get_surfaceArea());
		sum += proba[i];
	}

	for (int i = 0; i < _lights_objects.size(); i++)
	{
		if (rand*sum <= proba[i])
		{
			// We get a random direction (giving us a point on the Sphere)

			IAnimatable* light_object_i = _lights_objects[i];

			Vector light_center = light_object_i->get_center();
			Vector dir_center_light = (*intersection - light_center).normalize();

			r2 = Ray::getRandRayAngleUniform(light_center, _lights_objects[i]->get_surfaceArea(), dir_center_light);

			Vector randResultPoint = r2.get_center();

			Vector randResultDir = r2.get_dir().normalize();
			Vector randResultDirtoIntersection = (*intersection - randResultPoint).normalize();

			double d = (*intersection - randResultPoint).get_norm2();

			Vector nudged = *intersection + *n * _epsilon;

			if (computeShadows(nudged, randResultPoint, time) == false)
			{
				if (o->get_material()._phong == true)
				{

					Vector reflected = r.get_dir() - (*n * 2 * r.get_dir().scalar(*n));

					double exponent = o->get_material()._phong_exponent;

					double Phong = pow(reflected.scalar(randResultDirtoIntersection), exponent)*(exponent + 2.) / (2.*M_PI);

					Vector color_phonged = (color + o->get_material()._specular_color * (Phong - 1));
					
					/*if(n->scalar(randResultDirtoIntersection * (-1)) <= 0.)
						std::cout << "" << n->scalar(randResultDirtoIntersection * (-1)) << std::endl;*/

					intensity = _lights_objects[i]->get_material()._color * std::max(0., n->scalar(randResultDirtoIntersection * (-1))) * _lights_objects[i]->get_material()._emissivity
						* _lights_objects[i]->get_surfaceArea() * std::max(0., randResultDir.scalar(randResultDirtoIntersection))
						 * color_phonged / (std::max(0., (dir_center_light.scalar(randResultDir))) * d * 4 * M_PI);

				}
				else
				{
					intensity = _lights_objects[i]->get_material()._color * std::max(0., n->scalar(randResultDirtoIntersection * (-1))) * _lights_objects[i]->get_material()._emissivity
						* _lights_objects[i]->get_surfaceArea() * std::max(0., randResultDir.scalar(randResultDirtoIntersection))
						* color / (std::max(0., (dir_center_light.scalar(randResultDir))) * d * 4 * M_PI);
				}
			}
			break;
		}
	}
#endif
	return intensity;
}

Vector Scene::computeIntensityPointLight(Ray r, Vector * intersection, IObject * o, Vector * n, Vector color, int nombre_rebonds, double time, bool showEmissiveSurfaces) const
{
	Vector intensity = Vector(0., 0., 0.);

#ifdef POINT_LIGHT
	for (int k = 0; k < _lights.size(); k++)
	{
		Light* light = _lights[k];

		double d = (*intersection - light->get_center(time)).get_norm2();
		Vector l = (light->get_center(time) - *intersection).normalize();

		Vector nudged = *intersection + *n * _epsilon;
		if (computeShadows(nudged, light->get_center(time), time) == false)
		{
#ifdef INTENSITY
			intensity += light->get_intensity() * color * std::max(0., l.scalar(*n)) / (d * M_PI);
#else
			intensity += light->get_intensity() * _color / d;
#endif
		}
	}
#endif
	return intensity;
}

Vector Scene::computeIntensity(Ray r, Vector* intersection, IObject* o, Vector* n, Vector color, int nombre_rebonds, double time, bool showEmissiveSurfaces) const
{
	Vector intensity = Vector(0.,0.,0.);

#ifndef RUSSIAN_ROULETTE
	if (nombre_rebonds == 0)
		return Vector(0, 0, 0);
#endif

	Vector intensity_mirror = computeIntensityMirror(r, intersection, o, n, color, nombre_rebonds, time, showEmissiveSurfaces).max(0.);
	Vector intensity_transparent = computeIntensityTransparent(r, intersection, o, n, color, nombre_rebonds, time, showEmissiveSurfaces).max(0.);
	Vector intensity_indirect = computeIntensityIndirect(r, intersection, o, n, color, nombre_rebonds, time, showEmissiveSurfaces).max(0.);
	Vector intensity_direct = computeIntensityDirect(r, intersection, o, n, color, nombre_rebonds, time, showEmissiveSurfaces).max(0.);
	Vector intensity_emissive = computeIntensityEmissive(r, intersection, o, n, color, nombre_rebonds, time, showEmissiveSurfaces).max(0.);
	Vector intensity_point_light = computeIntensityPointLight(r, intersection, o, n, color, nombre_rebonds, time, showEmissiveSurfaces).max(0.);

	intensity = intensity_mirror + intensity_transparent + intensity_indirect + intensity_direct + intensity_emissive + intensity_point_light;

	return intensity;
}
