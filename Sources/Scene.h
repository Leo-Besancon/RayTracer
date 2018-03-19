#pragma once

#include "Light.h"
#include "IObject.h"

class Scene
{
	// Lights
	std::vector<Light*> _lights;
	std::vector<IObject*> _lights_objects;

	// Objects
	std::vector<IObject*> _objects;

	double _epsilon;
	double _n_air;

	double computeFresnel(Ray r, double n_air, double n_object, Vector * n) const;
	bool computeShadows(Vector nudged, Vector center, double time) const;
	Vector computeIntensityMirror(Ray r, Vector * intersection, IObject * o, Vector * n, Vector color, int nombre_rebonds, double time, bool showEmissiveSurfaces, bool fromTransparent = false) const;
	Vector computeIntensityTransparent(Ray r, Vector * intersection, IObject * o, Vector * n, Vector color, int nombre_rebonds, double time, bool showEmissiveSurfaces) const;
	Vector computeIntensityIndirect(Ray r, Vector * intersection, IObject * o, Vector * n, Vector color, int nombre_rebonds, double time, bool showEmissiveSurfaces) const;
	Vector computeIntensityEmissive(Ray r, Vector * intersection, IObject * o, Vector * n, Vector color, int nombre_rebonds, double time, bool showEmissiveSurfaces) const;
	Vector computeIntensityDirect(Ray r, Vector * intersection, IObject * o, Vector * n, Vector color, int nombre_rebonds, double time, bool showEmissiveSurfaces) const;
	Vector computeIntensityPointLight(Ray r, Vector * intersection, IObject * o, Vector * n, Vector color, int nombre_rebonds, double time, bool showEmissiveSurfaces) const;

public:
	Scene(double epsilon, double n_air) : _epsilon(epsilon), _n_air(n_air) { };
	~Scene();

	double get_epsilon() const { return _epsilon; };
	void set_epsilon(double epsilon) { _epsilon = epsilon; };

	double get_n_air() const { return _n_air; };
	void set_n_air(double n_air) { _n_air = n_air; };

	std::vector<Light*> get_lights() const { return _lights; };
	std::vector<IObject*> get_light_objects() const { return _lights_objects; };
	std::vector<IObject*> get_objects() const { return _objects; };

	void addLight(Light* l) { _lights.push_back(l); };
	void addLightObject(IObject* o) { _lights_objects.push_back(o); };
	void addObject(IObject* o) { _objects.push_back(o); };

	std::tuple<Vector*, Vector*, IObject*, Vector> computeIntersection(Ray r, double time) const;
	Vector computeIntensity(Ray r, Vector* intersection, IObject* o, Vector* n, Vector color, int nombre_rebonds, double time, bool showEmissiveSurfaces = true) const;

};

