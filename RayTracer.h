#pragma once

#include "Scene.h"
#include "Camera.h"
#include "Sphere.h"
#include "Mesh.h"

class RayTracer
{
	Scene* _scene;
	Camera* _camera;

	double _total_time;
	double _exposition_time_per_frame;
	int _nb_frames;
	int _nb_rays;
	int _nb_bounces;
	int _nb_box;
	double _gamma;


	std::vector<unsigned char> _img;

	void render_one_frame(int frame);
	void addAnimation(IAnimatable* animatable, Animation animation) { animatable->addAnimation(animation); };
	void save_image(const char * filename, const unsigned char * tableau, int w, int h);

public:
	RayTracer();
	~RayTracer();

	void createCamera();
	void createScene();
	void render_all_frames();
};

