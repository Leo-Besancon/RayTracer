#include "stdafx.h"
#include "RayTracer.h"

RayTracer::RayTracer()
{
	_total_time = 100.;
	_nb_frames = 1;
	_nb_rays = 300;
	_nb_bounces = 3;
	_nb_box = 2;

#ifdef MOTION_BLUR
	_exposition_time_per_frame = 1.;
#else
	_exposition_time_per_frame = 0.;
#endif

#ifdef GAMMA
	_gamma = 2.2;
#else
	_gamma = 1.0;
#endif
}

RayTracer::~RayTracer()
{
	delete _scene;
	delete _camera;
}

void RayTracer::createCamera()
{
	if (_camera != nullptr)
		delete _camera;

	int H = 300;
	int W = 300;
	_camera = new Camera(Vector(0, 0, 55), Vector(0, 0, -1), Vector(0, 1, 0), 60., 35., H, W);

	_img = std::vector<unsigned char>(W*H * 3, 0);
}

void RayTracer::createScene()
{
#ifdef MAIN_LOGS
	std::cout << "     -- Start Scene creation -- " << std::endl;
#endif

	if (_scene != nullptr)
		delete _scene;

	_scene = new Scene(0.001, 1.);
	
	// Create the Objects and Lights

	Light* l1 = new Light(Vector(-30, 5, 50), Vector(200000000, 200000000, 200000000) * 0.4);
	Light* l2 = new Light(Vector(10, 20, 40), Vector(200000000, 200000000, 200000000) * 0.4);

	double sphere_light_radius = 10.;
	Sphere* slight = new Sphere(Material::createEmissive(Vector(1, 1, 1), 2000000000.0 / (4.0 * M_PI * sphere_light_radius * sphere_light_radius)), Vector(-30, 5, 50), sphere_light_radius);

	Sphere* s1 = new Sphere(Material::createPhong(Vector(1., 1., 1.), Vector(0.2,0.2,0.2), 20), Vector(0, 5, 15), 15);
	//Sphere* s1 = new Sphere(Material::createMirror(Vector(1., 1., 1.)), Vector(0, 5, 15), 15);
	//Sphere* s1 = new Sphere(Material::createDiffuse(Vector(1., 1., 1.)), Vector(0, 5, 15), 15);
	Sphere* s2 = new Sphere(Material::createDiffuse(Vector(1, 0, 0)), Vector(0, 1000, 0), 940);
	Sphere* s3 = new Sphere(Material::createDiffuse(Vector(0, 1, 0)), Vector(0, 0, -1000), 940);
	Sphere* s4 = new Sphere(Material::createDiffuse(Vector(1, 0, 1)), Vector(0, 0, 1000), 940);
	Sphere* s5 = new Sphere(Material::createDiffuse(Vector(0, 0, 1)), Vector(0, -1000, 0), 990);

	Sphere* s6 = new Sphere(Material::createDiffuse(Vector(1, 1, 0.5)), Vector(0, 10, 35), 5);


#ifdef MAIN_LOGS
	std::cout << "         -- Start Mesh creation -- " << std::endl;
#endif
	Mesh* m = new Mesh(Material(), "T-Rex Model", "Texture.bmp");
	//Mesh* m = new Mesh(Material::createMirror(Vector(1,1,1)), "T-Rex Model", "Texture.bmp");
#ifdef MAIN_LOGS
	std::cout << "         -- End Mesh creation -- " << std::endl;
#endif

	m->normalize();
	m->scale(20);

	m->translate(Vector(0, 5, 15));
	
#ifdef MAIN_LOGS
	std::cout << "         -- Start Compute Bounding Boxes-- " << std::endl;
#endif
	m->computeBoxes(_nb_box);
#ifdef MAIN_LOGS
	std::cout << "         -- End Compute Bounding Boxes-- " << std::endl;
#endif

	//_camera->addAnimation(Animation(-0.01, 4., Vector(0, 0, 0), 1., 0., Vector(0., 0., 0.), 60., Vector(0, 5, 15), 0., Vector(0, 0, 0)));
	m->addAnimation(Animation(-.01, 30., Vector(0, 0, 0), 1., 0., Vector(0., 0., 0.), 360, Vector(0, 5, 15), 0., Vector(0, 0, 0)));

	//scene->addLight(l1);
	//scene->addLight(l2);
	_scene->addObject(slight);
	_scene->addLightObject(slight);

	_scene->addObject(s1);
	_scene->addObject(s2);
	_scene->addObject(s3);
	_scene->addObject(s4);
	_scene->addObject(s5);

	//scene->addObject(s6);
	//_scene->addObject(m);

#ifdef MAIN_LOGS
	std::cout << "     -- End Scene creation -- " << std::endl;
#endif
}

void RayTracer::render_all_frames()
{
	if (_camera == nullptr || _scene == nullptr)
	{
#ifdef MAIN_LOGS
		std::cout << "     -- ERROR : RayTracer needs both a Camera and a Scene ! -- " << std::endl;
#endif
		return;
	}

#ifdef MAIN_LOGS
	std::cout << "     -- Start Image Loop -- " << std::endl;
#endif

	for (int frame = 0; frame < _nb_frames; frame++)
	{
		render_one_frame(frame);
	}

#ifdef MAIN_LOGS
	std::cout << "     -- End Image Loop -- " << std::endl;
#endif
}

void RayTracer::render_one_frame(int frame)
{

#ifdef MAIN_LOGS
	std::cout << "          -- Start Image " << frame + 1 << " / " << _nb_frames << " -- " << std::endl;
	std::cout << "\r                Approx. progress : 0 %";
	int cur_iter = 0;
#endif

#ifdef MT	
#pragma omp parallel for schedule(dynamic,1)
#endif
	for (int i = 0; i < _camera->get_H(); i++)
	{
		for (int j = 0; j < _camera->get_W(); j++)
		{
			Vector intensity = Vector(0, 0, 0);

			Vector* intersection = nullptr;
			Vector* n = nullptr;
			IObject* o = nullptr;
			Vector color;

#ifdef MONTE_CARLO
			for (int N = 0; N < _nb_rays; N++)
			{
				double time = _total_time * frame * (1. / _nb_frames) + _exposition_time_per_frame * distrib(engine);

				Ray r;

#ifdef DoF
				r = Ray::getAADoFRay(i, j, _camera->get_H(), _camera->get_W(), _camera->get_depth(), _camera->get_center(), _camera->get_focal(), _camera->get_up(), _camera->get_direction());
#else
#ifdef AA
				r = Ray::getAARay(i, j, _camera->get_H(), _camera->get_W(), _camera->get_depth(), _camera->get_center(), _camera->get_up(), _camera->get_direction());
#else
				r = Ray::getNonAARay(i, j, _camera->get_H(), _camera->get_W(), _camera->get_depth(), _camera->get_center(), _camera->get_up(), _camera->get_direction());
#endif
#endif
				r.normalize();

				r.apply_animations(_camera->get_animations(), time);

				std::tuple<Vector*, Vector*, IObject*, Vector> tuple = _scene->computeIntersection(r, time);
				intersection = std::get<0>(tuple);
				n = std::get<1>(tuple);
				o = std::get<2>(tuple);
				color = std::get<3>(tuple);

				if (o != nullptr)
				{
					Vector intensity_this_ray = _scene->computeIntensity(r, intersection, o, n, color, _nb_bounces, time);
					intensity = intensity + intensity_this_ray;
				}
				delete n;
				delete intersection;
			}

			intensity = intensity * (1. / _nb_rays);
#else
			r = Ray::getNonAARay(i, j, H, W, scene->get_camera()->get_depth(), scene->get_camera()->get_center(), scene->get_camera()->get_up(), scene->get_camera()->get_direction());
			r.normalize();

			r.apply_animations(_camera->get_animations(), time);

			std::tuple<Vector*, Vector*, IObject*, Vector> tuple = scene->computeIntersection(r, time);
			intersection = std::get<0>(tuple);
			n = std::get<1>(tuple);
			o = std::get<2>(tuple);
			color = std::get<3>(tuple);

			intensity = scene->computeIntensity(r, intersection, o, n, color, _nb_bounces, time);
			delete n;
			delete intersection;
#endif

#ifndef COLOR
			double x = intensity.get_x();
			double y = intensity.get_y();
			double z = intensity.get_z();
			double moy = (x + y + z) / 3;
			intensity = Vector(moy,moy,moy);

#endif

#ifdef INTENSITY
			_img[((_camera->get_H() - 1 - i)* _camera->get_W() + j) * 3] = (unsigned char)std::min(255., pow(intensity.get_x(), 1. / _gamma));
			_img[((_camera->get_H() - 1 - i)* _camera->get_W() + j) * 3 + 1] = (unsigned char)std::min(255., pow(intensity.get_y(), 1. / _gamma));
			_img[((_camera->get_H() - 1 - i)* _camera->get_W() + j) * 3 + 2] = (unsigned char)std::min(255., pow(intensity.get_z(), 1. / _gamma));
#else
			if (std::min(255., pow(intensity.get_x(), 1. / _gamma)) >= 1.)
				img[((_camera->get_H() - 1 - i)* _camera->get_W() + j) * 3] = 255.;
			if (std::min(255., pow(intensity.get_y(), 1. / _gamma)) >= 1.)
				img[((_camera->get_H() - 1 - i)* _camera->get_W() + j) * 3 + 1] = 255.;
			if (std::min(255., pow(intensity.get_z(), 1. / _gamma)) >= 1.)
				img[((_camera->get_H() - 1 - i)* _camera->get_W() + j) * 3 + 2] = 255.;
#endif
		}
#ifdef MAIN_LOGS
#pragma omp critical
		{
			std::cout << "\r                Approx. progress : " << (int)(cur_iter++ * 100 / _camera->get_H()) << " % ";

		}
#endif
	}
#ifdef MAIN_LOGS
	std::cout << "\r                Approx. progress : 100 %" << std::endl;
	std::cout << "          -- End Image " << frame + 1 << " / " << _nb_frames << " -- " << std::endl;
#endif

	// Saving the resulting BMP
	std::ostringstream out_filename("");
	out_filename << "render_" << frame << ".bmp";
	save_image(out_filename.str().c_str(), &_img[0], _camera->get_W(), _camera->get_H());
}

void RayTracer::save_image(const char* filename, const unsigned char* tableau, int w, int h) {

	FILE *f;

	int filesize = 54 + 3 * w*h;

	unsigned char bmpfileheader[14] = { 'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0 };
	unsigned char bmpinfoheader[40] = { 40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0 };
	unsigned char bmppad[3] = { 0,0,0 };

	bmpfileheader[2] = (unsigned char)(filesize);
	bmpfileheader[3] = (unsigned char)(filesize >> 8);
	bmpfileheader[4] = (unsigned char)(filesize >> 16);
	bmpfileheader[5] = (unsigned char)(filesize >> 24);

	bmpinfoheader[4] = (unsigned char)(w);
	bmpinfoheader[5] = (unsigned char)(w >> 8);
	bmpinfoheader[6] = (unsigned char)(w >> 16);
	bmpinfoheader[7] = (unsigned char)(w >> 24);
	bmpinfoheader[8] = (unsigned char)(h);
	bmpinfoheader[9] = (unsigned char)(h >> 8);
	bmpinfoheader[10] = (unsigned char)(h >> 16);
	bmpinfoheader[11] = (unsigned char)(h >> 24);

	f = fopen(filename, "wb");
	fwrite(bmpfileheader, 1, 14, f);
	fwrite(bmpinfoheader, 1, 40, f);
	unsigned char *row = new unsigned char[w * 3];
	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++) {
			row[j * 3] = tableau[(w*(h - i - 1) * 3) + j * 3 + 2];
			row[j * 3 + 1] = tableau[(w*(h - i - 1) * 3) + j * 3 + 1];
			row[j * 3 + 2] = tableau[(w*(h - i - 1) * 3) + j * 3];
		}
		fwrite(row, 3, w, f);
		fwrite(bmppad, 1, (4 - (w * 3) % 4) % 4, f);
	}
	fclose(f);
	delete[] row;
}
