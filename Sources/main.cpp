// main.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "RayTracer.h"
#include "Sphere.h"
#include "Mesh.h"

int main()
{
#ifdef _WIN32
#ifdef LOW_PRIORITY
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
#endif

#ifdef MAIN_LOGS
	std::cout << " -- Start Program -- " << std::endl;
#endif

	RayTracer* raytracer = new RayTracer();

	raytracer->createCamera();

	raytracer->createScene();

	raytracer->render_all_frames();

	delete raytracer;

#ifdef MAIN_LOGS
	std::cout << " -- End Program -- " << std::endl;
#endif

	return 0;
}
