#include "stdafx.h"
#include "IAnimatable.h"
#include "Ray.h"

Vector IAnimatable::get_center(double time) const
{
	Ray ret = Ray(get_center(), Vector(1, 0, 0));

	ret.apply_animations(_animations, time);

	return ret.get_center();
}