// stdafx.h : fichier Include pour les fichiers Include système standard,
// ou les fichiers Include spécifiques aux projets qui sont utilisés fréquemment,
// et sont rarement modifiés
//

#pragma once

#define _CRT_SECURE_NO_WARNINGS

#include "targetver.h"

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif
#include <stdio.h>
#include <tchar.h>

#include <cmath>
#include <algorithm> 

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <limits>
#include <tuple>
#include <random>

#include "Vector.h"

#define M_PI 3.14159265358979323846

// Config parameters to add / remove functionalities easily

//#define LOGS // Mainly for Mesh parsing
#define MAIN_LOGS
#define LOW_PRIORITY
#define GAMMA
#define COLOR
#define INTENSITY
#define SHADOWS
#define MT
#define MONTE_CARLO
#define AA
#define FRESNEL
//#define DoF
#define INDIRECT_LIGHT
#define DIRECT_LIGHT
#define EMISSIVE_LIGHT
//#define MOTION_BLUR
//#define RUSSIAN_ROULETTE
//#define POINT_LIGHT

static std::default_random_engine engine;
static std::uniform_real_distribution<double> distrib(0, 1);
