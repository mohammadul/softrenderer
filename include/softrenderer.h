/**
 * @file softrenderer.h
 * @author Sk. Mohammadul Haque
 * @version 0.1.0.0
 * @copyright
 * Copyright (c) 2019 Sk. Mohammadul Haque.
 * @brief This header file contains declarations of all functions of softrenderer.
 */

#ifndef SOFTRENDERER_H_INCLUDED
#define SOFTRENDERER_H_INCLUDED


/*! \mainpage SoftRenderer
 *
 * \section intro_sec Introduction
 * SoftRenderer is a simple software-based rendering library written in C++.
 *
 * \section build_sec Build
 * To build the whole project, either Code::Blocks or Visual Studio 2012 (or later) is required.
 *
 * \author Sk. Mohammadul Haque
 * \copyright Copyright (c) 2019 Sk. Mohammadul Haque.
 */

#include <meshlib.h>
#include <vector>

#ifdef _WIN32
#  ifdef SOFTRENDERER_API_EXPORTS
#    define SOFTRENDERER_API __declspec(dllexport)
#  else
#    define SOFTRENDERER_API __declspec(dllimport)
#  endif
#else
#  define SOFTRENDERER_API extern
#endif


SOFTRENDERER_API int sr_mesh_render_persp(MESH m, std::vector<FLOATDATA>& depthmap, std::vector<FLOATDATA>& scalarmap, FLOATDATA f, INTDATA h, INTDATA w, bool doscalar, FLOATDATA cutoff = 0.0, FLOATDATA eps = 0.0001);
SOFTRENDERER_API int sr_mesh_render_ortho(MESH m, std::vector<FLOATDATA>& depthmap, std::vector<FLOATDATA>& scalarmap, FLOATDATA g, INTDATA h, INTDATA w, bool doscalar, FLOATDATA cutoff = 0.0, FLOATDATA eps = 0.0001);

#endif // SOFTRENDERER_H_INCLUDED


