/*
 * This file is part of ALVAR, A Library for Virtual and Augmented Reality.
 *
 * Copyright 2007-2012 VTT Technical Research Centre of Finland
 *
 * Contact: VTT Augmented Reality Team <alvar.info@vtt.fi>
 *          <http://www.vtt.fi/multimedia/alvar.html>
 *
 * ALVAR is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 2.1 of the License, or (at your option)
 * any later version.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with ALVAR; if not, see
 * <http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html>.
 */

#ifndef UTIL_H
#define UTIL_H

/**
 * \file Util.h
 *
 * \brief This file implements generic utility functions and a serialization
 * interface.
 */

// 2-10-16  Started cleaning up this module with better interface, removed some of the crazy inline returns, etc.  Lots to do.  TTC!!

/*****************************************************************************
 ******************************  I N C L U D E  ******************************
 ****************************************************************************/

#include "Alvar.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
//#include <cxcore.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
//#include <cv.h>
#include <cmath> //for abs
#include <map>
#include <type_traits>

namespace alvar {

const double PI = 3.14159265358979323846;

/**
* \brief Returns the sign of a number.
*/
template<class C> inline
int ALVAR_EXPORT Sign(const C& v)
   {
   return ((v < 0) ? -1 : 1);
   }

/**
 * \brief Simple \e Point class meant to be inherited from OpenCV point-classes. For example: Point<cv::Point2d>
 */
template<class C, typename D = int>
struct ALVAR_EXPORT Point : public C
   {
//   static_assert(std::is_base_of<cv::Point_<class T>, C>::value, "Point requires class with base cv::Point_<>");

   public:
      /**
    * \brief Additional value can be related to the point.
    */
      D val;

      Point(int vx = 0, int vy = 0)
         {
         C::x = vx;
         C::y = vy;
         }

      Point(double vx, double vy)
         {
         C::x = vx;
         C::y = vy;
         }
   };

/** 
  *  \brief The default integer point type.
*/
using PointInt = ALVAR_EXPORT Point<cv::Point2i>;

/**
  *  \brief The default double point type.
*/
using PointDouble = ALVAR_EXPORT Point<cv::Point2d>;

/** \brief Returns the squared distance of two points. 
  * \param p1	First point.
  * \param p2	Second point.
  * \return Squared distance.
*/

template<class PointType>
double PointSquaredDistance(PointType p1, PointType p2)
   {
   return ((p1.x - p2.x) * (p1.x - p2.x)) +
         ((p1.y - p2.y) * (p1.y - p2.y));
   }

/** 
  * \brief Limits a number to between two values. 
  * \param val		Input value.
  * \param min_val	Minimum value for the result.
  * \param max_val	Maximum value for the result.
  * \return			Resulting value that is between \e min_val and \e max_val.
*/
template<typename T>
T Limit(T val, T min_val, T max_val)
   {
   return (std::max(min_val, std::min(max_val, val)));
   }

} // namespace alvar

#endif
