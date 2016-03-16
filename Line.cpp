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

#include "Line.h"

using namespace std;

namespace alvar {

using namespace std;

/*****************************************************************************
 ***  class Line
 ****************************************************************************/

PointDouble Line::Intersection(const Line& l1, const Line& l2)
   {
   double vx = l1.s.x;
   double vy = l1.s.y;
   double ux = l2.s.x;
   double uy = l2.s.y;
   double wx = l2.c.x - l1.c.x;
   double wy = l2.c.y - l1.c.y;

   double s =0.0;
   double px = 0.0;
   double py = 0.0;
   double tmp = vx * uy - vy * ux;
   if (tmp == 0.0)
      tmp = 1.0;

   //if(/*tmp <= 1.f && tmp >= -1.f && */tmp != 0.f && ang > 0.1)
   {
   s = (vy * wx - vx * wy) / tmp;
   px = l2.c.x + s * ux;
   py = l2.c.y + s * uy;
   }

   return (PointDouble(px, py));

   }

int Line::FitLines(std::vector<Line> &lines,
             const vector<int>& corners,
             const vector<PointInt>& edge,
             cv::Mat& /* grey */) // grey image for future sub pixel accuracy
   {
   lines.clear();

   for (size_t j = 0 ; j < corners.size() ; ++j)
      {
      int start, end;
      int size = (int)edge.size();

      start = corners[j];

      if (j < (corners.size() - 1))
         {
         end  = corners[j + 1];
         } // end if
      else
         {
         end = corners[0];
         } // end else

      int len = end - start + 1;
      if (start >= end)
         {
         len += size;
         } // end if

      // OpenCV routine...
      // Change to std::vector. TTC!!
      cv::Mat line_data(1, len, CV_32FC2);
      for (int i = 0 ; i < len ; ++i)
         {
         int ind = i + start;
         if (ind >= size)
            ind = ind - size;

         double px = double(edge[ind].x);
         double py = double(edge[ind].y);
         line_data.at<cv::Point2f>(0, i) = cv::Point2f(px, py);
         }

      cv::Vec4f params;
      cv::fitLine(line_data, params, CV_DIST_L2, 0.0, 0.01, 0.01);
      lines.push_back(Line(params));
      }

   return (lines.size());

   }

} // namespace alvar
