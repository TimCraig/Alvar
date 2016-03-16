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

#include "Pose.h"

using namespace std;

namespace alvar {

using namespace std;

/*****************************************************************************
 ***  class Pose
 ****************************************************************************/

Pose::Pose() : rotation(), translation_mat(4, 1)
   {
   translation_mat = 0.0;
   translation_mat(3, 0) =  1.0;

   return;

   }

Pose::Pose(cv::Mat& tra, cv::Mat_<double>& rot, Rotation::RotationType t) : rotation(rot, t), translation_mat(4, 1)
   {
   // Homogeneous
   translation_mat(3, 0) = 1.0;

   // Fill in translation part
   translation_mat(0, 0) = tra.at<double>(0, 0);
   translation_mat(1, 0) = tra.at<double>(1, 0);
   translation_mat(2, 0) = tra.at<double>(2, 0);

   return;

   }

Pose::Pose(cv::Mat_<double>& mat) : rotation(mat, Rotation::RotationType::MAT), translation_mat(4, 1)
   {
   translation_mat = 0.0;
   translation_mat(3, 0) = 1.0;

   // Fill in translation part
   if (mat.cols == 4)
      {
      translation_mat(0, 0) = mat.at<double>(0, 3);
      translation_mat(1, 0) = mat.at<double>(1, 3);
      translation_mat(2, 0) = mat.at<double>(2, 3);
      }

   return;

   }

Pose::Pose(const Pose& src) : rotation(src.rotation)
   {
   translation_mat = src.translation_mat.clone();

   return;

   }

Pose& Pose::operator = (const Pose& rhs)
   {
   if (this != &rhs)
      {
      rotation = rhs.rotation;
      translation_mat = rhs.translation_mat.clone();
      } // end if

   return (*this);

   }

void Pose::Reset()
   {
   rotation.Reset();
   translation_mat = 0.0;
   translation_mat(3, 0) = 1.0;

   return;

   }

void Pose::SetMatrix(const cv::Mat& mat)
   {
   cv::Mat_<double> tmp;
   for(int i = 0 ; i < 3 ; ++i)
      for(int j = 0 ; j < 3 ; ++j)
         tmp(i, j) = mat.at<double>(i, j);

   rotation.FromMat9(tmp);
   if (mat.cols == 4)
      {
      translation_mat(0, 0) = mat.at<double>(0, 3);
      translation_mat(1, 0) = mat.at<double>(1, 3);
      translation_mat(2, 0) = mat.at<double>(2, 3);
      translation_mat(3, 0) = 1.0;
      }

   return;

   }

cv::Mat_<double> Pose::GetHomogeneousMatrix() const
   {
   cv::Mat_<double> mat = rotation.GetHomogeneousMatrix();

   mat(0, 3) = translation_mat(0, 0);
   mat(1, 3) = translation_mat(1, 0);
   mat(2, 3) = translation_mat(2, 0);


   return (mat);

   }

// Need to fix this!  Not currently using TTC!!
void Pose::GetMatrixGL(double /* gl */[16], bool /* mirror */)
   {
#if 0
   if (mirror)
      Mirror(false, true, true);

   cv::Mat gl_mat(4, 4, CV_64F, gl);
   GetMatrix(gl_mat);
   cv::Mat T = gl_mat.t();
   T.copyTo(gl_mat);

   if (mirror)
      Mirror(false, true, true);
#endif
   return;

   }

void Pose::SetMatrixGL(double gl[16], bool mirror)
{
	double gll[16];
   memcpy(gll, gl, sizeof(double) * 16);
   cv::Mat gl_mat(4, 4, CV_64F, gll);
   cv::Mat T = gl_mat.t();
   T.copyTo(gl_mat);
   SetMatrix(gl_mat);

   if (mirror)
      Mirror(false, true, true);

   return;

}

void Pose::Transpose()
   {
   SetMatrix(GetHomogeneousMatrix().t());

   return;

   }

void Pose::Invert()
   {
   cv::Mat tmp_mat = GetHomogeneousMatrix();
   tmp_mat = tmp_mat.inv();
   SetMatrix(tmp_mat);

   return;

   }

void Pose::Mirror(bool x, bool y, bool z)
   {
   cv::Mat tmp_mat = GetRotationMatrix();
   Rotation::MirrorMat(tmp_mat, x, y, z);
   SetMatrix(tmp_mat);

   return;

   }

void Pose::SetTranslation(const cv::Mat& tra)
   {
   translation_mat(0, 0) = tra.at<double>(0, 0);
   translation_mat(1, 0) = tra.at<double>(1, 0);
   translation_mat(2, 0) = tra.at<double>(2, 0);
   translation_mat(3, 0) = 1.0;

   return;

   }

void Pose::SetTranslation(const double* tra)
   {
   translation_mat(0, 0) = tra[0];
   translation_mat(1, 0) = tra[1];
   translation_mat(2, 0) = tra[2];
   translation_mat(3, 0) = 1.0;

   return;

   }

void Pose::SetTranslation(const double x, const double y, const double z)
   {
   translation_mat(0, 0) = x;
   translation_mat(1, 0) = y;
   translation_mat(2, 0) = z;
   translation_mat(3, 0) = 1.0;

   return;

   }

void Pose::GetTranslation(cv::Mat& tra) const
   {
   tra.at<double>(0, 0) = translation_mat(0, 0);
   tra.at<double>(1, 0) = translation_mat(1, 0);
   tra.at<double>(2, 0) = translation_mat(2, 0);
   if (tra.rows == 4)
      tra.at<double>(3, 0) = 1.0;

   return;

   }

void Pose::Output() const
   {
   rotation.Output();

   cout << translation_mat.at<double>(0, 0) << ","
        << translation_mat.at<double>(1, 0)<< ","
        << translation_mat.at<double>(2, 0)
        << endl;

   return;

   }

} // namespace alvar
