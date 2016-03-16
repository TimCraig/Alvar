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

#include "Rotation.h"

#include <opencv2/calib3d/calib3d.hpp>

using namespace std;

namespace alvar {
using namespace std;

/*****************************************************************************
 ***  class Rotation
 ****************************************************************************/

Rotation::Rotation()
   {
	Reset();

   return;

   }

Rotation::Rotation(cv::Mat_<double>& mat, RotationType t)
   {
	Reset();

	switch (t)
      {
      case RotationType::QUAT :
         SetQuaternion(mat);
			break;

      case RotationType::MAT :
         SetMatrix(mat);
			break;

      case RotationType::EUL :
         SetEuler(mat);
			break;

      case RotationType::ROD :
         SetRodrigues(mat);
			break;
      }

   return;

   }

// static
void Rotation::MirrorMat(cv::Mat& mat, bool x, bool y, bool z)
   {
//   cv::Mat mat_mul = cv::Mat::eye(mat.size(), mat.type());
   cv::Mat_<double> mat_mul = cv::Mat::eye(3, 3, CV_64F);

   if (x)
      mat_mul(0, 0) = -1.0;
   if (y)
      mat_mul(1, 1) = -1.0;
   if (z)
      mat_mul(2, 2) = -1.0;

   // Check order is correct  TTC!!
   mat = mat_mul.mul(mat);

   return;

   }
	
void Rotation::Mirror(bool x, bool y, bool z)
   {
   cv::Mat tmp_mat = GetMatrix();
   MirrorMat(tmp_mat, x, y, z);
   SetMatrix(tmp_mat);

   return;

   }

void Rotation::Reset()
   {
   quaternion[0] = 1.0;
   quaternion[1] = 0.0;
   quaternion[2] = 0.0;
   quaternion[3] = 0.0;

   return;

   }

// static
cv::Mat_<double> Rotation::Mat9ToRod(const cv::Mat_<double>& mat)
   {
   assert((mat.rows == 3) && (mat.cols == 3));
   cv::Mat_<double> rod(3, 1);

   cv::Rodrigues(mat, rod);

   return (rod);

   }

// static
cv::Mat_<double> Rotation::RodToMat9(const cv::Mat_<double>& rod)
   {
   assert((rod.rows == 3) && (rod.cols == 1));

   cv::Mat_<double> mat(3, 3);

   cv::Rodrigues(rod, mat);

   return (mat);

   }

// Conjugate the quaternion
void Rotation::Conjugate()
   {
   quaternion[1] = -quaternion[1];
   quaternion[2] = -quaternion[2];
   quaternion[3] = -quaternion[3];

   return;

   }

void Rotation::Normalize()
   {
   double length = sqrt(quaternion[0] * quaternion[0] +
         quaternion[1] * quaternion[1] + quaternion[2] * quaternion[2] +
         quaternion[3] * quaternion[3]);
	
   if (length != 0.0)
      {
      for (size_t i = 0 ; i < 4 ; ++i)
         {
         quaternion[i] = quaternion[i] / length;
         }
      }

   return;

   }

// static
void Rotation::Mul(const Rotation& q1, const Rotation& q2, Rotation& result)
   {
   result[0] = q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3];
   result[1] = q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2];
   result[2] = q1[0] * q2[2] + q1[2] * q2[0] + q1[3] * q2[1] - q1[1] * q2[3];
   result[3] = q1[0] * q2[3] + q1[3] * q2[0] + q1[1] * q2[2] - q1[2] * q2[1];

   result.Normalize();

   return;

   }

cv::Mat_<double> Rotation::ToMat9() const
   {
   cv::Mat_<double> mat(3, 3);

   double xx = quaternion[1] * quaternion[1];
   double xy = quaternion[1] * quaternion[2];
   double xz = quaternion[1] * quaternion[3];
   double xw = quaternion[1] * quaternion[0];

   double yy = quaternion[2] * quaternion[2];
   double yz = quaternion[2] * quaternion[3];
   double yw = quaternion[2] * quaternion[0];

   double zz = quaternion[3] * quaternion[3];
   double zw = quaternion[3] * quaternion[0];

   mat(0, 0) = 1 - 2 * ( yy + zz );
   mat(0, 1) =     2 * ( xy - zw );
   mat(0, 2) =     2 * ( xz + yw );
   mat(1, 0) =     2 * ( xy + zw );
   mat(1, 1) = 1 - 2 * ( xx + zz );
   mat(1, 2) =     2 * ( yz - xw );
   mat(2, 0) =     2 * ( xz - yw );
   mat(2, 1) =     2 * ( yz + xw );
   mat(2, 2) = 1 - 2 * ( xx + yy );

   return (mat);

   }

cv::Mat_<double> Rotation::ToMat16() const
   {
   cv::Mat_<double> mat(4, 4);

   double xx = quaternion[1] * quaternion[1];
   double xy = quaternion[1] * quaternion[2];
   double xz = quaternion[1] * quaternion[3];
   double xw = quaternion[1] * quaternion[0];

   double yy = quaternion[2] * quaternion[2];
   double yz = quaternion[2] * quaternion[3];
   double yw = quaternion[2] * quaternion[0];

   double zz = quaternion[3] * quaternion[3];
   double zw = quaternion[3] * quaternion[0];

   mat(0, 0) = 1 - 2 * ( yy + zz );
   mat(0, 1) =     2 * ( xy - zw );
   mat(0, 2) =     2 * ( xz + yw );
   mat(0, 3) = 0.0;
   mat(1, 0) =     2 * ( xy + zw );
   mat(1, 1) = 1 - 2 * ( xx + zz );
   mat(1, 2) =     2 * ( yz - xw );
   mat(1, 3) = 0.0;
   mat(2, 0) =     2 * ( xz - yw );
   mat(2, 1) =     2 * ( yz + xw );
   mat(2, 2) = 1 - 2 * ( xx + yy );
   mat(2, 3) = 0.0;
   mat(3, 0) = 0.0;
   mat(3, 1) = 0.0;
   mat(3, 2) = 0.0;
   mat(3, 3) = 1.0;

   return (mat);

   }


cv::Mat_<double> Rotation::ToEul() const
   {
   cv::Mat_<double> eul(3, 1);

   // Waste, not going to fix it now  TTC!!
   double qw = quaternion[0];
   double qx = quaternion[1];
   double qy = quaternion[2];
   double qz = quaternion[3];

   double heading = 0.0, bank = 0.0, attitude = 0.0;

   if ((2 * qx * qy + 2 * qz * qw) == 1.0)
      {
      heading = 2 * atan2(qx, qw);
      bank = 0;
      }
   else if ((2 * qx * qy + 2 * qz * qw) == -1.0)
      {
      heading = -2 * atan2(qx, qw);
      bank = 0.0;
      }
   else
      {
      heading = atan2(2 * qy * qw - 2 * qx * qz,
            1 - 2 * qy * qy - 2 * qz * qz);
      bank = atan2(2 * qx * qw - 2 * qy * qz,
            1 - 2 * qx * qx - 2 * qz * qz);
      }

   attitude = asin(2 * qx * qy + 2 * qz * qw);

   heading  = 180.0 * heading  / PI;
   attitude = 180.0 * attitude / PI;
   bank     = 180.0 * bank     / PI;

   eul(0, 0) = heading;
   eul(1, 0) = attitude;
   eul(2, 0) = bank;

   return (eul);

   }


void Rotation::FromMat9(const cv::Mat_<double>& mat9)
   {
   assert((mat9.rows == 3) && (mat9.cols == 3));

   quaternion[0] = sqrt(max(0.0, 1 + mat9(0, 0) + mat9(1, 1) + mat9(2, 2))) / 2.0;  // w
   quaternion[1] = sqrt(max(0.0, 1 + mat9(0, 0) - mat9(1, 1) - mat9(2, 2))) / 2.0;  // x
   quaternion[2] = sqrt(max(0.0, 1 - mat9(0, 0) + mat9(1, 1) - mat9(2, 2))) / 2.0;  // y
   quaternion[3] = sqrt(max(0.0, 1 - mat9(0, 0) - mat9(1, 1) + mat9(2, 2))) / 2.0;  // z

   quaternion[1] = quaternion[1] * Sign(mat9(2, 1) - mat9(1, 2)); // x
   quaternion[2] = quaternion[2] * Sign(mat9(0, 2) - mat9(2, 0)); // y
   quaternion[3] = quaternion[3] * Sign(mat9(1, 0) - mat9(0, 1)); // z

   Normalize();

   return;

   }


void Rotation::FromEul(const cv::Mat_<double>& eul)
   {
   double heading  = PI * eul(0, 0) / 180.0;
   double attitude = PI * eul(1, 0) / 180.0;
   double bank     = PI * eul(2, 0) / 180.0;

   double c1 = cos(heading / 2.0);
   double s1 = sin(heading / 2.0);
   double c2 = cos(attitude / 2.0);
   double s2 = sin(attitude / 2.0);
   double c3 = cos(bank / 2.0);
   double s3 = sin(bank / 2.0);
   double c1c2 = c1 * c2;
   double s1s2 = s1 * s2;

   quaternion[0] = c1c2 * c3  - s1s2 * s3;
   quaternion[1] = c1c2 * s3  + s1s2 * c3;
   quaternion[2] = s1 * c2 * c3 + c1 * s2 * s3;
   quaternion[3] = c1 * s2 * c3 - s1 * c2 * s3;

   Normalize();

   return;

   }

void Rotation::SetQuaternion(const cv::Mat_<double>& mat)
   {
   quaternion[0] = mat(0, 0);
   quaternion[1] = mat(1, 0);
   quaternion[2] = mat(2, 0);
   quaternion[3] = mat(3, 0);

   Normalize();

   return;

   }

void Rotation::SetQuaternion(const Rotation& rot)
   {
   quaternion = rot.quaternion;
   Normalize();

   return;

   }


void Rotation::SetEuler(const cv::Mat_<double>& mat)
   {
   FromEul(mat);

   return;

   }

void Rotation::SetMatrix(const cv::Mat& mat)
   {
   FromMat9(mat);

   return;

   }

Rotation& Rotation::operator += (const Rotation& r)
   {
   Rotation res;
   Mul(*this, r, res);
   *this = res;

   //x += v.x;
   //y += v.y;
   //z += v.z;

   return (*this);

   }

Rotation operator + (const Rotation& r1, const Rotation& r2)
   {
   Rotation ret = r1;
   ret += r2;

   return (ret);

   }

void Rotation::Output() const
   {
   cout << quaternion[0] << "," << quaternion[1] << "," << quaternion[2] <<
         "," << quaternion[3] << "|";

   return;

   }

} // namespace alvar
