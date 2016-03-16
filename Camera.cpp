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

#include "Camera.h"
//#include "FileFormatUtils.h"
#include <memory>

using namespace std;

namespace alvar {

using namespace std;

/*****************************************************************************
 ***  class ProjPoints
 ****************************************************************************/

void ProjPoints::Reset()
   {
   object_points.clear();
   image_points.clear();
   point_counts.clear();

   return;

   }

// TODO: Does it matter what the etalon_square_size is???
bool ProjPoints::AddPointsUsingChessboard(cv::Mat& image,
      double etalon_square_size, int etalon_rows, int etalon_columns,
      bool visualize)
   {
   if (image.cols == 0)
      return false;

   cv::Mat gray;
   std::vector<cv::Point2f> corners;

   if (image.channels() == 1)
      gray = image.clone();
   else
      cv::cvtColor(image, gray, CV_RGB2GRAY);

   width = image.cols;
   height = image.rows;

   int pattern_was_found = cv::findChessboardCorners(gray,
         cv::Size(etalon_rows, etalon_columns), corners,
         CV_CALIB_CB_ADAPTIVE_THRESH | CV_CALIB_CB_FAST_CHECK |
         CV_CALIB_CB_NORMALIZE_IMAGE);

   if (pattern_was_found != 0)
      {
      cv::cornerSubPix(gray, corners, cv::Size(5, 5), cv::Size(-1, -1),
            cv::TermCriteria(CV_TERMCRIT_ITER | CV_TERMCRIT_EPS, 10, 0.01f));
      for (size_t i = 0 ; i < corners.size() ; i++)
         {
         cv::Point3d po(etalon_square_size * (i % etalon_rows), etalon_square_size * (i / etalon_rows), 0.0);
         cv::Point2d pi(corners[i].x, corners[i].y);

         object_points.push_back(po);
         image_points.push_back(pi);
         }

      point_counts.push_back(corners.size());
      }

   if (visualize)
      {
      cv::drawChessboardCorners(image, cv::Size(etalon_rows, etalon_columns),
            corners, false /*pattern_was_found*/);
      }

   return (pattern_was_found != 0);

   }

bool ProjPoints::AddPointsUsingMarkers(vector<PointDouble>& marker_corners,
      vector<PointDouble>& marker_corners_img, cv::Mat& image)
   {
   width = image.cols;
   height = image.rows;
   if ((marker_corners.size() == marker_corners_img.size()) &&
         (marker_corners.size() == 4))
      {
      for (size_t p = 0 ; p < marker_corners.size() ; p++)
         {
         cv::Point3d po(marker_corners[p].x, marker_corners[p].y, 0.0);
         cv::Point2d pi(marker_corners_img[p]);

         object_points.push_back(po);
         image_points.push_back(pi);
         }

      point_counts.push_back(marker_corners.size());
      }

   return (true);

   }

/*****************************************************************************
 ***  class Camera
 ****************************************************************************/

Camera::Camera() : calib_K(3, 3), calib_D(4, 1)
   {
   calib_K = 0.0;
   calib_D = 0.0;

   calib_K(0, 0) = 550; // Just some focal length by default
   calib_K(1, 1) = 550; // Just some focal length by default
   calib_K(0, 2) = 320;
   calib_K(1, 2) = 240;
   calib_K(2, 2) = 1;

   calib_x_res = 640;
   calib_y_res = 480;
   x_res = 640;
   y_res = 480;

   return;

   }

void Camera::SetSimpleCalib(int _x_res, int _y_res, double f_fac)
   {
   calib_K = 0.0;
   calib_D = 0.0;

   calib_K(0, 0) = _x_res * f_fac; // Just some focal length by default
   calib_K(1, 1) = _x_res * f_fac; // Just some focal length by default
   calib_K(0, 2) = _x_res / 2;
   calib_K(1, 2) = _y_res / 2;
   calib_K(2, 2) = 1;

   calib_x_res = _x_res;
   calib_y_res = _y_res;

   return;

   }

bool Camera::LoadCalibOpenCV(const char* calibfile)
   {
   bool bRet = false;

   cv::FileStorage fs(calibfile, cv::FileStorage::READ);

   if (fs.isOpened())
      {
      cv::FileNode root_node = fs["calibration"];

		// K Intrinsic
      root_node["intrinsic_matrix"] >> calib_K;

		// D Distortion
      root_node["distortion"] >> calib_D;

		// Resolution
      root_node["width"] >> calib_x_res;
      root_node["height"] >> calib_y_res;

      bRet = true;
	}

   return (bRet);

   }

bool Camera::SetCalib(const char* calibfile, int _x_res, int _y_res, FILE_FORMAT format)
   {
   x_res = _x_res;
   y_res = _y_res;

   bool success = (calibfile != nullptr);
   if (success)
      {
      switch (format)
         {
         case FILE_FORMAT::FILE_FORMAT_OPENCV:
         case FILE_FORMAT::FILE_FORMAT_DEFAULT:
            success = LoadCalibOpenCV(calibfile);
            break;

         default:
            success = false;
            // TODO: throw exception?
            break;
         };

      if (success)
         {
         // Scale matrix in case of different resolution calibration.
         // The OpenCV documentation says:
         // - If an image from camera is up-sampled/down-sampled by some factor, all intrinsic camera parameters
         //   (fx, fy, cx and cy) should be scaled (multiplied/divided, respectively) by the same factor.
         // - The distortion coefficients remain the same regardless of the captured image resolution.
         if (calib_x_res != x_res)
            {
            double factor = static_cast<double>(x_res) / calib_x_res;
            calib_K(0, 0) *= factor;
            calib_K(0, 2) *= factor;
            }

         if (calib_y_res != y_res)
            {
            double factor = static_cast<double>(y_res) / calib_y_res;
            calib_K(1, 1) *= factor;
            calib_K(1, 2) *= factor;
            } // end if
         }
      } // end if

   return (success);

   }

bool Camera::SaveCalibOpenCV(const char* calibfile) const
   {
   cv::FileStorage fs(calibfile, cv::FileStorage::WRITE);

   bool ret = fs.isOpened();
   if (ret)
      {
      fs << "calibration" << "{"
         << "intrinsic_matrix" << calib_K
         << "distortion" << calib_D
         << "width" << calib_x_res
         << "height" << calib_y_res
         << "}";
      }

   return (ret);

   }

bool Camera::SaveCalib(const char* calibfile, FILE_FORMAT format) const
   {
   bool ret = false;

   if (calibfile != nullptr)
      {
      switch (format)
         {
         case FILE_FORMAT::FILE_FORMAT_OPENCV:
         case FILE_FORMAT::FILE_FORMAT_DEFAULT:
            ret = SaveCalibOpenCV(calibfile);
            break;

         default:
            break;
         }
      } // end if

   return (ret);

   }

// Needs to be brought up to date.
// I'm not using so don't feel like messing with it now!  TTC!!
void Camera::Calibrate(ProjPoints& pp)
   {
   CvMat *object_points = cvCreateMat(static_cast<int>(pp.object_points.size()), 1, CV_32FC3);
   CvMat *image_points = cvCreateMat(static_cast<int>(pp.image_points.size()), 1, CV_32FC2);
   const CvMat point_counts = cvMat(static_cast<int>(pp.point_counts.size()), 1, CV_32SC1, &pp.point_counts[0]);
   for (size_t i = 0 ; i < pp.object_points.size() ; i++)
      {
      object_points->data.fl[i*3+0] = static_cast<float>(pp.object_points[i].x);
      object_points->data.fl[i*3+1] = static_cast<float>(pp.object_points[i].y);
      object_points->data.fl[i*3+2] = static_cast<float>(pp.object_points[i].z);
      image_points->data.fl[i*2+0]  = static_cast<float>(pp.image_points[i].x);
      image_points->data.fl[i*2+1]  = static_cast<float>(pp.image_points[i].y);
      }

   CvMat CalibK = calib_K;
   CvMat CalibD = calib_D;
   cvCalibrateCamera2(object_points, image_points, &point_counts,
                      cvSize(pp.width, pp.height),
                      &CalibK, &CalibD, 0, 0, CV_CALIB_USE_INTRINSIC_GUESS);

   calib_x_res = pp.width;
   calib_y_res = pp.height;

   return;

   }

void Camera::SetRes(int _x_res, int _y_res)
   {
   x_res = _x_res;
   y_res = _y_res;

   // Scale matrix in case of different resolution calibration.
   // The OpenCV documentation says:
   // - If an image from camera is up-sampled/down-sampled by some factor, all intrinsic camera parameters
   //   (fx, fy, cx and cy) should be scaled (multiplied/divided, respectively) by the same factor.
   // - The distortion coefficients remain the same regardless of the captured image resolution.
   if (calib_x_res != x_res)
      {
      double scale = static_cast<double>(x_res) / calib_x_res;
      calib_K(0, 0) *= scale;
      calib_K(0, 2) *= scale;
      }

   if (calib_y_res != y_res)
      {
      double scale = static_cast<double>(y_res) / calib_y_res;
      calib_K(1, 1) *= scale;
      calib_K(1, 2) *= scale;
      } // end if

   return;

   }

// TODO: Better approach for this...
// Note, the proj_matrix[8] is now negated. This is due to the fact
// that with OpenCV and OpenGL projection matrices both y and z
// should be mirrored. All other components are 
void Camera::GetOpenglProjectionMatrix(double proj_matrix[16], const int width, const int height,
      const float far_clip /*= 1000.0f*/, const float near_clip /*= 0.1f*/) const
   {
   proj_matrix[0]	= 2.0 * calib_K(0, 0) / width;
   proj_matrix[1]	= 0.0;
   proj_matrix[2]	= 0.0;
   proj_matrix[3]	= 0.0;
   proj_matrix[4]  = 2.0 * calib_K(0, 1) / width; // skew
   proj_matrix[5]	= 2.0 * calib_K(1, 1) / height;
   proj_matrix[6]	= 0.0;
   proj_matrix[7]	= 0.0;
   //proj_matrix[8]	= (2.0f * calib_K_data[0][2] / float(width)) - 1.0f;
   proj_matrix[8]	= -(2.0 * calib_K(0, 2) / width) + 1.0;
   proj_matrix[9]	= (2.0 * calib_K(1, 2) / height) - 1.0;
   proj_matrix[10]	= -(far_clip + near_clip) / (far_clip - near_clip);
   proj_matrix[11]	= -1.0;
   proj_matrix[12]	= 0.0;
   proj_matrix[13]	= 0.0;
   proj_matrix[14]	= -2.0 * far_clip * near_clip / (far_clip - near_clip);
   proj_matrix[15]	= 0.0;

   return;

   }

void Camera::SetOpenglProjectionMatrix(double proj_matrix[16], const int width, const int height)
   {
   x_res = width;
   y_res = height;
   calib_x_res = width;
   calib_y_res = height;

   calib_K(0, 0) = proj_matrix[0] * width / 2.0;
   calib_K(0, 1) = proj_matrix[4] * width / 2.0;
   calib_K(1, 1) = proj_matrix[5] * height / 2.0;
   calib_K(0, 2) = (-proj_matrix[8] + 1.0) * width / 2.0; // Is this ok?
   calib_K(1, 2) = (proj_matrix[9] + 1.0) * height / 2.0;
   calib_K(2, 2) = 1.0;

   return;

   }

// A lot of duplicate code here  TTC!!
void Camera::Undistort(PointDouble& point) const
   {
   // focal length
   double ifx = 1./ calib_K(0, 0);
   double ify = 1./ calib_K(1, 1);

   // principal point
   double cx = calib_K(0, 2);
   double cy = calib_K(1, 2);

   // distortion coeffs
   const double* k = calib_D.ptr<double>();

   // compensate distortion iteratively
   double x = (point.x - cx) * ifx;
   double y = (point.y - cy) * ify;
   double x0 = x;
   double y0 = y;
   for (int j = 0 ; j < 5 ; j++)
      {
      double r2 = x * x + y * y;
      double icdist = 1.0 / (1.0  + k[0] * r2 + k[1] * r2 * r2);
      double delta_x = 2.0 * k[2] * x * y + k[3] * (r2 + 2 * x * x);
      double delta_y = k[2] * (r2 + 2 * y * y) + 2 * k[3] * x * y;
      x = (x0 - delta_x) * icdist;
      y = (y0 - delta_y) * icdist;
      }

   // apply compensation
   point.x = x / ifx + cx;
   point.y = y / ify + cy;

   return;

   }

void Camera::Undistort(vector<PointDouble>& points) const
   {
   // focal length
   double ifx = 1.0 / calib_K(0, 0);
   double ify = 1.0 / calib_K(1, 1);

   // principal point
   double cx = calib_K(0, 2);
   double cy = calib_K(1, 2);

   // distortion coeffs
   const double* k = calib_D.ptr<double>();

   for (size_t i = 0 ; i < points.size() ; i++)
      {
      // compensate distortion iteratively
      double x = (points[i].x - cx) * ifx;
      double y = (points[i].y - cy) * ify;
      double x0 = x;
      double y0 = y;
      for (int j = 0 ; j < 5 ; j++)
         {
         double r2 = x * x + y * y;
         double icdist = 1.0 / (1.0  + k[0] * r2 + k[1] * r2 * r2);
         double delta_x = 2.0 * k[2] * x * y + k[3] * (r2 + 2 * x * x);
         double delta_y = k[2] * (r2 + 2 * y * y) + 2 * k[3] * x * y;
         x = (x0 - delta_x) * icdist;
         y = (y0 - delta_y) * icdist;
         }

      // apply compensation
      points[i].x = x / ifx + cx;
      points[i].y = y / ify + cy;
      }

   return;

   }

void Camera::Undistort(cv::Point2f& point) const
   {
   // focal length
   double ifx = 1./ calib_K(0, 0);
   double ify = 1./ calib_K(1, 1);

   // principal point
   double cx = calib_K(0, 2);
   double cy = calib_K(1, 2);

   // distortion coeffs
   const double* k = calib_D.ptr<double>();

   // compensate distortion iteratively
   double x = (point.x - cx) * ifx;
   double y = (point.y - cy) * ify;
   double x0 = x;
   double y0 = y;
   for (int j = 0 ; j < 5 ; j++)
      {
      double r2 = x * x + y * y;
      double icdist = 1.0 / (1.0  + k[0] * r2 + k[1] * r2 * r2);
      double delta_x = 2.0 * k[2] * x * y + k[3] * (r2 + 2 * x * x);
      double delta_y = k[2] * (r2 + 2 * y * y) + 2 * k[3] * x * y;
      x = (x0 - delta_x) * icdist;
      y = (y0 - delta_y) * icdist;
      }

   // apply compensation
   point.x = static_cast<float>(x / ifx + cx);
   point.y = static_cast<float>(y / ify + cy);

   return;

   }

void Camera::Distort(vector<PointDouble>& points) const
   {
   // cx, cy
   double u0 = calib_K(0, 2);
   double v0 = calib_K(1, 2);
   double fx = calib_K(0, 0);
   double fy = calib_K(1, 1);
   double _fx = 1.0 / fx;
   double _fy = 1.0 / fy;
   const double* k = calib_D.ptr<double>();

   double k1 = k[0];
   double k2 = k[1];
   double p1 = k[2];
   double p2 = k[3];

   for (auto& pt : points)
      {
      // Distort
      double y = (pt.y - v0) * _fy;
      double y2 = y * y;
      double _2p1y = 2 * p1 * y;
      double _3p1y2 = 3 * p1 * y2;
      double p2y2 = p2 * y2;

      double x = (pt.x - u0) * _fx;
      double x2 = x * x;
      double r2 = x2 + y2;
      double d = 1.0 + (k1 + k2 * r2) * r2;

      pt.x = fx * (x * (d + _2p1y) + p2y2 + (3 * p2) * x2) + u0;
      pt.y = fy * (y * (d + (2 * p2) * x) + _3p1y2 + p1 * x2) + v0;
      }

   return;

   }

void Camera::Distort(PointDouble& point) const
   {
   // cx, cy
   double u0 = calib_K(0, 2);
   double v0 = calib_K(1, 2);
   double fx = calib_K(0, 0);
   double fy = calib_K(1, 1);
   double _fx = 1.0 / fx;
   double _fy = 1.0 / fy;
   const double* k = calib_D.ptr<double>();

   double k1 = k[0];
   double k2 = k[1];
   double p1 = k[2];
   double p2 = k[3];

   // Distort
   double y = (point.y - v0) * _fy;
   double y2 = y * y;
   double _2p1y = 2 * p1 * y;
   double _3p1y2 = 3 * p1 * y2;
   double p2y2 = p2 * y2;

   double x = (point.x - u0) * _fx;
   double x2 = x * x;
   double r2 = x2 + y2;
   double d = 1.0 + (k1 + k2 * r2) * r2;

   point.x = fx * (x * (d + _2p1y) + p2y2 + (3 * p2) * x2) + u0;
   point.y = fy * (y * (d + (2 * p2) * x) + _3p1y2 + p1 * x2) + v0;

   return;

   }

void Camera::Distort(cv::Point2f& point) const
   {
   // cx, cy
   double u0 = calib_K(0, 2);
   double v0 = calib_K(1, 2);
   double fx = calib_K(0, 0);
   double fy = calib_K(1, 1);
   double _fx = 1.0 / fx;
   double _fy = 1.0 / fy;
   const double* k = calib_D.ptr<double>();

   double k1 = k[0];
   double k2 = k[1];
   double p1 = k[2];
   double p2 = k[3];

   // Distort
   double y = (point.y - v0) * _fy;
   double y2 = y * y;
   double _2p1y = 2 * p1 * y;
   double _3p1y2 = 3 * p1 * y2;
   double p2y2 = p2 * y2;

   double x = (point.x - u0) * _fx;
   double x2 = x * x;
   double r2 = x2 + y2;
   double d = 1.0 + (k1 + k2 * r2) * r2;

   point.x = static_cast<float>(fx * (x * (d + _2p1y) + p2y2 + (3 * p2) * x2) + u0);
   point.y = static_cast<float>(fy * (y * (d + (2 * p2) * x) + _3p1y2 + p1 * x2) + v0);

   return;

   }

bool Camera::CalcExteriorOrientation(vector<cv::Point3d>& pw, vector<cv::Point2d>& pi, Pose& pose) const
   {
   cv::Mat_<double> ext_rodriques_mat(3, 1);
   cv::Mat ext_translate_mat(3, 1, CV_64F);

   bool ret = cv::solvePnP(pw, pi, calib_K, calib_D, ext_rodriques_mat, ext_translate_mat);

   if (ret)
      {
      pose.SetRodrigues(ext_rodriques_mat);
      pose.SetTranslation(ext_translate_mat);
      } // end if

   return (ret);

   }

bool Camera::CalcExteriorOrientation(vector<PointDouble>& pw,
      vector<PointDouble>& pi, cv::Mat& rodriques, cv::Mat& tra) const
   {
   //assert(pw.size() == pi.size());
   size_t size = pi.size();

   vector<cv::Point3d> pw3(size);
   vector<cv::Point2d> pi2(size);

   for (size_t i = 0 ; i < size ; i++)
      {
      pw3[i].x = pw[i].x;
      pw3[i].y = pw[i].y;
      pw3[i].z = 0;

      pi2[i] = pi[i];
      }

   return (CalcExteriorOrientation(pw3, pi2, rodriques, tra));

   }

// Currently the only version called.  TTC!!!
bool Camera::CalcExteriorOrientation(vector<PointDouble>& pw, vector<PointDouble>& pi, Pose& pose) const
   {
   cv::Mat ext_translate_mat(3, 1, CV_64F);
   cv::Mat_<double> ext_rodriques_mat(3, 1);

   bool ret = CalcExteriorOrientation(pw, pi, ext_rodriques_mat, ext_translate_mat);

   if (ret)
      {
      pose.SetRodrigues(ext_rodriques_mat);
      pose.SetTranslation(ext_translate_mat);
      } // end if

   return (ret);

   }

bool Camera::CalcExteriorOrientation(const cv::Mat& object_points, cv::Mat& image_points, Pose& pose) const
   {
   cv::Mat_<double> ext_rodriques_mat(3, 1);
   cv::Mat ext_translate_mat(3, 1, CV_64F);

   bool ret = CalcExteriorOrientation(object_points, image_points, ext_rodriques_mat, ext_translate_mat);

   if (ret)
      {
      pose.SetRodrigues(ext_rodriques_mat);
      pose.SetTranslation(ext_translate_mat);
      } // end if

   return (ret);

   }

void Camera::ProjectPoints(const std::vector<cv::Point3d>& pw, const Pose& pose, std::vector<cv::Point2d>& pi) const
   {
   cv::Mat ext_translate_mat(3, 1, CV_64F);
   pose.GetTranslation(ext_translate_mat);

   cv::projectPoints(pw, pose.GetRodrigues(), ext_translate_mat, calib_K, calib_D, pi);

   return;

   }

void Camera::ProjectPoints(const cv::Mat& object_points, const Pose& pose, cv::Mat& image_points) const
   {
   cv::Mat ext_translate_mat(3, 1, CV_64F);
   pose.GetTranslation(ext_translate_mat);
   cv::Mat_<double> ext_rodriques_mat = pose.GetRodrigues();

   cv::projectPoints(object_points, ext_rodriques_mat, ext_translate_mat, calib_K, calib_D, image_points);

   return;

   }

void Camera::ProjectPoints(const cv::Mat& object_points, double gl[16], cv::Mat& image_points) const
   {
   double glm[4][4] =
   {
      gl[0], gl[4], gl[8],  gl[12],
      gl[1], gl[5], gl[9],  gl[13],
      gl[2], gl[6], gl[10], gl[14],
      gl[3], gl[7], gl[11], gl[15],
   };

   cv::Mat glm_mat(4, 4, CV_64F, glm);

   // For some reason we need to mirror both y and z ???
   cv::Mat_<double> cv_mul(4, 4);
   cv_mul = cv::Mat::eye(4, 4, CV_64F);
   cv_mul(1, 1) = -1.0;
   cv_mul(2, 2) = -1.0;
   glm_mat = cv_mul * glm_mat;

   // Rotation
   Rotation r;
   r.SetMatrix(glm_mat);
   cv::Mat_<double> rod_mat = r.GetRodrigues();
   // Translation
   double tra[3] = { glm[0][3], glm[1][3], glm[2][3] };
   cv::Mat tra_mat(3, 1, CV_64F, tra);
   // Project points
   ProjectPoints(object_points, rod_mat, tra_mat, image_points);

   return;

   }

void Camera::ProjectPoint(const cv::Point3d& pw, const Pose& pose, cv::Point2d& pi) const
   {
   std::vector<cv::Point3d> object_points(1, pw);
   std::vector<cv::Point2d> image_points(1, pi);

   ProjectPoints(object_points, pose, image_points);

   pi = image_points[0];

   return;

   }
#if 0
void Camera::ProjectPoint(const cv::Point3f& pw, const Pose& pose, cv::Point2f& pi) const
   {
   cv::Mat_<cv::Point3f> object_points(1, 1);
   object_points(0, 0) = pw;
   cv::Mat_<cv::Point2f> image_points(1, 1);
   image_points(0, 0) = cv::Point2f(0.0f, 0.0f);
   ProjectPoints(object_points, pose, image_points);
   pi = image_points(0, 0);

   return;

   }
#endif

/*****************************************************************************
 ***  class Homography
 ****************************************************************************/

void Homography::Find(const vector<PointDouble>& pw, const vector<PointDouble>& pi)
   {
   assert(pw.size() == pi.size());
   size_t size = pi.size();

   std::vector<cv::Point2f> srcpts(size);
   std::vector<cv::Point2f> dstpts(size);

   for (size_t i = 0 ; i < size ; ++i)
      {
      srcpts[i] = pw[i];
      dstpts[i] = pi[i];
      }

   H = cv::findHomography(srcpts, dstpts);

   return;

   }

void Homography::ProjectPoints(const vector<PointDouble>& from, vector<PointDouble>& to) const
   {
   std::vector<cv::Point2d> srcpts;
   srcpts.reserve(from.size());
   std::vector<cv::Point2d> dstpts;

   for (auto& fpt : from)
      {
      srcpts.push_back(fpt);
      } // end for

   cv::perspectiveTransform(srcpts, dstpts, H);

   to.clear();
   for (auto& dstpt : dstpts)
      {
      PointDouble dpt(dstpt.x, dstpt.y);
      to.push_back(dpt);
      } // end for

   return;

   }

} // namespace alvar
