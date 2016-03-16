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

/*****************************************************************************
 ******************************  I N C L U D E  ******************************
 ****************************************************************************/

#include "Alvar.h"
#include "Marker.h"
#include <opencv2/highgui/highgui.hpp>

namespace alvar {

using namespace std;

/*****************************************************************************
 ***  class Marker
 ****************************************************************************/

Marker::Marker(double _edge_length, int _res, double _margin)
   {
   margin_error = 0;
   decode_error = 0;
   track_error = 0;
   SetMarkerSize(_edge_length, _res, _margin);

   return;

   }

Marker::Marker(const Marker& m) : marker_content(m.marker_content)
   {
   SetMarkerSize(m.edge_length, m.res, m.margin);

   pose = m.pose;
   margin_error = m.margin_error;
   decode_error = m.decode_error;
   track_error = m.track_error;
   //   cvCopy(m.marker_content, marker_content);

   marker_corners.resize(m.marker_corners.size());
   copy(m.marker_corners.begin(), m.marker_corners.end(), marker_corners.begin());

   marker_points.resize(m.marker_points.size());
   copy(m.marker_points.begin(), m.marker_points.end(), marker_points.begin());

   marker_corners_img.resize(m.marker_corners_img.size());
   copy(m.marker_corners_img.begin(), m.marker_corners_img.end(), marker_corners_img.begin());

#ifdef VISUALIZE_MARKER_POINTS
   marker_allpoints_img.resize(m.marker_allpoints_img.size());
   copy(m.marker_allpoints_img.begin(), m.marker_allpoints_img.end(), marker_allpoints_img.begin());
#endif

   return;

   }

void Marker::VisualizeMarkerPose(cv::Mat& image, Camera& /* cam */,
      std::vector<cv::Point2d>& visualize2d_points, cv::Scalar color) const
   {
	// Cube
   for (size_t i = 0 ; i < 4 ; i++)
      {
      cv::line(image, visualize2d_points[i],
            visualize2d_points[(i + 1) % 4], color);
      cv::line(image, visualize2d_points[i],
            visualize2d_points[4 + i], color);
      cv::line(image, visualize2d_points[4 + i],
            visualize2d_points[4 + (( i+ 1) % 4)], color);
      }

	// Coordinates
   cv::line(image, visualize2d_points[8], visualize2d_points[9], cv::Scalar(0, 0, 255));
   cv::line(image, visualize2d_points[8], visualize2d_points[10], cv::Scalar(0, 255, 0));
   cv::line(image, visualize2d_points[8], visualize2d_points[11], cv::Scalar(255, 0, 0));

   return;

   }

void Marker::VisualizeMarkerContent(cv::Mat& image, Camera& /* cam */, cv::Point2d datatext_point, cv::Point2d content_point) const
{
#ifdef VISUALIZE_MARKER_POINTS
   for (const auto& pt : marker_allpoints_img)
      {
      if (pt.val == 0)
         cv::circle(image, pt, 1, cv::Scalar(0, 255,0));
      else if (pt.val == 255)
         cv::circle(image, pt, 1, cv::Scalar(0, 0, 255));
      else
         cv::circle(image, pt, 2, cv::Scalar(0, 255, 255));
      }
#endif

   // Marker data
   std::stringstream val;
   val << int(GetId());
   cv::putText(image, val.str(), datatext_point, cv::FONT_HERSHEY_PLAIN, 0.5,
         cv::Scalar(0, 255, 255));

   // MarkerContent
   int xc = int(content_point.x);
   int yc = int(content_point.y);
   for (int j = 0 ; j < res * 3 ; j++)
      {
      for (int i = 0 ; i < res * 3 ; i++)
         {
         int x = xc + i;
         int y = yc + j;
         if ((x >= 0) && (x < image.cols) &&
             (y >= 0) && (y < image.rows))
            {
            if (marker_content.at<uchar>(j / 3, i / 3))
               {
               image.at<cv::Vec3b>(y, x) = cv::Vec3b(255, 255, 255);
               }
            else
               {
               image.at<cv::Vec3b>(y, x) = cv::Vec3b(0, 0, 0);
               }
            }
         }
      }

   return;

   }

void Marker::VisualizeMarkerError(cv::Mat& image, Camera& /* cam */, cv::Point2d errortext_point) const
   {
	std::stringstream val;
   if (GetError(MARGIN_ERROR | DECODE_ERROR) > 0)
      {
		val.str("");
      val << int(GetError(MARGIN_ERROR) * 100) << "% ";
      val << int(GetError(DECODE_ERROR) * 100) << "% ";
      cv::putText(image, val.str(), errortext_point, cv::FONT_HERSHEY_PLAIN, 0.5, cv::Scalar(0, 0, 255));
      }
   else if (GetError(TRACK_ERROR) > 0.01)
      {
		val.str("");
      val << int(GetError(TRACK_ERROR) * 100) << "%";
      cv::putText(image, val.str(), errortext_point, cv::FONT_HERSHEY_PLAIN, 0.5, cv::Scalar(0, 0, 128));
      }

   return;

   }

void MarkerData::VisualizeMarkerContent(cv::Mat& image, Camera& /* cam */,
      cv::Point2d datatext_point, cv::Point2d /* content_point */) const
   {
#ifdef VISUALIZE_MARKER_POINTS
   for (const auto& pt : marker_allpoints_img)
      {
      if (pt.val == 0)
         cv::circle(image, pt, 1, cv::Scalar(0, 255, 0));
      else if (pt.val == 255)
         cv::circle(image, pt, 1, cv::Scalar(0, 0, 255));
      else
         cv::circle(image, pt, 2, cv::Scalar(0, 255, 255));
      }
#endif

	// Marker data
   std::stringstream val;
   cv::Scalar rgb(0, 255, 255);
   if (content_type == MarkerContentType::MARKER_CONTENT_TYPE_NUMBER)
      {
      val << int(GetId());
      }
   else
      {
      if (content_type == MarkerContentType::MARKER_CONTENT_TYPE_FILE)
         rgb = cv::Scalar(255, 0, 0);
      if (content_type == MarkerContentType::MARKER_CONTENT_TYPE_HTTP)
         rgb = cv::Scalar(255, 0, 255);

      val << data.str;
      }

   cv::putText(image, val.str(), datatext_point, cv::FONT_HERSHEY_PLAIN, 2.0, rgb);

   return;

   }

void Marker::Visualize(cv::Mat& image, Camera& cam, cv::Scalar color) const
   {
   std::vector<cv::Point3d> visualize3d_points =
      {
      // Cube vertices
      {-(edge_length / 2), -(edge_length / 2), 0},
      {-(edge_length / 2),  (edge_length / 2), 0},
      { (edge_length / 2),  (edge_length / 2), 0},
      { (edge_length / 2), -(edge_length / 2), 0},
      {-(edge_length / 2), -(edge_length / 2), edge_length},
      {-(edge_length / 2),  (edge_length / 2), edge_length},
      { (edge_length / 2),  (edge_length / 2), edge_length},
      { (edge_length / 2), -(edge_length / 2), edge_length},
      // Coodinate axes
      {0, 0, 0},
      {edge_length, 0, 0},
      {0, edge_length, 0},
      {0, 0, edge_length}
      };

   std::vector<cv::Point2d> visualize2d_points(12);

   cam.ProjectPoints(visualize3d_points, pose, visualize2d_points);

	VisualizeMarkerPose(image, cam, visualize2d_points, color);
   VisualizeMarkerContent(image, cam, visualize2d_points[0], visualize2d_points[8]);
   VisualizeMarkerError(image, cam, visualize2d_points[2]);

   return;

   }

void Marker::CompareCorners(vector<PointDouble>& _marker_corners_img, int* orientation, double* error)
   {
   vector<PointDouble >::iterator corners_new = _marker_corners_img.begin();
   vector<PointDouble >::const_iterator corners_old = marker_corners_img.begin();
   vector<double> errors(4);
   for (int i = 0 ; i < 4 ; i++)
      {
      errors[0] += PointSquaredDistance(marker_corners_img[i], _marker_corners_img[i]);
      errors[1] += PointSquaredDistance(marker_corners_img[i], _marker_corners_img[(i + 1) % 4]);
      errors[2] += PointSquaredDistance(marker_corners_img[i], _marker_corners_img[(i + 2) % 4]);
      errors[3] += PointSquaredDistance(marker_corners_img[i], _marker_corners_img[(i + 3) % 4]);
      }

   *orientation = min_element(errors.begin(), errors.end()) - errors.begin();
   *error = sqrt(errors[*orientation] / 4);
   *error /= sqrt(max(PointSquaredDistance(marker_corners_img[0], marker_corners_img[2]),
         PointSquaredDistance(marker_corners_img[1], marker_corners_img[3])));

   return;

   }

#if 0  // TTC!!
void Marker::CompareContent(vector<PointDouble>& _marker_corners_img, cv::Mat& gray, Camera *cam, int *orientation) const {
	// TODO: Note, that to use this method you need to straighten the content
	// TODO: This method can be used with image based trackingt

}
#endif

bool Marker::UpdateContentBasic(vector<PointDouble>& _marker_corners_img, cv::Mat& gray, Camera& cam,
      int /* frame_no */ /*= 0*/)
   {
   vector<PointDouble> marker_corners_img_undist;
   marker_corners_img_undist.resize(_marker_corners_img.size());
   copy(_marker_corners_img.begin(), _marker_corners_img.end(), marker_corners_img_undist.begin());

   // Figure out the marker point position in the image
   Homography H;
   vector<PointDouble> marker_points_img(marker_points.size());
   marker_points_img.resize(marker_points.size());
   cam.Undistort(marker_corners_img_undist);
   H.Find(marker_corners, marker_corners_img_undist);
   H.ProjectPoints(marker_points, marker_points_img);
   cam.Distort(marker_points_img);

   // Read the content
   int x, y;
   double min = 255.0, max = 0.0;
   for (int j = 0 ; j < marker_content.rows ; j++)
      {
      for (int i = 0 ; i < marker_content.cols ; i++)
         {
         x = static_cast<int>(0.5 + Limit(marker_points_img[(j * marker_content.cols) + i].x, 1.0,
               static_cast<double>(gray.cols - 2)));
         y = static_cast<int>(0.5 + Limit(marker_points_img[(j * marker_content.cols) + i].y, 1.0,
             static_cast<double>(gray.rows - 2)));

         marker_points_img[(j * marker_content.cols) + i].val = gray.at<uchar>(y, x);
         /*
         // Use median of 5 neighbor pixels
         vector<int> vals;
         vals.clear();
         vals.push_back();
         vals.push_back((int)cvGetReal2D(gray, y-1, x));
         vals.push_back((int)cvGetReal2D(gray, y, x-1));
         vals.push_back((int)cvGetReal2D(gray, y+1, x));
         vals.push_back((int)cvGetReal2D(gray, y, x+1));
         nth_element(vals.begin(), vals.begin()+2, vals.end());
         tmp = vals[2];
         */

         marker_content.at<uchar>(j, i) = marker_points_img[(j * marker_content.cols) + i].val;
         if (marker_points_img[(j * marker_content.cols) + i].val > max)
            max = marker_points_img[(j * marker_content.cols) + i].val;
         if (marker_points_img[(j * marker_content.cols) + i].val < min)
            min = marker_points_img[(j * marker_content.cols) + i].val;
         }
      }

   // Take few additional points from border and just
   // outside the border to make the right thresholding
   vector<PointDouble> marker_margin_w_img(marker_margin_w.size());
   vector<PointDouble> marker_margin_b_img(marker_margin_b.size());
   H.ProjectPoints(marker_margin_w, marker_margin_w_img);
   H.ProjectPoints(marker_margin_b, marker_margin_b_img);
   cam.Distort(marker_margin_w_img);
   cam.Distort(marker_margin_b_img);

   min = max = 0; // Now min and max values are averages over black and white border pixels.
   for (auto& marker_margin : marker_margin_w_img)
      {
      x = static_cast<int>(0.5 + Limit(marker_margin.x, 0.0, static_cast<double>(gray.cols - 1)));
      y = static_cast<int>(0.5 + Limit(marker_margin.y, 0.0, static_cast<double>(gray.rows - 1)));
      marker_margin.val = gray.at<uchar>(y, x);
      max += marker_margin.val;
      //if(marker_margin_w_img[i].val > max) max = marker_margin_w_img[i].val;
      //if(marker_margin_w_img[i].val < min) min = marker_margin_w_img[i].val;
      }

   for (auto& marker_margin : marker_margin_b_img)
      {
      x = static_cast<int>(0.5 + Limit(marker_margin.x, 0.0, static_cast<double>(gray.cols - 1)));
      y = static_cast<int>(0.5 + Limit(marker_margin.y, 0.0, static_cast<double>(gray.rows - 1)));
      marker_margin.val = gray.at<uchar>(y, x);
      min += marker_margin.val;
      //if(marker_margin_b_img[i].val > max) max = marker_margin_b_img[i].val;
      //if(marker_margin_b_img[i].val < min) min = marker_margin_b_img[i].val;
      }

   max /= marker_margin_w_img.size();
   min /= marker_margin_b_img.size();

   // Threshold the marker content
   cv::threshold(marker_content, marker_content, (max + min) / 2.0, 255, cv::THRESH_BINARY);

   // Count erroneous margin nodes
   int erroneous = 0;
   int total = 0;
   for (auto& marker_margin : marker_margin_w_img)
      {
      if (marker_margin.val < (max + min) / 2.0)
         erroneous++;
      total++;
      }

   for (auto& marker_margin : marker_margin_b_img)
      {
      if (marker_margin.val > (max + min) / 2.0)
         erroneous++;
      total++;
      }

   margin_error = static_cast<double>(erroneous) / total;
   //	track_error;  // TTC!!

#ifdef VISUALIZE_MARKER_POINTS
   // Now we fill also this temporary debug table for visualizing marker code reading
   // TODO: this whole vector is only for debug purposes
   marker_allpoints_img.clear();
   for (auto& marker_margin : marker_margin_w_img)
      {
      PointDouble p = marker_margin;
      if (p.val < (max + min) / 2.0)
         p.val = 255; // error
      else
         p.val = 0; // ok

      marker_allpoints_img.push_back(p);
      }

   for (auto& marker_margin : marker_margin_b_img)
      {
      PointDouble p = marker_margin;
      if (p.val > (max + min) / 2.0)
         p.val = 255; // error
      else
         p.val = 0; // ok

      marker_allpoints_img.push_back(p);
      }

   for (auto& marker_point : marker_points_img)
      {
      PointDouble p = marker_point;
      p.val = 128; // Unknown?
      marker_allpoints_img.push_back(p);
      }
#endif

   return true;
   }

void Marker::UpdatePose(vector<PointDouble>& _marker_corners_img, Camera& cam, int orientation,
      int /* frame_no */ /* =0 */, bool update_pose /* =true */)
   {
   marker_corners_img.resize(_marker_corners_img.size());
   copy(_marker_corners_img.begin(), _marker_corners_img.end(), marker_corners_img.begin());

   // Calculate exterior orientation
   if (orientation > 0)
      std::rotate(marker_corners_img.begin(), marker_corners_img.begin() + orientation, marker_corners_img.end());

   if (update_pose)
      cam.CalcExteriorOrientation(marker_corners, marker_corners_img, pose);

   return;

   }

bool Marker::DecodeContent(int& orientation)
   {
   orientation = 0;
   decode_error = 0;

   return (true);

   }

void Marker::SaveMarkerImage(const char* filename, int save_res) const
   {
   if (save_res == 0)
      {
      // TODO: More intelligent deduction of a minimum save_res
      save_res = static_cast<int>((res + margin + margin) * 12);
      }

   double scale = save_res / (res + margin + margin);

   // Create a black image
   cv::Mat img(cv::Size(save_res, save_res), CV_8U, cv::Scalar(0));

   // Resize anc copy the marker to the black image at the correct position
   int submat_pos = static_cast<int>(margin * scale);
   int submat_size = static_cast<int>(res * scale);
   cv::Mat submat(img, cv::Rect(submat_pos, submat_pos, submat_size, submat_size));

   int resize = static_cast<int>(res * scale + 0.5);
   cv::Mat img_content;
   cv::resize(marker_content, img_content, cv::Size(resize, resize), cv::INTER_NEAREST);
   img_content.copyTo(submat);

   cv::imwrite(filename, img);

   return;

   }

void Marker::ScaleMarkerToImage(cv::Mat& image) const
   {
   const int multiplier = 96;

   // Create a black image
   int image_size = static_cast<int>(multiplier * (res + margin + margin) + 0.5);
   cv::Mat img(cv::Size(image_size, image_size), CV_8U, cv::Scalar(0));

   // Resize anc copy the marker to the black image at the correct position
   cv::Mat img_content;
   int submat_pos = static_cast<int>(multiplier * margin + 0.5);
   int submat_size = static_cast<int>(multiplier * res + 0.5);
   cv::Mat submat(img, cv::Rect(submat_pos, submat_pos, submat_size, submat_size));
   cv::resize(marker_content, img_content, cv::Size(submat_size, submat_size), cv::INTER_NEAREST);
   img_content.copyTo(submat);

   image = img;

   return;

   }

void Marker::SetMarkerSize(double _edge_length, int _res, double _margin)
   {
   // TODO: Is this right place for default marker size?
   edge_length = (_edge_length ? _edge_length : 1);
   res = _res; //(_res?_res:10);
   margin = ((_margin != 0.0) ? _margin : 1);

   double x_min = -0.5 * edge_length;
   double y_min = -0.5 * edge_length;
   double x_max = 0.5 * edge_length;
   double y_max = 0.5 * edge_length;

   double cx_min = (x_min * res) / (res + margin + margin);
   double cy_min = (y_min * res) / (res + margin + margin);
   // 	double cx_max = (x_max * res)/(res + margin + margin);  // TTC!!
   double cy_max = (y_max * res) / (res + margin + margin);
   double step = edge_length / (res + margin + margin);

   // marker_corners
   marker_corners_img.resize(4);

   // Same order as the detected corners
   marker_corners.clear();
   marker_corners.push_back(PointDouble(x_min, y_min));
   marker_corners.push_back(PointDouble(x_max, y_min));
   marker_corners.push_back(PointDouble(x_max, y_max));
   marker_corners.push_back(PointDouble(x_min, y_max));

   // Rest can be done only if we have existing resolution
   if (res > 0)
      {
      // marker_points
      marker_points.clear();
      for (int j = 0 ; j < res ; ++j)
         {
         for (int i = 0 ; i < res ; ++i)
            {
            PointDouble  pt;
            pt.y = cy_max - (step * j) - (step / 2);
            pt.x = cx_min + (step * i) + (step / 2);
            marker_points.push_back(pt);
            }
         }

      // Samples to be used in margins
      // TODO: Now this works only if the "margin" is without decimals
      // TODO: This should be made a lot cleaner
      marker_margin_w.clear();
      marker_margin_b.clear();
      for (int j = -1 ; j <= margin - 1 ; j++)
         {
         PointDouble  pt;
         // Sides
         for (int i = 0 ; i < res ; i++)
            {
            pt.x = cx_min + step * i + step / 2;
            pt.y = y_min + step * j + step / 2;
            if (j < 0)
               marker_margin_w.push_back(pt);
            else
               marker_margin_b.push_back(pt);

            pt.y = y_max - step * j - step / 2;

            if (j < 0)
               marker_margin_w.push_back(pt);
            else
               marker_margin_b.push_back(pt);

            pt.y = cy_min + step * i + step / 2;
            pt.x = x_min + step * j + step / 2;

            if (j < 0)
               marker_margin_w.push_back(pt);
            else
               marker_margin_b.push_back(pt);

            pt.x = x_max - step * j - step / 2;

            if (j < 0)
               marker_margin_w.push_back(pt);
            else
               marker_margin_b.push_back(pt);
            }

         // Corners
         for (int i = -1 ; i <= margin - 1 ; i++)
            {
            pt.x = x_min + step * i + step / 2;
            pt.y = y_min + step * j + step / 2;
            if ((j < 0) || (i < 0))
               marker_margin_w.push_back(pt);
            else
               marker_margin_b.push_back(pt);

            pt.x = x_min + step * i + step / 2;
            pt.y = y_max - step * j - step / 2;
            if ((j < 0) || (i < 0))
               marker_margin_w.push_back(pt);
            else
               marker_margin_b.push_back(pt);

            pt.x = x_max - step * i - step / 2;
            pt.y = y_max - step * j - step / 2;
            if ((j < 0) || (i < 0))
               marker_margin_w.push_back(pt);
            else
               marker_margin_b.push_back(pt);

            pt.x = x_max - step * i - step / 2;
            pt.y = y_min + step * j + step / 2;
            if ((j < 0) || (i < 0))
               marker_margin_w.push_back(pt);
            else
               marker_margin_b.push_back(pt);
            }
         }

      /*
   for(int j = -margin-1; j < res+margin+margin+2; ++j) {
      for(int i = 0; i < res+margin+margin+2; ++i) {
         PointDouble  pt;
         pt.y = y_min - step/2 + step*j;
         pt.x = x_min - step/2 + step*i;
         if ((pt.x < x_min) || (pt.y < y_min) ||
            (pt.x > x_max) || (pt.y > y_max))
         {
            marker_margin_w.push_back(pt);
         }
         else
         if ((pt.x < cx_min) || (pt.y < cy_min) ||
            (pt.x > cx_max) || (pt.y > cy_max))
         {
            marker_margin_b.push_back(pt);
         }
      }
   }
   */
      /*
   //double step = edge_length / (res + margin + margin);
   for(int j = 0; j < res+margin+margin+2; ++j) {
      for(int i = 0; i < res+margin+margin+2; ++i) {
         PointDouble  pt;
         pt.y = y_min - step/2 + step*j;
         pt.x = x_min - step/2 + step*i;
         if ((pt.x < x_min) || (pt.y < y_min) ||
            (pt.x > x_max) || (pt.y > y_max))
         {
            marker_margin_w.push_back(pt);
         }
         else
         if ((pt.x < cx_min) || (pt.y < cy_min) ||
            (pt.x > cx_max) || (pt.y > cy_max))
         {
            marker_margin_b.push_back(pt);
         }
      }
   }
   */
      /*
   marker_margin_w.clear();
   marker_margin_b.clear();
   for (double y=y_min-(step/2); y<y_max+(step/2); y+=step) {
      for (double x=x_min-(step/2); x<x_max+(step/2); x+=step) {
         PointDouble pt(x, y);
         if ((x < x_min) || (y < y_min) ||
            (x > x_max) || (y > y_max))
         {
            marker_margin_w.push_back(pt);
         }
         else
         if ((x < cx_min) || (y < cy_min) ||
            (x > cx_max) || (y > cy_max))
         {
            marker_margin_b.push_back(pt);
         }
      }
   }
   */
      /*
   marker_points.clear();
   marker_margin_w.clear();
   marker_margin_b.clear();
   for(int j = 0; j < res+margin+margin+2; ++j) {
      for(int i = 0; i < res+margin+margin+2; ++i) {
         PointDouble  pt;
      }
   }
   */

      // marker content
      cv::Mat temp(res, res, CV_8U, cv::Scalar(255));
      marker_content = temp;
      } // end if

   return;

   }


/*****************************************************************************
 ***  class MarkerArtoolkit
 ****************************************************************************/

bool MarkerArtoolkit::DecodeContent(int& orientation)
   {
   int a = marker_content.at<uchar>(0, 0);
   int b = marker_content.at<uchar>(res - 1, 0);
   int c = marker_content.at<uchar>(res - 1, res - 1);
   int d = marker_content.at<uchar>(0, res - 1);

   if (!a && !b && c)
      orientation = 0;
   else if (!b && !c && d)
      orientation = 1;
   else if (!c && !d && a)
      orientation = 2;
   else if (!d && !a && b)
      orientation = 3;
   else
      return (false);

   Bitset bs;
   bs.clear();
   for (int j = 0 ; j < res ; j++)
      {
      for (int i = 0 ; i < res ; i++)
         {
         if (orientation == 0)
            {
            if ((j == 0)     && (i == 0))
               continue;
            if ((j == res - 1) && (i == 0))
               continue;
            if ((j == res - 1) && (i == res - 1))
               continue;
            if (marker_content.at<uchar>(j, i))
               bs.push_back(false);
            else bs.push_back(true);
            }
         else if (orientation == 1)
            {
            if (((res - i - 1) == res - 1) && (j == 0))
               continue;
            if (((res - i - 1) == res - 1) && (j == res - 1))
               continue;
            if (((res - i - 1) == 0)     && (j == res - 1))
               continue;
            if (marker_content.at<uchar>(res-i-1, j))
               bs.push_back(false);
            else
               bs.push_back(true);
            }
         else if (orientation == 2)
            {
            if (((res - j - 1) == res - 1) && ((res - i - 1) == res - 1))
               continue;
            if (((res - j - 1) == 0)     && ((res - i - 1) == res - 1))
               continue;
            if (((res - j - 1) == 0)     && ((res - i - 1) == 0))
               continue;
            if (marker_content.at<uchar>(res - j - 1, res - i - 1))
               bs.push_back(false);
            else
               bs.push_back(true);
            }
         else if (orientation == 3)
            {
            if ((i == 0)     && ((res - j - 1) == res - 1))
               continue;
            if ((i == 0)     && ((res - j - 1) == 0))
               continue;
            if ((i == res - 1) && ((res - j - 1) == 0))
               continue;
            if (marker_content.at<uchar>(i, res - j - 1))
               bs.push_back(false);
            else
               bs.push_back(true);
            }
         }
      }

   id = bs.ulong();

   return (true);

   }

void MarkerArtoolkit::SetContent(unsigned long _id)
   {
   // Fill in the content values
   margin_error = 0;
   decode_error = 0;
   id = _id;
   Bitset bs;
   bs.push_back_meaningful(_id);

   for (int j = res - 1 ; j >= 0 ; j--)
      {
      for (int i = res - 1 ; i >= 0 ; i--)
         {
         if ((j == 0) && (i == 0))
            marker_content.at<uchar>(j, i) = 0;
         else if ((j == res-1) && (i == 0))
            marker_content.at<uchar>(j, i) = 0;
         else if ((j == res - 1) && (i == res - 1))
            marker_content.at<uchar>(j, i) = 255;
         else
            {
            if (bs.Length() && bs.pop_back())
               marker_content.at<uchar>(j, i) = 0;
            else
               marker_content.at<uchar>(j, i) = 255;
            }
         }
      }

   return;

   }

/*****************************************************************************
 ***  class MarkerData
 ****************************************************************************/

void MarkerData::DecodeOrientation(int& error, int& total, int& orientation)
   {
   vector<int> errors(4, 0);
   int color = 255;

   // Resolution identification
   int j = res / 2;
   for (int i = 0 ; i < res ; i++)
      {
      total++;
      if (marker_content.at<uchar>(j, i) != color)
         errors[0]++;
      if (marker_content.at<uchar>(i, j) != color)
         errors[1]++;

      color = (color ? 0 : 255);
      }

   errors[2] = errors[0];
   errors[3] = errors[1];

   // Orientation identification
   int i = res / 2;
   for (j = (res / 2) - 2 ; j <= (res / 2) + 2 ; j++)
      {
      if (j < (res / 2))
         {
         total++;
         if (marker_content.at<uchar>(j, i) != 0)
            errors[0]++;
         if (marker_content.at<uchar>(i, j) != 0)
            errors[1]++;
         if (marker_content.at<uchar>(j, i) != 255)
            errors[2]++;
         if (marker_content.at<uchar>(i, j) != 255)
            errors[3]++;
         }
      else if (j > (res / 2))
         {
         total++;
         if (marker_content.at<uchar>(j, i) != 255)
            errors[0]++;
         if (marker_content.at<uchar>(i, j) != 255)
            errors[1]++;
         if (marker_content.at<uchar>(j, i) != 0)
            errors[2]++;
         if (marker_content.at<uchar>(i, j) != 0)
            errors[3]++;
         }
      }

   orientation = min_element(errors.begin(), errors.end()) - errors.begin();
   error = errors[orientation];
   //*orientation = 0; // ttehop

   return;

   }

bool MarkerData::DetectResolution(vector<PointDouble>& _marker_corners_img, cv::Mat& gray, Camera& cam)
   {
   vector<PointDouble> marker_corners_img_undist;
   marker_corners_img_undist.resize(_marker_corners_img.size());
   copy(_marker_corners_img.begin(), _marker_corners_img.end(), marker_corners_img_undist.begin());

   // line_points
   std::vector<PointDouble> line_points;
   PointDouble pt;
   line_points.clear();
   pt.x = 0;
   pt.y = 0;
   line_points.push_back(pt);

   pt.x = -0.5 * edge_length;
   pt.y = 0;
   line_points.push_back(pt);

   pt.x = +0.5 * edge_length;
   pt.y = 0;
   line_points.push_back(pt);

   pt.x = 0;
   pt.y = -0.5 * edge_length;
   line_points.push_back(pt);

   pt.x = 0;
   pt.y = +0.5 * edge_length;
   line_points.push_back(pt);

   // Figure out the marker point position in the image
   // TODO: Note that line iterator cannot iterate outside image
   //       therefore we need to distort the endpoints and iterate straight lines.
   //       Right way would be to iterate undistorted lines and distort line points.
   Homography H;
   vector<PointDouble> line_points_img(line_points.size());
   line_points_img.resize(line_points.size());

   cam.Undistort(marker_corners_img_undist);

   H.Find(marker_corners, marker_corners_img_undist);
   H.ProjectPoints(line_points, line_points_img);

   cam.Distort(line_points_img);

   // Now we have undistorted line end points
   // Find lines and then distort the whole line
   int white_count[4] = { 0, 0, 0, 0 }; // white counts for lines 1->0, 2->0, 3->0, 4->0

   cv::Point2i pt2 = line_points_img[0];
   if ((pt2.x < 0) || (pt2.y < 0) || (pt2.x >= gray.cols) || (pt2.y >= gray.rows))
      {
      return (false);
      }

   bool white = true;
   for (int i = 0 ; i < 4 ; i++)
      {
      cv::Point2i pt1 = line_points_img[i + 1];
      if ((pt1.x < 0) || (pt1.y < 0) || (pt1.x >= gray.cols) || (pt1.y >= gray.rows))
         {
         return (false);
         }

      cv::LineIterator iterator(gray, pt1, pt2, 8, false);
      std::vector<uchar> vals;
      for (int ii = 0 ; ii < iterator.count ; ii++)
         {
         iterator++;
         vals.push_back(*(*iterator));
         }

      uchar vmin = *(std::min_element(vals.begin(), vals.end()));
      uchar vmax = *(std::max_element(vals.begin(), vals.end()));
      uchar thresh = (vmin + vmax) / 2;
      white = true;
      int bc = 0, wc = 0, N = 2;
      for (uchar& val : vals)
         {
         // change the color status if we had
         // N subsequent pixels of the other color
         if (val < thresh)
            {
            bc++;
            wc = 0;
            }
         else
            {
            wc++;
            bc = 0;
            }

         if (white && (bc >= N))
            {
            white = false;
            }
         else if (!white && (wc >= N))
            {
            white = true;
            white_count[i]++;
            }
         }
      }

   bool ret = true;
   if ((white_count[0] + white_count[1]) == (white_count[2] + white_count[3]))
      {
      ret = false;
      }
   else if ((white_count[0] + white_count[1]) > (white_count[2] + white_count[3]))
      {
      if ((white_count[0] != white_count[1]) || (white_count[0] < 2))
         {
         ret = false;
         } // end if
      else
         {
         int nof_whites = white_count[0] * 2 - (white ? 1 : 0); // 'white' contains middle color
         int new_res = 2 * nof_whites - 1;
         SetMarkerSize(edge_length, new_res, margin);
         } // end else
      }
   else
      {
      if ((white_count[2] != white_count[3]) || (white_count[2] < 2) || (((white_count[2] % 2) == 0) != white))
         {
         ret = false;
         } // end if
      else
         {
         int nof_whites = white_count[2] * 2 - (white ? 1 : 0);
         int new_res = 2 * nof_whites - 1;
         SetMarkerSize(edge_length, new_res, margin);
         } // end else
      }

   return (ret);

   }

bool MarkerData::UpdateContent(vector<PointDouble>& _marker_corners_img, cv::Mat& gray, Camera& cam, int frame_no /*= 0*/)
   {
   bool ret = false;

   if (res == 0)
      {
      ret = DetectResolution(_marker_corners_img, gray, cam);
      if (ret)
         {
         ret = UpdateContentBasic(_marker_corners_img, gray, cam, frame_no);
         } // end if
      } // end if
   else
      {
      ret = UpdateContentBasic(_marker_corners_img, gray, cam, frame_no);
      } // end else

   return (ret);

   }

int MarkerData::DecodeCode(int orientation, BitsetExt& bs, int& erroneous, int& total,
                           unsigned char& content_type)
   {
   // TODO: The orientation isn't fully understood?
   //for (int j = res-1; j >= 0; j--) {
   for (int j = 0 ; j < res ; j++)
      {
      for (int i = 0 ; i < res ; i++)
         {
         // TODO: Does this work ok for larger markers?
         if ((orientation == 0) || (orientation == 2))
            {
            if (j == res / 2)
               continue;
            if ((i == res / 2) && (j >= (res / 2) - 2) && (j <= (res / 2) + 2))
               continue;
            }
         else
            {
            if (i == res / 2)
               continue;
            if ((j == res / 2) && (i >= (res / 2) - 2) && (i <= (res / 2) + 2))
               continue;
            }

         int color = 0;
         if (orientation == 0)
            color = marker_content.at<uchar>(j, i);
         else if (orientation == 1)
            color = marker_content.at<uchar>(res - i - 1, j);
         else if (orientation == 2)
            color = marker_content.at<uchar>(res - j - 1, res - i - 1);
         else if (orientation == 3)
            color = marker_content.at<uchar>(i, res - j - 1);

         bs.push_back((color == 0));

         total++;
         }
      }

   unsigned char flags = 0;
   int errors = 0;

   // if we have larger than 16-bit code, then we have a header; 16-bit code has a
   // hamming(8,4) coded number
   if (bs.Length() > 16)
      {
      // read header (8-bit hamming(8,4) -> 4-bit flags)
      BitsetExt header;

      for (int i = 0 ; i < HEADER_SIZE ; i++)
         header.push_back(bs.pop_front());

      errors = header.hamming_dec(8);
      if (errors == -1)
         {
         //OutputDebugString("header decoding failed!!!!!\n");
         return (errors);
         }

      flags = header.uchar();
      }
   else
      {
      flags &= MarkerContentType::MARKER_CONTENT_TYPE_NUMBER;
      }

   // check which hamming we are using
   //bs->Output(cout); cout<<endl;
   if (flags & 0x8)
      errors = bs.hamming_dec(16);
   else
      errors = bs.hamming_dec(8);

   content_type = flags & 0x7;

   if (errors > 0)
      erroneous += errors;

   return (errors);

   }

void MarkerData::Read6bitStr(const BitsetExt& bs, char* s, size_t s_max_len)
   {
   const deque<bool> bits = bs.GetBits();
   deque<bool>::const_iterator iter;
   size_t len = 0;
   int bitpos = 5;
   unsigned long c = 0;

   for (iter = bits.begin() ; iter != bits.end() ; iter++)
      {
      if (*iter)
         c |= (0x01 << bitpos);

      bitpos--;

      if (bitpos < 0)
         {
         if (c == 000)                      s[len] = ':';
         else if ((c >= 001) && (c <= 032)) s[len] = 'a' + (char)c - 1;
         else if ((c >= 033) && (c <= 044)) s[len] = '0' + (char)c - 1;
         else if (c == 045)                 s[len] = '+';
         else if (c == 046)                 s[len] = '-';
         else if (c == 047)                 s[len] = '*';
         else if (c == 050)                 s[len] = '/';
         else if (c == 051)                 s[len] = '(';
         else if (c == 052)                 s[len] = ')';
         else if (c == 053)                 s[len] = '$';
         else if (c == 054)                 s[len] = '=';
         else if (c == 055)                 s[len] = ' ';
         else if (c == 056)                 s[len] = ',';
         else if (c == 057)                 s[len] = '.';
         else if (c == 060)                 s[len] = '#';
         else if (c == 061)                 s[len] = '[';
         else if (c == 062)                 s[len] = ']';
         else if (c == 063)                 s[len] = '%';
         else if (c == 064)                 s[len] = '\"';
         else if (c == 065)                 s[len] = '_';
         else if (c == 066)                 s[len] = '!';
         else if (c == 067)                 s[len] = '&';
         else if (c == 070)                 s[len] = '\'';
         else if (c == 071)                 s[len] = '?';
         else if (c == 072)                 s[len] = '<';
         else if (c == 073)                 s[len] = '>';
         else if (c == 074)                 s[len] = '@';
         else if (c == 075)                 s[len] = '\\';
         else if (c == 076)                 s[len] = '^';
         else if (c == 077)                 s[len] = ';';
         else s[len] = '?';

         len++;
         if (len >= (s_max_len - 1))
            break;

         bitpos = 5;
         c = 0;
         }
      }

   s[len] = 0;

   return;

   }

bool MarkerData::DecodeContent(int& orientation)
   {
   bool ret = true;

   orientation = 0;

   BitsetExt bs;
   int erroneous = 0;
   int total = 0;
   DecodeOrientation(erroneous, total, orientation);

   int err = DecodeCode(orientation, bs, erroneous, total, content_type);
   if (err != -1)
      {
      if (content_type == MarkerContentType::MARKER_CONTENT_TYPE_NUMBER)
         {
         data.id = bs.ulong();
         }
      else
         {
         Read6bitStr(bs, data.str, MAX_MARKER_STRING_LEN);
         }

      decode_error = static_cast<double>(erroneous) / total;
      }
   else
      {
      // couldn't fix
      decode_error = DBL_MAX;
      ret = false;
      } // end else

   return (ret);

   }

void MarkerData::Add6bitStr(BitsetExt& bs, char* s)
   {
   while (*s)
      {
      unsigned char c = (unsigned char)*s;
      if (c == ':')                      bs.push_back((unsigned char)0, 6);
      else if ((c >= 'A') && (c <= 'Z')) bs.push_back((unsigned char)(001 + c - 'A'), 6);
      else if ((c >= 'a') && (c <= 'z')) bs.push_back((unsigned char)(001 + c - 'a'), 6);
      else if ((c >= '0') && (c <= '9')) bs.push_back((unsigned char)(033 + c - '0'), 6);
      else if (c == '+')                 bs.push_back((unsigned char)045, 6);
      else if (c == '-')                 bs.push_back((unsigned char)046, 6);
      else if (c == '*')                 bs.push_back((unsigned char)047, 6);
      else if (c == '/')                 bs.push_back((unsigned char)050, 6);
      else if (c == '(')                 bs.push_back((unsigned char)051, 6);
      else if (c == ')')                 bs.push_back((unsigned char)052, 6);
      else if (c == '$')                 bs.push_back((unsigned char)053, 6);
      else if (c == '=')                 bs.push_back((unsigned char)054, 6);
      else if (c == ' ')                 bs.push_back((unsigned char)055, 6);
      else if (c == ',')                 bs.push_back((unsigned char)056, 6);
      else if (c == '.')                 bs.push_back((unsigned char)057, 6);
      else if (c == '#')                 bs.push_back((unsigned char)060, 6);
      else if (c == '[')                 bs.push_back((unsigned char)061, 6);
      else if (c == ']')                 bs.push_back((unsigned char)062, 6);
      else if (c == '%')                 bs.push_back((unsigned char)063, 6);
      else if (c == '\"')                bs.push_back((unsigned char)064, 6);
      else if (c == '_')                 bs.push_back((unsigned char)065, 6);
      else if (c == '!')                 bs.push_back((unsigned char)066, 6);
      else if (c == '&')                 bs.push_back((unsigned char)067, 6);
      else if (c == '\'')                bs.push_back((unsigned char)070, 6);
      else if (c == '?')                 bs.push_back((unsigned char)071, 6);
      else if (c == '<')                 bs.push_back((unsigned char)072, 6);
      else if (c == '>')                 bs.push_back((unsigned char)073, 6);
      else if (c == '@')                 bs.push_back((unsigned char)074, 6);
      else if (c == '\\')                bs.push_back((unsigned char)075, 6);
      else if (c == '^')                 bs.push_back((unsigned char)076, 6);
      else if (c == ';')                 bs.push_back((unsigned char)077, 6);
      else                               bs.push_back((unsigned char)071, 6);

      s++;
      }

   return;

   }

int MarkerData::UsableDataBits(int marker_res, int hamming)
   {
   int bits = 0;

   // Marker resolution must be 5 or greater and odd
   if ((marker_res >= 5) && ((marker_res % 2)) == 1)
      {
      bits = marker_res * marker_res;
      if (marker_res > 5)
         bits -= 8; // With larger resolutions we reserve 8 bits for hamming(8,4) encoded 4 flags

      bits -= marker_res;            // center line indicating the resolution
      bits -= 4;                     // the four pixels indicating the orientation
      int tail = bits % hamming;
      if (tail < 3)
         bits -= tail;    // hamming can't use tail pixels if there is only 2 or 1 of them
      } // end if

   return (bits);

   }

void MarkerData::SetContent(MarkerContentType _content_type, unsigned long _id, const char* _str,
      bool force_strong_hamming, bool verbose)
   {
   // Fill in the content values
   content_type = _content_type;
   margin_error = 0;
   decode_error = 0;

   if (content_type == MarkerContentType::MARKER_CONTENT_TYPE_NUMBER)
      {
      data.id = _id;
      }
   else
      {
      strncpy(data.str, _str, MAX_MARKER_STRING_LEN - 1);
      }

   // Encode
   const int max_marker_res = 127;
   BitsetExt bs_flags(verbose);
   BitsetExt bs_data(verbose);
   int enc_bits;  // How many encoded bits fits in the marker
   int data_bits; // How many data bits fit inside the encoded bits
   int hamming;   // Do we use 8-bit or 16-bit hamming?
   if (content_type == MarkerContentType::MARKER_CONTENT_TYPE_NUMBER)
      {
      bs_data.push_back_meaningful(data.id);

      bool loop = true;
      for (res = 5 ; (res < max_marker_res) && loop ; res += 2)
         {
         hamming = 8;
         enc_bits = UsableDataBits(res, hamming);
         data_bits = BitsetExt::count_hamming_dec_len(hamming, enc_bits);
         if (data_bits < bs_data.Length())
            {
            if ((res > 5) && !force_strong_hamming)
               {
               hamming = 16;
               enc_bits = UsableDataBits(res, hamming);
               data_bits = BitsetExt::count_hamming_dec_len(hamming, enc_bits);
               if (data_bits < bs_data.Length())
                  {
                  loop = false;
                  } // end if
               }
            } // end if
         else
            {
            loop = false;
            } // end else
         }  // end for

      bs_data.fill_zeros_left(data_bits);
      bs_data.hamming_enc(hamming);
      if (verbose)
         {
         cout << "Using hamming(" << hamming << ") for " << res << "x" << res << " marker" << endl;
         cout << bs_data.Length() << " bits are filled into " << data_bits;
         cout << " bits, and encoded into " << enc_bits << " bits" << endl;
         cout << "data src: ";
         bs_data.Output(cout);
         cout << endl;
         cout << "data enc: ";
         bs_data.Output(cout);
         cout << endl;
         }

      if (res > 5)
         {
         if (hamming == 16)
            bs_flags.push_back(true);
         else
            bs_flags.push_back(false);

         bs_flags.push_back(0UL, 3);
         bs_flags.hamming_enc(8);
         if (verbose)
            {
            cout << "flags src: ";
            bs_flags.Output(cout);
            cout << endl;
            cout << "flags enc: ";
            bs_flags.Output(cout);
            cout << endl;
            }
         }
      }
   else
      {
      Add6bitStr(bs_data, data.str);
      bool loop = true;
      for (res = 7 ; (res < max_marker_res) && loop ; res +=2)
         {
         hamming = 8;
         enc_bits = UsableDataBits(res, hamming);
         data_bits = BitsetExt::count_hamming_dec_len(hamming, enc_bits);
         if (data_bits < bs_data.Length())
            {
            if (!force_strong_hamming)
               {
               hamming = 16;
               enc_bits = UsableDataBits(res, hamming);
               data_bits = BitsetExt::count_hamming_dec_len(hamming, enc_bits);
               loop = (data_bits < bs_data.Length());
               }
            } // end if
         else
            {
            loop = false;
            } // end else
         }

      while (bs_data.Length() < ((data_bits / 6) * 6))
         {
         bs_data.push_back((unsigned char)055, 6); // Add space
         }

      while (bs_data.Length() < data_bits)
         {
         bs_data.push_back(false); // Add 0
         }

      bs_data.hamming_enc(hamming);
      if (hamming == 16)
         bs_flags.push_back(true);
      else
         bs_flags.push_back(false);

      if (content_type == MarkerContentType::MARKER_CONTENT_TYPE_STRING)
         bs_flags.push_back(1UL, 3);
      else if (content_type == MarkerContentType::MARKER_CONTENT_TYPE_FILE)
         bs_flags.push_back(2UL, 3);
      else if (content_type == MarkerContentType::MARKER_CONTENT_TYPE_HTTP)
         bs_flags.push_back(3UL, 3);
      bs_flags.hamming_enc(8);

      if (verbose)
         {
         cout << "Using hamming(" << hamming << ") for " << res << "x" << res << " marker" << endl;
         cout << bs_data.Length() << " bits are filled into " << data_bits;
         cout << " bits, and encoded into " << enc_bits << " bits" << endl;
         cout << "data src: ";
         bs_data.Output(cout);
         cout << endl;
         cout << "data enc: ";
         bs_data.Output(cout);
         cout << endl;
         cout << "flags src: ";
         bs_flags.Output(cout);
         cout << endl;
         cout << "flags enc: ";
         bs_flags.Output(cout);
         cout << endl;
         }
      }

   // Fill in the marker content
   deque<bool> bs(bs_flags.GetBits());
   bs.insert(bs.end(), bs_data.GetBits().begin(), bs_data.GetBits().end());
   deque<bool>::const_iterator iter = bs.begin();
   SetMarkerSize(edge_length, res, margin);
   marker_content = cv::Scalar(255);
   for (int j = 0 ; j < res ; j++)
      {
      for (int i = 0 ; i < res ; i++)
         {
         if (j == (res / 2))
            {
            if ((i % 2) != 0)
               {
               marker_content.at<uchar>(j, i) = 0;
               }
            }
         else if ((i == (res / 2)) && (j < (res / 2)) && (j >= ((res / 2) - 2)))
            {
            marker_content.at<uchar>(j, i) = 0;
            }
         else if ((i == (res / 2)) && (j > (res / 2)) && (j <= ((res / 2) + 2)))
            {
            marker_content.at<uchar>(j, i) = 255;
            }
         else
            {
            if (iter != bs.end())
               {
               if (*iter)
                  marker_content.at<uchar>(j, i) = 0;
               iter++;
               }
            }
         }
      }

   return;

   }

} // end namespace alvar
