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

#include "ConnectedComponents.h"
#include "Draw.h"
#include <cassert>

using namespace std;

namespace alvar {

using namespace std;

/*****************************************************************************
 ***  class Labeling
 ****************************************************************************/

Labeling::Labeling()
   {
   thresh_param1 = 31;
   thresh_param2 = 5;

   return;

   }

bool Labeling::CheckBorder(const std::vector<cv::Point2i>& contour, int width, int height) const
   {
   bool ret = true;
   for (auto& pt : contour)
      {
      if ((pt.x <= 1) || (pt.x >= width - 2) || (pt.y <= 1) || (pt.y >= height - 2))
         ret = false;
      }

   return (ret);

   }

/*****************************************************************************
 ***  class LabelingCvSeq
 ****************************************************************************/

LabelingCvSeq::LabelingCvSeq() : _n_blobs(0), _min_edge(20), _min_area(25)
   {
   SetOptions();

   return;

   }

void LabelingCvSeq::SetOptions(bool _detect_pose_grayscale)
   {
   detect_pose_grayscale = _detect_pose_grayscale;

   return;

   }

void LabelingCvSeq::LabelSquares(cv::Mat& image, Camera* pCamera /* = nullptr */,
      bool visualize /* = false */)
   {
   // Convert grayscale and threshold
   if (image.channels() == 4)
      {
      cv::cvtColor(image, gray, CV_RGBA2GRAY);
      } // end if
   else if(image.channels() == 3)
      {
      cv::cvtColor(image, gray, CV_RGB2GRAY);
      }  // end else if
   else if(image.channels() == 1)
      {
      gray = image.clone();
      }  // end else if
   else
      {
      cerr << "Unsupported image format" << endl;
      // Maybe want to notify the caller???  TTC!!
      }

   cv::adaptiveThreshold(gray, bw, 255, cv::ADAPTIVE_THRESH_MEAN_C,
         cv::THRESH_BINARY_INV, thresh_param1, thresh_param2);

   std::vector<std::vector<cv::Point2i>> squares;
   std::vector<std::vector<cv::Point2i>> square_contours;

   // Find all the closed contours in the image
   std::vector<std::vector<cv::Point2i>> contours;
   cv::findContours(bw, contours, cv::RETR_LIST, cv::CHAIN_APPROX_NONE);

   // Find contours that are possibly squares
   std::vector<std::vector<cv::Point2i>>::iterator contour = contours.begin();
   while (contour != contours.end())
      {
      if (static_cast<int>(contour->size()) >= _min_edge)
         {
         std::vector<cv::Point2i> poly;
         cv::approxPolyDP(*contour, poly, contour->size() * 0.035, true);

         if ((poly.size() == 4) && CheckBorder(poly, image.cols, image.rows) &&
            (fabs(cv::contourArea(poly, false)) > _min_area) &&
            cv::isContourConvex(poly))
            {
            squares.push_back(poly);
            square_contours.push_back(*contour);
            } // end if
         } // end if

      ++contour;
      } // end while

   _n_blobs = squares.size();
   blob_corners.resize(_n_blobs);

   // For every detected 4-corner blob
   for (int i = 0 ; i < _n_blobs ; ++i)
      {
      vector<Line> fitted_lines(4);
      blob_corners[i].resize(4);
      std::vector<cv::Point2i>& sq = squares[i];
      std::vector<cv::Point2i>& square_contour = square_contours[i];

      for (int j = 0 ; j < 4 ; ++j)
         {
         cv::Point2i& pt0 = sq[j];
         cv::Point2i& pt1 = sq[(j + 1) % 4];
         int k0 = -1, k1 = -1;
         for (size_t k = 0 ; k < square_contour.size() ; k++)
            {
            cv::Point2i& pt2 = square_contour[k];
            if ((pt0.x == pt2.x) && (pt0.y == pt2.y))
               k0 = k;

            if ((pt1.x == pt2.x) && (pt1.y == pt2.y))
               k1 = k;
            }

         int len = 1;
         if (k1 >= k0)
            len = k1 - k0 - 1; // neither k0 nor k1 are included
         else
            len = square_contour.size() - k0 + k1 - 1;

         if (len == 0)
            len = 1;

         cv::Mat_<cv::Point2f> line_data(1, len);
         for (int l = 0 ; l < len ; l++)
            {
            int ll = (k0 + l + 1) % square_contour.size();
            cv::Point2f pp = square_contour[ll];

            // Undistort
            if (pCamera != nullptr)
               {
               pCamera->Undistort(pp);
               } // end if

            line_data(0, l) = pp;
            }

         // Fit edge and put to vector of edges

         // TODO: The detect_pose_grayscale is still under work...
         /*
            if (detect_pose_grayscale &&
                (pt0->x > 3) && (pt0->y > 3) &&
                (pt0->x < (gray->width-4)) &&
                (pt0->y < (gray->height-4)))
            {
                // ttehop: Grayscale experiment
                FitLineGray(line_data, params, gray);
            }
            */

         cv::Vec4f params;
         cv::fitLine(line_data, params, CV_DIST_L2, 0, 0.01, 0.01);
         Line line(params);
         if (visualize)
            DrawLine(image, line);

         fitted_lines[j] = line;
         }

      // Calculated four intersection points
      for (size_t j = 0 ; j < 4 ; ++j)
         {
         PointDouble intc = Line::Intersection(fitted_lines[j], fitted_lines[(j + 1) % 4]);

         // TODO: Instead, test OpenCV find corner in sub-pix...
         //cv::Point2f pt = cv::Point2f(intc.x, intc.y);
         //cvFindCornerSubPix(gray, &pt,
         //                   1, cvSize(3,3), cvSize(-1,-1),
         //                   cvTermCriteria(
         //                   CV_TERMCRIT_ITER+CV_TERMCRIT_EPS,10,1e-4));

         // TODO: Now there is a wierd systematic 0.5 pixel error that is fixed here...
         //intc.x += 0.5;
         //intc.y += 0.5;

         if (pCamera != nullptr)
            {
            pCamera->Distort(intc);
            } // end if

         // TODO: Should we make this always counter-clockwise or clockwise?
         /*
            if (image->origin && j == 1) blob_corners[i][3] = intc;
            else if (image->origin && j == 3) blob_corners[i][1] = intc;
            else blob_corners[i][j] = intc;
            */

         blob_corners[i][j] = intc;
         }

      if (visualize)
         {
         std::vector<cv::Scalar> color =
            {
            {255, 255, 255},
            {0, 0, 255},
            {0, 255, 0},
            {255, 0, 0}
            };

         for(size_t j = 0 ; j < 4 ; ++j)
            {
            cv::circle(image, blob_corners[i][j], 5, color[j]);
            }
         }
      }

   return;

   }

#if 0
// Not called currently TTC!!
void LabelingCvSeq::LabelImage(cv::Mat& image,
      std::vector<std::vector<cv::Point2i>>& squares, int min_size, bool approx)
   {
   squares.clear();

   // Convert grayscale and threshold
   if(image.channels() == 4)
      cv::cvtColor(image, gray, CV_RGBA2GRAY);
   else if(image.channels() == 3)
      cv::cvtColor(image, gray, CV_RGB2GRAY);
   else if(image.channels() == 1)
      {
      gray = image.clone();
      cerr << "Unsupported image format" << endl;
      }

   cv::adaptiveThreshold(gray, bw, 255, cv::ADAPTIVE_THRESH_MEAN_C, cv::THRESH_BINARY_INV, thresh_param1, thresh_param2);

   std::vector<std::vector<cv::Point2i>> contours;
   cv::findContours(bw, contours, cv::RETR_LIST, cv::CHAIN_APPROX_NONE);

   std::vector<std::vector<cv::Point2i>>::iterator contour = contours.begin();
   while (contour != contours.end())
      {
      if (static_cast<int>(contour->size()) >= min_size)
         {
         if (approx)
            {
            std::vector<cv::Point2i> poly;
            cv::approxPolyDP(*contour, poly,
                  contour->size() * 0.02, 0 ); // TODO: Parameters?

            if (isContourConvex(poly))
               {
               squares.push_back(poly);
               }
            }
         else
            {
            squares.push_back(*contour);
            }
         }
      }

   return;

   }
#endif

//#define SHOW_DEBUG
#ifdef SHOW_DEBUG
#include "highgui.h"
#endif

#if 0
// TODO: This should be in LabelingCvSeq ???
// As near as I can tell, this function does squat in this form!!  TTC!!
// And it's only call is commented out.  TTC!!
void FitLineGray(cv::Mat_<cv::Point2f>& line_data, float /* params */ [4], cv::Mat& gray)
   {
	// this very simple approach works...
	/*
	float *cx = &(params[2]);
	float *cy = &(params[3]);
	float *sx = &(params[0]);
	float *sy = &(params[1]);
   cv::Point2f *p1 = (cv::Point2f*)CV_MAT_ELEM_PTR_FAST(*line_data, 0, 0, sizeof(cv::Point2f));
   cv::Point2f *p2 = (cv::Point2f*)CV_MAT_ELEM_PTR_FAST(*line_data, 0, line_data->cols-1, sizeof(cv::Point2f));
	*cx = p1->x; *cy = p1->y;
	*sx = p2->x - p1->x; *sy = p2->y - p1->y;
	return;
	*/

#ifdef SHOW_DEBUG
	IplImage *tmp = cvCreateImage(cvSize(gray->width, gray->height), IPL_DEPTH_8U, 3);
	IplImage *tmp2 = cvCreateImage(cvSize(gray->width*5, gray->height*5), IPL_DEPTH_8U, 3);
	cvCvtColor(gray, tmp, CV_GRAY2RGB);
	cvResize(tmp, tmp2, CV_INTER_NN);
#endif

	// Discover 1st the line normal direction
   cv::Point2f p1 = line_data(0, 0);
   cv::Point2f p2 = line_data(0, line_data.cols - 1);
   double dx = +(p2.y - p1.y);
   double dy = -(p2.x - p1.x);
   if ((dx == 0) && (dy == 0))
      // Fails but doesn't notify TTC!!
      return;
   else if (dx == 0)
      {
      dy /= dy;
      }
   else if (dy == 0)
      {
      dx /= dx;
      }
   else if (abs(dx) > abs(dy))
      {
      dy /= dx;
      dx /= dx;
      }
   else
      {
      dx /= dy;
      dy /= dy;
      }

	// Build normal search table
   const int win_size = 5;
   const int win_mid = win_size / 2;
   const int diff_win_size = win_size - 1;
	double xx[win_size], yy[win_size];
	double dxx[diff_win_size], dyy[diff_win_size];
   xx[win_mid] = 0;
   yy[win_mid] = 0;
   for (int i = 1 ; i <= win_size / 2 ; i++)
      {
      xx[win_mid + i] = std::round(i * dx);
		xx[win_mid - i] = -xx[win_mid + i];
      yy[win_mid + i] = std::round(i * dy);
		yy[win_mid - i] = -yy[win_mid + i];
      }

   for (int i = 0 ; i < diff_win_size ; i++)
      {
      dxx[i] = (xx[i] + xx[i + 1]) / 2;
      dyy[i] = (yy[i] + yy[i + 1]) / 2;
      }

	// Adjust the points
   for (int l = 0 ; l < line_data.cols ; l++)
      {
      cv::Point2f& p = line_data(0, l);

      double dx = 0, dy = 0, ww = 0;
      for (int i = 0 ; i < diff_win_size ; i++)
         {
         double c1 = gray.at<uchar>(int(p.y + yy[i]), int(p.x + xx[i]));
         double c2 = gray.at<uchar>(int(p.y + yy[i + 1]), int(p.x + xx[i + 1]));
#ifdef SHOW_DEBUG
         cv::circle(tmp2, cv::Point2i((p->x+xx[i])*5+2,(p->y+yy[i])*5+2), 0, CV_RGB(0,0,255));
         cv::circle(tmp2, cv::Point2i((p->x+xx[i+1])*5+2,(p->y+yy[i+1])*5+2), 0, CV_RGB(0,0,255));
#endif
         double w = std::abs(c1 - c2);
         dx += dxx[i] * w;
         dy += dyy[i] * w;
         ww += w;
         }

      if (ww > 0)
         {
         dx /= ww;
         dy /= ww;
         }

#ifdef SHOW_DEBUG
      cvLine(tmp2, cv::Point2i(p->x*5+2,p->y*5+2), cv::Point2i((p->x+dx)*5+2, (p->y+dy)*5+2), CV_RGB(0,255,0));
      p->x += float(dx);
      p->y += float(dy);
      cvCircle(tmp2, cv::Point2i(p->x*5+2,p->y*5+2), 0, CV_RGB(255,0,0));
#else
      p.x += float(dx);
      p.y += float(dy);
#endif
      }

#ifdef SHOW_DEBUG
	cvNamedWindow("tmp");
   cvShowImage("tmp", tmp2);
	cvWaitKey(0);
	cvReleaseImage(&tmp);
	cvReleaseImage(&tmp2);
#endif

   return;

   }
#endif
} // namespace alvar
