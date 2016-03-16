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

#include "Draw.h"
#include <cassert>

using namespace std;

namespace alvar {
using namespace std;

void DrawPoints(cv::Mat& image, const vector<cv::Point2i>& points, cv::Scalar color)
   {
   for (auto& point : points)
      cv::line(image, point, point, color);

   return;

   }

void DrawLine(cv::Mat& image, const Line& line, cv::Scalar color, double len /* = 200.0 */)
   {
   cv::Point2d Center(line.GetCenter());
   cv::Point2d Direction(line.GetDirection());
   cv::Point2i p1 = (Center + ((len / 2) * Direction));
   cv::Point2i p2 = (Center - ((len / 2) * Direction));

   cv::line(image, p1, p2, color);

   return;

   }

#if 0
// Not currently used
void DrawPoints(cv::Mat& image, const CvSeq* contour, cv::Scalar color)
{
	for(int i = 0; i < contour->total; ++i)
	{
      cv::Point2i* pt = (cv::Point2i*)cvGetSeqElem( contour, i);
      cv::line(image, *pt, *pt, color);
	}
}

void DrawCircles(cv::Mat& image, const CvSeq* contour, int radius, cv::Scalar color)
   {
   for (int i = 0 ; i < contour->total ; ++i)
      {
      cv::Point2i* pt = (cv::Point2i*)cvGetSeqElem(contour, i);
      cv::circle(image, *pt, radius, color);
      }

   return;

   }

void DrawLines(cv::Mat& image, const CvSeq* contour, cv::Scalar color)
   {
   if(contour->total >= 2)
      {
      for (size_t i = 0 ; i < contour->total ; ++i)`
         {
         cv::Point2i* pt1 = (cv::Point2i*)cvGetSeqElem( contour, i);
         cv::Point2i* pt2 = (cv::Point2i*)cvGetSeqElem( contour, (i+1)%(contour->total));
         cv::line(image, *pt1, *pt2, color);
         }
      }

   return;

   }
#endif

#if 0
void DrawCVEllipse(cv::Mat& image, const cv::RotatedRect& ellipse,
      cv::Scalar color, int thickness /* = 1 */, int lineType /* = 8 */)
   {
   //cout<<center.x<<" "<<center.y<<" "<<ellipse.size.width/2<<" "<<ellipse.size.height/2<<" "<<ellipse.angle<<endl;
   cv::ellipse(image, ellipse, color, thickness, lineType);

   return;

   }
#endif
#if 0
// Not currently used TTC!!
void BuildHideTexture(cv::Mat& image, cv::Mat& hide_texture,
      Camera& cam, double gl_modelview[16],
      PointDouble topleft, PointDouble botright)
   {
   double kx = 1.0;
   double ky = 1.0;

   double width = abs(botright.x - topleft.x);
   double height = abs(botright.y - topleft.y);

   double objx;
   double objy;

   unsigned char pixels[8][3];
   //	unsigned char color[3]={0,0,0}; // TTC!!

   int i = 0,j = 0 /* ,t=0 */; // TTC!!
   double ox, oy, ya, yb, xc, xd, offset;
   double sizex = width / 4, size2x = width / 2;
   double sizey = height / 4, size2y = height / 2;

   // Calculate extended coordinates of detected marker (+ border)
   objx = width / 2 * kx;
   objy = height / 2 * ky;

   //cout<<hide_texture->width<<","<<hide_texture->height<<endl;

   double l2r = 2 * width * kx;
   //	double l2s=2*height*ky; // TTC!!
   double lr = width * kx;
   double ls = height * ky;
   double r, s;
   double xstep = 2 * objx / hide_texture.cols;
   double ystep = 2 * objy / hide_texture.rows;
   for(i = 0 ; i < hide_texture.cols ; i++)
      {
      ox = -objx + i * xstep;
      offset = fmod((objx-ox), size2x);
      if (offset < sizex)
         xc = objx + offset;
      else
         xc = objx + size2x - offset;
      offset = fmod((objx+ox), size2x);
      if (offset < sizex)
         xd = -objx - offset;
      else
         xd = -objx-size2x+offset;

      r = (ox + objx);
      for (j = 0 ; j < hide_texture.rows ; j++)
         {
         oy = -objy+j*ystep;
         offset = fmod((objy-oy), size2y);
         if(	offset < sizey)
            ya = objy + offset;
         else
            ya = objy+size2y-offset;
         offset = fmod((oy+objy), size2y);
         if( offset < sizey)
            yb = -objy - offset;
         else
            yb = -objy-size2y+offset;
         s=(oy+objy);

         double points3d[4][3] =
         {
            ox, ya, 0,
            ox, yb, 0,
            xc, oy, 0,
            xd, oy, 0,
         };
         double points2d[4][2];
         cv::Mat points3d_mat(4, 3, CV_64F, points3d);
         cv::Mat points2d_mat(4, 2, CV_64F, points2d);
         cam->ProjectPoints(points3d_mat, gl_modelview, points2d_mat);
         int kuvanx4 = (int)Limit(points2d[0][0], 0, image.cols - 1);
         int kuvany4 = (int)Limit(points2d[0][1], 0, image.rows - 1);
         int kuvanx5 = (int)Limit(points2d[1][0], 0, image.cols - 1);
         int kuvany5 = (int)Limit(points2d[1][1], 0, image.rows - 1);
         int kuvanx6 = (int)Limit(points2d[2][0], 0, image.cols - 1);
         int kuvany6 = (int)Limit(points2d[2][1], 0, image.rows - 1);
         int kuvanx7 = (int)Limit(points2d[3][0], 0, image.cols - 1);
         int kuvany7 = (int)Limit(points2d[3][1], 0, image.rows - 1);

         pixels[4][0] = (unsigned char)cvGet2D(image, kuvany4, kuvanx4).val[0];
         pixels[4][1] = (unsigned char)cvGet2D(image, kuvany4, kuvanx4).val[1];
         pixels[4][2] = (unsigned char)cvGet2D(image, kuvany4, kuvanx4).val[2];
         pixels[5][0] = (unsigned char)cvGet2D(image, kuvany5, kuvanx5).val[0];
         pixels[5][1] = (unsigned char)cvGet2D(image, kuvany5, kuvanx5).val[1];
         pixels[5][2] = (unsigned char)cvGet2D(image, kuvany5, kuvanx5).val[2];
         pixels[6][0] = (unsigned char)cvGet2D(image, kuvany6, kuvanx6).val[0];
         pixels[6][1] = (unsigned char)cvGet2D(image, kuvany6, kuvanx6).val[1];
         pixels[6][2] = (unsigned char)cvGet2D(image, kuvany6, kuvanx6).val[2];
         pixels[7][0] = (unsigned char)cvGet2D(image, kuvany7, kuvanx7).val[0];
         pixels[7][1] = (unsigned char)cvGet2D(image, kuvany7, kuvanx7).val[1];
         pixels[7][2] = (unsigned char)cvGet2D(image, kuvany7, kuvanx7).val[2];

         // make the borders of the texture partly transparent
         int opaque;
         const int w = 1;
         if ((i < w) | (j < w) | (i > hide_texture.cols - w) | (j > hide_texture.cols - w))
            opaque=60;
         else if ((i<2*w)|(j<2*w)|(i>hide_texture.cols -2*w)|(j>hide_texture.cols -2*w))
            opaque=100;
         else if ((i<3*w)|(j<3*w)|(i>hide_texture.cols -3*w)|(j>hide_texture.cols -3*w))
            opaque=140;
         else if ((i<4*w)|(j<4*w)|(i>hide_texture.cols -4*w)|(j>hide_texture.cols -4*w))
            opaque=200;
         else
            opaque=255;

         cvSet2D(hide_texture, j, i, cv::Scalar(
                    (((lr-r)*pixels[7][0] + r*pixels[6][0]+ s* pixels[4][0] + (ls-s)* pixels[5][0])/l2r),
               (((lr-r)*pixels[7][1] + r*pixels[6][1]+ s* pixels[4][1] + (ls-s)* pixels[5][1])/l2r),
               (((lr-r)*pixels[7][2] + r*pixels[6][2]+ s* pixels[4][2] + (ls-s)* pixels[5][2])/l2r),
               opaque
               ));
         }
      }
   }


void DrawTexture(cv::Mat& image, cv::Mat& texture, Camera& cam,
      double gl_modelview[16], PointDouble topleft, PointDouble botright)
   {
   double width = abs(botright.x - topleft.x);
   double height = abs(botright.y - topleft.y);
   double objx = width /  2;
   double objy = height / 2;

   // Project corners
   double points3d[4][3] =
      {
      -objx, -objy, 0,
      -objx, objy, 0,
      objx,  objy, 0,
      objx, -objy, 0,
      };

   double points2d[4][2];
   cv::Mat points3d_mat(4, 3, CV_64F, points3d);
   cv::Mat points2d_mat(4, 2, CV_64F, points2d);
   cam->ProjectPoints(points3d_mat, gl_modelview, points2d_mat);

   // Warp texture and mask using the perspective that is based on the corners
   double map[9];
   cv::Mat map_mat(3, 3, CV_64F, map);
   cv::Point2f src[4] =
      {
      { 0, 0 },
      { 0, float(texture.height - 1) },
      { float(texture.width - 1), float(texture.height - 1) },
      { float(texture.width - 1), 0 },
      };

   cv::Point2f dst[4] =
   {
      { float(points2d[0][0]), float(points2d[0][1]) },
      { float(points2d[1][0]), float(points2d[1][1]) },
      { float(points2d[2][0]), float(points2d[2][1]) },
      { float(points2d[3][0]), float(points2d[3][1]) },
   };

   cvGetPerspectiveTransform(src, dst, &map_mat);
   cv::Mat img = image.clone();
   cv::Mat img2 = image.clone();
   cv::Mat_<uchar> mask(cv::Size(image.width, image.height));
   cv::Mat_<uchar> mask2(cv::Size(image.width, image.height));
   img = 0;
   img2 = 0;
   mask = 0;
   mask2 = 0;
   for (int j = 0 ; j < texture.height ; j++)
      {
      //ttesis: why must we copy the texture first?
      for (int i = 0 ; i < texture.width ; i++)
         {
         cv::Scalar s = cvGet2D(texture, j, i);
         cvSet2D(img, j, i, s);
         if ((i > 0) && (j > 0) && (i < (texture.width - 1)) && (j < (texture.height - 1)))
            mask(j, i) = 1; //ttesis: why are edges not included?
         }
      }
   cvWarpPerspective(img, img2, &map_mat);
   cvWarpPerspective(mask, mask2, &map_mat, 0);

   cvCopy(img2, image, mask2);

   return;

   }
#endif

} // namespace alvar
