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

#ifndef CONNECTEDCOMPONENTS_H
#define CONNECTEDCOMPONENTS_H

/**
 * \file ConnectedComponents.h
 *
 * \brief This file implements connected component labeling.
 */

#include "Alvar.h"
#include "Util.h"
#include "Line.h"
#include "Camera.h"

namespace alvar {


/*****************************************************************************
 *
 ***  class Labeling
 *
 * Base class for labeling connected components from binary image.
 *
 *****************************************************************************/

class ALVAR_EXPORT Labeling
   {
   public :
      /**
       * \brief Connected components labeling methods.
      */
      enum class LabelingMethod
         {
         CVSEQ
         };

      /**
       * \brief Two alternatives for thresholding the gray image. ADAPT (adaptive threshold) is only supported currently.
      */
      enum class ThresholdMethod
         {
         THRESH,
         ADAPT
         };

      /** Constructor */
      Labeling();

      Labeling(const Labeling& src) = delete;

      /** Destructor*/
      virtual ~Labeling() = default;

      Labeling& operator=(const Labeling& rhs) = delete;

      /**
    * \brief Grayscale image that is thresholded for labeling.
   */
      cv::Mat gray;
      /**
    * \brief Binary image that is then labeled.
   */
      cv::Mat bw;

      /**
    * \brief Vector of 4-length vectors where the corners of detected blobs are stored.
   */
      std::vector<std::vector<PointDouble>> blob_corners;

      /**
    * \brief Labels image and filters blobs to obtain square-shaped objects from the scene.
   */
      virtual void LabelSquares(cv::Mat& image, Camera* pCamera = nullptr, bool visualize = false) = 0;

      bool CheckBorder(const std::vector<cv::Point2i>& contour, int width, int height) const;

      void SetThreshParams(int param1, int param2)
         {
         thresh_param1 = param1;
         thresh_param2 = param2;

         return;
         }

   protected :
      int thresh_param1;
      int thresh_param2;

   private:

   }; // end of class Labeling


/*****************************************************************************
 *
 ***  class LabelingCvSeq
 *
 * Labeling class that uses OpenCV routines to find connected components.
 *
 *****************************************************************************/

class ALVAR_EXPORT LabelingCvSeq : public Labeling
   {
   public:
      LabelingCvSeq();

      LabelingCvSeq(const LabelingCvSeq& src) = delete;

      ~LabelingCvSeq() = default;

      LabelingCvSeq& operator=(const LabelingCvSeq& rhs) = delete;

      void SetOptions(bool _detect_pose_grayscale = false);

      virtual void LabelSquares(cv::Mat& image, Camera* pCamera = nullptr, bool visualize = false) override;

      // TODO: Releases memory inside, cannot return CvSeq*
      void LabelImage(cv::Mat& image, std::vector<std::vector<cv::Point2i>>& squares,
                      int min_size, bool approx = false);

   protected :
      int _n_blobs;
      int _min_edge;
      int _min_area;
      bool detect_pose_grayscale;

   private:      

   }; // end of class LabelingCvSeq

} // namespace alvar

#endif
