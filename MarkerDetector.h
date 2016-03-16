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
 *
 * Originally, this class was divided into two parts, a regular class for low level
 * implementation and a template class templated on the particular marker class.  It
 * "appeared" they thought about detecting different marker types polymorphically but
 * that's not what was done.  To simplify things, I eliminated the regular base class
 * and eliminated dealing with dynamically allocating marker instances and some ugly
 * casting.  Some protected methods that simply wrapped operations on the marker vectors
 * were eliminated in favor of calling the vector methods directly which seems much
 * cleaner.  TTC (3/14/16 Happy Pi Day!)
 *
 */

#ifndef MARKER_DETECTOR_H
#define MARKER_DETECTOR_H

/**
 * \file MarkerDetector.h
 *
 * \brief This file implements a generic marker detector.
 */

/*****************************************************************************
 ******************************  I N C L U D E  ******************************
 ****************************************************************************/

#include "Alvar.h"
#include "Util.h"
#include "ConnectedComponents.h"
#include "Draw.h"
#include "Camera.h"
#include "Marker.h"
#include "Rotation.h"
#include "Line.h"
#include <algorithm>
using std::rotate;
#include <list>
#include <vector>
#include <map>
#include <cassert>
#include <type_traits>

namespace alvar {


/*****************************************************************************
 *
 ***  class MarkerDetector
 *
 * MarkerDetector for detecting markers of type M
 * M Class that extends Marker
 *
*****************************************************************************/

// MARKER must be derived from Marker.
template<class MARKER>
class MarkerDetector
   {
   static_assert(std::is_base_of<Marker, MARKER>::value, "MarkerDetector requires class with base Marker");

   public:
      /** Constructor */
      MarkerDetector()
         {
         SetMarkerSize();
         SetOptions();

         labeling = nullptr;
         markers = new std::vector<MARKER>;
         track_markers = new std::vector<MARKER>;

         return;

         }

      MarkerDetector(const MarkerDetector& src) = delete;

      /** Destructor */
      ~MarkerDetector()
         {
         delete labeling;
         delete markers;
         delete track_markers;

         return;
         }

      MarkerDetector& operator=(const MarkerDetector& rhs) = delete;

      /** Clear the markers that are tracked */
      void TrackMarkersReset()
         {
         track_markers->clear();
         return;
         }

      void TrackMarkerAdd(int id, PointDouble corners[4]);

      /** Set the default marker size to be used for all markers unless
      * \param _edge_length Length of the marker's edge in whatever units you are using (e.g. cm)
      * \param _res The marker content resolution in pixels. By default we use 5x5 markers. If you use 0 with
      * \e MarkerData, the marker resolution is detected automatically.
      * \param _margin The marker margin resolution in pixels (The actual captured marker image has pixel resolution
      * of _margin+_res+_margin)
      *
      * \note The default marker content resolution (_res) of 5 can only detect marker ids from 0 to 255. For larger
      * marker ids, you need to increase the marker content resolution accordingly.
      */
      void SetMarkerSize(double _edge_length = 1, int _res = 5, double _margin = 2)
         {
         edge_length = _edge_length;
         res = _res;
         margin = _margin;
         map_edge_length.clear(); // TODO: Should we clear these here?

         return;

         }

      /** Set marker size for specified marker id. This needs to be called after setting the default marker size.
   * \param id The specified marker id
   * \param _edge_length Length of the marker's edge in whatever units you are using (e.g. cm)
   */
      void SetMarkerSizeForId(unsigned long id, double _edge_length = 1)
         {
         map_edge_length[id] = _edge_length;
         return;
         }

      /** Set marker size for specified marker id. This needs to be called after setting the default marker size.
   * \param _detect_pose_grayscale Do we detect marker pose using grayscale optimization?
   */
      void SetOptions(bool _detect_pose_grayscale = false)
         {
         detect_pose_grayscale = _detect_pose_grayscale;
         return;
         }

      int Detect(cv::Mat& image,
                 Camera& cam,
                 bool track = false,
                 bool visualize = false,
                 double max_new_marker_error = 0.08,
                 double max_track_error = 0.2,
                 Labeling::LabelingMethod labeling_method = Labeling::LabelingMethod::CVSEQ,
                 bool update_pose = true);

      int DetectAdditional(cv::Mat& image, Camera& cam, bool visualize = false, double max_track_error = 0.2);

      const std::vector<MARKER>* getMarkers() const
         {
         return (markers);
         }

      const std::vector<MARKER>* getTrackMarkers() const
         {
         return (track_markers);
         }

   protected:
      std::vector<MARKER>* markers;
      std::vector<MARKER>* track_markers;
      Labeling* labeling;

      std::map<unsigned long, double> map_edge_length;
      double edge_length;
      int res;
      double margin;
      bool detect_pose_grayscale;

      void swap_marker_tables()
         {
         std::swap(markers, track_markers);

         return;

         }

   private:

   }; // end of class MarkerDetector


/******************************************************************************
*
***  MarkerDetector::TrackMarkerAdd
*
* Add markers to be tracked
* Sometimes application or e.g. the \e MultiMarker implementation knows
* more about marker locations. Then this method can be used after \e Detect
* to indicate where additional trackable markers could be found. The
* \e DetectAdditional is called for tracking these.
*
******************************************************************************/

template<typename MARKER>
void MarkerDetector<MARKER>::TrackMarkerAdd(int id, PointDouble corners[4])
   {
   MARKER mn(edge_length, res, margin);
   if (map_edge_length.find(id) != map_edge_length.end())
      {
      mn.SetMarkerSize(map_edge_length[id], res, margin);
      }

   mn.SetId(id);
   mn.marker_corners_img.clear();
   mn.marker_corners_img.push_back(corners[0]);
   mn.marker_corners_img.push_back(corners[1]);
   mn.marker_corners_img.push_back(corners[2]);
   mn.marker_corners_img.push_back(corners[3]);
   track_markers->push_back(mn);

   return;

   } // end of method MarkerDetector::TrackMarkerAdd

/******************************************************************************
*
***  MarkerDetector::Detect
*
* Detect Markers from image
*
* The coordinates are little tricky. Here is a short summary.
*
* - Image (top-left origin).
* - The marker corners in the image are searched in sub-pixel accuracy in
*   counter-clockwise order starting from "lower-left" corner.
* - The corresponding marker corners and marker points are in marker coordinates
*   (x is to east, y is to north, and z is up from the marker)
* - The marker points are read from inside the margins starting from top-left
*   and reading the bits first left-to-right one line at a time.
*
******************************************************************************/

template<typename MARKER>
int MarkerDetector<MARKER>::Detect(cv::Mat& image, Camera& cam, bool track,
      bool visualize, double max_new_marker_error, double max_track_error,
      Labeling::LabelingMethod labeling_method, bool update_pose)
   {
   double error = -1;

   // Swap marker tables
   swap_marker_tables();
   markers->clear();

   switch (labeling_method)
      {
      case Labeling::LabelingMethod::CVSEQ :
         if (labeling == nullptr)
            labeling = new LabelingCvSeq;

         ((LabelingCvSeq*) labeling)->SetOptions(detect_pose_grayscale);
         break;
      }

   labeling->LabelSquares(image, &cam, visualize);
   std::vector<std::vector<PointDouble>>& blob_corners = labeling->blob_corners;
   cv::Mat& gray = labeling->gray;

   int orientation;

   // When tracking we find the best matching blob and test if it is near enough?
   if (track)
      {
      for (auto& mn : *track_markers)
         {
         if (mn.GetError(Marker::DECODE_ERROR | Marker::MARGIN_ERROR) == 0)
            {
            // We track only perfectly decoded markers
            int track_i = -1;
            int track_orientation = 0;
            double track_error = 1.0e200;
            for (size_t i = 0 ; i < blob_corners.size() ; ++i)
               {
               if (!blob_corners[i].empty())
                  {
                  mn.CompareCorners(blob_corners[i], &orientation, &error);
                  if (error < track_error)
                     {
                     track_i = i;
                     track_orientation = orientation;
                     track_error = error;
                     }
                  } // end if
               }

            if (track_error <= max_track_error)
               {
               mn.SetError(Marker::DECODE_ERROR, 0);
               mn.SetError(Marker::MARGIN_ERROR, 0);
               mn.SetError(Marker::TRACK_ERROR, track_error);
               mn.UpdatePose(blob_corners[track_i], cam, track_orientation, update_pose);
               markers->push_back(mn);
               blob_corners[track_i].clear(); // We don't want to handle this again...
               if (visualize)
                  mn.Visualize(image, cam, cv::Scalar(0, 255, 255));
               }
            } // end if
         }
      }

   // Now we go through the rest of the blobs -- in case there are new markers...
   for (auto& blob_corner : blob_corners)
      {
      if (!blob_corner.empty())
         {
         MARKER mn(edge_length, res, margin);
         if (mn.UpdateContent(blob_corner, gray, cam) &&
             mn.DecodeContent(orientation) &&
             (mn.GetError(Marker::MARGIN_ERROR | Marker::DECODE_ERROR) <= max_new_marker_error))
            {
            if (map_edge_length.find(mn.GetId()) != map_edge_length.end())
               {
               mn.SetMarkerSize(map_edge_length[mn.GetId()], res, margin);
               }

            mn.UpdatePose(blob_corner, cam, orientation, update_pose);
            markers->push_back(mn);

            if (visualize)
               mn.Visualize(image, cam, cv::Scalar(0, 0, 255));
            }
         } // end if
      }

   return (static_cast<int>(markers->size()));

   } // end of method MarkerDetector::Detect

/******************************************************************************
*
***  MarkerDetector::DetectAdditional
*
******************************************************************************/

template<typename MARKER>
int MarkerDetector<MARKER>::DetectAdditional(cv::Mat& image, Camera& cam, bool visualize, double max_track_error)
   {
   int count = -1;

   if (labeling != nullptr)
      {
      double error = -1;
      int orientation;
      count = 0;
      std::vector<std::vector<PointDouble>>& blob_corners = labeling->blob_corners;

      for (auto& mn : *track_markers)
         {
         if (mn.GetError(Marker::DECODE_ERROR | Marker::MARGIN_ERROR) == 0)
            {
            // We track only perfectly decoded markers
            int track_i = -1;
            int track_orientation = 0;
            double track_error = 1.0e200;

            for (size_t i ; blob_corners.size() ; i++)
               {
               if (!blob_corners[i].empty())
                  {
                  mn->CompareCorners(blob_corners[i], &orientation, &error);
                  if (error < track_error)
                     {
                     track_i = i;
                     track_orientation = orientation;
                     track_error = error;
                     }
                  } // end if
               }

            if (track_error <= max_track_error)
               {
               mn.SetError(Marker::DECODE_ERROR, 0);
               mn.SetError(Marker::MARGIN_ERROR, 0);
               mn.SetError(Marker::TRACK_ERROR, track_error);
               mn.UpdatePose(blob_corners[track_i], cam, track_orientation);
               markers->push_back(mn);
               count++;
               blob_corners[track_i].clear(); // We don't want to handle this again...

               if (visualize)
                  {
                  mn.Visualize(image, cam, cv::Scalar(255, 0, 255));
                  }
               }
            } // end if
         }
      } // end if

   return (count);

   } // end of method MarkerDetector::DetectAdditional

} // namespace alvar

#endif
