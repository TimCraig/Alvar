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

#ifndef POSE_H
#define POSE_H

/*****************************************************************************
 ******************************  I N C L U D E  ******************************
 ****************************************************************************/

#include "Alvar.h"
#include "Rotation.h"

/**
 * \file Pose.h
 *
 * \brief This file implements a pose.
 */

namespace alvar {

/*****************************************************************************
 *
 ***  class Pose
 *
 * Pose representation
 *
 * The rotation part of the transformation is handled by \e Rotation .
 * The translation part is stored internally using homogeneous 4-vector.
 *
 * Internally in ALVAR we assume coordinate system where
 * 'x' is right, 'y' is down, and 'z' is forward. However
 * the \e SetMatrixGL and \e GetMatrixGL change the pose
 * to support coordinates commonly used in OpenGL:
 * 'x' is right, 'y' is up, and 'z' is backward.
 *
 *****************************************************************************/

class ALVAR_EXPORT Pose
   {
   // Note, although we are using homogeneous coordinates x, y, z, w --  w is now mostly ignored
   public:      
      /** \e Constructor */
      Pose();

      ~Pose() = default;

      /** \e Constructor using the given translation and rotation elements
    *  \param tra Column vector containing three translation elements
    *  \param rot Handled using the \e Rotation class
    *  \param t   Handled using the \e Rotation class
    */
      Pose(cv::Mat& tra, cv::Mat_<double>& rot, Rotation::RotationType t);
      /** \e Constructor with 3x3, 3x4 or 4x4 matrix representation
    *  \param mat A 3x3 rotation matrix or 3x4 / 4x4 transformation matrix
    */
      Pose(cv::Mat_<double>& mat);

      /** \e Copy constructor */
      Pose(const Pose& src);

      /** Assignment operator for copying \e Pose class */
      Pose& operator  = (const Pose& rhs);

      /** \e Reset the pose */
      void Reset();

      /** Set the transformation from the given matrix \e mat
    *  \param mat A 3x3 rotation matrix or 3x4 / 4x4 transformation matrix
    */
      void SetMatrix(const cv::Mat& mat);
      /** \brief Set the \e Pose using OpenGL's transposed format.
    *  Note, that by default this also mirrors both the y- and z-axis (see \e Camera and \e Pose for more information)
    *  \param gl OpenGL 4x4 transformation matrix elements in column-order
    */
      void SetMatrixGL(double gl[16], bool mirror = true);

      void SetRodrigues(cv::Mat_<double>& rod)
         {
         rotation.SetRodrigues(rod);
         return;
         }

      cv::Mat_<double>  GetRodrigues() const
         {
         return (rotation.GetRodrigues());
         }

      /** Get the transformation into the given matrix \e mat
    *  \param mat A 3x3 rotation matrix
    */
      cv::Mat_<double> GetRotationMatrix() const
         {
         return (rotation.GetMatrix());
         }

      /** Get the transformation into the given matrix \e mat
    *  \param mat A 3x3 rotation matrix
    */
      cv::Mat_<double> GetHomogeneousMatrix() const;

      /** \brief Get the transformation matrix representation of the \e Pose using OpenGL's transposed format.
    *  Note, that by default this also mirrors both the y- and z-axis (see \e Camera and \e Pose for more information)
    *  \param gl OpenGL 4x4 transformation matrix elements in column-order
    */
      void GetMatrixGL(double gl[16], bool mirror=true);
      /** \e Transpose the transformation */
      void Transpose();

      /** Invert the pose */
      void Invert();

      /** \e Mirror the \e Pose
    *  \param x If true mirror the x-coordinates
    *  \param y If true mirror the y-coordinates
    *  \param z If true mirror the z-coordinates
    */
      void Mirror(bool x, bool y, bool z);
      /** Set the translation part for the \e Pose
    *  \param tra Column vector containing three translation elements
    */
      void SetTranslation(const cv::Mat& tra);
      /** Set the translation part for the \e Pose
    *  \param tra Array containing three translation elements
    */
      void SetTranslation(const double* tra);
      /** Set the translation part for the \e Pose */
      void SetTranslation(const double x, const double y, const double z);
      /** Get the translation part from the \e Pose
    *  \param tra Column vector where the three translation elements are filled in
    */
      void GetTranslation(cv::Mat& tra) const;

      /** \e Output for debugging purposes */
      void Output() const;

   protected:
      Rotation rotation;
      cv::Mat_<double>  translation_mat;

   }; // end of class Pose

} // namespace alvar

#endif
