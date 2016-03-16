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

#ifndef ROTATION_H
#define ROTATION_H

#include <iostream>

/**
 * \file Rotation.h
 *
 * \brief This file implements a parametrized rotation.
 */

#include "Alvar.h"
#include "Util.h"

namespace alvar {

// Use Boost::quaternion?  All this seems pretty archaic.  Revisit when working for cleanup.  TTC!!

/*****************************************************************************
 *
 ***  class Rotation
 *
 * Rotation structure and transformations between different parameterizations.
 *
 * Rotation is represented as a quaternion.  Rather than use quaterions to
 * perform the rotation math, this class is used to convert between
 * quaternions and rotation matrices in various forms.
 *
 *****************************************************************************/

class ALVAR_EXPORT Rotation
   {
   public:
      using Quaternion = std::array<double, 4>;

      /**
   * \brief Rotation can be represented in four ways: quaternion (QUAT), matrix (MAT), euler angles (EUL) and
   * exponential map (ROD).
   * Only used in constructor Rotation(cv::Mat& data, RotationType t)
   */
      enum class RotationType { QUAT, MAT, EUL, ROD };

      Rotation();

      /**
   * \brief Constructor.
   * \param data	Rotation data stored in cv::Mat. With RotationType::MAT both 3x3 and 4x4 matrices are allowed.
   * \param t		Rotation type that corresponds to data.
   */
      Rotation(cv::Mat_<double>& mat, RotationType t);

      Rotation(const Rotation& src) : quaternion(src.quaternion)
         {
         return;
         }

      ~Rotation() = default;

      Rotation& operator = (const Rotation& rhs)
         {
         if (this != &rhs)
            {
            quaternion = rhs.quaternion;
            } // end if
         return (*this);
         }

      Rotation& operator += (const Rotation& v);
      //Rotation& operator +  () {return *this;}

      double& operator[](size_t i)
         {
         return (quaternion[i]);
         }

      double operator[](size_t i) const
         {
         return (quaternion[i]);
         }

      void Transpose();

      /**
   * \brief Simple function to mirror a rotation matrix in different directions.
   * \param mat	Matrix to be mirrored.
   * \param x
   * \param y
   * \param z
   */
      static void MirrorMat(cv::Mat& mat, bool x, bool y, bool z);

      /**
   * \brief Mirrors the rotation in selected directions.
   * \param x
   * \param y
   * \param z
   */
      void Mirror(bool x, bool y, bool z);

      /**
   * \brief Resets the rotation into identity.
   */
      void Reset();

      /**
   * \brief Converts 3x3 rotation matrix into Rodrigues representation.
   * \param mat	3x3 rotation matrix.
   * \param rod	Resulting 3x1 rotation vector.
   */
      static cv::Mat_<double> Mat9ToRod(const cv::Mat_<double>& mat);

      /**
   * \brief Converts 3x1 rotation vector into 3x3 rotation matrix using Rodrigues' formula.
   * \param rod 3x1 rotation vector.
   * \param Resulting 3x3 rotation matrix.
   */
      static cv::Mat_<double> RodToMat9(const cv::Mat_<double>& rod);

      /**
   * \brief Inverts unit quaternion.
   * \param q	Unit quaternion to be inverted.
   * \param qi	Resulting quaternion.
   */
      void Conjugate();
      void Conjugate(Rotation& conj) const
         {
         conj = *this;
         conj.Conjugate();

         return;
         }

      /**
   * \brief Normalizes a quaternion.
   * \param q	Quaternion to be normalized.
   */
      void Normalize();

      /**
   * \brief Quaternion multiplication.
   * \param q1
   * \param q2
   * \param result
   */
      static void Mul(const Rotation& q1, const Rotation& q2, Rotation& result);

      //%  The quaternion has to be normalized!!!
      /**
   * \brief Converts a rotation described by a quaternion into 3x3 rotation matrix.
   * \param quat	Rotation in quaternion form.
   * \param mat	Corresponding 3x3 rotation matrix.
   */
      cv::Mat_<double> ToMat9() const;

      /**
   * \brief Converts a rotation described by a quaternion into 4x4 OpenGL-like transformation matrix.
   * \param quat	Rotation in quaternion form.
   * \param mat	Resulting 4x4 transformation matrix.
   */
      cv::Mat_<double>  ToMat16() const;

      /**
   * \brief Converts a rotation described by a quaternion into Euler angles.
   * \param q		Rotation in quaternion form.
   * \param eul	Resulting Euler angles.
   */
      cv::Mat_<double> ToEul() const;

      /**
   * \brief Converts a 3x3 rotation martix into quaternion form.
   * \param mat	3x3 rotation matrix.
   * \param quat	Resulting quaternion.
   */
      void FromMat9(const cv::Mat_<double>& mat9);

      /**
   * \brief Converts a rotation described by Euler angles into quaternion form.
   * \param eul	Rotation in Euler angles.
   * \param quat	Resulting quaternion.
   */
      void FromEul(const cv::Mat_<double>& eul);

      /**
   * \brief Sets the rotation from given quaternion.
   * \param mat	Input quaternion (4x1 cv::Mat).
   */
      void SetQuaternion(const cv::Mat_<double>& mat);
      void SetQuaternion(const Rotation& rot);

      /**
   * \brief Sets the rotation from given Euler angles.
   * \param mat	Input Euler angles (3x1 cv::Mat).
   */
      void SetEuler(const cv::Mat_<double>& eul);

      /**
   * \brief Sets the rotation from given rotation vector.
   * \param mat	Input rotation vector (3x1 cv::Mat).
   */
      void SetRodrigues(cv::Mat_<double>& rod)
         {
         assert((rod.rows == 3) && (rod.cols == 1));
         // Not sure why we have to go to a Mat9 here.  Just need to compute the components of the quaternion.
         FromMat9(RodToMat9(rod));
         return;
         }

      /**
   * \brief Sets the rotation from given rotation matrix. 3x3 and 4x4 matrices are allowed.
   * \param mat	Input rotation matrix (3x3 or 4x4 cv::Mat).
   */
      void SetMatrix(const cv::Mat&  mat);

      /**
   * \brief Returns the rotation in 3x3 matrix form.
   * \param mat	The rotation is stored here.
   */
      cv::Mat_<double> GetMatrix() const
         {
         return (ToMat9());
         }

      /**
   * \brief Returns the rotation in 4x4 homogeneous matrix form (Translation is zero).
   * \param mat	The rotation is stored here.
   */
      cv::Mat_<double> GetHomogeneousMatrix() const
         {
         return (ToMat16());
         }

      /**
   * \brief Returns the rotation in rotation vector form.
   * \param mat	The rotation is stored here (3x1 cv::Mat).
   */
      cv::Mat_<double> GetRodrigues() const
         {
         return (Rotation::Mat9ToRod(ToMat9()));
         }

      /**
   * \brief Returns the rotation in Euler angles.
   * \param mat	The rotation is stored here (3x1 cv::Mat).
   */
      cv::Mat_<double> GetEuler() const
         {
         return (ToEul());
         }

      /**
   * \brief Returns the rotation in quaternion form.
   * \param mat	The rotation is stored here (4x1 cv::Mat).
   */
      void GetQuaternion(Quaternion& quat) const
         {
         quat = quaternion;
         return;
         }

      // Dump the quaternion values to cout
      void Output() const;

   protected:
      // Storage for the 4 elements
      Quaternion quaternion;

   }; // end of class Rotation


} // namespace alvar

#endif
