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

/*
 * This file contains version information removed from the original Alvar.h
 * and consolidated into a class.
 *
 * Rational: I simply got tired of seeing all the unused variable warnings
 * and why (potentially) have copies of all these strings everywhere.  Just
 * neater.  TTC (12/31/15)
*/

#ifndef ALVARVERSION_H
#define ALVARVERSION_H


/**
 * \brief Main ALVAR namespace.
 */
namespace alvar {

/*****************************************************************************
 *
 ***  class AlvarVersion
 *
 *****************************************************************************/

class AlvarVersion
   {
   public:
      AlvarVersion() = default;

      AlvarVersion(const AlvarVersion& src) = delete;

      ~AlvarVersion() = default;

      AlvarVersion& operator=(const AlvarVersion& rhs) = delete;

      static void alvarInfo();

      /**
       * \brief Major version number.
       */
      const int ALVAR_VERSION_MAJOR = 3;

      /**
       * \brief Minor version number.
       */
      const int ALVAR_VERSION_MINOR = 0;

      /**
       * \brief Patch version number.
       */
      const int ALVAR_VERSION_PATCH = 0;

      /**
       * \brief Tag version string.
       *
       * The tag contains alpha, beta and release candidate versions.
       */
      static const char* ALVAR_VERSION_TAG;

      /**
       * \brief Revision version string.
       *
       * The revision contains an identifier from the source control system.
       */
      static const char* ALVAR_VERSION_REVISION;

      /**
       * \brief Entire version string.
       */
      static const char* ALVAR_VERSION;

      /**
       * \brief Entire version string without dots.
       */
      static const char* ALVAR_VERSION_NODOTS;

      /**
       * \brief Date the library was built.
       */
      static const char* ALVAR_DATE;

      /**
       * \brief System the library was built on.
       */
      static const char* ALVAR_SYSTEM;

   protected:

   private:

   }; // end of class AlvarVersion

}

#endif
