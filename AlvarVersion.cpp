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
 * Implementation of AlvarVersion.  Moved static initialization here from
 * the old Alvar.h.  See AlvarVersion.h for additional rationale.
 *
 * TTC (12/31/15)
 *
 * Probably need to change this to reflect my fork of the original code if this is released.
*/

/*****************************************************************************
 ******************************  I N C L U D E  ******************************
 ****************************************************************************/

#include "Alvar.h"
#include "AlvarVersion.h"

#include <iostream>

namespace alvar {

/**
 * \brief Tag version string.
 *
 * The tag contains alpha, beta and release candidate versions.
 */
const char* AlvarVersion::ALVAR_VERSION_TAG = "3.0.0D";

/**
 * \brief Revision version string.
 *
 * The revision contains an identifier from the source control system.
 */
const char* AlvarVersion::ALVAR_VERSION_REVISION = "";

/**
 * \brief Entire version string.
 */
const char* AlvarVersion::ALVAR_VERSION = "3.0.0D";

/**
 * \brief Entire version string without dots.
 */
const char* AlvarVersion::ALVAR_VERSION_NODOTS = "300D";

/**
 * \brief Date the library was built.
 */
const char* AlvarVersion::ALVAR_DATE = "2015-15-12";

/**
 * \brief System the library was built on.
 */
const char* AlvarVersion::ALVAR_SYSTEM = "Windows 6.2 AMD64";

void AlvarVersion::alvarInfo()
{
    std::cerr << "ALVAR Druai Consulting Fork" << ALVAR_VERSION << " - A Library for Virtual and Augmented Reality"
         << std::endl;
    std::cerr << "Copyright 2007-2012 VTT Technical Research Centre of Finland" << std::endl;
    std::cerr << "Licensed under the GNU Lesser General Public License" << std::endl;
    std::cerr << "Built on " << ALVAR_DATE << " for " << ALVAR_SYSTEM << std::endl;
    std::cerr << std::endl;
}

#if 0
struct AlvarLoader
   {
   AlvarLoader()
      {
		alvarInfo();
      return;
      }

   } alvarBasicLoader;

#endif

} // namespace alvar
