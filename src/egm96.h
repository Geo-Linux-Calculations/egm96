/*
 * Copyright (c) 2006 D.Ineiev <ineiev@yahoo.co.uk>
 * Copyright (c) 2020 Emeric Grange <emeric.grange@gmail.com>
 *
 * This software is provided 'as-is', without any express or implied warranty.
 * In no event will the authors be held liable for any damages arising from
 * the use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 *
 * 1. The origin of this software must not be misrepresented; you must not
 *    claim that you wrote the original software. If you use this software
 *    in a product, an acknowledgment in the product documentation would be
 *    appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 *    misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 */

#ifndef EGM96_H
#define EGM96_H
/* ************************************************************************** */

#define _coeffs    (65341) //!< Size of correction and harmonic coefficients arrays (361*181)
#define _nmax      (360)   //!< Maximum degree and orders of harmonic coefficients.
#define _361       (361)
#define EARTH_A    (6378137.0)
#define EARTH_GM   (0.3986004418e15)
#define EARTH_E2   (0.00669437999013)
#define EARTH_GEQT (9.7803253359)
#define EARTH_K    (0.00193185265246)
#define WGS84_J2   (0.108262982131e-2)
#define WGS84_J4   (-.237091120053e-05)
#define WGS84_J6   (0.608346498882e-8)
#define WGS84_J8   (-0.142681087920e-10)
#define WGS84_J10  (0.121439275882e-13)
#define WGS84_ACM  (0.01)
#define WGS84_DLT  (0.53)

void egm96_init_arrays(char* egmname, char* corname);
void egm96_write_arrays(char* genname);
/*!
 * \brief Compute the geoid undulation from the EGM96 potential coefficient model, for a given latitude and longitude.
 * \param latitude: Latitude (in degrees).
 * \param longitude: Longitude (in degrees).
 * \return The geoid undulation / altitude offset (in meters).
 */
double egm96_compute_altitude_offset(double latitude, double longitude);

/* ************************************************************************** */
#endif // EGM96_H
