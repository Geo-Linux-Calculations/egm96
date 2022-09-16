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

/*
 * This program is designed for the calculation of a geoid undulation at a point
 * whose latitude and longitude is specified.
 *
 * This program is designed to be used with the constants of EGM96 and those of
 * the WGS84(g873) system. The undulation will refer to the WGS84 ellipsoid.
 *
 * It's designed to use the potential coefficient model EGM96 and a set of
 * spherical harmonic coefficients of a correction term.
 * The correction term is composed of several different components, the primary
 * one being the conversion of a height anomaly to a geoid undulation.
 * The principles of this procedure were initially described in the paper:
 * - use of potential coefficient models for geoid undulation determination using
 *   a spherical harmonic representation of the height anomaly/geoid undulation
 *   difference by R.H. Rapp, Journal of Geodesy, 1996.
 *
 * This program is a modification of the program described in the following report:
 * - a fortran program for the computation of gravimetric quantities from high
 * degree spherical harmonic expansions, Richard H. Rapp, report 334, Department
 * of Geodetic Science and Surveying, the Ohio State University, Columbus, 1982
 */

#include <stdio.h>
#include <math.h>
#include "egm96.h"

//! EGM96 correction and harmonic coefficients
static double egm96_data[_coeffs+1][4]={0};

/* ************************************************************************** */

double egm96_hundu(double p[_coeffs+1],
             double sinml[_361+1], double cosml[_361+1],
             double gr, double re)
{
    // WGS 84 gravitational constant in m³/s² (mass of Earth’s atmosphere included)
    const double GM = EARTH_GM;
    // WGS 84 datum surface equatorial radius
    const double ae = EARTH_A;

    double ar = ae/re, arn = ar, ac = 0.0, a = 0.0;
    double sum, sumc, temp, tempc, ms, mc, target;
    unsigned k = 3, n, m;
    for (n = 2; n <= _nmax; n++)
    {
        arn *= ar;
        k++;
        sum = p[k] * egm96_data[k][2];
        sumc = p[k] * egm96_data[k][0];

        for (m = 1; m <= n; m++)
        {
            k++;
            mc = cosml[m];
            ms = sinml[m];
            tempc = egm96_data[k][0] * mc + egm96_data[k][1] * ms;
            temp  = egm96_data[k][2] * mc + egm96_data[k][3] * ms;
            sumc += (p[k] * tempc);
            sum  += (p[k] * temp);
        }
        ac += sumc;
        a += (sum * arn);
    }
    ac += egm96_data[1][0] + (p[2]*egm96_data[2][0]) + (p[3] * (egm96_data[3][0]*cosml[1] + egm96_data[3][1]*sinml[1]));

    // Add haco = ac/100 to convert height anomaly on the ellipsoid to the undulation
    // Add -0.53m to make undulation refer to the WGS84 ellipsoid
    target = ((a * GM) / (gr * re)) + (ac * WGS84_ACM) - WGS84_DLT;

    return target;
}

void egm96_dscml(double rlon, double sinml[_361+1], double cosml[_361+1])
{
    double a = sin(rlon);
    double b = cos(rlon);
    double b2 = b + b;
    unsigned m;

    sinml[1] = a;
    cosml[1] = b;
    sinml[2] = b2 * a;
    cosml[2] = b2 * b - 1;

    for (m = 3; m <= _nmax; m++)
    {
        sinml[m] = b2 * sinml[m-1] - sinml[m-2];
        cosml[m] = b2 * cosml[m-1] - cosml[m-2];
    }
}

/*!
 * \param m: order.
 * \param theta: Colatitude (radians).
 * \param rleg: Normalized legendre function.
 *
 * This subroutine computes all normalized legendre function in 'rleg'.
 * The dimensions of array 'rleg' must be at least equal to nmax+1.
 * All calculations are in double precision.
 *
 * Original programmer: Oscar L. Colombo, Dept. of Geodetic Science the Ohio State University, August 1980.
 * ineiev: I removed the derivatives, for they are never computed here.
 */
void egm96_legfdn(unsigned m, double theta, double rleg[_361+1])
{
    static double drts[1301], dirt[1301], cothet, sithet, rlnn[_361+1];
    static int ir; // TODO 'ir' must be set to zero before the first call to this sub.

    unsigned nmax1 = _nmax + 1;
    unsigned nmax2p = (2 * _nmax) + 1;
    unsigned m1 = m + 1;
    unsigned m2 = m + 2;
    unsigned m3 = m + 3;
    unsigned n, n1, n2;

    if (!ir)
    {
        ir = 1;
        for (n = 1; n <= nmax2p; n++)
        {
            drts[n] = sqrt(n);
            dirt[n] = 1.0 / drts[n];
        }
    }

    cothet = cos(theta);
    sithet = sin(theta);

    // compute the legendre functions
    rlnn[1] = 1;
    rlnn[2] = sithet * drts[3];
    for (n1 = 3; n1 <= m1; n1++)
    {
        n = n1 - 1;
        n2 = 2 * n;
        rlnn[n1] = drts[n2 + 1] * dirt[n2] * sithet * rlnn[n];
    }

    switch (m)
    {
        case 1:
            rleg[2] = rlnn[2];
            rleg[3] = drts[5] * cothet * rleg[2];
            break;
        case 0:
            rleg[1] = 1;
            rleg[2] = cothet * drts[3];
            break;
    }
    rleg[m1] = rlnn[m1];

    if (m2 <= nmax1)
    {
        rleg[m2] = drts[m1*2 + 1] * cothet * rleg[m1];
        if (m3 <= nmax1)
        {
            for (n1 = m3; n1 <= nmax1; n1++)
            {
                n = n1 - 1;
                if ((!m && n < 2) || (m == 1 && n < 3)) continue;
                n2 = 2 * n;
                rleg[n1] = drts[n2+1] * dirt[n+m] * dirt[n-m] * (drts[n2-1] * cothet * rleg[n1-1] - drts[n+m-1] * drts[n-m-1] * dirt[n2-3] * rleg[n1-2]);
            }
        }
    }
}

/*!
 * \param lat: Latitude in radians.
 * \param lon: Longitude in radians.
 * \param re: Geocentric radius.
 * \param rlat: Geocentric latitude.
 * \param gr: Normal gravity (m/sec²).
 *
 * This subroutine computes geocentric distance to the point, the geocentric
 * latitude, and an approximate value of normal gravity at the point based the
 * constants of the WGS84(g873) system are used.
 */
void egm96_radgra(double lat, double lon, double *rlat, double *gr, double *re)
{
    const double a = EARTH_A;
    const double e2 = EARTH_E2;
    const double geqt = EARTH_GEQT;
    const double k = EARTH_K;
    double lts = sin(lat), ltc = cos(lat);
    double lns = sin(lon), lnc = cos(lon);
    double t1 = lts * lts;
    double n = a / sqrt(1.0 - (e2 * t1));
    double t2 = n * ltc;
    double x = t2 * lnc;
    double y = t2 * lns;
    double z = (n * (1.0 - e2)) * lts;
    double x2 = x * x, y2 = y * y, z2 = z * z;

    *re = sqrt(x2 + y2 + z2);                              // compute the geocentric radius
    *rlat = atan(z / sqrt(x2 + y2));                       // compute the geocentric latitude
    *gr = geqt * (1.0 + (k * t1)) / sqrt(1.0 - (e2 * t1)); // compute normal gravity (m/sec²)
}

/*!
 * \brief Compute the geoid undulation from the EGM96 potential coefficient model, for a given latitude and longitude.
 * \param lat: Latitude in radians.
 * \param lon: Longitude in radians.
 * \return The geoid undulation / altitude offset (in meters).
 */
double egm96_undulation(double lat, double lon)
{
    double p[_coeffs+1], sinml[_361+1], cosml[_361+1], rleg[_361+1];
    double rlat, gr, re, target;
    unsigned int j, m, i, n, nmax1 = _nmax + 1;

    // compute the geocentric latitude, geocentric radius, normal gravity
    egm96_radgra(lat, lon, &rlat, &gr, &re);
    rlat = (M_PI / 2) - rlat;

    for (j = 1; j <= nmax1; j++)
    {
        m = j - 1;
        egm96_legfdn(m, rlat, rleg);
        for (i = j ; i <= nmax1; i++)
        {
            n = (((i - 1) * i) / 2) + m + 1;
            p[n] = rleg[i];
        }
     }
     egm96_dscml(lon, sinml, cosml);
     target = egm96_hundu(p, sinml, cosml, gr, re);

     return target;
}

/* ************************************************************************** */

double egm96_compute_altitude_offset(double lat, double lon)
{
    const double rad = (180.0 / M_PI);
    double target = egm96_undulation(lat/rad, lon/rad);
    return target;
}

/* ************************************************************************** */

/*!
 * \param f_12: EGM96 coefficients file.
 *
 * The even degree zonal coefficients given below were computed for the WGS84(g873)
 * system of constants and are identical to those values used in the NIMA gridding procedure.
 * Computed using subroutine grs written by N.K. PAVLIS.
 */
void egm96_dhcsin(FILE *f_12)
{
    double c, s, ec, es;

    const double j2 = WGS84_J2;
    const double j4 = WGS84_J4;
    const double j6 = WGS84_J6;
    const double j8 = WGS84_J8;
    const double j10 = WGS84_J10;

    unsigned n, m = (((_nmax+1) * (_nmax+2)) / 2);
    for (n = 1; n <= m; n++)
    {
        egm96_data[n][2] = egm96_data[n][3] = 0;
    }

    while (fscanf(f_12, "%i %i %lf %lf %lf %lf", &n, &m, &c, &s, &ec, &es) == 6)
    {
        if (n > _nmax) continue;
        n = ((n * (n + 1)) / 2) + m + 1;
        egm96_data[n][2] = c;
        egm96_data[n][3] = s;
    }
    egm96_data[4][2] += j2 / sqrt(5);
    egm96_data[11][2] += j4 / 3.0;
    egm96_data[22][2] += j6 / sqrt(13);
    egm96_data[37][2] += j8 / sqrt(17);
    egm96_data[56][2] += j10 / sqrt(21);
}

/* ************************************************************************** */

void egm96_init_arrays(char* egmname, char* corname)
{
    FILE *f_1, *f_12;
    int ig, n, m;
    double t1, t2;
    unsigned int i;

    f_1 = fopen(corname, "rb");   // correction coefficient file: modified with 'sed -e"s/D/e/g"' to be read with fscanf
    f_12 = fopen(egmname, "rb");    // potential coefficient file

    if (f_1 && f_12)
    {
        for (i = 1; i <= _coeffs; i++)
            egm96_data[i][0] = egm96_data[i][1] = 0;

        while (fscanf(f_1, "%i %i %lg %lg", &n, &m, &t1, &t2) == 4)
        {
            ig = ((n * (n+1)) / 2) + m + 1;
            egm96_data[ig][0] = t1;
            egm96_data[ig][1] = t2;
        }

        // the correction coefficients are now read in
        // the potential coefficients are now read in,
        // and the reference even degree zonal harmonic coefficients removed to degree 6
        egm96_dhcsin(f_12);

        fclose(f_1);
        fclose(f_12);
    }
}

/*!
 * \brief Write precomputed EGM96 correction and harmonic coefficients to egm96_data.h
 */
void egm96_write_arrays(char* genname)
{
    FILE *precomp_out;
    unsigned int i;

    precomp_out = fopen(genname, "wb");
    if (precomp_out)
    {
        fprintf(precomp_out, "#ifndef EGM96_DATA_H\n");
        fprintf(precomp_out, "#define EGM96_DATA_H\n\n");
        fprintf(precomp_out, "//! Precomputed EGM96 correction and harmonic coefficients\n");
        fprintf(precomp_out, "static const double egm96_data[65342][4] = {\n");

        for (i = 0; i <= _coeffs; i++)
        {
            fprintf(precomp_out, "{%g,%g,%g,%g},\n", egm96_data[i][0], egm96_data[i][1], egm96_data[i][2], egm96_data[i][3]);
        }

        fprintf(precomp_out, "};\n\n");
        fprintf(precomp_out, "#endif // EGM96_DATA_H\n");

        fclose(precomp_out);
    }
}

/* ************************************************************************** */
