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

/* ************************************************************************** */

/*!
 * \brief Main function.
 * \return 0 if success.
 *
 * The input files consist of:
 * - correction coefficient set ("CORRCOEF") => unit = 1
 * - potential coefficient set ("EGM96") => unit = 12
 * - points at which to compute ("INPUT.dat") => unit = 14
 * The output file is:
 * - computed geoid heights ("OUTPUT.dat") => unit = 20
 * - precomputed egm96_data.h (to use with the library)
 */
int main(int argc, char *argv[])
{
    FILE *f_14, *f_20;
    double flat, flon, u;
    if (argc < 3)
    {
        printf("EGM96\n");
        printf("Usage:\n");
        printf(" %s EGM96-file COR-file [INPUT-file] [OUTPUT-file]\n", argv[0]);
    }
    else
    {
        egm96_init_arrays(argv[1], argv[2]);

        if (argc > 3)
            f_14 = fopen(argv[3], "rb");
        else
            f_14 = stdin;
        if (argc > 4)
            f_20 = fopen(argv[4], "wb");
        else
            f_20 = stdout;
        if (f_14 && f_20)
        {
            // read geodetic latitude,longitude at point undulation is wanted
            while (2 == fscanf(f_14, "%lg %lg", &flat, &flon))
            {
                // compute the geocentric latitude, geocentric radius, normal gravity
                u = egm96_compute_altitude_offset(flat, flon);

                // u is the geoid undulation from the EGM96 potential coefficient model
                // including the height anomaly to geoid undulation correction term
                // and a correction term to have the undulations refer to the
                // WGS84 ellipsoid. the geoid undulation unit is meters.

                fprintf(f_20, "%14.7f %14.7f %16.7f\n", flat, flon, u);
            }

            fclose(f_14);
            fclose(f_20);
        }
    }

    return 0;
}

/* ************************************************************************** */
