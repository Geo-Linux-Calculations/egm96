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
    if (argc < 3)
    {
        printf("EGM96 generation head-file\n");
        printf("Usage:\n");
        printf(" %s EGM96-file COR-file [HEAD-file]\n", argv[0]);
    }
    else
    {
        egm96_init_arrays(argv[1], argv[2]);

        if (argc > 3)
            egm96_write_arrays(argv[3]);
        else
            egm96_write_arrays("EGM96_data.h");
    }

    return 0;
}

/* ************************************************************************** */
