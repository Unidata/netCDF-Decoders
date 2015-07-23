/*
 *   Copyright 1995, University Corporation for Atmospheric Research.
 *   See ../COPYRIGHT file for copying and redistribution conditions.
 */
/* $Id: get_prod.h,v 1.3 2002/12/05 19:46:31 rkambic Exp $ */

/* GRIB product delimiters */

#ifndef _GET_PROD_H_
#define _GET_PROD_H_

#include <stdio.h>

typedef struct prod {		/* a raw product with ID */
    unsigned int len;
    unsigned char *bytes;
    char *id;			/* WMO header, if any, or manufactured ID */
} prod;

enum PROD_MARK {READ_ERR, FOUND_START, FOUND_END, NOT_FOUND, FOUND_EOF };

#ifdef __cplusplus
extern "C" int get_prod (FILE *stream, int timeout, prod* prodp);
extern "C" int add_get_prod (FILE *stream, int timeout, prod* prodp);
#elif defined(__STDC__)
extern int get_prod (FILE *stream, int timeout, prod* prodp);
extern int add_get_prod (FILE *stream, int timeout, prod* prodp);
#else
extern int get_prod ( /* FILE *stream, int timeout, prod* prodp */ );
extern int add_get_prod ( /* FILE *stream, int timeout, prod* prodp */ );
#endif

#endif /* !_GET_PROD_H_ */
