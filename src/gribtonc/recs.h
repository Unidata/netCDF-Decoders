/*
 *   Copyright 1996, University Corporation for Atmospheric Research.
 *   See ../COPYRIGHT file for copying and redistribution conditions.
 */
/* $Id: recs.h,v 1.2 2003/03/25 19:49:53 rkambic Exp $ */

/* netCDF-specific structures and functions for record dimension */

#ifndef RECS_H_
#define RECS_H_

typedef struct rectimes {
    long nrecs;			/* current number of records */
    long size;			/* current size of table, >= nrecs */
    double *reftimes;		/* nth record has reftime[n] */
    double *valtimes;		/* nth record has valtime[n] */
} rectimes;

#ifdef __cplusplus
extern "C" int new_recs(ncfile*);
extern "C" int free_recs(rectimes*);
extern "C" long getrec(ncfile* nc, double reftime, double valtime, humtime *htp);
#elif defined(__STDC__)
extern int new_recs(ncfile* nc);
extern int free_recs(rectimes*);
extern long getrec(ncfile* nc, double reftime, double valtime,
		   humtime *htp);
#else
extern int new_recs(/* ncfile* nc*/);
extern int free_recs(/* rectimes* */);
extern long getrec(/* ncfile* nc, double reftime, double valtime, humtime *htp */);
#endif

#endif /* RECS_H_ */
