/*
 *   Copyright 1996, University Corporation for Atmospheric Research.
 *   See ../COPYRIGHT file for copying and redistribution conditions.
 */
/* $Id: recs.c,v 1.7 2004/06/01 22:20:50 rkambic Exp $ */

#include <netcdf.h>
#include <udunits.h>
#include "ulog.h"
#include "timeunits.h"
#include "nc.h"
#include "recs.h"
#include "nuwg.h"
#include "emalloc.h"

/*
 * Initializes (reftime,valtime) pair table from open netcdf file.
 * Returns -1 on failure.
 */
int
new_recs(nc)
    ncfile *nc;
{
    int ncid = nc->ncid;
    int recid = nc->recid;
    long nrecs;
    long size;
    char *ncname = nc->ncname;

    /* get number of records */
    if (ncdiminq(ncid, recid, (char *)0, &nrecs) == -1) {
	uerror("%s: can't get number of records", ncname);
	return -1;
    }
    nc->rt = (rectimes *) emalloc(sizeof(rectimes));
    nc->rt->nrecs = nrecs;

    if (nrecs == 0) {
#define RECS_INIT_SIZE 64	/* initial size of growable record table */
	size = RECS_INIT_SIZE;
    } else {
	size = 2*nrecs;
    }
    nc->rt->size = size;
    nc->rt->reftimes = (double *) emalloc(size * sizeof(double));
    nc->rt->valtimes = (double *) emalloc(size * sizeof(double));
    
    if (nrecs > 0) {
	int nerrs = 0;
	int irec;
	long start[1];
	long count[1];
	start[0] = 0;
	count[0] = nrecs;
	if (ncvarget(ncid, nc->reftimeid, start, count,
		     (void *)nc->rt->reftimes) == -1) {
	    uerror("%s: can't get reftimes", nc->ncname);
	    nerrs++;
	}
	if (ncvarget(ncid, nc->valtimeid, start, count,
		     (void *)nc->rt->valtimes ) == -1) {
	    uerror("%s: can't get valtimes", nc->ncname);
	    nerrs++;
	}
	if(nerrs)
	    return -1;
    }
    return 0;
}


/*
 * Frees memory used by (reftime,valtime) pair table from open netcdf file.
 */
int
free_recs(rp)
    rectimes *rp;
{
    if(rp) {
	if(rp->reftimes)
	    free(rp->reftimes);
	if(rp->valtimes)
	    free(rp->valtimes);
	free(rp);
    }
}


/*
 * Returns record number determined by (reftime,valtime) pair.  If
 * (reftime,valtime) is in table, returns corresponding record.  If
 * (reftime,valtime) is not in table, increments record count, enters
 * (reftime,valtime) for new record in table, updates reftime, valtime in
 * netCDF file, updates human-comprehensible time variables in netCDF file
 * (if any), and returns record.  Returns -1 on failure.
 */
long
getrec(nc, reftime, valtime, htp)
    ncfile *nc;
    double reftime;
    double valtime;
    humtime *htp;
{
    extern int nogrow;  /* grow/no grow unlimited dimension beyond CDL define */

    humtime hftp;
    int year, month, day, hour, minute;
    float second;
    double *dp;
    double *reftimes = nc->rt->reftimes;
    double *valtimes = nc->rt->valtimes;
    char datetime[100];
    char forecasttime[100];
    int ncid = nc->ncid;
    int i;

    if (nc->rt->nrecs == 0) {
	int nerrs = 0;
	int dimvaloffsetsid;
	long dimvaloffsets;
	double *valoffsets;
        utUnit valunits;
	long start[1];
	long count[1];
	long ix[2];
	long ccount[2];
	char datetime[100];
	char forecasttime[100];
	char* utunitsstr;

	/* get dimension and value of valtime_offset */
	dimvaloffsetsid = nuwg_getdim(ncid, DIM_VALOFFSETS );
	if (ncdiminq(ncid, dimvaloffsetsid, (char *)0, &dimvaloffsets) == -1) {
	    uerror("%s: can't get dimvaloffsets", nc->ncname);
	    nerrs++;
	}
/*
	uerror("dimvaloffsetsid = %d", dimvaloffsetsid);
	uerror("dimvaloffsets = %d", dimvaloffsets);
	    uerror("reftimes[i] %d=%f", i, reftimes[i] );
	    uerror("valtimes[i] %d=%f", i, valtimes[i] );
	    uerror("valoffsets %d=%f", i, valoffsets[i] );
	    utPrint( nc->vars[nc->reftimeid]->bunitp, &utunitsstr );
*/
        valoffsets = (double *) emalloc(dimvaloffsets * sizeof(double)); 
	start[0] = 0;
	count[0] = dimvaloffsets; 
	if (ncvarget(ncid, nc->valoffsetid, start, count,
		     (void *)valoffsets) == -1) {
	    uerror("%s: can't get valoffsets", nc->ncname);
	    nerrs++;
	}
	if(nerrs)
	    return -1;

	utCalendar(reftime, nc->vars[nc->reftimeid]->bunitp,
	    &year, &month, &day, &hour, &minute, &second);
	sprintf(datetime, "%.4d-%.2d-%.2d %.2d:%.2d:%.2dZ",
		    year, month, day, hour, minute, second);

	/* update humtime record variables, if any, from htp */
	if (nc->datetimeid > 0) {
	    sprintf(datetime, "%.4d-%.2d-%.2d %.2d:%.2d:%.2dZ",
		    htp->year, htp->month, htp->day,
		    htp->hour, htp->minute, htp->second);
	    ccount[1] = strlen(datetime)+1; /* include terminating null */
	}
	/* populate reftimes, valtimes, and update datetime in ncfile */
	ix[1] = 0;
	ccount[0] = 1;
    	for (i=0; i < dimvaloffsets; i++) { /* initialize reftimes & valtimes */
	    reftimes[i] = reftime;
	    valtimes[i] = reftime + valoffsets[i];
	    ix[0] = i;
	    ncvarput(ncid, nc->datetimeid, ix, ccount, datetime);
	    /* calculate and enter forecast time */
/*
            utCopy( nc->vars[nc->valtimeid]->bunitp, valunitsp );
	    get_units(nc, nc->valtimeid, valunitsp );
	    utCalendar(valtimes[i], nc->vars[nc->valtimeid]->bunitp,
            utCopy( nc->vars[nc->valtimeid]->bunitp, &valunits );
	    uerror("recs forecast time = %s", forecasttime );
	    uerror("\nrecs check reftime %f",  reftime );
	    utCalendar(reftime, nc->vars[nc->reftimeid]->bunitp,
	    &year, &month, &day, &hour, 
	    &minute, &second);
	    sprintf(datetime, "%.4d-%.2d-%.2d %.2d:%.2d:%.2dZ",
		    year, month, day,
		    hour, minute, second);
	    uerror("recs reftime time check = %s", datetime );
*/
	    utCalendar(valtimes[i], nc->vars[nc->valtimeid]->bunitp,
	    	&year, &month, &day, &hour, &minute, &second);

	    sprintf(forecasttime, "%.4d-%.2d-%.2d %.2d:%.2d:%.2dZ",
		    year, month, day, hour, minute, second);
	    ncvarput(ncid, nc->forecasttimeid, ix, ccount, forecasttime);
	}

    	/* Update reftimes and valtimes in ncfile */
	if (ncvarput(ncid, nc->reftimeid, start, count, reftimes) == -1 ||
	    ncvarput(ncid, nc->valtimeid, start, count, valtimes) == -1) {
	    uerror("%s: failed to add new reftime, valtime", nc->ncname);
	    return -1;
	}
	nc->rt->nrecs = dimvaloffsets;
    }

    /* First look in table of existing records */
    for (i=0; i < nc->rt->nrecs; i++) {
	if (reftime == reftimes[i] && valtime == valtimes[i])
	    return i;
    }

    if ( nogrow ) {
	uerror("%s: No grow mode, skipping new valoffset = %f", nc->ncname, htp->valoffset);
	return -1;
    }

    /* Didn't find it, so create a new record */
    if (nc->rt->nrecs+1 == nc->rt->size) {/* no room for another record,
					     double size of table before
					     adding */
	nc->rt->size *= 2;
	nc->rt->reftimes = (double *) erealloc(nc->rt->reftimes,
		     nc->rt->size * sizeof(double));
	nc->rt->valtimes = (double *) erealloc(nc->rt->valtimes,
		     nc->rt->size * sizeof(double));
	reftimes = nc->rt->reftimes;
	valtimes = nc->rt->valtimes;
    }
    reftimes[nc->rt->nrecs] = reftime;
    valtimes[nc->rt->nrecs] = valtime;
				
    {				/* Update reftimes and valtimes in ncfile */
	long ix[2];
	long count[2];

	ix[0] = nc->rt->nrecs;
	count[0] = 1;

	if (ncvarput1(ncid, nc->reftimeid, ix, &reftime) == -1 ||
	    ncvarput1(ncid, nc->valtimeid, ix, &valtime) == -1) {
	    uerror("%s: failed to add new reftime, valtime", nc->ncname);
	    return -1;
	}
	/* update humtime record variables, if any, from htp */
	ix[1] = 0;
	if (nc->datetimeid > 0) {
	    sprintf(datetime, "%.4d-%.2d-%.2d %.2d:%.2d:%.2dZ",
		    htp->year, htp->month, htp->day,
		    htp->hour, htp->minute, htp->second);
	    count[1] = strlen(datetime)+1; /* include terminating null */
	    ncvarput(ncid, nc->datetimeid, ix, count, datetime);
	}
	if (nc->forecasttimeid > 0) {
	    utCalendar(valtime, nc->vars[nc->valtimeid]->bunitp,
	    	&year, &month, &day, &hour, &minute, &second);

	    sprintf(forecasttime, "%.4d-%.2d-%.2d %.2d:%.2d:%.2dZ",
		    year, month, day, hour, minute, second);
	    ncvarput(ncid, nc->forecasttimeid, ix, count, forecasttime);
	}
	uerror("%s: adding new valoffset = %f", nc->ncname, htp->valoffset);
/*
	if (nc->valoffsetid > 0) {
	    ncvarput1(ncid, nc->valoffsetid, ix, &htp->valoffset);
	}
*/
    }
    nc->rt->nrecs++;
    return nc->rt->nrecs-1;
}
