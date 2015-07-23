/*
 *   Copyright 1996 University Corporation for Atmospheric Research
 *   See ../COPYRIGHT file for copying and redistribution conditions.
 */
/* $Id: timeunits.c,v 1.6 1995/12/07 15:31:00 russ Exp $ */

#include "ulog.h"
#include "timeunits.h"
#include "product_data.h"

/*
 * Name used in printed dumps
 */
char *
tunitsname(ii)
    int ii;
{
    switch (ii) {
      case TUNITS_MIN : return "Minute";
      case TUNITS_HOUR : return "Hour";
      case TUNITS_DAY : return "Day";
      case TUNITS_MONTH : return "Month";
      case TUNITS_YEAR : return "Year";
      case TUNITS_DECADE : return "Decade";
      case TUNITS_NORMAL : return "Normal (30 Years)";
      case TUNITS_CENTURY : return "Century";
      case TUNITS_3HR : return "3 hours";
      case TUNITS_6HR : return "6 hours";
      case TUNITS_12HR : return "12 hours";
      case TUNITS_SECOND : return "Second";
    }
    /* default */
    return "Unknown" ;
}


/*
 * Name used in units conversions with udunits library
 */
char *
tunits(ii)
    int ii;
{
    switch (ii) {
      case TUNITS_MIN : return "minute";
      case TUNITS_HOUR : return "hour";
      case TUNITS_DAY : return "day";
      case TUNITS_MONTH : return "year/12"; /* ??? */
      case TUNITS_YEAR : return "year";
      case TUNITS_DECADE : return "10 year";
      case TUNITS_NORMAL : return "30 year";
      case TUNITS_CENTURY : return "100 year";
      case TUNITS_3HR : return "3 hour";
      case TUNITS_6HR : return "6 hour";
      case TUNITS_12HR : return "12 hour";
      case TUNITS_SECOND : return "second";
    }
    /* default */
    return "Unknown" ;
}


char *
triname(ii)			/* Return time range indicator */
    int ii;
{
    switch (ii) {
      case TRI_P1: return	"Reference Time + P1" ;
      case TRI_IAP: return	"Initialized Analysis Product (P1=0)" ;
      case TRI_P12: return	"Valid from P1 to P2" ;
      case TRI_Ave: return	"Average from P1 to P2" ;
      case TRI_Acc: return	"Accumulation from P1 to P2" ;
      case TRI_Dif: return	"Difference from P2 to P1" ;
      case TRI_LP1: return	"Reference Time + Long P1" ;
      case TRI_AvgN: return	"Special average Algorithm 0" ;
      case TRI_AccN: return	"Special accumulation Algorithm 0" ;
      case TRI_AvgN1: return	"Special average Algorithm 1" ;
      case TRI_AccN1: return	"Special accumulation Algorithm 1" ;
      case TRI_AvgN2: return	"Special average Algorithm 2" ;
      case TRI_VarN: return 	"Temporal (co)variance";
      case TRI_SdN: return 	"Standard deviation";
      case TRI_AvgN3: return	"Special average Algorithm 3" ;
      case TRI_AccN3: return	"Special accumulation Algorithm 3" ;
    }
    /* default */
    return "Unknown" ;
}


char *
trisuffix(ii)			/* Time range indicator suffix */
    int ii;
{
    switch (ii) {
      case TRI_P1: return	"" ;
      case TRI_IAP: return	"init_times" ;
      case TRI_P12: return	"valid_times" ;
      case TRI_Ave: return	"average_times" ;
      case TRI_Acc: return	"accum_times" ;
      case TRI_Dif: return	"diff_times" ;
      case TRI_LP1: return	"" ;
      case TRI_AvgN: return	"average0" ;
      case TRI_AccN: return	"accum0" ;
      case TRI_AvgN1: return	"average1" ;
      case TRI_AccN1: return	"accum1" ;
      case TRI_AvgN2: return	"average2" ;
      case TRI_VarN: return 	"var";
      case TRI_SdN: return 	"stdev";
      case TRI_AvgN3: return	"average3" ;
      case TRI_AccN3: return	"accum3" ;
    }
    /* default */
    return "Unknown" ;
}


int
trinum(ii)			/* Number of time range indicators */
    int ii;
{
    switch (ii) {
      case TRI_P1:
      case TRI_IAP:
      case TRI_LP1:
	  return 1;

      case TRI_P12:
      case TRI_Ave:
      case TRI_Acc:
      case TRI_Dif:
      case TRI_AvgN:
      case TRI_AccN:
      case TRI_AvgN1:
      case TRI_AccN1:
      case TRI_AvgN2:
      case TRI_VarN:
      case TRI_SdN:
      case TRI_AvgN3:
      case TRI_AccN3:
	  return 2;
    }
    /* default */
    return 0 ;
}


/*
 * Return product valid time increment from product reference time.  This is
 * not well-defined for some accumulations and averages, so for these we
 * just return 0.
 */
int
frcst_time(pp)
    product_data *pp;
{
    switch(pp->tr_flg) {
    case TRI_P1:
    case TRI_IAP:
	return pp->tr[0];
    case TRI_LP1:
	return 256 * pp->tr[0] + pp->tr[1];

    case TRI_P12:
    case TRI_Ave:
	return (pp->tr[0] + pp->tr[1]) * 0.5;
    case TRI_Acc:
    case TRI_Dif:
	return pp->tr[1];
    case TRI_AvgN:
    case TRI_AccN:
    case TRI_AvgN1:
    case TRI_AccN1:
    case TRI_AvgN2:
    case TRI_VarN:
    case TRI_SdN:	
    case TRI_AvgN3:
    case TRI_AccN3:
	return 0;
    }
				/*     default: */
    uerror("frcst_time: unknown time range flag %d", pp->tr_flg);
    return 0;
}
