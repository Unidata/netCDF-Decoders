/*
 *   Copyright 1996 University Corporation for Atmospheric Research
 *   See ../COPYRIGHT file for copying and redistribution conditions.
 */
/* $Id: centers.c,v 1.2 1996/01/05 17:45:38 russ Exp $ */

#include "centers.h"

char *
centername(int center)
{
    switch(center){
    case CENTER_NMC:	
	return "US Weather Service - National Met. Center";
    case CENTER_NWSTG:
	return "US Weather Service - NWS Telecomms Gateway";
    case CENTER_NWSFS:
	return "US Weather Service - Field Stations";
    case CENTER_JMA:
	return "Japanese Meteorological Agency - Tokyo";
    case CENTER_NHC:
	return "US National Hurricane Center, Miami";
    case CENTER_CMS:
	return "Canadian Meteorological Service - Montreal";
    case CENTER_USAF:
	return "US Air Force - Global Weather Center";
    case CENTER_FNOC:
	return "US Navy  - Fleet Numerical Oceanography Center";
    case CENTER_FSL:
	return "NOAA Forecast Systems Lab, Boulder CO";
    case CENTER_NCAR:
	return "National Center for Atmospheric Research (NCAR), Boulder CO";
    case CENTER_UKMET:
	return "U.K. Met Office - Bracknell";
    case CENTER_FR:
	return "French Weather Service - Toulouse";
    case CENTER_ESA:
	return "European Space Agency (ESA)";
    case CENTER_ECMWF:
	return "European Center for Medium-Range Weather Forecasts - Reading";
    case CENTER_NL:
	return "DeBilt, Netherlands";
    }
    /* default */
    return "unknown";
}


char *
subcentername(int center, int sub)
{
    switch(center){
    case CENTER_NMC:	
	switch(sub) {
	case SUBCENTER_NCEP_REANA:
	    return "NCEP Re-Analysis Project";
	case SUBCENTER_NCEP_EP:
	    return "NCEP Ensemble Products";
	case SUBCENTER_NCEP_CO:
	    return "NCEP Central Operations";
	case SUBCENTER_NCEP_EMC:
	    return "NCEP Environmental Modelling Center";
	case SUBCENTER_NCEP_HPC:
	    return "NCEP Hydrometeorological Prediction Center";
	case SUBCENTER_NCEP_MPC:
	    return "NCEP Marine Prediction Center";
	case SUBCENTER_NCEP_CPC:
	    return "NCEP Climate Prediction Center";
	case SUBCENTER_NCEP_APC:
	    return "NCEP Aviation Weather Center";
	case SUBCENTER_NCEP_SPC:
	    return "NCEP Storm Prediction Center";
	case SUBCENTER_NCEP_TPC:
	    return "NCEP Tropical Prediction Center";
	default:
	    break;
	}
	break;
    case CENTER_NWSFS:
	switch(sub) {
	case SUBCENTER_NWSFS_ABRFC:
	    return "ABRFC - Arkansas-Red River RFC, Tulsa OK";
	case SUBCENTER_NWSFS_AKFC:
	    return "Alaska RFC, Anchorage, AK";
	case SUBCENTER_NWSFS_CBRFC:
	    return "CBRFC - Colorado Basin RFC, Salt Lake City, UT";
	case SUBCENTER_NWSFS_CNRFC:
	    return "CNRFC - California-Nevada RFC, Sacramento, CA";
	case SUBCENTER_NWSFS_LMRFC:
	    return "LMRFC - Lower Mississippi RFC, Slidel, LA";
	case SUBCENTER_NWSFS_MARFC:
	    return "MARFC - Mid Atlantic RFC, State College, PA";
	case SUBCENTER_NWSFS_MBRFC:
	    return "MBRFC - Missouri Basin RFC, Kansas City, MO";
	case SUBCENTER_NWSFS_NCRFC:
	    return "NCRFC - North Central RFC, Minneapolis, MN";
	case SUBCENTER_NWSFS_NERFC:
	    return "NERFC - Northeast RFC, Hartford, CT";
	case SUBCENTER_NWSFS_NWRFC:
	    return "NWRFC - Northwest RFC, Portland, OR";
	case SUBCENTER_NWSFS_OHRFC:
	    return "OHRFC - Ohio Basin RFC, Cincinnati, OH";
	case SUBCENTER_NWSFS_SERFC:
	    return "SERFC - Southeast RFC, Atlanta, GA";
	case SUBCENTER_NWSFS_WGRFC:
	    return "WGRFC - West Gulf RFC, Fort Worth, TX";
	case SUBCENTER_NWSFS_OUN:
	    return "OUN - Norman OK WFO";
	default:
	    break;
	}
	break;
    default:
	break;
    }
    return "unknown";
}
