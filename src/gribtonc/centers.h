/*
 *   Copyright 1995, University Corporation for Atmospheric Research.
 *   See ../COPYRIGHT file for copying and redistribution conditions.
 */
/* $Id: centers.h,v 1.2 1995/05/18 22:03:30 russ Exp $ */

/* National/International Originating Centers (Assigned by the WMO) */

#ifndef CENTERS_H_
#define CENTERS_H_
#define CENTER_NMC	7	/* US Weather Service - National Met. Center */
#define CENTER_NWSTG	8	/* US Weather Service - NWS Telecomms Gateway */
#define CENTER_NWSFS	9	/* US Weather Service - Field Stations */

#define CENTER_JMA	34	/* Japanese Meteorological Agency - Tokyo */
#define CENTER_NHC	52	/* US National Hurricane Center, Miami */
#define CENTER_CMS	54	/* Canadian Meteorological Service - Montreal */
#define CENTER_USAF	57	/* US Air Force - Global Weather Center */
#define CENTER_FNOC	58	/* US Navy - Fleet Numerical Oceanography Center */
#define CENTER_FSL	59	/* NOAA Forecast Systems Lab, Boulder CO */
#define CENTER_NCAR	60	/* National Center for Atmospheric Research (NCAR), Boulder CO */
#define CENTER_UKMET	74	/* U.K. Met Office - Bracknell */
#define CENTER_FR	85	/* French Weather Service - Toulouse */
#define CENTER_ESA	97	/* European Space Agency (ESA) */
#define CENTER_ECMWF	98	/* European Center for Medium-Range Weather
				   Forecasts - Reading */
#define CENTER_NL	99	/* DeBilt, Netherlands */

/* National Sub-Centers (Assigned by the Nation) */
#define SUBCENTER_NCEP_REANA      1	/* NCEP Re-Analysis Project */
#define SUBCENTER_NCEP_EP	  2	/* NCEP Ensemble Products */
#define SUBCENTER_NCEP_CO	  3	/* NCEP Central Operations */
#define SUBCENTER_NCEP_EMC	  4	/* NCEP Environmental Modelling Center */
#define SUBCENTER_NCEP_HPC	  5	/* NCEP Hydrometeorological Prediction Center */
#define SUBCENTER_NCEP_MPC	  6	/* NCEP Marine Prediction Center */
#define SUBCENTER_NCEP_CPC	  7	/* NCEP Climate Prediction Center */
#define SUBCENTER_NCEP_APC	  8	/* NCEP Aviation Weather Center */
#define SUBCENTER_NCEP_SPC	  9	/* NCEP Storm Prediction Center */
#define SUBCENTER_NCEP_TPC	 10	/* NCEP Tropical Prediction Center */

#define SUBCENTER_NWSFS_ABRFC	150	/* ABRFC - Arkansas-Red River RFC, Tulsa OK */
#define SUBCENTER_NWSFS_AKFC	151	/* Alaska RFC, Anchorage, AK */
#define SUBCENTER_NWSFS_CBRFC	152	/* CBRFC - Colorado Basin RFC, Salt Lake City, UT */
#define SUBCENTER_NWSFS_CNRFC	153	/* CNRFC - California-Nevada RFC, Sacramento, CA */
#define SUBCENTER_NWSFS_LMRFC	154	/* LMRFC - Lower Mississippi RFC, Slidel, LA */
#define SUBCENTER_NWSFS_MARFC	155	/* MARFC - Mid Atlantic RFC, State College, PA */
#define SUBCENTER_NWSFS_MBRFC	156	/* MBRFC - Missouri Basin RFC, Kansas City, MO */
#define SUBCENTER_NWSFS_NCRFC	157	/* NCRFC - North Central RFC, Minneapolis, MN */
#define SUBCENTER_NWSFS_NERFC	158	/* NERFC - Northeast RFC, Hartford, CT */
#define SUBCENTER_NWSFS_NWRFC	159	/* NWRFC - Northwest RFC, Portland, OR */
#define SUBCENTER_NWSFS_OHRFC	160	/* OHRFC - Ohio Basin RFC, Cincinnati, OH */
#define SUBCENTER_NWSFS_SERFC	161	/* SERFC - Southeast RFC, Atlanta, GA */
#define SUBCENTER_NWSFS_WGRFC	162	/* WGRFC - West Gulf RFC, Fort Worth, TX */
#define SUBCENTER_NWSFS_OUN	170	/* OUN - Norman OK WFO */

#ifdef __cplusplus
extern "C" char* centername(int center);
extern "C" char* subcentername(int center, int sub);
#elif defined(__STDC__)
extern char* centername(int center);
extern char* subcentername(int center, int sub);
#else
extern char* centername(/* int center */);
extern char* subcentername(/* int center, int sub */);
#endif

#endif /* CENTERS_H_ */
