/*
 *   Copyright 1995, University Corporation for Atmospheric Research.
 *   See ../COPYRIGHT file for copying and redistribution conditions.
 */
/* $Id: levels.h,v 1.14 2004/12/21 20:53:52 rkambic Exp $ */

/* GRIB levels */

#ifndef LEVELS_H_
#define LEVELS_H_

#define LEVEL_SURFACE	1	/* surface of the Earth */
#define LEVEL_CLOUD_BASE	2 /* cloud base level */
#define LEVEL_CLOUD_TOP	3	/* cloud top level */
#define LEVEL_ISOTHERM	4	/* 0 degree isotherm level */
#define LEVEL_ADIABAT	5	/* adiabatic condensation level */
#define LEVEL_MAX_WIND	6	/* maximium wind speed level */
#define LEVEL_TROP	7	/* at the tropopause */
#define LEVEL_TOP	8	/* nominal top of atmosphere */
#define LEVEL_SEABOT	9	/* sea bottom */
#define LEVEL_TMPL	20	/* temperature in 1/100 K (2 octets) */
#define LEVEL_ISOBARIC	100	/* isobaric level */
#define LEVEL_LISO	101	/* layer between two isobaric levels */
#define LEVEL_MEAN_SEA	102	/* mean sea level */
#define LEVEL_FH	103	/* fixed height level */
#define LEVEL_LFHM	104	/* layer between 2 height levels above MSL */
#define LEVEL_FHG	105	/* fixed height above ground */
#define LEVEL_LFHG	106	/* layer between 2 height levels above ground */
#define LEVEL_SIGMA	107	/* sigma level */
#define LEVEL_LS	108	/* layer between 2 sigma levels */
#define LEVEL_HY	109	/* Hybrid level */
#define LEVEL_LHY	110	/* Layer between 2 hybrid levels */
#define LEVEL_BLS	111	/* Depth below land surface */
#define LEVEL_LBLS	112	/* Layer between 2 depths below land surface */
#define LEVEL_ISEN	113	/* Isentropic (theta) level */
#define LEVEL_LISEN	114	/* Layer between 2 isentropic (theta) levels */
#define LEVEL_PDG	115	/* level at specified pressure difference from ground */
#define LEVEL_LPDG	116	/* layer between levels at specif. pressure diffs from ground */
#define LEVEL_PV	117	/* potential vorticity */
#define LEVEL_ETAL	119	/* ETA value in 1/10000 */
#define LEVEL_LETA	120	/* layer between two ETA levels */
#define LEVEL_LISH	121	/* layer between 2 isobaric surfaces (high precision) */
#define LEVEL_FHGH	125	/* height level above ground (high precision) */
#define LEVEL_LSH	128	/* layer between 2 sigma levels (high precision) */
#define LEVEL_LISM	141	/* layer between 2 isobaric surfaces (mixed precision) */
#define LEVEL_DBS	160	/* depth below sea level */
#define LEVEL_ATM	200	/* entire atmosphere considered as a single layer */
#define LEVEL_OCEAN	201	/* entire ocean considered as a single layer */

#define LEVEL_HTFL	204	/* Highest tropospheric freezing level */
#define LEVEL_GCBL	206	/* Grid scale cloud bottom level */
#define LEVEL_GCTL	207	/* Grid scale cloud top level */
#define LEVEL_BCBL	209	/* Boundary layer cloud bottom level */
#define LEVEL_BCTL	210	/* Boundary layer cloud top level */
#define LEVEL_BCY 	211	/* Boundary layer cloud layer */
#define LEVEL_LCBL	212	/* Low cloud bottom level */
#define LEVEL_LCTL	213	/* Low cloud top level */
#define LEVEL_LCY	214	/* Low cloud layer */
#define LEVEL_CEIL	215	/* Cloud ceiling */
#define LEVEL_PBLRI	220	/* Planetary Boundary Layer */
#define LEVEL_MCBL	222	/* Middle cloud bottom level */
#define LEVEL_MCTL	223	/* Middle cloud top level */
#define LEVEL_MCY	224	/* Middle cloud layer */
#define LEVEL_HCBL	232	/* High cloud bottom level */
#define LEVEL_HCTL	233	/* High cloud top level */
#define LEVEL_HCY	234	/* Highcloud layer */
#define LEVEL_OITL	235	/* Ocean Isotherm Level */
#define LEVEL_OLYR	236	/* Layer between two depths below ocean surface */
#define LEVEL_OBML	237	/* Bottom of Ocean Mixed Layer */
#define LEVEL_OBIL	238	/* Bottom of Ocean Isothermal Layer */
#define LEVEL_S26CY	239	/* Layer Ocean Surface and 26C Ocean Isothermal Level */
#define LEVEL_OMXL	240	/* Ocean Mixed Layer */
#define LEVEL_CCBL	242	/* convective cloud bottom level */
#define LEVEL_CCTL	243	/* convective cloud top level */
#define LEVEL_CCY 	244	/* convective cloud layer */
#define LEVEL_LLTW 	245	/* lowest level of the wet bulb zero */
#define LEVEL_MTHE	246	/* maximum equivalent potential temperature level */
#define LEVEL_EHLT	247	/* equilibrium level */
#define LEVEL_SCBL	248	/* Shallow convective cloud bottom level */
#define LEVEL_SCTL	249	/* Shallow convective cloud top level */
#define LEVEL_DCBL	251	/* Deep convective cloud bottom level */
#define LEVEL_DCTL	252	/* Deep convective cloud top level */
#define LEVEL_LBLSW	253	/* Lowest bottom level of supercooled liquid water layer */
#define LEVEL_HTLSW	254	/* Highest top level of supercooled liquid water layer */
#define LEVEL_FL	9999	/* FAA flight level */


#ifdef __cplusplus
extern "C" double mblev(int* levels);
extern "C" char* levelname(int level);
extern "C" char* levelsuffix(int level);
extern "C" char* levelunits(int level);
Ixtern "C" long level_index(double level, float* levels, long nlevels);
extern "C" long layer_index(double top, double bot, float* tops,
			    float* bots, long nlayers);
extern "C" int level1(int flag, int* levels);
extern "C" int level2(int flag, int* levels);
#elif defined(__STDC__)
extern double mblev(int* levels);
extern char* levelname(int level);
extern char* levelsuffix(int level);
extern char* levelunits(int level);
extern long level_index(double level, float* levels, long nlevels);
extern long layer_index(double top, double bot, float* tops,
			    float* bots, long nlayers);
extern int level1(int flag, int* levels);
extern int level2(int flag, int* levels);
#else
extern double mblev ( /* int* levels */ );
#define const
extern char* levelname ( /* int level */ );
extern char* levelsuffix ( /* int level */ );
extern char* levelunits ( /* int level */ );
extern long level_index ( /* double level, float* levels, long nlevels */ );
extern long layer_index( /* double top, double bot, float* tops,
			    float* bots, long nlayers */ );
extern int level1( /* int flag, int* levels */ );
extern int level2( /* int flag, int* levels */ );
#endif

#endif /* LEVELS_H_ */
