/*
 *   Copyright 1995, University Corporation for Atmospheric Research
 *   See ../COPYRIGHT file for copying and redistribution conditions.
 */
/* $Id: gdes.c,v 1.21 2005/06/01 17:03:34 rkambic Exp $ */

#include "test_noserc.h"

#include <stdlib.h>                        /* for malloc(), free(), ... */
#include <stdio.h>

#include "ulog.h"
#include "emalloc.h"
#include "gdes.h"
#include "gbytem.h"
#include "grib1.h"
#include "centers.h"
#include "quasi.h"

#ifdef __STDC__
static int fill_ll(grid_ll *raw, gdes *gd);
static int fill_rll(grid_rll *raw, gdes *gd);
static int fill_sll(grid_sll * raw, gdes * gd);
static int fill_srll(grid_srll * raw, gdes * gd);
static int fill_gau(grid_gau * raw, gdes * gd);
static int fill_rgau(grid_rgau * raw, gdes * gd);
static int fill_sgau(grid_sgau * raw, gdes * gd);
static int fill_srgau(grid_srgau * raw, gdes * gd);
static int fill_sph(grid_sph * raw, gdes * gd);
static int fill_rsph(grid_rsph * raw, gdes * gd);
static int fill_ssph(grid_ssph * raw, gdes * gd);
static int fill_srsph(grid_srsph * raw, gdes * gd);
static int fill_mercator(grid_mercator * raw, gdes * gd);
static int fill_polars(grid_polars * raw, gdes * gd);
static int fill_gnomon(grid_gnomon * raw, gdes * gd);
static int fill_lambert(grid_lambert * raw, gdes * gd);
static int fill_olambert(grid_olambert * raw, gdes * gd);
static int fill_albers(grid_albers * raw, gdes * gd);
static int fill_spacev(grid_spacev * raw, gdes * gd);
static void copy_gdes_ll(gdes_ll* gd, gdes_ll* g_ll);
static void copy_gdes_gau(gdes_gau* gd, gdes_gau* g_gau);
static void copy_gdes_sph(gdes_sph* gd, gdes_sph* g_sph);
static void copy_gdes_mercat( gdes_mercator* gd,gdes_mercator* g_mercat );
static void copy_gdes_polars( gdes_polars* gd,gdes_polars* g_polars );
static void copy_gdes_lambert( gdes_lambert* gd, gdes_lambert* g_lambert);
static void copy_gdes_spacev( gdes_spacev * gd, gdes_spacev* g_spacev);
static gdes* gds_to_gdes(gds * gdsp);
static void nmc_21_24(gdes *g);
static gdes* nmc_21(void);
static gdes* nmc_22(void);
static gdes* nmc_23(void);
static gdes* nmc_24(void);
static void nmc_25_26(gdes *g);
static gdes* nmc_25(void);
static gdes* nmc_26(void);
static gdes* nmc_50(void);
static void nmc_61_64(gdes *g);
static gdes* nmc_61(void);
static gdes* nmc_62(void);
static gdes* nmc_63(void);
static gdes* nmc_64(void);
static gdes* nmc_1(void);
static gdes* nmc_2(void);
static gdes* nmc_3(void);
static gdes* nmc_5(void);
static gdes* nmc_6(void);
static gdes* nmc_27(void);
static gdes* nmc_28(void);
static gdes* nmc_100(void);
static gdes* nmc_101(void);
static gdes* nmc_104(void);
static gdes* nmc_105(void);
static gdes* nmc_207(void);
static gdes* nmc_211(void);
static gdes* nmc_212(void);
static gdes* synth_gdes(int centerid, int gridid);
static int cmp_int_start( int* v1, int n1, int* v2, int n2 );
static int cmp_float_start( float* v1, int n1, float* v2, int n2 );
static int cmp_int_set_base( int* v1, int n1, int* v2, int n2 );
static int cmp_float_set_base( float* v1, int n1, float* v2, int n2 );
static int cmp_int_set( int* v1, int n1, int* v2, int n2 );
static int gengrid_cmp( gdes* gd1, gdes* gd2 );
static int rotstr_cmp(  rotated* rot1, rotated* rot2, 
                         stretched* strch1, stretched* strch2 );
static int gdes_ll_cmp( gdes_ll* ll1, gdes_ll* ll2 );
static int gdes_gau_cmp( gdes_gau* gau1, gdes_gau* gau2 );
static int gdes_sph_cmp( gdes_sph* sph1, gdes_sph* sph2 );
static int gdes_polars_cmp( gdes_polars* pol1, gdes_polars* pol2 );
static int gdes_mercator_cmp( gdes_mercator* m1, gdes_mercator* m2 );
static int gdes_lambert_cmp( gdes_lambert* lam1, gdes_lambert* lam2 );
static int gdes_spacev_cmp( gdes_spacev* spa1, gdes_spacev* spa2 );

#endif


static int
fill_ll(raw,gd)
    grid_ll *raw;
    gdes *gd;
{
    gdes_ll *cooked = &gd->grid.ll;

    cooked->ni = g2i(raw->ni);
    cooked->nj = g2i(raw->nj);
    cooked->la1 = g3si(raw->la1)*.001; /* millidegrees to degrees */
    cooked->lo1 = g3si(raw->lo1)*.001;
    cooked->la2 = g3si(raw->la2)*.001;
    cooked->lo2 = g3si(raw->lo2)*.001;
    cooked->di = g2i(raw->di)*.001;
    cooked->dj = g2i(raw->dj)*.001;
    cooked->rot=0;
    cooked->strch=0;

    if(cooked->ni == G2I_MISSING && g2i(raw->di) == G2I_MISSING) {
        cooked->ni = GDES_INT_MISSING;
        cooked->di = GDES_FLOAT_MISSING;
        gd->quasi = QUASI_ROWS;  /* A quasi-regular grid with varying rows */
    } else if (cooked->nj == G2I_MISSING && g2i(raw->dj) == G2I_MISSING) {
        cooked->nj = GDES_INT_MISSING;
        cooked->dj = GDES_FLOAT_MISSING;
        gd->quasi = QUASI_COLS;   /* A quasi-regular grid with varying cols */
    }
    gd->ncols = cooked->ni;
    gd->nrows = cooked->nj;
    gd->npts = gd->ncols*gd->nrows;
    gd->res_flags = raw->res_flags;
    gd->scan_mode = raw->scan_mode;
    return 0;
}

static int
fill_rll(raw,gd)
    grid_rll *raw;
    gdes *gd;
{
    gdes_ll *cooked = &gd->grid.ll;
    cooked->ni = g2i(raw->ni);
    cooked->nj = g2i(raw->nj);
    cooked->la1 = g3si(raw->la1)*.001; /* millidegrees to degrees */
    cooked->lo1 = g3si(raw->lo1)*.001;
    cooked->la2 = g3si(raw->la2)*.001;
    cooked->lo2 = g3si(raw->lo2)*.001;
    cooked->di = g2i(raw->di)*.001;
    cooked->dj = g2i(raw->dj)*.001;
    cooked->rot=(rotated *)emalloc(sizeof(rotated));
    cooked->rot->lat = g3si(raw->lapole)*.001;
    cooked->rot->lon = g3si(raw->lopole)*.001;
    cooked->rot->angle = g4f(raw->angrot);
    cooked->strch=0;

    if(cooked->ni == G2I_MISSING && g2i(raw->di) == G2I_MISSING) {
        cooked->ni = GDES_INT_MISSING;
        cooked->di = GDES_FLOAT_MISSING;
        gd->quasi = QUASI_ROWS;   /* A quasi-regular grid with varying rows */
    } else if (cooked->nj == G2I_MISSING && g2i(raw->dj) == G2I_MISSING) {
        cooked->nj = GDES_INT_MISSING;
        cooked->dj = GDES_FLOAT_MISSING;
        gd->quasi = QUASI_COLS;   /* A quasi-regular grid with varying cols */
    }
    gd->ncols = cooked->ni;
    gd->nrows = cooked->nj;
    gd->npts = gd->ncols*gd->nrows;
    gd->res_flags = raw->res_flags;
    gd->scan_mode = raw->scan_mode;
    return 0;
}

static int
fill_sll(raw,gd)
    grid_sll *raw;
    gdes *gd;
{
    gdes_ll *cooked = &gd->grid.ll;
    cooked->ni = g2i(raw->ni);
    cooked->nj = g2i(raw->nj);
    cooked->la1 = g3si(raw->la1)*.001; /* millidegrees to degrees */
    cooked->lo1 = g3si(raw->lo1)*.001;
    cooked->la2 = g3si(raw->la2)*.001;
    cooked->lo2 = g3si(raw->lo2)*.001;
    cooked->di = g2i(raw->di)*.001;
    cooked->dj = g2i(raw->dj)*.001;
    cooked->rot=0;
    cooked->strch=(stretched *)emalloc(sizeof(stretched));
    cooked->strch->lat = g3si(raw->lastr)*.001;
    cooked->strch->lon = g3si(raw->lostr)*.001;
    cooked->strch->factor = g4f(raw->stretch);

    if(cooked->ni == G2I_MISSING && g2i(raw->di) == G2I_MISSING) {
        cooked->ni = GDES_INT_MISSING;
        cooked->di = GDES_FLOAT_MISSING;
        gd->quasi = QUASI_ROWS;  /* A quasi-regular grid with varying rows */
    } else if (cooked->nj == G2I_MISSING && g2i(raw->dj) == G2I_MISSING) {
        cooked->nj = GDES_INT_MISSING;
        cooked->dj = GDES_FLOAT_MISSING;
        gd->quasi = QUASI_COLS; /* A quasi-regular grid with varying cols */
    }
    gd->ncols = cooked->ni;
    gd->nrows = cooked->nj;
    gd->npts = gd->ncols*gd->nrows;
    gd->res_flags = raw->res_flags;
    gd->scan_mode = raw->scan_mode;
    return 0;
}

static int
fill_srll(raw,gd)
    grid_srll *raw;
    gdes *gd;
{
    gdes_ll *cooked = &gd->grid.ll;
    cooked->ni = g2i(raw->ni);
    cooked->nj = g2i(raw->nj);
    cooked->la1 = g3si(raw->la1)*.001; /* millidegrees to degrees */
    cooked->lo1 = g3si(raw->lo1)*.001;
    cooked->la2 = g3si(raw->la2)*.001;
    cooked->lo2 = g3si(raw->lo2)*.001;
    cooked->di = g2i(raw->di)*.001;
    cooked->dj = g2i(raw->dj)*.001;
    cooked->rot=(rotated *)emalloc(sizeof(rotated));
    cooked->rot->lat = g3si(raw->lapole)*.001;
    cooked->rot->lon = g3si(raw->lopole)*.001;
    cooked->rot->angle = g4f(raw->angrot);
    cooked->strch=(stretched *)emalloc(sizeof(stretched));
    cooked->strch->lat = g3si(raw->lastr)*.001;
    cooked->strch->lon = g3si(raw->lostr)*.001;
    cooked->strch->factor = g4f(raw->stretch);

    if(cooked->ni == G2I_MISSING && g2i(raw->di) == G2I_MISSING) {
        cooked->ni = GDES_INT_MISSING;
        cooked->di = GDES_FLOAT_MISSING;
        gd->quasi = QUASI_ROWS;  /* A quasi-regular grid with varying rows */
    } else if (cooked->nj == G2I_MISSING && g2i(raw->dj) == G2I_MISSING) {
        cooked->nj = GDES_INT_MISSING;
        cooked->dj = GDES_FLOAT_MISSING;
        gd->quasi = QUASI_COLS;  /* A quasi-regular grid with varying cols */
    }
    gd->ncols = cooked->ni;
    gd->nrows = cooked->nj;
    gd->npts = gd->ncols*gd->nrows;
    gd->res_flags = raw->res_flags;
    gd->scan_mode = raw->scan_mode;
    return 0;
}

static int
fill_gau(raw,gd)
    grid_gau *raw;
    gdes *gd;
{
    gdes_gau *cooked = &gd->grid.gau;
    cooked->ni = g2i(raw->ni);
    cooked->nj = g2i(raw->nj);
    cooked->la1 = g3si(raw->la1)*.001; /* millidegrees to degrees */
    cooked->lo1 = g3si(raw->lo1)*.001;
    cooked->la2 = g3si(raw->la2)*.001;
    cooked->lo2 = g3si(raw->lo2)*.001;
    cooked->di = g2i(raw->di)*.001;
    cooked->n = g2i(raw->n);
    cooked->rot=0;
    cooked->strch=0;

    if(cooked->ni == G2I_MISSING && g2i(raw->di) == G2I_MISSING) {
        cooked->ni = GDES_INT_MISSING;
        cooked->di = GDES_FLOAT_MISSING;
        gd->quasi = QUASI_ROWS;  /* A quasi-regular grid with varying rows */
    }
    gd->ncols = cooked->ni;
    gd->nrows = cooked->nj;
    gd->npts = gd->ncols*gd->nrows;
    gd->res_flags = raw->res_flags;
    gd->scan_mode = raw->scan_mode;
    return 0;
}

static int
fill_rgau(raw,gd)
    grid_rgau *raw;
    gdes *gd;
{
    gdes_gau *cooked = &gd->grid.gau;
    cooked->ni = g2i(raw->ni);
    cooked->nj = g2i(raw->nj);
    cooked->la1 = g3si(raw->la1)*.001; /* millidegrees to degrees */
    cooked->lo1 = g3si(raw->lo1)*.001;
    cooked->la2 = g3si(raw->la2)*.001;
    cooked->lo2 = g3si(raw->lo2)*.001;
    cooked->di = g2i(raw->di)*.001;
    cooked->n = g2i(raw->n);
    cooked->rot=(rotated *)emalloc(sizeof(rotated));
    cooked->rot->lat = g3si(raw->lapole)*.001;
    cooked->rot->lon = g3si(raw->lopole)*.001;
    cooked->rot->angle = g4f(raw->angrot);
    cooked->strch = 0;

    if(cooked->ni == G2I_MISSING && g2i(raw->di) == G2I_MISSING) {
        cooked->ni = GDES_INT_MISSING;
        cooked->di = GDES_FLOAT_MISSING;
        gd->quasi = QUASI_ROWS; /* A quasi-regular grid with varying rows */
    }
    gd->ncols = cooked->ni;
    gd->nrows = cooked->nj;
    gd->npts = gd->ncols*gd->nrows;
    gd->res_flags = raw->res_flags;
    gd->scan_mode = raw->scan_mode;
    return 0;
}

static int
fill_sgau(raw,gd)
    grid_sgau *raw;
    gdes *gd;
{
    gdes_gau *cooked = &gd->grid.gau;
    cooked->ni = g2i(raw->ni);
    cooked->nj = g2i(raw->nj);
    cooked->la1 = g3si(raw->la1)*.001; /* millidegrees to degrees */
    cooked->lo1 = g3si(raw->lo1)*.001;
    cooked->la2 = g3si(raw->la2)*.001;
    cooked->lo2 = g3si(raw->lo2)*.001;
    cooked->di = g2i(raw->di)*.001;
    cooked->n = g2i(raw->n);
    cooked->rot=0;
    cooked->strch=(stretched *)emalloc(sizeof(stretched));
    cooked->strch->lat = g3si(raw->lastr)*.001;
    cooked->strch->lon = g3si(raw->lostr)*.001;
    cooked->strch->factor = g4f(raw->stretch);

    if(cooked->ni == G2I_MISSING && g2i(raw->di) == G2I_MISSING) {
        cooked->ni = GDES_INT_MISSING;
        cooked->di = GDES_FLOAT_MISSING;
        gd->quasi = QUASI_ROWS;   /* A quasi-regular grid with varying rows */
    }
    gd->ncols = cooked->ni;
    gd->nrows = cooked->nj;
    gd->npts = gd->ncols*gd->nrows;
    gd->res_flags = raw->res_flags;
    gd->scan_mode = raw->scan_mode;
    return 0;
}

static int
fill_srgau(raw,gd)
    grid_srgau *raw;
    gdes *gd;
{
    gdes_gau *cooked = &gd->grid.gau;
    cooked->ni = g2i(raw->ni);
    cooked->nj = g2i(raw->nj);
    cooked->la1 = g3si(raw->la1)*.001; /* millidegrees to degrees */
    cooked->lo1 = g3si(raw->lo1)*.001;
    cooked->la2 = g3si(raw->la2)*.001;
    cooked->lo2 = g3si(raw->lo2)*.001;
    cooked->di = g2i(raw->di)*.001;
    cooked->n = g2i(raw->n);
    cooked->rot=(rotated *)emalloc(sizeof(rotated));
    cooked->rot->lat = g3si(raw->lapole)*.001;
    cooked->rot->lon = g3si(raw->lopole)*.001;
    cooked->rot->angle = g4f(raw->angrot);
    cooked->strch=(stretched *)emalloc(sizeof(stretched));
    cooked->strch->lat = g3si(raw->lastr)*.001;
    cooked->strch->lon = g3si(raw->lostr)*.001;
    cooked->strch->factor = g4f(raw->stretch);

    if(cooked->ni == G2I_MISSING && g2i(raw->di) == G2I_MISSING) {
        cooked->ni = GDES_INT_MISSING;
        cooked->di = GDES_FLOAT_MISSING;
        gd->quasi = QUASI_ROWS; /* A quasi-regular grid with varying rows */
    }
    gd->ncols = cooked->ni;
    gd->nrows = cooked->nj;
    gd->npts = gd->ncols*gd->nrows;
    gd->res_flags = raw->res_flags;
    gd->scan_mode = raw->scan_mode;
    return 0;
}

static int
fill_sph(raw,gd)
    grid_sph *raw;
    gdes *gd;
{
    gdes_sph *cooked = &gd->grid.sph;
    cooked->j = g2i(raw->j);
    cooked->k = g2i(raw->k);
    cooked->m = g2i(raw->m);
    cooked->type = g1i(raw->type);
    cooked->mode = g1i(raw->mode);
    cooked->rot=0;
    cooked->strch=0;
    gd->ncols = cooked->j;        /* *** probably not right ? *** */
    gd->nrows = cooked->k;
    gd->npts = gd->ncols*gd->nrows;
    gd->res_flags = 0;
    gd->scan_mode = 0;
    return 0;
}

static int
fill_rsph(raw,gd)
    grid_rsph *raw;
    gdes *gd;
{
    gdes_sph *cooked = &gd->grid.sph;
    cooked->j = g2i(raw->j);
    cooked->k = g2i(raw->k);
    cooked->m = g2i(raw->m);
    cooked->type = g1i(raw->type);
    cooked->mode = g1i(raw->mode);
    cooked->rot=(rotated *)emalloc(sizeof(rotated));
    cooked->rot->lat = g3si(raw->lapole)*.001;
    cooked->rot->lon = g3si(raw->lopole)*.001;
    cooked->rot->angle = g4f(raw->angrot);
    cooked->strch=0;
    gd->ncols = cooked->j;     /* *** probably not right *** */
    gd->nrows = cooked->k;
    gd->npts = gd->ncols*gd->nrows;
    gd->res_flags = 0;
    gd->scan_mode = 0;
    return 0;
}

static int
fill_ssph(raw,gd)
    grid_ssph *raw;
    gdes *gd;
{
    gdes_sph *cooked = &gd->grid.sph;
    cooked->j = g2i(raw->j);
    cooked->k = g2i(raw->k);
    cooked->m = g2i(raw->m);
    cooked->type = g1i(raw->type);
    cooked->mode = g1i(raw->mode);
    cooked->rot=0;
    cooked->strch=(stretched *)emalloc(sizeof(stretched));
    cooked->strch->lat = g3si(raw->lastr)*.001;
    cooked->strch->lon = g3si(raw->lostr)*.001;
    cooked->strch->factor = g4f(raw->stretch);
    gd->ncols = cooked->j;        /* *** probably not right *** */
    gd->nrows = cooked->k;
    gd->npts = gd->ncols*gd->nrows;
    gd->res_flags = 0;
    gd->scan_mode = 0;
    return 0;
}

static int
fill_srsph(raw,gd)
    grid_srsph *raw;
    gdes *gd;
{
    gdes_sph *cooked = &gd->grid.sph;
    cooked->j = g2i(raw->j);
    cooked->k = g2i(raw->k);
    cooked->m = g2i(raw->m);
    cooked->type = g1i(raw->type);
    cooked->mode = g1i(raw->mode);
    cooked->rot=(rotated *)emalloc(sizeof(rotated));
    cooked->rot->lat = g3si(raw->lapole)*.001;
    cooked->rot->lon = g3si(raw->lopole)*.001;
    cooked->rot->angle = g4f(raw->angrot);
    cooked->strch=(stretched *)emalloc(sizeof(stretched));
    cooked->strch->lat = g3si(raw->lastr)*.001;
    cooked->strch->lon = g3si(raw->lostr)*.001;
    cooked->strch->factor = g4f(raw->stretch);
    gd->ncols = cooked->j;        /* *** probably not right *** */
    gd->nrows = cooked->k;
    gd->npts = gd->ncols*gd->nrows;
    gd->res_flags = 0;
    gd->scan_mode = 0;
    return 0;
}

static int
fill_mercator(raw,gd)
    grid_mercator *raw;
    gdes *gd;
{
    gdes_mercator *cooked = &gd->grid.mercator;
    cooked->ni = g2i(raw->ni);
    cooked->nj = g2i(raw->nj);
    cooked->la1 = g3si(raw->la1)*.001; /* millidegrees to degrees */
    cooked->lo1 = g3si(raw->lo1)*.001;
    cooked->la2 = g3si(raw->la2)*.001;
    cooked->lo2 = g3si(raw->lo2)*.001;
    cooked->latin = g3si(raw->latin)*.001;
    cooked->di = g3i(raw->di);        /* in meters */
    cooked->dj = g3i(raw->dj);        /* in meters */
    gd->ncols = cooked->ni;
    gd->nrows = cooked->nj;
    gd->npts = gd->ncols*gd->nrows;
    gd->res_flags = raw->res_flags;
    gd->scan_mode = raw->scan_mode;
    return 0;
}

static int
fill_polars(raw,gd)
    grid_polars *raw;
    gdes *gd;
{
    gdes_polars *cooked = &gd->grid.polars;
    cooked->nx = g2i(raw->nx);
    cooked->ny = g2i(raw->ny);
    cooked->la1 = g3si(raw->la1)*.001; /* millidegrees to degrees */
    cooked->lo1 = g3si(raw->lo1)*.001;
    cooked->lov = g3si(raw->lov)*.001;
    cooked->dx = g3i(raw->dx);        /* in meters */
    cooked->dy = g3i(raw->dy);        /* in meters */
    cooked->pole = ((raw->pole & 0x80) == 0x80);
    gd->ncols = cooked->nx;
    gd->nrows = cooked->ny;
    gd->npts = gd->ncols*gd->nrows;
    gd->res_flags = raw->res_flags;
    gd->scan_mode = raw->scan_mode;
    return 0;
}

static int
fill_gnomon(raw,gd)
    grid_gnomon *raw;
    gdes *gd;
{
    gdes_polars *cooked = &gd->grid.polars;
    cooked->nx = g2i(raw->nx);
    cooked->ny = g2i(raw->ny);
    cooked->la1 = g3si(raw->la1)*.001; /* millidegrees to degrees */
    cooked->lo1 = g3si(raw->lo1)*.001;
    cooked->lov = g3si(raw->lov)*.001;
    cooked->dx = g3i(raw->dx);        /* in meters */
    cooked->dy = g3i(raw->dy);        /* in meters */
    cooked->pole = ((raw->pole & 0x80) == 0x80);
    gd->ncols = cooked->nx;
    gd->nrows = cooked->ny;
    gd->npts = gd->ncols*gd->nrows;
    gd->res_flags = raw->res_flags;
    gd->scan_mode = raw->scan_mode;
    return 0;
}

static int
fill_lambert(raw,gd)
    grid_lambert *raw;
    gdes *gd;
{
    gdes_lambert *cooked = &gd->grid.lambert;
    cooked->nx = g2i(raw->nx);
    cooked->ny = g2i(raw->ny);
    cooked->la1 = g3si(raw->la1)*.001; /* millidegrees to degrees */
    cooked->lo1 = g3si(raw->lo1)*.001;
    cooked->lov = g3si(raw->lov)*.001;
    cooked->dx = g3i(raw->dx);        /* in meters */
    cooked->dy = g3i(raw->dy);        /* in meters */
    cooked->pole = ((raw->pole & 0x80) == 0x80);
    cooked->centers = ((raw->pole & 0x40) == 0x40) + 1;
    cooked->latin1 = g3si(raw->latin1)*.001;
    cooked->latin2 = g3si(raw->latin2)*.001;
    cooked->splat = g3si(raw->splat)*.001;
    cooked->splon = g3si(raw->splon)*.001;
    gd->ncols = cooked->nx;
    gd->nrows = cooked->ny;
    gd->npts = gd->ncols*gd->nrows;
    gd->res_flags = raw->res_flags;
    gd->scan_mode = raw->scan_mode;
    return 0;
}

static int
fill_olambert(raw,gd)
    grid_olambert *raw;
    gdes *gd;
{
    gdes_lambert *cooked = &gd->grid.lambert;
    cooked->nx = g2i(raw->nx);
    cooked->ny = g2i(raw->ny);
    cooked->la1 = g3si(raw->la1)*.001; /* millidegrees to degrees */
    cooked->lo1 = g3si(raw->lo1)*.001;
    cooked->lov = g3si(raw->lov)*.001;
    cooked->dx = g3i(raw->dx);        /* in meters */
    cooked->dy = g3i(raw->dy);        /* in meters */
    cooked->pole = ((raw->pole & 0x80) == 0x80);
    cooked->centers = ((raw->pole & 0x40) == 0x40) + 1;
    cooked->latin1 = g3si(raw->latin1)*.001;
    cooked->latin2 = g3si(raw->latin2)*.001;
    cooked->splat = g3si(raw->splat)*.001;
    cooked->splon = g3si(raw->splon)*.001;
    gd->ncols = cooked->nx;
    gd->nrows = cooked->ny;
    gd->npts = gd->ncols*gd->nrows;
    gd->res_flags = raw->res_flags;
    gd->scan_mode = raw->scan_mode;
    return 0;
}

static int
fill_albers(raw,gd)
    grid_albers *raw;
    gdes *gd;
{
    gdes_lambert *cooked = &gd->grid.lambert;
    cooked->nx = g2i(raw->nx);
    cooked->ny = g2i(raw->ny);
    cooked->la1 = g3si(raw->la1)*.001; /* millidegrees to degrees */
    cooked->lo1 = g3si(raw->lo1)*.001;
    cooked->lov = g3si(raw->lov)*.001;
    cooked->dx = g3i(raw->dx);        /* in meters */
    cooked->dy = g3i(raw->dy);        /* in meters */
    cooked->pole = ((raw->pole & 0x80) == 0x80);
    cooked->centers = ((raw->pole & 0x40) == 0x40) + 1;
    cooked->latin1 = g3si(raw->latin1)*.001;
    cooked->latin2 = g3si(raw->latin2)*.001;
    cooked->splat = g3si(raw->splat)*.001;
    cooked->splon = g3si(raw->splon)*.001;
    gd->ncols = cooked->nx;
    gd->nrows = cooked->ny;
    gd->npts = gd->ncols*gd->nrows;
    gd->res_flags = raw->res_flags;
    gd->scan_mode = raw->scan_mode;
    return 0;
}

static int
fill_spacev(raw,gd)
    grid_spacev *raw;
    gdes *gd;
{
    gdes_spacev *cooked = &gd->grid.spacev;
    cooked->nx = g2i(raw->nx);
    cooked->ny = g2i(raw->ny);
    cooked->lap = g3si(raw->lap)*.001; /* millidegrees to degrees */
    cooked->lop = g3si(raw->lop)*.001;
    cooked->dx = g3i(raw->dx);        /* in grid-lengths */
    cooked->dy = g3i(raw->dy);        /* in grid-lengths */
    cooked->xp = g2i(raw->xp);
    cooked->yp = g2i(raw->yp);
    cooked->orient = g3si(raw->orient)*.001; /* millidegrees to degrees */
    cooked->nr = g3i(raw->nr);
    cooked->xo = g2i(raw->xo);
    cooked->yo = g2i(raw->yo);
    gd->ncols = cooked->nx;
    gd->nrows = cooked->ny;
    gd->npts = gd->ncols*gd->nrows;
    gd->res_flags = raw->res_flags;
    gd->scan_mode = raw->scan_mode;
    return 0;
}

static void copy_gdes_ll(gdes_ll* gd, gdes_ll* g_ll){
/* This function take a deep copy from gd to g_ll */ 
    *g_ll=*gd;

    if( gd->rot!=0 ){
        g_ll->rot=(rotated*)malloc(sizeof(rotated));
        *(g_ll->rot)=*(gd->rot);
    }

    if( gd->strch!=0 ){
        g_ll->strch=(stretched*)malloc(sizeof(stretched));
        *(g_ll->strch)=*(gd->strch);
    }
    
return;
} /* end of  copy_gdes_ll */

static void copy_gdes_gau(gdes_gau* gd, gdes_gau* g_gau){
/* This function take a deep copy from gd to g_gau */ 
    *g_gau=*gd;

    if( gd->rot!=0 ){
        g_gau->rot=(rotated*)malloc(sizeof(rotated));
        *(g_gau->rot)=*(gd->rot);
    }

    if( gd->strch!=0 ){
        g_gau->strch=(stretched*)malloc(sizeof(stretched));
        *(g_gau->strch)=*(gd->strch);
    }
    
return;
}

static void copy_gdes_sph(gdes_sph* gd, gdes_sph* g_sph){
/* This function take a deep copy from gd to g_sph */ 
    *g_sph=*gd;

    if( gd->rot!=0 ){
        g_sph->rot=(rotated*)malloc(sizeof(rotated));
        *(g_sph->rot)=*(gd->rot);
    }

    if( gd->strch!=0 ){
        g_sph->strch=(stretched*)malloc(sizeof(stretched));
        *(g_sph->strch)=*(gd->strch);
    }
    
return;
}

/* This function take a deep copy from gd to g_mercat */ 
static void copy_gdes_mercat( gdes_mercator* gd,gdes_mercator* g_mercat ){
    *g_mercat=*gd;return;
}

/* This function take a deep copy from gd to g_polars */ 
static void copy_gdes_polars( gdes_polars* gd,gdes_polars* g_polars ){
    *g_polars=*gd;return;
}

/* This function take a deep copy from gd to g_lambert */ 
static void copy_gdes_lambert( gdes_lambert* gd, gdes_lambert* g_lambert){
    *g_lambert=*gd;return;
}

/* This function take a deep copy from gd to g_spacev */ 
static void copy_gdes_spacev( gdes_spacev * gd, gdes_spacev* g_spacev){
    *g_spacev=*gd;return;
}        

/*
 * Make a gdes from an existing raw gds.
 * User should call free_gdes() when done to free.
 */
static gdes*
gds_to_gdes(gdsp)
    gds *gdsp;
{
    int nv = g1i(gdsp->nv);
    int pv = g1i(gdsp->pv);
    int type = g1i(gdsp->type);
    int nerrs = 0;
    gdes *ret = (gdes *) emalloc(sizeof(gdes));

    ret->type = type;
    ret->quasi = QUASI_RECT;        /* ordinary rectangular grid is default */
    ret->nv = 0;
    ret->vc = 0;
    ret->lc = 0;
    ret->keep = 0;
    if(nv != 0 && nv != G1I_MISSING) { /* we have vert. coordinates */
        g4flt *fp = (g4flt *) ((char *)gdsp + pv);
        int i;
        
        ret->nv = nv;
        ret->vc = (float *) emalloc(ret->nv * sizeof(float));
        /* unpack the vertical coords into floats */
        for (i = 0; i < ret->nv; i++) {
            ret->vc[i] = g4f(*fp++);
        }
    }

    switch(type) {
    case GRID_LL:
        nerrs = fill_ll(&gdsp->grid.ll, ret);
        break;
    case GRID_RLL:
        nerrs = fill_rll(&gdsp->grid.rll, ret);
        break;
    case GRID_SLL:
        nerrs = fill_sll(&gdsp->grid.sll, ret);
        break;
    case GRID_SRLL:
        nerrs = fill_srll(&gdsp->grid.srll, ret);
        break;
    case GRID_GAU:
        nerrs = fill_gau(&gdsp->grid.gau, ret);
        break;
    case GRID_RGAU:
        nerrs = fill_rgau(&gdsp->grid.rgau, ret);
        break;
    case GRID_SGAU:
        nerrs = fill_sgau(&gdsp->grid.sgau, ret);
        break;
    case GRID_SRGAU:
        nerrs = fill_srgau(&gdsp->grid.srgau, ret);
        break;
    case GRID_SPH:
        nerrs = fill_sph(&gdsp->grid.sph, ret);
        break;
    case GRID_RSPH:
        nerrs = fill_rsph(&gdsp->grid.rsph, ret);
        break;
    case GRID_SSPH:
        nerrs = fill_ssph(&gdsp->grid.ssph, ret);
        break;
    case GRID_SRSPH:
        nerrs = fill_srsph(&gdsp->grid.srsph, ret);
        break;
    case GRID_MERCAT:
        nerrs = fill_mercator(&gdsp->grid.mercator, ret);
        break;
    case GRID_POLARS:
        nerrs = fill_polars(&gdsp->grid.polars, ret);
        break;
    case GRID_GNOMON:
        nerrs = fill_gnomon(&gdsp->grid.gnomon, ret);
        break;
    case GRID_LAMBERT:
        nerrs = fill_lambert(&gdsp->grid.lambert, ret);
        break;
    case GRID_ALBERS:
        nerrs = fill_albers(&gdsp->grid.albers, ret);
        break;
    case GRID_OLAMBERT:
        nerrs = fill_olambert(&gdsp->grid.olambert, ret);
        break;
    case GRID_SPACEV:
        nerrs = fill_spacev(&gdsp->grid.spacev, ret);
        break;
    case GRID_UTM:
    case GRID_SIMPOL:
    case GRID_MILLER:
    default:
        nerrs++;
    }
    if (nerrs) {
        free_gdes(ret);
        ret = 0;
        return ret;
    }

    if(ret->quasi == QUASI_ROWS) {
        ret->ncols = 1;
        ret->lc = (int *)emalloc((1 + ret->nrows) * sizeof(int));
        {        /* unpack list of row indexes */
            g2int *ip = (g2int *) ((char *)gdsp + g1i(gdsp->pv)-1 + 4*ret->nv);
            int i;
            int n=0;
            int maxlc = 0;
            for (i = 0; i < ret->nrows; i++) {
                int nl;
                nl = g2i(*ip++);
                if (nl > maxlc)
                    maxlc = nl;
                ret->lc[i] = n;
                n += nl;
            }
            ret->lc[i] = n;
            ret->npts = n;
            ret->maxlc = maxlc;
        }
    } else if(ret->quasi == QUASI_COLS) {
        ret->nrows = 1;
        ret->lc = (int *)emalloc((1 + ret->ncols) * sizeof(int));
        {        /* unpack list of col indexes */
            g2int *ip = (g2int *) ((char *)gdsp + g1i(gdsp->pv)-1 + 4*ret->nv);
            int i;
            int n=0;
            int maxlc = 0;
            for (i = 0; i < ret->ncols; i++) {
                int nl;
                nl = g2i(*ip++);
                if (nl > maxlc)
                    maxlc = nl;
                ret->lc[i] = n;
                n += nl;
            }
            ret->lc[i] = n;
            ret->npts = n;
            ret->maxlc = maxlc;
        }
    }
    if (nerrs) {
        free_gdes(ret);
        ret = 0;
    }
    return ret;
}

static void
nmc_21_24(gdes *g)
{
    g->grid.ll.ni = 37;
    g->grid.ll.nj = 37;
    g->grid.ll.di = 5.0;
    g->grid.ll.dj = 2.5;
    g->grid.ll.rot = 0;
    g->grid.ll.strch = 0;
    g->type = GRID_LL;
    g->ncols = 37;
    g->nrows = 37;
    g->npts = g->ncols*g->nrows;
    g->res_flags = RESCMP_DIRINC | RESCMP_UVRES;
    g->scan_mode = SCAN_J_PLUS;
    g->nv = 0;
    g->vc = 0;
    g->quasi = QUASI_RECT;
    g->lc = 0;
    g->keep = 1;
}

static gdes*
nmc_21()
{
    static gdes_ll *ll;
    static gdes g;

    nmc_21_24(&g);
    ll = &g.grid.ll;
    ll->la1 = 0;
    ll->lo1 = 0;
    ll->la2 = 90;
    ll->lo2 = 180;
    return &g;
}

static gdes*
nmc_22()
{
    static gdes_ll *ll;
    static gdes g;

    ll = &g.grid.ll;
    nmc_21_24(&g);
    ll->la1 = 0;
    ll->lo1 = -180;
    ll->la2 = 90;
    ll->lo2 = 0;
    return &g;
}

static gdes*
nmc_23()
{
    static gdes_ll *ll;
    static gdes g;

    ll = &g.grid.ll;
    nmc_21_24(&g);
    ll->la1 = -90;
    ll->lo1 = 0;
    ll->la2 = 0;
    ll->lo2 = 180;
    return &g;
}

static gdes*
nmc_24()
{
    static gdes_ll *ll;
    static gdes g;

    ll = &g.grid.ll;
    nmc_21_24(&g);
    ll->la1 = -90;
    ll->lo1 = -180;
    ll->la2 = 0;
    ll->lo2 = 0;
    return &g;
}

static void
nmc_25_26(gdes *g)
{
    g->grid.ll.ni = 72;
    g->grid.ll.nj = 19;
    g->grid.ll.di = 5.0;
    g->grid.ll.dj = 5.0;
    g->grid.ll.rot = 0;
    g->grid.ll.strch = 0;
    g->type = GRID_LL;
    g->ncols = 72;
    g->nrows = 19;
    g->npts = g->ncols*g->nrows;
    g->res_flags = RESCMP_DIRINC | RESCMP_UVRES;
    g->scan_mode = SCAN_J_PLUS;
    g->nv = 0;
    g->vc = 0;
    g->quasi = QUASI_RECT;
    g->lc = 0;
    g->keep = 1;
}

static gdes*
nmc_25()
{
    static gdes_ll *ll;
    static gdes g;

    ll = &g.grid.ll;
    nmc_25_26(&g);
    ll->la1 = 0;
    ll->lo1 = 0;
    ll->la2 = 90;
    ll->lo2 = 355;
    return &g;
}

static gdes*
nmc_26()
{
    static gdes_ll *ll;
    static gdes g;

    ll = &g.grid.ll;
    nmc_25_26(&g);
    ll->la1 = -90;
    ll->lo1 = 0;
    ll->la2 = 0;
    ll->lo2 = 355;
    return &g;
}

static gdes*
nmc_50()
{
    static gdes_ll *ll;
    static gdes g;

    ll = &g.grid.ll;

    ll->ni = 36;
    ll->nj = 33;
    ll->di = 2.5;
    ll->dj = 1.25;
    ll->la1 = 20;
    ll->lo1 = -140;
    ll->la2 = 60;
    ll->lo2 = -52.5;
    ll->rot = 0;
    ll->strch = 0;
    g.type = GRID_LL;
    g.ncols = 36;
    g.nrows = 33;
    g.npts = g.ncols*g.nrows;
    g.res_flags = RESCMP_DIRINC | RESCMP_UVRES;
    g.scan_mode = SCAN_J_PLUS;
    g.nv = 0;
    g.vc = 0;
    g.quasi = QUASI_RECT;
    g.lc = 0;
    g.keep = 1;
    return &g;
}

static void
nmc_61_64(gdes *g)
{
    g->grid.ll.ni = 91;
    g->grid.ll.nj = 46;
    g->grid.ll.di = 2.0;
    g->grid.ll.dj = 2.0;
    g->grid.ll.rot = 0;
    g->grid.ll.strch = 0;
    g->type = GRID_LL;
    g->ncols = 91;
    g->nrows = 46;
    g->npts = g->ncols*g->nrows;
    g->res_flags = RESCMP_DIRINC | RESCMP_UVRES;
    g->scan_mode = SCAN_J_PLUS;
    g->nv = 0;
    g->vc = 0;
    g->quasi = QUASI_RECT;
    g->lc = 0;
    g->keep = 1;
}

static gdes*
nmc_61()
{
    static gdes_ll *ll;
    static gdes g;

    ll = &g.grid.ll;
    nmc_61_64(&g);
    ll->la1 = 0;
    ll->lo1 = 0;
    ll->la2 = 90;
    ll->lo2 = 180;
    return &g;
}

static gdes*
nmc_62()
{
    static gdes_ll *ll;
    static gdes g;

    ll = &g.grid.ll;
    nmc_61_64(&g);
    ll->la1 = 0;
    ll->lo1 = -180;
    ll->la2 = 90;
    ll->lo2 = 0;
    return &g;
}

static gdes*
nmc_63()
{
    static gdes_ll *ll;
    static gdes g;

    ll = &g.grid.ll;
    nmc_61_64(&g);
    ll->la1 = -90;
    ll->lo1 = 0;
    ll->la2 = 0;
    ll->lo2 = 180;
    return &g;
}

static gdes*
nmc_64()
{
    static gdes_ll *ll;
    static gdes g;

    ll = &g.grid.ll;
    nmc_61_64(&g);
    ll->la1 = -90;
    ll->lo1 = -180;
    ll->la2 = 0;
    ll->lo2 = 0;
    return &g;
}

static gdes*
nmc_1()                                /* Tropical strip (Mercator) */
{
    static gdes_mercator *mercator;
    static gdes g;

    mercator = &g.grid.mercator;

    mercator->ni = 73;
    mercator->nj = 23;
    mercator->di = 513669;
    mercator->dj = 513669;
    mercator->la1 = -48.09;
    mercator->lo1 = 0;
    mercator->la2 = 48.09;
    mercator->lo2 = 360.;
    mercator->latin = 22.5;

    g.type = GRID_MERCAT;
    g.ncols = 73;
    g.nrows = 23;
    g.npts = g.ncols*g.nrows;
    g.res_flags = RESCMP_DIRINC | RESCMP_UVRES;
    g.scan_mode = SCAN_J_PLUS;
    g.nv = 0;
    g.vc = 0;
    g.quasi = QUASI_RECT;
    g.lc = 0;
    g.keep = 1;
    return &g;
}

static gdes*
nmc_2()
{
    static gdes_ll *ll;
    static gdes g;

    ll = &g.grid.ll;

    ll->ni = 144;
    ll->nj = 73;
    ll->di = 2.5;
    ll->dj = 2.5;
    ll->la1 = 90;
    ll->lo1 = 0;
    ll->la2 = -90;
    ll->lo2 = 355;
    ll->rot = 0;
    ll->strch = 0;
    g.type = GRID_LL;
    g.ncols = 144;
    g.nrows = 73;
    g.npts = g.ncols*g.nrows;
    g.res_flags = RESCMP_DIRINC | RESCMP_UVRES;
    g.scan_mode = 0;
    g.nv = 0;
    g.vc = 0;
    g.quasi = QUASI_RECT;
    g.lc = 0;
    g.keep = 1;
    return &g;
}

static gdes*
nmc_3()
{
    static gdes_ll *ll;
    static gdes g;

    ll = &g.grid.ll;

    ll->ni = 360;
    ll->nj = 181;
    ll->di = 1.0;
    ll->dj = 1.0;
    ll->la1 = 90;
    ll->lo1 = 0;
    ll->la2 = -90;
    ll->lo2 = 359;
    ll->rot = 0;
    ll->strch = 0;
    g.type = GRID_LL;
    g.ncols = ll->ni;
    g.nrows = ll->nj;
    g.npts = g.ncols*g.nrows;
    g.res_flags = RESCMP_DIRINC ;
    g.scan_mode = 0;
    g.nv = 0;
    g.vc = 0;
    g.quasi = QUASI_RECT;
    g.lc = 0;
    g.keep = 1;
    return &g;
}

static gdes*
nmc_5()         /* N. Hemisphere polar stereographic oriented 105W */
{
    static gdes_polars *polars;
    static gdes g;

    polars = &g.grid.polars;

    polars->nx = 53;
    polars->ny = 57;
    polars->la1 = 7.64713;
    polars->lo1 = -133.443;
    polars->lov = -105;
    polars->dx = 109500.;
    polars->dy = 109500.;
    polars->pole = 0;
    g.type = GRID_POLARS;
    g.ncols = polars->nx;
    g.nrows = polars->ny;
    g.npts = g.ncols*g.nrows;
    g.res_flags = RESCMP_UVRES ;
    g.scan_mode = SCAN_J_PLUS ;
    g.nv = 0;
    g.vc = 0;
    g.quasi = QUASI_RECT;
    g.lc = 0;
    g.keep = 1;
    return &g;
}

static gdes*
nmc_6()            /* N. Hemisphere polar stereographic oriented 105W */
{
    static gdes_polars *polars;
    static gdes g;

    polars = &g.grid.polars;

    polars->nx = 53;
    polars->ny = 45;
    polars->la1 = 7.64713;
    polars->lo1 = -133.443;
    polars->lov = -105;
    polars->dx = 109500.;
    polars->dy = 109500.;
    polars->pole = 0;
    g.type = GRID_POLARS;
    g.ncols = polars->nx;
    g.nrows = polars->ny;
    g.npts = g.ncols*g.nrows;
    g.res_flags = RESCMP_UVRES ;
    g.scan_mode = SCAN_J_PLUS ;
    g.nv = 0;
    g.vc = 0;
    g.quasi = QUASI_RECT;
    g.lc = 0;
    g.keep = 1;
    return &g;
}

static gdes*
nmc_27()       /* N. Hemisphere polar stereographic oriented 80W */
{
    static gdes_polars *polars;
    static gdes g;

    polars = &g.grid.polars;

    polars->nx = 65;
    polars->ny = 65;
    polars->la1 = -20.8255;
    polars->lo1 = -125.;
    polars->lov = -80;
    polars->dx = 381000.;
    polars->dy = 381000.;
    polars->pole = 0;
    g.type = GRID_POLARS;
    g.ncols = polars->nx;
    g.nrows = polars->ny;
    g.npts = g.ncols*g.nrows;
    g.res_flags = RESCMP_UVRES ;
    g.scan_mode = SCAN_J_PLUS ;
    g.nv = 0;
    g.vc = 0;
    g.quasi = QUASI_RECT;
    g.lc = 0;
    g.keep = 1;
    return &g;
}

static gdes*
nmc_28()            /* S. Hemisphere polar stereographic oriented 100E */
{
    static gdes_polars *polars;
    static gdes g;

    polars = &g.grid.polars;

    polars->nx = 65;
    polars->ny = 65;
    polars->la1 = 20.8255;
    polars->lo1 = 145.;
    polars->lov = 100.;
    polars->dx = 381000.;
    polars->dy = 381000.;
    polars->pole = 1;
    g.type = GRID_POLARS;
    g.ncols = polars->nx;
    g.nrows = polars->ny;
    g.npts = g.ncols*g.nrows;
    g.res_flags = RESCMP_UVRES ;
    g.scan_mode = SCAN_J_PLUS ;
    g.nv = 0;
    g.vc = 0;
    g.quasi = QUASI_RECT;
    g.lc = 0;
    g.keep = 1;
    return &g;
}

static gdes*
nmc_100()            /* N. Hemisphere polar stereographic oriented 105W */
{
    static gdes_polars *polars;
    static gdes g;

    polars = &g.grid.polars;

    polars->nx = 83;
    polars->ny = 83;
    polars->la1 = 17.1101;
    polars->lo1 = -129.296;
    polars->lov = -105;
    polars->dx = 91452.;
    polars->dy = 91452.;
    polars->pole = 0;
    g.type = GRID_POLARS;
    g.ncols = polars->nx;
    g.nrows = polars->ny;
    g.npts = g.ncols*g.nrows;
    g.res_flags = RESCMP_UVRES ;
    g.scan_mode = SCAN_J_PLUS ;
    g.nv = 0;
    g.vc = 0;
    g.quasi = QUASI_RECT;
    g.lc = 0;
    g.keep = 1;
    return &g;
}

static gdes*
nmc_101()              /* N. Hemisphere polar stereographic oriented 105W */
{
    static gdes_polars *polars;
    static gdes g;

    polars = &g.grid.polars;

    polars->nx = 113;
    polars->ny = 91;
    polars->la1 = 10.52797;
    polars->lo1 = -137.146;
    polars->lov = -105;
    polars->dx = 91452.;
    polars->dy = 91452.;
    polars->pole = 0;
    g.type = GRID_POLARS;
    g.ncols = polars->nx;
    g.nrows = polars->ny;
    g.npts = g.ncols*g.nrows;
    g.res_flags = RESCMP_UVRES ;
    g.scan_mode = SCAN_J_PLUS ;
    g.nv = 0;
    g.vc = 0;
    g.quasi = QUASI_RECT;
    g.lc = 0;
    g.keep = 1;
    return &g;
}

static gdes*
nmc_104()            /* N. Hemisphere polar stereographic oriented 105W */
{
    static gdes_polars *polars;
    static gdes g;

    polars = &g.grid.polars;

    polars->nx = 147;
    polars->ny = 110;
    polars->la1 = -0.268327;
    polars->lo1 = -139.475;
    polars->lov = -105;
    polars->dx = 90754.64;
    polars->dy = 90754.64;
    polars->pole = 0;
    g.type = GRID_POLARS;
    g.ncols = polars->nx;
    g.nrows = polars->ny;
    g.npts = g.ncols*g.nrows;
    g.res_flags = RESCMP_UVRES ;
    g.scan_mode = SCAN_J_PLUS ;
    g.nv = 0;
    g.vc = 0;
    g.quasi = QUASI_RECT;
    g.lc = 0;
    g.keep = 1;
    return &g;
}

static gdes*
nmc_105()             /* N. Hemisphere polar stereographic oriented 105W */
{
    static gdes_polars *polars;
    static gdes g;

    polars = &g.grid.polars;

    polars->nx = 83;
    polars->ny = 83;
    polars->la1 = 17.529;
    polars->lo1 = -129.296;
    polars->lov = -105;
    polars->dx = 90754.64;
    polars->dy = 90754.64;
    polars->pole = 0;
    g.type = GRID_POLARS;
    g.ncols = polars->nx;
    g.nrows = polars->ny;
    g.npts = g.ncols*g.nrows;
    g.res_flags = RESCMP_UVRES ;
    g.scan_mode = SCAN_J_PLUS ;
    g.nv = 0;
    g.vc = 0;
    g.quasi = QUASI_RECT;
    g.lc = 0;
    g.keep = 1;
    return &g;
}

static gdes*
nmc_207()               /* Regional - Alaska (polar stereographic) */
{
    static gdes_polars *polars;
    static gdes g;

    polars = &g.grid.polars;

    polars->nx = 49;
    polars->ny = 35;
    polars->la1 = 42.085;
    polars->lo1 = -175.641;
    polars->lov = -150;
    polars->dx = 95250;
    polars->dy = 95250;
    polars->pole = 0;
    g.type = GRID_POLARS;
    g.ncols = polars->nx;
    g.nrows = polars->ny;
    g.npts = g.ncols*g.nrows;
    g.res_flags = RESCMP_UVRES ;
    g.scan_mode = SCAN_J_PLUS ;
    g.nv = 0;
    g.vc = 0;
    g.quasi = QUASI_RECT;
    g.lc = 0;
    g.keep = 1;
    return &g;
}

static gdes*
nmc_211()                        /* Regional - CONUS (Lambert Conformal) */
{
    static gdes_lambert *lambert;
    static gdes g;

    lambert = &g.grid.lambert;

    lambert->nx = 93;
    lambert->ny = 65;
    lambert->la1 = 12.190;
    lambert->lo1 = -133.459;
    lambert->lov = -95;
    lambert->dx = 81270.5;
    lambert->dy = 81270.5;
    lambert->pole = 0;
    lambert->centers = 1;        /* not bipolar */
    lambert->latin1 = 25.0;
    lambert->latin2 = 25.0;        /* tangent cone */
    g.type = GRID_LAMBERT;
    g.ncols = lambert->nx;
    g.nrows = lambert->ny;
    g.npts = g.ncols*g.nrows;
    g.res_flags = RESCMP_UVRES ;
    g.scan_mode = SCAN_J_PLUS ;
    g.nv = 0;
    g.vc = 0;
    g.quasi = QUASI_RECT;
    g.lc = 0;
    g.keep = 1;
    return &g;
}

static gdes*
nmc_212()                        /* Regional - CONUS (Lambert Conformal) */
{
    static gdes_lambert *lambert;
    static gdes g;

    lambert = &g.grid.lambert;

    lambert->nx = 185;
    lambert->ny = 129;
    lambert->la1 = 12.190;
    lambert->lo1 = -133.459;
    lambert->lov = -95;
    lambert->dx = 40635;
    lambert->dy = 40635;
    lambert->pole = 0;
    lambert->centers = 1;        /* not bipolar */
    lambert->latin1 = 25.0;
    lambert->latin2 = 25.0;        /* tangent cone */
    g.type = GRID_LAMBERT;
    g.ncols = lambert->nx;
    g.nrows = lambert->ny;
    g.npts = g.ncols*g.nrows;
    g.res_flags = RESCMP_UVRES ;
    g.scan_mode = SCAN_J_PLUS ;
    g.nv = 0;
    g.vc = 0;
    g.quasi = QUASI_RECT;
    g.lc = 0;
    g.keep = 1;
    return &g;
}

/*
 * Sythesize a gdes from the center ID and grid ID.
 * Returns 0 if failed.
 */
static gdes*
synth_gdes(centerid, gridid)
    int centerid;
    int gridid;
{
    /* If it's an international exchange grid, it doesn't matter what center
       it's from */
    switch(gridid) {
    case 21: return nmc_21();
    case 22: return nmc_22();
    case 23: return nmc_23();
    case 24: return nmc_24();
    case 25: return nmc_25();
    case 26: return nmc_26();
    case 50: return nmc_50();
    case 61: return nmc_61();
    case 62: return nmc_62();
    case 63: return nmc_63();
    case 64: return nmc_64();
    }

    switch(centerid) {
    case CENTER_NMC:
        switch(gridid) {

         case  1: return nmc_1();
         case  2: return nmc_2();
         case  3: return nmc_3(); /* and so on ..., but we don't see these on
                                     HRS, so we'll finish them later */
/*         case  4: return nmc_4(); */
         case  5: return nmc_5();
         case  6: return nmc_6();
         case  27: return nmc_27();
         case  28: return nmc_28();
/*         case  29: return nmc_29(); */
/*         case  30: return nmc_30(); */
/*         case  33: return nmc_33(); */
/*         case  34: return nmc_34(); */
/*         case  45: return nmc_45(); */
/*         case  55: return nmc_55(); */
/*         case  56: return nmc_56(); */
/*         case  67: return nmc_67(); */
/*         case  68: return nmc_68(); */
/*         case  69: return nmc_69(); */
/*         case  70: return nmc_70(); */
/*         case  71: return nmc_71(); */
/*         case  72: return nmc_72(); */
/*         case  73: return nmc_73(); */
/*         case  74: return nmc_74(); */
/*         case  75: return nmc_75(); */
/*         case  76: return nmc_76(); */
/*         case  77: return nmc_77(); */
/*         case  85: return nmc_85(); */
/*         case  86: return nmc_86(); */
/*         case  87: return nmc_87(); */
/*         case  90: return nmc_90(); */
/*         case  91: return nmc_91(); */
/*         case  92: return nmc_92(); */
/*         case  93: return nmc_93(); */
/*         case  94: return nmc_94(); */
/*         case  95: return nmc_95(); */
/*         case  98: return nmc_98(); */
        case  100: return nmc_100();
        case  101: return nmc_101();
/*         case  103: return nmc_103(); */
        case  104: return nmc_104();
        case  105: return nmc_105();
/*         case  106: return nmc_106(); */
/*         case  107: return nmc_107(); */
/*         case  126: return nmc_126(); */
/*         case  201: return nmc_201(); */
/*         case  202: return nmc_202(); */
/*         case  203: return nmc_203(); */
/*         case  204: return nmc_204(); */
/*         case  205: return nmc_205(); */
/*         case  206: return nmc_206(); */
        case  207: return nmc_207();
/*         case  208: return nmc_208(); */
/*         case  209: return nmc_209(); */
/*         case  210: return nmc_210(); */
        case  211: return nmc_211();
        case  212: return nmc_212();
/*         case  213: return nmc_213(); */
/*         case  214: return nmc_214(); */
        default:
            uerror("unrecognized NMC grid id %d", gridid );
            return 0;
        }
    case CENTER_FNOC:
        switch (gridid) {
/*         case  220: return fnoc_220(); */
/*         case  221: return fnoc_221(); */
/*         case  223: return fnoc_223(); */
        default:
            uerror("unrecognized FNOC grid id %d", gridid );
            return 0;
        }
    case CENTER_ECMWF:
        switch (gridid) {        /* These come with Grid Description Sections */
/*         case 1: return ecmwf_1(); */
/*         case 2: return ecmwf_2(); */
/*         case 3: return ecmwf_3(); */
/*         case 4: return ecmwf_4(); */
/*         case 5: return ecmwf_5(); */
/*         case 6: return ecmwf_6(); */
/*         case 7: return ecmwf_7(); */
/*         case 8: return ecmwf_8(); */
/*         case 9: return ecmwf_9(); */
/*         case 10: return ecmwf_10(); */
/*         case 11: return ecmwf_11(); */
/*         case 12: return ecmwf_12(); */
/*         case 13: return ecmwf_13(); */
/*         case 14: return ecmwf_14(); */
/*         case 15: return ecmwf_15(); */
/*         case 16: return ecmwf_16(); */
        default:
            uerror("unrecognized ECMWF grid id %d", gridid );
            return 0;
        }
    default:
        uerror("unrecognized (center,grid) combination: (%d,%d)",
               centerid, gridid);
        return 0;
    }
}

/*
 * Make a gdes from raw gds or pds.  Returns 0 on failure.
 */
gdes*
make_gdes(gb)
    grib1 *gb;
{
    if (gb->gdsp) {                /* If GDS exists, use it */
        TEST_NOSERC("make_gdes",0 );
        return gds_to_gdes(gb->gdsp);
    } else if(gb->pdsp->grid == NONCATALOGED_GRID) {
        uerror("grid id = 255, but no Grid Description Section");
        return 0;
    } else {      /* Otherwise, manufacture from PDS center, model */
        TEST_NOSERC("make_gdes",1 );
        return synth_gdes(g1i(gb->pdsp->center), g1i(gb->pdsp->grid));
    }
}

/*
 * Free gdes, unless statically allocated
 */
void
free_gdes(gd)
    gdes *gd;
{

    if (gd) {
        if(gd->keep)
            return;
        if (gd->vc)
            free(gd->vc);
        if (gd->lc)
            free(gd->lc);
        switch(gd->type) {        /* free type-specific stuff */
        case GRID_LL:
            break;
        case GRID_RLL:
            if(gd->grid.ll.rot)
                free(gd->grid.ll.rot);
            break;
        case GRID_SLL:
            if(gd->grid.ll.strch)
                free(gd->grid.ll.strch);
            break;
        case GRID_SRLL:
            if(gd->grid.ll.rot)
                free(gd->grid.ll.rot);
            if(gd->grid.ll.strch)
                free(gd->grid.ll.strch);
            break;
        case GRID_GAU:
            break;
        case GRID_RGAU:
            if(gd->grid.gau.rot)
                free(gd->grid.gau.rot);
            break;
        case GRID_SGAU:
            if(gd->grid.gau.strch)
                free(gd->grid.gau.strch);
            break;
        case GRID_SRGAU:
            if(gd->grid.gau.rot)
                free(gd->grid.gau.rot);
            if(gd->grid.gau.strch)
                free(gd->grid.gau.strch);
            break;
        case GRID_SPH:
            break;
        case GRID_RSPH:
            if(gd->grid.sph.rot)
                free(gd->grid.sph.rot);
            break;
        case GRID_SSPH:
            if(gd->grid.sph.strch)
                free(gd->grid.sph.strch);
            break;
        case GRID_SRSPH:
            if(gd->grid.sph.rot)
                free(gd->grid.sph.rot);
            if(gd->grid.sph.strch)
                free(gd->grid.sph.strch);
            break;
        case GRID_MERCAT:
        case GRID_POLARS:
        case GRID_GNOMON:
        case GRID_LAMBERT:
        case GRID_ALBERS:
        case GRID_OLAMBERT:
        case GRID_SPACEV:
        case GRID_UTM:
        case GRID_SIMPOL:
        case GRID_MILLER:
        default:
            break;
        }
        free(gd);
    }
}

/*
 * return string describing type of grid projection
 */
char*
gds_typename(int type)
{
        switch(type) {
        case GRID_LL:
            return "Latitude/Longitude";
        case GRID_RLL:
            return "Rotated latitude/longitude";
        case GRID_SLL:
            return "Stretched latitude/longitude";
        case GRID_SRLL:
            return "Stretched and rotated latitude/longitude";
        case GRID_GAU:
            return "Gaussian latitude/longitude";
        case GRID_RGAU:
            return "Rotated Gaussian latitude/longitude";
        case GRID_SGAU:
            return "Stretched Gaussian latitude/longitude";
        case GRID_SRGAU:
            return "Stretched and rotated Gaussian latitude/longitude";
        case GRID_SPH:
            return "Spherical harmonic coefficients";
        case GRID_RSPH:
            return "Rotated spherical harmonics";
        case GRID_SSPH:
            return "Stretched spherical harmonics";
        case GRID_SRSPH:
            return "Stretched and rotated spherical harmonics";
        case GRID_MERCAT:
            return "Mercator projection";
        case GRID_POLARS:
            return "Polar stereographic projection";
        case GRID_GNOMON:
            return "Gnomonic projection";
        case GRID_LAMBERT:
            return "Lambert conformal projection";
        case GRID_ALBERS:
            return "Albers equal-area projection";
        case GRID_OLAMBERT:
            return "Oblique Lambert projection";
        case GRID_SPACEV:
            return "Space view";
        case GRID_UTM:
            return "Universal Transverse Mercator (UTM) projection";
        case GRID_SIMPOL:
            return "Simple polyconic projection";
        case GRID_MILLER:
            return "Miller's cylindrical projection";
        default:
            return "Unknown GRIB GDS data representation type";
        }
}

/*
 * Dump gdes in text form
 */
void
print_gdes(gd)
    gdes *gd;
{
    printf("   %24s : %d (%s)\n","GDS representation type",
           gd->type, gds_typename(gd->type));
    printf("   %24s : %d\n","Number of columns",
           gd->ncols);        
    printf("   %24s : %d\n","Number of rows",
           gd->nrows);        
    printf("   %24s : %d\n","Number of points",
           gd->npts);
    printf("   %24s : ","Kind of grid");
    switch(gd->quasi) {
    case QUASI_RECT:
        printf("rectangular\n");
        break;
    case QUASI_ROWS:
        printf("quasi-regular (varying rows)\n");
        printf("   %24s : ","Row lengths");
        {
            int ii, *ip;
            for (ii=0, ip=gd->lc; ii < gd->nrows; ii++) {
                printf("%d ", *(ip+1) - *ip);
                ip++;
                if(ii%16 == 15 && ii < gd->nrows-1)
                    printf("\n   %24s   ", "");
            }
            printf("\n");
        }
        break;
    case QUASI_COLS:
        printf("quasi-regular (varying columns)\n");
        {
            int ii, *ip;
            for (ii=0, ip=gd->lc; ii < gd->ncols; ii++) {
                printf("%d ", *(ip+1) - *ip);
                ip++;
                if(ii%16 == 15 && ii < gd->ncols-1)
                    printf("\n   %24s   ", "");
            }
            printf("\n");
        }
        break;
    default:
        printf("invalid code for quasi-regularity, %d\n", gd->quasi);
    }
    printf("   %24s : %#x\n", "GDS res/comp flag", gd->res_flags);
    printf("   %24s : %d\n", "GDS scan mode flag", gd->scan_mode);
/*
    printf("   %24s : %o\n", "GDS Octal scan mode flag", gd->scan_mode);
*/
    printf("   %24s : %d\n", "GDS no. of vert. coords", gd->nv);
    switch (gd->type) {
    case GRID_LL:                /* fall through */
    case GRID_RLL:                /* fall through */
    case GRID_SLL:                /* fall through */
    case GRID_SRLL:
    {
        gdes_ll *gg = &gd->grid.ll;
        printf("   %24s : %d\n", "GDS Ni", gg->ni);
        printf("   %24s : %d\n", "GDS Nj", gg->nj);
        printf("   %24s : %f\n", "GDS La1", gg->la1);
        printf("   %24s : %f\n", "GDS Lo1", gg->lo1);
        printf("   %24s : %f\n", "GDS La2", gg->la2);
        printf("   %24s : %f\n", "GDS Lo2", gg->lo2);
        printf("   %24s : %f\n", "GDS Di", gg->di);
        printf("   %24s : %f\n", "GDS Dj", gg->dj);
        if (gg->rot) {
            printf("   %24s : %f\n", "GDS Lat of S. pole of rotation",
                   gg->rot->lat);
            printf("   %24s : %f\n", "GDS Lon of S. pole of rotation",
                   gg->rot->lon);
            printf("   %24s : %f\n", "GDS Angle of rotation",
                   gg->rot->angle);
        }
        if (gg->strch) {
            printf("   %24s : %f\n", "GDS Lat of S. pole of stretching",
                   gg->strch->lat);
            printf("   %24s : %f\n", "GDS Lon of S. pole of stretching",
                   gg->strch->lon);
            printf("   %24s : %f\n", "GDS Stretching factor",
                   gg->strch->factor);
        }
    }
        break;
    case GRID_GAU:                /* fall through */
    case GRID_RGAU:                /* fall through */
    case GRID_SGAU:                /* fall through */
    case GRID_SRGAU:                /* fall through */
    {
        gdes_gau *gg = &gd->grid.gau;
        printf("   %24s : %d\n", "GDS Ni", gg->ni);
        printf("   %24s : %d\n", "GDS Nj", gg->nj);
        printf("   %24s : %f\n", "GDS La1", gg->la1);
        printf("   %24s : %f\n", "GDS Lo1", gg->lo1);
        printf("   %24s : %f\n", "GDS La2", gg->la2);
        printf("   %24s : %f\n", "GDS Lo2", gg->lo2);
        printf("   %24s : %f\n", "GDS Di", gg->di);
        printf("   %24s : %d\n", "GDS N", gg->n);
        if (gg->rot) {
            printf("   %24s : %f\n", "GDS Lat of S. pole of rotation",
                   gg->rot->lat);
            printf("   %24s : %f\n", "GDS Lon of S. pole of rotation",
                   gg->rot->lon);
            printf("   %24s : %f\n", "GDS Angle of rotation",
                   gg->rot->angle);
        }
        if (gg->strch) {
            printf("   %24s : %f\n", "GDS Lat of S. pole of stretching",
                   gg->strch->lat);
            printf("   %24s : %f\n", "GDS Lon of S. pole of stretching",
                   gg->strch->lon);
            printf("   %24s : %f\n", "GDS Stretching factor",
                   gg->strch->factor);
        }
    }
    break;
    case GRID_SPH:
    case GRID_RSPH:
    case GRID_SSPH:
    case GRID_SRSPH:
    {
        gdes_sph *gg = &gd->grid.sph;
        printf("   %24s : %d\n", "GDS j", gg->j);
        printf("   %24s : %d\n", "GDS k", gg->k);
        printf("   %24s : %d\n", "GDS m", gg->m);
        printf("   %24s : %d\n", "GDS type", gg->type);
        printf("   %24s : %d\n", "GDS mode", gg->mode);
        if (gg->rot) {
            printf("   %24s : %f\n", "GDS Lat of S. pole of rotation",
                   gg->rot->lat);
            printf("   %24s : %f\n", "GDS Lon of S. pole of rotation",
                   gg->rot->lon);
            printf("   %24s : %f\n", "GDS Angle of rotation",
                   gg->rot->angle);
        }
        if (gg->strch) {
            printf("   %24s : %f\n", "GDS Lat of S. pole of stretching",
                   gg->strch->lat);
            printf("   %24s : %f\n", "GDS Lon of S. pole of stretching",
                   gg->strch->lon);
            printf("   %24s : %f\n", "GDS Stretching factor",
                   gg->strch->factor);
        }
    }
    break;
    case GRID_MERCAT:
    {
        gdes_mercator *gg = &gd->grid.mercator;
        printf("   %24s : %d\n", "GDS Ni", gg->ni);
        printf("   %24s : %d\n", "GDS Nj", gg->nj);
        printf("   %24s : %f\n", "GDS La1", gg->la1);
        printf("   %24s : %f\n", "GDS Lo1", gg->lo1);
        printf("   %24s : %f\n", "GDS La2", gg->la2);
        printf("   %24s : %f\n", "GDS Lo2", gg->lo2);
        printf("   %24s : %f\n", "GDS Latin", gg->latin);
        printf("   %24s : %f\n", "GDS Di", gg->di);
        printf("   %24s : %f\n", "GDS Dj", gg->dj);
    }
    break;
    case GRID_GNOMON:                /* fall through */
    case GRID_POLARS:
    {
        gdes_polars *gg = &gd->grid.polars;
        printf("   %24s : %d\n", "GDS Nx", gg->nx);
        printf("   %24s : %d\n", "GDS Ny", gg->ny);
        printf("   %24s : %f\n", "GDS La1", gg->la1);
        printf("   %24s : %f\n", "GDS Lo1", gg->lo1);
        printf("   %24s : %f\n", "GDS Lov", gg->lov);
        printf("   %24s : %f\n", "GDS Dx", gg->dx);
        printf("   %24s : %f\n", "GDS Dy", gg->dy);
        printf("   %24s : %s\n", "GDS Pole in proj. plane",
               gg->pole == 0 ? "North" : "South");
    }
        break;
    case GRID_LAMBERT:
    {
        gdes_lambert *gg = &gd->grid.lambert;
        printf("   %24s : %d\n", "GDS Nx", gg->nx);
        printf("   %24s : %d\n", "GDS Ny", gg->ny);
        printf("   %24s : %f\n", "GDS La1", gg->la1);
        printf("   %24s : %f\n", "GDS Lo1", gg->lo1);
        printf("   %24s : %f\n", "GDS Lov", gg->lov);
        printf("   %24s : %f\n", "GDS Dx", gg->dx);
        printf("   %24s : %f\n", "GDS Dy", gg->dy);
        printf("   %24s : %s\n", "GDS Pole in proj. plane",
               gg->pole == 0 ? "North" : "South");
        printf("   %24s : %d\n", "GDS centers", gg->centers);
        printf("   %24s : %f\n", "GDS Latin1", gg->latin1);
        printf("   %24s : %f\n", "GDS Latin2", gg->latin2);
        printf("   %24s : %f\n", "GDS SpLat", gg->splat);
        printf("   %24s : %f\n", "GDS SpLon", gg->splon);
    }
        break;
    case GRID_SPACEV:
    {
        gdes_spacev *gg = &gd->grid.spacev;
        printf("   %24s : %d\n", "GDS Nx", gg->nx);
        printf("   %24s : %d\n", "GDS Ny", gg->ny);
        printf("   %24s : %f\n", "GDS Lap", gg->lap);
        printf("   %24s : %f\n", "GDS Lop", gg->lop);
        printf("   %24s : %f\n", "GDS dx", gg->dx);
        printf("   %24s : %f\n", "GDS dy", gg->dy);
        printf("   %24s : %f\n", "GDS Xp", gg->xp);
        printf("   %24s : %f\n", "GDS Yp", gg->yp);
        printf("   %24s : %f\n", "GDS Orientation", gg->orient);
        printf("   %24s : %f\n", "GDS Nr", gg->nr);
        printf("   %24s : %f\n", "GDS Xo", gg->xo);
        printf("   %24s : %f\n", "GDS Yo", gg->yo);
    }
    break;
    case GRID_ALBERS:
    case GRID_OLAMBERT:
    case GRID_UTM:
    case GRID_SIMPOL:
    case GRID_MILLER:
    default:
        break;
    }
}

void
print_netcdf_gdes_nav(gdes *gd){
/* 
This function is based on void print_gdes(gd) gdes *gd; and the code inside 
switch (gd->type) because  print_gdes is take as the basis for the 
categorisation of the different grids.
*/
/* This function prints the navigation variables to the standard output */
/* CALLED BY: print_cdl() in the file gribtocdl. */

printf("        // navigation variables all use nav dimension\n");
printf("\n");
printf("        char   nav_model(nav, nav_len) ;        // navigation parameterization\n");
printf("               nav_model:long_name = \"navigation model name\" ;\n");
printf("\n");
printf("        int    grid_type_code(nav) ;\n");
printf("               grid_type_code:long_name = \"GRIB-1 GDS data representation type\" ;\n");
printf("\n");
printf("        char   grid_type(nav, nav_len) ;\n");
printf("               grid_type:long_name = \"GRIB-1 grid type\" ;\n");
printf("\n");
printf("        char   grid_name(nav, nav_len) ;\n");
printf("               grid_name:long_name = \"grid name\" ;\n");
printf("\n");
printf("        int    grid_center(nav) ;\n");
printf("               grid_center:long_name = \"GRIB-1 originating center ID\" ;\n");
printf("\n");
printf("        int    grid_number(nav, ngrids) ;\n");
printf("               grid_number:long_name = \"GRIB-1 catalogued grid numbers\" ;\n");
printf("               grid_number:_FillValue = -9999 ;\n");
printf("\n");

switch (gd->type) {
    case GRID_LL:                /* fall through */
    case GRID_RLL:                /* fall through */
    case GRID_SLL:                /* fall through */
    case GRID_SRLL:
    {
gdes_ll *gg = &gd->grid.ll;
printf("        char   i_dim(nav, nav_len) ;\n");
printf("               i_dim:long_name = \"longitude dimension name\" ;\n");
printf("\n");
printf("        char   j_dim(nav, nav_len) ;\n");
printf("               j_dim:long_name = \"latitude dimension name\" ;\n");
printf("\n");
printf("        int    Ni(nav) ;\n");
printf("               Ni:long_name = \"number of points along a latitude circle\" ;\n");
printf("\n");
printf("        int    Nj(nav) ;\n");
printf("               Nj:long_name = \"number of points along a longitude circle\" ;\n");
printf("\n");
printf("        float  La1(nav) ;\n");
printf("               La1:long_name = \"latitude of first grid point\" ;\n");
printf("               La1:units = \"degrees_north\" ;\n");
printf("\n");
printf("        float  Lo1(nav) ;\n");
printf("               Lo1:long_name = \"longitude of first grid point\" ;\n");
printf("               Lo1:units = \"degrees_east\" ;\n");
printf("\n");
printf("        float  La2(nav) ;\n");
printf("               La2:long_name = \"latitude of last grid point\" ;\n");
printf("               La2:units = \"degrees_north\" ;\n");
printf("\n");
printf("        float  Lo2(nav) ;\n");
printf("               Lo2:long_name = \"longitude of last grid point\" ;\n");
printf("               Lo2:units = \"degrees_east\" ;\n");
printf("\n");
printf("        float  Di(nav) ;\n");
printf("               Di:long_name = \"longitudinal direction increment\" ;\n");
printf("               Di:units = \"degrees\" ;\n");
printf("\n");
printf("        float  Dj(nav) ;\n");
printf("               Dj:long_name = \"latitudinal direction increment\" ;\n");
printf("               Dj:units = \"degrees\" ;\n");

if( gg->rot ){

printf("\n");
printf("        float  RotAngle(nav) ;\n");
printf("               RotAngle:long_name = \"Angle of rotation\" ;\n");
printf("               RotAngle:units = \"degrees\" ;\n");
printf("\n");
printf("        float  RotLat(nav) ;\n");
printf("               RotLat:long_name = \"Lat of S. pole of rotation\" ;\n");
printf("               RotLat:units = \"degrees\" ;\n");
printf("\n");
printf("        float  RotLon(nav) ;\n");
printf("               RotLon:long_name = \"Lon of S. pole of rotation\" ;\n");
printf("               RotLon:units = \"degrees\" ;\n");
}

if( gg->strch ){ 

printf("\n");
printf("        float  StrchFactor(nav) ;\n");
printf("               StrchFactor:long_name = \"Stretching factor\" ;\n");
printf("               StrchFactor:units = \"    \" ;\n");
printf("\n");
printf("        float  StrchLat(nav) ;\n");
printf("               StrchLat:long_name = \"Lat of S. pole of stretching\" ;\n");
printf("               StrchLat:units = \"     \" ;\n");
printf("\n");
printf("        float  StrchLon(nav) ;\n");
printf("               StrchLon:long_name = \"Lon of S. pole of stretching\" ;\n");
printf("               StrchLon:units = \"      \" ;\n");
}
    }
        break;
    case GRID_GAU:                /* fall through */
    case GRID_RGAU:                /* fall through */
    case GRID_SGAU:                /* fall through */
    case GRID_SRGAU:                /* fall through */
    {
gdes_gau *gg = &gd->grid.gau;
printf("        char   i_dim(nav, nav_len) ;\n");
printf("               i_dim:long_name = \"longitude dimension name\" ;\n");
printf("\n");
printf("        char   j_dim(nav, nav_len) ;\n");
printf("               j_dim:long_name = \"latitude dimension name\" ;\n");
printf("\n");
printf("        int    Ni(nav) ;\n");
printf("               Ni:long_name = \"number of points along a latitude circle\" ;\n");
printf("\n");
printf("        int    Nj(nav) ;\n");
printf("               Nj:long_name =        \"number of points along a longitude circle\" ;\n");
printf("\n");
printf("        float  La1(nav) ;\n");
printf("               La1:long_name = \"latitude of first grid point\" ;\n");
printf("               La1:units = \"degrees_north\" ;\n");
printf("\n");
printf("        float  Lo1(nav) ;\n");
printf("               Lo1:long_name = \"longitude of first grid point\" ;\n");
printf("               Lo1:units = \"degrees_east\" ;\n");
printf("\n");
printf("        float  La2(nav) ;\n");
printf("               La2:long_name = \"latitude of last grid point\" ;\n");
printf("               La2:units = \"degrees_north\" ;\n");
printf("\n");
printf("        float  Lo2(nav) ;\n");
printf("               Lo2:long_name = \"longitude of last grid point\" ;\n");
printf("               Lo2:units = \"degrees_east\" ;\n");
printf("\n");
printf("        float  Di(nav) ;\n");
printf("               Di:long_name = \"longitudinal direction increment\" ;\n");
printf("               Di:units = \"degrees\" ;\n");
printf("\n");
printf("        float  N(nav) ;\n");
printf("               N:long_name = \"# of parallels between a pole and equator\" ;\n");
printf("               N:units = \"\" ;\n");
if( gg->rot ){
printf("\n");
printf("        float  RotAngle(nav) ;\n");
printf("               RotAngle:long_name = \"Angle of rotation\" ;\n");
printf("               RotAngle:units = \"degrees\" ;\n");
printf("\n");
printf("        float  RotLat(nav) ;\n");
printf("               RotLat:long_name = \"Lat of S. pole of rotation\" ;\n");
printf("               RotLat:units = \"degrees\" ;\n");
printf("\n");
printf("        float  RotLon(nav) ;\n");
printf("               RotLon:long_name = \"Lon of S. pole of rotation\" ;\n");
printf("               RotLon:units = \"degrees\" ;\n");
}

if( gg->strch ){ 

printf("\n");
printf("        float  StrchFactor(nav) ;\n");
printf("               StrchFactor:long_name = \"Stretching factor\" ;\n");
printf("               StrchFactor:units = \"    \" ;\n");
printf("\n");
printf("        float  StrchLat(nav) ;\n");
printf("               StrchLat:long_name = \"Lat of S. pole of stretching\" ;\n");
printf("               StrchLat:units = \"     \" ;\n");
printf("\n");
printf("        float  StrchLon(nav) ;\n");
printf("               StrchLon:long_name = \"Lon of S. pole of stretching\" ;\n");
printf("               StrchLon:units = \"     \" ;\n");
}        
    }
    break;
    case GRID_SPH:
    case GRID_RSPH:
    case GRID_SSPH:
    case GRID_SRSPH:
    {
        gdes_sph *gg = &gd->grid.sph;
        printf(" %s = %d\n", "j", gg->j);
        printf(" %s = %d\n", "k", gg->k);
        printf(" %s = %d\n", "m", gg->m);
        printf(" %s = %d\n", "type", gg->type);
        printf(" %s = %d\n", "mode", gg->mode);

        if( gg->rot ){
printf("\n");
printf("        float  RotAngle(nav) ;\n");
printf("               RotAngle:long_name = \"Angle of rotation\" ;\n");
printf("               RotAngle:units = \"degrees\" ;\n");
printf("\n");
printf("        float  RotLat(nav) ;\n");
printf("               RotLat:long_name = \"Lat of S. pole of rotation\" ;\n");
printf("               RotLat:units = \"degrees\" ;\n");
printf("\n");
printf("        float  RotLon(nav) ;\n");
printf("               RotLon:long_name = \"Lon of S. pole of rotation\" ;\n");
printf("               RotLon:units = \"degrees\" ;\n");
}

if( gg->strch ){ 

printf("\n");
printf("        float  StrchFactor(nav) ;\n");
printf("               StrchFactor:long_name = \"Stretching factor\" ;\n");
printf("               StrchFactor:units = \"    \" ;\n");
printf("\n");
printf("        float  StrchLat(nav) ;\n");
printf("               StrchLat:long_name = \"Lat of S. pole of stretching\" ;\n");
printf("               StrchLat:units = \"     \" ;\n");
printf("\n");
printf("        float  StrchLon(nav) ;\n");
printf("               StrchLon:long_name = \"Lon of S. pole of stretching\" ;\n");
printf("               StrchLon:units = \"     \" ;\n");
}                
    }
    break;
    case GRID_MERCAT:
    {
printf("        char   i_dim(nav, nav_len) ;\n");
printf("               i_dim:long_name = \"longitude dimension name\" ;\n");
printf("\n");
printf("        char   j_dim(nav, nav_len) ;\n");
printf("               j_dim:long_name = \"latitude dimension name\" ;\n");
printf("\n");
printf("        int    Ni(nav) ;\n");
printf("               Ni:long_name = \"number of points along a latitude circle\" ;\n");
printf("\n");
printf("        int    Nj(nav) ;\n");
printf("               Nj:long_name =        \"number of points along a longitude circle\" ;\n");
printf("\n");
printf("        float  La1(nav) ;\n");
printf("               La1:long_name = \"latitude of first grid point\" ;\n");
printf("               La1:units = \"degrees_north\" ;\n");
printf("\n");
printf("        float  Lo1(nav) ;\n");
printf("               Lo1:long_name = \"longitude of first grid point\" ;\n");
printf("               Lo1:units = \"degrees_east\" ;\n");
printf("\n");
printf("        float  La2(nav) ;\n");
printf("               La2:long_name = \"latitude of last grid point\" ;\n");
printf("               La2:units = \"degrees_north\" ;\n");
printf("\n");
printf("        float  Lo2(nav) ;\n");
printf("               Lo2:long_name = \"longitude of last grid point\" ;\n");
printf("               Lo2:units = \"degrees_east\" ;\n");
printf("\n");
printf("        float  Di(nav) ;\n");
printf("               Di:long_name = \"longitudinal direction increment\" ;\n");
printf("               Di:units = \"degrees\" ;\n");
printf("\n");
printf("        float  Latin(nav) ;\n");
printf("               Latin:long_name = \"XXX\" ;\n");
printf("               Latin:units = \" \" ;\n");        
    }
    break;
    case GRID_GNOMON:                /* fall through */
    case GRID_POLARS:
    {
printf("        char   x_dim(nav, nav_len) ;\n");
printf("               x_dim:long_name = \"x dimension name\" ;\n");
printf("\n");
printf("        char   y_dim(nav, nav_len) ;\n");
printf("               y_dim:long_name = \"y dimension name\" ;\n");
printf("\n");
printf("        long   Nx(nav) ;\n");
printf("               Nx:long_name = \"number of points along x-axis\" ;\n");
printf("\n");
printf("        long   Ny(nav) ;\n");
printf("               Ny:long_name =  \"number of points along y-axis\" ;\n");
printf("\n");
printf("        float  La1(nav) ;\n");
printf("               La1:long_name = \"latitude of first grid point\" ;\n");
printf("               La1:units = \"degrees_north\" ;\n");
printf("\n");
printf("        float  Lo1(nav) ;\n");
printf("               Lo1:long_name = \"longitude of first grid point\" ;\n");
printf("               Lo1:units = \"degrees_east\" ;\n");
printf("\n");
printf("        float  Lov(nav) ;\n");
printf("               Lov:long_name = \"orientation of the grid\" ;\n");
printf("               Lov:units = \"degrees_east\" ;\n");
printf("\n");
printf("        float  Dx(nav) ;\n");
printf("               Dx:long_name = \"x-direction grid length\" ;\n");
printf("               Dx:units = \"m\" ;\n");
printf("\n");
printf("        float  Dy(nav) ;\n");
printf("               Dy:long_name = \"y-direction grid length\" ;\n");
printf("               Dy:units = \"m\" ;\n");
printf("\n");
printf("        byte   ProjFlag(nav) ;\n");
printf("               ProjFlag:long_name = \"projection center flag\" ;\n");
printf("\n");
    }
        break;
    case GRID_LAMBERT:
    {
printf("        char   x_dim(nav, nav_len) ;\n");
printf("               x_dim:long_name = \"x dimension name\" ;\n");
printf("\n");
printf("        char   y_dim(nav, nav_len) ;\n");
printf("               y_dim:long_name = \"y dimension name\" ;\n");
printf("\n");
printf("        long   Nx(nav) ;\n");
printf("               Nx:long_name = \"number of points along x-axis\" ;\n");
printf("\n");
printf("        long   Ny(nav) ;\n");
printf("               Ny:long_name =  \"number of points along y-axis\" ;\n");
printf("\n");
printf("        float  La1(nav) ;\n");
printf("               La1:long_name = \"latitude of first grid point\" ;\n");
printf("               La1:units = \"degrees_north\" ;\n");
printf("\n");
printf("        float  Lo1(nav) ;\n");
printf("               Lo1:long_name = \"longitude of first grid point\" ;\n");
printf("               Lo1:units = \"degrees_east\" ;\n");
printf("\n");
printf("        float  Lov(nav) ;\n");
printf("               Lov:long_name = \"orientation of the grid\" ;\n");
printf("               Lov:units = \"degrees_east\" ;\n");
printf("\n");
printf("        float  Dx(nav) ;\n");
printf("               Dx:long_name = \"x-direction grid length\" ;\n");
printf("               Dx:units = \"m\" ;\n");
printf("\n");
printf("        float  Dy(nav) ;\n");
printf("               Dy:long_name = \"y-direction grid length\" ;\n");
printf("               Dy:units = \"m\" ;\n");
printf("\n");
printf("        byte   ProjFlag(nav) ;\n");
printf("               ProjFlag:long_name = \"projection center flag\" ;\n");
printf("\n");
printf("        float  Latin1(nav) ;\n");
printf("               Latin1:long_name = \"first intersecting latitude\" ;\n");
printf("               Latin1:units = \"degrees_north\" ;\n");
printf("\n");
printf("        float  Latin2(nav) ;\n");
printf("               Latin2:long_name = \"second intersecting latitude\" ;\n");
printf("               Latin2:units = \"degrees_north\" ;\n");
printf("\n");
printf("        float  SpLat(nav) ;\n");
printf("               SpLat:long_name = \"latitude of the southern pole\" ;\n");
printf("               SpLat:units = \"degrees_north\" ;\n");
printf("\n");
printf("        float  SpLon(nav) ;\n");
printf("               SpLon:long_name = \"longitude of the southern pole\" ;\n");
printf("               SpLon:units = \"degrees_east\" ;\n");

    }
        break;
    case GRID_SPACEV:
    {
        gdes_spacev *gg = &gd->grid.spacev;
        printf(" %s = %d\n", "Nx", gg->nx);
        printf(" %s = %d\n", "Ny", gg->ny);
        printf(" %s = %f\n", "Lap", gg->lap);
        printf(" %s = %f\n", "Lop", gg->lop);
        printf(" %s = %f\n", "dx", gg->dx);
        printf(" %s = %f\n", "dy", gg->dy);
        printf(" %s = %f\n", "Xp", gg->xp);
        printf(" %s = %f\n", "Yp", gg->yp);
        printf(" %s = %f\n", "Orientation", gg->orient);
        printf(" %s = %f\n", "Nr", gg->nr);
        printf(" %s = %f\n", "Xo", gg->xo);
        printf(" %s = %f\n", "Yo", gg->yo);
    }
    break;
    case GRID_ALBERS:
    case GRID_OLAMBERT:
    case GRID_UTM:
    case GRID_SIMPOL:
    case GRID_MILLER:
    default:
        break;
    }
    printf("\n");
    printf("        byte   ResCompFlag(nav) ;\n");
    printf("               ResCompFlag:long_name = \"resolution and component flags\" ;\n");
    printf("\n");
    printf("        // end of navigation variables\n\n");

    return;
} /* end of print_netcdf_gdes_nav */

void
print_netcdf_gdes_type(gdes *gd)
/* 
This function is based on void print_gdes(gd) gdes *gd; and the code inside 
switch (gd->type)
*/
/* "GDS " is replaced by "" */
/* This function prints the content of the navigation variables to the standard output */
/* CALLED BY: print_cdl() in the file gribtocdl. */
{
/*
   printf("   %24s : %d\n", "GDS scan mode flag", gd->scan_mode);
   printf("   %24s : %d\n", "GDS no. of vert. coords", gd->nv);
*/

switch (gd->type) {
    case GRID_LL:                /* fall through */
    case GRID_RLL:                /* fall through */
    case GRID_SLL:                /* fall through */
    case GRID_SRLL:
    {
        gdes_ll *gg = &gd->grid.ll;
        /* updated 18.12.2001 by ab 
           Adjust grid so that it is always South-North [-90:+90]
           and West-East [-180:180] or [0:360] */
        /* adjust numbers so that they are increasing */
        /* Assume longitude is allways west-east (sic) */
        while (gg->lo2 < gg->lo1) gg->lo2 += 360;
        /* Just swap latitudes if wrong direction */

        if (gg->la1 > gg->la2 ) {
          float tlat=gg->la1;
          gg->la1=gg->la2;
          gg->la2=tlat;
        }
        printf(" %s = \"%s\" ;\n", "i_dim", "lon");
        printf(" %s = \"%s\" ;\n", "j_dim", "lat");
        printf(" %s = %d ;\n", "Ni", gg->ni);
        printf(" %s = %d ;\n", "Nj", gg->nj);
        printf(" %s = %f ;\n", "La1", gg->la1);
        printf(" %s = %f ;\n", "Lo1", gg->lo1);
        printf(" %s = %f ;\n", "La2", gg->la2);
        printf(" %s = %f ;\n", "Lo2", gg->lo2);
        printf(" %s = %f ;\n", "Di", gg->di);
        printf(" %s = %f ;\n", "Dj", gg->dj);
        if (gg->rot) {
            printf(" %s = %f ;\n", "RotLat",gg->rot->lat);
            printf(" %s = %f ;\n", "RotLon",gg->rot->lon);
            printf(" %s = %f ;\n", "RotAngle",gg->rot->angle);
        }
        if (gg->strch) {
            printf(" %s = %f ;\n", "StrchLat", gg->strch->lat);
            printf(" %s = %f ;\n", "StrchLon", gg->strch->lon);
            printf(" %s = %f ;\n", "StrchFactor", gg->strch->factor);
        }
    }
        break;
    case GRID_GAU:                /* fall through */
    case GRID_RGAU:                /* fall through */
    case GRID_SGAU:                /* fall through */
    case GRID_SRGAU:                /* fall through */
    {
        gdes_gau *gg = &gd->grid.gau;

        if (gg->la1 > gg->la2 ) {
          float tlat=gg->la1;
          gg->la1=gg->la2;
          gg->la2=tlat;
        }
        printf(" %s = %d ;\n", "Ni", gg->ni);
        printf(" %s = %d ;\n", "Nj", gg->nj);
        printf(" %s = %f ;\n", "La1", gg->la1);
        printf(" %s = %f ;\n", "Lo1", gg->lo1);
        printf(" %s = %f ;\n", "La2", gg->la2);
        printf(" %s = %f ;\n", "Lo2", gg->lo2);
        printf(" %s = %f ;\n", "Di", gg->di);
        printf(" %s = %d ;\n", "N", gg->n);
        if (gg->rot) {
            printf(" %s = %f ;\n", "RotLat",gg->rot->lat);
            printf(" %s = %f ;\n", "RotLon",gg->rot->lon);
            printf(" %s = %f ;\n", "RotAngle",gg->rot->angle);
        }
        if (gg->strch) {
            printf(" %s = %f ;\n", "StrchLat", gg->strch->lat);
            printf(" %s = %f ;\n", "StrchLon", gg->strch->lon);
            printf(" %s = %f ;\n", "StrchFactor", gg->strch->factor);
        }
    }
    break;
    case GRID_SPH:
    case GRID_RSPH:
    case GRID_SSPH:
    case GRID_SRSPH:
    {
        gdes_sph *gg = &gd->grid.sph;
        printf(" %s = %d ;\n", "j", gg->j);
        printf(" %s = %d ;\n", "k", gg->k);
        printf(" %s = %d ;\n", "m", gg->m);
        printf(" %s = %d ;\n", "type", gg->type);
        printf(" %s = %d ;\n", "mode", gg->mode);
        if (gg->rot) {
            printf(" %s = %f ;\n", "RotLat",gg->rot->lat);
            printf(" %s = %f ;\n", "RotLon",gg->rot->lon);
            printf(" %s = %f ;\n", "RotAngle",gg->rot->angle);
        }
        if (gg->strch) {
            printf(" %s = %f ;\n", "StrchLat", gg->strch->lat);
            printf(" %s = %f ;\n", "StrchLon", gg->strch->lon);
            printf(" %s = %f ;\n", "StrchFactor", gg->strch->factor);
        }
    }
    break;
    case GRID_MERCAT:
    {
        gdes_mercator *gg = &gd->grid.mercator;
        printf(" %s = %d ;\n", "Ni", gg->ni);
        printf(" %s = %d ;\n", "Nj", gg->nj);
        printf(" %s = %f ;\n", "La1", gg->la1);
        printf(" %s = %f ;\n", "Lo1", gg->lo1);
        printf(" %s = %f ;\n", "La2", gg->la2);
        printf(" %s = %f ;\n", "Lo2", gg->lo2);
        printf(" %s = %f ;\n", "Latin", gg->latin);
        printf(" %s = %f ;\n", "Di", gg->di);
        printf(" %s = %f ;\n", "Dj", gg->dj);
    }
    break;
    case GRID_GNOMON:                /* fall through */
    case GRID_POLARS:
    {
        gdes_polars *gg = &gd->grid.polars;
        printf(" x_dim = \"x\" ;\n");
        printf(" y_dim = \"y\" ;\n");
        printf(" %s = %d ;\n", "Nx", gg->nx);
        printf(" %s = %d ;\n", "Ny", gg->ny);
        printf(" %s = %f ;\n", "La1", gg->la1);
        printf(" %s = %f ;\n", "Lo1", gg->lo1);
        printf(" %s = %f ;\n", "Lov", gg->lov);
        printf(" %s = %f ;\n", "Dx", gg->dx);
        printf(" %s = %f ;\n", "Dy", gg->dy);
        printf(" // %s = %s ;\n", "Pole in proj. plane",
               gg->pole == 0 ? "North" : "South");
        printf(" %s = %d ;\n", "ProjFlag", gg->pole);
    }
        break;
    case GRID_LAMBERT:
    {
        gdes_lambert *gg = &gd->grid.lambert;
        printf(" x_dim = \"x\" ;\n");
        printf(" y_dim = \"y\" ;\n");
        printf(" %s = %d ;\n", "Nx", gg->nx);
        printf(" %s = %d ;\n", "Ny", gg->ny);
        printf(" %s = %f ;\n", "La1", gg->la1);
        printf(" %s = %f ;\n", "Lo1", gg->lo1);
        printf(" %s = %f ;\n", "Lov", gg->lov);
        printf(" %s = %f ;\n", "Dx", gg->dx);
        printf(" %s = %f ;\n", "Dy", gg->dy);
        printf(" %s = %f ;\n", "Latin1", gg->latin1);
        printf(" %s = %f ;\n", "Latin2", gg->latin2);
        printf(" %s = %f ;\n", "SpLat", gg->splat);
        printf(" %s = %f ;\n", "SpLon", gg->splon);
        printf(" // %s = %s ;\n", "Pole in proj. plane",
               gg->pole == 0 ? "North" : "South");
        printf(" %s = %d ;\n", "ProjFlag", gg->pole);
    }
        break;
    case GRID_SPACEV:
    {
        gdes_spacev *gg = &gd->grid.spacev;
        printf(" %s = %d ;\n", "Nx", gg->nx);
        printf(" %s = %d ;\n", "Ny", gg->ny);
        printf(" %s = %f ;\n", "Lap", gg->lap);
        printf(" %s = %f ;\n", "Lop", gg->lop);
        printf(" %s = %f ;\n", "dx", gg->dx);
        printf(" %s = %f ;\n", "dy", gg->dy);
        printf(" %s = %f ;\n", "Xp", gg->xp);
        printf(" %s = %f ;\n", "Yp", gg->yp);
        printf(" %s = %f ;\n", "Orientation", gg->orient);
        printf(" %s = %f ;\n", "Nr", gg->nr);
        printf(" %s = %f ;\n", "Xo", gg->xo);
        printf(" %s = %f ;\n", "Yo", gg->yo);
    }
    break;
    case GRID_ALBERS:
    case GRID_OLAMBERT:
    case GRID_UTM:
    case GRID_SIMPOL:
    case GRID_MILLER:
    default:
        break;
    }
} /* end of print_netcdf_gdes_type */

gdes* copy_gdes( gdes* gd, gdes *ret ){
/* This function take a deep copy from gd to ret*/ 

    int nerrs = 0;

    ret->type  = gd->type;
    ret->ncols = gd->ncols;
    ret->nrows = gd->nrows;
    ret->npts  = gd->npts;
    ret->res_flags=gd->res_flags;
    ret->scan_mode=gd->scan_mode;
    ret->maxlc = gd->maxlc;

    ret->quasi = gd->quasi;        /* ordinary rectangular grid is default */
    ret->nv = gd->nv;
    ret->vc = 0;
    ret->lc = 0;
    ret->keep = gd->keep;

    if(gd->nv != 0 ){
        
        int i;
        ret->vc = (float *) emalloc(ret->nv * sizeof(float));
        /* unpack the vertical coords into floats */
        for (i = 0; i < ret->nv; i++) {
            ret->vc[i] = gd->vc[i];
        }
    }

    switch (gd->type) {
    case GRID_LL:                /* fall through */
    case GRID_RLL:                /* fall through */
    case GRID_SLL:                /* fall through */
    case GRID_SRLL:
    {
        copy_gdes_ll(&(gd->grid.ll),&(ret->grid.ll));
    }
        break;
    case GRID_GAU:                /* fall through */
    case GRID_RGAU:                /* fall through */
    case GRID_SGAU:                /* fall through */
    case GRID_SRGAU:                /* fall through */
    {
        copy_gdes_gau(&(gd->grid.gau),&(ret->grid.gau));  
    }
    break;
    case GRID_SPH:
    case GRID_RSPH:
    case GRID_SSPH:
    case GRID_SRSPH:
    {
        copy_gdes_sph(&(gd->grid.sph),&(ret->grid.sph));
      
    }
    break;
    case GRID_MERCAT:
    {
       copy_gdes_mercat(&(gd->grid.mercator),&(ret->grid.mercator));
    }
    break;
    case GRID_GNOMON:                /* fall through */
    case GRID_POLARS:
    {
        copy_gdes_polars( &(gd->grid.polars), &(ret->grid.polars) );
    }
        break;
    case GRID_LAMBERT:
    {
        copy_gdes_lambert( &(gd->grid.lambert), &(ret->grid.lambert) );        
    }
        break;
    case GRID_SPACEV:
    {
        /*ret.grid.spacev=copy_gdes_spacev(gd);*/
    }
    break;
    case GRID_ALBERS:
    case GRID_OLAMBERT:
    case GRID_UTM:
    case GRID_SIMPOL:
    case GRID_MILLER:   
    default:
        nerrs++;
    }

    if (nerrs) {
        free_gdes(ret);
        ret = 0;
        return ret;
    }

    if(ret->quasi == QUASI_ROWS) {
        ret->ncols = 1;
        ret->lc = (int *)emalloc((1 + ret->nrows) * sizeof(int));
        {        /* unpack list of row indexes */
            int i;
            for (i = 0; i < ret->nrows+1; i++) {
                ret->lc[i] = gd->lc[i];
            }
         
        }
    } else if(ret->quasi == QUASI_COLS) {
        ret->nrows = 1;
        ret->lc = (int *)emalloc((1 + ret->ncols) * sizeof(int));
        {        /* unpack list of col indexes */
            int i;
            for (i = 0; i < ret->ncols+1; i++) {
                ret->lc[i] = gd->lc[i];
            }
            
        }
    }
    if (nerrs) {
        free_gdes(ret);
        ret = 0;
    }
    return ret;
} /* end of  copy_gdes */

/* extern int testNoserc; */
#define CMP( A )  if(0) printf("   %2s \n", A)

static int cmp_int_start( int* v1, int n1, int* v2, int n2 ){
/* the first half of the function cmp_int_set and  cmp_int_set_vector(not implemented) */
/* this function must be called before cmp_int_set_base */

     /* both pointers equal zero */
    if( v1 == NULL && v2 == NULL ){CMP("CMP zero pointers"); return 1;}

    /* one pointer is equal to zero, and the other is not and 
       the values for the other are greater than zero */
    if( (v1 == NULL && v2 != NULL && n2 > 0) ||  
		(v2 == NULL && v1 != NULL  && n1 > 0)){
         CMP("CMP v NULL"); 
         return 0;}

    /* both pointers are different from zero */
    if( v1 != NULL && v2 != NULL ){
        /* cmp number of items in array */
        if( n1 != n2 ){ CMP("CMP n 0"); return 0;}
        if( n1 == n2 == 0 ){ CMP("CMP n 1"); return 1;}
    }
    
    return 1;
}

static int cmp_float_start( float* v1, int n1, float* v2, int n2 ){
/* the first half of the function cmp_float_set and  cmp_float_set_vector(not implemented) */
/* this function must be called before cmp_float_set_base */

     /* both pointers equal zero */
    if( v1 == NULL && v2 == NULL ){CMP("CMP zero pointers"); return 1;}

    /* one pointer is equal to zero, and the other is not and 
       the values for the other are greater than zero */
    if( (v1 == NULL && v2 != NULL && n2 > 0) ||  (v2 == NULL && v1 != NULL  && n1 > 0)){
         CMP("CMP v NULL"); 
         return 0;}

    /* both pointers are different from zero */
    if( v1 != NULL && v2 != NULL ){
        /* cmp number of items in array */
        if( n1 != n2 ){ CMP("CMP n 0"); return 0;}
        if( n1 == n2 == 0 ){ CMP("CMP n 1"); return 1;}
    }
    
    return 1;
}

static int cmp_int_set_base( int* v1, int n1, int* v2, int n2 ){
/* the second half of the function cmp_int_set */
/* this function must be called after cmp_int_start*/
/* This function compares two sets stored as arrays */
        
    if( n1 == n2 && n1 > 0 &&  n2 > 0 && v1 != NULL && v2 != NULL ){
        int b;
        int i;
        int j;
        for( i=0; i<n1; i++ ){/* the arrays are sets and not vectors */
            b=0;
            for( j=0; j<n2; j++ ){
                if( v1[i] == v2[j] ){ b=1; break;}
            }
            if( b==0 ){
                CMP("CMP v");
                return 0;
            }
        }
    }

      return 1;
}

static int cmp_float_set_base( float* v1, int n1, float* v2, int n2 ){
/* the second half of the function cmp_float_set */
/* this function must be called after cmp_float_start*/
/* This function compares two sets stored as arrays */
        
    if( n1 == n2 && n1 > 0 &&  n2 > 0 && v1 != NULL && v2 != NULL ){
        int b;
        int i;
        int j;
        for( i=0; i<n1; i++ ){/* the arrays are sets and not vectors */
            b=0;
            for( j=0; j<n2; j++ ){
                if( v1[i] == v2[j] ){ b=1; break;}
            }
            if( b==0 ){
                CMP("CMP v");
                return 0;
            }
        }
    }

    return 1;
}

static int cmp_int_set( int* v1, int n1, int* v2, int n2 ){
/* This function compares two sets stored as arrays */
    if(cmp_int_start( v1, n1, v2, n2 ))
        return 0;

    return cmp_int_set_base( v1, n1, v2, n2 );
}

cdl_gdes_data* get_gdes_data(gdes* gd){
/* This function is based on void print_gdes(gd) gdes *gd; and the code inside switch (gd->type)*/
/* This function takes out and returns data from the union gd.grid */
/* CALLED BY: print_cdl() in the file gribtocdl. */

 cdl_gdes_data* lgdat=(cdl_gdes_data*)emalloc(sizeof(cdl_gdes_data));

 lgdat->ni=-1;
 lgdat->nj=-1;  
 lgdat->la1=-1.;
 lgdat->lo1=-1.;
 lgdat->la2=-1.;
 lgdat->lo2=-1.;
 lgdat->di=-1.;
 lgdat->dj=-1.;
 /*printf(" %s \n", gds_typename(gd->type));*/

switch (gd->type) {
    case GRID_LL:                /* fall through */
    case GRID_RLL:                /* fall through */
    case GRID_SLL:                /* fall through */
    case GRID_SRLL:
    {        
        gdes_ll *gg = &gd->grid.ll;
        lgdat->ni=gg->ni;  /*printf("   %24s : %d\n", "GDS Ni", lgdat->ni);*/
        lgdat->nj=gg->nj;  /*printf("   %24s : %d\n", "GDS Nj", lgdat->nj);*/
        /*Nj=gg->nj;
        lgdat->nj=Nj;*/
        lgdat->la1=gg->la1; /* printf("   %24s : %f\n", "GDS La1", gg->la1);*/
        lgdat->lo1=gg->lo1; /* printf("   %24s : %f\n", "GDS Lo1", gg->lo1);*/
        lgdat->la2=gg->la2; /* printf("   %24s : %f\n", "GDS La2", gg->la2);*/
        lgdat->lo2=gg->lo2; /*printf("   %24s : %f\n", "GDS Lo2", gg->lo2);*/
        lgdat->di=gg->di;  /*printf("   %24s : %f\n", "GDS Di", gg->di);*/
        lgdat->dj=gg->dj;  /*printf("   %24s : %f\n", "GDS Dj", gg->dj); */
        gg=0;
        return lgdat;
    }
    break;
    case GRID_GAU:                /* fall through */
    case GRID_RGAU:                /* fall through */
    case GRID_SGAU:                /* fall through */
    case GRID_SRGAU:                /* fall through */
    {
        gdes_gau *gg = &gd->grid.gau;
        lgdat->ni=gg->ni;/*printf("   %24s : %d\n", "GDS Ni", gg->ni);*/
        lgdat->nj=gg->nj;/*printf("   %24s : %d\n", "GDS Nj", gg->nj);*/
        lgdat->la1=gg->la1; /* printf("   %24s : %f\n", "GDS La1", gg->la1);*/
        lgdat->lo1=gg->lo1; /* printf("   %24s : %f\n", "GDS Lo1", gg->lo1);*/
        lgdat->la2=gg->la2; /* printf("   %24s : %f\n", "GDS La2", gg->la2);*/
        lgdat->lo2=gg->lo2; /*printf("   %24s : %f\n", "GDS Lo2", gg->lo2);*/
        lgdat->di=gg->di;  /*printf("   %24s : %f\n", "GDS Di", gg->di);*/
                           /*printf("   %24s : %d\n", "GDS n", gg->n);*/ 
    }
    break;
    case GRID_SPH:
    case GRID_RSPH:
    case GRID_SSPH:
    case GRID_SRSPH:
    {
        /* gdes_sph *gg = &gd->grid.sph; */
        lgdat->ni=-1;
        /*lgdat->nj=-1;*/
        /* printf("   %24s : %d\n", "GDS j", gg->j);
        printf("   %24s : %d\n", "GDS k", gg->k);
        printf("   %24s : %d\n", "GDS m", gg->m);
        printf("   %24s : %d\n", "GDS type", gg->type);
        printf("   %24s : %d\n", "GDS mode", gg->mode); */

    }
    break;
    case GRID_MERCAT:
    {
        gdes_mercator *gg = &gd->grid.mercator;    
        lgdat->ni=gg->ni;/*printf("   %24s : %d\n", "GDS Ni", gg->ni);*/
        /*lgdat->nj=gg->nj;*//*printf("   %24s : %d\n", "GDS Nj", gg->nj);*/                             
        lgdat->la1=gg->la1; /* printf("   %24s : %f\n", "GDS La1", gg->la1);*/
        lgdat->lo1=gg->lo1; /* printf("   %24s : %f\n", "GDS Lo1", gg->lo1);*/
        lgdat->la2=gg->la2; /* printf("   %24s : %f\n", "GDS La2", gg->la2);*/
        lgdat->lo2=gg->lo2; /*printf("   %24s : %f\n", "GDS Lo2", gg->lo2);*/
        lgdat->di=gg->di;  /*printf("   %24s : %f\n", "GDS Di", gg->di);*/
        /*printf("   %24s : %f\n", "GDS Latin", gg->latin);*/
    }
    break;
    case GRID_GNOMON:                /* fall through */
    case GRID_POLARS:
    {
        gdes_polars *gg = &gd->grid.polars;
        lgdat->ni=-1;
        /*lgdat->nj=-1;*/
        lgdat->la1=gg->la1; /* printf("   %24s : %f\n", "GDS La1", gg->la1);*/
        lgdat->lo1=gg->lo1; /* printf("   %24s : %f\n", "GDS Lo1", gg->lo1);*/
     /* printf("   %24s : %d\n", "GDS Nx", gg->nx);
        printf("   %24s : %d\n", "GDS Ny", gg->ny);
        printf("   %24s : %f\n", "GDS La1", gg->la1);
        printf("   %24s : %f\n", "GDS Lo1", gg->lo1);
        printf("   %24s : %f\n", "GDS Lov", gg->lov);
        printf("   %24s : %f\n", "GDS Dx", gg->dx);
        printf("   %24s : %f\n", "GDS Dy", gg->dy);
        printf("   %24s : %s\n", "GDS Pole in proj. plane",
               gg->pole == 0 ? "North" : "South"); */
    }
        break;
    case GRID_LAMBERT:
    {
        gdes_lambert *gg = &gd->grid.lambert;
        lgdat->ni=-1;
        /*lgdat->nj=-1;*/
        lgdat->la1=gg->la1; /* printf("   %24s : %f\n", "GDS La1", gg->la1);*/
        lgdat->lo1=gg->lo1; /* printf("   %24s : %f\n", "GDS Lo1", gg->lo1);*/
        /*printf("   %24s : %d\n", "GDS Nx", gg->nx);
        printf("   %24s : %d\n", "GDS Ny", gg->ny);
        printf("   %24s : %f\n", "GDS La1", gg->la1);
        printf("   %24s : %f\n", "GDS Lo1", gg->lo1);
        printf("   %24s : %f\n", "GDS Lov", gg->lov);
        printf("   %24s : %f\n", "GDS Dx", gg->dx);
        printf("   %24s : %f\n", "GDS Dy", gg->dy);
        printf("   %24s : %s\n", "GDS Pole in proj. plane",
               gg->pole == 0 ? "North" : "South");
        printf("   %24s : %d\n", "GDS centers", gg->centers);
        printf("   %24s : %f\n", "GDS Latin1", gg->latin1);
        printf("   %24s : %f\n", "GDS Latin2", gg->latin2);
        printf("   %24s : %f\n", "GDS SpLat", gg->splat);
        printf("   %24s : %f\n", "GDS SpLon", gg->splon);*/
    }
        break;
    case GRID_SPACEV:
    {
      /* gdes_spacev *gg = &gd->grid.spacev; */
        lgdat->ni=-1;
        /*lgdat->nj=-1;*/
        /* printf("   %24s : %d\n", "GDS Nx", gg->nx);
        printf("   %24s : %d\n", "GDS Ny", gg->ny);
        printf("   %24s : %f\n", "GDS Lap", gg->lap);
        printf("   %24s : %f\n", "GDS Lop", gg->lop);
        printf("   %24s : %f\n", "GDS dx", gg->dx);
        printf("   %24s : %f\n", "GDS dy", gg->dy);
        printf("   %24s : %f\n", "GDS Xp", gg->xp);
        printf("   %24s : %f\n", "GDS Yp", gg->yp);
        printf("   %24s : %f\n", "GDS Orientation", gg->orient);
        printf("   %24s : %d\n", "GDS Nr", gg->nr);
        printf("   %24s : %f\n", "GDS Xo", gg->xo);
        printf("   %24s : %f\n", "GDS Yo", gg->yo); */
    }
    break;
    case GRID_ALBERS:
    case GRID_OLAMBERT:
    case GRID_UTM:
    case GRID_SIMPOL:
    case GRID_MILLER:
    default:
        break;
 }
  
   return lgdat;   
} /* end of get_gdes_data */

int gdes_cmp( gdes* gd1, gdes* gd2){    
/* This function compares two gdes objects to see if they are equal */    
    TEST_NOSERC("gdes_cmp called: ", 0 );
    if( gd1->type != gd2->type ){ CMP("CMP gd1->type");return 0;}
    if( gd1->ncols != gd2->ncols ){ CMP("CMP gd1->ncols");return 0;}
    if( gd1->nrows != gd2->nrows ){ CMP("CMP gd1->nrows"); return 0;}
    if( gd1->npts != gd2->npts ){ CMP("CMP gd1->npts"); return 0;}
    if( gd1->res_flags != gd2->res_flags ){ CMP("CMP gd1->res_flags"); return 0;}
    if( gd1->scan_mode != gd2->scan_mode ){ CMP("CMP gd1->scan_mode"); return 0;}
    if( gd1->quasi != gd2->quasi ){ CMP("CMP gd1->quasi"); return 0;}

    if(!cmp_float_start( gd1->vc, gd1->nv, gd2->vc, gd2->nv )){ CMP("CMP nv"); return 0;}
    if(!cmp_int_start( gd1->lc, gd1->maxlc, gd2->lc, gd2->maxlc )){ CMP("CMP maxlc"); return 0;}

    /* if( gd1->keep != gd2->keep ){CMP("CMP keep");return 0;}// This comparison is no longer done */
    /* because keep only stores the storage status */

    if(!cmp_float_set_base( gd1->vc, gd1->nv, gd2->vc, gd2->nv )){ CMP("CMP vc"); return 0;}
    if(!cmp_int_set_base( gd1->lc, gd1->maxlc, gd2->lc, gd2->maxlc )){ CMP("CMP lc"); return 0;}
    
    if( gengrid_cmp( gd1 , gd2 )){
        TEST_NOSERC("gdes_cmp", 1 );
        return 1;
    }
    else{
        CMP("CMP gdes_cmp");
        TEST_NOSERC("gdes_cmp", 0 );
        return 0;
    }    

    return 1;
} /* end of gdes_cmp */

static int gengrid_cmp( gdes* gd1, gdes* gd2 ){
/* This function compares the gengrid unions in two gdes objects to see if they are equal */    
gengrid* gr1=&(gd1->grid);
gengrid* gr2=&(gd2->grid);

if( gd1->type != gd2->type ) return 0;
 
switch (gd1->type) {
    case GRID_LL:                /* fall through */
    case GRID_RLL:                /* fall through */
    case GRID_SLL:                /* fall through */
    case GRID_SRLL:
    {
        if( !gdes_ll_cmp( &(gr1->ll), &(gr2->ll) ))  
                return 0;
        else
        return 1;
    }
        break;
    case GRID_GAU:                /* fall through */
    case GRID_RGAU:                /* fall through */
    case GRID_SGAU:                /* fall through */
    case GRID_SRGAU:                /* fall through */
    {
        if( !gdes_gau_cmp( &(gr1->gau), &(gr2->gau) ))
        return 0;
        else
        return 1;
    }
    break;
    case GRID_SPH:
    case GRID_RSPH:
    case GRID_SSPH:
    case GRID_SRSPH:
    {
        if( !gdes_sph_cmp( &(gr1->sph), &(gr2->sph) ))
        return 0;
        else
        return 1; 
    }
    break;
    case GRID_MERCAT:
    {
        
        if( !gdes_mercator_cmp( &(gr1->mercator), &(gr2->mercator) ))
        return 0;
        else
        return 1; 
    }
    break;
    case GRID_GNOMON:                /* fall through */
    case GRID_POLARS:
    {
        if( !gdes_polars_cmp( &(gr1->polars), &(gr2->polars) ))
        return 0;
        else
        return 1; 
    }
        break;
    case GRID_LAMBERT:
    {
        if( !gdes_lambert_cmp( &(gr1->lambert), &(gr2->lambert) ))
        return 0;
        else
        return 1;
        
    }
        break;
    case GRID_SPACEV:
    {
        if( !gdes_spacev_cmp( &(gr1->spacev), &(gr2->spacev) ))
        return 0;
        else
        return 1;
    }
    break;
    case GRID_ALBERS:
    case GRID_OLAMBERT:
    case GRID_UTM:
    case GRID_SIMPOL:
    case GRID_MILLER:
    default:
        return 1;
        break;
    }
        return 1;
} /* end of gdes_cmp */
    
static int rotstr_cmp(  rotated* rot1, rotated* rot2, 
                 stretched* strch1, stretched* strch2 ){
/* This function compares two rotated and  two stretched objects to 
   see if both comparisons are equal */ 
    if( rot1 == 0  && rot2 != 0 ) return 0; 
    if( rot1 != 0  && rot2 == 0 ) return 0;
    if( strch1 == 0  && strch2 != 0 ) return 0; 
    if( strch1 != 0  && strch2 == 0 ) return 0; 
    if( rot1 != 0 && rot2 != 0 ){
            if( rot1->lat != rot2->lat) return 0;                        
            if( rot1->lon != rot2->lon) return 0;                
        if( rot1->angle != rot2->angle) return 0;
    }
    if( strch1 != 0 && strch2 != 0 ){
            if( strch1->lat != strch2->lat) return 0;                        
            if( strch1->lon != strch2->lon) return 0;                
        if( strch1->factor != strch2->factor) return 0;
    }
  
    return 1;
}

static int gdes_ll_cmp( gdes_ll* ll1, gdes_ll* ll2 ){
/* This function compares two gdes_ll objects to see if they are equal */     
    if( ll1->ni != ll2->ni ) return 0;
    if( ll1->nj != ll2->nj ) return 0;
    if( ll1->la1 != ll2->la1 ) return 0;
    if( ll1->lo1 != ll2->lo1 ) return 0;
    if( ll1->la2 != ll2->la2 ) return 0;
    if( ll1->lo2 != ll2->lo2 ) return 0;
    if( ll1->di != ll2->di ) return 0;
    if( ll1->dj != ll2->dj ) return 0;
    return rotstr_cmp( ll1->rot, ll2->rot, ll1->strch, ll2->strch ); 
}

static int gdes_gau_cmp( gdes_gau* gau1, gdes_gau* gau2 ){
/* This function compares two gdes_gau objects to see if they are equal */
    if( gau1->ni != gau2->ni ) return 0;
    if( gau1->nj != gau2->nj ) return 0;
    if( gau1->la1 != gau2->la1 ) return 0;
    if( gau1->lo1 != gau2->lo1 ) return 0;
    if( gau1->la2 != gau2->la2 ) return 0;
    if( gau1->lo2 != gau2->lo2 ) return 0;
    if( gau1->di != gau2->di ) return 0;
    if( gau1->n != gau2->n ) return 0;
    return rotstr_cmp( gau1->rot, gau2->rot, gau1->strch, gau2->strch );
}

static int gdes_sph_cmp( gdes_sph* sph1, gdes_sph* sph2 ){
/* This function compares two gdes_sph objects to see if they are equal */
    if( sph1->j != sph2->j ) return 0;
    if( sph1->k != sph2->k ) return 0;
    if( sph1->m != sph2->m ) return 0;
    if( sph1->type != sph2->type ) return 0;
    if( sph1->mode != sph2->mode ) return 0;
    return rotstr_cmp( sph1->rot, sph2->rot, sph1->strch, sph2->strch );
}

static int gdes_polars_cmp( gdes_polars* pol1, gdes_polars* pol2 ){
/* This function compares two gdes_polars objects to see if they are equal */
    if( pol1->nx != pol2->nx ) return 0;                        
    if( pol1->ny != pol2->ny ) return 0;                        
    if( pol1->la1 != pol2->la1 ) return 0;                        
    if( pol1->lo1 != pol2->lo1 ) return 0;
    if( pol1->lov != pol2->lov ) return 0;
    if( pol1->dx != pol2->dx ) return 0;
    if( pol1->dy != pol2->dy ) return 0;
    if( pol1->pole != pol2->pole ) return 0;
    return 1;
}

static int gdes_mercator_cmp( gdes_mercator* m1, gdes_mercator* m2 ){
/* This function compares two gdes_mercator objects to see if they are equal */
    if( m1->ni != m2->ni ) return 0;                        
    if( m1->nj != m2->nj ) return 0;                        
    if( m1->la1 != m2->la1 ) return 0;                        
    if( m1->lo1 != m2->lo1 ) return 0;
    if( m1->la2 != m2->la2 ) return 0;
    if( m1->lo2 != m2->lo2 ) return 0;
    if( m1->latin != m2->latin ) return 0;
    if( m1->di != m2->di ) return 0;
    if( m1->dj != m2->dj ) return 0;
    
    return 1;
}

static int gdes_lambert_cmp( gdes_lambert* lam1, gdes_lambert* lam2 ){
/* This function compares two gdes_lambert objects to see if they are equal */
    if( lam1->nx != lam2->nx ) return 0;
    if( lam1->ny != lam2->ny ) return 0;
    if( lam1->la1 != lam2->la1 ) return 0;
    if( lam1->lo1 != lam2->lo1 ) return 0;
    if( lam1->lov != lam2->lov ) return 0;
    if( lam1->dx != lam2->dx ) return 0;
    if( lam1->dy != lam2->dy ) return 0;
    if( lam1->pole != lam2->pole ) return 0;
    if( lam1->centers != lam2->centers ) return 0;
    if( lam1->latin1 != lam2->latin1 ) return 0;
    if( lam1->latin2 != lam2->latin2 ) return 0;
    if( lam1->splat != lam2->splat ) return 0;
    if( lam1->splon != lam2->splon ) return 0;
    return 1;
}

static int gdes_spacev_cmp( gdes_spacev* spa1, gdes_spacev* spa2 ){
/* This function compares two gdes_spacev objects to see if they are equal */
    if( spa1->nx != spa2->nx ) return 0; 
    if( spa1->ny != spa2->ny ) return 0;
    if( spa1->lap != spa2->lap ) return 0;
    if( spa1->lop != spa2->lop ) return 0;
    if( spa1->dx != spa2->dx ) return 0;
    if( spa1->dy != spa2->dy ) return 0;
    if( spa1->xp != spa2->xp ) return 0;
    if( spa1->yp != spa2->yp ) return 0;
    if( spa1->orient != spa2->orient ) return 0;
    if( spa1->nr != spa2->nr ) return 0;
    if( spa1->xo != spa2->xo ) return 0;
    if( spa1->yo != spa2->yo ) return 0;
    return 1;
}
