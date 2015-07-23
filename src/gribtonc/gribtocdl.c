/*
 * Noserc, based on University Corporation for Atmospheric Research gribdump.c.
 */

/* Remember to update this constant for each change in the program */
/*#define PROGVER "gribtocdl 1.1 - 18.12.2001" */
#define PROGVER "gribtocdl 1.4 - 12.12.2002"

/* 
 * decodes GRIB products into cdl structure ASCII
 */

#include <stdio.h>
#include <signal.h>
#include <errno.h>
#include <unistd.h>
#include <limits.h>		/* _POSIX_PATH_MAX */
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "ulog.h"
#include "mkdirs_open.h"
#include "centers.h"
#include "models.h"
#include "params.h"
#include "levels.h"
#include "timeunits.h"
#include "grib1.h"
#include "product_data.h"
#include "quasi.h"
#include "user_param.h"

#ifdef NO_ATEXIT
#include <atexit.h>
#endif

#define NR 250
#define NR_GDES 25

#include "test_noserc.h" 
int  testNoserc; 

/* Keep the original filename, if any - added 18.12.2001 by ab */
char* grib_file_name;

/* This struct is used to store the different levels in the struct gribdata */
typedef struct level{
int level_flg;       
char* clevelname;    
char* clevelsuffix;
int nflevel;          
float flevel[2][NR]; 
int   level[2][NR];
}level;

typedef struct parameter{
int param;
int nlevel;
int level_flg[NR];
}parameter;

typedef struct lgrib{
    char delim[4];	/* "GRIB", indicator block */ 
    int edition; /* 0 or 1, so far */
    int center;	/* NMC is 7 */
    int subcenter;	/* 0 if none */
}lgrib;

/* This struct is the datatype of the object m_gb where 
   the values from grib header information are 
   stored if not identical values already 
   exist in the object */

typedef struct gribdata{
/* this struct is based on the struct product_data in the file product_data.h */
int nlgrib;
lgrib* llgrib;
int nmodel;
int* model;  /* allocated by the originating center */
int ngrid;
int* grid;   /* Grid ID, may be center-specific */
int nparam_vers;
int* param_vers;/* parameter table version */
int nparam;          /*1*/
parameter* param;   /*1*//* Indicator of parameter */
int nlevel;          
level* gblevel;   
int nreftime;
char** creftime;
int ntunit;
int* tunit;
int ntr_flg;
int* tr_flg;/**/
int ntr;
int* tr;/* has the valtime offsets */
int nhas_gds;
int* has_gds;
int nnpts;
int* npts;
int ngd;
gdes* gd; 
}gribdata;

typedef struct list_gribdata{
int ngb;
gribdata* gb;
}list_gribdata;

static gribdata m_gb;
static list_gribdata m_lgb;

#ifdef __STDC__
static int AddInt( int* arr, int narr, int val );
static int AddCharp( char** arr, int narr, char* val );
static void init_gribdata( gribdata* gb, int nr_gdes );
static void print_gribdata(int ngdes);
static void add_grib(struct product_data* gp, int prec, int ngdes );
static int  do_add(FILE * fp, FILE * ep, int timeout, int prec, quas*);
static void cleanup(void);
static void signal_handler(int sig);
static void set_sigactions(void);
static void usage(char* av0);
static char* parameter_long_name( int i );
static char* level_long_name( int lev );
static void print_time_cdl();
static void print_levels(gribdata* l_gb);
static void print_parameters(gribdata* l_gb, char* idx1, char* idx2);
static int meta_level(int flag);
static int is_ixg(int id);
static product_data* grib_decode(prod*, quas*);
static void print_cdl(); 
int main(int ac, char** av);
#endif

unsigned long num_wmo_messages; /* for statistics on exit */
unsigned long num_gribs_decoded;
unsigned long num_gribs_written;
unsigned long num_degribbed;

/*
 * Timeout in seconds.  If no input is received for this interval, process
 * closes output file and exits.
 */
#define DEFAULT_TIMEOUT  600

/*
 * Default precision used for printing floats.
 */
#define DEFAULT_PRECISION	7

/*
 * Called at exit.
 * This callback routine registered by atexit().
 */
static void
cleanup()
{
    uinfo("Exiting") ;
    uinfo("%lu WMO msgs, %lu GRIBs decoded, %lu written ",
	  num_wmo_messages, num_gribs_decoded, num_gribs_written);
    (void) closeulog();
}

/*
 * Called upon receipt of signals.
 * This callback routine registered in set_sigactions() .
 */
static void
signal_handler(sig)
     int sig ;
{
#ifdef SVR3SIGNALS
    /* 
     * Some systems reset handler to SIG_DFL upon entry to handler.
     * In that case, we reregister our handler.
     */
    (void) signal(sig, signal_handler) ;
#endif
    switch(sig) {
      case SIGHUP :
	udebug("SIGHUP") ;
	return ;
      case SIGINT :
	unotice("Interrupt") ;
	exit(0) ;
      case SIGTERM :
	udebug("SIGTERM") ;
	exit(0) ;
      case SIGUSR1 :
	udebug("SIGUSR1") ;
	return ;
      case SIGUSR2 :
	if (toggleulogpri(LOG_INFO))
	  unotice("Going verbose") ;
	else
	  unotice("Going silent") ;
	return ;
      case SIGPIPE :
	unotice("SIGPIPE") ;
	exit(0) ;
    }
    udebug("signal_handler: unhandled signal: %d", sig) ;
}

static void
set_sigactions()
{
#ifndef NO_POSIXSIGNALS
    struct sigaction sigact ;
    
    sigact.sa_handler = signal_handler;
    sigemptyset(&sigact.sa_mask) ;
    sigact.sa_flags = 0 ;
    
    (void) sigaction(SIGHUP, &sigact, NULL) ;
    (void) sigaction(SIGINT, &sigact, NULL) ;
    (void) sigaction(SIGTERM, &sigact, NULL) ;
    (void) sigaction(SIGUSR1, &sigact, NULL) ;
    (void) sigaction(SIGUSR2, &sigact, NULL) ;
    (void) sigaction(SIGPIPE, &sigact, NULL) ;
#else
    (void) signal(SIGHUP, signal_handler) ;
    (void) signal(SIGINT, signal_handler) ;
    (void) signal(SIGTERM, signal_handler) ;
    (void) signal(SIGUSR1, signal_handler) ;
    (void) signal(SIGUSR2, signal_handler) ;
    (void) signal(SIGPIPE, signal_handler) ;
#endif
}


static void
usage(av0)
     char *av0 ; /*  id string */
{
    (void)fprintf(stderr,
		  "Usage: %s [options] [files ...]\t\nOptions:\n", av0);
    (void)fprintf(stderr,
		  "\t-v           Verbose, report decoding steps\n") ;
    (void)fprintf(stderr,
		  "\t-g tabfile   Use grib-netcdf translation defined in this file\n") ;
    (void)fprintf(stderr,
		  "\t-t timeout   If no input, exit after \"timeout\" seconds (default %d)\n",
		  DEFAULT_TIMEOUT) ;
    (void)fprintf(stderr,
		  "\t-q meth      Specs for expanding quasi-regular grids,\n\t\t\te.g. -q \"lin,dlat=2.5,dlon=5.0\"\n") ;
    exit(1);
}

int headers_only;
int brief_output;
int gdes_analyse;
int printcdl;

static int
is_ixg(id)			/* is an international exchange grid? */
     int id;
{
    switch (id) {
      case 21:
      case 22:
      case 23:
      case 24:
      case 25:
      case 26:
      case 37:
      case 38:
      case 39:
      case 40:
      case 41:
      case 42:
      case 43:
      case 44:
      case 50:
      case 61:
      case 62:
      case 63:
      case 64:
	return 1;
    }
    /* default: */
    return 0;
}

/*****************************************************/
static void print_gribdata( int ngdes ){
/* This function prints to screen from the objects m_lgb.gb[k] for all k, or the object m_gb */ 
 int i=0;
 int j=0;
 int k=0;
 int max=1;
 gribdata* l_gb;

 if(ngdes){ 
    max = m_lgb.ngb;
 }
 else 
    l_gb = &m_gb;

 for( k=0; k < max; k++ ){
   if(ngdes){ 
       l_gb = &m_lgb.gb[k];
   }
   printf("   %24s : %d\n","nmodel", l_gb->nmodel);
   for( i=0; i< l_gb->nmodel; i++ ){
     for( j=0; j<l_gb->nlgrib; j++ ){
	 printf("   %24s : %d %s\n","model", l_gb->model[i], modelname(l_gb->llgrib[j].center,l_gb->model[i]) );
     }
   }
   printf("   %24s : %d\n","ngrid", l_gb->ngrid);
   for( i=0; i< l_gb->ngrid; i++ )
     printf("   %24s : %d   %d \n","grid", l_gb->grid[i], is_ixg(l_gb->grid[i]));
   printf("   %24s : %d\n","nnpts", l_gb->nnpts);
   for( i=0; i< l_gb->nnpts; i++ )
     printf("   %24s : %d\n","npts", l_gb->npts[i]);
   printf("   %24s : %d\n","nparam", l_gb->nparam);
   for( i=0; i< l_gb->nparam; i++ )
     printf("   %24s : %d  %s  %s\n","param", l_gb->param[i].param, grib_pname(l_gb->param[i].param), grib_units(l_gb->param[i].param));

   printf("   %24s : %d\n","nlevel", l_gb->nlevel);
   for( i=0; i< l_gb->nlevel; i++ ){
      printf("   %24s : %d  %s : %s/%s \n","level_flg",l_gb->gblevel[i].level_flg,
	     "clevel",l_gb->gblevel[i].clevelname,l_gb->gblevel[i].clevelsuffix );
      printf("   %24s : %d\n","nflevel", l_gb->gblevel[i].nflevel);

      for( j=0; j<l_gb->gblevel[i].nflevel; j++){
	 int level[2];

	 if( l_gb->gblevel[i].flevel[0][j] >-1 )  
	    printf("   %24s : %f  ","flevel", l_gb->gblevel[i].flevel[0][j] );
	 if( l_gb->gblevel[i].flevel[1][j] > -1 )
	    printf("  %f", l_gb->gblevel[i].flevel[1][j] );

	 level[0] = l_gb->gblevel[i].level[0][j];
	 level[1] = l_gb->gblevel[i].level[1][j];
	 printf(" ** %4d %4d", level1(l_gb->gblevel[i].level_flg, level),
	       level2(l_gb->gblevel[i].level_flg, level) ); 
	 printf("\n");	
      }
   }

   printf("   %24s : %d\n","nreftime", l_gb->nreftime);
   for( i=0; i< l_gb->nreftime; i++ )
     printf("   %24s : %s\n","creftime", l_gb->creftime[i]);

   printf("   %24s : %d\n","ntunit", l_gb->ntunit);
   for( i=0; i< l_gb->ntunit; i++ )
     printf("   %24s : %s\n","tunitsname", tunitsname(l_gb->tunit[i]));
   printf("   %24s : %d\n","ntr_flg", l_gb->ntr_flg);
   for( i=0; i< l_gb->ntr_flg; i++ )
     printf("   %24s : %s\n","Time Range Indicator",  triname(l_gb->tr_flg[i]));
   printf("   %24s : %d\n","ntr", l_gb->ntr);
   for( i=0; i< l_gb->ntr; i++ )
     printf("   %24s : %d\n","ValTime_offset",  l_gb->tr[i]);
   printf("   %24s : %d\n","nhas_gds", l_gb->nhas_gds);
   for( i=0; i< l_gb->nhas_gds; i++ )
      printf("   %24s : %d\n","has_gds", l_gb->has_gds[i]);    
  printf("   %24s : %d\n","nedition==nlgrib", l_gb->nlgrib); 
  for( i=0; i< l_gb->nlgrib; i++ )
     printf("   %24s : %d\n","edition", l_gb->llgrib[i].edition);
   printf("   %24s : %d\n","nparam_vers", l_gb->nparam_vers);
   for( i=0; i< l_gb->nparam_vers; i++ )
     printf("   %24s : %d\n","param_vers", l_gb->param_vers[i]); 
   printf("   %24s : %d\n","ngd", l_gb->ngd);
   for( i=0; i< l_gb->ngd; i++ )
     print_gdes(&(l_gb->gd[i]));
 printf("****************************************** \n");
 printf("****************************************** \n");

 }/*end for k */
} /* end of print_gribdata */

static int meta_level(flag)
/* This function is directly based on the function level1 in the file levels.c*/
/* This function categorises the different levels in different meta levels. */
    int flag;			/* GRIB level flag */
{
    switch(flag){
    case LEVEL_SURFACE: 
    case LEVEL_CLOUD_BASE: 
    case LEVEL_CLOUD_TOP: 
    case LEVEL_ISOTHERM: 
    case LEVEL_ADIABAT: 
    case LEVEL_MAX_WIND: 
    case LEVEL_TROP: 
    case LEVEL_TOP:
    case LEVEL_SEABOT:
    case LEVEL_MEAN_SEA: 
    case LEVEL_ATM:
    case LEVEL_OCEAN:
    case LEVEL_HTFL:
    case LEVEL_GCBL:
    case LEVEL_GCTL:
    case LEVEL_BCBL:
    case LEVEL_BCTL:
    case LEVEL_BCY:
    case LEVEL_LCBL:
    case LEVEL_LCTL:
    case LEVEL_LCY:
    case LEVEL_CEIL:
    case LEVEL_PBLRI:
    case LEVEL_MCBL:
    case LEVEL_MCTL:
    case LEVEL_MCY:
    case LEVEL_HCBL:
    case LEVEL_HCTL:
    case LEVEL_HCY:
    case LEVEL_OITL:
    case LEVEL_OBML:
    case LEVEL_OBIL:
    case LEVEL_S26CY:
    case LEVEL_OMXL:
    case LEVEL_CCBL:
    case LEVEL_CCTL:
    case LEVEL_CCY:
    case LEVEL_LLTW:
    case LEVEL_MTHE:
    case LEVEL_EHLT:
    case LEVEL_SCBL:
    case LEVEL_SCTL:
    case LEVEL_DCBL:
    case LEVEL_DCTL:
    case LEVEL_LBLSW:
    case LEVEL_HTLSW:
    case LEVEL_FL:
	return 1;

    case LEVEL_TMPL:
    case LEVEL_ISOBARIC: 
    case LEVEL_FH: 
    case LEVEL_FHG: 
    case LEVEL_SIGMA: 
    case LEVEL_HY:
    case LEVEL_BLS: 
    case LEVEL_ISEN:
    case LEVEL_PDG:
    case LEVEL_PV:
    case LEVEL_ETAL:
    case LEVEL_FHGH:
    case LEVEL_DBS: 
	return 2;	/* 2-octet level */

    case LEVEL_LISO: 
    case LEVEL_LFHM: 
    case LEVEL_LFHG: 
    case LEVEL_LS: 
    case LEVEL_LHY:
    case LEVEL_LBLS: 
    case LEVEL_LISEN:
    case LEVEL_LPDG:
    case LEVEL_LETA: 
    case LEVEL_LISH: 
    case LEVEL_LSH: 
    case LEVEL_LISM: 
    case LEVEL_OLYR:
	return 3;		/* 1-octet level */
    }
    /* default */
    printf("ERROR Level=%d\n", flag);
    return -1;
} /* end meta_level */

static void print_time_cdl(){
/* This function prints the time variables to the standard output */
/* CALLED BY: print_cdl()*/
printf("        double reftime(record) ;	// reference time of the model\n");
printf("               reftime:long_name = \"reference time\" ;\n");
printf("               reftime:units = \"hours since 1992-1-1\" ;\n");
printf("\n");
printf("        double valtime(record) ;	// forecast time (\"valid\" time)\n");
printf("               valtime:long_name = \"valid time\" ;\n");
printf("               valtime:units = \"hours since 1992-1-1\" ;\n");
printf("\n");
printf("        :record = \"reftime, valtime\" ;	// \"dimension attribute\" -- means\n");
printf("                                        // (reftime, valtime) uniquely\n");
printf("                                        // determine record\n");
printf("\n");
printf("        char   datetime(record, time_len) ; // derived from reftime\n");
printf("               datetime:long_name = \"reference date and time\" ;\n");
printf("               // units YYYY-MM-DD hh:mm:ssZ  (ISO 8601)\n");
printf("\n");
printf("        double valtime_offset(valtime_offset) ; // valtime - reftime\n");
printf("               valtime_offset:long_name = \"hours from reference time\" ;\n");
printf("               valtime_offset:units = \"hours\" ;\n");
printf("\n");
printf("        char   forecasttime(record, time_len) ; // derived from valtime\n");
printf("               forecasttime:long_name = \"forecast date and time\" ;\n");
printf("               // units YYYY-MM-DD hh:mm:ssZ  (ISO 8601)\n");
printf("\n");
}

static char* parameter_long_name( int i ){
/* this function is based upon the definitions and parameters in params.h */
/* this function returns the longname to the different parameters based on the parameterflag */

  /* First check if anything was user-defined */
  char *p;
  p=user_plongname(i);
  if (p!=NULL) return p;

 switch (i) {    
    case PARM_RESERVED   : return "Reserved";
    case PARM_PRESSURE	 : return "Pressure";
    case PARM_PMSL	 : return "Pressure reduced to MSL";
    case PARM_PTND	 : return "Pressure tendency";
    case PARM_ICAHT	 : return "ICAO standard atmosphere reference height";
    case PARM_GPT	 : return "Geopotential";
    case PARM_GPT_HGT	 : return "Geopotential height";
    case PARM_GEOM_HGT	 : return "Geometric height";
    case PARM_HSTDV	 : return "Standard deviation of height";
    case PARM_TOZNE	 : return "Total ozone";
    case PARM_TEMP	 : return "Temperature";
    case PARM_VTEMP	 : return "Virtual temperature";
    case PARM_POT_TEMP	 : return "Potential temperature";
    case PARM_APOT_TEMP	 : return "Pseudo-adiabatic potential temperature";
    case PARM_MAX_TEMP	 : return "Maximum temperature";
    case PARM_MIN_TEMP	 : return "Minimum temperature";
    case PARM_DP_TEMP	 : return "Dew point temperature";
    case PARM_DP_DEP	 : return "Dew point depression (or deficit)";
    case PARM_LAPSE	 : return "Lapse rate";

    case PARM_VIS	 : return "Visibility";

    case PARM_RAD1	 : return "Radar spectra, direction and frequency";
    case PARM_RAD2	 : return "Radar spectra, direction and radial num";
    case PARM_RAD3	 : return "Radar spectra, radial num and radial num";
    case PARM_PLI	 : return "Parcel lifted index (to 500 hPa)";

    case PARM_TANOM	 : return "Temperature anomaly";
    case PARM_PANOM	 : return "Pressure anomaly";
    case PARM_ZANOM	 : return "Geopotential height anomaly";
    case PARM_WAV1	 : return "Wave spectra, direction and frequency";
    case PARM_WAV2	 : return "Wave spectra, direction and radial num";
    case PARM_WAV3	 : return "Wave spectra, radial num and radial num";
    case PARM_WND_DIR	 : return "Wind direction";
    case PARM_WND_SPEED	 : return "Wind speed";
    case PARM_U_WIND	 : return "u-component of wind";
    case PARM_V_WIND	 : return "v-component of wind";
    case PARM_STRM_FUNC	 : return "Stream function";
    case PARM_VPOT	 : return "Velocity potential";

    case PARM_MNTSF	 : return "Montgomery stream function";

    case PARM_SIG_VEL	 : return "Sigma coordinate of vertical velocity";
    case PARM_VERT_VEL	 : return "Pressure vertical velocity";
    case PARM_GEOM_VEL	 : return "Geometric vertical velocity";
    case PARM_ABS_VOR	 : return "Absolute vorticity";
    case PARM_ABS_DIV	 : return "Absolute divergence";
    case PARM_REL_VOR	 : return "Relative vorticity";
    case PARM_REL_DIV	 : return "Relative divergence";
    case PARM_U_SHR	 : return "Vertical u-component shear";
    case PARM_V_SHR	 : return "Vertical v-component shear";
    case PARM_CRNT_DIR	 : return "Direction of current";
    case PARM_CRNT_SPD	 : return "Speed of current";
    case PARM_U_CRNT	 : return "u-component of current";
    case PARM_V_CRNT	 : return "v-component of current";
    case PARM_SPEC_HUM 	 : return "Specific humidity";
    case PARM_REL_HUM	 : return "Relative humidity";
    case PARM_HUM_MIX	 : return "Humidity mixing ratio";
    case PARM_PR_WATER	 : return "Precipitable water";
    case PARM_VAP_PR	 : return "Vapor pressure";
    case PARM_SAT_DEF	 : return "Saturation deficit";
    case PARM_EVAP	 : return "Evaporation";

    case PARM_C_ICE	 : return "Cloud ice";

    case PARM_PRECIP_RT	 : return "Precipitation rate";
    case PARM_THND_PROB	 : return "Thunderstorm probability";
    case PARM_PRECIP_TOT : return "Total precipitation";
    case PARM_PRECIP_LS	 : return "Large scale precipitation";
    case PARM_PRECIP_CN	 : return "Convective precipitation";
    case PARM_SNOW_RT	 : return "Snowfall rate water equivalent";
    case PARM_SNOW_WAT	 : return "Water equivalent of accumulated snow depth";
    case PARM_SNOW	 : return "Snow depth";
    case PARM_MIXED_DPTH : return "Mixed layer depth";
    case PARM_TT_DEPTH	 : return "Transient thermocline depth";
    case PARM_MT_DEPTH	 : return "Main thermocline depth";
    case PARM_MTD_ANOM	 : return "Main thermocline anomaly";
    case PARM_CLOUD	 : return "Total cloud cover";
    case PARM_CLOUD_CN	 : return "Convective cloud cover";
    case PARM_CLOUD_LOW	 : return "Low cloud cover";
    case PARM_CLOUD_MED	 : return "Medium cloud cover";
    case PARM_CLOUD_HI	 : return "High cloud cover";
    case PARM_CLOUD_WAT	 : return "Cloud water";
    case PARM_BLI	 : return "Best lifted index (to 500 hPa)";

    case PARM_SNO_C	 : return "Convective snow";
    case PARM_SNO_L	 : return "Large scale snow";

    case PARM_SEA_TEMP	 : return "Sea temperature";
    case PARM_LAND_MASK	 : return "Land-sea mask (1=land; 0=sea)";
    case PARM_SEA_MEAN	 : return "Deviation of sea level from mean";
    case PARM_SRF_RN	 : return "Surface roughness";
    case PARM_ALBEDO	 : return "Albedo";
    case PARM_SOIL_TEMP	 : return "Soil temperature";
    case PARM_SOIL_MST 	 : return "Soil moisture content";
    case PARM_VEG	 : return "Vegetation";
    case PARM_SAL	 : return "Salinity";
    case PARM_DENS	 : return "Density";

    case PARM_WATR	 : return "Water runoff";

    case PARM_ICE_CONC	 : return "Ice concentration (ice=l; no ice=O)";
    case PARM_ICE_THICK	 : return "Ice thickness";
    case PARM_ICE_DIR	 : return "Direction of ice drift";
    case PARM_ICE_SPD	 : return "Speed of ice drift";
    case PARM_ICE_U	 : return "u-component of ice drift";
    case PARM_ICE_V	 : return "v-component of ice drift";
    case PARM_ICE_GROWTH : return "Ice growth";
    case PARM_ICE_DIV	 : return "Ice divergence";

    case PARM_SNO_M	 : return "Snow melt";

    case PARM_WAVE_HGT	 : return "Significant height of combined wind waves and swell";
    case PARM_SEA_DIR	 : return "Direction of wind waves";
    case PARM_SEA_HGT	 : return "Significant height of wind waves";
    case PARM_SEA_PER	 : return "Mean period of wind waves";
    case PARM_SWELL_DIR	 : return "Direction of swell waves";
    case PARM_SWELL_HGT	 : return "Significant height of swell waves";
    case PARM_SWELL_PER	 : return "Mean period of swell waves";
    case PARM_WAVE_DIR	 : return "Primary wave direction";
    case PARM_WAVE_PER	 : return "Primary wave mean period";
    case PARM_WAVE2_DIR	 : return "Secondary wave direction";
    case PARM_WAVE2_PER	 : return "Secondary wave mean period";
    case PARM_RDN_SWSRF	 : return "Net shortwave radiation (surface)";
    case PARM_RDN_LWSRF	 : return "Net longwave radiation (surface)";
    case PARM_RDN_SWTOP	 : return "Net shortwave radiation (top of atmos.)";
    case PARM_RDN_LWTOP	 : return "Net longwave radiation (top of atmos.)";
    case PARM_RDN_LW	 : return "Long wave radiation";
    case PARM_RDN_SW	 : return "Short wave radiation";
    case PARM_RDN_GLBL	 : return "Global radiation";
    case PARM_BRTMP	 : return "Brightness temperature";
    case PARM_LWRAD	 : return "Long-wave radiation";
    case PARM_SWRAD	 : return "Short-wave radiation";

    case PARM_LAT_HT	 : return "Latent heat flux";
    case PARM_SEN_HT	 : return "Sensible heat flux";
    case PARM_BL_DISS	 : return "Boundary layer dissipation";

    case PARM_U_FLX	 : return "Momentum flux, u component";
    case PARM_V_FLX	 : return "Momentum flux, v component";
    case PARM_WMIXE	 : return "Wind mixing energy";

    case PARM_IMAGE	 : return "Image data";
    case PARM_MSLSA	 : return "Mean sea level pressure (std. atms. reduction)";
    case PARM_PM	 : return "Mean sea level pressure (MAPS system reduction)";
    case PARM_MSLET	 : return "Mean sea level pressure (ETA model reduction)";
    case PARM_LIFT_INDX	 : return "Surface lifted index";
    case PARM_LIFT_INDX4 : return "Best (4 layer) lifted index";
    case PARM_K_INDX	 : return "K index";
    case PARM_SW_INDX	 : return "Sweat index";
    case PARM_HM_DIV	 : return "Horizontal moisture divergence";
    case PARM_VERT_SSHR	 : return "Vertical speed shear";
    case PARM_TSLSA	 : return "3-hr pressure tendency (std. atms. reduction)";
    case PARM_BVF_2	 : return "Brunt-vaisala frequency (squared)";
    case PARM_PV_MW	 : return "Potential vorticity (density weighted)";
    case PARM_CRAIN	 : return "Categorical rain  (yes=1; no=0)";
    case PARM_CFRZRN	 : return "Categorical freezing rain  (yes=1; no=0)";
    case PARM_CICEPL	 : return "Categorical ice pellets  (yes=1; no=0)";
    case PARM_CSNOW	 : return "Categorical snow  (yes=1; no=0)";
    case PARM_SOILW	 : return "Volumetric soil moisture content";
    case PARM_PEVPR	 : return "Potential evaporation rate";
    case PARM_CWORK	 : return "Cloud workfunction";
    case PARM_U_GWD	 : return "Zonal flux of gravity wave stress";
    case PARM_V_GWD	 : return "Meridional flux of gravity wave stress";
    case PARM_PV	 : return "Potential vorticity";
    case PARM_COVMZ	 : return "Covariance between meridional and zonal components of wind";
    case PARM_COVTZ	 : return "Covariance between temperature and zonal component of wind";
    case PARM_COVTM	 : return "Covariance between temperature and meridional component of wind";
    case PARM_CLWMR	 : return "Cloud water";
    case PARM_O3MR	 : return "Ozone mixing ratio";

    case PARM_GFLUX	 : return "Ground heat flux";
    case PARM_CIN	 : return "Convective inhibition";
    case PARM_CAPE	 : return "Convective available potential energy";
    case PARM_TKE	 : return "Turbulent kinetic energy";
    case PARM_CONDP	 : return "Condensation pressure of parcel lifted from indicated surface";
    case PARM_CSUSF	 : return "Clear sky upward solar flux";
    case PARM_CSDSF	 : return "Clear sky downward solar flux";
    case PARM_CSULF	 : return "Clear sky upward long wave flux";
    case PARM_CSDLF	 : return "Clear sky downward long wave flux";
    case PARM_CFNSF	 : return "Cloud forcing net solar flux";
    case PARM_CFNLF	 : return "Cloud forcing net long wave flux";
    case PARM_VBDSF	 : return "Visible beam downward solar flux";
    case PARM_VDDSF	 : return "Visible diffuse downward solar flux";
    case PARM_NBDSF	 : return "Near IR beam downward solar flux";
    case PARM_NDDSF	 : return "Near IR diffuse downward solar flux";
    case PARM_RWMR	 : return "Rain water mixing ratio";
    case PARM_SNMR	 : return "Snow water mixing ratio";
    case PARM_M_FLX	 : return "Momentum flux";
    case PARM_LMH	 : return "Mass point model surface";
    case PARM_LMV	 : return "Velocity point model surface";
    case PARM_MLYNO	 : return "Model layer number (from bottom up)";
    case PARM_NLAT    	 : return "Latitude (-90 to +90)";
    case PARM_ELON	 : return "East longitude (0-360)";
    case PARM_ICMR	 : return "Ice water mixing ratio";
    case PARM_GRMR	 : return "Graupel water mixing ratio";
    case PARM_GUST	 : return "Surface wind gust";
    case PARM_LPS_X	 : return "x-gradient of log pressure";
    case PARM_LPS_Y	 : return "y-gradient of log pressure";
    case PARM_HGT_X	 : return "x-gradient of height";
    case PARM_HGT_Y	 : return "y-gradient of height";
    case PARM_TPFI	 : return "Turbulence potential forecast index";
    case PARM_TIPD	 : return "Total icing potential diagnostic";
    case PARM_LTNG	 : return "Lightning";

    case PARM_VPTMP	 : return "Virtual potential temperature";
    case PARM_HLCY	 : return "Storm relative helicity";
    case PARM_PROB	 : return "Probability from ensemble";
    case PARM_PROBN      : return "Probability from ensemble normalized with respect to climate expectancy";
    case PARM_POP	 : return "Probability of precipitation";
    case PARM_CPOFP	 : return "Probability of frozen precipitation";
    case PARM_CPOZP	 : return "Probability of freezing precipitation";
    case PARM_USTM	 : return "u-component of storm motion";
    case PARM_VSTM	 : return "v-component of storm motion";
    case PARM_NCIP	 : return "Number concentration for ice particles";
    case PARM_EVBS	 : return "Direct evaporation from bare soil";
    case PARM_EVCW	 : return "Canopy water evaporation";
    case PARM_NOICE_WAT	 : return "Ice-free water surface";
    case PARM_CWDI	 : return "Convective weather detection index";

    case PARM_DSWRF	 : return "Downward short wave rad. flux";
    case PARM_DLWRF	 : return "Downward long wave rad. flux";
    case PARM_UVPI	 : return "Ultra violet potential index";
    case PARM_MSTR_AVL	 : return "Moisture availability";
    case PARM_XCHG_COF	 : return "Exchange coefficient";
    case PARM_NMIX_LYRS	 : return "Number of mixed layers next to surface";
    case PARM_TRANS	 : return "Transpiration";

    case PARM_USWRF	 : return "Upward short wave rad. flux";
    case PARM_ULWRF	 : return "Upward long wave rad. flux";

    case PARM_CLOUD_NCN	 : return "Amount of non-convective cloud";

    case PARM_CPRAT	 : return "Convective precipitation rate";
    case PARM_TTDIA	 : return "Temperature tendency by all physics";

    case PARM_RDN_TTND	 : return "Temperature tendency by all radiation";

    case PARM_TTPHY	 : return "Temperature tendency by non-radiation physics";
    case PARM_PREIX	 : return "precip.index(0.0-1.00)(see note)";
    case PARM_TSD1D	 : return "Std. Dev. of IR T over 1x1 deg area";

    case PARM_LN_PRES	 : return "Natural log of surface pressure";
    case PARM_HPBL	 : return "Planetary boundary layer height";
    case PARM_GPT_HGT5	 : return "5-wave geopotential height";

    case PARM_C_WAT	 : return "Plant canopy surface water";
    case PARM_SOTYP	 : return "Soil type (as in Zobler)";
    case PARM_VGTYP	 : return "Vegitation type (as in SiB)";
    case PARM_BMIXL	 : return "Blackadar's mixing length scale";
    case PARM_AMIXL	 : return "Asymptotic mixing length scale";
    case PARM_PEVAP	 : return "Potential evaporation";
    case PARM_SNOHF	 : return "Snow phase-change heat flux";
    case PARM_GPT_HGTA5	 : return "5-wave geopotential height anomaly";
    case PARM_MFLUX	 : return "Convective cloud mass flux";
    case PARM_DTRF	 : return "Downward total radiation flux";
    case PARM_UTRF	 : return "Upward total radiation flux";
    case PARM_BGRUN	 : return "Baseflow-groundwater runoff";
    case PARM_SSRUN	 : return "Storm surface runoff";
    case PARM_SIPD	 : return "Supercooled Large Droplet (SLD) icing potential diagnostic ";
    case PARM_OZTOT	 : return "Total ozone";
    case PARM_SNO_CVR	 : return "Snow cover";
    case PARM_SNO_T	 : return "Snow temperature";
    case PARM_COVTW	 : return "Covariance between temperature and vertical component of the wind";
    case PARM_LRGHR	 : return "Large scale condensation heat rate";
    case PARM_CNVHR	 : return "Deep convective heating rate";
    case PARM_CNVMR	 : return "Deep convective moistening rate";
    case PARM_SHAHR	 : return "Shallow convective heating rate";
    case PARM_SHAMR	 : return "Shallow convective moistening rate";
    case PARM_VDFHR	 : return "Vertical diffusion heating rate";
    case PARM_VDFUA	 : return "Vertical diffusion zonal acceleration";
    case PARM_VDFVA	 : return "Vertical diffusion meridional accel";
    case PARM_VDFMR	 : return "Vertical diffusion moistening rate";
    case PARM_SWHR	 : return "Solar radiative heating rate";
    case PARM_LWHR	 : return "Long wave radiative heating rate";
    case PARM_CD	 : return "Drag coefficient";
    case PARM_FRICV	 : return "Friction velocity";
    case PARM_RI	 : return "Richardson number";

    case PARM_MISSING	 : return "Missing";
 }
 return "Parameter not found";
} /* end of parameter_long_name */

static char* level_long_name( int lev ){
/* this function is based upon the function levelsuffix in levels.c */
/* this function returns the longname to the different levels based on the levelflag */
    switch (lev){
    case LEVEL_SURFACE: 
	return "Surface";
    case LEVEL_CLOUD_BASE: 
	return "Cloud Base";
    case LEVEL_CLOUD_TOP: 
	return "Cloud Top";
    case LEVEL_ISOTHERM: 
	return "0 Isotherm";
    case LEVEL_ADIABAT: 
	return "Adiabatic Condensation";
    case LEVEL_MAX_WIND: 
	return "Maximum Wind";
    case LEVEL_TROP: 
	return "Tropopause";
    case LEVEL_TOP:
	return "Top of Atmosphere";
    case LEVEL_SEABOT:
	return "Sea Bottom";
    case LEVEL_TMPL:
	return "Temperature in 1/100 K";
    case LEVEL_ISOBARIC: 
	return "Isobaric";
    case LEVEL_LISO: 
	return "Layer Between Two Isobaric";
    case LEVEL_MEAN_SEA: 
	return "Mean Sea";
    case LEVEL_FH: 
	return "Fixed Height";
    case LEVEL_LFHM: 
	return "Layer Between Two Heights Above MSL";
    case LEVEL_FHG: 
	return "Fixed Height Above Ground";
    case LEVEL_LFHG: 
	return "Layer Between Two Fixed Heights Above Ground";
    case LEVEL_SIGMA: 
	return "Sigma";
    case LEVEL_LS: 
	return "Layer Between Two Sigma";
    case LEVEL_HY:
	return "Hybrid level";
    case LEVEL_LHY:
	return "Layer between 2 hybrid levels";
    case LEVEL_BLS: 
	return "Below Land Surface";
    case LEVEL_LBLS: 
	return "Layer Between Two Depths Below Land Surface";
    case LEVEL_ISEN:
	return "Isentropic (theta) level";
    case LEVEL_LISEN:
	return "Layer between 2 isentropic (theta) levels";
    case LEVEL_PDG:
	return "level at specified pressure difference from ground to level";
    case LEVEL_LPDG:
	return "layer between 2 levels at specified pressure differences from ground to levels";
    case LEVEL_PV:
	return "potential vorticity";
    case LEVEL_ETAL:
	return "ETA level";
    case LEVEL_LETA: 
	return "Layer between two ETA levels";
    case LEVEL_LISH: 
	return "Layer Between Two Isobaric Surfaces, High Precision";
    case LEVEL_FHGH:
	return "Height level above ground (high precision)";
    case LEVEL_LSH: 
	return "Layer Between Two Sigma Levels, High Precision";
    case LEVEL_LISM: 
	return "Layer Between Two Isobaric Surfaces, Mixed Precision";
    case LEVEL_DBS: 
	return "Depth Below Sea";
    case LEVEL_ATM:
	return "Entire atmosphere considered as a single layer";
    case LEVEL_OCEAN:
	return "Entire ocean considered as a single layer";
    case LEVEL_HTFL:
	return "Highest tropospheric freezing level";
    case LEVEL_GCBL:
	return "Grid scale cloud bottom level";
    case LEVEL_GCTL:
	return "Grid scale cloud top level";
    case LEVEL_BCBL:
	return "Boundary layer cloud bottom level";
    case LEVEL_BCTL:
	return "Boundary layer cloud top level";
    case LEVEL_BCY:
	return "Boundary layer cloud layer";
    case LEVEL_LCBL:
	return "Low cloud bottom level";
    case LEVEL_LCTL:
	return "Low cloud top level";
    case LEVEL_LCY:
	return "Low cloud layer";
    case LEVEL_CEIL: 
        return "Cloud ceiling";
    case LEVEL_PBLRI: 
        return "Planetary Boundary Layer";
    case LEVEL_MCBL:
	return "Middle cloud bottom level";
    case LEVEL_MCTL:
	return "Middle cloud top level";
    case LEVEL_MCY:
	return "Middle cloud layer";
    case LEVEL_HCBL:
	return "High cloud bottom level";
    case LEVEL_HCTL:
	return "High cloud top level";
    case LEVEL_HCY:
	return "High cloud layer";
    case LEVEL_OITL: 
        return "Ocean Isotherm Level";
    case LEVEL_OLYR: 
        return "Layer between two depths below ocean surface";
    case LEVEL_OBML: 
        return "Bottom of Ocean Mixed Layer";
    case LEVEL_OBIL: 
        return "Bottom of Ocean Isothermal Layer";
    case LEVEL_S26CY: 
        return "Layer Ocean Surface and 26C Ocean Isothermal Level";
    case LEVEL_OMXL: 
        return "Ocean Mixed Layer";
    case LEVEL_CCBL:
	return "Convective cloud bottom level";
    case LEVEL_CCTL:
	return "Convective cloud top level";
    case LEVEL_CCY:
	return "Convective cloud layer";
    case LEVEL_LLTW:
	return "Lowest level of the wet bulb zero";
    case LEVEL_MTHE:
 	return "maximum equivalent potential temperature level";
    case LEVEL_EHLT:
 	return "equilibrium level";
    case LEVEL_SCBL:
	return "Shallow convective cloud bottom level";
    case LEVEL_SCTL:
	return "Shallow convective cloud top level";
    case LEVEL_DCBL:
	return "Deep convective cloud bottom level";
    case LEVEL_DCTL:
	return "Deep convective cloud top level";
    case LEVEL_LBLSW: 
        return "Lowest bottom level of supercooled liquid water layer";
    case LEVEL_HTLSW: 
        return "Highest top level of supercooled liquid water layer";
    case LEVEL_FL:
 	return "flight_level";
    }
				/* default: */
    uerror("bad level flag: %d", lev);
    return "";
} /* end of level_long_name */

static void  print_levels(gribdata* l_gb){
/* INPUT: gribdata* l_gb */
/* This function prints the level variables to the standard output */
/* CALLED BY: print_cdl()*/
int i;
for( i=0; i<l_gb->nlevel; i++ ){
    int meta = meta_level(l_gb->gblevel[i].level_flg);

    if(meta==2){
	switch(l_gb->gblevel[i].level_flg){
	   case LEVEL_ISOBARIC:
	   {
	       
	       printf("        float  level(level) ;\n");
	       printf("               level:long_name = \"isobaric level\" ;\n");
	       printf("               level:units = \"hectopascals\" ;\n");
	       printf("\n");
	   }
	   break;
	   default:
	   {
	       char* suf=levelsuffix(l_gb->gblevel[i].level_flg);
	       char* unit=levelunits(l_gb->gblevel[i].level_flg);
	       printf("        float  %s(%s) ;\n", suf,suf);
               printf("               %s:long_name = ",suf);
               printf("\"%s\" ;\n",level_long_name(l_gb->gblevel[i].level_flg));
               printf("               %s:units = \"%s\" ;\n", suf, unit );
	       printf("\n");
	   }
	}
    }
    if(meta==3){ 
       char* suf=levelsuffix(l_gb->gblevel[i].level_flg);
       char* unit=levelunits(l_gb->gblevel[i].level_flg);
       printf("        :%s = \"%s_bot, %s_top\" ; // (\"%s_bot, %s_top\") uniquely\n", suf,suf ,suf,suf,suf);
       printf("                                       // determines %s\n", suf);
       printf("\n");
       printf("        float  %s_bot(%s) ;\n", suf, suf );
       printf("               %s_bot:long_name = \"bottom level of ",suf);
       printf("%s \" ;\n",level_long_name( l_gb->gblevel[i].level_flg));
       printf("               %s_bot:units = \"%s\" ;\n", suf, unit );
       printf("\n");
       printf("        float  %s_top(%s) ;\n",suf,suf);
       printf("               %s_top:long_name = \"top level of ",suf); 
       printf("%s \" ;\n",level_long_name( l_gb->gblevel[i].level_flg));
       printf("               %s_top:units = \"%s\" ;\n",suf,unit );
       printf("\n");
    }
    if(meta==-1){
       printf("ERROR");
    }
} /* end for number of levels */

} /* end print_levels */

static void print_parameters( gribdata* l_gb, char* idx1, char* idx2 ){
/* INPUT: gribdata* l_gb */
/* NB! see the parallel in the  print_levels function */
/* This function prints the parameter variables to the standard output */
/* CALLED BY: print_cdl()*/
 int i;
 int j;
 char* at = " at ";

 for( i=0; i< l_gb->nparam; i++ ){
     int parm = l_gb->param[i].param ;
     char* param= grib_pname(l_gb->param[i].param);
/* used for debug
printf("parameter number = %d\n", parm );
*/
     for( j=0; j< l_gb->param[i].nlevel; j++ ){
	 char* datatype="float";
	 int level = l_gb->param[i].level_flg[j] ;
	 char* suf=levelsuffix(l_gb->param[i].level_flg[j]);

	 /* level == 1 (surface)*/
         if( meta_level(l_gb->param[i].level_flg[j]) == 1 ) {
	    switch( l_gb->param[i].level_flg[j] ){
	      case LEVEL_SURFACE:
	      {
		if( sfcparam( parm )) {
	            printf("%13s  %s(record,%s,%s) ;\n", 
			datatype, param, idx1, idx2 );
		} else if( parm == 7 && level == 1 ) { /* Geopotential height */
	            printf("%13s  %s_%s(%s,%s) ;\n", 
			datatype, param, suf, idx1, idx2 );
		} else {
	            printf("%13s  %s_%s(record,%s,%s) ;\n", 
			datatype, param, suf, idx1, idx2 );
		}
	      }
	      break;
	      case LEVEL_MEAN_SEA:
	      {
		if( mslparam( parm )) {
	            printf("%13s  %s(record,%s,%s) ;\n", 
			datatype, param, idx1, idx2 );
		} else {
	            printf("%13s  %s_%s(record,%s,%s) ;\n", 
			datatype, param, suf, idx1, idx2 );
		}
	      }
	      break;
	      default:
	      {  
	         printf("%13s  %s_%s(record,%s,%s) ;\n", 
			datatype, param, suf, idx1, idx2 );
	      }
              break;
	    }
         } else { /* level != 1 */
	    switch( l_gb->param[i].level_flg[j] ){
	      case LEVEL_ISOBARIC:
	      {
	        printf("%13s  %s(record,level,%s,%s) ;\n", 
			datatype, param, idx1, idx2 );
	      }
	      break;
	      case LEVEL_PV:
	      {
	        printf("%13s  %s_%s(record,%s,%s) ;\n", 
			datatype, param, suf, idx1, idx2 );
	      }
	      break;
	      case LEVEL_LISO:
	      {
		if( lisoparam( parm )) {
	            printf("%13s  %s(record,%s,%s,%s) ;\n", 
			datatype, param, suf, idx1, idx2 );
		} else {
	            printf("%13s  %s_%s(record,%s,%s,%s) ;\n", 
			datatype, param, suf, suf, idx1, idx2 );
		}
	      }
	      break;
	      default:
	      {  
	        printf("%13s  %s_%s(record,%s,%s,%s) ;\n", 
			datatype, param,suf,suf, idx1, idx2);
	      }
              break;
	    }
	 } /* end of level != 1 */	 

	 switch( l_gb->param[i].level_flg[j] ){
	      case LEVEL_SURFACE:
	      {
		if( sfcparam( parm )) {
	            printf("%13s  %s:long_name = \"%s%s%s\" ;\n"," ", 
			param, parameter_long_name( l_gb->param[i].param ),
			at, level_long_name( level) );
	            printf("%13s  %s:GRIB_parameter_number = %d ;\n"," ", 
			param, parm);
	            printf("%13s  %s:GRIB_level_flag = %d ;\n"," ", 
			param, level);
	            printf("%13s  %s:units = \"%s\" ;\n"," ", 
			param, grib_units(l_gb->param[i].param));
	            printf("%13s  %s:_FillValue = -9999.f ;\n"," ", param);
	            printf("%13s  %s:navigation = \"nav\" ;\n"," ", param);
		} else {
	            printf("%13s  %s_%s:long_name = \"%s%s%s\" ;\n"," ", 
			param, suf, parameter_long_name( l_gb->param[i].param),
			at, level_long_name( level) );
	            printf("%13s  %s_%s:GRIB_parameter_number = %d ;\n"," ", 
			param, suf, parm);
	            printf("%13s  %s_%s:GRIB_level_flag = %d ;\n"," ", 
			param, suf, level);
	            printf("%13s  %s_%s:units = \"%s\" ;\n"," ", 
			param, suf, grib_units(l_gb->param[i].param));
	            printf("%13s  %s_%s:_FillValue = -9999.f ;\n"," ", 
			param, suf);
	            printf("%13s  %s_%s:navigation = \"nav\" ;\n"," ", 
			param, suf );
		}
	      }
	      break;
	      case LEVEL_MEAN_SEA:
	      {
		if( mslparam( parm )) {
	            printf("%13s  %s:long_name = \"%s%s%s\" ;\n"," ", 
			param, parameter_long_name( l_gb->param[i].param ),
			at, level_long_name( level) );
	            printf("%13s  %s:GRIB_parameter_number = %d ;\n"," ", 
			param, parm);
	            printf("%13s  %s:GRIB_level_flag = %d ;\n"," ", 
			param, level);
	            printf("%13s  %s:units = \"%s\" ;\n"," ", 
			param, grib_units(l_gb->param[i].param));
	            printf("%13s  %s:_FillValue = -9999.f ;\n"," ", param);
	            printf("%13s  %s:navigation = \"nav\" ;\n"," ", param);
		} else {
	            printf("%13s  %s_%s:long_name = \"%s%s%s\" ;\n"," ", 
			param, suf, parameter_long_name( l_gb->param[i].param),
			at, level_long_name( level) );
	            printf("%13s  %s_%s:GRIB_parameter_number = %d ;\n"," ", 
			param, suf, parm);
	            printf("%13s  %s_%s:GRIB_level_flag = %d ;\n"," ", 
			param, suf, level);
	            printf("%13s  %s_%s:units = \"%s\" ;\n"," ", 
			param, suf, grib_units(l_gb->param[i].param));
	            printf("%13s  %s_%s:_FillValue = -9999.f ;\n"," ", 
			param, suf);
	            printf("%13s  %s_%s:navigation = \"nav\" ;\n"," ", 
			param, suf );
		}
	      }
	      break;
	      case LEVEL_ISOBARIC:
	      {
	        printf("%13s  %s:long_name = \"%s%s%s\" ;\n"," ", 
			param, parameter_long_name( l_gb->param[i].param ),
			at, level_long_name( level) );
	        printf("%13s  %s:GRIB_parameter_number = %d ;\n"," ", 
			param, parm);
	        printf("%13s  %s:GRIB_level_flag = %d ;\n"," ", 
			param, level);
	        printf("%13s  %s:units = \"%s\" ;\n"," ", 
			param, grib_units(l_gb->param[i].param));
	        printf("%13s  %s:_FillValue = -9999.f ;\n"," ", param);
	        printf("%13s  %s:navigation = \"nav\" ;\n"," ", param);
	      }
	      break;
	      case LEVEL_LISO:
	      {
		if( lisoparam( parm )) {
	            printf("%13s  %s:long_name = \"%s%s%s\" ;\n"," ", 
			param, parameter_long_name( l_gb->param[i].param ),
			at, level_long_name( level) );
	            printf("%13s  %s:GRIB_parameter_number = %d ;\n"," ", 
			param, parm);
	            printf("%13s  %s:GRIB_level_flag = %d ;\n"," ", 
			param, level);
	            printf("%13s  %s:units = \"%s\" ;\n"," ", 
			param, grib_units(l_gb->param[i].param));
	            printf("%13s  %s:_FillValue = -9999.f ;\n"," ", param);
	            printf("%13s  %s:navigation = \"nav\" ;\n"," ", param);
		} else {
	            printf("%13s  %s_%s:long_name = \"%s%s%s\" ;\n"," ", 
			param, suf, parameter_long_name( l_gb->param[i].param),
			at, level_long_name( level) );
	            printf("%13s  %s_%s:GRIB_parameter_number = %d ;\n"," ", 
			param, suf, parm);
	            printf("%13s  %s_%s:GRIB_level_flag = %d ;\n"," ", 
			param, suf, level);
	            printf("%13s  %s_%s:units = \"%s\" ;\n"," ", 
			param, suf, grib_units(l_gb->param[i].param));
	            printf("%13s  %s_%s:_FillValue = -9999.f ;\n"," ", 
			param, suf);
	            printf("%13s  %s_%s:navigation = \"nav\" ;\n"," ", 
			param, suf );
		}
	      }
	      break;
	      default:
	      {
	         printf("%13s  %s_%s:long_name = \"%s%s%s\" ;\n"," ", 
			param, suf, parameter_long_name( l_gb->param[i].param),
			at, level_long_name( level) );
	         printf("%13s  %s_%s:GRIB_parameter_number = %d ;\n"," ", 
			param, suf, parm);
	         printf("%13s  %s_%s:GRIB_level_flag = %d ;\n"," ", 
			param, suf, level);
	         printf("%13s  %s_%s:units = \"%s\" ;\n"," ", 
			param, suf, grib_units(l_gb->param[i].param));
	         printf("%13s  %s_%s:_FillValue = -9999.f ;\n"," ", 
			param, suf);
	         printf("%13s  %s_%s:navigation = \"nav\" ;\n"," ", 
			param, suf );
	      }
	      break;
        } 
	printf("\n");
     }/* end for j   */
  } /* end for gb->nparam */
} /* end print_parameters */

static void init_gribdata( gribdata* gb, int nr_gdes ){
/* This function innitialises a struct of type gribdata */
   /* printf("sizeof(char*) = %d", sizeof(char*)); */
   gb->nlgrib=0;
   gb->llgrib=(lgrib*)malloc(sizeof(lgrib)*nr_gdes);
   gb->nmodel=0;
   gb->model=(int*)malloc(sizeof(int)*NR);
   gb->ngrid=0;
   gb->grid=(int*)malloc(sizeof(int)*NR);
   gb->nparam_vers=0;
   gb->param_vers=(int*)malloc(sizeof(int)*NR);
   gb->nparam=0;
   gb->param=(parameter*)malloc(sizeof(parameter)*NR);
   gb->nlevel=0;
   gb->gblevel=(level*)malloc(sizeof(level)*NR);
   gb->nreftime=0;
   gb->creftime=(char**)malloc(sizeof(char*)*NR);
   gb->ntunit=0;
   gb->tunit=(int*)malloc(sizeof(int)*NR);
   gb->ntr_flg=0;
   gb->tr_flg=(int*)malloc(sizeof(int)*NR);
   gb->ntr=0;
   gb->tr=(int*)malloc(sizeof(int)*NR);
   gb->nhas_gds=0;
   gb->has_gds=(int*)malloc(sizeof(int)*NR);
   gb->nnpts=0;
   gb->npts=(int*)malloc(sizeof(int)*NR);
   gb->ngd=0;
   gb->gd=(gdes*)malloc(sizeof(gdes)*nr_gdes);
   return;
}

static void
add_grib(gp, prec, ngdes )
     struct product_data *gp;
     int prec;			/* precision to use for printing float vals */
     int ngdes;                 /* is one if one orders the data by gdes */ 
{
/* This function is originally based on the function print_grib, but instead 
   of printing the values to the terminal screen they are stored in the object 
   struct static gribdata m_gb or m_lgb.gb if not identical values already 
   are stored in the struct. The INPUT parameter ngdes control whether the 
   data is ordered in a structure after gdes type(m_lgb.gb) or not( m_gb). 
*/
    char* ps;
    int b=0;
    int i=0;
    int j=0;
    int index=0;
    gribdata* l_gb=NULL;

    float lflevel[2];
    lflevel[0]= -1;
    lflevel[1]= -1;

    if(ngdes){
	int index_lgb=0;
	b=0;
	
	for( i=0; i< m_lgb.ngb; i++ ){
	    /*m_lgb.gb[i].gd[0].keep=0; */
	    if(gdes_cmp( gp->gd, m_lgb.gb[i].gd )){
		b=1; index_lgb=i; break;
	    }
	   /* m_lgb.gb[i].gd[0].keep=1;*/
	}

	if(!b){
	    if( !(m_lgb.ngb%NR_GDES) && m_lgb.ngb !=0 ){
		/* reallocate more memory*/	        
		if(!(m_lgb.gb=(gribdata*)realloc(m_lgb.gb,
			    (sizeof(gribdata)*( m_lgb.ngb + NR_GDES)))))
		    uerror("Error: can't alloc for parameter array: %d", 
			m_lgb.ngb );
	    }/* end reallocate more memory*/

	    init_gribdata( &m_lgb.gb[m_lgb.ngb],1 );
	    index_lgb = m_lgb.ngb;
	    m_lgb.gb[m_lgb.ngb].ngd=1;
	    gp->gd->keep=1;
	    /* This protects the memory against free_gdes called by the 
			      function free_product_data*/
	    m_lgb.gb[m_lgb.ngb].gd[0] = *(gp->gd); 
	    m_lgb.ngb++;	
	}

	l_gb = &m_lgb.gb[index_lgb];

    }
    else{
	 b=0;

	 for( i=0; i< m_gb.ngd; i++ ){
	     /*m_gb.gd[i].keep=0;*/
	     if(gdes_cmp( gp->gd, &(m_gb.gd[i]) )){
		 b=1;break;
	     }
	     /*m_gb.gd[i].keep=1;*/
	 }

	 if(!b){
	    if( !(m_gb.ngd%NR_GDES) && m_gb.ngd !=0 ){
		/* reallocate more memory*/	        
		if(!(m_gb.gd=(gdes*)realloc(m_gb.gd,
			    (sizeof(gribdata)*( m_gb.ngd + NR_GDES)))))
		    uerror("Error: can't alloc for parameter array: %d", 
			m_gb.ngd );
	    }/* end reallocate more memory*/
	     
	    gp->gd->keep=1;
		/* This protects the memory against free_gdes called by the 
			       function free_product_data*/
	     /*copy_gdes( gp->gd, &(m_gb.gd[i]) );*/
	     m_gb.gd[m_gb.ngd] = *(gp->gd);
	     m_gb.ngd++;
	 }
	 l_gb = &m_gb;
    }
   
    if( !(l_gb->nlgrib%NR_GDES) ){	        
	if(!(l_gb->llgrib=(lgrib*)realloc(l_gb->llgrib,
			    (sizeof(char*)*(l_gb->nlgrib + NR_GDES)))))
        uerror("Error: can't alloc for creftime array: %d", l_gb->nlgrib);
    }
    b=0;
    for( i=0; i< l_gb->nreftime; i++ ){
	if( (!strcmp(l_gb->llgrib[i].delim,gp->delim)) && 
		l_gb->llgrib[i].edition==gp->edition &&
	    l_gb->llgrib[i].center==gp->center && 
		l_gb->llgrib[i].subcenter==gp->subcenter ){
	    b=1;break;
	}
    }
    if(!b){
	strcpy(l_gb->llgrib[l_gb->nlgrib].delim,gp->delim);
	l_gb->llgrib[l_gb->nlgrib].edition=gp->edition;
	l_gb->llgrib[l_gb->nlgrib].center=gp->center;
	l_gb->llgrib[l_gb->nlgrib].subcenter=gp->subcenter;
	l_gb->nlgrib++;
    }
    l_gb->nmodel=AddInt( l_gb->model, l_gb->nmodel, gp->model);
    l_gb->ngrid=AddInt( l_gb->grid, l_gb->ngrid, gp->grid);
    l_gb->nnpts=AddInt( l_gb->npts, l_gb->nnpts, gp->npts);
    if( !(l_gb->nparam%NR) && l_gb->nparam !=0 ){	        
	if(!(l_gb->param=(parameter*)realloc(l_gb->param,
			    (sizeof(parameter)*(l_gb->nparam + NR)))))
        uerror("Error: can't alloc for parameter array: %d", l_gb->nparam);
    }

    l_gb->param[l_gb->nparam].param=gp->param;
    
    b=0;
    index=0;
    for( i=0; i< l_gb->nparam; i++ ){
	if( l_gb->param[l_gb->nparam].param == l_gb->param[i].param ){
	   b=1;
	   index=i;
	   break;
	}
    }

    if(b){
	int b_parlev=0; 
	/* here this variable means boolean control for the level */ 
        for( i=0; i< l_gb->param[index].nlevel; i++ ){
	    if( gp->level_flg == l_gb->param[index].level_flg[i] ){
		b_parlev=1;break;
	    }
	}
	if(!b_parlev){
	    l_gb->param[index].level_flg[l_gb->param[index].nlevel]=gp->level_flg;
	    l_gb->param[index].nlevel++;
	}
	    
    }
    else{ /*if(!b)*/
	
	for( j=0; j<NR; j++ ){  
	    l_gb->param[l_gb->nparam].level_flg[j]=-1;/* initialise values */
	}
	l_gb->param[l_gb->nparam].nlevel = 1; 
	/* the new parameter has at least one level */
        l_gb->param[l_gb->nparam].level_flg[0] = gp->level_flg;
	l_gb->nparam++;
    } 
   
   if( !(l_gb->nlevel%NR) && l_gb->nlevel!=0 ){	        
	if(!(l_gb->gblevel=(level*)realloc(l_gb->gblevel,
			    (sizeof(level)*(l_gb->nlevel + NR)))))
        uerror("Error: can't alloc for level array: %d", l_gb->nlevel);
    } 

    l_gb->gblevel[l_gb->nlevel].level_flg = gp->level_flg;            
    
    b=0;
    index=0;
    /* index is the index of position in the array 'struct level gblevel[NR]' */
    for( i=0; i< l_gb->nlevel; i++ ){
	if( l_gb->gblevel[l_gb->nlevel].level_flg == l_gb->gblevel[i].level_flg ){
	   b=1;
	   index=i;
	   break;
	}
    }

    if(!b){
        l_gb->gblevel[l_gb->nlevel].level_flg=gp->level_flg;
        l_gb->gblevel[l_gb->nlevel].clevelname = levelname(gp->level_flg);
	l_gb->gblevel[l_gb->nlevel].clevelsuffix = levelsuffix(gp->level_flg);
        index=l_gb->nlevel;
        l_gb->gblevel[index].nflevel=0; /* initialise values */
	for( j=0; j<NR; j++ ){  /* initialise values */
	    l_gb->gblevel[index].flevel[0][j]=-1;
	    l_gb->gblevel[index].flevel[1][j]=-1;
	    l_gb->gblevel[index].level[0][j]=-1;
	    l_gb->gblevel[index].level[1][j]=-1;
	}
	    
        l_gb->nlevel++;
    }    

    switch(meta_level(gp->level_flg)) {
    case 1:
	break;
    case 2:      
	lflevel[0] = level1(gp->level_flg,gp->level);
	/* lflevel[1] = -1; */
    case 3:
	lflevel[0] = level1( gp->level_flg, gp->level );
        lflevel[1] = level2( gp->level_flg, gp->level );
	break;
    }
      
    b=0;
    for( i=0; i< l_gb->gblevel[index].nflevel; i++ ){
	if( ( lflevel[0] == l_gb->gblevel[index].flevel[0][i] ) && 
            ( lflevel[1] == l_gb->gblevel[index].flevel[1][i] ) &&
	    ( gp->level[0] == l_gb->gblevel[index].level[0][i] ) && 
            ( gp->level[1] == l_gb->gblevel[index].level[1][i] )
	    ){
	    b=1;break;
	}
    }

    if(!b){ 
	l_gb->gblevel[index].flevel[0][l_gb->gblevel[index].nflevel] = 
		lflevel[0];
        l_gb->gblevel[index].flevel[1][l_gb->gblevel[index].nflevel] = 
		lflevel[1];
        l_gb->gblevel[index].level[0][ l_gb->gblevel[index].nflevel] = 
		gp->level[0];
	l_gb->gblevel[index].level[1][ l_gb->gblevel[index].nflevel] = 
		gp->level[1];
   	l_gb->gblevel[index].nflevel++;
    }

    ps=(char*)malloc(80);
    sprintf( ps," %04d/%02d/%02d:%02d:%02d", 
				/* century 21 doesn't start until 2001 */
	   gp->year+(gp->century -1)*100,
	   gp->month, gp->day, gp->hour, 
	   gp->minute);


    if( !(l_gb->nreftime%NR) ){	        
	if(!(l_gb->creftime=(char**)realloc(l_gb->creftime,
			    (sizeof(char*)*(l_gb->nreftime + NR)))))
        uerror("Error: can't alloc for creftime array: %d", l_gb->nreftime);

	for( i=l_gb->nreftime; i<(l_gb->nreftime + NR); i++ )
		l_gb->creftime[i]=NULL;
    }

    l_gb->creftime[l_gb->nreftime]=ps;
    
    b=0;
    for( i=0; i< l_gb->nreftime; i++ ){
	if(!strcmp( l_gb->creftime[l_gb->nreftime], l_gb->creftime[i] )){
	   b=1;
	}
    }
    if(!b){
	l_gb->nreftime++;
    }

    l_gb->ntr_flg= AddInt( l_gb->tr_flg, l_gb->ntr_flg, gp->tr_flg );
/*
    Following line assumes (gp->tunit == TUNITS_HOUR) NCEP models are all hours
*/
    l_gb->ntr = AddInt( l_gb->tr, l_gb->ntr, frcst_time( gp ) );

    l_gb->nhas_gds= AddInt( l_gb->has_gds, l_gb->nhas_gds, gp->has_gds );
    l_gb->nparam_vers= AddInt( l_gb->param_vers, l_gb->nparam_vers,gp->param_vers);
} /* end of add_grib */ 

static int AddInt( int* arr, int narr, int val ){
/* This is an utility function local to this file used by the function 
   add_grib. */ 
       
    int i=0;
    int b=0;

    if( !(narr%NR) && narr !=0 ){	        
	if(!(arr=(int*)realloc(arr,(sizeof(int)*(narr + NR)))))
        uerror("Error: can't alloc for creftime array: %d", narr);

	for( i=narr; i<(narr + NR); i++ )
		arr[i]=0;
    }

    
    for( i=0; i< narr; i++ ){
	if( val == arr[i] ){
	   b=1; break;
	}
    }

    if(!b){
	arr[narr]= val;
	narr++;
    }

  return narr;
} /* end of AddInt */ 

static int AddCharp( char** arr, int narr, char* val ){
/* This is an utility function local to this file used by the function 
   add_grib. */ 
       
    int i=0;
    int b=0;

    if( !(narr%NR) && narr !=0 ){	        
	if(!(arr=(char**)realloc(arr,(sizeof(char*)*(narr + NR)))))
        uerror("Error: can't alloc for creftime array: %d", narr);

	for( i=narr; i<(narr + NR); i++ )
		arr[i]=NULL;
    }

    for( i=0; i< narr; i++ ){
	if(!strcmp( arr[i], val )){
	   b=1; break;
	}
    }

    if(!b){
	arr[narr]= val;
	narr++;
    }

  return narr; 
} /* end of AddCharp */

/*
 * Parse raw product bytes into product_data structure.  Returns 0 if
 * failed.  User should call free_product_data() on result when done with
 * it.  Also expands quasi-regular grids to full rectangular grids if qmeth
 * is non-null.
 */
static product_data *
grib_decode(prodp, quasp)
     prod *prodp ;		/* input raw product bytes */
     quas *quasp ;		/* if non-null, method used to expand
				   quasi-regular "grids" */
{
    grib1 *gp = new_grib1(prodp); /* overlay raw bits on raw grib1 structure */
    product_data *pdp;

    if (gp == 0)
	return 0;

    pdp = new_product_data(gp);	/* compute cooked product structure, with GDS
				   (manufactured, if necessary) and bytemap */
    free_grib1(gp);

    if (!pdp)
	return 0;

    if (pdp->gd->quasi && quasp) {
	int ret = expand_quasi(quasp, pdp) ; /* Changes *pdp */
	if (!ret) {
	    uerror("%s: can't expand quasi-regular grid", pdp->header);
	    return 0;
	}
    }
    return pdp;
} /* end of grib_decode */

static int
do_add (fp, ep, timeout, prec, quasp)
/* This function is closely based on the function do_dump in the file 
   gribdump.c and gets, decodes and stores the header information from the 
   gribobject 
*/
    FILE *fp;			/* input */
    FILE *ep;			/* if non-null, where to append bad GRIBs */
    int timeout;
    int prec;			/* precision to use for printing floats */
    quas *quasp;		/* if non-null, quasi-reg expansion method */
{
    struct prod the_prod;
    struct product_data *gribp;	/* UPC structure */
    
    num_wmo_messages = 0;
    num_gribs_decoded = 0;

    if (brief_output)
	printf(" cnt mdl grd prm    lvlf lev1 lev2  trf tr0 tr1 bits bms gds   npts header\n");

    TEST_NOSERC("tralla",1000);

    while(1){
        int bytes= add_get_prod(fp, timeout, &the_prod);

        if( bytes == -1)
	    return -1;

	if (bytes == 0){
	  TEST_NOSERC( " bytes == 0", num_gribs_written );
	  continue;/*continue;*/
        }
	else
	  num_wmo_messages++;
	
	/* decode it into a grib1 product structure */
	gribp = grib_decode(&the_prod, quasp);
	if (gribp == 0) {
	    TEST_NOSERC( "num_gribs_written gribp == 0", num_gribs_written );
	    if (ep) {
		if (fwrite(the_prod.bytes, the_prod.len, 1, ep) == 0) {
		    serror("writing bad GRIB to error file");
		}
	    }
	    continue;/*continue;*/
	}
	num_gribs_decoded++;

	if( gdes_analyse )
	    add_grib(gribp,prec, 1);
	else
	    add_grib(gribp, prec, 0);
	free_product_data(gribp);
	num_gribs_written++;

	TEST_NOSERC("num_gribs_written", num_gribs_written );
    }
} /* do_add */

void rtll ( tlm0d, tph0d, tlmd, tphd, almd, aphd ) 
double tlm0d, tph0d, tlmd, tphd, *almd, *aphd ;
{

/* function [almd,aphd]=rtll(tlm0d,tph0d,tlmd,tphd);
 *-------------------------------------------------------------------------
 * RTLL transforms rotated coordinates (tlmd,tphd) into
 * ordinary geographic coordinates (almd,aphd). i/o decimal degrees (DD)
 * tlm0d, tph0d: lon & lat of the center of rotation ("North Pole") in DD.
*/
	double PI, dtr, ctph0, stph0, stph, ctph, ctlm, stlm, aph, cph ;
	PI = atan2( 1, 1 ) * 4 ;
	dtr = PI / 180 ;

	ctph0 = cos( tph0d * dtr );
	stph0 = sin( tph0d * dtr );

	stph = sin( tphd * dtr );
	ctph = cos( tphd * dtr );

	ctlm = cos( tlmd * dtr );
	stlm = sin( tlmd * dtr );

	aph = asin( stph0 * ctph * ctlm + ctph0 * stph );
	cph = cos( aph );

	*almd = tlm0d + asin( stlm * ctph / cph ) / dtr ;
	*aphd = aph / dtr ;

}

static void print_cdl(){
/* 
   This function prints the cdl to standard output (the screen),
   the form of the cdl is based on an analysis of the grib code
   that stores data in the gribdata object m_gb 
*/ 

 int i=0;
 int j=0;
 int k=0;
 float ilon;
 float jlat;
 gribdata* l_gb=NULL;
 cdl_gdes_data* g_data=NULL;
 gdes* lgdes=NULL;
 int ngrids;
 char* idx1 ;
 char* idx2 ;
 time_t t;
 char time_now[35];

 l_gb = &m_gb;
 ngrids=m_gb.ngd;
 lgdes=&(l_gb->gd[0]);
 g_data= get_gdes_data(lgdes);

/* start printing output cdl */
 printf("netcdf %s{ \n", "Replace_with_model_name" );
 printf("\n");
 printf("\n");
 printf("dimensions:\n");
 printf("        record = UNLIMITED ;   // (reference time, forecast time)\n");
 switch (lgdes->type) {
    case GRID_POLARS:
    {
	gdes_polars *gg = &lgdes->grid.polars ;
	printf("        x = %d ;\n",gg->nx);
	printf("        y = %d ;\n",gg->ny);
	idx1 = "y" ;
	idx2 = "x" ;
    }
    break;
    case GRID_LAMBERT:
    {
	gdes_lambert *gg = &lgdes->grid.lambert;
	printf("        x = %d ;\n",gg->nx);
	printf("        y = %d ;\n",gg->ny);
	idx1 = "y" ;
	idx2 = "x" ;
    }
    break;
    case GRID_GAU:
    {
	gdes_gau *gg = &lgdes->grid.gau;
	printf("        lat = %d ;            // latitude\n",gg->nj);
	printf("        lon = %d ;            // longitude\n",gg->ni); 
	/* kludge cuz gg->n gets overwritten, needed for writing out lat later*/
	g_data->dj = gg->n; 
	idx1 = "lat" ;
	idx2 = "lon" ;
    }
    break;
    default:
    {
	printf("        lat = %d ;            // latitude\n",g_data->nj);
	printf("        lon = %d ;            // longitude\n",g_data->ni);
	idx1 = "lat" ;
	idx2 = "lon" ;
    }
    break;
}

/* print dimensions */
 for( i=0; i<l_gb->nlevel; i++ ){
    int meta = meta_level(l_gb->gblevel[i].level_flg);

    if(meta==2){
	switch(l_gb->gblevel[i].level_flg){
	   case LEVEL_ISOBARIC:
	   {
	       printf("        level = %d ;           // isobaric levels\n", 
			l_gb->gblevel[i].nflevel);
	   }
	   break;
	   default:
	   {
	       char* suf=levelsuffix(l_gb->gblevel[i].level_flg);
	       printf("        %s = %d ;              // %s\n", 
		suf,l_gb->gblevel[i].nflevel, 
		level_long_name( l_gb->gblevel[i].level_flg) );
	   }
	}
    }
    if(meta==3){ 
       char* suf=levelsuffix(l_gb->gblevel[i].level_flg);
       printf("        %s = %d ;              // %s\n", suf,
	l_gb->gblevel[i].nflevel, level_long_name( l_gb->gblevel[i].level_flg));
    }
    if(meta==-1){
       printf("ERROR");
    }
 }

 printf("        time_len = 21 ;        // string length for datetime strings\n");
 printf("        valtime_offset = %d ;   // number of offset times\n",
	l_gb->ntr );
 printf("        nmodels = %d ;          // number of models\n", l_gb->nmodel);
 printf("        ngrids = %d ;           // number of grids\n", ngrids );
 printf("        nav = 1 ;              // for navigation\n");
 printf("        nav_len = 100 ;        // max string length for navigation strings\n");
 printf("\n");
 printf("variables:\n");
 printf("\n");

 print_time_cdl();

 print_levels(l_gb);

/* print lat/lon info is needed */
 printf("\n");
 if(  ! strcmp( idx1, "lat" ) ) {
    printf("        // The following lat and lon coordinate variables are redundant,\n");
    printf("        // since the navigation variables provide the necessary information.\n");
    printf("        // The extra information is included here for human readability.\n");
    printf("\n");
    printf("        float  lat(lat) ;\n");
    printf("               lat:long_name = \"latitude\" ;\n");
    printf("               lat:units = \"degrees_north\" ;\n");
    printf("\n");
    printf("        float  lon(lon) ;\n");
    printf("               lon:long_name = \"longitude\" ;\n");
    printf("               lon:units = \"degrees_east\" ;\n");
    printf("\n");
}

 printf("        long   model_id(nmodels) ;\n");
 printf("               model_id:long_name = \"generating process ID number\" ;\n");
 printf("\n");
 
 print_netcdf_gdes_nav(lgdes);

 print_parameters(l_gb, idx1, idx2 );

 printf("\n");
 printf("// global attributes\n");
 /* updated traceability 18.12.2001 by ab */
 t = time(NULL);
 strftime(time_now,35,"%Y-%m-%d %T",localtime(&t));
 printf("               :history = \"%s - created by %s\" ; \n",
	time_now, "gribtocdl" );
 printf("               :title = \"%s\" ;\n", "Enter model definition here" );
 printf("               :Conventions = \"NUWG\" ;\n");
 printf("               :GRIB_reference = \"Office Note 388 GRIB\" ;\n");
 printf("               :GRIB_URL = \"http://www.nco.ncep.noaa.gov/pmb/docs/on388/\" ;\n");
 printf("               :version = 0.0 ;\n");
 printf("\n");
 printf("data:\n");
 printf("\n"); 

 /* write out levels in proper order */
 for( i=0; i<l_gb->nlevel; i++ ){
    int meta = meta_level(l_gb->gblevel[i].level_flg);

    if(meta==2){
	switch(l_gb->gblevel[i].level_flg){
	case LEVEL_ISOBARIC:
	{
	    float a = -9999 ;
	    float b = -9999 ;
	    int c = 1 ;
	    k=0;
	    printf(" level = ");    
	    /* sort levels, starting with highest */
	    for( j=0; j<l_gb->gblevel[i].nflevel ; j++ ){
		if( l_gb->gblevel[i].flevel[0][j] > a ) 
			a = l_gb->gblevel[i].flevel[0][j] ;
	    }
	    printf("%.1f", a);
	    while( c < l_gb->gblevel[i].nflevel ) {
		printf(", " );
	        if(!(++k%10)){ printf("\n      ");}
		for( j=0; j<l_gb->gblevel[i].nflevel ; j++ ){
		    if( l_gb->gblevel[i].flevel[0][j] < a && 
		    	 b < l_gb->gblevel[i].flevel[0][j] ) 
			b = l_gb->gblevel[i].flevel[0][j] ;
		}
		printf("%.1f", b );
		a = b;
		b = -9999 ;
		c++ ;
	    }
	    printf(" ;\n") ;
	}
	break;
	default:
	{
	    float a = 999999 ;
	    float b = 999999 ;
	    int c = 1 ;
	    char* suf=levelsuffix(l_gb->gblevel[i].level_flg);
	    printf(" %s = ", suf );   
	    k=0;

	    /* sort levels, starting with lowest */
	    for( j=0; j<l_gb->gblevel[i].nflevel ; j++ ){
		if( l_gb->gblevel[i].flevel[0][j] < a ) 
			a = l_gb->gblevel[i].flevel[0][j] ;
	    }
	    printf("%.1f", a);
	    while( c < l_gb->gblevel[i].nflevel ) {
		printf(", " );
	        if(!(++k%10)){ printf("\n      ");}
		for( j=0; j<l_gb->gblevel[i].nflevel ; j++ ){
		    if( l_gb->gblevel[i].flevel[0][j] > a && 
		    	 b > l_gb->gblevel[i].flevel[0][j] ) 
			b = l_gb->gblevel[i].flevel[0][j] ;
		}
		printf("%.1f", b );
		a = b;
		b = 999999 ;
		c++ ;
	    }
	    printf(" ;\n") ;
	}
	break;
	}
    }
    if(meta==3){
	switch(l_gb->gblevel[i].level_flg){
	default:
	{
	char* suf=levelsuffix(l_gb->gblevel[i].level_flg);
        printf(" %s_bot = ", suf );
	k=0;
	for( j=0; j<l_gb->gblevel[i].nflevel -1; j++ ){    
	       printf("%.1f, ", l_gb->gblevel[i].flevel[1][j]);
	       if(!(++k%8)){ printf("\n");printf("      ");}
	}
	printf("%.1f ;\n", l_gb->gblevel[i].flevel[1][l_gb->gblevel[i].nflevel-1]);
	k=0;
	printf(" %s_top = ", suf );
	for( j=0; j<l_gb->gblevel[i].nflevel - 1; j++ ){  
	       printf("%.1f, ", l_gb->gblevel[i].flevel[0][j]);
	       if(!(++k%8)){ printf("\n");printf("      ");}
	}
        printf("%.1f ;\n", l_gb->gblevel[i].flevel[0][l_gb->gblevel[i].nflevel-1]);            
	}
	break;
	} /* end switch on level_flg */
    } /* end meta == 3 */  
    if(meta==-1){
       printf("ERROR");
    }
  }/* end for( i=0; i<l_gb->nlevel; i++ ) */

 /* Updated 18.12.2001 by ab - print all model numbers */
 printf(" model_id = ");
 for( i=0; i<l_gb->nmodel-1; i++ ) {
    printf("%d, ", l_gb->model[i]);
 }
 printf("%d ;\n",l_gb->model[l_gb->nmodel-1]);

 /* write out offsets in ascending order */
 {
 int l, h, val;

 printf(" valtime_offset = "); 
 l = -1;
 for( j=0; j<l_gb->ntr; j++ ){
 	h = 9999;
 	for( i=0; i<l_gb->ntr ; i++ ){
    		if( l_gb->tr[ i ] > l && l_gb->tr[ i ] < h ) {
			val = l_gb->tr[ i ];
			h = l_gb->tr[ i ];
		}
 	}
	l = val;
	if( j == (l_gb->ntr -1 )){
    		printf( "%d ;\n\n", val );
	} else {
    		printf( "%d, ", val );
	}
 }
 }
 
 printf("\n");
 printf(" // Navigation\n");
 printf(" nav_model = \"GRIB1\" ;\n");
 printf(" grid_type_code = %d ;\n",l_gb->gd[0].type);
 printf(" grid_type = \"%s\" ;\n", gds_typename(l_gb->gd[0].type));
 printf(" grid_name = \" \" ;\n");
 if(l_gb->nlgrib>1) printf("ERROR, klarer ikke nlgrib>1");
 printf(" grid_center = %d ;\n",l_gb->llgrib[0].center);
 /* updated 18.12.2001 by ab - grid numbers are retrieved */
 printf(" grid_number = ");
 i = 0;
 while( i < ngrids ) {
    if( i > 0 )
	 printf(", ") ;
    printf("%d",l_gb->grid[i]);
    i++ ;
 }
 printf(" ;\n") ;

/*
 printf(" %s = %o ;\n", "ScanMode", l_gb->gd[0].scan_mode );
 printf(" %s = %o ;\n", "ScanMode", (l_gb->gd[0].scan_mode& 0x40) == 0 ? 1: 0 );
*/
 print_netcdf_gdes_type(&(l_gb->gd[0]));

 printf(" %s = %d ;\n", "ResCompFlag", l_gb->gd[0].res_flags );
 printf("\n");

 /* exit here before writing out lat and lon for certain grids*/
 switch (lgdes->type) {
    case GRID_POLARS:
    case GRID_LAMBERT:
    {
	printf("}\n");
	return;
    }
    break;
    default:
    {
	int count ;

	/* updated 18.12.2001 by ab 
           Adjust grid so that it is always South-North [-90:+90]
           and West-East [-180:180] or [0:360] */
        /* adjust numbers so that they are increasing */
	/* Assume longitude is always west-east (sic) */
        while (g_data->lo2 < g_data->lo1) g_data->lo2 += 360;
        /* Just swap latitudes if wrong direction */
	if (g_data->la1 > g_data->la2) {
	  float tlat=g_data->la1;
	  g_data->la1=g_data->la2;
	  g_data->la2=tlat;
	}

	if( lgdes->type == GRID_RLL ) {
		double pole_lat = 90 -32.5 ;
		double pole_lon = 10 ;
		double incrementdj, incrementdi ;
		double lat[ 273 ], lon[ 235 ] ;
		double almd, aphd ;
        	incrementdi = g_data->di / 1000 ;
        	incrementdj = g_data->dj / 1000 ;
		printf("g_data->di=%6.2f\n", g_data->di );  
		printf("g_data->dj=%6.2f\n", g_data->dj );  
		lat[ 0 ] = g_data->la1 ;
		lon[ 0 ] = g_data->lo1 ;
/*
		g_data->la2 = -15.563 ;
		g_data->lo2 = 24.5 ;
*/
		for( j = 1; j <= 272; j++ ) {
			lat[ j ] = lat[ j -1 ] + incrementdj ;
			for( i = 1; i <= 234; i++ ) {
				lon[ i ] = lon[ i -1 ] + incrementdj ;
/*
				printf("lon[ %d ] =%8.6f, lat[ %d ] =%8.6f\n" ,
					 i, lon[ i ], j, lat[ j ]) ;
*/
				rtll( pole_lon, pole_lat, lon[ i ], lat[ j ],
					&almd, &aphd );
				printf(" lon_true[ i, j ]=%10.6f, lat_true[ i, j ] =%10.6f\n", almd, aphd ) ;
/*
				rtll( pole_lon, pole_lat, lon[ i ], lat[ j ] );
				print "lon_true[ i, j ] =lon_true[ i, j ], " ,
			 	"lat_true[ i, j ] =lat_true[ i, j ]\n" ;
*/
			}
		}
		return ;

    	} else if( lgdes->type == GRID_GAU ) {
		/* create increment from g_data->dj which is really gg->n */
		g_data->dj = abs(g_data->la1) / (g_data->dj -1);
    	}

	/* Write longitude grid points, count is number to write */
	count = (int) ( (g_data->lo2 - g_data->lo1) / g_data->di ) + 1.5 ;

        ilon=g_data->lo1;
    
        printf(" lon =");
        k=0;
	for( i = 1; i <  count ; i++ ){
	    printf("%6.2f,",ilon );  
            if(!(++k%8)){ printf("\n");printf("      ");}
	    ilon += g_data->di ;
        }
        printf("%6.2f ;\n",ilon);
        printf("\n");

        /* Write latitude grid points, count is number to write */
    
	count = (int) ( (g_data->la2 - g_data->la1) / g_data->dj ) + 1.5 ;

       	jlat=g_data->la1;
       	printf(" lat =");
       	k=0;
	for( i = 1; i <  count ; i++ ){
    		printf("%6.2f,",jlat );
       		if(!(++k%8)){ printf("\n");printf("      ");}
    		jlat+=  g_data->dj ;
       	}
       	printf("%6.2f ;\n",jlat);

       	printf("\n}\n");

    } /* end of default */
    break;
} /* end of switch */
  return;
} /* end of print_cdl */


int
main(ac,av)
     int ac ;
     char *av[] ;
{
    /* This main function is similar to the main in the gribdump.c */
    int counter=0;
    extern int optind;
    char *logfname = "-" ;	/* log file name, default is stdout */
    FILE *fp = stdin;		/* default input */
    FILE *ep = 0;		/* file handle for bad GRIBS output, when
				   -e badfname used */
    int prec = DEFAULT_PRECISION; /* precision to use for printing floats */
    int timeo = DEFAULT_TIMEOUT ; /* timeout */
    quas *quasp = 0;		/* default, don't expand quasi-regular grids */
    /*
     * register exit handler
     */
    init_gribdata(&m_gb, NR_GDES );
    m_lgb.ngb=0;
    m_lgb.gb=(gribdata*)malloc(sizeof(gribdata)*NR_GDES);
    /* allocate memory for the list of gribdata*/
   
    /*TEST_NOSERC("HVA SKJER?",0 );*/

    if(atexit(cleanup) != 0)
    {
	  serror("atexit") ;
	  exit(1) ;
    }
    
    /*
     * set up signal handlers
     */
    set_sigactions() ;
    
    {
	extern int opterr;
	extern char *optarg;
	int ch;
	int logmask = (LOG_MASK(LOG_ERR) | LOG_MASK(LOG_NOTICE)) ;
	int errflg = 0;
	
	opterr = 1;
	
	brief_output = 0;
        headers_only = 0;
	testNoserc=0;
	gdes_analyse=0;
        /* Default is to generate CDL */
        printcdl=1;
        
        /* New option "-g gribtonctable" implemented May 2002 by Arild Burud*/
	while ((ch = getopt(ac, av, "vxhwyzbl:t:e:p:q:g:")) != EOF) {
	    switch (ch) {
	    case 'v': /* verbose */
		logmask |= LOG_MASK(LOG_INFO) ;
		break;
	    case 'x': /* ? */
		logmask |= LOG_MASK(LOG_DEBUG) ;
		break;
	    case 'g': /* gribtonctable */
	       fprintf(stderr, "Conversion table in file %s\n",optarg) ; 
                user_makeparamtable(optarg);
/*
                user_printparamtable(); 
*/
		break;
	    case 'l': /* logfile */
		logfname = optarg ;
		break;
	    case 't': /* timeout */
		timeo = atoi(optarg) ;
		if(timeo < 1) {
		    fprintf(stderr, "%s: invalid timeout %s",
			    av[0], optarg) ;
		    errflg++;
		}
		break;
	    case 'h': /* headers */
		headers_only = 1;
		break;
	    case 'b': /* brief */
		brief_output = 1;
		headers_only = 1;
		break;
	    case 'e': /* errorfile */
		ep = fopen(optarg, "w");
		if(!ep) {
		    serror("can't open %s", optarg);
		    errflg++;
		}
		break;
	    case 'p': /* precision */
		prec = atoi(optarg) ;
		if(prec < 1 || prec > 99) {
		    fprintf(stderr, "% s: invalid precision %s",
			    av[0], optarg) ;
		    errflg++;
		}
		break;
	    case 'q': /* grid expansion */
		quasp = qmeth_parse(optarg);
		if(!quasp) {
		    fprintf(stderr,
			    "%s: invalid quasi-regular expansion method %s",
			    av[0], optarg) ;
		    errflg++;
		}
		break;
	    case 'z': /* debug */
		testNoserc=1;
		break;
	    case 'y': /* debug */
		gdes_analyse=1;
		break;
	    case 'w': /* produce CDL */
		printcdl=1;
		break; 
	    case '?': /* print help */
		errflg++;
		break;
	    }
	}

	if(errflg)
	  usage(av[0]);	
	(void) setulogmask(logmask) ;
    }
        
    /*
     * initialize logger
     */
    (void) openulog(ubasename(av[0]),
		    (LOG_CONS|LOG_PID), LOG_LOCAL0, logfname) ;
    uinfo("Starting Up") ;

    do {
	if (optind < ac) {	/* arguments still left, use as filenames */
	    fp = fopen(av[optind++], "r");
	    if(!fp) {
		serror("can't open %s", av[optind-1]);
		return 1;
	    }
            grib_file_name=av[optind-1]; /* added 18.12.2001 by a.b. */
	}
	do_add(fp, ep, timeo, prec, quasp);

	if(testNoserc){
	    counter++;
	    printf("   %24s : %d\n","counter", counter ); 
	}  
    } while (optind < ac);

    if( gdes_analyse )
	print_gribdata(1);
    else if(printcdl)
	print_cdl();
    else
	print_gribdata(0);

    
    return 0;
} /* end of main */
