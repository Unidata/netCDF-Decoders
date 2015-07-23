/*
 * Here are methods to include user-defined mapping of
 * Grib codes to NetCDF codes, by reading from a user-supplied
 * file.
 * Written by Arild Burud/NoSerC, May 2002 
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#define MAXLINE 1024    /* Longest line to read from file */

struct user_params {
  char *name;           /* Short netCDF name of parameter */
  char *long_name;      /* Long netCDF name of parameter */
  char *units;          /* Units as understood by udunits(3) */
};

static struct user_params user_param_table[256] = {
  { NULL, NULL, NULL } };

static int user_table_initialized = 0; /* Set to 1 when users table read */

/*
 * Initialize the whole table to NULL
 */

void
user_inittable(void)
{
  int i;
  /* printf("user_inittable entry\n"); */
  for (i=0; i<256; i++) {
    user_param_table[i].name=NULL;
    user_param_table[i].long_name=NULL;
    user_param_table[i].units=NULL;
  }
  user_table_initialized=1;
  /* printf("user_inittable exit\n"); */
}

/*
 * Just a simple method to determine comment lines in userfile
 * A comment line does not start with a digit...
 * Allow white space before digits 
 */

int
user_commentline(line)
     char *line;
{
  char *p=line;
  while (*p!='\0') {
    if (isspace(*p)) p++;
    else {
      if (isdigit(*p)) return 0;
      else return 1;
    }
  }
  return 1; /* blank line == comment line */
}

/*
 * A simple stripper for trailing spaces
 */
char *
user_strip_trailer(str)
     char str[];
{
  int ln;
  ln=strlen(str)-1;
  while (ln>=0) {
    if (isspace(str[ln])) {
      str[ln--]='\0';
    }
    else return str;
  }
  return NULL;
}

/*
 * Open user file and read into table
 */

void
user_makeparamtable(userfile)
     char *userfile;
{
  char line[MAXLINE];
  FILE *fp;
  char *p;
  int gcode;
  char *gname;
  char *glongname;
  char *gunits;
  char *blank=" ";
  int nl=0;
  /* printf("user_makeparamtable(%s) entry\n",userfile); */
  if (user_table_initialized == 0) user_inittable();
  /* open user file */
  fp=fopen(userfile,"r");
  if (!fp) {
    fprintf(stderr,"Error: can't open %s\n",userfile);
    return;
  }
  /* start reading it */
  while (fgets(line,MAXLINE,fp)!=NULL) {
    nl++;
    if (!user_commentline(line)) {
      /* parse line into code, name, longname and units */
      p=strtok(line,";");
      gname=user_strip_trailer(strtok(NULL,";"));
      glongname=user_strip_trailer(strtok(NULL,";"));
      gunits=user_strip_trailer(strtok(NULL,";"));
      /* just check that we have something to store */
      /* 14.8.2002 by AB, now allow none/null units */
      if (p==NULL || gname==NULL || glongname==NULL){
        fprintf(stderr,"Error in %s, line %d not processed correctly\n",userfile,nl);
      }
      else {
	if (gunits==NULL) gunits=blank; /* units=none  (14.8.2002 by AB)*/
	if (   !strcmp(gunits,"none") 
            || !strcmp(gunits,"None") 
            || !strcmp(gunits,"NONE")) gunits=blank; 
	gcode=atoi(p);
	if (gcode<0 || gcode>=256) {
          fprintf(stderr,"Error in %s, invalid parameter %d on line %d\n",userfile,gcode,nl);
          return;
	} 
	user_param_table[gcode].name=strdup(gname);
	user_param_table[gcode].long_name=strdup(glongname);
	user_param_table[gcode].units=strdup(gunits);
      }
    }
  }
  fclose(fp);
  /* printf("user_makeparamtable() exit\n"); */
}

/*
 * Just a little debug printout of the stored table
 */

void
user_printparamtable(void)
{
  int i;
  printf(" // User defined parameter table\n");
  printf(" // Grib number;netCDF name;netCDF long name;Units;\n");
  for (i=0; i<256; i++) {
    if (user_param_table[i].name!=NULL) {
      printf("// %d;%s;%s;%s;\n",
             i,
             user_param_table[i].name,
             user_param_table[i].long_name,
	     user_param_table[i].units);
    }
  }
  printf(" // End of table\n");
}

/*
 * With a given grib code, return netCDF name
 */

char *
user_pname(param)
     int param;
{
  /* printf("user_pname(%d) entry\n",param); */
  if (user_table_initialized == 0) user_inittable();
  if (param < 0 || param >= 256) return NULL;
  return user_param_table[param].name;
}

/*
 * With a given grib code, return netCDF long parameter name
 */

char *
user_plongname(param)
     int param;
{
  /* printf("user_plongname(%d) entry\n",param); */
  if (user_table_initialized == 0) user_inittable();
  if (param < 0 || param >= 256) return NULL;
  return user_param_table[param].long_name;
}

/*
 * With a given netCDF short parameter name, return the Grib code
 * Use linear search in table
 */

int
user_gribpcode(pname)
     char *pname;
{
  int i;
  /* printf("user_gribpcode(%s) entry\n",pname); */
  if (user_table_initialized == 0) user_inittable();
  for (i=0; i<256; i++) {
    if (user_param_table[i].name!=NULL) {
      if (strcmp(pname,user_param_table[i].name)==0)return i;
    }
  }
  return -1; /* not found */
}

/*
 * Given the parameter code, return units string
 */

char *
user_gribunits(param)
     int param;
{
  /* printf("user_gribunits(%d) entry\n",param); */
  if (user_table_initialized == 0) user_inittable();
  if (param < 0 || param >= 256) return NULL;
  return user_param_table[param].units;
}


