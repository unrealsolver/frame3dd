/*
 FRAME3DD:
 Static and dynamic structural analysis of 2D and 3D frames and trusses with
 elastic and geometric stiffness.
 ---------------------------------------------------------------------------
 http://frame3dd.sourceforge.net/
 ---------------------------------------------------------------------------
 Copyright (C) 1992-2009  Henri P. Gavin

 FRAME3DD is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 FRAME3DD is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with FRAME3DD.  If not, see <http://www.gnu.org/licenses/>.
*//**
	@file
	Input/output routines for FRAME.

	@note The file format for FRAME is defined in doc/user_manual.html.
*/

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>

#include "common.h"
#include "frame3dd_io.h"
#include "coordtrans.h"
#include "lu_dcmp.h"
#include "nrutil.h"

/* #define MASSDATA_DEBUG */

/* forward decls */

static void getline_no_comment(
	FILE *fp,    /**< pointer to the file from which to read */
	char *s,     /**< pointer to the string to which to write */
	int lim      /**< the longest anticipated line length  */
);


/*------------------------------------------------------------------------------
PARSE_OPTIONS -  parse command line options		
command line options over-ride values in the input data file	 	
04 Mar 2009, 22 Sep 2009
------------------------------------------------------------------------------*/
void parse_options (
	int argc, char *argv[], 
	char IN_file[], char OUT_file[], 
	int *shear_flag,
	int *geom_flag,
	int *anlyz_flag,
	double *exagg_flag,
	int *lump_flag,
	int *modal_flag,
	double *tol_flag, 
	double *shift_flag,
	float *pan_flag,
	int *write_matrix,
	int *condense_flag,
	int *verbose,
	int *debug
){

	char	option;
	int	sfrv;	/* *scanf return value	*/

	/* default values */

	*shear_flag = *geom_flag  = *anlyz_flag = *lump_flag = *modal_flag = -1;
	*exagg_flag = *tol_flag = *shift_flag = -1.0;
	*pan_flag = *condense_flag = -1.0;
	*write_matrix = 0;
	*debug = 0; *verbose = 1;

	/* set up file names for the the input data and the output data */

	switch ( argc ) {
	 case 1: {
 		fprintf(stderr,"\n Frame3DD version: %s\n", VERSION);
		fprintf(stderr," Analysis of 2D and 3D structural frames with elastic and geometric stiffness.\n");
		fprintf(stderr," http://frame3dd.sourceforge.net\n\n");
		fprintf (stderr," Please enter the  input data file name: ");
		sfrv=scanf("%s", IN_file );
		if (sfrv != 1) sferr("IN_file");
		fprintf (stderr," Please enter the output data file name: ");
		sfrv=scanf("%s", OUT_file );
		if (sfrv != 1) sferr("OUT_file");
		return;
	 }
	 case 3: {
		if ( argv[1][0] != '-' ) {
			strcpy(  IN_file , argv[1] );
			strcpy( OUT_file , argv[2] );
			return;
		}
	 }
	}

	while ((option=getopt(argc,argv, "i:o:hvaqcdws:g:e:l:m:t:f:p:r:")) != -1){
		switch ( option ) {
			case 'i':		/* input data file name */
				strcpy(IN_file,optarg);
				break;
			case 'o':		/* output data file name */
				strcpy(OUT_file,optarg);
				break;
			case 'h':		/* help	*/
				display_help();
				exit(0);
			case 'v':		/*version */
				display_version();
				exit(0);
			case 'a':		/*version */
				display_version_about();
				exit(0);
			case 'q':		/* quiet */
				*verbose = 0;
				break;
			case 'c':		/* data check only */
				*anlyz_flag = 0;
				break;
			case 'd':		/* debug */
				*debug = 1;
				break;
			case 'w':		/* write stiffness and mass */
				*write_matrix = 1;
				break;
			case 's':		/* shear deformation */
				if (strcmp(optarg,"Off")==0)
					*shear_flag = 0;
				else if (strcmp(optarg,"On")==0)
					*shear_flag = 1;
				else {
					printf(" frame3dd command-line error"); 
					printf(": argument to -s option"); 
					printf(" should be either On or Off\n");
					exit(1);
				}
				break;
			case 'g':		/* geometric stiffness */
				if (strcmp(optarg,"Off")==0)
					*geom_flag = 0;
				else if (strcmp(optarg,"On")==0)
					*geom_flag = 1;
				else {
					printf(" frame3dd command-line error"); 
					printf(": argument to -g option"); 
					printf(" should be either On or Off\n");
					exit(1);
				}
				break;
			case 'e':		/* exaggeration factor */
				*exagg_flag = atof(optarg);
				break;
			case 'l':		/* lumped or consistent mass */
				if (strcmp(optarg,"Off")==0)
					*lump_flag = 0;
				else if (strcmp(optarg,"On")==0)
					*lump_flag = 1;
				else {
					printf(" frame3dd command-line error"); 
					printf(": argument to -l option"); 
					printf(" should be either On or Off\n");
					exit(1);
				}
				break;
			case 'm':		/* modal analysis method */
				if (strcmp(optarg,"J")==0)
					*modal_flag = 1;
				else if (strcmp(optarg,"S")==0)
					*modal_flag = 2;
				else {
					printf(" frame3dd command-line error"); 
					printf(": argument to -m option"); 
					printf(" should be either J or S\n\n");
					exit(1);
				}
				break;
			case 't':		/* modal analysis tolerence */
				*tol_flag = atof(optarg);
				if (*tol_flag == 0.0) {
					printf(" frame3dd command-line error"); 
					printf(": argument to -t option"); 
					printf(" should be a number.\n\n");
					exit(1);
				}
				break;
			case 'f':		/* modal analysis freq. shift */
				*shift_flag = atof(optarg);
				if (*shift_flag == 0.0) {
					printf(" frame3dd command-line error"); 
					printf(": argument to -f option"); 
					printf(" should be a number.\n\n");
					exit(1);
				}
				break;
			case 'p':		/* pan rate	*/
				*pan_flag = atof(optarg);
				if (*pan_flag == 0.0) {
					printf(" frame3dd command-line error"); 
					printf(": argument to -p option"); 
					printf(" should be a number.\n");
					exit(1);
				}
				break;
			case 'r':		/* matrix condensation method */
				*condense_flag = atoi(optarg);
				if (*condense_flag < 0 || *condense_flag > 3) {
					printf(" frame3dd command-line error"); 
					printf(": argument to -r option"); 
					printf(" should be 0, 1, or 2.\n\n");
					exit(1);
				}
				break;
			case '?':
				printf("  Missing argument or Unknown option: -%c\n\n", option );
				display_help();
				exit(1);
		}
	}

	if ( strcmp(IN_file,"\0") == 0 ) {
		fprintf (stderr," Please enter the  input data file name: ");
		sfrv=scanf("%s", IN_file );
		if (sfrv != 1) sferr("IN_file");
		fprintf (stderr," Please enter the output data file name: ");
		sfrv=scanf("%s", OUT_file );
		if (sfrv != 1) sferr("OUT_file");
	}
	if ( strcmp(IN_file,"\0") != 0 && strcmp(OUT_file,"\0") == 0 ) {
		strcpy( OUT_file, IN_file );
		strcat( OUT_file, ".out" );
	}
}


/*------------------------------------------------------------------------------
DISPLAY_HELP -  display help information to stderr	
04 Mar 2009, 22 Sep 2009
------------------------------------------------------------------------------*/
void display_help()
{
 fprintf(stderr,"\n Frame3DD version: %s\n", VERSION);
 fprintf(stderr," Analysis of 2D and 3D structural frames with elastic and geometric stiffness.\n");
 fprintf(stderr," http://frame3dd.sourceforge.net\n\n");
/* fprintf(stderr,"  Usage: frame3dd -i<input> -o<output> [-hvcq] [-s<On|Off>] [-g<On|Off>] [-e<value>] [-l<On|Off>] [-f<value>] [-m J|S] [-t<value>] [-p<value>] \n");
 */
 fprintf(stderr,"  Frame3DD may be run with interactive prompting for file names by typing ...\n");
 fprintf(stderr,"       frame3dd \n\n");
 fprintf(stderr,"  Frame3DD may be run without command-line options by typing ...\n");
 fprintf(stderr,"       frame3dd <InFile> <OutFile> \n\n");

 fprintf(stderr,"  Frame3DD may be run with command-line options by typing ...\n");
 fprintf(stderr,"       frame3dd -i <InFile> -o <OutFile> [OPTIONS] \n\n");

 fprintf(stderr," ... where [OPTIONS] over-rides values in the input data file and includes\n");
 fprintf(stderr,"     one or more of the following:\n\n");

 fprintf(stderr," -------------------------------------------------------------------------\n");
 fprintf(stderr,"  -i  <InFile>  the  input data file name --- described in the manual\n");
 fprintf(stderr,"  -o <OutFile>  the output data file name\n");
 fprintf(stderr,"  -h            print this help message and exit\n");
 fprintf(stderr,"  -v            display program version, website, brief help info and exit\n");
 fprintf(stderr,"  -a            display program version, website and exit\n");
 fprintf(stderr,"  -c            data check only - the output data reviews the input data\n");
 fprintf(stderr,"  -w            write stiffness and mass matrices to files named Ks Kd Md\n");
 fprintf(stderr,"  -q            suppress screen output except for warning messages\n");
 fprintf(stderr,"  -s  On|Off    On: include shear deformation or Off: neglect ...\n");
 fprintf(stderr,"  -g  On|Off    On: include geometric stiffness or Off: neglect ...\n");
 fprintf(stderr,"  -e <value>    level of deformation exaggeration for Gnuplot output\n");
 fprintf(stderr,"  -l  On|Off    On: lumped mass matrix or Off: consistent mass matrix\n");
 fprintf(stderr,"  -f <value>    modal frequency shift for unrestrained structures\n");
 fprintf(stderr,"  -m   J|S      modal analysis method: J=Jacobi-Subspace or S=Stodola\n");
 fprintf(stderr,"  -t <value>    convergence tolerance for modal analysis\n");
 fprintf(stderr,"  -p <value>    pan rate for mode shape animation\n");
 fprintf(stderr,"  -r <value>    matrix condensation method: 0, 1, 2, or 3 \n");
 fprintf(stderr," -------------------------------------------------------------------------\n");


}


/*------------------------------------------------------------------------------
DISPLAY_USAGE -  display usage information to stderr	
04 Mar 2009
------------------------------------------------------------------------------*/
void display_usage()
{
 fprintf(stderr,"\n Frame3DD version: %s\n", VERSION);
 fprintf(stderr," Analysis of 2D and 3D structural frames with elastic and geometric stiffness.\n");
 fprintf(stderr," http://frame3dd.sourceforge.net\n\n");
/* fprintf(stderr,"  Usage: frame3dd -i<input> -o<output> [-hvcq] [-s<On|Off>] [-g<On|Off>] [-e<value>] [-l<On|Off>] [-f<value>] [-m J|S] [-t<value>] [-p<value>] \n");
 */
 fprintf(stderr,"  Usage: frame3dd -i <input> -o <output> [OPTIONS] \n\n");

 fprintf(stderr,"  Type ...   frame3dd -h   ... for additional help information.\n\n");

}

/*------------------------------------------------------------------------------
DISPLAY_VERSION_HELP -  display version, website, and brief help info. to stderr
04 Mar 2009
------------------------------------------------------------------------------*/
void display_version()
{
 fprintf(stderr,"\n Frame3DD version: %s\n", VERSION);
 fprintf(stderr," Analysis of 2D and 3D structural frames with elastic and geometric stiffness.\n");
 fprintf(stderr," http://frame3dd.sourceforge.net\n\n");

 fprintf(stderr,"  Usage: frame3dd -i <input> -o <output> [OPTIONS] \n\n");

 fprintf(stderr,"  Type ...   frame3dd -h   ... for additional help information.\n\n");
}


/*------------------------------------------------------------------------------
DISPLAY_VERSION_ABOUT-  display version and website to stderr for 
running as a background process 
22 Sep 2009
Contributed by Barry Sanford, barry.sanford@trimjoist.com
------------------------------------------------------------------------------*/
void display_version_about()
{
 fprintf(stderr," Frame3DD version: %s\n", VERSION);
 fprintf(stderr," Analysis of 2D and 3D structural frames with elastic and geometric stiffness\n");
 fprintf(stderr," http://frame3dd.sourceforge.net\n");
}


/*------------------------------------------------------------------------------
READ_JOINT_DATA  -  read joint location data		
04 Jan 2009
------------------------------------------------------------------------------*/
void read_joint_data( FILE *fp, int nJ, vec3 *xyz, float *r )
{
	int	i, j,
		sfrv;		/* *scanf result	*/

	for (i=1;i<=nJ;i++) {		/* read joint coordinates	*/
		sfrv=fscanf(fp, "%d", &j );
		if (sfrv != 1) sferr("joint number in joint data");
		if ( j <= 0 || j > nJ ) {
		    fprintf(stderr,"\nERROR: in joint coordinate data, joint number out of range\n");
		    fprintf(stderr,"(joint id is %d <= 0 or > %d)\n", j, nJ);
		    exit(1);
		}
		sfrv=fscanf(fp, "%lf %lf %lf %f", &xyz[j].x, &xyz[j].y, &xyz[j].z, &r[j]);
		if (sfrv != 4) sferr("joint coordinates in joint data");
		/* fprintf(stderr,"\nj = %d, pos = (%lf, %lf, %lf), r = %f", j, xyz[j].x, xyz[j].y, xyz[j].z, r[j]); */
		r[j] = fabs(r[j]);
	}
	return;
}


/*------------------------------------------------------------------------------
READ_FRAME_ELEMENT_DATA  -  read frame element property data		
04 Jan 2009
------------------------------------------------------------------------------*/
void read_frame_element_data(
	FILE *fp,
	int nJ, int nE, vec3 *xyz, float *r,
	double *L, double *Le,
	int *J1, int *J2,
	float *Ax, float *Asy, float *Asz,
	float *J, float *Iy, float *Iz, float *E, float *G, float *p, float *d
){
	int	j1, j2, i, b;
	int	sfrv;	/* *scanf return value */

	for (i=1;i<=nE;i++) {		/* read frame element properties */
		sfrv=fscanf(fp, "%d", &b );
		if (sfrv != 1) sferr("frame element number in element data");
		if ( b <= 0 || b > nE ) {
		    fprintf(stderr,"  error in frame element property data: Element number out of range  ");
		    fprintf(stderr,"  Frame element number: %d  \n", b);
		    exit(1);
		}
		sfrv=fscanf(fp, "%d %d", &J1[b], &J2[b] );
		if (sfrv != 2) sferr("joint numbers in frame element data");
		if ( J1[b] <= 0 || J1[b] > nJ || J2[b] <= 0 || J2[b] > nJ ) {
		    fprintf(stderr,"  error in frame element property data: joint number out of range  ");
		    fprintf(stderr,"  Frame element number: %d \n", b);
		    exit(1);
		}
		sfrv=fscanf(fp, "%f %f %f", &Ax[b], &Asy[b], &Asz[b] );
		if (sfrv != 3) sferr("section areas in frame element data");
		sfrv=fscanf(fp, "%f %f %f", &J[b],  &Iy[b],  &Iz[b] );
		if (sfrv != 3) sferr("section inertias in frame element data");
		sfrv=fscanf(fp, "%f %f", &E[b], &G[b] );
		if (sfrv != 2) sferr("material moduli in frame element data");
		sfrv=fscanf(fp, "%f", &p[b]);
		if (sfrv != 1) sferr("roll angle in frame element data");

		p[b] = p[b]*PI/180.0;	/* convert from degrees to radians */

		sfrv=fscanf(fp, "%f",  &d[b]);
		if (sfrv != 1) sferr("mass density in frame element data");

		if ( Ax[b] < 0 || Asy[b] < 0 || Asz[b] < 0 ||
		      J[b] < 0 ||  Iy[b] < 0 ||  Iz[b] < 0	) {
		 fprintf(stderr,"  error in frame element property data: section property < 0  ");
		 fprintf(stderr,"  Frame element number: %d  \n", b);
		 exit(1);
		}
		if ( Ax[b] == 0 ) {
		 fprintf(stderr,"  error in frame element property data: cross section area is zero   ");
		 fprintf(stderr,"  Frame element number: %d  \n", b);
		 exit(1);
		}
		if ( (Asy[b] == 0 || Asz[b] == 0) && G[b] == 0 ) {
		 fprintf(stderr,"  error in frame element property data: a shear area and shear modulus are zero   ");
		 fprintf(stderr,"  Frame element number: %d  \n", b);
		 exit(1);
		}
		if ( J[b] == 0 ) {
		 fprintf(stderr,"  error in frame element property data: torsional moment of inertia is zero   ");
		 fprintf(stderr,"  Frame element number: %d  \n", b);
		 exit(1);
		}
		if ( Iy[b] == 0 || Iz[b] == 0 ) {
		 fprintf(stderr,"  error: cross section bending moment of inertia is zero   ");
		 fprintf(stderr,"  Frame element number : %d  \n", b);
		 exit(1);
		}
		if ( E[b] <= 0 || G[b] <= 0 ) {
		 fprintf(stderr,"  error : material elastic modulus E or G is not positive   ");
		 fprintf(stderr,"  Frame element number: %d  \n", b);
		 exit(1);
		}
		if ( d[b] <= 0 ) {
		 fprintf(stderr,"  error : mass density d is not positive   ");
		 fprintf(stderr,"  Frame element number: %d  \n", b);
		 exit(1);
		}
	}

	for (b=1;b<=nE;b++) {		/* calculate frame element lengths */
		j1 = J1[b];
		j2 = J2[b];

#define SQ(X) ((X)*(X))
		L[b] =	SQ( xyz[j2].x - xyz[j1].x ) +
			SQ( xyz[j2].y - xyz[j1].y ) +
			SQ( xyz[j2].z - xyz[j1].z );
#undef SQ

		L[b] = sqrt( L[b] );
		Le[b] = L[b] - r[j1] - r[j2];
		if ( j1 == j2 || L[b] == 0.0 ) {
		   fprintf(stderr,
			" Frame elements must start and stop at different joints\n");
		   fprintf(stderr,
			" frame element %d  J1= %d J2= %d L= %e\n", b, j1,j2, L[b] );
		   fprintf(stderr,
			"  Perhaps frame element number %d has not been specified.\n",i);
		   fprintf(stderr,
			"  or perhaps the Input Data file is missing expected data.\n");
		   exit(1);
		}
		if ( Le[b] <= 0.0 ) {
		   fprintf(stderr, " Joint radii are too large.\n");
		   fprintf(stderr,
			" frame element %d  J1= %d J2= %d L= %e \n", b, j1,j2, L[b] );
		   fprintf(stderr,
			" r1= %e r2= %e Le= %e \n", r[j1], r[j2], Le[b] );
		   exit(1);
		}
	}

	return;
}


/*------------------------------------------------------------------------------
READ_RUN_DATA  -  read information for analysis
29 Dec 2008
------------------------------------------------------------------------------*/
void read_run_data (
	FILE	*fp, 
	char	*OUT_file,	/* output data file name */
	int	*shear,
	int	shear_flag,
	int	*geom,
	int	geom_flag,
	char	*meshpath,
	char	*plotpath,
	double	*exagg_static,
	double	exagg_flag,
	int	*anlyz,
	int	anlyz_flag,
	int	debug
){
	int	full_len=0, len=0, i;
	char	base_file[96] = "EMPTY_BASE";
	char	mesh_file[96] = "EMPTY_MESH";
	int	sfrv;	/* *scanf return value */

	strcpy(base_file,OUT_file);	
	while ( base_file[len++] != '\0' )
		/* the length of the base_file */ ;
	full_len = len;
	while ( base_file[len--] != '.' && len > 0 )
		/* find the last '.' in base_file */ ;
	if ( len == 0 )	len = full_len;
	base_file[++len] = '\0';	/* end base_file at the last '.' */

	strcpy(plotpath,base_file);
	strcat(plotpath,".plt");

	while ( base_file[len] != '/' && base_file[len] != '\\' && len > 0 )
		len--;	/* find the last '/' or '\' in base_file */ 
	i = 0;
	while ( base_file[len] != '\0' )
		mesh_file[i++] = base_file[len++];
	mesh_file[i] = '\0';
	strcat(mesh_file,"-msh");
	output_path(mesh_file,meshpath,FRAME3DD_PATHMAX,NULL);

	if ( debug) {
		fprintf(stderr,"OUT_FILE = %s \n", OUT_file);
		fprintf(stderr,"BASE_FILE = %s \n", base_file);
		fprintf(stderr,"PLOTPATH = %s \n", plotpath);
		fprintf(stderr,"MESHPATH = %s \n", meshpath);
	}

	sfrv=fscanf( fp, "%d %d %lf", shear, geom,  exagg_static );
	if (sfrv != 3) sferr("shear, geom or exagg_static variables");

	if (*shear != 0 && *shear != 1) {
	    fprintf(stderr," Rember to specify shear deformations");
	    fprintf(stderr," with a 0 or a 1 after the frame element property info.\n");
	    exit(1);
	}

	if (*geom != 0 && *geom != 1) {
	    fprintf(stderr," Rember to specify geometric stiffness");
	    fprintf(stderr," with a 0 or a 1 after the frame element property info.\n");
	    exit(1);
	}

	if ( *exagg_static < 0.0 ) {
	    fprintf(stderr," Remember to specify an exageration");
	    fprintf(stderr," factor greater than zero\n");
	    exit(1);
	}

	/* over-ride values from input data file with command-line options */
	if ( shear_flag != -1   )	*shear = shear_flag;
	if ( geom_flag  != -1   )	*geom = geom_flag;
	if ( exagg_flag != -1.0 )	*exagg_static = exagg_flag;
	if ( anlyz_flag != -1.0 )	*anlyz = anlyz_flag;


	return;
}


/*-----------------------------------------------------------------------------
FRAME3DD_GETLINE -  get line into a character string. from K&R        03feb94
-----------------------------------------------------------------------------*/
void frame3dd_getline (
FILE	*fp,
char    *s,
int     lim
){
    int     c=0, i=0;

    while (--lim > 0 && (c=getc(fp)) != EOF && c != '\n' )
	s[i++] = c;
/*      if (c == '\n')  s[i++] = c;	*/
    s[i] = '\0';
    return;
}


/* platform-dependent path sperator character ... */

#if defined(WIN32) || defined(DJGPP)
static const char sep = '\\';
#else
static const char sep = '/';
#endif

/*----------------------------------------------------------------------------
TEMP_DIR
return platform-specific temp file location -- 
John Pye, Feb 2009
----------------------------------------------------------------------------*/
static const char *temp_dir(){
#if defined(WIN32) || defined(DJGPP)
	char *tmp;
	tmp = getenv("TEMP");
	if ( tmp==NULL ) {
		fprintf(stderr,
"ERROR: Environment Variables %%TEMP%% and %%FRAME3DD_OUTDIR%% are not set.  "
"At least one of these variables must be set so that Frame3DD knows where to "
"write its temporary files.  Set one of these variable, then re-run Frame3DD."
"The Frame3DD on-line documentation provides help on this issue.");
		exit(1);
	} 
#else
	const char *tmp = "/tmp";	/* Linux, Unix, OS X	*/
#endif
	return tmp;
}


/*------------------------------------------------------------------------------
OUTPUT_PATH
return path for output files using either current directory, or FRAME3DD_OUTDIR
if specified. -- 
John Pye, Feb 2009.
------------------------------------------------------------------------------*/
void output_path(const char *fname, char fullpath[], const int len, const char *default_outdir) {
	assert(fname!=NULL);
	int res;
	if ( fname[0]==sep ) {
		/* absolute output path specified */
//		res = snprintf(fullpath,len,"%s",fname);
		res = sprintf(fullpath,"%s",fname);
	} else {
		/* fprintf(stderr,"Generating output path for file '%s'\n",fname); */
		const char *outdir;
		outdir = getenv("FRAME3DD_OUTDIR");
		if (outdir==NULL) {
			if (default_outdir==NULL) {
				outdir = temp_dir();
			} else {
				outdir = default_outdir;
			}
		}
//		res = snprintf(fullpath,len,"%s%c%s",outdir,sep,fname);
		res = sprintf(fullpath,"%s%c%s",outdir,sep,fname);
	}
	if(res > len){
		fprintf(stderr,"ERROR: unable to construct output filename: overflow.\n");
		exit(1);
	}
	/* fprintf(stderr,"Output file path generated: %s\n",fullpath); */
}


/*-----------------------------------------------------------------------------
PARSE_INPUT                                                             7may03
 remove comments from the input file, and write a 'clean' input file
-----------------------------------------------------------------------------*/
void parse_input(FILE *fp, const char *tpath){
	FILE	*fpc;		/* cleaned inout/output file pointer	*/
	char	line[256];

	if ((fpc = fopen (tpath, "w")) == NULL) {
		fprintf (stderr,"ERROR: cannot open file '%s'\n",tpath);
		exit(1);
	}

	do {
		getline_no_comment(fp, line, 256);
		fprintf(fpc, "%s \n", line );
	} while ( line[0] != '_' && line[0] != EOF );

	fclose(fpc);

}


/*-----------------------------------------------------------------------------
GETLINE_NO_COMMENT                                                
 get a line into a character string. from K&R
 get the line only up to one of the following characters:  \n  %  #  ;  ? 
 ignore all comma (,) characters, ignore all double quote (") characters
09 Feb 2009
-----------------------------------------------------------------------------*/
void getline_no_comment(
	FILE *fp,   /**< pointer to the file from which to read */
	char *s,    /**< pointer to the string to which to write */
	int lim    /**< the longest anticipated line length  */
){
	int     c=0, i=0;

	while (--lim > 0 && (c=getc(fp)) != EOF && c != '\n' && c != '%' &&
		c != '#' && c != ';' && c != '?' ) {
		if (c != ',' && c != '"' )
			s[i++] = c;
		else
			s[i++] = ' ';
	/*      if (c == '\n')  s[i++] = c;     */
	}
	s[i] = '\0';
	if (c != '\n')
		while (--lim > 0 && (c=getc(fp)) != EOF && c != '\n')
		/* read the rest of the line, otherwise do nothing */ ;

	if ( c == EOF ) s[0] = EOF;

	return;
}


/*------------------------------------------------------------------------------
READ_REACTION_DATA - Read fixed joint displacement boundary conditions
29 Dec 2009
------------------------------------------------------------------------------*/
void read_reaction_data (
	FILE *fp, int DoF, int nJ, int *nR, int *R, int *sumR, int verbose
){
	int	i,j,l;
	int	sfrv;	/* *scanf return value */

	for (i=1; i<=DoF; i++)	R[i] = 0;

	sfrv=fscanf(fp,"%d", nR );	/* read restrained degrees of freedom */
	if (sfrv != 1) sferr("number of reactions in reaction data");
	if ( verbose ) {
		printf(" number of joints with reactions ");
		dots(stdout,20);
		printf(" nR = %3d ", *nR );
	}
	if ( *nR < 0 || *nR > DoF/6 ) {
		fprintf(stderr," number of joints with reactions ");
		dots(stderr,20);
		fprintf(stderr," nR = %3d ", *nR );
		fprintf(stderr,"  error: valid ranges for nR is 0 ... %d \n", DoF/6 );
		exit(1);
	}

	for (i=1; i <= *nR; i++) {
	    sfrv=fscanf(fp,"%d", &j);
	    if (sfrv != 1) sferr("joint number in reaction data");
	    for (l=5; l >=0; l--) {

		sfrv=fscanf(fp,"%d", &R[6*j-l] );
		if (sfrv != 1) sferr("reaction value in reaction data");

		if ( j > nJ ) {
		    fprintf(stderr,"  error in reaction data: joint number %d is greater than the number of joints, %d \n", j, nJ );
		    exit(1);
		}
		if ( R[6*j-l] != 0 && R[6*j-l] != 1 ) {
		    fprintf(stderr,"  error in reaction data: Reaction data must be 0 or 1\n");
		    fprintf(stderr,"  Data for joint %d, DoF %d is %d\n",
							j, 6-l, R[6*j-l] );
		    exit(1);
		}
	    }
	    *sumR = 0;
	    for (l=5; l >=0; l--) 	*sumR += R[6*j-l];
	    if ( *sumR == 0 ) {
		fprintf(stderr,"  error: joint %3d has no reactions\n", j);
		fprintf(stderr,"  Remove joint %3d from the list of reactions\n", j);
		fprintf(stderr,"  and set nR to %3d \n", *nR-1);
		exit(1);
	    }
	}
	*sumR=0;	for (i=1;i<=DoF;i++)	*sumR += R[i];
	if ( *sumR < 4 ) {
	    fprintf(stderr,"  warning:  Un-restrained structure\n");
	    fprintf(stderr,"  %d imposed reactions.\n", *sumR );
	    fprintf(stderr,"  At least 4 reactions are required to support static loads.\n");
	    /*	exit(1); */
	}
	if ( *sumR >= DoF ) {
	    fprintf(stderr,"  error in reaction data:  Fully restrained structure\n");
	    fprintf(stderr,"  %d imposed reactions >= %d degrees of freedom\n",
								*sumR, DoF );
	    exit(1);
	}

	return;
}


/*------------------------------------------------------------------------------
READ_AND_ASSEMBLE_LOADS  -
read load information data, assemble un-restrained load vectors
09 Sep 2008
------------------------------------------------------------------------------*/
void read_and_assemble_loads(
		FILE *fp,
		int nJ, int nE, int nL, int DoF,
		vec3 *xyz,
		double *L, double *Le,
		int *J1, int *J2,
		float *Ax, float *Asy, float *Asz,
		float *Iy, float *Iz, float *E, float *G,
		float *p,
		float *d, float *gX, float *gY, float *gZ, 
		int *R,
		int shear,
		int *nF, int *nU, int *nW, int *nP, int *nT, int *nD,
		double **Q,
		double **F_mech, double **F_temp,
		float ***U, float ***W, float ***P, float ***T, float **Dp,
		double ***feF_mech, double ***feF_temp, 
		int verbose
){
	float	hy, hz;			/* section dimensions in local coords */

	float	x1,x2, w1,w2;
	double	Ln, R1o, R2o, f01, f02; 

	double	Nx1, Vy1, Vz1, Mx1=0.0, My1=0.0, Mz1=0.0, /* fixed end forces */
		Nx2, Vy2, Vz2, Mx2=0.0, My2=0.0, Mz2=0.0,
		Ksy, Ksz, 		/* shear deformatn coefficients	*/
		a, b,			/* point load locations */
		t1, t2, t3, t4, t5, t6, t7, t8, t9;	/* 3D coord Xfrm coeffs */
	int	i,j,l, lc, n, j1, j2;
	int	sfrv;	/* *scanf return value */

	for (j=1; j<=DoF; j++)
		for (lc=1; lc <= nL; lc++)
			F_mech[lc][j] = F_temp[lc][j] = 0.0;
	for (i=1; i<=12; i++)
		for (n=1; n<=nE; n++)
			for (lc=1; lc <= nL; lc++)
				feF_mech[lc][n][i] = feF_temp[lc][n][i] = 0.0;

	for (i=1; i<=DoF; i++)	for (lc=1; lc<=nL; lc++) Dp[lc][i] = 0.0;

	for (i=1;i<=nE;i++)	for(j=1;j<=12;j++)	Q[i][j] = 0.0;

	for (lc = 1; lc <= nL; lc++) {		/* begin load-case loop */

	  if ( verbose ) printf(" load case %d of %d: \n", lc, nL );

	  /* gravity loads applied uniformly to all frame elements	*/
	  sfrv=fscanf(fp,"%f %f %f", &gX[lc], &gY[lc], &gZ[lc] );
	  if (sfrv != 3) sferr("gX gY gZ values in load data");

	  for (n=1; n<=nE; n++) {

		j1 = J1[n];	j2 = J2[n];

		coord_trans ( xyz, L[n], j1, j2,
			&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p[n] );

		feF_mech[lc][n][1]  = d[n]*Ax[n]*L[n]*gX[lc] / 2.0;
		feF_mech[lc][n][2]  = d[n]*Ax[n]*L[n]*gY[lc] / 2.0;
		feF_mech[lc][n][3]  = d[n]*Ax[n]*L[n]*gZ[lc] / 2.0;

		feF_mech[lc][n][4]  = d[n]*Ax[n]*L[n]*L[n] / 12.0 *
			( (-t4*t8+t5*t7)*gY[lc] + (-t4*t9+t6*t7)*gZ[lc] );
		feF_mech[lc][n][5]  = d[n]*Ax[n]*L[n]*L[n] / 12.0 *
			( (-t5*t7+t4*t8)*gX[lc] + (-t5*t9+t6*t8)*gZ[lc] );
		feF_mech[lc][n][6]  = d[n]*Ax[n]*L[n]*L[n] / 12.0 *
			( (-t6*t7+t4*t9)*gX[lc] + (-t6*t8+t5*t9)*gY[lc] );

		feF_mech[lc][n][7]  = d[n]*Ax[n]*L[n]*gX[lc] / 2.0;
		feF_mech[lc][n][8]  = d[n]*Ax[n]*L[n]*gY[lc] / 2.0;
		feF_mech[lc][n][9]  = d[n]*Ax[n]*L[n]*gZ[lc] / 2.0;

		feF_mech[lc][n][10] = d[n]*Ax[n]*L[n]*L[n] / 12.0 *
			( ( t4*t8-t5*t7)*gY[lc] + ( t4*t9-t6*t7)*gZ[lc] );
		feF_mech[lc][n][11] = d[n]*Ax[n]*L[n]*L[n] / 12.0 *
			( ( t5*t7-t4*t8)*gX[lc] + ( t5*t9-t6*t8)*gZ[lc] );
		feF_mech[lc][n][12] = d[n]*Ax[n]*L[n]*L[n] / 12.0 *
			( ( t6*t7-t4*t9)*gX[lc] + ( t6*t8-t5*t9)*gY[lc] );

		/* debugging  
		printf("n=%d ", n);
		for (l=1;l<=12;l++) {
			if (feF_mech[lc][n][l] != 0)
			   printf(" feF %d = %9.2e ", l, feF_mech[lc][n][l] );
		}
		printf("\n");   */
	  }					/* end gravity loads */


	  sfrv=fscanf(fp,"%d", &nF[lc] );	/* joint point loads	*/
	  if (sfrv != 1) sferr("nF value in load data");
	  if ( verbose ) {
		printf("  number of loaded joints ");
	  	dots(stdout,27);
		printf(" nF = %3d\n", nF[lc]);
	  }
	  for (i=1; i <= nF[lc]; i++) {	/* ! global structural coordinates ! */
		sfrv=fscanf(fp,"%d", &j);
		if (sfrv != 1) sferr("joint value in point load data");
		if ( j < 1 || j > nJ ) {
		    fprintf(stderr,"  error in joint load data: joint number out of range  ");
		    fprintf(stderr,"  Joint: %d  \n", j);
		    fprintf(stderr,"  Perhaps you did not specify %d joint loads \n", nF[lc] );
		   fprintf(stderr,
			"  or perhaps the Input Data file is missing expected data.\n");
		    exit(1);
		}

		for (l=5; l>=0; l--) {
			sfrv=fscanf(fp,"%lf", &F_mech[lc][6*j-l] );
			if (sfrv != 1) sferr("force value in point load data");
		}

		if ( F_mech[lc][6*j-5]==0 && F_mech[lc][6*j-4]==0 && F_mech[lc][6*j-3]==0 && F_mech[lc][6*j-2]==0 && F_mech[lc][6*j-1]==0 && F_mech[lc][6*j]==0 )
		    fprintf(stderr,"   warning: All joint loads applied at joint %d  are zero\n", j );
	  }					/* end joint point loads  */

	  sfrv=fscanf(fp,"%d", &nU[lc] );	/* uniformly distributed loads */
	  if (sfrv != 1) sferr("nU value in uniform load data");
	  if ( verbose ) {
		printf("  number of uniformly distributed loads ");
	  	dots(stdout,13);
	  	printf(" nU = %3d\n", nU[lc]);
	  }
	  if ( nU[lc] < 0 || nU[lc] > nE ) {
		fprintf(stderr,"  number of uniformly distributed loads ");
	  	dots(stderr,13);
	  	fprintf(stderr," nU = %3d\n", nU[lc]);
		fprintf(stderr,"  error: valid ranges for nU is 0 ... %d \n", nE );
		exit(1);
	  }
	  for (i=1; i <= nU[lc]; i++) {	/* ! local element coordinates ! */
		sfrv=fscanf(fp,"%d", &n );
	  	if (sfrv != 1) sferr("frame element number in uniform load data");
		if ( n < 1 || n > nE ) {
		    fprintf(stderr,"  error in uniform distributed loads: element number %d is out of range\n",n);
		    exit(1);
		}
		U[lc][i][1] = (double) n;
		for (l=2; l<=4; l++) {
			sfrv=fscanf(fp,"%f", &U[lc][i][l] );
	  		if (sfrv != 1) sferr("load value in uniform load data");
		}

		if ( U[lc][i][2]==0 && U[lc][i][3]==0 && U[lc][i][4]==0 )
		    fprintf(stderr,"   warning: All distributed loads applied to frame element %d  are zero\n", n );

		Nx1 = Nx2 = U[lc][i][2]*Le[n] / 2.0;
		Vy1 = Vy2 = U[lc][i][3]*Le[n] / 2.0;
		Vz1 = Vz2 = U[lc][i][4]*Le[n] / 2.0;
		Mx1 = Mx2 = 0.0;
		My1 = -U[lc][i][4]*Le[n]*Le[n] / 12.0;	My2 = -My1;
		Mz1 =  U[lc][i][3]*Le[n]*Le[n] / 12.0;	Mz2 = -Mz1;

		/* debugging 
		printf("n=%d Vy=%9.2e Vz=%9.2e My=%9.2e Mz=%9.2e\n",
						n, Vy1,Vz1, My1,Mz1 );	*/

		j1 = J1[n];	j2 = J2[n];

		coord_trans ( xyz, L[n], j1, j2,
			&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p[n] );

		/* debugging
		printf("t1=%5.2f t2=%5.2f t3=%5.2f \n", t1, t2, t3 );
		printf("t4=%5.2f t5=%5.2f t6=%5.2f \n", t4, t5, t6 );
		printf("t7=%5.2f t8=%5.2f t9=%5.2f \n", t7, t8, t9 ); */

		/* {F} = [T]'{Q} */
		feF_mech[lc][n][1]  += ( Nx1*t1 + Vy1*t4 + Vz1*t7 );
		feF_mech[lc][n][2]  += ( Nx1*t2 + Vy1*t5 + Vz1*t8 );
		feF_mech[lc][n][3]  += ( Nx1*t3 + Vy1*t6 + Vz1*t9 );
		feF_mech[lc][n][4]  += ( Mx1*t1 + My1*t4 + Mz1*t7 );
		feF_mech[lc][n][5]  += ( Mx1*t2 + My1*t5 + Mz1*t8 );
		feF_mech[lc][n][6]  += ( Mx1*t3 + My1*t6 + Mz1*t9 );

		feF_mech[lc][n][7]  += ( Nx2*t1 + Vy2*t4 + Vz2*t7 );
		feF_mech[lc][n][8]  += ( Nx2*t2 + Vy2*t5 + Vz2*t8 );
		feF_mech[lc][n][9]  += ( Nx2*t3 + Vy2*t6 + Vz2*t9 );
		feF_mech[lc][n][10] += ( Mx2*t1 + My2*t4 + Mz2*t7 );
		feF_mech[lc][n][11] += ( Mx2*t2 + My2*t5 + Mz2*t8 );
		feF_mech[lc][n][12] += ( Mx2*t3 + My2*t6 + Mz2*t9 );

		/* debugging  
		printf("n=%d ", n);
		for (l=1;l<=12;l++) {
			if (feF_mech[lc][n][l] != 0)
			   printf(" feF %d = %9.2e ", l, feF_mech[lc][n][l] );
		}
		printf("\n"); */ 
	  }				/* end uniformly distributed loads */

	  sfrv=fscanf(fp,"%d", &nW[lc] ); /* trapezoidally distributed loads */
	  if (sfrv != 1) sferr("nW value in load data");
	  if ( verbose ) {
		printf("  number of trapezoidally distributed loads ");
	  	dots(stdout,9);
	  	printf(" nW = %3d\n", nW[lc]);
	  }
	  if ( nW[lc] < 0 || nW[lc] > 10*nE ) {
		fprintf(stderr,"  error: valid ranges for nW is 0 ... %d \n", 10*nE );
		exit(1);
	  }
	  for (i=1; i <= nW[lc]; i++) {	/* ! local element coordinates ! */
		sfrv=fscanf(fp,"%d", &n );
	  	if (sfrv != 1) sferr("frame element number in uniform load data");
		if ( n < 1 || n > nE ) {
		    fprintf(stderr,"  error in uniform distributed loads: element number %d is out of range\n",n);
		    exit(1);
		}
		W[lc][i][1] = (double) n;
		for (l=2; l<=13; l++) {
			sfrv=fscanf(fp,"%f", &W[lc][i][l] );
			if (sfrv != 1) sferr("value in trapezoidal load data");
		}

		Ln = L[n];

		/* error checking */

		if ( W[lc][i][ 4]==0 && W[lc][i][ 5]==0 &&
		     W[lc][i][ 7]==0 && W[lc][i][ 8]==0 &&
		     W[lc][i][12]==0 && W[lc][i][13]==0 ) {
		  fprintf(stderr,"   warning: All trapezoidal loads applied to frame element %d  are zero\n", n );
		  fprintf(stderr,"     load case: %d , element %d , load %d\n ", lc, n, i );
		}

		if ( W[lc][i][ 2] < 0 ) {
		  fprintf(stderr,"   error in x-axis trapezoidal loads, ");
		  fprintf(stderr,"    load case: %d , element %d , load %d\n ", lc, n, i );
		  fprintf(stderr,"    starting location = %f < 0\n",W[lc][i][2]);
		  exit(1);
		}
		if ( W[lc][i][ 2] > W[lc][i][3] ) {
		  fprintf(stderr,"   error in x-axis trapezoidal loads, ");
		  fprintf(stderr,"    load case: %d , element %d , load %d\n ", lc, n, i );
		  fprintf(stderr,"    starting location = %f > ending location = %f \n",W[lc][i][2], W[lc][i][3] );
		  exit(1);
		}
		if ( W[lc][i][ 3] > Ln ) {
		  fprintf(stderr,"   error in x-axis trapezoidal loads, ");
		  fprintf(stderr,"    load case: %d , element %d , load %d\n ", lc, n, i );
		  fprintf(stderr,"    ending location = %f > L (%f) \n",W[lc][i][3], Ln );
		  exit(1);
		}
		if ( W[lc][i][ 6] < 0 ) {
		  fprintf(stderr,"   error in y-axis trapezoidal loads, ");
		  fprintf(stderr,"    load case: %d , element %d , load %d\n ", lc, n, i );
		  fprintf(stderr,"    starting location = %f < 0\n",W[lc][i][6]);
		  exit(1);
		}
		if ( W[lc][i][ 6] > W[lc][i][7] ) {
		  fprintf(stderr,"   error in y-axis trapezoidal loads, ");
		  fprintf(stderr,"    load case: %d , element %d , load %d\n ", lc, n, i );
		  fprintf(stderr,"    starting location = %f > ending location = %f \n",W[lc][i][6], W[lc][i][7] );
		  exit(1);
		}
		if ( W[lc][i][ 7] > Ln ) {
		  fprintf(stderr,"   error in y-axis trapezoidal loads, ");
		  fprintf(stderr,"    load case: %d , element %d , load %d\n ", lc, n, i );
		  fprintf(stderr,"    ending location = %f > L (%f) \n",W[lc][i][7],Ln );
		  exit(1);
		}
		if ( W[lc][i][10] < 0 ) {
		  fprintf(stderr,"   error in z-axis trapezoidal loads, ");
		  fprintf(stderr,"    load case: %d , element %d , load %d\n ", lc, n, i );
		  fprintf(stderr,"    starting location = %f < 0\n",W[lc][i][10]);
		  exit(1);
		}
		if ( W[lc][i][10] > W[lc][i][11] ) {
		  fprintf(stderr,"   error in z-axis trapezoidal loads, ");
		  fprintf(stderr,"    load case: %d , element %d , load %d\n ", lc, n, i );
		  fprintf(stderr,"    starting location = %f > ending location = %f \n",W[lc][i][10], W[lc][i][11] );
		  exit(1);
		}
		if ( W[lc][i][11] > Ln ) {
		  fprintf(stderr,"   error in z-axis trapezoidal loads, ");
		  fprintf(stderr,"    load case: %d , element %d , load %d\n ", lc, n, i );
		  fprintf(stderr,"    ending location = %f > L (%f) \n",W[lc][i][11], Ln );
		  exit(1);
		}

		if ( shear ) {
			Ksy = (12.0*E[n]*Iz[n]) / (G[n]*Asy[n]*Le[n]*Le[n]);
			Ksz = (12.0*E[n]*Iy[n]) / (G[n]*Asz[n]*Le[n]*Le[n]);
		} else	Ksy = Ksz = 0.0;

		/* x-axis trapezoidal loads (along the frame element length) */
		x1 =  W[lc][i][2]; x2 =  W[lc][i][3];
		w1 =  W[lc][i][4]; w2 =  W[lc][i][5];

		Nx1 = ( 3.0*(w1+w2)*Ln*(x2-x1) - (2.0*w2+w1)*x2*x2 + (w2-w1)*x2*x1 + (2.0*w1+w2)*x1*x1 ) / (6.0*Ln);
		Nx2 = ( -(2.0*w1+w2)*x1*x1 + (2.0*w2+w1)*x2*x2  - (w1-w2)*x1*x2 ) / ( 6.0*Ln );

		/* y-axis trapezoidal loads (across the frame element length) */
		x1 =  W[lc][i][6];  x2 = W[lc][i][7];
		w1 =  W[lc][i][8]; w2 =  W[lc][i][9];

		R1o = ( (2.0*w1+w2)*x1*x1 - (w1+2.0*w2)*x2*x2 + 
			 3.0*(w1+w2)*Ln*(x2-x1) - (w1-w2)*x1*x2 ) / (6.0*Ln);
		R2o = ( (w1+2.0*w2)*x2*x2 + (w1-w2)*x1*x2 - 
			(2.0*w1+w2)*x1*x1 ) / (6.0*Ln);

		f01 = (  3.0*(w2+4.0*w1)*x1*x1*x1*x1 -  3.0*(w1+4.0*w2)*x2*x2*x2*x2
		      - 15.0*(w2+3.0*w1)*Ln*x1*x1*x1 + 15.0*(w1+3.0*w2)*Ln*x2*x2*x2
		      -  3.0*(w1-w2)*x1*x2*(x1*x1 + x2*x2)
		      + 20.0*(w2+2.0*w1)*Ln*Ln*x1*x1 - 20.0*(w1+2.0*w2)*Ln*Ln*x2*x2
		      + 15.0*(w1-w2)*Ln*x1*x2*(x1+x2)
		      -  3.0*(w1-w2)*x1*x1*x2*x2 - 20.0*(w1-w2)*Ln*Ln*x1*x2 ) / 360.0;

		f02 = (  3.0*(w2+4.0*w1)*x1*x1*x1*x1 - 3.0*(w1+4.0*w2)*x2*x2*x2*x2
		      -  3.0*(w1-w2)*x1*x2*(x1*x1+x2*x2)
		      - 10.0*(w2+2.0*w1)*Ln*Ln*x1*x1 + 10.0*(w1+2.0*w2)*Ln*Ln*x2*x2
		      -  3.0*(w1-w2)*x1*x1*x2*x2 + 10.0*(w1-w2)*Ln*Ln*x1*x2 ) / 360.0;

		Mz1 = -( 4.0*f01 + 2.0*f02 + Ksy*(f01 - f02) ) / ( Ln*Ln*(1.0+Ksy) );
		Mz2 = -( 2.0*f01 + 4.0*f02 - Ksy*(f01 - f02) ) / ( Ln*Ln*(1.0+Ksy) );

		Vy1 =  R1o + Mz1/Ln + Mz2/Ln;
		Vy2 =  R2o - Mz1/Ln - Mz2/Ln;

		/* z-axis trapezoidal loads (across the frame element length) */
		x1 =  W[lc][i][10]; x2 =  W[lc][i][11];
		w1 =  W[lc][i][12]; w2 =  W[lc][i][13];

		R1o = ( (2.0*w1+w2)*x1*x1 - (w1+2.0*w2)*x2*x2 + 
			 3.0*(w1+w2)*Ln*(x2-x1) - (w1-w2)*x1*x2 ) / (6.0*Ln);
		R2o = ( (w1+2.0*w2)*x2*x2 + (w1-w2)*x1*x2 -
			(2.0*w1+w2)*x1*x1 ) / (6.0*Ln);

		f01 = (  3.0*(w2+4.0*w1)*x1*x1*x1*x1 -  3.0*(w1+4.0*w2)*x2*x2*x2*x2
		      - 15.0*(w2+3.0*w1)*Ln*x1*x1*x1 + 15.0*(w1+3.0*w2)*Ln*x2*x2*x2
		      -  3.0*(w1-w2)*x1*x2*(x1*x1 + x2*x2)
		      + 20.0*(w2+2.0*w1)*Ln*Ln*x1*x1 - 20.0*(w1+2.0*w2)*Ln*Ln*x2*x2
		      + 15.0*(w1-w2)*Ln*x1*x2*(x1+x2)
		      -  3.0*(w1-w2)*x1*x1*x2*x2 - 20.0*(w1-w2)*Ln*Ln*x1*x2 ) / 360.0;

		f02 = (  3.0*(w2+4.0*w1)*x1*x1*x1*x1 - 3.0*(w1+4.0*w2)*x2*x2*x2*x2
		      -  3.0*(w1-w2)*x1*x2*(x1*x1+x2*x2)
		      - 10.0*(w2+2.0*w1)*Ln*Ln*x1*x1 + 10.0*(w1+2.0*w2)*Ln*Ln*x2*x2
		      -  3.0*(w1-w2)*x1*x1*x2*x2 + 10.0*(w1-w2)*Ln*Ln*x1*x2 ) / 360.0;

		My1 = ( 4.0*f01 + 2.0*f02 + Ksz*(f01 - f02) ) / ( Ln*Ln*(1.0+Ksz) );
		My2 = ( 2.0*f01 + 4.0*f02 - Ksz*(f01 - f02) ) / ( Ln*Ln*(1.0+Ksz) );

		Vz1 =  R1o - My1/Ln - My2/Ln;
		Vz2 =  R2o + My1/Ln + My2/Ln;

		/* debugging 
		printf("n=%d\n Nx1=%9.3f\n Nx2=%9.3f\n Vy1=%9.3f\n Vy2=%9.3f\n Vz1=%9.3f\n Vz2=%9.3f\n My1=%9.3f\n My2=%9.3f\n Mz1=%9.3f\n Mz2=%9.3f\n",
				n, Nx1,Nx2,Vy1,Vy2,Vz1,Vz2, My1,My2,Mz1,Mz2 );	*/

		j1 = J1[n];	j2 = J2[n];

		coord_trans ( xyz, Ln, j1, j2,
			&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p[n] );

		/* debugging
		printf("t1=%5.2f t2=%5.2f t3=%5.2f \n", t1, t2, t3 );
		printf("t4=%5.2f t5=%5.2f t6=%5.2f \n", t4, t5, t6 );
		printf("t7=%5.2f t8=%5.2f t9=%5.2f \n", t7, t8, t9 ); */

		/* {F} = [T]'{Q} */
		feF_mech[lc][n][1]  += ( Nx1*t1 + Vy1*t4 + Vz1*t7 );
		feF_mech[lc][n][2]  += ( Nx1*t2 + Vy1*t5 + Vz1*t8 );
		feF_mech[lc][n][3]  += ( Nx1*t3 + Vy1*t6 + Vz1*t9 );
		feF_mech[lc][n][4]  += ( Mx1*t1 + My1*t4 + Mz1*t7 );
		feF_mech[lc][n][5]  += ( Mx1*t2 + My1*t5 + Mz1*t8 );
		feF_mech[lc][n][6]  += ( Mx1*t3 + My1*t6 + Mz1*t9 );

		feF_mech[lc][n][7]  += ( Nx2*t1 + Vy2*t4 + Vz2*t7 );
		feF_mech[lc][n][8]  += ( Nx2*t2 + Vy2*t5 + Vz2*t8 );
		feF_mech[lc][n][9]  += ( Nx2*t3 + Vy2*t6 + Vz2*t9 );
		feF_mech[lc][n][10] += ( Mx2*t1 + My2*t4 + Mz2*t7 );
		feF_mech[lc][n][11] += ( Mx2*t2 + My2*t5 + Mz2*t8 );
		feF_mech[lc][n][12] += ( Mx2*t3 + My2*t6 + Mz2*t9 );

		/* debugging
		for (l=1;l<=13;l++) printf(" %9.2e ", W[lc][i][l] );
		printf("\n"); 
		printf("n=%d ", n);
		for (l=1;l<=12;l++) {
			if (feF_mech[lc][n][l] != 0)
			   printf(" feF %d = %9.3f ", l, feF_mech[lc][n][l] );
		}
		printf("\n"); */
	  }			/* end trapezoidally distributed loads */

	  sfrv=fscanf(fp,"%d", &nP[lc] );	/* element point loads	*/
	  if (sfrv != 1) sferr("nP value load data");
	  if ( verbose ) {
	  	printf("  number of concentrated frame element point loads ");
	  	dots(stdout,2);
	  	printf(" nP = %3d\n", nP[lc]);
	  }
	  if ( nP[lc] < 0 || nP[lc] > nE ) {
	  	fprintf(stderr,"  number of concentrated frame element point loads ");
	  	dots(stderr,3);
	  	fprintf(stderr," nP = %3d\n", nP[lc]);
		fprintf(stderr,"  error: valid ranges for nP is 0 ... %d \n", nE );
		exit(1);
	  }
	  for (i=1; i <= nP[lc]; i++) {	/* ! local element coordinates ! */
		sfrv=fscanf(fp,"%d", &n );
		if (sfrv != 1) sferr("frame element number value point load data");
		if ( n < 1 || n > nE ) {
		    fprintf(stderr,"   error in internal point loads: frame element number %d is out of range\n",n);
		    exit(1);
		}
		P[lc][i][1] = (double) n;
		for (l=2; l<=5; l++) { 
			sfrv=fscanf(fp,"%f", &P[lc][i][l] );
			if (sfrv != 1) sferr("value in point load data");
		}
		a = P[lc][i][5];	b = L[n] - a;

		if ( a < 0 || L[n] < a || b < 0 || L[n] < b ) {
		    fprintf(stderr,"  error in point load data: Point load coord. out of range\n");
		    fprintf(stderr,"  Frame element number: %d  L: %lf  load coord.: %lf\n",
							n, L[n], P[lc][i][5] );
		    exit(1);
		}

		if ( shear ) {
			Ksy = (12.0*E[n]*Iz[n]) / (G[n]*Asy[n]*Le[n]*Le[n]);
			Ksz = (12.0*E[n]*Iy[n]) / (G[n]*Asz[n]*Le[n]*Le[n]);
		} else	Ksy = Ksz = 0.0;

		Ln = L[n];

		Nx1 = P[lc][i][2]*a/Ln;
		Nx2 = P[lc][i][2]*b/Ln;

		Vy1 = (1./(1.+Ksz))    * P[lc][i][3]*b*b*(3.*a + b) / ( Ln*Ln*Ln ) +
			(Ksz/(1.+Ksz)) * P[lc][i][3]*b/Ln;
		Vy2 = (1./(1.+Ksz))    * P[lc][i][3]*a*a*(3.*b + a) / ( Ln*Ln*Ln ) +
			(Ksz/(1.+Ksz)) * P[lc][i][3]*a/Ln;

		Vz1 = (1./(1.+Ksy))    * P[lc][i][4]*b*b*(3.*a + b) / ( Ln*Ln*Ln ) +
			(Ksy/(1.+Ksy)) * P[lc][i][4]*b/Ln;
		Vz2 = (1./(1.+Ksy))    * P[lc][i][4]*a*a*(3.*b + a) / ( Ln*Ln*Ln ) +
			(Ksy/(1.+Ksy)) * P[lc][i][4]*a/Ln;

		Mx1 = Mx2 = 0.0;

		My1 = -(1./(1.+Ksy))  * P[lc][i][4]*a*b*b / ( Ln*Ln ) -
			(Ksy/(1.+Ksy))* P[lc][i][4]*a*b   / (2.*Ln);
		My2 =  (1./(1.+Ksy))  * P[lc][i][4]*a*a*b / ( Ln*Ln ) +
			(Ksy/(1.+Ksy))* P[lc][i][4]*a*b   / (2.*Ln);

		Mz1 =  (1./(1.+Ksz))  * P[lc][i][3]*a*b*b / ( Ln*Ln ) +
			(Ksz/(1.+Ksz))* P[lc][i][3]*a*b   / (2.*Ln);
		Mz2 = -(1./(1.+Ksz))  * P[lc][i][3]*a*a*b / ( Ln*Ln ) -
			(Ksz/(1.+Ksz))* P[lc][i][3]*a*b   / (2.*Ln);

		j1 = J1[n];	j2 = J2[n];

		coord_trans ( xyz, Ln, j1, j2,
			&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p[n] );

		/* {F} = [T]'{Q} */
		feF_mech[lc][n][1]  += ( Nx1*t1 + Vy1*t4 + Vz1*t7 );
		feF_mech[lc][n][2]  += ( Nx1*t2 + Vy1*t5 + Vz1*t8 );
		feF_mech[lc][n][3]  += ( Nx1*t3 + Vy1*t6 + Vz1*t9 );
		feF_mech[lc][n][4]  += ( Mx1*t1 + My1*t4 + Mz1*t7 );
		feF_mech[lc][n][5]  += ( Mx1*t2 + My1*t5 + Mz1*t8 );
		feF_mech[lc][n][6]  += ( Mx1*t3 + My1*t6 + Mz1*t9 );

		feF_mech[lc][n][7]  += ( Nx2*t1 + Vy2*t4 + Vz2*t7 );
		feF_mech[lc][n][8]  += ( Nx2*t2 + Vy2*t5 + Vz2*t8 );
		feF_mech[lc][n][9]  += ( Nx2*t3 + Vy2*t6 + Vz2*t9 );
		feF_mech[lc][n][10] += ( Mx2*t1 + My2*t4 + Mz2*t7 );
		feF_mech[lc][n][11] += ( Mx2*t2 + My2*t5 + Mz2*t8 );
		feF_mech[lc][n][12] += ( Mx2*t3 + My2*t6 + Mz2*t9 );
	  }					/* end element point loads */

	  sfrv=fscanf(fp,"%d", &nT[lc] );	/* thermal loads	*/
	  if (sfrv != 1) sferr("nT value in load data");
	  if ( verbose ) {
	  	printf("  number of temperature changes ");
	  	dots(stdout,21);
	  	printf(" nT = %3d\n", nT[lc] );
	  }
	  if ( nT[lc] < 0 || nT[lc] > nE ) {
	  	fprintf(stderr,"  number of temperature changes ");
	  	dots(stderr,21);
	  	fprintf(stderr," nT = %3d\n", nT[lc] );
		fprintf(stderr,"  error: valid ranges for nT is 0 ... %d \n", nE );
		exit(1);
	  }
	  for (i=1; i <= nT[lc]; i++) {	/* ! local element coordinates ! */
		sfrv=fscanf(fp,"%d", &n );
		if (sfrv != 1) sferr("frame element number in temperature load data");
		if ( n < 1 || n > nE ) {
		    fprintf(stderr,"  error in temperature loads: frame element number %d is out of range\n",n);
		    exit(1);
		}
		T[lc][i][1] = (double) n;
		for (l=2; l<=8; l++) {
			sfrv=fscanf(fp,"%f", &T[lc][i][l] );
			if (sfrv != 1) sferr("value in temperature load data");
		}
		a  = T[lc][i][2];
		hy = T[lc][i][3];
		hz = T[lc][i][4];

		if ( hy < 0 || hz < 0 ) {
		    fprintf(stderr,"  error in thermal load data: section dimension < 0\n");
		    fprintf(stderr,"  Frame element number: %d  hy: %f  hz: %f\n", n,hy,hz);
		    exit(1);
		}

		Nx2 = (a/4.0)*( T[lc][i][5]+T[lc][i][6]+T[lc][i][7]+T[lc][i][8])*E[n]*Ax[n];
		Nx1 = -Nx2;
		Vy1 = Vy2 = Vz1 = Vz2 = 0.0;
		Mx1 = Mx2 = 0.0;
		My1 =  (a/hz)*(T[lc][i][8]-T[lc][i][7])*E[n]*Iy[n];
		My2 = -My1;
		Mz1 =  (a/hy)*(T[lc][i][5]-T[lc][i][6])*E[n]*Iz[n];
		Mz2 = -Mz1;

		j1 = J1[n];	j2 = J2[n];

		coord_trans ( xyz, L[n], j1, j2,
			&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p[n] );

		/* {F} = [T]'{Q} */
		feF_temp[lc][n][1]  += ( Nx1*t1 + Vy1*t4 + Vz1*t7 );
		feF_temp[lc][n][2]  += ( Nx1*t2 + Vy1*t5 + Vz1*t8 );
		feF_temp[lc][n][3]  += ( Nx1*t3 + Vy1*t6 + Vz1*t9 );
		feF_temp[lc][n][4]  += ( Mx1*t1 + My1*t4 + Mz1*t7 );
		feF_temp[lc][n][5]  += ( Mx1*t2 + My1*t5 + Mz1*t8 );
		feF_temp[lc][n][6]  += ( Mx1*t3 + My1*t6 + Mz1*t9 );

		feF_temp[lc][n][7]  += ( Nx2*t1 + Vy2*t4 + Vz2*t7 );
		feF_temp[lc][n][8]  += ( Nx2*t2 + Vy2*t5 + Vz2*t8 );
		feF_temp[lc][n][9]  += ( Nx2*t3 + Vy2*t6 + Vz2*t9 );
		feF_temp[lc][n][10] += ( Mx2*t1 + My2*t4 + Mz2*t7 );
		feF_temp[lc][n][11] += ( Mx2*t2 + My2*t5 + Mz2*t8 );
		feF_temp[lc][n][12] += ( Mx2*t3 + My2*t6 + Mz2*t9 );
	  }				/* end thermal loads	*/

	  /* debugging ...  check feF's prior to asembly 
	  for (n=1; n<=nE; n++) {	
		printf("n=%d ", n);
		for (l=1;l<=12;l++) {
			if (feF_mech[lc][n][l] != 0)
			   printf(" feF %d = %9.2e ", l, feF_mech[lc][n][l] );
		}
		printf("\n"); 
	  } */

	  for (n=1; n<=nE; n++) {
	     j1 = J1[n];	j2 = J2[n];
	     for (i=1; i<= 6; i++) F_mech[lc][6*j1- 6+i] += feF_mech[lc][n][i];
	     for (i=7; i<=12; i++) F_mech[lc][6*j2-12+i] += feF_mech[lc][n][i];
	     for (i=1; i<= 6; i++) F_temp[lc][6*j1- 6+i] += feF_temp[lc][n][i];
	     for (i=7; i<=12; i++) F_temp[lc][6*j2-12+i] += feF_temp[lc][n][i];
	  }

	  sfrv=fscanf(fp,"%d", &nD[lc] ); /* read prescribed displacements */
	  if (sfrv != 1) sferr("nD value in load data");
	  if ( verbose ) {
	  	printf("  number of prescribed displacements ");
	  	dots(stdout,16);
	  	printf(" nD = %3d\n", nD[lc] );
	  }
	  for (i=1; i <= nD[lc]; i++) {
		sfrv=fscanf(fp,"%d", &j);
		if (sfrv != 1) sferr("joint number value in prescribed displacement data");
		for (l=5; l >=0; l--) {
			sfrv=fscanf(fp,"%f", &Dp[lc][6*j-l] );
			if (sfrv != 1) sferr("prescribed displacement value");
			if ( R[6*j-l] == 0 && Dp[lc][6*j-l] != 0.0 ) {
			    printf(" Initial displacements can be prescribed");
			    printf(" only at restrained coordinates\n");
			    printf(" joint: %d  dof: %d  R: %d\n",
							j, 6-l, R[6*j-l] );
			    exit(1);
			}
		}
	  }

	}					/* end load-case loop */

	return;
}


/*------------------------------------------------------------------------------
READ_MASS_DATA  -  read element densities and extra inertial mass data	16aug01
------------------------------------------------------------------------------*/
void read_mass_data(
		FILE *fp,
		char *OUT_file, 
		int nJ, int nE, int *nI, int *nX, 
		float *d, float *BMs,
		float *JMs, float *JMx, float *JMy, float *JMz,
		double *L, float *Ax,
		double *total_mass, double *struct_mass,
		int *nM, int *Mmethod, int modal_flag, 
		int *lump, int lump_flag, 
		double *tol, double tol_flag, double *shift, double shift_flag,
		double *exagg_modal, 
		char modepath[],
		int anim[], float *pan, float pan_flag,
		int verbose, int debug
){
/*	double	ms = 0.0; */
	int	i,j, jnt, m, b, nA;
	int	full_len=0, len=0;
	int	sfrv;	/* *scanf return value	*/

	char	base_file[96] = "EMPTY_BASE";
	char	mode_file[96] = "EMPTY_MODE";

	*total_mass = *struct_mass = 0.0;

	sfrv=fscanf ( fp, "%d", nM );
	if (sfrv != 1) sferr("nM value in mass data");

	if ( verbose ) {
		printf(" number of dynamic modes ");
		dots(stdout,28);
		printf(" nM = %3d\n", *nM);
	}

	if ( *nM < 1 || sfrv != 1 ) {
		*nM = 0;
		return;
	}

	sfrv=fscanf( fp, "%d", Mmethod );
	if (sfrv != 1) sferr("Mmethod value in mass data");
	if ( modal_flag != -1 )	*Mmethod = modal_flag;

	if ( verbose ) {
		printf(" modal analysis method ");
		dots(stdout,30);
		printf(" %3d ",*Mmethod);
		if ( *Mmethod == 1 ) printf(" (Subspace-Jacobi)\n");
		if ( *Mmethod == 2 ) printf(" (Stodola)\n");
	}


#ifdef MASSDATA_DEBUG
	FILE	*mf;				// mass data file
	mf = fopen("MassData.txt","w");		// open mass data file
	if ((mf = fopen ("MassData.txt", "w")) == NULL) {
	  fprintf (stderr," error: cannot open file 'MassData.txt'\n");
	  exit(1);
	}
	fprintf(mf,"%% structural mass data \n");
	fprintf(mf,"%% element\tAx\t\tlength\t\tdensity\t\tmass \n");
#endif

	sfrv=fscanf( fp, "%d", lump );
	if (sfrv != 1) sferr("lump value in mass data");
	sfrv=fscanf( fp, "%lf", tol );
	if (sfrv != 1) sferr("tol value in mass data");
	sfrv=fscanf( fp, "%lf", shift );
	if (sfrv != 1) sferr("shift value in mass data");
	sfrv=fscanf( fp, "%lf", exagg_modal );
	if (sfrv != 1) sferr("exagg_modal value in mass data");

	if (  lump_flag != -1   )	*lump = lump_flag;
	if (  tol_flag  != -1.0 )	*tol  = tol_flag;
	if ( shift_flag != -1.0 )	*shift = shift_flag;


	/* number of joints with extra inertias */
	sfrv=fscanf(fp,"%d", nI );
	if (sfrv != 1) sferr("nI value in mass data");
	if ( verbose ) {
		printf(" number of joints with extra lumped inertia ");
		dots(stdout,9);
		printf(" nI = %3d\n",*nI);
	}
	for (j=1; j <= *nI; j++) {
		sfrv=fscanf(fp, "%d", &jnt );
		if (sfrv != 1) sferr("joint value in extra joint mass data");
		if ( jnt < 1 || jnt > nJ ) {
	    		fprintf(stderr,"  error in joint mass data: joint number out of range  ");
	    		fprintf(stderr,"  Joint: %d  \n", jnt );
	    		fprintf(stderr,"  Perhaps you did not specify %d extra masses \n", *nI );
			fprintf(stderr,
			"  or perhaps the Input Data file is missing expected data.\n");
	    		exit(1);
		}
		sfrv=fscanf(fp, "%f %f %f %f",
			&JMs[jnt], &JMx[jnt], &JMy[jnt], &JMz[jnt] );
		if (sfrv != 4) sferr("joint inertia in extra mass data");
		*total_mass += JMs[jnt];

		if ( JMs[jnt]==0 && JMx[jnt]==0 && JMy[jnt]==0 && JMz[jnt]==0 )
	    	fprintf(stderr,"  warning: All extra joint inertia at joint %d  are zero\n", jnt );
	}

	/* number of frame elements with extra beam mass */
	sfrv=fscanf(fp,"%d", nX );
	if (sfrv != 1) sferr("nX value in mass data");
	if ( verbose ) {
		printf(" number of frame elements with extra mass ");
		dots(stdout,11);
		printf(" nX = %3d\n",*nX);
		if (sfrv != 1) sferr("element value in extra element mass data");
	}
	for (m=1; m <= *nX; m++) {
		sfrv=fscanf(fp, "%d", &b );
		if (sfrv != 1) sferr("element number in extra element mass data");
		if ( b < 1 || b > nE ) {
	    		fprintf(stderr,"  error in element mass data: element number out of range  ");
	    		fprintf(stderr,"  Element: %d  \n", b );
	    		fprintf(stderr,"  Perhaps you did not specify %d extra masses \n", *nX );
			fprintf(stderr,
			"  or perhaps the Input Data file is missing expected data.\n");
	    		exit(1);
		}
		sfrv=fscanf(fp, "%f", &BMs[b] );
		if (sfrv != 1) sferr("extra element mass value in mass data");
	}


	/* calculate the total mass and the structural mass */
	for (b=1; b <= nE; b++) {
		*total_mass  += d[b]*Ax[b]*L[b] + BMs[b];
		*struct_mass += d[b]*Ax[b]*L[b];
#ifdef MASSDATA_DEBUG
		fprintf(mf," %4d\t\t%12.5e\t%12.5e\t%12.5e\t%12.5e \n",
		 b, Ax[b], L[b], d[b], d[b]*Ax[b]*L[b] );
#endif
	}

#ifdef MASSDATA_DEBUG
	fclose(mf);
#endif

	for (m=1;m<=nE;m++) {			/* check inertia data	*/
	    if ( d[m] < 0.0 || BMs[m] < 0.0 || d[m]+BMs[m] <= 0.0 ) {
		fprintf(stderr,"  error: Non-positive mass or density\n");
		fprintf(stderr,"  d[%d]= %f  BMs[%d]= %f\n",m,d[m],m,BMs[m]);
		exit(1);
	    }
	}
/*	for (m=1;m<=nE;m++) ms += BMs[m]; // consistent mass doesn't agree  */
/*	if ( ms > 0.0 )	    *lump = 1;    // with concentrated masses, BMs  */

	if ( verbose ) {
		printf(" structural mass ");
		dots(stdout,36);
		printf("  %12.4e\n",*struct_mass);
		printf(" total mass ");
		dots(stdout,41);
		printf("  %12.4e\n",*total_mass);
	}
	sfrv=fscanf ( fp, "%d", &nA );
	if (sfrv != 1) sferr("nA value in mode animation data");
	if ( verbose ) {
		printf(" number of modes to be animated ");
		dots(stdout,21);
		printf(" nA = %3d\n",nA);
	}
	if (nA > 20)
	  fprintf(stderr," nA = %d, only 20 or fewer modes may be animated\n", nA );
	for ( m = 0; m < 20; m++ )	anim[m] = 0;
	for ( m = 0; m < nA; m++ ) {
		sfrv=fscanf ( fp, "%d", &anim[m] );
		if (sfrv != 1) sferr("mode number in mode animation data");
	}

	sfrv=fscanf ( fp, "%f", pan );
	if (sfrv != 1) sferr("pan value in mode animation data");
	if ( pan_flag != -1.0 )	*pan = pan_flag;


	strcpy(base_file,OUT_file);	
	while ( base_file[len++] != '\0' )
		/* the length of the base_file */ ;
	full_len = len;

	while ( base_file[len--] != '.' && len > 0 )
		/* find the last '.' in base_file */ ;
	if ( len == 0 )	len = full_len;
	base_file[++len] = '\0';	/* end base_file at the last '.' */

	while ( base_file[len] != '/' && base_file[len] != '\\' && len > 0 )
		len--;	/* find the last '/' or '\' in base_file */ 
	i = 0;
	while ( base_file[len] != '\0' )
		mode_file[i++] = base_file[len++];
	mode_file[i] = '\0';
	strcat(mode_file,"-m");
	output_path(mode_file,modepath,FRAME3DD_PATHMAX,NULL);

	return;
}


/*------------------------------------------------------------------------------
READ_CONDENSE   -  read matrix condensation information 	        30aug01
------------------------------------------------------------------------------*/
void read_condensation_data (
		FILE *fp,
		int nJ, int nM,
		int *nC, int *Cdof,
		int *Cmethod, int condense_flag, int *q, int *m, int verbose
){
	int	i,j,k,  **qm;
	int	sfrv;	/* *scanf return value */

	*Cmethod = *nC = *Cdof = 0;

	if ( (sfrv=fscanf ( fp, "%d", Cmethod )) != 1 )   {
		*Cmethod = *nC = *Cdof = 0;
		return;
	}

	if ( condense_flag != -1 )	*Cmethod = condense_flag;

	if ( *Cmethod <= 0 )  {
		*Cmethod = *nC = *Cdof = 0;
		return;
	}

	if ( *Cmethod > 3 ) *Cmethod = 1;	/* default */
	if ( verbose ) {
		printf(" condensation method ");
		dots(stdout,32);
		printf(" %d ", *Cmethod );
		if ( *Cmethod == 1 )	printf(" (static only) \n");
		if ( *Cmethod == 2 )	printf(" (Guyan) \n");
		if ( *Cmethod == 3 )	printf(" (dynamic) \n");
	}

	if ( (sfrv=fscanf ( fp, "%d", nC )) != 1 )  {
		*Cmethod = *nC = *Cdof = 0;
		return;
	}

	if ( verbose ) {
		printf(" number of joints with condensed DoF's ");
		dots(stdout,14);
		printf(" nC = %3d\n", *nC );
	}

	if ( (*nC) > nJ ) {
	  fprintf(stderr," error in matrix condensation data: \n");
	  fprintf(stderr,"  error: nC > nJ ... nC=%d; nJ=%d;\n",*nC,nJ);
	  fprintf(stderr,"  The number of joints with condensed DoF's ");
	  fprintf(stderr,"may not exceed the total number of joints.\n");
	  exit(1);
	}

	qm = imatrix( 1, *nC, 1,7 );

	for ( i=1; i <= *nC; i++) {
	 sfrv=fscanf( fp, "%d %d %d %d %d %d %d",
	 &qm[i][1],
	 &qm[i][2], &qm[i][3], &qm[i][4], &qm[i][5], &qm[i][6], &qm[i][7]);
	 if (sfrv != 7) sferr("DoF numbers in condensation data");
	 if ( qm[i][1] < 1 || qm[i][1] > nJ ) {		/* error check */
	  fprintf(stderr," error in matrix condensation data: ");
	  fprintf(stderr," condensed joint number out of range\n");
	  fprintf(stderr,"  cj[%d] = %d  ... nJ = %d  \n", i, qm[i][1], nJ );
	  exit(1);
	 }
	}

	for (i=1; i <= *nC; i++)  for (j=2; j<=7; j++)  if (qm[i][j]) (*Cdof)++;

	k=1;
	for (i=1; i <= *nC; i++) {
		for (j=2; j<=7; j++) {
			if (qm[i][j]) {
				q[k] = 6*(qm[i][1]-1) + j-1;
				++k;
			}
		}
	}

	for (i=1; i<= *Cdof; i++) {
	 sfrv=fscanf( fp, "%d", &m[i] );
	 if (sfrv != 1) sferr("mode number in condensation data");
	 if ( (m[i] < 0 || m[i] > nM) && *Cmethod == 3 ) {
	  fprintf(stderr," error in matrix condensation data: \n");
	  fprintf(stderr,"  error: m[%d] = %d \n",i,m[i]);
	  fprintf(stderr,"  The condensed mode number must be between ");
	  fprintf(stderr,"  1 and %d (modes).\n", nM);
	  exit(1);
	 }
	}

	free_imatrix(qm,1, *nC, 1,7);
	return;
}


/*------------------------------------------------------------------------------
WRITE_INPUT_DATA  -  save input data					07nov02
------------------------------------------------------------------------------*/
void write_input_data(
	FILE *fp,
	char *title, int nJ, int nE, int nL,
	int *nD, int nR,
	int *nF, int *nU, int *nW, int *nP, int *nT,
	vec3 *xyz, float *r,
	int *J1, int *J2,
	float *Ax, float *Asy, float *Asz, float *J, float *Iy, float *Iz,
	float *E, float *G, float *p, float *d, 
	float *gX, float *gY, float *gZ, 
	double **F, float **Dp,
	int *R,
	float ***U, float ***W, float ***P, float ***T,
	int shear, int anlyz, int geom
){
	int	i,j,n, lc;
	time_t  now;		/* modern time variable type    (DJGPP) */

	(void) time(&now);

	for (i=1; i<=80; i++)	fprintf(fp,"_");
  	fprintf(fp,"\nFrame3DD version: %s ", VERSION );
	fprintf(fp,"              http://frame3dd.sf.net/\n");
	fprintf(fp,"GPL Copyright (C) 1992-2009, Henri P. Gavin \n");
	fprintf(fp,"Frame3DD is distributed in the hope that it will be useful");
	fprintf(fp," but with no warranty.\n");
	fprintf(fp,"For details see the GNU Public Licence:");
	fprintf(fp," http://www.fsf.org/copyleft/gpl.html\n");
	for (i=1; i<=80; i++)	fprintf(fp,"_"); fprintf(fp,"\n\n");
	fprintf(fp,"%s\n",title);
	fprintf(fp, "%s", ctime(&now) );

	for (i=1; i<=80; i++)	fprintf(fp,"_");	fprintf(fp,"\n");

	fprintf(fp,"In 2D problems the Y-axis is vertical.  ");
#if Zvert
	fprintf(fp,"In 3D problems the Z-axis is vertical.\n");
#else
	fprintf(fp,"In 3D problems the Y-axis is vertical.\n");
#endif

	for (i=1; i<=80; i++)	fprintf(fp,"_");	fprintf(fp,"\n");

	fprintf(fp,"%5d JOINTS         ", nJ ); 
	fprintf(fp,"%5d FIXED JOINTS   ", nR );
	fprintf(fp,"%5d FRAME ELEMENTS ", nE ); 
	fprintf(fp,"%3d LOAD CASES   \n", nL );

	for (i=1; i<=80; i++)	fprintf(fp,"_"); fprintf(fp,"\n");

	fprintf(fp,"J O I N T   D A T A     ");
	fprintf(fp,"                                    R E S T R A I N T S\n");
	fprintf(fp,"  Joint      X              Y              Z");
	fprintf(fp,"         radius  Fx Fy Fz Mx My Mz\n");
	for (i=1; i<=nJ; i++) {
	 j = 6*(i-1);
	 fprintf(fp,"%5d %14.6f %14.6f %14.6f %8.3f  %2d %2d %2d %2d %2d %2d\n",
		i, xyz[i].x, xyz[i].y, xyz[i].z, r[i],
			R[j+1], R[j+2], R[j+3], R[j+4], R[j+5], R[j+6] );
	}
	fprintf(fp,"F R A M E   E L E M E N T   D A T A\t\t\t\t\t(local)\n");
	fprintf(fp,"  Elmnt  J1    J2     Ax   Asy   Asz    ");
	fprintf(fp,"Jxx     Iyy     Izz       E       G roll  density\n");
	for (i=1; i<= nE; i++) {
		fprintf(fp,"%5d %5d %5d %6.1f %5.1f %5.1f",
					i, J1[i],J2[i], Ax[i], Asy[i], Asz[i] );
		fprintf(fp," %6.1f %7.1f %7.1f %8.1f %7.1f %3.0f %8.2e\n",
			J[i], Iy[i], Iz[i], E[i], G[i], p[i]*180.0/PI, d[i] );
	}
	if ( shear )	fprintf(fp,"  Include shear deformations.\n");
	else		fprintf(fp,"  Neglect shear deformations.\n");
	if ( geom )	fprintf(fp,"  Include geometric stiffness.\n");
	else		fprintf(fp,"  Neglect geometric stiffness.\n");

	for (lc = 1; lc <= nL; lc++) {		/* start load case loop */

	  fprintf(fp,"\nL O A D   C A S E   %d   O F   %d  ... \n\n", lc,nL);
	  fprintf(fp,"   Gravity X = ");
	  if (gX[lc] == 0) fprintf(fp," 0.0 "); else fprintf(fp," %.3f ", gX[lc]);
	  fprintf(fp,"   Gravity Y = ");
	  if (gY[lc] == 0) fprintf(fp," 0.0 "); else fprintf(fp," %.3f ", gY[lc]);
	  fprintf(fp,"   Gravity Z = ");
	  if (gZ[lc] == 0) fprintf(fp," 0.0 "); else fprintf(fp," %.3f ", gZ[lc]);
	  fprintf(fp,"\n");
	  fprintf(fp," %3d concentrated loads\n", nF[lc] );
	  fprintf(fp," %3d uniformly distributed loads\n", nU[lc]);
	  fprintf(fp," %3d trapezoidally distributed loads\n", nW[lc]);
	  fprintf(fp," %3d concentrated point loads\n", nP[lc] );
	  fprintf(fp," %3d temperature loads\n", nT[lc] );
	  fprintf(fp," %3d prescribed displacements\n", nD[lc] );
	  if ( nF[lc] > 0 || nU[lc] > 0 || nW[lc] > 0 || nP[lc] > 0 || nT[lc] > 0 ) {
	    fprintf(fp," J O I N T   L O A D S");
	    fprintf(fp,"  +  E Q U I V A L E N T   J O I N T   L O A D S  (global)\n");
	    fprintf(fp,"  Joint       Fx          Fy          Fz");
	    fprintf(fp,"          Mxx         Myy         Mzz\n");
	    for (j=1; j<=nJ; j++) {
		i = 6*(j-1);
		if ( F[lc][i+1]!=0.0 || F[lc][i+2]!=0.0 || F[lc][i+3]!=0.0 ||
		     F[lc][i+4]!=0.0 || F[lc][i+5]!=0.0 || F[lc][i+6]!=0.0 ) {
			fprintf(fp, " %5d", j);
			for (i=5; i>=0; i--) fprintf(fp, " %11.3f", F[lc][6*j-i] );
			fprintf(fp, "\n");
		}
	    }
	  }

	  if ( nU[lc] > 0 ) {
	    fprintf(fp," U N I F O R M   L O A D S");
	    fprintf(fp,"\t\t\t\t\t\t(local)\n");
	    fprintf(fp,"  Elmnt       Ux               Uy               Uz\n");
	    for (n=1; n<=nU[lc]; n++) {
		fprintf(fp, " %5d", (int) (U[lc][n][1]) );
		for (i=2; i<=4; i++) fprintf(fp, " %16.8f", U[lc][n][i] );
		fprintf(fp, "\n");
	    }
	  }

	  if ( nW[lc] > 0 ) {
	    fprintf(fp," T R A P E Z O I D A L   L O A D S");
	    fprintf(fp,"\t\t\t\t\t(local)\n");
	    fprintf(fp,"  Elmnt       x1               x2               W1               W2\n");
	    for (n=1; n<=nW[lc]; n++) {
		fprintf(fp, " %5d", (int) (W[lc][n][1]) );
		for (i=2; i<=5; i++) fprintf(fp, " %16.8f", W[lc][n][i] );
		fprintf(fp, "  (x)\n");
		fprintf(fp, " %5d", (int) (W[lc][n][1]) );
		for (i=6; i<=9; i++) fprintf(fp, " %16.8f", W[lc][n][i] );
		fprintf(fp, "  (y)\n");
		fprintf(fp, " %5d", (int) (W[lc][n][1]) );
		for (i=10; i<=13; i++) fprintf(fp, " %16.8f", W[lc][n][i] );
		fprintf(fp, "  (z)\n");
	    }
	  }

	  if ( nP[lc] > 0 ) {
	    fprintf(fp," C O N C E T R A T E D   P O I N T   L O A D S");
	    fprintf(fp,"\t\t\t\t(local)\n");
	    fprintf(fp,"  Elmnt       Px          Py          Pz          x\n");
	    for (n=1; n<=nP[lc]; n++) {
		fprintf(fp, " %5d", (int) (P[lc][n][1]) );
		for (i=2; i<=5; i++) fprintf(fp, " %11.3f", P[lc][n][i] );
		fprintf(fp, "\n");
	    }
	  }

	  if ( nT[lc] > 0 ) {
	    fprintf(fp," T E M P E R A T U R E   C H A N G E S");
	    fprintf(fp,"\t\t\t\t\t(local)\n");
	    fprintf(fp,"  Elmnt     coef      hy        hz");
	    fprintf(fp,"        Ty+       Ty-       Tz+       Tz-\n");
	    for (n=1; n<=nT[lc]; n++) {
		fprintf(fp, " %5d", (int) (T[lc][n][1]) );
		fprintf(fp, " %9.2e", T[lc][n][2] );
		for (i=3; i<=8; i++) fprintf(fp, " %9.3f", T[lc][n][i] );
		fprintf(fp, "\n");
	    }
	  }

	  if ( nD[lc] > 0 ) {
	    fprintf(fp,"\n P R E S C R I B E D   D I S P L A C E M E N T S");
	    fprintf(fp,"                        (global)\n");
	    fprintf(fp,"  Joint       Dx          Dy          Dz");
	    fprintf(fp,"          Dxx         Dyy         Dzz\n");
	    for (j=1; j<=nJ; j++) {
		i = 6*(j-1);
		if ( Dp[lc][i+1]!=0.0 || Dp[lc][i+2]!=0.0 || Dp[lc][i+3]!=0.0 ||
		     Dp[lc][i+4]!=0.0 || Dp[lc][i+5]!=0.0 || Dp[lc][i+6]!=0.0 ){
			fprintf(fp, " %5d", j);
			for (i=5; i>=0; i--) fprintf(fp, " %11.3f",
							Dp[lc][6*j-i] );
			fprintf(fp, "\n");
		}
	    }
	  }

	}					/* end load case loop	*/

	if (anlyz) {
	 fprintf(fp,"\nE L A S T I C   S T I F F N E S S   A N A L Y S I S");
	 fprintf(fp,"   via  L D L'  decomposition\n\n");
	}
	else		fprintf(fp,"D A T A   C H E C K   O N L Y\n");
	fflush(fp);
	return;
}


/*------------------------------------------------------------------------------
WRITE_STATIC_RESULTS -  save joint displacements and frame element end forces
09 Sep 2008
------------------------------------------------------------------------------*/
void write_static_results (
		FILE *fp,
		int nJ, int nE, int nL, int lc, int DoF,
		int *J1, int *J2,
		double *F, double *D, int *R, double **Q,
		double err, int ok
){
	double	disp;
	int	i,j,n;

	if ( ok < 0 ) {
	 fprintf(fp,"  * The Stiffness Matrix is not positive-definite *\n");
	 fprintf(fp,"    Check that all six rigid-body translations are restrained\n");
	 fprintf(fp,"    If geometric stiffness is included, reduce the loads.\n");
/*	 return; */
	}

	fprintf(fp,"\nL O A D   C A S E   %d   O F   %d  ... \n\n", lc, nL);

	fprintf(fp,"J O I N T   D I S P L A C E M E N T S");
	fprintf(fp,"\t\t\t\t\t(global)\n");
	fprintf(fp,"  Joint   X-dsp       Y-dsp       Z-dsp");
	fprintf(fp,"       X-rot       Y-rot       Z-rot\n");
	for (j=1; j<= nJ; j++) {
	    disp = 0.0;
	    for ( i=5; i>=0; i-- ) disp += fabs( D[6*j-i] );
	    if ( disp > 0.0 ) {
		fprintf(fp," %5d", j);
		for ( i=5; i>=0; i-- ) {
			if ( fabs(D[6*j-i]) < 1.e-8 )
				fprintf (fp, "    0.0     ");
			else    fprintf (fp, " %11.6f",  D[6*j-i] );
		}
		fprintf(fp,"\n");
	    }
	}
	fprintf(fp,"F R A M E   E L E M E N T   E N D   F O R C E S");
	fprintf(fp,"\t\t\t\t(local)\n");
	fprintf(fp,"  Elmnt  Joint      Nx          Vy         Vz");
	fprintf(fp,"        Txx        Myy        Mzz\n");
	for (n=1; n<= nE; n++) {
		fprintf(fp," %5d  %5d", n, J1[n]);
		if ( fabs(Q[n][1]) < 0.0001 )
			fprintf (fp, "      0.0   ");
		else    fprintf (fp, " %10.3f", Q[n][1] );
		if ( Q[n][1] >=  0.0001 ) fprintf(fp, "c");
		if ( Q[n][1] <= -0.0001 ) fprintf(fp, "t");
		for (i=2; i<=6; i++) {
			if ( fabs(Q[n][i]) < 0.0001 )
				fprintf (fp, "      0.0  ");
			else    fprintf (fp, " %10.3f", Q[n][i] );
		}
		fprintf(fp,"\n");
		fprintf(fp," %5d  %5d", n, J2[n]);
		if ( fabs(Q[n][7]) < 0.0001 )
			fprintf (fp, "      0.0   ");
		else    fprintf (fp, " %10.3f", Q[n][7] );
		if ( Q[n][7] >=  0.0001 ) fprintf(fp, "t");
		if ( Q[n][7] <= -0.0001 ) fprintf(fp, "c");
		for (i=8; i<=12; i++) {
			if ( fabs(Q[n][i]) < 0.0001 )
				fprintf (fp, "      0.0  ");
			else    fprintf (fp, " %10.3f", Q[n][i] );
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"R E A C T I O N S\t\t\t\t\t\t\t(global)\n");
	fprintf(fp,"  Joint       Fx          Fy          Fz");
	fprintf(fp,"         Mxx         Myy         Mzz\n");
	for (j=1; j<=nJ; j++) {
		i = 6*(j-1);
		if ( R[i+1] || R[i+2] || R[i+3] ||
		     R[i+4] || R[i+5] || R[i+6] ) {
			fprintf(fp, " %5d", j);
			for (i=5; i>=0; i--) {
				if ( !R[6*j-i] || fabs(F[6*j-i]) < 0.0001 )
					fprintf (fp, "       0.0  ");
				else    fprintf (fp, " %11.3f", -F[6*j-i] );
			}
			fprintf(fp, "\n");
		}
	}
	fprintf(fp,"R M S   E Q U I L I B R I U M    E R R O R: %9.3e\n", err );
	fflush(fp);
	return;
}


/*------------------------------------------------------------------------------
WRITE_STATIC_CSV -  save joint displacements and frame element end forces
31 Dec 2008
------------------------------------------------------------------------------*/
void write_static_csv(
		char *OUT_file,
		char *title,
		int nJ, int nE, int nL, int lc, int DoF,
		int *J1, int *J2,
		double *F, double *D, int *R, double **Q,
		double err, int ok
){

	FILE	*fpcsv;
	int	i,j,n;
	char	*wa;
	char	CSV_file[128];
	time_t  now;		/* modern time variable type    (DJGPP) */

	(void) time(&now);

	i=0;
	j=0;
	while (i<128) {
		CSV_file[j] = OUT_file[i];
		if ( CSV_file[j] == '+' ||
		     CSV_file[j] == '-' ||
		     CSV_file[j] == '*' ||
		     CSV_file[j] == '^' ||
		     CSV_file[j] == '.' ||
		     CSV_file[j] == '\0') {
			CSV_file[j] = '_';
			break;
		}
		i++;
		j++;
	}
	CSV_file[++j] = '\0';
	strcat(CSV_file,"out.CSV");


	wa  = "a";
	if (lc == 1) wa = "w";

	if ((fpcsv = fopen (CSV_file, wa)) == NULL) {
	  fprintf (stderr," error: cannot open CSV file %s\n", CSV_file);
	  exit(1);
	}


	if ( lc == 1 ) {
  	 fprintf(fpcsv,"\" Frame3DD version: %s ", VERSION );
	 fprintf(fpcsv,"              http://frame3dd.sf.net/\"\n");
	 fprintf(fpcsv,"\"GPL Copyright (C) 1992-2009, Henri P. Gavin \"\n");
	 fprintf(fpcsv,"\"Frame3DD is distributed in the hope that it will be useful");
	 fprintf(fpcsv," but with no warranty.\"\n");
	 fprintf(fpcsv,"\"For details see the GNU Public Licence:");
	 fprintf(fpcsv," http://www.fsf.org/copyleft/gpl.html\"\n");
	 fprintf(fpcsv,"\" %s \"\n",title);
	 fprintf(fpcsv,"\" %s \"\n", ctime(&now) );

	 fprintf(fpcsv,"\" .CSV formatted results of Frame3DD analysis \"\n");
	 fprintf(fpcsv,"\n , Load Case , Displacements , End Forces , Reactions \n");
	 for (i = 1; i <= nL; i++) {
	 	fprintf(fpcsv," First Row , %d , %d , %d , %d  \n",
			i,
			15+(i-1)*(nJ*2+nE*2+10) + 2*nL,
			17+(i-1)*(nJ*2+nE*2+10) + 2*nL+ nJ,
			19+(i-1)*(nJ*2+nE*2+10) + 2*nL+ nJ + 2*nE );
	 	fprintf(fpcsv," Last Row , %d , %d , %d , %d  \n",
			i,
			15+(i-1)*(nJ*2+nE*2+10) + 2*nL + nJ - 1,
			17+(i-1)*(nJ*2+nE*2+10) + 2*nL + nJ + 2*nE - 1,
			19+(i-1)*(nJ*2+nE*2+10) + 2*nL + 2*nJ + 2*nE - 1 );
	 }

	}


	if ( ok < 0 ) {
	 fprintf(fpcsv,"\"  * The Stiffness Matrix is not positive-definite * \"\n");
	 fprintf(fpcsv,"\" Check that all six rigid-body translations are restrained\"\n");
	 fprintf(fpcsv,"\" If geometric stiffness is included, reduce the loads.\"\n");
/*	 return; */
	}


	fprintf(fpcsv,"\n\"L O A D   C A S E   %d   O F   %d  ... \"\n\n", lc, nL);

	fprintf(fpcsv,"\"J O I N T   D I S P L A C E M E N T S");
	fprintf(fpcsv,"  (global)\"\n");
	fprintf(fpcsv,"Joint ,  X-dsp   ,   Y-dsp  ,    Z-dsp");
	fprintf(fpcsv," ,     X-rot  ,    Y-rot   ,   Z-rot\n");
	for (j=1; j<= nJ; j++) {
		fprintf(fpcsv," %5d,", j);
		for ( i=5; i>=0; i-- ) {
			if ( fabs(D[6*j-i]) < 1.e-8 )
				fprintf (fpcsv, "    0.0,    ");
			else    fprintf (fpcsv, " %12.5e,",  D[6*j-i] );
		}
		fprintf(fpcsv,"\n");
	}
	fprintf(fpcsv,"\"F R A M E   E L E M E N T   E N D   F O R C E S");
	fprintf(fpcsv,"  (local)\"\n");
	fprintf(fpcsv,"Elmnt , Joint ,    Nx     ,    Vy   ,     Vz");
	fprintf(fpcsv,"   ,     Txx   ,    Myy  ,     Mzz\n");
	for (n=1; n<= nE; n++) {
		fprintf(fpcsv," %5d, %5d,", n, J1[n]);
		if ( fabs(Q[n][1]) < 0.0001 )
			fprintf (fpcsv, "      0.0,  ");
		else    fprintf (fpcsv, " %12.5e,", Q[n][1] );
		for (i=2; i<=6; i++) {
			if ( fabs(Q[n][i]) < 0.0001 )
				fprintf (fpcsv, "      0.0, ");
			else    fprintf (fpcsv, " %12.5e,", Q[n][i] );
		}
		fprintf(fpcsv,"\n");
		fprintf(fpcsv," %5d, %5d,", n, J2[n]);
		if ( fabs(Q[n][7]) < 0.0001 )
			fprintf (fpcsv, "      0.0,  ");
		else    fprintf (fpcsv, " %12.5e,", Q[n][7] );
		for (i=8; i<=12; i++) {
			if ( fabs(Q[n][i]) < 0.0001 )
				fprintf (fpcsv, "      0.0, ");
			else    fprintf (fpcsv, " %12.5e,", Q[n][i] );
		}
		fprintf(fpcsv,"\n");
	}
	fprintf(fpcsv,"\"R E A C T I O N S  (global)\"\n");
	fprintf(fpcsv," Joint  ,    Fx      ,   Fy   ,      Fz");
	fprintf(fpcsv,"   ,     Mxx    ,    Myy    ,    Mzz\n");
	for (j=1; j<=nJ; j++) {
		i = 6*(j-1);
		fprintf(fpcsv, " %5d,", j);
		for (i=5; i>=0; i--) {
			if ( !R[6*j-i] || fabs(F[6*j-i]) < 0.0001 )
				fprintf (fpcsv, "       0.0, ");
			else    fprintf (fpcsv, " %12.5e,", -F[6*j-i] );
		}
		fprintf(fpcsv, "\n");
	}
	fprintf(fpcsv,"\"R M S   E Q U I L I B R I U M    E R R O R:\", %9.3e\n", err );

	fclose(fpcsv);

	return;
}

/*------------------------------------------------------------------------------
WRITE_STATIC_MFILE -  						
save joint displacements and frame element end forces in an m-file
this function interacts with frame_3dd.m, an m-file interface to frame3dd
09 Sep 2008
------------------------------------------------------------------------------*/
void write_static_mfile (
		char *OUT_file, char *title,
		int nJ, int nE, int nL, int lc, int DoF,
		int *J1, int *J2,
		double *F, double *D, int *R, double **Q,
		double err, int ok
){
	FILE	*fpm;
	int	i,j,n;
	char	*wa;
	char	M_file[128];
	time_t  now;	/* modern time variable type    (DJGPP) */

	(void) time(&now);

	i=0;
	j=0;
	while (i<128) {
		M_file[j] = OUT_file[i];
		if ( M_file[j] == '+' ||
		     M_file[j] == '-' ||
		     M_file[j] == '*' ||
		     M_file[j] == '^' ||
		     M_file[j] == '.' ||
		     M_file[j] == '\0') {
			M_file[j] = '_';
			break;
		}
		i++;
		j++;
	}
	M_file[++j] = '\0';
	strcat(M_file,"out.m");

	wa  = "a";
	if (lc == 1) wa = "w";

	if ((fpm = fopen (M_file, wa)) == NULL) {
	  fprintf (stderr," error: cannot open file %s\n", M_file );
	  exit(1);
	}

	if ( lc == 1 ) {
  	 fprintf(fpm,"%% Frame3DD version: %s ", VERSION );
	 fprintf(fpm,"              http://frame3dd.sf.net/\n");
	 fprintf(fpm,"%%GPL Copyright (C) 1992-2009, Henri P. Gavin \n");
	 fprintf(fpm,"%%Frame3DD is distributed in the hope that it will be useful");
	 fprintf(fpm," but with no warranty.\n");
	 fprintf(fpm,"%%For details see the GNU Public Licence:");
	 fprintf(fpm," http://www.fsf.org/copyleft/gpl.html\n");
	 fprintf(fpm,"%% %s\n",title);
	 fprintf(fpm, "%% %s", ctime(&now) );

	 fprintf(fpm,"%% m-file formatted results of frame3dd analysis\n");
	 fprintf(fpm,"%% to be read by frame_3dd.m\n");
	}


	if ( ok < 0 ) {
	 fprintf(fpm,"%%  The Stiffness Matrix is not positive-definite *\n");
	 fprintf(fpm,"%%  Check that all six rigid-body translations are restrained\n");
	 fprintf(fpm,"%%  If geometric stiffness is included, reduce the loads.\n");
/*	 return; */
	}

	fprintf(fpm,"\n%% L O A D   C A S E   %d   O F   %d  ... \n\n", lc, nL);

	fprintf(fpm,"%% J O I N T   D I S P L A C E M E N T S");
	fprintf(fpm,"\t\t(global)\n");
	fprintf(fpm,"%%\tX-dsp\t\tY-dsp\t\tZ-dsp\t\tX-rot\t\tY-rot\t\tZ-rot\n");
	fprintf(fpm,"D%d=[",lc);
	for (j=1; j<= nJ; j++) {
		for ( i=5; i>=0; i-- ) {
			if ( fabs(D[6*j-i]) < 1.e-8 )
				fprintf (fpm, "\t0.0\t");
			else    fprintf (fpm, "\t%13.6e",  D[6*j-i] );
		}
		if ( j < nJ )	fprintf(fpm," ; \n");
		else		fprintf(fpm," ]'; \n\n");
	}

	fprintf(fpm,"%% F R A M E   E L E M E N T   E N D   F O R C E S");
	fprintf(fpm,"\t\t(local)\n");
	fprintf(fpm,"%%\tNx_1\t\tVy_1\t\tVz_1\t\tTxx_1\t\tMyy_1\t\tMzz_1\t");
	fprintf(fpm,"  \tNx_2\t\tVy_2\t\tVz_2\t\tTxx_2\t\tMyy_2\t\tMzz_2\n");
	fprintf(fpm,"F%d=[",lc);
	for (n=1; n<= nE; n++) {
		if ( fabs(Q[n][1]) < 0.0001 )
			fprintf (fpm, "\t0.0\t");
		else    fprintf (fpm, "\t%13.6e", Q[n][1] );
		for (i=2; i<=6; i++) {
			if ( fabs(Q[n][i]) < 0.0001 )
				fprintf (fpm, "\t0.0\t");
			else    fprintf (fpm, "\t%13.6e", Q[n][i] );
		}
		if ( fabs(Q[n][7]) < 0.0001 )
			fprintf (fpm, "\t0.0\t");
		else    fprintf (fpm, "\t%13.6e", Q[n][7] );
		for (i=8; i<=12; i++) {
			if ( fabs(Q[n][i]) < 0.0001 )
				fprintf (fpm, "\t0.0\t");
			else    fprintf (fpm, "\t%13.6e", Q[n][i] );
		}
		if ( n < nE )	fprintf(fpm," ; \n");
		else		fprintf(fpm," ]'; \n\n");
	}

	fprintf(fpm,"%% R E A C T I O N S\t\t\t\t(global)\n");
	fprintf(fpm,"%%\tFx\t\tFy\t\tFz\t\tMxx\t\tMyy\t\tMzz\n");
	fprintf(fpm,"R%d=[",lc);
	for (j=1; j<=nJ; j++) {
		i = 6*(j-1);
		for (i=5; i>=0; i--) {
			if ( !R[6*j-i] || fabs(F[6*j-i]) < 0.0001 )
				fprintf (fpm, "\t0.0\t");
			else    fprintf (fpm, "\t%13.6e", -F[6*j-i] );
		}
		if ( j < nJ )	fprintf(fpm," ; \n");
		else		fprintf(fpm," ]'; \n\n");
	}

	fprintf(fpm,"%% R M S   E Q U I L I B R I U M    E R R O R: %9.3e\n", err );
	fprintf(fpm,"\n\n  load Ks \n\n");

	fclose(fpm);

	return;
}


/*------------------------------------------------------------------------------
WRITE_MODAL_RESULTS -  save modal frequencies and mode shapes	
16 Aug 2001
------------------------------------------------------------------------------*/
void write_modal_results(
		FILE *fp,
		int nJ, int nE, int nI, int DoF,
		double **M, double *f, double **V,
		double total_mass, double struct_mass,
		int iter, int sumR, int nM,
		double shift, int lump, double tol, int ok
){
	int	i, j, k, m, num_modes;
	double	mpfX, mpfY, mpfZ,	/* mode participation factors	*/
		*msX, *msY, *msZ;
	double	fs;

	msX = dvector(1,DoF);
	msY = dvector(1,DoF);
	msZ = dvector(1,DoF);

	for (i=1; i<=DoF; i++) {
		msX[i] = msY[i] = msZ[i] = 0.0;
		for (j=1; j<=DoF; j+=6) msX[i] += M[i][j];
		for (j=2; j<=DoF; j+=6) msY[i] += M[i][j];
		for (j=3; j<=DoF; j+=6) msZ[i] += M[i][j];
	}

	if ( (DoF - sumR) > nM )	num_modes = nM;
	else	num_modes = DoF - sumR;

	fprintf(fp,"\nM O D A L   A N A L Y S I S   R E S U L T S\n");
	fprintf(fp,"  Total Mass:  %e   ", total_mass );
	fprintf(fp,"  Structural Mass:  %e \n", struct_mass );
	fprintf(fp,"J O I N T   M A S S E S");
	fprintf(fp,"\t(diagonal of the mass matrix)\t\t\t(global)\n");
	fprintf(fp,"  Joint X-mass      Y-mass      Z-mass");
	fprintf(fp,"      X-inrta     Y-inrta     Z-inrta\n");
	for (j=1; j <= nJ; j++) {
		k = 6*(j-1);
		fprintf(fp," %5d", j);
		for ( i=1; i<=6; i++ )
			fprintf (fp, " %11.5e", M[k+i][k+i] );
		fprintf(fp,"\n");
	}
	if ( lump )	fprintf(fp,"  Lump masses at joints.\n");
	else		fprintf(fp,"  Use consistent mass matrix.\n");
	fprintf(fp,"N A T U R A L   F R E Q U E N C I E S   & \n");
	fprintf(fp,"M A S S   N O R M A L I Z E D   M O D E   S H A P E S \n");
	fprintf(fp," convergence tolerance: %.3e \n", tol);
	for (m=1; m<=num_modes; m++) {
	    mpfX = 0.0;	for (i=1; i<=DoF; i++)    mpfX += V[i][m]*msX[i];
	    mpfY = 0.0;	for (i=1; i<=DoF; i++)    mpfY += V[i][m]*msY[i];
	    mpfZ = 0.0;	for (i=1; i<=DoF; i++)    mpfZ += V[i][m]*msZ[i];
	    fprintf(fp,"  MODE %5d:   f= %lf Hz,  T= %lf sec\n",m,f[m],1./f[m]);
	    fprintf(fp,"\t\tX- modal participation factor = %12.4e \n", mpfX);
	    fprintf(fp,"\t\tY- modal participation factor = %12.4e \n", mpfY);
	    fprintf(fp,"\t\tZ- modal participation factor = %12.4e \n", mpfZ);

	    fprintf(fp,"  Joint   X-dsp       Y-dsp       Z-dsp");
	    fprintf(fp,"       X-rot       Y-rot       Z-rot\n");
	    for (j=1; j<= nJ; j++) {
		fprintf(fp," %5d", j);
		for ( i=5; i>=0; i-- )	fprintf (fp, " %11.3e", V[6*j-i][m] );
		fprintf(fp,"\n");
	    }
	}

	fprintf(fp,"M A T R I X    I T E R A T I O N S: %d\n", iter );

	fs = sqrt(4.0*PI*PI*f[nM]*f[nM] + tol) / (2.0*PI);

	fprintf(fp,"There are %d modes below %f Hz.", -ok, fs );
	if ( -ok > nM ) {
		fprintf(fp," ... %d modes were not found.\n", -ok-nM );
		fprintf(fp," Try increasing the number of modes in \n");
		fprintf(fp," order to get the missing modes below %f Hz.\n",fs);
	} else  fprintf(fp," ... All %d modes were found.\n", nM );


	free_dvector(msX,1,DoF);
	free_dvector(msY,1,DoF);
	free_dvector(msZ,1,DoF);
	fflush(fp);
	return;
}


/*------------------------------------------------------------------------------
STATIC_MESH  -
create mesh data of deformed and undeformed mesh, use gnuplot	22 Feb 1999
useful gnuplot options: set noxtics noytics noztics noborder view nokey
------------------------------------------------------------------------------*/
void static_mesh(
		char IN_file[], char meshpath[], char plotpath[],
		char *title, int nJ, int nE, int nL, int lc, int DoF,
		vec3 *xyz, double *L,
		int *J1, int *J2, float *p, double *D,
		double exagg_static, int anlyz
){
	FILE	*fpmfx, *fpm;
	double	mx, my, mz;	/* coordinates of the frame element number labels */
	int	j1, j2, i, j, m, X=0, Y=0, Z=0;
	char	meshfl[128], D3 = '#';
	time_t  now;		/* modern time variable type    (DJGPP) */

	(void) time(&now);

	sprintf( meshfl, "%sf.%03d", meshpath, lc );

	if ((fpmfx = fopen (meshfl, "w")) == NULL) {
		printf (" error: cannot open meshpath: %s\n", meshfl );
		exit(1);
	}

	if ((fpm = fopen (meshpath, "w")) == NULL) {
		printf (" error: cannot open meshpath: %s\n", meshpath );
		exit(1);
	}

	if (!anlyz) exagg_static = 0.0;



	fprintf(fpm,"# FRAME3DD ANALYSIS RESULTS  http://frame3dd.sf.net/\n");
	fprintf(fpm,"# %s\n", title );
	fprintf(fpm,"# L O A D  C A S E   %d  of   %d \n", lc, nL );
	fprintf(fpm,"# %s", ctime(&now) );
	fprintf(fpm,"# M E S H   D A T A   (global coordinates)");
	fprintf(fpm," deflection exaggeration: %.1f\n", exagg_static );
	fprintf(fpm,"# Joint      X           Y           Z");
	fprintf(fpm,"          X-dsp       Y-dsp       Z-dsp\n");

	fprintf(fpmfx,"# FRAME3DD ANALYSIS RESULTS  http://frame3dd.sf.net/\n");
	fprintf(fpmfx,"# %s\n", title );
	fprintf(fpmfx,"# %s", ctime(&now) );
	fprintf(fpmfx,"# F L E X E D   M E S H   D A T A ");
	fprintf(fpmfx,"  deflection exaggeration: %.1f\n", exagg_static );
	fprintf(fpmfx,"#       X-dsp        Y-dsp        Z-dsp\n");

	for (m=1; m<=nE; m++) {

		bent_beam ( fpmfx, J1[m],J2[m],xyz, L[m],p[m],D, exagg_static );

		j = J1[m];	i = 6*(j-1);
		fprintf (fpm,"%5d %11.3e %11.3e %11.3e", j, xyz[j].x,xyz[j].y,xyz[j].z);
		fprintf (fpm," %11.3e %11.3e %11.3e\n",
			xyz[j].x + exagg_static*D[i+1],
			xyz[j].y + exagg_static*D[i+2],
			xyz[j].z + exagg_static*D[i+3]
		);

		j = J2[m];	i = 6*(j-1);
		fprintf (fpm,"%5d %11.3e %11.3e %11.3e", j, xyz[j].x,xyz[j].y,xyz[j].z);
		fprintf (fpm," %11.3e %11.3e %11.3e\n",
			xyz[j].x + exagg_static*D[i+1],
			xyz[j].y + exagg_static*D[i+2],
			xyz[j].z + exagg_static*D[i+3]
		);
		fprintf(fpm,"\n\n");
	}

	for ( j=1; j<=nJ; j++ ) {
		if (xyz[j].x != 0.0) X=1;	/* check for three-dimensional frame */
		if (xyz[j].y != 0.0) Y=1;
		if (xyz[j].z != 0.0) Z=1;
	}
	if ( X && Y && Z ) D3 = ' ';

	fclose(fpmfx);
	fclose(fpm);

	if (lc == 1) {
	    if ((fpm = fopen (plotpath, "w")) == NULL) {
		printf (" error: cannot open plot file: %s\n", plotpath);
		exit(1);
	    }
	} else {
	    if ((fpm = fopen (plotpath, "a")) == NULL) {
		printf (" error: cannot open plot file: %s\n", plotpath);
		exit(1);
	    }
	}

	if (lc == 1) {		/* first load case */

	 fprintf(fpm,"# FRAME3DD ANALYSIS RESULTS  http://frame3dd.sf.net/\n");
	 fprintf(fpm,"# %s\n", title );
	 fprintf(fpm,"# %s", ctime(&now) );
	 fprintf(fpm,"# M E S H   A N N O T A T I O N   F I L E \n");

	 fprintf(fpm,"set title \"%s\\n", title );
	 fprintf(fpm,"analysis file: %s ", IN_file );
	 fprintf(fpm,"  deflection exaggeration: %.1f ", exagg_static );
	 fprintf(fpm,"  load case %d of %d \"\n", lc, nL );

	 fprintf(fpm,"set autoscale\n");
	 fprintf(fpm,"set noborder\n");
	 fprintf(fpm,"set pointsize 1.0\n");
	 fprintf(fpm,"set xtics; set ytics; set ztics; \n");
	 fprintf(fpm,"set nozeroaxis\n");
	 fprintf(fpm,"set nokey\n");
	 fprintf(fpm,"set nolabel\n");

 	 fprintf(fpm,"# NODE NUMBER LABELS\n");
	 for (j=1; j<=nJ; j++)
		fprintf(fpm,"set label ' %d' at %12.4e, %12.4e, %12.4e\n",
					j, xyz[j].x,xyz[j].y,xyz[j].z );

	 fprintf(fpm,"# ELEMENT NUMBER LABELS\n");
	 for (m=1; m<=nE; m++) {
		j1 = J1[m];	j2 = J2[m];
		mx = 0.5 * ( xyz[j1].x + xyz[j2].x );
		my = 0.5 * ( xyz[j1].y + xyz[j2].y );
		mz = 0.5 * ( xyz[j1].z + xyz[j2].z );
		fprintf(fpm,"set label ' %d' at %12.4e, %12.4e, %12.4e\n",
								m, mx, my, mz );
	 }
	 fprintf(fpm,"plot '%s' u 2:3 t 'undeformed mesh' w lp ", meshpath);
	 if (!anlyz) fprintf(fpm,"lw 2 lt 1 pt 6 \n");
	 else fprintf(fpm,"lw 1 lt 5 pt 6, '%s' u 1:2 t 'load case %d of %d' w l lw 2 lt 3\n", meshfl, lc, nL );

	 fprintf(fpm,"%c set parametric\n", D3 );
	 fprintf(fpm,"%c set view 60, 70, 1 \n", D3 );
	 fprintf(fpm,"%c set nokey\n", D3 );
	 fprintf(fpm,"%c set xlabel 'x'\n", D3 );
	 fprintf(fpm,"%c set ylabel 'y'\n", D3 );
	 fprintf(fpm,"%c set zlabel 'z'\n", D3 );
/*	 fprintf(fpm,"%c set nolabel\n", D3 );	*/
	 fprintf(fpm,"%c splot '%s' u 2:3:4 t 'load case %d of %d' w lp ",
							D3, meshpath, lc, nL );
	 if (!anlyz) fprintf(fpm," lw 2 lt 1 pt 6 \n");
	 else fprintf(fpm," lw 1 lt 5 pt 6, '%s' u 1:2:3 t 'load case %d of %d' w l lw 2 lt 3\n",meshfl, lc, nL );

	} else { 		/* additional load cases */

	 fprintf(fpm,"pause -1\n");

	 fprintf(fpm,"set title \"%s\\n", title );
	 fprintf(fpm,"analysis file: %s ", IN_file );
	 fprintf(fpm,"  deflection exaggeration: %.1f ", exagg_static );
	 fprintf(fpm,"  load case %d of %d \"\n", lc, nL );

	 fprintf(fpm,"plot '%s' u 2:3 t 'undeformed mesh' w lp ", meshpath);
	 if (!anlyz) fprintf(fpm,"lw 2 lt 1 pt 6 \n");
	 else fprintf(fpm,"lw 1 lt 5 pt 6, '%s' u 1:2 t 'load case %d of %d' w l lw 2 lt 3\n", meshfl, lc, nL );

	 fprintf(fpm,"%c set parametric\n", D3 );
	 fprintf(fpm,"%c set view 60, 70, 1 \n", D3 );
	 fprintf(fpm,"%c set nokey\n", D3 );
	 fprintf(fpm,"%c set xlabel 'x'\n", D3 );
	 fprintf(fpm,"%c set ylabel 'y'\n", D3 );
	 fprintf(fpm,"%c set zlabel 'z'\n", D3 );
/*	 fprintf(fpm,"%c set nolabel\n", D3 );	*/
	 fprintf(fpm,"%c splot '%s' u 2:3:4 t 'undeformed mesh' w lp ",
								D3, meshpath );
	 if (!anlyz) fprintf(fpm," lw 2 lt 1 pt 6 \n");
	 else fprintf(fpm," lw 1 lt 5 pt 6, '%s' u 1:2:3 t 'load case %d of %d' w l lw 2 lt 3\n",meshfl, lc, nL );
	}

	fclose(fpm);

	return;
}


/*------------------------------------------------------------------------------
MODAL_MESH  -  create mesh data of the mode-shape meshes, use gnuplot	19oct98
	 useful gnuplot options: set noxtics noytics noztics noborder view nokey
------------------------------------------------------------------------------*/
void modal_mesh(
		char IN_file[], char meshpath[], char modepath[],
		char plotpath[], char *title,
		int nJ, int nE, int DoF, int nM,
		vec3 *xyz, double *L,
		int *J1, int *J2, float *p,
		double **M, double *f, double **V,
		double exagg_modal, int anlyz
){
	FILE	*fpm;
	double mpfX, mpfY, mpfZ;	/* mode participation factors	*/
	double *msX, *msY, *msZ;
	double *v;		/* a mode-shape vector */

	int	i, j, m,n, X=0, Y=0, Z=0;
	char	D3 = '#', modefl[128];


	msX = dvector(1,DoF);
	msY = dvector(1,DoF);
	msZ = dvector(1,DoF);
	v   = dvector(1,DoF);

	for (i=1; i<=DoF; i++) {	/* modal participation factors */
		msX[i] = msY[i] = msZ[i] = 0.0;
		for (j=1; j<=DoF; j+=6) msX[i] += M[i][j];
		for (j=2; j<=DoF; j+=6) msY[i] += M[i][j];
		for (j=3; j<=DoF; j+=6) msZ[i] += M[i][j];
	}

	if (!anlyz) exagg_modal = 0.0;

	for (m=1; m<=nM; m++) {

		sprintf( modefl,"%s-%02d-", modepath, m );

		if ((fpm = fopen (modefl, "w")) == NULL) {
			printf (" error: cannot open modal mesh file: %s\n", modefl);
			exit(1);
		}

		fprintf(fpm,"# FRAME3DD ANALYSIS RESULTS  http://frame3dd.sf.net/\n");
		fprintf(fpm,"# %s\n", title );
		fprintf(fpm,"# M O D E   S H A P E   D A T A   F O R   M O D E");
		fprintf(fpm,"   %d\t(global coordinates)\n", m );
		fprintf(fpm,"# deflection exaggeration: %.1f\n\n", exagg_modal );
		mpfX = 0.0;	for (i=1; i<=DoF; i++)    mpfX += V[i][m]*msX[i];
		mpfY = 0.0;	for (i=1; i<=DoF; i++)    mpfY += V[i][m]*msY[i];
		mpfZ = 0.0;	for (i=1; i<=DoF; i++)    mpfZ += V[i][m]*msZ[i];
		fprintf(fpm,"# MODE %5d:   f= %lf Hz, T= %lf sec\n", m,f[m],1./f[m]);
		fprintf(fpm,"#\t\tX- modal participation factor = %12.4e \n", mpfX);
		fprintf(fpm,"#\t\tY- modal participation factor = %12.4e \n", mpfY);
		fprintf(fpm,"#\t\tZ- modal participation factor = %12.4e \n", mpfZ);

		for(i=1; i<=DoF; i++)	v[i] = V[i][m];

		fprintf(fpm,"#      X-dsp       Y-dsp       Z-dsp\n\n");

		for(n=1; n<=nE; n++)
			bent_beam ( fpm, J1[n], J2[n], xyz, L[n], p[n], v, exagg_modal );

		for ( j=1; j<=nJ; j++ ) {
			if (xyz[j].x != 0.0) X=1;	/* check for three-dimensional frame */
			if (xyz[j].y != 0.0) Y=1;
			if (xyz[j].z != 0.0) Z=1;
		}

		if ( X && Y && Z ) D3 = ' ';

		fclose(fpm);

		if ((fpm = fopen (plotpath, "a")) == NULL) {
			printf (" error: cannot append plot file: %s\n",plotpath);
			exit(1);
		}
		fprintf(fpm,"pause -1\n");
		fprintf(fpm,"set nolabel\n");
		fprintf(fpm,"set title '%s     mode %d     %lf Hz'\n",IN_file,m,f[m]);
		fprintf(fpm,"plot '%s' u 2:3 t 'undeformed mesh' w l ", meshpath );
		if (!anlyz) fprintf(fpm," lw 2 lt 1 \n");
		else fprintf(fpm," lw 1 lt 5 , '%s' u 1:2 t 'mode-shape %d' w l lw 2 lt 3\n",
								modefl, m );
		fprintf(fpm,"%c pause -1\n", D3 );
		fprintf(fpm,"%c set nokey\n", D3 );
		fprintf(fpm,"%c splot '%s' u 2:3:4 t 'undeformed mesh' w l ",
								D3, meshpath);
		if (!anlyz) fprintf(fpm," lw 2 lt 1 \n");
		else fprintf(fpm," lw 1 lt 5 , '%s' u 1:2:3 t 'mode-shape %d' w l lw 2 lt 3\n",
								modefl, m );

		fclose(fpm);

	}

	free_dvector(msX,1,DoF);
	free_dvector(msY,1,DoF);
	free_dvector(msZ,1,DoF);
	free_dvector(v,1,DoF);
}

/*------------------------------------------------------------------------------
ANIMATE -  create mesh data of animated mode-shape meshes, use gnuplot	16dec98
	 useful gnuplot options: set noxtics noytics noztics noborder view nokey
	 mpeg movie example:   % convert mesh_file-03-f-*.ps mode-03.mpeg
	 ... requires ImageMagick and mpeg2vidcodec packages
------------------------------------------------------------------------------*/
void animate(
	char IN_file[], char meshpath[], char modepath[], char plotpath[],
	char *title,
	int anim[],
	int nJ, int nE, int DoF, int nM,
	vec3 *xyz, double *L, float *p,
	int *J1, int *J2, double *f, double **V,
	double exagg_modal,
	float pan
){
	FILE	*fpm;

	double x_min = 0, x_max = 0,
		y_min = 0, y_max = 0,
		z_min = 0, z_max = 0,
		rot_x_init = 70,	/* inital x-rotation in 3D animation */
		rot_x_final = 60,	/* final  x-rotation in 3D animation */
		rot_z_init = 100,	/* inital z-rotation in 3D animation */
		rot_z_final = 120,	/* final  z-rotation in 3D animation */
		zoom_init = 1.1,	/* inital zoom scale in 3D animation */
		zoom_final = 1.1,	/* final  zoom scale in 3D animation */
		frames = 25,	/* number of frames in animation	*/
		ex=10,		/* an exageration factor, for animation */
		*v;

	int	fr, i,j, m,n, X=0, Y=0, Z=0, c, CYCLES=3,
		frame_number = 0,
		total_frames;	/* total number of frames in animation */

	char	D3 = '#',
		Movie = '#',	/* use '#' for no-movie  -OR-  ' ' for movie */
		modefl[128], framefl[128];

	for (j=1; j<=nJ; j++) {		/* check for three-dimensional frame */
		if (xyz[j].x != 0.0) X=1;
		if (xyz[j].y != 0.0) Y=1;
		if (xyz[j].z != 0.0) Z=1;
		if (xyz[j].x < x_min ) x_min = xyz[j].x;
		if (xyz[j].y < y_min ) y_min = xyz[j].y;
		if (xyz[j].z < z_min ) z_min = xyz[j].z;
		if ( x_max < xyz[j].x ) x_max = xyz[j].x;
		if ( y_max < xyz[j].y ) y_max = xyz[j].y;
		if ( z_max < xyz[j].z ) z_max = xyz[j].z;
	}
	if ( X && Y && Z ) D3 = ' ';


	if ((fpm = fopen (plotpath, "a")) == NULL) {
		printf (" error: cannot append plot file: %s\n",plotpath);
		exit(1);
	}
	i = 0;
	while ( (m = anim[i]) != 0 && i < 20) {
	 if ( i==0 ) {

	   fprintf(fpm,"\n# --- M O D E   S H A P E   A N I M A T I O N ---\n");
	   fprintf(fpm,"set noborder\n");
	   fprintf(fpm,"set nozeroaxis\n");
	   fprintf(fpm,"set autoscale\n");
	   fprintf(fpm,"set noxtics; set noytics; set noztics; \n");
	   fprintf(fpm,"set nokey\n");
	   if ( x_min != x_max )
		fprintf(fpm,"set xrange [ %lf : %lf ] \n",
	 		x_min-0.2*(x_max-x_min), x_max+0.2*(x_max-x_min) );
	   else fprintf(fpm,"set xrange [ %lf : %lf ] \n",
			x_min-exagg_modal, x_max+exagg_modal );
	   if (y_min != y_max)
		fprintf(fpm,"set yrange [ %lf : %lf ] \n",
	 		y_min-0.2*(y_max-y_min), y_max+0.2*(y_max-y_min) );
	   else fprintf(fpm,"set yrange [ %lf : %lf ] \n",
			y_min-exagg_modal, y_max+exagg_modal );
	   if (z_min != z_max)
	   	fprintf(fpm,"set zrange [ %lf : %lf ] \n",
			z_min-0.2*(z_max-z_min), z_max+0.2*(z_max-z_min) );
	   else fprintf(fpm,"set zrange [ %lf : %lf ] \n",
			z_min-exagg_modal, z_max+exagg_modal );

	   fprintf(fpm,"%c set parametric\n", D3 );
	   fprintf(fpm,"%c set view 60, 70, 1 \n", D3 );
	   fprintf(fpm,"%c set xlabel \n", D3 );
	   fprintf(fpm,"%c set ylabel \n", D3 );
	   fprintf(fpm,"%c set zlabel \n", D3 );
	   fprintf(fpm,"%c set nolabel \n", D3 );
	 }

	 fprintf(fpm,"pause -1 \n");
	 fprintf(fpm,"set title '%s     mode %d      %lf Hz'\n",IN_file,m,f[m]);

	 frame_number = 0;
	 total_frames = 2*CYCLES*frames;
	 for ( c=1; c <= CYCLES; c++ ) {
	  for ( fr=0; fr<=frames; fr++ ) {

	    sprintf(modefl,"%s-%02d.%03d", modepath, m, fr  );
	    sprintf(framefl,"%s-%02d-f-%03d.ps", modepath, m, fr  );

	    if ( D3 == '#' ) {
		fprintf(fpm,"plot '%s' u 2:3 w l lw 1 lt 5, ", meshpath );
	 	fprintf(fpm," '%s' u 1:2 w l lw 2 lt 3 ;", modefl );
	    } else {
	      if ( pan != 0.0 )
		fprintf(fpm,"%c set view %5.1f, %5.1f, %4.2f \n", D3,
		rot_x_init + pan*(rot_x_final-rot_x_init)*frame_number/total_frames,
		rot_z_init + pan*(rot_z_final-rot_z_init)*frame_number/total_frames,
		zoom_init + pan*(zoom_final-zoom_init)*frame_number/total_frames );
	      fprintf(fpm,"%c splot '%s' u 2:3:4 w l lw 1 lt 5, ",D3,meshpath);
	      fprintf(fpm," '%s' u 1:2:3 w l lw 2 lt 3;", modefl );
	    }
	    if ( fr==0 && c==1 )	fprintf(fpm,"  pause 1.5 \n");
	    else			fprintf(fpm,"  pause 0.05 \n");
	    fprintf(fpm,"%c  load 'saveplot';\n",Movie);
	    fprintf(fpm,"%c  !mv my-plot.ps %s\n", Movie, framefl );
	  }
	  for ( fr = frames-1; fr > 0; fr-- ) {

	    sprintf(modefl,"%s-%02d.%03d", modepath, m, fr  );
	    sprintf(framefl,"%s-%02d-f-%03d.ps", modepath, m, fr  );

	    if ( D3 == '#' ) {
	 	fprintf(fpm,"plot '%s' u 2:3 w l lw 1 lt 5, ", meshpath );
		fprintf(fpm," '%s' u 1:2 w l lw 2 lt 3;", modefl );
	    } else {
	      if ( pan != 0.0 )
		fprintf(fpm,"%c set view %5.1f, %5.1f, %4.2f \n", D3,
		rot_x_init + pan*(rot_x_final-rot_x_init)*frame_number/total_frames,
		rot_z_init + pan*(rot_z_final-rot_z_init)*frame_number/total_frames,
		zoom_init + pan*(zoom_final-zoom_init)*frame_number/total_frames );
	      fprintf(fpm,"%c splot '%s' u 2:3:4 w l lw 1 lt 5, ",D3,meshpath);
	      fprintf(fpm," '%s' u 1:2:3 w l lw 2 lt 3;", modefl );
	    }
	    fprintf(fpm,"  pause 0.05 \n");
	    fprintf(fpm,"%c  load 'saveplot';\n",Movie);
	    fprintf(fpm,"%c  !mv my-plot.ps %s\n", Movie, framefl );
	  }
	 }
	 fr = 0;

	 sprintf(modefl,"%s-%02d.%03d", modepath, m, fr  );

	 if ( D3 == '#' ) {
	 	fprintf(fpm,"plot '%s' u 2:3 w l lw 2 lt 5, ", meshpath );
		fprintf(fpm," '%s' u 1:2 w l lw 3 lt 3 \n", modefl );
	 } else {
		fprintf(fpm,"%c splot '%s' u 2:3:4 w l lw 2 lt 5, ",D3,meshpath);
		fprintf(fpm," '%s' u 1:2:3 w l lw 3 lt 3 \n", modefl );
	 }

	 i++;
	}
	fclose(fpm);

	v = dvector(1,DoF);

	i = 0;
	while ( (m = anim[i]) != 0 ) {
	  for ( fr=0; fr<=frames; fr++ ) {

	    sprintf(modefl,"%s-%02d.%03d", modepath, m, fr  );

	    if ((fpm = fopen (modefl, "w")) == NULL) {
		printf (" error: cannot open modal mesh file: %s\n", modefl);
		exit(1);
	    }

	    ex = exagg_modal*cos( PI*fr/frames );

	    fprintf(fpm,"# FRAME3DD ANALYSIS RESULTS  http://frame3dd.sf.net/\n");
	    fprintf(fpm,"# %s\n", title );
	    fprintf(fpm,"# A N I M A T E D   M O D E   S H A P E   D A T A \n");
	    fprintf(fpm,"# deflection exaggeration: %.1f\n", ex );
	    fprintf(fpm,"# MODE %5d: f= %lf Hz  T= %lf sec\n\n",m,f[m],1./f[m]);

	    for (j=1; j<=DoF; j++)	v[j] = V[j][m];		/* mode "m" */

	    fprintf(fpm,"#      X-dsp       Y-dsp       Z-dsp\n\n");

	    for (n=1; n<=nE; n++)

		bent_beam ( fpm, J1[n], J2[n], xyz, L[n], p[n], v, ex );

	    fclose(fpm);
	  }
	  i++;
	}
	free_dvector(v,1,DoF);
	return;
}


/*------------------------------------------------------------------------------
BENT_BEAM  -  computes cubic deflection functions from end deflections
and end rotations.  Saves deflected shapes to a file.  These bent shapes
are exact for mode-shapes, and for frames loaded at their joints.
15 May 2009
------------------------------------------------------------------------------*/
void bent_beam(
	FILE *fp, int j1, int j2, vec3 *xyz,
	double L, float p, double *D, double exagg
){
	double	t1, t2, t3, t4, t5, t6, t7, t8, t9, 	/* coord xfmn	*/
		u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12,
		*a, *b, **A,
		s, v, w, dx, dy, dz;
	int	i1, i2, pd;

	A = dmatrix(1,4,1,4);
	a = dvector(1,4);
	b = dvector(1,4);

	coord_trans ( xyz, L, j1, j2,
			&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p );

	i1 = 6*(j1-1);	i2 = 6*(j2-1);

		/* compute end deflections in local coordinates */

	u1  = exagg*(t1*D[i1+1] + t2*D[i1+2] + t3*D[i1+3]);
	u2  = exagg*(t4*D[i1+1] + t5*D[i1+2] + t6*D[i1+3]);
	u3  = exagg*(t7*D[i1+1] + t8*D[i1+2] + t9*D[i1+3]);

	u4  = exagg*(t1*D[i1+4] + t2*D[i1+5] + t3*D[i1+6]);
	u5  = exagg*(t4*D[i1+4] + t5*D[i1+5] + t6*D[i1+6]);
	u6  = exagg*(t7*D[i1+4] + t8*D[i1+5] + t9*D[i1+6]);

	u7  = exagg*(t1*D[i2+1] + t2*D[i2+2] + t3*D[i2+3]);
	u8  = exagg*(t4*D[i2+1] + t5*D[i2+2] + t6*D[i2+3]);
	u9  = exagg*(t7*D[i2+1] + t8*D[i2+2] + t9*D[i2+3]);

	u10 = exagg*(t1*D[i2+4] + t2*D[i2+5] + t3*D[i2+6]);
	u11 = exagg*(t4*D[i2+4] + t5*D[i2+5] + t6*D[i2+6]);
	u12 = exagg*(t7*D[i2+4] + t8*D[i2+5] + t9*D[i2+6]);

		/* curve-fitting problem for a cubic polynomial */

	a[1] = u2;		b[1] = u3;
	a[2] = u8;   		b[2] = u9;
	a[3] = u6;		b[3] = -u5;
	a[4] = u12;		b[4] = -u11;
/*
	a[3] = exagg*tan(u6);	b[3] = exagg*tan(-u5);
	a[4] = exagg*tan(u12);	b[4] = exagg*tan(-u11);
*/

	u7 += L;
	A[1][1] = 1.0;   A[1][2] = u1;   A[1][3] = u1*u1;   A[1][4] = u1*u1*u1;
	A[2][1] = 1.0;   A[2][2] = u7;   A[2][3] = u7*u7;   A[2][4] = u7*u7*u7;
	A[3][1] = 0.0;   A[3][2] = 1.;   A[3][3] = 2.*u1;   A[3][4] = 3.*u1*u1;
	A[4][1] = 0.0;   A[4][2] = 1.;   A[4][3] = 2.*u7;   A[4][4] = 3.*u7*u7;
	u7 -= L;

	lu_dcmp ( A, 4, a, 1, 1, &pd );		/* solve for cubic coef's */

	if (!pd) {
	 fprintf(stderr," j1 = %d  j2 = %d  L = %e  u7 = %e \n", j1,j2,L,u7);
	 exit(1);
	}

	lu_dcmp ( A, 4, b, 0, 1, &pd );		/* solve for cubic coef's */

	for ( s = u1; s <= 1.01*(L+u7); s += (L+u7-u1) / 10.0 ) {

			/* deformed shape in local coordinates */
		v = a[1] + a[2]*s + a[3]*s*s + a[4]*s*s*s;
		w = b[1] + b[2]*s + b[3]*s*s + b[4]*s*s*s;

			/* deformed shape in global coordinates */
		dx = t1*s + t4*v + t7*w;
		dy = t2*s + t5*v + t8*w;
		dz = t3*s + t6*v + t9*w;

		fprintf (fp," %12.4e %12.4e %12.4e\n",
			xyz[j1].x + dx, xyz[j1].y + dy, xyz[j1].z + dz
		);
	}
	fprintf(fp,"\n\n");

	free_dmatrix(A,1,4,1,4);
	free_dvector(a,1,4);
	free_dvector(b,1,4);

	return;
}


/*------------------------------------------------------------------------------
SFERR  -  Display error message upon an erronous *scanf operation
------------------------------------------------------------------------------*/
void sferr ( char s[] ) {
	fprintf(stderr," >> Input Data file error while reading %s\n",s);
	return;
}


/*------------------------------------------------------------------------------
MY_ITOA  -  Convert an integer n to charcters in s, from K&R, 1978,   p. 59-60
... specialized for portability between GNU GCC and DJGPP GCC
------------------------------------------------------------------------------*/
void my_itoa ( int n, char s[], int k ) {
	int	c, i, j, sign;

	if ((sign = n) < 0) 		/* record sign */
		n = -n;			/* make n positive */
	i = 0;
	do {				/* generate digits in reverse order */
		s[i++] = n % 10 + '0';	/* get next digit */
	} while ((n /= 10) > 0);	/* delete it */
	for (;i<k;)	s[i++] = '0';	/* add leading '0' */
	if (sign < 0)
		s[i++] = '-';
	s[i] = '\0';
					/* reverse order of string s */
	j = 0;
	while ( s[j] != '\0' )	j++;	/* j is length of s - 1 */
	--j;

	for (i = 0; i < j; i++, j--) {
		c = s[i];
		s[i] = s[j];
		s[j] = c;
	}
	return;
}


/*------------------------------------------------------------------------------
GET_FILE_EXT  -  get the file extension,
		return 1 if the extension is ".csv"
		return 2 if the extension is ".fmm"
		return 0 otherwise
------------------------------------------------------------------------------*/
int get_file_ext( char *filename, char *ext )
{
	int	i=0, full_len=0, len=0;

	while ( filename[len++] != '\0' ) /* the length of file filename */ ;
	full_len = len;
	while ( filename[len--] != '.' && len > 0 ) /* the last '.' in filename */ ;
	if ( len == 0 )	len = full_len;
	++len;

	for ( i=0; len < full_len; i++,len++ ) ext[i] = tolower(filename[len]);

	/* test */
//	printf(" filename '%s' has length %d and extension = '%s' \n",
//							filename, len, ext);
//	printf(" Is .CSV? ... = %d \n", !strcmp(ext,".csv") );

	if ( !strcmp(ext,".csv") ) return (1);
	if ( !strcmp(ext,".fmm") ) return (2);
	return(0);
}


/*------------------------------------------------------------------------------
DOTS  -  print a set of dots (periods)
------------------------------------------------------------------------------*/
void dots ( FILE *fp, int n ) {
	int i;
	for (i=1; i<=n; i++)	fprintf(fp,".");
}


/*------------------------------------------------------------------------------
EVALUATE -  displays a randomly-generated goodbye message.  
------------------------------------------------------------------------------*/
void evaluate ( float error ) {
	int r;

	r = rand() % 10;

	if ( error < 1e-5 ) {

	    switch ( r ) {
		case 0: printf(" ...  awesome!  \n"); break; 
		case 1: printf(" ...  excellent!  \n"); break; 
		case 2: printf(" ...  woo-hoo!  \n"); break; 
		case 3: printf(" ...  yipee!  \n"); break; 
		case 4: printf(" ...  hoo-ray!  \n"); break; 
		case 5: printf(" ...  outstanding!  \n"); break; 
		case 6: printf(" ...  job well done!  \n"); break; 
		case 7: printf(" ...  priceless!  \n"); break; 
		case 8: printf(" ...  very nice!  \n"); break; 
		case 9: printf(" ...  very good!  \n"); break; 
	    }
	    return;
	}
	
	if ( error < 1e-3 ) {

	    switch ( r ) {
		case 0: printf(" ...  certainly acceptable!  \n"); break; 
		case 1: printf(" ...  bling!  \n"); break; 
		case 2: printf(" ...  that will certainly do!  \n"); break; 
		case 3: printf(" ...  not at all shabby!  \n"); break; 
		case 4: printf(" ...  quite reasonable!  \n"); break; 
		case 5: printf(" ...  and there's nothing wrong with that!  \n"); break; 
		case 6: printf(" ...  up to snuff!  \n"); break; 
		case 7: printf(" ...  bully!  \n"); break; 
		case 8: printf(" ...  nice!  \n"); break; 
		case 9: printf(" ...  good!  \n"); break; 
	    }
	    return;
	}

	if ( error < 1e-2 ) {

	    switch ( r ) {
		case 0: printf(" ...  adequate.  \n"); break; 
		case 1: printf(" ...  passable.  \n"); break; 
		case 2: printf(" ...  all right.  \n"); break; 
		case 3: printf(" ...  ok.  \n"); break; 
		case 4: printf(" ...  not bad.  \n"); break; 
		case 5: printf(" ...  fine.  \n"); break; 
		case 6: printf(" ...  fair. \n"); break; 
		case 7: printf(" ...  respectable.  \n"); break; 
		case 8: printf(" ...  tolerable.  \n"); break; 
		case 9: printf(" ...  just ok.  \n"); break; 
	    }
	    return;
	}

	if ( error > 1e-2 ) {

	    switch ( r ) {
		case 0: printf(" ...  abominable.  \n"); break; 
		case 1: printf(" ...  cruddy.  \n"); break; 
		case 2: printf(" ...  atrocious.  \n"); break; 
		case 3: printf(" ...  not ok.  \n"); break; 
		case 4: printf(" ...  garbage.  \n"); break; 
		case 5: printf(" ...  crappy.  \n"); break; 
		case 6: printf(" ...  oh noooo. \n"); break; 
		case 7: printf(" ...  abominable.  \n"); break; 
		case 8: printf(" ...  bummer.  \n"); break; 
		case 9: printf(" ...  awful.  \n"); break; 
	    }
	    return;
	}

}

