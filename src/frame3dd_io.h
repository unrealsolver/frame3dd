/*
 This file is part of FRAME3DD:
 Static and dynamic structural analysis of 2D and 3D frames and trusses with
 elastic and geometric stiffness.
 ---------------------------------------------------------------------------
 http://frame3dd.sourceforge.net/
 ---------------------------------------------------------------------------
 Copyright (C) 1992-2014  Henri P. Gavin

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
	Frame input/output functions (FRAME3DD native format)
*/

#include <time.h>
#include <stdio.h>
#include <unistd.h>	/* getopt for parsing command-line options	*/

#include "common.h"
#include "compat_types.h"
#include "microstran/vec3.h"
#include "types.h"


/**
 PARSE_OPTIONS -  parse command line options			     04mar09
 command line options over-ride values in the input data file
*/
RuntimeArgs parse_options (
	int argc, char *argv[],
	char IN_file[], char OUT_file[]
);

/**
 DISPLAY_HELP -  display help information to stderr
04mar09
*/
void display_help();

/**
 DISPLAY_USAGE -  display usage information to stderr
04mar09
*/
void display_usage();

/**
 DISPLAY_VERSION -  display version, website, and some help info to stderr
04mar09
*/
void display_version();

/**
 DISPLAY_VERSION_BACKGROUND -  display version and website to stderr
22sep09
*/
void display_version_about();

/**
FRAME3DD_GETLINE
	get line into a character string. from K&R.
	@NOTE this is different from the GNU 'getline'.
*/
int frame3dd_getline( FILE *fp, char *s, int lim );


/**
	Re-write input file without comments or commas
	and output to file with path tpath

	@param fp input file pointer
	@param tpath output clean file path (temp file)
*/
void parse_input(FILE *fp, const char *tpath);


/**
 * Read number of nodes from the input file
 */
void read_node_number(
	FILE *fp,	/**< input data file pointer			*/
	int *nN,
	const int verbose
);

/**
	Read node coordinate data
*/
void read_node_data (
	FILE *fp,	/**< input data file pointer			*/
	Frame *frame,	/**< Refacotred data				*/
	InputScope *scope /**< 'Compatibility' input data scope		*/
);

/**
	Read frame element property data
*/
void read_frame_element_data (
	FILE *fp,	/**< input data file pointer			*/
	Frame *frame,	/**< New frame data struct */
	InputScope *scope
);


/**
	Read data controlling certain aspects of the analysis
*/
void read_run_data (
	FILE *fp,	/**< input data file pointer			*/
	char OUT_file[], /**< output data file name			*/
	char meshpath[],/**< file name for mesh data output		*/
	char plotpath[],/**< file name for Gnuplot script		*/
	char infcpath[],/**< file name for internal force data		*/
	RunOptions *run_options, /**< RunOptions object          	*/
	RuntimeArgs args,  /**< Command line options			*/
	InputScope *scope
);


/**
	Read fixed node displacement boundary conditions
*/
void read_reaction_data(
	FILE *fp,	/**< input data file pointer			*/
	int verbose,	/**< 1: copious screen output; 0: none		*/
	Frame *frame,
	InputScope *scope
);

void read_element_number(
	FILE *fp,	/**< input data file pointer			*/
	int *nE,
	const int verbose
);

/**
	read load information data, form un-restrained load vector
*/
void read_and_assemble_loads(
	FILE *fp,	/**< input data file pointer			*/
	int verbose,		/**< 1: copious output to screen, 0: none */
	Frame *frame,
	LoadCases *load_cases,	/**< Load cases array  */
	InputScope *scope	/**< Input data struct */
);


/**
	read member densities and extra inertial mass data
*/
void read_mass_data(
	FILE *fp,	/**< input data file pointer			*/
	char *OUT_file,	/**< input output data file name 		*/
	char modepath[], /**< filename for mode shape data for plotting	*/
	RuntimeArgs args, /**< Command line options		*/
	InputScope *scope /**< Input data combined */
);


/*
 * READ_CONDENSATION_DATA - read matrix condensation information
 */
void read_condensation_data(
	FILE *fp,	/**< input data file pointer			*/
	int verbose,	/**< 1: copious output to screen, 0: none	*/
	RuntimeArgs args, /**< Command line options		*/
	InputScope *scope /**< Input data combined */
);


/*
 * WRITE_INPUT DATA - write input data to a file
 */
void write_input_data(
	FILE *fp,	/**< input data file pointer			*/
	char *title,
	int *nD, int nR,
	int *nF, int *nU, int *nW, int *nP, int *nT,
	float *p, float *d,
	double **Ft, double **Fm, float **Dp,
	int *R,
	float ***U, float ***W, float ***P, float ***T,
	int shear, int anlyz, int geom,
	Frame *frame,
	LoadCases *load_cases
);


/*
 * WRITE_STATIC_RESULTS
 * save node displacements, reactions, and element end forces in a text file
 * 9sep08, 2014-05-15
 */
void write_static_results(
	FILE *fp,
	int lc, int DoF,
	int *N1, int *N2,
	double *F, double *D, double *R, int *r, double **Q,
	double err, int ok, int axial_sign,
	Frame *frame, LoadCases *load_cases
);


/*
 * CSV_filename - return the file name for the .CSV file and
 * whether the file should be written or appended (wa)
 * 1 Nov 2015
*/
void CSV_filename( char CSV_file[], char wa[], char OUT_file[], int lc );


/*
 * WRITE_STATIC_CSV
 * save node displacements, reactions, and  element end forces in a .CSV file
 * 31dec08, 2015-05-15
 */
void write_static_csv(
	char *OUT_file, char *title,
	int lc, int DoF,
	int *N1, int *N2,
	double *F, double *D, double *R, int *r, double **Q,
	double err, int ok,
	Frame *frame, LoadCases *load_cases
);



/*
 * WRITE_STATIC_MFILE
 * save node displacements, reactions, and element end forces in an m-file
 * 9sep08, 2015-05-15
 */
void write_static_mfile(
	char *OUT_file, char *title,
	int lc, int DoF,
	int *N1, int *N2,
	double *F, double *D, double *R, int *r, double **Q,
	double err, int ok,
	Frame *frame, LoadCases *load_cases
);


/*
 * PEAK_INTERNAL_FORCES
 *	calculate frame element internal forces, Nx, Vy, Vz, Tx, My, Mz
 *	calculate frame element local displacements, Rx, Dx, Dy, Dz
 *	return the peak (maximum absolute) values
 *	18 jun 2013
 */
void peak_internal_forces (
        int lc,         // load case number
        vec3 *xyz,      // node locations
        double **Q, double *L, int *N1, int *N2,
        float *p,
        float *d, float gX, float gY, float gZ,
        int nU, float **U, int nW, float **W, int nP, float **P,
        double *D, int shear,
	float dx,	// x-axis increment along frame element

        // vectors of peak forces, moments, displacements and slopes
	// for each frame element, for load case "lc"
        double **pkNx, double **pkVy, double **pkVz,
        double **pkTx, double **pkMy, double **pkMz,
        double **pkDx, double **pkDy, double **pkDz,
        double **pkRx, double **pkSy, double **pkSz,
	Frame *frame, LoadCases *load_cases
);


/*
 * WRITE_INTERNAL_FORCES
 *	calculate frame element internal forces, Nx, Vy, Vz, Tx, My, Mz
 *	calculate frame element local displacements, Rx, Dx, Dy, Dz
 *	write internal forces and local displacements to an output data file
 *	4jan10
 */
void write_internal_forces(
	char *OUT_file, /**< output data filename                       */
	FILE *fp,	/**< pointer to output data file		*/
	char infcpath[],/**< interior force data file			*/
	int lc,		/**< load case number				*/
	char title[],	/**< title of the analysis			*/
	float dx,	/**< increment distance along local x axis      */
	vec3 *xyz,	/**< XYZ locations of each node                */
	double **Q,	/**< frame element end forces                   */
	double *L,	/**< length of each frame element               */
	int *N1, int *N2, /**< node connectivity                       */
	float *p,	/**< roll angle, radians                        */
	float *d,	/**< mass density                               */
	float gX, float gY, float gZ,	/**< gravitational acceleration */
	int nU,		/**< number of uniformly-distributed loads	*/
	float **U,	/**< uniformly distributed load data            */
	int nW,		/**< number of trapezoidally-distributed loads	*/
	float **W,	/**< trapezoidally distributed load data        */
	int nP,		/**< number of internal point loads		*/
	float **P,	/**< internal point load data                   */
	double *D,	/**< node displacements                        */
	int shear,	/**< shear deformation flag                     */
	double error,	/**< RMS equilibrium error			*/
	Frame *frame,
	LoadCases *load_cases
);


/*
 * WRITE_MODAL_RESULTS
 *	save modal frequencies and mode shapes			16aug01
 */
void write_modal_results(
	FILE *fp,
	int nI, int DoF,
	double **M, double *f, double **V,
	double total_mass, double struct_mass,
	int iter, int sumR, int nM,
	double shift, int lump, double tol, int ok, Frame *frame
);


/*
 * SFERR - display *scanf errors from output from *scanf function
 */
void sferr( char *s );


/*
 * GET_FILE_EXT
 *	return the file extension, including the period (.)
 *        return 1 if the extension is ".csv"
 *        return 2 if the extension is ".fmm"
 *        return 0 otherwise
 */
int get_file_ext( char *filename, char *ext );


/*
 * TEMP_FILE_LOCATION
 *	Return location for temporary 'clean' file
 *
 *	@param fname name of the file, excluding directory.
 *	@param fullpath (returned) full path to temporary file.
 *	@param len available length for string being returned.
 */
void temp_file_location(const char *fname, char fullpath[], const int len);


/*
 * OUTPUT_PATH
 *	Return path to an output file.
 *
 *	If the fname starts with a path separator ('\' on Windows, else '/') then
 *	it is assumed to be an absolute path, and will be returned unchanged.
 *
 *	If the fname doesn't start with a path separator, it is assumed to be a
 *	relative path. It will be appended to the contents of FRAME3DD_OUTDIR if
 *	that path is defined.
 *
 *	If FRAME3DD_OUTDIR is not defined, it will be appended to the contents of
 *	the parameter default_outdir. If the parameter default_outdir is given a
 *	NULL value, then it will be replaced by the appropriate temp file location
 *	(see temp_file_location and temp_dir in frame3dd_io.c).
 *
 *	@param fname name of the file, excluding directory.
 *	@param fullpath (returned) full path to temporary file.
 *	@param len available length for string being returned.
 *	@param default_outdir default output directory in case where no env var.
 */
void output_path(const char *fname, char fullpath[], const int len, const char *default_outdir);


/*
 * DOTS - print a set of dots (periods)
 */
void dots ( FILE *fp, int n );


/*
 *  EVALUATE -  displays a randomly-generated evaluation message.
 */
void evaluate (  float error, float rms_resid, float tol );

