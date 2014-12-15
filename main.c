#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>

#include "bool.h"
#include "pdb.h"
#include "fit.h"
#include "proginfo.h"

static void coord_extract (const struct coordinates *c, struct coordinates *o, int from, int to);
static void gnuplot_out (const char *name_gp, int xfrom, int xto, int wfrom, int wto, double topz);

//-----------------------------------------------------------------------------

//static void mod2_CAcoord (const struct model *m, struct coordinates *c);
//static void mod2_ALLcoord (const struct model *m, struct coordinates *c, int from, int to);

static void	findhinges (bool_t isquiet, struct model *model_a, struct model *model_b, 
					int botwindow, int topwindow, bool_t corrected, bool_t autoz, double topz, FILE *outf );

/*
|
|	GENERAL OPTIONS
|
\*--------------------------------------------------------------*/

#include "myopt.h"

const char *license_str = "\n"
"   Copyright (c) 2014 Miguel A. Ballicora\n"
"   Rmsdscan is program for detecting hinges in protein structures\n"
"\n"
"   Rmsdscan is free software: you can redistribute it and/or modify\n"
"   it under the terms of the GNU General Public License as published by\n"
"   the Free Software Foundation, either version 3 of the License, or\n"
"   (at your option) any later version.\n"
"\n"
"   Rmsdscan is distributed in the hope that it will be useful,\n"
"   but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
"   GNU General Public License for more details.\n"
"\n"
"   You should have received a copy of the GNU General Public License\n"
"   along with Rmsdscan.  If not, see <http://www.gnu.org/licenses/>.\n"
;

static void parameter_error(void);
static void example (void);
static void usage (void);

/* VARIABLES */

	static const char *copyright_str = 
		"Copyright (c) 2014 Miguel A. Ballicora\n"
		"There is NO WARRANTY of any kind\n"
		;

	static const char *intro_str =
		"Program to detect hinges by comparing two conformers\n"
		;

	const char *example_options = 
		"-a file1.pdb -b file2.pdb -w 21 -o out.csv";

	static const char *example_str =
		"  - Processes file1.pdb and file2.pdb with a window size of 21 residues\n"
		"  - Outputs a list of hinge scores for each residue in out.csv\n"
		;

	static const char *help_str =
		" -h          print this help\n"
		" -H          print just the switches\n"
		" -v          print version number and exit\n"
		" -L          display the license information\n"
		" -q          quiet mode (no screen progress updates)\n"
		" -w <n>      window size (always odd). Min size if -W is provided\n"
		" -W <n>      multi scan, from -w to -W window sizes (always odd numbers)\n"
		" -Z <n>      top hinge score for plotting, default is the maximum observed\n"
		" -c          corrected hinge scores\n"
		" -a <file>   first input file in pdb format\n"
		" -b <file>   second input file in pdb format\n"
		" -o <file>   output file (text format), goes to the screen if not present\n"
		"\n"
		;

	const char *usage_options = 
		"[-OPTION]";
		;
	/*	 ....5....|....5....|....5....|....5....|....5....|....5....|....5....|....5....|*/
		

const char *OPTION_LIST = "hHvLqa:b:o:w:W:Z:c";

//=============================================================================

static void parameter_error(void) {	printf ("Error in parameters, switch unknown\n"); return;}

static void
example (void)
{
	printf ("\n"
		"quick example: %s %s\n"
		"%s"
		, proginfo_name()
		, example_options
		, example_str);
	return;
}

static void
usage (void)
{
	printf ("\n"
		"usage: %s %s\n"
		"%s"
		, proginfo_name()
		, usage_options
		, help_str);
}

//===================================================================================================

static void
collectmypdb	(const char *name_i, struct model *pmodel_input);

int main(int argc, char *argv[])
{
	FILE *outf;
	struct model MIA, MIB;
	int window = 31;
	int topwindow = 31;
	bool_t multiwin = FALSE;
	bool_t corrected = FALSE;
	bool_t autoz = TRUE;
	double topz = 0;	
	
	const char *inputa = "";
	const char *inputb = "";
	const char *textstr= "";
	const char *inputf;

	int op = 0;

	bool_t QUIET_MODE;
	int version_mode, help_mode, switch_mode, license_mode, input_mode;


	/* defaults */
	version_mode = FALSE;
	license_mode = FALSE;
	help_mode    = FALSE;
	switch_mode  = FALSE;
	input_mode   = FALSE;
	QUIET_MODE   = FALSE;
	inputa       = NULL;
	inputb       = NULL;
	textstr      = NULL;
	corrected    = FALSE;
	
	while (END_OF_OPTIONS != (op = options (argc, argv, OPTION_LIST))) {
		switch (op) {
			case 'v':	version_mode = TRUE; 	break;
			case 'L':	version_mode = TRUE; 	
						license_mode = TRUE;
						break;
			case 'h':	help_mode = TRUE;		break;
			case 'H':	switch_mode = TRUE;		break;
			case 'a': 	input_mode = TRUE;
					 	inputa = opt_arg;
						break;
			case 'b': 	input_mode = TRUE;
					 	inputb = opt_arg;
						break;
			case 'o': 	textstr = opt_arg;
						break;
			case 'w': 	if (1 != sscanf(opt_arg,"%d", &window) || window < 0) {
							fprintf(stderr, "wrong simulation parameter\n");
							exit(EXIT_FAILURE);
						} else if ((((unsigned)window)&1) == 0) {
							fprintf(stderr, "-w <window> needs to be an odd number\n");
							exit(EXIT_FAILURE);								
						}
						break;
			case 'W': 	if (1 != sscanf(opt_arg,"%d", &topwindow) || topwindow < 0) {
							fprintf(stderr, "wrong simulation parameter\n");
							exit(EXIT_FAILURE);
						} else if ((((unsigned)window)&1) == 0) {
							fprintf(stderr, "-W <n> needs to be an odd number\n");
							exit(EXIT_FAILURE);								
						}
						multiwin = TRUE;
						break;
			case 'Z': 	if (1 != sscanf(opt_arg,"%lf", &topz) || topz < 0) {
							fprintf(stderr, "wrong maximum hinge score provided\n");
							exit(EXIT_FAILURE);
						}
						autoz = FALSE;
						break;						
			case 'q':	QUIET_MODE = TRUE;	break;
			case 'c':	corrected = TRUE;	break;
			case '?': 	parameter_error();
						exit(EXIT_FAILURE);
						break;
			default:	fprintf (stderr, "ERROR: %d\n", op);
						exit(EXIT_FAILURE);
						break;
		}		
	}

	/*----------------------------------*\
	|	Return version
	\*----------------------------------*/
	if (version_mode) {
		printf ("%s %s\n",proginfo_name(),proginfo_version());
		if (license_mode)
 			printf ("%s\n", license_str);
		return EXIT_SUCCESS;
	}
	if (argc < 2) {
		fprintf (stderr, "%s %s\n",proginfo_name(),proginfo_version());
		fprintf (stderr, "%s", copyright_str);
		fprintf (stderr, "for help type:\n%s -h\n\n", proginfo_name());
		exit (EXIT_FAILURE);
	}
	if (help_mode) {
		printf ("\n%s", intro_str);
		example();
		usage();
		printf ("%s\n", copyright_str);
		exit (EXIT_SUCCESS);
	}
	if (switch_mode && !help_mode) {
		usage();
		exit (EXIT_SUCCESS);
	}
	if ((argc - opt_index) > 1) {
		/* too many parameters */
		fprintf (stderr, "ERROR: Extra parameters present\n");
		fprintf (stderr, "Make sure to surround parameters with \"quotes\" if they contain spaces\n\n");
		exit(EXIT_FAILURE);
	}
	if (input_mode && argc != opt_index) {
		fprintf (stderr, "Extra parameters present\n");
		fprintf (stderr, "Make sure to surround parameters with \"quotes\" if they contain spaces\n\n");
		exit(EXIT_FAILURE);
	}
	if (!input_mode && argc == opt_index) {
		fprintf (stderr, "Need file name to proceed\n\n");
		exit(EXIT_FAILURE);
	}
	
	/* get rest, should be only one at this point */
	while (opt_index < argc ) {
		inputf = argv[opt_index++];
		printf ("Extra parameter in command line: %s\n",inputf);
	}

	multiwin = multiwin && topwindow > window;
	if (!multiwin) topwindow = window;

	/*========*/

	collectmypdb	(inputa, &MIA);
	collectmypdb	(inputb, &MIB);

	if (!QUIET_MODE)
		printf ("output to = %s\n",textstr==NULL? "stdout": textstr);

	if (NULL != textstr && NULL != (outf = fopen(textstr, "w"))) {
		findhinges (QUIET_MODE, &MIA, &MIB, window, topwindow, corrected, autoz, topz, outf);
		fclose(outf);
	} else {
		findhinges (QUIET_MODE, &MIA, &MIB, window, topwindow, corrected, autoz, topz, stdout);
	}
	
    return EXIT_SUCCESS;
}

//=============================================================================

static void
collectmypdb	(const char *name_i, struct model *pmodel_input)
{
	FILE *fi;
	if (NULL != (fi = fopen(name_i, "r"))) {
		modelload(fi, pmodel_input);
		fclose(fi);
	} else {
		printf("problems with %s\n", name_i);
		exit(EXIT_FAILURE);
	}
	return;
}

//-----------------------------------------------------------------------------

static
void coord_extract (const struct coordinates *c, struct coordinates *o, int from, int to)
{
	int i, j, end;
	int n = c->n;
	
	end = to+1 > n? n: to+1;
		
	for (i = from, j = 0; i < end; i++) {
			o->atm[j][0] = c->atm[i][0];
			o->atm[j][1] = c->atm[i][1];
			o->atm[j][2] = c->atm[i][2];
  			j++;
 	}
 	o->n = j;
}

static void mod2_CAcoord(const struct model *m, struct coordinates *c)
{
	int i, j;
	int n = m->n;
	
	for (i = 0, j = 0; i < n ; i++) {
		
		if (m->a[i].atmlabel[0] == 'C' && m->a[i].atmlabel[1] == 'A') {
  			c->atm[j][0] = m->a[i].x;
  			c->atm[j][1] = m->a[i].y;
  			c->atm[j][2] = m->a[i].z;
  			j++;
    	} 
    	if (j >= MAXATOMS) {
    		fprintf (stderr, "Number of atoms exceeded limit\n");
    		break;
		}
 	}
 	c->n = j;
}

#if 0
static void mod2_ALLcoord(const struct model *m, struct coordinates *c, int from, int to)
{
	int i, j;
	int n = m->n;

	for (i = 0, j = 0; i < n ; i++) {
		
		if (   m->a[i].atmnum >= from 
			&& m->a[i].atmnum <= to
			&& m->a[i].atmlabel[0] != 'H'
			) {
  			c->atm[j][0] = m->a[i].x;
  			c->atm[j][1] = m->a[i].y;
  			c->atm[j][2] = m->a[i].z;
  			j++;
    	} 
    	if (j >= MAXATOMS) {
    		fprintf (stderr, "Number of atoms exceeded limit\n");
    		break;
		}
 	}
 	c->n = j;
}
#endif

//==================================================================================

static void
findhinges (bool_t isquiet, struct model *model_a, struct model *model_b, int botwindow, int topwindow, bool_t corrected, bool_t autoz, double topz, FILE *outf)
{
	int reference_window;
	double rmsd;
	double zmax;
	int n_slices, fr, to, av, j;
	int shift = 0;
	int window;
	bool_t multi = topwindow > botwindow;
	bool_t auto_zscale = autoz;

	struct coordinates SCA, SCB, SCA_all, SCB_all;   
	struct transrot tr;

	mod2_CAcoord (model_a, &SCA_all);
	mod2_CAcoord (model_b, &SCB_all);

	// validation
	if (SCA_all.n != SCB_all.n) {
		fprintf (stderr, "Number of total residues do not match:\n"); 
		fprintf (stderr, "model 1: %d residues\n", SCA_all.n); 
		fprintf (stderr, "model 2: %d residues\n", SCB_all.n); 			
		exit(0);
	} 

	reference_window = botwindow;
	zmax = 0;
	
	for (window = botwindow; window < topwindow+1; window += 2) {

		//--------------------------------------------------
		n_slices = (SCA_all.n - window + 1);

		shift = model_get_first_residue_number (model_a);

		if (!isquiet) {
			if (multi) {
				printf ("\nwindow = %d\n",window); 
			}
			printf ("alpha carbons = %d\n", SCA_all.n); 
			printf ("atoms in model = %d\n", model_a->n); 
			printf ("number of slices = %d\n", n_slices); 
			printf ("shift to first residue numbers = %d\n",shift); 
			if (corrected) 
				printf ("window reference for correction = %d\n",reference_window); 
		}

		// head
		if (multi) {
			for (j = 0; j < (window-1)/2; j++) {
				// output a blank score
				fprintf (outf, "%4d ", 0);
			}
		}
		
		for (j = 0; j < n_slices; j++) {

			fr = j;
			to = fr + window - 1;
			av = (fr+to)/2;

			coord_extract (&SCA_all, &SCA, fr, to);
			coord_extract (&SCB_all, &SCB, fr, to);

			//fprintf (stderr,"from=%d, to=%d\n",fr,to);
			//fprintf (stderr,"n transferred=%d, %d\n",SCA.n, SCB.n);

			if (SCA.n == 0 || SCB.n == 0) {
				fprintf (stderr, "Warning: File could be empty\n"); exit(0);
			} 

			if (SCA.n != SCB.n) {
				fprintf (stderr, "Number of residues do not match for slice: %d\n",j); 
				fprintf (stderr, "model 1: %d residues\n", SCA.n); 
				fprintf (stderr, "model 2: %d residues\n", SCB.n); 			
				exit(0);
			} 
		
			rmsd = fit (&SCA, &SCB, &tr);

			// debug code
			if (0) {
				if (av == 514) {
					FILE *fo;
					model_transrot(&tr, model_b);
					if (NULL != (fo = fopen("btmpout_B.pdb","w"))) {
						fprintpdb(fo, model_b); 
						fclose (fo);
					}
					if (NULL != (fo = fopen("btmpout_A.pdb","w"))) {
						fprintpdb(fo, model_a); 
						fclose (fo);
					}
				}
			}

			if (!multi)	{
				fprintf(outf,"%d, %lf\n",av+shift,rmsd);
			} else{
				double r = (double)window/(double)reference_window;
				double corr_x = 1.0 * rmsd / (1.0 + 0.5 * log(r));				
				double out_x = corrected? corr_x: rmsd;

				zmax = zmax < out_x? out_x: zmax;
			
				fprintf (outf, "%4d ", (int)(1000.0 * out_x));
			}
		}

		// Tail
		if (multi) {
			for (j = 0; j < (window-1)/2; j++) {
				fprintf (outf, "%4d ", 0);
			}
			fprintf (outf, "\n");
		}

		printf ("\n");
		
	} // end loop windows
	
	if (!isquiet && multi) {
		printf ("Maximum Hinge Score: %lf\n", zmax);
	}				
	
	gnuplot_out ("gp.plt", shift, SCA_all.n+shift-1, botwindow, topwindow, auto_zscale? zmax: topz);

}



static const char *block1[9] =
{"reset"
,"set terminal svg"
,"set output 'out.svg'"
,"unset key"
,""
,"set style line 11 lc rgb '#808080' lt 1"
,"set border 3 front ls 11"
,"set tics nomirror out scale 0.75"
,NULL
};


static const char *block3[7] = 
{"set xlabel 'Residue'"
,"set ylabel 'Window'"
,""
,"set colorbox noborder"
,"set cblabel \"Hinge Score\""
,""
,NULL
};

static const  char *block4[3] = 
{"show cbtics"
,""
,NULL
};

static const char *palette[37] =
{"set palette defined(\\"
,"0       0.2314  0.2980  0.7529,\\"
,"0.03125 0.2667  0.3529  0.8000,\\"
,"0.0625  0.3020  0.4078  0.8431,\\"
,"0.09375 0.3412  0.4588  0.8824,\\"
,"0.125   0.3843  0.5098  0.9176,\\"
,"0.15625 0.4235  0.5569  0.9451,\\"
,"0.1875  0.4667  0.6039  0.9686,\\"
,"0.21875 0.5098  0.6471  0.9843,\\"
,"0.25    0.5529  0.6902  0.9961,\\"
,"0.28125 0.5961  0.7255  1.0000,\\"
,"0.3125  0.6392  0.7608  1.0000,\\"
,"0.34375 0.6824  0.7882  0.9922,\\"
,"0.375   0.7216  0.8157  0.9765,\\"
,"0.40625 0.7608  0.8353  0.9569,\\"
,"0.4375  0.8000  0.8510  0.9333,\\"
,"0.46875 0.8353  0.8588  0.9020,\\"
,"0.5     0.8667  0.8667  0.8667,\\"
,"0.53125 0.8980  0.8471  0.8196,\\"
,"0.5625  0.9255  0.8275  0.7725,\\"
,"0.59375 0.9451  0.8000  0.7255,\\"
,"0.625   0.9608  0.7686  0.6784,\\"
,"0.65625 0.9686  0.7333  0.6275,\\"
,"0.6875  0.9686  0.6941  0.5804,\\"
,"0.71875 0.9686  0.6510  0.5294,\\"
,"0.75    0.9569  0.6039  0.4824,\\"
,"0.78125 0.9451  0.5529  0.4353,\\"
,"0.8125  0.9255  0.4980  0.3882,\\"
,"0.84375 0.8980  0.4392  0.3451,\\"
,"0.875   0.8706  0.3765  0.3020,\\"
,"0.90625 0.8353  0.3137  0.2588,\\"
,"0.9375  0.7961  0.2431  0.2196,\\"
,"0.96875 0.7529  0.1569  0.1843,\\"
,"1       0.7059  0.0157  0.1490\\"
,")"
,""
,NULL
};


static void
gnuplot_out (const char *name_gp, int xfrom, int xto, int wfrom, int wto, double topz)
{
	FILE *fi;
	if (NULL != (fi = fopen(name_gp, "w"))) {
		// output here
		const char **s;
		int i;

		s = block1;			
		for (i = 0; s[i] != NULL; i++) {
			fprintf (fi, "%s\n", s[i]);
		}

		fprintf (fi,"set xrange [%d:%d]\n",xfrom-1, xto+1);
		fprintf (fi,"set yrange [%d:%d]\n",wfrom-1, wto+1);
		
		s = block3;			
		for (i = 0; s[i] != NULL; i++) {
			fprintf (fi, "%s\n", s[i]);
		}
		
		fprintf (fi,"set cbrange [0:%lf]\n",topz);		
				
		s = block4;			
		for (i = 0; s[i] != NULL; i++) {
			fprintf (fi, "%s\n", s[i]);
		}		
		
		s = palette;			
		for (i = 0; s[i] != NULL; i++) {
			fprintf (fi, "%s\n", s[i]);
		}

		fprintf (fi,"plot 'out.txt' u (($1+(%d))/1.0):(%d+2*($2)):($3/1000.0) matrix with image\n",xfrom,wfrom);
		
		fclose(fi);
	} else {
		fprintf(stderr, "problems to open for output: %s\n", name_gp);
		exit(EXIT_FAILURE);
	}
	return;
}
