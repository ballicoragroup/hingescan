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

//-----------------------------------------------------------------------------

//static void mod2_CAcoord (const struct model *m, struct coordinates *c);
//static void mod2_ALLcoord (const struct model *m, struct coordinates *c, int from, int to);
static void	findhinges (struct model *model_a, struct model *model_b, int window, FILE *outf );



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
		"Program to detect hinges\n"
		;

	const char *example_options = 
		"example__options________here";

	static const char *example_str =
		"  - Processes <> to...........\n"
		"  - and then ______________________________\n"
		;

	static const char *help_str =
		" -h          print this help\n"
		" -H          print just the switches\n"
		" -v          print version number and exit\n"
		" -L          display the license information\n"
		" -q          quiet mode (no screen progress updates)\n"
		" -a <file>   first input file in pdb format\n"
		" -b <file>   second input file in pdb format\n"
		" -o <file>   output file (text format), goes to the screen if not present\n"
		"\n"
		;

	const char *usage_options = 
		"[-OPTION]";
		;
	/*	 ....5....|....5....|....5....|....5....|....5....|....5....|....5....|....5....|*/
		

const char *OPTION_LIST = "hHvLqa:b:o:";

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
	int window = 81;
	const char *inputa = "data/ns-c.pdb";
	const char *inputb = "data/ns-o.pdb";
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
			case 'q':	QUIET_MODE = TRUE;	break;
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
	if (QUIET_MODE) {
		printf ("Quiet mode (-q) version not implemented\n");
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

	/*==== SET INPUT ====*/

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	collectmypdb	(inputa, &MIA);
	collectmypdb	(inputb, &MIB);

printf ("textstr=%s\n",textstr==NULL? "NULL": textstr);

	if (NULL != textstr && NULL != (outf = fopen(textstr, "w"))) {
		findhinges (&MIA, &MIB, window, outf);
		fclose(outf);
	} else {
		findhinges (&MIA, &MIB, window, stdout);
	}
    return EXIT_SUCCESS;
}

//=============================================================================
//-----------------------------------------------------------------------------

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
findhinges (struct model *model_a, struct model *model_b, int window, FILE *outf)
{
	double rmsd;
	int n_slices, fr, to, av, j;

	struct coordinates SCA, SCB, SCA_all, SCB_all;   
	struct transrot tr;

//

	mod2_CAcoord (model_a, &SCA_all);
	mod2_CAcoord (model_b, &SCB_all);

	n_slices = (SCA_all.n - window + 1);

	printf ("n_slices=%d, MIA.n=%d\n",n_slices,model_a->n); 
	printf ("window=%d\n",window); 

	for (j = 0; j < n_slices; j++) {

		fr = j;
		to = fr + window - 1;
		av = (fr+to)/2;

		coord_extract (&SCA_all, &SCA, fr, to);
		coord_extract (&SCB_all, &SCB, fr, to);

		//fprintf (stderr,"from=%d, to=%d\n",fr,to);
		//fprintf (stderr,"n transferred=%d, %d\n",SCA.n, SCB.n);

		if (SCA.n == 0 || SCB.n == 0 || SCA.n != SCB.n) {
			fprintf (stderr, "Warning: File could be empty\n"); exit(0);
		} 

		rmsd = fit (&SCA, &SCB, &tr);
				
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
				
		fprintf(outf,"%d, %lf\n",av+16,rmsd);
	}

}

