/*
	Hingescan is program for detecting hinges in protein structures
    Copyright 2015 Miguel A. Ballicora

    This file is part of Hingescan.

    Hingescan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Ordo is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Hingescan.  If not, see <http://www.gnu.org/licenses/>.
*/

#if !defined(H_PDB)
#define H_PDB
/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
#include <stdio.h>
#include "fit.h"

enum {MAXLINE = 255, MAXRECORDS = 10000};
enum {ATOM = 1, HETATM = 2 };

struct atom {
	int type;
	int atmnum;
	char atmlabel[6];
	char alternative;
	char reslabel[6];
 	char chain;
	int  resnumber;
	double x;
	double y;
	double z;
	double o;
	double b;
};

struct model {
    int n;
    struct atom a[MAXRECORDS];
};    

void modelload(FILE *fi, struct model *pmodel);
void fprintatom (FILE *f, struct atom *p);
void model_rotate (const struct rotdata *r, struct model *m);
void fprintpdb (FILE *f, struct model *model);
void model_transrot (const struct transrot *tr, struct model *m);

int model_get_first_residue_number (struct model *m);

/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
#endif

