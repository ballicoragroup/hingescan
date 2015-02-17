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

#if !defined(H_FIT)
#define H_FIT
/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

enum {NPARAM = 6, NAXIS = 3};
enum {MAXATOMS = 10000};

struct coordinates {
	int n;
	double atm [MAXATOMS] [NAXIS];
};

struct rotdata {
	double mat [NAXIS] [NAXIS];
	double tns [NAXIS];
};

struct transrot {
	double m [NAXIS+1] [NAXIS+1];
};

/*-----------------------------------------------*/

extern void masscenter(const struct coordinates *tp, double *c);

extern double sqdev(const struct coordinates *tp, 
					const struct coordinates *mb, 
     				const double *p);

extern double fit	( const struct coordinates *t
					, const struct coordinates *m
					, struct transrot *out);

void do_transrot (const struct transrot *inp, double *io);     				
     				


void atmtns (const double *delta, double *v);
/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
#endif

