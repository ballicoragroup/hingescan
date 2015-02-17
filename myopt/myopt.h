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


#if !defined(H_MYOPT)
#define H_MYOPT
/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

enum char_options {END_OF_OPTIONS = -1};

extern int 			opt_index;      
extern const char *	opt_arg;    

extern int 			options(int argc, char *argv[], const char *legal);

/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
#endif

