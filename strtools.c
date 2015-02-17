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

#include <string.h>
#include "strtools.h"

char *
trim (char *s)
{
	char *b, *e, *p, *q;
	
	e = s + strlen(s);
	b = s;

   	while (e > b && *(e-1)==' ') e--;
    while (e > b && *b    ==' ') b++;
    
    for (p = s, q = b; q < e; p++, q++) {
    	*p = *q;
    }
    *p = '\0';
    return s;
}


char *
mygets (char *s, int n, FILE *iop)
{
	register int c = EOF;
	register char *cs;

	if (!(n > 0))
		return NULL;

	cs = s;
	while (--n > 0 && (c = getc(iop)) != EOF) {
		*cs++ = (char) c;
		if (c == '\n')
			break;
	}
	*cs = '\0';
	return (c == EOF && cs == s) ? NULL : s;
}

