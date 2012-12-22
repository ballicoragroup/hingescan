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

