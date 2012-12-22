#include <string.h>
#include <ctype.h>
#include "pdb.h"
#include "fit.h"
#include "strtools.h"
#include "bool.h"

void atomadd(const char *s, struct model *m);
bool_t ishead (char *s, char *t);
bool_t isEND (char *s);


void modelload(FILE *fi, struct model *pmodel)
{

	char s[MAXLINE+1];
	
	pmodel->n = 0; /* reset model */
	
	while (NULL != mygets(s, MAXLINE, fi)) {
		if (ishead(s, "ATOM  ")) {
			atomadd(s, pmodel);
  		}
		if (ishead(s, "HETATM")) {
			atomadd(s, pmodel);
  		}  		
  		if (isEND (s)) {
  			break;
    	}
 	}

}


void atomadd(const char *s, struct model *m)
{
	enum {MAXVARS = 12};
	int cols[MAXVARS+1] = {1, 7, 12, 17, 18, 21, 23, 31, 39, 47, 55, 61, 67};

	char inp[MAXLINE+1];
	int v, i, k, sz;
	struct atom at;
 	struct atom at_none = {0, 0, "", '_', "", '_', -1, 0.0, 0.0, 0.0, 0.0, 0.0}; 
	
	at = at_none;
	
	sz = strlen(s);
	
	for (v = 0; v < MAXVARS; v++) {
		int from = cols[v    ] - 1;
		int top  = cols[v + 1] - 1;
		k = 0;
		for (i = from; i < sz && i < top; i++) {
 			inp[k++] = s[i];
  		} 
  		inp[k] = '\0';
 		
  		switch (v) {
  		    case  0: if (ishead(inp, "ATOM  ")) { 
        					at.type = ATOM; 
   					} else if (ishead(inp, "HETATM")) {
                			at.type = HETATM;
     			    };
                   	break;	 
  		    case  1: sscanf(inp, "%d", &at.atmnum); 			break;
            case  2: strncpy(at.atmlabel, inp, 6); 	
            			trim (at.atmlabel);
                		break;
            case  3: at.alternative = inp[0]; 					break;
            case  4: strncpy(at.reslabel, inp, 6);
            			trim(at.reslabel);
           				break;
            case  5: at.chain = inp[1]; 						break;
            case  6: sscanf(inp, "%d" , &at.resnumber); 		break;               
            case  7: sscanf(inp, "%lf", &at.x); 				break;
            case  8: sscanf(inp, "%lf", &at.y); 				break;
            case  9: sscanf(inp, "%lf", &at.z); 				break;
            case 10: sscanf(inp, "%lf", &at.o); 				break;
            case 11: sscanf(inp, "%lf", &at.b); 				break;  		
            default: 											break;
        }    
  	}
  	if (m->n < MAXRECORDS) {
   		m->a[m->n] = at;
   		m->n++;
    } 
}


void fprintpdb (FILE *f, struct model *model)
{
	int i;
	for (i = 0; i < model->n; i ++) {
    	fprintatom(f, &model->a[i]);
	}  
}


void fprintatom (FILE *f, struct atom *p)
{
	if (p->type == ATOM)
		fprintf(f, "%s","ATOM  ");
	else if (p->type == HETATM)
		fprintf(f, "%s","HETATM");
 	else
 		fprintf(f, "%s","XXXXXX");
 		
   	fprintf(f, "%5d  %-3s%c%-3s %c%4d",
    		p->atmnum, p->atmlabel, p->alternative, 
    		p->reslabel, p->chain, p->resnumber);
  
    fprintf(f, "    "
    "%8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n", 
    p->x, p->y, p->z, p->o, p->b); 	
}



char *s_atom (struct atom *p, char *buffer_out, int max)
{
	char *atm;
	char buffer[256];
	char *s = buffer;
 	char *e = buffer;
	char *q = buffer_out;
	int accum;
	*s = '\0';
	
	if (p->type == ATOM) {
		atm = "ATOM  ";
	} else if (p->type == HETATM) {
		atm = "HETATM";		
	} else { 
		atm = "XXXXXX";
	}
	e += sprintf(e, "%s", atm);
   	e += sprintf(e, "%6d:%3s%c:%3s:%c:%5d",
    		p->atmnum, 
      		p->atmlabel, 
        	p->alternative, 
         	p->reslabel, 
          	p->chain, 
           	p->resnumber);
  						
	accum = 0;
	while (accum < (max -1) && *s) {
 		*q++ = *s++; accum++;
 	}
    *q = '\0';
    
    return buffer_out;
}


bool_t
ishead (char *s, char *t)
{
	return 	s[0]==t[0] && s[1]==t[1] && s[2]==t[2] && 
            s[3]==t[3] && s[4]==t[4] && s[5]==t[5];
}


bool_t
isEND (char *s)
{
	return 	s[0]=='E' && s[1]=='N' && s[2]=='D' && 
            (s[3]=='\0' || s[3]==EOF || isspace(s[3]));
}

void model_rotate (const struct rotdata *r, struct model *m)
{
	int d, i, k, j;
	double v[NAXIS], w[NAXIS];
	double M[NAXIS] [NAXIS];

  	for (i = 0; i < NAXIS; i++) 
  	  	for (j = 0; j < NAXIS; j++)     	
      	  	M[i][j] = r->mat[i][j]; 	

    for (d = 0; d < m->n; d++) {
   
        v[0] = m->a[d].x; 
        v[1] = m->a[d].y;
        v[2] = m->a[d].z;
        
	    /*atmbck (r->ctr, v);*/
	
 	    /* matrix multiply vector v, result in w */
  	    for (i = 0; i < NAXIS; i++) {
  	  		double tmp = 0.0;
       		for (k = 0; k < NAXIS; k++ ) {
            	tmp += M[i][k] * v[k];
            } 		
            w[i] = tmp;
      	} 	
    
    	/*atmtns (r->ctr, w);*/
    	atmtns (r->tns, w);
        m->a[d].x = w[0]; 
        m->a[d].y = w[1];
        m->a[d].z = w[2];         
     }	
}


void model_transrot (const struct transrot *tr, struct model *m)
{
	int d;
 	double v[NAXIS];

    for (d = 0; d < m->n; d++) {
   
        v[0] = m->a[d].x; 
        v[1] = m->a[d].y;
        v[2] = m->a[d].z;
        
        do_transrot(tr, v);
        
        m->a[d].x = v[0]; 
        m->a[d].y = v[1];
        m->a[d].z = v[2];         
     }	
}
