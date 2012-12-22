#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bool.h"
#include "pdb.h"
#include "fit.h"
#include <math.h>

void mod2_CAcoord(const struct model *m, struct coordinates *c);

struct model model [2];    

double dist(struct atom *a, struct atom *b) ;

int main(int argc, char *argv[])
{

	double p[NPARAM];
	double res;
	int i;
	int Target;
	FILE *fi, *fb, *fo;
	char *namei, *nameb, *nameo;
	double D_cutoff = 4.0;

	struct coordinates mdA, mdB;
 

    if (argc < 4) {
    	printf("Not enough parameters\n");
    	printf("Usage: %s [pdbfile1] [pdbfile2]\n", "rmsd");
        exit(EXIT_FAILURE);	
    } 
    
 	namei = argv[1];
 	nameb = argv[2];
 	nameo = argv[3];  

      
	if (NULL != (fi = fopen(namei, "r"))) {
 		if (NULL != (fb = fopen(nameb, "r"))) {
 			if (NULL != (fo = fopen(nameo, "w"))) { 		
	 			int i;
	 			struct atom *pa;
		 		struct rotdata r;
		 		struct transrot tr;
           			
     			modelload(fi, &model[0]);
     			modelload(fb, &model[1]);
     			
        		mod2_CAcoord (&model[0], &mdA);
     			mod2_CAcoord (&model[1], &mdB);
     			
#if 0
     			for (i = 0; i < NPARAM; i++) 
        			p[i] = 0;
     			masscenter(&mdB, &p[6]);
     			res = sqrt(sqdev(&mdA, &mdB, p));
				p[0] = 3.1459;
				getrotation(p, &r);
     			model_rotate(&r, &model[1]);     			
           		printf("%d, %d, %.2lf\n",mdA.n, mdB.n, res); 
     			fprintpdb(fo, &model[1]);
     			system("PAUSE");                        		
#endif

#if 1
     			for (i = 0; i < NPARAM; i++) 
        			p[i] = 0;
				printf("starting RMSD: %.4lf\n",sqrt(sqdev(&mdA, &mdB, p)));
     			fit   (&mdA, &mdB, &tr);
 			
     			model_transrot(&tr, &model[1]);
     			
     			/*model_rotate(&r, &model[1]);*/   
     			fprintpdb(fo, &model[1]);             			
    			//system("PAUSE");
#endif
     			
                          			
#if 0
     			for (i = 0; i < 6; i ++) {
        			printf("%.3lf, %.3lf, %.3lf\n",
           				mdA.atm[i][0],mdA.atm[i][1],mdA.atm[i][2]);
     			}
     			printf("\n");
     			for (i = 0; i < 6; i ++) {
        			printf("%.3lf, %.3lf, %.3lf\n",
           				mdB.atm[i][0],mdB.atm[i][1],mdB.atm[i][2]);
     			}
     			printf("\n");
#endif
				fclose(fo);
			} else {
				printf("problems with %s\n", nameo);			
   			}
   			
   			fclose(fb);
   		} else {
   			printf("problems with %s\n", nameb);
     	}
     	
   	   	fclose(fi);
   	} else {
		printf("problems with %s\n", namei);
    }

    return 0;
}




double dist(struct atom *a, struct atom *b) 
{
	double dx, dy, dz, t;
	dx = a->x - b->x;
	dy = a->y - b->y;
	dz = a->z - b->z;
	
	t = dx*dx + dy*dy + dz*dz;
	
	return sqrt(t);
}

void mod2_CAcoord(const struct model *m, struct coordinates *c)
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
    	if (j >= MAXATOMS)
    		break;
	
 	}
 	c->n = j;
}


