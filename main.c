#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include "bool.h"
#include "pdb.h"
#include "fit.h"
#include <math.h>


static void coord_extract (const struct coordinates *c, struct coordinates *o, int from, int to);

//-----------------------------------------------------------------------------

void mod2_CAcoord(const struct model *m, struct coordinates *c);
void mod2_ALLcoord(const struct model *m, struct coordinates *c, int from, int to);

static void	findhinges (struct model *model_a, struct model *model_b, int window );

//=============================================================================

static void
collectmypdb	(const char *name_i, struct model *pmodel_input);


int main(int argc, char *argv[])
{
	struct model MIA, MIB;
	int window = 81;
	const char *name_A = "data/ns-c.pdb";
	const char *name_B = "data/ns-o.pdb";


    if (argc < 2) {
    	printf("Not enough parameters\n");
    	printf("Usage: %s [pdblist] ...\n", "rmsd");
        exit(EXIT_FAILURE);	
    } 

	collectmypdb	(name_A, &MIA);
	collectmypdb	(name_B, &MIB);

	findhinges (&MIA, &MIB, window);

    return EXIT_SUCCESS;
}

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

//=============================================================================

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
    	if (j >= MAXATOMS) {
    		fprintf (stderr, "Number of atoms exceeded limit\n");
    		break;
		}
 	}
 	c->n = j;
}

void mod2_ALLcoord(const struct model *m, struct coordinates *c, int from, int to)
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

//==================================================================================

static void
findhinges (struct model *model_a, struct model *model_b, int window )
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
				
		printf("%d, %lf\n",av+16,rmsd);
	}

}

