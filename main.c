#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include "bool.h"
#include "pdb.h"
#include "fit.h"
#include <math.h>

#define MAXCOOR 10000

struct coordinates_set {
	struct coordinates coor[MAXCOOR];
	char *label[MAXCOOR]; // Label names
	int n;

};

struct coordinates_set Coor_set;

//--------------------------------------------------------------------

void mod2_CAcoord(const struct model *m, struct coordinates *c);
void mod2_ALLcoord(const struct model *m, struct coordinates *c, int from, int to);

//struct model model [2];    

int ATOM_FROM;
int ATOM_TO;

char Folder_line[1024];
char *Folder_name;


char *Fn[MAXCOOR]; // File names
int N_files = 0;

char Buffer_names[MAXCOOR*100];
char *Endbuffer = Buffer_names;

//=================================================

static char *trimCRLF (char *x)
{
	size_t y = strlen(x)-1;
	x [y] = '\0';		
	return x;
}

static char *
make_label_alloc(char *prefix, char *name)
{
	char *ori = Endbuffer;

	assert(Endbuffer && prefix && name);

	strcpy(Endbuffer,prefix);
	Endbuffer += strlen(prefix);
	strcpy(Endbuffer,name);
	Endbuffer += strlen(name);
	Endbuffer += 1;
	return ori;
}

//=================================================

static char *
process_setfile	( const char *name_source
				, /*@OUT@*/	char *infolder
				, char *fnam[]
				, int maxcoor
				, int *fnam_n
				, char *endbuffer
				);

static void
coord_collect	( char *folder
				, char *inputfiles[]
				, int n_input
				, int atom_from
				, int atom_to
				, char *prefix
				, struct coordinates_set *cs);

static int
triangle_key(int i, int j, int side_size)
{
	int waste, key;
	assert ( j >= i);
	waste = i * (i+1) / 2;
	key   = i * side_size + j;
	return key - waste;
};

int main(int argc, char *argv[])
{
	double rmsd, rmsd0;
	double *rbuf;
	int r;

    if (argc < 4) {
    	printf("Not enough parameters\n");
    	printf("Usage: %s [atom from] [atom to] [pdblist] [pdblist]\n", "rmsd");
        exit(EXIT_FAILURE);	
    } 
    
	Folder_name = "";
	Coor_set.n = 0;
	Endbuffer = Buffer_names;

	if (	1!=sscanf (argv[1],"%d",&ATOM_FROM)
		||	1!=sscanf (argv[2],"%d",&ATOM_TO)
		) {
		fprintf(stderr, "Error in input parameters\n");
		exit(EXIT_FAILURE);
	}


	Endbuffer = process_setfile	( argv[3]
								, Folder_line
								, Fn
								, MAXCOOR
								, &N_files
								, Endbuffer
								);

	Folder_name = Folder_line;

	printf ("File names read=%d\n", N_files);
	assert(Endbuffer);

	coord_collect	( Folder_name
					, Fn
					, N_files
					, ATOM_FROM
					, ATOM_TO
					, "A"
					, &Coor_set);
//

	Endbuffer = process_setfile	( argv[4]
								, Folder_line
								, Fn
								, MAXCOOR
								, &N_files
								, Endbuffer
								);

	Folder_name = Folder_line;

	printf ("File names read=%d\n", N_files);
	assert(Endbuffer);

	coord_collect	( Folder_name
					, Fn
					, N_files
					, ATOM_FROM
					, ATOM_TO
					, "B"
					, &Coor_set);

printf ("Total elements: %d\n", Coor_set.n);

	r = 0;
	rbuf = malloc (sizeof(double) * Coor_set.n * Coor_set.n);

	if (rbuf) {

		FILE *ofile;
		int i,j,n;
		struct transrot tr;
		char buffname [1024] = "infile";
		char *oname = buffname;

		n = Coor_set.n; 

		// calculation, store in triangular buffer
		for (i = 0; i < n; i++) {

			fprintf (stderr,"reference: %d\n",i);
	
			for (j = i; j < n; j++) {
				assert (Coor_set.coor[i].n > 0 && Coor_set.coor[j].n > 0);
				assert (j >= i);
				assert (r == triangle_key(i,j,n));
				rmsd = fit (&Coor_set.coor[i], &Coor_set.coor[j], &tr);
				assert(j != i || rmsd < 0.0001 || 0==fprintf(stderr,"Label=%s, self rmsd=%.4lf\n",Coor_set.label[i],rmsd));
				rbuf[r++] = rmsd;
			}
		}

		if (NULL != (ofile = fopen(oname, "w"))) {

			fprintf(ofile,"%d\n",n);
			for (i = 0; i < n; i++) {

				fprintf(ofile, "%-10s",Coor_set.label[i]);
		
				for (j = 0; j < n; j++) {
					assert (Coor_set.coor[i].n > 0 && Coor_set.coor[j].n);
					rmsd = rbuf [triangle_key(i<j?i:j,  i<j?j:i, n)];

//					rmsd0 = fit (&Coor_set.coor[i], &Coor_set.coor[j], &tr, i==j && i==628);
//					assert((rmsd0 - rmsd) < 0.0001 && (rmsd - rmsd0) < 0.0001);

					fprintf(ofile," %.4lf",rmsd);
				}
				fprintf(ofile,"\n");
			}
			fclose(ofile);
		} else {
				printf("problems opening %s\n", oname);
				exit(EXIT_FAILURE);
		}

		free(rbuf);

	} else {
		fprintf(stderr, "memory not enough\n");
		exit(EXIT_FAILURE);
	}

    return EXIT_SUCCESS;
}

//----------------------------------------------------------------------------

static void
coord_collect	( char *folder
				, char *inputfiles[]
				, int n_input
				, int atom_from
				, int atom_to
				, char *prefix
				, struct coordinates_set *cs)
{
	FILE *fi;
	int i;
	char *name_line;
	char name_i[1024];
	struct model model_input;
	int nc = cs->n;
		
	for (i = 0; i < n_input; i++) {
		name_line = inputfiles[i];

		name_i[0] = '\0';
		strcpy(name_i, folder);
		strcat(name_i, name_line);

		if (NULL != (fi = fopen(name_i, "r"))) {
			struct coordinates *pma = &(cs->coor[nc]);

			modelload(fi, &model_input);
			mod2_ALLcoord (&model_input, pma, atom_from, atom_to);

			if (pma->n == 0) {
				fprintf (stderr, "Warning: File %s could be empty\n", name_i);
			} else {
				// make a label for this file
				cs->label[nc] = make_label_alloc (prefix,name_line);
				nc++;
			}

			fclose(fi);
		} else {
			printf("problems with %s\n", name_i);
			exit(EXIT_FAILURE);
		}
	}

	printf ("Files read=%d\n", nc);

	cs->n = nc;

	return;
}

//========================================================================


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
    	if (j >= MAXATOMS)
    		break;
	
 	}
 	c->n = j;

}

//==================================================================================

static char *
process_setfile	( const char *name_source
				, /*@OUT@*/	char *infolder
				, char *fnam[]
				, int maxcoor
				, int *fnam_n
				, char *endbuffer
				)
{
	FILE *fs;
	char name_line[1024];
	int N_fnam = 0;

	if (NULL != (fs = fopen(name_source, "r"))) {

		printf ("open: %s\n",name_source);

		if (fgets(infolder, 1024, fs)) {
			trimCRLF(infolder);

			while (N_fnam < maxcoor && fgets(name_line, 1024, fs)) {
				char *j;
				trimCRLF(name_line);
				for (j = name_line; *j && isspace(*j); j++) {;}
				if (*j == '\0') {
					continue;
				}
				if (name_line[strlen(name_line)-1] == '~') {
					printf ("Ignored=%s\n",name_line);
					continue;
				}
				fnam[N_fnam++] = endbuffer;
				strcpy(endbuffer,name_line);
				endbuffer += strlen(name_line) + 1;

			}

			*endbuffer = '\0';

		}
		fclose(fs);	
	} else {
		printf("problems with %s\n", name_source);
		exit(EXIT_FAILURE);
	}

	*fnam_n = N_fnam;
	return endbuffer;
}


