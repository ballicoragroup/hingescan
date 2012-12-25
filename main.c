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
struct coordinates Coor[MAXCOOR];
int N_coor = 0;

void mod2_CAcoord(const struct model *m, struct coordinates *c);
void mod2_ALLcoord(const struct model *m, struct coordinates *c, int from, int to);

struct model model [2];    

int ATOM_FROM;
int ATOM_TO;

char Folder_line[1024];
char *Folder_name;


char *Fn[MAXCOOR];
int N_files = 0;

char buffer_names[MAXCOOR*100];
char *endbuffer = buffer_names;


static char *trimCRLF (char *x)
{
	size_t y = strlen(x)-1;
	x [y] = '\0';		
	return x;
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

int main(int argc, char *argv[])
{
	double rmsd;
	FILE *fi;
	char *namei;

    if (argc < 4) {
    	printf("Not enough parameters\n");
    	printf("Usage: %s [pdblist] [atom from] [atom to]\n", "rmsd");
        exit(EXIT_FAILURE);	
    } 
    
 	namei = argv[1];
	Folder_name = "";

	if (	1!=sscanf (argv[2],"%d",&ATOM_FROM)
		||	1!=sscanf (argv[3],"%d",&ATOM_TO)
		) {
		fprintf(stderr, "Error in input parameters\n");
		exit(EXIT_FAILURE);
	}

#if 0
{
	char *name_source = argv[1];
	FILE *fs;
	char name_line[1024];


	if (NULL != (fs = fopen(name_source, "r"))) {

		printf ("open: %s\n",name_source);

		if (fgets(folder_line, 1024, fs)) {
			trimCRLF(folder_line);
			folder_name	= folder_line;

			printf ("folder_name=%s\n",folder_name);

			while (fgets(name_line, 1024, fs)) {
				char *j;
				trimCRLF(name_line);
				for (j = name_line; *j && isspace(*j); j++) {;}
				if (*j == '\0') {
					continue;
				}
				Fn[N_files++] = endbuffer;
				strcpy(endbuffer,name_line);
				endbuffer += strlen(name_line) + 1;

			}

			endbuffer = '\0';

		}
		// {int i;	for (i = 0; i < N_files; i++) {printf ("%s\n",Fn[i]);}}

		fclose(fs);	
	} else {
		printf("problems with %s\n", name_source);
	}

	printf ("File names read=%d\n", N_files);
}
#endif

{
	char *name_source;

	name_source = argv[1];

	endbuffer = process_setfile	( name_source
								, Folder_line
								, Fn
								, MAXCOOR
								, &N_files
								, endbuffer
								);

	Folder_name = Folder_line;

	printf ("File names read=%d\n", N_files);


}
//

{
	int i;
	char *name_line;
	char name_i[1024];
	struct model model_input;

		
	for (i = 0; i < N_files; i++) {
		name_line = Fn[i];

		name_i[0] = '\0';
		strcpy(name_i, Folder_name);
		strcat(name_i, name_line);

		if (NULL != (fi = fopen(name_i, "r"))) {
			struct coordinates *pma = &Coor[N_coor];

			modelload(fi, &model_input);
			mod2_ALLcoord (&model_input, pma, ATOM_FROM, ATOM_TO);

			if (pma->n == 0) {
				fprintf (stderr, "Warning: File %s could be empty\n", name_i);
			} else {
				N_coor++;
			}

			fclose(fi);
		} else {
			printf("problems with %s\n", namei);
		}
	}

	printf ("Files read=%d\n", N_coor);
}

{
	int i;
	struct transrot tr;
	for (i = 0; i < N_coor; i++) {
		assert (Coor[i].n > 0);
		rmsd = fit (&Coor[0], &Coor[i], &tr);
		//printf("RMSD [%d,%d]: %.4lf\n",0,i,rmsd);
	}
}



{
	FILE *ofile;
	int i,j,n;
	struct transrot tr;
	n = N_coor; 

	if (NULL != (ofile = fopen("infile", "w"))) {

		fprintf(ofile,"%d\n",n);
		for (i = 0; i < n; i++) {

			fprintf (stderr,"reference: %d\n",i);
			fprintf(ofile, "%-10d",i);
	
			for (j = 0; j < n; j++) {
				assert (Coor[i].n > 0 && Coor[j].n);
				rmsd = fit (&Coor[i], &Coor[j], &tr);
				fprintf(ofile," %.4lf",rmsd);
			}
			fprintf(ofile,"\n");
		}
		fclose(ofile);
	}
}

    return EXIT_SUCCESS;
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

void mod2_ALLcoord(const struct model *m, struct coordinates *c, int from, int to)
{
	int i, j;
	int n = m->n;
	
	for (i = 0, j = 0; i < n ; i++) {
		
		if (m->a[i].atmnum >= from && m->a[i].atmnum <= to) {
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
				fnam[N_fnam++] = endbuffer;
				strcpy(endbuffer,name_line);
				endbuffer += strlen(name_line) + 1;

			}

			endbuffer = '\0';

		}
		// {int i;	for (i = 0; i < N_fnam; i++) {printf ("%s\n",fnam[i]);}}

		fclose(fs);	
	} else {
		printf("problems with %s\n", name_source);
	}

	printf ("File names read=%d\n", N_fnam);

	*fnam_n = N_fnam;

	return endbuffer;
}


