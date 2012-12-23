#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
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

#define ATOM_FROM  1
#define ATOM_TO   43

char *folder_name;




char *Fn[MAXCOOR];
int N_files = 0;

char buffer_names[MAXCOOR*100];
char *endbuffer = buffer_names;





int main(int argc, char *argv[])
{
	double rmsd;
	FILE *fi;
	char *namei;

    if (argc < 2) {
    	printf("Not enough parameters\n");
    	printf("Usage: %s [pdbfile1] [pdbfile2]\n", "rmsd");
        exit(EXIT_FAILURE);	
    } 
    
 	namei = argv[1];

	folder_name = "/home/miguel/Dropbox/arg32/md/ATP_frames_ALA/";





{
	char *name_source = argv[1];
	FILE *fs;
	char name_line[1024];

	if (NULL != (fs = fopen(name_source, "r"))) {

		printf ("open: %s\n",name_source);

		while (fgets(name_line, 1024, fs)) {
			char *j;
			//FIXME trim  blanks and spaces
			for (j = name_line; *j && isspace(*j); j++) {;}
			if (*j == '\0') continue;

			name_line[strlen(name_line)-1] = '\0';


			Fn[N_files++] = endbuffer;
			strcpy(endbuffer,name_line);
			endbuffer += strlen(name_line) + 1;

		}

		endbuffer = '\0';

		{int i;
		for (i = 0; i < N_files; i++) {
			printf ("%s\n",Fn[i]);
		}
		}
		fclose(fs);	
	} else {
		printf("problems with %s\n", name_source);
	}

	printf ("Files read=%d\n", N_coor);
}

//



{
	int i;
	char *name_source = argv[1];
	FILE *fs;
	char *name_line;
	char name_i[1024];
	struct model model_input;

	if (NULL != (fs = fopen(name_source, "r"))) {

		printf ("open: %s\n",name_source);


		
		for (i = 0; i < N_files; i++) {
			name_line = Fn[i];

			name_i[0] = '\0';
			strcpy(name_i, folder_name);
			strcat(name_i, name_line);

			if (NULL != (fi = fopen(name_i, "r"))) {
				struct coordinates *pma = &Coor[N_coor++];

				printf ("read: %s\n",name_i);

				modelload(fi, &model_input);
				mod2_ALLcoord (&model_input, pma, ATOM_FROM, ATOM_TO);
				fclose(fi);
			} else {
				printf("problems with %s\n", namei);
			}
		}
		fclose(fs);	
	} else {
		printf("problems with %s\n", name_source);
	}

	printf ("Files read=%d\n", N_coor);
}

{
	int i;
	struct transrot tr;
	for (i = 1; i < N_coor; i++) {
		rmsd = fit (&Coor[0], &Coor[i], &tr);
		printf("RMSD [%d,%d]: %.4lf\n",0,i,rmsd);
	}

	exit(0);
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


