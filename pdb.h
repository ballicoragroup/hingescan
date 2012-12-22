#if !defined(H_PDB)
#define H_PDB
/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
#include <stdio.h>
#include "fit.h"

enum {MAXLINE = 255, MAXRECORDS = 10000};
enum {ATOM = 1, HETATM = 2 };

struct atom {
	int type;
	int atmnum;
	char atmlabel[6];
	char alternative;
	char reslabel[6];
 	char chain;
	int  resnumber;
	double x;
	double y;
	double z;
	double o;
	double b;
};

struct model {
    int n;
    struct atom a[MAXRECORDS];
};    

void modelload(FILE *fi, struct model *pmodel);
void fprintatom (FILE *f, struct atom *p);
void model_rotate (const struct rotdata *r, struct model *m);
void fprintpdb (FILE *f, struct model *model);
void model_transrot (const struct transrot *tr, struct model *m);

/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
#endif

