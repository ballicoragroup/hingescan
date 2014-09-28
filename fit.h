#if !defined(H_FIT)
#define H_FIT
/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

enum {NPARAM = 6, NAXIS = 3};
enum {MAXATOMS = 10000};

struct coordinates {
	int n;
	double atm [MAXATOMS] [NAXIS];
};

struct rotdata {
	double mat [NAXIS] [NAXIS];
	double tns [NAXIS];
};

struct transrot {
	double m [NAXIS+1] [NAXIS+1];
};

/*-----------------------------------------------*/

extern void masscenter(const struct coordinates *tp, double *c);

extern double sqdev(const struct coordinates *tp, 
					const struct coordinates *mb, 
     				const double *p);

extern double fit	( const struct coordinates *t
					, const struct coordinates *m
					, struct transrot *out);

void do_transrot (const struct transrot *inp, double *io);     				
     				


void atmtns (const double *delta, double *v);
/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
#endif

