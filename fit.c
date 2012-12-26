#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "fit.h"

#if 0
#define INSPECT
#endif

struct coordinates t_buffer, m_buffer;

/* 
	P[] 
	0: angle X
	1: angle Y
	2: angle Z
	3: trans X
	4: trans Y
	5: trans Z
	6: center X
	7: center Y
	8: center Z
*/


void coor_rotate (const struct coordinates *ori, const struct rotdata *r, struct coordinates *to);

void getrotation (const double *p, struct rotdata *r);
void atmcpy (const double *in, double *out);
void atmtns (const double *delta, double *v);
void atmbck (const double *delta, double *v);
double vector_sqdist(const double *v, const double *w);
double iterfit (const struct coordinates *t, 
				const struct coordinates *m, double *s_io, double *p_io);
double iterfitQ(const struct coordinates *t, 
				const struct coordinates *m, double *f_io, double *p_io);				
void prmcpy(const double *a, double *b);
void prmsca(const double *d, double s, double *b);
void showvector6(char *s, double *v);

void get_transrot (struct rotdata *r, struct transrot *out);

void compose_r_t (struct transrot *io, const double *v);
void compose_t_r (const double *v, struct transrot *io);

void transrotcpy (const struct transrot *inp, struct transrot *out);


void showrot (const struct rotdata *inp);
void showtransrot (const struct transrot *inp);

/*-----------------------------------------------------*/

void getrotation (const double *p, struct rotdata *r)
{
	double angleX, angleY, angleZ;
	double X[NAXIS][NAXIS];
 	double Y[NAXIS][NAXIS];
  	double Z[NAXIS][NAXIS];
   	double R[NAXIS][NAXIS];
	double cx, cy, cz, sx, sy, sz;
	int i, j, k;

	angleX = p[0];
	angleY = p[1];
	angleZ = p[2];
	
	cx = cos(angleX);
	sx = sin(angleX);
	cy = cos(angleY);
	sy = sin(angleY);
 	cz = cos(angleZ);
	sz = sin(angleZ);
 
  	Z[0][0] =  cz;
  	Z[0][1] = -sz;
  	Z[0][2] = 0.0;
  	Z[1][0] =  sz;
  	Z[1][1] =  cz;
  	Z[1][2] = 0.0; 
  	Z[2][0] = 0.0;
  	Z[2][1] = 0.0;
  	Z[2][2] = 1.0;
  	
  	Y[0][0] =  cy;
  	Y[0][1] = 0.0;
  	Y[0][2] =  sy;
  	Y[1][0] = 0.0;
  	Y[1][1] = 1.0;
  	Y[1][2] = 0.0; 
  	Y[2][0] = -sy;
  	Y[2][1] = 0.0;
  	Y[2][2] =  cy;     	
  	
  	X[0][0] = 1.0;
  	X[0][1] = 0.0;
  	X[0][2] = 0.0;
  	X[1][0] = 0.0;
  	X[1][1] =  cx;
  	X[1][2] = -sx; 
  	X[2][0] = 0.0;
  	X[2][1] =  sx;
  	X[2][2] =  cx;   	
  	
  	/* Multiply R = Z * X */
  	for (i = 0; i < NAXIS; i++) {
  	  	for (j = 0; j < NAXIS; j++) {    	
  	  		R[i][j] = 0.0;
      	}
    }
  	
  	for (i = 0; i < NAXIS; i++) {
  	  	for (j = 0; j < NAXIS; j++) {    	
  	  		double tmp = 0.0;
       		for (k = 0; k < NAXIS; k++ ) {
            	tmp += Z[i][k] * X[k][j];
            } 		
            R[i][j] = tmp;
      	}
    }

  	/* Multiply Z = Y * R */    
  	for (i = 0; i < NAXIS; i++) {
  	  	for (j = 0; j < NAXIS; j++) {    	
  	  		Z[i][j] = 0.0;
      	}
    }    
  	for (i = 0; i < NAXIS; i++) {
  	  	for (j = 0; j < NAXIS; j++) {    	
  	  		double tmp = 0.0;
       		for (k = 0; k < NAXIS; k++ ) {
            	tmp += Y[i][k] * R[k][j];
            } 		
            Z[i][j] = tmp;
      	}
    }        
  	
  	/* OUTPUT */
  	for (i = 0; i < NAXIS; i++) {
  	  	for (j = 0; j < NAXIS; j++) {    	
  	  		r->mat[i][j] = Z[i][j];
      	}
    } 
    r->tns[0] = p[3];
    r->tns[1] = p[4];
    r->tns[2] = p[5];
                               	
}



double sqdev(const struct coordinates *tp, const struct coordinates *mb, const double *p)
{
	int d, i, k, j, n;
	double v[NAXIS], w[NAXIS];
	double accum;
	struct rotdata r;
	double M[NAXIS][NAXIS];

    assert(mb->n == tp-> n && mb->n > 0);
              	
	getrotation(p, &r);	
  	
   	for (i = 0; i < NAXIS; i++) 
  	  	for (j = 0; j < NAXIS; j++)     	
      	  	M[i][j] = r.mat[i][j]; 
      	  	
	accum = 0.0;
	
    n = mb->n;
	if (n > tp->n)
        n = tp->n;
             	
    for (d = 0; d < n; d++) {
    	
     	atmcpy(mb->atm[d], v);
   	
     	/* center ATOM to rotate */
	    /*atmbck (r.ctr, v);*/
	    
 	    /* matrix multiply vector v, result in w */
  	    for (i = 0; i < NAXIS; i++) {
  	  		double tmp = 0.0;
       		for (k = 0; k < NAXIS; k++ ) {
            	tmp += M[i][k] * v[k];
            } 		
            w[i] = tmp;
      	}		
      	/* center and translate */
		/*atmtns (r.ctr, w);*/
		
		atmtns (r.tns, w);
  		/*------------*/
  
	    accum += vector_sqdist (tp->atm[d], w);
	
	}	
	return accum/n;
	
}

void atmcpy (const double *in, double *out)
{
	out[0] = in[0];
	out[1] = in[1];
 	out[2] = in[2];	
}

void atmtns (const double *delta, double *v)
{
	v[0] += delta[0];
	v[1] += delta[1];
	v[2] += delta[2]; 	
} 

void atmbck (const double *delta, double *v)
{
	v[0] -= delta[0];
	v[1] -= delta[1];
	v[2] -= delta[2]; 	
} 

double vector_sqdist(const double *v, const double *w)
{	
	double dx, dy, dz;
	dx = v[0] - w[0];
	dy = v[1] - w[1];
	dz = v[2] - w[2];
	return dx*dx + dy*dy + dz*dz;
}


static void coor_masscenter(const struct coordinates *tp, double *c)
{
	int d, n;
	double ac_x, ac_y, ac_z;
  	
  	n = tp->n;
  	ac_x = ac_y = ac_z = 0;
   	for (d = 0; d < n; d++) {
    	ac_x += tp->atm[d][0];
    	ac_y += tp->atm[d][1];    
    	ac_z += tp->atm[d][2];    
    }
    c[0] = ac_x / n;
    c[1] = ac_y / n;
    c[2] = ac_z / n;
}

void coor_rotate (const struct coordinates *ori, const struct rotdata *r, struct coordinates *to)
{
	int d, i, k, j;
	double v[NAXIS], w[NAXIS];
//	const double *v_ori;
	double M[NAXIS][NAXIS];

  	for (i = 0; i < NAXIS; i++) 
  	  	for (j = 0; j < NAXIS; j++)     	
      	  	M[i][j] = r->mat[i][j]; 	

    to->n = ori->n;

    for (d = 0; d < ori->n; d++) {
   
       	atmcpy (ori->atm[d], v); 	
	    /*atmbck (r->ctr, v);*/
	
 	    /* matrix multiply vector v, result in w */
  	    for (i = 0; i < NAXIS; i++) {
  	  		double tmp = 0.0;
       		for (k = 0; k < NAXIS; k++ ) {
            	tmp += M[i][k] * v[k];
            } 		
            w[i] = tmp;
      	} 	
    
    	/* atmtns (r->ctr, w);*/
    	atmtns (r->tns, w);
        atmcpy (w, to->atm[d]);
     }	

}


extern double
fit (const struct coordinates *t, const struct coordinates *m, struct transrot *out)
{	
	double p[NPARAM] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	
	const struct coordinates *ptbuf;
	const struct coordinates *pmbuf;
	
	double z, z0;
	double scale;
//	double qfactor;
	int i;
	
	double ct[3], cm[3];
	double t_mv[NPARAM] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double m_mv[NPARAM] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};	
	struct rotdata t_rot, m_rot;
	
	struct rotdata r;
	struct transrot tr;
	/* search for optimum parameter list to match t and m */

	scale = 1.0;
	
#if 1
	coor_masscenter(t, ct);
	coor_masscenter(m, cm);
	
	
	t_mv[3] = -ct[0];
	t_mv[4] = -ct[1];
	t_mv[5] = -ct[2];
	m_mv[3] = -cm[0];
	m_mv[4] = -cm[1];
	m_mv[5] = -cm[2];

	
	getrotation(t_mv, &t_rot);
	getrotation(m_mv, &m_rot);	
	
	coor_rotate(t, &t_rot, &t_buffer);
	coor_rotate(m, &m_rot, &m_buffer);
	
	ptbuf = /*t;*/   &t_buffer; 
	pmbuf = /* m;*/  &m_buffer;
#else
	ptbuf = t;   
	pmbuf = m;
#endif
	
	

	
#if 1
	z0 = iterfit(ptbuf, pmbuf, &scale, p);

	for (i = 0; i < 2000; i++) {

		z = iterfit(ptbuf, pmbuf, &scale, p);

		if (i != 0) {
			if (sqrt(z0) < 0.0000001) {
				z = z > z0? z0: z; 
				break;
			}
  			if (
				z < z0 
				&&	((z0-z)/z0) < 0.0000000001) {

				break;	
			}
  			z0 = z;
		} 
		else {	z0 = z;}

	}
#else
	
	for (i = 0, qfactor = 1.0; i < 1000 && qfactor > 0.000001; i++) {
		z = iterfitQ(ptbuf, pmbuf, &qfactor, p);
	}	
	
#endif

	getrotation (p, &r);
	get_transrot (&r, &tr);

	/*showtransrot(&tr);*/
	/*showrot(&r);*/
	
	for (i = 0; i < NAXIS; i++)
		cm[i] = -cm[i];
		
	compose_r_t(&tr, cm);
	compose_t_r(ct, &tr);
	/*showtransrot(&tr);*/
	transrotcpy (&tr, out);

	return sqrt(z);
}


double iterfit (const struct coordinates *t, 
				const struct coordinates *m, double *s_io, double *p_io)
{
	enum {SUBLIMIT = 100};
	enum {LOW = 0, MED = 1, HIG = 2, NSEEDS = 3};
	int best;
	double sc;
	double sci  [NSEEDS] = {0.05, 0.10, 0.20};
	double adj  [NSEEDS] = {0.8, 1.0, 1.2};
	double fi   [NSEEDS];
	double pi   [NSEEDS] [NPARAM];
	double pori [NPARAM];
	double dp   [NPARAM];
	double pS   [NPARAM];
	double pN   [NPARAM];
	double p	[NPARAM];
 	double fS, fN; 	
 	double sS, sN; 	
	int i, c, counter;
	double po;
	double finc, fdec;
	double delta[NPARAM] ={	0.00000001, 0.00000001, 0.00000001,
        					0.00000001, 0.00000001, 0.00000001 };

double fS_ori;

	/************/


	for (i = 0; i < NSEEDS; i++) 
 		sci [i] = *s_io * adj[i];
 	
 	prmcpy(p_io, pori);                     /* preserves original values */
 	prmcpy(p_io, p);

  	
	for (c = 0; c < NPARAM; c++) {
		po = p[c];
  		p[c] = po + delta[c];		
  		finc = sqdev (t, m, p);             /* incremented */
  		p[c] = po - delta[c];
  		fdec = sqdev (t, m, p);             /* decreased */
    	p[c] = po;                          /* restores paramenter */
     	dp[c] = (finc - fdec)/(2*delta[c]); /* derivative */ 		
     	dp[c] = -dp[c];                     /* nega  value of the derivative */
	}

	prmcpy(pori, pS);
 	prmsca(dp, 0, pS); 		/* apply dp to pS, scaled by 0. No changes here...*/
 	fS = sqdev (t, m, pS); 	/* fS = starting deviation */
	
	fS_ori = fS;
 
#ifdef INSPECT
 	showvector6("dp calculated in iter fit", dp);
	showvector6("p input in iterfit", pori);
 	printf("START: %.5lf\n", sqrt(fS));
	system("PAUSE");
#endif	

	for (i = 0; i < NSEEDS; i++) {
 		prmcpy(pori, pi[i]);
 		prmsca(dp, sci[i], pi[i]);
 		fi[i] = sqdev (t, m, pi[i]); 
  	}

#ifdef INSPECT     
    printf("\nLMH\n");
	for (i = 0; i < NSEEDS; i++) {
		printf("%9.5lf:  %8.5lf\n", sci[i], sqrt(fi[i]));
  	}     
#endif
     
	if (fS < fi[LOW] && fS < fi[MED] && fS < fi[HIG]) {
		/* The starting point is the best, so, it will search*/
		/* for much lower sN factors until it finds one that is */
		/* lower than the starting point */

			sN /= 2;
			prmcpy(pori, pN);
			prmsca(dp, sN, pN);               
			fN = sqdev (t, m, pN);	

		for (i = 0, sN = sci[LOW]; i < SUBLIMIT && fS < fN; i++) {
			sN /= 2;
			prmcpy(pori, pN);
			prmsca(dp, sN, pN);               
			fN = sqdev (t, m, pN);			
		}

assert(fN <= fS);
assert(sqdev (t, m, pN) <= sqdev (t, m, pS));

	 	sc = 0.8; /* adjusting factor */
	 	prmcpy(pN, pS);
   		fS   = fN;
   	   	sS   = sN;	 	
    }
    else 
   	if (fi[MED] <= fi[LOW] && fi[MED] <= fi[HIG]) {
   		prmcpy(pi[MED], p_io);
   		*s_io = sci[MED];   		

assert(fi[MED] <= fS_ori);
assert(sqdev (t, m, pi[MED]) <= sqdev (t, m, pS));

   		return fi[MED];
    } else if (fi[HIG] <= fi[LOW]) {
    	best = HIG;
    	sc   =adj[best];
    	prmcpy(pi[best], pS);
   		fS   = fi[best];
   		sS   =sci[best];
    } else {
    	best = LOW;
   	    sc   =adj[best];
   	    prmcpy(pi[best], pS);
   	    fS   = fi[best];
   	    sS   =sci[best];    	
    }	

assert(fS <= fS_ori);
        	
	prmcpy(pS, pN);                           /* pN starts with pS values */
	fN = fS;
	sN = sS;
	
 	sN *= sc; 	
	prmcpy(pori, pN);
  	prmsca(dp, sN, pN);                       /* scale param Next */
	fN = sqdev (t, m, pN);

	counter = 0; 

#ifdef INSPECT	
	printf("counting\n");
	printf("%5d:  %9.5lf: %8.5lf\n", counter, sN, sqrt(fN));
#endif



	  		
    while (fN < fS && counter++ < SUBLIMIT) {
 	   		prmcpy(pN, pS);                   /* pS advances      */
 	   		fS = fN;
 	   		sS = sN;
 	   		
 	   		sN *= sc; 	
	 	   	prmcpy(pori, pN);
  		 	prmsca(dp, sN, pN);               /* scale param Next */
	  		fN = sqdev (t, m, pN);
  
#ifdef INSPECT  
	printf("%5d:  %9.5lf: %8.5lf\n", counter, sN, sqrt(fN));
#endif
             	
	}	

assert(fN >= fS || 0==fprintf(stderr,"counter=%d\n, fN=%lf, fS%lf",counter,fN,fS));
assert(sqdev (t, m, pN) >= sqdev (t, m, pS));

#ifdef INSPECT
	printf("end\n");
	printf("%5d:  %9.5lf: %8.5lf\n", counter, sS, sqrt(fS));	
 	printf("\n");
 	system("PAUSE");
#endif

	assert(fS <= fS_ori);
	assert(sqdev (t, m, pS) <= sqdev (t, m, pori)); 

	if (fS <= fS_ori) {
		prmcpy(pS, p_io);
		*s_io = sS;
		return fS;
	} else {
	 	prmcpy(pori, p_io);                     /* preserves original values */
		return fS_ori;
	}
}


void prmcpy(const double *a, double *b)
{
	int i;
	for (i = 0; i < NPARAM; i++) {
		b[i] = a[i];
	}
}


void prmsca(const double *d, double s, double *b)
{
	int i;
	for (i = 0; i < NPARAM; i++) {
		b[i] += s * d[i];
	}
}


void showvector6(char *s, double *v)
{
	printf(	"%s\n"
 			"%15.12lf %15.12lf %15.12lf\n"
			"%15.12lf %15.12lf %15.12lf\n\n"	
   			, s
      		, v[0], v[1], v[2]
   			, v[3], v[4], v[5]);	      	
}


double iterfitQ(const struct coordinates *t, 
				const struct coordinates *m, double *f_io, double *p_io)
{
#define QLIMIT 0.000001	
	enum {SUBLIMIT = 100};
	
	double pori [NPARAM];
	double dp   [NPARAM];
	double pS   [NPARAM];
	double pN   [NPARAM];
	double p	[NPARAM];
 	double fS, fN; 	
	int i, c, iter;
	double qfactor;
	double po;
	double q[NPARAM];
	double Q[NPARAM] ={		1.0, 1.0, 1.0,
        					5.0, 5.0, 5.0 };
	double delta[NPARAM] ={	0.00000001, 0.00000001, 0.00000001,
        					0.00000001, 0.00000001, 0.00000001 };      
    double finc, fdec;  					
	/************/
      
    qfactor = *f_io;  
    	
 	prmcpy(p_io, pori);                     	/* preserves original values */
 	prmcpy(p_io, p);
	
  	/*starting values without modification */
	prmcpy(p, pS);
 	fS = sqdev (t, m, pS); 

fN = fS; //in case there is no iterations;

 	for (iter = 0; iter < SUBLIMIT; iter++) {
	
  		/* obtain derivative dp */
  		for (c = 0; c < NPARAM; c++) {
  			po = p[c];
  			p[c] = po + delta[c];		
  			finc = sqdev (t, m, p);             /* incremented */
  			p[c] = po - delta[c];
  			fdec = sqdev (t, m, p);             /* decreased */
    		p[c] = po;                          /* restores paramenter */
     		dp[c] = (finc - fdec)/(2*delta[c]); /* derivative */ 		
     		dp[c] = -dp[c];                     /* nega value of the deriv */
     	}

        while (qfactor > QLIMIT) {
 			for (i = 0; i < NPARAM; i++) {
 				if (dp[i] < 0) {
 					q[i] = -Q[i] * qfactor;
				} else {
 					q[i] =  Q[i] * qfactor;
				}
  			}

	  		prmcpy(p, pN);
 			prmsca(q, 1, pN);
 			fN = sqdev (t, m, pN); 	 
  
#ifdef INSPECT 
	if (iter%10 == 0) 
	printf("i=%5d q=%9.5lf rmsd=%8.5lf\n", iter, qfactor, sqrt(fN));

#endif
 			if (fN < fS) {
 				break;
 	 		} else {
 				qfactor /= 3;
 	 		}
        }

 		prmcpy(pN, p);
  		prmcpy(pN, pS);
  		fS = fN;
  	}

	prmcpy(pS, p_io);
	*f_io = qfactor;
	return fS;
}

void get_transrot (struct rotdata *r, struct transrot *out) 
{
	int i, j, k;
	struct transrot a, b, c;
	int WIDTH = NAXIS + 1;
	
	for (i = 0; i < WIDTH; i++) {
		for (j = 0; j < WIDTH; j ++) {
			a.m[i][j] = 0;
   			b.m[i][j] = 0;  		
    	} 
 	}
	for (i = 0; i < WIDTH; i++) {
		a.m[i][i] = 1.0;
		b.m[i][i] = 1.0;		
 	}
	for (i = 0; i < NAXIS; i++) {
		for (j = 0; j < NAXIS; j ++) {
			a.m[i][j] = r->mat[i][j];  		
    	} 
 	} 	
	for (i = 0; i < NAXIS; i++) {
		b.m[i][NAXIS] = r->tns[i];
 	} 
	
  	/* Multiply c = b * a */
  	for (i = 0; i < WIDTH; i++) {
  	  	for (j = 0; j < WIDTH; j++) {    	
  	  		c.m[i][j] = 0.0;
      	}
    }
  	
  	for (i = 0; i < WIDTH; i++) {
  	  	for (j = 0; j < WIDTH; j++) {    	
  	  		double tmp = 0.0;
       		for (k = 0; k < WIDTH; k++ ) {
            	tmp += b.m[i][k] * a.m[k][j];
            } 		
            c.m[i][j] = tmp;
      	}
    }  	
    
    /* OUTPUT */
	for (i = 0; i < WIDTH; i++) {
		for (j = 0; j < WIDTH; j ++) {
			out->m[i][j] = c.m[i][j];
  		
    	} 
 	}
}

void compose_r_t (struct transrot *io, const double *v)
{
	int i, j, k;
	struct transrot a, b, c;
	int WIDTH = NAXIS + 1;

	for (i = 0; i < WIDTH; i++) {
		for (j = 0; j < WIDTH; j ++) {
   			a.m[i][j] = io->m[i][j];  		
    	} 
 	}
	
	for (i = 0; i < WIDTH; i++) {
		for (j = 0; j < WIDTH; j ++) {
   			b.m[i][j] = 0;  		
    	} 
 	}
	for (i = 0; i < WIDTH; i++) {
		b.m[i][i] = 1.0;		
 	}
	for (i = 0; i < NAXIS; i++) {
		b.m[i][NAXIS] = v[i];
 	} 
	
  	/* Multiply c = a * b */
  	for (i = 0; i < WIDTH; i++) {
  	  	for (j = 0; j < WIDTH; j++) {    	
  	  		c.m[i][j] = 0.0;
      	}
    }
  	
  	for (i = 0; i < WIDTH; i++) {
  	  	for (j = 0; j < WIDTH; j++) {    	
  	  		double tmp = 0.0;
       		for (k = 0; k < WIDTH; k++ ) {
            	tmp += a.m[i][k] * b.m[k][j];
            } 		
            c.m[i][j] = tmp;
      	}
    }  	
    
    /* OUTPUT */
	for (i = 0; i < WIDTH; i++) {
		for (j = 0; j < WIDTH; j ++) {
			io->m[i][j] = c.m[i][j];
  		
    	} 
 	}
	
}

void compose_t_r (const double *v, struct transrot *io)
{
	int i, j, k;
	struct transrot a, b, c;
	int WIDTH = NAXIS + 1;

	for (i = 0; i < WIDTH; i++) {
		for (j = 0; j < WIDTH; j ++) {
   			a.m[i][j] = io->m[i][j];  		
    	} 
 	}
	
	for (i = 0; i < WIDTH; i++) {
		for (j = 0; j < WIDTH; j ++) {
   			b.m[i][j] = 0;  		
    	} 
 	}
	for (i = 0; i < WIDTH; i++) {
		b.m[i][i] = 1.0;		
 	}
	for (i = 0; i < NAXIS; i++) {
		b.m[i][NAXIS] = v[i];
 	} 
	
  	/* Multiply c = a * b */
  	for (i = 0; i < WIDTH; i++) {
  	  	for (j = 0; j < WIDTH; j++) {    	
  	  		c.m[i][j] = 0.0;
      	}
    }
  	
  	for (i = 0; i < WIDTH; i++) {
  	  	for (j = 0; j < WIDTH; j++) {    	
  	  		double tmp = 0.0;
       		for (k = 0; k < WIDTH; k++ ) {
            	tmp += b.m[i][k] * a.m[k][j];
            } 		
            c.m[i][j] = tmp;
      	}
    }  	
    
    /* OUTPUT */
	for (i = 0; i < WIDTH; i++) {
		for (j = 0; j < WIDTH; j ++) {
			io->m[i][j] = c.m[i][j];
  		
    	} 
 	}
	
}

void do_transrot (const struct transrot *inp, double *io)
{
	double v[NAXIS+1];
	double w[NAXIS+1];	
	int i;
//	int j;
	int k;
//	struct transrot a, b, c;
	int WIDTH = NAXIS + 1;

	for (i = 0; i < NAXIS; i++) {
		v[i] = io[i];
	}
	v[NAXIS] = 1.0;


 	    /* matrix multiply vector v, result in w */
  	    for (i = 0; i < WIDTH; i++) {
  	  		double tmp = 0.0;
       		for (k = 0; k < WIDTH; k++ ) {
            	tmp += inp->m[i][k] * v[k];
            } 		
            w[i] = tmp;
      	} 	

  

    /* OUTPUT */
	for (i = 0; i < NAXIS; i++) {
		io[i] = w[i];
 	}
	
}

void transrotcpy (const struct transrot *inp, struct transrot *out)
{
	int i, j;
	int WIDTH = NAXIS + 1;

    
    /* OUTPUT */
	for (i = 0; i < WIDTH; i++) {
		for (j = 0; j < WIDTH; j ++) {
			out->m[i][j] = inp->m[i][j];
    	} 
 	}	
}

void showrot (const struct rotdata *inp)
{
	int i, j;
//	int WIDTH = NAXIS + 1;

    for (i = 0; i < NAXIS; i++) {
       		for (j = 0; j < NAXIS; j++ ) {
       			printf("%.5lf  ",inp->mat[i][j]);
            } 
            printf("\n");		
  	} 	
  	printf("%.5lf  %.5lf  %.5lf  ", inp->tns[0], inp->tns[1], inp->tns[2]);
   	printf("\n");
}

void showtransrot (const struct transrot *inp)
{
	int i, j;
	int WIDTH = NAXIS + 1;

    for (i = 0; i < WIDTH; i++) {
       		for (j = 0; j < WIDTH; j++ ) {
       			printf("%.5lf  ",inp->m[i][j]);
            } 
            printf("\n");		
  	} 	
   	printf("\n");
}
