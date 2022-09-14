#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

/** 
I commented out the following lines so that -Wall option of
gcc would not crib about a defined variable not being used.
I might as well have removed it. But, if I ever happen to use
Numerical Recipes code, it might be needed.
**/

/*
static float sqrarg;
#define SQR(a) ((sqrarg =(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
*/

#if defined (__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

void nrerror (char error_text[]);
float *vector(long nl, long nh);
int *ivector (long nl, long nh);
unsigned char *cvector (long nl, long nh);
unsigned long *lvector (long nl, long nh);
double *dvector (long nl, long nh);
float **matrix (long nrl, long nrh, long ncl, long nch);
double **dmatrix (long nrl, long nrh, long ncl, long nch);
int **imatrix (long nrl, long nrh, long ncl, long nch);

void free_vector (float *v, long nl, long nh);
void free_ivector (int *v, long nl, long nh);
void free_cvector (unsigned char *v, long nl, long nh);
void free_lvector (unsigned long *v, long nl, long nh);
void free_dvector (double *v, long nl, long nh);
void free_matrix (float **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix (double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix (int **m, long nrl, long nrh, long ncl, long nch);

#else /* ANSI */
/* traditional K&R */

void nrerror ();

/* etc, etc, etc.
// since we are not likely to use the earlier, and outdated,
// K&R version of 'C', this stuff has been left out.
*/
#endif /* ANSI */

#endif /* _NR_UTILS_H_ */


