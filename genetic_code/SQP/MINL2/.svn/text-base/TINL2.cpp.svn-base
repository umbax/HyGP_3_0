/* TINL2.F -- translated by f2c (version 20050501).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/
#include <iostream>  // basic i/o commands

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"



/* Table of constant values */

static integer c__0 = 0;
static integer c__1 = 1;
static integer n = 2;  // n: number of unknowns (parameters)   c__2
static integer m = 3;  // m: number of functions (terms of the sum)  c__3
static integer c__195 = 195;


/*      TEST OF MINL2     21.11.1991 */


int MAIN__()   //<----refers to the master.cpp!
{
    extern /* Subroutine */ int opti_(integer *);

    //opti_(&c__0);  // to check the gradient
    opti_(&c__1);  // to perform optimization 
    cout << "\nDone...\n";
	return 0;
}


int opti_(integer *method)
{
    /* Format strings */
    static char fmt_10[] = "(\002 INPUT ERROR. PARAMETER NUMBER \002,i1,\002 IS OUTSIDE ITS RANGE.\002)";
    static char fmt_21[] = "(\002 TEST OF GRADIENTS \002//\002 MAXIMUM FORWA\
RD  DIFFERENCE: \002,d8.2,\002 AT FUNCTION NO \002,i1,\002 AND VARIABLE NO\
 \002,i1/\002 MAXIMUM BACKWARD DIFFERENCE: \002,d8.2,\002 AT FUNCTION NO \
\002,i1,\002 AND VARIABLE NO \002,i1/\002 MAXIMUM CENTRAL  DIFFERENCE: \002,\
d8.2,\002 AT FUNCTION NO \002,i1,\002 AND VARIABLE NO \002,i1)";
    static char fmt_22[] = "(//\002 MAXIMUM ELEMENT IN DF: \002,d8.2)";
    static char fmt_30[] = "(\002 UPPER LIMIT FOR FUNCTION EVALUATIONS EXCEEDED.\002/)";
    static char fmt_31[] = "(23x,\002 SOLUTION: \002,d18.10/d52.10//\002 NUM\
BER OF CALLS OF FDF: \002,i4//\002 FUNCTION VALUES AT THE SOLUTION: \002,d18.10,2(/d52.10))";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();

    /* Local variables */
    static integer i__, j, k;
    static doublereal w[195];    // size of the workspace: don't touch it!
	static doublereal x[2];    // size of the parameters vector: get it as input
	static doublereal dx;  //<----------------------------------------------------
    extern /* Subroutine */ int fdf_(...);   // 
    static doublereal eps;
    extern /* Subroutine */ int minl2_(U_fp, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *,
	     integer *);
    static integer index[8]	/* was [4][2] */;  //?????
    static logical optim;
    static integer icontr, maxfun;

    /* Fortran I/O blocks */
    static cilist io___8 = { 0, 6, 0, fmt_10, 0 };
    static cilist io___11 = { 0, 6, 0, fmt_21, 0 };
    static cilist io___12 = { 0, 6, 0, fmt_22, 0 };
    static cilist io___13 = { 0, 6, 0, fmt_30, 0 };
    static cilist io___14 = { 0, 6, 0, fmt_31, 0 };



/*     SET PAPAMETERS */
    eps = 1e-10;
    maxfun = 25;
/*     SET INITIAL GUESS */
    x[0] = 1.;  //initial guess of x1
    x[1] = 1.; 	//initial guess of x2
    icontr = *method; // <---------------------------------------------------------
    optim = icontr > 0;
    if (! optim) {
	dx = .001;
    }

    minl2_(	(U_fp)fdf_,
		&n,	//n : no. of unknown parameters - dimension of x
		&m,	// m : no. of functions, terms of the sum
		x,	// x : initial guess of the parameters 
		&dx,	// dx : damping factor (10E-4 if initial guess close to solution, otherwise 1)
		&eps,	// eps : accuracy
		&maxfun,// maxfun : max no. of calls to the function
		w,	// workspace
		&c__195,// dimension of the workspace: don't touch it	
		&icontr); // icontr : used to control the function

    if (icontr < 0) {
	goto L100;
    }
    if (! optim) {
	goto L200;
    }
    goto L300;

/*     PARAMETER OUTSIDE RANGE */
L100:
    s_wsfe(&io___8);
    i__1 = -icontr;
    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
    e_wsfe();
    goto L400;

/*     RESULTS FROM GRADIENT TEST */
L200:
    for (k = 2; k <= 4; ++k) {
	index[k - 1] = (integer) w[k * 2];
	index[k + 3] = (integer) w[(k << 1) + 1];
/* L20: */
    }
    s_wsfe(&io___11);
    for (k = 2; k <= 4; ++k) {
	do_fio(&c__1, (char *)&w[k - 1], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&index[k - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&index[k + 3], (ftnlen)sizeof(integer));
    }
    e_wsfe();
    s_wsfe(&io___12);
    do_fio(&c__1, (char *)&w[0], (ftnlen)sizeof(doublereal));
    e_wsfe();
    goto L400;

/*     RESULTS FROM OPTIMIZATION */
L300:
    if (icontr == 2) {
		s_wsfe(&io___13);
		e_wsfe();
    }
	// only this part is executed if icontr =1
    s_wsfe(&io___14);
    // this cycle controls the independent variables x
	for (i__ = 1; i__ <= 2; ++i__) {       //2 is the number of independent variables
		do_fio(&c__1, (char *)&x[i__ - 1], (ftnlen)sizeof(doublereal));
    }
    do_fio(&c__1, (char *)&maxfun, (ftnlen)sizeof(integer));
	  // this cycle controls the dependent functions f   
	for (j = 1; j <= 3; ++j) {      //3 is the number of functions the F sum is made of
		do_fio(&c__1, (char *)&w[j - 1], (ftnlen)sizeof(doublereal));
    }
    e_wsfe(); //this function prints "SOLUTION"

L400:
	printf ("\nSolution: %f  %f      Solution size n: %i", x[0], x[1],  sizeof(x));
	printf ("\nValues of the single f : %E %E %E", w[0],  w[1], w[2] );

    return 0;

}



// Subroutine fdf : this has to be defined by the user
// n : no. of unknown parameters - dimension of x (in)
// m : no. of functions, terms of the sum (in)
// x : initial guess of the parameters (in)
// df : vector containing the values of the derivatives in x (out)
// f : vector containing the values of the functions in x (out)
 int fdf_(integer *n, integer *m, doublereal *x, doublereal *df, doublereal *f)
{
    /* System generated locals */
    int df_dim1, df_offset;
    double d__1;

    /* Parameter adjustments */
    --x;
    --f;
    df_dim1 = *m;
    df_offset = 1 + df_dim1;
    df -= df_offset;

    /* Function Body */
    // functions
    f[1] = 1.5 - x[1] * (1. - x[2]);	//F1			
    f[2] = 2.25 - x[1] * (1. - x[2] * x[2]);	//F2
    f[3] = 2.625 - x[1] * (1. - x[2] * (x[2] * x[2]));	//F3
    
    // derivatives
    df[df_dim1 + 1] = x[2] - 1.;  	// dF1/dx1
    df[(df_dim1 << 1) + 1] = x[1];	// dF1/dx2       << = shift left?
    df[df_dim1 + 2] = x[2] * x[2] - 1.;		// dF2/dx1
    df[(df_dim1 << 1) + 2] = x[1] * 2. * x[2];	// dF2/dx2
    df[df_dim1 + 3] = x[2] * (x[2] * x[2]) - 1.;	// dF3/dx1
    df[(df_dim1 << 1) + 3] = x[1] * 3. * (x[2] * x[2]);	// dF3/dx2

    return 0;
} 


#ifdef __cplusplus
	}
#endif
