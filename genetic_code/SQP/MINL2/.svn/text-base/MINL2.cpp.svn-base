/* MINL2.F -- translated by f2c (version 20050501).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#ifdef __cplusplus
extern "C" {
#endif
//#include "f2c.h"
/* Table of constant values */
static integer COMMENT = 1;
static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__5 = 5;
static doublereal c_b76 = 1.;
static doublereal c_b95 = 0.;
static integer c__0 = 0;
static integer c__2 = 2;
static doublereal c_b147 = -1.;
static integer c__4 = 4;
static doublereal c_b235 = .33333333333333331;
static integer c__6 = 6;
static integer c__27 = 27;
// in the following definition, P and ntree are inputs added to avoid global variables. They are passed recursively
// through all the functions called by minl2_ to reach the final destination fdf_ and then allow to call the 
// Population member functions from there. To know which are these intermediate functions find "P and ntree added".
// Could exist a simpler way... static variables as defined above? I am trying this, search for "P and tree removed" instead 
// of "P and ntree added" to find where those variables were added
/* Subroutine */ int minl2_(U_fp fdf, integer *n, integer *m, doublereal *x, 
	doublereal *dx, doublereal *eps, integer *maxfun, doublereal *w, 
	integer *iw, integer *icontr) //P and ntree removed
{	
	if (COMMENT) {
	printf("\n\nminl2_ : Entered minl2_");
	//printf("\nminl2_ : P = %i, \nminl2_ : ntree = %i", P, ntree);
 	}   
	
	
	/* Initialized data */

    static integer uiparm[9] = { 0,0,0,0,0,0,0,0,0 };

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle();

    /* Local variables */
    static doublereal d__[4];
    static integer i__, id[6], na, nf, il, iv[1000], nr, nw, na1, nf1, nf2, 
	    maxf;
    extern /* Subroutine */ int calcj_(...), calcr_(...);
    static logical optim;
    extern /* Subroutine */ int nl2sol_(integer *, integer *, doublereal *, 
	    U_fp, U_fp, integer *, doublereal *, integer *, doublereal *, 
	    U_fp), checkd_(U_fp, integer *, integer *, doublereal *, //P and ntree removed: nl2sol_ :,Population *, integer)
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *), dfault_(integer *, //P and ntree removed: checkd_ :,Population *, integer)
	    doublereal *);
    static logical conprt;
    static real ivprop;

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 6, 0, 0, 0 };
    static cilist io___19 = { 0, 6, 0, 0, 0 };
    static cilist io___20 = { 0, 6, 0, 0, 0 };
    static cilist io___21 = { 0, 6, 0, 0, 0 };
    static cilist io___22 = { 0, 6, 0, 0, 0 };
    static cilist io___23 = { 0, 6, 0, 0, 0 };
    static cilist io___24 = { 0, 6, 0, 0, 0 };
    static cilist io___25 = { 0, 6, 0, 0, 0 };



/*   UNCONSTRAINED MINIMIZATION OF THE L2-NORM OF A VECTOR */
/*   FUNCTION. */
/*   THIS SUBROUTINE USES A SLIGHTLY MODIFIED VERSION OF THE SUBROUTINE */
/*   NL2SOL (AVAILABLE FROM NETLIB) AS DESCRIBED IN: */
/*   DENNIS, J.E., GAY, D.M., AND WELSCH, R.E. (1981), */
/*   'AN ADAPTIVE NONLINEAR LEAST-SQUARES ALGORITHM', ACM TRANS. MATH. */
/*   SOFTWARE, VOL. 7, NO. 3., PP. 348-368. */
/*   DENNIS, J.E., GAY, D.M., AND WELSCH, R.E. (1981), */
/*   'NL2SOL - AN ADAPTIVE NONLINEAR LEAST-SQUARES ALGORITHM', */
/*   ACM TRANS. MATH. SOFTWARE, VOL. 7, NO. 3., PP. 369-383. */

/*   FOR A PROGRAM DESCRIPTION SEE: */
/*   KAJ MADSEN,OLE TINGLEFF, PER CHRISTIAN HANSEN & WOJCIECH OWCZARZ: */
/*   'ROBUST SUBROUTINES FOR NON-LINEAR OPTIMIZATION', REPORT NI 90-06, */
/*   INSTITUTE FOR NUMERICAL ANALYSIS, TECHNICAL UNIVERSITY OF DENMARK, */
/*   APRIL 1990. */

/*     IV INTEGER ARRAY OF LENGTH = N + 60 */


/*     W REAL ARRAY OF LENGTH = 93 + M(N + 3) + N(3N + 33)/2 */

    /* Parameter adjustments */
    --x;
    --w;

    /* Function Body */
    uiparm[2] = 0;
    conprt = FALSE_;
    if (conprt) {
	s_wsle(&io___3);
	do_lio(&c__9, &c__1, " ====== MINL2  START ======", (ftnlen)27);
	e_wsle();
    }
    optim = *icontr > 0;
    il = *m * (*n + 3) + 93 + *n * (*n * 3 + 33) / 2 + *n * *m + *m;
    if (*n <= 0) {
	*icontr = -2;
    }
    if (*m <= 0) {
	*icontr = -3;
    }
    if (optim) {
	goto L2;
    }
    if (*dx == 0.) {
	*icontr = -5;
    }
L2:
    if (*eps <= 0.) {
	*icontr = -6;
    }
    if (*maxfun <= 0) {
	*icontr = -7;
    }
    if (*iw < il) {
	*icontr = -9;
    }
    if (*icontr < 0) {
	return 0;
    }
/*   ****   EXIT, ERROR IN THE PARAMETERS   *** */

/* L6: */
    if (optim) {
	goto L9;
    }
/*   ****   CHECK OF THE GRADIENTS          *** */
    na = 1;
    na1 = na + *n * *m;
    nf = na1 + *n * *m;
    nf1 = nf + *m;
    nf2 = nf1 + *m;
    checkd_((U_fp)fdf, n, m, &x[1], dx, d__, id, &w[na], &w[na1], &w[nf], &w[
	    nf1], &w[nf2]);  //P and ntree removed
    for (i__ = 1; i__ <= 4; ++i__) {
/* L7: */
	w[i__] = d__[i__ - 1];
    }
    for (i__ = 1; i__ <= 6; ++i__) {
/* L8: */
	w[i__ + 4] = (doublereal) id[i__ - 1];
    }
    return 0;
/*   ****   EXIT, GRADIENTS CHECKED         *** */

L9:
/*   ****   THE OPTIMIZATION PROPER         *** */
    nr = 1;
    nw = *n * *m + *m + 1;
    maxf = *maxfun;
    iv[0] = 0;
    dfault_(iv, &w[nw]);

/*    SET MAX CALL FUNCTION */

    ivprop = (real) iv[16] / (real) (*maxfun);
    iv[16] = *maxfun;

/*    SET MAX ITERATION LIMIT */

    iv[17] = (integer)(iv[17] / ivprop) + 1;      //added (integer) manually

/* SET PRINT FLAGS FOR NL2SOL OF */

    if (*icontr > 1) {
	goto L1217;
    }
    iv[18] = 0;
    iv[19] = 0;
    iv[20] = 0;
    iv[21] = 0;
    iv[22] = 0;
    iv[23] = 0;
L1217:

    w[nw + 33] = *eps;
    iv[0] = 12;
    for (i__ = 1; i__ <= 147; ++i__) {
	if (conprt && i__ <= 62) {
	    s_wsle(&io___19);
	    do_lio(&c__3, &c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&iv[i__ - 1], (ftnlen)sizeof(integer)
		    );
	    do_lio(&c__5, &c__1, (char *)&w[i__ + nw], (ftnlen)sizeof(
		    doublereal));
	    e_wsle();
	}
	if (conprt && i__ > 62) {
	    s_wsle(&io___20);
	    do_lio(&c__3, &c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_lio(&c__9, &c__1, "---", (ftnlen)3);
	    do_lio(&c__5, &c__1, (char *)&w[i__ + nw], (ftnlen)sizeof(
		    doublereal));
	    e_wsle();
	}
/* L2408: */
    }
    if (conprt) {
	s_wsle(&io___21);
	do_lio(&c__9, &c__1, " ====== ENTRY NL2SOL =====", (ftnlen)26);
	e_wsle();
    }
    nl2sol_(m, n, &x[1], (U_fp)calcr_, (U_fp)calcj_, iv, &w[nw], uiparm, &w[
	    nr], (U_fp)fdf);  //P and ntree removed
    *icontr = 0;
    if (conprt) {
	s_wsle(&io___22);
	do_lio(&c__9, &c__1, " ====== EXIT  NL2SOL =====", (ftnlen)26);
	e_wsle();
    }
    *maxfun = uiparm[2];
    for (i__ = 1; i__ <= 147; ++i__) {
	if (conprt && i__ <= 62) {
	    s_wsle(&io___23);
	    do_lio(&c__3, &c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&iv[i__ - 1], (ftnlen)sizeof(integer)
		    );
	    do_lio(&c__5, &c__1, (char *)&w[i__ + nw], (ftnlen)sizeof(
		    doublereal));
	    e_wsle();
	}
	if (conprt && i__ > 62) {
	    s_wsle(&io___24);
	    do_lio(&c__3, &c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_lio(&c__9, &c__1, "---", (ftnlen)3);
	    do_lio(&c__5, &c__1, (char *)&w[i__ + nw], (ftnlen)sizeof(
		    doublereal));
	    e_wsle();
	}
/* L2406: */
    }
    if (*maxfun >= maxf) {
	*icontr = 2;
    }
    if (conprt) {
	s_wsle(&io___25);
	do_lio(&c__9, &c__1, " ====== MINL2  EXIT ======", (ftnlen)26);
	e_wsle();
    }
    return 0;
/*   **** NORMAL EXIT             **** */
} /* minl2_ */


/* Subroutine */ int checkd_(S_fp fdf, integer *n, integer *m, doublereal *x, 
	doublereal *stepl, doublereal *diff, integer *indx, doublereal *a, 
	doublereal *a1, doublereal *f, doublereal *f1, doublereal *f2) //P and ntree removed
{
    /* System generated locals */
    integer a_dim1, a_offset, a1_dim1, a1_offset, i__1, i__2;

    /* Local variables */
    static doublereal c__;
    static integer i__, j;
    static doublereal s, db, dc, df;
    static integer ib, ic, jc, jb, jf;
    static logical chb, chf;
    static doublereal aji;
    static logical chm;
    static integer iif;
    static doublereal amax;
    extern /* Subroutine */ int diferr_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, logical *, logical *, logical *);

    /* Parameter adjustments */
    --x;
    --f2;
    --f1;
    --f;
    a1_dim1 = *m;
    a1_offset = 1 + a1_dim1;
    a1 -= a1_offset;
    a_dim1 = *m;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --diff;
    --indx;

    /* Function Body */
    amax = 0.;
    dc = 0.;
    df = 0.;
    db = 0.;
    (*fdf)(n, m, &x[1], &a[a_offset], &f[1]); //P and ntree removed
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	c__ = x[i__];
	x[i__] = c__ + *stepl;
	(*fdf)(n, m, &x[1], &a1[a1_offset], &f1[1]); //P and ntree removed
	x[i__] = c__ - *stepl * .5;
	(*fdf)(n, m, &x[1], &a1[a1_offset], &f2[1]); //P and ntree removed
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    aji = a[j + i__ * a_dim1];
	    s = abs(aji);
	    if (amax < s) {
		amax = s;
	    }
	    diferr_(&f[j], &f1[j], &f2[j], &aji, stepl, &dc, &df, &db, &chm, &
		    chf, &chb);
	    if (! chm) {
		goto L1;
	    }
	    jc = j;
	    ic = i__;
L1:
	    if (! chf) {
		goto L2;
	    }
	    jf = j;
	    iif = i__;
L2:
	    if (! chb) {
		goto L3;
	    }
	    jb = j;
	    ib = i__;
L3:
	    ;
	}
/* L4: */
	x[i__] = c__;
    }
    diff[1] = amax;
    diff[2] = df;
    diff[3] = db;
    diff[4] = dc;
    indx[1] = jf;
    indx[2] = iif;
    indx[3] = jb;
    indx[4] = ib;
    indx[5] = jc;
    indx[6] = ic;
    return 0;
} /* checkd_ */


/* Subroutine */ int diferr_(doublereal *ym, doublereal *yf, doublereal *yb, 
	doublereal *diff, doublereal *stepl, doublereal *erm, doublereal *erf,
	 doublereal *erb, logical *chm, logical *chf, logical *chb)
{
    static doublereal s;

    s = (*yf - *ym) / *stepl - *diff;
    *chf = abs(s) > *erf;
    if (*chf) {
	*erf = abs(s);
    }
    s = (*ym - *yb) * 2. / *stepl - *diff;
    *chb = abs(s) > *erb;
    if (*chb) {
	*erb = abs(s);
    }
    s = (*ym + *yf / 3. - *yb * 4. / 3.) / *stepl - *diff;
    *chm = abs(s) > *erm;
    if (*chm) {
	*erm = abs(s);
    }
    return 0;
} /* diferr_ */


/* Subroutine */ int calcr_(integer *n, integer *p, doublereal *x, integer *
	nf, doublereal *r__, integer *uiparm, doublereal *urparm, S_fp ufparm) //P and ntree removed
{
	if (COMMENT) {
	printf("\ncalcr_ : Entered calcr_");
	//printf("\ncalcr_ : P = %i \ncalcr_ : ntree = %i", P, ntree);
	}
    /* Parameter adjustments */
    --r__;
    --urparm;
    --x;
    --uiparm;

    /* Function Body */
    (*ufparm)(p, n, &x[1], &urparm[*n + 1], &r__[1]); //P and ntree removed
    uiparm[4] = 1;
    ++uiparm[3];
    return 0;
} /* calcr_ */


/* Subroutine */ int calcj_(integer *n, integer *p, doublereal *x, integer *
	nf, doublereal *j, integer *uiparm, doublereal *urparm, S_fp ufparm, Population *P, integer ntree)  
{	
	if (COMMENT) {
	printf("\n\ncalcj_ : Entered calcj_");
	printf("\ncalcj_ : P = %i, \ncalcj_ : ntree = %i", P, ntree);
	}	    
	/* System generated locals */
    integer i__1;

    /* Local variables */
    static integer kk;

    /* Parameter adjustments */
    --urparm;
    --j;
    --x;
    --uiparm;

    /* Function Body */
    if (uiparm[4] == 1) {
	uiparm[4] = 0;
	i__1 = *n * *p;
	for (kk = 1; kk <= i__1; ++kk) {
	    j[kk] = urparm[*n + kk];
/* L10: */
	}
    } else {
	(*ufparm)(p, n, &x[1], &j[1], &urparm[1]); //P and ntree removed 
	++uiparm[3];
    }
    return 0;
} /* calcj_ */


integer imdcon_(integer *k)
{
    /* Initialized data */

    static integer mdcon[3] = { 6,8,5 };

    /* System generated locals */
    integer ret_val;



/*  ***  RETURN INTEGER MACHINE-DEPENDENT CONSTANTS  *** */

/*     ***  K = 1 MEANS RETURN STANDARD OUTPUT UNIT NUMBER.   *** */
/*     ***  K = 2 MEANS RETURN ALTERNATE OUTPUT UNIT NUMBER.  *** */
/*     ***  K = 3 MEANS RETURN  INPUT UNIT NUMBER.            *** */
/*          (NOTE -- K = 2, 3 ARE USED ONLY BY TEST PROGRAMS.) */


    ret_val = mdcon[(0 + (0 + (*k - 1 << 2))) / 4];
    return ret_val;
/*  ***  LAST CARD OF IMDCON FOLLOWS  *** */
} /* imdcon_ */

doublereal rmdcon_(integer *k)
{
    /* Initialized data */

    static doublereal one001 = 1.001;
    static doublereal pt999 = .999;
    static doublereal big = 1.79769e308;
    static doublereal eta = 9.999e-307;
    static doublereal machep = 2.22e-16;

    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double sqrt(doublereal);


/*  ***  RETURN MACHINE DEPENDENT CONSTANTS USED BY NL2SOL  *** */

/* +++  COMMENTS BELOW CONTAIN DATA STATEMENTS FOR VARIOUS MACHINES.  +++ */
/* +++  TO CONVERT TO ANOTHER MACHINE, PLACE A C IN COLUMN 1 OF THE   +++ */
/* +++  DATA STATEMENT LINE(S) THAT CORRESPOND TO THE CURRENT MACHINE +++ */
/* +++  AND REMOVE THE C FROM COLUMN 1 OF THE DATA STATEMENT LINE(S)  +++ */
/* +++  THAT CORRESPOND TO THE NEW MACHINE.                           +++ */


/*  ***  THE CONSTANT RETURNED DEPENDS ON K... */

/*  ***        K = 1... SMALLEST POS. ETA SUCH THAT -ETA EXISTS. */
/*  ***        K = 2... SQUARE ROOT OF 1.001*ETA. */
/*  ***        K = 3... UNIT ROUNDOFF = SMALLEST POS. NO. MACHEP SUCH */
/*  ***                 THAT 1 + MACHEP .GT. 1 .AND. 1 - MACHEP .LT. 1. */
/*  ***        K = 4... SQUARE ROOT OF 0.999*MACHEP. */
/*  ***        K = 5... SQUARE ROOT OF 0.999*BIG (SEE K = 6). */
/*  ***        K = 6... LARGEST MACHINE NO. BIG SUCH THAT -BIG EXISTS. */

/* /+ */
/* / */

/*  +++  IEEE : INTEL 80X86, TITAN note ETA's real value is 2.2D-308 +++ */

/*  +++  IBM 360, IBM 370, OR XEROX  +++ */

/*     DATA BIG/Z7FFFFFFFFFFFFFFF/, ETA/Z0010000000000000/, */
/*    1     MACHEP/Z3410000000000000/ */

/*  +++  DATA GENERAL  +++ */

/*     DATA BIG/0.7237005577D+76/, ETA/0.5397605347D-78/, */
/*    1     MACHEP/2.22044605D-16/ */

/*  +++  DEC 11  +++ */

/*     DATA BIG/1.7D+38/, ETA/2.938735878D-39/, MACHEP/2.775557562D-17/ */

/*  +++  HP3000  +++ */

/*     DATA BIG/1.157920892D+77/, ETA/8.636168556D-78/, */
/*    1     MACHEP/5.551115124D-17/ */

/*  +++  HONEYWELL  +++ */

/*     DATA BIG/1.69D+38/, ETA/5.9D-39/, MACHEP/2.1680435D-19/ */

/*  +++  DEC10  +++ */

/*     DATA BIG/"377777100000000000000000/, */
/*    1     ETA/"002400400000000000000000/, */
/*    2     MACHEP/"104400000000000000000000/ */

/*  +++  BURROUGHS  +++ */

/*     DATA BIG/O0777777777777777,O7777777777777777/, */
/*    1     ETA/O1771000000000000,O7770000000000000/, */
/*    2     MACHEP/O1451000000000000,O0000000000000000/ */

/*  +++  CONTROL DATA  +++ */


/*     DATA BIG/37767777777777777777B,37167777777777777777B/, */
/*    1     ETA/00014000000000000000B,00000000000000000000B/, */
/*    2     MACHEP/15614000000000000000B,15010000000000000000B/ */

/*  +++  PRIME  +++ */

/*     DATA BIG/1.0D+9786/, ETA/1.0D-9860/, MACHEP/1.4210855D-14/ */

/*  +++  UNIVAC  +++ */

/*     DATA BIG/8.988D+307/, ETA/1.2D-308/, MACHEP/1.734723476D-18/ */

/*  +++  VAX  +++ */

/*     DATA BIG/1.7D+38/, ETA/2.939D-39/, MACHEP/1.3877788D-17/ */

/* -------------------------------  BODY  -------------------------------- */

    switch (*k) {
	case 1:  goto L10;
	case 2:  goto L20;
	case 3:  goto L30;
	case 4:  goto L40;
	case 5:  goto L50;
	case 6:  goto L60;
    }

L10:
    ret_val = eta;
    goto L999;

L20:
    ret_val = sqrt(one001 * eta);
    goto L999;

L30:
    ret_val = machep;
    goto L999;

L40:
    ret_val = sqrt(pt999 * machep);
    goto L999;

L50:
    ret_val = sqrt(pt999 * big);
    goto L999;

L60:
    ret_val = big;

L999:
    return ret_val;
/*  ***  LAST CARD OF RMDCON FOLLOWS  *** */
} /* rmdcon_ */


/* Subroutine */ int nl2sol_(integer *n, integer *p, doublereal *x, S_fp 
	calcr, S_fp calcj, integer *iv, doublereal *v, integer *uiparm, 
	doublereal *urparm, U_fp ufparm) //P and ntree removed
{	
	if (COMMENT) {
	printf("\nnl2sol_ : entered nl2sol");
	//printf("\nnl2sol_ : P = %i \nnl2sol_ : ntree = %i", P, ntree);
	}    
	/* System generated locals */
    integer i__1;

    /* Local variables */
    static integer d1, j1, r1, nf;
    extern /* Subroutine */ int nl2itr_(doublereal *, integer *, doublereal *,
	     integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *);
    static logical strted;
    extern /* Subroutine */ int itsmry_(doublereal *, integer *, integer *, 
	    doublereal *, doublereal *);


/*  ***  MINIMIZE NONLINEAR SUM OF SQUARES USING ANALYTIC JACOBIAN  *** */
/*  ***  (NL2SOL VERSION 2.2)                                       *** */
/*  ***                                                             *** */
/*  ***  THIS SUBROUTINE WAS SLIGHTLY MODIFIED IN ORDER BE CON-     *** */
/*  ***  SISTENS WITH THE OTHER SUBROUTINES IN THIS PACKAGE         *** */

/*     DIMENSION IV(60+P),  V(93 + N*P + 3*N + P*(3*P+33)/2) */
/*     DIMENSION UIPARM(*), URPARM(*) */

/*  ***  PURPOSE  *** */

/*        GIVEN A P-VECTOR X OF PARAMETERS, CALCR COMPUTES AN N-VECTOR */
/*     R = R(X) OF RESIDUALS CORRESPONDING TO X.  (R(X) PROBABLY ARISES */
/*     FROM A NONLINEAR MODEL INVOLVING P PARAMETERS AND N OBSERVATIONS.) */
/*     THIS ROUTINE INTERACTS WITH NL2ITR TO SEEK A PARAMETER VECTOR X */
/*     THAT MINIMIZES THE SUM OF THE SQUARES OF (THE COMPONENTS OF) R(X), */
/*     I.E., THAT MINIMIZES THE SUM-OF-SQUARES FUNCTION */
/*     F(X) = (R(X)**T) * R(X) / 2.  R(X) IS ASSUMED TO BE A TWICE CON- */
/*     TINUOUSLY DIFFERENTIABLE FUNCTION OF X. */

/* --------------------------  PARAMETER USAGE  -------------------------- */

/* N........ (INPUT) THE NUMBER OF OBSERVATIONS, I.E., THE NUMBER OF */
/*                  COMPONENTS IN R(X).  N MUST BE .GE. P. */
/* P........ (INPUT) THE NUMBER OF PARAMETERS (COMPONENTS IN X).  P MUST */
/*                  BE POSITIVE. */
/* X........ (INPUT/OUTPUT).  ON INPUT, X IS AN INITIAL GUESS AT THE */
/*                  DESIRED PARAMETER ESTIMATE.  ON OUTPUT, X CONTAINS */
/*                  THE BEST PARAMETER ESTIMATE FOUND. */
/* CALCR.... (INPUT) A SUBROUTINE WHICH, GIVEN X, COMPUTES R(X).  CALCR */
/*                  MUST BE DECLARED EXTERNAL IN THE CALLING PROGRAM. */
/*                  IT IS INVOKED BY */
/*                       CALL CALCR(N,P,X,NF,R,UIPARM,URPARM,UFPARM) */
/*                  WHEN CALCR IS CALLED, NF IS THE INVOCATION COUNT */
/*                  FOR CALCR.  IT IS INCLUDED FOR POSSIBLE USE WITH */
/*                  CALCJ.  IF X IS OUT OF BOUNDS (E.G. IF IT WOULD */
/*                  CAUSE OVERFLOW IN COMPUTING R(X)), THEN CALCR SHOULD */
/*                  SET NF TO 0.  THIS WILL CAUSE A SHORTER STEP TO BE */
/*                  ATTEMPTED.  THE OTHER PARAMETERS ARE AS DESCRIBED */
/*                  ABOVE AND BELOW.  CALCR SHOULD NOT CHANGE N, P, OR X. */
/* CALCJ.... (INPUT) A SUBROUTINE WHICH, GIVEN X, COMPUTES THE JACOBIAN */
/*                  MATRIX J OF R AT X, I.E., THE N BY P MATRIX WHOSE */
/*                  (I,K) ENTRY IS THE PARTIAL DERIVATIVE OF THE I-TH */
/*                  COMPONENT OF R WITH RESPECT TO X(K).  CALCJ MUST BE */
/*                  DECLARED EXTERNAL IN THE CALLING PROGRAM.  IT IS */
/*                  INVOKED BY */
/*                       CALL CALCJ(N,P,X,NF,J,UIPARM,URPARM,UFPARM) */
/*                  NF IS THE INVOCATION COUNT FOR CALCR AT THE TIME */
/*                  R(X) WAS EVALUATED.  THE X PASSED TO CALCJ IS */
/*                  USUALLY THE ONE PASSED TO CALCR ON EITHER ITS MOST */
/*                  RECENT INVOCATION OR THE ONE PRIOR TO IT.  IF CALCR */
/*                  SAVES INTERMEDIATE RESULTS FOR USE BY CALCJ, THEN IT */
/*                  IS POSSIBLE TO TELL FROM NF WHETHER THEY ARE VALID */
/*                  FOR THE CURRENT X (OR WHICH COPY IS VALID IF TWO */
/*                  COPIES ARE KEPT).  IF J CANNOT BE COMPUTED AT X, */
/*                  THEN CALCJ SHOULD SET NF TO 0.  IN THIS CASE, NL2SOL */
/*                  WILL RETURN WITH IV(1) = 15.  THE OTHER PARAMETERS */
/*                  TO CALCJ ARE AS DESCRIBED ABOVE AND BELOW.  CALCJ */
/*                  SHOULD NOT CHANGE N, P, OR X. */
/* IV....... (INPUT/OUTPUT) AN INTEGER VALUE ARRAY OF LENGTH AT LEAST */
/*                  60 + P THAT HELPS CONTROL THE NL2SOL ALGORITHM AND */
/*                  THAT IS USED TO STORE VARIOUS INTERMEDIATE QUANTI- */
/*                  TIES.  OF PARTICULAR INTEREST ARE THE INITIALIZATION/ */
/*                  RETURN CODE IV(1) AND THE ENTRIES IN IV THAT CONTROL */
/*                  PRINTING AND LIMIT THE NUMBER OF ITERATIONS AND FUNC- */
/*                  TION EVALUATIONS.  SEE THE SECTION ON IV INPUT */
/*                  VALUES BELOW. */
/* V........ (INPUT/OUTPUT) A FLOATING-POINT VALUE ARRAY OF LENGTH AT */
/*                  LEAST 93 + N*P + 3*N + P*(3*P+33)/2 THAT HELPS CON- */
/*                  TROL THE NL2SOL ALGORITHM AND THAT IS USED TO STORE */
/*                  VARIOUS INTERMEDIATE QUANTITIES.  OF PARTICULAR IN- */
/*                  TEREST ARE THE ENTRIES IN V THAT LIMIT THE LENGTH OF */
/*                  THE FIRST STEP ATTEMPTED (LMAX0), SPECIFY CONVER- */
/*                  GENCE TOLERANCES (AFCTOL, RFCTOL, XCTOL, XFTOL), */
/*                  AND HELP CHOOSE THE STEP SIZE USED IN COMPUTING THE */
/*                  COVARIANCE MATRIX (DELTA0).  SEE THE SECTION ON */
/*                  (SELECTED) V INPUT VALUES BELOW. */
/* UIPARM... (INPUT) USER INTEGER PARAMETER ARRAY PASSED WITHOUT CHANGE */
/*                  TO CALCR AND CALCJ. */
/* URPARM... (INPUT) USER FLOATING-POINT PARAMETER ARRAY PASSED WITHOUT */
/*                  CHANGE TO CALCR AND CALCJ. */
/* UFPARM... (INPUT) USER EXTERNAL SUBROUTINE OR FUNCTION PASSED WITHOUT */
/*                  CHANGE TO CALCR AND CALCJ. */

/*  ***  IV INPUT VALUES (FROM SUBROUTINE DFAULT)  *** */

/* IV(1)...  ON INPUT, IV(1) SHOULD HAVE A VALUE BETWEEN 0 AND 12...... */
/*             0 AND 12 MEAN THIS IS A FRESH START.  0 MEANS THAT */
/*             DFAULT(IV, V) IS TO BE CALLED TO PROVIDE ALL DEFAULT */
/*             VALUES TO IV AND V.  12 (THE VALUE THAT DFAULT ASSIGNS TO */
/*             IV(1)) MEANS THE CALLER HAS ALREADY CALLED DFAULT(IV, V) */
/*             AND HAS POSSIBLY CHANGED SOME IV AND/OR V ENTRIES TO NON- */
/*             DEFAULT VALUES.  DEFAULT = 12. */
/* IV(COVPRT)... IV(14) = 1 MEANS PRINT A COVARIANCE MATRIX AT THE SOLU- */
/*             TION.  (THIS MATRIX IS COMPUTED JUST BEFORE A RETURN WITH */
/*             IV(1) = 3, 4, 5, 6.) */
/*             IV(COVPRT) = 0 MEANS SKIP THIS PRINTING.  DEFAULT = 1. */
/* IV(COVREQ)... IV(15) = NONZERO MEANS COMPUTE A COVARIANCE MATRIX */
/*             JUST BEFORE A RETURN WITH IV(1) = 3, 4, 5, 6.  IN */
/*             THIS CASE, AN APPROXIMATE COVARIANCE MATRIX IS OBTAINED */
/*             IN ONE OF SEVERAL WAYS.  LET K = ABS(IV(COVREQ)) AND LET */
/*             SCALE = 2*F(X)/MAX(1,N-P),  WHERE 2*F(X) IS THE RESIDUAL */
/*             SUM OF SQUARES.  IF K = 1 OR 2, THEN A FINITE-DIFFERENCE */
/*             HESSIAN APPROXIMATION H IS OBTAINED.  IF H IS POSITIVE */
/*             DEFINITE (OR, FOR K = 3, IF THE JACOBIAN MATRIX J AT X */
/*             IS NONSINGULAR), THEN ONE OF THE FOLLOWING IS COMPUTED... */
/*                  K = 1....  SCALE * H**-1 * (J**T * J) * H**-1. */
/*                  K = 2....  SCALE * H**-1. */
/*                  K = 3....  SCALE * (J**T * J)**-1. */
/*             (J**T IS THE TRANSPOSE OF J, WHILE **-1 MEANS INVERSE.) */
/*             IF IV(COVREQ) IS POSITIVE, THEN BOTH FUNCTION AND GRAD- */
/*             IENT VALUES (CALLS ON CALCR AND CALCJ) ARE USED IN COM- */
/*             PUTING H (WITH STEP SIZES DETERMINED USING V(DELTA0) -- */
/*             SEE BELOW), WHILE IF IV(COVREQ) IS NEGATIVE, THEN ONLY */
/*             FUNCTION VALUES (CALLS ON CALCR) ARE USED (WITH STEP */
/*             SIZES DETERMINED USING V(DLTFDC) -- SEE BELOW).  IF */
/*             IV(COVREQ) = 0, THEN NO ATTEMPT IS MADE TO COMPUTE A CO- */
/*             VARIANCE MATRIX (UNLESS IV(COVPRT) = 1, IN WHICH CASE */
/*             IV(COVREQ) = 1 IS ASSUMED).  SEE IV(COVMAT) BELOW. */
/*             DEFAULT = 1. */
/* IV(DTYPE).... IV(16) TELLS HOW THE SCALE VECTOR D (SEE REF. 1) SHOULD */
/*             BE CHOSEN.  IV(DTYPE) .GE. 1 MEANS CHOOSE D AS DESCRIBED */
/*             BELOW WITH V(DFAC).  IV(DTYPE) .LE. 0 MEANS THE CALLER */
/*             HAS CHOSEN D AND HAS STORED IT IN V STARTING AT */
/*             V(94 + 2*N + P*(3*P + 31)/2).  DEFAULT = 1. */
/* IV(INITS).... IV(25) TELLS HOW THE S MATRIX (SEE REF. 1) SHOULD BE */
/*             INITIALIZED.  0 MEANS INITIALIZE S TO 0 (AND START WITH */
/*             THE GAUSS-NEWTON MODEL).  1 AND 2 MEAN THAT THE CALLER */
/*             HAS STORED THE LOWER TRIANGLE OF THE INITIAL S ROWWISE IN */
/*             V STARTING AT V(87+2*P).  IV(INITS) = 1 MEANS START WITH */
/*             THE GAUSS-NEWTON MODEL, WHILE IV(INITS) = 2 MEANS START */
/*             WITH THE AUGMENTED MODEL (SEE REF. 1).  DEFAULT = 0. */
/* IV(MXFCAL)... IV(17) GIVES THE MAXIMUM NUMBER OF FUNCTION EVALUATIONS */
/*             (CALLS ON CALCR, EXCLUDING THOSE USED TO COMPUTE THE CO- */
/*             VARIANCE MATRIX) ALLOWED.  IF THIS NUMBER DOES NOT SUF- */
/*             FICE, THEN NL2SOL RETURNS WITH IV(1) = 9.  DEFAULT = 200. */
/* IV(MXITER)... IV(18) GIVES THE MAXIMUM NUMBER OF ITERATIONS ALLOWED. */
/*             IT ALSO INDIRECTLY LIMITS THE NUMBER OF GRADIENT EVALUA- */
/*             TIONS (CALLS ON CALCJ, EXCLUDING THOSE USED TO COMPUTE */
/*             THE COVARIANCE MATRIX) TO IV(MXITER) + 1.  IF IV(MXITER) */
/*             ITERATIONS DO NOT SUFFICE, THEN NL2SOL RETURNS WITH */
/*             IV(1) = 10.  DEFAULT = 150. */
/* IV(OUTLEV)... IV(19) CONTROLS THE NUMBER AND LENGTH OF ITERATION SUM- */
/*             MARY LINES PRINTED (BY ITSMRY).  IV(OUTLEV) = 0 MEANS DO */
/*             NOT PRINT ANY SUMMARY LINES.  OTHERWISE, PRINT A SUMMARY */
/*             LINE AFTER EACH ABS(IV(OUTLEV)) ITERATIONS.  IF IV(OUTLEV) */
/*             IS POSITIVE, THEN SUMMARY LINES OF LENGTH 117 (PLUS CARRI- */
/*             AGE CONTROL) ARE PRINTED, INCLUDING THE FOLLOWING...  THE */
/*             ITERATION AND FUNCTION EVALUATION COUNTS, CURRENT FUNC- */
/*             TION VALUE (V(F) = HALF THE SUM OF SQUARES), RELATIVE */
/*             DIFFERENCE IN FUNCTION VALUES ACHIEVED BY THE LATEST STEP */
/*             (I.E., RELDF = (F0-V(F))/F0, WHERE F0 IS THE FUNCTION */
/*             VALUE FROM THE PREVIOUS ITERATION), THE RELATIVE FUNCTION */
/*             REDUCTION PREDICTED FOR THE STEP JUST TAKEN (I.E., */
/*             PRELDF = V(PREDUC) / F0, WHERE V(PREDUC) IS DESCRIBED */
/*             BELOW), THE SCALED RELATIVE CHANGE IN X (SEE V(RELDX) */
/*             BELOW), THE MODELS USED IN THE CURRENT ITERATION (G = */
/*             GAUSS-NEWTON, S=AUGMENTED), THE MARQUARDT PARAMETER */
/*             STPPAR USED IN COMPUTING THE LAST STEP, THE SIZING FACTOR */
/*             USED IN UPDATING S, THE 2-NORM OF THE SCALE VECTOR D */
/*             TIMES THE STEP JUST TAKEN (SEE REF. 1), AND NPRELDF, I.E., */
/*             V(NREDUC)/F0, WHERE V(NREDUC) IS DESCRIBED BELOW -- IF */
/*             NPRELDF IS POSITIVE, THEN IT IS THE RELATIVE FUNCTION */
/*             REDUCTION PREDICTED FOR A NEWTON STEP (ONE WITH */
/*             STPPAR = 0).  IF NPRELDF IS ZERO, EITHER THE GRADIENT */
/*             VANISHES (AS DOES PRELDF) OR ELSE THE AUGMENTED MODEL */
/*             IS BEING USED AND ITS HESSIAN IS INDEFINITE (WITH PRELDF */
/*             POSITIVE).  IF NPRELDF IS NEGATIVE, THEN IT IS THE NEGA- */
/*             OF THE RELATIVE FUNCTION REDUCTION PREDICTED FOR A STEP */
/*             COMPUTED WITH STEP BOUND V(LMAX0) FOR USE IN TESTING FOR */
/*             SINGULAR CONVERGENCE. */
/*                  IF IV(OUTLEV) IS NEGATIVE, THEN LINES OF MAXIMUM */
/*             LENGTH 79 (OR 55 IS IV(COVPRT) = 0) ARE PRINTED, INCLUD- */
/*             ING ONLY THE FIRST 6 ITEMS LISTED ABOVE (THROUGH RELDX). */
/*             DEFAULT = 1. */
/* IV(PARPRT)... IV(20) = 1 MEANS PRINT ANY NONDEFAULT V VALUES ON A */
/*             FRESH START OR ANY CHANGED V VALUES ON A RESTART. */
/*             IV(PARPRT) = 0 MEANS SKIP THIS PRINTING.  DEFAULT = 1. */
/* IV(PRUNIT)... IV(21) IS THE OUTPUT UNIT NUMBER ON WHICH ALL PRINTING */
/*             IS DONE.  IV(PRUNIT) = 0 MEANS SUPPRESS ALL PRINTING. */
/*             (SETTING IV(PRUNIT) TO 0 IS THE ONLY WAY TO SUPPRESS THE */
/*             ONE-LINE TERMINATION REASON MESSAGE PRINTED BY ITSMRY.) */
/*             DEFAULT = STANDARD OUTPUT UNIT (UNIT 6 ON MOST SYSTEMS). */
/* IV(SOLPRT)... IV(22) = 1 MEANS PRINT OUT THE VALUE OF X RETURNED (AS */
/*             WELL AS THE CORRESPONDING GRADIENT AND SCALE VECTOR D). */
/*             IV(SOLPRT) = 0 MEANS SKIP THIS PRINTING.  DEFAULT = 1. */
/* IV(STATPR)... IV(23) = 1 MEANS PRINT SUMMARY STATISTICS UPON RETURN- */
/*             ING.  THESE CONSIST OF THE FUNCTION VALUE (HALF THE SUM */
/*             OF SQUARES) AT X, V(RELDX) (SEE BELOW), THE NUMBER OF */
/*             FUNCTION AND GRADIENT EVALUATIONS (CALLS ON CALCR AND */
/*             CALCJ RESPECTIVELY, EXCLUDING ANY CALLS USED TO COMPUTE */
/*             THE COVARIANCE), THE RELATIVE FUNCTION REDUCTIONS PREDICT- */
/*             ED FOR THE LAST STEP TAKEN AND FOR A NEWTON STEP (OR PER- */
/*             HAPS A STEP BOUNDED BY V(LMAX0) -- SEE THE DESCRIPTIONS */
/*             OF PRELDF AND NPRELDF UNDER IV(OUTLEV) ABOVE), AND (IF AN */
/*             ATTEMPT WAS MADE TO COMPUTE THE COVARIANCE) THE NUMBER OF */
/*             CALLS ON CALCR AND CALCJ USED IN TRYING TO COMPUTE THE */
/*             COVARIANCE.  IV(STATPR) = 0 MEANS SKIP THIS PRINTING. */
/*             DEFAULT = 1. */
/* IV(X0PRT).... IV(24) = 1 MEANS PRINT THE INITIAL X AND SCALE VECTOR D */
/*             (ON A FRESH START ONLY).  IV(X0PRT) = 0 MEANS SKIP THIS */
/*             PRINTING.  DEFAULT = 1. */

/*  ***  (SELECTED) IV OUTPUT VALUES  *** */
//legenda
/* IV(1)........ ON OUTPUT, IV(1) IS A RETURN CODE.... */
/*             3 = X-CONVERGENCE.  THE SCALED RELATIVE DIFFERENCE BE- */
/*                  TWEEN THE CURRENT PARAMETER VECTOR X AND A LOCALLY */
/*                  OPTIMAL PARAMETER VECTOR IS VERY LIKELY AT MOST */
/*                  V(XCTOL). */
/*             4 = RELATIVE FUNCTION CONVERGENCE.  THE RELATIVE DIFFER- */
/*                  ENCE BETWEEN THE CURRENT FUNCTION VALUE AND ITS LO- */
/*                  CALLY OPTIMAL VALUE IS VERY LIKELY AT MOST V(RFCTOL). */
/*             5 = BOTH X- AND RELATIVE FUNCTION CONVERGENCE (I.E., THE */
/*                  CONDITIONS FOR IV(1) = 3 AND IV(1) = 4 BOTH HOLD). */
/*             6 = ABSOLUTE FUNCTION CONVERGENCE.  THE CURRENT FUNCTION */
/*                  VALUE IS AT MOST V(AFCTOL) IN ABSOLUTE VALUE. */
/*             7 = SINGULAR CONVERGENCE.  THE HESSIAN NEAR THE CURRENT */
/*                  ITERATE APPEARS TO BE SINGULAR OR NEARLY SO, AND A */
/*                  STEP OF LENGTH AT MOST V(LMAX0) IS UNLIKELY TO YIELD */
/*                  A RELATIVE FUNCTION DECREASE OF MORE THAN V(RFCTOL). */
/*             8 = FALSE CONVERGENCE.  THE ITERATES APPEAR TO BE CONVERG- */
/*                  ING TO A NONCRITICAL POINT.  THIS MAY MEAN THAT THE */
/*                  CONVERGENCE TOLERANCES (V(AFCTOL), V(RFCTOL), */
/*                  V(XCTOL)) ARE TOO SMALL FOR THE ACCURACY TO WHICH */
/*                  THE FUNCTION AND GRADIENT ARE BEING COMPUTED, THAT */
/*                  THERE IS AN ERROR IN COMPUTING THE GRADIENT, OR THAT */
/*                  THE FUNCTION OR GRADIENT IS DISCONTINUOUS NEAR X. */
/*             9 = FUNCTION EVALUATION LIMIT REACHED WITHOUT OTHER CON- */
/*                  VERGENCE (SEE IV(MXFCAL)). */
/*            10 = ITERATION LIMIT REACHED WITHOUT OTHER CONVERGENCE */
/*                  (SEE IV(MXITER)). */
/*            11 = STOPX RETURNED .TRUE. (EXTERNAL INTERRUPT).  SEE THE */
/*                  USAGE NOTES BELOW. */
/*            13 = F(X) CANNOT BE COMPUTED AT THE INITIAL X. */
/*            14 = BAD PARAMETERS PASSED TO ASSESS (WHICH SHOULD NOT */
/*                  OCCUR). */
/*            15 = THE JACOBIAN COULD NOT BE COMPUTED AT X (SEE CALCJ */
/*                  ABOVE). */
/*            16 = N OR P (OR PARAMETER NN TO NL2ITR) OUT OF RANGE -- */
/*                  P .LE. 0 OR N .LT. P OR NN .LT. N. */
/*            17 = RESTART ATTEMPTED WITH N OR P (OR PAR. NN TO NL2ITR) */
/*                  CHANGED. */
/*            18 = IV(INITS) IS OUT OF RANGE. */
/*            19...45 = V(IV(1)) IS OUT OF RANGE. */
/*            50 = IV(1) WAS OUT OF RANGE. */
/*            87...(86+P) = JTOL(IV(1)-86) (I.E., V(IV(1)) IS NOT */
/*                  POSITIVE (SEE V(DFAC) BELOW). */
/* IV(COVMAT)... IV(26) TELLS WHETHER A COVARIANCE MATRIX WAS COMPUTED. */
/*             IF (IV(COVMAT) IS POSITIVE, THEN THE LOWER TRIANGLE OF */
/*             THE COVARIANCE MATRIX IS STORED ROWWISE IN V STARTING AT */
/*             V(IV(COVMAT)).  IF IV(COVMAT) = 0, THEN NO ATTEMPT WAS */
/*             MADE TO COMPUTE THE COVARIANCE.  IF IV(COVMAT) = -1, */
/*             THEN THE FINITE-DIFFERENCE HESSIAN WAS INDEFINITE.  AND */
/*             AND IF IV(COVMAT) = -2, THEN A SUCCESSFUL FINITE-DIFFER- */
/*             ENCING STEP COULD NOT BE FOUND FOR SOME COMPONENT OF X */
/*             (I.E., CALCR SET NF TO 0 FOR EACH OF TWO TRIAL STEPS). */
/*             NOTE THAT IV(COVMAT) IS RESET TO 0 AFTER EACH SUCCESSFUL */
/*             STEP, SO IF SUCH A STEP IS TAKEN AFTER A RESTART, THEN */
/*             THE COVARIANCE MATRIX WILL BE RECOMPUTED. */
/* IV(D)........ IV(27) IS THE STARTING SUBSCRIPT IN V OF THE CURRENT */
/*             SCALE VECTOR D. */
/* IV(G)........ IV(28) IS THE STARTING SUBSCRIPT IN V OF THE CURRENT */
/*             LEAST-SQUARES GRADIENT VECTOR (J**T)*R. */
/* IV(NFCALL)... IV(6) IS THE NUMBER OF CALLS SO FAR MADE ON CALCR (I.E., */
/*             FUNCTION EVALUATIONS, INCLUDING THOSE USED IN COMPUTING */
/*             THE COVARIANCE). */
/* IV(NFCOV).... IV(40) IS THE NUMBER OF CALLS MADE ON CALCR WHEN */
/*             TRYING TO COMPUTE COVARIANCE MATRICES. */
/* IV(NGCALL)... IV(30) IS THE NUMBER OF GRADIENT EVALUATIONS (CALLS ON */
/*             CALCJ) SO FAR DONE (INCLUDING THOSE USED FOR COMPUTING */
/*             THE COVARIANCE). */
/* IV(NGCOV).... IV(41) IS THE NUMBER OF CALLS MADE ON CALCJ WHEN */
/*             TRYING TO COMPUTE COVARIANCE MATRICES. */
/* IV(NITER).... IV(31) IS THE NUMBER OF ITERATIONS PERFORMED. */
/* IV(R)........ IV(50) IS THE STARTING SUBSCRIPT IN V OF THE RESIDUAL */
/*             VECTOR R CORRESPONDING TO X. */

/*  ***  (SELECTED) V INPUT VALUES (FROM SUBROUTINE DFAULT)  *** */

/* V(AFCTOL)... V(31) IS THE ABSOLUTE FUNCTION CONVERGENCE TOLERANCE. */
/*             IF NL2SOL FINDS A POINT WHERE THE FUNCTION VALUE (HALF */
/*             THE SUM OF SQUARES) IS LESS THAN V(AFCTOL), AND IF NL2SOL */
/*             DOES NOT RETURN WITH IV(1) = 3, 4, OR 5, THEN IT RETURNS */
/*             WITH IV(1) = 6.  DEFAULT = MAX(10**-20, MACHEP**2), WHERE */
/*             MACHEP IS THE UNIT ROUNDOFF. */
/* V(DELTA0)... V(44) IS A FACTOR USED IN CHOOSING THE FINITE-DIFFERENCE */
/*             STEP SIZE USED IN COMPUTING THE COVARIANCE MATRIX WHEN */
/*             IV(COVREQ) = 1 OR 2.  FOR COMPONENT I, STEP SIZE */
/*                  V(DELTA0) * MAX(ABS(X(I)), 1/D(I)) * SIGN(X(I)) */
/*             IS USED, WHERE D IS THE CURRENT SCALE VECTOR (SEE REF. 1). */
/*             (IF THIS STEP RESULTS IN CALCR SETTING NF TO 0, THEN -0.5 */
/*             TIMES THIS STEP IS ALSO TRIED.)  DEFAULT = MACHEP**0.5, */
/*             WHERE MACHEP IS THE UNIT ROUNDOFF. */
/* V(DFAC)..... V(41) AND THE D0 AND JTOL ARRAYS (SEE V(D0INIT) AND */
/*             V(JTINIT)) ARE USED IN UPDATING THE SCALE VECTOR D WHEN */
/*             IV(DTYPE) .GT. 0.  (D IS INITIALIZED ACCORDING TO */
/*             V(DINIT).)  LET D1(I) = */
/*               MAX(SQRT(JCNORM(I)**2 + MAX(S(I,I),0)), V(DFAC)*D(I)), */
/*             WHERE JCNORM(I) IS THE 2-NORM OF THE I-TH COLUMN OF THE */
/*             CURRENT JACOBIAN MATRIX AND S IS THE S MATRIX OF REF. 1. */
/*             IF IV(DTYPE) = 1, THEN D(I) IS SET TO D1(I) UNLESS */
/*             D1(I) .LT. JTOL(I), IN WHICH CASE D(I) IS SET TO */
/*                                MAX(D0(I), JTOL(I)). */
/*             IF IV(DTYPE) .GE. 2, THEN D IS UPDATED DURING THE FIRST */
/*             ITERATION AS FOR IV(DTYPE) = 1 (AFTER ANY INITIALIZATION */
/*             DUE TO V(DINIT)) AND IS LEFT UNCHANGED THEREAFTER. */
/*             DEFAULT = 0.6. */
/* V(DINIT).... V(38), IF NONNEGATIVE, IS THE VALUE TO WHICH THE SCALE */
/*             VECTOR D IS INITIALIZED.  DEFAULT = 0. */
/* V(DLTFDC)... V(40) HELPS CHOOSE THE STEP SIZE USED WHEN COMPUTING THE */
/*             COVARIANCE MATRIX WHEN IV(COVREQ) = -1 OR -2.  FOR */
/*             DIFFERENCES INVOLVING X(I), THE STEP SIZE FIRST TRIED IS */
/*                       V(DLTFDC) * MAX(ABS(X(I)), 1/D(I)), */
/*             WHERE D IS THE CURRENT SCALE VECTOR (SEE REF. 1).  (IF */
/*             THIS STEP IS TOO BIG THE FIRST TIME IT IS TRIED, I.E., IF */
/*             CALCR SETS NF TO 0, THEN -0.5 TIMES THIS STEP IS ALSO */
/*             TRIED.)  DEFAULT = MACHEP**(1/3), WHERE MACHEP IS THE */
/*             UNIT ROUNDOFF. */
/* V(D0INIT)... V(37), IF POSITIVE, IS THE VALUE TO WHICH ALL COMPONENTS */
/*             OF THE D0 VECTOR (SEE V(DFAC)) ARE INITIALIZED.  IF */
/*             V(DFAC) = 0, THEN IT IS ASSUMED THAT THE CALLER HAS */
/*             STORED D0 IN V STARTING AT V(P+87).  DEFAULT = 1.0. */
/* V(JTINIT)... V(39), IF POSITIVE, IS THE VALUE TO WHICH ALL COMPONENTS */
/*             OF THE JTOL ARRAY (SEE V(DFAC)) ARE INITIALIZED.  IF */
/*             V(JTINIT) = 0, THEN IT IS ASSUMED THAT THE CALLER HAS */
/*             STORED JTOL IN V STARTING AT V(87).  DEFAULT = 10**-6. */
/* V(LMAX0).... V(35) GIVES THE MAXIMUM 2-NORM ALLOWED FOR D TIMES THE */
/*             VERY FIRST STEP THAT NL2SOL ATTEMPTS.  IT IS ALSO USED */
/*             IN TESTING FOR SINGULAR CONVERGENCE -- IF THE FUNCTION */
/*             REDUCTION PREDICTED FOR A STEP OF LENGTH BOUNDED BY */
/*             V(LMAX0) IS AT MOST V(RFCTOL) * ABS(F0), WHERE  F0  IS */
/*             THE FUNCTION VALUE AT THE START OF THE CURRENT ITERATION, */
/*             AND IF NL2SOL DOES NOT RETURN WITH IV(1) = 3, 4, 5, OR 6, */
/*             THEN IT RETURNS WITH IV(1) = 7.    DEFAULT = 100. */
/* V(RFCTOL)... V(32) IS THE RELATIVE FUNCTION CONVERGENCE TOLERANCE. */
/*             IF THE CURRENT MODEL PREDICTS A MAXIMUM POSSIBLE FUNCTION */
/*             REDUCTION (SEE V(NREDUC)) OF AT MOST V(RFCTOL)*ABS(F0) AT */
/*             THE START OF THE CURRENT ITERATION, WHERE  F0  IS THE */
/*             THEN CURRENT FUNCTION VALUE, AND IF THE LAST STEP ATTEMPT- */
/*             ED ACHIEVED NO MORE THAN TWICE THE PREDICTED FUNCTION */
/*             DECREASE, THEN NL2SOL RETURNS WITH IV(1) = 4 (OR 5). */
/*             DEFAULT = MAX(10**-10, MACHEP**(2/3)), WHERE MACHEP IS */
/*             THE UNIT ROUNDOFF. */
/* V(TUNER1)... V(26) HELPS DECIDE WHEN TO CHECK FOR FALSE CONVERGENCE */
/*             AND TO CONSIDER SWITCHING MODELS.  THIS IS DONE IF THE */
/*             ACTUAL FUNCTION DECREASE FROM THE CURRENT STEP IS NO MORE */
/*             THAN V(TUNER1) TIMES ITS PREDICTED VALUE.  DEFAULT = 0.1. */
/* V(XCTOL).... V(33) IS THE X-CONVERGENCE TOLERANCE.  IF A NEWTON STEP */
/*             (SEE V(NREDUC)) IS TRIED THAT HAS V(RELDX) .LE. V(XCTOL) */
/*             AND IF THIS STEP YIELDS AT MOST TWICE THE PREDICTED FUNC- */
/*             TION DECREASE, THEN NL2SOL RETURNS WITH IV(1) = 3 (OR 5). */
/*             (SEE THE DESCRIPTION OF V(RELDX) BELOW.) */
/*             DEFAULT = MACHEP**0.5, WHERE MACHEP IS THE UNIT ROUNDOFF. */
/* V(XFTOL).... V(34) IS THE FALSE CONVERGENCE TOLERANCE.  IF A STEP IS */
/*             TRIED THAT GIVES NO MORE THAN V(TUNER1) TIMES THE PREDICT- */
/*             ED FUNCTION DECREASE AND THAT HAS V(RELDX) .LE. V(XFTOL), */
/*             AND IF NL2SOL DOES NOT RETURN WITH IV(1) = 3, 4, 5, 6, OR */
/*             7, THEN IT RETURNS WITH IV(1) = 8.  (SEE THE DESCRIPTION */
/*             OF V(RELDX) BELOW.)  DEFAULT = 100*MACHEP, WHERE */
/*             MACHEP IS THE UNIT ROUNDOFF. */
/* V(*)........ DFAULT SUPPLIES TO V A NUMBER OF TUNING CONSTANTS, WITH */
/*             WHICH IT SHOULD ORDINARILY BE UNNECESSARY TO TINKER.  SEE */
/*             VERSION 2.2 OF THE NL2SOL USAGE SUMMARY (WHICH IS AN */
/*             APPENDIX TO REF. 1). */

/*  ***  (SELECTED) V OUTPUT VALUES  *** */

/* V(DGNORM)... V(1) IS THE 2-NORM OF (D**-1)*G, WHERE G IS THE MOST RE- */
/*             CENTLY COMPUTED GRADIENT AND D IS THE CORRESPONDING SCALE */
/*             VECTOR. */
/* V(DSTNRM)... V(2) IS THE 2-NORM OF D*STEP, WHERE STEP IS THE MOST RE- */
/*             CENTLY COMPUTED STEP AND D IS THE CURRENT SCALE VECTOR. */
/* V(F)........ V(10) IS THE CURRENT FUNCTION VALUE (HALF THE SUM OF */
/*             SQUARES). */
/* V(F0)....... V(13) IS THE FUNCTION VALUE AT THE START OF THE CURRENT */
/*             ITERATION. */
/* V(NREDUC)... V(6), IF POSITIVE, IS THE MAXIMUM FUNCTION REDUCTION */
/*             POSSIBLE ACCORDING TO THE CURRENT MODEL, I.E., THE FUNC- */
/*             TION REDUCTION PREDICTED FOR A NEWTON STEP (I.E., */
/*             STEP = -H**-1 * G,  WHERE  G = (J**T) * R  IS THE CURRENT */
/*             GRADIENT AND H IS THE CURRENT HESSIAN APPROXIMATION -- */
/*             H = (J**T)*J  FOR THE GAUSS-NEWTON MODEL AND */
/*             H = (J**T)*J + S  FOR THE AUGMENTED MODEL). */
/*                  V(NREDUC) = ZERO MEANS H IS NOT POSITIVE DEFINITE. */
/*                  IF V(NREDUC) IS NEGATIVE, THEN IT IS THE NEGATIVE OF */
/*             THE FUNCTION REDUCTION PREDICTED FOR A STEP COMPUTED WITH */
/*             A STEP BOUND OF V(LMAX0) FOR USE IN TESTING FOR SINGULAR */
/*             CONVERGENCE. */
/* V(PREDUC)... V(7) IS THE FUNCTION REDUCTION PREDICTED (BY THE CURRENT */
/*             QUADRATIC MODEL) FOR THE CURRENT STEP.  THIS (DIVIDED BY */
/*             V(F0)) IS USED IN TESTING FOR RELATIVE FUNCTION */
/*             CONVERGENCE. */
/* V(RELDX).... V(17) IS THE SCALED RELATIVE CHANGE IN X CAUSED BY THE */
/*             CURRENT STEP, COMPUTED AS */
/*                  MAX(ABS(D(I)*(X(I)-X0(I)), 1 .LE. I .LE. P) / */
/*                     MAX(D(I)*(ABS(X(I))+ABS(X0(I))), 1 .LE. I .LE. P), */
/*             WHERE X = X0 + STEP. */

/* -------------------------------  NOTES  ------------------------------- */

/*  ***  ALGORITHM NOTES  *** */

/*        SEE REF. 1 FOR A DESCRIPTION OF THE ALGORITHM USED. */
/*        ON PROBLEMS WHICH ARE NATURALLY WELL SCALED, BETTER PERFORM- */
/*     ANCE MAY BE OBTAINED BY SETTING V(D0INIT) = 1.0 AND IV(DTYPE) = 0, */
/*     WHICH WILL CAUSE THE SCALE VECTOR D TO BE SET TO ALL ONES. */

/*  ***  USAGE NOTES  *** */

/*        AFTER A RETURN WITH IV(1) .LE. 11, IT IS POSSIBLE TO RESTART, */
/*     I.E., TO CHANGE SOME OF THE IV AND V INPUT VALUES DESCRIBED ABOVE */
/*     AND CONTINUE THE ALGORITHM FROM THE POINT WHERE IT WAS INTERRUPT- */
/*     ED.  IV(1) SHOULD NOT BE CHANGED, NOR SHOULD ANY ENTRIES OF IV */
/*     AND V OTHER THAN THE INPUT VALUES (THOSE SUPPLIED BY DFAULT). */
/*        THOSE WHO DO NOT WISH TO WRITE A CALCJ WHICH COMPUTES THE JA- */
/*     COBIAN MATRIX ANALYTICALLY SHOULD CALL NL2SNO RATHER THAN NL2SOL. */
/*     NL2SNO USES FINITE DIFFERENCES TO COMPUTE AN APPROXIMATE JACOBIAN. */
/*        THOSE WHO WOULD PREFER TO PROVIDE R AND J (THE RESIDUAL AND */
/*     JACOBIAN) BY REVERSE COMMUNICATION RATHER THAN BY WRITING SUBROU- */
/*     TINES CALCR AND CALCJ MAY CALL ON NL2ITR DIRECTLY.  SEE THE COM- */
/*     MENTS AT THE BEGINNING OF NL2ITR. */
/*        THOSE WHO USE NL2SOL INTERACTIVELY MAY WISH TO SUPPLY THEIR */
/*     OWN STOPX FUNCTION, WHICH SHOULD RETURN .TRUE. IF THE BREAK KEY */
/*     HAS BEEN PRESSED SINCE STOPX WAS LAST INVOKED.  THIS MAKES IT POS- */
/*     SIBLE TO EXTERNALLY INTERRUPT NL2SOL (WHICH WILL RETURN WITH */
/*     IV(1) = 11 IF STOPX RETURNS .TRUE.). */
/*        STORAGE FOR J IS ALLOCATED AT THE END OF V.  THUS THE CALLER */
/*     MAY MAKE V LONGER THAN SPECIFIED ABOVE AND MAY ALLOW CALCJ TO USE */
/*     ELEMENTS OF J BEYOND THE FIRST N*P AS SCRATCH STORAGE. */

/*  ***  PORTABILITY NOTES  *** */

/*        THE NL2SOL DISTRIBUTION TAPE CONTAINS BOTH SINGLE- AND DOUBLE- */
/*     PRECISION VERSIONS OF THE NL2SOL SOURCE CODE, SO IT SHOULD BE UN- */
/*     NECESSARY TO CHANGE PRECISIONS. */
/*        ONLY THE FUNCTIONS IMDCON AND RMDCON CONTAIN MACHINE-DEPENDENT */
/*     CONSTANTS.  TO CHANGE FROM ONE MACHINE TO ANOTHER, IT SHOULD */
/*     SUFFICE TO CHANGE THE (FEW) RELEVANT LINES IN THESE FUNCTIONS. */
/*        INTRINSIC FUNCTIONS ARE EXPLICITLY DECLARED.  ON CERTAIN COM- */
/*     PUTERS (E.G. UNIVAC), IT MAY BE NECESSARY TO COMMENT OUT THESE */
/*     DECLARATIONS.  SO THAT THIS MAY BE DONE AUTOMATICALLY BY A SIMPLE */
/*     PROGRAM, SUCH DECLARATIONS ARE PRECEDED BY A COMMENT HAVING C/+ */
/*     IN COLUMNS 1-3 AND BLANKS IN COLUMNS 4-72 AND ARE FOLLOWED BY */
/*     A COMMENT HAVING C/ IN COLUMNS 1 AND 2 AND BLANKS IN COLUMNS 3-72. */
/*        THE NL2SOL SOURCE CODE IS EXPRESSED IN 1966 ANSI STANDARD */
/*     FORTRAN.  IT MAY BE CONVERTED TO FORTRAN 77 BY */
/*     COMMENTING OUT ALL LINES THAT FALL BETWEEN A LINE HAVING C/6 IN */
/*     COLUMNS 1-3 AND A LINE HAVING C/7 IN COLUMNS 1-3 AND BY REMOVING */
/*     (I.E., REPLACING BY A BLANK) THE C IN COLUMN 1 OF THE LINES THAT */
/*     FOLLOW THE C/7 LINE AND PRECEED A LINE HAVING C/ IN COLUMNS 1-2 */
/*     AND BLANKS IN COLUMNS 3-72.  THESE CHANGES CONVERT SOME DATA */
/*     STATEMENTS INTO PARAMETER STATEMENTS, CONVERT SOME VARIABLES FROM */
/*     REAL TO CHARACTER*4, AND MAKE THE DATA STATEMENTS THAT INITIALIZE */
/*     THESE VARIABLES USE CHARACTER STRINGS DELIMITED BY PRIMES INSTEAD */
/*     OF HOLLERITH CONSTANTS.  (SUCH VARIABLES AND DATA STATEMENTS */
/*     APPEAR ONLY IN MODULES ITSMRY AND PARCHK.  PARAMETER STATEMENTS */
/*     APPEAR NEARLY EVERYWHERE.) */

/*  ***  REFERENCES  *** */

/* 1.  DENNIS, J.E., GAY, D.M., AND WELSCH, R.E. (1981), AN ADAPTIVE */
/*             NONLINEAR LEAST-SQUARES ALGORITHM, ACM TRANS. MATH. */
/*             SOFTWARE, VOL. 7, NO. 3. */


/*  ***  GENERAL  *** */

/*     CODED BY DAVID M. GAY (WINTER 1979 - WINTER 1980). */
/*     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH */
/*     SUPPORTED BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS */
/*     MCS-7600324, DCR75-10143, 76-14311DSS, MCS76-11989, AND */
/*     MCS-7906671. */

/* ----------------------------  DECLARATIONS  --------------------------- */

/* ITSMRY... PRINTS ITERATION SUMMARY AND INFO ABOUT INITIAL AND FINAL X. */
/* NL2ITR... REVERSE-COMMUNICATION ROUTINE THAT CARRIES OUT NL2SOL ALGO- */
/*             RITHM. */


/*  ***  SUBSCRIPTS FOR IV AND V  *** */


/*  ***  IV SUBSCRIPT VALUES  *** */

/* /6 */
/*     DATA NFCALL/6/, NFGCAL/7/, TOOBIG/2/ */
/* /7 */
/* / */

/*  ***  V SUBSCRIPT VALUES  *** */

/* /6 */
/*     DATA D/27/, J/33/, R/50/ */
/* /7 */
/* / */

/* +++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++ */

    /* Parameter adjustments */
    --x;
    --iv;
    --v;
    --uiparm;
    --urparm;

    /* Function Body */
    d1 = (*n << 1) + 94 + *p * (*p * 3 + 31) / 2;
    iv[27] = d1;
    r1 = d1 + *p;
    iv[50] = r1;
    j1 = r1 + *n;
    iv[33] = j1;
    strted = TRUE_;
    if (iv[1] != 0 && iv[1] != 12) {
	goto L40;
    }
    strted = FALSE_;
    iv[6] = 1;
    iv[7] = 1;

L10:
    nf = iv[6];
    (*calcr)(n, p, &x[1], &nf, &v[r1], &uiparm[1], &urparm[1], (U_fp)ufparm); //P and ntree removed
    if (strted) {
	goto L20;
    }
    if (nf > 0) {
	goto L30;
    }
    iv[1] = 13;
    goto L60;

L20:
    if (nf <= 0) {
	iv[2] = 1;
    }
    goto L40;

L30:
    (*calcj)(n, p, &x[1], &iv[7], &v[j1], &uiparm[1], &urparm[1], (U_fp)
	    ufparm);   //P and ntree removed
    if (iv[7] == 0) {
	goto L50;
    }
    strted = TRUE_;

L40:
    nl2itr_(&v[d1], &iv[1], &v[j1], n, n, p, &v[r1], &v[1], &x[1]);
    if ((i__1 = iv[1] - 2) < 0) {
	goto L10;
    } else if (i__1 == 0) {
	goto L30;
    } else {
	goto L999;
    }

L50:
    iv[1] = 15;
L60:
    itsmry_(&v[d1], &iv[1], p, &v[1], &x[1]);

L999:
    return 0;
/*  ***  LAST CARD OF NL2SOL FOLLOWS  *** */
} /* nl2sol_ */

/* Subroutine */ int nl2sno_(integer *n, integer *p, doublereal *x, S_fp 
	calcr, integer *iv, doublereal *v, integer *uiparm, doublereal *
	urparm, U_fp ufparm)
{
    /* Initialized data */

    static doublereal hlim = 0.;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal h__;
    static integer i__, k, d1, j1, r1, dk, nf, rn;
    static doublereal xk;
    static integer j1k, dinit;
    extern /* Subroutine */ int nl2itr_(doublereal *, integer *, doublereal *,
	     integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *), dfault_(integer *, doublereal *);
    extern doublereal rmdcon_(integer *);
    static logical strted;
    extern /* Subroutine */ int vscopy_(integer *, doublereal *, doublereal *)
	    , itsmry_(doublereal *, integer *, integer *, doublereal *, 
	    doublereal *);


/*  ***  LIKE NL2SOL, BUT WITHOUT CALCJ -- MINIMIZE NONLINEAR SUM OF  *** */
/*  ***  SQUARES USING FINITE-DIFFERENCE JACOBIAN APPROXIMATIONS      *** */
/*  ***  (NL2SOL VERSION 2.2)  *** */

/*     DIMENSION IV(60+P),  V(93 + N*P + 3*N + P*(3*P+33)/2) */

/* -----------------------------  DISCUSSION  ---------------------------- */

/*        THE PARAMETERS FOR NL2SNO ARE THE SAME AS THOSE FOR NL2SOL */
/*     (WHICH SEE), EXCEPT THAT CALCJ IS OMITTED.  INSTEAD OF CALLING */
/*     CALCJ TO OBTAIN THE JACOBIAN MATRIX OF R AT X, NL2SNO COMPUTES */
/*     AN APPROXIMATION TO IT BY FINITE (FORWARD) DIFFERENCES -- SEE */
/*     V(DLTFDJ) BELOW.  NL2SNO USES FUNCTION VALUES ONLY WHEN COMPUT- */
/*     THE COVARIANCE MATRIX (RATHER THAN THE FUNCTIONS AND GRADIENTS */
/*     THAT NL2SOL MAY USE).  TO DO SO, NL2SNO SETS IV(COVREQ) TO -1 IF */
/*     IV(COVPRT) = 1 WITH IV(COVREQ) = 0 AND TO MINUS ITS ABSOLUTE */
/*     VALUE OTHERWISE.  THUS V(DELTA0) IS NEVER REFERENCED AND ONLY */
/*     V(DLTFDC) MATTERS -- SEE NL2SOL FOR A DESCRIPTION OF V(DLTFDC). */
/*        THE NUMBER OF EXTRA CALLS ON CALCR USED IN COMPUTING THE JACO- */
/*     BIAN APPROXIMATION ARE NOT INCLUDED IN THE FUNCTION EVALUATION */
/*     COUNT IV(NFCALL) AND ARE NOT OTHERWISE REPORTED. */

/* V(DLTFDJ)... V(36) HELPS CHOOSE THE STEP SIZE USED WHEN COMPUTING THE */
/*             FINITE-DIFFERENCE JACOBIAN MATRIX.  FOR DIFFERENCES IN- */
/*             VOLVING X(I), THE STEP SIZE FIRST TRIED IS */
/*                       V(DLTFDJ) * MAX(ABS(X(I)), 1/D(I)), */
/*             WHERE D IS THE CURRENT SCALE VECTOR (SEE REF. 1).  (IF */
/*             THIS STEP IS TOO BIG, I.E., IF CALCR SETS NF TO 0, THEN */
/*             SMALLER STEPS ARE TRIED UNTIL THE STEP SIZE IS SHRUNK BE- */
/*             LOW 1000 * MACHEP, WHERE MACHEP IS THE UNIT ROUNDOFF. */
/*             DEFAULT = MACHEP**0.5. */

/*  ***  REFERENCES  *** */

/* 1.  DENNIS, J.E., GAY, D.M., AND WELSCH, R.E. (1981), AN ADAPTIVE */
/*             NONLINEAR LEAST-SQUARES ALGORITHM, ACM TRANS. MATH. */
/*             SOFTWARE, VOL. 7, NO. 3. */

/*  ***  GENERAL  *** */

/*     CODED BY DAVID M. GAY. */
/*     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH */
/*     SUPPORTED BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS */
/*     MCS-7600324, DCR75-10143, 76-14311DSS, MCS76-11989, AND */
/*     MCS-7906671. */

/* +++++++++++++++++++++++++++  DECLARATIONS  ++++++++++++++++++++++++++++ */

/*  ***  INTRINSIC FUNCTIONS  *** */
/* /+ */
/* / */
/*  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  *** */


/* DFAULT... SUPPLIES DEFAULT PARAMETER VALUES. */
/* ITSMRY... PRINTS ITERATION SUMMARY AND INFO ABOUT INITIAL AND FINAL X. */
/* NL2ITR... REVERSE-COMMUNICATION ROUTINE THAT CARRIES OUT NL2SOL ALGO- */
/*             RITHM. */
/* RMDCON... RETURNS MACHINE-DEPENDENT CONSTANTS. */
/* VSCOPY... SETS ALL ELEMENTS OF A VECTOR TO A SCALAR. */


/*  ***  SUBSCRIPTS FOR IV AND V  *** */


/* /6 */
/*     DATA HFAC/1.D+3/, NEGPT5/-0.5D+0/, ONE/1.D+0/, ZERO/0.D+0/ */
/* /7 */
/* / */

/*  ***  IV SUBSCRIPT VALUES  *** */

/* /6 */
/*     DATA COVPRT/14/, COVREQ/15/, D/27/, DTYPE/16/, J/33/, */
/*    1     NFCALL/6/, NFGCAL/7/, R/50/, TOOBIG/2/ */
/* /7 */
/* / */

/*  ***  V SUBSCRIPT VALUES  *** */

/* /6 */
/*     DATA DLTFDJ/36/, DINIT/38/ */
/* /7 */
/* / */
    /* Parameter adjustments */
    --x;
    --iv;
    --v;
    --uiparm;
    --urparm;

    /* Function Body */

/* +++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++ */

    d1 = 94 + 2 * *n + *p * (3 * *p + 31) / 2;
    iv[27] = d1;
    r1 = d1 + *p;
    iv[50] = r1;
    j1 = r1 + *n;
    iv[33] = j1;
    rn = j1 - 1;
    if (iv[1] == 0) {
	dfault_(&iv[1], &v[1]);
    }
    iv[15] = -abs(iv[15]);
    if (iv[14] != 0 && iv[15] == 0) {
	iv[15] = -1;
    }
    strted = TRUE_;
    if (iv[1] != 12) {
	goto L80;
    }
    strted = FALSE_;
    iv[6] = 1;
    iv[7] = 1;
/*        ***  INITIALIZE SCALE VECTOR D TO ONES FOR COMPUTING */
/*        ***  INITIAL JACOBIAN. */
    if (iv[16] > 0) {
	vscopy_(p, &v[d1], &c_b76);
    }
    if (v[dinit] > 0.) {
	vscopy_(p, &v[d1], &v[dinit]);
    }

L10:
    nf = iv[6];
    (*calcr)(n, p, &x[1], &nf, &v[r1], &uiparm[1], &urparm[1], (U_fp)ufparm);//, P, ntree); //P and ntree added
    if (strted) {
	goto L20;
    }
    if (nf > 0) {
	goto L30;
    }
    iv[1] = 13;
    goto L90;

L20:
    if (nf <= 0) {
	iv[2] = 1;
    }
    goto L80;

/*  ***  COMPUTE FINITE-DIFFERENCE JACOBIAN  *** */

L30:
    j1k = j1;
    dk = d1;
    i__1 = *p;
    for (k = 1; k <= i__1; ++k) {
	xk = x[k];
/* Computing MAX */
	d__1 = abs(xk), d__2 = 1. / v[dk];
	h__ = v[36] * max(d__1,d__2);
	++dk;
L40:
	x[k] = xk + h__;
	nf = iv[7];
	(*calcr)(n, p, &x[1], &nf, &v[j1k], &uiparm[1], &urparm[1], (U_fp)
		ufparm);//, P, ntree);  // P and ntree added
	if (nf > 0) {
	    goto L50;
	}
	if (hlim == 0.) {
	    hlim = rmdcon_(&c__3) * 1e3;
	}
/*             ***  HLIM = HFAC TIMES THE UNIT ROUNDOFF  *** */
	h__ *= -.5;
	if (abs(h__) >= hlim) {
	    goto L40;
	}
	iv[1] = 15;
	goto L90;
L50:
	x[k] = xk;
	i__2 = rn;
	for (i__ = r1; i__ <= i__2; ++i__) {
	    v[j1k] = (v[j1k] - v[i__]) / h__;
	    ++j1k;
/* L60: */
	}
/* L70: */
    }

    strted = TRUE_;

L80:
    nl2itr_(&v[d1], &iv[1], &v[j1], n, n, p, &v[r1], &v[1], &x[1]);
    if ((i__1 = iv[1] - 2) < 0) {
	goto L10;
    } else if (i__1 == 0) {
	goto L30;
    } else {
	goto L999;
    }

L90:
    itsmry_(&v[d1], &iv[1], p, &v[1], &x[1]);

L999:
    return 0;
/*  ***  LAST CARD OF NL2SNO FOLLOWS  *** */
} /* nl2sno_ */

/* Subroutine */ int nl2itr_(doublereal *d__, integer *iv, doublereal *j, 
	integer *n, integer *nn, integer *p, doublereal *r__, doublereal *v, 
	doublereal *x)
{
    /* System generated locals */
    integer j_dim1, j_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static doublereal e;
    static integer i__, k, l, m;
    static doublereal t;
    static integer g1, h0, h1, s1;
    static doublereal t1;
    static integer w1, g01, x01, rd0, im1, rd1, km1, rdk, ipk, smh, dig1, 
	    lky1, qtr1, pp1o2;
    static doublereal rdof1;
    static integer lmat1, temp1, temp2, ipiv1, step1, ipivi, ipivk, dummy, 
	    sstep;
    extern /* Subroutine */ int vcopy_(integer *, doublereal *, doublereal *),
	     vaxpy_(integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    extern logical stopx_(integer *);
    static integer rsave1;
    extern doublereal v2norm_(integer *, doublereal *);
    extern /* Subroutine */ int parchk_(integer *, integer *, integer *, 
	    integer *, doublereal *), covclc_(integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *), qrfact_(integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, doublereal *), dupdat_(doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *);
    extern doublereal dotprd_(integer *, doublereal *, doublereal *);
    extern /* Subroutine */ int assess_(doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), lmstep_(doublereal *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static integer stpmod;
    extern /* Subroutine */ int qapply_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *), slupdt_(doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *);
    static integer lstgst;
    extern /* Subroutine */ int gqtstp_(doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *), rptmul_(integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), slvmul_(integer *, doublereal *, 
	    doublereal *, doublereal *), vscopy_(integer *, doublereal *, 
	    doublereal *), itsmry_(doublereal *, integer *, integer *, 
	    doublereal *, doublereal *);
    static doublereal sttsst;


/*  ***  CARRY OUT NL2SOL (NONLINEAR LEAST-SQUARES) ITERATIONS  *** */
/*  ***  (NL2SOL VERSION 2.2)  *** */

/*  ***  PARAMETER DECLARATIONS  *** */

/*     DIMENSION IV(60+P), V(93 + 2*N + P*(3*P+31)/2) */


/* --------------------------  PARAMETER USAGE  -------------------------- */

/* D.... SCALE VECTOR. */
/* IV... INTEGER VALUE ARRAY. */
/* J.... N BY P JACOBIAN MATRIX (LEAD DIMENSION NN). */
/* N.... NUMBER OF OBSERVATIONS (COMPONENTS IN R). */
/* NN... LEAD DIMENSION OF J. */
/* P.... NUMBER OF PARAMETERS (COMPONENTS IN X). */
/* R.... RESIDUAL VECTOR. */
/* V.... FLOATING-POINT VALUE ARRAY. */
/* X.... PARAMETER VECTOR. */

/*  ***  DISCUSSION  *** */

/*        PARAMETERS IV, N, P, V, AND X ARE THE SAME AS THE CORRESPOND- */
/*     ING ONES TO NL2SOL (WHICH SEE), EXCEPT THAT V CAN BE SHORTER */
/*     (SINCE THE PART OF V THAT NL2SOL USES FOR STORING D, J, AND R IS */
/*     NOT NEEDED).  MOREOVER, COMPARED WITH NL2SOL, IV(1) MAY HAVE THE */
/*     TWO ADDITIONAL OUTPUT VALUES 1 AND 2, WHICH ARE EXPLAINED BELOW, */
/*     AS IS THE USE OF IV(TOOBIG) AND IV(NFGCAL).  THE VALUES IV(D), */
/*     IV(J), AND IV(R), WHICH ARE OUTPUT VALUES FROM NL2SOL (AND */
/*     NL2SNO), ARE NOT REFERENCED BY NL2ITR OR THE SUBROUTINES IT CALLS. */
/*        ON A FRESH START, I.E., A CALL ON NL2ITR WITH IV(1) = 0 OR 12, */
/*     NL2ITR ASSUMES THAT R = R(X), THE RESIDUAL AT X, AND J = J(X), */
/*     THE CORRESPONDING JACOBIAN MATRIX OF R AT X. */

/* IV(1) = 1 MEANS THE CALLER SHOULD SET R TO R(X), THE RESIDUAL AT X, */
/*             AND CALL NL2ITR AGAIN, HAVING CHANGED NONE OF THE OTHER */
/*             PARAMETERS.  AN EXCEPTION OCCURS IF R CANNOT BE EVALUATED */
/*             AT X (E.G. IF R WOULD OVERFLOW), WHICH MAY HAPPEN BECAUSE */
/*             OF AN OVERSIZED STEP.  IN THIS CASE THE CALLER SHOULD SET */
/*             IV(TOOBIG) = IV(2) TO 1, WHICH WILL CAUSE NL2ITR TO IG- */
/*             NORE R AND TRY A SMALLER STEP.  THE PARAMETER NF THAT */
/*             NL2SOL PASSES TO CALCR (FOR POSSIBLE USE BY CALCJ) IS A */
/*             COPY OF IV(NFCALL) = IV(6). */
/* IV(1) = 2 MEANS THE CALLER SHOULD SET J TO J(X), THE JACOBIAN MATRIX */
/*             OF R AT X, AND CALL NL2ITR AGAIN.  THE CALLER MAY CHANGE */
/*             D AT THIS TIME, BUT SHOULD NOT CHANGE ANY OF THE OTHER */
/*             PARAMETERS.  THE PARAMETER NF THAT NL2SOL PASSES TO */
/*             CALCJ IS IV(NFGCAL) = IV(7).  IF J CANNOT BE EVALUATED */
/*             AT X, THEN THE CALLER MAY SET IV(NFGCAL) TO 0, IN WHICH */
/*             CASE NL2ITR WILL RETURN WITH IV(1) = 15. */

/*  ***  GENERAL  *** */

/*     CODED BY DAVID M. GAY. */
/*     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH */
/*     SUPPORTED BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS */

/*     MCS-7600324, DCR75-10143, 76-14311DSS, MCS76-11989, AND */
/*     MCS-7906671. */
/*        (SEE NL2SOL FOR REFERENCES.) */

/* +++++++++++++++++++++++++++  DECLARATIONS  ++++++++++++++++++++++++++++ */

/*  ***  LOCAL VARIABLES  *** */


/*     ***  CONSTANTS  *** */


/*  ***  INTRINSIC FUNCTIONS  *** */
/* /+ */
/* / */
/*  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  *** */


/* ASSESS... ASSESSES CANDIDATE STEP. */
/* COVCLC... COMPUTES COVARIANCE MATRIX. */
/* DOTPRD... RETURNS INNER PRODUCT OF TWO VECTORS. */
/* DUPDAT... UPDATES SCALE VECTOR D. */
/* GQTSTP... COMPUTES GOLDFELD-QUANDT-TROTTER STEP (AUGMENTED MODEL). */
/* ITSMRY... PRINTS ITERATION SUMMARY AND INFO ABOUT INITIAL AND FINAL X. */
/* LMSTEP... COMPUTES LEVENBERG-MARQUARDT STEP (GAUSS-NEWTON MODEL). */
/* PARCHK... CHECKS VALIDITY OF INPUT IV AND V VALUES. */
/* QAPPLY... APPLIES ORTHOGONAL MATRIX Q FROM QRFACT TO A VECTOR. */
/* QRFACT... COMPUTES QR DECOMPOSITION OF A MATRIX VIA HOUSEHOLDER TRANS. */
/* RPTMUL... MULTIPLIES VECTOR BY THE R MATRIX (AND/OR ITS TRANSPOSE) */
/*             STORED BY QRFACT. */
/* SLUPDT... PERFORMS QUASI-NEWTON UPDATE ON COMPACTLY STORED LOWER TRI- */
/*             ANGLE OF A SYMMETRIC MATRIX. */
/* STOPX.... RETURNS .TRUE. IF THE BREAK KEY HAS BEEN PRESSED. */
/* VAXPY.... COMPUTES SCALAR TIMES ONE VECTOR PLUS ANOTHER. */
/* VCOPY.... COPIES ONE VECTOR TO ANOTHER. */
/* VSCOPY... SETS ALL ELEMENTS OF A VECTOR TO A SCALAR. */
/* V2NORM... RETURNS THE 2-NORM OF A VECTOR. */

/*  ***  SUBSCRIPTS FOR IV AND V  *** */


/*  ***  IV SUBSCRIPT VALUES  *** */

/* /6 */
/*     DATA CNVCOD/34/, COVMAT/26/, COVPRT/14/, */
/*    1     COVREQ/15/, DIG/43/, DTYPE/16/, G/28/, H/44/, */
/*    2     IERR/32/, INITS/25/, IPIVOT/61/, IPIV0/60/, */
/*    3     IRC/3/, KAGQT/35/, KALM/36/, LKY/37/, LMAT/58/, */
/*    4     MODE/38/, MODEL/5/, MXFCAL/17/, MXITER/18/, */
/*    5     NFCALL/6/, NFGCAL/7/, NFCOV/40/, NGCOV/41/, */
/*    6     NGCALL/30/, NITER/31/, QTR/49/, */
/*    7     RADINC/8/, RD/51/, RESTOR/9/, RSAVE/52/, S/53/, */
/*    8     STEP/55/, STGLIM/11/, STLSTG/56/, SUSED/57/, */
/*    9     SWITCH/12/, TOOBIG/2/, W/59/, XIRC/13/, X0/60/ */
/* /7 */
/* / */

/*  ***  V SUBSCRIPT VALUES  *** */

/* /6 */
/*     DATA COSMIN/43/, DGNORM/1/, DINIT/38/, DSTNRM/2/, */
/*    1     D0INIT/37/, F/10/, FDIF/11/, FUZZ/45/, */
/*    2     F0/13/, GTSTEP/4/, INCFAC/23/, */
/*    3     JTINIT/39/, JTOL1/87/, LMAX0/35/, */
/*    4     NVSAVE/9/, PHMXFC/21/, PREDUC/7/, */
/*    5     RADFAC/16/, RADIUS/8/, RAD0/9/, RLIMIT/42/, */
/*    6     SIZE/47/, STPPAR/5/, TUNER4/29/, TUNER5/30/, */
/*    7     VSAVE1/78/, WSCALE/48/ */
/* /7 */
/* / */


/* /6 */
/*     DATA HALF/0.5D+0/, NEGONE/-1.D+0/, ONE/1.D+0/, ZERO/0.D+0/ */
/* /7 */
/* / */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

    /* Parameter adjustments */
    --iv;
    --r__;
    --x;
    j_dim1 = *nn;
    j_offset = 1 + j_dim1;
    j -= j_offset;
    --d__;
    --v;

    /* Function Body */
    i__ = iv[1];
    if (i__ == 1) {
	goto L20;
    }
    if (i__ == 2) {
	goto L50;
    }

/*  ***  CHECK VALIDITY OF IV AND V INPUT VALUES  *** */

/*     ***  NOTE -- IF IV(1) = 0, THEN PARCHK CALLS DFAULT(IV, V)  *** */
    parchk_(&iv[1], n, nn, p, &v[1]);   //call to parchk_
    i__ = iv[1] - 2;
    if (i__ > 10) {
	goto L999;
    }
    switch (i__) {
	case 1:  goto L350;
	case 2:  goto L350;
	case 3:  goto L350;
	case 4:  goto L350;
	case 5:  goto L350;
	case 6:  goto L350;
	case 7:  goto L195;
	case 8:  goto L160;
	case 9:  goto L195;
	case 10:  goto L10;
    }

/*  ***  INITIALIZATION AND STORAGE ALLOCATION  *** */

L10:
    iv[31] = 0;
    iv[6] = 1;
    iv[30] = 1;
    iv[7] = 1;
    iv[38] = -1;
    iv[11] = 2;
    iv[2] = 0;
    iv[34] = 0;
    iv[26] = 0;
    iv[40] = 0;
    iv[41] = 0;
    iv[36] = -1;
    iv[8] = 0;
    iv[53] = (*p << 1) + 87;
    pp1o2 = *p * (*p + 1) / 2;
    iv[60] = iv[53] + pp1o2;
    iv[55] = iv[60] + *p;
    iv[56] = iv[55] + *p;
    iv[43] = iv[56] + *p;
    iv[28] = iv[43] + *p;
    iv[37] = iv[28] + *p;
    iv[51] = iv[37] + *p;
    iv[52] = iv[51] + *p;
    iv[49] = iv[52] + *n;
    iv[44] = iv[49] + *n;
    iv[59] = iv[44] + pp1o2;
    iv[58] = iv[59] + (*p << 2) + 7;
/*     +++ LENGTH OF W = P*(P+9)/2 + 7.  LMAT IS CONTAINED IN W. */
    if (v[38] >= 0.) {
	vscopy_(p, &d__[1], &v[38]);
    }
    if (v[39] > 0.) {
	vscopy_(p, &v[87], &v[39]);
    }
    i__ = *p + 87;
    if (v[37] > 0.) {
	vscopy_(p, &v[i__], &v[37]);
    }
    v[9] = 0.;
    v[5] = 0.;
    v[8] = v[35] / (v[21] + 1.);

/*  ***  SET INITIAL MODEL AND S MATRIX  *** */

    iv[5] = 1;
    if (iv[25] == 2) {
	iv[5] = 2;
    }
    s1 = iv[53];
    if (iv[25] == 0) {
	vscopy_(&pp1o2, &v[s1], &c_b95);
    }

/*  ***  COMPUTE FUNCTION VALUE (HALF THE SUM OF SQUARES)  *** */

L20:
    t = v2norm_(n, &r__[1]);
    if (t > v[42]) {
	iv[2] = 1;
    }
    if (iv[2] != 0) {
	goto L30;
    }
/* Computing 2nd power */
    d__1 = t;
    v[10] = d__1 * d__1 * .5;
L30:
    if (iv[38] < 0) {
	goto L40;
    } else if (iv[38] == 0) {
	goto L350;
    } else {
	goto L730;
    }

L40:
    if (iv[2] == 0) {
	goto L60;
    }
    iv[1] = 13;
    goto L900;

/*  ***  MAKE SURE JACOBIAN COULD BE COMPUTED  *** */

L50:
    if (iv[7] != 0) {
	goto L60;
    }
    iv[1] = 15;
    goto L900;

/*  ***  COMPUTE GRADIENT  *** */

L60:
    iv[36] = -1;
    g1 = iv[28];
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[g1] = dotprd_(n, &r__[1], &j[i__ * j_dim1 + 1]);
	++g1;
/* L70: */
    }
    if (iv[38] > 0) {
	goto L710;
    }

/*  ***  UPDATE D AND MAKE COPIES OF R FOR POSSIBLE USE LATER  *** */

    if (iv[16] > 0) {
	dupdat_(&d__[1], &iv[1], &j[j_offset], n, nn, p, &v[1]);
    }
    rsave1 = iv[52];
    vcopy_(n, &v[rsave1], &r__[1]);
    qtr1 = iv[49];
    vcopy_(n, &v[qtr1], &r__[1]);

/*  ***  COMPUTE  D**-1 * GRADIENT  *** */

    g1 = iv[28];
    dig1 = iv[43];
    k = dig1;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[k] = v[g1] / d__[i__];
	++k;
	++g1;
/* L80: */
    }
    v[1] = v2norm_(p, &v[dig1]);

    if (iv[34] != 0) {
	goto L700;
    }
    if (iv[38] == 0) {
	goto L570;
    }
    iv[38] = 0;


/* -----------------------------  MAIN LOOP  ----------------------------- */


/*  ***  PRINT ITERATION SUMMARY, CHECK ITERATION LIMIT  *** */

L150:
    itsmry_(&d__[1], &iv[1], p, &v[1], &x[1]);
L160:
    k = iv[31];
    if (k < iv[18]) {
	goto L170;
    }
    iv[1] = 10;
    goto L900;
L170:
    iv[31] = k + 1;

/*  ***  UPDATE RADIUS  *** */

    if (k == 0) {
	goto L185;
    }
    step1 = iv[55];
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[step1] = d__[i__] * v[step1];
	++step1;
/* L180: */
    }
    step1 = iv[55];
    v[8] = v[16] * v2norm_(p, &v[step1]);

/*  ***  INITIALIZE FOR START OF NEXT ITERATION  *** */

L185:
    x01 = iv[60];
    v[13] = v[10];
    iv[35] = -1;
    iv[3] = 4;
    iv[44] = -abs(iv[44]);
    iv[57] = iv[5];

/*     ***  COPY X TO X0  *** */

    vcopy_(p, &v[x01], &x[1]);

/*  ***  CHECK STOPX AND FUNCTION EVALUATION LIMIT  *** */

L190:
    if (! stopx_(&dummy)) {
	goto L200;
    }
    iv[1] = 11;
    goto L205;

/*     ***  COME HERE WHEN RESTARTING AFTER FUNC. EVAL. LIMIT OR STOPX. */

L195:
    if (v[10] >= v[13]) {
	goto L200;
    }
    v[16] = 1.;
    k = iv[31];
    goto L170;

L200:
    if (iv[6] < iv[17] + iv[40]) {
	goto L210;
    }
    iv[1] = 9;
L205:
    if (v[10] >= v[13]) {
	goto L900;
    }

/*        ***  IN CASE OF STOPX OR FUNCTION EVALUATION LIMIT WITH */
/*        ***  IMPROVED V(F), EVALUATE THE GRADIENT AT X. */

    iv[34] = iv[1];
    goto L560;

/* . . . . . . . . . . . . .  COMPUTE CANDIDATE STEP  . . . . . . . . . . */

L210:
    step1 = iv[55];
    w1 = iv[59];
    if (iv[5] == 2) {
	goto L240;
    }

/*  ***  COMPUTE LEVENBERG-MARQUARDT STEP  *** */

    qtr1 = iv[49];
    if (iv[36] >= 0) {
	goto L215;
    }
    rd1 = iv[51];
    if (-1 == iv[36]) {
	qrfact_(nn, n, p, &j[j_offset], &v[rd1], &iv[61], &iv[32], &c__0, &v[
		w1]);
    }
    qapply_(nn, n, p, &j[j_offset], &v[qtr1], &iv[32]);
L215:
    h1 = iv[44];
    if (h1 > 0) {
	goto L230;
    }

/*        ***  COPY R MATRIX TO H  *** */

    h1 = -h1;
    iv[44] = h1;
    k = h1;
    rd1 = iv[51];
    v[k] = v[rd1];
    if (*p == 1) {
	goto L230;
    }
    i__1 = *p;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i__2 = i__ - 1;
	vcopy_(&i__2, &v[k + 1], &j[i__ * j_dim1 + 1]);
	k += i__;
	++rd1;
	v[k] = v[rd1];
/* L220: */
    }

L230:
    g1 = iv[28];
    lmstep_(&d__[1], &v[g1], &iv[32], &iv[61], &iv[36], p, &v[qtr1], &v[h1], &
	    v[step1], &v[1], &v[w1]);
    goto L310;

/*  ***  COMPUTE GOLDFELD-QUANDT-TROTTER STEP (AUGMENTED MODEL)  *** */

L240:
    if (iv[44] > 0) {
	goto L300;
    }

/*     ***  SET H TO  D**-1 * ( (J**T)*J + S) ) * D**-1.  *** */

    h1 = -iv[44];
    iv[44] = h1;
    s1 = iv[53];
    if (-1 != iv[36]) {
	goto L270;
    }

/*        ***  J IS IN ITS ORIGINAL FORM  *** */

    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t = 1. / d__[i__];
	i__2 = i__;
	for (k = 1; k <= i__2; ++k) {
	    v[h1] = t * (dotprd_(n, &j[i__ * j_dim1 + 1], &j[k * j_dim1 + 1]) 
		    + v[s1]) / d__[k];
	    ++h1;
	    ++s1;
/* L250: */
	}
/* L260: */
    }
    goto L300;

/*  ***  LMSTEP HAS APPLIED QRFACT TO J  *** */

L270:
    smh = s1 - h1;
    h0 = h1 - 1;
    ipiv1 = iv[61];
    t1 = 1. / d__[ipiv1];
    rd0 = iv[51] - 1;
    rdof1 = v[rd0 + 1];
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l = i__ + 60;
	ipivi = iv[l];
	h1 = h0 + ipivi * (ipivi - 1) / 2;
	l = h1 + ipivi;
	m = l + smh;
/*             ***  V(L) = H(IPIVOT(I), IPIVOT(I))  *** */
/*             ***  V(M) = S(IPIVOT(I), IPIVOT(I))  *** */
	t = 1. / d__[ipivi];
	rdk = rd0 + i__;
/* Computing 2nd power */
	d__1 = v[rdk];
	e = d__1 * d__1;
	if (i__ > 1) {
	    i__2 = i__ - 1;
	    e += dotprd_(&i__2, &j[i__ * j_dim1 + 1], &j[i__ * j_dim1 + 1]);
	}
/* Computing 2nd power */
	d__1 = t;
	v[l] = (e + v[m]) * (d__1 * d__1);
	if (i__ == 1) {
	    goto L290;
	}
	l = h1 + ipiv1;
	if (ipivi < ipiv1) {
	    l += (ipiv1 - ipivi) * (ipiv1 + ipivi - 3) / 2;
	}
	m = l + smh;
/*             ***  V(L) = H(IPIVOT(I), IPIVOT(1))  *** */
/*             ***  V(M) = S(IPIVOT(I), IPIVOT(1))  *** */
	v[l] = t * (rdof1 * j[i__ * j_dim1 + 1] + v[m]) * t1;
	if (i__ == 2) {
	    goto L290;
	}
	im1 = i__ - 1;
	i__2 = im1;
	for (k = 2; k <= i__2; ++k) {
	    ipk = k + 60;
	    ipivk = iv[ipk];
	    l = h1 + ipivk;
	    if (ipivi < ipivk) {
		l += (ipivk - ipivi) * (ipivk + ipivi - 3) / 2;
	    }
	    m = l + smh;
/*                  ***  V(L) = H(IPIVOT(I), IPIVOT(K))  *** */
/*                  ***  V(M) = S(IPIVOT(I), IPIVOT(K))  *** */
	    km1 = k - 1;
	    rdk = rd0 + k;
	    v[l] = t * (dotprd_(&km1, &j[i__ * j_dim1 + 1], &j[k * j_dim1 + 1]
		    ) + v[rdk] * j[k + i__ * j_dim1] + v[m]) / d__[ipivk];
/* L280: */
	}
L290:
	;
    }

/*  ***  COMPUTE ACTUAL GOLDFELD-QUANDT-TROTTER STEP  *** */

L300:
    h1 = iv[44];
    dig1 = iv[43];
    lmat1 = iv[58];
    gqtstp_(&d__[1], &v[dig1], &v[h1], &iv[35], &v[lmat1], p, &v[step1], &v[1]
	    , &v[w1]);


/*  ***  COMPUTE R(X0 + STEP)  *** */

L310:
    if (iv[3] == 6) {
	goto L350;
    }
    x01 = iv[60];
    step1 = iv[55];
    vaxpy_(p, &x[1], &c_b76, &v[step1], &v[x01]);
    ++iv[6];
    iv[1] = 1;
    iv[2] = 0;
    goto L999;

/* . . . . . . . . . . . . .  ASSESS CANDIDATE STEP  . . . . . . . . . . . */

L350:
    step1 = iv[55];
    lstgst = iv[56];
    x01 = iv[60];
    assess_(&d__[1], &iv[1], p, &v[step1], &v[lstgst], &v[1], &x[1], &v[x01]);

/*  ***  IF NECESSARY, SWITCH MODELS AND/OR RESTORE R  *** */

    if (iv[12] == 0) {
	goto L360;
    }
    iv[44] = -abs(iv[44]);
    iv[57] += 2;
    vcopy_(&c__9, &v[1], &v[78]);
L360:
    if (iv[9] == 0) {
	goto L390;
    }
    rsave1 = iv[52];
    vcopy_(n, &r__[1], &v[rsave1]);
L390:
    l = iv[3] - 4;
    stpmod = iv[5];
    if (l > 0) {
	switch (l) {
	    case 1:  goto L410;
	    case 2:  goto L440;
	    case 3:  goto L450;
	    case 4:  goto L450;
	    case 5:  goto L450;
	    case 6:  goto L450;
	    case 7:  goto L450;
	    case 8:  goto L450;
	    case 9:  goto L640;
	    case 10:  goto L570;
	}
    }

/*  ***  DECIDE WHETHER TO CHANGE MODELS  *** */

    e = v[7] - v[11];
    sstep = iv[37];
    s1 = iv[53];
    slvmul_(p, &v[sstep], &v[s1], &v[step1]);
    sttsst = dotprd_(p, &v[step1], &v[sstep]) * .5;
    if (iv[5] == 1) {
	sttsst = -sttsst;
    }
    if ((d__1 = e + sttsst, abs(d__1)) * v[45] >= abs(e)) {
	goto L400;
    }

/*     ***  SWITCH MODELS  *** */

    iv[5] = 3 - iv[5];
    if (iv[5] == 1) {
	iv[35] = -1;
    }
    if (iv[5] == 2 && iv[36] > 0) {
	iv[36] = 0;
    }
    if (-2 < l) {
	goto L480;
    }
    iv[44] = -abs(iv[44]);
    iv[57] += 2;
    vcopy_(&c__9, &v[78], &v[1]);
    goto L420;

L400:
    if (-3 < l) {
	goto L480;
    }

/*     ***  RECOMPUTE STEP WITH DECREASED RADIUS  *** */

    v[8] = v[16] * v[2];
    goto L190;

/*  ***  RECOMPUTE STEP, SAVING V VALUES AND R IF NECESSARY  *** */

L410:
    v[8] = v[16] * v[2];
L420:
    if (v[10] >= v[13]) {
	goto L190;
    }
    rsave1 = iv[52];
    vcopy_(n, &v[rsave1], &r__[1]);
    goto L190;

/*  ***  COMPUTE STEP OF LENGTH V(LMAX0) FOR SINGULAR CONVERGENCE TEST */

L440:
    v[8] = v[35];
    goto L210;

/*  ***  CONVERGENCE OR FALSE CONVERGENCE  *** */

L450:
    iv[34] = l;
    if (v[10] >= v[13]) {
	goto L700;
    }
    if (iv[13] == 14) {
	goto L700;
    }
    iv[13] = 14;

/* . . . . . . . . . . . .  PROCESS ACCEPTABLE STEP  . . . . . . . . . . . */

L480:
    iv[26] = 0;

/*  ***  SET  LKY = (J(X0)**T) * R(X)  *** */

    lky1 = iv[37];
    if (iv[36] >= 0) {
	goto L500;
    }

/*     ***  JACOBIAN HAS NOT BEEN MODIFIED  *** */

    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[lky1] = dotprd_(n, &j[i__ * j_dim1 + 1], &r__[1]);
	++lky1;
/* L490: */
    }
    goto L510;

/*  ***  QRFACT HAS BEEN APPLIED TO J.  STORE COPY OF R IN QTR AND  *** */
/*  ***  APPLY Q TO IT.                                             *** */

L500:
    qtr1 = iv[49];
    vcopy_(n, &v[qtr1], &r__[1]);
    qapply_(nn, n, p, &j[j_offset], &v[qtr1], &iv[32]);

/*  ***  MULTIPLY TOP P-VECTOR IN QTR BY PERMUTED UPPER TRIANGLE    *** */
/*  ***  STORED BY QRFACT IN J AND RD.                              *** */

    rd1 = iv[51];
    temp1 = iv[56];
    rptmul_(&c__3, &iv[61], &j[j_offset], nn, p, &v[rd1], &v[qtr1], &v[lky1], 
	    &v[temp1]);

/*  ***  SEE WHETHER TO SET V(RADFAC) BY GRADIENT TESTS  *** */

L510:
    if (iv[3] != 3) {
	goto L560;
    }
    step1 = iv[55];
    temp1 = iv[56];
    temp2 = iv[60];

/*     ***  SET  TEMP1 = HESSIAN * STEP  FOR USE IN GRADIENT TESTS  *** */

    if (stpmod == 2) {
	goto L530;
    }

/*        ***  STEP COMPUTED USING GAUSS-NEWTON MODEL  *** */
/*        ***  -- QRFACT HAS BEEN APPLIED TO J         *** */

    rd1 = iv[51];
    rptmul_(&c__2, &iv[61], &j[j_offset], nn, p, &v[rd1], &v[step1], &v[temp1]
	    , &v[temp2]);
    goto L560;

/*     ***  STEP COMPUTED USING AUGMENTED MODEL  *** */

L530:
    h1 = iv[44];
    k = temp2;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[k] = d__[i__] * v[step1];
	++k;
	++step1;
/* L540: */
    }
    slvmul_(p, &v[temp1], &v[h1], &v[temp2]);
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[temp1] = d__[i__] * v[temp1];
	++temp1;
/* L550: */
    }

/*  ***  SAVE OLD GRADIENT AND COMPUTE NEW ONE  *** */

L560:
    ++iv[30];
    g1 = iv[28];
    g01 = iv[59];
    vcopy_(p, &v[g01], &v[g1]);
    iv[1] = 2;
    goto L999;

/*  ***  INITIALIZATIONS -- G0 = G - G0, ETC.  *** */

L570:
    g01 = iv[59];
    g1 = iv[28];
    vaxpy_(p, &v[g01], &c_b147, &v[g01], &v[g1]);
    step1 = iv[55];
    temp1 = iv[56];
    temp2 = iv[60];
    if (iv[3] != 3) {
	goto L600;
    }

/*  ***  SET V(RADFAC) BY GRADIENT TESTS  *** */

/*     ***  SET  TEMP1 = D**-1 * (HESSIAN * STEP  +  (G(X0) - G(X)))  *** */

    k = temp1;
    l = g01;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[k] = (v[k] - v[l]) / d__[i__];
	++k;
	++l;
/* L580: */
    }

/*        ***  DO GRADIENT TESTS  *** */

    if (v2norm_(p, &v[temp1]) <= v[1] * v[29]) {
	goto L590;
    }
    if (dotprd_(p, &v[g1], &v[step1]) >= v[4] * v[30]) {
	goto L600;
    }
L590:
    v[16] = v[23];

/*  ***  FINISH COMPUTING LKY = ((J(X) - J(X0))**T) * R  *** */

/*     ***  CURRENTLY LKY = (J(X0)**T) * R  *** */

L600:
    lky1 = iv[37];
    vaxpy_(p, &v[lky1], &c_b147, &v[lky1], &v[g1]);

/*  ***  DETERMINE SIZING FACTOR V(SIZE)  *** */

/*     ***  SET TEMP1 = S * STEP  *** */
    s1 = iv[53];
    slvmul_(p, &v[temp1], &v[s1], &v[step1]);

    t1 = (d__1 = dotprd_(p, &v[step1], &v[temp1]), abs(d__1));
    t = (d__1 = dotprd_(p, &v[step1], &v[lky1]), abs(d__1));
    v[47] = 1.;
    if (t < t1) {
	v[47] = t / t1;
    }

/*  ***  UPDATE S  *** */

    slupdt_(&v[s1], &v[43], p, &v[47], &v[step1], &v[temp1], &v[temp2], &v[
	    g01], &v[48], &v[lky1]);
    iv[1] = 2;
    goto L150;

/* . . . . . . . . . . . . . .  MISC. DETAILS  . . . . . . . . . . . . . . */

/*  ***  BAD PARAMETERS TO ASSESS  *** */

L640:
    iv[1] = 14;
    goto L900;

/*  ***  CONVERGENCE OBTAINED -- COMPUTE COVARIANCE MATRIX IF DESIRED *** */

L700:
    if (iv[15] == 0 && iv[14] == 0) {
	goto L760;
    }
    if (iv[26] != 0) {
	goto L760;
    }
    if (iv[34] >= 7) {
	goto L760;
    }
    iv[38] = 0;
L710:
    covclc_(&i__, &d__[1], &iv[1], &j[j_offset], n, nn, p, &r__[1], &v[1], &x[
	    1]);
    switch (i__) {
	case 1:  goto L720;
	case 2:  goto L720;
	case 3:  goto L740;
	case 4:  goto L750;
    }
L720:
    ++iv[40];
    ++iv[6];
    iv[9] = i__;
    iv[1] = 1;
    goto L999;

L730:
    if (iv[9] == 1 || iv[2] != 0) {
	goto L710;
    }
    iv[7] = iv[6];
L740:
    ++iv[41];
    ++iv[30];
    iv[1] = 2;
    goto L999;

L750:
    iv[38] = 0;
    if (iv[31] == 0) {
	iv[38] = -1;
    }

L760:
    iv[1] = iv[34];
    iv[34] = 0;

/*  ***  PRINT SUMMARY OF FINAL ITERATION AND OTHER REQUESTED ITEMS  *** */

L900:
    itsmry_(&d__[1], &iv[1], p, &v[1], &x[1]);

L999:
    return 0;

/*  ***  LAST CARD OF NL2ITR FOLLOWS  *** */
} /* nl2itr_ */

/* Subroutine */ int assess_(doublereal *d__, integer *iv, integer *p, 
	doublereal *step, doublereal *stlstg, doublereal *v, doublereal *x, 
	doublereal *x0)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, nfc;
    static doublereal gts, emax, xmax, rfac1;
    static logical goodx;
    extern /* Subroutine */ int vcopy_(integer *, doublereal *, doublereal *);
    static doublereal reldx1;
    extern doublereal reldst_(integer *, doublereal *, doublereal *, 
	    doublereal *);


/*  ***  ASSESS CANDIDATE STEP (NL2SOL VERSION 2.2)  *** */


/*  ***  PURPOSE  *** */

/*        THIS SUBROUTINE IS CALLED BY AN UNCONSTRAINED MINIMIZATION */
/*     ROUTINE TO ASSESS THE NEXT CANDIDATE STEP.  IT MAY RECOMMEND ONE */
/*     OF SEVERAL COURSES OF ACTION, SUCH AS ACCEPTING THE STEP, RECOM- */
/*     PUTING IT USING THE SAME OR A NEW QUADRATIC MODEL, OR HALTING DUE */
/*     TO CONVERGENCE OR FALSE CONVERGENCE.  SEE THE RETURN CODE LISTING */
/*     BELOW. */

/* --------------------------  PARAMETER USAGE  -------------------------- */

/*     IV (I/O) INTEGER PARAMETER AND SCRATCH VECTOR -- SEE DESCRIPTION */
/*             BELOW OF IV VALUES REFERENCED. */
/*      D (IN)  SCALE VECTOR USED IN COMPUTING V(RELDX) -- SEE BELOW. */
/*      P (IN)  NUMBER OF PARAMETERS BEING OPTIMIZED. */
/*   STEP (I/O) ON INPUT, STEP IS THE STEP TO BE ASSESSED.  IT IS UN- */
/*             CHANGED ON OUTPUT UNLESS A PREVIOUS STEP ACHIEVED A */
/*             BETTER OBJECTIVE FUNCTION REDUCTION, IN WHICH CASE STLSTG */
/*             WILL HAVE BEEN COPIED TO STEP. */
/* STLSTG (I/O) WHEN ASSESS RECOMMENDS RECOMPUTING STEP EVEN THOUGH THE */
/*             CURRENT (OR A PREVIOUS) STEP YIELDS AN OBJECTIVE FUNC- */
/*             TION DECREASE, IT SAVES IN STLSTG THE STEP THAT GAVE THE */
/*             BEST FUNCTION REDUCTION SEEN SO FAR (IN THE CURRENT ITERA- */
/*             TION).  IF THE RECOMPUTED STEP YIELDS A LARGER FUNCTION */
/*             VALUE, THEN STEP IS RESTORED FROM STLSTG AND */
/*             X = X0 + STEP IS RECOMPUTED. */
/*      V (I/O) REAL PARAMETER AND SCRATCH VECTOR -- SEE DESCRIPTION */
/*             BELOW OF V VALUES REFERENCED. */
/*      X (I/O) ON INPUT, X = X0 + STEP IS THE POINT AT WHICH THE OBJEC- */
/*             TIVE FUNCTION HAS JUST BEEN EVALUATED.  IF AN EARLIER */
/*             STEP YIELDED A BIGGER FUNCTION DECREASE, THEN X IS */
/*             RESTORED TO THE CORRESPONDING EARLIER VALUE.  OTHERWISE, */
/*             IF THE CURRENT STEP DOES NOT GIVE ANY FUNCTION DECREASE, */
/*             THEN X IS RESTORED TO X0. */
/*     X0 (IN)  INITIAL OBJECTIVE FUNCTION PARAMETER VECTOR (AT THE */
/*             START OF THE CURRENT ITERATION). */

/*  ***  IV VALUES REFERENCED  *** */

/*    IV(IRC) (I/O) ON INPUT FOR THE FIRST STEP TRIED IN A NEW ITERATION, */
/*             IV(IRC) SHOULD BE SET TO 3 OR 4 (THE VALUE TO WHICH IT IS */
/*             SET WHEN STEP IS DEFINITELY TO BE ACCEPTED).  ON INPUT */
/*             AFTER STEP HAS BEEN RECOMPUTED, IV(IRC) SHOULD BE */
/*             UNCHANGED SINCE THE PREVIOUS RETURN OF ASSESS. */
/*                ON OUTPUT, IV(IRC) IS A RETURN CODE HAVING ONE OF THE */
/*             FOLLOWING VALUES... */
/*                  1 = SWITCH MODELS OR TRY SMALLER STEP. */
/*                  2 = SWITCH MODELS OR ACCEPT STEP. */
/*                  3 = ACCEPT STEP AND DETERMINE V(RADFAC) BY GRADIENT */
/*                       TESTS. */
/*                  4 = ACCEPT STEP, V(RADFAC) HAS BEEN DETERMINED. */
/*                  5 = RECOMPUTE STEP (USING THE SAME MODEL). */
/*                  6 = RECOMPUTE STEP WITH RADIUS = V(LMAX0) BUT DO NOT */
/*                       EVAULATE THE OBJECTIVE FUNCTION. */
/*                  7 = X-CONVERGENCE (SEE V(XCTOL)). */
/*                  8 = RELATIVE FUNCTION CONVERGENCE (SEE V(RFCTOL)). */
/*                  9 = BOTH X- AND RELATIVE FUNCTION CONVERGENCE. */
/*                 10 = ABSOLUTE FUNCTION CONVERGENCE (SEE V(AFCTOL)). */
/*                 11 = SINGULAR CONVERGENCE (SEE V(LMAX0)). */
/*                 12 = FALSE CONVERGENCE (SEE V(XFTOL)). */
/*                 13 = IV(IRC) WAS OUT OF RANGE ON INPUT. */
/*             RETURN CODE I HAS PRECDENCE OVER I+1 FOR I = 9, 10, 11. */
/* IV(MLSTGD) (I/O) SAVED VALUE OF IV(MODEL). */
/*  IV(MODEL) (I/O) ON INPUT, IV(MODEL) SHOULD BE AN INTEGER IDENTIFYING */
/*             THE CURRENT QUADRATIC MODEL OF THE OBJECTIVE FUNCTION. */
/*             IF A PREVIOUS STEP YIELDED A BETTER FUNCTION REDUCTION, */
/*             THEN IV(MODEL) WILL BE SET TO IV(MLSTGD) ON OUTPUT. */
/* IV(NFCALL) (IN)  INVOCATION COUNT FOR THE OBJECTIVE FUNCTION. */
/* IV(NFGCAL) (I/O) VALUE OF IV(NFCALL) AT STEP THAT GAVE THE BIGGEST */
/*             FUNCTION REDUCTION THIS ITERATION.  IV(NFGCAL) REMAINS */
/*             UNCHANGED UNTIL A FUNCTION REDUCTION IS OBTAINED. */
/* IV(RADINC) (I/O) THE NUMBER OF RADIUS INCREASES (OR MINUS THE NUMBER */
/*             OF DECREASES) SO FAR THIS ITERATION. */
/* IV(RESTOR) (OUT) SET TO 0 UNLESS X AND V(F) HAVE BEEN RESTORED, IN */
/*             WHICH CASE ASSESS SETS IV(RESTOR) = 1. */
/*  IV(STAGE) (I/O) COUNT OF THE NUMBER OF MODELS TRIED SO FAR IN THE */
/*             CURRENT ITERATION. */
/* IV(STGLIM) (IN)  MAXIMUM NUMBER OF MODELS TO CONSIDER. */
/* IV(SWITCH) (OUT) SET TO 0 UNLESS A NEW MODEL IS BEING TRIED AND IT */
/*             GIVES A SMALLER FUNCTION VALUE THAN THE PREVIOUS MODEL, */
/*             IN WHICH CASE ASSESS SETS IV(SWITCH) = 1. */
/* IV(TOOBIG) (IN)  IS NONZERO IF STEP WAS TOO BIG (E.G. IF IT CAUSED */
/*             OVERFLOW). */
/*   IV(XIRC) (I/O) VALUE THAT IV(IRC) WOULD HAVE IN THE ABSENCE OF */
/*             CONVERGENCE, FALSE CONVERGENCE, AND OVERSIZED STEPS. */

/*  ***  V VALUES REFERENCED  *** */

/* V(AFCTOL) (IN)  ABSOLUTE FUNCTION CONVERGENCE TOLERANCE.  IF THE */
/*             ABSOLUTE VALUE OF THE CURRENT FUNCTION VALUE V(F) IS LESS */
/*             THAN V(AFCTOL), THEN ASSESS RETURNS WITH IV(IRC) = 10. */
/* V(DECFAC) (IN)  FACTOR BY WHICH TO DECREASE RADIUS WHEN IV(TOOBIG) IS */
/*             NONZERO. */
/* V(DSTNRM) (IN)  THE 2-NORM OF D*STEP. */
/* V(DSTSAV) (I/O) VALUE OF V(DSTNRM) ON SAVED STEP. */
/*   V(DST0) (IN)  THE 2-NORM OF D TIMES THE NEWTON STEP (WHEN DEFINED, */
/*             I.E., FOR V(NREDUC) .GE. 0). */
/*      V(F) (I/O) ON BOTH INPUT AND OUTPUT, V(F) IS THE OBJECTIVE FUNC- */
/*             TION VALUE AT X.  IF X IS RESTORED TO A PREVIOUS VALUE, */
/*             THEN V(F) IS RESTORED TO THE CORRESPONDING VALUE. */
/*   V(FDIF) (OUT) THE FUNCTION REDUCTION V(F0) - V(F) (FOR THE OUTPUT */
/*             VALUE OF V(F) IF AN EARLIER STEP GAVE A BIGGER FUNCTION */
/*             DECREASE, AND FOR THE INPUT VALUE OF V(F) OTHERWISE). */
/* V(FLSTGD) (I/O) SAVED VALUE OF V(F). */
/*     V(F0) (IN)  OBJECTIVE FUNCTION VALUE AT START OF ITERATION. */
/* V(GTSLST) (I/O) VALUE OF V(GTSTEP) ON SAVED STEP. */
/* V(GTSTEP) (IN)  INNER PRODUCT BETWEEN STEP AND GRADIENT. */
/* V(INCFAC) (IN)  MINIMUM FACTOR BY WHICH TO INCREASE RADIUS. */
/*  V(LMAX0) (IN)  MAXIMUM REASONABLE STEP SIZE (AND INITIAL STEP BOUND). */
/*             IF THE ACTUAL FUNCTION DECREASE IS NO MORE THAN TWICE */
/*             WHAT WAS PREDICTED, IF A RETURN WITH IV(IRC) = 7, 8, 9, */
/*             OR 10 DOES NOT OCCUR, IF V(DSTNRM) .GT. V(LMAX0), AND IF */
/*             V(PREDUC) .LE. V(RFCTOL) * ABS(V(F0)), THEN ASSESS RE- */
/*             TURNS WITH IV(IRC) = 11.  IF SO DOING APPEARS WORTHWHILE, */
/*             THEN ASSESS REPEATS THIS TEST WITH V(PREDUC) COMPUTED FOR */
/*             A STEP OF LENGTH V(LMAX0) (BY A RETURN WITH IV(IRC) = 6). */
/* V(NREDUC) (I/O)  FUNCTION REDUCTION PREDICTED BY QUADRATIC MODEL FOR */
/*             NEWTON STEP.  IF ASSESS IS CALLED WITH IV(IRC) = 6, I.E., */
/*             IF V(PREDUC) HAS BEEN COMPUTED WITH RADIUS = V(LMAX0) FOR */
/*             USE IN THE SINGULAR CONVERVENCE TEST, THEN V(NREDUC) IS */
/*             SET TO -V(PREDUC) BEFORE THE LATTER IS RESTORED. */
/* V(PLSTGD) (I/O) VALUE OF V(PREDUC) ON SAVED STEP. */
/* V(PREDUC) (I/O) FUNCTION REDUCTION PREDICTED BY QUADRATIC MODEL FOR */
/*             CURRENT STEP. */
/* V(RADFAC) (OUT) FACTOR TO BE USED IN DETERMINING THE NEW RADIUS, */
/*             WHICH SHOULD BE V(RADFAC)*DST, WHERE  DST  IS EITHER THE */
/*             OUTPUT VALUE OF V(DSTNRM) OR THE 2-NORM OF */
/*             DIAG(NEWD)*STEP  FOR THE OUTPUT VALUE OF STEP AND THE */
/*             UPDATED VERSION, NEWD, OF THE SCALE VECTOR D.  FOR */
/*             IV(IRC) = 3, V(RADFAC) = 1.0 IS RETURNED. */
/* V(RDFCMN) (IN)  MINIMUM VALUE FOR V(RADFAC) IN TERMS OF THE INPUT */
/*             VALUE OF V(DSTNRM) -- SUGGESTED VALUE = 0.1. */
/* V(RDFCMX) (IN)  MAXIMUM VALUE FOR V(RADFAC) -- SUGGESTED VALUE = 4.0. */
/*  V(RELDX) (OUT) SCALED RELATIVE CHANGE IN X CAUSED BY STEP, COMPUTED */
/*             BY FUNCTION  RELDST  AS */
/*                 MAX (D(I)*ABS(X(I)-X0(I)), 1 .LE. I .LE. P) / */
/*                    MAX (D(I)*(ABS(X(I))+ABS(X0(I))), 1 .LE. I .LE. P). */
/*             IF AN ACCEPTABLE STEP IS RETURNED, THEN V(RELDX) IS COM- */
/*             PUTED USING THE OUTPUT (POSSIBLY RESTORED) VALUES OF X */
/*             AND STEP.  OTHERWISE IT IS COMPUTED USING THE INPUT */
/*             VALUES. */
/* V(RFCTOL) (IN)  RELATIVE FUNCTION CONVERGENCE TOLERANCE.  IF THE */
/*             ACTUAL FUNCTION REDUCTION IS AT MOST TWICE WHAT WAS PRE- */
/*             DICTED AND  V(NREDUC) .LE. V(RFCTOL)*ABS(V(F0)),  THEN */
/*             ASSESS RETURNS WITH IV(IRC) = 8 OR 9.  SEE ALSO V(LMAX0). */
/* V(STPPAR) (IN)  MARQUARDT PARAMETER -- 0 MEANS FULL NEWTON STEP. */
/* V(TUNER1) (IN)  TUNING CONSTANT USED TO DECIDE IF THE FUNCTION */
/*             REDUCTION WAS MUCH LESS THAN EXPECTED.  SUGGESTED */
/*             VALUE = 0.1. */
/* V(TUNER2) (IN)  TUNING CONSTANT USED TO DECIDE IF THE FUNCTION */
/*             REDUCTION WAS LARGE ENOUGH TO ACCEPT STEP.  SUGGESTED */
/*             VALUE = 10**-4. */
/* V(TUNER3) (IN)  TUNING CONSTANT USED TO DECIDE IF THE RADIUS */
/*             SHOULD BE INCREASED.  SUGGESTED VALUE = 0.75. */
/*  V(XCTOL) (IN)  X-CONVERGENCE CRITERION.  IF STEP IS A NEWTON STEP */
/*             (V(STPPAR) = 0) HAVING V(RELDX) .LE. V(XCTOL) AND GIVING */
/*             AT MOST TWICE THE PREDICTED FUNCTION DECREASE, THEN */
/*             ASSESS RETURNS IV(IRC) = 7 OR 9. */
/*  V(XFTOL) (IN)  FALSE CONVERGENCE TOLERANCE.  IF STEP GAVE NO OR ONLY */
/*             A SMALL FUNCTION DECREASE AND V(RELDX) .LE. V(XFTOL), */
/*             THEN ASSESS RETURNS WITH IV(IRC) = 12. */

/* -------------------------------  NOTES  ------------------------------- */

/*  ***  APPLICATION AND USAGE RESTRICTIONS  *** */

/*        THIS ROUTINE IS CALLED AS PART OF THE NL2SOL (NONLINEAR */
/*     LEAST-SQUARES) PACKAGE.  IT MAY BE USED IN ANY UNCONSTRAINED */
/*     MINIMIZATION SOLVER THAT USES DOGLEG, GOLDFELD-QUANDT-TROTTER, */
/*     OR LEVENBERG-MARQUARDT STEPS. */

/*  ***  ALGORITHM NOTES  *** */

/*        SEE (1) FOR FURTHER DISCUSSION OF THE ASSESSING AND MODEL */
/*     SWITCHING STRATEGIES.  WHILE NL2SOL CONSIDERS ONLY TWO MODELS, */
/*     ASSESS IS DESIGNED TO HANDLE ANY NUMBER OF MODELS. */

/*  ***  USAGE NOTES  *** */

/*        ON THE FIRST CALL OF AN ITERATION, ONLY THE I/O VARIABLES */
/*     STEP, X, IV(IRC), IV(MODEL), V(F), V(DSTNRM), V(GTSTEP), AND */
/*     V(PREDUC) NEED HAVE BEEN INITIALIZED.  BETWEEN CALLS, NO I/O */
/*     VALUES EXECPT STEP, X, IV(MODEL), V(F) AND THE STOPPING TOLER- */
/*     ANCES SHOULD BE CHANGED. */
/*        AFTER A RETURN FOR CONVERGENCE OR FALSE CONVERGENCE, ONE CAN */
/*     CHANGE THE STOPPING TOLERANCES AND CALL ASSESS AGAIN, IN WHICH */
/*     CASE THE STOPPING TESTS WILL BE REPEATED. */

/*  ***  REFERENCES  *** */

/*     (1) DENNIS, J.E., JR., GAY, D.M., AND WELSCH, R.E. (1981), */
/*        AN ADAPTIVE NONLINEAR LEAST-SQUARES ALGORITHM, */
/*        ACM TRANS. MATH. SOFTWARE, VOL. 7, NO. 3. */

/*     (2) POWELL, M.J.D. (1970)  A FORTRAN SUBROUTINE FOR SOLVING */
/*        SYSTEMS OF NONLINEAR ALGEBRAIC EQUATIONS, IN NUMERICAL */
/*        METHODS FOR NONLINEAR ALGEBRAIC EQUATIONS, EDITED BY */
/*        P. RABINOWITZ, GORDON AND BREACH, LONDON. */

/*  ***  HISTORY  *** */

/*        JOHN DENNIS DESIGNED MUCH OF THIS ROUTINE, STARTING WITH */
/*     IDEAS IN (2). ROY WELSCH SUGGESTED THE MODEL SWITCHING STRATEGY. */
/*        DAVID GAY AND STEPHEN PETERS CAST THIS SUBROUTINE INTO A MORE */
/*     PORTABLE FORM (WINTER 1977), AND DAVID GAY CAST IT INTO ITS */
/*     PRESENT FORM (FALL 1978). */

/*  ***  GENERAL  *** */

/*     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH */
/*     SUPPORTED BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS */
/*     MCS-7600324, DCR75-10143, 76-14311DSS, MCS76-11989, AND */
/*     MCS-7906671. */

/* ------------------------  EXTERNAL QUANTITIES  ------------------------ */

/*  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  *** */


/* VCOPY.... COPIES ONE VECTOR TO ANOTHER. */

/*  ***  INTRINSIC FUNCTIONS  *** */
/* /+ */
/* / */
/*  ***  NO COMMON BLOCKS  *** */

/* --------------------------  LOCAL VARIABLES  -------------------------- */


/*  ***  SUBSCRIPTS FOR IV AND V  *** */


/*  ***  DATA INITIALIZATIONS  *** */

/* /6 */
/*     DATA HALF/0.5D+0/, ONE/1.D+0/, TWO/2.D+0/, ZERO/0.D+0/ */
/* /7 */
/* / */

/* /6 */
/*     DATA IRC/3/, MLSTGD/4/, MODEL/5/, NFCALL/6/, */
/*    1     NFGCAL/7/, RADINC/8/, RESTOR/9/, STAGE/10/, */
/*    2     STGLIM/11/, SWITCH/12/, TOOBIG/2/, XIRC/13/ */
/* /7 */
/* / */
/* /6 */
/*     DATA AFCTOL/31/, DECFAC/22/, DSTNRM/2/, DST0/3/, */
/*    1     DSTSAV/18/, F/10/, FDIF/11/, FLSTGD/12/, F0/13/, */
/*    2     GTSLST/14/, GTSTEP/4/, INCFAC/23/, */
/*    3     LMAX0/35/, NREDUC/6/, PLSTGD/15/, PREDUC/7/, */
/*    4     RADFAC/16/, RDFCMN/24/, RDFCMX/25/, */
/*    5     RELDX/17/, RFCTOL/32/, STPPAR/5/, TUNER1/26/, */
/*    6     TUNER2/27/, TUNER3/28/, XCTOL/33/, XFTOL/34/ */
/* /7 */
/* / */

/* +++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++ */

    /* Parameter adjustments */
    --iv;
    --x0;
    --x;
    --stlstg;
    --step;
    --d__;
    --v;

    /* Function Body */
    nfc = iv[6];
    iv[12] = 0;
    iv[9] = 0;
    rfac1 = 1.;
    goodx = TRUE_;
    i__ = iv[3];
    if (i__ >= 1 && i__ <= 12) {
	switch (i__) {
	    case 1:  goto L20;
	    case 2:  goto L30;
	    case 3:  goto L10;
	    case 4:  goto L10;
	    case 5:  goto L40;
	    case 6:  goto L360;
	    case 7:  goto L290;
	    case 8:  goto L290;
	    case 9:  goto L290;
	    case 10:  goto L290;
	    case 11:  goto L290;
	    case 12:  goto L140;
	}
    }
    iv[3] = 13;
    goto L999;

/*  ***  INITIALIZE FOR NEW ITERATION  *** */

L10:
    iv[10] = 1;
    iv[8] = 0;
    v[12] = v[13];
    if (iv[2] == 0) {
	goto L90;
    }
    iv[10] = -1;
    iv[13] = i__;
    goto L60;

/*  ***  STEP WAS RECOMPUTED WITH NEW MODEL OR SMALLER RADIUS  *** */
/*  ***  FIRST DECIDE WHICH  *** */

L20:
    if (iv[5] != iv[4]) {
	goto L30;
    }
/*        ***  OLD MODEL RETAINED, SMALLER RADIUS TRIED  *** */
/*        ***  DO NOT CONSIDER ANY MORE NEW MODELS THIS ITERATION  *** */
    iv[10] = iv[11];
    iv[8] = -1;
    goto L90;

/*  ***  A NEW MODEL IS BEING TRIED.  DECIDE WHETHER TO KEEP IT.  *** */

L30:
    ++iv[10];

/*     ***  NOW WE ADD THE POSSIBILTIY THAT STEP WAS RECOMPUTED WITH  *** */
/*     ***  THE SAME MODEL, PERHAPS BECAUSE OF AN OVERSIZED STEP.     *** */

L40:
    if (iv[10] > 0) {
	goto L50;
    }

/*        ***  STEP WAS RECOMPUTED BECAUSE IT WAS TOO BIG.  *** */

    if (iv[2] != 0) {
	goto L60;
    }

/*        ***  RESTORE IV(STAGE) AND PICK UP WHERE WE LEFT OFF.  *** */

    iv[10] = -iv[10];
    i__ = iv[13];
    switch (i__) {
	case 1:  goto L20;
	case 2:  goto L30;
	case 3:  goto L90;
	case 4:  goto L90;
	case 5:  goto L70;
    }

L50:
    if (iv[2] == 0) {
	goto L70;
    }

/*  ***  HANDLE OVERSIZE STEP  *** */

    if (iv[8] > 0) {
	goto L80;
    }
    iv[10] = -iv[10];
    iv[13] = iv[3];

L60:
    v[16] = v[22];
    --iv[8];
    iv[3] = 5;
    goto L999;

L70:
    if (v[10] < v[12]) {
	goto L90;
    }

/*     *** THE NEW STEP IS A LOSER.  RESTORE OLD MODEL.  *** */

    if (iv[5] == iv[4]) {
	goto L80;
    }
    iv[5] = iv[4];
    iv[12] = 1;

/*     ***  RESTORE STEP, ETC. ONLY IF A PREVIOUS STEP DECREASED V(F). */

L80:
    if (v[12] >= v[13]) {
	goto L90;
    }
    iv[9] = 1;
    v[10] = v[12];
    v[7] = v[15];
    v[4] = v[14];
    if (iv[12] == 0) {
	rfac1 = v[2] / v[18];
    }
    v[2] = v[18];
    nfc = iv[7];
    goodx = FALSE_;


/*  ***  COMPUTE RELATIVE CHANGE IN X BY CURRENT STEP  *** */

L90:
    reldx1 = reldst_(p, &d__[1], &x[1], &x0[1]);

/*  ***  RESTORE X AND STEP IF NECESSARY  *** */

    if (goodx) {
	goto L105;
    }
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	step[i__] = stlstg[i__];
	x[i__] = x0[i__] + stlstg[i__];
/* L100: */
    }

L105:
    v[11] = v[13] - v[10];
    if (v[11] > v[27] * v[7]) {
	goto L120;
    }

/*        ***  NO (OR ONLY A TRIVIAL) FUNCTION DECREASE */
/*        ***  -- SO TRY NEW MODEL OR SMALLER RADIUS */

    v[17] = reldx1;
    if (v[10] < v[13]) {
	goto L110;
    }
    iv[4] = iv[5];
    v[12] = v[10];
    v[10] = v[13];
    vcopy_(p, &x[1], &x0[1]);
    iv[9] = 1;
    goto L115;
L110:
    iv[7] = nfc;
L115:
    iv[3] = 1;
    if (iv[10] < iv[11]) {
	goto L130;
    }
    iv[3] = 5;
    --iv[8];
    goto L130;

/*  ***  NONTRIVIAL FUNCTION DECREASE ACHIEVED  *** */

L120:
    iv[7] = nfc;
    rfac1 = 1.;
    if (goodx) {
	v[17] = reldx1;
    }
    v[18] = v[2];
    if (v[11] > v[7] * v[26]) {
	goto L200;
    }

/*  ***  DECREASE WAS MUCH LESS THAN PREDICTED -- EITHER CHANGE MODELS */
/*  ***  OR ACCEPT STEP WITH DECREASED RADIUS. */

    if (iv[10] >= iv[11]) {
	goto L125;
    }
/*        ***  CONSIDER SWITCHING MODELS  *** */
    iv[3] = 2;
    goto L130;

/*     ***  ACCEPT STEP WITH DECREASED RADIUS  *** */

L125:
    iv[3] = 4;

/*  ***  SET V(RADFAC) TO FLETCHER*S DECREASE FACTOR  *** */

L130:
    iv[13] = iv[3];
    emax = v[4] + v[11];
    v[16] = rfac1 * .5;
    if (emax < v[4]) {
/* Computing MAX */
	d__1 = v[24], d__2 = v[4] * .5 / emax;
	v[16] = rfac1 * max(d__1,d__2);
    }

/*  ***  DO FALSE CONVERGENCE TEST  *** */

L140:
    if (v[17] <= v[34]) {
	goto L160;
    }
    iv[3] = iv[13];
    if (v[10] < v[13]) {
	goto L230;
    }
    goto L300;

L160:
    iv[3] = 12;
    goto L310;

/*  ***  HANDLE GOOD FUNCTION DECREASE  *** */

L200:
    if (v[11] < -v[28] * v[4]) {
	goto L260;
    }

/*     ***  INCREASING RADIUS LOOKS WORTHWHILE.  SEE IF WE JUST */
/*     ***  RECOMPUTED STEP WITH A DECREASED RADIUS OR RESTORED STEP */
/*     ***  AFTER RECOMPUTING IT WITH A LARGER RADIUS. */

    if (iv[8] < 0) {
	goto L260;
    }
    if (iv[9] == 1) {
	goto L260;
    }

/*        ***  WE DID NOT.  TRY A LONGER STEP UNLESS THIS WAS A NEWTON */
/*        ***  STEP. */

    v[16] = v[25];
    gts = v[4];
    if (v[11] < (.5 / v[16] - 1.) * gts) {
/* Computing MAX */
	d__1 = v[23], d__2 = gts * .5 / (gts + v[11]);
	v[16] = max(d__1,d__2);
    }
    iv[3] = 4;
    if (v[5] == 0.) {
	goto L300;
    }
/*             ***  STEP WAS NOT A NEWTON STEP.  RECOMPUTE IT WITH */
/*             ***  A LARGER RADIUS. */
    iv[3] = 5;
    ++iv[8];

/*  ***  SAVE VALUES CORRESPONDING TO GOOD STEP  *** */

L230:
    v[12] = v[10];
    iv[4] = iv[5];
    vcopy_(p, &stlstg[1], &step[1]);
    v[18] = v[2];
    iv[7] = nfc;
    v[15] = v[7];
    v[14] = v[4];
    goto L300;

/*  ***  ACCEPT STEP WITH RADIUS UNCHANGED  *** */

L260:
    v[16] = 1.;
    iv[3] = 3;
    goto L300;

/*  ***  COME HERE FOR A RESTART AFTER CONVERGENCE  *** */

L290:
    iv[3] = iv[13];
    if (v[18] >= 0.) {
	goto L310;
    }
    iv[3] = 12;
    goto L310;

/*  ***  PERFORM CONVERGENCE TESTS  *** */

L300:
    iv[13] = iv[3];
L310:
    if (abs(v[10]) < v[31]) {
	iv[3] = 10;
    }
    if (v[11] * .5 > v[7]) {
	goto L999;
    }
    emax = v[32] * abs(v[13]);
    if (v[2] > v[35] && v[7] <= emax) {
	iv[3] = 11;
    }
    if (v[3] < 0.) {
	goto L320;
    }
    i__ = 0;
    if (v[6] > 0. && v[6] <= emax || v[6] == 0. && v[7] == 0.) {
	i__ = 2;
    }
    if (v[5] == 0. && v[17] <= v[33] && goodx) {
	++i__;
    }
    if (i__ > 0) {
	iv[3] = i__ + 6;
    }

/*  ***  CONSIDER RECOMPUTING STEP OF LENGTH V(LMAX0) FOR SINGULAR */
/*  ***  CONVERGENCE TEST. */

L320:
    if ((i__1 = iv[3] - 3, abs(i__1)) > 2 && iv[3] != 12) {
	goto L999;
    }
    if (v[2] > v[35]) {
	goto L330;
    }
    if (v[7] >= emax) {
	goto L999;
    }
    if (v[3] <= 0.) {
	goto L340;
    }
    if (v[3] * .5 <= v[35]) {
	goto L999;
    }
    goto L340;
L330:
    if (v[2] * .5 <= v[35]) {
	goto L999;
    }
    xmax = v[35] / v[2];
    if (xmax * (2. - xmax) * v[7] >= emax) {
	goto L999;
    }
L340:
    if (v[6] < 0.) {
	goto L370;
    }

/*  ***  RECOMPUTE V(PREDUC) FOR USE IN SINGULAR CONVERGENCE TEST  *** */

    v[14] = v[4];
    v[18] = v[2];
    if (iv[3] == 12) {
	v[18] = -v[18];
    }
    v[15] = v[7];
    iv[3] = 6;
    vcopy_(p, &stlstg[1], &step[1]);
    goto L999;

/*  ***  PERFORM SINGULAR CONVERGENCE TEST WITH RECOMPUTED V(PREDUC)  *** */

L360:
    v[4] = v[14];
    v[2] = abs(v[18]);
    vcopy_(p, &step[1], &stlstg[1]);
    iv[3] = iv[13];
    if (v[18] <= 0.) {
	iv[3] = 12;
    }
    v[6] = -v[7];
    v[7] = v[15];
L370:
    if (-v[6] <= v[32] * abs(v[13])) {
	iv[3] = 11;
    }

L999:
    return 0;

/*  ***  LAST CARD OF ASSESS FOLLOWS  *** */
} /* assess_ */

/* Subroutine */ int covclc_(integer *covirc, doublereal *d__, integer *iv, 
	doublereal *j, integer *n, integer *nn, integer *p, doublereal *r__, 
	doublereal *v, doublereal *x)
{
    /* System generated locals */
    integer j_dim1, j_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, k, l, m;
    static doublereal t;
    static integer g1, w0, w1, hc, gp, kl;
    static doublereal wk;
    static integer wl, rd1, ip1, mm1;
    static doublereal del;
    static integer hmi, irc, hpi, hpm, cov, stp0, qtr1, kind, mm1o2, pp1o2, 
	    stpi, stpm;
    static logical havej;
    static integer ipivi, ipivk;
    extern /* Subroutine */ int vcopy_(integer *, doublereal *, doublereal *),
	     lsqrt_(integer *, integer *, doublereal *, doublereal *, integer 
	    *);
    static integer gsave1;
    extern /* Subroutine */ int qrfact_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    doublereal *), ltsqar_(integer *, doublereal *, doublereal *), 
	    livmul_(integer *, doublereal *, doublereal *, doublereal *), 
	    linvrt_(integer *, doublereal *, doublereal *), litvmu_(integer *,
	     doublereal *, doublereal *, doublereal *), vscopy_(integer *, 
	    doublereal *, doublereal *);


/*  ***  COMPUTE COVARIANCE MATRIX FOR NL2ITR (NL2SOL VERSION 2.2)  *** */

/*  ***  LET K = IABS(IV(COVREQ).  FOR K .LE. 2, A FINITE-DIFFERENCE */
/*  ***  HESSIAN H IS COMPUTED (USING FUNC. AND GRAD. VALUES IF */
/*  ***  IV(COVREQ) IS NONNEGATIVE, AND USING ONLY FUNC. VALUES IF */
/*  ***  IV(COVREQ) IS NEGATIVE).  FOR SCALE = 2*F(X) / MAX(1, N-P), */
/*  ***  WHERE 2*F(X) IS THE RESIDUAL SUM OF SQUARES, COVCLC COMPUTES... */
/*  ***             K = 0 OR 1...  SCALE * H**-1 * (J**T * J) * H**-1. */
/*  ***             K = 2...  SCALE * H**-1. */
/*  ***             K .GE. 3...  SCALE * (J**T * J)**-1. */

/*  ***  PARAMETER DECLARATIONS  *** */

/*     DIMENSION IV(*), V(*) */

/*  ***  LOCAL VARIABLES  *** */


/*  ***  INTRINSIC FUNCTIONS  *** */
/* /+ */
/* / */
/*  ***  EXTERNAL SUBROUTINES  *** */


/* LINVRT... INVERT LOWER TRIANGULAR MATRIX. */
/* LITVMU... APPLY INVERSE-TRANSPOSE OF COMPACT LOWER TRIANG. MATRIX. */
/* LIVMUL... APPLY INVERSE OF COMPACT LOWER TRIANG. MATRIX. */
/* LSQRT.... COMPUTE CHOLESKY FACTOR OF (LOWER TRINAG. OF) A SYM. MATRIX. */
/* LTSQAR... GIVEN LOWER TRIANG. MATRIX L, COMPUTE (L**T)*L. */
/* QRFACT... COMPUTE QR DECOMPOSITION OF A MATRIX. */
/* VCOPY.... COPY ONE VECTOR TO ANOTHER. */
/* VSCOPY... SET ALL ELEMENTS OF A VECTOR TO A SCALAR. */

/*  ***  SUBSCRIPTS FOR IV AND V  *** */


/* /6 */
/*     DATA HALF/0.5D+0/, NEGPT5/-0.5D+0/, ONE/1.D+0/, TWO/2.D+0/, */
/*    1     ZERO/0.D+0/ */
/* /7 */
/* / */

/* /6 */
/*     DATA COVMAT/26/, COVREQ/15/, DELTA/50/, DELTA0/44/, */
/*    1     DLTFDC/40/, F/10/, FX/46/, G/28/, H/44/, IERR/32/, */
/*    2     IPIVOT/61/, IPIV0/60/, KAGQT/35/, KALM/36/, */
/*    3     LMAT/58/, MODE/38/, NFGCAL/7/, QTR/49/, */
/*    4     RD/51/, RSAVE/52/, SAVEI/54/, SWITCH/12/, */
/*    5     TOOBIG/2/, W/59/, XMSAVE/49/ */
/* /7 */
/* / */

/* +++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++ */

    /* Parameter adjustments */
    --iv;
    --r__;
    --x;
    j_dim1 = *nn;
    j_offset = 1 + j_dim1;
    j -= j_offset;
    --d__;
    --v;

    /* Function Body */
    *covirc = 4;
    kind = iv[15];
    m = iv[38];
    if (m > 0) {
	goto L10;
    }
    iv[35] = -1;
    if (iv[36] > 0) {
	iv[36] = 0;
    }
    if (abs(kind) >= 3) {
	goto L300;
    }
    v[46] = v[10];
    k = iv[52];
    vcopy_(n, &v[k], &r__[1]);
L10:
    if (m > *p) {
	goto L200;
    }
    if (kind < 0) {
	goto L100;
    }

/*  ***  COMPUTE FINITE-DIFFERENCE HESSIAN USING BOTH FUNCTION AND */
/*  ***  GRADIENT VALUES. */

    gsave1 = iv[59] + *p;
    g1 = iv[28];
    if (m > 0) {
	goto L15;
    }
/*        ***  FIRST CALL ON COVCLC.  SET GSAVE = G, TAKE FIRST STEP  *** */
    vcopy_(p, &v[gsave1], &v[g1]);
    iv[12] = iv[7];
    goto L80;

L15:
    del = v[50];
    x[m] = v[49];
    if (iv[2] == 0) {
	goto L30;
    }

/*     ***  HANDLE OVERSIZE V(DELTA)  *** */

    if (del * x[m] > 0.) {
	goto L20;
    }
/*             ***  WE ALREADY TRIED SHRINKING V(DELTA), SO QUIT  *** */
    iv[26] = -2;
    goto L190;

/*        ***  TRY SHRINKING V(DELTA)  *** */
L20:
    del *= -.5;
    goto L90;

L30:
    cov = iv[58];
    gp = g1 + *p - 1;

/*  ***  SET  G = (G - GSAVE)/DEL  *** */

    i__1 = gp;
    for (i__ = g1; i__ <= i__1; ++i__) {
	v[i__] = (v[i__] - v[gsave1]) / del;
	++gsave1;
/* L40: */
    }

/*  ***  ADD G AS NEW COL. TO FINITE-DIFF. HESSIAN MATRIX  *** */

    k = cov + m * (m - 1) / 2;
    l = k + m - 2;
    if (m == 1) {
	goto L60;
    }

/*  ***  SET  H(I,M) = 0.5 * (H(I,M) + G(I))  FOR I = 1 TO M-1  *** */

    i__1 = l;
    for (i__ = k; i__ <= i__1; ++i__) {
	v[i__] = (v[i__] + v[g1]) * .5;
	++g1;
/* L50: */
    }

/*  ***  ADD  H(I,M) = G(I)  FOR I = M TO P  *** */

L60:
    ++l;
    i__1 = *p;
    for (i__ = m; i__ <= i__1; ++i__) {
	v[l] = v[g1];
	l += i__;
	++g1;
/* L70: */
    }

L80:
    ++m;
    iv[38] = m;
    if (m > *p) {
	goto L190;
    }

/*  ***  CHOOSE NEXT FINITE-DIFFERENCE STEP, RETURN TO GET G THERE  *** */

/* Computing MAX */
    d__2 = 1. / d__[m], d__3 = (d__1 = x[m], abs(d__1));
    del = v[44] * max(d__2,d__3);
    if (x[m] < 0.) {
	del = -del;
    }
    v[49] = x[m];
L90:
    x[m] += del;
    v[50] = del;
    *covirc = 2;
    goto L999;

/*  ***  COMPUTE FINITE-DIFFERENCE HESSIAN USING FUNCTION VALUES ONLY. */

L100:
    stp0 = iv[59] + *p - 1;
    mm1 = m - 1;
    mm1o2 = m * mm1 / 2;
    if (m > 0) {
	goto L105;
    }
/*        ***  FIRST CALL ON COVCLC.  *** */
    iv[54] = 0;
    goto L180;

L105:
    i__ = iv[54];
    if (i__ > 0) {
	goto L160;
    }
    if (iv[2] == 0) {
	goto L120;
    }

/*     ***  HANDLE OVERSIZE STEP  *** */

    stpm = stp0 + m;
    del = v[stpm];
    if (del * x[49] > 0.) {
	goto L110;
    }
/*             ***  WE ALREADY TRIED SHRINKING THE STEP, SO QUIT  *** */
    iv[26] = -2;
    goto L999;

/*        ***  TRY SHRINKING THE STEP  *** */
L110:
    del *= -.5;
    x[m] = x[49] + del;
    v[stpm] = del;
    *covirc = 1;
    goto L999;

/*  ***  SAVE F(X + STP(M)*E(M)) IN H(P,M)  *** */

L120:
    pp1o2 = *p * (*p - 1) / 2;
    cov = iv[58];
    hpm = cov + pp1o2 + mm1;
    v[hpm] = v[10];

/*  ***  START COMPUTING ROW M OF THE FINITE-DIFFERENCE HESSIAN H.  *** */

    hmi = cov + mm1o2;
    if (mm1 == 0) {
	goto L140;
    }
    hpi = cov + pp1o2;
    i__1 = mm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[hmi] = v[46] - (v[10] + v[hpi]);
	++hmi;
	++hpi;
/* L130: */
    }
L140:
    v[hmi] = v[10] - v[46] * 2.;

/*  ***  COMPUTE FUNCTION VALUES NEEDED TO COMPLETE ROW M OF H.  *** */

    i__ = 1;

L150:
    iv[54] = i__;
    stpi = stp0 + i__;
    v[50] = x[i__];
    x[i__] += v[stpi];
    if (i__ == m) {
	x[i__] = v[49] - v[stpi];
    }
    *covirc = 1;
    goto L999;

L160:
    x[i__] = v[50];
    if (iv[2] == 0) {
	goto L170;
    }
/*        ***  PUNT IN THE EVENT OF AN OVERSIZE STEP  *** */
    iv[26] = -2;
    goto L999;

/*  ***  FINISH COMPUTING H(M,I)  *** */

L170:
    stpi = stp0 + i__;
    hmi = cov + mm1o2 + i__ - 1;
    stpm = stp0 + m;
    v[hmi] = (v[hmi] + v[10]) / (v[stpi] * v[stpm]);
    ++i__;
    if (i__ <= m) {
	goto L150;
    }
    iv[54] = 0;
    x[m] = v[49];

L180:
    ++m;
    iv[38] = m;
    if (m > *p) {
	goto L190;
    }

/*  ***  PREPARE TO COMPUTE ROW M OF THE FINITE-DIFFERENCE HESSIAN H. */
/*  ***  COMPUTE M-TH STEP SIZE STP(M), THEN RETURN TO OBTAIN */
/*  ***  F(X + STP(M)*E(M)), WHERE E(M) = M-TH STD. UNIT VECTOR. */

/* Computing MAX */
    d__2 = 1. / d__[m], d__3 = (d__1 = x[m], abs(d__1));
    del = v[40] * max(d__2,d__3);
    if (x[m] < 0.) {
	del = -del;
    }
    v[49] = x[m];
    x[m] += del;
    stpm = stp0 + m;
    v[stpm] = del;
    *covirc = 1;
    goto L999;

/*  ***  RESTORE R, V(F), ETC.  *** */

L190:
    k = iv[52];
    vcopy_(n, &r__[1], &v[k]);
    v[10] = v[46];
    if (kind < 0) {
	goto L200;
    }
    iv[7] = iv[12];
    qtr1 = iv[49];
    vcopy_(n, &v[qtr1], &r__[1]);
    if (iv[26] < 0) {
	goto L999;
    }
    *covirc = 3;
    goto L999;

L200:
    cov = iv[58];

/*  ***  THE COMPLETE FINITE-DIFF. HESSIAN IS NOW STORED AT V(COV).   *** */
/*  ***  USE IT TO COMPUTE THE REQUESTED COVARIANCE MATRIX.           *** */

/*     ***  COMPUTE CHOLESKY FACTOR C OF H = C*(C**T)  *** */
/*     ***  AND STORE IT AT V(HC).  *** */

    hc = cov;
    if (abs(kind) == 2) {
	goto L210;
    }
    hc = abs(iv[44]);
    iv[44] = -hc;
L210:
    lsqrt_(&c__1, p, &v[hc], &v[cov], &irc);
    iv[26] = -1;
    if (irc != 0) {
	goto L999;
    }

    w1 = iv[59] + *p;
    if (abs(kind) > 1) {
	goto L350;
    }

/*  ***  COVARIANCE = SCALE * H**-1 * (J**T * J) * H**-1  *** */

    i__1 = *p * (*p + 1) / 2;
    vscopy_(&i__1, &v[cov], &c_b95);
    havej = iv[36] == -1;
/*     ***  HAVEJ = .TRUE. MEANS J IS IN ITS ORIGINAL FORM, WHILE */
/*     ***  HAVEJ = .FALSE. MEANS QRFACT HAS BEEN APPLIED TO J. */

    m = *p;
    if (havej) {
	m = *n;
    }
    w0 = w1 - 1;
    rd1 = iv[51];
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (havej) {
	    goto L240;
	}

/*        ***  SET W = IPIVOT * (ROW I OF R MATRIX FROM QRFACT).  *** */

	vscopy_(p, &v[w1], &c_b95);
	ipivi = i__ + 60;
	l = w0 + iv[ipivi];
	v[l] = v[rd1];
	++rd1;
	if (i__ == *p) {
	    goto L260;
	}
	ip1 = i__ + 1;
	i__2 = *p;
	for (k = ip1; k <= i__2; ++k) {
	    ipivk = k + 60;
	    l = w0 + iv[ipivk];
	    v[l] = j[i__ + k * j_dim1];
/* L230: */
	}
	goto L260;

/*        ***  SET W = (ROW I OF J).  *** */

L240:
	l = w0;
	i__2 = *p;
	for (k = 1; k <= i__2; ++k) {
	    ++l;
	    v[l] = j[i__ + k * j_dim1];
/* L250: */
	}

/*        ***  SET W = H**-1 * W.  *** */

L260:
	livmul_(p, &v[w1], &v[hc], &v[w1]);
	litvmu_(p, &v[w1], &v[hc], &v[w1]);

/*        ***  ADD  W * W**T  TO COVARIANCE MATRIX.  *** */

	kl = cov;
	i__2 = *p;
	for (k = 1; k <= i__2; ++k) {
	    l = w0 + k;
	    wk = v[l];
	    i__3 = k;
	    for (l = 1; l <= i__3; ++l) {
		wl = w0 + l;
		v[kl] += wk * v[wl];
		++kl;
/* L270: */
	    }
/* L280: */
	}
/* L290: */
    }
    goto L380;

/*  ***  COVARIANCE = SCALE * (J**T * J)**-1.  *** */

L300:
    rd1 = iv[51];
    if (iv[36] != -1) {
	goto L310;
    }

/*        ***  APPLY QRFACT TO J  *** */

    qtr1 = iv[49];
    vcopy_(n, &v[qtr1], &r__[1]);
    w1 = iv[59] + *p;
    qrfact_(nn, n, p, &j[j_offset], &v[rd1], &iv[61], &iv[32], &c__0, &v[w1]);
    iv[36] = -2;
L310:
    iv[26] = -1;
    if (iv[32] != 0) {
	goto L999;
    }
    cov = iv[58];
    hc = abs(iv[44]);
    iv[44] = -hc;

/*     ***  SET HC = (R MATRIX FROM QRFACT).  *** */

    l = hc;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ > 1) {
	    i__2 = i__ - 1;
	    vcopy_(&i__2, &v[l], &j[i__ * j_dim1 + 1]);
	}
	l = l + i__ - 1;
	v[l] = v[rd1];
	++l;
	++rd1;
/* L340: */
    }

/*  ***  THE CHOLESKY FACTOR C OF THE UNSCALED INVERSE COVARIANCE MATRIX */
/*  ***  (OR PERMUTATION THEREOF) IS STORED AT V(HC). */

/*  ***  SET C = C**-1. */

L350:
    linvrt_(p, &v[hc], &v[hc]);

/*  ***  SET C = C**T * C. */

    ltsqar_(p, &v[hc], &v[hc]);

    if (hc == cov) {
	goto L380;
    }

/*     ***  C = PERMUTED, UNSCALED COVARIANCE. */
/*     ***  SET COV = IPIVOT * C * IPIVOT**T. */

    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	m = i__ + 60;
	ipivi = iv[m];
	kl = cov - 1 + ipivi * (ipivi - 1) / 2;
	i__2 = i__;
	for (k = 1; k <= i__2; ++k) {
	    m = k + 60;
	    ipivk = iv[m];
	    l = kl + ipivk;
	    if (ipivk > ipivi) {
		l += (ipivk - ipivi) * (ipivk + ipivi - 3) / 2;
	    }
	    v[l] = v[hc];
	    ++hc;
/* L360: */
	}
/* L370: */
    }

L380:
    iv[26] = cov;

/*  ***  APPLY SCALE FACTOR = (RESID. SUM OF SQUARES) / MAX(1,N-P). */

/* Computing MAX */
    i__1 = 1, i__2 = *n - *p;
    t = v[10] / ((real) max(i__1,i__2) * .5);
    k = cov - 1 + *p * (*p + 1) / 2;
    i__1 = k;
    for (i__ = cov; i__ <= i__1; ++i__) {
/* L390: */
	v[i__] = t * v[i__];
    }

L999:
    return 0;
/*  ***  LAST CARD OF COVCLC FOLLOWS  *** */
} /* covclc_ */

/* Subroutine */ int dfault_(integer *iv, doublereal *v)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static doublereal machep;
    extern integer imdcon_(integer *);
    extern doublereal rmdcon_(integer *);
    static doublereal mepcrt, sqteps;


/*  ***  SUPPLY NL2SOL (VERSION 2.2) DEFAULT VALUES TO IV AND V  *** */

/* /+ */
/* / */


/*  ***  SUBSCRIPTS FOR IV AND V  *** */


/* /6 */
/*     DATA ONE/1.D+0/, THREE/3.D+0/ */
/* /7 */
/* / */

/*  ***  IV SUBSCRIPT VALUES  *** */

/* /6 */
/*     DATA COVPRT/14/, COVREQ/15/, DTYPE/16/, INITS/25/, */
/*    1     MXFCAL/17/, MXITER/18/, OUTLEV/19/, */
/*    2     PARPRT/20/, PRUNIT/21/, SOLPRT/22/, */
/*    3     STATPR/23/, X0PRT/24/ */
/* /7 */
/* / */

/*  ***  V SUBSCRIPT VALUES  *** */

/* /6 */
/*     DATA AFCTOL/31/, COSMIN/43/, DECFAC/22/, */
/*    1     DELTA0/44/, DFAC/41/, DINIT/38/, DLTFDC/40/, */
/*    2     DLTFDJ/36/, D0INIT/37/, EPSLON/19/, FUZZ/45/, */
/*    3     INCFAC/23/, JTINIT/39/, LMAX0/35/, PHMNFC/20/, */
/*    4     PHMXFC/21/, RDFCMN/24/, RDFCMX/25/, */
/*    5     RFCTOL/32/, RLIMIT/42/, TUNER1/26/, */
/*    6     TUNER2/27/, TUNER3/28/, TUNER4/29/, */
/*    7     TUNER5/30/, XCTOL/33/, XFTOL/34/ */
/* /7 */
/* / */

/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --v;
    --iv;

    /* Function Body */
    iv[1] = 12;
    iv[14] = 1;
    iv[15] = 1;
    iv[16] = 1;
    iv[25] = 0;
    iv[17] = 200;
    iv[18] = 150;
    iv[19] = 1;
    iv[20] = 1;
    iv[21] = imdcon_(&c__1);
    iv[22] = 1;
    iv[23] = 1;
    iv[24] = 1;

    machep = rmdcon_(&c__3);
    v[31] = 1e-20;
    if (machep > 1e-10) {
/* Computing 2nd power */
	d__1 = machep;
	v[31] = d__1 * d__1;
    }
/* Computing MAX */
    d__1 = 1e-6, d__2 = machep * 100.;
    v[43] = max(d__1,d__2);
    v[22] = .5;
    sqteps = rmdcon_(&c__4);
    v[44] = sqteps;
    v[41] = .6;
    v[38] = 0.;
    mepcrt = pow_dd(&machep, &c_b235);
    v[40] = mepcrt;
    v[36] = sqteps;
    v[37] = 1.;
    v[19] = .1;
    v[45] = 1.5;
    v[23] = 2.;
    v[39] = 1e-6;
    v[35] = 100.;
    v[20] = -.1;
    v[21] = .1;
    v[24] = .1;
    v[25] = 4.;
/* Computing MAX */
/* Computing 2nd power */
    d__3 = mepcrt;
    d__1 = 1e-10, d__2 = d__3 * d__3;
    v[32] = max(d__1,d__2);
    v[42] = rmdcon_(&c__5);
    v[26] = .1;
    v[27] = 1e-4;
    v[28] = .75;
    v[29] = .5;
    v[30] = .75;
    v[33] = sqteps;
    v[34] = machep * 100.;

/* L999: */
    return 0;
/*  ***  LAST CARD OF DFAULT FOLLOWS  *** */
} /* dfault_ */

doublereal dotprd_(integer *p, doublereal *x, doublereal *y)
{
    /* Initialized data */

    static doublereal sqteta = 0.;

    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2, d__3, d__4;

    /* Local variables */
    static integer i__;
    static doublereal t;
    extern doublereal rmdcon_(integer *);


/*  ***  RETURN THE INNER PRODUCT OF THE P-VECTORS X AND Y.  *** */


/* /+ */
/* / */

/*  ***  RMDCON(2) RETURNS A MACHINE-DEPENDENT CONSTANT, SQTETA, WHICH */
/*  ***  IS SLIGHTLY LARGER THAN THE SMALLEST POSITIVE NUMBER THAT */
/*  ***  CAN BE SQUARED WITHOUT UNDERFLOWING. */

/* /6 */
/*     DATA ONE/1.D+0/, SQTETA/0.D+0/, ZERO/0.D+0/ */
/* /7 */
    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
/* / */

    ret_val = 0.;
    if (*p <= 0) {
	goto L999;
    }
    if (sqteta == 0.) {
	sqteta = rmdcon_(&c__2);
    }
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	d__3 = (d__1 = x[i__], abs(d__1)), d__4 = (d__2 = y[i__], abs(d__2));
	t = max(d__3,d__4);
	if (t > 1.) {
	    goto L10;
	}
	if (t < sqteta) {
	    goto L20;
	}
	t = x[i__] / sqteta * y[i__];
	if (abs(t) < sqteta) {
	    goto L20;
	}
L10:
	ret_val += x[i__] * y[i__];
L20:
	;
    }

L999:
    return ret_val;
/*  ***  LAST CARD OF DOTPRD FOLLOWS  *** */
} /* dotprd_ */

/* Subroutine */ int dupdat_(doublereal *d__, integer *iv, doublereal *j, 
	integer *n, integer *nn, integer *p, doublereal *v)
{
    /* System generated locals */
    integer j_dim1, j_offset, i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal t;
    static integer d0, s1;
    static doublereal sii, vdfac;
    static integer jtoli;
    extern doublereal v2norm_(integer *, doublereal *);


/*  ***  UPDATE SCALE VECTOR D FOR NL2ITR (NL2SOL VERSION 2.2)  *** */

/*  ***  PARAMETER DECLARATIONS  *** */

/*     DIMENSION IV(*), V(*) */

/*  ***  LOCAL VARIABLES  *** */


/*     ***  CONSTANTS  *** */

/*  ***  INTRINSIC FUNCTIONS  *** */
/* /+ */
/* / */
/*  ***  EXTERNAL FUNCTION  *** */


/*  ***  SUBSCRIPTS FOR IV AND V  *** */

/* /6 */
/*     DATA DFAC/41/, DTYPE/16/, JTOL0/86/, NITER/31/, S/53/ */
/* /7 */
/* / */

/* /6 */
/*     DATA ZERO/0.D+0/ */
/* /7 */
/* / */

/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --iv;
    j_dim1 = *nn;
    j_offset = 1 + j_dim1;
    j -= j_offset;
    --d__;
    --v;

    /* Function Body */
    i__ = iv[16];
    if (i__ == 1) {
	goto L20;
    }
    if (iv[31] > 0) {
	goto L999;
    }

L20:
    vdfac = v[41];
    d0 = *p + 86;
    s1 = iv[53] - 1;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s1 += i__;
	sii = v[s1];
	t = v2norm_(n, &j[i__ * j_dim1 + 1]);
	if (sii > 0.) {
	    t = sqrt(t * t + sii);
	}
	jtoli = i__ + 86;
	++d0;
	if (t < v[jtoli]) {
/* Computing MAX */
	    d__1 = v[d0], d__2 = v[jtoli];
	    t = max(d__1,d__2);
	}
/* Computing MAX */
	d__1 = vdfac * d__[i__];
	d__[i__] = max(d__1,t);
/* L30: */
    }

L999:
    return 0;
/*  ***  LAST CARD OF DUPDAT FOLLOWS  *** */
} /* dupdat_ */

/* Subroutine */ int gqtstp_(doublereal *d__, doublereal *dig, doublereal *
	dihdi, integer *ka, doublereal *l, integer *p, doublereal *step, 
	doublereal *v, doublereal *w)
{
    /* Initialized data */

    static doublereal dgxfac = 0.;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, q;
    static doublereal t;
    static integer x, k1, q0;
    static doublereal t1;
    static integer x0;
    static doublereal lk, si, sk, uk, wi, sw;
    static integer im1, lk0, uk0;
    static doublereal aki, akk, rad;
    static integer inc, irc;
    static doublereal phi, dst;
    static integer diag, emin, emax;
    static doublereal root;
    static integer diag0;
    static doublereal epso6, delta;
    static integer kalim;
    extern /* Subroutine */ int lsqrt_(integer *, integer *, doublereal *, 
	    doublereal *, integer *);
    extern doublereal v2norm_(integer *, doublereal *);
    static doublereal alphak, psifac;
    static integer dggdmx;
    static doublereal oldphi;
    extern doublereal rmdcon_(integer *);
    static doublereal phimin, phimax;
    static integer phipin;
    extern doublereal dotprd_(integer *, doublereal *, doublereal *);
    static integer dstsav;
    extern /* Subroutine */ int livmul_(integer *, doublereal *, doublereal *,
	     doublereal *);
    extern doublereal lsvmin_(integer *, doublereal *, doublereal *, 
	    doublereal *);
    extern /* Subroutine */ int litvmu_(integer *, doublereal *, doublereal *,
	     doublereal *);
    static logical restrt;
    static doublereal twopsi;


/*  *** COMPUTE GOLDFELD-QUANDT-TROTTER STEP BY MORE-HEBDEN TECHNIQUE *** */
/*  ***  (NL2SOL VERSION 2.2)  *** */

/*  ***  PARAMETER DECLARATIONS  *** */

/*     DIMENSION DIHDI(P*(P+1)/2), L(P*(P+1)/2), W(4*P+7) */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*  ***  PURPOSE  *** */

/*        GIVEN THE (COMPACTLY STORED) LOWER TRIANGLE OF A SCALED */
/*     HESSIAN (APPROXIMATION) AND A NONZERO SCALED GRADIENT VECTOR, */
/*     THIS SUBROUTINE COMPUTES A GOLDFELD-QUANDT-TROTTER STEP OF */
/*     APPROXIMATE LENGTH V(RADIUS) BY THE MORE-HEBDEN TECHNIQUE.  IN */
/*     OTHER WORDS, STEP IS COMPUTED TO (APPROXIMATELY) MINIMIZE */
/*     PSI(STEP) = (G**T)*STEP + 0.5*(STEP**T)*H*STEP  SUCH THAT THE */
/*     2-NORM OF D*STEP IS AT MOST (APPROXIMATELY) V(RADIUS), WHERE */
/*     G  IS THE GRADIENT,  H  IS THE HESSIAN, AND  D  IS A DIAGONAL */
/*     SCALE MATRIX WHOSE DIAGONAL IS STORED IN THE PARAMETER D. */
/*     (GQTSTP ASSUMES  DIG = D**-1 * G  AND  DIHDI = D**-1 * H * D**-1.) */
/*     IF G = 0, HOWEVER, STEP = 0 IS RETURNED (EVEN AT A SADDLE POINT). */

/*  ***  PARAMETER DESCRIPTION  *** */

/*     D (IN)  = THE SCALE VECTOR, I.E. THE DIAGONAL OF THE SCALE */
/*              MATRIX  D  MENTIONED ABOVE UNDER PURPOSE. */
/*   DIG (IN)  = THE SCALED GRADIENT VECTOR, D**-1 * G.  IF G = 0, THEN */
/*              STEP = 0  AND  V(STPPAR) = 0  ARE RETURNED. */
/* DIHDI (IN)  = LOWER TRIANGLE OF THE SCALED HESSIAN (APPROXIMATION), */
/*              I.E., D**-1 * H * D**-1, STORED COMPACTLY BY ROWS., I.E., */
/*              IN THE ORDER (1,1), (2,1), (2,2), (3,1), (3,2), ETC. */
/*    KA (I/O) = THE NUMBER OF HEBDEN ITERATIONS (SO FAR) TAKEN TO DETER- */
/*              MINE STEP.  KA .LT. 0 ON INPUT MEANS THIS IS THE FIRST */
/*              ATTEMPT TO DETERMINE STEP (FOR THE PRESENT DIG AND DIHDI) */
/*              -- KA IS INITIALIZED TO 0 IN THIS CASE.  OUTPUT WITH */
/*              KA = 0  (OR V(STPPAR) = 0)  MEANS  STEP = -(H**-1)*G. */
/*     L (I/O) = WORKSPACE OF LENGTH P*(P+1)/2 FOR CHOLESKY FACTORS. */
/*     P (IN)  = NUMBER OF PARAMETERS -- THE HESSIAN IS A  P X P  MATRIX. */
/*  STEP (I/O) = THE STEP COMPUTED. */
/*     V (I/O) CONTAINS VARIOUS CONSTANTS AND VARIABLES DESCRIBED BELOW. */
/*     W (I/O) = WORKSPACE OF LENGTH 4*P + 6. */

/*  ***  ENTRIES IN V  *** */

/* V(DGNORM) (I/O) = 2-NORM OF (D**-1)*G. */
/* V(DSTNRM) (OUTPUT) = 2-NORM OF D*STEP. */
/* V(DST0)   (I/O) = 2-NORM OF D*(H**-1)*G (FOR POS. DEF. H ONLY), OR */
/*             OVERESTIMATE OF SMALLEST EIGENVALUE OF (D**-1)*H*(D**-1). */
/* V(EPSLON) (IN)  = MAX. REL. ERROR ALLOWED FOR PSI(STEP).  FOR THE */
/*             STEP RETURNED, PSI(STEP) WILL EXCEED ITS OPTIMAL VALUE */
/*             BY LESS THAN -V(EPSLON)*PSI(STEP).  SUGGESTED VALUE = 0.1. */
/* V(GTSTEP) (OUT) = INNER PRODUCT BETWEEN G AND STEP. */
/* V(NREDUC) (OUT) = PSI(-(H**-1)*G) = PSI(NEWTON STEP)  (FOR POS. DEF. */
/*             H ONLY -- V(NREDUC) IS SET TO ZERO OTHERWISE). */
/* V(PHMNFC) (IN)  = TOL. (TOGETHER WITH V(PHMXFC)) FOR ACCEPTING STEP */
/*             (MORE*S SIGMA).  THE ERROR V(DSTNRM) - V(RADIUS) MUST LIE */
/*             BETWEEN V(PHMNFC)*V(RADIUS) AND V(PHMXFC)*V(RADIUS). */
/* V(PHMXFC) (IN)  (SEE V(PHMNFC).) */
/*             SUGGESTED VALUES -- V(PHMNFC) = -0.25, V(PHMXFC) = 0.5. */
/* V(PREDUC) (OUT) = PSI(STEP) = PREDICTED OBJ. FUNC. REDUCTION FOR STEP. */
/* V(RADIUS) (IN)  = RADIUS OF CURRENT (SCALED) TRUST REGION. */
/* V(RAD0)   (I/O) = VALUE OF V(RADIUS) FROM PREVIOUS CALL. */
/* V(STPPAR) (I/O) IS NORMALLY THE MARQUARDT PARAMETER, I.E. THE ALPHA */
/*             DESCRIBED BELOW UNDER ALGORITHM NOTES.  IF H + ALPHA*D**2 */
/*             (SEE ALGORITHM NOTES) IS (NEARLY) SINGULAR, HOWEVER, */
/*             THEN V(STPPAR) = -ALPHA. */

/*  ***  USAGE NOTES  *** */

/*     IF IT IS DESIRED TO RECOMPUTE STEP USING A DIFFERENT VALUE OF */
/*     V(RADIUS), THEN THIS ROUTINE MAY BE RESTARTED BY CALLING IT */
/*     WITH ALL PARAMETERS UNCHANGED EXCEPT V(RADIUS).  (THIS EXPLAINS */
/*     WHY STEP AND W ARE LISTED AS I/O).  ON AN INTIIAL CALL (ONE WITH */
/*     KA .LT. 0), STEP AND W NEED NOT BE INITIALIZED AND ONLY COMPO- */
/*     NENTS V(EPSLON), V(STPPAR), V(PHMNFC), V(PHMXFC), V(RADIUS), AND */
/*     V(RAD0) OF V MUST BE INITIALIZED.  TO COMPUTE STEP FROM A SADDLE */
/*     POINT (WHERE THE TRUE GRADIENT VANISHES AND H HAS A NEGATIVE */
/*     EIGENVALUE), A NONZERO G WITH SMALL COMPONENTS SHOULD BE PASSED. */

/*  ***  APPLICATION AND USAGE RESTRICTIONS  *** */

/*     THIS ROUTINE IS CALLED AS PART OF THE NL2SOL (NONLINEAR LEAST- */
/*     SQUARES) PACKAGE (REF. 1), BUT IT COULD BE USED IN SOLVING ANY */
/*     UNCONSTRAINED MINIMIZATION PROBLEM. */

/*  ***  ALGORITHM NOTES  *** */

/*        THE DESIRED G-Q-T STEP (REF. 2, 3, 4) SATISFIES */
/*     (H + ALPHA*D**2)*STEP = -G  FOR SOME NONNEGATIVE ALPHA SUCH THAT */
/*     H + ALPHA*D**2 IS POSITIVE SEMIDEFINITE.  ALPHA AND STEP ARE */
/*     COMPUTED BY A SCHEME ANALOGOUS TO THE ONE DESCRIBED IN REF. 5. */
/*     ESTIMATES OF THE SMALLEST AND LARGEST EIGENVALUES OF THE HESSIAN */
/*     ARE OBTAINED FROM THE GERSCHGORIN CIRCLE THEOREM ENHANCED BY A */
/*     SIMPLE FORM OF THE SCALING DESCRIBED IN REF. 6.  CASES IN WHICH */
/*     H + ALPHA*D**2 IS NEARLY (OR EXACTLY) SINGULAR ARE HANDLED BY */
/*     THE TECHNIQUE DISCUSSED IN REF. 2.  IN THESE CASES, A STEP OF */
/*     (EXACT) LENGTH V(RADIUS) IS RETURNED FOR WHICH PSI(STEP) EXCEEDS */
/*     ITS OPTIMAL VALUE BY LESS THAN -V(EPSLON)*PSI(STEP). */

/*  ***  FUNCTIONS AND SUBROUTINES CALLED  *** */

/* DOTPRD - RETURNS INNER PRODUCT OF TWO VECTORS. */
/* LITVMU - APPLIES INVERSE-TRANSPOSE OF COMPACT LOWER TRIANG. MATRIX. */
/* LIVMUL - APPLIES INVERSE OF COMPACT LOWER TRIANG. MATRIX. */
/* LSQRT  - FINDS CHOLESKY FACTOR (OF COMPACTLY STORED LOWER TRIANG.). */
/* LSVMIN - RETURNS APPROX. TO MIN. SING. VALUE OF LOWER TRIANG. MATRIX. */
/* RMDCON - RETURNS MACHINE-DEPENDENT CONSTANTS. */
/* V2NORM - RETURNS 2-NORM OF A VECTOR. */

/*  ***  REFERENCES  *** */

/* 1.  DENNIS, J.E., GAY, D.M., AND WELSCH, R.E. (1981), AN ADAPTIVE */
/*             NONLINEAR LEAST-SQUARES ALGORITHM, ACM TRANS. MATH. */
/*             SOFTWARE, VOL. 7, NO. 3. */
/* 2.  GAY, D.M. (1981), COMPUTING OPTIMAL LOCALLY CONSTRAINED STEPS, */
/*             SIAM J. SCI. STATIST. COMPUTING, VOL. 2, NO. 2, PP. */
/*             186-197. */
/* 3.  GOLDFELD, S.M., QUANDT, R.E., AND TROTTER, H.F. (1966), */
/*             MAXIMIZATION BY QUADRATIC HILL-CLIMBING, ECONOMETRICA 34, */
/*             PP. 541-551. */
/* 4.  HEBDEN, M.D. (1973), AN ALGORITHM FOR MINIMIZATION USING EXACT */
/*             SECOND DERIVATIVES, REPORT T.P. 515, THEORETICAL PHYSICS */
/*             DIV., A.E.R.E. HARWELL, OXON., ENGLAND. */
/* 5.  MORE, J.J. (1978), THE LEVENBERG-MARQUARDT ALGORITHM, IMPLEMEN- */
/*             TATION AND THEORY, PP.105-116 OF SPRINGER LECTURE NOTES */
/*             IN MATHEMATICS NO. 630, EDITED BY G.A. WATSON, SPRINGER- */
/*             VERLAG, BERLIN AND NEW YORK. */
/* 6.  VARGA, R.S. (1965), MINIMAL GERSCHGORIN SETS, PACIFIC J. MATH. 15, */
/*             PP. 719-729. */

/*  ***  GENERAL  *** */

/*     CODED BY DAVID M. GAY. */
/*     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH */
/*     SUPPORTED BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS */
/*     MCS-7600324, DCR75-10143, 76-14311DSS, MCS76-11989, AND */
/*     MCS-7906671. */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*  ***  LOCAL VARIABLES  *** */


/*     ***  CONSTANTS  *** */

/*  ***  INTRINSIC FUNCTIONS  *** */
/* /+ */
/* / */
/*  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  *** */


/*  ***  SUBSCRIPTS FOR V  *** */

/* /6 */
/*     DATA DGNORM/1/, DSTNRM/2/, DST0/3/, EPSLON/19/, */
/*    1     GTSTEP/4/, NREDUC/6/, PHMNFC/20/, */
/*    2     PHMXFC/21/, PREDUC/7/, RADIUS/8/, */
/*    3     RAD0/9/, STPPAR/5/ */
/* /7 */
/* / */

/* /6 */
/*     DATA EPSFAC/50.0D+0/, FOUR/4.0D+0/, HALF/0.5D+0/, */
/*    1     KAPPA/2.0D+0/, NEGONE/-1.0D+0/, ONE/1.0D+0/, P001/1.0D-3/, */
/*    2     SIX/6.0D+0/, THREE/3.0D+0/, TWO/2.0D+0/, ZERO/0.0D+0/ */
/* /7 */
/* / */
    /* Parameter adjustments */
    --dihdi;
    --l;
    --step;
    --dig;
    --d__;
    --v;
    --w;

    /* Function Body */

/*  ***  BODY  *** */

/*     ***  STORE LARGEST ABS. ENTRY IN (D**-1)*H*(D**-1) AT W(DGGDMX). */
    dggdmx = *p + 1;
/*     ***  STORE GERSCHGORIN OVER- AND UNDERESTIMATES OF THE LARGEST */
/*     ***  AND SMALLEST EIGENVALUES OF (D**-1)*H*(D**-1) AT W(EMAX) */
/*     ***  AND W(EMIN) RESPECTIVELY. */
    emax = dggdmx + 1;
    emin = emax + 1;
/*     ***  FOR USE IN RECOMPUTING STEP, THE FINAL VALUES OF LK, UK, DST, */
/*     ***  AND THE INVERSE DERIVATIVE OF MORE*S PHI AT 0 (FOR POS. DEF. */
/*     ***  H) ARE STORED IN W(LK0), W(UK0), W(DSTSAV), AND W(PHIPIN) */
/*     ***  RESPECTIVELY. */
    lk0 = emin + 1;
    phipin = lk0 + 1;
    uk0 = phipin + 1;
    dstsav = uk0 + 1;
/*     ***  STORE DIAG OF (D**-1)*H*(D**-1) IN W(DIAG),...,W(DIAG0+P). */
    diag0 = dstsav;
    diag = diag0 + 1;
/*     ***  STORE -D*STEP IN W(Q),...,W(Q0+P). */
    q0 = diag0 + *p;
    q = q0 + 1;
    rad = v[8];
/*     ***  PHITOL = MAX. ERROR ALLOWED IN DST = V(DSTNRM) = 2-NORM OF */
/*     ***  D*STEP. */
    phimax = v[21] * rad;
    phimin = v[20] * rad;
/*     ***  EPSO6 AND PSIFAC ARE USED IN CHECKING FOR THE SPECIAL CASE */
/*     ***  OF (NEARLY) SINGULAR H + ALPHA*D**2 (SEE REF. 2). */
/* Computing 2nd power */
    d__1 = rad;
    psifac = v[19] * 2. / (((v[20] + 1.) * 4. * 3. + 2. + 2.) * 3. * (d__1 * 
	    d__1));
/*     ***  OLDPHI IS USED TO DETECT LIMITS OF NUMERICAL ACCURACY.  IF */
/*     ***  WE RECOMPUTE STEP AND IT DOES NOT CHANGE, THEN WE ACCEPT IT. */
    oldphi = 0.;
    epso6 = v[19] / 6.;
    irc = 0;
    restrt = FALSE_;
    kalim = *ka + 50;

/*  ***  START OR RESTART, DEPENDING ON KA  *** */

    if (*ka >= 0) {
	goto L310;
    }

/*  ***  FRESH START  *** */

    k = 0;
    uk = -1.;
    *ka = 0;
    kalim = 50;

/*     ***  STORE DIAG(DIHDI) IN W(DIAG0+1),...,W(DIAG0+P)  *** */

    j = 0;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j += i__;
	k1 = diag0 + i__;
	w[k1] = dihdi[j];
/* L20: */
    }

/*     ***  DETERMINE W(DGGDMX), THE LARGEST ELEMENT OF DIHDI  *** */

    t1 = 0.;
    j = *p * (*p + 1) / 2;
    i__1 = j;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t = (d__1 = dihdi[i__], abs(d__1));
	if (t1 < t) {
	    t1 = t;
	}
/* L30: */
    }
    w[dggdmx] = t1;

/*  ***  TRY ALPHA = 0  *** */

L40:
    lsqrt_(&c__1, p, &l[1], &dihdi[1], &irc);
    if (irc == 0) {
	goto L60;
    }
/*        ***  INDEF. H -- UNDERESTIMATE SMALLEST EIGENVALUE, USE THIS */
/*        ***  ESTIMATE TO INITIALIZE LOWER BOUND LK ON ALPHA. */
    j = irc * (irc + 1) / 2;
    t = l[j];
    l[j] = 1.;
    i__1 = irc;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L50: */
	w[i__] = 0.;
    }
    w[irc] = 1.;
    litvmu_(&irc, &w[1], &l[1], &w[1]);
    t1 = v2norm_(&irc, &w[1]);
    lk = -t / t1 / t1;
    v[3] = -lk;
    if (restrt) {
	goto L210;
    }
    v[6] = 0.;
    goto L70;

/*     ***  POSITIVE DEFINITE H -- COMPUTE UNMODIFIED NEWTON STEP.  *** */
L60:
    lk = 0.;
    livmul_(p, &w[q], &l[1], &dig[1]);
    v[6] = dotprd_(p, &w[q], &w[q]) * .5;
    litvmu_(p, &w[q], &l[1], &w[q]);
    dst = v2norm_(p, &w[q]);
    v[3] = dst;
    phi = dst - rad;
    if (phi <= phimax) {
	goto L280;
    }
    if (restrt) {
	goto L210;
    }

/*  ***  PREPARE TO COMPUTE GERSCHGORIN ESTIMATES OF LARGEST (AND */
/*  ***  SMALLEST) EIGENVALUES.  *** */

L70:
    v[1] = v2norm_(p, &dig[1]);
    if (v[1] == 0.) {
	goto L450;
    }
    k = 0;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wi = 0.;
	if (i__ == 1) {
	    goto L90;
	}
	im1 = i__ - 1;
	i__2 = im1;
	for (j = 1; j <= i__2; ++j) {
	    ++k;
	    t = (d__1 = dihdi[k], abs(d__1));
	    wi += t;
	    w[j] += t;
/* L80: */
	}
L90:
	w[i__] = wi;
	++k;
/* L100: */
    }

/*  ***  (UNDER-)ESTIMATE SMALLEST EIGENVALUE OF (D**-1)*H*(D**-1)  *** */

    k = 1;
    t1 = w[diag] - w[1];
    if (*p <= 1) {
	goto L120;
    }
    i__1 = *p;
    for (i__ = 2; i__ <= i__1; ++i__) {
	j = diag0 + i__;
	t = w[j] - w[i__];
	if (t >= t1) {
	    goto L110;
	}
	t1 = t;
	k = i__;
L110:
	;
    }

L120:
    sk = w[k];
    j = diag0 + k;
    akk = w[j];
    k1 = k * (k - 1) / 2 + 1;
    inc = 1;
    t = 0.;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ == k) {
	    goto L130;
	}
	aki = (d__1 = dihdi[k1], abs(d__1));
	si = w[i__];
	j = diag0 + i__;
	t1 = (akk - w[j] + si - aki) * .5;
	t1 += sqrt(t1 * t1 + sk * aki);
	if (t < t1) {
	    t = t1;
	}
	if (i__ < k) {
	    goto L140;
	}
L130:
	inc = i__;
L140:
	k1 += inc;
/* L150: */
    }

    w[emin] = akk - t;
    uk = v[1] / rad - w[emin];

/*  ***  COMPUTE GERSCHGORIN (OVER-)ESTIMATE OF LARGEST EIGENVALUE  *** */

    k = 1;
    t1 = w[diag] + w[1];
    if (*p <= 1) {
	goto L170;
    }
    i__1 = *p;
    for (i__ = 2; i__ <= i__1; ++i__) {
	j = diag0 + i__;
	t = w[j] + w[i__];
	if (t <= t1) {
	    goto L160;
	}
	t1 = t;
	k = i__;
L160:
	;
    }

L170:
    sk = w[k];
    j = diag0 + k;
    akk = w[j];
    k1 = k * (k - 1) / 2 + 1;
    inc = 1;
    t = 0.;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ == k) {
	    goto L180;
	}
	aki = (d__1 = dihdi[k1], abs(d__1));
	si = w[i__];
	j = diag0 + i__;
	t1 = (w[j] + si - aki - akk) * .5;
	t1 += sqrt(t1 * t1 + sk * aki);
	if (t < t1) {
	    t = t1;
	}
	if (i__ < k) {
	    goto L190;
	}
L180:
	inc = i__;
L190:
	k1 += inc;
/* L200: */
    }

    w[emax] = akk + t;
/* Computing MAX */
    d__1 = lk, d__2 = v[1] / rad - w[emax];
    lk = max(d__1,d__2);

/*     ***  ALPHAK = CURRENT VALUE OF ALPHA (SEE ALG. NOTES ABOVE).  WE */
/*     ***  USE MORE*S SCHEME FOR INITIALIZING IT. */
    alphak = abs(v[5]) * v[9] / rad;

    if (irc != 0) {
	goto L210;
    }

/*  ***  COMPUTE L0 FOR POSITIVE DEFINITE H  *** */

    livmul_(p, &w[1], &l[1], &w[q]);
    t = v2norm_(p, &w[1]);
    w[phipin] = dst / t / t;
/* Computing MAX */
    d__1 = lk, d__2 = phi * w[phipin];
    lk = max(d__1,d__2);

/*  ***  SAFEGUARD ALPHAK AND ADD ALPHAK*I TO (D**-1)*H*(D**-1)  *** */

L210:
    ++(*ka);
    if (-v[3] >= alphak || alphak < lk || alphak >= uk) {
/* Computing MAX */
	d__1 = .001, d__2 = sqrt(lk / uk);
	alphak = uk * max(d__1,d__2);
    }
    k = 0;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k += i__;
	j = diag0 + i__;
	dihdi[k] = w[j] + alphak;
/* L220: */
    }

/*  ***  TRY COMPUTING CHOLESKY DECOMPOSITION  *** */

    lsqrt_(&c__1, p, &l[1], &dihdi[1], &irc);
    if (irc == 0) {
	goto L250;
    }

/*  ***  (D**-1)*H*(D**-1) + ALPHAK*I  IS INDEFINITE -- OVERESTIMATE */
/*  ***  SMALLEST EIGENVALUE FOR USE IN UPDATING LK  *** */

    j = irc * (irc + 1) / 2;
    t = l[j];
    l[j] = 1.;
    i__1 = irc;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L230: */
	w[i__] = 0.;
    }
    w[irc] = 1.;
    litvmu_(&irc, &w[1], &l[1], &w[1]);
    t1 = v2norm_(&irc, &w[1]);
    lk = alphak - t / t1 / t1;
    v[3] = -lk;
    goto L210;

/*  ***  ALPHAK MAKES (D**-1)*H*(D**-1) POSITIVE DEFINITE. */
/*  ***  COMPUTE Q = -D*STEP, CHECK FOR CONVERGENCE.  *** */

L250:
    livmul_(p, &w[q], &l[1], &dig[1]);
    litvmu_(p, &w[q], &l[1], &w[q]);
    dst = v2norm_(p, &w[q]);
    phi = dst - rad;
    if (phi <= phimax && phi >= phimin) {
	goto L290;
    }
    if (phi == oldphi) {
	goto L290;
    }
    oldphi = phi;
    if (phi > 0.) {
	goto L260;
    }
/*        ***  CHECK FOR THE SPECIAL CASE OF  H + ALPHA*D**2  (NEARLY) */
/*        ***  SINGULAR.  DELTA IS .GE. THE SMALLEST EIGENVALUE OF */
/*        ***  (D**-1)*H*(D**-1) + ALPHAK*I. */
    if (v[3] > 0.) {
	goto L260;
    }
    delta = alphak + v[3];
    twopsi = alphak * dst * dst + dotprd_(p, &dig[1], &w[q]);
    if (delta < psifac * twopsi) {
	goto L270;
    }

/*  ***  UNACCEPTABLE ALPHAK -- UPDATE LK, UK, ALPHAK  *** */

L260:
    if (*ka >= kalim) {
	goto L290;
    }
    livmul_(p, &w[1], &l[1], &w[q]);
    t1 = v2norm_(p, &w[1]);
/*     ***  THE FOLLOWING DMIN1 IS NECESSARY BECAUSE OF RESTARTS  *** */
    if (phi < 0.) {
	uk = min(uk,alphak);
    }
    alphak += phi / t1 * (dst / t1) * (dst / rad);
    lk = max(lk,alphak);
    goto L210;

/*  ***  DECIDE HOW TO HANDLE (NEARLY) SINGULAR H + ALPHA*D**2  *** */

/*     ***  IF NOT YET AVAILABLE, OBTAIN MACHINE DEPENDENT VALUE DGXFAC. */
L270:
    if (dgxfac == 0.) {
	dgxfac = rmdcon_(&c__3) * 50.;
    }

/*     ***  NOW DECIDE.  *** */
    if (delta > dgxfac * w[dggdmx]) {
	goto L350;
    }
/*        ***  DELTA IS SO SMALL WE CANNOT HANDLE THE SPECIAL CASE IN */
/*        ***  THE AVAILABLE ARITHMETIC.  ACCEPT STEP AS IT IS. */
    goto L290;

/*  ***  ACCEPTABLE STEP ON FIRST TRY  *** */

L280:
    alphak = 0.;

/*  ***  SUCCESSFUL STEP IN GENERAL.  COMPUTE STEP = -(D**-1)*Q  *** */

L290:
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = q0 + i__;
	step[i__] = -w[j] / d__[i__];
/* L300: */
    }
    v[4] = -dotprd_(p, &dig[1], &w[q]);
    v[7] = (abs(alphak) * dst * dst - v[4]) * .5;
    goto L430;


/*  ***  RESTART WITH NEW RADIUS  *** */

L310:
    if (v[3] <= 0. || v[3] - rad > phimax) {
	goto L330;
    }

/*     ***  PREPARE TO RETURN NEWTON STEP  *** */

    restrt = TRUE_;
    ++(*ka);
    k = 0;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k += i__;
	j = diag0 + i__;
	dihdi[k] = w[j];
/* L320: */
    }
    uk = -1.;
    goto L40;

L330:
    if (*ka == 0) {
	goto L60;
    }

    dst = w[dstsav];
    alphak = abs(v[5]);
    phi = dst - rad;
    t = v[1] / rad;
    if (rad > v[9]) {
	goto L340;
    }

/*        ***  SMALLER RADIUS  *** */
    uk = t - w[emin];
    lk = 0.;
    if (alphak > 0.) {
	lk = w[lk0];
    }
/* Computing MAX */
    d__1 = lk, d__2 = t - w[emax];
    lk = max(d__1,d__2);
    if (v[3] > 0.) {
/* Computing MAX */
	d__1 = lk, d__2 = (v[3] - rad) * w[phipin];
	lk = max(d__1,d__2);
    }
    goto L260;

/*     ***  BIGGER RADIUS  *** */
L340:
    uk = t - w[emin];
    if (alphak > 0.) {
/* Computing MIN */
	d__1 = uk, d__2 = w[uk0];
	uk = min(d__1,d__2);
    }
/* Computing MAX */
    d__1 = 0., d__2 = -v[3], d__1 = max(d__1,d__2), d__2 = t - w[emax];
    lk = max(d__1,d__2);
    if (v[3] > 0.) {
/* Computing MAX */
	d__1 = lk, d__2 = (v[3] - rad) * w[phipin];
	lk = max(d__1,d__2);
    }
    goto L260;

/*  ***  HANDLE (NEARLY) SINGULAR H + ALPHA*D**2  *** */

/*     ***  NEGATE ALPHAK TO INDICATE SPECIAL CASE  *** */
L350:
    alphak = -alphak;
/*     ***  ALLOCATE STORAGE FOR SCRATCH VECTOR X  *** */
    x0 = q0 + *p;
    x = x0 + 1;

/*  ***  USE INVERSE POWER METHOD WITH START FROM LSVMIN TO OBTAIN */
/*  ***  APPROXIMATE EIGENVECTOR CORRESPONDING TO SMALLEST EIGENVALUE */
/*  ***  OF (D**-1)*H*(D**-1). */

    delta *= 2.;
    t = lsvmin_(p, &l[1], &w[x], &w[1]);

    k = 0;
/*     ***  NORMALIZE W  *** */
L360:
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L370: */
	w[i__] = t * w[i__];
    }
/*     ***  COMPLETE CURRENT INV. POWER ITER. -- REPLACE W BY (L**-T)*W. */
    litvmu_(p, &w[1], &l[1], &w[1]);
    t1 = 1. / v2norm_(p, &w[1]);
    t = t1 * t;
    if (t <= delta) {
	goto L390;
    }
    if (k > 30) {
	goto L290;
    }
    ++k;
/*     ***  START NEXT INV. POWER ITER. BY STORING NORMALIZED W IN X. */
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = x0 + i__;
	w[j] = t1 * w[i__];
/* L380: */
    }
/*     ***  COMPUTE W = (L**-1)*X. */
    livmul_(p, &w[1], &l[1], &w[x]);
    t = 1. / v2norm_(p, &w[1]);
    goto L360;

L390:
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L400: */
	w[i__] = t1 * w[i__];
    }

/*  ***  NOW W IS THE DESIRED APPROXIMATE (UNIT) EIGENVECTOR AND */
/*  ***  T*X = ((D**-1)*H*(D**-1) + ALPHAK*I)*W. */

    sw = dotprd_(p, &w[q], &w[1]);
    t1 = (rad + dst) * (rad - dst);
    root = sqrt(sw * sw + t1);
    if (sw < 0.) {
	root = -root;
    }
    si = t1 / (sw + root);
/*     ***  ACCEPT CURRENT STEP IF ADDING SI*W WOULD LEAD TO A */
/*     ***  FURTHER RELATIVE REDUCTION IN PSI OF LESS THAN V(EPSLON)/3. */
    v[7] = twopsi * .5;
    t1 = 0.;
    t = si * (alphak * sw - si * .5 * (alphak + t * dotprd_(p, &w[x], &w[1])))
	    ;
    if (t < epso6 * twopsi) {
	goto L410;
    }
    v[7] += t;
    dst = rad;
    t1 = -si;
L410:
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = q0 + i__;
	w[j] = t1 * w[i__] - w[j];
	step[i__] = w[j] / d__[i__];
/* L420: */
    }
    v[4] = dotprd_(p, &dig[1], &w[q]);

/*  ***  SAVE VALUES FOR USE IN A POSSIBLE RESTART  *** */

L430:
    v[2] = dst;
    v[5] = alphak;
    w[lk0] = lk;
    w[uk0] = uk;
    v[9] = rad;
    w[dstsav] = dst;

/*     ***  RESTORE DIAGONAL OF DIHDI  *** */

    j = 0;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j += i__;
	k = diag0 + i__;
	dihdi[j] = w[k];
/* L440: */
    }
    goto L999;

/*  ***  SPECIAL CASE -- G = 0  *** */

L450:
    v[5] = 0.;
    v[7] = 0.;
    v[2] = 0.;
    v[4] = 0.;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L460: */
	step[i__] = 0.;
    }

L999:
    return 0;

/*  ***  LAST CARD OF GQTSTP FOLLOWS  *** */
} /* gqtstp_ */

/* Subroutine */ int itsmry_(doublereal *d__, integer *iv, integer *p, 
	doublereal *v, doublereal *x)
{
    /* Initialized data */

    static char model1[4*6+1] = "                  G   S ";
    static char model2[4*6+1] = " G   S  G-S S-G -S-G-G-S";

    /* Format strings */
    static char fmt_1010[] = "(\0020   IT    NF\002,6x,\002F\002,8x,\002RE\
LDF\002,6x,\002PRELDF\002,5x,\002RELDX\002)";
    static char fmt_1017[] = "(1x,i5,i6,4d11.3,a3,a4,4d11.3)";
    static char fmt_1015[] = "(\0020   IT    NF\002,6x,\002F\002,8x,\002RE\
LDF\002,6x,\002PRELDF\002,5x,\002RELDX\002,4x,\002MODEL    STPPAR\002,6x,\
\002SIZE\002,6x,\002D*STEP\002,5x,\002NPRELDF\002)";
    static char fmt_1030[] = "(\0020***** X-CONVERGENCE *****\002)";
    static char fmt_1035[] = "(\0020***** RELATIVE FUNCTION CONVERGENCE **\
***\002)";
    static char fmt_1040[] = "(\0020***** X- AND RELATIVE FUNCTION CONVERGEN\
CE *****\002)";
    static char fmt_1045[] = "(\0020***** ABSOLUTE FUNCTION CONVERGENCE **\
***\002)";
    static char fmt_1050[] = "(\0020***** SINGULAR CONVERGENCE *****\002)";
    static char fmt_1060[] = "(\0020***** FALSE CONVERGENCE *****\002)";
    static char fmt_1070[] = "(\0020***** FUNCTION EVALUATION LIMIT *****\
\002)";
    static char fmt_1080[] = "(\0020***** ITERATION LIMIT *****\002)";
    static char fmt_1090[] = "(\0020***** STOPX *****\002)";
    static char fmt_1100[] = "(\0020***** INITIAL SUM OF SQUARES OVERFLOWS *\
****\002)";
    static char fmt_1120[] = "(\0020***** BAD PARAMETERS TO ASSESS *****\002)"
	    ;
    static char fmt_1130[] = "(\0020***** J COULD NOT BE COMPUTED *****\002)";
    static char fmt_1140[] = "(\0020***** IV(1) =\002,i5,\002 *****\002)";
    static char fmt_1150[] = "(\0020    I     INITIAL X(I)\002,7x,\002D(I\
)\002//(1x,i5,d17.6,d14.3))";
    static char fmt_1160[] = "(\0020    0     1\002,d11.3,11x,d11.3)";
    static char fmt_1180[] = "(\0020FUNCTION\002,d17.6,\002   RELDX\002,d20.\
6/\002 FUNC. EVALS\002,i8,9x,\002GRAD. EVALS\002,i8/\002 PRELDF\002,d19.6,3x,\
\002NPRELDF\002,d18.6)";
    static char fmt_1185[] = "(\0020\002,i4,\002 EXTRA FUNC. EVALS FOR COVAR\
IANCE.\002)";
    static char fmt_1186[] = "(1x,i4,\002 EXTRA GRAD. EVALS FOR COVARIANCE\
.\002)";
    static char fmt_1190[] = "(\0020    I      FINAL X(I)\002,8x,\002D(I)\
\002,10x,\002G(I)\002/)";
    static char fmt_1200[] = "(1x,i5,d17.6,2d14.3)";
    static char fmt_1220[] = "(\0020++++++ INDEFINITE COVARIANCE MATRIX ++++\
++\002)";
    static char fmt_1225[] = "(\0020++++++ OVERSIZE STEPS IN COMPUTING COVAR\
IANCE +++++\002)";
    static char fmt_1230[] = "(\0020++++++ COVARIANCE MATRIX NOT COMPUTED ++\
++++\002)";
    static char fmt_1241[] = "(\0020COVARIANCE = SCALE * H**-1 * (J**T * J) \
* H**-1\002/)";
    static char fmt_1242[] = "(\0020COVARIANCE = SCALE * H**-1\002/)";
    static char fmt_1243[] = "(\0020COVARIANCE = SCALE * (J**T * J)**-1\002/)"
	    ;
    static char fmt_1250[] = "(\002 ROW\002,i3,2x,9d12.4/(9x,9d12.4))";
    static char fmt_1270[] = "(\002 ROW\002,i3,2x,5d12.4/(9x,5d12.4))";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j, m, g1, i1, ii, nf, ng, ol, pu, iv1, cov1;
    static doublereal oldf, reldf, nreldf, preldf;

    /* Fortran I/O blocks */
    static cilist io___218 = { 0, 0, 0, fmt_1010, 0 };
    static cilist io___219 = { 0, 0, 0, fmt_1017, 0 };
    static cilist io___220 = { 0, 0, 0, fmt_1015, 0 };
    static cilist io___223 = { 0, 0, 0, fmt_1017, 0 };
    static cilist io___224 = { 0, 0, 0, fmt_1030, 0 };
    static cilist io___225 = { 0, 0, 0, fmt_1035, 0 };
    static cilist io___226 = { 0, 0, 0, fmt_1040, 0 };
    static cilist io___227 = { 0, 0, 0, fmt_1045, 0 };
    static cilist io___228 = { 0, 0, 0, fmt_1050, 0 };
    static cilist io___229 = { 0, 0, 0, fmt_1060, 0 };
    static cilist io___230 = { 0, 0, 0, fmt_1070, 0 };
    static cilist io___231 = { 0, 0, 0, fmt_1080, 0 };
    static cilist io___232 = { 0, 0, 0, fmt_1090, 0 };
    static cilist io___233 = { 0, 0, 0, fmt_1100, 0 };
    static cilist io___234 = { 0, 0, 0, fmt_1120, 0 };
    static cilist io___235 = { 0, 0, 0, fmt_1130, 0 };
    static cilist io___236 = { 0, 0, 0, fmt_1140, 0 };
    static cilist io___237 = { 0, 0, 0, fmt_1150, 0 };
    static cilist io___239 = { 0, 0, 0, fmt_1010, 0 };
    static cilist io___240 = { 0, 0, 0, fmt_1015, 0 };
    static cilist io___241 = { 0, 0, 0, fmt_1160, 0 };
    static cilist io___243 = { 0, 0, 0, fmt_1180, 0 };
    static cilist io___244 = { 0, 0, 0, fmt_1185, 0 };
    static cilist io___245 = { 0, 0, 0, fmt_1186, 0 };
    static cilist io___247 = { 0, 0, 0, fmt_1190, 0 };
    static cilist io___248 = { 0, 0, 0, fmt_1200, 0 };
    static cilist io___250 = { 0, 0, 0, fmt_1220, 0 };
    static cilist io___251 = { 0, 0, 0, fmt_1225, 0 };
    static cilist io___252 = { 0, 0, 0, fmt_1230, 0 };
    static cilist io___253 = { 0, 0, 0, fmt_1241, 0 };
    static cilist io___254 = { 0, 0, 0, fmt_1242, 0 };
    static cilist io___255 = { 0, 0, 0, fmt_1243, 0 };
    static cilist io___258 = { 0, 0, 0, fmt_1250, 0 };
    static cilist io___260 = { 0, 0, 0, fmt_1270, 0 };



/*  ***  PRINT NL2SOL (VERSION 2.2) ITERATION SUMMARY  *** */

/*  ***  PARAMETER DECLARATIONS  *** */

/*     DIMENSION IV(*), V(*) */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*  ***  LOCAL VARIABLES  *** */

/* /6 */
/*     REAL MODEL1(6), MODEL2(6) */
/* /7 */
/* / */

/*  ***  INTRINSIC FUNCTIONS  *** */
/* /+ */
/* / */
/*  ***  NO EXTERNAL FUNCTIONS OR SUBROUTINES  *** */

/*  ***  SUBSCRIPTS FOR IV AND V  *** */


/*  ***  IV SUBSCRIPT VALUES  *** */

/* /6 */
/*     DATA COVMAT/26/, COVPRT/14/, G/28/, COVREQ/15/, */
/*    1     NEEDHD/39/, NFCALL/6/, NFCOV/40/, NGCOV/41/, */
/*    2     NGCALL/30/, NITER/31/, OUTLEV/19/, PRNTIT/48/, */
/*    3     PRUNIT/21/, SOLPRT/22/, STATPR/23/, SUSED/57/, */
/*    4     X0PRT/24/ */
/* /7 */
/* / */

/*  ***  V SUBSCRIPT VALUES  *** */

/* /6 */
/*     DATA DSTNRM/2/, F/10/, F0/13/, FDIF/11/, NREDUC/6/, */
/*    1     PREDUC/7/, RELDX/17/, SIZE/47/, STPPAR/5/ */
/* /7 */
/* / */

/* /6 */
/*     DATA ZERO/0.D+0/ */
/* /7 */
/* / */
/* /6 */
/*     DATA MODEL1(1)/4H    /, MODEL1(2)/4H    /, MODEL1(3)/4H    /, */
/*    1     MODEL1(4)/4H    /, MODEL1(5)/4H  G /, MODEL1(6)/4H  S /, */
/*    2     MODEL2(1)/4H G  /, MODEL2(2)/4H S  /, MODEL2(3)/4HG-S /, */
/*    3     MODEL2(4)/4HS-G /, MODEL2(5)/4H-S-G/, MODEL2(6)/4H-G-S/ */
/* /7 */
    /* Parameter adjustments */
    --iv;
    --x;
    --d__;
    --v;

    /* Function Body */
/* / */

/* ----------------------------------------------------------------------- */

    pu = iv[21];
    if (pu == 0) {
	goto L999;
    }
    iv1 = iv[1];
    ol = iv[19];
    if (iv1 < 2 || iv1 > 15) {
	goto L140;
    }
    if (ol == 0) {
	goto L20;
    }
    if (iv1 >= 12) {
	goto L20;
    }
    if (iv1 >= 10 && iv[48] == 0) {
	goto L20;
    }
    if (iv1 > 2) {
	goto L10;
    }
    ++iv[48];
    if (iv[48] < abs(ol)) {
	goto L999;
    }
L10:
    nf = iv[6] - abs(iv[40]);
    iv[48] = 0;
    reldf = 0.;
    preldf = 0.;
    oldf = v[13];
    if (oldf <= 0.) {
	goto L12;
    }
    reldf = v[11] / oldf;
    preldf = v[7] / oldf;
L12:
    if (ol > 0) {
	goto L15;
    }

/*        ***  PRINT SHORT SUMMARY LINE  *** */

    if (iv[39] == 1) {
	io___218.ciunit = pu;
	s_wsfe(&io___218);
	e_wsfe();
    }
    iv[39] = 0;
    io___219.ciunit = pu;
    s_wsfe(&io___219);
    do_fio(&c__1, (char *)&iv[31], (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nf, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&v[10], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&reldf, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&preldf, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&v[17], (ftnlen)sizeof(doublereal));
    e_wsfe();
    goto L20;

/*     ***  PRINT LONG SUMMARY LINE  *** */

L15:
    if (iv[39] == 1) {
	io___220.ciunit = pu;
	s_wsfe(&io___220);
	e_wsfe();
    }
    iv[39] = 0;
    m = iv[57];
    nreldf = 0.;
    if (oldf > 0.) {
	nreldf = v[6] / oldf;
    }
    io___223.ciunit = pu;
    s_wsfe(&io___223);
    do_fio(&c__1, (char *)&iv[31], (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nf, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&v[10], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&reldf, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&preldf, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&v[17], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, model1 + (m - 1 << 2), (ftnlen)4);
    do_fio(&c__1, model2 + (m - 1 << 2), (ftnlen)4);
    do_fio(&c__1, (char *)&v[5], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&v[47], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&v[2], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&nreldf, (ftnlen)sizeof(doublereal));
    e_wsfe();

L20:
    switch (iv1) {
	case 1:  goto L999;
	case 2:  goto L999;
	case 3:  goto L30;
	case 4:  goto L35;
	case 5:  goto L40;
	case 6:  goto L45;
	case 7:  goto L50;
	case 8:  goto L60;
	case 9:  goto L70;
	case 10:  goto L80;
	case 11:  goto L90;
	case 12:  goto L150;
	case 13:  goto L110;
	case 14:  goto L120;
	case 15:  goto L130;
    }

L30:
    io___224.ciunit = pu;
    s_wsfe(&io___224);
    e_wsfe();
    goto L180;

L35:
    io___225.ciunit = pu;
    s_wsfe(&io___225);
    e_wsfe();
    goto L180;

L40:
    io___226.ciunit = pu;
    s_wsfe(&io___226);
    e_wsfe();
    goto L180;

L45:
    io___227.ciunit = pu;
    s_wsfe(&io___227);
    e_wsfe();
    goto L180;

L50:
    io___228.ciunit = pu;
    s_wsfe(&io___228);
    e_wsfe();
    goto L180;

L60:
    io___229.ciunit = pu;
    s_wsfe(&io___229);
    e_wsfe();
    goto L180;

L70:
    io___230.ciunit = pu;
    s_wsfe(&io___230);
    e_wsfe();
    goto L180;

L80:
    io___231.ciunit = pu;
    s_wsfe(&io___231);
    e_wsfe();
    goto L180;

L90:
    io___232.ciunit = pu;
    s_wsfe(&io___232);
    e_wsfe();
    goto L180;

L110:
    io___233.ciunit = pu;
    s_wsfe(&io___233);
    e_wsfe();

    goto L150;

L120:
    io___234.ciunit = pu;
    s_wsfe(&io___234);
    e_wsfe();
    goto L999;

L130:
    io___235.ciunit = pu;
    s_wsfe(&io___235);
    e_wsfe();
    if (iv[31] > 0) {
	goto L190;
    }
    goto L150;

L140:
    io___236.ciunit = pu;
    s_wsfe(&io___236);
    do_fio(&c__1, (char *)&iv1, (ftnlen)sizeof(integer));
    e_wsfe();
    goto L999;

/*  ***  INITIAL CALL ON ITSMRY  *** */

L150:
    if (iv[24] != 0) {
	io___237.ciunit = pu;
	s_wsfe(&io___237);
	i__1 = *p;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&d__[i__], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
    if (iv1 >= 13) {
	goto L999;
    }
    iv[39] = 0;
    iv[48] = 0;
    if (ol == 0) {
	goto L999;
    }
    if (ol < 0) {
	io___239.ciunit = pu;
	s_wsfe(&io___239);
	e_wsfe();
    }
    if (ol > 0) {
	io___240.ciunit = pu;
	s_wsfe(&io___240);
	e_wsfe();
    }
    io___241.ciunit = pu;
    s_wsfe(&io___241);
    do_fio(&c__1, (char *)&v[10], (ftnlen)sizeof(doublereal));
    e_wsfe();
    goto L999;

/*  ***  PRINT VARIOUS INFORMATION REQUESTED ON SOLUTION  *** */

L180:
    iv[39] = 1;
    if (iv[23] == 0) {
	goto L190;
    }
    oldf = v[13];
    preldf = 0.;
    nreldf = 0.;
    if (oldf <= 0.) {
	goto L185;
    }
    preldf = v[7] / oldf;
    nreldf = v[6] / oldf;
L185:
    nf = iv[6] - iv[40];
    ng = iv[30] - iv[41];
    io___243.ciunit = pu;
    s_wsfe(&io___243);
    do_fio(&c__1, (char *)&v[10], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&v[17], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&nf, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ng, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&preldf, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&nreldf, (ftnlen)sizeof(doublereal));
    e_wsfe();

    if (iv[40] > 0) {
	io___244.ciunit = pu;
	s_wsfe(&io___244);
	do_fio(&c__1, (char *)&iv[40], (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (iv[41] > 0) {
	io___245.ciunit = pu;
	s_wsfe(&io___245);
	do_fio(&c__1, (char *)&iv[41], (ftnlen)sizeof(integer));
	e_wsfe();
    }

L190:
    if (iv[22] == 0) {
	goto L210;
    }
    iv[39] = 1;
    g1 = iv[28];
    io___247.ciunit = pu;
    s_wsfe(&io___247);
    e_wsfe();
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	io___248.ciunit = pu;
	s_wsfe(&io___248);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&d__[i__], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&v[g1], (ftnlen)sizeof(doublereal));
	e_wsfe();
	++g1;
/* L200: */
    }

L210:
    if (iv[14] == 0) {
	goto L999;
    }
    cov1 = iv[26];
    iv[39] = 1;
    if (cov1 < 0) {
	goto L220;
    } else if (cov1 == 0) {
	goto L230;
    } else {
	goto L240;
    }
L220:
    if (-1 == cov1) {
	io___250.ciunit = pu;
	s_wsfe(&io___250);
	e_wsfe();
    }
    if (-2 == cov1) {
	io___251.ciunit = pu;
	s_wsfe(&io___251);
	e_wsfe();
    }
    goto L999;

L230:
    io___252.ciunit = pu;
    s_wsfe(&io___252);
    e_wsfe();
    goto L999;

L240:
    i__ = abs(iv[15]);
    if (i__ <= 1) {
	io___253.ciunit = pu;
	s_wsfe(&io___253);
	e_wsfe();
    }
    if (i__ == 2) {
	io___254.ciunit = pu;
	s_wsfe(&io___254);
	e_wsfe();
    }
    if (i__ >= 3) {
	io___255.ciunit = pu;
	s_wsfe(&io___255);
	e_wsfe();
    }
    ii = cov1 - 1;
    if (ol <= 0) {
	goto L260;
    }
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i1 = ii + 1;
	ii += i__;
	io___258.ciunit = pu;
	s_wsfe(&io___258);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	i__2 = ii;
	for (j = i1; j <= i__2; ++j) {
	    do_fio(&c__1, (char *)&v[j], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
/* L250: */
    }
    goto L999;

L260:
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i1 = ii + 1;
	ii += i__;
	io___260.ciunit = pu;
	s_wsfe(&io___260);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	i__2 = ii;
	for (j = i1; j <= i__2; ++j) {
	    do_fio(&c__1, (char *)&v[j], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
/* L270: */
    }

L999:
    return 0;
/*  ***  LAST CARD OF ITSMRY FOLLOWS  *** */
} /* itsmry_ */

/* Subroutine */ int linvrt_(integer *n, doublereal *lin, doublereal *l)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, k;
    static doublereal t;
    static integer j0, j1, k0, ii, jj, im1, np1;


/*  ***  COMPUTE  LIN = L**-1,  BOTH  N X N  LOWER TRIANG. STORED   *** */
/*  ***  COMPACTLY BY ROWS.  LIN AND L MAY SHARE THE SAME STORAGE.  *** */

/*  ***  PARAMETERS  *** */

/*     DIMENSION L(N*(N+1)/2), LIN(N*(N+1)/2) */

/*  ***  LOCAL VARIABLES  *** */

/* /6 */
/*     DATA ONE/1.D+0/, ZERO/0.D+0/ */
/* /7 */
/* / */

/*  ***  BODY  *** */

    /* Parameter adjustments */
    --l;
    --lin;

    /* Function Body */
    np1 = *n + 1;
    j0 = *n * np1 / 2;
    i__1 = *n;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = np1 - ii;
	lin[j0] = 1. / l[j0];
	if (i__ <= 1) {
	    goto L999;
	}
	j1 = j0;
	im1 = i__ - 1;
	i__2 = im1;
	for (jj = 1; jj <= i__2; ++jj) {
	    t = 0.;
	    j0 = j1;
	    k0 = j1 - jj;
	    i__3 = jj;
	    for (k = 1; k <= i__3; ++k) {
		t -= l[k0] * lin[j0];
		--j0;
		k0 = k0 + k - i__;
/* L10: */
	    }
	    lin[j0] = t / l[k0];
/* L20: */
	}
	--j0;
/* L30: */
    }
L999:
    return 0;
/*  ***  LAST CARD OF LINVRT FOLLOWS  *** */
} /* linvrt_ */

/* Subroutine */ int litvmu_(integer *n, doublereal *x, doublereal *l, 
	doublereal *y)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, i0, ii, ij;
    static doublereal xi;
    static integer im1, np1;


/*  ***  SOLVE  (L**T)*X = Y,  WHERE  L  IS AN  N X N  LOWER TRIANGULAR */
/*  ***  MATRIX STORED COMPACTLY BY ROWS.  X AND Y MAY OCCUPY THE SAME */
/*  ***  STORAGE.  *** */

/* /6 */
/*     DATA ZERO/0.D+0/ */
/* /7 */
/* / */

    /* Parameter adjustments */
    --y;
    --x;
    --l;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	x[i__] = y[i__];
    }
    np1 = *n + 1;
    i0 = *n * (*n + 1) / 2;
    i__1 = *n;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = np1 - ii;
	xi = x[i__] / l[i0];
	x[i__] = xi;
	if (i__ <= 1) {
	    goto L999;
	}
	i0 -= i__;
	if (xi == 0.) {
	    goto L30;
	}
	im1 = i__ - 1;
	i__2 = im1;
	for (j = 1; j <= i__2; ++j) {
	    ij = i0 + j;
	    x[j] -= xi * l[ij];
/* L20: */
	}
L30:
	;
    }
L999:
    return 0;
/*  ***  LAST CARD OF LITVMU FOLLOWS  *** */
} /* litvmu_ */

/* Subroutine */ int livmul_(integer *n, doublereal *x, doublereal *l, 
	doublereal *y)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k;
    static doublereal t;
    extern doublereal dotprd_(integer *, doublereal *, doublereal *);


/*  ***  SOLVE  L*X = Y, WHERE  L  IS AN  N X N  LOWER TRIANGULAR */
/*  ***  MATRIX STORED COMPACTLY BY ROWS.  X AND Y MAY OCCUPY THE SAME */
/*  ***  STORAGE.  *** */

/* /6 */
/*     DATA ZERO/0.D+0/ */
/* /7 */
/* / */

    /* Parameter adjustments */
    --y;
    --x;
    --l;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (y[k] != 0.) {
	    goto L20;
	}
	x[k] = 0.;
/* L10: */
    }
    goto L999;
L20:
    j = k * (k + 1) / 2;
    x[k] = y[k] / l[j];
    if (k >= *n) {
	goto L999;
    }
    ++k;
    i__1 = *n;
    for (i__ = k; i__ <= i__1; ++i__) {
	i__2 = i__ - 1;
	t = dotprd_(&i__2, &l[j + 1], &x[1]);
	j += i__;
	x[i__] = (y[i__] - t) / l[j];
/* L30: */
    }
L999:
    return 0;
/*  ***  LAST CARD OF LIVMUL FOLLOWS  *** */
} /* livmul_ */

/* Subroutine */ int lmstep_(doublereal *d__, doublereal *g, integer *ierr, 
	integer *ipivot, integer *ka, integer *p, doublereal *qtr, doublereal 
	*r__, doublereal *step, doublereal *v, doublereal *w)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal a, b;
    static integer i__, k, l;
    static doublereal t, d1, d2;
    static integer i1, j1;
    static doublereal lk, si, sj, uk, wl;
    static integer lk0, ip1, uk0;
    static doublereal adi, rad, phi;
    static integer res;
    static doublereal dst;
    static integer res0, pp1o2, rmat;
    static doublereal dtol;
    static integer rmat0, kalim;
    extern /* Subroutine */ int vcopy_(integer *, doublereal *, doublereal *);
    extern doublereal v2norm_(integer *, doublereal *);
    static doublereal alphak, dfacsq, psifac, oldphi, phimin, phimax;
    static integer phipin;
    extern doublereal dotprd_(integer *, doublereal *, doublereal *);
    static integer dstsav;
    static doublereal sqrtak;
    extern /* Subroutine */ int livmul_(integer *, doublereal *, doublereal *,
	     doublereal *), litvmu_(integer *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal twopsi;


/*  ***  COMPUTE LEVENBERG-MARQUARDT STEP USING MORE-HEBDEN TECHNIQUE  ** */
/*  ***  NL2SOL VERSION 2.2.  *** */

/*  ***  PARAMETER DECLARATIONS  *** */

/*     DIMENSION W(P*(P+5)/2 + 4) */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*  ***  PURPOSE  *** */

/*        GIVEN THE R MATRIX FROM THE QR DECOMPOSITION OF A JACOBIAN */
/*     MATRIX, J, AS WELL AS Q-TRANSPOSE TIMES THE CORRESPONDING */
/*     RESIDUAL VECTOR, RESID, THIS SUBROUTINE COMPUTES A LEVENBERG- */
/*     MARQUARDT STEP OF APPROXIMATE LENGTH V(RADIUS) BY THE MORE- */
/*     TECHNIQUE. */

/*  ***  PARAMETER DESCRIPTION  *** */

/*      D (IN)  = THE SCALE VECTOR. */
/*      G (IN)  = THE GRADIENT VECTOR (J**T)*R. */
/*   IERR (I/O) = RETURN CODE FROM QRFACT OR QRFGS -- 0 MEANS R HAS */
/*             FULL RANK. */
/* IPIVOT (I/O) = PERMUTATION ARRAY FROM QRFACT OR QRFGS, WHICH COMPUTE */
/*             QR DECOMPOSITIONS WITH COLUMN PIVOTING. */
/*     KA (I/O).  KA .LT. 0 ON INPUT MEANS THIS IS THE FIRST CALL ON */
/*             LMSTEP FOR THE CURRENT R AND QTR.  ON OUTPUT KA CON- */
/*             TAINS THE NUMBER OF HEBDEN ITERATIONS NEEDED TO DETERMINE */
/*             STEP.  KA = 0 MEANS A GAUSS-NEWTON STEP. */
/*      P (IN)  = NUMBER OF PARAMETERS. */
/*    QTR (IN)  = (Q**T)*RESID = Q-TRANSPOSE TIMES THE RESIDUAL VECTOR. */
/*      R (IN)  = THE R MATRIX, STORED COMPACTLY BY COLUMNS. */
/*   STEP (OUT) = THE LEVENBERG-MARQUARDT STEP COMPUTED. */
/*      V (I/O) CONTAINS VARIOUS CONSTANTS AND VARIABLES DESCRIBED BELOW. */
/*      W (I/O) = WORKSPACE OF LENGTH P*(P+5)/2 + 4. */

/*  ***  ENTRIES IN V  *** */

/* V(DGNORM) (I/O) = 2-NORM OF (D**-1)*G. */
/* V(DSTNRM) (I/O) = 2-NORM OF D*STEP. */
/* V(DST0)   (I/O) = 2-NORM OF GAUSS-NEWTON STEP (FOR NONSING. J). */
/* V(EPSLON) (IN) = MAX. REL. ERROR ALLOWED IN TWONORM(R)**2 MINUS */
/*             TWONORM(R - J*STEP)**2.  (SEE ALGORITHM NOTES BELOW.) */
/* V(GTSTEP) (OUT) = INNER PRODUCT BETWEEN G AND STEP. */
/* V(NREDUC) (OUT) = HALF THE REDUCTION IN THE SUM OF SQUARES PREDICTED */
/*             FOR A GAUSS-NEWTON STEP. */
/* V(PHMNFC) (IN)  = TOL. (TOGETHER WITH V(PHMXFC)) FOR ACCEPTING STEP */
/*             (MORE*S SIGMA).  THE ERROR V(DSTNRM) - V(RADIUS) MUST LIE */
/*             BETWEEN V(PHMNFC)*V(RADIUS) AND V(PHMXFC)*V(RADIUS). */
/* V(PHMXFC) (IN)  (SEE V(PHMNFC).) */
/* V(PREDUC) (OUT) = HALF THE REDUCTION IN THE SUM OF SQUARES PREDICTED */
/*             BY THE STEP RETURNED. */
/* V(RADIUS) (IN)  = RADIUS OF CURRENT (SCALED) TRUST REGION. */
/* V(RAD0)   (I/O) = VALUE OF V(RADIUS) FROM PREVIOUS CALL. */
/* V(STPPAR) (I/O) = MARQUARDT PARAMETER (OR ITS NEGATIVE IF THE SPECIAL */
/*             CASE MENTIONED BELOW IN THE ALGORITHM NOTES OCCURS). */

/* NOTE -- SEE DATA STATEMENT BELOW FOR VALUES OF ABOVE SUBSCRIPTS. */

/*  ***  USAGE NOTES  *** */

/*     IF IT IS DESIRED TO RECOMPUTE STEP USING A DIFFERENT VALUE OF */
/*     V(RADIUS), THEN THIS ROUTINE MAY BE RESTARTED BY CALLING IT */
/*     WITH ALL PARAMETERS UNCHANGED EXCEPT V(RADIUS).  (THIS EXPLAINS */
/*     WHY MANY PARAMETERS ARE LISTED AS I/O).  ON AN INTIIAL CALL (ONE */
/*     WITH KA = -1), THE CALLER NEED ONLY HAVE INITIALIZED D, G, KA, P, */
/*     QTR, R, V(EPSLON), V(PHMNFC), V(PHMXFC), V(RADIUS), AND V(RAD0). */

/*  ***  APPLICATION AND USAGE RESTRICTIONS  *** */

/*     THIS ROUTINE IS CALLED AS PART OF THE NL2SOL (NONLINEAR LEAST- */
/*     SQUARES) PACKAGE (REF. 1). */

/*  ***  ALGORITHM NOTES  *** */

/*     THIS CODE IMPLEMENTS THE STEP COMPUTATION SCHEME DESCRIBED IN */
/*     REFS. 2 AND 4.  FAST GIVENS TRANSFORMATIONS (SEE REF. 3, PP. 60- */
/*     62) ARE USED TO COMPUTE STEP WITH A NONZERO MARQUARDT PARAMETER. */
/*        A SPECIAL CASE OCCURS IF J IS (NEARLY) SINGULAR AND V(RADIUS) */
/*     IS SUFFICIENTLY LARGE.  IN THIS CASE THE STEP RETURNED IS SUCH */
/*     THAT  TWONORM(R)**2 - TWONORM(R - J*STEP)**2  DIFFERS FROM ITS */
/*     OPTIMAL VALUE BY LESS THAN V(EPSLON) TIMES THIS OPTIMAL VALUE, */
/*     WHERE J AND R DENOTE THE ORIGINAL JACOBIAN AND RESIDUAL.  (SEE */
/*     REF. 2 FOR MORE DETAILS.) */

/*  ***  FUNCTIONS AND SUBROUTINES CALLED  *** */

/* DOTPRD - RETURNS INNER PRODUCT OF TWO VECTORS. */
/* LITVMU - APPLY INVERSE-TRANSPOSE OF COMPACT LOWER TRIANG. MATRIX. */
/* LIVMUL - APPLY INVERSE OF COMPACT LOWER TRIANG. MATRIX. */
/* VCOPY  - COPIES ONE VECTOR TO ANOTHER. */
/* V2NORM - RETURNS 2-NORM OF A VECTOR. */

/*  ***  REFERENCES  *** */

/* 1.  DENNIS, J.E., GAY, D.M., AND WELSCH, R.E. (1981), AN ADAPTIVE */
/*             NONLINEAR LEAST-SQUARES ALGORITHM, ACM TRANS. MATH. */
/*             SOFTWARE, VOL. 7, NO. 3. */
/* 2.  GAY, D.M. (1981), COMPUTING OPTIMAL LOCALLY CONSTRAINED STEPS, */
/*             SIAM J. SCI. STATIST. COMPUTING, VOL. 2, NO. 2, PP. */
/*             186-197. */
/* 3.  LAWSON, C.L., AND HANSON, R.J. (1974), SOLVING LEAST SQUARES */
/*             PROBLEMS, PRENTICE-HALL, ENGLEWOOD CLIFFS, N.J. */
/* 4.  MORE, J.J. (1978), THE LEVENBERG-MARQUARDT ALGORITHM, IMPLEMEN- */
/*             TATION AND THEORY, PP.105-116 OF SPRINGER LECTURE NOTES */
/*             IN MATHEMATICS NO. 630, EDITED BY G.A. WATSON, SPRINGER- */
/*             VERLAG, BERLIN AND NEW YORK. */

/*  ***  GENERAL  *** */

/*     CODED BY DAVID M. GAY. */
/*     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH */
/*     SUPPORTED BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS */
/*     MCS-7600324, DCR75-10143, 76-14311DSS, MCS76-11989, AND */
/*     MCS-7906671. */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*  ***  LOCAL VARIABLES  *** */


/*     ***  CONSTANTS  *** */

/*  ***  INTRINSIC FUNCTIONS  *** */
/* /+ */
/* / */
/*  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  *** */


/*  ***  SUBSCRIPTS FOR V  *** */

/* /6 */
/*     DATA DGNORM/1/, DSTNRM/2/, DST0/3/, EPSLON/19/, */
/*    1     GTSTEP/4/, NREDUC/6/, PHMNFC/20/, */
/*    2     PHMXFC/21/, PREDUC/7/, RADIUS/8/, */
/*    3     RAD0/9/, STPPAR/5/ */
/* /7 */
/* / */

/* /6 */
/*     DATA DFAC/256.D+0/, EIGHT/8.D+0/, HALF/0.5D+0/, NEGONE/-1.D+0/, */
/*    1     ONE/1.D+0/, P001/1.D-3/, THREE/3.D+0/, TTOL/2.5D+0/, */
/*    2     ZERO/0.D+0/ */
/* /7 */
/* / */

/*  ***  BODY  *** */

/*     ***  FOR USE IN RECOMPUTING STEP, THE FINAL VALUES OF LK AND UK, */
/*     ***  THE INVERSE DERIVATIVE OF MORE*S PHI AT 0 (FOR NONSING. J) */
/*     ***  AND THE VALUE RETURNED AS V(DSTNRM) ARE STORED AT W(LK0), */
/*     ***  W(UK0), W(PHIPIN), AND W(DSTSAV) RESPECTIVELY. */
    /* Parameter adjustments */
    --step;
    --qtr;
    --ipivot;
    --g;
    --d__;
    --r__;
    --v;
    --w;

    /* Function Body */
    lk0 = *p + 1;
    phipin = lk0 + 1;
    uk0 = phipin + 1;
    dstsav = uk0 + 1;
    rmat0 = dstsav;
/*     ***  A COPY OF THE R-MATRIX FROM THE QR DECOMPOSITION OF J IS */
/*     ***  STORED IN W STARTING AT W(RMAT), AND A COPY OF THE RESIDUAL */
/*     ***  VECTOR IS STORED IN W STARTING AT W(RES).  THE LOOPS BELOW */
/*     ***  THAT UPDATE THE QR DECOMP. FOR A NONZERO MARQUARDT PARAMETER */
/*     ***  WORK ON THESE COPIES. */
    rmat = rmat0 + 1;
    pp1o2 = *p * (*p + 1) / 2;
    res0 = pp1o2 + rmat0;
    res = res0 + 1;
    rad = v[8];
    if (rad > 0.) {
/* Computing 2nd power */
	d__1 = rad;
	psifac = v[19] / (((v[20] + 1.) * 8. + 3.) * (d__1 * d__1));
    }
    phimax = v[21] * rad;
    phimin = v[20] * rad;
/*     ***  DTOL, DFAC, AND DFACSQ ARE USED IN RESCALING THE FAST GIVENS */
/*     ***  REPRESENTATION OF THE UPDATED QR DECOMPOSITION. */
    dtol = .00390625;
    dfacsq = 65536.;
/*     ***  OLDPHI IS USED TO DETECT LIMITS OF NUMERICAL ACCURACY.  IF */
/*     ***  WE RECOMPUTE STEP AND IT DOES NOT CHANGE, THEN WE ACCEPT IT. */
    oldphi = 0.;
    lk = 0.;
    uk = 0.;
    kalim = *ka + 12;

/*  ***  START OR RESTART, DEPENDING ON KA  *** */

    if (*ka < 0) {
	goto L10;
    } else if (*ka == 0) {
	goto L20;
    } else {
	goto L370;
    }

/*  ***  FRESH START -- COMPUTE V(NREDUC)  *** */

L10:
    *ka = 0;
    kalim = 12;
    k = *p;
    if (*ierr != 0) {
	k = abs(*ierr) - 1;
    }
    v[6] = dotprd_(&k, &qtr[1], &qtr[1]) * .5;

/*  ***  SET UP TO TRY INITIAL GAUSS-NEWTON STEP  *** */

L20:
    v[3] = -1.;
    if (*ierr != 0) {
	goto L90;
    }

/*  ***  COMPUTE GAUSS-NEWTON STEP  *** */

/*     ***  NOTE -- THE R-MATRIX IS STORED COMPACTLY BY COLUMNS IN */
/*     ***  R(1), R(2), R(3), ...  IT IS THE TRANSPOSE OF A */
/*     ***  LOWER TRIANGULAR MATRIX STORED COMPACTLY BY ROWS, AND WE */
/*     ***  TREAT IT AS SUCH WHEN USING LITVMU AND LIVMUL. */
    litvmu_(p, &w[1], &r__[1], &qtr[1]);
/*     ***  TEMPORARILY STORE PERMUTED -D*STEP IN STEP. */
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j1 = ipivot[i__];
	step[i__] = d__[j1] * w[i__];
/* L60: */
    }
    dst = v2norm_(p, &step[1]);
    v[3] = dst;
    phi = dst - rad;
    if (phi <= phimax) {
	goto L410;
    }
/*     ***  IF THIS IS A RESTART, GO TO 110  *** */
    if (*ka > 0) {
	goto L110;
    }

/*  ***  GAUSS-NEWTON STEP WAS UNACCEPTABLE.  COMPUTE L0  *** */

    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j1 = ipivot[i__];
	step[i__] = d__[j1] * (step[i__] / dst);
/* L70: */
    }
    livmul_(p, &step[1], &r__[1], &step[1]);
    t = 1. / v2norm_(p, &step[1]);
    w[phipin] = t / dst * t;
    lk = phi * w[phipin];

/*  ***  COMPUTE U0  *** */

L90:
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L100: */
	w[i__] = g[i__] / d__[i__];
    }
    v[1] = v2norm_(p, &w[1]);
    uk = v[1] / rad;
    if (uk <= 0.) {
	goto L390;
    }

/*     ***  ALPHAK WILL BE USED AS THE CURRENT MARQUARDT PARAMETER.  WE */
/*     ***  USE MORE*S SCHEME FOR INITIALIZING IT. */
    alphak = abs(v[5]) * v[9] / rad;


/*  ***  TOP OF LOOP -- INCREMENT KA, COPY R TO RMAT, QTR TO RES  *** */

L110:
    ++(*ka);
    vcopy_(&pp1o2, &w[rmat], &r__[1]);
    vcopy_(p, &w[res], &qtr[1]);

/*  ***  SAFEGUARD ALPHAK AND INITIALIZE FAST GIVENS SCALE VECTOR.  *** */

    if (alphak <= 0. || alphak < lk || alphak >= uk) {
/* Computing MAX */
	d__1 = .001, d__2 = sqrt(lk / uk);
	alphak = uk * max(d__1,d__2);
    }
    sqrtak = sqrt(alphak);
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L120: */
	w[i__] = 1.;
    }

/*  ***  ADD ALPHAK*D AND UPDATE QR DECOMP. USING FAST GIVENS TRANS.  *** */

    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*        ***  GENERATE, APPLY 1ST GIVENS TRANS. FOR ROW I OF ALPHAK*D. */
/*        ***  (USE STEP TO STORE TEMPORARY ROW)  *** */
	l = i__ * (i__ + 1) / 2 + rmat0;
	wl = w[l];
	d2 = 1.;
	d1 = w[i__];
	j1 = ipivot[i__];
	adi = sqrtak * d__[j1];
	if (adi >= abs(wl)) {
	    goto L150;
	}
L130:
	a = adi / wl;
	b = d2 * a / d1;
	t = a * b + 1.;
	if (t > 2.5) {
	    goto L150;
	}
	w[i__] = d1 / t;
	d2 /= t;
	w[l] = t * wl;
	a = -a;
	i__2 = *p;
	for (j1 = i__; j1 <= i__2; ++j1) {
	    l += j1;
	    step[j1] = a * w[l];
/* L140: */
	}
	goto L170;

L150:
	b = wl / adi;
	a = d1 * b / d2;
	t = a * b + 1.;
	if (t > 2.5) {
	    goto L130;
	}
	w[i__] = d2 / t;
	d2 = d1 / t;
	w[l] = t * adi;
	i__2 = *p;
	for (j1 = i__; j1 <= i__2; ++j1) {
	    l += j1;
	    wl = w[l];
	    step[j1] = -wl;
	    w[l] = a * wl;
/* L160: */
	}

L170:
	if (i__ == *p) {
	    goto L280;
	}

/*        ***  NOW USE GIVENS TRANS. TO ZERO ELEMENTS OF TEMP. ROW  *** */

	ip1 = i__ + 1;
	i__2 = *p;
	for (i1 = ip1; i1 <= i__2; ++i1) {
	    l = i1 * (i1 + 1) / 2 + rmat0;
	    wl = w[l];
	    si = step[i1 - 1];
	    d1 = w[i1];

/*             ***  RESCALE ROW I1 IF NECESSARY  *** */

	    if (d1 >= dtol) {
		goto L190;
	    }
	    d1 *= dfacsq;
	    wl /= 256.;
	    k = l;
	    i__3 = *p;
	    for (j1 = i1; j1 <= i__3; ++j1) {
		k += j1;
		w[k] /= 256.;
/* L180: */
	    }

/*             ***  USE GIVENS TRANS. TO ZERO NEXT ELEMENT OF TEMP. ROW */

L190:
	    if (abs(si) > abs(wl)) {
		goto L220;
	    }
	    if (si == 0.) {
		goto L260;
	    }
L200:
	    a = si / wl;
	    b = d2 * a / d1;
	    t = a * b + 1.;
	    if (t > 2.5) {
		goto L220;
	    }
	    w[l] = t * wl;
	    w[i1] = d1 / t;
	    d2 /= t;
	    i__3 = *p;
	    for (j1 = i1; j1 <= i__3; ++j1) {
		l += j1;
		wl = w[l];
		sj = step[j1];
		w[l] = wl + b * sj;
		step[j1] = sj - a * wl;
/* L210: */
	    }
	    goto L240;

L220:
	    b = wl / si;
	    a = d1 * b / d2;
	    t = a * b + 1.;
	    if (t > 2.5) {
		goto L200;
	    }
	    w[i1] = d2 / t;
	    d2 = d1 / t;
	    w[l] = t * si;
	    i__3 = *p;
	    for (j1 = i1; j1 <= i__3; ++j1) {
		l += j1;
		wl = w[l];
		sj = step[j1];
		w[l] = a * wl + sj;
		step[j1] = b * sj - wl;
/* L230: */
	    }

/*             ***  RESCALE TEMP. ROW IF NECESSARY  *** */

L240:
	    if (d2 >= dtol) {
		goto L260;
	    }
	    d2 *= dfacsq;
	    i__3 = *p;
	    for (k = i1; k <= i__3; ++k) {
/* L250: */
		step[k] /= 256.;
	    }
L260:
	    ;
	}
/* L270: */
    }

/*  ***  COMPUTE STEP  *** */

L280:
    litvmu_(p, &w[res], &w[rmat], &w[res]);
/*     ***  RECOVER STEP AND STORE PERMUTED -D*STEP AT W(RES)  *** */
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j1 = ipivot[i__];
	k = res0 + i__;
	t = w[k];
	step[j1] = -t;
	w[k] = t * d__[j1];
/* L290: */
    }
    dst = v2norm_(p, &w[res]);
    phi = dst - rad;
    if (phi <= phimax && phi >= phimin) {
	goto L430;
    }
    if (oldphi == phi) {
	goto L430;
    }
    oldphi = phi;

/*  ***  CHECK FOR (AND HANDLE) SPECIAL CASE  *** */

    if (phi > 0.) {
	goto L310;
    }
    if (*ka >= kalim) {
	goto L430;
    }
    twopsi = alphak * dst * dst - dotprd_(p, &step[1], &g[1]);
    if (alphak >= twopsi * psifac) {
	goto L310;
    }
    v[5] = -alphak;
    goto L440;

/*  ***  UNACCEPTABLE STEP -- UPDATE LK, UK, ALPHAK, AND TRY AGAIN  *** */

L300:
    if (phi < 0.) {
	uk = min(uk,alphak);
    }
    goto L320;
L310:
    if (phi < 0.) {
	uk = alphak;
    }
L320:
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j1 = ipivot[i__];
	k = res0 + i__;
	step[i__] = d__[j1] * (w[k] / dst);
/* L330: */
    }
    livmul_(p, &step[1], &w[rmat], &step[1]);
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L340: */
	step[i__] /= sqrt(w[i__]);
    }
    t = 1. / v2norm_(p, &step[1]);
    alphak += t * phi * t / rad;
    lk = max(lk,alphak);
    goto L110;

/*  ***  RESTART  *** */

L370:
    lk = w[lk0];
    uk = w[uk0];
    if (v[3] > 0. && v[3] - rad <= phimax) {
	goto L20;
    }
    alphak = abs(v[5]);
    dst = w[dstsav];
    phi = dst - rad;
    t = v[1] / rad;
    if (rad > v[9]) {
	goto L380;
    }

/*        ***  SMALLER RADIUS  *** */
    uk = t;
    if (alphak <= 0.) {
	lk = 0.;
    }
    if (v[3] > 0.) {
/* Computing MAX */
	d__1 = lk, d__2 = (v[3] - rad) * w[phipin];
	lk = max(d__1,d__2);
    }
    goto L300;

/*     ***  BIGGER RADIUS  *** */
L380:
    if (alphak <= 0. || uk > t) {
	uk = t;
    }
    lk = 0.;
    if (v[3] > 0.) {
/* Computing MAX */
	d__1 = lk, d__2 = (v[3] - rad) * w[phipin];
	lk = max(d__1,d__2);
    }
    goto L300;

/*  ***  SPECIAL CASE -- RAD .LE. 0 OR (G = 0 AND J IS SINGULAR)  *** */

L390:
    v[5] = 0.;
    dst = 0.;
    lk = 0.;
    uk = 0.;
    v[4] = 0.;
    v[7] = 0.;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L400: */
	step[i__] = 0.;
    }
    goto L450;

/*  ***  ACCEPTABLE GAUSS-NEWTON STEP -- RECOVER STEP FROM W  *** */

L410:
    alphak = 0.;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j1 = ipivot[i__];
	step[j1] = -w[i__];
/* L420: */
    }

/*  ***  SAVE VALUES FOR USE IN A POSSIBLE RESTART  *** */

L430:
    v[5] = alphak;
L440:
    v[4] = dotprd_(p, &step[1], &g[1]);
    v[7] = (alphak * dst * dst - v[4]) * .5;
L450:
    v[2] = dst;
    w[dstsav] = dst;
    w[lk0] = lk;
    w[uk0] = uk;
    v[9] = rad;

/* L999: */
    return 0;

/*  ***  LAST CARD OF LMSTEP FOLLOWS  *** */
} /* lmstep_ */

/* Subroutine */ int lsqrt_(integer *n1, integer *n, doublereal *l, 
	doublereal *a, integer *irc)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal t;
    static integer i0, j0, ij, ik, jk;
    static doublereal td;
    static integer im1, jm1;


/*  ***  COMPUTE ROWS N1 THROUGH N OF THE CHOLESKY FACTOR  L  OF */
/*  ***  A = L*(L**T),  WHERE  L  AND THE LOWER TRIANGLE OF  A  ARE BOTH */
/*  ***  STORED COMPACTLY BY ROWS (AND MAY OCCUPY THE SAME STORAGE). */
/*  ***  IRC = 0 MEANS ALL WENT WELL.  IRC = J MEANS THE LEADING */
/*  ***  PRINCIPAL  J X J  SUBMATRIX OF  A  IS NOT POSITIVE DEFINITE -- */
/*  ***  AND  L(J*(J+1)/2)  CONTAINS THE (NONPOS.) REDUCED J-TH DIAGONAL. */

/*  ***  PARAMETERS  *** */

/*     DIMENSION L(N*(N+1)/2), A(N*(N+1)/2) */

/*  ***  LOCAL VARIABLES  *** */


/*  ***  INTRINSIC FUNCTIONS  *** */
/* /+ */
/* / */
/* /6 */
/*     DATA ZERO/0.D+0/ */
/* /7 */
/* / */

/*  ***  BODY  *** */

    /* Parameter adjustments */
    --a;
    --l;

    /* Function Body */
    i0 = *n1 * (*n1 - 1) / 2;
    i__1 = *n;
    for (i__ = *n1; i__ <= i__1; ++i__) {
	td = 0.;
	if (i__ == 1) {
	    goto L40;
	}
	j0 = 0;
	im1 = i__ - 1;
	i__2 = im1;
	for (j = 1; j <= i__2; ++j) {
	    t = 0.;
	    if (j == 1) {
		goto L20;
	    }
	    jm1 = j - 1;
	    i__3 = jm1;
	    for (k = 1; k <= i__3; ++k) {
		ik = i0 + k;
		jk = j0 + k;
		t += l[ik] * l[jk];
/* L10: */
	    }
L20:
	    ij = i0 + j;
	    j0 += j;
	    t = (a[ij] - t) / l[j0];
	    l[ij] = t;
	    td += t * t;
/* L30: */
	}
L40:
	i0 += i__;
	t = a[i0] - td;
	if (t <= 0.) {
	    goto L60;
	}
	l[i0] = sqrt(t);
/* L50: */
    }

    *irc = 0;
    goto L999;

L60:
    l[i0] = t;
    *irc = i__;

L999:
    return 0;

/*  ***  LAST CARD OF LSQRT  *** */
} /* lsqrt_ */

doublereal lsvmin_(integer *p, doublereal *l, doublereal *x, doublereal *y)
{
    /* Initialized data */

    static integer ix = 2;

    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1;

    /* Local variables */
    static doublereal b;
    static integer i__, j;
    static doublereal t;
    static integer j0, ii, ji, jj, jm1, jjj;
    static doublereal psj, splus, xplus;
    extern doublereal v2norm_(integer *, doublereal *);
    static integer pplus1;
    static doublereal sminus, xminus;


/*  ***  ESTIMATE SMALLEST SING. VALUE OF PACKED LOWER TRIANG. MATRIX L */

/*  ***  PARAMETER DECLARATIONS  *** */

/*     DIMENSION L(P*(P+1)/2) */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*  ***  PURPOSE  *** */

/*     THIS FUNCTION RETURNS A GOOD OVER-ESTIMATE OF THE SMALLEST */
/*     SINGULAR VALUE OF THE PACKED LOWER TRIANGULAR MATRIX L. */

/*  ***  PARAMETER DESCRIPTION  *** */

/*  P (IN)  = THE ORDER OF L.  L IS A  P X P  LOWER TRIANGULAR MATRIX. */
/*  L (IN)  = ARRAY HOLDING THE ELEMENTS OF  L  IN ROW ORDER, I.E. */
/*             L(1,1), L(2,1), L(2,2), L(3,1), L(3,2), L(3,3), ETC. */
/*  X (OUT) IF LSVMIN RETURNS A POSITIVE VALUE, THEN X IS A NORMALIZED */
/*             APPROXIMATE LEFT SINGULAR VECTOR CORRESPONDING TO THE */
/*             SMALLEST SINGULAR VALUE.  THIS APPROXIMATION MAY BE VERY */
/*             CRUDE.  IF LSVMIN RETURNS ZERO, THEN SOME COMPONENTS OF X */
/*             ARE ZERO AND THE REST RETAIN THEIR INPUT VALUES. */
/*  Y (OUT) IF LSVMIN RETURNS A POSITIVE VALUE, THEN Y = (L**-1)*X IS AN */
/*             UNNORMALIZED APPROXIMATE RIGHT SINGULAR VECTOR CORRESPOND- */
/*             ING TO THE SMALLEST SINGULAR VALUE.  THIS APPROXIMATION */
/*             MAY BE CRUDE.  IF LSVMIN RETURNS ZERO, THEN Y RETAINS ITS */
/*             INPUT VALUE.  THE CALLER MAY PASS THE SAME VECTOR FOR X */
/*             AND Y (NONSTANDARD FORTRAN USAGE), IN WHICH CASE Y OVER- */
/*             WRITES X (FOR NONZERO LSVMIN RETURNS). */

/*  ***  APPLICATION AND USAGE RESTRICTIONS  *** */

/*     THERE ARE NO USAGE RESTRICTIONS. */

/*  ***  ALGORITHM NOTES  *** */

/*     THE ALGORITHM IS BASED ON (1), WITH THE ADDITIONAL PROVISION THAT */
/*     LSVMIN = 0 IS RETURNED IF THE SMALLEST DIAGONAL ELEMENT OF L */
/*     (IN MAGNITUDE) IS NOT MORE THAN THE UNIT ROUNDOFF TIMES THE */
/*     LARGEST.  THE ALGORITHM USES A RANDOM NUMBER GENERATOR PROPOSED */
/*     IN (4), WHICH PASSES THE SPECTRAL TEST WITH FLYING COLORS -- SEE */
/*     (2) AND (3). */

/*  ***  SUBROUTINES AND FUNCTIONS CALLED  *** */

/*        V2NORM - FUNCTION, RETURNS THE 2-NORM OF A VECTOR. */

/*  ***  REFERENCES  *** */

/*     (1) CLINE, A., MOLER, C., STEWART, G., AND WILKINSON, J.H.(1977), */
/*         AN ESTIMATE FOR THE CONDITION NUMBER OF A MATRIX, REPORT */
/*         TM-310, APPLIED MATH. DIV., ARGONNE NATIONAL LABORATORY. */

/*     (2) HOAGLIN, D.C. (1976), THEORETICAL PROPERTIES OF CONGRUENTIAL */
/*         RANDOM-NUMBER GENERATORS --  AN EMPIRICAL VIEW, */
/*         MEMORANDUM NS-340, DEPT. OF STATISTICS, HARVARD UNIV. */

/*     (3) KNUTH, D.E. (1969), THE ART OF COMPUTER PROGRAMMING, VOL. 2 */
/*         (SEMINUMERICAL ALGORITHMS), ADDISON-WESLEY, READING, MASS. */

/*     (4) SMITH, C.S. (1971), MULTIPLICATIVE PSEUDO-RANDOM NUMBER */
/*         GENERATORS WITH PRIME MODULUS, J. ASSOC. COMPUT. MACH. 18, */
/*         PP. 586-593. */

/*  ***  HISTORY  *** */

/*     DESIGNED AND CODED BY DAVID M. GAY (WINTER 1977/SUMMER 1978). */

/*  ***  GENERAL  *** */

/*     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH */
/*     SUPPORTED BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS */
/*     MCS-7600324, DCR75-10143, 76-14311DSS, AND MCS76-11989. */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*  ***  LOCAL VARIABLES  *** */


/*  ***  CONSTANTS  *** */


/*  ***  INTRINSIC FUNCTIONS  *** */
/* /+ */
/* / */
/*  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  *** */


/* /6 */
/*     DATA HALF/0.5D+0/, ONE/1.D+0/, R9973/9973.D+0/, ZERO/0.D+0/ */
/* /7 */
/*     SAVE IX */
/* / */
    /* Parameter adjustments */
    --y;
    --x;
    --l;

    /* Function Body */

/*  ***  BODY  *** */

/*  ***  FIRST CHECK WHETHER TO RETURN LSVMIN = 0 AND INITIALIZE X  *** */

    ii = 0;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = 0.;
	ii += i__;
	if (l[ii] == 0.) {
	    goto L300;
	}
/* L10: */
    }
    if (ix % 9973 == 0) {
	ix = 2;
    }
    pplus1 = *p + 1;

/*  ***  SOLVE (L**T)*X = B, WHERE THE COMPONENTS OF B HAVE RANDOMLY */
/*  ***  CHOSEN MAGNITUDES IN (.5,1) WITH SIGNS CHOSEN TO MAKE X LARGE. */

/*     DO J = P TO 1 BY -1... */
    i__1 = *p;
    for (jjj = 1; jjj <= i__1; ++jjj) {
	j = pplus1 - jjj;
/*       ***  DETERMINE X(J) IN THIS ITERATION. NOTE FOR I = 1,2,...,J */
/*       ***  THAT X(I) HOLDS THE CURRENT PARTIAL SUM FOR ROW I. */
	ix = ix * 3432 % 9973;
	b = ((real) ix / 9973. + 1.) * .5;
	xplus = b - x[j];
	xminus = -b - x[j];
	splus = abs(xplus);
	sminus = abs(xminus);
	jm1 = j - 1;
	j0 = j * jm1 / 2;
	jj = j0 + j;
	xplus /= l[jj];
	xminus /= l[jj];
	if (jm1 == 0) {
	    goto L30;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ji = j0 + i__;
	    splus += (d__1 = x[i__] + l[ji] * xplus, abs(d__1));
	    sminus += (d__1 = x[i__] + l[ji] * xminus, abs(d__1));
/* L20: */
	}
L30:
	if (sminus > splus) {
	    xplus = xminus;
	}
	x[j] = xplus;
/*       ***  UPDATE PARTIAL SUMS  *** */
	if (jm1 == 0) {
	    goto L100;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ji = j0 + i__;
	    x[i__] += l[ji] * xplus;
/* L40: */
	}
L100:
	;
    }

/*  ***  NORMALIZE X  *** */

    t = 1. / v2norm_(p, &x[1]);
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L110: */
	x[i__] = t * x[i__];
    }

/*  ***  SOLVE L*Y = X AND RETURN SVMIN = 1/TWONORM(Y)  *** */

    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	psj = 0.;
	jm1 = j - 1;
	j0 = j * jm1 / 2;
	if (jm1 == 0) {
	    goto L130;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ji = j0 + i__;
	    psj += l[ji] * y[i__];
/* L120: */
	}
L130:
	jj = j0 + j;
	y[j] = (x[j] - psj) / l[jj];
/* L200: */
    }

    ret_val = 1. / v2norm_(p, &y[1]);
    goto L999;

L300:
    ret_val = 0.;
L999:
    return ret_val;
/*  ***  LAST CARD OF LSVMIN FOLLOWS  *** */
} /* lsvmin_ */

/* Subroutine */ int ltsqar_(integer *n, doublereal *a, doublereal *l)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, m, i1, ii;
    static doublereal lj, lii;
    static integer iim1;


/*  ***  SET A TO LOWER TRIANGLE OF (L**T) * L  *** */

/*  ***  L = N X N LOWER TRIANG. MATRIX STORED ROWWISE.  *** */
/*  ***  A IS ALSO STORED ROWWISE AND MAY SHARE STORAGE WITH L.  *** */

/*     DIMENSION A(N*(N+1)/2), L(N*(N+1)/2) */


    /* Parameter adjustments */
    --l;
    --a;

    /* Function Body */
    ii = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i1 = ii + 1;
	ii += i__;
	m = 1;
	if (i__ == 1) {
	    goto L30;
	}
	iim1 = ii - 1;
	i__2 = iim1;
	for (j = i1; j <= i__2; ++j) {
	    lj = l[j];
	    i__3 = j;
	    for (k = i1; k <= i__3; ++k) {
		a[m] += lj * l[k];
		++m;
/* L10: */
	    }
/* L20: */
	}
L30:
	lii = l[ii];
	i__2 = ii;
	for (j = i1; j <= i__2; ++j) {
/* L40: */
	    a[j] = lii * l[j];
	}
/* L50: */
    }

/* L999: */
    return 0;
/*  ***  LAST CARD OF LTSQAR FOLLOWS  *** */
} /* ltsqar_ */




//     PARAMETER CHECK ----
/* Subroutine */ int parchk_(integer *iv, integer *n, integer *nn, integer *p,
	 doublereal *v)
{
	if (COMMENT)
		printf("\nparchk : entered");
	/* Initialized data */

    static doublereal big = 0.;
    static doublereal tiny = 1.;
    static char vn[4*2*27+1] = "EPSLON..PHMNFC..PHMXFC..DECFAC..INCFAC..RDFC\
MN..RDFCMX..TUNER1..TUNER2..TUNER3..TUNER4..TUNER5..AFCTOL..RFCTOL..XCTOL...\
XFTOL...LMAX0...DLTFDJ..D0INIT..DINIT...JTINIT..DLTFDC..DFAC....RLIMIT..COSM\
IN..DELTA0..FUZZ....";
    static doublereal vm[27] = { .001,-.99,.001,.01,1.2,.01,1.2,0.,0.,.001,
	    -1.,0.0,0.0,0.0,0.,0.,0.0,0.0,0.,-10.,0.,0.0,0.,1e10,0.0,0.0,1.01 
	    };
    static doublereal vx[27] = { .9,-.001,10.,.8,100.,.8,100.,.5,.5,1.,1.,0.0,
	    0.0,.1,1.,1.,0.0,1.,0.0,0.0,0.0,1.,1.,0.0,1.,1.,100. };
    static char cngd[4*3+1] = "---CHANGED V";
    static char dflt[4*3+1] = "NONDEFAULT V";

    /* Format strings */
    static char fmt_10[] = "(\0020///// BAD NN, N, OR P... NN =\002,i5,\002,\
 N =\002,i5,\002, P =\002,i5)";
    static char fmt_40[] = "(\0020///// (NN,N,P) CHANGED FROM (\002,i5,\002\
,\002,i5,\002,\002,i3,\002) TO (\002,i5,\002,\002,i5,\002,\002,i3,\002).\002)"
	    ;
    static char fmt_60[] = "(\0020/////  IV(1) =\002,i5,\002 SHOULD BE BETWE\
EN 0 AND 12.\002)";
    static char fmt_100[] = "(\0020/////  INITS... IV(25) =\002,i4,\002 SHOU\
LD BE BETWEEN 0\002,\002 AND 2.\002)";
    static char fmt_120[] = "(\0020/////  \002,2a4,\002.. V(\002,i2,\002) \
=\002,d11.3,\002 SHOULD\002,\002 BE BETWEEN\002,d11.3,\002 AND\002,d11.3)";
    static char fmt_150[] = "(\0020///// JTOL(\002,i3,\002) = V(\002,i3,\002\
) =\002,d11.3,\002 SHOULD BE POSITIVE.\002)";
    static char fmt_190[] = "(\0020NONDEFAULT VALUES....\002/\002 INITS.....\
 IV(25) =\002,i3)";
    static char fmt_215[] = "(\0020\002,3a4,\002ALUES....\002/)";
    static char fmt_205[] = "(\002 DTYPE..... IV(16) =\002,i3)";
    static char fmt_220[] = "(1x,2a4,\002.. V(\002,i2,\002) =\002,d15.7)";
    static char fmt_250[] = "(\0020(INITIAL) JTOL ARRAY...\002/(1x,6d12.3))";
    static char fmt_270[] = "(\0020(INITIAL) D0 ARRAY...\002/1x,6d12.3)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, k, l, m;
    static doublereal vk;
    static integer pu, iv1;
    static char which[4*3];
    static integer jtolp;
    extern /* Subroutine */ int vcopy_(integer *, doublereal *, doublereal *);
    static doublereal machep;
    extern /* Subroutine */ int dfault_(integer *, doublereal *);
    extern doublereal rmdcon_(integer *);

    /* Fortran I/O blocks */
    static cilist io___369 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___372 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___373 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___376 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___379 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___381 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___382 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___383 = { 0, 0, 0, fmt_215, 0 };
    static cilist io___384 = { 0, 0, 0, fmt_205, 0 };
    static cilist io___386 = { 0, 0, 0, fmt_215, 0 };
    static cilist io___387 = { 0, 0, 0, fmt_220, 0 };
    static cilist io___388 = { 0, 0, 0, fmt_250, 0 };
    static cilist io___389 = { 0, 0, 0, fmt_270, 0 };



/*  ***  CHECK NL2SOL (VERSION 2.2) PARAMETERS, PRINT CHANGED VALUES  *** */

/*     DIMENSION IV(*), V(*) */

/* DFAULT -- SUPPLIES DFAULT PARAMETER VALUES. */
/* RMDCON -- RETURNS MACHINE-DEPENDENT CONSTANTS. */
/* VCOPY  -- COPIES ONE VECTOR TO ANOTHER. */

/*  ***  LOCAL VARIABLES  *** */

/* /6 */
/*     REAL CNGD(3), DFLT(3), VN(2,27), WHICH(3) */
/* /7 */
/* / */

/*  ***  IV AND V SUBSCRIPTS  *** */


/* /6 */
/*     DATA NVDFLT/27/, ZERO/0.D+0/ */
/* /7 */
/* / */

/* /6 */
/*     DATA DTYPE/16/, DTYPE0/29/, D0INIT/37/, EPSLON/19/, */
/*    1     INITS/25/, JTINIT/39/, JTOL0/86/, JTOL1/87/, */
/*    2     OLDN/45/, OLDNN/46/, OLDP/47/, PARPRT/20/, */
/*    3     PARSV1/51/, PRUNIT/21/ */
/* /7 */
/* / */

    /* Parameter adjustments */
    --v;
    --iv;

    /* Function Body */
/* /6 */
/*     DATA VN(1,1),VN(2,1)/4HEPSL,4HON../ */
/*     DATA VN(1,2),VN(2,2)/4HPHMN,4HFC../ */
/*     DATA VN(1,3),VN(2,3)/4HPHMX,4HFC../ */
/*     DATA VN(1,4),VN(2,4)/4HDECF,4HAC../ */
/*     DATA VN(1,5),VN(2,5)/4HINCF,4HAC../ */
/*     DATA VN(1,6),VN(2,6)/4HRDFC,4HMN../ */
/*     DATA VN(1,7),VN(2,7)/4HRDFC,4HMX../ */
/*     DATA VN(1,8),VN(2,8)/4HTUNE,4HR1../ */
/*     DATA VN(1,9),VN(2,9)/4HTUNE,4HR2../ */
/*     DATA VN(1,10),VN(2,10)/4HTUNE,4HR3../ */
/*     DATA VN(1,11),VN(2,11)/4HTUNE,4HR4../ */
/*     DATA VN(1,12),VN(2,12)/4HTUNE,4HR5../ */
/*     DATA VN(1,13),VN(2,13)/4HAFCT,4HOL../ */
/*     DATA VN(1,14),VN(2,14)/4HRFCT,4HOL../ */
/*     DATA VN(1,15),VN(2,15)/4HXCTO,4HL.../ */
/*     DATA VN(1,16),VN(2,16)/4HXFTO,4HL.../ */
/*     DATA VN(1,17),VN(2,17)/4HLMAX,4H0.../ */
/*     DATA VN(1,18),VN(2,18)/4HDLTF,4HDJ../ */
/*     DATA VN(1,19),VN(2,19)/4HD0IN,4HIT../ */
/*     DATA VN(1,20),VN(2,20)/4HDINI,4HT.../ */
/*     DATA VN(1,21),VN(2,21)/4HJTIN,4HIT../ */
/*     DATA VN(1,22),VN(2,22)/4HDLTF,4HDC../ */
/*     DATA VN(1,23),VN(2,23)/4HDFAC,4H..../ */
/*     DATA VN(1,24),VN(2,24)/4HRLIM,4HIT../ */
/*     DATA VN(1,25),VN(2,25)/4HCOSM,4HIN../ */
/*     DATA VN(1,26),VN(2,26)/4HDELT,4HA0../ */
/*     DATA VN(1,27),VN(2,27)/4HFUZZ,4H..../ */
/* /7 */
/* / */


/* /6 */
/*     DATA CNGD(1),CNGD(2),CNGD(3)/4H---C,4HHANG,4HED V/, */
/*    1     DFLT(1),DFLT(2),DFLT(3)/4HNOND,4HEFAU,4HLT V/ */
/* /7 */
/* / */

/* ....................................................................... */

    if (iv[1] == 0) {
	dfault_(&iv[1], &v[1]);
    }
    pu = iv[21];
    iv1 = iv[1];
    if (iv1 != 12) {
	goto L30;
    }
    if (*nn >= *n && *n >= *p && *p >= 1) {
	goto L20;
    }
    iv[1] = 16;
    if (pu != 0) {
	io___369.ciunit = pu;
	s_wsfe(&io___369);
	do_fio(&c__1, (char *)&(*nn), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*p), (ftnlen)sizeof(integer));
	e_wsfe();
    }
    goto L999;
L20:
    k = iv[21];
    dfault_(&iv[21], &v[33]);
    iv[21] = k;
    iv[29] = iv[36];
    iv[45] = *n;
    iv[46] = *nn;
    iv[47] = *p;
    s_copy(which, dflt, (ftnlen)4, (ftnlen)4);
    s_copy(which + 4, dflt + 4, (ftnlen)4, (ftnlen)4);
    s_copy(which + 8, dflt + 8, (ftnlen)4, (ftnlen)4);
    goto L80;
L30:
    if (*n == iv[45] && *nn == iv[46] && *p == iv[47]) {
	goto L50;
    }
    iv[1] = 17;
    if (pu != 0) {
	io___372.ciunit = pu;
	s_wsfe(&io___372);
	do_fio(&c__1, (char *)&iv[46], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&iv[45], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&iv[47], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*nn), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*p), (ftnlen)sizeof(integer));
	e_wsfe();
    }
    goto L999;

L50:
    if (iv1 <= 11 && iv1 >= 1) {
	goto L70;
    }
    iv[1] = 50;
    if (pu != 0) {
	io___373.ciunit = pu;
	s_wsfe(&io___373);
	do_fio(&c__1, (char *)&iv1, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    goto L999;

L70:
    s_copy(which, cngd, (ftnlen)4, (ftnlen)4);
    s_copy(which + 4, cngd + 4, (ftnlen)4, (ftnlen)4);
    s_copy(which + 8, cngd + 8, (ftnlen)4, (ftnlen)4);

L80:
    if (big > tiny) {
	goto L90;
    }
    tiny = rmdcon_(&c__1);
    machep = rmdcon_(&c__3);
    big = rmdcon_(&c__6);
    vm[11] = machep;
    vx[11] = big;
    vm[12] = tiny;
    vx[12] = big;
    vm[13] = machep;
    vm[16] = tiny;
    vx[16] = big;
    vm[17] = machep;
    vx[18] = big;
    vx[19] = big;
    vx[20] = big;
    vm[21] = machep;
    vx[23] = rmdcon_(&c__5);
    vm[24] = machep;
    vm[25] = machep;
L90:
    m = 0;
    if (iv[25] >= 0 && iv[25] <= 2) {
	goto L110;
    }
    m = 18;
    if (pu != 0) {
	io___376.ciunit = pu;
	s_wsfe(&io___376);
	do_fio(&c__1, (char *)&iv[25], (ftnlen)sizeof(integer));
	e_wsfe();
    }
L110:
    k = 19;
    for (i__ = 1; i__ <= 27; ++i__) {
	vk = v[k];
	if (vk >= vm[i__ - 1] && vk <= vx[i__ - 1]) {
	    goto L130;
	}
	m = k;
	if (pu != 0) {
	    io___379.ciunit = pu;
	    s_wsfe(&io___379);
	    do_fio(&c__1, vn + ((i__ << 1) - 2 << 2), (ftnlen)4);
	    do_fio(&c__1, vn + ((i__ << 1) - 1 << 2), (ftnlen)4);
	    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&vk, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&vm[i__ - 1], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&vx[i__ - 1], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
L130:
	++k;
/* L140: */
    }

    if (iv1 == 12 && v[39] > 0.) {
	goto L170;
    }

/*  ***  CHECK JTOL VALUES  *** */

    jtolp = *p + 86;
    i__1 = jtolp;
    for (i__ = 87; i__ <= i__1; ++i__) {
	if (v[i__] > 0.) {
	    goto L160;
	}
	k = i__ - 86;
	if (pu != 0) {
	    io___381.ciunit = pu;
	    s_wsfe(&io___381);
	    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&v[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	m = i__;
L160:
	;
    }

L170:
    if (m == 0) {
	goto L180;
    }
    iv[1] = m;
    goto L999;

L180:
    if (pu == 0 || iv[20] == 0) {
	goto L999;
    }
    if (iv1 != 12 || iv[25] == 0) {
	goto L200;
    }
    m = 1;
    io___382.ciunit = pu;
    s_wsfe(&io___382);
    do_fio(&c__1, (char *)&iv[25], (ftnlen)sizeof(integer));
    e_wsfe();
L200:
    if (iv[16] == iv[29]) {
	goto L210;
    }
    if (m == 0) {
	io___383.ciunit = pu;
	s_wsfe(&io___383);
	do_fio(&c__3, which, (ftnlen)4);
	e_wsfe();
    }
    m = 1;
    io___384.ciunit = pu;
    s_wsfe(&io___384);
    do_fio(&c__1, (char *)&iv[16], (ftnlen)sizeof(integer));
    e_wsfe();
L210:
    k = 19;
    l = 51;
    for (i__ = 1; i__ <= 27; ++i__) {
	if (v[k] == v[l]) {
	    goto L230;
	}
	if (m == 0) {
	    io___386.ciunit = pu;
	    s_wsfe(&io___386);
	    do_fio(&c__3, which, (ftnlen)4);
	    e_wsfe();
	}
	m = 1;
	io___387.ciunit = pu;
	s_wsfe(&io___387);
	do_fio(&c__1, vn + ((i__ << 1) - 2 << 2), (ftnlen)4);
	do_fio(&c__1, vn + ((i__ << 1) - 1 << 2), (ftnlen)4);
	do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&v[k], (ftnlen)sizeof(doublereal));
	e_wsfe();
L230:
	++k;
	++l;
/* L240: */
    }
    iv[29] = iv[16];
    vcopy_(&c__27, &v[51], &v[19]);
    if (iv1 != 12) {
	goto L999;
    }
    if (v[39] > 0.) {
	goto L260;
    }
    jtolp = *p + 86;
    io___388.ciunit = pu;
    s_wsfe(&io___388);
    i__1 = jtolp;
    for (i__ = 87; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&v[i__], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
L260:
    if (v[37] > 0.) {
	goto L999;
    }
    k = *p + 87;
    l = k + *p - 1;
    io___389.ciunit = pu;
    s_wsfe(&io___389);
    i__1 = l;
    for (i__ = k; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&v[i__], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();

L999:
    return 0;
/*  ***  LAST CARD OF PARCHK FOLLOWS  *** */
	if (COMMENT)
		printf("parchk: leaving...");
} /* parchk_ */




/* Subroutine */ int qapply_(integer *nn, integer *n, integer *p, doublereal *
	j, doublereal *r__, integer *ierr)
{
    /* System generated locals */
    integer j_dim1, j_offset, i__1, i__2;

    /* Local variables */
    static integer i__, k, l;
    static doublereal t;
    static integer nl1;
    extern doublereal dotprd_(integer *, doublereal *, doublereal *);

/*     *****PARAMETERS. */

/*     .................................................................. */
/*     .................................................................. */

/*     *****PURPOSE. */
/*     THIS SUBROUTINE APPLIES TO R THE ORTHOGONAL TRANSFORMATIONS */
/*     STORED IN J BY QRFACT */

/*     *****PARAMETER DESCRIPTION. */
/*     ON INPUT. */

/*        NN IS THE ROW DIMENSION OF THE MATRIX J AS DECLARED IN */
/*             THE CALLING PROGRAM DIMENSION STATEMENT */

/*        N IS THE NUMBER OF ROWS OF J AND THE SIZE OF THE VECTOR R */

/*        P IS THE NUMBER OF COLUMNS OF J AND THE SIZE OF SIGMA */

/*        J CONTAINS ON AND BELOW ITS DIAGONAL THE COLUMN VECTORS */
/*             U WHICH DETERMINE THE HOUSEHOLDER TRANSFORMATIONS */
/*             IDENT - U*U.TRANSPOSE */

/*        R IS THE RIGHT HAND SIDE VECTOR TO WHICH THE ORTHOGONAL */
/*             TRANSFORMATIONS WILL BE APPLIED */

/*        IERR IF NON-ZERO INDICATES THAT NOT ALL THE TRANSFORMATIONS */
/*             WERE SUCCESSFULLY DETERMINED AND ONLY THE FIRST */
/*             ABS(IERR) - 1 TRANSFORMATIONS WILL BE USED */

/*     ON OUTPUT. */

/*        R HAS BEEN OVERWRITTEN BY ITS TRANSFORMED IMAGE */

/*     *****APPLICATION AND USAGE RESTRICTIONS. */
/*     NONE */

/*     *****ALGORITHM NOTES. */
/*     THE VECTORS U WHICH DETERMINE THE HOUSEHOLDER TRANSFORMATIONS */
/*     ARE NORMALIZED SO THAT THEIR 2-NORM SQUARED IS 2.  THE USE OF */
/*     THESE TRANSFORMATIONS HERE IS IN THE SPIRIT OF (1). */

/*     *****SUBROUTINES AND FUNCTIONS CALLED. */

/*     DOTPRD - FUNCTION, RETURNS THE INNER PRODUCT OF VECTORS */

/*     *****REFERENCES. */
/*     (1) BUSINGER, P. A., AND GOLUB, G. H. (1965), LINEAR LEAST SQUARES */
/*        SOLUTIONS BY HOUSEHOLDER TRANSFORMATIONS, NUMER. MATH. 7, */
/*        PP. 269-276. */

/*     *****HISTORY. */
/*     DESIGNED BY DAVID M. GAY, CODED BY STEPHEN C. PETERS (WINTER 1977) */

/*     *****GENERAL. */

/*     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH */
/*     SUPPORTED BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS */
/*     MCS-7600324, DCR75-10143, 76-14311DSS, AND MCS76-11989. */

/*     .................................................................. */
/*     .................................................................. */

/*     *****LOCAL VARIABLES. */
/*     *****INTRINSIC FUNCTIONS. */
/* /+ */
/* / */
/*     *****FUNCTIONS. */

    /* Parameter adjustments */
    --r__;
    j_dim1 = *nn;
    j_offset = 1 + j_dim1;
    j -= j_offset;

    /* Function Body */
    k = *p;
    if (*ierr != 0) {
	k = abs(*ierr) - 1;
    }
    if (k == 0) {
	goto L999;
    }

    i__1 = k;
    for (l = 1; l <= i__1; ++l) {
	nl1 = *n - l + 1;
	t = -dotprd_(&nl1, &j[l + l * j_dim1], &r__[l]);

	i__2 = *n;
	for (i__ = l; i__ <= i__2; ++i__) {
/* L10: */
	    r__[i__] += t * j[i__ + l * j_dim1];
	}
/* L20: */
    }
L999:
    return 0;
/*     .... LAST CARD OF QAPPLY ......................................... */
} /* qapply_ */

/* Subroutine */ int qrfact_(integer *nm, integer *m, integer *n, doublereal *
	qr, doublereal *alpha, integer *ipivot, integer *ierr, integer *
	nopivk, doublereal *sum)
{
    /* Initialized data */

    static doublereal rktol = 0.;
    static doublereal ufeta = 0.;

    /* System generated locals */
    integer qr_dim1, qr_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, k1, mk1;
    static doublereal beta;
    static integer jbar;
    static doublereal temp, qrkk, sumj, sigma;
    static integer minum;
    extern /* Subroutine */ int vaxpy_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal rktol1;
    extern doublereal v2norm_(integer *, doublereal *);
    static doublereal alphak;
    extern doublereal rmdcon_(integer *), dotprd_(integer *, doublereal *, 
	    doublereal *);
    static doublereal qrkmax;
    extern /* Subroutine */ int vscopy_(integer *, doublereal *, doublereal *)
	    ;


/*  ***  COMPUTE THE QR DECOMPOSITION OF THE MATRIX STORED IN QR  *** */

/*     *****PARAMETERS. */
/*     *****LOCAL VARIABLES. */
/*     *****FUNCTIONS. */
/* /+ */
/* / */
/* DOTPRD... RETURNS INNER PRODUCT OF TWO VECTORS. */
/* RMDCON... RETURNS MACHINE-DEPENDENT CONSTANTS. */
/* VAXPY... COMPUTES SCALAR TIMES ONE VECTOR PLUS ANOTHER. */
/* VSCOPY... SETS ALL ELEMENTS OF A VECTOR TO A SCALAR. */
/* V2NORM... RETURNS THE 2-NORM OF A VECTOR. */

/*     *****CONSTANTS. */
/* /6 */
/*     DATA ONE/1.0D+0/, P01/0.01D+0/, P99/0.99D+0/, ZERO/0.0D+0/ */
/* /7 */
/* / */


/*     .................................................................. */
/*     .................................................................. */


/*     *****PURPOSE. */

/*     THIS SUBROUTINE DOES A QR-DECOMPOSITION ON THE M X N MATRIX QR, */
/*        WITH AN OPTIONALLY MODIFIED COLUMN PIVOTING, AND RETURNS THE */
/*        UPPER TRIANGULAR R-MATRIX, AS WELL AS THE ORTHOGONAL VECTORS */
/*        USED IN THE TRANSFORMATIONS. */

/*     *****PARAMETER DESCRIPTION. */
/*     ON INPUT. */

/*        NM MUST BE SET TO THE ROW DIMENSION OF THE TWO DIMENSIONAL */
/*             ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM */
/*             DIMENSION STATEMENT. */

/*        M MUST BE SET TO THE NUMBER OF ROWS IN THE MATRIX. */

/*        N MUST BE SET TO THE NUMBER OF COLUMNS IN THE MATRIX. */

/*        QR CONTAINS THE REAL RECTANGULAR MATRIX TO BE DECOMPOSED. */

/*     NOPIVK IS USED TO CONTROL PIVOTTING.  COLUMNS 1 THROUGH */
/*        NOPIVK WILL REMAIN FIXED IN POSITION. */

/*        SUM IS USED FOR TEMPORARY STORAGE FOR THE SUBROUTINE. */

/*     ON OUTPUT. */

/*        QR CONTAINS THE NON-DIAGONAL ELEMENTS OF THE R-MATRIX */
/*             IN THE STRICT UPPER TRIANGLE. THE VECTORS U, WHICH */
/*             DEFINE THE HOUSEHOLDER TRANSFORMATIONS   I - U*U-TRANSP, */
/*             ARE IN THE COLUMNS OF THE LOWER TRIANGLE. THESE VECTORS U */
/*             ARE SCALED SO THAT THE SQUARE OF THEIR 2-NORM IS 2.0. */

/*        ALPHA CONTAINS THE DIAGONAL ELEMENTS OF THE R-MATRIX. */

/*        IPIVOT REFLECTS THE COLUMN PIVOTING PERFORMED ON THE INPUT */
/*             MATRIX TO ACCOMPLISH THE DECOMPOSITION. THE J-TH */
/*             ELEMENT OF IPIVOT GIVES THE COLUMN OF THE ORIGINAL */
/*             MATRIX WHICH WAS PIVOTED INTO COLUMN J DURING THE */
/*             DECOMPOSITION. */

/*        IERR IS SET TO. */
/*             0 FOR NORMAL RETURN, */
/*             K IF NO NON-ZERO PIVOT COULD BE FOUND FOR THE K-TH */
/*                  TRANSFORMATION, OR */
/*             -K FOR AN ERROR EXIT ON THE K-TH THANSFORMATION. */
/*             IF AN ERROR EXIT WAS TAKEN, THE FIRST (K - 1) */
/*             TRANSFORMATIONS ARE CORRECT. */


/*     *****APPLICATIONS AND USAGE RESTRICTIONS. */
/*     THIS MAY BE USED WHEN SOLVING LINEAR LEAST-SQUARES PROBLEMS -- */
/*     SEE SUBROUTINE QR1 OF ROSEPACK.  IT IS CALLED FOR THIS PURPOSE */
/*     BY LLSQST IN THE NL2SOL (NONLINEAR LEAST-SQUARES) PACKAGE. */

/*     *****ALGORITHM NOTES. */
/*     THIS VERSION OF QRFACT TRIES TO ELIMINATE THE OCCURRENCE OF */
/*     UNDERFLOWS DURING THE ACCUMULATION OF INNER PRODUCTS.  RKTOL1 */
/*     IS CHOSEN BELOW SO AS TO INSURE THAT DISCARDED TERMS HAVE NO */
/*     EFFECT ON THE COMPUTED TWO-NORMS. */

/*     ADAPTED FROM THE ALGOL ROUTINE SOLVE (1). */

/*     *****REFERENCES. */
/*     (1)     BUSINGER,P. AND GOLUB,G.H., LINEAR LEAST SQUARES */
/*     SOLUTIONS BY HOUSHOLDER TRANSFORMATIONS, IN WILKINSON,J.H. */
/*     AND REINSCH,C.(EDS.), HANDBOOK FOR AUTOMATIC COMPUTATION, */
/*     VOLUME II. LINEAR ALGEBRA, SPRINGER-VERLAG, 111-118 (1971). */
/*     PREPUBLISHED IN NUMER.MATH. 7, 269-276 (1965). */

/*     *****HISTORY. */
/*     THIS AMOUNTS TO THE SUBROUTINE QR1 OF ROSEPACK WITH RKTOL1 USED */
/*     IN PLACE OF RKTOL BELOW, WITH V2NORM USED TO INITIALIZE (AND */
/*     SOMETIMES UPDATE) THE SUM ARRAY, AND WITH CALLS ON DOTPRD AND */
/*     VAXPY IN PLACE OF SOME LOOPS. */

/*     *****GENERAL. */

/*     DEVELOPMENT OF THIS PROGRAM SUPPORTED IN PART BY */
/*     NATIONAL SCIENCE FOUNDATION GRANT GJ-1154X3 AND */
/*     NATIONAL SCIENCE FOUNDATION GRANT DCR75-08802 */
/*     TO NATIONAL BUREAU OF ECONOMIC RESEARCH, INC. */



/*     .................................................................. */
/*     .................................................................. */


/*     ..........  UFETA IS THE SMALLEST POSITIVE FLOATING POINT NUMBER */
/*        S.T. UFETA AND -UFETA CAN BOTH BE REPRESENTED. */

/*     ..........  RKTOL IS THE SQUARE ROOT OF THE RELATIVE PRECISION */
/*        OF FLOATING POINT ARITHMETIC (MACHEP). */
    /* Parameter adjustments */
    --sum;
    --ipivot;
    --alpha;
    qr_dim1 = *nm;
    qr_offset = 1 + qr_dim1;
    qr -= qr_offset;

    /* Function Body */
/*     *****BODY OF PROGRAM. */
    if (ufeta > 0.) {
	goto L10;
    }
    ufeta = rmdcon_(&c__1);
    rktol = rmdcon_(&c__4);
L10:
    *ierr = 0;
    rktol1 = rktol * .01;

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	sum[j] = v2norm_(m, &qr[j * qr_dim1 + 1]);
	ipivot[j] = j;
/* L20: */
    }

    minum = min(*m,*n);

    i__1 = minum;
    for (k = 1; k <= i__1; ++k) {
	mk1 = *m - k + 1;
/*        ..........K-TH HOUSEHOLDER TRANSFORMATION.......... */
	sigma = 0.;
	jbar = 0;
/*        ..........FIND LARGEST COLUMN SUM.......... */
	if (k <= *nopivk) {
	    goto L50;
	}
	i__2 = *n;
	for (j = k; j <= i__2; ++j) {
	    if (sigma >= sum[j]) {
		goto L30;
	    }
	    sigma = sum[j];
	    jbar = j;
L30:
	    ;
	}

	if (jbar == 0) {
	    goto L220;
	}
	if (jbar == k) {
	    goto L50;
	}
/*        ..........COLUMN INTERCHANGE.......... */
	i__ = ipivot[k];
	ipivot[k] = ipivot[jbar];
	ipivot[jbar] = i__;
	sum[jbar] = sum[k];
	sum[k] = sigma;

	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sigma = qr[i__ + k * qr_dim1];
	    qr[i__ + k * qr_dim1] = qr[i__ + jbar * qr_dim1];
	    qr[i__ + jbar * qr_dim1] = sigma;
/* L40: */
	}
/*        ..........END OF COLUMN INTERCHANGE.......... */
L50:
/*        ..........  SECOND INNER PRODUCT  .......... */
	qrkmax = 0.;

	i__2 = *m;
	for (i__ = k; i__ <= i__2; ++i__) {
	    if ((d__1 = qr[i__ + k * qr_dim1], abs(d__1)) > qrkmax) {
		qrkmax = (d__2 = qr[i__ + k * qr_dim1], abs(d__2));
	    }
/* L60: */
	}

	if (qrkmax < ufeta) {
	    goto L210;
	}
	alphak = v2norm_(&mk1, &qr[k + k * qr_dim1]) / qrkmax;
/* Computing 2nd power */
	d__1 = alphak;
	sigma = d__1 * d__1;

/*        ..........  END SECOND INNER PRODUCT  .......... */
	qrkk = qr[k + k * qr_dim1];
	if (qrkk >= 0.) {
	    alphak = -alphak;
	}
	alpha[k] = alphak * qrkmax;
	beta = qrkmax * sqrt(sigma - qrkk * alphak / qrkmax);
	qr[k + k * qr_dim1] = qrkk - alpha[k];
	i__2 = *m;
	for (i__ = k; i__ <= i__2; ++i__) {
/* L65: */
	    qr[i__ + k * qr_dim1] /= beta;
	}
	k1 = k + 1;
	if (k1 > *n) {
	    goto L120;
	}

	i__2 = *n;
	for (j = k1; j <= i__2; ++j) {
	    temp = -dotprd_(&mk1, &qr[k + k * qr_dim1], &qr[k + j * qr_dim1]);

/*             ***  SET QR(I,J) = QR(I,J) + TEMP*QR(I,K), I = K,...,M. */

	    vaxpy_(&mk1, &qr[k + j * qr_dim1], &temp, &qr[k + k * qr_dim1], &
		    qr[k + j * qr_dim1]);

	    if (k1 > *m) {
		goto L110;
	    }
	    sumj = sum[j];
	    if (sumj < ufeta) {
		goto L110;
	    }
	    temp = (d__1 = qr[k + j * qr_dim1] / sumj, abs(d__1));
	    if (temp < rktol1) {
		goto L110;
	    }
	    if (temp >= .99) {
		goto L90;
	    }
/* Computing 2nd power */
	    d__1 = temp;
	    sum[j] = sumj * sqrt(1. - d__1 * d__1);
	    goto L110;
L90:
	    i__3 = *m - k;
	    sum[j] = v2norm_(&i__3, &qr[k1 + j * qr_dim1]);
L110:
	    ;
	}
/*        ..........END OF K-TH HOUSEHOLDER TRANSFORMATION.......... */
L120:
	;
    }

    goto L999;
/*     ..........ERROR EXIT ON K-TH TRANSFORMATION.......... */
L210:
    *ierr = -k;
    goto L230;
/*     ..........NO NON-ZERO ACCEPTABLE PIVOT FOUND.......... */
L220:
    *ierr = k;
L230:
    i__1 = *n;
    for (i__ = k; i__ <= i__1; ++i__) {
	alpha[i__] = 0.;
	if (i__ > k) {
	    i__2 = i__ - k;
	    vscopy_(&i__2, &qr[k + i__ * qr_dim1], &c_b95);
	}
/* L240: */
    }
/*     ..........RETURN TO CALLER.......... */
L999:
    return 0;
/*     ..........LAST CARD OF QRFACT.......... */
} /* qrfact_ */

doublereal reldst_(integer *p, doublereal *d__, doublereal *x, doublereal *x0)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2;

    /* Local variables */
    static integer i__;
    static doublereal t, emax, xmax;


/*  ***  COMPUTE AND RETURN RELATIVE DIFFERENCE BETWEEN X AND X0  *** */
/*  ***  NL2SOL VERSION 2.2  *** */

/* /+ */
/* / */
/* /6 */
/*     DATA ZERO/0.D+0/ */
/* /7 */
/* / */

    /* Parameter adjustments */
    --x0;
    --x;
    --d__;

    /* Function Body */
    emax = 0.;
    xmax = 0.;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t = (d__1 = d__[i__] * (x[i__] - x0[i__]), abs(d__1));
	if (emax < t) {
	    emax = t;
	}
	t = d__[i__] * ((d__1 = x[i__], abs(d__1)) + (d__2 = x0[i__], abs(
		d__2)));
	if (xmax < t) {
	    xmax = t;
	}
/* L10: */
    }
    ret_val = 0.;
    if (xmax > 0.) {
	ret_val = emax / xmax;
    }
/* L999: */
    return ret_val;
/*  ***  LAST CARD OF RELDST FOLLOWS  *** */
} /* reldst_ */

/* Subroutine */ int rptmul_(integer *func, integer *ipivot, doublereal *j, 
	integer *nn, integer *p, doublereal *rd, doublereal *x, doublereal *y,
	 doublereal *z__)
{
    /* System generated locals */
    integer j_dim1, j_offset, i__1, i__2;

    /* Local variables */
    static integer i__, k;
    static doublereal zk;
    static integer im1, km1;
    extern doublereal dotprd_(integer *, doublereal *, doublereal *);


/*  ***  FUNC = 1... SET  Y = RMAT * (PERM**T) * X. */
/*  ***  FUNC = 2... SET  Y = PERM * (RMAT**T) * RMAT * (PERM**T) * X. */
/*  ***  FUNC = 3... SET  Y = PERM * (RMAT**T) X. */


/*  ***  PERM = MATRIX WHOSE I-TH COL. IS THE IPIVOT(I)-TH UNIT VECTOR. */
/*  ***  RMAT IS THE UPPER TRIANGULAR MATRIX WHOSE STRICT UPPER TRIANGLE */
/*  ***       IS STORED IN  J  AND WHOSE DIAGONAL IS STORED IN RD. */
/*  ***  Z IS A SCRATCH VECTOR. */
/*  ***  X AND Y MAY SHARE STORAGE. */


/*  ***  LOCAL VARIABLES  *** */


/*  ***  EXTERNAL FUNCTION  *** */


/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --z__;
    --y;
    --x;
    --rd;
    j_dim1 = *nn;
    j_offset = 1 + j_dim1;
    j -= j_offset;
    --ipivot;

    /* Function Body */
    if (*func > 2) {
	goto L50;
    }

/*  ***  FIRST SET  Z = (PERM**T) * X  *** */

    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = ipivot[i__];
	z__[i__] = x[k];
/* L10: */
    }

/*  ***  NOW SET  Y = RMAT * Z  *** */

    y[1] = z__[1] * rd[1];
    if (*p <= 1) {
	goto L40;
    }
    i__1 = *p;
    for (k = 2; k <= i__1; ++k) {
	km1 = k - 1;
	zk = z__[k];
	i__2 = km1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L20: */
	    y[i__] += j[i__ + k * j_dim1] * zk;
	}
	y[k] = zk * rd[k];
/* L30: */
    }

L40:
    if (*func <= 1) {
	goto L999;
    }
    goto L70;

L50:
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L60: */
	y[i__] = x[i__];
    }

/*  ***  SET  Z = (RMAT**T) * Y  *** */

L70:
    z__[1] = y[1] * rd[1];
    if (*p == 1) {
	goto L90;
    }
    i__1 = *p;
    for (i__ = 2; i__ <= i__1; ++i__) {
	im1 = i__ - 1;
	z__[i__] = y[i__] * rd[i__] + dotprd_(&im1, &j[i__ * j_dim1 + 1], &y[
		1]);
/* L80: */
    }

/*  ***  NOW SET  Y = PERM * Z  *** */

L90:
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = ipivot[i__];
	y[k] = z__[i__];
/* L100: */
    }

L999:
    return 0;
/*  ***  LAST CARD OF RPTMUL FOLLOWS  *** */
} /* rptmul_ */

/* Subroutine */ int slupdt_(doublereal *a, doublereal *cosmin, integer *p, 
	doublereal *size, doublereal *step, doublereal *u, doublereal *w, 
	doublereal *wchmtd, doublereal *wscale, doublereal *y)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j, k;
    static doublereal t, ui, wi;
    extern doublereal v2norm_(integer *, doublereal *);
    static doublereal denmin;
    extern doublereal dotprd_(integer *, doublereal *, doublereal *);
    static doublereal sdotwm;
    extern /* Subroutine */ int slvmul_(integer *, doublereal *, doublereal *,
	     doublereal *);


/*  ***  UPDATE SYMMETRIC  A  SO THAT  A * STEP = Y  *** */
/*  ***  (LOWER TRIANGLE OF  A  STORED ROWWISE       *** */

/*  ***  PARAMETER DECLARATIONS  *** */

/*     DIMENSION A(P*(P+1)/2) */

/*  ***  LOCAL VARIABLES  *** */


/*     ***  CONSTANTS  *** */

/*  ***  INTRINSIC FUNCTIONS  *** */
/* /+ */
/* / */
/*  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  *** */


/* /6 */
/*     DATA HALF/0.5D+0/, ONE/1.D+0/, ZERO/0.D+0/ */
/* /7 */
/* / */

/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --a;
    --y;
    --wchmtd;
    --w;
    --u;
    --step;

    /* Function Body */
    sdotwm = dotprd_(p, &step[1], &wchmtd[1]);
    denmin = *cosmin * v2norm_(p, &step[1]) * v2norm_(p, &wchmtd[1]);
    *wscale = 1.;
    if (denmin != 0.) {
/* Computing MIN */
	d__2 = 1., d__3 = (d__1 = sdotwm / denmin, abs(d__1));
	*wscale = min(d__2,d__3);
    }
    t = 0.;
    if (sdotwm != 0.) {
	t = *wscale / sdotwm;
    }
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	w[i__] = t * wchmtd[i__];
    }
    slvmul_(p, &u[1], &a[1], &step[1]);
    t = (*size * dotprd_(p, &step[1], &u[1]) - dotprd_(p, &step[1], &y[1])) * 
	    .5;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
	u[i__] = t * w[i__] + y[i__] - *size * u[i__];
    }

/*  ***  SET  A = A + U*(W**T) + W*(U**T)  *** */

    k = 1;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ui = u[i__];
	wi = w[i__];
	i__2 = i__;
	for (j = 1; j <= i__2; ++j) {
	    a[k] = *size * a[k] + ui * w[j] + wi * u[j];
	    ++k;
/* L30: */
	}
/* L40: */
    }

/* L999: */
    return 0;
/*  ***  LAST CARD OF SLUPDT FOLLOWS  *** */
} /* slupdt_ */

/* Subroutine */ int slvmul_(integer *p, doublereal *y, doublereal *s, 
	doublereal *x)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k;
    static doublereal xi;
    static integer im1;
    extern doublereal dotprd_(integer *, doublereal *, doublereal *);


/*  ***  SET  Y = S * X,  S = P X P SYMMETRIC MATRIX.  *** */
/*  ***  LOWER TRIANGLE OF  S  STORED ROWWISE.         *** */

/*  ***  PARAMETER DECLARATIONS  *** */

/*     DIMENSION S(P*(P+1)/2) */

/*  ***  LOCAL VARIABLES  *** */


/*  ***  NO INTRINSIC FUNCTIONS  *** */

/*  ***  EXTERNAL FUNCTION  *** */


/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --x;
    --y;
    --s;

    /* Function Body */
    j = 1;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] = dotprd_(&i__, &s[j], &x[1]);
	j += i__;
/* L10: */
    }

    if (*p <= 1) {
	goto L999;
    }
    j = 1;
    i__1 = *p;
    for (i__ = 2; i__ <= i__1; ++i__) {
	xi = x[i__];
	im1 = i__ - 1;
	++j;
	i__2 = im1;
	for (k = 1; k <= i__2; ++k) {
	    y[k] += s[j] * xi;
	    ++j;
/* L30: */
	}
/* L40: */
    }

L999:
    return 0;
/*  ***  LAST CARD OF SLVMUL FOLLOWS  *** */
} /* slvmul_ */

logical stopx_(integer *idummy)
{
    /* System generated locals */
    logical ret_val;

/*     *****PARAMETERS... */

/*     .................................................................. */

/*     *****PURPOSE... */
/*     THIS FUNCTION MAY SERVE AS THE STOPX (ASYNCHRONOUS INTERRUPTION) */
/*     FUNCTION FOR THE NL2SOL (NONLINEAR LEAST-SQUARES) PACKAGE AT */
/*     THOSE INSTALLATIONS WHICH DO NOT WISH TO IMPLEMENT A */
/*     DYNAMIC STOPX. */

/*     *****ALGORITHM NOTES... */
/*     AT INSTALLATIONS WHERE THE NL2SOL SYSTEM IS USED */
/*     INTERACTIVELY, THIS DUMMY STOPX SHOULD BE REPLACED BY A */
/*     FUNCTION THAT RETURNS .TRUE. IF AND ONLY IF THE INTERRUPT */
/*     (BREAK) KEY HAS BEEN PRESSED SINCE THE LAST CALL ON STOPX. */

/*     .................................................................. */

    ret_val = FALSE_;
    return ret_val;
} /* stopx_ */

/* Subroutine */ int vaxpy_(integer *p, doublereal *w, doublereal *a, 
	doublereal *x, doublereal *y)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;


/*  ***  SET W = A*X + Y  --  W, X, Y = P-VECTORS, A = SCALAR  *** */



    /* Parameter adjustments */
    --y;
    --x;
    --w;

    /* Function Body */
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	w[i__] = *a * x[i__] + y[i__];
    }
    return 0;
} /* vaxpy_ */

/* Subroutine */ int vcopy_(integer *p, doublereal *y, doublereal *x)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;


/*  ***  SET Y = X, WHERE X AND Y ARE P-VECTORS  *** */



    /* Parameter adjustments */
    --x;
    --y;

    /* Function Body */
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	y[i__] = x[i__];
    }
    return 0;
} /* vcopy_ */

/* Subroutine */ int vscopy_(integer *p, doublereal *y, doublereal *s)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;


/*  ***  SET P-VECTOR Y TO SCALAR S  *** */



    /* Parameter adjustments */
    --y;

    /* Function Body */
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	y[i__] = *s;
    }
    return 0;
} /* vscopy_ */

doublereal v2norm_(integer *p, doublereal *x)
{
    /* Initialized data */

    static doublereal sqteta = 0.;

    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal r__, t, xi, scale;
    extern doublereal rmdcon_(integer *);


/*  ***  RETURN THE 2-NORM OF THE P-VECTOR X, TAKING  *** */
/*  ***  CARE TO AVOID THE MOST LIKELY UNDERFLOWS.    *** */


/* /+ */
/* / */

/* /6 */
/*     DATA ONE/1.D+0/, ZERO/0.D+0/ */
/* /7 */
/* / */
    /* Parameter adjustments */
    --x;

    /* Function Body */

    if (*p > 0) {
	goto L10;
    }
    ret_val = 0.;
    goto L999;
L10:
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (x[i__] != 0.) {
	    goto L30;
	}
/* L20: */
    }
    ret_val = 0.;
    goto L999;

L30:
    scale = (d__1 = x[i__], abs(d__1));
    if (i__ < *p) {
	goto L40;
    }
    ret_val = scale;
    goto L999;
L40:
    t = 1.;
    if (sqteta == 0.) {
	sqteta = rmdcon_(&c__2);
    }

/*     ***  SQTETA IS (SLIGHTLY LARGER THAN) THE SQUARE ROOT OF THE */
/*     ***  SMALLEST POSITIVE FLOATING POINT NUMBER ON THE MACHINE. */
/*     ***  THE TESTS INVOLVING SQTETA ARE DONE TO PREVENT UNDERFLOWS. */

    j = i__ + 1;
    i__1 = *p;
    for (i__ = j; i__ <= i__1; ++i__) {
	xi = (d__1 = x[i__], abs(d__1));
	if (xi > scale) {
	    goto L50;
	}
	r__ = xi / scale;
	if (r__ > sqteta) {
	    t += r__ * r__;
	}
	goto L60;
L50:
	r__ = scale / xi;
	if (r__ <= sqteta) {
	    r__ = 0.;
	}
	t = t * r__ * r__ + 1.;
	scale = xi;
L60:
	;
    }

    ret_val = scale * sqrt(t);
L999:
    return ret_val;
/*  ***  LAST CARD OF V2NORM FOLLOWS  *** */
} /* v2norm_ */

#ifdef __cplusplus
	}
#endif
