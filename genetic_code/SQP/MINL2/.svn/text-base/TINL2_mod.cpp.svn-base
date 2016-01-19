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





// these three lines tells the compiler that the following code is written in C (porting C into C++)

#ifdef __cplusplus
extern "C" {
#endif
static Population *Pop; //fdf should get Pop without passing it through MINL2- it works!
static int ntree;  //try  - the same as before - it works!

int opti_cpp(Population *c_Pop, int *method, int n_unknown, int n_fun, int c_ntree, double *x)
{	
	Pop = c_Pop;
	ntree = c_ntree;
	Val t;
	char* expr;

	int COMMENT = 0; //1 if comment on the screen, 0 if silent...
	//if (COMMENT) 
		printf("\n\nEntered opti_");
	
	integer c__1 = 1;    
	integer n = integer(n_unknown);  // n: number of unknowns (parameters) 
	integer m = integer(n_fun);  // m: number of functions (terms of the sum)
	integer c__195 = m*(n+3)+93+n*(n*3+33)/2+n*m+m;//size if the workspace array IW (see MINL2, first lines...)
	// in fortran version is used c__195 = IW=M*(2*N+4)+(N/2+1)*(3*N+33)+93)
	if (COMMENT) {
		//verifies the values
		*(v_list[0]->var) = z_test[0];
		t = (Val)(Pop->tree_value(ntree));
		printf ("\nValue of the tree in z1=-10 :  %f",t);
		printf("\nopti_ : ntree = %i", ntree);
		printf("\nopti_ : P = %i", Pop);
		printf("\nopti_ : n = %i", n); 
		printf("\nopti_ : m = %i", m);
		printf("\nopti_ : c__195 = %i", c__195);
	}

	/* Format strings */
    char fmt_10[] = "(\002 INPUT ERROR. PARAMETER NUMBER \002,i1,\002 IS OUTSIDE ITS RANGE.\002)";
    char fmt_21[] = "(\002 TEST OF GRADIENTS \002//\002 MAXIMUM FORWA\
RD  DIFFERENCE: \002,d8.2,\002 AT FUNCTION NO \002,i1,\002 AND VARIABLE NO\
 \002,i1/\002 MAXIMUM BACKWARD DIFFERENCE: \002,d8.2,\002 AT FUNCTION NO \
\002,i1,\002 AND VARIABLE NO \002,i1/\002 MAXIMUM CENTRAL  DIFFERENCE: \002,\
d8.2,\002 AT FUNCTION NO \002,i1,\002 AND VARIABLE NO \002,i1)";
    char fmt_22[] = "(//\002 MAXIMUM ELEMENT IN DF: \002,d8.2)";
    char fmt_30[] = "(\002 UPPER LIMIT FOR FUNCTION EVALUATIONS EXCEEDED.\002/)";
    char fmt_31[] = "(23x,\002 SOLUTION: \002,d18.10/d52.10//\002 NUM\
BER OF CALLS OF FDF: \002,i4//\002 FUNCTION VALUES AT THE SOLUTION: \002,d18.10,2(/d52.10))";

	 /* Fortran I/O blocks */
   	cilist io___8 = { 0, 6, 0, fmt_10, 0 };
   	cilist io___11 = { 0, 6, 0, fmt_21, 0 };
   	cilist io___12 = { 0, 6, 0, fmt_22, 0 };
   	cilist io___13 = { 0, 6, 0, fmt_30, 0 };


    /* System generated locals */
    integer i__1;

    // Builtin functions
   	integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();

    // Local variables
  	integer i__, j, k;  
	
	/*
	// Array x (parameters)  
	doublereal * x ;    
	x = new (nothrow) doublereal [n];
	// verify if the the "new" has been succesful
	if (x == 0)
		cout << "\nError : memory can't be allocated (new (nothrow) doublereal x [n]  failed)";
	else 
		if (COMMENT) 
			printf("\nopti_ : new (nothrow) doublereal [n] - OK : pointer = %i", x);
	*/

	// Array w (workspace)  
	doublereal * w;
	w = new (nothrow) doublereal [c__195];
	if (w == 0) 
		cout << "\nError : memory can't be allocated (new (nothrow) doublereal [c_195] )";
	else
		if (COMMENT) 
			printf("\nopti_ : new (nothrow) doublereal [c_195] - OK : pointer = %i", w);
	
	doublereal dx;
   	doublereal eps;
 	integer index[8]	/* was [4][2] */;  //?????
 	logical optim;  
	integer icontr, maxfun; 

	// Initial guess of parameters
	if (COMMENT) printf("\nopti_ : ");
	for (int i=0; i<n; i++) {
		//x[i]=10.;	
		if (COMMENT) printf (" x[%i] = %f", i, x[i]);
	}

	if (COMMENT) {
		expr = Pop->print(ntree);
		printf ("\nopti_ : Tree before updating parameters %i :   %s",ntree,expr);
		//verifies the value
		*(v_list[0]->var) = z_test[1];
		t = (Val)(Pop->tree_value(ntree));
		printf ("\nValue of the tree in z1=-9 :  %f",t);
	}

	

	//     SET PARAMETERS
    eps = 1e-10;   //tolerance
    maxfun = 50;  //max number of function calls   25;
    icontr = *method; 
    optim = icontr > 0;
	if (! optim) {dx = .001;}	//original choice 	dx = .001;}
	
	
	// call to minl2_ (in MINL2.cpp, translated from fortran)
	printf("\nopti_ : call to MINL2");
    minl2_(	(U_fp)(fdf_),					
		&n,	//n : no. of unknown parameters - dimension of x
		&m,	// m : no. of functions, terms of the sum
		x,	// x : initial guess of the parameters 
		&dx,	// dx : damping factor (10E-4 if initial guess close to solution, otherwise 1)
		&eps,	// eps : accuracy
		&maxfun,// maxfun : max no. of calls to the function
		w,	// workspace
		&c__195,// dimension of the workspace: don't touch it	
		&icontr); // icontr : used to control the function
		//Pop,		// Pop : pointer to class Population
		//ntree); // ntree: reference number of the tree
	

	printf("\n");

    if (icontr < 0) {
	goto L100;
    }
    if (! optim) {
	goto L200;
    }
    goto L300;

/*     PARAMETER OUTSIDE RANGE */
L100:
    printf("COMPUTATION DID NOT START");
	s_wsfe(&io___8);
    i__1 = -icontr;
    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	printf("\nICONTR = %i", icontr);    
	e_wsfe();
    goto L400;

/*     RESULTS FROM GRADIENT TEST */
L200:
    printf("\nResults from gradient test");
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
		printf("ITERATION STOPPED (TOO MANY ITERATIONS NEEDED OR ROUNDING ERRORS DOMINATE)");
		s_wsfe(&io___13);
		e_wsfe();
    }
    // this cycle controls the independent variables x
	printf ("RESULTS FROM THE OPTIMIZATION");
	for (int i = 0; i < n; i++) {       //n is the number of independent variables
		printf("\nopti_ : Value of the %i -th solution: x(%i) = %f", i, i, x[i] );
	}
    printf("\nopti_ : Number of calls: %i", maxfun);
	
	// this cycle controls the dependent functions f   
	if (COMMENT) 
		for (int j = 0; j < m; j++) {      //m is the number of functions the F sum is made of
			printf("\nopti_ :  g-t(%i) = w(%i) = %E", j, j, w[j]);
		}
	
	//total error
	doublereal sum;   
	for (j=0; j<m; j++) 
		sum=sum+w[j]*w[j];	
	//value of .5*S(fi^2)
	printf("\nopti_ : Value of half the sum of the squares   F: %E",  .5*sum);
 	
L400:
	printf("\nSize of w = %i", sizeof w);
	printf(" ICONTR = %i", icontr);
	if (!icontr)
		printf("Successful call");
	

	
	delete[] w;
    return 0;

}

#ifdef __cplusplus
	}
#endif



#ifdef __cplusplus     //if compiled with g++ it says the compiler to consider the code as C code
extern "C" {
#endif
// Subroutine fdf : this has to be defined by the user. It is called many times by MINL2 (see number of calls).
// THE FOLLOWING INPUT CANNOT BE CHANGED!!
// n : no. of unknown parameters - dimension of x (in)
// m : no. of functions, terms of the sum (in)
// x : initial guess of the parameters (in)
// df : vector containing the values of the derivatives in x (out)
// f : vector containing the values of the functions in x (out)
int fdf_(int *n, int *m, double *x, double *df, double *f) //Population *Pop, int ntree
{	
	char *expr;
	int COMMENT =0; //1 if comments on the screen, 0 if silent...
	
	//if (COMMENT) {
		printf("\n\nEntered fdf_ ");  
		// just to check ntree and pointer to Population are right
		printf("\nfdf_ : ntree = %i", ntree);
		printf("\nfdf_ : Pop = %i", Pop);
	//}

	// System generated locals 
    int df_dim1, df_offset;
    double d__1;
	
    // Parameter adjustments 
    --x;
    --f;
    df_dim1 = *m;   //
    df_offset = 1 + df_dim1;
    df -= df_offset;

	// variables used in differentiation
	double x_curr;
	double der_ij;
	double sum_sq_error = 0.;    
	
	//Function Body : the first element on the array is 1!!!! NOT 0!
    // update the fitness function term with the new parameters values
	double *x_shift = new (nothrow) double [*n]; 
	if (x_shift == 0)
		if (COMMENT)
			cout << "\nError : memory can't be allocated (new (nothrow) double x_shift [*n]  failed)";
	else 
		printf("\nfdf_ : double* x_shift = new (nothrow) doublereal [n] - OK : pointer = %i", x);
	
	// update the original tree with the newly computed parameters 
	// (shift of the entries of the array required when passing parameters (inside fdf it starts from 1))
	for (int i=0; i<*n; i++)
		x_shift[i] = x[i+1];
	Pop->update(ntree, x_shift, *n);
	// verify (just a check that the parameters have been updated in the actual tree)
	expr = Pop->print(ntree);
	printf ("\nfdf_ : Tree %i :   %s",ntree,expr);
	
	if (COMMENT) {
		//verifies the value
		*(v_list[0]->var) = z_test[1];
		Val tr = (Val)(Pop->tree_value(ntree));
		printf ("\nfdf_ : Value of the tree in z1=-9 :  %f",tr);
	}

	// functions
	Val g, t;
	for (int i=1; i<*m+1; i++) { // counter for the function
		
		// functions	
		g =(g_obj(z_test[i-1]));
		*(v_list[0]->var) = z_test[i-1];   //assigns  the values to the independent variables (NOT parameters!) 
		t =(Val)(Pop->tree_value(ntree));  
		f[i] = g-t;
		
		if (COMMENT) {
		//check values on the screen
		expr = Pop->print(ntree);
		printf("\nfdf_ : function number %i  corrisponding to test point   z[%i] = %f",i,i,z_test[i-1]);
		printf ("\nfdf_ : Tree %i :   %s",ntree,expr);
		printf("\nfdf_ : g = %E   t = %E   g-t = %E", g, t, f[i]);
		}

		//compute a term of fitness value
		sum_sq_error = sum_sq_error + .5*(f[i]*f[i]); //refers to the precedent set of parameters 
		
		// derivatives		
		for (int j=1; j<*n+1; j++)   { // counter for the parameter
			x_curr = x[j];
			der_ij = Pop->jacobian_ij(x_curr, ntree, j);
			if (j==1)
				df[df_dim1 + i] = der_ij;   // derivative respect to the first variable
			if (j>1)
				df[(df_dim1 << (j-1)) +i] = der_ij;   //<------------------------
			
			if (COMMENT) printf ("\nfdf_ : df %i / dx %i = %f ", i, j, der_ij);  
		}
	}
	
	

	if (COMMENT) {
		printf("\nfdf_ : Fitness = %E", sum_sq_error);
		for (int h=1; h<*n+1; h++)	
			printf("  x[%i]=%f", h, x[h]);
	}
	delete [] x_shift;
    return 0;
}


#ifdef __cplusplus
	}
#endif



		
