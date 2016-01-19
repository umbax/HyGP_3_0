#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
extern "C" {
int  fdftinl2_c__ (double *f,  double *df, double *x, int *m, int *n)  
{ 

	int COMMENT = 1; //1 if comments on the screen, 0 if silent...
	static int times = 1;
	int i, mloc, nloc, nvar;
  	double z1;
	char *expr;
	Val g, t;
	int ntree = Pop->ntree_fdf_c;   //really ugly way to read a variable...Find another one!
	
	if (COMMENT) {
	cout << "\n\nfdf_c : ntree = " << ntree << "   (call n. " << times << ")";
	cout << "\nfdf_c : Pop = " << Pop ;
	cout << "\nfdf_c : *n = " << *n ;
	cout << "\nfdf_c : *m = " << *m ;
	}

	mloc =*m;  // no. of function to minimize (related to no. of test point or fitness cases or records)
  	nloc = *n; // no. of parameters to optimize
	z1=0;
	nvar = Pop->num_vars;

	// System generated locals 
    int df_dim1, df_offset;
    double d__1;
	
    // Parameter adjustments - OK
  	--x;  //?
 	--f;   //?
   df_dim1 = mloc;   //number of rows of df
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
	// x : array from 1 (FORTRAN)
	// x_shift : array from 0 (C++)
	// (shift of the entries of the array required when passing parameters (inside fdf it starts from 1))
	for (int i=0; i<nloc; i++) {   
		x_shift[i] = x[i+1];
		if (COMMENT) {
		printf("\nx[%i] = %f ", i,x[i]);
		printf("\nx_shift[%i] = %f", i, x_shift[i]);
		}
	}
	if (COMMENT) printf("\nx[%i] = %f ", nloc,x[nloc]);
	Pop->update(ntree, x_shift, nloc);

	// verify (just a check that the parameters have been updated in the actual tree)
	expr = Pop->print(ntree,NULL);
	if (COMMENT)
		cout << "\nfdf_ c : Tree " << ntree << " : " << expr << endl;
	
	// functions
	for (int i=0; i<mloc; i++) { // counter for the function (each i is a fitness case... mloc = n_test_cases_tune)
		
		// functions	
		g  = Pop->data_tune[i][nvar];    // value of the output
		for (int j=0; j<nvar; j++) {			// assign the right value to all the variables for the i-th fitness case
			*(v_list[j]->var)=Pop->data_tune[i][j];   //*((Pop->v_list[j])->var)  a primo membro?
		}
		t =(Val)(Pop->tree_value(ntree));  //compute the value of the tree for the point data_tune[i][j]
		f[i+1] = g-t;
		
		if (COMMENT) {
		//check values on the screen
		expr = Pop->print(ntree, NULL);   //trick not to write Pop->complete_trees,  that is private...
		printf("\nfdf_c: fdf_c called %i times",times);
		printf("\nfdf_ c: function number %i  corresponding to output n.%i = %f",i,i,Pop->data[i][nvar]);
		printf ("\nfdf_c : Tree %i :   %s",ntree,expr);
		printf("\nfdf_ c: g = %E   t = %E   g-t = %E", g, t, f[i+1]);
		}


		//compute a term of fitness value
		sum_sq_error = sum_sq_error +  .5*(f[i+1]*f[i+1]); //refers to the precedent set of parameters 
		
		// derivatives - (columns of the same row i+1)
		for (int j=1; j<nloc+1; j++)   { // counter for the parameter    OK  from 1 to nloc+1
			x_curr = x[j];
			der_ij = Pop->jacobian_ij(x_curr, ntree, j);

			// unfortunately matrix notation df[i][j] can't be used because the matrix fdf was defined in FORTRAN...
			// REMEMBER!!! FORTRAN allocate memory in arrays column after column!
			df[((df_dim1)*(j)) + i + 1] = der_ij;   
			
			
	
			if (COMMENT){
				printf ("\nfdf_ : der_ij(%i)(%i) = %f ", i+1, j, der_ij); 
				printf ("   df[%i, %i] = %f ", i+1, j, df[(df_dim1*(j)) + i + 1] );  
			}
		}
	}
	
	

	if (COMMENT) {
		printf("\nfdf_ : Fitness = %E", sum_sq_error);
		for (int h=1; h<nloc+1; h++)	
			printf("  x[%i]=%f", h, x[h]);
	}
	cout << endl;
	
	times++;




	delete [] x_shift;
    return (1);
}
}
