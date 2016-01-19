#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
extern "C" {
int  fdf_c__ (double *f,  double *x, int *m, int *n)  
{ 
	
	int COMMENT = 0; //if 1 comments, if 0 silent...	
  	static int times = 1;
	int i, mloc, nloc;
  	double z1 = 0.;
	double pen, e_arg;	
	char *expr;
	int nvar = Pop->parameters.nvar;
	//----- can be simplified --------
   	int ntree = Pop->ntree_fdf_c;   //really ugly way to read a variable...Find another one!
   	Binary_Node *tree = Pop->complete_trees[ntree];   //address of the complete tree to be optimised
	//---------------------------------------
   	Val g,t;
	double sum_sq_error = 0.;    
	
	mloc =*m - tree->n_pulsations;  // no. of "distance" terms to minimise (related to no. fitness cases or records)
	nloc = *n; // no. of parameters to optimise

	if (COMMENT) {
		cout << "\n\nfdf_c : ntree = " << ntree << "   (call n. " << times << ")";
		cout << "\nfdf_c : Pop = " << Pop ;
		cout << "\nfdf_c : *n = " << *n ;
		cout << "\nfdf_c : *m = " << *m ;
		cout << "\nfdf_c__ : n_puls = " << tree->n_pulsations;	
		cout << "\n fdf_c__ : mloc = " << mloc; 
	}
	
	//---------------------------------------------------------------------------------------------------------------------
    // fetch the parameters from x 
	// (shift really needed?? in FORTRAN the first element on the array is 1!!!! NOT 0!)
	//---------------------------------------------------------------------------------------------------------------------
	double *x_shift = new (nothrow) double [nloc]; 
	if (x_shift == 0) {
			cerr << "\nfdf_c__ : ERROR : memory can't be allocated (new (nothrow) double x_shift [*n]  failed)";
			exit (-1);
	}
	else
		if (COMMENT)
			cout << "\nfdf_ : double* x_shift = new (nothrow) doublereal [n] - OK : pointer = " << x;
	
	// update the original tree with the newly computed parameters 
	// (shift of the entries of the array required when passing parameters (INSIDE fdf it starts from 1))
	for (int i=0; i<nloc ; i++) {
		x_shift[i] = x[i];
		if (COMMENT) {
		printf("\nx[%i] = %f ", i,x[i]);
		printf("\nx_shift[%i] = %f", i, x_shift[i]);
		}
	}
	Pop->update_complete_tree(tree, x_shift, nloc);

	// verify (just a check that the parameters have been updated in the actual tree)
	if (COMMENT) {
		expr = Pop->print(ntree, NULL);   //trick not to write Pop->complete_trees,  that is private...
		cout << "\nfdf_ c : Tree " << ntree << " : " << expr << endl;
		cout << "\nData_tune : output" << endl;
		delete [] expr;
	}

	//----------------------------------------------------------------------------------------------------------------------------------
	// DEFINITION OF THE SINGLE TERM OF THE SUM (difference between tree output and actual output)
	//----------------------------------------------------------------------------------------------------------------------------------	
	for (i = 0; i <  mloc; i++) {				// each i is a fitness case... mloc = n_test_cases_tune    
	// functions	
		g = Pop->problem.data_tuning[i][nvar];    // value of the output
		
		if (COMMENT) 
			cout << "\n" << i << ") ";
		for (int j=0; j<nvar; j++) {			// assign the right value to all the variables for the i-th fitness case
			(*Pop->problem.v_list[j]).value = Pop->problem.data_tuning[i][j];
			if (COMMENT)
				cout << (*Pop->problem.v_list[j]).value << "  ";
		}
	
		// compute tree value for point data_tune[i][j]
		t =(Val)(Pop->tree_value(tree, NULL));
		
		//------------------------------------------------
		// definition of the term used by MI0L2 		
		//------------------------------------------------
		f[i] = g-t;  //so far : it works perfectly
/*
		if (abs(g)<1.0e-12)   //for normalised RMSE
			f[i] = g-t;
		else 
			f[i] = abs((g-t)/g);
//*/		
		if (COMMENT) {
		//check values on the screen
			expr = Pop->print(ntree, NULL);   //trick not to write Pop->complete_trees,  that is private...
			printf("\nfdf_ c: function number %i  corrisponding to output %i = %f",i,i,Pop->problem.data_tuning[i][Pop->problem.get_n_cols() - 1]);
			printf ("\nfdf_c : Tree %i :   %s",ntree,expr);
			printf("\nfdf_ c: g = %E   t = %E   g-t = %E", g, t, f[i]);
		}

		//compute a term of fitness value
		sum_sq_error = sum_sq_error + .5*(f[i]*f[i]); //refers to the precedent set of parameters 
	
    }

	//-------------------------------------------------------------------------
	// definition of extra penalisations (i.e. for pulsation)
	//-------------------------------------------------------------------------
	for (int i=0; i< tree->n_pulsations; i++) {
		// compare the pulsation to the corresponding omega_lim		
		//if (abs(x[Pop->index_puls[i]]) <= Pop->omega_lim[i])  //omega_lim=2pi*N_oscillations/Zk_ampl		
		if (abs(x[tree->index_puls[i]]) <= (Pop->problem.v_list[tree->index_var[i]])->omega_lim)  //omega_lim=2pi*N_oscillations/Zk_ampl
			f[mloc+i] = 0.0;
		else {
			//e_arg = abs(x[Pop->index_puls[i]]) - Pop->omega_lim[i];
			e_arg = abs(x[tree->index_puls[i]]) - (Pop->problem.v_list[tree->index_var[i]])->omega_lim;
			pen = exp(e_arg*e_arg) - 1.0;    //in ASMO UK paper I forgot - 1.0!!!

			if (pen<=MAX_VAL)
				f[mloc+i] = pen;
			else
				f[mloc+i] = MAX_VAL;
		}
		if (COMMENT) {	
			cout << "\n" << mloc+i << ") penalisation term = " << f[mloc+i] << endl;
			cout << "variable : " << (Pop->problem.v_list[tree->index_var[i]])->name << endl;
			cout << "omega_lim = " << (*Pop->problem.v_list[tree->index_var[i]]).omega_lim;
			cout << " pulsation = " <<  x[tree->index_puls[i]] << endl;
		}

	}

	// add other terms here (say, weight decay...)

	////

	times++;
    delete[] x_shift;

	return(1);
}
}  //extern C ...
