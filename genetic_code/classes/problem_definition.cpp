// Copyright 2016 Dr Umberto Armani
//  
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//  
//      http://www.apache.org/licenses/LICENSE-2.0
//  
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


/*
 * problem_definition.cpp
 *
 *  Created on: Jan 31, 2011
 *      Author: cnua
 */

using namespace std;

#include "./problem_definition.h"


// ProblemDefinition constructor
ProblemDefinition::ProblemDefinition(void)
{
	// private members
	data = NULL;			//original data provided by the user (NOT to be touched)
	n_data = -1;			//total number of rows in data
	n_var = -1;			// number of variables
	n_cols = -1;			//number of columns in data (n_var+1)
	n_folds=0;
	folds_table=NULL;
	points_per_fold=NULL;
	validation_fold=0;

	// public members
	//folds = NULL; // Val***
	data_tuning = NULL;
	n_tuning = -1;
	data_validation = NULL;
	n_validation = -1;

	// statistics of corresponding output (target) of the WHOLE data set (DATA) contained in the main input file (refers to the whole data set - tuning+evaluation)
	sum_output = 0.0;	// sum of the squares of each target value
	y_ave = 0.0;		// average value of input (target) data
	Sy = 0.0;			// SStot total sum of squares of (observed data - average observed data) // defined in read_input_file function (read_file_new.cpp)
	y_var = 0.0;		// variance of target data
	y_max = 0.0;			// max value of target
	index_max=-1;
	y_min = 0.0;			// min value of target
	index_min=-1;
	first_acf_root_input = 0.0; // closest root to 0 of autocorrelation function of input data - only for n_var=1
	tot_variation_input = 0.0; // total variation

	// "Nyquist" variable
	Ny_omega_max = 0.0;

	// statistics of corresponding output (target) of the TEST data set (TEST DATA SET)
	data_test = NULL;
	n_test = -1;
	sum_output_test = 0.0;
	y_ave_test = 0.0;
	Sy_test = 0.0;
	y_var_test = 0.0;
	y_test_max = 0.0;
	y_test_min = 0.0;

	// inequality constraints on values (order 0)
	data_inequality0 = NULL;
	n_inequality0 = -1;
	constraints0 = NULL;

	// inequality constraints on first derivatives (order 1)
	data_inequality1 = NULL;
	n_inequality1 = -1;
	constraints1 = NULL;

	// variables initialisation
	variables_initialised = 0;
	symbol = "not_initialised";      //letter used for the variables
	v_list = NULL;    //list of variables addresses

	//Andrey's idea to avoid dynamic allocation of functions list
	// in the future use a vector
	for (int k=0; k<15; k++ ) {
		//dummy_uni[k] = NULL;
		u_func_list[k] = NULL;
	}
	num_u_funcs = -1;

	//Andrey's idea to avoid dynamic allocation of functions list
	// in the future use a vector
	for (int k=0; k<15; k++ ) {
		//dummy_bin[k] = NULL;
		b_func_list[k] = NULL;
	}
	num_b_funcs = -1;
	division = NULL;

	// input data autocorrelation
	delay_max=0;
	r_k=NULL;
}

//// ProblemDefinition copy constructor (used for parallelisation - see firstprivate)
//ProblemDefinition::ProblemDefinition(const ProblemDefinition& p)
//{
//	int COMMENT = 0;
//
//# pragma omp critical
//	{
//	// here the instructions to make a deep copy
//	if (COMMENT) cout << "\n\ncopy constructor entered" << endl;
//
//	// COPY PRIVATE MEMBERS
//	// copy n_data, n_var, n_cols, data
//	n_data = p.n_data;			//total number of rows in data
//	n_var = p.n_var;			// number of variables
//	n_cols = p.n_cols;			//number of columns in data (n_var+1)
//	data = new Val*[n_data];
//	if (data==NULL)  {
//		cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of data failed!! Input data can't be imported" << endl;
//		exit(-1);
//	}
//	for (int i=0; i<n_data; i++) {
//		data[i] = new Val[n_var+1];
//		if (data[i]==NULL)  {
//			cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of data[" << i << "]  failed!! Input data can't be imported" << endl;
//			exit(-1);
//		}
//	}
//	for (int i=0; i<n_data; i++)
//		for (int j=0; j<n_var+1; j++)
//			data[i][j]= p.data[i][j];
//
//	// copy n_folds, validation_fold, folds_table (matrix n_datax2), points_per_fold (n_foldsx1)
//	n_folds = p.n_folds;
//	validation_fold = p.validation_fold;
//	folds_table = new int*[n_data];
//	if (folds_table==NULL)  {
//		cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of folds_table failed!!" << endl;
//		exit(-1);
//	}
//	for (int i=0; i<n_data; i++) {
//		folds_table[i] = new int[2];
//		if (folds_table[i]==NULL)  {
//			cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of folds_table[" << i << "]  failed!!" << endl;
//			exit(-1);
//		}
//	}
//	for (int i=0; i<n_data; i++)
//		for (int j=0; j<2; j++)
//			folds_table[i][j]= p.folds_table[i][j];
//
//	points_per_fold = new int[n_folds];
//	if (points_per_fold==NULL)  {
//		cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of points_per_fold failed!!" << endl;
//		exit(-1);
//	}
//	for (int i=0; i<n_folds; i++)
//		points_per_fold[i]= p.points_per_fold[i];
//
//
//
//
//	// COPY PUBLIC MEMBERS
//
//	// data_tuning, validation and test
//	// data tuning
//	n_tuning = p.n_tuning;
//	data_tuning = new Val*[n_tuning];
//	if (data_tuning==NULL)  {
//		cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of data_tuning failed!! Input data can't be imported" << endl;
//		exit(-1);
//	}
//	for (int i=0; i<n_tuning; i++) {
//		data_tuning[i] = new Val[n_var+1];
//		if (data_tuning[i]==NULL)  {
//			cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of data_tuning[" << i << "]  failed!! Input data can't be imported" << endl;
//			exit(-1);
//		}
//	}
//	for (int i=0; i<n_tuning; i++)
//		for (int j=0; j<n_var+1; j++)
//			data_tuning[i][j]= p.data_tuning[i][j];
//
//	// a little statistics on tuning data corresponding output
//	sum_output = p.sum_output;
//	y_ave = p.y_ave;
//	Sy = p.Sy;
//	y_var = p.y_var;
//
//	// data validation
//	n_validation = p.n_validation;
//	data_validation = new Val*[n_validation];
//	if (data_validation==NULL)  {
//		cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of data_validation failed!! Input data can't be imported" << endl;
//		exit(-1);
//	}
//	for (int i=0; i<n_validation; i++) {
//		data_validation[i] = new Val[n_var+1];
//		if (data_validation[i]==NULL)  {
//			cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of data_validation[" << i << "]  failed!! Input data can't be imported" << endl;
//			exit(-1);
//		}
//	}
//	for (int i=0; i<n_validation; i++)
//			for (int j=0; j<n_var+1; j++)
//				data_validation[i][j]= p.data_validation[i][j];
//
//	// data test
//	n_test = p.n_test;
//	data_test = new Val* [n_test];
//	if (data_test==NULL)  {
//		cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of data_test failed!! Test data can't be imported" << endl;
//		exit(-1);
//	}
//	for (int i=0; i<n_test; i++) {
//		data_test[i] = new Val[n_var+1];
//		if (data_test[i]==NULL)  {
//			cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of data_test[" << i << "]  failed!! Input data can't be imported" << endl;
//			exit(-1);
//		}
//	}
//	for (int i=0; i<n_test; i++)
//			for (int j=0; j<n_var+1; j++)
//				data_test[i][j]= p.data_test[i][j];
//
//	// a little statistics on testing data corresponding output
//	sum_output_test = p.sum_output_test;
//	y_ave_test = p.y_ave_test;
//	Sy_test = p.Sy_test;
//
//
//	// inequality constraints on values (order 0)
//	n_inequality0 = p.n_inequality0;
//	//char* constraints0;   no need to be copied
//	data_inequality0= new Val*[n_inequality0];
//	if (data_inequality0==NULL)  {
//		cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of data_inequality0 failed!! Input data can't be imported" << endl;
//		exit(-1);
//	}
//	for (int i=0; i<n_inequality0; i++) {
//		data_inequality0[i] = new Val[n_var+1];
//		if (data_inequality0[i]==NULL)  {
//			cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of data_inequality0[" << i << "]  failed!! Input data can't be imported" << endl;
//			exit(-1);
//		}
//	}
//	for (int i=0; i<n_inequality0; i++)
//		for (int j=0; j<n_var+1; j++)
//			data_inequality0[i][j]= p.data_inequality0[i][j];
//
//	// inequality constraints on first derivatives (order 1)
//	n_inequality1 = p.n_inequality1;
//	//char* constraints1;   no need to be copied
//	data_inequality1= new Val*[n_inequality1];
//	if (data_inequality1==NULL)  {
//		cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of data_inequality1 failed!! Input data can't be imported" << endl;
//		exit(-1);
//	}
//	for (int i=0; i<n_inequality1; i++) {
//		data_inequality1[i] = new Val[n_var+1];
//		if (data_inequality1[i]==NULL)  {
//			cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of data_inequality1[" << i << "]  failed!! Input data can't be imported" << endl;
//			exit(-1);
//		}
//	}
//	for (int i=0; i<n_inequality1; i++)
//			for (int j=0; j<n_var+1; j++)
//				data_inequality1[i][j]= p.data_inequality1[i][j];
//
//	// symbol
//	int len = (p.symbol).length();
//	char* expr = new char [len+1];
//	strcpy(expr,p.symbol.c_str());      //letter used for the variables
//	symbol.assign(expr);
//	delete [] expr;
//// so far so good...
//
//
//	v_list = new Variable *[n_var]; //array of pointers to Variable
//	if (v_list == NULL) {
//		cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of v_list failed. Exit." << endl;
//		exit(-1);
//	}
//
//	for (int k=0; k<n_var; k++)
//		v_list[k] = NULL;    //list of variables addresses
//
//
//
//	//Andrey's idea to avoid dynamic allocation of functions list
//	// in the future use a vector
//	for (int k=0; k<15; k++)
//			dummy_uni[k] = p.dummy_uni[k];
//	num_u_funcs = p.num_u_funcs;
//	for (int k=0; k<num_u_funcs; k++)
//		u_func_list[k] = (p.u_func_list)[k];
//
//	//Andrey's idea to avoid dynamic allocation of functions list
//	// in the future use a vector
//	for (int k=0; k<7; k++)
//		dummy_bin[k] = p.dummy_bin[k];
//	num_b_funcs = p.num_b_funcs;
//	for (int k=0; k<num_b_funcs; k++)
//			b_func_list[k] = (p.b_func_list)[k];
//
//	division = p.division;
//
//
//	variables_initialised = p.variables_initialised;
//
//	if (COMMENT) cout << "\nProblemDefinition copy constructor exit";
//	}
//}


// ProblemDefinition destructor
ProblemDefinition::~ProblemDefinition(void)
{
	// delete all the dynamically allocated variables!!!
	// still to understand why this destructor is called:
	// - after Population::get_tree_derivative_given_norm_vector

	cout << "ProblemDefinition::~ProblemDefinition(void) : destructor called" << endl;

	// delete folds_table array
	if (folds_table!=NULL) {
			for (int i = 0; i < n_data; ++i) delete[] folds_table[i];
			delete[] folds_table;
	}

	// delete points_per_fold array
	if (points_per_fold!=NULL)
		delete[] points_per_fold;

	// folds ??

	// array storing autocorrelation values
	if (n_var==1) {
		delete[] r_k;
	}

	cout << "ProblemDefinition::~ProblemDefinition(void) : exit" << endl;

}


//function to initialise varibles (see v_list)
void ProblemDefinition::initialise_variables(Variable**p_Z, double max_n_periods)
{
	int COMMENT = 1;

	// if variables already initialised skip this function
	if (variables_initialised) {
		cout << "\nProblemDefinition::initialise_variables : variables already initialised.";
		cout << "\nLeaving function.";
		return;
	}

	// the following is executed if variables have not been initialised previously

	// check if data has been defined
	if (data==NULL) {
		cerr << "\nProblemDefinition::initialise_variables() : ERROR! data not defined (data=NULL) ! Exit" << endl;
		exit(-1);
		return;
	}

	// check for errors in n_var
	if (n_var<1) {
		cerr << "\nProblemDefinition::initialise_variables() : ERROR! Cannot initialise variables as n_var = " << n_var << endl;
		exit(-1);
		return;
	}

	// DYNAMIC DECLARATION OF VARIABLE SET: declare and initialise nvar Variable and assign their names
	*p_Z = new Variable[n_var];  // array of Variables : its scope goes beyond the read_file function as it is allocated using new
	if (*p_Z==NULL) {
		cerr << "ProblemDefinition::initialise_variables() : ERROR ! Can't allocate *p_Z = new Variable[n_var]";
		exit(-1);
	}

	// in the future use ProblemDefinition constructor to create v_list...
	v_list = new Variable *[n_var]; //array of pointers to Variable
	if (v_list == NULL) {
		cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of v_list failed. Exit." << endl;
		exit(-1);
	}

	symbol = "Z";      //letter used for the variables - see problem definition


	// Nyquist constraint on omega for sin/cos terms
	// computed on independent variable values for n_var = 1: in the future create a "Nyquist" array so to store a different value for each variable depending on sampling frequency
	if (n_var==1) {
		cout << "ProblemDefinition::compute_inputdata_stats() - Ny_omega_max" << endl;
		Ny_omega_max = 0.5*3.141592653/((data[n_data-1][0]-data[0][0])/(n_data-1.0));
		cout << "ProblemDefinition::compute_inputdata_stats() - Ny_omega_max = " << Ny_omega_max;
	}

	for (int k=0; k<n_var; k++ ) {
		// NAME
		// establish first the link between the pointers' list and the Variable array
		v_list[k] = &((*p_Z)[k]);
		char num_field[10];
		sprintf(num_field, "%d", k+1);
		string num_field_s = num_field;
		string var_name;
		var_name = symbol + num_field_s;
		strcpy( v_list[k]->name, var_name.c_str());

		// RANGE and OMEGA_LIM
		Val max = data[0][k];
		Val min = data[0][k];
		for (int i=0; i<n_data; i++)	{
			if (data[i][k] > max)
				max = data[i][k];
			if (data[i][k] < min)
				min = data[i][k];
		}
		v_list[k]->lower_b = min;
		v_list[k]->upper_b = max;
		v_list[k]->range = max-min;
		v_list[k]->omega_lim = 0.5*3.141592653/((max-min)/(n_data-1.0)); // Nyquist value  //=2.0*PI*(max_n_periods)/(max-min);  // historical value of omega_lim (see PhD thesis)

		// print variable's status
		if (COMMENT) v_list[k]->show_status();
		//cout << "Have a go?" << endl;
		//cin.get();
	}


	// change variables_initialised state
	variables_initialised = 1;


}


// function to compute input data statistics (mainly statistical properties of target values set)
void ProblemDefinition::compute_inputdata_stats(void)
{
	// check that the data has been correctly imported
	if (n_var==-1) {  // number of variables
		cerr << "\nProblemDefinition::compute_inputdata_stats : ERROR : n_var=-1, stats cannot be computed!!! Data correctly imported?";
		exit (-1);
	}
	if (n_data==-1) {  //total number of rows (fitness cases) in data
		cerr << "\nProblemDefinition::compute_inputdata_stats : ERROR : n_data=-1, stats cannot be computed!!! Data correctly imported?";
		exit (-1);
	}
	if (data==NULL) {  //original data provided by the user (NOT to be touched)
		cerr << "\nProblemDefinition::compute_inputdata_stats : ERROR : data=NULL, stats cannot be computed!!! Data correctly imported?";
		exit (-1);
	}

	// compute statistical properties of DATA target values
	//--------------------------------------------------------

	// sum_output, y_ave, Y_min and y_max
	Val y_ave_temp = .0;
	Val a=data[0][n_var];
	Val a_max=data[0][n_var];
	int k_max=0;
	Val a_min=data[0][n_var];
	int k_min=0;
	sum_output=.0;
	for (int k=0; k < n_data; k++) {
		a=data[k][n_var];
		// sum of target squares
		sum_output = sum_output + a*a;
		// sum of the target values
		y_ave_temp = y_ave_temp + a;
		// min and max
		if (a>=a_max) {
			a_max=a;
			k_max=k;
		}
		if (a<=a_min) {
			a_min=a;
			k_min=k;
		}
		// total variation
		if (n_var==1) {
			if (k>0) tot_variation_input=tot_variation_input+fabs(data[k][n_var]-data[k-1][n_var]);
		}
	}
	// mean target value
	y_ave = y_ave_temp/(double)(n_data);
	// min target value
	y_min=a_min;
	index_max=k_max;
	// max target value
	y_max=a_max;
	index_min=k_min;

	// Sy=sum((output- average_output)^2)=(n-1)*output variance on the whole dataset
	// variance y_var
	Sy=.0;
	for (int k=0; k < n_data; k++) Sy = Sy + (data[k][n_var] - y_ave)*(data[k][n_var] - y_ave);
	y_var = Sy/(double)(n_data-1);

	// autocorrelation r_k of original signal (only for 1D case) on building data set (also r_k on test data set should be computed...)
	if (n_var==1) {
		delay_max=(int)floor(n_data/2); // 9/1/21 value to be discussed, this is just temporary. Also used in Population::fitness_func
		r_k= new Val[delay_max]; 	//array storing autocorrelation values
		for (int i=0; i<delay_max; i++) {
			r_k[i]=0.0;
		}
		cout << "delay_max = " << delay_max << endl;
		cout << "k	r_k" << endl;

		int first_acf_root_input_found=0;
		for (int delay=0; delay<delay_max; delay++) {
			//cout << "\ndelay=" << delay;
			for (int k=0; k < n_data-delay; k++) {
				r_k[delay] = r_k[delay] + (data[k][n_var]-y_ave)*(data[k+delay][n_var]-y_ave);
				//cout << k << "  c_k[" << delay << "] = " << c_k[delay] << endl;
			}
			r_k[delay]= r_k[delay]/Sy;
			//cout << "r_k[" << delay << "]= " << r_k[delay] << endl;

			// 20/2/21 search for first point at which ACF halves -------------- no longer first root of autocorrelation function (the one closest to 0)
			// if r_k[delay-1] r_k[delay]>0 are both positive or negative do nothing: r_k[delay-1]*r_k[delay]>0
			if (first_acf_root_input_found==0) {
				if (fabs(r_k[delay]-0.5)<1.0E-12) {
					first_acf_root_input=data[delay][0];  //mind! Only for n_var=1!
					first_acf_root_input_found=1;
				} else {
					if (delay>0) {
						if ( (r_k[delay-1]-0.5>0) && (r_k[delay]-0.5<0) ) {
							first_acf_root_input=(data[delay-1][0]+data[delay][0])/2.0; //mind! Only for n_var=1 and delay>0!
							first_acf_root_input_found=1;
						}
					}
				}
			}
		}
	} // end if n_var==1


	// compute statistical properties of TEST DATA target values
	//------------------------------------------------------------
	Val y_ave_test_temp = .0;
	a=data_test[0][n_var];
	a_max=data_test[0][n_var];
	a_min=data_test[0][n_var];
	sum_output_test = .0;
	for (int k=0; k < n_test; k++) {
		a=data_test[k][n_var];
		// sum of target squares
		sum_output_test = sum_output_test + a*a;
		// sum of the target values
		y_ave_test_temp = y_ave_test_temp + a;
		// min and max
		if (a>=a_max) a_max=a;
		if (a<=a_min) a_min=a;
	}
	// mean target value
	y_ave_test = y_ave_test_temp/(double)(n_test);
	// min target value
	y_test_min=a_min;
	// max target value
	y_test_max=a_max;


	// Sy=sum((output- average_output)^2)=(n-1)*output variance on the whole dataset
	// variance y_var
	Sy_test = .0;
	for (int k=0; k < n_test; k++) Sy_test = Sy_test + (data_test[k][n_var] - y_ave_test)*(data_test[k][n_var] - y_ave_test);
	y_var_test=Sy_test/(double)(n_test-1);

}


// show data
void ProblemDefinition::show_data(void)
{
	cout << "\n\ndata[][]" << endl;
	cout << setw(3) << " Z"; // << symbol;
	for (int j=0; j< n_cols - 2;  j++)
		cout << left << setw(22) << j+1 << "Z";
	cout << left << setw(22) << n_cols-1;
	cout << left << setw(22) << "TARGET";
	cout << endl;
	for (int i=0; i< n_data; i++)  {//m must be equal to ndata
		cout << left << setw(3) << i;
		for (int j=0; j< n_cols;  j++) {
			cout <<  scientific <<  left << setw(22)  << data[i][j];
		}
	cout << endl;
	}

	// stats
	cout << "\nsum output :  sum_output = " << sum_output;
	cout << "\naverage output : y_ave = " << y_ave;
	cout << "\n(output variance)*(n-1) : Sy = " << Sy << endl;

}

// show data_tuning
void ProblemDefinition::show_data_tuning(void) {
	cout << "\n\ndata_tuning[][]" << endl;
	cout << setw(3) << " Z"; // << symbol;
	for (int j=0; j< n_cols - 2;  j++)
		cout << left << setw(22) << j+1 << "Z";
	cout << left << setw(22) << n_cols-1;
	cout << left << setw(22) << "TARGET";
	cout << endl;
	for (int i=0; i< n_tuning; i++)  {//m must be equal to ntuning
		cout << left << setw(3) << i;
		for (int j=0; j< n_cols;  j++) {
			cout <<  scientific <<  left << setw(22)  << data_tuning[i][j];
		}
	cout << endl;
	}
}

// show data_validation
void ProblemDefinition::show_data_validation(void) {
	cout << "\n\ndata_validation[][]" << endl;
	cout << setw(3) << " Z"; // << symbol;
	for (int j=0; j< n_cols - 2;  j++)
		cout << left << setw(22) << j+1 << "Z";
	cout << left << setw(22) << n_cols-1;
	cout << left << setw(22) << "TARGET";
	cout << endl;
	for (int i=0; i< n_validation; i++)  {//m must be equal to n_validation
		cout << left << setw(3) << i;
		for (int j=0; j< n_cols;  j++) {
			cout <<  scientific <<  left << setw(22)  << data_validation[i][j];
		}
	cout << endl;
	}
}


// show data_test
void ProblemDefinition::show_data_test(void) {
	cout << "\n\ndata_test[][]" << endl;
	cout << setw(3) << " Z"; // << symbol;
	for (int j=0; j< n_cols - 2;  j++)
		cout << left << setw(22) << j+1 << "Z";
	cout << left << setw(22) << n_cols-1;
	cout << left << setw(22) << "TARGET";
	cout << endl;
	for (int i=0; i< n_test; i++)  {//m must be equal to n_validation
		cout << left << setw(3) << i;
		for (int j=0; j< n_cols;  j++) {
			cout <<  scientific <<  left << setw(22)  << data_test[i][j];
		}
	cout << endl;
	}

	// stats
	cout << "\nsum output :  sum_output_test = " << sum_output_test;
	cout << "\naverage output : y_ave_test = " << y_ave_test;
	cout << "\n(output variance)*(n-1) : Sy_test = " << Sy_test << endl;
}



//show data_inequality0
void ProblemDefinition::show_data_inequality0(void)
{
	// text field width
	int s = 18;

	cout << "\n\ndata_inequality0[][]" << endl;
	cout << "n_inequality0 = " << n_inequality0 << endl;
	cout << setw(3) << " " << symbol;
	for (int j=0; j< n_cols - 2;  j++)
		cout << left << setw(s-1) << j+1 << "Z";
	cout << left << setw(s-1) << n_cols-1;
	cout << left << setw(s-1) << "Constraint Value";
	cout << left << setw(s-1) << "Relationship";
	cout << endl;
	for (int i=0; i< n_inequality0; i++)  {//m must be equal to n_inequality0
		cout << left << setw(3) << i;
		for (int j=0; j< n_cols;  j++) {
			cout <<  scientific <<  left << setw(s)  << data_inequality0[i][j];
		}
		cout << constraints0[i];
		cout << endl;
	}
}


//show data_inequality1
void ProblemDefinition::show_data_inequality1(void)
{
	// text field width
	int s = 18;

	cout << "\n\ndata_inequality1[][]" << endl;
	cout << "n_inequality1 = " << n_inequality1 << endl;
	cout << setw(3) << " " << symbol;
	for (int j=0; j< n_var -1 ;  j++)
		cout << left << setw(s-1) << j+1 << symbol;
	cout << left << setw(s-1) << n_var;
	cout << "n";
	for (int j=0; j< n_var - 1;  j++)
			cout << left << setw(s-1) << j+1 << "n";
	cout << left << setw(s-1) << n_var;
	cout << left << setw(s-1) << "Partial der.";
	cout << left << setw(s-1) << "Relationship";
	cout << endl;
	for (int i=0; i< n_inequality1; i++)  {//m must be equal to n_inequality0
		cout << left << setw(3) << i;
		for (int j=0; j< 2*n_var+1;  j++) {
			cout <<  scientific <<  left << setw(s)  << data_inequality1[i][j];
		}

		cout << constraints1[i];
		cout << endl;
	}
}


// show binary functions
void ProblemDefinition::show_binary_functions(void)
{
	cout << "\n\nnum_b_funcs = " << num_b_funcs;
	cout << " Binary functions: ";
	int k;
	if (num_b_funcs>0) {
		for (k=0; k< num_b_funcs-1; k++) {
			cout << b_func_list[k]->sign << ", ";
		}
		cout << b_func_list[k]->sign;
	}
	else
		cout << "No binary functions used";
}


// show unary functions
void ProblemDefinition::show_unary_functions(void)
{
	cout << "\nnum_u_funcs = " << num_u_funcs;
	cout << " Unary functions: ";
	int k;
	if (num_u_funcs>0) {
		for (k=0; k < num_u_funcs-1; k++) {
			cout << u_func_list[k]->pre_sign << u_func_list[k]->post_sign << ", ";
		}
		cout << u_func_list[k]->pre_sign << u_func_list[k]->post_sign;
		cout << endl;
	}
	else
		cout << "No unary functions used" << endl;
}

// show lots of data
void ProblemDefinition::show_all(void)
{
	cout << "\n\nProblemDefinition::show_all(void)" << endl;
	show_binary_functions();
	show_unary_functions();
	show_data();    //original whole data set (tuning+evaluation) from input_file

	show_data_tuning();
	show_data_validation();
	show_data_test();

	show_data_inequality0();
	show_data_inequality1();
}


void ProblemDefinition::set_folds_table(void)
{
	// allocate folds_table
	folds_table = new int*[n_data];
	for (int i = 0; i < n_data; ++i)
		folds_table[i] = new int[2];

	// allocate points_per_fold
	if (n_folds>0) points_per_fold = new int[n_folds];

}



int ProblemDefinition::get_fold_from_row(int i)
{

	if (i>n_data-1) {
		cerr << "ProblemDefinition::get_fold_from_row(int) error : i>(n_data-1)";
		exit(EXIT_FAILURE);
	}

	// return the no. of the fold the data row i belongs to
	return folds_table[i][1];
}



int ProblemDefinition::get_points_per_fold(int i)
{
	if (i>(n_folds-1)) {
		cerr << "ProblemDefinition::get_points_per_fold(int) error : i>(n_folds-1)";
		exit(EXIT_FAILURE);
	}

	return points_per_fold[i];
}


// function to randomly split the whole input data in k data sets using the k-fold technique
// used to perform Cross validation using the PRESS predictor
// (http://scikit-learn.org/stable/modules/cross_validation.html)
// 26/5/2017 now it just defines an association table storing the fold number of each data row
void ProblemDefinition::kfold_split(int split_switch)
{
	cout << "\n\nProblemDefinition::kfold_split(int)" << endl;
	cout << "n_folds = " << n_folds << endl;

	if (split_switch==0) {
		cout << "\n\nCrossvalidation (K-folds method) not requested, p_folds_table allocated so to allow tuning on whole data set." << endl;

		for (int i=0; i<n_data; i++) {
			folds_table[i][0] = i;  	// column 0: data row number
			folds_table[i][1] = -1;  // column 1: fold number (number starts from 0)
		}
		set_validation_fold(-1);

		return;
	}

	// in the absence of explicit methods to generate/populate data_tuning for crossvalidation, set the pointer to NULL?


	// check that n_folds is not larger than n_data (if equal, the method turns into leave one out)
	if (n_folds > n_data) {
		// output error message
	    cerr << "ProblemDefinition::kfold_split : ERROR ! n_folds > n_data !!!!\n";
	    exit(-1);
	}


	// compute number of records per fold and at the same time initialise Val*** folds (array no, row no, col no)
	//folds = new Val**[n_folds]; USE VECTORS OR JUST POINTERS!
	int min_points = int((n_data-n_data%n_folds)/n_folds);
	for (int i=0; i<n_folds; i++) {
		points_per_fold[i] = min_points;
		if ((i+1)<=(n_data%n_folds)) points_per_fold[i]++;
		cout << "\nPoints_per_fold[" << i << "] = " << points_per_fold[i];
	}

	// Initialise folds_table (allocated in read_file_new line 568) to create an association table: data row number and corresponding fold number
	// scan the whole data set "data" and randomly assign each row to a different fold
	for (int i=0; i<n_data; i++) {
		folds_table[i][0] = i;  	// column 0: data row number
		folds_table[i][1] = i%n_folds;  // column 1: fold number (number starts from 0)
	}


	// show array
	cout << "\n\np_folds_table[]" << endl;
	//cout << setw(3) << " Z"; // << symbol;
	cout << left << setw(4) << "Row";
	cout << left << setw(10) << "data_row";
	cout << left << setw(10) << "fold_id_no";
	cout << left << setw(15) << "Pointed data" << endl;
	for (int i=0; i< n_data; i++)  {//m must be equal to ndata
		cout << left << setw(4) << i;
		cout <<  left << setw(10)  << folds_table[i][0] << left << setw(10) << folds_table[i][1];
		cout << left << setw(15) << data[folds_table[i][0]][0] << endl;
	}

	// print the folds to file for external use - or let reporter do it...

}








// show statistics on output
//void ProblemDefinition::show_output_statistics(void)
//{
//	cout << "\nsum output :  sum_output = " << sum_output;
//	cout << "\naverage output : y_ave = " << y_ave;
//	cout << "\n(output variance)*(n-1) : Sy = " << Sy << endl;
//}
