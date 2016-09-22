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

// ProblemDefinition constructor
ProblemDefinition::ProblemDefinition(void)
{
	// private members
	data = NULL;			//original data provided by the user (NOT to be touched)
	n_data = -1;			//total number of rows in data
	n_var = -1;			// number of variables
	n_cols = -1;			//number of columns in data (n_var+1)

	// public members
	data_tuning = NULL;
	n_tuning = -1;
	data_evaluation = NULL;
	n_evaluation = -1;
	// a little statistics on output (refers to the whole data set though - tuning+evaluation)
	sum_output = -1.0;
	y_ave = -1.0;
	Sy = -1.0;

	data_test = NULL;
	n_test = -1;
	sum_output_test = -1.0;
	y_ave_test = -1.0;
	Sy_test = -1;



	// inequality constraints on values (order 0)
	data_inequality0 = NULL;
	n_inequality0 = -1;
	constraints0 = ""; //previous version
	//strcpy(constraints0, "" );
	// inequality constraints on first derivatives (order 1)
	data_inequality1 = NULL;
	n_inequality1 = -1;
	constraints1 = ""; //previous version
	//strcpy(constraints1, "" );


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

}

// ProblemDefinition copy constructor (used for parallelisation - see firstprivate)
ProblemDefinition::ProblemDefinition(const ProblemDefinition& p)
{
	int COMMENT = 0;

# pragma omp critical
	{
	// here the instructions to make a deep copy
	if (COMMENT) cout << "\n\ncopy constructor entered" << endl;

	// copy private members
	n_data = p.n_data;			//total number of rows in data
	n_var = p.n_var;			// number of variables
	n_cols = p.n_cols;			//number of columns in data (n_var+1)
	data = new Val*[n_data];
	if (data==NULL)  {
		cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of data failed!! Input data can't be imported" << endl;
		exit(-1);
	}
	for (int i=0; i<n_data; i++) {
		data[i] = new Val[n_var+1];
		if (data[i]==NULL)  {
			cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of data[" << i << "]  failed!! Input data can't be imported" << endl;
			exit(-1);
		}
	}
	for (int i=0; i<n_data; i++)
		for (int j=0; j<n_var+1; j++)
			data[i][j]= p.data[i][j];



	// copy public members

	// data_tuning and evaluation
	// data tuning
	n_tuning = p.n_tuning;
	data_tuning = new Val*[n_tuning];
	if (data_tuning==NULL)  {
		cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of data_tuning failed!! Input data can't be imported" << endl;
		exit(-1);
	}
	for (int i=0; i<n_tuning; i++) {
		data_tuning[i] = new Val[n_var+1];
		if (data_tuning[i]==NULL)  {
			cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of data_tuning[" << i << "]  failed!! Input data can't be imported" << endl;
			exit(-1);
		}
	}
	for (int i=0; i<n_tuning; i++)
		for (int j=0; j<n_var+1; j++)
			data_tuning[i][j]= p.data_tuning[i][j];



	// data evaluation
	n_evaluation = p.n_evaluation;
	data_evaluation = new Val*[n_evaluation];
	if (data_evaluation==NULL)  {
		cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of data_evaluation failed!! Input data can't be imported" << endl;
		exit(-1);
	}
	for (int i=0; i<n_evaluation; i++) {
		data_evaluation[i] = new Val[n_var+1];
		if (data_evaluation[i]==NULL)  {
			cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of data_evaluation[" << i << "]  failed!! Input data can't be imported" << endl;
			exit(-1);
		}
	}
	for (int i=0; i<n_evaluation; i++)
			for (int j=0; j<n_var+1; j++)
				data_evaluation[i][j]= p.data_evaluation[i][j];


	// inequality constraints on values (order 0)
	n_inequality0 = p.n_inequality0;
	//char* constraints0;   no need to be copied
	data_inequality0= new Val*[n_inequality0];
	if (data_inequality0==NULL)  {
		cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of data_inequality0 failed!! Input data can't be imported" << endl;
		exit(-1);
	}
	for (int i=0; i<n_inequality0; i++) {
		data_inequality0[i] = new Val[n_var+1];
		if (data_inequality0[i]==NULL)  {
			cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of data_inequality0[" << i << "]  failed!! Input data can't be imported" << endl;
			exit(-1);
		}
	}
	for (int i=0; i<n_inequality0; i++)
		for (int j=0; j<n_var+1; j++)
			data_inequality0[i][j]= p.data_inequality0[i][j];

	// inequality constraints on first derivatives (order 1)
	n_inequality1 = p.n_inequality1;
	//char* constraints1;   no need to be copied
	data_inequality1= new Val*[n_inequality1];
	if (data_inequality1==NULL)  {
		cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of data_inequality1 failed!! Input data can't be imported" << endl;
		exit(-1);
	}
	for (int i=0; i<n_inequality1; i++) {
		data_inequality1[i] = new Val[n_var+1];
		if (data_inequality1[i]==NULL)  {
			cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of data_inequality1[" << i << "]  failed!! Input data can't be imported" << endl;
			exit(-1);
		}
	}
	for (int i=0; i<n_inequality1; i++)
			for (int j=0; j<n_var+1; j++)
				data_inequality1[i][j]= p.data_inequality1[i][j];

	// symbol
	int len = (p.symbol).length();
	char* expr = new char [len+1];
	strcpy(expr,p.symbol.c_str());      //letter used for the variables
	symbol.assign(expr);
	delete [] expr;
// so far so good...


	v_list = new Variable *[n_var]; //array of pointers to Variable
	if (v_list == NULL) {
		cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of v_list failed. Exit." << endl;
		exit(-1);
	}

	for (int k=0; k<n_var; k++)
		v_list[k] = NULL;    //list of variables addresses



	//Andrey's idea to avoid dynamic allocation of functions list
	// in the future use a vector
	for (int k=0; k<15; k++)
			dummy_uni[k] = p.dummy_uni[k];
	num_u_funcs = p.num_u_funcs;
	for (int k=0; k<num_u_funcs; k++)
		u_func_list[k] = (p.u_func_list)[k];

	//Andrey's idea to avoid dynamic allocation of functions list
	// in the future use a vector
	for (int k=0; k<7; k++)
		dummy_bin[k] = p.dummy_bin[k];
	num_b_funcs = p.num_b_funcs;
	for (int k=0; k<num_b_funcs; k++)
			b_func_list[k] = (p.b_func_list)[k];

	division = p.division;


	// a little statistics on output
	sum_output = p.sum_output;
	y_ave = p.y_ave;
	Sy = p.Sy;

	variables_initialised = p.variables_initialised;

	if (COMMENT) cout << "\nProblemDefinition copy constructor exit";
	}
}


// ProblemDefinition destructor
ProblemDefinition::~ProblemDefinition(void)
{


}


//function to initialise varibles (see v_list)
void ProblemDefinition::initialise_variables(Variable**p_Z, double max_n_periods)
{
	int COMMENT = 1;

	// if variables already initialised exit...
	if (variables_initialised) {
		cout << "\nProblemDefinition::initialise_variables : variables already initialised.";
		cout << "\nLeaving function.";
		return;
	}


	// here variables has not been initialised previously

	// DYNAMIC DECLARATION OF VARIABLE SET: declare and initialise nvar Variable and assign their names
	//Variable *
	*p_Z = new Variable[n_var];  // array of Variables : its scope goes beyond the read_file function as it is allocated using new


	// in the future use ProblemDefinition constructor to create v_list...
	v_list = new Variable *[n_var]; //array of pointers to Variable
	if (v_list == NULL) {
		cerr << "\nERROR: ProblemDefinition copy constructor : dynamic allocation of v_list failed. Exit." << endl;
		exit(-1);
	}

	symbol = "Z";      //letter used for the variables - see problem definition

	for (int k=0; k<n_var; k++ ) {
		// NAME
		// establish first the link between the pointers' list and the Variable array
		v_list[k] = &((*p_Z)[k]);
		char num_field[10];
		sprintf(num_field, "%d", k+1);
		string num_field_s = num_field;
		string var_name;
		var_name = symbol + num_field_s;
		strcpy( (*(v_list[k])).name, var_name.c_str());

		// RANGE and OMEGA_LIM
		Val max = data[0][k];
		Val min = data[0][k];
		for (int i=0; i<n_data; i++)	{
			if (data[i][k] > max)
				max = data[i][k];
			if (data[i][k] < min)
				min = data[i][k];
		}
		(*(v_list[k])).lower_b = min;
		(*(v_list[k])).upper_b = max;
		(*(v_list[k])).range = max-min;
		(*(v_list[k])).omega_lim = 2.0*PI*(max_n_periods)/(max-min);

		// print variable's status
		if (COMMENT) (*(v_list[k])).show_status();
		//cout << "Have a go?" << endl;
		//cin.get();
	}

	// change variables_initialised state
	variables_initialised = 1;


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

// show data_evaluation
void ProblemDefinition::show_data_evaluation(void) {
	cout << "\n\ndata_evaluation[][]" << endl;
	cout << setw(3) << " Z"; // << symbol;
	for (int j=0; j< n_cols - 2;  j++)
		cout << left << setw(22) << j+1 << "Z";
	cout << left << setw(22) << n_cols-1;
	cout << left << setw(22) << "TARGET";
	cout << endl;
	for (int i=0; i< n_evaluation; i++)  {//m must be equal to n_evaluation
		cout << left << setw(3) << i;
		for (int j=0; j< n_cols;  j++) {
			cout <<  scientific <<  left << setw(22)  << data_evaluation[i][j];
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
	for (int i=0; i< n_test; i++)  {//m must be equal to n_evaluation
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
			cout << (*(b_func_list[k])).sign << ", ";
		}
		cout << (*(b_func_list[k])).sign;
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
			cout << (*u_func_list[k]).sign << ", ";
		}
		cout << (*u_func_list[k]).sign;
		cout << endl;
	}
	else
		cout << "No unary functions used" << endl;
}


// show statistics on output
//void ProblemDefinition::show_output_statistics(void)
//{
//	cout << "\nsum output :  sum_output = " << sum_output;
//	cout << "\naverage output : y_ave = " << y_ave;
//	cout << "\n(output variance)*(n-1) : Sy = " << Sy << endl;
//}