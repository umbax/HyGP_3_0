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


#ifndef PROBLEMDEFINITION_H_
#define PROBLEMDEFINITION_H_

// dependencies
#include <iostream>  // basic i/o commands: cout, cin, scientific, fixed, cerr
#include <cstdlib>   // NULL, exit, EXIT_FAILURE
#include <iomanip>  // manipulators (text format): setw
#include <cstdlib>   // NULL, exit, EXIT_FAILURE
#include <cstdio>   // sprintf
#include <string>    // to manipulate strings (C)
#include <cstring>   // to manipulate strings (C) (strcpy,)
#include <vector>
#include <cmath>	// floor


#include "../modules/Val_type.h"
#include "../modules/variable.h"
#include "../modules/func_primitives_prototypes.h"

// class containing the data used during the evolution

class ProblemDefinition
{
	private:
		// original imported data set
		Val** data;			// original data provided by the user in the main input file (NOT to be touched)
		int n_data;			// total number of rows in data (nfitcases in RunParameters)
		int n_var;			// number of variables
		int n_cols;			// number of columns in data (n_var+1)
		// cross validation (PRESS evaluation) purposes
		int n_folds;		// number of folds for crossvalidation
		int** folds_table; // association table between fold and data row (n_datax2), used in model tuning
		int* points_per_fold; // array of int containing the number of points (records) per fold (n_foldsx1) - so the n. of validation points
		int validation_fold; // current fold (int) used for validation
		
	public:

		// constructor
		ProblemDefinition(void);
		// copy constructor
		ProblemDefinition(const ProblemDefinition&);
		// destructor
		~ProblemDefinition(void);

		// a little statistics of corresponding output of the WHOLE data set (DATA) contained in the main input file
		// (the one containing HyGP hyperparameters)
		Val sum_output; // sum of the squares of each target value
		Val y_ave; 		// average value of input data
		Val Sy;  		// SStot total sum of squares of (observed data - average observed data) // defined in read_input_file function (read_file_new.cpp)
		Val y_var; 		// variance of input data
		Val y_max;		// max value of target
		Val y_min;		// min value of target

		// folds for crossvalidation
		//Val*** folds; // pointer to a 3d array (no of fold x row x column) - initialised in kfold_split

		// data tuning and evaluation (tuning = building data, evaluation = cross validation data)
		Val** data_tuning;		// BUILDING/TUNING DATA SET : not defined explicitly
		int  n_tuning;
		Val** data_validation;	// CROSS VALIDATION DATA SET (coincides with building data if split = 0)
		int n_validation;	 //ex n_evaluation;

		// TEST DATA SET
		Val** data_test;
		int n_test;
		Val	sum_output_test;
		Val	y_ave_test;
		Val Sy_test;
		Val y_var_test;
		Val y_test_min;
		Val y_test_max;

		// inequality constraints on values (order 0)
		Val** data_inequality0;
		int n_inequality0;
		char* constraints0;
		// inequality constraints on first derivatives (order 1)
		Val** data_inequality1;
		int n_inequality1;
		char* constraints1;

		string symbol;      //letter used for the variables
		Variable** v_list;    //list of variables addresses
		//Variable *Z;  //experiment to bring Variable* Z in ProblemDefinition
		
		//Andrey's idea to avoid dynamic allocation of functions list
		// in the future use a vector
		Unary_Func dummy_uni[15];  
		Unary_Func *u_func_list[15];   // max of binary functions allowed is 15!
		int num_u_funcs;
		
		//Andrey's idea to avoid dynamic allocation of functions list
		// in the future use a vector
		Binary_Func dummy_bin[7];  
		Binary_Func *b_func_list[7];   // max of binary functions allowed is 7!
		Binary_Func *division;
		int num_b_funcs;

		//setter methods for private members
		void set_data(Val **ext_data) {data = ext_data;};
		void set_n_data(int tot_rows) {n_data = tot_rows;};
		void set_n_var(int n_variables) {n_var = n_variables;};
		void set_n_cols(int tot_cols) {n_cols = tot_cols;};
		void set_n_folds(int n_f) {n_folds=n_f;};
		void set_folds_table(void); // used in read_file_new.cpp line 568
		void set_validation_fold(int i) {validation_fold=i;};  //add feasibility checks!!!!
		// getter methods for private members
		Val get_data(int row, int column) {return data[row][column];};
		Val** get_data_address(void) {return data;};
		int get_n_data(void) {return n_data;};
		int get_n_var(void) {return n_var;};
		int get_n_cols(void) {return n_cols;};
		int get_n_folds(void) {return n_folds;};
		int get_points_per_fold(int i);
		int get_fold_from_row(int i);
		int get_validation_fold(void) {return validation_fold;};
		
		//function to initialise variables (see v_list)
		int variables_initialised;
		void initialise_variables(Variable** , double);

		// function to compute input data statistics
		void compute_inputdata_stats(void);

		// function to display data on the screen
		void show_data(void);
		void show_data_tuning(void);
		void show_data_validation(void);
		void show_data_test(void);

		void show_data_inequality0(void);
		void show_data_inequality1(void);
		void show_binary_functions(void);
		void show_unary_functions(void);
		void show_all(void);
		//void show_output_statistics(void);

		// function to split the whole input data in k data sets using the k-fold technique
		// used to perform Cross validation using the PRESS predictor
		// (http://scikit-learn.org/stable/modules/cross_validation.html)
		void kfold_split(int);
};

#endif /* PROBLEMDEFINITION_H_ */
