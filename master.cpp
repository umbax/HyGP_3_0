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


#include <iostream>  // basic i/o commands: cout, cin, scientific, fixed, cerr
#include <iomanip>  // manipulators (text format): setw
#include <string>    // to manipulate strings (C++)
#include <cstring>   // to manipulate strings (C) (strcpy, strcmp)
#include <cstdlib>   // NULL, exit, EXIT_FAILURE
#include <cmath>	// pow, sqrt
#include <algorithm> // max_element, min_element
#include <exception>

/// compatibility issues (Microsoft tempting to modify commands?):commented if under linux, uncommented under WINDOWS
//#define strdup _strdup

using namespace std;

// headers
#include "./genetic_code/classes/run_parameters.h"
#include "./genetic_code/classes/problem_definition.h"
#include "./genetic_code/read_input/read_file_new.h" // read test data not implemented yet in OpenMP version
#include "./genetic_code/input_checks/input_check.h"
#include "./genetic_code/classes/reporter.h"
#include "./genetic_code/classes/class_POPULATION.h"



// global variables - population purposes
//const double PI = 4.0*atan(1.0); // redefined in Val_type.h
//double MAX_VAL = 1.8e+19;   // redefined in Val_type.h
//double MIN_VAL = 1.8e-19;		// redefined in Val_type.h
// try! Simple way to make fdf_c able to detect the values of a few variables. Think about inserting fdf_c in class Population!!!
Population *Pop; // try! Simple way to make fdf_c able to get the values of the tree
int VERBOSE = 1;     //set to 1 if you want to print on the screen all the comments! 0 for a "clean" and wordless execution...

// use: ./gp location/name_input location/name_output

int main (int argc, char *argv[])
{
	// check the number of arguments
	if (argc<3) {
		cerr << "\nERROR!!! Too few arguments!!!"	;
		cerr << "\nUSAGE: >> ./gp  location/input_file  existing_directory_output";
		cerr << "\nExample: >> ./gp ./input/input_file.txt ./output\n";
		exit(-1);
	}
	if (argc>4) {
		cerr << "\nERROR!!! Too many arguments!!!"	;
		cerr << "\nUSAGE: >> ./gp  location/input_file  location/test_data  existing_directory_output";
		cerr << "\nExample: >> ./gp ./input/input_file.txt ./input/test_data_file.txt ./output\n";
		exit(-1);
	}
	
	string FILE_INPUT;
	string FILE_TEST_DATA;
	string DIR_OUTPUT;
	FILE_TEST_DATA = "Not defined";
	FILE_INPUT=argv[1];
	if (argc == 3) {
		DIR_OUTPUT=argv[2];
	} else {
		FILE_TEST_DATA=argv[2];
		DIR_OUTPUT=argv[3];
	}


 	cout << "\nNAME_INPUT = " << FILE_INPUT;
 	cout << "\nNAME_TEST_DATA = " << FILE_TEST_DATA;
 	cout << "\nNAME_OUTPUT = " << DIR_OUTPUT;
	


	//presentation
	printf("\n\n			****************************			");
	printf("\n			 FIRST EXPERIMENTS IN GP   			");
	printf("\n														           ");
	printf("\n			 by Umberto Armani		            ");
	printf("\n			****************************			");
	cout << "\n\nAvailable operations:";
	cout << "\nBINARY: add , sub, mult, sdiv, spow";
	cout << "\nUNARY: shift, neg, square, cube, exp, nxp, sin, cos, inv, abs, log, sinh, cosh, tanh, rectwave, hfreqsin";
	cout << "\n(list updated 7/11/2020)";
	cout << "\n\nATTENTION! If you use more than 15 unary operations or 7 binary operations";
	cout << "\nincrease the size of u_func_list or b_func_list (see problem_definition.h)";
	cout << "\n(FUTURE solution : use vectors to read primitives...)";
	cout << endl;

	//set the number of decimal figures.
	cout.precision(5);  

	// other classes instantiations
	Variable* Z;    //to be included in ProblemDefinition...
	RunParameters Mparam;  		// to distinguish it from parameters object in Population
	ProblemDefinition Mprobl;	// to distinguish it from problem object in Population
	Reporter pop_reporter;
	
	// read inputs
	read_input_file(FILE_INPUT, &Mparam, &Mprobl);  // also initialise parameters and problem objects
	if (argc==4) read_test_data(FILE_TEST_DATA, &Mparam, &Mprobl);

	// check for errors on the input parameter between parenthesis
	input_check(&Mparam, &Mprobl);
	
	// compute additional attributes or set up structures in Problem Definition (variables, for example)
	// now all is done in read_input_file (but it's too messy there...)
	Mprobl.compute_inputdata_stats();
	
	// print to file input data statistics
	pop_reporter.inputdatastats2file(&Mprobl, DIR_OUTPUT);

	// show imported data (input parameters and input data)
	Mparam.show();  // run hyperparameters
	//Mprobl.show_all();    // input data matrix and other data

		
	// to stop the execution
	//cout << "Have a go?" << endl;
	//cin.get();
	
	// random value generator
	if (Mparam.seed<0) {
        // if seed = -1 use random seed (time)
        time_t *tp = NULL; // seed the random value generator
	    srand((unsigned int) time (tp));    //to be used normally
	    Mparam.seed = time(tp); // attention! two runs that starts within 1 second have the same seed, so are identical!
        cout << "\n\nSEED=-1 in input file: seed randomly generated = " << Mparam.seed << endl;
    } else {
        // if seed >0 use the value given in input file
        srand(Mparam.seed);
        cout << "\n\nused seed = " << Mparam.seed << endl;
    }
	

	// PARALLELISATION. FROM HERE ....

	time_t start, finish;
	double delta_t;
    time(&start);	
    int elapsed_time;
    
	//---------------------------------------------------------------------
	// CREATE THE POPULATION (constructor called)
	//---------------------------------------------------------------------
    // initialise variables
    Mprobl.initialise_variables(&Z, Mparam.max_n_periods);
    //for (int k=0; k<Mprobl.get_n_var(); k++)
    //	cout << "\n v_list" << k << " = " << Mprobl.v_list[k];

	// do you really need new? no, but don't touch it for now
	Population *P = new Population(&Mparam, &Mprobl);
	if (!P) {
		cerr << "\nmain : Error creating population!!!\n";
		exit(-1); 
	}
	
	// pause execution
	//cin.get(); // problem with int_rand in Population constructor

	//as it's hard to pass Pop as a parameter to fdf_c__ through fortran functions, treat it as a global variable
	Pop = P;


	int n =Pop->parameters->nvar;
	cout << "\nn= " << n;

	// split the whole input dataset in k folds for cross validation
	Mprobl.kfold_split(Mparam.crossvalidation);  // test with 3 folds



	///// INITIAL GENERATION (0) ///////////////////////////////////////
	
	// trees before parameters insertion		
	printf ("\nInitialization of the population (generation 0)\n");
	P->print_population_without_parameters(0);
	

	/// split the data set in tuning set (data_tuning) and validation set (data_validation).
	// See SPLIT and VALIDATING_LINES in input file
	// this function will also allow to increase the number of fitness cases during the run...
	// 23/5/2017 rewrite the split function to implement correctly the CROSSVALIDATION and the PRESS error calculation
	// IMPORTANT! Check that the correct Sy is used in computing RMSE (note taken on 10/5/2014)
	//P->split_data(&Mparam, &Mprobl, 0,1);
	

	/////////// OTHER GENERATIONS ///////////////////////////////////////
	int check_end = 0;
	int last_gen=0;
	for (int i=0; i<Mparam.G+1; i++) {
		
		if (i) {   //skip generation 0
			// split the data for the current generation (this function will also allow to increase the number of fitness cases during the run...)
			//P->split_data(i,G,split); // not used...

			// GENETIC OPERATORS: sorting, reproduction, crossover, mutation, pruning
			P->new_spawn(Mparam, Mprobl, Mparam.nfitcases,i);

		}

		// evaluate fitness function (in hybrid/memetic GP parameters are added and tuned first, then the evaluation is performed)
		P->evaluate(i,Mparam.G);

		// sort according to F (aggregate fitness, not RMSE!) : VITAL! Both populations must be sorted, trees[] and complete_trees[]
		P->sort(i,tree_comp_F);

		// update the best individual - structure and complete tree (for PARAMETER INHERITANCE)
		P->update_ext_archive();

		// compute elapsed time
		elapsed_time = (int)(P->compute_time(start, finish, &delta_t));

		if (VERBOSE) {
			// print elapsed time
			cout << "\nElapsed time: " << elapsed_time << " sec";  //total seconds

			// print out the best member - population WITHOUT parameters
			P->print_population_without_parameters(i);
	
			// print out the best member - population WITH parameters
			P->print_population_with_parameters(i);
		}

		// compute statistical data relating to population (vital if data is shared through populations)
		// IT's REALLY IMPORTANT that this function is executed after evaluate and sort, 
		// as evaluate uses statistical data referring to the previous generation
		P->compute_statistics();

		// evaluate termination condition
		check_end=P->terminate(Mparam.threshold);
		last_gen = i;   


		// PRINT TO FILE OPERATIONS (in case of crash, data is saved)  -------------------------
		cout << "\nPrint generation results to file...";
		// write evolution statistical data to file "data_gp.txt"
		pop_reporter.stats2file(&Mparam, P, DIR_OUTPUT, i, check_end);
		// write the target points (training set), best individual corresponding values and residuals to file "points_gp.txt"
		pop_reporter.points2file(&Mparam, &Mprobl, P, DIR_OUTPUT, i, check_end, start, finish, delta_t, Mparam.seed);
		// write best individual's expression as per aggregate value F (!!!) on training data set to "best_gp.txt"
		pop_reporter.update_best2file_build(P, DIR_OUTPUT, i, check_end);
		// write list of the best-so-far individuals (see elite or archive - first repr_tot individuals) on training data set to "latest_archive.txt"
		pop_reporter.archive2file_build(P, DIR_OUTPUT, i, check_end);
		// write no of tree evaluations at each generation to "n_tree_evaluations.txt"
		pop_reporter.n_tree_eval2file(P, DIR_OUTPUT, i, check_end);
		// write related to adaptive approach (eps_neutral, the constructive, destructive and neutral genetic operations rates, etc) to "adaptation_data.txt"
		pop_reporter.adaptive_gen_ops_data2file(P, DIR_OUTPUT, i, check_end);
		cout << "OK";
		// -------------------------------------------------------------------------------------
		
		if (check_end)
			break;
	
		// update genetic operators rates (adaptive approach - can be turned on and off inside the function)
		if (i)
			P->adapt_genetic_operators_rates();
  
	}
	// end evolution


	//termination criterion (successful run) satisfied
	if (check_end) {
		cout << "Termination criterion satisfied (RMSE < " << Mparam.threshold << ")." << endl;
		cout << "Possible solution: " << endl; 
		P->print_population_with_parameters(last_gen);
		cout << "Check latest_archive.txt for solutions\n" << endl;

	}
	else {
		P->print_population_with_parameters(last_gen);
	}

	// just for test
	//P->get_tree_derivative_given_norm_vector(Mprobl, P->complete_trees[0]);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	//---------------------------------------------------------------------
	// END OPERATIONS
	//---------------------------------------------------------------------
	// write node statistics to file
	pop_reporter.node_stats2file(&Mparam, P, DIR_OUTPUT);
	
	// evaluate fitness (RMSE and R2) on test data set (only if test data has been provided)
	if (argc==4) {
		// show data_test
		cout << "problem->show_data_validation() : show current data_test :" << endl;
		P->problem->show_data_test();
		// evaluate complete individuals on test data set provided by the user
		P->evaluate_complete_trees(); // SET CORRECTLY Mprobl.data_test, n_test, Sy_test after implementing function to read test data set
		// sort according to error (RMSE) - better not to use it to keep order and to recognise performance on building and test data sets...
		//P->sort(last_gen,tree_comp_fitness); // non va: ordina in ordine decrescente e alcune volte pone a 0 RMSE e R2. Perché?
		// find and print best individual on test data set as per RMSE (!!!!) to file "best_gp_TEST.txt"
		pop_reporter.best2file_test(P, DIR_OUTPUT, last_gen);
		// print archive evaluated on the test data set to file
		pop_reporter.archive2file_test(P, DIR_OUTPUT, last_gen);  // insert a function to order in rmse decreasing order, leaving however the name of the run
	}

	// free memory allocated to Population, as not used anymore (in the future declare statically P, as you will always need one population...)
	delete P;

	// PARALLELISATION ... TO HERE.

	cout << "\n\n END" << endl;

	return 0;
}


// DIFFERENT OPTIMISERS
//#include "./genetic code/SQP/MINL2.cpp"    - only for optimizer translated in C++
//#include "./genetic code/SQP/TINL2_mod.cpp"   - only for optimizer translated in C++

//for Andrey's method - TI0L2 and MI0L2 - IT WORKS PERFECTLY
#include "./genetic_code/SQP/MI0L2_c/fdf_c.cpp"
//for Umberto's method - TINL2 and MINL2 - WORKS?
//#include "./genetic code/SQP/fdfTINL2_c.cpp"  


