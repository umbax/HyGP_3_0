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

#include <iostream>  // basic i/o commands
#include <fstream>  // file I/O
#include <iomanip>  // manipulators (text format)
#include <string>    // to manipulate strings (C)
#include <cstring>   // to manipulate strings (C)
#include <string>    //  string class (C++)
#include <cstdlib>   // NULL
#include <ctime>	 // to work with variables of type time_t and random number generator "srandom"  
#include <cmath>
#include <algorithm>
#include <vector>
#include <exception>
#include <sstream>    // for stringstream

///// compatibility issues (Microsoft tempting to modify commands?)
// commented if under linux, uncommented under WINDOWS
//#define strdup _strdup

using namespace std;

// global variables - numeric settings
const double PI = 4.0*atan(1.0);
double MAX_VAL = 1.8e+19;
double MIN_VAL = 1.8e-19;

// HEADERS
// do not change the order of the following include commands! 
//Otherwise you may have undeclared objects //
#include "./genetic_code/modules/nodes_types.h"   	//just a header, no source file
#include "./genetic_code/modules/Val_type.h"		//just a header, no source file
#include "./genetic_code/modules/variable.h"
#include "./genetic_code/modules/func_primitives_prototypes.h"
#include "./genetic_code/classes/run_parameters.h"
#include "./genetic_code/classes/problem_definition.h"
#include "./genetic_code/read_input/show_loaded_data.h"
#include "./genetic_code/read_input/read_file_new.h"
#include "./genetic_code/read_input/read_test_data.h"
#include "./genetic_code/input_checks/input_check.h"
#include "./genetic_code/classes/class_NODE_base.h"
#include "./genetic_code/tree_functions/tree_operations.h"
#include "./genetic_code/tree_functions/vector_derivative_functions.h"
#include "./genetic_code/modules/primitives.cpp"
#include "./genetic_code/classes/class_BINARY_NODE.h"
#include "./genetic_code/classes/class_UNARY_NODE.h"
#include "./genetic_code/classes/class_TERMINAL_VAR.h"
#include "./genetic_code/classes/class_TERMINAL_CONST.h"
#include "./genetic_code/classes/class_POPULATION.h"
#include "./genetic_code/classes/reporter.h"

extern "C" {extern void opti_(int*, double*, int*, int*, int*);}  //for Andrey's method - TI0L2 and MI0L2 - IT WORKS PERFECTLY
extern "C" {int  fdf_c__ (double *,  double *, int *,int *); }

//extern "C" {extern void optitinl2_(int*, double*, int*, int*,int*);}  //for Umberto's method - TINL2 and MINL2 - WORKS? WORK STILL NEEDED... 
//extern "C" {int  fdf_ctinl2__ (double*,  double*, double*, int*, int*); }


// global variables - population purposes
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
		cerr << "\nUSAGE: >> ./gp  location/input_file  location/test_data existing_directory_output";
		cerr << "\nExample: >> ./gp ./input/input_file.txt ./input/test_data_file.txt ./output\n";
		exit(-1);
	}
	
	string FILE_INPUT, FILE_TEST_DATA, DIR_OUTPUT;
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
	cout << "\nUNARY: shift, neg, square, cube, exp, nxp, sin, cos, inv, abs, log, sinh, cosh, tanh";
	cout << "\n(list updated 20/10/2010)";
	cout << "\n\nATTENTION! If you use more than 15 unary operations or 7 binary operations";
	cout << "\nincrease the size of u_func_list or b_func_list (see problem_definition.h)";
	cout << "\n(FUTURE solution : use vectors to read primitives...)";
	cout << endl;

	//set the number of decimal figures.
	cout.precision(5);  

	// other classes instantiations
	Variable* Z;    //to be included in ProblemDefinition...
	RunParameters parameters;
	ProblemDefinition problem;
	Reporter pop_reporter;

	//int n_point_evals=0;   ????
	
	// read inputs
	read_input_file(FILE_INPUT, &parameters, &problem);
	if (argc==4) read_test_data(FILE_TEST_DATA, &parameters, &problem);

	// check for errors on the input parameter between parenthesis
	input_check(&parameters, &problem);
	
	// compute additional attributes or set up structures in Problem Definition (variables, for example)
	// now all is done in read_input_file (but it's too messy there...)
	
	// show the results
	//show_loaded_data(&parameters, &problem);
		
	// to stop the execution
	//cout << "Have a go?" << endl;
	//cin.get();
	
	// random value generator
	if (parameters.seed<0) {
        // if seed = -1 use random seed (time)
        time_t *tp = NULL; // seed the random value generator
	    srand((unsigned int) time (tp));    //to be used normally
	    parameters.seed = time(tp); // attention! two runs that starts within 1 second have the same seed, so are identical!
        cout << "\n\nSEED=-1 in input file: seed randomly generated = " << parameters.seed << endl;   	
    } else {
        // if seed >0 use the value given in input file
        srand(parameters.seed); 
        cout << "\n\nused seed = " << parameters.seed << endl;   	   
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
    problem.initialise_variables(&Z, parameters.max_n_periods);
    //for (int k=0; k<problem.get_n_var(); k++)
    //	cout << "\n v_list" << k << " = " << problem.v_list[k];

	// do you really need new? no, but don't touch it for now
	Population *P = new Population(&parameters, &problem);
	if (!P) {
		cerr << "\nmain : Error creating population!!!\n";
		exit(-1); 
	}
	
	//as it's hard to pass Pop as a parameter to fdf_c__ through fortran functions, treat it as a global variable
	Pop = P;


	///// INITIAL GENERATION (0) ///////////////////////////////////////
	
	// trees before parameters insertion		
	printf ("\nInitialization of the population (generation 0)\n");
	P->print_population_without_parameters(0);
	

	/// split the data set in tuning set (data_tuning) and validation set (data_validation).
	// See SPLIT and VALIDATING_LINES in input file
	// this function will also allow to increase the number of fitness cases during the run...
	// IMPORTANT! Check that the correct Sy is used in computing RMSE (note taken on 10/5/2014)
	P->split_data(&parameters, &problem, 0,1);
	
	/////////// OTHER GENERATIONS ///////////////////////////////////////
	int check_end = 0;
	int last_gen=0;
	for (int i=0; i<parameters.G+1; i++) {
		
		if (i) {   //skip generation 0
			// split the data for the current generation (this function will also allow to increase the number of fitness cases during the run...)
			//P->split_data(i,G,split); // not used...
/*
			if ((i%6)==0) {   //6
				// KILLING and FILLING
				P->kill_and_fill(&problem);
				//cin.get();	
			}
			else
//*/
				// GENETIC OPERATORS: sorting, reproduction, crossover, mutation, pruning
				P->new_spawn(parameters, problem, parameters.nfitcases,i);

			// print population WITHOUT parameters after genetic operations
			//if (COMMENT) {
			//	printf("\n\n***********Generation %d after genetic operations (not sorted, new trees marked with f9.999999E+99)************\n", i-1);
			//	P->print_population_without_parameters(i-1);
			//	printf("\n***********************************************************************\n");
			//}
		}

		// evaluate fitness function (in structural GP parameters are added and tuned first, then the evaluation is performed)
		P->evaluate(i,parameters.G);  
		

			//----------------------------------------------------------------------------------
			// extfitness - keeping the individual with least fitness value -------
			// sort according to fitness (error)
			//P->sort(i,tree_comp_fitness);

			///printf("\n\n***********Generation %d after genetic operations (not sorted, new trees marked with f9.999999E+99)************\n", i-1);
			//P->print_population_without_parameters(i-1);
			//printf("\n***********************************************************************\n");

			// update the best individual - structure and complete tree (for PARAMETER INHERITANCE)
			//P->update_ext_archive();
			//---------------------------------------------------------------------------------------------


		//// temporary
		//cout << "\n+++++++++++++++ POPULATION EVALUATED BUT NOT SORTED+++++++++++++";
		//P->print_population_without_parameters(i);
		//cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
		// sort according to F
		// VITAL! Both populations must be sorted, trees[] and complete_trees[]...
		P->sort(i,tree_comp_F);
	 	
			///printf("\n\n***********Generation %d after genetic operations (not sorted, new trees marked with f9.999999E+99)************\n", i-1);
			///P->print_population_without_parameters(i-1);
			///printf("\n***********************************************************************\n");

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
		check_end=P->terminate(parameters.threshold);
		last_gen = i;   

		// for the split data set, re-tune and re-evaluate the individuals on the merged data set
	 	if (parameters.split) {
			if ((check_end) || (i==parameters.G)) {
				cout << "\nBest Individual re-tuning and re-evaluation on the whole dataset" << endl;
				P->split_data(&parameters, &problem,last_gen,last_gen);
				P->evaluate(i,parameters.G);  
			}
		}

		// PRINT TO FILE OPERATIONS (in case of crash, data is saved)  -------------------------
		// print to file (stream of data to file opened and closed within the function)
		pop_reporter.stats2file(&parameters, P, DIR_OUTPUT, i, check_end);
		// write the test-points (training set), and the corresponding values of g_obj and the best individual (only one)
		pop_reporter.points2file(&parameters, &problem, P, DIR_OUTPUT, i, check_end, start, finish, delta_t, parameters.seed);
		// write best individual's expression
		pop_reporter.update_best2file_build(P, DIR_OUTPUT, i, check_end);
		// write the list of the best-so-far individuals (see elite or archive) - truncation!
		pop_reporter.archive2file_build(P, DIR_OUTPUT, i, check_end);
		// write no of tree evaluations
		pop_reporter.n_tree_eval2file(P, DIR_OUTPUT, i, check_end);
		// write adaptive approach data
		pop_reporter.adaptive_gen_ops_data2file(P, DIR_OUTPUT, i, check_end);
		// -------------------------------------------------------------------------------------
		
		if (check_end)
			break;
	
		// update genetic operators rates (adaptive approach)
		if (i)
			P->adapt_genetic_operators_rates();
  
	}
	// end evolution


	//termination criterion (successful run) satisfied
	if (check_end) {
		cout << "Termination criterion satisfied (RMSE < " << parameters.threshold << ")." << endl;
		cout << "Possible solution: " << endl; 
		P->print_population_with_parameters(last_gen);
		cout << "Check latest_archive.txt for solutions\n" << endl;

	}
	else {
		P->print_population_with_parameters(last_gen);
	}

	// just for test
	//P->get_tree_derivative_given_norm_vector(problem, P->complete_trees[0]);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	//---------------------------------------------------------------------
	// END OPERATIONS
	//---------------------------------------------------------------------
	// write node statistics to file
	pop_reporter.node_stats2file(&parameters, P, DIR_OUTPUT);
	
	// evaluate fitness (RMSE and R2) on test data set (only if test data has been provided)
	if (argc==4) {
		// evaluate complete individuals on test data set provided by the user
		P->evaluate_complete_trees(); // SET CORRECTLY problem.data_test, n_test, Sy_test after implementing function to read test data set
		// sort according to error (RMSE) - better not to use it to keep order and to recognise performance on building and test data sets...
		//P->sort(last_gen,tree_comp_fitness); // non va: ordina in ordine decrescente e alcune volte pone a 0 RMSE e R2. Perché?
		// find and print best individual on test data set to file
		pop_reporter.best2file_test(P, DIR_OUTPUT, last_gen);
		// print archive evaluated on the test data set to file
		pop_reporter.archive2file_test(P, DIR_OUTPUT, last_gen);  // insert a function to order in rmse decreasing order, leaving however the name of the run
	}

	// free memory allocated to Population, as not used anymore
	delete P;

	// PARALLELISATION ... TO HERE.

	cout << "\n\n END" << endl;

	return 0;
}

// to get rid of the following source files inclusions, compile their source files separately and link them in makefile!!!!
#include "./genetic_code/classes/class_BINARY_NODE.cpp"
#include "./genetic_code/classes/class_UNARY_NODE.cpp"
#include "./genetic_code/classes/class_TERMINAL_VAR.cpp"
#include "./genetic_code/classes/class_TERMINAL_CONST.cpp"
#include "./genetic_code/read_input/read_file_new.cpp"
#include "./genetic_code/read_input/read_test_data.cpp"
#include "./genetic_code/read_input/show_loaded_data.cpp"
//#include "./genetic code/SQP/MINL2.cpp"    - only for optimizer translated in C++
//#include "./genetic code/SQP/TINL2_mod.cpp"   - only for optimizer translated in C++
#include "./genetic_code/classes/class_POPULATION.cpp"
#include "./genetic_code/tree_functions/tree_operations.cpp"
//for Andrey's method - TI0L2 and MI0L2 - IT WORKS PERFECTLY
#include "./genetic_code/SQP/MI0L2_c/fdf_c.cpp"
//for Umberto's method - TINL2 and MINL2 - WORKS?
//#include "./genetic code/SQP/fdfTINL2_c.cpp"  
#include "genetic_code/tree_functions/vector_derivative_functions.cpp"
#include "./genetic_code/classes/run_parameters.cpp"
#include "./genetic_code/classes/problem_definition.cpp"
#include "./genetic_code/classes/reporter.cpp"
