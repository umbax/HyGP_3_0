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
 * parallel_master.cpp
 *
 * this is the parallelised version of master.cpp
 * A whole experiment made of n_runs is launched when parallel_master is executed.
 * The parallelisation is done using OpenMP.
 * Each run is allocated a thread:
 * ideally, the number of runs should match the number of cores
 *
 * parallel use:
 * ./gp  location/name_input location/name_output_folder  n_runs p
 *
 * sequential use:
 * ./gp  location/name_input location/name_output_folder  n_runs s
 *
 *  Created on: Mar 1, 2011
 *      Author: Umberto Armani
 */

// use: ./parallel_gp_exp  location/name_input location/name_output  n_runs

#include <iostream>  // basic i/o commands
#include <fstream>  // file I/O
#include <iomanip>  // manipulators (text format)
#include <string>    // to manipulate strings (C)
#include <cstring>    // to manipulate strings (C)
#include <string>    //  string class (C++)
#include <cstdlib>
#include <ctime>	 // to work with variables of type time_t and random number generator "srandom"
#include <cmath>
#include <algorithm>
#include <vector>
#include <exception>
#include <sstream>    // for stringstream

#include <sys/stat.h>
#include <sys/types.h>
#include <omp.h>   // to be added only in mpi computers

using namespace std;

// global variables - numeric settings
const double PI = 4.0*atan(1.0);
//double MAX_VAL = 1.8e+19;   // can be defined in Val_type.h
//double MIN_VAL = 1.8e-19;		// can be defined in Val_type.h

// HEADERS
// do not change the order of the following include commands!
//Otherwise you may have undeclared objects //
#include "./genetic_code/I_O_functions/split_file_name.cpp"
#include "./genetic_code/I_O_functions/copy_file.cpp"

#include "./genetic_code/modules/nodes_types.h"
#include "./genetic_code/modules/Val_type.h"
#include "./genetic_code/modules/variable.h"
#include "./genetic_code/modules/func_primitives_prototypes.h"
#include "./genetic_code/classes/run_parameters.h"
#include "./genetic_code/classes/problem_definition.h"
#include "./genetic_code/read_input/show_loaded_data.h"
#include "./genetic_code/read_input/read_file_new.h"
#include "./genetic_code/input_checks/input_check.h"
#include "./genetic_code/classes/class_NODE_base.h"
#include "./genetic_code/tree_functions/tree_operations.h"
#include "./genetic_code/tree_functions/vector_derivative_functions.h"
#include "./genetic_code/modules/primitives.h"
#include "./genetic_code/classes/class_BINARY_NODE.h"
#include "./genetic_code/classes/class_UNARY_NODE.h"
#include "./genetic_code/classes/class_TERMINAL_VAR.h"
#include "./genetic_code/classes/class_TERMINAL_CONST.h"
#include "./genetic_code/classes/class_POPULATION.h"
#include "./genetic_code/classes/reporter.h"
#include "./genetic_code/single_run.h"
#include "./genetic_code/run_seeding/seed_generator.h"

extern "C" {extern void opti_(int*, double*, int*, int*, int*);}  //for Andrey's method - TI0L2 and MI0L2 - IT WORKS PERFECTLY
extern "C" {int  fdf_c__ (double *,  double *, int *,int *); }


// primitives.h : list of primitives, made of functions declarations but also of variables (that are global).
// In the future it may be a good idea to use threadprivate to protect these functions
// and variables in each thread. Otherwise, use pragma omp critical

// global variables - population purposes
//-------------------------------------------------------------------------------------------------------------------------
// list of Populations (one for each run)
Population *Pop;// = NULL; //  Simple way to make fdf_c able to get the values of the tree
# pragma omp threadprivate(Pop)

int VERBOSE = 0;     //set to 1 if you want to print on the screen all the comments! 0 for a "clean" and wordless execution...
//-------------------------------------------------------------------------------------------------------------------------



int main (int argc, char *argv[])
{
	// check the number of arguments
	if (argc!=5) {
		cerr << "\nERROR!!! Too few/ too many arguments!!!"	;
		cerr << "\nUSAGE: >> ./gp  location/name_input  existing_directory_output n_runs  mode";
		cerr << "\nExample: >> ./gp ./input/input_file.txt ./output 10 p\n";
		exit(-1);
	}
	if (strcmp(argv[4],"p") && strcmp(argv[4],"s")) {
		cerr << "\nERROR : mode must be p (parallel) or s (sequential)" << endl;
		exit(-1);
	}



	cout << "\n\n FIRST SIMPLE PARALLEL GP IMPLEMENTATION\n";

	// input arguments and environment variables
	string INPUT_STRING=argv[1];
	string INPUT_DIR, INPUT_FILE;
	split_file_name(INPUT_STRING, INPUT_DIR, INPUT_FILE);

	string OUTPUT_STRING=argv[2];
	string DESTINATION, OUTPUT_DIR;
	split_file_name(OUTPUT_STRING, DESTINATION, OUTPUT_DIR);

	int N_RUNS = 0;
	sscanf(argv[3], "%d", &N_RUNS);

	int PARALLEL_EX;
	if (!strcmp(argv[4],"p"))
		PARALLEL_EX = 1;
	else
		PARALLEL_EX = 0;

	cout << "\nINPUT_STRING = " << INPUT_STRING;
 	cout << "\nINPUT_DIR = " <<  INPUT_DIR;
 	cout << "\nINPUT_FILE = " <<  INPUT_FILE;

 	cout << "\nOUTPUT_STRING = " << OUTPUT_STRING;
 	cout << "\nDESTINATION = " << DESTINATION;
 	cout << "\nOUTPUT_DIR = " << OUTPUT_DIR;

 	cout << "\nN_RUNS = " << N_RUNS << endl;

 	// create directory OUTPUT_DIR
 	mkdir(OUTPUT_STRING.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
	// copy input file in DESTINATION/OUTPUT_DIR and reassign INPUT_FILE
 	string INPUT_FILE_PATH =  OUTPUT_STRING + "/input_file.txt";
 	copy_file(INPUT_STRING, INPUT_FILE_PATH);


 	// MANAGE INPUT FILE and DATA ACQUISITION
	// other classes instantiations
	RunParameters parameters;
	ProblemDefinition problem;
	Reporter pop_reporter;

	// read input and check for errors in loaded data
	read_input_file(INPUT_FILE_PATH, &parameters, &problem);
	// check for errors in loaded data
	input_check(&parameters, &problem); //pass by reference, otherwise copy constructor is called...
	// show the results
	show_loaded_data(&parameters, &problem);


	// RANDOM SEEDS GENERATOR
	cout << "\nSetting seeds for random number generator...";
	int* SEED;
	seed_generator(&SEED, N_RUNS, parameters.seed);


	// TIME VARIABLE DECLARATION
	time_t start, finish;
	double delta_t;

	//SETTING PRECISION (number of decimal figures)
	cout.precision(5);

	// PARALLELISATION. FROM HERE ....
	#pragma omp parallel  if (PARALLEL_EX) num_threads(N_RUNS) \
				default(none) \
				shared(cout, cin) \
				private(finish, delta_t) \
				firstprivate(start,SEED, parameters, problem, pop_reporter,\
									N_RUNS,OUTPUT_STRING)  //,MAX_VAL,MIN_VAL)
				// firstprivate created copies of the variable using copy constructor for ProblemDefinition
	{
		# pragma omp for
		for (int cur_run=0; cur_run < N_RUNS; cur_run++) {

			int id = omp_get_thread_num();
			int tot_threads =  omp_get_num_threads();

			# pragma omp critical
			{
				cout << "\n Single run : id = " << id << " tot_threads = " << tot_threads;
				//cout << "\nData passed to the run:";
				//show_loaded_data(&parameters, &problem);
				//cout << "\nPop = " << Pop;
			}


			// in the future pass parameters and problem by reference (address) so a copy constructor is not called...
			int r = single_run(id,cur_run, SEED, OUTPUT_STRING,\
					&parameters, &problem, pop_reporter,\
					start, finish, delta_t);


		// end #pragma omp for
		}

	// PARALLELISATION ... TO HERE.
	// end #pragma omp parallel
	}

	// let all the threads catch up...
	# pragma omp barrier

	// operations to collect data ? Now done by ./posteriori
	// for now just create a file to tell that the execution has finished
	string SUCCESS_FILE =  OUTPUT_STRING + "/end_OK.txt";
	ofstream fout;
	fout.open(SUCCESS_FILE.c_str(), ios_base::out | ios_base::trunc);
	fout << N_RUNS << " runs have been successfully completed." << endl;
	fout << "End.";
	fout.close();

	cout << "\n\n END" << endl;

 	return 0;
}


#include "./genetic_code/read_input/read_file_new.cpp"
#include "./genetic_code/read_input/show_loaded_data.cpp"
#include "./genetic_code/classes/class_POPULATION.cpp"
#include "./genetic_code/tree_functions/tree_operations.cpp"
//for Andrey's method - TI0L2 and MI0L2 - IT WORKS PERFECTLY
#include "./genetic_code/SQP/MI0L2_c/fdf_c.cpp"
#include "./genetic_code/tree_functions/vector_derivative_functions.cpp"
#include "./genetic_code/classes/run_parameters.cpp"
#include "./genetic_code/classes/problem_definition.cpp"
#include "./genetic_code/classes/reporter.cpp"
#include "./genetic_code/single_run.cpp"
#include "./genetic_code/run_seeding/seed_generator.cpp"
