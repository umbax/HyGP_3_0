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


#ifndef REPORTER_H_
#define REPORTER_H_

#include <iomanip>  // manipulators (text format): setw
//#include <sstream>    // for stringstream
#include <fstream>  // file I/O

#include "../nodes/nodes_types.h"    //just a header, no source file
#include "./run_parameters.h"
#include "./problem_definition.h"
#include "./class_POPULATION.h"


//class to manage file operations



class Reporter
{
	public:
		
		// functions to print data
		//function to print generation statistics to file  - SINGLE RUN
		void stats2file(RunParameters*, Population *, string, int, int);
		//function to print points of the best tree to file - SINGLE RUN
		void points2file(RunParameters*, ProblemDefinition*, Population* , string, int, int, time_t, time_t, double, int);
		// function to print the expression of the best individual to file - SINGLE RUN
		void update_best2file_build(Population *, string, int, int);
		// function that finds and prints the best individual on the test data set
		void best2file_test(Population *, string, int);
		//function to print all individuals in the archive with their score evaluated on training (building) data set - WHOLE RUN
		void archive2file_build(Population *, string, int, int);
		//function to print all individuals in the archive with their score evaluated on test (validation) data set - WHOLE RUN
		void archive2file_test(Population *, string, int);
		//function to print the node selection statistics to file "node_selection.txt"
		void node_stats2file(RunParameters *, Population *, string);
		// function to print the no. of tree evaluations to file
		void n_tree_eval2file(Population *, string , int , int);
		// function to print parameters related to the behaviuour of genetic operators (to assess adaptive approach)
		void adaptive_gen_ops_data2file(Population *, string, int, int);
};

#endif  /* REPORTER_H_ */
