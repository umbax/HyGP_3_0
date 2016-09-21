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



#ifndef RUNPARAMETERS_H_
#define RUNPARAMETERS_H_

using namespace std;

// dependencies
#include "../modules/Val_type.h"

// class containing the parameters defining the evolution
class RunParameters
{
	public:
			int seed;
			int nvar;				// number of variables
			Val minrand;	
			Val maxrand;	
			Val step;
			double max_n_periods;	// maximum number of complete oscillations in each direction
			int nfitcases;			// total number of fitness cases (no of rows of the sample matrix)
			int method;				// method for generating the initial pop. -> 1=no limits, 2=FULL, 3=GROW, 4=RAMPED
			int depth_max;		// maximum depth of newly generated tree (must be > 0 or > depth_min)
			int depth_min;		// minimum depth of newly generated tree (must be > 0)  --- only for RAMPED method
			int depth_lim;			// limit on depth during the evolution (no trees will ever reach a depth bigger than depth_lim)
			double p_FULL;			// % of trees built with FULL method (p_FULL + p_GROW = 1)    --- only for RAMPED method
			double repr_rate;		// reproduction rate
			double cross_rate;	 	// crossover rate
			double mut_rate;		// mutation rate 
			double comp_rate;		//percentage of the individuals created composing two trees, with respect to the individual killed
			double new_rate;		//percentage of the new individuals randomly initialized, with respect to the individual killed
			int M;							// size of the population - number of trees 
			int G;							// number of generations 
			bool normalised;		// 1 if normalised RMSE is used; 0 for normal RMSE (see Population::fitness_func)
			bool minmax;				// 1 if MinMax F function, 0 if usual weighted function
			double w_complexity;		//weight given to the second normalized objective (measure of complexity of the tree)
			double w_n_corrections;	//weight given to the number of times protected division "corrects" the result (division by zero)
			double w_size;					 //weight given to the number of nodes a tree is made of
			double w_factorisation;  //penalisation for lack of factorisation
			int n_guesses;			 // no. of guesses for the SQP optimization
			int split;							// 1 if data set is split each generation, 0 if not
			int validating_lines;		//number of rows containing validating points (from the upper row of the matrix) 
			double threshold;				//value of fitness (RMSE) under which the evolution stops
			int n_inequality0;  	 	//total number of order 0 inequality constraints (number of points constrained) 
			double w_pen_ord0;  	//weight for the penalisation on order 0 inequality constraint
			int n_inequality1;    	//total number of order 1 inequality constraints (number of points constrained) 
			double w_pen_ord1;    //weight for the penalisation on order 1 inequality constraint
			
			// function to show data
			void show(void);
};

#endif  /* RUNPARAMETERS_H_ */
