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
 * run_parameters.cpp
 *
 *  Created on: Jan 31, 2011
 *      Author: cnua
 */

#include "./run_parameters.h"


// RunParameters constructor
RunParameters::RunParameters(void)
{
	// all parameters initialised with 0 are loaded from input file (see ::read_file_new.cpp)
	seed=0;
	nvar=0;				// number of variables
	minrand=0.0;		// lower bound of the interval in which random numbers are selected from
	maxrand=0.0;		// upper bound of the interval in which random numbers are selected from
	step=0.01;			// delta used to define the resolution of the interval in which random numbers are selected from
	max_n_periods=0.0;	// maximum number of complete oscillations in each input space dimension
	nfitcases=0;		// total number of fitness cases (no of rows of the sample matrix)
	method=0;			// method for generating the initial pop. -> 1=no limits, 2=FULL, 3=GROW, 4=RAMPED
	depth_max=0;		// maximum depth of newly generated tree (must be > 0 or > depth_min)
	depth_min=0;		// minimum depth of newly generated tree (must be > 0)  --- only for RAMPED method
	depth_lim=0;		// limit on depth during the evolution (no trees will ever reach a depth bigger than depth_lim)
	p_FULL=0.0;			// % of trees built with FULL method (p_FULL + p_GROW = 1)    --- only for RAMPED method
	repr_rate=0.0;		// reproduction rate
	cross_rate=0.0;	 	// crossover rate
	mut_rate=0.0;		// mutation rate
	comp_rate=0.0;		//percentage of the individuals created composing two trees, with respect to the individual killed
	new_rate=0.0;		//percentage of the new individuals randomly initialized, with respect to the individual killed
	M=0;				// size of the population - number of trees
	G=0;				// number of generations
	bounded = 0;

	normalised=0;		// 7/11/20 PARAMETER NOT USED ANY LONGER: 1 if normalised RMSE is used; 0 for normal RMSE (see Population::fitness_func)
	minmax=0;			// 7/11/20 PARAMETER NOT USED ANY LONGER: 1 if MinMax F function, 0 if usual weighted function

	pso_nparticles=0;
	pso_niterations=0;

	strat_statp=0;			// Strategy to handle statistical moments as objectives - see Population::aggregate_F
	w_strat_statp=0.0;		// Numerical weight of strategy to handle statistical moments as objectives - see Population::aggregate_F
	w_ACF=0.0;				// weight related to autocorrelation finction (ACF) properties - only used for 1D problems
	w_tvariation=0.0;			// weight related to total variation
	w_complexity=0.0;		//weight given to the second normalized objective (measure of complexity of the tree)
	w_n_corrections=0.0;	//weight given to the number of times protected division "corrects" the result (division by zero)
	w_size=0.0;				//weight given to the number of nodes a tree is made of
	w_factorisation=0.0;  	//penalisation for lack of factorisation
	n_guesses=0;			// no. of guesses for the SQP optimization
	crossvalidation=0;		// 1 if data set is split each generation, 0 if not
	folds_n=0;				//number of rows containing validating points (from the upper row of the matrix)
	threshold=0.0;				//value of fitness (RMSE) under which the evolution stops
	n_inequality0=0;  	 	//total number of order 0 inequality constraints (number of points constrained)
	w_pen_ord0=0.0;  	//weight for the penalisation on order 0 inequality constraint
	n_inequality1=0;    	//total number of order 1 inequality constraints (number of points constrained)
	w_pen_ord1=0.0;    //weight for the penalisation on order 1 inequality constraint
}


void RunParameters::show(void)
{
	cout << "\n\nRunParameters::show(void) : data imported in RunParameters object: " << endl;
    cout << "seed = " << seed << endl;
	cout << "nvar = " << nvar << endl;
  	cout << "minrand= " << minrand << endl;
	cout << "maxrand= " << maxrand << endl;
	cout << "max_n_periods= " << max_n_periods << endl;
	cout << "pso_nparticles= " << pso_nparticles << endl;		// "step= " << step << endl;
	cout << "pso_niterations= " << pso_niterations << endl;
	cout << "nfitcases = " << nfitcases << endl;
	cout << "method = " << method << endl;
	cout << "depth_max = " << depth_max << endl;
 	cout << "depth_min = " << depth_min << endl;
	cout << "depth_lim = " << depth_lim << endl;
   	cout << "p_FULL = " << p_FULL << endl;
	cout << "repr_rate = " <<  repr_rate << endl;
	cout << "cross_rate = " <<  cross_rate << endl;
	cout << "mut_rate = " <<  mut_rate << endl;
	cout << "comp_rate = " << comp_rate << endl;
	cout << "new_rate = " << new_rate << endl;
	cout << "M = " << M << endl;
	cout << "G = " << G << endl;
	cout << "bounded = " << bounded << endl;
	cout << "strat_statp = " << strat_statp << endl;
	cout << "w_strat_statp = " << w_strat_statp << endl;
	cout << "w_ACF = " << w_ACF << endl;
	cout << "w_tvariation = " << w_tvariation << endl;
	cout << "w_complexity = " << w_complexity << endl;
	cout << "w_n_corrections = " << w_n_corrections << endl;
	cout << "w_size = " << w_size << endl;
	cout << "n_guesses = " << n_guesses << endl;
	cout << "crossvalidation = " << crossvalidation << endl;
	cout << "folds_n = " << folds_n << endl;
	cout << "threshold = " << threshold << endl;
	cout << "n_inequality0 = " << n_inequality0 << endl;
	cout << "w_pen_ord0 = " << w_pen_ord0 << endl;
	cout << "n_inequality1 = " << n_inequality1 << endl;
	cout << "w_pen_ord1 = " << w_pen_ord1 << endl;
}
