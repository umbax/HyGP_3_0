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

#include <iostream>
#include <stdlib.h>   //exit

using namespace std;

// dependencies
#include "../classes/run_parameters.h"
#include "../classes/problem_definition.h"


void input_check(RunParameters *pr, ProblemDefinition *pb)
{
	if ((pr->method<1) || (pr->method>4)) {
		cerr << "\nError: initialization : it MUST be  1<=method<=4 !!!\n";
		exit(-1); 
	}

	if (pr->method==4) {
		if (pr->depth_min<1) {
			cerr << "\nError: RAMPED initialization : depth_min MUST be > 1 !!!\n";
			exit(-1); 
		}
		if (pr->depth_min>pr->depth_max) {
			cerr << "\nError: RAMPED initialization : depth_max MUST be > depth_min !!!\n";
			exit(-1); 
		}
	}
 
	if ((pr->method==2) || (pr->method==3)) {
		if (pr->depth_max<1) {
			cerr << "\nError: initialization : depth_max MUST be > 0 !!!\n";
			exit(-1);
		}
	}
	
	if (abs((pr->repr_rate+pr->cross_rate+pr->mut_rate)-1.0) > 1.0e-5) {
		cerr << "\nError: repr_rate+cross_rate+mut_rate = " << pr->repr_rate+pr->cross_rate+pr->mut_rate << " !!!";
		cerr << "\nError: the sum (repr_rate+cross_rate+mut_rate) MUST be equal to 1 !!!\n";
		exit(-1); 
	}
	
	// check on population size
	if (pr->M<0) {
		cerr << "\nError: the size of the population (M) is negative!!!\n";
		exit(-1); 
	}

	// check on number of generations
	if (pr->G<0) {
		cerr << "\nError: the number of generations (G) is negative!!!\n";
		exit(-1); 
	}

	// check on penalisation weights
	// complexity (a2)
	if ((pr->w_complexity<0) || (pr->w_complexity>1)) {
		cerr << "\nError: the weight on complexity (W_COMPLEXITY) has to be between 0 and 1!!!\n";
		exit(-1); 
	}
    // w_n_corrections (a3)
	if ((pr->w_n_corrections<0) || (pr->w_n_corrections>1)) {
			cerr << "\nError: the weight on corrections (W_N_CORRECTIONS) has to be between 0 and 1!!!\n";
			exit(-1);
	}
	// w_size (a4)
	if ((pr->w_size<0) || (pr->w_size>1)) {
			cerr << "\nError: the weight on size (W_SIZE) has to be between 0 and 1!!!\n";
			exit(-1);
	}
	// w_pen_ord0 (a5)
	if ((pr->w_pen_ord0<0) || (pr->w_pen_ord0>1)) {
			cerr << "\nError: the weight on 0-order inequality constraints (W_PEN_ORD_0) has to be between 0 and 1!!!\n";
			exit(-1);
	}
	// w_pen_ord1 (a6)
	if ((pr->w_pen_ord1<0) || (pr->w_pen_ord1>1)) {
			cerr << "\nError: the weight on 1-order inequality constraints (W_PEN_ORD_1) has to be between 0 and 1!!!\n";
			exit(-1);
	}
	// w_factorisation
	if (pr->w_factorisation>0) {
		cout << "\nW_FACTORISATION = " << pr->w_factorisation << " > 0 : Fitness function : Factorisation approach enabled\n";
	} else {
		cout << "\nW_FACTORISATION = " << pr->w_factorisation << " <= 0 : Fitness function : Standard approach enabled\n";
	}

	// check sum of weights
	double tot = pr->w_complexity + pr->w_n_corrections + pr->w_size + pr->w_pen_ord0+ pr->w_pen_ord1;
	if (tot>1) {
		cerr << "\nError: the sum of the weights is bigger than 1!! This way RMSE does not guide the search!!!\n";
		exit(-1);
	}

	// check on number of guesses
	if (pr->n_guesses<0) {
		cerr << "\nError: the number of initial random guesses (N_GUESSES) is negative!!!\n";
		exit(-1); 
	}

	if ((pr->split<0) || (pr->split>1)) {
		cerr << "\nError: split can be enabled (SPLIT = 1) or disabled (SPLIT = 0). No other values!!!\n";
		exit(-1); 
	}
	
	if (pr->split==1) {
		if ((pr->validating_lines>=pr->nfitcases) || (pr->validating_lines<1)) {
			cerr << "\nError: the validating entries must be fewer than the total no of fitness cases!!!\n";
			cerr << "0 < validating_lines < nfitcases\n";
			exit(-1); 
		}
	}	

	if ((pb->num_b_funcs==0) && (pb->num_u_funcs==0)) {
		cerr << "\nError: NO available functions!!! num_b_funcs = " << pb->num_b_funcs << ", num_u_funcs = " << pb->num_u_funcs << "\n";
		exit(-1); 
	}
	

}