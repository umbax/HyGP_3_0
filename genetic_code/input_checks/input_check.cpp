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


#include "./input_check.h"



void input_check(RunParameters *pr, ProblemDefinition *pb)
{
	// checks on method
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
	
	if (fabs((pr->repr_rate+pr->cross_rate+pr->mut_rate)-1.0) > 1.0e-5) {
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

	// check strategy for statistical properties
	if ( (pr->strat_statp!=4) && (pr->strat_statp!=6) && (pr->strat_statp!=7) && (pr->strat_statp!=8) && (pr->strat_statp!=9) && (pr->strat_statp!=11) && (pr->strat_statp!=13) && (pr->strat_statp!=14) && (pr->strat_statp!=15) && (pr->strat_statp!=16)) {
		cerr << "\nError: the strategies currently available are STRAT_STATP= 4, STRAT_STATP= 6, STRAT_STATP= 7, STRAT_STATP= 8, STRAT_STATP= 11, STRAT_STATP= 13, STRAT_STATP= 14, STRAT_STATP= 15, , STRAT_STATP= 16!!!\n";
		exit(-1);
	}

	// ATTENTION pr->strat_statp==12 only possible for nvar==1

	// check on penalisation weights
	// statistical properties objective  (a8)
	if ((pr->w_strat_statp<0) || (pr->w_strat_statp>1.0)) {
		cerr << "\nError: the weight linked to statistical properties objective (W_STRAT_STATP) has to be between 0 and 1!!!\n";
		exit(-1);
	}
	// complexity (a2)
	if ((pr->w_complexity<0) || (pr->w_complexity>1.0)) {
		cerr << "\nError: the weight on complexity (W_COMPLEXITY) has to be between 0 and 1!!!\n";
		exit(-1); 
	}
    // w_n_corrections (a3)
	if ((pr->w_n_corrections<0) || (pr->w_n_corrections>1.0)) {
			cerr << "\nError: the weight on corrections (W_N_CORRECTIONS) has to be between 0 and 1!!!\n";
			exit(-1);
	}
	// w_size (a4)
	if ((pr->w_size<0) || (pr->w_size>1.0)) {
			cerr << "\nError: the weight on size (W_SIZE) has to be between 0 and 1!!!\n";
			exit(-1);
	}
	// w_pen_ord0 (a5)
	if ((pr->w_pen_ord0<0) || (pr->w_pen_ord0>1.0)) {
			cerr << "\nError: the weight on 0-order inequality constraints (W_PEN_ORD_0) has to be between 0 and 1!!!\n";
			exit(-1);
	}
	// w_pen_ord1 (a6)
	if ((pr->w_pen_ord1<0) || (pr->w_pen_ord1>1.0)) {
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

	// check on crossvalidation
	if ((pr->crossvalidation<0) || (pr->crossvalidation>1)) {
		cerr << "\nError: Crossvalidation can be enabled (CROSSVALIDATION = 1) or disabled (CROSSVALIDATION = 0). No other values!!!\n";
		exit(-1); 
	}
	
	if (pr->crossvalidation==1) {
		if ((pr->folds_n>pr->nfitcases) || (pr->folds_n<2)) {
			cerr << "\nError: the number of folds cannot be larger than the total no of fitness cases nor smaller than 2!!!\n";
			cerr << "0 < validating_lines < nfitcases\n";
			exit(-1); 
		}
	}	

	// check on loaded primitives
	if ((pb->num_b_funcs==0) && (pb->num_u_funcs==0)) {
		cerr << "\nError: NO available functions!!! num_b_funcs = " << pb->num_b_funcs << ", num_u_funcs = " << pb->num_u_funcs << "\n";
		exit(-1); 
	}
	

}
