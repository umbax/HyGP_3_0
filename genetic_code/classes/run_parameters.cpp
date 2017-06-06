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


void RunParameters::show(void)
{
	cout << "\n\nRunParameters::show(void) : data imported in RunParameters object: " << endl;
    cout << "seed = " << seed << endl;
	cout << "nvar = " << nvar << endl;
  	cout << "minrand= " << minrand << endl;
	cout << "maxrand= " << maxrand << endl;
	cout << "max_n_periods= " << max_n_periods << endl;
	cout << "step= " << step << endl;
	cout << "nfitcases = " << nfitcases << endl;
	cout << "method = " << method << endl;
	cout << "depth_max = " << depth_max << endl;
 	cout <<  "depth_min = " << depth_min << endl;
	cout <<  "depth_lim = " << depth_lim << endl;
   	cout << "p_FULL = " << p_FULL << endl;
	cout << "repr_rate = " <<  repr_rate << endl;
	cout << "cross_rate = " <<  cross_rate << endl;
	cout << "mut_rate = " <<  mut_rate << endl;
	cout << "comp_rate = " << comp_rate << endl;
	cout << "new_rate = " << new_rate << endl;
	cout << "M = " << M << endl;
	cout << "G = " << G << endl;
	cout << "normalised_err = " << normalised << endl;
	cout << "minmax = " << minmax << endl;
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
