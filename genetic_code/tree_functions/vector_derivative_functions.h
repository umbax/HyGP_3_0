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


#ifndef VECTOR_DERIVATIVE_FUNCTIONS_H_
#define VECTOR_DERIVATIVE_FUNCTIONS_H_

// included dependencies
#include <iostream>  // basic i/o commands: cout, cin, scientific, fixed, cerr
#include <cstdlib>   // NULL, exit, EXIT_FAILURE
#include <algorithm> // max_element, min_element
#include <cmath>	// pow, sqrt


int int_rand(int);
void basic_stat_analysis(double*,double*,int);
void vect_sub(double*, double*, int, double*);
void vect_sum(double*, double*, int, double*);
void vect_div_scal(double*, double, int, double*);
void vect_mult_scal(double*, double, int, double*);
double point_dist(double*, double*, int);
double vect_magn(double*, int);
int point_comp(const void *, const void *);
double fin_differences(double, double, double);
void get_deriv_to_closest_neighbour(double **, int, int, double **, int);

#endif  /* VECTOR_DERIVATIVE_FUNCTIONS_H_ */
