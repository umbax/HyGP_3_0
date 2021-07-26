/*
 * hypso_launcher.h
 *
 *  Created on: Nov 8, 2016
 *      Author: umba
 */

#ifndef GENETIC_CODE_HYPSO_HYPSO_LAUNCHER_H_
#define GENETIC_CODE_HYPSO_HYPSO_LAUNCHER_H_

#include <iostream>  // basic i/o commands: cout, cin, scientific, fixed
#include <string>    // to manipulate strings (C)
#include <cstring>   // to manipulate strings (C)
#include <cmath>
#include <stdio.h>
#include <cstdlib>   // NULL, exit, EXIT_FAILURE

#include "../nodes/nodes_types.h"
#include "../nodes/nodes.h"
#include "../classes/class_POPULATION.h"


double** zeros(int, int);
double* zeros(int);

int regenerate_on_boundary(double**, int, double**);
int regenerate_from_swarm_best(double**, int, int, double**, double**);

//int psominimize(Binary_Node *, double (*)(double*, int), int , double* , double* , int , int , double** , double* );
int psominimize(Binary_Node*, Population*, int , double* , double* , int , int , double** , double* );

void hypso_launcher(Population*, Binary_Node*, int, double*);


#endif /* GENETIC_CODE_HYPSO_HYPSO_LAUNCHER_H_ */
