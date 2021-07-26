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

#include <iostream>  // basic i/o commands: cout, cin, scientific, fixed, cerr
#include <iomanip>  // manipulators (text format): setw
#include <string>    // to manipulate strings (C)
#include <cstring>   // to manipulate strings (C) (strcpy, strcmp)
#include <cstdlib>   // NULL, exit, EXIT_FAILURE
#include <cmath>	// pow, sqrt, abs
#include <algorithm> // max_element, min_element
#include <exception>
#include <vector>

//#include "../tree_functions/vector_derivative_functions.h"
#include "../tree_functions/tree_operations.h"
#include "../nodes/nodes.h"
#include "../modules/func_primitives_prototypes.h"
#include "../modules/primitives.h"
#include "../modules/variable.h"
#include "../classes/problem_definition.h"
#include "../classes/run_parameters.h"
#include "../classes/class_POPULATION.h"

//#include "./hypso_launcher.h"

using namespace std;

Population *Pop; // try! Simple way to make fdf_c able to get the values of the tree

int main ()
{
	// test HyPSO algorithm to get the target model y=sin(x)+cos(5*x)
	// define data
	int n_var=1;
	int rows=101;  // number of building records (building data set size)
	int cols=2;
	Val** p_data;
	p_data = new Val*[rows];
	double u_bound=6.40;
	for (int k=0; k<rows; k++) {
		p_data[k]=new Val[cols];
		//Z1
		p_data[k][0]=(u_bound/(rows-1))*k;
		cout << "\np_data[" << k << "][0]= " << p_data[k][0];  // x
		// target
		p_data[k][1]=sin(p_data[k][0])+cos(5*p_data[k][0]);  // y=sin(x)+cos(5*x)
		cout << "  p_data[" << k << "][1]= " << p_data[k][1];
	}
	cin.get();
	//Val** data_valid;



	// define variables and ProblemDefinition
	Variable* Z;
	// define ProblemDefinition
	ProblemDefinition pb;
	pb.set_data(p_data);
	pb.set_n_var(n_var);
	pb.set_n_cols(cols);	//number of columns in data (n_var+1)
	pb.set_n_data(rows);	//total number of rows in data
	pb.initialise_variables(&Z, 2);
	pb.b_func_list[0]=&Add;         // check dummy_bin and the use of pointers... there must be something wrong
	pb.b_func_list[1]=&Sub;
	pb.b_func_list[2]=&Mult;
	pb.b_func_list[3]=&SDiv;
	pb.num_b_funcs = 4; // no. of binary operators
	pb.u_func_list[0]=&Sin;
	pb.u_func_list[1]=&Cos;
	pb.num_u_funcs = 2; // no. of unary operators
	// Sy=sum((output- average_output)^2)=(n-1)*output variance on the whole dataset
	pb.sum_output=0.0;
	for (int i=0; i<pb.get_n_data(); i++) {
		pb.sum_output=pb.sum_output+pb.get_data(i,1);   // sum of each target value
	}
	cout << "\npb.sum_output= " << pb.sum_output << endl;
	pb.y_ave=pb.sum_output/pb.get_n_data(); 	//average value of input data
	pb.Sy=0.0;
	for (int i=0; i<pb.get_n_data(); i++) {
		pb.Sy=pb.Sy + pow((pb.get_data(i,1)-pb.y_ave),2.0);   // SStot total sum of squares of (observed data - average observed data)
	}
	cout << "\npb.Sy= " << pb.Sy << endl;
	pb.data_validation=p_data;
	pb.n_validation=rows;

	// define Parameters
	RunParameters pr;
	pr.seed = 1;
	pr.nvar = 1;
	pr.minrand = 0;
	pr.maxrand = 3000;
	pr.step = .001;
	pr.max_n_periods = 2;
	pr.nfitcases = 1;
	pr.method = 4;
	pr.depth_max = 3;
	pr.depth_min = 2;
	pr.depth_lim = 7;
	pr.p_FULL = 0.5;
	pr.repr_rate = .2;
	pr.cross_rate = .4;
	pr.mut_rate = .4;
	pr.comp_rate = 0;
	pr.new_rate = 0;
	pr.M = 1;  //4
	pr.G = 1;  //2
	pr.normalised=0;
	pr.minmax=0;
	pr.w_complexity = 0;
	pr.w_n_corrections = 0;
	pr.w_size = 0;
	pr.w_factorisation = 0;  // switch between standard approach (<0) and factorisation bonus (>0)


	// define Population
	Population *P = new Population(&pr, &pb);     //Population P(&pr, &pb);
	//as it's hard to pass Pop as a parameter to fdf_c__ through fortran functions, treat it as a global variable
	Pop = P;


	// define tree sin(3*Z1)+cos(2*Z1)
	Binary_Node T7_tree(NULL, &Add);
	T7_tree.n_corrections=0;
	// D1 (depth 1)
	Unary_Node T7_D1_P1(&T7_tree, &Sin);  // left child, depth 1 : sin()
	T7_tree.set_left((Node *)&T7_D1_P1);
	Unary_Node T7_D1_P2(&T7_tree, &Cos);  // right child, depth 1 : cos()
	T7_tree.set_right((Node *)&T7_D1_P2);
	// D2 (depth 2)
	Binary_Node T7_D2_P1(&T7_D1_P1,&Mult);  // child of sin() : 3*Z1
	T7_D1_P1.set_child((Node *)&T7_D2_P1);
	Binary_Node T7_D2_P2(&T7_D1_P2, &Mult);	// child of cos() : 2*Z1
	T7_D1_P2.set_child((Node *)&T7_D2_P2);
	// D3 (depth 3)
	Terminal_Const T7_D3_P1(&T7_D2_P1, 3.0);		// left child of mult in sin: 3
	T7_D2_P1.set_left((Node *)&T7_D3_P1);
	Terminal_Var T7_D3_P2(&T7_D2_P1,pb.v_list[0]);  // right child of mult in sin: Z1
	T7_D2_P1.set_right((Node *)&T7_D3_P2);

	Terminal_Const T7_D3_P3(&T7_D2_P2, 2.0);		// left child of mult in cos: 2
	T7_D2_P2.set_left((Node *)&T7_D3_P3);
	Terminal_Var T7_D3_P4(&T7_D2_P2,pb.v_list[0]);  // right child of mult in cos: Z1
	T7_D2_P2.set_right((Node *)&T7_D3_P4);

	// extract numerical coefficients
	int n_param = T7_tree.find_parameters();
	cout << "\nFind_parameters(): n_param=" << n_param << endl;
	double x[2]={T7_D3_P1.value(0), T7_D3_P3.value(0)};  // use in the future p_par

	cout << "\nTree : ";
	P->print_individual((Node *)&T7_tree);  // till here ok
	cout << "\nNumerical coefficients (original values) : ";
	cout << "x[0]=" << x[0] << " x[1]=" << x[1] << endl;

	cin.get();

	// define HyPSO parameters
	int spacedim=2;  //2 as 2 are the parameters to be optimised (x is an array of size 2)

	// double Population::pso_objfunction(double* x, int n_param, Binary_Node *ntree)
	cout << "\nOriginal tree: ";
	P->print_individual((Node *)&T7_tree);
	P->hypso_launcher(&T7_tree, spacedim, x);
	// print optimised model
	cout << "\nTarget model: y=sin(x)+cos(5*x)";
	cout << "\nOptimised tree: ";
	P->print_individual((Node *)&T7_tree);
	cout << endl;
	cout << endl;

}

//for Andrey's method - TI0L2 and MI0L2 - IT WORKS PERFECTLY
#include "../SQP/MI0L2_c/fdf_c.cpp"
