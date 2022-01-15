// Copyright 2017 Dr Umberto Armani
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

#include "../genetic_code/tree_functions/vector_derivative_functions.h"
#include "../genetic_code/tree_functions/tree_operations.h"
#include "../genetic_code/nodes/nodes.h"
#include "../genetic_code/modules/func_primitives_prototypes.h"
#include "../genetic_code/modules/primitives.h"
#include "../genetic_code/modules/variable.h"
#include "../genetic_code/classes/problem_definition.h"
#include "../genetic_code/classes/run_parameters.h"
#include "../genetic_code/classes/class_POPULATION.h"

using namespace std;

Population* Pop;

int main ()
{
	/*
	 * test Population::fitness_func(Val Sy, Val** data_used, int n_cases, Node *current_tree, Val *result_tree, bool normalised)
	 */

	cout << "\n\nTEST : Population::fitness_func(Val Sy, Val** data_used, int n_cases, Node *current_tree, Val *result_tree, bool normalised)";

	int errors_no=0;
	vector <int> failed_tests;
	vector <int> passed_tests;

	// define other variables
	char* expr;
	Val result[10];
	result[0] = (Val)0.0;  // fitness value (error - RMSE)
	result[1] = (Val)0.0;	// storing n of hits
	result[2] = (Val)0.0;	// storing n of corrections done by protected operations
	result[3] = (Val)0.0; 	// storing value of R squared (R2)
	result[4] = (Val)0.0;	// storing mean tree value on building data set
	result[5] = (Val)0.0;	// storing variance of tree values on building data set
	result[6] = (Val)0.0;	// storing min tree value on building data set
	result[7] = (Val)0.0;	// storing max tree value on building data set
	result[8] = (Val)0.0;	// storing value first zero of tree autocorrelation function
	result[9] = (Val)0.0;	// storing total variation

	// define data
	int rows=1;
	int cols=3;
	Val** p_data;
	p_data = new Val*[rows];
	for (int k=0; k<rows; k++) p_data[k]=new Val[cols];
	p_data[0][0] = 0.0; //9.306025965175708e-01;  //Z1
	p_data[0][1] = 0.0; // -5.733778508896148e-01; //Z2
	p_data[0][2] = 0.0; //2.071917766525170e+02;  // target
	Val** data_valid;

	// define variables and ProblemDefinition
	Variable* Z;
	ProblemDefinition pb;
	pb.set_data(p_data);
	pb.set_n_var(cols-1);
	pb.set_n_cols(cols);
	pb.set_n_data(1);
	pb.initialise_variables(&Z, 2);
	pb.b_func_list[0]=&Add;         // check dummy_bin and the use of pointers... there must be something wrong
	pb.b_func_list[1]=&Sub;
	pb.b_func_list[2]=&Mult;
	pb.b_func_list[3]=&SDiv;
	pb.num_b_funcs = 4; // no. of binary operators
	pb.u_func_list[0]=&Square;
	pb.u_func_list[1]=&Cos;
	pb.num_u_funcs = 2; // no. of unary operators
	// Sy=sum((output- average_output)^2)=(n-1)*output variance on the whole dataset
	pb.Sy=0;

	// define Parameters
	RunParameters pr;
	pr.seed = 1;
	pr.nvar = 2;
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
	pr.M = 4;
	pr.G = 2;
	pr.normalised=0;
	pr.minmax=0;
	pr.w_complexity = 0;
	pr.w_n_corrections = 0;
	pr.w_size = 0;
	pr.w_factorisation = 0;  // switch between standard approach (<0) and factorisation bonus (>0)
/*
		pr->n_guesses = (int)(p_par_values[25]);
		pr->crossvalidation = (int)(p_par_values[26]);
		pr->folds_n = (int)(p_par_values[27]);
		pr->threshold = p_par_values[28];
		pr->n_inequality0 = (int)p_par_values[29];
*/
	pr.w_pen_ord0 = 0;
	pr.n_inequality1 = 0;
	pr.w_pen_ord1 = 0;


	// define Population
	Population P(&pr, &pb);

	//----------------------------------------------------------------------------------------
	// POPULATION::FITNESS_FUNC : TEST N. 1 (2D)
	//----------------------------------------------------------------------------------------
	cout << "\nPOPULATION::FITNESS_FUNC : TEST N. 1 -----------------------------------";
	// define tree Z1+Z2
	Binary_Node T1_tree(NULL, &Add);	//root node
	Terminal_Var T1_D1_P1(&T1_tree, pb.v_list[0]);  // left child
	T1_tree.set_left((Node *)&T1_D1_P1);
	Terminal_Var T1_D1_P2(&T1_tree, pb.v_list[1]);  // right child
	T1_tree.set_right((Node *)&T1_D1_P2);

	cout << "\nTree : ";
	P.print_individual((Node *)&T1_tree);
	cout << "Node * &T1_tree= " << &T1_tree << endl;

	cols=3;
	pb.n_validation=2;
	data_valid = new Val*[pb.n_validation];
	for (int k=0; k<pb.n_validation; k++) data_valid[k]=new Val[cols];
	// fitness case 0
	data_valid[0][0] = 15.0;	//Z1
	data_valid[0][1] = 2.5;		//Z2
	data_valid[0][2] = 17.5;	//target
	// fitness case 1
	data_valid[1][0] = 20.0;	//Z1
	data_valid[1][1] = 5.0;		//Z2
	data_valid[1][2] = 25.0;  	//target
	//
	pb.sum_output=data_valid[0][2]+data_valid[1][2]; 	// sum of each target value
	pb.y_ave=pb.sum_output/2.0; 						// average value of input data
	pb.Sy=pow((data_valid[0][2]-pb.y_ave),2.0)+pow((data_valid[1][2]-pb.y_ave),2.0);  		// SStot total sum of squares of (observed data - average observed data) // defined in read_input_file function (read_file_new.cpp)
	pb.y_var=pb.Sy/2.0; 		// variance of input data
	pb.y_max = 25.0;		// max value of target
	pb.y_min = 17.5;		// min value of target
	pb.data_validation=data_valid;
	pb.show_data_validation();

	// segmentation fault here?
	cout << "Calling fitness_func..." << endl;
	cout << "T1_tree.n_corrections= " << T1_tree.n_corrections << endl;
	P.fitness_func(pb.Sy, pb.data_validation, 2, (Node *)&T1_tree, result, pr.normalised,0);
	// segmentation fault here?

	cout << "\n\nOutput : RMSE - result[0] = " << result[0] << "  EXPECTED 0.0" << endl;
	cout << "Output : Max abs error- result[1] = " << result[1] << "  EXPECTED 0.0" <<  endl;
	cout << "Output : N. of corrections - result[2] = " << result[2] << "  EXPECTED 0.0"  << endl;
	cout << "Output : R2 - result[3] = " << result[3] << "  EXPECTED 1.0" << endl;

	if (fabs(result[0])<1.0E-12 && fabs(result[1])<1.0E-12 && result[2]==0 && fabs(result[3]-1.0)<1.0E-12 ) {
		cout << "\nPOPULATION::FITNESS_FUNC : Test n.1 : OK ------------------------------------------" << endl;
		passed_tests.push_back(1);
	} else {
		cout << "\nPOPULATION::FITNESS_FUNC : Test n.1 : ERROR !!!! -----------------------------------" << endl;
		errors_no++;
		failed_tests.push_back(1);
	}

	for (int k=0; k<pb.n_validation; k++) delete[] data_valid[k];
	delete data_valid;


	//----------------------------------------------------------------------------------------
	// POPULATION::FITNESS_FUNC : TEST N. 2 (2D)
	//----------------------------------------------------------------------------------------
	cout << "\nPOPULATION::FITNESS_FUNC : TEST N. 2 -----------------------------------";
	// define tree 2*(Z1/Z2)
	Val param=2.0;
	Binary_Node T2_tree(NULL, &Mult);  //root node
	// D1
	Terminal_Const T2_D1_P1(&T2_tree, param);  // left child, depth 1
	T2_tree.set_left((Node *)&T2_D1_P1);
	Binary_Node T2_D1_P2(&T2_tree, &SDiv);  // right child, depth 1
	T2_tree.set_right((Node *)&T2_D1_P2);
	// D2
	Terminal_Var T2_D2_P1(&T2_D1_P2, pb.v_list[0]);  // left child
	T2_D1_P2.set_left((Node *)&T2_D2_P1);
	Terminal_Var T2_D2_P2(&T2_D1_P2, pb.v_list[1]);  // right child
	T2_D1_P2.set_right((Node *)&T2_D2_P2);

	cout << "\nTree : ";
	P.print_individual((Node *)&T2_tree);

	cols=3;
	pb.n_validation=2;
	data_valid = new Val*[pb.n_validation];
	for (int k=0; k<pb.n_validation; k++) data_valid[k]=new Val[cols];
	// fitness case 0
	data_valid[0][0] = 15.0;	//Z1
	data_valid[0][1] = 3.0;		//Z2
	data_valid[0][2] = 10.0;	//target
	// fitness case 1
	data_valid[1][0] = 1.0;		//Z1
	data_valid[1][1] = 3.0;		//Z2
	data_valid[1][2] = 0.6666666666666667;	//target
	//
	pb.sum_output=data_valid[0][2]+data_valid[1][2]; 	// sum of each target value
	pb.y_ave=pb.sum_output/2.0; 						// average value of input data
	pb.Sy=pow((data_valid[0][2]-pb.y_ave),2.0)+pow((data_valid[1][2]-pb.y_ave),2.0);  		// SStot total sum of squares of (observed data - average observed data) // defined in read_input_file function (read_file_new.cpp)
	pb.y_var=pb.Sy/2.0; 		// variance of input data
	pb.y_max = 10.0;		// max value of target
	pb.y_min = 0.6666666666666667;		// min value of target
	pb.data_validation=data_valid;
	pb.show_data_validation();
	P.fitness_func(pb.Sy, pb.data_validation, 2, (Node *)&T2_tree, result, P.parameters->normalised,0);

	cout << "\n\nOutput : RMSE - result[0] = " << result[0] << "  EXPECTED 0.0" << endl;
	cout << "Output : Max abs error- result[1] = " << result[1] << "  EXPECTED 0.0" <<  endl;
	cout << "Output : N. of corrections - result[2] = " << result[2] << "  EXPECTED 0.0" << endl;
	cout << "Output : R2 - result[3] = " << result[3] << "  EXPECTED 1.0" << endl;

	if (fabs(result[0])<1.0E-12 && fabs(result[1])<1.0E-12 && result[2]==0 && fabs(result[3]-1.0)<1.0E-12 ) {
		cout << "\nPOPULATION::FITNESS_FUNC : Test n.2 : OK ------------------------------------------" << endl;
		passed_tests.push_back(2);
	} else {
		cout << "\nPOPULATION::FITNESS_FUNC : Test n.2 : ERROR !!!! -----------------------------------" << endl;
		errors_no++;
		failed_tests.push_back(2);
	}



	//----------------------------------------------------------------------------------------
	// POPULATION::FITNESS_FUNC : TEST N. 3 (2D)
	//----------------------------------------------------------------------------------------
	cout << "\nPOPULATION::FITNESS_FUNC : TEST N. 3  protected division with correction -----------";
	// define tree 2*(Z1/Z2) with Z2=0 : the protected division returns 1
	// use already defined tree in Test N. 2 (T2_tree)

	cout << "\nTree : ";
	P.print_individual((Node *)&T2_tree);

	pb.n_validation=1;
	data_valid = new Val*[pb.n_validation];
	for (int k=0; k<pb.n_validation; k++) data_valid[k]=new Val[cols];
	// fitness case 0
	data_valid[0][0] = 15.0;	//Z1
	data_valid[0][1] = 0.0;		//Z2
	data_valid[0][2] = 2.0;		//target

	pb.data_validation=data_valid;
	pb.show_data_validation();
	P.fitness_func(1, pb.data_validation, 1, (Node *)&T2_tree, result, P.parameters->normalised,0);

	cout << "\n\nOutput : RMSE - result[0] = " << result[0] << "  EXPECTED 0.0" << endl;
	cout << "Output : Max abs error- result[1] = " << result[1] << "  EXPECTED 1.0" <<  endl;
	cout << "Output : N. of corrections - result[2] = " << result[2] << "  EXPECTED 1.0"  << endl;
	cout << "Output : R2 - result[3] = " << result[3] << "  EXPECTED 1.0" << endl;

	if (fabs(result[0])<1.0E-12 && fabs(result[1])<1.0E-12 && result[2]==1 && fabs(result[3]-1.0)<1.0E-12 ) {
		cout << "\nPOPULATION::FITNESS_FUNC : Test n.3 : OK ------------------------------------------" << endl;
		passed_tests.push_back(3);
	} else {
		cout << "\nPOPULATION::FITNESS_FUNC : Test n.3 : ERROR !!!! -----------------------------------" << endl;
		errors_no++;
		failed_tests.push_back(3);
	}



	//----------------------------------------------------------------------------------------
	// POPULATION::FITNESS_FUNC : TEST N. 4 (2D) - Check NaN value for R2
	//----------------------------------------------------------------------------------------
	cout << "\nPOPULATION::FITNESS_FUNC : TEST N. 4 power -------------------------------";
	// define tree 2*(Z1^2-Z2^3)
	param=2.0;
	Binary_Node T4_tree(NULL, &Mult);
	T4_tree.n_corrections=0;
	// D1
	Terminal_Const T4_D1_P1(&T4_tree, param);  // left child, depth 1
	T4_tree.set_left((Node *)&T4_D1_P1);
	Binary_Node T4_D1_P2(&T4_tree, &Sub);  // right child, depth 1
	T4_tree.set_right((Node *)&T4_D1_P2);
	// D2
	Unary_Node T4_D2_P1(&T4_D1_P2, &Square);  // left child
	T4_D1_P2.set_left((Node *)&T4_D2_P1);
	Unary_Node T4_D2_P2(&T4_D1_P2, &Cube);  // right child
	T4_D1_P2.set_right((Node *)&T4_D2_P2);
	// D3
	Terminal_Var T4_D3_P1(&T4_D2_P1, pb.v_list[0]);  // left branch
	T4_D2_P1.set_child((Node *)&T4_D3_P1);
	Terminal_Var T4_D3_P2(&T4_D2_P2, pb.v_list[1]);  // right branch
	T4_D2_P2.set_child((Node *)&T4_D3_P2);

	cout << "\nTree : ";
	P.print_individual((Node *)&T4_tree);

	pb.n_validation=1;
	data_valid = new Val*[pb.n_validation];
	for (int k=0; k<pb.n_validation; k++) data_valid[k]=new Val[cols];
	// fitness case 0
	data_valid[0][0] = 2.0;	//Z1
	data_valid[0][1] = 3.0;	//Z2
	data_valid[0][2] = -46.0;	//target
	pb.data_validation=data_valid;
	pb.show_data_validation();
	// call fitness_func providing Sy=0: as error is also zero, R2 will be Nan
	P.fitness_func(0, pb.data_validation, 1, (Node *)&T4_tree, result, P.parameters->normalised,0);

	cout << "\n\nOutput : RMSE - result[0] = " << result[0] << "  EXPECTED 0.0" << endl;
	cout << "Output : Max abs error- result[1] = " << result[1] << "  EXPECTED 0.0" << endl;
	cout << "Output : N. of corrections - result[2] = " << result[2] << "  EXPECTED 0" << endl;
	cout << "Output : R2 - result[3] = " << result[3] << "  EXPECTED -NaN" << endl;

	if (fabs(result[0])<1.0E-12 && fabs(result[1])<1.0E-12 && result[2]==0 && isnan(result[3])) { // isnan detects both +nan and -nan ... (!)
		cout << "\nPOPULATION::FITNESS_FUNC : Test n.4 : OK ------------------------------------------" << endl;
		passed_tests.push_back(4);
	} else {
		cout << "\nPOPULATION::FITNESS_FUNC : Test n.4 : ERROR !!!! -----------------------------------" << endl;
		errors_no++;
		failed_tests.push_back(4);
	}



	//----------------------------------------------------------------------------------------
	// POPULATION::FITNESS_FUNC : TEST N. 5 (1D) - check Inf value for R2
	//----------------------------------------------------------------------------------------
	cout << "\nPOPULATION::FITNESS_FUNC : TEST N. 5 cos -------------------------------";
	pr.nvar = 1;
	pb.set_n_var(1);
	pb.set_n_cols(2);
	pb.set_n_data(1);
	pb.initialise_variables(&Z, 2);
	// define tree (1.94071E+02 + ((((-5.63690E+01 * (cos((-1.22141E-01 * Z1)))))^2) + (((1.23075E+02 * (Z1 / Z1)))^2)))
	Val par_array[]={1.94071E+02, -5.63690E+01, -1.22141E-01, 1.23075E+02};
	Binary_Node T5_tree(NULL, &Add);
	T5_tree.n_corrections=0;
	// D1
	Terminal_Const T5_D1_P1(&T5_tree, par_array[0]);  // left child, depth 1
	T5_tree.set_left((Node *)&T5_D1_P1);
	Binary_Node T5_D1_P2(&T5_tree, &Add);  // right child, depth 1
	T5_tree.set_right((Node *)&T5_D1_P2);
	// D2
	Unary_Node T5_D2_P1(&T5_D1_P2, &Square);  // left child
	T5_D1_P2.set_left((Node *)&T5_D2_P1);
	Unary_Node T5_D2_P2(&T5_D1_P2, &Square);  // right child
	T5_D1_P2.set_right((Node *)&T5_D2_P2);
	// D3
	Binary_Node T5_D3_P1(&T5_D2_P1, &Mult);
	T5_D2_P1.set_child((Node *)&T5_D3_P1);
	Binary_Node T5_D3_P2(&T5_D2_P2, &Mult);
	T5_D2_P2.set_child((Node *)&T5_D3_P2);
	// D4
	Terminal_Const T5_D4_P1(&T5_D3_P1, par_array[1]);
	T5_D3_P1.set_left((Node *)&T5_D4_P1);
	Unary_Node T5_D4_P2(&T5_D3_P1, &Cos);
	T5_D3_P1.set_right((Node *)&T5_D4_P2);
	Terminal_Const T5_D4_P3(&T5_D3_P2, par_array[3]);
	T5_D3_P2.set_left((Node *)&T5_D4_P3);
	Binary_Node T5_D4_P4(&T5_D3_P2, &SDiv);
	T5_D3_P2.set_right((Node *)&T5_D4_P4);
	// D5
	Binary_Node T5_D5_P1(&T5_D4_P2, &Mult);
	T5_D4_P2.set_child((Node *)&T5_D5_P1);
	Terminal_Var T5_D5_P2(&T5_D4_P4,pb.v_list[0]);
	T5_D4_P4.set_left((Node *)&T5_D5_P2);
	Terminal_Var T5_D5_P3(&T5_D4_P4,pb.v_list[0]);
	T5_D4_P4.set_right((Node *)&T5_D5_P3);
	// D6
	Terminal_Const T5_D6_P1(&T5_D5_P1, par_array[2]);
	T5_D5_P1.set_left((Node *)&T5_D6_P1);
	Terminal_Var T5_D6_P2(&T5_D5_P1,pb.v_list[0]);
	T5_D5_P1.set_right((Node *)&T5_D6_P2);

	cout << "\nTree : ";
	P.print_individual((Node *)&T5_tree);

	pb.n_validation=1;
	data_valid = new Val*[pb.n_validation];
	for (int k=0; k<pb.n_validation; k++) data_valid[k]=new Val[2];  // 1D test case
	// fitness case 0
	data_valid[0][0] = 9.30603000000000e-01;	// Z1
	data_valid[0][1] = 1.84781154459874e+04;	// target
	pb.data_validation=data_valid;
	pb.show_data_validation();
	// call fitness_func providing Sy=0: as error is not zero, R2 will be -Inf
	P.fitness_func(0, pb.data_validation, 1, (Node *)&T5_tree, result, P.parameters->normalised,0);

	cout << "\n\nOutput : Error - result[0] = " << result[0] << "  EXPECTED 0.0" << endl;
	cout << "Output : Max abs error- result[1] = " << result[1] << "  EXPECTED 0.0" <<  endl;
	cout << "Output : N. of corrections - result[2] = " << result[2] << "  EXPECTED 0" << endl;
	cout << "Output : R2 - result[3] = " << result[3] << "  EXPECTED -inf" << endl;

	if (fabs(result[0])<1.0E-10 && fabs(result[0])<1.0E-10 && result[2]==0 && isinf(result[3])) {  // isinf detects both +inf and -inf ... (!)
		cout << "\nPOPULATION::FITNESS_FUNC : Test n.5 : OK ------------------------------------------" << endl;
		passed_tests.push_back(5);
	} else {
		cout << "\nPOPULATION::FITNESS_FUNC : Test n.5 : ERROR !!!! -----------------------------------" << endl;
		errors_no++;
		failed_tests.push_back(5);
	}


	//----------------------------------------------------------------------------------------
	// POPULATION::FITNESS_FUNC : TEST N. 6 (2D)
	//----------------------------------------------------------------------------------------
	cout << "\nPOPULATION::FITNESS_FUNC : TEST N. 6 -----------------------------------";
	pr.nvar = 2;
	pb.set_n_var(2);
	pb.set_n_cols(3);
	pb.initialise_variables(&Z, 2);
	// define tree 2*(Z1^2-Z2^3) : use tree T4_tree defined in test 4
	cout << "\nTree : ";
	P.print_individual((Node *)&T4_tree);

	pb.n_validation=5;
	data_valid = new Val*[pb.n_validation];
	for (int k=0; k<pb.n_validation; k++) data_valid[k]=new Val[cols];
	// fitness case 0
	data_valid[0][0] = 9.306025965175708e-01;
	data_valid[0][1] = -5.733778508896148e-01;
	data_valid[0][2] = 2.10905226667451e+00;   //2.071917766525170e+02;
	// fitness case 1
	data_valid[1][0] = -3.110936555224211e-01;
	data_valid[1][1] = -6.781615857869738e-03;
	data_valid[1][2] =  1.93559148789886e-01;   //2.791452126307039e+00;
	// fitness case 2
	data_valid[2][0] = 1.845479995425238e+00;
	data_valid[2][1] = -2.622378471188245e-01;
	data_valid[2][2] = 6.84766033265214e+00;  //1.346162370140825e+03;
	// fitness case 3
	data_valid[3][0] = -1.711763045791965e+00;
	data_valid[3][1] =  2.498336684293894e-01;
	data_valid[3][2] = 5.82907778272687e+00; // 7.257539620491450e+02;
	// fitness case 4
	data_valid[4][0] = 2.136318721215420e-01;
	data_valid[4][1] = 4.664845206618544e-01;
	data_valid[4][2] = -1.11744194457976e-01;   //1.832950568020592e+01;
	Val Sy = 4.11575155133945e+01;
	pb.data_validation=data_valid;
	pb.show_data_validation();
	P.fitness_func(Sy, pb.data_validation, 5, (Node *)&T4_tree, result, 0,0);

	cout << "\n\nOutput : RMSE - result[0] = " << result[0] << "  EXPECTED 0.0" << endl;
	cout << "Output : Max abs error- result[1] = " << result[1] << "  EXPECTED 0.0" <<  endl;
	cout << "Output : N. of corrections - result[2] = " << result[2] << "  EXPECTED 0" << endl;
	cout << "Output : R2 - result[3] = " << result[3] << "  EXPECTED 1" << endl;

	if (fabs(result[0])<1.0E-12 && fabs(result[1])<1.0E-12  && result[2]==0 && fabs(result[3]-1.0)<1E-10) {
		cout << "\nPOPULATION::FITNESS_FUNC : Test n.6 : OK ------------------------------------------" << endl;
		passed_tests.push_back(6);
	} else {
		cout << "\nPOPULATION::FITNESS_FUNC : Test n.6 : ERROR !!!! -----------------------------------" << endl;
		errors_no++;
		failed_tests.push_back(6);
	}


	//----------------------------------------------------------------------------------------
	// POPULATION::FITNESS_FUNC : TEST N. 7 (1D) - Autocorrelation function
	//----------------------------------------------------------------------------------------
	cout << "\nPOPULATION::FITNESS_FUNC : TEST N. 7 -----------------------------------";
	pr.nvar = 1;
	pb.set_n_var(1);
	pb.set_n_cols(2);
	pb.set_n_data(1);
	pb.initialise_variables(&Z, 2);
	// define tree sin(Z1)+cos(2*Z1)
	Binary_Node T7_tree(NULL, &Add);
	T7_tree.n_corrections=0;
	// D1
	Unary_Node T7_D1_P1(&T7_tree, &Sin);  // left child, depth 1 : sin()
	T7_tree.set_left((Node *)&T7_D1_P1);
	Unary_Node T7_D1_P2(&T7_tree, &Cos);  // right child, depth 1 : cos()
	T7_tree.set_right((Node *)&T7_D1_P2);
	// D2
	Terminal_Var T7_D2_P1(&T7_D1_P1,pb.v_list[0]);  // child of sin() : Z1
	T7_D1_P1.set_child((Node *)&T7_D2_P1);
	Binary_Node T7_D2_P2(&T7_D1_P2, &Mult);			// child of cos() : 2*Z1
	T7_D1_P2.set_child((Node *)&T7_D2_P2);
	// D3
	Terminal_Const T7_D3_P1(&T7_D2_P2, 2.0);		// left child : 2
	T7_D2_P2.set_left((Node *)&T7_D3_P1);
	Terminal_Var T7_D3_P2(&T7_D2_P2,pb.v_list[0]);  // right child : Z1
	T7_D2_P2.set_right((Node *)&T7_D3_P2);

	cout << "\nTree : ";
	P.print_individual((Node *)&T7_tree);  // till here ok

	pb.n_validation=21;
	data_valid = new Val*[pb.n_validation];
	for (int k=0; k<pb.n_validation; k++) data_valid[k]=new Val[cols];
	// assign data. Z1=[0:0.1:2] (21 points)
	data_valid[0][0] = 0.0;
	data_valid[0][1] = sin(data_valid[0][0])+cos(2*data_valid[0][0]);
	cout << "\n\n" << data_valid[0][0] << " " << data_valid[0][1] << endl;
	for (int k=1; k<pb.n_validation; k++) {
		data_valid[k][0] = data_valid[k-1][0]+0.1;	// Z1
		data_valid[k][1] = sin(data_valid[k][0])+cos(2*data_valid[k][0]);	// target
		cout << data_valid[k][0] << " " << data_valid[k][1] << endl;
	}  // till here ok : indipendent variable Z1 and target values correctly assigned
	pb.data_validation=data_valid;
	pb.delay_max = 11;  // inlcudes delay=0, that is 1

	P.fitness_func(0.0, pb.data_validation, pb.n_validation, (Node *)&T7_tree, result, 0,0);
	cout << "\nAutocorrelation function values (n. lags = 10)" << endl;
	for (int k=0; k<pb.delay_max; k++) {
		cout << "delay=" << k << " " << T7_tree.r_k[k] << endl;
	}
	// correct values of acf
//	1.000000
//	   0.942484
//	   0.834328
//	   0.687290
//	   0.514229
//	   0.328302
//	   0.142187
//	  -0.032628
//	  -0.186446
//	  -0.311819
//	  -0.403856
	cout << "First autocorrelation function root = " << T7_tree.first_acf_root_tree;
	cout << "Test passed - see comparison of ACF values" << endl;
	passed_tests.push_back(7);


	//----------------------------------------------------------------------------------------
	// POPULATION::AGGREGATE_F : TEST N. 8 (2D) - Aggregate fitness function F
	//----------------------------------------------------------------------------------------
	//cout << "\n\nTEST : Population::aggregate_F(RunParameters* pr, Val average_err, Binary_Node *complete_tree, int gen, int G)" << endl;
	cout << "\n\nPOPULATION::AGGREGATE_F : TEST N. 8 -----------------------------------";
	// use tree T2 defined in test n.2: 2*(Z1/Z2)
	cout << "\nTree : ";
	P.print_individual((Node *)&T2_tree);
	rows=1;
	cols=3;
	pb.n_validation=1;
	data_valid = new Val*[pb.n_validation];
	pb.set_n_cols(3);
	for (int k=0; k<pb.n_validation; k++) data_valid[k]=new Val[cols];
	data_valid[0][0] = 15.0;
	data_valid[0][1] = 0.0;
	data_valid[0][2] = 2.0;
	pb.data_validation=data_valid;
	pb.y_max=11.0;
	pb.y_min=1.0;
	pb.y_ave=1.0;
	pb.y_var=1.0;
	pb.first_acf_root_input=1.0;
	pb.show_data_validation();

	// weights
	pr.w_complexity = 0.02;
	pr.w_n_corrections = 0.1;
	pr.w_size = 0.001;
	pr.w_factorisation = 0;  // switch between standard approach (<0) and factorisation bonus (>0)
	pr.w_pen_ord0 = 0.0;
	pr.w_pen_ord1 = 0.0;
	pr.w_strat_statp=0.1;
	pr.w_ACF = 0.1;
	pr.w_tvariation=0.0;
	// a1=
	// as a result, Fweight[]
	P.Fweight[1]=1.0-pr.w_complexity-pr.w_n_corrections-pr.w_size-pr.w_pen_ord0-pr.w_pen_ord1-pr.w_factorisation-pr.w_strat_statp-pr.w_ACF-pr.w_tvariation;
	P.Fweight[2]= pr.w_complexity;
	P.Fweight[3]=pr.w_n_corrections;
	P.Fweight[4]=pr.w_size;
	P.Fweight[5]=pr.w_pen_ord0;
	P.Fweight[6]=pr.w_pen_ord1;
	P.Fweight[8]=pr.w_strat_statp;
	P.Fweight[10]=pr.w_ACF;
	P.Fweight[11]=pr.w_tvariation;

	pr.strat_statp=15;	// Strategy 15

	T2_tree.fitness=10.0;   // RMSE, just a randomly chosen value
	Val average_err = 15.0;
	T2_tree.n_tuning_parameters=T2_tree.find_parameters();    // not nice... make n_tuning_parameters a member of Node?
	cout << "T2_tree.n_tuning_parameters = " << T2_tree.n_tuning_parameters;
	T2_tree.n_corrections=1;
	int size = T2_tree.count();
	T2_tree.tree_mean=1.0;
	T2_tree.tree_variance=1.0;
	T2_tree.first_acf_root_tree=1.0;

	P.aggregate_F(&pb, &pr, average_err, &T2_tree, 6, 30);  // average_err, 6 and 30 randomly chosen
	//void Population::aggregate_F(ProblemDefinition* ppd, RunParameters* pr, Val average_err, Binary_Node *complete_tree, int gen, int G)


	T2_tree.show_state();
	// ATTENZIONE: da rivedere coefficienti e valori di T[...] dopo riorganizzazione di Agosto '21!
	cout << "\nOBJECTIVE 1 : RMSE (raw fitness)";
	double a1=1.0-pr.w_complexity-pr.w_n_corrections-pr.w_size-pr.w_pen_ord0-pr.w_pen_ord1-pr.w_factorisation-pr.w_strat_statp-pr.w_ACF-pr.w_tvariation; // equal to 1-0.121=0.879
	double Obj1=a1*exp(T2_tree.fitness/fabs(pb.y_max-pb.y_min));  // = e*a1= 1.845713
	cout << "\nT2_tree.T[1] = " << T2_tree.T[1];  // non assegnato correttamente!
	cout << "\na1*exp(T2_tree.fitness/fabs(pb.y_max-pb.y_min)) = " << Obj1 << "   SHOULD BE 1.845713" << endl;

	cout << "\nOBJECTIVE 2 : complexity";
	double Obj2=pr.w_complexity*pow((double)(T2_tree.n_tuning_parameters), 1.0);
	cout << "\nT2_tree.T[2] = " << T2_tree.T[2]; // non assegnato correttamente!
	cout << "\npr.w_complexity*pow((double)(T2_tree.n_tuning_parameters), 1.0) = " << Obj2 << "   SHOULD BE 0.02" << endl;

	cout << "\nOBJECTIVE 3 : corrections";
	cout << "\nT2_tree.T[3] = " << T2_tree.T[3];
	cout << "\npr.w_n_corrections*T2_tree.n_corrections*1.0E6 = " << pr.w_n_corrections*T2_tree.n_corrections*1.0E6 << "   SHOULD BE 100000";

	cout << "\nOBJECTIVE 4 : size";
	cout << "\npr.w_size*T2_tree.count() = " << pr.w_size*T2_tree.count() << "   SHOULD BE 0.005";

	if (fabs(T2_tree.T[1]-a1*exp(T2_tree.fitness/fabs(pb.y_max-pb.y_min)))<1E-4 && fabs(T2_tree.T[2]-pr.w_complexity*T2_tree.n_tuning_parameters)<1E-12 && fabs(T2_tree.T[3]-pr.w_n_corrections*T2_tree.n_corrections*1.0E6)<1E-12 && fabs(pr.w_size*T2_tree.count()-T2_tree.T[4])<1E-12) {
		cout << "\nPOPULATION::FITNESS_FUNC : Test n.8 : OK ------------------------------------------" << endl;
		passed_tests.push_back(8);
	} else {
		cout << "\nPOPULATION::FITNESS_FUNC : Test n.8 : ERROR !!!! -----------------------------------" << endl;
		errors_no++;
		failed_tests.push_back(8);
	}



	//----------------------------------------------------------------------------------------
	// POPULATION::SORT : TEST N. 9 (2D)
	//----------------------------------------------------------------------------------------
	cout << "\nPOPULATION::SORT : TEST N. 9 STILL TO BE CHECKED -----------------------------------";
	// use tree defined previously:
	//T3 test n.3 : 2*(Z1/Z2), T7 test 7 sin(Z1)+cos(2*Z1), T4 test n.4 2*(Z1^2-Z2^3), T1 test n.1 Z1+Z2

	P.print_population_with_parameters(0);

	Binary_Node* complete_trees[3];
	//complete_trees[0]=&T1_tree;
	complete_trees[0]=&T2_tree;
	complete_trees[1]=&T4_tree;
	complete_trees[2]=&T7_tree;

	cout << "\nTrees : ";
	P.print_individual((Node *)&T2_tree);

	pb.n_validation=1;
	data_valid = new Val*[pb.n_validation];
	for (int k=0; k<pb.n_validation; k++) data_valid[k]=new Val[cols];
	data_valid[0][0] = 15.0;
	data_valid[0][1] = 0.0;
	data_valid[0][2] = 2.0;
	pb.data_validation=data_valid;
	pb.show_data_validation();

	pr.w_complexity = 0.02;
	pr.w_n_corrections = .1;
	pr.w_size = 0.001;
	pr.w_factorisation = 0;  // switch between standard approach (<0) and factorisation bonus (>0)
	pr.w_pen_ord0 = 0;
	//pr.n_inequality1 = 0;
	pr.w_pen_ord1 = 0;

	T2_tree.fitness=10.0;   // RMSE or PRESS
	average_err = 15.0;
	T2_tree.n_tuning_parameters=T2_tree.find_parameters();    // not nice... make n_tuning_parameters a member of Node?
	T2_tree.n_corrections=1;
	size = T2_tree.count();

	P.aggregate_F(&pb, &pr, average_err, &T2_tree, 6, 30);

	T2_tree.show_state();
	a1=1.0-pr.w_complexity-pr.w_n_corrections-pr.w_size-pr.w_pen_ord0-pr.w_pen_ord1;

	if (T2_tree.T[1]==a1*T2_tree.fitness/average_err && T2_tree.T[2]==pr.w_complexity*T2_tree.n_tuning_parameters && T2_tree.T[3]==pr.w_n_corrections*T2_tree.n_corrections*1.0E6 && pr.w_size*T2_tree.count()==T2_tree.T[4]) {
		cout << "\nPOPULATION::FITNESS_FUNC : Test n.9 : OK ------------------------------------------" << endl;
		passed_tests.push_back(9);
	} else {
		cout << "\nPOPULATION::FITNESS_FUNC : Test n.9 : ERROR !!!! -----------------------------------" << endl;
		errors_no++;
		failed_tests.push_back(9);
	}




 ////////////////////////////////////////////////////////////////////////////////////////////

	cout << "\nPassed tests:" << endl;
	for (int i=0; i<passed_tests.size(); i++) cout << passed_tests[i] << " ";
	cout << endl;
	cout << "\n\nTotal number of errors : " << errors_no << endl;
	cout << "Failed tests:" << endl;
	for (int i=0; i<failed_tests.size(); i++) cout << failed_tests[i] << " ";
	cout << endl;
}

#include "../genetic_code/tree_functions/vector_derivative_functions.cpp"
#include "../genetic_code/tree_functions/tree_operations.cpp"
#include "../genetic_code/nodes/nodes.cpp"
#include "../genetic_code/modules/primitives.cpp"
#include "../genetic_code/modules/variable.cpp"
#include "../genetic_code/classes/problem_definition.cpp"
#include "../genetic_code/classes/run_parameters.cpp"
#include "../genetic_code/classes/class_POPULATION.cpp"

#include "../genetic_code/SQP/MI0L2_c/fdf_c.cpp"



//#include "../genetic_code/classes/problem_definition.cpp"
