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


#include "./reporter.h"


// function to print evaluated points of the best individual (shape) to file "inputdata_stats.txt"
void Reporter::inputdatastats2file(ProblemDefinition *pb, const char* DIR_OUTPUT_char)
{
	string dir_output,file,r,s;
	const char *expr1;
	// string to char conversion
	dir_output=DIR_OUTPUT_char;
	file = "inputdata_stats.txt";
	r = "/";
	s =  dir_output+r+file;
	expr1 = s.c_str();
	Val target=0.0;
	Val tree=0.0;
	Val residual=0.0;


	//declare stream
	ofstream fout;  //NOT static, so every time (each new generation) a new "fout" is opened! This let the function overwrite the previous file (otherwise too big!)

	// open file: output stream
	fout.open(expr1, ios_base::out | ios_base::trunc);    //it doesn't append...

	static int w_c1 = 4;
	static int w_c = 14;
	// Features of the BUILDING data set
	//first line
	fout <<"# Input data (whole - see input file) statistics" << endl;
	fout << pb->get_n_data() << " # Number of record/fitness cases: n_data " << endl;
	fout << pb->get_n_var() << " # Number of input variables: n_var " << endl;
	fout << pb->sum_output << " # Sum of target values: sum_output  " << pb->Ny_omega_max << " # Ny_omega_max" << endl;
	fout << pb->y_ave << " # Average target value: y_ave " << endl;
	fout << pb->Sy << " # SStot Total sum of squares of (target- y_ave) : Sy " << endl;
	fout << pb->y_var << " # Target variance: y_var" << endl;
	fout << pb->y_max << " # Max value of target: y_max" << endl;
	fout << pb->y_min << " # Min value of target: y_min" << endl;
	fout << pb->first_acf_root_input << " # First autocorrelation function root" << endl;
	fout << pb->tot_variation_input << " # Total variation: tot_variation_input " << endl;
	fout << "# First 5 record/fitness cases with target= " << endl;
	for (int i=0; i<5; i++) {
		for (int j=0; j<pb->get_n_var()+1; j++) {
			fout << "# " << pb->get_data(i, j) << " ";
		}
		fout << endl;
	}
	// save autocorrelation values (each row corresponds to a delay value)
	fout << "# Autocorrelation function values (ProblemDefinition::compute_inputdata_stats()):" << endl;
	for (int i=0; i<pb->delay_max; i++) {
		fout << "# " << pb->r_k[i] << endl;
	}
	fout << "###" << endl;


	// Features of the TEST data set
	if (pb->data_test) {
		fout <<"# Test data statistics" << endl;
		fout << pb->n_test << " # Number of record/fitness cases: n_test " << endl;
		//fout << "Number of input variables: n_var = " << pb-> << endl; //n of input variables not checked in test data?
		fout << pb->sum_output_test << " # Sum of target values: sum_output_test " << endl;
		fout << pb->y_ave_test << " # Average target value: y_ave_test = " << endl;
		fout << pb->Sy_test << " # SStot Total sum of squares of (target- y_ave) : Sy_test" << endl;
		fout << pb->y_var_test << " # Target variance: y_var_test" << endl;
		fout << pb->y_test_max << " # Max value of target: y_max_test" << endl;
		fout << pb->y_test_min << " # Min value of target: y_min_test" << endl;
		fout << "# First 5 record/fitness cases with target:" << endl;
		for (int i=0; i<5; i++) {
			for (int j=0; j<pb->get_n_var()+1; j++) {
				fout << "# " << pb->data_test[i][j] << " ";
			}
			fout << endl;
		}
	} else {
		fout <<"# No Test data imported" << endl;
	}


	// close stream (at each call is open and then closed)
	fout.close();
}


// function to print statistical data to file "data_gp.txt"
void Reporter::stats2file(RunParameters *pr, Population *P, string DIR_OUTPUT, int gi, int check_end)
{ 

	string file,r,s;
	const char *expr1;
	// string to char conversion
	file = "data_gp.txt";
	r = "/";
	s =  DIR_OUTPUT+r+file;
	expr1 = s.c_str();
		
	ofstream fout;
	
	// width of the columns (c1 -> column 1), in characters
	static int w_c1 = 6; //6;
	static int w_c = 10; //14;

	// done only at the first time (generation 0)
	if (gi==0) {	
		// open file for writing, truncating
		fout.open(expr1, ios_base::out | ios_base::trunc);    //it doesn't append...
		//first line
		fout << "# Fitness, Size and Depth statistics" << endl;
		fout << "# M = " << pr->M << endl;
		fout << "# G = " << pr->G << endl;
		fout << "# Method = " << pr->method << endl;
		fout << left <<setw(w_c1) << "# Gen" ; 
		fout << left << setw(w_c) <<"Fit_min" << setw(w_c) << "Fit_ave" << setw(w_c) << "Fit_max" << setw(w_c) << "Fit_var";
		fout << left << setw(w_c) << "S_min" << setw(w_c) << "S_ave" << setw(w_c) << "S_max";
		fout << left << setw(w_c) << "D_min" << setw(w_c) << "D_ave" << setw(w_c) << "D_max";
		fout << left << setw(w_c) << "F_min" << setw(w_c) << "F_ave" << setw(w_c) << "F_max" << setw(w_c) << "F_var";
		fout << endl;
	} else {
		// open file for writing, appending
		fout.open(expr1, ios_base::out | ios_base::app);
	}

	//fitness data: print to file
	fout << gi << " " << scientific << P->Fit_min << " " << scientific << P->Fit_ave << " " <<  scientific << P->Fit_max  << " " <<  scientific << P->Fit_var;

	//size data: print to file
	fout << " " << scientific << P->S_min << " " << scientific << P->S_ave << " " << scientific << P->S_max;

	//depth statistics: print to file
	fout << " " << scientific << P->D_min << " " << scientific <<  P->D_ave << " " << scientific << P->D_max;
	
	//fitness data: print to file
	fout << " " << scientific << P->F_min << " " << scientific << P->F_ave << " " <<  scientific << P->F_max  << " " <<  scientific << P->F_var << endl;

	// close stream (at each call is open and then closed)
	fout.close();
}

// function to print evaluated points of the best individual (P->complete_trees[0]) to file "points_gp.txt"
void Reporter::points2file(RunParameters *pr, ProblemDefinition *pb, Population* P, string DIR_OUTPUT, int gi, int check_end, time_t start, time_t finish, double delta_t,int seed)
{

	string file, r,s;
	const char *expr1;
	// string to char conversion
	file = "points_gp.txt";
	r = "/";
	s =  DIR_OUTPUT+r+file;
	expr1 = s.c_str();
	Val target=0.0;
	Val tree=0.0;
	Val residual=0.0;

	
	//declare stream
	ofstream fout;  //NOT static, so every time (each new generation) a new "fout" is opened! This let the function overwrite the previous file (otherwise too big!)
	
	// open file: output stream
	fout.open(expr1, ios_base::out | ios_base::trunc);    //it doesn't append...
	

	static int w_c1 = 4;
	static int w_c = 14;
	//first line
	fout <<"# Training set and corresponding values of objective function (target) and best individual (tree) " << endl;
	fout << "M= " << pr->M << endl;
	fout << "G= " << pr->G << endl;
	fout << "Method= " << pr->method << endl;
	fout << "Generation= " << gi << endl; 
	fout << "Elapsed_time(sec)= " << delta_t << endl;
	fout << "Used_seed= " << seed << endl;
	fout << "Tree mean value on training data set= " << P->complete_trees[0]->tree_mean << endl;
	fout << "Tree variance on training data set= " << P->complete_trees[0]->tree_variance << endl;
	fout << "Tree min value:" << P->complete_trees[0]->tree_min <<endl;
	fout << "Tree max value:" << P->complete_trees[0]->tree_max <<endl;
	fout << "First autocorrelation function root:" << P->complete_trees[0]->first_acf_root_tree << endl;
	fout << "Tree total variation: " << P->complete_trees[0]->tot_variation_tree << endl;
	//fout << "###" << endl;
	char *expr;
	expr = P->complete_trees[0]->print();
	fout << "Tree expression:" << endl;
	fout << expr << endl;
	fout << left << setw(w_c1) << " " << setw(w_c) << "target" << setw(w_c) << "tree" << setw(w_c) << "residual (tree-target)" << endl;
	
	// this re-evaluation cycle can be saved if tree values are saved in Population::fitness_func() and here just recalled
	for (int i=0; i< pr->nfitcases; i++) {
		//print # of column
		fout << left << setw(w_c1) << i;
	
		//update the variables' value (using the whole data matrix)
		for (int j=0; j<pr->nvar; j++) {			// assign the right value to all the variables for the i-th fitness case
			(*(pb->v_list[j])).value = pb->get_data(i,j);
		}
		
		//print data relating to the first individual in tree (the best one): target, tree, residual=tree-target
		target=pb->get_data(i,pr->nvar);
		tree=P->complete_trees[0]->value(NULL);
		residual=tree-target;
		fout <<  " " << scientific << setprecision(10) << target << " " << scientific << tree << " " << scientific << residual << endl;
	}

	//fout << endl;

//	// save autocorrelation values (each row corresponds to a delay value)
//	// ATTENTION! All the attributes that are not saved in the tree after tuning (see ::evaluate)
//	// correspond to the last set of coefficient optimised, not necessarily to the best set and so to the best tree!!
//	fout << "Autocorrelation function values:" << endl;
//	for (int i=0; i<pb->delay_max; i++) {
//		fout << (P->complete_trees[0])->r_k[i] << endl;
//	}


	// close stream (at each call is open and then closed)
	fout.close();
}


// function that writes to file the best (complete) individual as per aggregate function F (not RMSE or R2) on training data set
void Reporter::update_best2file_build(Population *P, string DIR_OUTPUT, int gi, int check_end)
{ 
	string file, r,s;
	const char *expr1;
	char *expr; 
	// string to char conversion
	file = "best_gp.txt";
	r = "/";
	s =  DIR_OUTPUT+r+file;
	expr1 = s.c_str();

	ofstream fout;

	// done only at the first time (generation 0)
	if (gi==0) {	
		// open file for writing, truncating
		fout.open(expr1, ios_base::out | ios_base::trunc);
		fout << "# Gen" <<" F" << " Fitness(RMSE)" << " Corrections" << " R2(adim.)" << " Mean" << " Var" << " Min" << " Max" << " First_ACF" << " Tot_variation" << " Maxabserror";
		fout <<  " Expression"  << endl;
	} else {
		// open file for writing, appending
		fout.open(expr1, ios_base::out | ios_base::app);
	}

	// fetch the expression of the best individual (minimum F)
	expr = P->print(0, P->complete_trees);  //0 is the position in the array "complete_trees" of the individual with minimum F, not RMSE or R2!
	
	// print state variables and expression of the best individual (all these features refer to tree evaluated on building data set)
	fout << gi << " ";  // generation
	fout << scientific << P->complete_trees[0]->F;  // aggregate fitness F
	fout << " " << scientific << P->complete_trees[0]->fitness; // RMSE
	fout << " " << P->complete_trees[0]->n_corrections; // corrections
	fout << " " << scientific << P->complete_trees[0]->R2;  // R2
	fout << " " << scientific << P->complete_trees[0]->tree_mean;  // mean
	fout << " " << scientific << P->complete_trees[0]->tree_variance;  // variance
	fout << " " << scientific << P->complete_trees[0]->tree_min; // min
	fout << " " << scientific << P->complete_trees[0]->tree_max; // max
	fout << " " << scientific << P->complete_trees[0]->first_acf_root_tree; // first ACF root
	fout << " " << scientific << P->complete_trees[0]->tot_variation_tree; // total variation
	fout << " " << scientific << P->complete_trees[0]->maxabserror; // max absolute error on building data set
	fout << "  \"" << expr << "\"" << endl;  

	// free memory
	delete [] expr;

	// close stream (at each call is open and then closed)
	fout.close();
}


// function that finds and prints the best individual as per RMSE (!!!!) on the test data set
void Reporter::best2file_test(Population *P, string DIR_OUTPUT, int gi)
{
	int repr_tot = P->get_repr_tot();
	string file, r,s;
	const char *expr1;
	char *expr;
	// string to char conversion
	file = "best_gp_TEST.txt";
	r = "/";
	s =  DIR_OUTPUT+r+file;
	expr1 = s.c_str();

	// find the best individual as per RMSE on the test data set (it might not be the one in position 0...//expr = P->print(0, P->complete_trees);  //0 is the position in the array "complete_trees" of the individual with minimum F, not fitness!
	double fitness_test_max = 999999.0e+10; // too low? See MAX_VAL
	int i_best_test = -1;
	int i;
	// there must be at least one individual in the archive
	fitness_test_max=P->complete_trees[0]->fitness_test;   // ATTENTION: fitness is RMSE error (or PRESS), different from F
	i_best_test=0;
	// check if there is another individual in the archive that is better than the first one
	for (i=1; i<repr_tot; i++) {
		if (P->complete_trees[i]->fitness_test < fitness_test_max) { // && (P->complete_trees[i]->n_corrections_test==0)) {
			i_best_test = i;
			fitness_test_max = P->complete_trees[i]->fitness_test;
		}
	}

	// open file for writing, truncating
	ofstream fout;
	fout.open(expr1, ios_base::out | ios_base::trunc);
	// follow the same data pattern used to report results on the tuning/building data set
	fout << "#" << "Tree with lowest RMSE on test data set: i = " << i_best_test << endl;  //size of the "archive"
	fout << "# Gen" <<" F" << " Fitness(RMSE)" << " Corrections" << " R2(adim.)" << " Mean" << " Var" << " Min" << " Max" << " First_ACF" << " Tot_variation" << " Maxabserror";
	fout <<  " Expression"  << endl;

	// print state variables and expression of the best individual
	if (i_best_test>=0) { // safety check on value of i_best_test not to leave best_gp_TEST empty
		expr = P->print(i_best_test, P->complete_trees);
		// print state variables and expression of the best individual (all these features refer to tree evaluated on building data set)
		fout << gi << " ";  // generation
		fout << scientific << "N/A";  // aggregate fitness F on test data set not computed yet.. is it useful?
		fout << " " << scientific << P->complete_trees[i_best_test]->fitness_test; // RMSE
		fout << " " << P->complete_trees[i_best_test]->n_corrections_test; // corrections
		fout << " " << scientific << P->complete_trees[i_best_test]->R2_test;  // R2
		fout << " " << scientific << P->complete_trees[i_best_test]->tree_mean_test;  // mean
		fout << " " << scientific << P->complete_trees[i_best_test]->tree_variance_test;  // variance
		fout << " " << scientific << P->complete_trees[i_best_test]->tree_min_test; // min
		fout << " " << scientific << P->complete_trees[i_best_test]->tree_max_test; // max
		fout << " " << scientific << P->complete_trees[i_best_test]->first_acf_root_tree_test; // first ACF root
		fout << " " << scientific << P->complete_trees[i_best_test]->tot_variation_tree_test; // total variation
		fout << " " << scientific << P->complete_trees[i_best_test]->maxabserror_test; // max absolute error on building data set
		fout << "  \"" << expr << "\"" << endl;
	}
	else {
		fout << "Reporter::best2file_test not able to select the best individual evaluated on the test data set : repr_tot=" << repr_tot << endl;
	}

	// free memory
	if (expr!=NULL) delete [] expr;

	// close stream (at each call is open and then closed)
	fout.close();
}



// function that prints to file all the individuals in the archive (first repr_tot) together
// with their scores (F, RMSE or fitness, R2, Hits and expression) evaluated on the
// TRAINING (BUILDING) data set
void Reporter::archive2file_build(Population *P, string DIR_OUTPUT, int gi, int check_end)
{ 
	string file, r,s;
	const char *expr1;
	char *expr; 
	// string to char conversion
	file = "latest_archive.txt";
	r = "/";
	s =  DIR_OUTPUT+r+file;
	expr1 = s.c_str();
	int repr_tot = P->get_repr_tot();

	// open file: output stream (only the last generation objectives are recorded!)
	ofstream fout;				//NOT static, so every time a new "fout" is opened! This let the function overwrite the previous file.
	fout.open(expr1, ios_base::out | ios_base::trunc);    

	// print header (references to the objectives that have to be compared globally)
	fout << "#" << " R " << repr_tot << endl;  //size of the "archive" 
	fout << "# " <<  "Last generation " << " F " << " Fitness(RMSE) " << " R2(adim.)" << " Hits " << " Expression " << endl;

	// print to file all the repr_tot individuals in the archive, evaluated on the training/building data set
	for   (int i=0; i<repr_tot; i++) { //the first repr_tot individuals are the best-so-far. Check: They must have previously been sorted!
		expr = P->print(i, P->complete_trees);  

		fout << gi << " ";  // generation
		fout << scientific << P->complete_trees[i]->F;  // aggregate fitness F
		fout << " " << scientific << P->complete_trees[i]->fitness; // RMSE
		fout << " " << P->complete_trees[i]->n_corrections; // corrections
		fout << " " << scientific << P->complete_trees[i]->R2;  // R2
		fout << " " << scientific << P->complete_trees[i]->tree_mean;  // mean
		fout << " " << scientific << P->complete_trees[i]->tree_variance;  // variance
		fout << " " << scientific << P->complete_trees[i]->tree_min; // min
		fout << " " << scientific << P->complete_trees[i]->tree_max; // max
		fout << " " << scientific << P->complete_trees[i]->first_acf_root_tree; // first ACF root
		fout << " " << scientific << P->complete_trees[i]->tot_variation_tree; // total variation

		fout << "  \"" << expr << "\"" << endl;
		
		// free memory
		delete [] expr;	
	}


	// close stream (at each call is open and then closed)
	fout.close();

}


// function that evaluate on the TEST (VALIDATION) data set all the individuals
// in the archive (first repr_tot) together and print them with their score (only RMSE or fitness)
// to file
void Reporter::archive2file_test(Population *P, string DIR_OUTPUT, int gi)
{

	int repr_tot = P->get_repr_tot();

	// print to file
	string file, r,s;
	const char *expr1;
	char *expr;
	// string to char conversion
	file = "final_archive_TEST.txt";
	r = "/";
	s =  DIR_OUTPUT+r+file;
	expr1 = s.c_str();

	// open file: output stream (only the last generation objectives are recorded!)
	ofstream fout;				//NOT static, so every time a new "fout" is opened! This let the function overwrite the previous file.
	fout.open(expr1, ios_base::out | ios_base::trunc);

	// print header (references to the objectives that have to be compared globally)
	fout << "#" << "Size of the archive : repr_tot = " << repr_tot << endl;  //size of the "archive"
	fout << "# Gen" <<" F" << " Fitness(RMSE)" << " Corrections" << " R2(adim.)" << " Mean" << " Var" << " Min" << " Max" << " First_ACF" << " Tot_variation" << " Maxabserror";
	fout <<  " Expression"  << endl;

	// print to file all the repr_tot individuals in the archive, evaluated on the test data set
	for (int i=0; i<repr_tot; i++) {
		expr = P->print(i, P->complete_trees);
		// print state variables and expression of the best individual (all these features refer to tree evaluated on building data set)
		fout << gi << " ";  // generation
		fout << scientific << "N/A";  // aggregate fitness F on test data set not computed yet.. is it useful?
		fout << " " << scientific << P->complete_trees[i]->fitness_test; // RMSE
		fout << " " << P->complete_trees[i]->n_corrections_test; // corrections
		fout << " " << scientific << P->complete_trees[i]->R2_test;  // R2
		fout << " " << scientific << P->complete_trees[i]->tree_mean_test;  // mean
		fout << " " << scientific << P->complete_trees[i]->tree_variance_test;  // variance
		fout << " " << scientific << P->complete_trees[i]->tree_min_test; // min
		fout << " " << scientific << P->complete_trees[i]->tree_max_test; // max
		fout << " " << scientific << P->complete_trees[i]->first_acf_root_tree_test; // first ACF root
		fout << " " << scientific << P->complete_trees[i]->tot_variation_tree_test; // total variation
		fout << " " << scientific << P->complete_trees[i]->maxabserror_test; // max absolute error on building data set
		fout << "  \"" << expr << "\"" << endl;

		// free memory
		delete [] expr;
	}

	// close stream
	fout.close();

}


void Reporter::node_stats2file(RunParameters *pr, Population *P, string DIR_OUTPUT)
// function that prints data relative to node selection frequencies to file
// this function is called at the end of the execution
// (so data is not written if the execution is interrupted...)
{
	string file, r,s;
	const char *expr1;
	// string to char conversion
	file = "node_selection.txt";
	r = "/";
	s =  DIR_OUTPUT+r+file;
	expr1 = s.c_str();
	// open the output stream
	ofstream fout;
	fout.open(expr1, ios_base::out | ios_base::trunc);    //it doesn't append...
	
	// write, in the order: total_nodes_selected selected_nodes_per_type[1 ... 4] selected_nodes_per_depth[0 .... d_lim]
	fout << "# Node selection statistics" << endl;
	fout << "# Total n. selected   |   N. Selected per type: Binary = " << NODE_BINARY << ", Unary = " << NODE_UNARY;
	fout << ", Var = " << NODE_VAR << ", Const = " << NODE_CONST << "      |       N. selected per depth [0 1 2   ... " << pr->depth_lim << " ]" << endl;
	fout << P->total_nodes_selected << "   ";
	for (int k=0; k<4; k++)
		fout << P->selected_nodes_per_type[k] << "  ";
	fout << "  ";
	for (int k=0; k<pr->depth_lim; k++)
		fout << P->selected_nodes_per_depth[k] << "  ";
	fout << P->selected_nodes_per_depth[pr->depth_lim] << endl;

	// close stream (at each call is open and then closed)
	fout.close(); 
}


// function to print to file the number of evaluations at each generation
void Reporter::n_tree_eval2file(Population *P, string DIR_OUTPUT, int gi, int check_end)
{
	string file, r,s;
	const char *expr1;
	// string to char conversion
	file = "n_tree_evaluations.txt";
	r = "/";
	s =  DIR_OUTPUT+r+file;
	expr1 = s.c_str();

	ofstream fout;

	// done only at the first time (generation 0)
	if (gi==0) {
		// open file for writing, truncating
		fout.open(expr1, ios_base::out | ios_base::trunc);
		fout << "# Gen    No. of tree evaluations"  << endl;
	} else {
		// open file for writing, appending
		fout.open(expr1, ios_base::out | ios_base::app);
	}

	// print to file generation number and no. of tree evaluations
	fout << gi << " " << P->get_n_tree_evaluations() << endl;

	// close stream (at each call is open and then closed)
	fout.close();
}


// function that writes to file all data related to adaptive approach for genetic operators
// (eps_neutral, learning_window, the constructive, destructive and
// neutral genetic operations rates, so that they can be plotted (write script in matlab).
void Reporter::adaptive_gen_ops_data2file(Population *P, string DIR_OUTPUT, int gi, int check_end)
{
	string file, r,s;
	const char *expr1;
	// string to char conversion
	file = "adaptation_data.txt";
	r = "/";
	s =  DIR_OUTPUT+r+file;
	expr1 = s.c_str();

	ofstream fout;

	// done only at the first time (generation 0)
	if (gi==0) {
		// open file for writing, truncating
		fout.open(expr1, ios_base::out | ios_base::trunc);
		fout << "# Parameters related to adaptive approach"  << endl;
		fout << "# eps_neutral = "  << P->eps_neutral << endl;
		fout << "# learning_window = "  << P->learning_window << endl;
		fout << "# Gen  Pop_size repr_rate cross_rate  mut_rate repr_tot reproduction_perf repr_av_delta cross_tot crossover_perf cross_av_delta tot_smut s_mutation_perf smut_av_delta tot_pmut p_mutation_perf pmut_av_delta"  << endl;

	} else {
		// open file for writing, appending
		fout.open(expr1, ios_base::out | ios_base::app);
	}

	// print to file generation no. and population size
	fout << gi << " " << P->get_size() <<  " ";
	// print repr_rate, cross_rate, mut_rate
	fout << scientific << P->get_repr_rate() << " " << P->get_cross_rate() << " " << P->get_mut_rate() <<  " ";
	// print reproduction performance
	fout << P->tot_repr << " ";
	for (int k=0; k<3; k++)
		fout << P->reproduction_perf[k] << " ";
	for (int k=0; k<3; k++)
			fout << P->repr_av_delta[k] << " ";
	// print crossover performance
	fout << P->tot_cross << " ";
	for (int k=0; k<3; k++)
		fout << P->crossover_perf[k] << " ";
	for (int k=0; k<3; k++)
		fout << P->cross_av_delta[k] << " ";
	// print subtree mutation performance
	fout << P->tot_smut << " ";
	for (int k=0; k<3; k++)
		fout << P->s_mutation_perf[k] << " ";
	for (int k=0; k<3; k++)
		fout << P->smut_av_delta[k] << " ";
	// print point mutation performance
	fout << P->tot_pmut << " ";
	for (int k=0; k<3; k++)
		fout << P->p_mutation_perf[k] << " ";
	for (int k=0; k<3; k++)
			fout << P->pmut_av_delta[k] << " ";
	// end line
	fout << endl;
	// close the stream to file if termination criterion met or at the last generation
	//if  ((gi== pr.G) || (check_end))

	// close stream (at each call is open and then closed)
	fout.close();

}


// function that prints coefficients and corresponding objective contributions to file
void Reporter::F_coefficients2file(Population *P, string DIR_OUTPUT, int gi, int check_end)
{
	string file, r,s;
	const char *expr1;
	// string to char conversion
	file = "F_coefficients.txt";
	r = "/";
	s =  DIR_OUTPUT+r+file;
	expr1 = s.c_str();

	ofstream fout;

	// done only at the first time (generation 0)
	if (gi==0) {
		// open file for writing, truncating
		fout.open(expr1, ios_base::out | ios_base::trunc);
		fout << "# Values of F weight coefficients throughout the evolution"  << endl;
		fout << "# Fweight[0]	sum of all weights (must be 1)" << endl;
		fout << "# Fweight[1]  	RMSE or PRESS value" << endl;
		fout << "# Fweight[2]  	complexity - no of tuning parameters - W_COMPLEXITY" << endl;
		fout << "# Fweight[3] 	no of corrections - W_N_CORRECTIONS" << endl;
		fout << "# Fweight[4] 	no of nodes - tree size - W_SIZE" << endl;
		fout << "# Fweight[5]	penalisation from inequality constraints order 0 - W_PEN_ORD0" << endl;
		fout << "# Fweight[6] 	penalisation from inequality constraints order 1 - W_PEN_ORD1" << endl;
		fout << "# Fweight[7] 	penalisation to increase factorisation (depth of first division) - W_FACTORISATION" << endl;
		fout << "# Fweight[8] 	penalisation on statistical properties of the tree - W_STRAT_STATP " << endl;
		fout << "# Fweight[9] 	penalisation of diverging trees (high level polynomials are present)" << endl;
		fout << "# Fweight[10]  penalisation linked to autocorrelation function (point at which ACF halves) - W_ACF" << endl;
		fout << "# Fweight[11]  penalisation linked to total variation - W_TVARIATION" << endl;
		fout << "# Gen Fweight[i] Fc_perc_ave[i] for i=1,...,11 Fweight[0]"  << endl;

	} else {
		// open file for writing, appending
		fout.open(expr1, ios_base::out | ios_base::app);
	}

	// print to file generation no. and population size
	fout << gi << "  " << fixed;
	for (int i=1; i<12; i++) {
		// save F weight value and average objective error in the archive
		fout << P->Fweight[i] << " " << P->Fc_perc_ave[i] << "  ";
	}

	fout << "  " << P->Fweight[0] << endl;  //sum of all weights, it must be 1

	// close the stream to file if termination criterion met or at the last generation
	//if  ((gi== pr.G) || (check_end))

	// close stream (at each call is open and then closed)
	fout.close();

}
