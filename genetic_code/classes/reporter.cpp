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


// function to print evaluated points of the best individual (shape) to file "points_gp.txt"
void Reporter::inputdatastats2file(ProblemDefinition *pb, string DIR_OUTPUT)
{
	string file, r,s;
	const char *expr1;
	// string to char conversion
	file = "inputdata_stats.txt";
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
	fout <<"# Input data (whole - see input file) statistics" << endl;
	fout << "Number of record/fitness cases: n_data = " << pb->get_n_data() << endl;
	fout << "Number of input variables: n_var = " << pb->get_n_var() << endl;
	fout << "Sum of target values: sum_output = "<< pb->sum_output << endl;
	fout << "Average target value: y_ave = " << pb->y_ave << endl;
	fout << "SStot Total sum of squares of (target- y_ave) : Sy = " << pb->Sy << endl;
	fout << "Target variance: y_var = " << pb->y_var << endl;
	fout << "Max value of target: y_max = " << pb->y_max << endl; // usa basic_stat_analysis.cpp
	fout << "Min value of target: y_min = " << pb->y_min << endl; // usa basic_stat_analysis.cpp

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

// function to print evaluated points of the best individual (shape) to file "points_gp.txt"
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
	fout << "Tree average value on training data set= " << P->complete_trees[0]->tree_mean << endl;
	fout << "Tree variance on training data set= " << P->complete_trees[0]->tree_variance << endl;
	
	fout << left <<setw(w_c1) << " " << setw(w_c) << "target"<< setw(w_c) << "tree" << setw(w_c) << "residual (tree-target)" << endl;
	
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
		fout <<  " " << scientific  << target << " " << scientific << tree << " " << scientific << residual << endl;
	}

	//fout << endl;
	// close stream (at each call is open and then closed)
	fout.close();
}


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
		fout << "# Gen" <<" F" << " Fitness(RMSE)" << " R2(adim.)" << " Hits " <<  "Expression"  << endl;
	} else {
		// open file for writing, appending
		fout.open(expr1, ios_base::out | ios_base::app);
	}

	// fetch the expression of the best individual
	expr = P->print(0, P->complete_trees);  //0 is the position in the array "complete_trees" of the individual with minimum F, not fitness!
	
	// print state variables and expression of the best individual
	fout << gi << " " << scientific << P->complete_trees[0]->F;
	fout << " " << scientific << P->complete_trees[0]->fitness;
	fout << " " << scientific << P->complete_trees[0]->R2;
	fout << " " << P->complete_trees[0]->hits;
	fout << "  \"" << expr << "\"" << endl;  

	// free memory
	delete [] expr;

	// close stream (at each call is open and then closed)
	fout.close();
}


// function that finds and prints the best individual on the test data set
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

	// find the best individual on the test data set
	double f_max = 999999.0e+10; // too low? See MAX_VAL
	int i_best_test = -1;
	int i;
	// there must be at least one individual in the archive
	f_max=P->complete_trees[0]->fitness_test;   // ATTENTION: with fitness is meant "error" (RMSE or PRESS), different from F
	i_best_test=0;
	// check if there is another individual in the archive that is better than the first one
	for (i=0; i<repr_tot-1; i++) {
		if (P->complete_trees[i]->fitness_test <= f_max) { // && (P->complete_trees[i]->n_corrections_test==0)) {
			i_best_test = i;
			f_max = P->complete_trees[i]->fitness_test;
		}
	}
	// safety check on value of i_best_test not to leave best_gp_TEST empty?


	ofstream fout;

	// open file for writing, truncating
	fout.open(expr1, ios_base::out | ios_base::trunc);
	// follow the same data pattern used to report results on the tuning/building data set
	// # gen  F    RMSE   R2(adim)  Hits  Expr
	fout << "# " <<  "Final gen " << "F " << "RMSE " << "R2(adim.) " << "Hits " << "Expression " << endl;


	// fetch the expression of the best individual
	//expr = P->print(0, P->complete_trees);  //0 is the position in the array "complete_trees" of the individual with minimum F, not fitness!

	// print state variables and expression of the best individual
	if (i_best_test>=0) {
		expr = P->print(i_best_test, P->complete_trees);
		fout << gi;  // Generation
		fout << " " << scientific << "not_implemented_yet_Reporter::best2file_test";// P->complete_trees[i_best_test]->F_test;   // F... SHOULD BE F_TEST!!!! Not implemented yet
		fout << " " << scientific << P->complete_trees[i_best_test]->fitness_test;  // ATTENTION: with fitness is meant "error" (RMSE or PRESS), different from F
		fout << " " << scientific << P->complete_trees[i_best_test]->R2_test; // R2
		fout << " "  << P->complete_trees[i_best_test]->hits_test; // Hits
		fout << "  \"" << expr << "\"" << endl;  // Expression
	}
	else {
		fout << "Reporter::best2file_test not able to select the best individual evaluated on the test data set : repr_tot=" << repr_tot << endl;
	}

	// free memory
	if (expr!=NULL) delete [] expr;

	// close stream (at each call is open and then closed)
	fout.close();
}
//



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
		fout << gi << " " << scientific << P->complete_trees[i]->F;
		fout << " " << scientific << P->complete_trees[i]->fitness;
		fout << " " << scientific << P->complete_trees[i]->R2;
		fout << " "  << P->complete_trees[i]->hits;
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
	fout << "# " <<  "Final generation " << " Fitness(RMSE) " << " R2(adim.)" << "Hits" << " Expression " << endl;

	// print to file all the repr_tot individuals in the archive, evaluated on the training/building data set
	for (int i=0; i<repr_tot; i++) {
		expr = P->print(i, P->complete_trees);
		fout << gi << " " << scientific << P->complete_trees[i]->fitness_test;
		fout << " " << scientific << P->complete_trees[i]->R2_test;
		fout << " "  << P->complete_trees[i]->hits_test;
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


// function that writes to file all data related to adaptive approach
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
