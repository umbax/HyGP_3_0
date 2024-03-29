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

#ifndef CLASS_POPULATION_H_
#define CLASS_POPULATION_H_


#include <iostream>  // basic i/o commands: cout, cin, scientific, fixed
#include <string>    // to manipulate strings (C)
#include <cstring>   // to manipulate strings (C)
#include <cmath>
#include <stdio.h>
#include <math.h>       /* isinf, sqrt */
#include <cstdlib>   // NULL, exit, EXIT_FAILURE
#include <algorithm> // max_element

#include "../nodes/nodes_types.h"
#include "../nodes/nodes.h"
#include "./run_parameters.h"
#include "./problem_definition.h"
#include "../tree_functions/vector_derivative_functions.h"
#include "../tree_functions/tree_operations.h"


extern "C" {extern void opti_(int*, double*, int*, int*, int*);}  //for Andrey's method - TI0L2 and MI0L2 - IT WORKS PERFECTLY
extern "C" {int fdf_c__ (double *,  double *, int *,int *); }

// the following two lines are a test on for a different optimiser (TINL2 and MINL2 - Umberto's method)
//extern "C" {extern void optitinl2_(int*, double*, int*, int*,int*);}  //for Umberto's method - TINL2 and MINL2 - WORKS? WORK STILL NEEDED...
//extern "C" {int  fdf_ctinl2__ (double*,  double*, double*, int*, int*); }



class Population {
 private:
	Binary_Node **trees;        // list of trees without parameters at the current generation
    int size;                   // size of the population (number of individuals, which is constant throughout the generations)
	Binary_Node *p_new_subtree;  // pointer to the new subtree created each time subtree mutation is called

	// TREE SIZE AND DEPTH
	int depth;				//current depth of one node in the initial tree
	int depth_max;		//higher than 0 (root node is binary)
	int full;					//	1 or 0 (for RAMPED method only) 
	int ntotf;					// quantity of trees built with FULL method (for RAMPED method only) 
	int n_full;				// number of trees built with FULL method so far (for RAMPED method only) 
	int n_grow;				// number of trees built with GROW method so far (for RAMPED method only)
	
	//GENETIC OPERATIONS' PARAMETERS
	double repr_rate;
	double cross_rate;
	double mut_rate;

	int comp_tot;	// no of archive individuals that have been copied. Updated each generation
	int new_tot;	// no of archive individuals that have been generated from scratch to fill void left by purged copies. Updated each generation
	
	int n_test_cases;   //number of fitness cases
	int n_test_cases_tune; // number of fitness cases used for tuning
 	int n_test_cases_fitness; // number of fitness cases used for actual fitness evaluation

	// computational cost of the run
 	int n_tree_evaluations;


	// functions to generate random subtrees from either binary or unary node
    // - in the future replace with "tree_generators"
	Node *new_node(Node *);
    void build_tree(Binary_Node *);
    void build_tree(Unary_Node *);
	//functions to implement the FULL, GROW method to generate the initial population
	int FULL_method(void);
	int GROW_method(void);


	// function to evaluate the depth of a node in a tree
	int evaluate_root_distance(Node *p_node);

	//function to insert parameters (constants) into a single individual
	void parameters_allocation(Binary_Node*,  Binary_Node**);
	
	// functions that implement selection methods
	int tournament(int, int, int);
	
	// functions that implement genetic operations
	int select_node(Binary_Node*);    // you can choose simple method (biased), or fair method (each depth has the same probability to be sampled)
	int potential_depth(Node *, Node *);
	void crossover(Binary_Node *,Binary_Node *,int, int, Binary_Node **,Binary_Node **, int*);
	int point_mutation(Binary_Node *, int);
	Binary_Node *generate_subtree(int);

	Val MAX_F_VALUE=1E+99;

  public:
  	Binary_Node **complete_trees;   //list of trees with parameters at the current generation	
  	
    // Population constructor 
  	//Population(RunParameters, ProblemDefinition);
	Population(RunParameters*, ProblemDefinition*);
	
	// Population destructor					
    ~Population(void);        
	
	ProblemDefinition* problem;
	RunParameters* parameters;
	
	// function to evaluate each individual (formerly a private function, turned to public 4/6/2017)
	// input: address to matrix of data to be used for evaluation, no of records, address of a tree (root node), address of array of result)
	void fitness_func(Val, Val**, int, Node*, Val*, bool, bool);
	// parameters and function to evaluate aggregate function F (see variable in Binary Node)
	void aggregate_F(ProblemDefinition*, RunParameters*, Val, Binary_Node*, int, int);

	int get_size(void);    // returns the population size
	double get_repr_rate(void);  // returns repr_rate
	double get_cross_rate(void);  // returns cross_rate
	double get_mut_rate(void);   // returns mut_rate
	double compute_time(time_t, time_t, double*);  		// function to compute the time needed to compute each generation
	void print_individual(Node *);   // print to screen a single individual (tree) given its address
	void print_population_with_parameters(int);  // print the entire population of trees with parameters (complete_trees[])
	void print_population_without_parameters(int);  // print the entire population of trees without parameters (trees[])
	void sort(int,int (*)(const void *, const void *));            // sort population by fitness
    void kill_and_fill(ProblemDefinition*);                //function to kill comp_rate+new_rate percentage of the population and to fill the vacuum with new individuals

    int composition[4];		// contains the number of individuals copied, freshly introduced, generated by crossover and mutation at each generation
    void new_spawn(RunParameters, int, int); // 2nd argument ProblemDefinition, 			//implemented by me: performs REPRODUCTION, CROSSOVER, MUTATION!
	void population_reproduction(Binary_Node **, Binary_Node **, int , bool, int*);  // performs reproduction with filter for copies at population level
	void population_crossover(Binary_Node **, Binary_Node **, int, int*);
	void population_mutation(Binary_Node **, Binary_Node **, int, int, int*);

	int terminate (double); // sets the number of hits on each member (added 16/10/08)
    char *print(int,Binary_Node**);             // returns the expression of a single tree
	
	// "getter" methods (return private members)
    int get_n_tree_evaluations(void) {return n_tree_evaluations;};
    int members(void) {return size;}; // returns the current size of population
    double fitness(int i) {return trees[i]->fitness;};
	// returns the number of hits of a member
    int no_hits(int i) {return trees[i]->hits;};
	// returns the number of nodes of a member
	int no_nodes(int i) {return trees[i]->count();};
	// returns the depth of a member
	int tree_depth(int i) {return trees[i]->calc_depth();};
	int get_repr_tot() {return (composition[0]+composition[1]);}; //returns the number of individuals copied and newly generated during reproduction
	

	// function to bypass unary node (used for parameter insertion) - NOT USED CURRENTLY
	Node *remove_unary_node(Unary_Node *u_node);
	
	void population_parameters_allocation(void);

	int ntree_fdf_c; // SQP parameter optimisation: simple way to make fdf_c able to get the values of the tree whose costants are being optimized

	void evaluate(int, int); // sets error, F and other attributes  of each tree
	void evaluate_complete_trees_on_test_dataset(void); // evaluates fitness, hits, corrections, R2 of complete trees on test data set
	void perform_editing(int i); // edits trees[i], individual without parameters
	int tuning_individual(int, Binary_Node *, Binary_Node *, int);      //finds, extracts the terminal_const variables and optimise them through SQP
	void tuning_individual_RMSE_single_guess(Binary_Node*, double*, int*, int*, Val*); // performs RMSE error evaluation on a single individual for a given set of parameters (tree coefficients)
	void tuning_individual_PRESS_single_guess(Binary_Node*, double*, int*, int*, Val*);  // performs PRESS error evaluation (crossvalidation) on a single individual for a given set of parameters (tree coefficients)
	double constraint_evaluation(Val**, int, char*, Variable **, int, Binary_Node *);

	// pulsation control
	void find_pulsations(Binary_Node *);     //function that find pulsations updating n_pulsations and index_puls 

	//function searching for first division in a tree
	void search_first_op(Binary_Node*, Node *, int);
	// function that searches for high level polynomials in complete trees (22/11/20)
	void search_high_level_polynomials(Binary_Node*, Node *);
	// function that purges diverging terms in trees without parameters (29/11/20)
	void purge_diverging_terms(Binary_Node*, Node *);


	// function which returns the value of a given tree (indicated by ntree, index of trees[])
	Val tree_value (Binary_Node*, int* );
	// function to update the values of the parameters in the SQP routine (called recursively by fdf_)
	void update_complete_tree(Binary_Node*, double *, int);	

	// function to compute fitness value and derivatives respect to the parameters to be tuned (see SQP - TINL2 and MINL2)
	double jacobian_ij (double, int, int);
	// function that computes the tree derivative in a given point along a direction specified by a second point (see matrix fun_der)
	void get_tree_derivative_given_points(Val **,Binary_Node*, double **, int, double *);
	// function that computes the tree derivative in a given point along a direction declared by a normalised vector
	double get_tree_derivative_given_norm_vector(ProblemDefinition, Binary_Node*);

	
	int identical_nodes(Node*, Node*);
	int identical_trees(Binary_Node*, Binary_Node*);   //function to check if two trees are identical
	void prune_tree(Binary_Node*);   //function that deletes from the tree the least important branch (radical Lamarck approach...)	

	// ARCHIVE STORING THE BEST TREE
	//vector <Node *> p_par_best;  //array containing the pointers to terminal const nodes in the best tree of the previous generation
	//void inherit_parameters(Binary_Node *, Binary_Node *); //function that fetch and store the parameters of the best individual
	Binary_Node* best_tree;
	Binary_Node* best_complete_tree;
	void update_ext_archive(void);
		
	// RANDOM VALUE GENERATOR
	Val constant_generation(Val, Val, Val);
	Val constant_generation(Val, Val, Val, Val, Val);

	// Aggregate F (fitness) weights and components (objectives)
	double F[12]; // entry [0] stores the weighted sum of components F[i] (not elegant.. as F is the sum of T[i])
	double Fweight[12]; // values of weights/coefficients of aggregate F - entry [0] has no useful purpose...
	void adjust_Fweight(int);
//	double F1, a1;  //used to store the main objective, RMSE or PRESS value
//	double F2, a2;  // second objective, related to complexity (no of tuning parameters) - W_COMPLEXITY
//	double F3, a3; 	// third objective (no of corrections) - W_N_CORRECTIONS
//	double F4, a4; 	// fourth objective (no of nodes - tree size) - W_SIZE
//	double F5, a5;	//penalisation from inequality constraints order 0
//	double F6, a6; 	//penalisation from inequality constraints order 1
//	double F7, a7; 	// penalisation to increase factorisation (depth of first division)
//	double F8, a8; 	// ADDED 11/8/20: penalisation on statistical properties of the tree (average and variance)
//	double F9, a9; 	// ADDED 22/11/20: penalisation of diverging trees (high level polynomials - that is not argument of any function - are present)
//	double F10, a10; // ADDED 23/1/21: penalisation linked to autocorrelation function (point at which ACF halves) - see ::fitness_func()
//	double F11, a11; // ADDED 6/2/21: penalisation linked to total variation

	// adaptive genetic operations rates
	// adaptive approach settings
	double sum_observed_var;
	double observed_var_thr;  //threshold for the observed variable used in adaptive approach
	double eps_neutral; // HYPERPARAMETER: fitness increase under which the genetic operations is considered neutral
	int learning_window; // length in number of generations of the window over which the neutral, constructive or destructive rates are computed
	int learning_on; // flag for learning phase
	double learning_var;
	double sum_learning_var;
	int window_counter;
	// counters of total individuals generated by genetic operations referring to current generation
	int tot_repr;
	int tot_cross;
	int tot_smut;
	int tot_pmut;
	// performance variables referring to current generation
	int reproduction_perf[3]; // 0: destructive, 1: neutral, 2: constructive
	int crossover_perf[3]; // 0: destructive, 1: neutral, 2: constructive
	int s_mutation_perf[3];  // 0: destructive, 1: neutral, 2: constructive
	int p_mutation_perf[3];  // 0: destructive, 1: neutral, 2: constructive
	double repr_av_delta[3]; // 0: destructive, 1: neutral, 2: constructive
	double cross_av_delta[3]; // 0: destructive, 1: neutral, 2: constructive
	double smut_av_delta[3]; // 0: destructive, 1: neutral, 2: constructive
	double pmut_av_delta[3]; // 0: destructive, 1: neutral, 2: constructive
	// performance variables referring to current generation
	int window_reproduction_perf[3]; // destructive, 1: neutral, 2: constructive
	int window_crossover_perf[3]; // destructive, 1: neutral, 2: constructive
	int window_s_mutation_perf[3];  // destructive, 1: neutral, 2: constructive
	int window_p_mutation_perf[3];  // destructive, 1: neutral, 2: constructive
	void measure_genetic_op_performance(Binary_Node*);
	void adapt_genetic_operators_rates(void);
	void adapt_genetic_operators_rates_notused(void); //actually not used in the execution...for tests
	void adapt_genetic_operators_rates_notused2(void); //actually not used in the execution...for tests

	// STATISTICS
	// ARCHIVE STATISTICS in the current generation
	Val Fit_min;
	Val Fit_ave;
	Val Fit_max;
	Val Fit_var;
	double S_min;
	double S_ave;
	double S_max;
	double D_min;
	double D_ave;
	double D_max;
	double F_min;
	double F_ave;
	double F_max;
	double F_var;
	double pen_ord1_ave;
	double T_perc_ave[12];	// array storing the average value of each single aggregate F component in the elite mulitplied by corresp. weight (used to implement the adaptive weight change - ::adjust_Fweight)
	double Fc_perc_ave[12]; // archive-averaged percentual contribution of each objective with respect to F (see Population::compute_statistics()) - Fc_perc_ave[0] not used
	double Fc_ave[12];    // average value of each single aggregate F component in the elite (see Population::compute_statistics()) - Fc_ave[0] not used
	void compute_statistics(void);


	// NODE SELECTION STATISTICS 
	int total_nodes_selected;
	int* selected_nodes_per_depth; 
	int selected_nodes_per_type[4];
	void compute_selected_nodes_statistics(Binary_Node* , int);

	// to allow use of an externally defined tree
	Binary_Node* external_tree;
	void assign_external_tree(Binary_Node* tree) {external_tree = tree;};
	Val ext_result[4];

	// HyPSO optimiser
	double pso_objfunction(double*, int, Binary_Node*); // HyPSO parameter optimisation
	double** zeros(int, int);
	double* zeros(int);
	int regenerate_on_boundary(double**, int, double**);
	int regenerate_from_swarm_best(double**, int, int, double**, double**);
	int psominimize(Binary_Node*, int, double*, double*, int, int, double**, double*);
	void hypso_launcher(RunParameters*, Binary_Node*, int, double*);


};


#endif /* CLASS_POPULATION_H_ */
