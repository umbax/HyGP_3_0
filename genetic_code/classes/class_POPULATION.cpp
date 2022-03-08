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


using namespace std;

#include "./class_POPULATION.h"


// class POPULATION - function definitions

// POPULATION CLASS CONSTRUCTOR
// takes the RunParameters class and the problem definition class
//Population::Population(RunParameters pr, ProblemDefinition pb)
Population::Population(RunParameters* pr, ProblemDefinition* pb)
{
	int COMMENT=1;

	// data in RunParameters and ProblemDefinition
    problem = pb;
    parameters = pr;

    // tree generation - shouldn't be general members, but local parameters (generation methods)
	ntotf = (int)(pr->p_FULL*pr->M);				// --- only for RAMPED method (see parameters)
	n_full = 0;
	n_grow = 0;
	
	// STATE variables - they may change by the generation
	size = pr->M;
	repr_rate = pr->repr_rate; // repr_rate;
	cross_rate = pr->cross_rate;
	mut_rate = pr->mut_rate;
	// number of individuals to be generated using reproduction, crossover, mutations
	composition[0] = (int)(floor(repr_rate*size));   //reproduction
	composition[1] = 0;								// new
	composition[2] = (int)(floor(cross_rate*size));	// crossover
	composition[3] = (int)(floor(mut_rate*size)); 	// mutation
	comp_tot = (int)(floor(pr->comp_rate*size));
	new_tot = (int)(floor(pr->new_rate*size));
	
	// adaptive approach settings
	sum_observed_var = .0;
	observed_var_thr = .0;
	eps_neutral = 1.0e-5;
	learning_window = 0; //updated in Population::adapt_genetic_operators
	learning_on = 1; // flag for learning phase : 1 if learning on, 0 if off
	learning_var = .0;
	sum_learning_var = .0;
	window_counter = 0;
	// number of individuals ACTUALLY generated using reproduction, crossover, mutations
	tot_repr= 0;
	tot_cross= 0;
	tot_smut= 0;
	tot_pmut= 0;
	for (int i=0; i<3; i++) {  // initialisation of genetic operations performance counters
		// for current generation
		reproduction_perf[i]=0;
		crossover_perf[i]=0;
		s_mutation_perf[i]=0;
		p_mutation_perf[i]=0;
		// for window
		window_reproduction_perf[i]=0;
		window_crossover_perf[i]=0;
		window_s_mutation_perf[i]=0;
		window_p_mutation_perf[i]=0;
	}

	// initial setting: tuning (building) and evaluation are done on the same dataset
	n_test_cases = pr->nfitcases;
	n_test_cases_tune = pr->nfitcases;
	n_test_cases_fitness = pr->nfitcases;
	
	ntree_fdf_c = -1;   //-1 is a fake value just to start

	n_tree_evaluations = 0;

	// aggregate fitness function (::aggregate_F) weights and components initialisation
	for (int i=0; i<12; i++) 	{
		F[i]=0.0; // entry [0] stores the weighted sum of components F[i]
		Fc_perc_ave[i]=0.0;
		Fc_ave[i]=0.0;
		T_perc_ave[i]=0.0;
	}
	Fc_perc_ave[0]=0.0;
	Fc_ave[0]=0.0;
	T_perc_ave[0]=0.0;
	Fweight[2]=pr->w_complexity;		// a2 -> number of tuning parameters (see strategy)
	Fweight[3]=pr->w_n_corrections;	// a3 -> number of corrections performed by protected operations
	Fweight[4]=pr->w_size;			// a4 -> model size (number of nodes)
	Fweight[5]=pr->w_pen_ord0;		// a5 -> penalisation of unsatisfied inequality constraint, order 0 (value)
	Fweight[6]=pr->w_pen_ord1;		// a6 -> penalisation of unsatisfied inequality constraint, order 1 (first derivative)
	Fweight[7]=pr->w_factorisation;   // a7 -> penalisation for lack of factorisation (depth of first division)
	Fweight[8]=pr->w_strat_statp; 	// a8 -> difference in mean and variance (see corresp. strategy) between original signal and evolved model
	Fweight[9]=0.0; //0.4; //0.3; //0.2; // a9 -> presence of high level polynomials that cause divergent behaviour
	Fweight[10]=pr->w_ACF; //1.0E-1;  //23/1/21 test  // a10 -> difference in point at which ACF halves for original signal and evolved model
	Fweight[11]=pr->w_tvariation; //0.0; //1.0E-1; //1.0E-1; // a11 -> difference in total variation between original signal and evolved model
	// the weight of the primary objective, RMSE error, is the residual to 1 of the sum of previous coefficients a2 to a8
	Fweight[1]=double(1.-Fweight[2]-Fweight[3]-Fweight[4]-Fweight[5]-Fweight[6]-Fweight[7]-Fweight[8]-Fweight[9]-Fweight[10]-Fweight[11]); //-a7); // The sum of all a_i coefficients must be 1!!
	Fweight[0]=0.0;
	for (int i=1; i<12; i++) 	{
		Fweight[0]=Fweight[0]+Fweight[i]; // entry [0] stores the sum of all weights Fweight[i]
	};


	// initialise population statistics data (with values that make no harm)
	Fit_min = 1.0;
	Fit_ave = 1.0;
	Fit_max = 1.0;
	Fit_var = 1.0;
	S_min = 1.0;
	S_ave = 1.0;
	S_max = 1.0;
	D_min = 1.0;
	D_ave = 1.0;
	D_max = 1.0;
	pen_ord1_ave = 1.0;
	
	external_tree=NULL;


	// initialise variables for node selection statistics -  to be put in RunStatistics!
	total_nodes_selected = 0;
	selected_nodes_per_depth = new int[parameters->depth_lim+1];
	if (selected_nodes_per_depth == NULL) {
		cerr << "\nError : not enough memory to create : selected_nodes_per_depth";
		exit(-1);
	}
	for (int k=0; k<parameters->depth_lim+1; k++) {
			selected_nodes_per_depth[k] = 0;
			//cout << "Selected_nodes_per_depth [" << k << "] = " <<  selected_nodes_per_depth[k] << endl;
	}
	for (int k=0; k<4; k++) { 
			selected_nodes_per_type[k] = 0;
			//cout << "Selected_nodes_per_type [" << k << "] = " <<  selected_nodes_per_type[k] << endl;
	}


	// external super archive - initialisation with NULL value
	best_tree = NULL;
	best_complete_tree = NULL;

    // initialise the pointer to new subtree for subtree mutation
	p_new_subtree = NULL;

	// allocate pointer to complete_trees array (containing the trees with parameters).
	// MIND: complete_trees are not generated! This is done in Population::evaluate
	// the root node is assumed to be a binary node
	complete_trees = new Binary_Node *[size];    //accessible from everywhere
	 if (!complete_trees) {
        // output error message
        cerr << "Population::Population : ERROR ! Can't allocate complete_trees pointer list!\n";
        exit(-1);
    }
	for (int k =0; k<size; k++)
		complete_trees[k] = NULL;

	// allocate the tree pointer list
	// the root node is assumed to be a binary node
    trees = new Binary_Node *[size];
    if (!trees) {
        // output error message
        cerr << "Population::Population : ERROR ! Can't allocate tree pointer list!\n";
        exit(-1);
    }
    for (int k =0; k<size; k++)
    		trees[k] = NULL;


    // go through the list and create new top level (root) nodes, which by definition have no parent (parent=NULL).
    // Currently root nodes can be binary functions only.
    // In the future allow any kind of node to be root node
    for (int i=0;i<size;i++) {
		//cout << "\n Population:Population CHECK: problem->num_b_funcs = " << problem->num_b_funcs << endl;
    	trees[i] = new Binary_Node(NULL,problem->b_func_list[int_rand(problem->num_b_funcs)]);    //int_rand(2)+2]);    // 2]); //gives multiplication!
    }

    // now go through the nodes and recursively build random parameterless trees
	if ((parameters->method==2) || (parameters->method==3)) // method for generating the initial pop. -> 1=no limits, 2=FULL, 3=GROW, 4=RAMPED
		depth_max = parameters->depth_max;

	for (int i=0;i<size;i++) {
		if (parameters->method==4) {						//just for RAMPED method
			depth_max = int_rand(parameters->depth_max - parameters->depth_min+1)+parameters->depth_min;
			full = int_rand(2);   //if 1 FULL, if 0 GROW
			if  (full) {
				if (n_full == ntotf)     
					full = 0;
				}
			else {
				if (n_grow ==(size-ntotf))
					full = 1;
			}
		}

		depth = 1;
		// generate the parameterless tree
		build_tree(trees[i]);

		if (parameters->method==4) {
			if (full)
				n_full++;
			else
				n_grow++;
		}
	};

	//depending on the method:
	switch (parameters->method) {
	
		case 1:	// no limits (original method)
			{
			printf ("\nPopulation::Population : INITIALIZATION with NO LIMITS method");
			}
			break;


		case 2:  // FULL method
			{
			printf ("\nPopulation::Population : INITIALIZATION with FULL method");
			}
			break;

		case 3:	//GROW method
			{
			printf ("\nPopulation::Population : INITIALIZATION with GROW method");
			}
			break;



		case 4:	//RAMPED method
			{
			cout << "\nPopulation::Population : INITIALIZATION with RAMPED method";
			cout << "\nPopulation size = " << size;
			cout << "\nNo. of trees expected to be built with FULL method: ntotf=p_FULL*M= " << ntotf;
			cout << "\nNo. of trees effectively built with FULL method: n_full = " << n_full;
			cout << "\nNo. of trees effectively built with GROW method: n_grow = " << n_grow;
			cout << "\nTotal No. of trees created: " << n_full + n_grow << endl;
			//check for errors for RAMPED method
			if ((n_full + n_grow)!=size) {
				cerr << "\nERROR!!! RAMPED method : the numbers of tree created is not equal to the size of the population!";
				exit(-1);
			}
			}
			break;
	}

	// here trees[] has been initialised and "contains" parameterless trees, complete_trees[] instead is an array of NULL values

	cout << "Population::Population : INITIALIZATION terminated. OK." << endl;
}



// Population destructor
Population::~Population(void)
{
	// delete
	delete[] selected_nodes_per_depth;

	// delete content
	for (int j=0;j<size;j++) {
		delete trees[j];
		//if (complete_trees[j]) delete complete_trees[j];
	}
	// delete containers
	delete[] trees;
	//if (complete_trees) delete[] complete_trees;
	
}



// this generates a random node, whose parent is n. it is mutually recursive
// with the build_tree's
// INPUT : new node's parent address
// OUTPUT : pointer to the newly created node 
Node *Population::new_node(Node *n)
{   
    static int c = 0;
    
    // pick a node type for new node
    // (MH - probably want to replace this with something that takes current
    // tree depth into account or terminal nodes probability to be chosen 
	// to limit the complexity of the initial trees)
	int type =-1;
	
	// just a try to force the use of a constant in case spow is the parent's function
	if (&Spow)   //check if the primitive is used... to avoid run-time errors 
		if (n->type == NODE_BINARY)
			if  (!strcmp( ((Binary_Node*)n)->get_func()->sign,"^") )
				if (problem->num_u_funcs)
					type = int_rand(4);
				else
					type = 2+int_rand(2);
	

	if (type < 0) {	
		switch (parameters->method) {

			case 1:
			{
				// no limits (original method)
				if (problem->num_u_funcs)
					type = int_rand(3);
				else
					type = 2*int_rand(2);
			}
			break;
	
			case 2:  // FULL method
			{
				type = FULL_method();
			}
			break;
	
			case 3:	//GROW method
			{
				type = GROW_method();
			}
			break;

			case 4:	//RAMPED method
			{
				if (full) {
					type = FULL_method();
				}
				else {
					type = GROW_method();	
				}
			}
			
			break;
		}
	}

    // pointer to the new node
    Node *newnode;

    // create the new node
    switch (type) {

      case NODE_BINARY:
        {
            // create it with a random function
			Binary_Node *nn= new Binary_Node(n,problem->b_func_list[int_rand(problem->num_b_funcs)]);   // int_rand(2)+2]); //<- solo mult and div          //
            // save the pointer
            newnode = (Node *)nn;
            
            // recursively build the new subtree
          depth++;	
			build_tree(nn);
        }
        break;

      case NODE_UNARY:
        {
            // create it with a random function
            Unary_Node *nn = new Unary_Node(n,problem->u_func_list[int_rand(problem->num_u_funcs)]);
            // save the pointer
            newnode = (Node *)nn;
            
            // recursively build the new subtree
			depth++;	
			build_tree(nn);
        }
        break;

      case NODE_CONST:
        {
            // create it with a random constant value
            Terminal_Const *nn = new Terminal_Const(n, (Val)1.0);   //const_rand());
            // save the pointer
            newnode = (Node *)nn;
			
        }
        break;

      case NODE_VAR:
        {
            // create it with a random variable
            Terminal_Var *nn = new Terminal_Var(n,problem->v_list[int_rand(parameters->nvar)]);
            // save the pointer
            newnode = (Node *)nn;
			
        }
        break;
    }

    // return this new node
	return newnode;
}





// function that generates a random parameterless subtree from a binary node base.
void Population::build_tree(Binary_Node *n)
{
	// make the left subtree
    Node *left = new_node((Node *)n);
	// assign it to this node's left
	n->set_left(left);

    // make the right subtree
	Node *right = new_node((Node *)n);
    // assign it to this node's right
    n->set_right(right);
	depth--;
}





// function that generates a random parameterless subtree from a unary node base
void Population::build_tree(Unary_Node *n)
{
    // make the child subtree
    Node *child = new_node((Node *)n);
    // assign it
    n->set_child(child);
	depth--;
}



int Population::get_size(void)
{
	return size;
}

// returns repr_rate
double Population::get_repr_rate(void)
{
	return repr_rate;
}

// returns cross_rate
double Population::get_cross_rate(void)
{
	return cross_rate;
}

// returns mut_rate
double Population::get_mut_rate(void)
{
	return mut_rate;
}


// selection of crossover point
// input: the number of the tree in the trees[] array
// output: the number of the node
int Population::select_node(Binary_Node* p_tree)//(int p1)
{
	int COMMENT = 0;  //1 comments, 0 silent...
	
	int n_nodes = 0; 
	int depth = 0;
	int selected_depth = 0;
	int selected_node = 0;
	vector <int> depth_nodes;  
	Node* p_node;	
	if (COMMENT) cout << "\n\nPopulation::select_node" << endl;
	//--------------------------------------------------------------------------------------------------------------
	// simple method (biased) : select a node randomly among all
	// selected_node = int_rand((p_tree->count()-1)) + 2;  //this exclude the root
 	// selected_node = int_rand(p_tree->count()) + 1;   //this includes the root, but gives problems due to the exchange of the root node during crossover
	//-----------------------------------------------------------------------------------------------------------

	// fair method: each depth has the same probability to be sampled
	// retrieve the depth of the current tree
	depth = p_tree->calc_depth();

	// randomly select a depth smaller or equal to the depth of the current tree
	//selected_depth = int_rand(depth)+1;   //root node (depth 0) is so excluded...
	selected_depth = int_rand(depth+1);  // root node (depth 0) is included
	if (COMMENT) cout << "Selected depth : " << selected_depth;

	// find the pointers to the nodes sharing the same depth (as the number is unknown a priori, use a vector...)
	if (selected_depth == 0)
		// if depth 0 is selected, there's only root node...
		depth_nodes.push_back(1);
	else {
		// selected depth different from 0...
		n_nodes = p_tree->count();
		for (int i=1; i<n_nodes+1; i++) {		//searching through the nodes, following the numeric order
			if (COMMENT) cout << "\nChecking node n. " << i << " depth...";
			p_node = p_tree->find(i);
			if (evaluate_root_distance(p_node) == selected_depth) {		
				depth_nodes.push_back(i);						// appends the pointer to terminal_const to the vector
			}
		}
	}	

	if (COMMENT) {
		cout << "\nPopulation::select_node" << endl;
		cout << "selected_depth = " << selected_depth << endl;
		cout << "depth_nodes.size() = " << (int)(depth_nodes.size()) << endl;
		for (int i = 0; i<(int)(depth_nodes.size()); i++) {
			cout << "depth_nodes[" << i << "] = " << (int)depth_nodes[i] << endl;
		}
	}

	// now the vector depth_nodes stores all the nodes having depth selected_depth.
	// choose randomly a node among them
	int k = int_rand((int)(depth_nodes.size()));
	selected_node = depth_nodes[k];
	
	if (COMMENT)
		cout << "\nPopulation::select_node :  selected node = " <<  selected_node << endl;

	// update node selection statistics
	compute_selected_nodes_statistics(p_tree, selected_node);

	// return selected node
	return selected_node;  
}


// function that returns the depth of the current node from the root
// (the function calc_depth begins the count from the leaves... )
// input: pointer to the node
// output: depth of the node (distance from root in number of steps)
int Population::evaluate_root_distance(Node *p_node)
{
	int COMMENT = 0; //1 comments, 0 silent

	int root_distance;
	Node *current_parent;
	if (COMMENT) {
		cout << "\nPopulation::evaluate_root_distance :  ";
		cout << " address of the node = " << p_node;
	}

	root_distance = 0;
	current_parent = p_node->get_parent();
	while (current_parent != NULL) {
		root_distance++;
		current_parent = current_parent->get_parent();
	}
	
	if (COMMENT) cout << "\nEvaluate_root_distance : root distance = " << root_distance;
	//return the distance of the given node from the root
	return root_distance;
}






//function that verify the potential depth of a tree after a subtree insertion (see subtree mutation or crossover)
//input: pointer to the root nodes of the subtree to be replaced, pointer to the subtree to be introduced
// output: depth of the single branch affected by the virtual substitution
int Population::potential_depth(Node *p_old_subtree, Node *p_new_subtree)
{
	int COMMENT = 0; //1 comments, 0 silent

	int root_distance;
	Node *current_parent;
	if (COMMENT)
		cout << "\nPopulation::potential_depth :  ";
	
	// compute root distance
	root_distance = evaluate_root_distance(p_old_subtree);

	// return new branch depth. Note that as the individual before the virtual substitution
	// was within the depth limit, the only way it can exceed the limit is that the new branch 
	// is exactly the branch responsible for breaking it. So the number returned would be the maximum 
	// depth among the branches, that is the depth of the tree...
	return (p_new_subtree->calc_depth() + root_distance);
}


// function to perform crossover
// input:
// Binary_Node *p1	pointer to parent 1
// Binary_Node *p2	pointer to parent 2
// int cp1			crossover point in parent 1
// int cp2			crossover point in parent 2
// Binary_Node **c1	pointer to pointer to child 1 (double ** because redeclaration is done)
//	Binary_Node **c2	pointer to pointer to child 2 (double ** because redeclaration is done)
//	int* cross_point		pointer to array containing information on crossover points (just to check)
// output:
// none
void Population::crossover (Binary_Node *p1, Binary_Node *p2, int cp1, int cp2, Binary_Node **c1,
               Binary_Node **c2, int* cross_point)
{
	int COMMENT = 0; //1 comments, 0 silent

	char *expr;

	static int n_times = 1;
	if (COMMENT) {
		cout << "\nEntered Population::crossover --- " << n_times << " times" << endl;
		cout << "p1 = " << p1<< ", p2 = " << p2 << endl;
		cout << "cp1 = " << cp1 << ", cp2 = " << cp2 << endl;
		cout << "c1 = " << c1 << ", c2 = " << c2 << endl;
	} 
	// copy the first parent
	if (COMMENT) 
		printf("\nPopulation::crossover : copy of the 1st parent...");
	*c1 = (Binary_Node *)tree_copy(p1,NULL);
	if (COMMENT) 
		printf("done");

    // copy the second
	if (COMMENT) 
		printf("\nPopulation::crossover : copy of the 2nd parent...");
    *c2 = (Binary_Node *)tree_copy(p2,NULL);
	if (COMMENT) 
		printf("done");
    

    // get the crossover sections from the parents (right now this will not pick the root node)
	// crossover point in the first parent
   cross_point[0] = cp1;
	// get the address cs1 of the first crossover point (node cp1) 
    Node *cs1 = p1->find(cp1);
	if (COMMENT) 
		printf("\nPopulation::crossover : 1st parent crossover point = %i ",cp1);
	

   	// crossover point in the second parent
   cross_point[1] = cp2;
	//  get the address cs2 of the second crossover point (node cp2) 
    Node *cs2 = p2->find(cp2);
	if (COMMENT) 
		printf("\nPopulation::crossover : 2nd parent crossover point = %i ",cp2);
   
	// ---------------------------------------------------------------------------------------------------------------------------
    // delete the selected subtree rooted in cp1 from the 1st child and copy in the 2nd parent's subtree rooted in cp2

    // find the node that will be deleted from the 1st child
    Node *delnode = (*c1)->find(cp1);
    // save its parent
    Node *parent = delnode->get_parent();
    // if the parent is a binary node, figure out if this was a left or
    // right child
    int left=0;
	if (parent)
    	if (parent->type == NODE_BINARY)
      		if (((Binary_Node *)parent)->get_left() == delnode)
        		left = 1;
      
    // delete the subtree
	delete delnode;
	if (!parent)
		*c1 = NULL; // to avoid problems...

    // make a copy of the subtree rooted in cp2 from parent 2
    Node *newnode1 = tree_copy(cs2,parent);

    // reset the parent's pointer
    if (parent) {
    	// if it was a unary pointer, its easy
    	if (parent->type == NODE_UNARY)
      		((Unary_Node *)parent)->set_child(newnode1);
    	else {
		 	// it must be a binary node then
      		if (left)
        		((Binary_Node *)parent)->set_left(newnode1);
      		else
        		((Binary_Node *)parent)->set_right(newnode1);  
		} 
	}
	else {
		*c1 =  (Binary_Node*)newnode1;  
	}
	// -------------------------------------------------------------------------------------------
	if (COMMENT) {
		cout << "\nChild 1 done" << endl;
		cout << "*c1 = " << *c1 << endl;
		cout << "parent = " << parent << " , newnode1 = " << newnode1 << endl;
	}

	// -------------------------------------------------------------------------------------------
    // now delete the cs from the 2nd child and link in a copy of the 1st
    // parent's cs  (begin cut-and-paste :)
    //

    // find the node that will be deleted from the 2nd child
    delnode = (*c2)->find(cp2);
    // save its parent
    parent = delnode->get_parent();
    // if the parent is a binary node, figure out if this was a left or
    // right child
    left=0;
	if (parent)
    	if (parent->type == NODE_BINARY)
      		if (((Binary_Node *)parent)->get_left() == delnode)
        		left = 1;
      
    // delete the subtree
    delete delnode;
	if (!parent)
		*c2 = NULL; // to avoid problems...

    // make a copy of the subtree rooted in cp1 from parent 1
    Node *newnode2 = tree_copy(cs1,parent);

    // reset the parent's pointer
    if (parent) {
    	// if it was a unary pointer, its easy
    	if (parent->type == NODE_UNARY)
      		((Unary_Node *)parent)->set_child(newnode2);
    	else {
			// it must be a binary node then
      		if (left)
        		((Binary_Node *)parent)->set_left(newnode2);
      		else
        		((Binary_Node *)parent)->set_right(newnode2); 
		}
	}
	else {
		*c2 = (Binary_Node*)newnode2; 
	}
	//------------------------------------------------------------------------------------------
	if (COMMENT) {
		cout << "\nChild 2 done" << endl;
		cout << "*c2 = " << *c2 << endl;
		cout << "parent = " << parent << " , newnode2 = " << newnode2 << endl;
	}

    // check the type of the children's root nodes
	if ((*c1)->type != NODE_BINARY) {
			if (COMMENT) {
				cout << "\nERROR! THE ROOT NODE OF *c1 IS NOT BINARY!!!!" << endl;
				cout << "A 3-NODE TREE IS BUILT MULTIPLYING THE VARIABLE BY 1" << endl;
			}
			insert_parameter((Node*)*c1, (Node**)c1, &Mult, (Val)1.) ;
	}
	if ((*c2)->type != NODE_BINARY) {
			if (COMMENT) {
				cout << "\nERROR! THE ROOT NODE OF *c2 IS NOT BINARY!!!!" << endl;
				cout << "A 3-NODE TREE IS BUILT MULTIPLYING THE VARIABLE BY 1" << endl;
			}
			insert_parameter((Node*)*c2, (Node **)c2, &Mult, (Val)1.) ;
	}

	// set the fitness of the two new offspring to a high number!
	// this declaration is
	(*c1)->fitness=9.999999E99;
	(*c2)->fitness=9.999999E99;

    n_times++; //counts the times crossover is executed
	if (COMMENT) {
		cout << "\n*c1 = " << *c1 << ", *c2 = " << *c2;

		cout << "\nPopulation::crossover : exit";	
	}
	
}



//point mutation function - in structural GP a terminal_var CANNOT be replaced by a terminal_const!!!
// terminal nodes are swapped with other terminal nodes
// functional nodes are redefined with other function of the
// same arity
// INPUT: pointer to the tree undergoing mutation (Binary_Node *p_tree), mutation point (int)
// OUTPUT: mutation node nm if ok, 0 if mutation not performed  
int Population::point_mutation(Binary_Node *tree, int nm)
{	
	int COMMENT = 0; // 1 if comments, 0 silent...
	Binary_Func *p_b_fun;
	

	// fetch the pointer to the node nm (Node *p_nm)
	Node *p_nm = tree->find(nm);  
	//get the pointer to parent 
	Node *par = p_nm->get_parent();
	//check the kind of node
	int t;
	if ((p_nm->type ==  NODE_BINARY) && (problem->num_b_funcs>1)) {
			Binary_Func *p_b_fun;
			Binary_Node *p_b_node = (Binary_Node*)p_nm;
			t = NODE_BINARY;
			// select randomly a function different from the current one and get the pointer 
			if (COMMENT) {
				cout << "\nCurrent f: p_nm->f = " << p_b_node->get_func() << " Sign = " << (p_b_node->get_func())->sign << endl; 
			}
			p_b_fun = p_b_node->get_func();
			while ((p_b_node->get_func()) == p_b_fun) {
				// select the pointer to the new function
				p_b_fun = problem->b_func_list[int_rand(problem->num_b_funcs)];
				if (COMMENT) cout << "\nChosen f: p_b_fun = " << p_b_fun << " Sign = " << p_b_fun->sign << endl; 
			}
			// call the function in Binary Node class to change the function of the node 
			((Binary_Node *)p_nm)->change_f(p_b_fun);
	}
		

	if ((p_nm->type ==  NODE_UNARY) && (problem->num_u_funcs>1)) {
			Unary_Func *p_u_fun;
			Unary_Node *p_u_node = (Unary_Node*)p_nm;
			t = NODE_UNARY;
			// select randomly a function different from the current one and get the pointer 
			if (COMMENT) {
				cout << "\nCurrent f: p_nm->f = " << p_u_node->get_func() << " Sign = " << (p_u_node->get_func())->pre_sign << (p_u_node->get_func())->post_sign;
			}
			p_u_fun = p_u_node->get_func();
			while ((p_u_node->get_func()) == p_u_fun) {
				// select the pointer to the new function
				p_u_fun = problem->u_func_list[int_rand(problem->num_u_funcs)];
				if (COMMENT)  cout << "\nChosen f: p_u_fun = " << p_u_fun << " Sign = " << p_u_fun->pre_sign << p_u_fun->post_sign;
			}
			// call the function in Unary Node class to change the function of the node 
			((Unary_Node *)p_nm)->change_f(p_u_fun); 
	}
			

	if ( ((p_nm->type ==  NODE_VAR) && (parameters->nvar>1)) || (p_nm->type ==  NODE_CONST) ) {
			
			if (p_nm->type ==  NODE_VAR) 
				t = NODE_VAR;
			else 
				t = NODE_CONST;
			
			//create new random terminal node (const or var)
			Node *new_node;
			//select which kind of node to create (function for that?)
			int kind = p_nm->type;  //int_rand(2)+2; //put 3 or 2 to create only terminal_const or terminal_var nodes, respectively. 
			///////////////////////// the whole thing can be optimized!!! Creating a new node is not needed, just replace the variable!
			// current effect: a terminal const remains the same, a terminal var is mutated with another random variable
			if (kind== NODE_CONST ) {
				// create a new terminal_const node
				Terminal_Const *newconst = new Terminal_Const((Node *)par,  1);
				new_node = (Node *)newconst;
			}
			if (kind==NODE_VAR) {
				// choose a variable different from the current one
				Variable *p_var;
				Terminal_Var *p_terminal_var = (Terminal_Var*)p_nm;
				if (COMMENT) {
					cout << "\nCurrent var: p_nm->var = " << p_terminal_var->val_p() << " Sign = " << (p_terminal_var->val_p())->name; 
				}
				p_var = p_terminal_var->val_p();
				while (p_var == p_terminal_var->val_p()) {
					// select the pointer to a new variable
					p_var = problem->v_list[int_rand(parameters->nvar)];
					if (COMMENT) {
						cout << "\nChosen var: p_var = " << p_var << " Sign = " << p_var->name; 
					}
				}
				// create a new terminal_var node
				Terminal_Var *newvar = new Terminal_Var((Node *)par, p_var);
				new_node = (Node *)newvar;
			}

			//set the parent's pointer to new node
			// parent node is unary
			if (par->type == 1)
				((Unary_Node *)par)->set_child(new_node);
			//parent is binary: nm is left or right child? Check...
			if	(par->type == 0) {
				if (((Binary_Node *)par)->get_left() == p_nm)
					// nm is left child
					((Binary_Node *)par)->set_left(new_node);
				else
					// nm is right child
					((Binary_Node *)par)->set_right(new_node);
			}
			
			// set nm parent to Null 
			p_nm->set_parent(NULL);
			//destroy node pointed by nm (already substituted)
			delete p_nm;
        }
	
	if (COMMENT) {
		cout<< "\nPopulation::mutation : node nm = "<< nm << " type = " << t << " pointer p_nm = " << p_nm;  
		cout <<	"\nPopulation::mutation : no. of nodes of the tree = " << tree->count();	
	}

	return 1;
}


//function to generate randomly a single subtree (or tree)
// input: void
// output: pointer to new subtree's root node (Binary_Node *)
Binary_Node *Population::generate_subtree(int max_depth_subtree)
{
	int COMMENT =0;
	Val value;
	char *expr;
	if (COMMENT) cout << "\nPopulation::generate_subtree: CREATION OF A NEW SUBTREE" << endl;
	
	// create the root of the new subtree (the first argument of the constructor is the pointer to parent, NULL, as the node is root)
	p_new_subtree = new Binary_Node(NULL,problem->b_func_list[int_rand(problem->num_b_funcs)]);    //int_rand(2)+2]);    // 2]); //gives multiplication!

 	// set the maximum depth of the new subtree ()FULL and GROW methods) or the range of the depth
	depth_max = max_depth_subtree;    
	//depth_max = int_rand(max_depth_subtree-d_min+1)+d_min;   //where d_min is taken from initialization parameters (see DEPTH_MIN in input_file) 
	full = int_rand(2);	//choose randomly if FULL (full=1) or GROW (full=0) mode

	// set the depth of the current node, root.
	depth = 1;

	// build recursively the new (sub)tree, without parameters
	build_tree(p_new_subtree);

	//the method used for initialization is the same as used for initial population generation  
	if (COMMENT) {
		switch (parameters->method) {
	
			case 1:	// no limits (original method)
				{
					printf ("\n INITIALIZATION with NO LIMITS method");
				}
				break;

			case 2:  // FULL method
				{
					printf ("\n INITIALIZATION with FULL method");
				}
				break;

			case 3:	//GROW method
				{
					printf ("\n INITIALIZATION with GROW method");
				}
				break;

			case 4:	//RAMPED method
				{
					printf ("\n INITIALIZATION with RAMPED method");
				}
			break;
		}
	}
	
	// print the new subtree before parameter insertion
	expr = p_new_subtree->print();  
	if (COMMENT) {
		cout << "\nSingle new subtree created: parameters not introduced yet" << endl;
		cout << "\nDepth = " << p_new_subtree->calc_depth() << " \nExpr =" << expr << endl;
		delete [] expr;
	}

	//insert parameters to the "structure" of the tree (all of them, according to the usual algorithm)
	//parameters_allocation(p_new_subtree, &p_new_subtree);    //ONLY THE STRUCTURE IS CREATED in structuralGP!!! 
/*//
	//insert a single parameter to the "structure" of the tree
	//IT'S A TRY! ADDING A SINGLE PARAMETER TO THE STRUCTURE (x^0...) 
	// choose a random value for the final ("shift") parameter
	value = constant_generation();
	// insert the external parameter with the chosen value (final parameter - shift of the total tree )
	insert_parameter((Node*)p_new_subtree, (Node **)&p_new_subtree, &Add, value);     
	if (COMMENT) 
		cout << "\nFree parameter " << value << " added (shift of the tree)";
//*/

	// print the expression after parameter insertion
	if (COMMENT) {
		expr = p_new_subtree->print();  
		cout << "\nSingle new subtree created: after parameter insertion" << endl;
		cout << "\nDepth = " << p_new_subtree->calc_depth() << " \nExpr =" << expr << endl;
		delete [] expr;
	}

	// return the pointer to the root of the just created new subtree.
	// As it is dynamically allocated, it remains allocated even outside of the function
	return p_new_subtree;
}





// function to prune the branches that are multiplied by the "least important" tuning parameters
void Population::prune_tree (Binary_Node *p_tree)
{
	int COMMENT =0;  //1 comments, 0 silent
	char* expr;
	double min, max, cur;
	int n_param, n_param_copy; 
	int min_pos, max_pos;

	if (COMMENT) {
		cout << "Population::prune_tree" << endl;
	}

	//retrieve the parameters of the selected tree
	// (p_par already exists and is updated, contains the addresses of the parameters' nodes)
	// count the number of parameters and return if there is only one
	n_param = p_tree->p_par.size();
	if (COMMENT)
		cout << "\nThe selected tree has " << n_param << " parameters" << endl;

	if (n_param==1) {
		cout << "\nOnly one parameter. Exit prune_tree()" << endl;
		return;
	}
	
	// Now we are sure that there is more than 1 parameter.
	// Copy the selected tree whose pointer is given as input
	// (it should be the best tree of the population, first individual of the archive - complete_trees[0] - so make sure that sorting has been done!)
	Binary_Node *copyp_tree = (Binary_Node*)tree_copy((Node*)p_tree,NULL);
	 
	if (COMMENT) {
		expr = copyp_tree->print();
		cout << "\nPRUNING";
		cout << "\nIndividual selected for pruning:\n " <<  expr << endl;
		delete [] expr;
	}
	
	// Retrieve the parameters of the copied tree (their addresses are stored in copy_tree->p_par)
	n_param_copy = copyp_tree->find_parameters(); 
	if (n_param != n_param_copy) {
		cout << "\nERROR! THE SELECTED TREE AND ITS COPY HAVE A DIFFERENT NUMBER OF PARAMETERS! EXIT..." << endl;	
		exit(-1);		
	}
	if (COMMENT) {
		cout << "\nThe selected tree has " << n_param << " parameters" << endl;
		cout << "Parameters' value: ";	
		for (int k=0; k<n_param_copy; k++) 
			cout << ((Terminal_Const *)copyp_tree->p_par[k])->value(NULL) << " ";
	}
	
	
	// find the maximum and the minimum absolute values of the parameters
	min = fabs( ((Terminal_Const *)copyp_tree->p_par[0])->value(NULL) );
	max = min;
	min_pos = 0;
	max_pos = 0;
	for (int k=1; k<n_param; k++) {
		cur = fabs(((Terminal_Const *)copyp_tree->p_par[k])->value(NULL)); //values of the parameters
		// update minimum
		if (cur<min) {
			min=cur;
			min_pos = k;
		}
		// update maximum
		if (cur>max) {
			max=cur;
			max_pos = k;
		}
	}

	if (COMMENT) {
 		cout << "\nMaximum absolute value: " << max << ", found in position " << max_pos;
		cout << "\nMinimum absolute value: " << min << ", found in position " << min_pos;
	}
	

	// create an array of the same size of p_par and store the importance of each parameter
	double *importance = new double[n_param_copy];
	// set importance: for now it's just the absolute value of the parameter...
	importance[0] = fabs(((Terminal_Const *)copyp_tree->p_par[0])->value(NULL));
	min = importance[0];
	min_pos = 0;
	
	for (int k=1; k<n_param_copy; k++) {
		// set importance: for now it's just the absolute value of the parameter...
		importance[k] = fabs(((Terminal_Const *)copyp_tree->p_par[k])->value(NULL));
		cur = importance[k];
		// update minimum
		if (cur<min) {
			min=cur;
			min_pos = k;
		}
	}


	if (COMMENT) {
		cout << "\nImportance : ";
		for (int k=0; k<n_param_copy; k++) {
			cout << importance[k] << " ";
		}
	}

	// retrieve the address of the least "important" constant node
	Node *p_least = copyp_tree->p_par[min_pos];

	if (COMMENT)
		cout << "\nThe least important parameter is " << ((Terminal_Var*)p_least)->value(NULL) << ", adddress " << p_least;	

	// get the address of the parent of the least important parameter
	Node *parent = p_least->get_parent();
	int left = 0;
	if (((Binary_Node*)parent)->get_left() == p_least) {
		left = 1;
		if (COMMENT) cout << "\nThe selected node is a left child";
	}
	else {
		left = 0;
		if (COMMENT) cout << "\nThe selected node is a right child";
	}
		

	// prune the block corresponding to the least "important" parameter
	// eliminate the subtree only if the parent is a multiplication or a division, and in this case the parameter must be left child 
	if (parent->type == NODE_BINARY)
		if ( (((Binary_Node*)parent)->get_func() == &Mult) ||    ((((Binary_Node*)parent)->get_func()== &SDiv) && left) )  {
			eliminate_subtree(parent, (Val)0.0);
			cout << "\n CUT!" << endl;
		}


	if (COMMENT) {
		cout << "Subtree rooted in " << parent << " has been successfully replaced by 0" << endl;
		expr = copyp_tree->print();
		cout << expr ; 
		delete [] expr;		
	}
	
	// free memory occupied by importance_list
	delete[] importance;

	//copy the pruned tree back into the new_trees

	if (COMMENT)
		cout << "\nPopulation:prune_tree  .... normal exit" << endl;
}




// FULL method to generate the initial population of parameterless trees
int Population::FULL_method(void)
{
	int t; //type of the node
	if (depth<depth_max) {
		if (problem->num_u_funcs)
			t = int_rand(2);    // NODE_BINARY 0    NODE_UNARY 1
		else
			t=0;				// NODE_BINARY 0
	}			
	else
		t = 2;	  // NODE_VAR 2

	return t;
}





// GROW method to generate the initial population of parameterless trees
int Population::GROW_method(void)
{
	int t;
	if (depth<depth_max) {
		if (problem->num_u_funcs)
			t = int_rand(3);	   // NODE_BINARY 0    NODE_UNARY 1    NODE_VAR 2
		else 
			t = 2*int_rand(2);     // NODE_BINARY 0    NODE_VAR 2
	}		
	else
		t = 2;

	return t;
}


double Population::compute_time(time_t start, time_t finish, double *p_delta_t)
{	
	int COMMENT =0;

	time(&finish);
   *p_delta_t = difftime(finish, start);


	if (COMMENT) {	
		cout << "\nElapsed time  (hrs:mins:secs) :  " ;
		cout << (int)(*p_delta_t/3600.0) << " : ";    //hours
		cout << (int)(fmod(*p_delta_t, 3600.0))/60 << " : ";  //minutes
		cout << (int)fmod((fmod(*p_delta_t, 3600.0)), 60.0); //seconds

		cout << "       (total secs " << *p_delta_t << ")";  //total seconds
	}

	return (double)(*p_delta_t);
}


void Population::print_individual(Node *tree)
{
	char *expr;
	expr = tree->print();
	cout << " " << expr << " ";
	delete [] expr;
	return;
}


void Population::print_population_without_parameters(int gen)
{
	// if COMMENT=1 prints all trees, if COMMENT=0 prints only the best tree (lowest F)
	int COMMENT = 1;

	char *expr;	
	cout <<"\n\n --------------------  GENERATION " << gen << " without parameters -------------------" << endl; 
	cout << " n F fitness hits n_nodes depth expr" << endl; 	
	for (int j=0;j<size;j++) {
		expr = print(j,trees);
		cout << j << "  " << scientific << trees[j]->F << " " << trees[j]->fitness;
		cout << "  " << trees[j]->hits << "  " << trees[j]->count();
		cout << "  " << trees[j]->calc_depth() << "  " << expr << endl;
		
		// free memory used for single expression	
		delete [] expr;		

		//print only the best parameterless tree if COMMENT=0
		if (!COMMENT) break; 
	}
	cout << "------------------------------------------------------------------------" << endl; 

}


void Population::print_population_with_parameters(int gen)
{
	// if COMMENT=1 prints all trees, if COMMENT=0 prints only the best tree (lowest F)
	int COMMENT =0;

	cout << "\nPopulation::print_population_with_parameters()" << endl;

	int k=0;
	while (k<size) {
		if (complete_trees[k]==NULL) {
			cout << "\ncomplete_trees[] not initialised. Exit function";
			return;
		}
		k++;
	}

	cout <<"\n\n --------------------  GENERATION " << gen << " with parameters -------------------" << endl; 
	cout << " n F RMSE R2 n_nodes depth expr" << endl;
	//for (int j=0;j<(int)(floor(parameters->repr_rate*parameters->M));j++) {  // to print the first repr_rate*size trees

	for (int j=0;j<size; j++) {  // to print all the trees in the population

		// usual approach to print trees' features
		cout << "\nComplete tree n. " << j;
		complete_trees[j]->show_state(); // to show all the data referring to the trees

		cout << endl;

		//prints only the best complete individual if COMMENT=0
		if (!COMMENT) break; 
	}

	cout << "------------------------------------------------------------------------\n" << endl; 
	cout << "Population::print_population_with_parameters: selected strategy = " << parameters->strat_statp << endl;
	
}

// sorts the generation according to F (not fitness=RMSE)
// (see tree_functions/function copy and sort tree.cpp/ tree_comp)
void Population::sort (int gen, int (*p_tree_comp)(const void *, const void *))
{
	int COMMENT =0;

	if (COMMENT)  {
		cout << "\nPopulation::sort";
		cout << "\nGeneration " << gen; 
	}
	// sort the population in trees[] in descending order
    qsort(trees,size,sizeof(Binary_Node *),p_tree_comp);

	// sort the population in complete_trees[] in descending order
    qsort(complete_trees,size,sizeof(Binary_Node *),p_tree_comp);
}



int Population::tournament(int first, int last, int tournament_size)
{
	int COMMENT =0; //1 comments, 0 silent...
	if (COMMENT)
		cout << "\n\nPopulation::tournament()" ;
 
	int s; 
	struct fighter {
		double quality;   // will be a val...
		int ID_number;
	};

	// allocate two "fighters": one is the current winner, the other is the current rival
	fighter competitor;
	fighter winner;
	
	// select randomly tournament_size fighters...and assign fitness and ID number (referring to trees array)
	winner.ID_number = int_rand(last-first+1) + first;

	winner.quality = trees[winner.ID_number]->F;     // selection based on F

	//winner.quality = trees[winner.ID_number]->fitness;    // selection based on error
	if (COMMENT) {
		cout << "\nFighter 1" << endl;
		cout << "  ID_number = " << winner.ID_number << endl; 
		cout << "  F = " << winner.quality << endl;
	}

	for (int i=1; i<tournament_size; i++) {
		//  record the individual ID and fitness (has been chosen to be a fighter)
		competitor.ID_number = int_rand(last-first+1) + first; //see position in trees array
		competitor.quality = trees[competitor.ID_number]->F;
		if (COMMENT) {
			cout << "Fighter " << i+1 << endl;
			cout << "  ID_number = " << competitor.ID_number << endl; 
			cout << "  F = " << competitor.quality << endl;
		}
		if (competitor.quality < winner.quality) {    //lower quality is better...	
			winner.quality = competitor.quality;
			winner.ID_number = competitor.ID_number;
		}
		
	}

	// winner declaration
	if (COMMENT) 
		cout << "Winner's ID_number -> " <<  winner.ID_number << endl;
	return winner.ID_number; 
}



void Population::kill_and_fill (ProblemDefinition *pb)
{
	int COMMENT = 0;
	int kill_tot = comp_tot+new_tot;
	cout << "\nPopulation::kill_and_fill" << endl;
	cout << "comp_tot = " << comp_tot << " , new_tot = " << new_tot << " , kill_tot = " << kill_tot << endl;
	if ((!comp_tot) && (!new_tot)) {
		cout << "Kill_and_fill not applied... Exit" << endl;
		return;	
	}
		
	
	char* expr;
	char* cursed_op = SDiv.sign; // previously " / ";
	int p1, p2, count, valid_composition;
	int tree_depth_max = 10000;
	Binary_Node* new_left_child;
	Binary_Node* new_right_child;
	Binary_Node** new_trees;


	//-----------------------------------------------------------------------------------------------------------------------------
	// INNOVATE : introduce new individuals in the population 
	//-----------------------------------------------------------------------------------------------------------------------------
	if (COMMENT) cout << "\n\nINNOVATE" << endl;
	
	if (!new_tot)
		cout << "new_tot = " << new_tot << ". Not applied." << endl;
	else {
		new_trees  = new Binary_Node *[new_tot];

		for (int k=0; k<new_tot; k++) {
			// generate a new tree with a maximum depth specified as input
			new_trees[k] = generate_subtree(3);     //if 1 ONLY NEW BINARY FUNCTIONS ARE INTRODUCED!! 
		}

		// COPY THE NEW TREES BACK INTO THE POPULATION
		//delete the memory allocated to the tree to be substituted
		for (int k=0; k<new_tot; k++) {
			delete trees[size-comp_tot-k-1];
			// copy the new tree to trees[size-k-1]
			trees[size-comp_tot-k-1] = (Binary_Node *)tree_copy(new_trees[k],NULL); 
			trees[size-comp_tot-k-1]->fitness = 444444.;
		
			if (COMMENT) {
				expr = trees[size-comp_tot-k-1]->print();
				cout << k+1 << ") new tree :  " << expr << endl;
				delete [] expr;		
			}
		}
	
		// free memory
		for (int k=0; k<new_tot; k++) delete new_trees[k];
		delete[] new_trees;
	}

	//-----------------------------------------------------------------------------------------------------------------------------
	// COMPOSE : take two trees and consider them as children of a randomly chosen binary function 
	//-----------------------------------------------------------------------------------------------------------------------------
	Binary_Func* bin_func;

	if (COMMENT) cout << "\nCOMPOSE" << endl;

	if (!comp_tot)
		cout << "comp_tot = " << comp_tot << ". Not applied." << endl;
	else {
		new_trees  = new Binary_Node *[comp_tot];
	
		for (int k=0; k<comp_tot; k++) {
			//-------------------------------------------------------------------------------
			// generate a new binary node (root node, randomly chosen between binary functions)
			// 1 - each binary operation can be chosen
			//  bin_func = problem.b_func_list[int_rand(problem.num_b_funcs)];
			// 2 - division cannot be chosen
			//while (!strcmp(bin_func->sign,cursed_op))
			//		bin_func = b_funcs[int_rand(num_b_funcs)];
			// 3 - the root node IS protected division (in order to ease the search for rational expressions - Kotanchek, RatPol2D)
			if (pb->division)
				bin_func = pb->division;
			else {
				cout << "\n\nComposition with division as root node : " << endl;
				cout << "Division not among primitives : composition skipped." << endl;
				// you should delete the arrays "new trees" previously allocated...
				return;
			};

			// -----------------------------------------------------------------------------
			new_trees[k] =  new Binary_Node(NULL,bin_func); 
		
			if (COMMENT) cout << "\nTree n. " << k << " . Selected function : " << bin_func->sign << endl;

			//select two parents randomly, checking the maximum depth constraint
			count = 0;
			valid_composition = 0;	
			while (!valid_composition) {
				// tournament selection doesn't make much sense, as the fitness of the composed individual is not related to the two chosen for composition...
				//p1 = tournament(0, composition[0]+composition[1]-1, 3);//(0, repr_tot-1, 3);		// = int_rand(repr_children); //TRUNCATION with T=repr_rate:parent 1 taken from the copied elite
				//p2 = tournament(composition[0]+composition[1], size-1, 3);				// = int_rand(size); // for TRUNCATION with T=1: parent 2 taken from the entire population
				p1 = int_rand(size);
				p2 = int_rand(size);
				tree_depth_max = max(trees[p1]->calc_depth(), trees[p2]->calc_depth()) + 1;  //+1 is for the new binary root node
				if (((tree_depth_max<=parameters->depth_lim) || (count==9)) && (p1!=p2)) {
					if (COMMENT) cout << "\nComposition accepted" << endl;
					valid_composition = 1;
				}
				count++;
			}
		
			if (COMMENT) {
				expr = trees[p1]->print();
				cout << "\nFirst tree : " << p1 << ") " << expr << endl;
				delete [] expr;		
				expr = trees[p2]->print();
				cout << "\nSecond tree : " << p2  << ") " << expr << endl;
				delete [] expr;		
			}
		
			// copy the selected individuals (not deleting the original ones!)
			new_left_child = (Binary_Node *)tree_copy(trees[p1],NULL); 
			new_right_child =  (Binary_Node *)tree_copy(trees[p2],NULL); 
		
			// make connection from new root to left and right child
			(new_trees[k])->set_left((Node*)new_left_child);
			(new_trees[k])->set_right((Node*)new_right_child);
		
			// make connection from children to new root
			((Node*)new_left_child)->set_parent((Node*)new_trees[k]);
			((Node*)new_right_child)->set_parent((Node*)new_trees[k]);

		}	

		
		// COPY THE NEW TREES BACK INTO THE POPULATION
		//delete the memory allocated to the tree to be substituted
		for (int k=0; k<comp_tot; k++) {
			delete trees[size-k-1];
			// copy the new tree to trees[size-k-1]
			trees[size-k-1] = (Binary_Node *)tree_copy(new_trees[k],NULL); 
			trees[size-k-1]->fitness = 555555.;
		
			if (COMMENT) {
				expr = trees[size-k-1]->print();
				cout << k+1 << ") Composed tree : " << expr << endl;
				delete [] expr;		
			}
		}
	
		// free memory
		for (int k=0; k<comp_tot; k++) delete new_trees[k];
		delete[] new_trees;
	}
}


// function that performs the reproduction of a portion of the population (elite)
// input: pointer to old population, pointer to new population, no. of trees to copy, flag to keep or eliminate copies, pointer to counter variable)
// output: none
void Population::population_reproduction(Binary_Node **trees, Binary_Node **new_trees, int n_repr, bool flag_no_copies, int *counter)
{
	int COMMENT =0;
	char *expr;

	// counter[0]: stores the number of trees belonging to old population copied to the new one
	// counter[1]: stores the number of trees generated from scratch to reach n_repr individuals in the new population

	 //REPRODUCTION (elite)

	int unique = 1;
	if (COMMENT) cout << "\nPopulation::population_reproduction() : start" << endl;
	if (n_repr) {
		// copy best tree of old population trees[0] to new_trees[0] (not identical to any other tree, as it is the first one)
		new_trees[0] = (Binary_Node *)tree_copy(trees[0],NULL);
		// update son_of and parent_fitness
		new_trees[0]->parent_fitness = trees[0]->fitness;
		new_trees[0]->son_of = 0;
		counter[0]++;

		//if (COMMENT) {
			expr = trees[0]->print();
			cout << "\ntrees[0] = " << expr;
			delete [] expr;
			expr = new_trees[0]->print();
			cout << "\n ==> new_trees[0] = " << expr;
			delete [] expr;
		//}

		int i;
		for (i=1; i<n_repr; i++) {   //for (i=1; i<n_repr; i++) {
			// TRUNCATION
			// copy the i-th individual in trees[] in i-th individual in new_trees[]
			unique = 1;
			if (COMMENT) {
				expr = trees[i]->print();
				cout << "\ntrees[" << i << "] = " << expr;
				delete [] expr;
			}
			//start NO COPIES---------
			//if flag_no_copies is true, copies are filtered out, that is only unique individuals are copied into the new population
			if (flag_no_copies) {
				for (int j=0; j< counter[0]; j++) {
					if (COMMENT) {
						expr = new_trees[j]->print();
						cout << "\n		new_trees[" << j << "] = " << expr;
						delete [] expr;
					}
					if (identical_trees(trees[i], new_trees[j])) {
						if (COMMENT) cout << "\nIdentical trees : trees["<< i << "] not copied!" << endl;
						unique = 0;
					}
					if (!unique) break;
				}
			}
			// end NO COPIES--------------------------------------

			if (unique) {
				new_trees[counter[0]] = (Binary_Node *)tree_copy(trees[i],NULL);
				// update son_of and parent_fitness
				new_trees[counter[0]]->parent_fitness = trees[i]->fitness;
				new_trees[counter[0]]->son_of = 0;
 				if (COMMENT) {
 					expr = new_trees[counter[0]]->print();
 					cout << "\n ==> new_trees[" << counter[0] << "] = " << expr;
 					delete [] expr;
 				}
 				counter[0]++;
			}
		}
 	// end individual copying (end if (n_repr) {)
	}
 
 	// FILL THE MISSING TREES (DELETED COPIES...) WITH RANDOMLY INITIALISED INDIVIDUALS
 	int n_missing = n_repr - counter[0] - counter[1];
 	if (COMMENT) {
 		if (!n_missing)
 			 cout << "\nArchive filled : no other individuals will be copied into" ;
 		else
 		if (COMMENT)
 			cout << "\nFill in the archive (filling the empty spaces...)" ;
 	}
 	while (counter[0] + counter[1] < n_repr) {
 		new_trees[counter[0] + counter[1]] = generate_subtree(2); //here depth = 2 : if depth = 1 ONLY NEW BINARY FUNCTIONS ARE INTRODUCED!!
 		// update son_of and parent_fitness
 		new_trees[counter[0] + counter[1]]->parent_fitness = trees[0]->fitness;
 		new_trees[counter[0] + counter[1]]->son_of = 0;

 		if (COMMENT) {
 			expr = new_trees[counter[0] + counter[1]]->print();
 			cout << "\nnew_trees[" << counter[0] + counter[1] << "] = " << expr;
 			delete [] expr;
 		}
 		counter[1]++;
 	}

 	// here n_repr = counter[0] + counter[1]  !
 	if (counter[0] + counter[1] != n_repr) {
		cerr << "\nPopulation::population_reproduction : ERROR!";
		cerr << "\nThe total number of individuals copied/generated ( "<< counter[0] + counter[1] << ") is smaller than n_repr ( " << n_repr << ")" << endl;
		exit(-1);
 	}

 	if (COMMENT) {
 		cout << "\nPopulation::population_reproduction : n_repr = " << n_repr << " counter[0] = " << counter[0] <<  " counter[1] = " << counter[1];
 		cout << "\nPopulation:: REPRODUCTION :";
 		cout << "\ncounter[0] + counter[1] = " << counter[0] + counter[1] << " (individuals actually copied - till now " << (double)((counter[0] + counter[1])*100/size) << " perc. of population)" << endl;
 		cout << "\nList of individuals copied: " << endl;
 		for (int k=0; k< counter[0] + counter[1]; k++) {
 			if (k==counter[0] + counter[1]-n_missing) cout << "\n --------  new " << n_missing << " individuals generated to fill the archive: ------------------";
 			expr = new_trees[k]->print();
 			cout << "\nnew_trees[" << k << "] = " << expr;
 			delete [] expr;
 		}
     }

}

// function that performs crossover
// input: pointer to old population, pointer to new population, no. of trees to cross, pointer to counter variable)
// output: none
void Population::population_crossover(Binary_Node **trees, Binary_Node **new_trees, int n_cross, int *counter)
 {
	 char *expr;
	 int COMMENT = 0;
	 if (COMMENT)
	 	cout << "\nCROSSOVER";
	 int cross_point[2];
	 int crossover_invalid =1;
	 int d1,d2;
	 int first;
	 int p1 = 0;
	 int p2 = 1;
	 int crossover_point1 = 1;  // if crossover can't be performed due to depth limit reasons, just swap trees...
	 int crossover_point2 = 1;  // if crossover can't be performed due to depth limit reasons, just swap trees...

	 for (int i=0;i<n_cross;i+=2) {  //two individuals are created at a time...
	   	if (COMMENT) cout << "\n\n Children " << i << " and " << i+1;
		while (crossover_invalid) {
	 		//select two parents randomly ( build a function for it!)
	 		p1 = tournament(0, counter[0]+counter[1]-1, 3);		//TRUNCATION with T=repr_rate:parent 1 taken from the copied elite. NON ADAPTIVE VERSION
	 		//p1 = tournament(0, (int)(floor(mut_rate*size))-1, 3); // ADAPTIVE VERSION
	 		p2 = tournament(0, size-1, 3);				// = int_rand(size); // for TRUNCATION with T=1: parent 2 taken from the entire population
	 		// select two nodes randomly (crossover points), one in each parent (use member variables crossover_point1 and crossover_point2)
	 		crossover_point1 = select_node(trees[p1]);
	 		crossover_point2 = select_node(trees[p2]);
	 		// check that the offspring generated by crossover don't exceed limits on depth (see depth_limit)
	 		// first argument is the pointer to the subtree to be replaced, the second the pointer to the subtree to be introduced
	 		d1 = potential_depth((Node*)(trees[p1]->find(crossover_point1)), (Node*)(trees[p2]->find(crossover_point2)));
	 		d2 = potential_depth((Node*)(trees[p2]->find(crossover_point2)), (Node*)(trees[p1]->find(crossover_point1)));
	 		if ((d1<=parameters->depth_lim) && (d2<=parameters->depth_lim)) {
	 			if (COMMENT)  cout << "\n Crossover is valid (offspring depths within the limit depth_lim): new max depths d1 =" << d1 << " , d2 = " << d2 << endl;
	 			crossover_invalid = 0;
	 		}
	 		else {
	 			if (COMMENT) cout << "\n Crossover is NOT valid for the parameters chosen (offspring depths exceed depth_lim): max depth d1 =" << d1 << " , d2 = " << d2 << endl;
	 		}
	 	}
	 	// here the offspring will not exceed depth limit: perform crossover
	 	//the first two parameters are pointers to the parents,
	 	//the last two param. are pointers to pointers to the children got from crossover and the last is an array with information)
	  	int pos = counter[0] + counter[1] + counter[2];
		crossover(trees[p1],trees[p2], crossover_point1, crossover_point2 ,&new_trees[pos],&new_trees[pos+1], cross_point);
		// update son_of and parent_fitness
		double best_parent_fitness = min(trees[p1]->fitness,trees[p2]->fitness);
		new_trees[pos]->parent_fitness = best_parent_fitness;
		new_trees[pos]->son_of = 1;
		new_trees[pos+1]->parent_fitness = best_parent_fitness;
		new_trees[pos+1]->son_of = 1;

		// set invalid flag back to invalid (important)
     	crossover_invalid = 1;
     	if (COMMENT) {
     		// parents
     		expr = print(p1,trees);
     		cout << "\nParent 1:";
     		cout << "\nFitness : " << trees[p1]->fitness << endl;
     		cout << expr << endl;
     		delete [] expr;
     		expr = print(p2,trees);
     		cout << "\nParent 2:";
     		cout << "\nFitness : " << trees[p2]->fitness << endl;
     		cout << expr << endl;
     		delete [] expr;
     		// children
     		expr = print(pos,new_trees);
			cout << "\nChild 1:";
			cout << "\nson_of : " << new_trees[pos]->son_of;
			cout << "\nParent_Fitness: " << new_trees[pos]->parent_fitness << endl;
			cout << expr << endl;
			delete [] expr;
			expr = print(pos+1,new_trees);
			cout << "\nChild 2:";
			cout << "\nson_of : " << new_trees[pos+1]->son_of;
			cout << "\nParent_Fitness: " << new_trees[pos+1]->parent_fitness << endl;
			cout << expr << endl;
			delete [] expr;
			printf("\nPopulation::new_spawn : CROSSOVER : tree %i and % i => p1=%i p2=%i cross_points %i %i", pos,pos+1,p1,p2,cross_point[0],cross_point[1]);
		}

		counter[2]+=2;
	}

	if (COMMENT) {
		cout << "\n\nPopulation::population_crossover : n_cross = " << n_cross << ",  counter = " << counter[2];
		cout << "\nPopulation::population_crossover : CROSSOVER :  = " << counter[2] << " individuals created through crossover (till now " << (double)(counter[2]*100/size) << " perc. of population)";
	}
	if (abs(n_cross-counter[2]) > 1) {
			cerr << "\nPopulation::population_crossover : ERROR!";
			cerr << "\n n_cross " << n_cross << ") and counter[2] (" << counter[2] << ") are too different!!";
			exit(-1);
	}
 }


 void Population::population_mutation(Binary_Node **trees, Binary_Node **new_trees,
 																	int n_mut, int gen, int *counter)
 {
	 char *expr;
	 int COMMENT =0;

	 if (COMMENT)
	 	cout << "\n\n\nMUTATION";
	 int dm;
	 int mutation_invalid = 1;
	 int pm = size - 1; //if no individuals are valid, the worst is chosen from trees...
	 int node_mutation = 1; // if no nodes are valid, the root is chosen (so instead of a mutation there is an introduction of a new tree...)
	 Node *p_node_mutation;

	 // DIFFERENT MUTATION ACCORDING TO THE NUMBER OF GENERATION
	 if (gen%2) {
	 		if (COMMENT)  cout << "\nGeneration " << gen << " is ODD: performing point mutation";
	 		// ----------------------------------------------------------------------
	 		// POINT MUTATION: as changes functions with same arity, if there are not unary functions, it is not able to introduce them!
	 		//-----------------------------------------------------------------------

	 		for (int i=0; i<(size-(counter[0]+counter[1]+counter[2])); i++) {     	// instead of i<(size-cross_tot-repr_tot), so that not more than "size" individuals are actually created
	 			while (mutation_invalid) { //££
	 				//select an individual for mutation from the ones copied by reproduction
	 				pm = tournament(0, counter[0]+counter[1]-1, 3);	 //TRUNCATION with T=repr_rate NON ADAPTIVE VERSION
          //pm = tournament(0, (int)floor(mut_rate*size)-1, 3);    // ADAPTIVE VERSION   
	 				//select a mutation point in the chosen individual
	 				node_mutation = select_node(trees[pm]);
	 				mutation_invalid = 0;
	 			}//££

	 			if (COMMENT) {
	 					cout << "\n\nPopulation::new_spawn : tree undergoing mutation = " << pm << ", mutation point = " << node_mutation << endl;
	 			}

	 			// copy the tree to mutate in the new population new_trees
	 			if (COMMENT) cout << "Copy the selected tree to new_trees[no_so_far] " << endl;
	 			int no_so_far = counter[0]+counter[1]+counter[2]+counter[3];
	 			new_trees[no_so_far] = (Binary_Node *)tree_copy(trees[pm],NULL);

	 			// do the mutation
	 			//point mutation. Arguments: pointer to the root, node number (//remove the forward slashes to have point mutation)
	 			//remember to set the fitness of the mutated tree to 999999E99!
	 			if (!point_mutation(new_trees[no_so_far], node_mutation))
	 				cout << "\nMutation on tree " << i << " not performed" << endl;
	 			// update son_of and parent_fitness
	 			new_trees[no_so_far]->parent_fitness = trees[pm]->fitness;
	 			new_trees[no_so_far]->son_of = 3;

	 			if (COMMENT) {
	 				cout << "\nCheck:";
	 				cout << "\n	Address of tree root = " << new_trees[no_so_far];
	 				cout << "\n	Address of address of tree root = " << &new_trees[no_so_far];
	 				cout << "\n	Address of the selected node of the old subtree = " << (Node*)(new_trees[no_so_far]->find(node_mutation));
	 				cout << "\n	Address of the parent of selected node of the old subtree = " << (Node*)(new_trees[no_so_far]->find(node_mutation))->get_parent();
	 				cout << endl;
	 				// parent
	 				expr = print(pm,trees);
	 				cout << "\nParent:";
	 				cout << "\nFitness : " << trees[pm]->fitness << endl;
	 				cout << expr << endl;
	 				delete [] expr;
	 				// children
	 				expr = print(no_so_far,new_trees);
	 				cout << "\nChild:";
	 				cout << "\nson_of : " << new_trees[no_so_far]->son_of;
	 				cout << "\nParent_Fitness: " << new_trees[no_so_far]->parent_fitness << endl;
	 				cout << expr << endl;
	 				delete [] expr;
	 			}


	 			// set the flag
	 			mutation_invalid = 1; //important!


	 			// update counter
	 			counter[3]++;

	 		}
	 	//-------------------------------------------------------------------------
	 	//-------------------------------------------------------------------------
	 	}
	 	else {
	 		if (COMMENT) cout << "\nGeneration " << gen << " is EVEN: performing subtree mutation";
	 		// ----------------------------------------------------------------------
	 		// SUBTREE MUTATION: it is able to introduce any kind of function (binary and unary) IF the new subtree has depth greater than 1!
	 		//-----------------------------------------------------------------------
	 		int max_depth_new_subtree, root_dist;
	 		for (int i=0; i<(size-(counter[0]+counter[1]+counter[2])); i++) {     	// instead of i<(size-cross_tot-repr_tot), so that not more than "size" individuals are actually created
	 			while (mutation_invalid) { //££
	 				if (COMMENT) cout << "\n\nSelection of tree and node for subtree mutation" << endl;
	 				//select an individual for mutation
	 				pm = tournament(0, counter[0]+counter[1]-1, 3);	//TRUNCATION with T=repr_rate NON ADAPTIVE VERSION
          //pm = tournament(0, (int)floor(mut_rate*size)-1, 3);    // ADAPTIVE VERSION   
	 				//select a mutation point in the chosen individual
	 				node_mutation = select_node(trees[pm]);
	 				p_node_mutation = trees[pm]->find(node_mutation);
	 				root_dist = evaluate_root_distance(p_node_mutation);
	 				max_depth_new_subtree = parameters->depth_lim - root_dist;
	 				if (COMMENT) {
	 					cout << "\n\nRoot distance of node " << node_mutation << " in tree " << pm << " = " << root_dist;
	 					cout << "\nMaximum new subtree depth = " << max_depth_new_subtree;
	 					cout << "\nDepth limit = " << parameters->depth_lim << endl;
	 				}
	 				// generate a new subtree with a maximum depth specified as input
	 				p_new_subtree = generate_subtree(3);     //if 1 ONLY NEW BINARY FUNCTIONS ARE INTRODUCED!!
	 				// check that the offspring generated by mutation don't exceed limits on depth (see depth_limit)
	 				// first argument is the pointer to the subtree to be replaced, the second the pointer to the subtree to be introduced
	 				dm = potential_depth((Node*)(trees[pm]->find(node_mutation)), (Node*)(p_new_subtree));
	 				if (dm<=parameters->depth_lim) {
	 					if (COMMENT) cout << "\n Subtree mutation is valid (offspring depths within the limit depth_lim): new max depth  dm =" << dm << endl;
	 					mutation_invalid = 0;
	 				}
	 				else {
	 						//if (COMMENT)
	 							cout << "\n Subtree mutation is NOT valid for the parameters chosen (offspring depth exceed depth_lim): potential max depth dm =" << dm << endl;
	 				}

	 			}//££

	 			if (COMMENT) {
	 					cout << "\n\nPopulation::new_spawn : tree undergoing mutation = " << pm << ", mutation point = " << node_mutation << endl;
	 			}

	 			// copy the tree to mutate in the new population new_trees
	 			if (COMMENT) cout << "Copy the selected tree to new_trees[no_so_far] " << endl;
	 			int no_so_far = counter[0]+counter[1]+counter[2]+counter[3];
	 			new_trees[no_so_far] = (Binary_Node *)tree_copy(trees[pm],NULL);

	 			// subtree mutation. Arguments: pointers to old and new subtrees' root nodes.
	 			//remember to set the fitness of the mutated tree to 999999E99!
	 			if (!insert_subtree((Node**)(&(new_trees[no_so_far])), (Node*)(new_trees[no_so_far]->find(node_mutation)), (Node*)p_new_subtree))
	 					if (COMMENT) cout << "\nMutation on " << i << "tree not performed" << endl;

	 			// update son_of and parent_fitness
	 			new_trees[no_so_far]->parent_fitness = trees[pm]->fitness;
	 			new_trees[no_so_far]->son_of = 2;

	 			if (COMMENT) {
	 				cout << "\nCheck:";
	 				cout << "\nAddress of tree root = " << new_trees[no_so_far];
	 				cout << "\nAddress of address of tree root = " << &new_trees[no_so_far];
	 				cout << "\nAddress of the selected node of the old subtree = " << (Node*)(new_trees[no_so_far]->find(node_mutation));
	 				cout << "\nAddress of the parent of selected node of the old subtree = " << (Node*)(new_trees[no_so_far]->find(node_mutation))->get_parent();
	 				// parent
	 				expr = print(pm,trees);
	 				cout << "\nParent:";
	 				cout << "\nFitness : " << trees[pm]->fitness << endl;
	 				cout << expr << endl;
	 				delete [] expr;
	 				// children
	 				expr = print(no_so_far,new_trees);
	 				cout << "\nChild:";
	 				cout << "\nson_of : " << new_trees[no_so_far]->son_of;
	 				cout << "\nParent_Fitness: " << new_trees[no_so_far]->parent_fitness << endl;
	 				cout << expr << endl;
	 				delete [] expr;
	 			}


	 			mutation_invalid = 1; //important!

	 			counter[3]++;

	 		}

	 	}

	 	// summarize results got from mutation (point or subtree)
	 	if (COMMENT) {
	 		printf("\nPopulation::new_spawn : n_mut = %i   offspring = %i", n_mut, counter[3]);
	 		printf("\nPopulation::new_spawn : MUTATION : %i individuals created (till now %.2f perc. of population)",counter[3], (double)(counter[3]*100/size));
	 	}

 }


void Population::new_spawn(RunParameters pr, int n_test_cases, int gen)  //ProblemDefinition pb
{

	int COMMENT =0;  //1 comments, 0 silent
	char *expr;

	cout << "\n\nPopulation::new_spawn" << endl;

	// set the number of individuals to copy/generate with each genetic operator
	int repr_tot = (int)(floor(repr_rate*size));   //originally pr.repr_rate
	if (repr_tot==0) repr_tot = 1;  // IMPORTANT: to ensure at least 1 individual in the elite
	int cross_tot = (int)(floor(cross_rate*size));  //originally pr.cross_rate
	int mut_tot = size - cross_tot - repr_tot;

	// allocate a new tree pointer list
   	Binary_Node **new_trees;
	new_trees = new Binary_Node *[size];
	// check
    if (!new_trees) {
       cerr << "\n\nPopulation::new_spawn : ERROR : can't allocate new_trees pointer list!\n";
       cerr << "\nExit";
       exit(-1);
    }

	
    // SORTING both trees and complete_trees first (with respect to F, not fitness=RMSE value)
    // do you really need to sort the population? This is already done in master.cpp after trees evaluation. So turned off.
    //sort(gen, tree_comp_F);

    // genetic operations
    // initialise the counter copied/generated by genetic operators
    // ([0] for copied, [1] for new, [2] for crossover, [3] for mutation)
    for (int k=0; k<4; k++) composition[k] = 0;

    //REPRODUCTION (elite) : results written to composition[0] (no. individuals copied) and composition[1] (no. individuals new)
    population_reproduction(trees, new_trees, repr_tot, true, composition);   // "true" removes the copies when performing reproduction

    //CROSSOVER : number of offspring generated written to composition[2]
    population_crossover(trees, new_trees, cross_tot, composition);

    //MUTATION: alternate point and subtree mutation. Number of individual mutated written to composition[3]
    if (mut_tot) population_mutation(trees, new_trees, mut_tot, gen, composition);

  
    // CHECK THAT NUMBER OF GENERATED/REPLICATED INDIVIDUALS = POPULATION SIZE
    int total_new_trees = 0;
    for (int k=0; k<4; k++) total_new_trees = total_new_trees + composition[k];

    cout << "\nTotal number of new trees generated/replicated : " << total_new_trees;
    cout << "\n  by reproduction : " << composition[0];
    cout << "\n  generated from scratch : " << composition[1];
    cout << "\n  by crossover : " << composition[2];
    cout << "\n  by mutation : " << composition[3];
  
  
    if (total_new_trees != size) {
    	cout << "\n\nPopulation::new_spawn : ERROR! total_new_trees ("<< total_new_trees <<") different from size (" << size << ")! Exit...\n";
    	exit(-1);
    }
 
    // PRUNING  (experimental...  still to be checked)
	// try to prune the best complete tree (remember that complete_trees[], was sorted before reproduction)
	//prune_tree(complete_trees[0]);    


	// COPY THE INDIVIDUAL WITHOUT PARAMETERS - UPDATE THE TREES ARRAY 
	//(it's not enough to make trees point to new_trees...this way deleting new_trees would delete also trees!!)
	if (COMMENT)  cout << "\n\nCOPY of the new_trees[] to trees[]" << endl;
	double v=0.;
	for (int i=0; i<size; i++) {
		// warning if the root node is not a binary node
		if (new_trees[i]->type != NODE_BINARY) {
			cout << "\nERROR! THE ROOT NODE OF new_trees[" << i << "] IS NOT BINARY!!!!" << endl;
			exit(-1);
		}
		// delete trees[i] to avoid increasing memory usage: REALLY IMPORTANT!!!
		delete trees[i]; //old_tree; 
		//
		trees[i] = (Binary_Node *)tree_copy(new_trees[i],NULL); 
		trees[i]->fitness = 999999E99;  //just to mark the new individuals...
		// in the future you may want to insert these lines in tree_copy:
		trees[i]->son_of = new_trees[i]->son_of;
		trees[i]->parent_fitness = new_trees[i]->parent_fitness;
	}

	//delete new_trees
	cout << "\nStart new_trees deletion process... ";
	int last_k=0;
	for (int k=0; k<size; k++) {
		delete new_trees[k];
		last_k = k;
	}
	cout << "Deleted up to new_trees[ " << last_k << " ] on " << size << " individuals";
	delete[] new_trees;
	
	//if (COMMENT)
		cout << "\nPopulation::new_spawn : exit" << endl;
}




// for hybrid/memetic GP evaluate() carries out more tasks than it does in classic (not hybrid) GP:
// EDITING - if needed
// 0 - copies the parameterless trees in a new array (complete_trees)
// 1 - inserts parameters in the trees (otherwise without)
// 2 - tunes newly generated complete trees, evaluating their fitness values
// 3 - evaluates fitness, adding extra penalisation terms or introducing new formulations
// 4 - sets the found fitness value to the corresponding parameterless tree
void Population::evaluate(int gen, int G)
{
	int COMMENT = 0; //1 comments, 0 silent...
	char *expr;
	Val F;				 //used to store the F aggregated value

	//if (COMMENT)
		cout << "\nPopulation :: evaluate()" << endl;

	// initialisation of genetic operations performance counters referring to current generation
	for (int i=0; i<3; i++) {
			// for current generation
			// counters of destructive/neutral/constructive genetic ops
			reproduction_perf[i]=0;
			crossover_perf[i]=0;
			s_mutation_perf[i]=0;
			p_mutation_perf[i]=0;
			// variables storing average delta of destructive/neutral/constructive genetic ops
			repr_av_delta[i]=.0;
			cross_av_delta[i]=.0;
			smut_av_delta[i]=.0;
			pmut_av_delta[i]=.0;
	}
	tot_repr= 0;
	tot_cross= 0;
	tot_smut= 0;
	tot_pmut= 0;
	


	// cycle through each individual in the population of size "size"
	cout << "Evaluating tree n. : " << endl;
	for (int i=0; i<size; i++) {

		cout << "\r" << i << " " << flush;  // \r deletes what was printed on screen the previous generation
		// print the individual without parameters
		expr = trees[i]->print();
		if (COMMENT)  {
			cout << "\n Individual  " << i << " :" << endl;
			cout << " Individual in trees[" << i << "]: " << endl;
			cout << " " << expr << endl;
			delete [] expr;
		}
		
		
		// EDITING (pag 266 PhD thesis): Edit trees to avoid undefined operations, for now only x/0
		// Doing this on the array of parameterless trees allow the modification to be transmitted throughout the evolution.
		// Strategy adopted: Ed2 - if a variable is found as divisor (right child of a division),
		// such variable is either replaced with 1 (Strategy : Ed) or might be summed with 1 (Strategy: Ed2)
		// perform_editing(i);
//			cout << "\n Individual  " << i << " :" << endl;
//			cout << " Individual in trees[" << i << "]: " << endl;
//			expr = trees[i]->print();
//			cout << " " << expr << endl;
//			delete [] expr;
		// 9/12/20 TREE EDITING to convert diverging functions to bounded functions
		// (single converging terms as variables are turned into constants)
		if (parameters->bounded==1) purge_diverging_terms(trees[i],(Node*)(trees[i]));
//			cout << "\nAfter purge_diverging_terms:";
//			cout << "\n Individual  " << i << " :" << endl;
//			cout << " Individual in trees[" << i << "]: " << endl;
//			expr = trees[i]->print();
//			cout << " " << expr << endl;
//			delete [] expr;
//			cin.get();

		// 0 - copy the parameterless trees in a new array (complete_trees)
		if (COMMENT) cout << "\n\nCopying individual trees["<< i << "] to complete_trees[" << i << "]..." << endl;
		// delete previous (old) complete tree if it exists
		if (complete_trees[i]) delete complete_trees[i];
		complete_trees[i] = (Binary_Node *)tree_copy(trees[i],NULL);
		if (COMMENT) {
			expr = complete_trees[i]->print();
			cout << " Individual in complete_trees[" << i << "]:" << endl;
			cout << " " << expr << endl;
			delete [] expr;
		}


		// 1 - insert parameters in the trees in complete_trees (otherwise without)
		if (COMMENT) cout << " Inserting parameters (constant=1) in individual complete_trees["<< i << "]" << endl;
		parameters_allocation(complete_trees[i],&complete_trees[i]);
		if (COMMENT) {
			expr = complete_trees[i]->print();
			cout << "Individual in complete_trees[" << i << "]" << endl;
			cout << " " << expr << endl;
			delete [] expr;
		}
		// find tree parameters: the tree with parameters is now built, so the number of parameters will not change until genetic operations
		complete_trees[i]->find_parameters();  //this function builds the vector storing the addresses of the parameters and their number
		//cout << "\nPopulation::evaluate : complete_trees["<<i<<"]->n_tuning_parameters = " << complete_trees[i]->n_tuning_parameters; //correct

		// 2 - tune newly generated complete trees - which means tuning all the numerical parameters in a tree
		// tuning is done on the TRAINING (BUILDING) SET and only if the number of numerical parameters is <= a given threshold
		// evaluation is done on the EVALUATION SET (they are the same set if CROSSVALIDATION is not enabled)
		// RMSE (error) and other tree's attributes are evaluated at this stage
		// (this is done also for trees copied as a result of reproduction - as the structure is copied, not the parameters!)
		if (COMMENT) cout << " Tuning parameters in individual complete_trees["<< i << "]" << endl;
		tuning_individual(parameters->n_guesses, trees[i], complete_trees[i], i);
		if (COMMENT) {
			cout << "Individual in complete_trees[" << i << "]" << endl;
			print_individual((Node *)complete_trees[i]);
			cout << " Number of parameters: " << complete_trees[i]->n_tuning_parameters;
		}
		// 22/11/20 insert here search for high level polynomials
		//search_high_level_polynomials(complete_trees[i],(Node*)(complete_trees[i]));
		//cin.get();

		// 3 - assign the found fitness value to the corresponding parameterless tree
		if (COMMENT) cout << " Updating fitness value of the corresponding parameterless tree trees["<< i << "]" << endl;
		trees[i]->fitness = complete_trees[i]->fitness;
		
		// 4 - zero-order constraint penalisation evaluation
		complete_trees[i]->pen_ord0 = constraint_evaluation(problem->data_inequality0, problem->n_inequality0, problem->constraints0,problem->v_list,parameters->nvar,complete_trees[i]);
		if (COMMENT) cout << "\ncomplete_trees[i]->pen_ord0 = " << complete_trees[i]->pen_ord0;

		// 5 - first-order constraint penalisation evaluation
		// STILL TO BE CHECKED!!! USE WITH CAUTION!!!
		// if (COMMENT) cout << " Evaluating first-order constraint penalisation evaluation" << endl;
		//complete_trees[i]->pen_ord1 = get_tree_derivative_given_norm_vector(problem, complete_trees[i]);
		complete_trees[i]->pen_ord1=0;

		// 6 - factorisation bonus : find the depth of the first factorising operation in the complete tree
		////if (problem.division)
		if (COMMENT) cout << " Finding depth of first factorising operation" << endl;
		if (parameters->w_factorisation) {
			search_first_op(complete_trees[i],(Node*)(complete_trees[i]),0); // <= IMPORTANT for FACTORISE!
		}

		// 7 - evaluate F (at last)
		// REMEMBER: during optimization fitness values (sheer RMSE) are evaluated...
		// IMPORTANT if you want to introduce penalization for complexity. Extra penalization terms must be added here to fitness
		if (COMMENT) cout << " Assigning F value to tree trees["<< i << "]" << endl;
		aggregate_F(problem, parameters, Fit_ave, complete_trees[i], gen, G);
		// assign the F value to the corresponding tree without parameters
		trees[i]->F = complete_trees[i]->F;

		// 8 - measure the performance of each genetic operator (counter and delta) and update window counter
		measure_genetic_op_performance(trees[i]);


		// final check
		if (COMMENT)  {
			expr = complete_trees[i]->print();
			cout << "complete_trees["<< i << "] : " << endl;
			cout << expr << endl;
			delete [] expr;
			cout << "Fitness: " << complete_trees[i]->fitness << " F: " << complete_trees[i]->F;
			cout << " hits: " << complete_trees[i]->hits << " tuning parameters: " << complete_trees[i]->n_tuning_parameters << endl;
		
			expr = trees[i]->print();
			cout << "trees["<< i << "] : " << endl;
			cout << expr << endl;
			delete [] expr;
			cout << "Fitness: " << trees[i]->fitness << " F: " << trees[i]->F << " hits: " << trees[i]->hits << " tuning parameters: " << trees[i]->n_tuning_parameters << endl;

			cout << "Individual tree["<< i << "] done" << endl << endl;
		}

	} //end of "for"

	// update parameters and counters referring to the whole generation
	// compute average error delta for each genetic operator
	for (int i=0; i<3; i++) {
		// for current generation
		if (reproduction_perf[i])
			repr_av_delta[i]=repr_av_delta[i]/(double)reproduction_perf[i];
		else
			repr_av_delta[i]=.0;
		if (crossover_perf[i])
			cross_av_delta[i]=cross_av_delta[i]/(double)crossover_perf[i];
		else
			cross_av_delta[i]=.0;
		if (s_mutation_perf[i])
			smut_av_delta[i]=smut_av_delta[i]/(double)s_mutation_perf[i];
		else
			smut_av_delta[i]=.0;
		if (p_mutation_perf[i])
			pmut_av_delta[i]=pmut_av_delta[i]/(double)p_mutation_perf[i];
		else
			pmut_av_delta[i]=.0;
	}
	

	// check counters for adaptive approach
	if (COMMENT)  {
		cout << "\n\nreproduction_perf[]=" << reproduction_perf[0] << ", " << reproduction_perf[1] << ", " << reproduction_perf[2] << ", " << " : tot_repr = " <<tot_repr;
		cout << "\ncrossover_perf[]=" << crossover_perf[0] << ", " << crossover_perf[1] << ", " << crossover_perf[2] << ", " << " : tot_cross = " << tot_cross;
		cout << "\ns_mutation_perf=" << s_mutation_perf[0] << ", " << s_mutation_perf[1] << ", " << s_mutation_perf[2] << ", " << " : tot_smut = " << tot_smut;
		cout << "\np_mutation_perf=" << p_mutation_perf[0] << ", " << p_mutation_perf[1] << ", " << p_mutation_perf[2] << ", " << " : tot_pmut = " << tot_pmut;
	}

	cout << "\nPopulation :: evaluate() : exit" << endl;

}

// function that evaluate RMSE (fitness), hits and R2 of complete trees belonging to the archive (repr_tot) on test data set
// 25/11/21 expanded with all the 10 objectives currently in use on building data set
void Population::evaluate_complete_trees_on_test_dataset()
{
	int COMMENT = 1;
	//-----------------------------------------------------------------
	Val result[10];
	for (int k=0; k<10; k++) result[k]=0.0;
//	result[0] = (Val)0.0;  // fitness value (error - RMSE)
//	result[1] = (Val)0.0;	// storing max absolute error
//	result[2] = (Val)0.0;	// storing n of corrections done by protected operations
//	result[3] = (Val)0.0; 	// storing value of R squared (R2)
//	result[4] = (Val)0.0;	// storing mean tree value on building data set
//	result[5] = (Val)0.0;	// storing variance of tree values on building data set
//	result[6] = (Val)0.0;	// storing min tree value on building data set
//	result[7] = (Val)0.0;	// storing max tree value on building data set
//	result[8] = (Val)0.0;	// storing inipd. variable at which the tree autocorrelation function reduces by half
//	result[9] = (Val)0.0;	// storing total variation

	int repr_tot = get_repr_tot();
	cout << "\n\nEVALUATION of the n. " << repr_tot << " first repr_tot complete trees on the TEST DATA SET";
	for (int i=0; i<repr_tot; i++) {
		// evaluate tree errors on test data set
		fitness_func(problem->Sy_test, problem->data_test, problem->n_test, complete_trees[i], result, parameters->normalised, parameters->crossvalidation);   //IMPORTANT: fitness evaluated on data_test !!!

		complete_trees[i]->n_corrections_test = (int)result[2];
		complete_trees[i]->fitness_test = result[0];  // basically RMSE, see fitness_func
		complete_trees[i]->maxabserror_test = result[1];
		complete_trees[i]->R2_test = result[3];
		complete_trees[i]->tree_mean_test = result[4];
		complete_trees[i]->tree_variance_test = result[5];
		complete_trees[i]->tree_min_test = result[6];
		complete_trees[i]->tree_max_test = result[7];
		complete_trees[i]->first_acf_root_tree_test = result[8];
		complete_trees[i]->tot_variation_tree_test = result[9];

		// evaluate F - the line below will not work on test data set features as aggregate_F uses features evaluated on building data set
		// aggregate_F(problem, parameters, Fit_ave, complete_trees[i], gen, G);

		if (COMMENT) cout << "\n complete_trees[" << i << "]->n_corrections_test = " << scientific << complete_trees[i]->n_corrections_test;
		if (COMMENT) cout << "\n complete_trees[" << i << "]->fitness_test = " << scientific << complete_trees[i]->fitness_test;
		if (COMMENT) cout << "\n complete_trees[" << i << "]->R2_test = " << scientific << complete_trees[i]->R2_test;
		if (COMMENT) cout << "\n Aggregate F not evaluated"<< endl;
	}
	//-----------------------------------------------------------------
	if (COMMENT) cout << "\n Exiting Population::evaluate_complete_trees()";
}



// fitness function definition (what is known as prediction error of the model)
// Input:
// - address of the data matrix to be used for evaluation of fitness
// - no of records or cases (no of rows to be used)
// - address of the complete tree to be evaluated (root node)
// - address of the array where to put results
// Output (stored indirectly in result_tree[] and referring to the specific dataset fed in):
// - result_tree[]
void Population::fitness_func(Val Sy, Val** data_used, int n_cases, Node *current_tree, Val *result_tree, bool normalised, bool crossvalid)
{

	// FOR THE FUTURE: change the way you import Sy: you might implement a class for storing data sets and all the related statistics
	int COMMENT =0;

	// initialisation of internal variables
	Val* treeval = new Val[n_cases];  // 10/1/21 turned Val treeval into a dynamic array
	Val a_max, a_min;
	Val error = 0.0;
	Val error_norm = 0.0;
	Val square_err = 0.0;
	Val sum_values = 0.0;
	Val sum_square_values = 0.0;
	Val SStot_tree = 0.0;
	Val result_tree_norm = 0.0;
	Val threshold = 0.0001;   // for counting hits
	Val tot_variation=0.0;
	Val maxabserror=0.0;
	int n_corrections = 0;
	((Binary_Node*)current_tree)->n_corrections=0;
	// autocorrelation function
	if (parameters->nvar==1) {
		((Binary_Node*)current_tree)->r_k_size=problem->delay_max;
		((Binary_Node*)current_tree)->r_k= new Val[problem->delay_max]; 	//array storing autocorrelation values: add delete in Destructor!!
		for (int i=0; i<(problem->delay_max); i++) {  // da spostare
			((Binary_Node*)current_tree)->r_k[i]=0.0; // da spostare
		}
	}
	Val first_acf_root=0;

	// initialization of output array
	result_tree[0] = (Val)0.0;  // fitness value (error - RMSE)
	result_tree[1] = (Val)0.0;	// maxabserror - 24/11/21 previously storing n of hits, from today storing maxabserror
	result_tree[2] = (Val)0.0;	// storing n of corrections done by protected operations
	result_tree[3] = (Val)0.0; 	// storing value of R squared (R2)
	result_tree[4] = (Val)0.0;	// storing mean tree value on building data set
	result_tree[5] = (Val)0.0;	// storing variance of tree values on building data set
	result_tree[6] = (Val)0.0;	// storing min tree value on building data set
	result_tree[7] = (Val)0.0;	// storing max tree value on building data set
	result_tree[8] = (Val)0.0;	// storing value first zero of tree autocorrelation function
	result_tree[9] = (Val)0.0;	// storing total variation

	if (COMMENT) {
		cout << "\nPopulation::fitness_func" << endl;
		cout << "Node *current_tree= " << current_tree << endl;
		cout << "Tree to be evaluated:\n" << endl;
		print_individual(current_tree);
		cout << "\nnum_vars = " << parameters->nvar << " n_cases = " << n_cases << endl;
		cout << " data_used: variables    given output    tree output" << endl; 
	}



	// cycle on the fitness cases, or points in the building/tuning data set (nfitcases)
	//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	for (int i=0; i<n_cases; i++) {
		if (COMMENT) cout << "\ni=" << i << ") ";
		
		// assign the right value to all the variables for the i-th fitness case
		for (int j=0; j<parameters->nvar; j++) {
			(problem->v_list[j])->value = data_used[i][j];
			if (COMMENT) cout << (problem->v_list[j])->value << "  ";   //var is field of Z defined at line 247 in master.cpp (it is not a member of Terminal_Var!!!)
		}

		// tree evaluation
		// (here n_corrections is evaluated) - ATTENTION if crossvalidation is enabled, corrections might be counted more than once per each tree!
		treeval[i] = tree_value((Binary_Node*)current_tree, &(((Binary_Node*)current_tree)->n_corrections));
		if (COMMENT) cout << " treeval[" << i<< "]" << treeval[i];

		// MAX and MIN of tree
		if (i==0) {
			a_min= treeval[i];
			a_max= treeval[i];
		} else {
			if (treeval[i]>a_max) a_max=treeval[i];
			if (treeval[i]<a_min) a_min=treeval[i];
		}

		// ERRORS
		// 	RMSE versions
		error =  (Val)(treeval[i] - data_used[i][parameters->nvar]);
		if (COMMENT) {
			cout << " target = " << data_used[i][parameters->nvar];
			cout << "\nerror predicted-target = " << error;
		}
		square_err = square_err + error*error;   //sum of the square of the errors - SSres (https://en.wikipedia.org/wiki/Coefficient_of_determination)
		sum_values = sum_values + treeval[i]; 				// sum of the tree values to compute mean
		sum_square_values = sum_square_values + treeval[i]*treeval[i];
		if (parameters->nvar==1) {
			if (i>0) tot_variation=tot_variation+fabs(treeval[i]-treeval[i-1]);
		}
		// maximum absolute error
		if (i==0) {
			maxabserror=fabs(error);
		} else {
			if (fabs(error)>maxabserror) maxabserror=fabs(error);
		}
		

		// normalised RMSE version
		if (normalised) { 
    		if (fabs(data_used[i][parameters->nvar])>1.0e-12)   //1.0e-12
				error_norm =  (Val)fabs((data_used[i][parameters->nvar] - treeval[i])/data_used[i][parameters->nvar]);
			else
				error_norm =  (Val)fabs(data_used[i][parameters->nvar] - treeval[i]);
    		//
			result_tree_norm = result_tree_norm + error_norm*error_norm;   //sum of the square of the errors - for RMSE *
		}
		//---------------------------------------
		

	}   // end cycle on the fitness cases
	//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

	//----------------------------------------------------------------------
	// OUTPUT
	//----------------------------------------------------------------------
	// RMSE error - result_tree[0]
	if (normalised==0) {
		if (crossvalid) {
			// root of sum of squared errors, for HyGP with crossvalidation capability
			result_tree[0] = sqrt((Val)square_err);
		} else {
			// common RMSE VERSION (for original HyGP, no crossvalidation)
			result_tree[0] = sqrt((Val)square_err/(Val)n_cases);
		}
	}
	// normalised RMSE version
	if (normalised==1) {
			result_tree[0] = sqrt(result_tree_norm/(Val)n_cases);
	}
	if (COMMENT) cout << "\nRMSE  = " << result_tree[0] << endl;

	// MAX ABSOLUTE ERROR - result_tree[1]
	result_tree[1]=maxabserror;					// before 24/11/21 it stored hits;
	if (COMMENT) cout << "Hits  = " << result_tree[1] << endl;

	// CORRECTIONS - result_tree[2]
	n_corrections=((Binary_Node*)current_tree)->n_corrections;
	result_tree[2] = n_corrections;	
	if (COMMENT) cout << "N. corrections  = " << result_tree[2] << endl;

	// R2 - result_tree[3]
	result_tree[3] = 1.0 - square_err/Sy;  // Sy (also referred to as SStot) is computed from the original, given data
	if (COMMENT) cout << "R2  = " << result_tree[3] << endl;

	// average value of the tree on the given data set - result_tree[4]
	result_tree[4] = sum_values/((Val)n_cases);
	if (COMMENT) cout << "Average value on building data set  = " << result_tree[4] << endl;

	// variance of tree values on building data set - result_tree[5] (see reformulation in http://datagenetics.com/blog/november22017/index.html)
	result_tree[5] = sum_square_values/((Val)n_cases)-result_tree[4]*result_tree[4];
	if (COMMENT) cout << "Variance on building data set  = " << result_tree[5] << endl;

	// min tree value on building data set - result_tree[6]
	result_tree[6] = a_min;
	if (COMMENT) cout << "Min value on building data set  = " << result_tree[6] << endl;

	// max tree value on building data set - result_tree[7]
	result_tree[7] = a_max;
	if (COMMENT) cout << "Max value on building data set  = " << result_tree[7] << endl;
	//--------------------------------------------------------------------------------------------------

	// AUTOCORRELATION r_k : only for 1D case and autocorrelation objective activated
	// 18/2/22 UNDER DEVELOPMENT: make it more efficient? So far 2x n_cases cycles are needed to complete the calculations...
	if ((parameters->nvar==1) && (parameters->w_ACF>0.0)) {
		// first compute SStot_tree = sum((y_tree- average(y_tree))^2) (https://en.wikipedia.org/wiki/Coefficient_of_determination)
		// good! Direct formula for SStot:
		SStot_tree = sum_square_values-sum_values*sum_values/(Val)n_cases;
		// then compute autocorrelation function r_k
		int ACF_first_root_found=0;
		int delay_max=(int)floor(n_cases/2);
		for (int delay=0; delay<delay_max; delay++) {
			for (int k=0; k < n_cases-delay; k++) {
				((Binary_Node*)current_tree)->r_k[delay] = ((Binary_Node*)current_tree)->r_k[delay] + (treeval[k]-result_tree[4])*(treeval[k+delay]-result_tree[4]);
			}
			((Binary_Node*)current_tree)->r_k[delay]= (((Binary_Node*)current_tree)->r_k[delay])/SStot_tree;


			// 20/2/21 search for first point at which ACF halves -------------- no longer first root of autocorrelation function (the one closest to 0)
			// if r_k[delay-1] r_k[delay]>0 are both positive or negative do nothing: r_k[delay-1]*r_k[delay]>0
			if (ACF_first_root_found==0) {
				if ( fabs(((Binary_Node*)current_tree)->r_k[delay]-0.5)<1.0E-12 ) {
					first_acf_root=data_used[delay][0];  //mind! Only for n_var=1!
					ACF_first_root_found=1;
				} else {
					if ((((Binary_Node*)current_tree)->r_k[delay-1]-0.5>0) && (((Binary_Node*)current_tree)->r_k[delay]-0.5<0)) {
						first_acf_root=(data_used[delay-1][0]+data_used[delay][0])/2.0; //mind! Only for n_var=1!
						ACF_first_root_found=1;
					}
				}
			}

		} // end delay cycle
	} // end AUTOCORRELATION CALCULATION


	// just to check autocorrelation function values
//	cout << "\nPopulation::fitness_func : Computed autocorrelation function:" << endl;
//	for (int delay=0; delay<problem->delay_max; delay++) {
//		cout << "r_k[" << delay << "]= " << ((Binary_Node*)current_tree)->r_k[delay] << endl;
//	}
//	cin.get();

	// first zero of autocorrelation function - result_tree[8]
	result_tree[8] = first_acf_root;

	// total variation - result_tree[9]
	result_tree[9] = tot_variation;

	delete[] treeval;  // why not storing the values in each root node, so when you have to save them to file in Reporter you do not have to recompute them?
}



// Function that defines the objective to minimise for PSO algorithm (coarse tree's parameters optimisation)
// Input:
// - array of tree's parameters
// - number of parameters (spacedim)
// - pointer to tree root node
// Output:
// - tree RMSE error (simple fitness, not aggregated one!)
double Population::pso_objfunction(double* x, int n_param, Binary_Node *ntree)
{
	int COMMENT=0;

	if (COMMENT) cout << "\nPopulation::pso_objfunction()" << endl;

	double fitness; // actually an error (RMSE). Be aware that this should be a Val type...
	Val result[10];
	result[0] = (Val)0.0;  // fitness value (error - RMSE)
	result[1] = (Val)0.0;	// storing max absolute error
	result[2] = (Val)0.0;	// storing n of corrections done by protected operations
	result[3] = (Val)0.0; 	// storing value of R squared (R2)
	result[4] = (Val)0.0;	// storing mean tree value on building data set
	result[5] = (Val)0.0;	// storing variance of tree values on building data set
	result[6] = (Val)0.0;	// storing min tree value on building data set
	result[7] = (Val)0.0;	// storing max tree value on building data set
	result[8] = (Val)0.0;	// storing value first zero of tree autocorrelation function
	result[9] = (Val)0.0;	// storing total variation

	// update tree parameters with input values x
	if (COMMENT) {
		cout << "n_param= " << n_param << endl;
		cout << "x is correctly initialised/passed to Population::pso_objfunction()?" << endl;
		for (int i=0; i< n_param; i++) {
			cout << "x[" << i << "]=" << x[i] << "  ";
		}
	}
	update_complete_tree(ntree, x, n_param);

	// evaluate tree fitness (RMSE error)
	fitness_func(problem->Sy, problem->data_validation, problem->n_validation, ntree, result, parameters->normalised, parameters->crossvalidation);   //IMPORTANT: fitness evaluated on data_validation !!!
	fitness = (double)result[0];
	//hits = (int)result[1];
	//n_corrections = (int)result[2];
	//R2 = result[3];

	return fitness;
}


double Population::constraint_evaluation(Val**data, int n_cases, char* constraints, Variable **var_list, int n_vars, Binary_Node *complete_tree)
{
	int COMMENT = 0;
	double leak = 0.0;
	Val value = 0.0;

	for (int i=0; i< n_cases; i++) {
		// assign the right value to all the variables for the i-th fitness case

		if (COMMENT) cout << "\n\nConstraint evaluation point:";
		for (int j=0; j< n_vars; j++) {
			(*var_list[j]).value = data[i][j];
			if (COMMENT)
				cout << (*var_list[j]).value << "  ";
		}
		// evaluate tree at given point
		value =  tree_value(complete_tree, NULL);

		// constraint evaluation
		if (constraints[i]=='>') {
			if (value <= data[i][n_vars])
				// 1st approach: just count the unfeasible points
				//leak++;
				// 2nd approach: measure the distance from the feasible region (dist)
				leak = leak + fabs(data[i][n_vars] - value);
		}
		if (constraints[i]=='<') {
			if (value >= data[i][n_vars])
				// 1st approach: just count the unfeasible points
				//leak++;
				// 2nd approach: measure the distance from the feasible region (dist)
				leak = leak + fabs(data[i][n_vars] - value);
		}

		if (COMMENT) {
			cout << "\n\ntree value = " << value;
			cout << "\n constraint : " << constraints[i];
			cout << "\n constraint value = " << data[i][n_vars];
			cout << "\nleak = " << leak << endl;
		}
	}

	if (COMMENT) cout << "\nleak total = " << leak;

	return leak;
}


// function to adjust adaptively the weights of aggregate fitness function F
// 5/12/21 currently applied only to 3 objectives: Fc_perc_ave[1], Fc_perc_ave[8], Fc_perc_ave[10]
void Population::adjust_Fweight(int gen) {
	int COMMENT = 1; //1 comments, 0 silent...

	Fweight[9]= 0.0;

//	// ADAPTIVE STRATEGY A1 - first idea
//	// objective 1: main objective, RMSE value Fweight[1]
//	if ((Fc_perc_ave[1]>5.0) && (Fweight[8]>0.05) && (Fweight[10]>0.05)) {
//		Fweight[8]= Fweight[8]-0.05;		// statistical properties
//		Fweight[10]= Fweight[10]-0.05;	// difference in point at which ACF halves
//		// The sum of all a_i coefficients must be 1!!
//		Fweight[1]=double(1.-Fweight[2]-Fweight[3]-Fweight[4]-Fweight[5]-Fweight[6]-Fweight[7]-Fweight[8]-Fweight[9]-Fweight[10]-Fweight[11]);
//		return;
//	}
//
//	// objective 2: penalisation on statistical properties of the tree (average and variance) Fweight[8]
//	if ((Fc_perc_ave[8]>5.0) && (Fweight[10]>0.05) && (Fweight[1]>0.05)) {
//		Fweight[8]= Fweight[8]+0.1;		// statistical properties
//		Fweight[10]= Fweight[10]-0.05;	// difference in point at which ACF halves
//		// The sum of all a_i coefficients must be 1!!
//		Fweight[1]=double(1.-Fweight[2]-Fweight[3]-Fweight[4]-Fweight[5]-Fweight[6]-Fweight[7]-Fweight[8]-Fweight[9]-Fweight[10]-Fweight[11]);
//		return;
//	}
//
//	// objective 3: penalisation linked to autocorrelation function (point at which ACF halves) Fweight[10]
//	if ((Fc_perc_ave[10]>5.0) && (Fweight[1]>0.05) && (Fweight[8]>0.05)) {
//		Fweight[8]= Fweight[8]-0.05;		// statistical properties
//		Fweight[10]= Fweight[10]+0.1;	// difference in point at which ACF halves
//		// The sum of all a_i coefficients must be 1!!
//		Fweight[1]=double(1.-Fweight[2]-Fweight[3]-Fweight[4]-Fweight[5]-Fweight[6]-Fweight[7]-Fweight[8]-Fweight[9]-Fweight[10]-Fweight[11]);
//		return;
//	}
//	//  Fweight 2, 3, 4, 5, 6, 7 are constant and equal to the values in the input file


//	// ADAPTIVE STRATEGY A2 - MINMAX
//	// objective 1: main objective, RMSE value Fweight[1]
//	if ((Fc_perc_ave[1]>Fc_perc_ave[8]) && (Fc_perc_ave[1]>Fc_perc_ave[10])) {
//		Fweight[1]=0.8;		// RMSE on building data set
//		Fweight[8]= 0.1;		// statistical properties
//		Fweight[10]= 0.0;	// difference in point at which ACF halves
//		// The sum of all a_i coefficients must be 1!!
//		return;
//	}
//
//	// objective 2: penalisation on statistical properties of the tree (average and variance) Fweight[8]
//	if ((Fc_perc_ave[8]>Fc_perc_ave[10]) && (Fc_perc_ave[8]>Fc_perc_ave[1])) {
//		Fweight[1]=0.1;
//		Fweight[8]= 0.8;		// statistical properties
//		Fweight[10]= 0.0;	// difference in point at which ACF halves
//		// The sum of all a_i coefficients must be 1!!
//		return;
//	}
//
//	// objective 3: penalisation linked to autocorrelation function (point at which ACF halves) Fweight[10]
//	if ((Fc_perc_ave[10]>Fc_perc_ave[1]) && (Fc_perc_ave[10]>Fc_perc_ave[8])) {
//		Fweight[1]=0.0;
//		Fweight[8]= 0.0;		// statistical properties
//		Fweight[10]= 1.0;	// difference in point at which ACF halves
//		// The sum of all a_i coefficients must be 1!!
//		return;
//	}
//  Fweight 2, 3, 4, 5, 6, 7 are constant and equal to the values in the input file

	// ADAPTIVE STRATEGY A3 - Sergey's adaptation of conjugate gradients
	// find largest objective in absolute value and update only the corresponding weight
	double Fc_ave_abs[12];
	int index_max_abs=0;
	for (int i=1; i<12; i++) Fc_ave_abs[i]=fabs(Fc_ave[i]);
	index_max_abs=distance(Fc_ave_abs, max_element(Fc_ave_abs, Fc_ave_abs + 12));  // defines the largest objective in absolute value

	double sum=0.0;
	for (int i=1; i<12; i++) {
		if (i!=index_max_abs) {
			sum=sum+Fweight[i]*Fc_ave[i];
		}
	}

	if (COMMENT) {
		for (int i=1; i<12; i++) cout << "\nFweight[ " << i << "] = " << Fweight[i] << "  Fc_ave[ " << i << "] = " << Fc_ave[i];
		cout << "\nMax abs value for i_max=" << index_max_abs;

	}

	Fweight[index_max_abs]=-sum/Fc_ave[index_max_abs]; // does it make sense to have negative weights!!?

	if (COMMENT) {
		cout << "\nNew weight for i_max: Fweight[ " << index_max_abs << "] = " << Fweight[index_max_abs];
	}

}


// function to compute the aggregate version of fitness (called F)
// input: number of the tree in trees[]  (and complete_trees[])
// output: value of the aggregate function F
void Population::aggregate_F(ProblemDefinition* ppd, RunParameters* pr, Val average_err, Binary_Node *complete_tree, int gen, int G)
{	

	// not correct value of n_tuning_parameters, shown 999999
	//cout << "\nPopulation::aggregate_F : complete_tree->n_tuning_parameters = " << complete_tree->n_tuning_parameters;

	int COMMENT = 0; //1 comments, 0 silent...
	char* expr;
	double f_5=0.0;

	//-------------------------------------------------------------
	// first objective: FITNESS (RMSE)
	//-------------------------------------------------------------
	//normalization: the fitness value (RMSE here) is divided by the average fitness in the last population
	//double nt = (double)n_test_cases;
	//F1 = (double)((complete_tree->fitness)*(complete_trees[i]->fitness)*nt/sum_output);   //Alvarez's version - too SMALL!
	//F1 = 1.-1./(1.+complete_tree->fitness); 
	//F1 = complete_tree->fitness/sum_output;	 //the term is constant, but does not give great results...
	if (pr->crossvalidation==0)
		// F1 = (complete_tree->fitness/average_err);	//used in the PhD Thesis: average_err is the av. fitness in the archive of the previous gen. THE TERM VARIES DURING THE EVOLUTION!!!
		F[1] = 100.0*(complete_tree->fitness/average_err);	//used for Mode 1 Ocean problem
		// F1 = pow(complete_tree->fitness,2);	// test 1/2/2021
	else
		F[1] = complete_tree->fitness;

	//---------------------------------------------------------------------
	// second objective: NUMBER OF TUNING PARAMETERS
	//---------------------------------------------------------------------
	// normalization: not necessary if the objective is already a pure number
	//F2 = (double)(pow(v,2.)/d_lim); //n of tuning parameters as a measure of complexity! 
	//F2 = complete_tree->calc_depth()/d_lim;  //depth as a measure of complexity (already pure number, dividing by d_lim was not necessary)!
	//F2 = (double)(complete_tree->n_corrections);
	//-------------------------------------------------------------
	// third objective: No of corrections (singularities)
	//-------------------------------------------------------------
	//-------------------------------------------------------------
	// fourth objective : SIZE
	//-------------------------------------------------------------
	//-------------------------------------------------------------
	// fifth objective : PENALISATION INEQUALITY CONSTRAINTS ORDER 0 (values)
	//-------------------------------------------------------------
	//F5 = complete_tree->pen_ord0;

	//-------------------------------------------------------------
	// sixth objective : PENALISATION INEQUALITY CONSTRAINTS ORDER 1 (derivatives)
	//-------------------------------------------------------------
	// megaexp1
	// F6 = (1000000.0)*(exp(complete_tree->pen_ord1)-1.0);
	// megaexp2
	// F6 = 1000000.0*exp(pow(complete_tree->pen_ord1,2.0))-1.0);
	// divave
	//F6 = (complete_tree->pen_ord1)/(pen_ord1_ave+1.0);//complete_tree->pen_ord1;
	// kdivave
	//F6 = 1000.0*(complete_tree->pen_ord1)/(pen_ord1_ave+1.0);
	// divave001
	F[6] = (complete_tree->pen_ord1)/(pen_ord1_ave+.001);
	// expdivave
	//F6 = exp((complete_tree->pen_ord1)/(pen_ord1_ave+1.0))-1.0;
	// expdivave2
	//F6 = exp(pow((complete_tree->pen_ord1)/(pen_ord1_ave+1.0),2.0))-1.0;

	//-------------------------------------------------------------
	// seventh objective : FACTORISATION (related to depth of first division)
	//-------------------------------------------------------------
	// factorisation bonus enabled only if w_factorisation > 0 (see input file)
	F[7] = 0.0;
	if (pr->w_factorisation>0) {
		// FACTORISE APPROACH (also called FACTORISATION BONUS)
		double d = (double)(complete_tree->depth_first_op);
		// OLD STUFF not tested enough //if (d == -1)
		// OLD STUFF not tested enough //	F7 = 10.0*complete_tree->count();//exp(2.0)-1.0; //  //0.0;
		// OLD STUFF not tested enough //else
		// OLD STUFF not tested enough //	F7 = 10.0*pow(d-1.0,2.0);   //1) exp(d)-1.0; //2) exp(d-1.0)-1.0; //3)10.0*d*d
		// "BONUS on ERROR": factorise approach (factorisation bonus)
		if (( d < 0.2*(double)(complete_tree->calc_depth()) ) && (d <= 5.0))
			//if ( d < 0.05*((double)(complete_tree->calc_depth())*(1.0 - F4)) )  // depth to size approach
			F[7]= 0.1;
		else
			F[7] = 1.0;
	}

	//-------------------------------------------------------------------------
	// eighth objective : statistical properties of the tree (mean and variance)
	//-------------------------------------------------------------------------
	// first attempt only valid for target function with zero mean and similar variance to input data
	//F8 = fabs(complete_tree->tree_mean) // Strategy 1
	//F8 = sqrt(fabs(ppd->y_var - complete_tree->tree_variance)); // Strategy 2
	//F8 = fabs(complete_tree->tree_mean) + sqrt(fabs(ppd->y_var - complete_tree->tree_variance)); // Strategy 3
	//F8 = sqrt(fabs(ppd->y_var - complete_tree->tree_variance))/(fabs(1+fabs(ppd->y_ave-complete_tree->tree_mean)));// Strategy 4
	//F8 = fabs(ppd->y_var - complete_tree->tree_variance)/pow(fabs(1+fabs(ppd->y_ave-complete_tree->tree_mean)),2); // Strategy 5
	// strategies below still under testing

	//F8 = pow(sqrt(fabs(ppd->y_var-complete_tree->tree_variance)),1)+pow(fabs(ppd->y_ave-complete_tree->tree_mean),1); //Strategy 8
	//(exp((fabs(ppd->y_ave-complete_tree->tree_mean)))-1.0);

	//-------------------------------------------------------------------------
	// 9th objective : presence of high level polynomials that cause divergent behaviour (extrapolation issues)
	//-------------------------------------------------------------------------
	//F9 = 1.0*(complete_tree->diverging);  //1.0*(complete_tree->diverging);

	//-------------------------------------------------------------------------
	// 10th objective : difference in autocorrelation function properties between input signal and model
	//-------------------------------------------------------------------------

	//-------------------------------------------------------------------------
	// 11th objective : difference in total variation
	//-------------------------------------------------------------------------
	//F11=pow(sqrt(fabs(ppd->tot_variation_input-complete_tree->tot_variation_tree)),3);

	// STRATEGY 4
	if ((pr->strat_statp)==4) {
		// fitness value (RMSE)
		F[1] = complete_tree->fitness/average_err;  //MIND that this expression of error does not allow to reduce monotonically the error by the generation (it is not absoluite, the average error reduces by generation...)
		// number of tuning parameters
		F[2] = pow((double)(complete_tree->n_tuning_parameters), 2.0);
		//No of corrections performed by protected operations
		F[3] = (double)(complete_tree->n_corrections)*1.0E6;
		//SIZE
		F[4] = (double)(complete_tree->count());
		f_5=pow(complete_tree->pen_ord0, 3.0);
		F[5]=1.0E6*(exp(f_5)-1.0);   // F[5] = 1.0E6*(exp(F[5]*F[5]*F[5])-1.0);
		F[6] = (complete_tree->pen_ord1)/(pen_ord1_ave+.001);
		F[7]=0.0;
		F[8] = sqrt(fabs(ppd->y_var-complete_tree->tree_variance))/(1.0+fabs(ppd->y_ave-complete_tree->tree_mean));
		F[9] = 0.0;
		F[10] = 0.0;
		F[11] = 0.0;
	}

	// STRATEGY 6
	if ((pr->strat_statp)==6) {
		// fitness value (RMSE)
		F[1] = complete_tree->fitness/average_err;
		// number of tuning parameters
		F[2] = pow((double)(complete_tree->n_tuning_parameters), 2.0);
		//No of corrections performed by protected operations
		F[3] = (double)(complete_tree->n_corrections)*1.0E6;
		//SIZE
		F[4] = (double)(complete_tree->count());
		f_5=pow(complete_tree->pen_ord0, 3.0);
		F[5]=1.0E6*(exp(f_5)-1.0);   // F[5] = 1.0E6*(exp(F[5]*F[5]*F[5])-1.0);
		F[7]=0.0;
		F[8] = pow(sqrt(fabs(ppd->y_var-complete_tree->tree_variance)),3)+pow(fabs(ppd->y_ave-complete_tree->tree_mean),3);
		F[9] = 0.0;
		F[10] = 0.0;
		F[11] = 0.0;
	}

	// STRATEGY 7
	if ((pr->strat_statp)==7) {
		// fitness value (RMSE)
		F[1] = complete_tree->fitness/average_err;
		// number of tuning parameters
		F[2] = pow((double)(complete_tree->n_tuning_parameters), 2.0);
		//No of corrections performed by protected operations
		F[3] = (double)(complete_tree->n_corrections)*1.0E6;
		//SIZE
		F[4] = (double)(complete_tree->count());
		f_5=pow(complete_tree->pen_ord0, 3.0);
		F[5]=1.0E6*(exp(f_5)-1.0);   // F[5] = 1.0E6*(exp(F[5]*F[5]*F[5])-1.0);
		F[7]=0.0;
		F[8] = pow(sqrt(fabs(ppd->y_var-complete_tree->tree_variance)),2)+pow(fabs(ppd->y_ave-complete_tree->tree_mean),2);//Strategy 7
		F[9] = 0.0;
		F[10] = 0.0;
		F[11] = 0.0;
	}

	//STRATEGY 8
	if ((pr->strat_statp)==8) {
		// fitness value (RMSE)
		F[1] = complete_tree->fitness/average_err;
		// number of tuning parameters
		F[2] = pow((double)(complete_tree->n_tuning_parameters), 2.0);
		//No of corrections performed by protected operations
		F[3] = (double)(complete_tree->n_corrections)*1.0E6;
		//SIZE
		F[4] = (double)(complete_tree->count());
		f_5=pow(complete_tree->pen_ord0, 3.0);
		F[5]=1.0E6*(exp(f_5)-1.0);   // F[5] = 1.0E6*(exp(F[5]*F[5]*F[5])-1.0);
		F[7]=0.0;
		F[8] = pow(sqrt(fabs(ppd->y_var-complete_tree->tree_variance)),3)+pow(fabs(ppd->y_ave-complete_tree->tree_mean),3)+pow(fabs(ppd->y_max-complete_tree->tree_max),3)+pow(fabs(ppd->y_min-complete_tree->tree_min),3);//Strategy 8
		F[9] = 0.0;
		F[10] = 0.0;
		F[11] = 0.0;
	}

	//STRATEGY 11
	if ((pr->strat_statp)==11) {
		//
		// all the other objectives?!
		//
		f_5=pow(complete_tree->pen_ord0, 3.0);
		F[5]=1.0E6*(exp(f_5)-1.0);   // F[5] = 1.0E6*(exp(F[5]*F[5]*F[5])-1.0);
		F[7]=0.0;
		F[8] = pow(sqrt(fabs(ppd->y_var-complete_tree->tree_variance)),3)+pow(fabs(ppd->y_ave-complete_tree->tree_mean),3);
		//F8 = pow(fabs(ppd->y_max-complete_tree->tree_max),3)+pow(fabs(ppd->y_min-complete_tree->tree_min),3);
		// test to monitor ratio between difference in peaks and standard deviation
		//F8= pow(fabs(ppd->y_max-complete_tree->tree_max)/(sqrt(ppd->y_var)),3)+pow(fabs(ppd->y_min-complete_tree->tree_min)/(sqrt(ppd->y_var)),3);
		// test to monitor ratio between difference in peaks and target function amplitude
		//F8= pow(100*fabs(ppd->y_max-complete_tree->tree_max)/(fabs(ppd->y_max-ppd->y_min)),3)+pow(100*fabs(ppd->y_min-complete_tree->tree_min)/(fabs(ppd->y_max-ppd->y_min)),3);
	}

	// STRATEGY 13 - all objectives
	if ((pr->strat_statp)==13) {
		// fitness value (RMSE)
		F[1] = exp(10.0*complete_tree->fitness/fabs(ppd->y_max-ppd->y_min)); //20/2/21 try exp(complete_tree->fitness/average_err);
		// number of tuning parameters squared
		F[2] = (complete_tree->n_tuning_parameters);  //0.0001*
		// F[2] = pow((double)(complete_tree->n_tuning_parameters), 2.0); // 12/12/21 currently in use, but very large
		// No of corrections performed by protected operations
		F[3] = (double)(complete_tree->n_corrections)*1.0E6;
		// SIZE (number of nodes)
		F[4] = (double)(complete_tree->count());
		//F4 = pow((double)(complete_tree->count()),2.0);
		// factorisation bonus enabled only if w_factorisation > 0 (see input file)
		f_5=pow(complete_tree->pen_ord0, 3.0);
		F[5]=1.0E6*(exp(f_5)-1.0);   // F[5] = 1.0E6*(exp(F[5]*F[5]*F[5])-1.0);
		F[6] = 0.0;     //(complete_tree->pen_ord1)/(pen_ord1_ave+.001);
		F[7] = 0.0;
		// statistical properties of the tree (mean and variance)
		//F[8] = exp(10.0*fabs(ppd->y_var-complete_tree->tree_variance)/ppd->y_var) + exp(10.0*fabs(ppd->y_ave-complete_tree->tree_mean)/(fabs(ppd->y_ave)+1.0)) + exp(10.0*fabs(ppd->y_max-complete_tree->tree_max)/fabs(ppd->y_max));
		F[8] = exp(10.0*fabs(ppd->y_var-complete_tree->tree_variance)/ppd->y_var) + exp(10.0*fabs(ppd->y_ave-complete_tree->tree_mean)/(fabs(ppd->y_ave)+1.0)) + exp(10.0*fabs(ppd->y_max-complete_tree->tree_max)/fabs(ppd->y_max))+ exp(10.0*fabs(ppd->y_min-complete_tree->tree_min)/fabs(ppd->y_min));

		// presence of high level polynomials that cause divergent behaviour
		F[9] = 0.0;
		// difference in point at which ACF halves
		//F10 = pow(sqrt(fabs(ppd->first_acf_root_input-complete_tree->first_acf_root_tree)),3);
		F[10] = exp(10.0*fabs(ppd->first_acf_root_input-complete_tree->first_acf_root_tree)/ppd->first_acf_root_input)-1.0;
		// difference in total variation
		F[11] = pow(fabs(ppd->tot_variation_input-complete_tree->tot_variation_tree)/ppd->tot_variation_input,3);
	}

	// STRATEGY 14 - RMSE is not used as objective
	if ((pr->strat_statp)==14) {
		// fitness value (RMSE)
		F[1] = 0.0;  // the aim of strategy 14 is to evolve a model with the same statistical properties of the original signal, not local accuracy is not a priority
		// number of tuning parameters
		F[2] = pow((double)(complete_tree->n_tuning_parameters), 2.0);
		//No of corrections performed by protected operations
		F[3] = (double)(complete_tree->n_corrections)*1.0E6;
		//SIZE
		F[4] = (double)(complete_tree->count());
		//F4 = pow((double)(complete_tree->count()),2.0);
		f_5=pow(complete_tree->pen_ord0, 3.0);
		F[5]=1.0E6*(exp(f_5)-1.0);   // F[5] = 1.0E6*(exp(F[5]*F[5]*F[5])-1.0);
		F[6] = (complete_tree->pen_ord1)/(pen_ord1_ave+.001);
		// factorisation bonus enabled only if w_factorisation > 0 (see input file)
		F[7] = 0.0;
		// statistical properties of the tree (mean and variance). 20/2/22 added exp(10.0* to "max" objective
		F[8] = exp(10.0*fabs(ppd->y_var-complete_tree->tree_variance)/ppd->y_var) + exp(10.0*fabs(ppd->y_ave-complete_tree->tree_mean)/(fabs(ppd->y_ave)+1)) + exp(10.0*fabs(ppd->y_max-complete_tree->tree_max)/fabs(ppd->y_max)); //+pow(fabs(ppd->y_ave-complete_tree->tree_mean)/fabs(ppd->y_max-ppd->y_min),3)
		//F[8] = fabs(ppd->y_max-complete_tree->tree_max)/fabs(ppd->y_max); //+pow(fabs(ppd->y_ave-complete_tree->tree_mean)/fabs(ppd->y_max-ppd->y_min),3)

		// presence of high level polynomials that cause divergent behaviour
		F[9] = 0.0;
		// difference in point at which ACF halves
		//F10 = pow(sqrt(fabs(ppd->first_acf_root_input-complete_tree->first_acf_root_tree)),3);
		F[10] = exp(10.0*fabs(ppd->first_acf_root_input-complete_tree->first_acf_root_tree)/ppd->first_acf_root_input)-1.0;
		// difference in total variation
		F[11] = pow(fabs(ppd->tot_variation_input-complete_tree->tot_variation_tree)/ppd->tot_variation_input,3);
	}

	// STRATEGY 15 : under test with SPOD data set
	if ((pr->strat_statp)==15) {
		// fitness value (RMSE)
		F[1] = exp(complete_tree->fitness/fabs(ppd->y_max-ppd->y_min)); //20/2/21 try exp(complete_tree->fitness/average_err);
		// complexity - number of tuning parameters
		F[2] = pow((double)(complete_tree->n_tuning_parameters), 1.0);
		// No of corrections performed by protected operations
		F[3] = (double)(complete_tree->n_corrections)*1.0E6;
		// SIZE
		F[4] = (double)(complete_tree->count());
		//F4 = pow((double)(complete_tree->count()),2.0);
		// penalisation 0 order
		f_5=pow(complete_tree->pen_ord0, 3.0);
		F[5]=1.0E6*(exp(f_5)-1.0);   // F[5] = 1.0E6*(exp(F[5]*F[5]*F[5])-1.0);
		// penalisation 1st order
		F[6] = (complete_tree->pen_ord1)/(pen_ord1_ave+.001);
		// factorisation bonus enabled only if w_factorisation > 0 (see input file)
		F[7] = 0.0;
		// statistical properties of the tree (mean, variance, max)
		//F8 = exp(fabs(ppd->y_var-complete_tree->tree_variance)/ppd->y_var) + exp(10.0*fabs(ppd->y_ave-complete_tree->tree_mean)/(fabs(ppd->y_ave)+1)) + fabs(ppd->y_max-complete_tree->tree_max)/fabs(ppd->y_max);
		F[8] = pow(sqrt(fabs(ppd->y_var-complete_tree->tree_variance)/ppd->y_var),2)+pow(fabs(ppd->y_ave-complete_tree->tree_mean)/(fabs(ppd->y_ave)+1),2);//Strategy 7
		// presence of high level polynomials that cause divergent behaviour
		F[9] = 0.0;
		// difference in point at which ACF halves
		//F10 = exp(fabs(ppd->first_acf_root_input-complete_tree->first_acf_root_tree)/ppd->first_acf_root_input)-1.0;
		F[10] = pow(sqrt(fabs(ppd->first_acf_root_input-complete_tree->first_acf_root_tree)/ppd->first_acf_root_input),2); // why sqrt? puoi toglierla!
		// difference in total variation
		F[11] = 0.0; //pow(fabs(ppd->tot_variation_input-complete_tree->tot_variation_tree)/ppd->tot_variation_input,3);
	}

	// STRATEGY 16- same objectives as in Strategy 13 but normalised with respect to sum of all objectives
	if ((pr->strat_statp)==16) {
		// fitness value (RMSE)
		F[1] = exp(10.0*complete_tree->fitness/fabs(ppd->y_max-ppd->y_min)); //20/2/21 try exp(complete_tree->fitness/average_err);
		// number of tuning parameters squared
		F[2] = 0.0001*(Val)(complete_tree->n_tuning_parameters);
		// F[2] = pow((double)(complete_tree->n_tuning_parameters), 2.0); // 12/12/21 currently in use, but very large
		// No of corrections performed by protected operations
		F[3] = (double)(complete_tree->n_corrections)*1.0E6;
		// SIZE (number of nodes)
		F[4] = (double)(complete_tree->count());
		//F4 = pow((double)(complete_tree->count()),2.0);
		// factorisation bonus enabled only if w_factorisation > 0 (see input file)
		f_5=pow(complete_tree->pen_ord0, 3.0);
		F[5]=1.0E6*(exp(f_5)-1.0);   // F[5] = 1.0E6*(exp(F[5]*F[5]*F[5])-1.0);
		F[6] = 0.0;     //(complete_tree->pen_ord1)/(pen_ord1_ave+.001);
		F[7] = 0.0;
		// statistical properties of the tree (mean and variance)
		F[8] = exp(10.0*fabs(ppd->y_var-complete_tree->tree_variance)/ppd->y_var) + exp(10.0*fabs(ppd->y_ave-complete_tree->tree_mean)/(fabs(ppd->y_ave)+1.0)) + fabs(ppd->y_max-complete_tree->tree_max)/fabs(ppd->y_max); //+pow(fabs(ppd->y_ave-complete_tree->tree_mean)/fabs(ppd->y_max-ppd->y_min),3)
		// presence of high level polynomials that cause divergent behaviour
		F[9] = 0.0;
		// difference in point at which ACF halves
		//F10 = pow(sqrt(fabs(ppd->first_acf_root_input-complete_tree->first_acf_root_tree)),3);
		F[10] = exp(10.0*fabs(ppd->first_acf_root_input-complete_tree->first_acf_root_tree)/ppd->first_acf_root_input)-1.0;
		// difference in total variation
		F[11] = pow(fabs(ppd->tot_variation_input-complete_tree->tot_variation_tree)/ppd->tot_variation_input,3);

		double Fsum=0.0;
		for (int i=1; i<12; i++) Fsum=Fsum+F[i];
		for (int i=1; i<12; i++) F[i]=F[i]/Fsum;
	}


	//------------------------------------------------------------
	// copy locally stored objectives values to tree attributed
	//------------------------------------------------------------
	for (int i=1; i<12; i++) complete_tree->Fc[i] = F[i];  // from 1 to 11
	// complete_tree->Fc[5] = F[5]; //megadistexp3
	// complete_tree->Fc[7] = F[7]; // not used in the standard approach


	//------------------------------------------------------------
	// calculation of contributions (objectives multiplied by weights) at tree level
	//------------------------------------------------------------
	for (int i=1; i<12; i++) complete_tree->T[i] = Fweight[i]*complete_tree->Fc[i];  // from 1 to 11

	//------------------------------------------------------------
	// fitness function definition
	//------------------------------------------------------------
	// dynamic coefficient
	// double dc;
	// if (G) 
	//	dc = (double)(gen/G);
	// else dc = 1.;

	// MINMAX approach - Strategy 11
	if (pr->strat_statp==11) {	// if (pr->minmax) {
		// build the list of elements among which you want to find the maximum
		double list[5];
		list[0] = complete_tree->T[1];
		list[1] = complete_tree->T[2];
		list[2] = complete_tree->T[3];
		list[3] = complete_tree->T[4];
		list[4] = complete_tree->T[8];

		complete_tree->F = *max_element(list, list+5);
	}

	// AGGREGATED (WEIGHTED SUM) approach - not strategy 11
	if (pr->strat_statp!=11) {
		// compute F as linear combination of fitness value and other objectives
		//F = a1*F1+a2*dc*dc*F2+a3*dc*dc*F3+a4*dc*dc*F4;   //dynamic weights
		//F = a1*F1+a2*F2+a3*(exp(F3*F3)-1)*1000000.0+a4*F4+(1000000.0)*(exp(F5*F5)-1.0)*a5;
		//F = a1*F1+a2*F2+a3*F3*1000000.0+a4*F4+(1000000.0)*(exp(F5*F5)-1.0)*a5;
		
		if (pr->w_factorisation<=0) {
			// STANDARD APPROACH (HyGP)
			complete_tree->F = 0.0;
			for (int i=1; i<12; i++) complete_tree->F = complete_tree->F + complete_tree->T[i];  // 11 objectives, from 1 to 11
			if (isinf(complete_tree->F)) complete_tree->F=MAX_F_VALUE;    // see Population.h
		}

		// adaptive approach
		//	complete_tree->F = complete_tree->T1 + (-1.0*((double)learning_on-1.0))*complete_tree->T2 + complete_tree->T3 + ((double)learning_on+1.0)*complete_tree->T4 + complete_tree->T5 + complete_tree->T6 + complete_tree->T7;

		if (pr->w_factorisation>0) {
			// FACTORISE APPROACH (also called FACTORISATION BONUS)
			// mind that "search_first_op" has to be enabled for factorisation bonus to work!
			complete_tree->F = F[7]*(complete_tree->T[1] + complete_tree->T[2] + complete_tree->T[3] + complete_tree->T[4] + complete_tree->T[5] + complete_tree->T[6]);
			if (isinf(complete_tree->F)) complete_tree->F=MAX_F_VALUE;    // see Population.h
		}
		// factorise2
		//complete_tree->F = F7*(complete_tree->T1 + complete_tree->T2 + complete_tree->T4 + complete_tree->T5 + complete_tree->T6) + complete_tree->T3 ;
	}



	if (COMMENT) {  // && (i==0)) {
		cout << "\nPopulation::aggregate_F" << endl;
		complete_tree->show_state();
		cout << "\nd_lim = " <<parameters->depth_lim;
		cout << "\ncomplete_tree->fitness = " << complete_tree->fitness;
		cout << "\nFit_ave = " << Fit_ave;
		cout << "\na1 = " << Fweight[1] << "  F1 = " << F[1] << "  a1*F1 = " << complete_tree->T[1];
		cout << "\na2 = " << Fweight[2] << "  F2 = " << F[2] << "  a2*F2 = " << complete_tree->T[2];
		cout << "\na3 = " << Fweight[3] << "  F3 = " << F[3] << "  a3*F3 = " << complete_tree->T[3];
		cout << "\na4 = " << Fweight[4] << "  F4 = " << F[4] << "  a4*F4 = " << complete_tree->T[4];
		cout << "\na5 = " << Fweight[5] << "  F5 = " << F[5] << "  a5*F5 = " << complete_tree->T[5];
		cout << "\na7 = " << Fweight[7] << "  F7 = " << F[7] << "  a7*F7 = " << complete_tree->T[7];
		cout << "\na8 = " << Fweight[8] << "  F8 = " << F[8] << "  a8*F8 = " << complete_tree->T[8];
		//for (int k=0; k<4; k++) cout << "\nlist[ " << k << " ] = " << list[k];
		cout << "\nF = " << complete_tree->F << endl;
	}

}


// uses the provided function to set the number of hits of each member
int Population::terminate(double THRESHOLD)
{
    // the evolution stops if there is an individual that scores n_test_cases hits...
    int check_end = 0;
	for (int i=0; i<1; i++) {          //if you want to check the whole population put size
		// terminate if number of hits = n test cases
		//if (trees[i]->hits==n_test_cases) {
		// terminate if fitness value (RMSE) is smaller than THRESHOLD
		//cout << "Fitness complete_trees[ " << i << "] = " <<  complete_trees[i]->fitness << endl;
		if ((complete_trees[i]->fitness)<THRESHOLD) {	
				check_end = 1;  
				// returns the number of the good individual
				return check_end;
		}
	}
	return check_end;
}




// returns the expression string for member n of the population contained in tree_array
char *Population::print (int n, Binary_Node** tree_array)
{
    // temporary trick to let fdf_c_ print individuals
	if (tree_array == NULL) 
		tree_array = complete_trees;

	// just call it and return it
    return (tree_array[n]->print());
}





// inserts parameters (terminal_const nodes with value=1) according to rules defined in Table 5.1 page 132 and at page 142 Armani's thesis
//input: pointer to tree root, pointer to pointer to tree root
void Population::parameters_allocation(Binary_Node *tree, Binary_Node**p_tree)
{
	int COMMENT =0;//1 for comments, 0 silent...
	Val value;
	if (COMMENT)
		cout << "\nParameters_allocation called";
	
	
		//INNER PARAMETERS (call function Binary_Node::check_allocation of the root node (binary_Node))
		tree->check_allocation(((Node **)&tree),10);    //10 is a fake number...
		if (COMMENT) 
			cout << "\nAllocation completed";
		
		// EXTERNAL PARAMETER
		// choose a random value for the final ("shift") parameter
		//value = constant_generation();
		// insert the external parameter with the chosen value (final parameter - shift of the total tree )
		insert_parameter((Node*)tree, (Node **)p_tree, &Add, (Val)1.0); 		//instead of value

		if (COMMENT) 
			cout << "\nFree parameter " << value << " added (shift of the tree)";

}


//function to add parameters to the whole population (recursive call to parameters_allocation)
void Population::population_parameters_allocation(void)  
{
	cout << "\nPopulation::population_parameters_allocation. Size = " << size << endl;
	for (int i=0; i<size; i++) {
		cout << "\n\nIndividual " << i << endl;
		parameters_allocation(trees[i],&trees[i]);
	}
}


// function that performs EDITING on the individuals without parameters
//PROBLEM : still at an early stage of development, as currently avoids only x/a for a=0,
// and a must be a single node to be recognised...
// input: number of the tree to be edited in trees[]
void Population::perform_editing(int i)
{
	int COMMENT = 0;
	// check if critical operations are used
	if (problem->division) { // add check if division is in the tree, so to skip in case this operation (editing)
		// gets here if division is among primitives

		if (COMMENT) {
			cout << "\n\n" << i << ") Tree before editing: ";
			print_individual((Node *)trees[i]);
		}

		edit_tree((Node*)trees[i], (Node**)&trees[i]);

		if (COMMENT) {
			cout << "\nDone";
			cout << "\n" << i << ") Tree after editing: ";
			print_individual((Node *)trees[i]);
		}

	}

}


// function to tune the numerical coefficients of a single tree (single or more initial guesses),
// based on RMSE computed on building data set (simple RMSE or PRESS approach)
// input: number of guesses, address of the parameterless tree, address of the complete tree, corresponding number of the tree
// output: integer (just to check)
int Population::tuning_individual(int n_guesses, Binary_Node *tree_no_par, Binary_Node *ntree, int tree_no)
{	
	int COMMENT = 0;
	int c;
	char* expr;
	int n_param, n_nodes, hits, hits_best;
	int IW;
	int n_tuning_parameters, n_tuning_parameters_best, n_corrections, n_corrections_best;
	bool param_inherited=false;
			
	double node_value; 
	double* x;	// array of tree parameters (numerical coefficients)
	double* x_best;   // best array of tree parameters (numerical coefficients)
	Val fitness, fitness_best;
	Val maxabserror, maxabserror_best;
	Val R2, R2_best;
	Val tree_mean, tree_mean_best;
	Val tree_variance, tree_variance_best;
	Val tree_min, tree_min_best;
	Val tree_max, tree_max_best;
	Val tree_acf_root, tree_acf_root_best;
	Val tree_tot_variation, tree_tot_variation_best;
	Val result[10];
	
	result[0] = (Val)0.0;  	//storing fitness value (error - RMSE)
	result[1] = (Val)0.0;	// storing maxabserror
	result[2] = (Val)0.0;	// storing n of corrections done by protected operations
	result[3] = (Val)0.0; 	// storing value of R squared (R2)
	result[4] = (Val)0.0;	// storing mean tree value on building data set
	result[5] = (Val)0.0;	// storing variance of tree values on building data set
	result[6] = (Val)0.0;	// storing min tree value on building data set
	result[7] = (Val)0.0;	// storing max tree value on building data set
	result[8] = (Val)0.0;	// storing first zero of autocorrelation function
	result[9] = (Val)0.0;	// storing total variation


	if (COMMENT) {
		cout << "\n\nPopulation::tuning_individual : TUNING TREE n. " << tree_no;
		cout << "\nTree before optimisation (values of parameters as they happen to be after tree creation)\n ";
		print_individual((Node *)ntree);
		cout << "\ntuning_individual called - tree no. " << tree_no;
	}
/*	
	//----------------------------------------------------------------------------------------------------------------------------------
	// FUNCTION TO ESTIMATE THE DERIVATIVE IN THE DIRECTION OF  THE CLOSEST NEIGHBOUR 
	// in the future call it from main (as the algorithm need to be run only once)	
	//----------------------------------------------------------------------------------------------------------------------------------
	// estimate the derivative of the unknown function in a randomly chosen fitness point (see DATA)
	// compute derivative - not tested for nested DoE (SPLIT = 1) !!!!!
	 // entry 0: value of the derivative of the unknown function in selected point,
	//  entry 1: row number of selected or initial point (where the derivative is computed)  
	//  entry 2: row number final point (see matrix data_tune)
	int n_der = 7;   //equal to VALIDATING_LINES when SPLIT=0
	cout << "n_der = " << n_der << endl; 
	double **fun_der = new double *[n_der];
	for (int k=0; k<n_der; k++)
		fun_der[k] = new double [3];  
	double *tree_der;
	get_fun_derivative(data_fitness, n_test_cases_fitness, num_vars, fun_der, n_der);  //to use this, SPLIT=0 necessarily!!!
	//get_fun_derivative(data, n_test_cases_tune, num_vars, fun_der);
	cout << "RESULT: number of points where derivatives are estimated: n_der = " << n_der << endl;
	for (int k=0; k<n_der; k++) 
	cout << "point " << k << ") df/ds = " << fun_der[k][0] << ", n. row initial point: " <<  (int)fun_der[k][1] << ", n. row final point: " <<  (int)fun_der[k][2] << endl;
	double tol = 3.0;
	int der_check;
	// -----------------------------------------------------------------------------------------------------------------------------	
	
*/	
	//--------------------------------------------------------------------------------------------------------------------------------
	// find the addresses of the parameters to be tuned (pulsations are recognised by find_pulsations)
	//--------------------------------------------------------------------------------------------------------------------------------
	n_param = ntree->n_tuning_parameters; // 16/2/22 removed =ntree->find_parameters();  // number of parameters found in the tree
	find_pulsations(ntree);    //update n_pulsations and index_puls (both Population members)
	n_nodes = ntree->count();

	if (COMMENT) {
		cout << "\n Tree n. " << tree_no << "\n : CHECK after find_pulsation function \n\n";
		print_individual((Node *)ntree);
		cout << "\nindex_puls = " << ntree->index_puls;
		cout << "\nParameters:" << endl;
		for (int k=0; k<n_param; k++)
			cout << " " << ((Terminal_Const *)ntree->p_par[k])->value(NULL);
		cout << "\nn_pulsations = " << ntree->n_pulsations << endl;
		if (ntree->n_pulsations)  {
			cout << "index_puls[] : " << endl;
			for (int k=0; k< ntree->n_pulsations; k++)
				cout << " " << ntree->index_puls[k];
			cout  << "\nPulsations:" << endl;	
			for (int k=0; k < ntree->n_pulsations; k++)
				cout << " " << (ntree->p_par[ntree->index_puls[k]])->value(NULL);
		}
	}


	
	//-----------------------------------------------------------------------------------------------------------
	// Creation and initialisation of the arrays x and x_best containing the parameters to be tuned 
	//-----------------------------------------------------------------------------------------------------------
	// Array x (parameters)  
	x = new (nothrow) double [n_param];
	if (x==NULL) {
		cerr << "\nPopulation::tuning_individual : ERROR ! Not enough memory to create x";
		exit(-1);
	}
	else 
		if (COMMENT) 
			cout << "\nPopulation::tuning_individual : array x successfully created : pointer = %i" << x;
	
	// Array x_best (best set of parameters)
	x_best = new (nothrow) double [n_param];
	if (x_best==NULL) {
		cerr << "\nPopulation::tuning_individual : ERROR! Not enough memory to create x_best ";
		exit (-1);
	}
	else 
		if (COMMENT) 
			cout << "\nPopulation::tuning_individual :  array x_best successfully created : pointer = %i" << x;
		
	// initialize parameters (if tuning is not successful, these are the values that show up)	
	for (int i=0; i<n_param; i++) {
		x[i] = 1.; 
		x_best[i] = 1.;		
	}

	//---------------------------------------------------------------------
	//initialization of state variables of best tree
	//---------------------------------------------------------------------	
	fitness_best = 999999E99;
	maxabserror_best = 999999E99;
	n_corrections_best = 999999;
	R2_best = 999999;
	tree_mean_best = 999999;
	tree_variance_best = 999999;
	tree_min_best=999999;
	tree_max_best=999999;
	tree_tot_variation_best=999999;
	//	result[0] = (Val)0.0;  	// fitness value (error - RMSE)
	//	result[1] = (Val)0.0;	// maxabserror
	//	result[2] = (Val)0.0;	// n of corrections done by protected operations
	//	result[3] = (Val)0.0; 	// value of R squared (R2)
	//	result[4] = (Val)0.0;	// mean tree value on building data set
	//	result[5] = (Val)0.0;	// variance of tree values on building data set
	//	result[6] = (Val)0.0;	// min tree value on building data set
	//	result[7] = (Val)0.0;	// max tree value on building data set
	//  result[8] = (Val)0.0;	// first zero of autocorrelation function
	//  result[9] = (Val)0.0;	// total variation
	
	//--------------------------------------------------------------------------------------------------------------------------------
	// HERE THE MULTIPLE GUESSES CYCLE STARTS 
	// n_guesses is the no of random initial guesses, imposed by the user 
	// (n_guesses is a member variable of Population -- see Population.h)
	//--------------------------------------------------------------------------------------------------------------------------------	
	int n_guess_ok=0;
	//int i_guess = 0;


	for (int i_guess=0; i_guess<n_guesses; i_guess++) {     
	//while (n_guess_ok < n_guesses) {
		// ----------------------------------------------------------------------------------------------------------------------------
		param_inherited=false;
		if (COMMENT) cout << "\nRANDOM GUESS n. " << i_guess << endl; 
		// if first guess and not initial generation (best_tree is not defined for generation 0) and the current tree is identical to best_tree in terms of structure
		if ((i_guess==0) && (best_tree) && (identical_trees(tree_no_par, best_tree))) { // why not leaving only i_guess==0 ?

		// PARAMETERS INHERITED if tree is identical to best tree. Why not inheriting coefficients for all individuals reproduced?
			if (COMMENT) { //start comment
				cout << "\nbest_tree and trees[" << ntree << "] are identical:";
				cout << "\nbest_tree : ";
				print_individual((Node *)best_tree);
				cout << "\ntree_no_par : ";
				print_individual((Node *)tree_no_par);
				cout << "\nparameters inherited...";
				cout << "\nbest_complete_tree : ";
				print_individual((Node *)best_complete_tree);
				cout << "\ncomplete_trees[" << ntree << "] : ";
				print_individual((Node *)ntree);
				Val node_value;
				cout << endl;
				for (int i=0; i<n_param; i++) {
					node_value = ((Terminal_Const *)best_complete_tree->p_par[i])->value(NULL);
					cout << "value [" << i << "] = " << node_value << endl;
				}
				cout << "\n  x inherited (size " << n_param << ").";
			} //end comment

			//fetch the parameters from best_complete_tree
			for (int i=0; i<n_param; i++) x[i] = ((Terminal_Const *)best_complete_tree->p_par[i])->value(NULL); //values of the parameters
			param_inherited=true;
		} // end if

		else {
			// PARAMETERS RANDOMLY GENERATED
			// initialization of x (random guess)
			int r=0;			
			for (int j=0; j<n_param; j++) {
				// check if x[j] is a pulsation		
				if (!(ntree->n_pulsations) || (r >= ntree->n_pulsations))
					x[j] = constant_generation(parameters->minrand, parameters->maxrand, parameters->step);
				else {
					// there's  at least one pulsation				
					if (j== ntree->index_puls[r])	{
						Val a = (problem->v_list[ntree->index_var[r]])->omega_lim;
						x[j] = constant_generation(parameters->minrand, parameters->maxrand, parameters->step, -a, a);
						//x[j] = constant_generation();	// more general, but relies on the right guess!					
						r++;
					}
					else 
						x[j] = constant_generation(parameters->minrand, parameters->maxrand, parameters->step);
				}
			}
			if (COMMENT)  cout << "\n  x randomly initialized (size " << n_param << ").";

		}

		// show parameters (inherited or randomly chosen)
		if (COMMENT) {
			if (param_inherited)
				cout << "\nTuning_individual : parameters inherited";
			else
				cout << "\nTuning_individual : parameters randomly generated";
			cout << "\nx = [";
			for (int j=0; j<n_param-1; j++) 
					cout << x[j] << " , ";
			cout << x[n_param-1] << "]" << endl;
		}
				
		if (COMMENT) {
			cout << "\n  address Population P = " << this; 
			cout << "\n Number of nodes = " << n_nodes;
			cout << "\n  Number of parameters: n = " << n_param;
			cout << "\n  Number of functions/test cases: m = " << n_test_cases_tune;	
			cout << "\n  Dimension of w = " << IW;
			cout << endl;
		}	  
   		
		// ntree_fdf_c is public member of Population
		ntree_fdf_c = tree_no;
		if (COMMENT) cout << "\n  ntree_fdf_c = " << ntree_fdf_c;


		// if parameters are inherited from best tree, then just reevaluate objectives without tuning
		if (param_inherited) {
			ntree->n_corrections=0;
			// update ntree with inherited parameters prior to fitness evaluation
			update_complete_tree(ntree, x, n_param);
			// 16/2/22 if parameters are inherited, why not also copying fitness function value? It makes no sense to recompute it!
			fitness_func(problem->Sy, problem->get_data_address(), problem->get_n_data(), ntree, result, parameters->normalised, parameters->crossvalidation);

			// update tree properties:
			fitness = result[0];  // RMSE error (CROSSVALIDATION=0) or PressRMSE error (CROSSVALIDATION=1)
			maxabserror = result[1];
			n_corrections = (int)result[2];
			R2 = result[3];
			tree_mean = result[4];
			tree_variance = result[5];
			tree_min=result[6];
			tree_max=result[7];
			tree_acf_root = result[8];
			tree_tot_variation = result[9];
		} else {
			// if parameters are randomly generated then tune the tree
			// tune and validate a single individual using PRESS crossvalidation, using randomly generated initial guesses for x (numerical coefficients)
			if (parameters->crossvalidation==1) tuning_individual_PRESS_single_guess(ntree, x, &n_param, &n_guess_ok, result);
			// tune and validate a single individual using simple RMSE approach, using randomly generated initial guesses for x (numerical coefficients)
			if (parameters->crossvalidation==0) tuning_individual_RMSE_single_guess(ntree, x, &n_param, &n_guess_ok, result);

			// as a result of tuning, result[] is used to update tree properties:
			fitness = result[0];  // RMSE error (CROSSVALIDATION=0) or PressRMSE error (CROSSVALIDATION=1)
			maxabserror = result[1];
			n_corrections = (int)result[2];
			R2 = result[3];
			tree_mean = result[4];
			tree_variance = result[5];
			tree_min=result[6];
			tree_max=result[7];
			tree_acf_root = result[8];
			tree_tot_variation = result[9];
			// also the array of coefficients x is returned...
		}

/*
		// SELECTION of the BEST SET of PARAMETERS (x_best update) 
		// -------------------------------------------
		// here you can check the derivatives, the correlation matrix, ecc, of the current individual with the best individual so far
		// the one that gives closer result to the derivatives, correlation matrix, ecc, computed from the original output data is chosen
		// -------------------------------------------
		// compute the derivative of the tree in the selected point and in the direction given by get_fun_derivative(...)
		tree_der = new double [n_der];	
		get_tree_derivative(complete_trees[ntree], fun_der, n_der, tree_der);
		cout << "\nRESULTS from get_tree_derivative: n_der = " << n_der << endl;
		for (int k=0; k<n_der; k++) 
			cout << "point " << k << ") derivative = " << tree_der[k] << endl;

		for (int j=0; j<n_der; j++) {
			cout << 1./tol << " ### " << tree_der[j]/fun_der[j][0] << " ### " << tol << endl;
			if ( (tree_der[j]/fun_der[j][0]>0.0) && (tree_der[j]/fun_der[j][0]<tol)) {
				der_check = 1;
				cout << "ok" << endl;
			}
			else {
				der_check =0;
				cout << "FAILED!" << endl;
			}
		}
		//-------------------------------------------------------------------------------------------------------
*/		
		//-----------------------------------------------------------------------------------------------------
		// UPDATE local best results if the current tuning cycle has given better results
		// Mind!!!! Coefficient tuning is based on RMSE, not on the aggregated function F (obvious, as you are not comparing trees with the same structure)!!!
		//-----------------------------------------------------------------------------------------------------
		//if ((fitness < fitness_best) && (method>=0) && (der_check)) {
		// what if the paramaters x is not acceptable (all NaN)? Implement a function to check! Here or in tuning_individual_PRESS_single_guess
		if (fitness < fitness_best)  {	   // single-objective comparison on validation error (measured on data_validation)- RMSE or PressRMSE (CROSSVALIDATION=1)
			fitness_best = fitness;
			maxabserror_best = maxabserror;
			n_corrections_best = n_corrections; 
			R2_best = R2;
			tree_mean_best = tree_mean;
			tree_variance_best = tree_variance;
			tree_min_best=tree_min;
			tree_max_best=tree_max;
			tree_acf_root_best = tree_acf_root;
			tree_tot_variation_best = tree_tot_variation;

			for (int j=0; j<n_param; j++)
				 x_best[j] = x[j];  //
		}
		
		//end of the single guess
		if (COMMENT) {
			expr = ntree->print();
			cout << "\nPopulation_tuning_individual : tree after optimisation n. " << i_guess << " (optimised set of parameters inserted)\n " << expr << endl;
			delete [] expr;			
			cout << "\nEnd of guess n. " << i_guess << endl;	
		}


	}	// here the cycle related to the number of random guesses on tree coefficients (i_guess<n_guesses) stops -------

	//---------------------------------------------------------------------------------------------------------------
	// after tuning cycles, update state variables of the (complete) tree with the best set of (tuned) parameters
	//---------------------------------------------------------------------------------------------------------------
	update_complete_tree(ntree, x_best, n_param);
	ntree->fitness = fitness_best;
	ntree->maxabserror = maxabserror_best;
	ntree->n_corrections = n_corrections_best;	
	ntree->R2 = R2_best;
	ntree->tree_mean = tree_mean_best;
	ntree->tree_variance = tree_variance_best;
	ntree->tree_min=tree_min_best;
	ntree->tree_max=tree_max_best;
	ntree->first_acf_root_tree = tree_acf_root_best;
	ntree->tot_variation_tree = tree_tot_variation_best;

	if (COMMENT) {	
		//  print the value of the tuned constants
		for (int i=0; i<n_param; i++) { 
			node_value = ((Terminal_Const *)ntree->p_par[i])->value(NULL);
			cout << "\nPopulation::tuning_individual : Pointer to " << i+1 << "-th const : " << ntree->p_par[i];
			cout << ", Population::tuning_individual : Constant value = " << node_value;
		}	
		expr = ntree->print();
		cout << "\n Population::tuning_individual : Tree after optimisation (best set of parameters inserted)\n " << expr << endl;
		// free memory	
		delete [] expr;	

		cout << "\nPopulation::tuning_individual : Tree n. " << tree_no << " TUNED (constants optimized)." << endl;
	}
	
	// uncomment if you use derivative estimates
	//for (int k=0; k<n_der; k++) 
	//	delete [] fun_der[k];
	//delete [] fun_der; 
	//delete [] tree_der;
	
	// free memory used for x and x_best
	delete[] x;
	delete[] x_best;
	
	// free array with indexes of pulsations (index_puls and index_var are allocated dynamically in Population::find_pulsations)
	if (ntree->n_pulsations) {
 		if (ntree->n_pulsations == 1) { 
			delete ntree->index_puls;
			delete ntree->index_var;
		}
		else {
			delete [] ntree->index_puls;
			delete [] ntree->index_var;
		}
	}

	return 1;    
}



//-/-/-/-/-/-/--/-/-/-/-/-/-/-/
// function that performs standard RMS error evaluation on a single individual, for a single guess of the individual parameters
// Inputs:
// - ntree: pointer to the tree to be tuned
// - x_original: pointer to array storing original values of the numerical coefficients (parameters)
// - p_n_param: pointer to variable with number of numerical coefficients
// - p_n_guess_ok: pointer to variable with number of guesses (check it...)
// - p_result: pointer to array in which the results of the numerical optimisation - tree objectives - will be stored as output
// Outputs:
// - None
void Population::tuning_individual_RMSE_single_guess(Binary_Node *ntree, double* x_original, int* p_n_param, int* p_n_guess_ok, Val* p_result)
{
	int COMMENT=0;
	int m_total=0;
	int IW=0;
	int method=0; // method used as input in opti:
				// for MINL2.cpp: (derivatives values are used)
				//1		->start minimization
				//2		->start minimization and print information during the iteration
				//=0	-> check gradient. No iteration
				// for MI0L2.cpp: (derivatives values are NOT used)
				//0		->start minimization with no gradient
				//1		->start minimization with gradient computed using function values (no actual differentiation)
	// for tuning (*p_n_param is the number of parameters found in the tree)
	int n_param_threshold = floor(0.5*problem->n_tuning);  // if the coefficients are too many, tuning is not performed...

	if (COMMENT) cout << "\nPopulation::tuning_individual_RMSE_single_guess : enter";


	// perform a last tuning and validation on the whole data set to get all other outputs right, including array of parameters x
	problem->n_tuning = problem->get_n_data();
	problem->n_validation = problem->get_n_data();
	m_total = problem->n_tuning + ntree->n_pulsations;

	//--------------------------------------------------------------------
	// set the size of vector W (see SQP_guide pag 13)
	//--------------------------------------------------------------------
	if (method>0)
		IW=10*(2*m_total*(*p_n_param+1)+*p_n_param*(*p_n_param+3));
	else
		IW=10*(2*m_total*(*p_n_param+2)+*p_n_param+10);


	//-----------------------------------------------------------
	// TREE TUNING on the whole data set
	// call the optimization functions: HYPSO (C++) and/or SQP (Fortran) (LINKS TO EXTERNAL FUNCTION)
	//-----------------------------------------------------------



	// PSO parameters optimisation (just try to run HyPSO for top get closer to global minima and provide good start guesses for SQP)
	// steps involved: call to HYPSO/pso_launcher, launch psominimize, internal call to objective function in model.cpp (tree evaluator and error definition with weights for pulsations...)
	// 8/11/21 test: PSO brought out of the control on number of coefficients
	// -------------------------------------------------------------------------------------------
	if (COMMENT) {
		cout << "\n\nPopulation::tuning_individual_RMSE_single_guess - x_original before HyPSO:" << endl;
		for (int i=0; i<*p_n_param; i++) cout << "x_original[" << i << "]=" << x_original[i] << endl;
		cout << "\n\nPopulation::tuning_individual_RMSE_single_guess : call hypso launcher";
	}
	hypso_launcher(parameters, ntree, *p_n_param, x_original); // method of Population
	if (COMMENT) {
		cout << "\n\nPopulation::tuning_individual_RMSE_single_guess - x_original after HyPSO:" << endl;
		for (int i=0; i<*p_n_param; i++) cout << "x_original[" << i << "]=" << x_original[i] << endl;
	}

	// 25/7/2016 : HERE YOU COULD ADD A CHECK TO SKIP TUNING (AND SO LEAVE CONSTANT NODES RANDOMLY INITIALISED)
	//             IF THE NUMBER OF PARAMETERS IS LARGER THAN SAY HALF THE QUANTITY OF THE TRAINING POINTS...
	//             SEE BISHOP 1996 ARTICLE
	// IF (CONDITION FOR TUNING is TRUE) THEN    - example n_param <= n_param_max = 0.5 * size training data set
	if (*p_n_param<=n_param_threshold) {

		// 8/11/21 PSO removed from here

		// SQP parameters optimisation
		// ---------------------------
		if (COMMENT) cout << "\nPopulation::tuning_individual_RMSE_single_guess : before optimisation ICONTR = method = " << method;
		// to perform optimization with TINL2_mod.cpp, MINL2.cpp through f2c: DOESN'T WORK
		//c = opti_cpp(this, &method, n_param, n_test_cases, ntree,x);

		// CALL FORTRAN FUNCTION: ALWAYS USE ADDRESSES OF VARIABLES, NOT VARIABLES DIRECTLY!
		// ALSO: DON'T USE CAPITAL LETTERS... otherwise problems during compiling ("undefined reference")
		//to perform optimization with TINL2.FOR  - Umberto's method (copied from Andrey's)
		//optitinl2_(&method, x, &n_param, &n_test_cases_tune, &IW);
		if (COMMENT) cout << "\n\nPopulation::tuning_individual_RMSE_single_guess : *p_n_param = " << *p_n_param << ", n_tuning = " << problem->n_tuning;

		//to perform optimization with TIOL2.FOR  - Andrey's method - WORKS PERFECTLY!
		if (COMMENT) cout << "\nPopulation::tuning_individual_RMSE_single_guess : Call opti_()";
		opti_(&method, x_original, p_n_param, &m_total, &IW);

		if (COMMENT) {
			cout << "\nPopulation::tuning_individual_RMSE_single_guess : after optimisation ICONTR = method = " << method;
			if (method==2) cout << "\n => UPPER LIMIT FOR FUNCTION EVALUATIONS EXCEEDED.";
			if (method==1) cout << "\n => SUM OF SQUARES FAILS TO DECREASE";
			if (method==0) cout << "\n => successful convergence";
			if (method<0) cout << "\n => COMPUTATION DID NOT START: see SQP_guide (pag 14) for details";
		}

	} else {
		// if tuning is not performed, the tree coefficients stay as they are (randomly generated or improved by PSO)
		if (COMMENT) cout << "\nPopulation::tuning_individual_RMSE_single_guess : *p_n_param>n_param_threshold : tuning not performed, tree coefficients remain randomly generated" << endl;
	}

	// update the tree with the new, optimised values of the parameters
	update_complete_tree(ntree, x_original, *p_n_param);


	//-----------------------------------------------------------
	// TREE VALIDATION on the whole data set
	//-----------------------------------------------------------
	if (method==0) (*p_n_guess_ok)++;
	if (COMMENT) {
		cout << "\nPopulation::tuning_individual_RMSE_single_guess : show current data_validation" << endl;
		cout << "problem->n_validation = " << problem->n_validation;
		//problem->show_data_validation();
	}
	ntree->n_corrections=0;
	fitness_func(problem->Sy, problem->get_data_address(), problem->get_n_data(), ntree, p_result, parameters->normalised, parameters->crossvalidation);



	//------------------------------------------------------------------------------------
	// RETURN RESULTS AND UPDATE INDIVIDUAL PARAMETERS (COEFFICIENTS) AFTER CROSSVALIDATION
	//------------------------------------------------------------------------------------

	// finally, define output values
	//p_result[0];  // error predictor if CROSSVALIDATION=1 ok, if not...?   result[0] : error (fitness value - RMSE)
	//	result[1] = (Val)0.0;	// storing n of hits
	//	result[2] = (Val)0.0;	// storing n of corrections done by protected operations
	//	result[3] = (Val)0.0; // storing value of R squared (R2)
	//	result[4] = (Val)0.0;	// storing mean tree value on building data set
	//	result[5] = (Val)0.0;	// storing variance of tree values on building data set
	//	result[6] = (Val)0.0;	// storing min tree value on building data set
	//	result[7] = (Val)0.0;	// storing max tree value on building data set
	// what if they are all NaN?! The individual expression cannot be printed... CHECK
	for (int i=0; i<(*p_n_param); i++) {
		if (isnan(x_original[i])) x_original[i]=1.0;   // correction in case of illegal parameter
	}


	if (COMMENT) {
		cout << "\nPopulation::tuning_individual_RMSE_single_guess : tuning and validation complete" << endl;
		cout << "Population::tuning_individual_RMSE_single_guess : parameters x chosen for insertion:" << endl;
		cout << "\nx = [";
		for (int j=0; j<*p_n_param-1; j++) cout << x_original[j] << " , ";
		cout << x_original[*p_n_param-1] << "]" << endl;
	}

	if (COMMENT) cout << "\n\nPopulation::tuning_individual_RMSE_single_guess : exit";

}



// function that performs the PRESS error evaluation (CROSSVALIDATION) on a single individual, for a single guess of the individual parameters
void Population::tuning_individual_PRESS_single_guess(Binary_Node *ntree, double* x_original, int* p_n_param, int* p_n_guess_ok, Val* p_result)
{
	int COMMENT=1;
	int m_total=0;
	int IW=0;
	int method; // method used as input in opti:
				// for MINL2.cpp: (derivatives values are used)
				//1		->start minimization
				//2		->start minimization and print information during the iteration
				//=0	-> check gradient. No iteration
				// for MI0L2.cpp: (derivatives values are NOT used)
				//0		->start minimization with no gradient
				//1		->start minimization with gradient computed using function values (no actual differentiation)

	// for tuning (*p_n_param is the number of parameters found in the tree)
	int n_param_threshold = floor(0.5*problem->n_tuning);  // if the coefficients are too many, tuning is not performed...
	double* x=new double[*p_n_param];
	double* x_out=new double[*p_n_param];
	Val* exv;
	double exv_tot=0.0;
	double PressRMS=0.0;

	if (COMMENT) cout << "\nPopulation::tuning_individual_PRESS_single_guess : enter";

	// array containing the crossvalidation errors (see Viana "Making the most out of surrogate models:" 2010, pag.4
	exv = new Val[problem->get_n_folds()];
	for (int i=0; i<(*p_n_param); i++) x_out[i]=0.0;

	// start crossvalidation cycle: at each iteration a erms is computed, which will be used to compute PRESS
	for (int vf=0; vf<problem->get_n_folds(); vf++) {   // vf stands for validation fold...

		// initialise method (start minimization with no gradient) to prevent changes from opti execution to change the aim of the optimisation
		// (ex:if method is changed by opti from 0 to 1 - SUM OF SQUARES FAIL TO DECREASE - this 1 is taken as input for the opti at the next iteration, stopping the optimisation
		method=0;

		// initialise array of coefficient x = x_original (so all tuning processes starts with the identical tree)
		for (int i=0; i<(*p_n_param); i++) x[i]=x_original[i];

		// show optimised parameters after all the tuning in the different building data sets
		if (COMMENT) {
			cout << "\nPopulation::tuning_individual_PRESS_single_guess : vf = " << vf << "/////////////" << endl;
			cout << "\nPopulation::tuning_individual_PRESS_single_guess : parameters BEFORE optimisation vf = " << vf << " :";
			cout << "\nx = [";
			for (int j=0; j<*p_n_param-1; j++) cout << x[j] << " , ";
			cout << x[*p_n_param-1] << "]" << endl;
		}

		// set validation fold
		problem->set_validation_fold(vf);
		problem->n_tuning = problem->get_n_data() - problem->get_points_per_fold(vf);

		// set the total number of summands in SQP error metric for tuning individuals (see Armani PhD thesis pag. 143)
		m_total = problem->n_tuning + ntree->n_pulsations;

		if (COMMENT) {
			cout << "\nPopulation::tuning_individual_PRESS_single_guess : *p_method = " << method << endl;
			cout << "Population::tuning_individual_PRESS_single_guess : problem->get_validation_fold(vf) = " << problem->get_validation_fold() << endl;
			cout << "Population::tuning_individual_PRESS_single_guess : problem->n_tuning = " << problem->n_tuning << endl;
			cout << "Population::tuning_individual_PRESS_single_guess : problem->data_tuning = " << problem->data_tuning << endl;
			cout << "Population::tuning_individual_PRESS_single_guess : problem->data = " << problem->get_data_address() << endl;
			cout << "Population::tuning_individual_PRESS_single_guess : problem->get_fold_from_row(0) = " << problem->get_fold_from_row(0) << endl;
			cout << "Population::tuning_individual_PRESS_single_guess : *p_n_param = " << *p_n_param << "   n_param_threshold = " << n_param_threshold << endl;
		}

		//--------------------------------------------------------------------
		// set the size of vector W (see SQP_guide pag 13)
		//--------------------------------------------------------------------
		if (method>0)
			IW=10*(2*m_total*(*p_n_param+1)+*p_n_param*(*p_n_param+3));
		else
			IW=10*(2*m_total*(*p_n_param+2)+*p_n_param+10);


		//-----------------------------------------------------------
		// TREE TUNING on tuning data set = all data except the current validation data set
		// call the optimization functions: HYPSO (C++) and/or SQP (Fortran) (LINKS TO EXTERNAL FUNCTION)
		//-----------------------------------------------------------

		// 25/7/2016 : HERE YOU COULD ADD A CHECK TO SKIP TUNING (AND SO LEAVE CONSTANT NODES RANDOMLY INITIALISED)
		//             IF THE NUMBER OF PARAMETERS IS LARGER THAN SAY HALF THE QUANTITY OF THE TRAINING POINTS...
		//             SEE BISHOP 1996 ARTICLE
		// IF (CONDITION FOR TUNING is TRUE) THEN    - example n_param <= n_param_max = 0.5 * size training data set
		if (*p_n_param<=n_param_threshold) {

			// insert here call to HyPSO (just try to run HyPSO for a couple of iterations to search for global minima)
			// steps involved: call to HYPSO/pso_launcher, launch psominimize, internal call to objective function in model.cpp (tree evaluator and error definition with weights for pulsations...)
			//double (*p_objfun)(double*, int, Binary_Node*);
			//p_objfun = &Population::pso_objfunction;  // method in Population
			hypso_launcher(parameters, ntree, *p_n_param, x_original); // method of Population

			// SQP parameters optimisation
			if (COMMENT)
				cout << "\nPopulation::tuning_individual_PRESS_single_guess : before optimisation ICONTR = method = " << method;
			// to perform optimization with TINL2_mod.cpp, MINL2.cpp through f2c: DOESN'T WORK
			//c = opti_cpp(this, &method, n_param, n_test_cases, ntree,x);

			// CALL FORTRAN FUNCTION: ALWAYS USE ADDRESSES OF VARIABLES, NOT VARIABLES DIRECTLY!
			// ALSO: DON'T USE CAPITAL LETTERS... otherwise problems during compiling ("undefined reference")
			//to perform optimization with TINL2.FOR  - Umberto's method (copied from Andrey's)
			//optitinl2_(&method, x, &n_param, &n_test_cases_tune, &IW);
			if (COMMENT) cout << "\n\nPopulation::tuning_individual_PRESS_single_guess : *p_n_param = " << *p_n_param << ", n_tuning = " << problem->n_tuning;
			if (1==1) {
				//to perform optimization with TIOL2.FOR  - Andrey's method - WORKS PERFECTLY!
				if (COMMENT) cout << "\nPopulation::tuning_individual : Call opti_()";
				opti_(&method, x, p_n_param, &m_total, &IW);
			} else {
				if (COMMENT) cout << "\nPopulation::tuning_individual_PRESS_single_guess : tuning through opti_() skipped ...";
			}


			if (COMMENT) {
				cout << "\nPopulation::tuning_individual_PRESS_single_guess : after optimisation ICONTR = method = " << method;
				if (method==2) cout << "\n => UPPER LIMIT FOR FUNCTION EVALUATIONS EXCEEDED.";
				if (method==1) cout << "\n => SUM OF SQUARES FAILS TO DECREASE";
				if (method==0) cout << "\n => successful convergence";
				if (method<0) cout << "\n => COMPUTATION DID NOT START: see SQP_guide (pag 14) for details";
			}

		} else {
			// if tuning is not performed, the tree coefficients stay as they are, randomly generated
			if (COMMENT) cout << "\nPopulation::tuning_individual_PRESS_single_guess : *p_n_param>n_param_threshold : tuning not performed, tree coefficients remain randomly generated" << endl;
		}


		// show optimised parameters after all the tuning in the different building data sets
		if (COMMENT) {
			cout << "\nPopulation::tuning_individual_PRESS_single_guess : parameters AFTER optimisation vf =" << vf << " :";
			cout << "\nx = [";
			for (int j=0; j<*p_n_param-1; j++) cout << x[j] << " , ";
			cout << x[*p_n_param-1] << "]" << endl;
		}


		//-----------------------------------------------------------------------------------------------------
		// TREE VALIDATION on the current validation data set vf
		// (it doesn't assign the values to the tree...wait for the final update)
		//-----------------------------------------------------------------------------------------------------
		// define data_validation, made of points in validation fold
		// you can either actually copy values, or use the association table and check each point during error calculation as done in tuning
		// qui
		int row;
		if (parameters->crossvalidation==1) {
			problem->data_validation = new Val*[problem->get_points_per_fold(vf)];
			for (int i=0; i<problem->get_points_per_fold(vf); i++)
				problem->data_validation[i]= new Val[problem->get_n_cols()];

			// fill data_validation with correct values through folds_table
			row=0;
			for (int i=0; i<problem->get_n_data(); i++) {
				if (problem->get_fold_from_row(i)==vf) {
					for (int j=0; j<problem->get_n_cols(); j++)
						problem->data_validation[row][j]=problem->get_data(i,j);
					row++;
				}
			}

			problem->n_validation = row;
		}  // end cycle crossvalidation=1


		if (COMMENT) {
			cout << "\n\nPopulation::tuning_individual_PRESS_single_guess : tree VALIDATION" << endl;
			cout << "Population::tuning_individual_PRESS_single_guess : show current data_validation" << endl;
			cout << "problem->n_validation = " << row;
			problem->show_data_validation();
		}

		// evaluate individual on validation data set
		fitness_func(problem->Sy, problem->data_validation, problem->n_validation, ntree, p_result, parameters->normalised, parameters->crossvalidation);   //IMPORTANT: fitness evaluated on data_validation !!!

		//-----------------------------------------------------------------------------------------------------
		// STORE Exv resulting from tuning and validation on different data sets
		//-----------------------------------------------------------------------------------------------------
		exv[vf]=p_result[0];
		cout << "\nPopulation::tuning_individual_PRESS_single_guess : exv[" << vf << "] = "<< exv[vf];

		// delete data_validation? OK, but then you have to copy DATA in data_validation in read_input otherwise the standard execution does not work (CROSSOVALIDATION=0)...

		for (int i=0; i<(*p_n_param); i++) x_out[i]=x_out[i]+x[i];

	} // end cycle in folds  -------------------------------------------


	//------------------------------------------------------------------------------------
	// RETURN RESULTS AND UPDATE INDIVIDUAL PARAMETERS (COEFFICIENTS) AFTER CROSSVALIDATION
	//------------------------------------------------------------------------------------
	// compute PRESS and assign it to result[0]
	for (int i=0; i<problem->get_n_folds(); i++) exv_tot = exv_tot + exv[i]*exv[i];
	PressRMS=sqrt(exv_tot/problem->get_n_data());

	if (COMMENT) cout << "\nPopulation::tuning_individual_PRESS_single_guess : PressRMS = " << PressRMS;

///
	// perform a last tuning and validation on the whole data set to get all other outputs right, including array of parameters x
	problem->n_tuning = problem->get_n_data();
	m_total = problem->n_tuning + ntree->n_pulsations;
	method=0;
	if (method>0)
		IW=10*(2*m_total*(*p_n_param+1)+*p_n_param*(*p_n_param+3));
	else
		IW=10*(2*m_total*(*p_n_param+2)+*p_n_param+10);
	for (int i=0; i<(*p_n_param); i++) x[i]=x_original[i];
	if (*p_n_param<=n_param_threshold) opti_(&method, x, p_n_param, &m_total, &IW);
	if (method==0) (*p_n_guess_ok)++;
	// evaluate individual on validation data set
	ntree->n_corrections=0;
	fitness_func(problem->Sy, problem->get_data_address(), problem->get_n_data(), ntree, p_result, parameters->normalised, parameters->crossvalidation);

///

	// finally, define output values
	p_result[0]=PressRMS;  // error predictor if CROSSVALIDATION=1 ok, if not...?   result[0] : error (RMSE)
	//p_result[1] // result[1] : n. of hits scored
	//p_result[2]	// result[2]: n. of corrections done by protected operations
	//p_result[3]	// result[3]: R2
	// define which array x of tree parameters... currently x obtained from the last tuning on the whole data set
	// what if they are all NaN?! The individual expression cannot be printed... CHECK
	for (int i=0; i<(*p_n_param); i++) {
		if (isnan(x[i])) x[i]=1.0;   // correction in case of illegal parameter
		x_original[i]=x[i];  // x is the array of parameters chosen as final
	}

	// update the tree with the new, optimised values of the parameters
	update_complete_tree(ntree, x_original, *p_n_param);


	if (COMMENT) {
		cout << "\nPopulation::tuning_individual_PRESS_single_guess : crossvalidation complete" << endl;
		cout << "Population::tuning_individual_PRESS_single_guess : parameters x chosen for insertion:" << endl;
		cout << "\nx = [";
		for (int j=0; j<*p_n_param-1; j++) cout << x_original[j] << " , ";
		cout << x_original[*p_n_param-1] << "]" << endl;
	}


	// delete
	delete[] exv;
	delete[] x, x_out;


	if (COMMENT) cout << "\n\nPopulation::tuning_individual_PRESS_single_guess : exit";

}



// function to evaluate a tree at population level (n_tree_evaluations is updated)
// (if you use value in Binary_Node n_tree_evaluations is NOT updated!)
Val Population::tree_value (Binary_Node* p_tree, int* p_n_corrections)
{
	Val v;
	// evaluate the tree
	// (the pointer between parentheses is NULL if you don't want to update n_corrections member)
	v = p_tree->value(p_n_corrections);    //value(NULL) as here there's no interest in counting the number of corrections (it is done only in fitness_func)

	// update evaluation counter
	n_tree_evaluations++;

	// returns the value of the tree # ntree (from 0 to size-1)
	return v;
}


// function to update the tree's parameters with the optimized ones (operates on complete_trees)
// (called recursively by fdf_) 
// Input:
// - Pointer to Tree root (Binary Node*)
// - Pointer to parameters array (tree new coefficients)
// - Size of parameters array
// Output: Null
void Population::update_complete_tree(Binary_Node *ctree, double *x, int n)
{	
	int COMMENT=0;
	if (COMMENT) cout << "\n\nPopulation::update_complete_tree" << endl;
	int u=0;

	if (ctree->p_par.size()==0) {
		cerr << "Error! Array of coefficients addresses not built! Exit..." << endl;
		exit(-1);
	}

	for (int i=0; i<n; i++) {
		 if (COMMENT) cout << "x[" << i << "]=" << x[i] << endl;
		 u=0;
		 if (COMMENT)cout << "p_par[" << i << "]=" << ((Terminal_Const *)ctree->p_par[i])->value(&u) << endl;
		 ((Terminal_Const *)ctree->p_par[i])->assign(x[i]);
		 if (COMMENT)cout << "p_par[" << i << "]=" << ((Terminal_Const *)ctree->p_par[i])->value(&u) << endl;
	}

	if (COMMENT) cout << "\n\nPopulation::update_complete_tree : exit" << endl;
}



// function which computes the partial derivative of a single term of
// the fitness function respect to a single parameter among those
// which have to be tuned
// input : actual value of the parameter, reference tree, parameter identificative number
// output : partial derivative of the fitness f. term respect to the parameter  
double Population::jacobian_ij (double x_curr, int ntree, int j)
{
// differentiation: centred method (df/dx = (f(x+dx) - f(x-dx)/2dx)
	
	double dx;
	double x_s, x_f;
	double f_s, f_f;
	double der_ij;
    Terminal_Const *p;

	dx = 0.00001;    //displacement: just a guess, now
	// value in x_j + dx_j
	x_s = x_curr + dx;
	p = (Terminal_Const *)(complete_trees[ntree]->p_par[j-1]);
	p->assign(x_s);
	f_s = tree_value(complete_trees[ntree], NULL);     //value(NULL) as here there's no interest in counting the number of corrections
	
	// value in x_j - dx_j
	x_f = x_curr - dx;
	p->assign(x_f);
	f_f = tree_value(complete_trees[ntree], NULL);    //value(NULL) as here there's no interest in counting the number of corrections
	
	

	// final value of dfi/dxj:
	der_ij = (f_s-f_f)/(2.*dx);     // - or not?
 	
	//set the parameters back to the initial value (if not, the parameter remain x_f, that is the actual value - dx !!!)
	p->assign(x_curr);

	return der_ij;
}


// function that computes the derivative of a tree in a given point (initial point).
// The direction is found declaring a second point (called final point)
// INPUT:
// data_used - matrix containing the points used to define initial and final points
// fun_der - matrix containing the number of the row corresponding to initial point (column 1) and final point (column 2) in matrix data_used
// n_der - number of derivatives requested
//
// OUTPUT:
// void
void  Population::get_tree_derivative_given_points(Val **data_used,Binary_Node *current_tree, double **fun_der, int n_der, double *tree_der)
{
	int COMMENT = 0;
	// step for derivative
	double dv = 1.0E-06;
	double tree_initial, tree_final;
	double *initial = new double [parameters->nvar];
	double *final = new double [parameters->nvar];
	double* conn_v = new double [parameters->nvar];
	cout << "\nPopulation::get_tree_derivative" << endl;
	
	for (int row=0; row<n_der; row++) {
		
		// extract nvar coordinates of initial point
		for (int k=0; k<parameters->nvar; k++)
			initial[k] = data_used[(int)(fun_der[row][1])][k];  //previous version (see get_deriv_to_closest_neighbour)
			// INSERT HERE THE MATRIX
		if (COMMENT)
			if (parameters->nvar==2) cout << "\n initial = (" << initial[0] << ", " << initial[1] << ")" << endl;
	
		//evaluate tree in initial point
		for (int k=0; k<parameters->nvar; k++) 			// assign the right value to all the variables for the i-th sample case
				(problem->v_list[k])->value = initial[k];
		tree_initial = tree_value(current_tree, NULL);    // NULL as n_corrections is not a matter here
		if (COMMENT) cout << "\nValue of the tree at the initial point : tree_initial = " << tree_initial << endl;

		// find final point
		// find connecting vector
		vect_sub(&(data_used[(int)(fun_der[row][2])][0]), initial, parameters->nvar, conn_v);  //ok till here
		if (COMMENT)
			if (parameters->nvar==2) cout << "Connecting vector : conn_v = (" << conn_v[0] << ", " << conn_v[1] << ")" << endl;

		// compute connecting vector magnitude
		double magn;
		magn = vect_magn(conn_v,parameters->nvar);
		if (COMMENT) cout << "magn = " << magn << endl;

		// normalise connecting vector
		vect_div_scal(conn_v, magn, parameters->nvar, conn_v);
		if (COMMENT)
			if (parameters->nvar==2) cout << "Normalised conn_v = (" << conn_v[0] << ", " << conn_v[1] << ")" << endl;
		magn = vect_magn(conn_v,parameters->nvar);
		if (COMMENT) cout << "magn = " << magn << endl;

		// multiply normalised connecting vector by an infinitesimal number...
		if (COMMENT) cout << "dv = " << dv << endl;
		vect_mult_scal(conn_v, dv, parameters->nvar, conn_v);
		if (COMMENT)
			if (parameters->nvar==2) cout << "Normalised conn_v multiplied by dv = (" << conn_v[0] << ", " << conn_v[1] << ")" << endl;
		magn = vect_magn(conn_v,parameters->nvar);
		if (COMMENT) cout << "magn = " << magn << endl;

		// finally, get final point
		vect_sum(initial, conn_v, parameters->nvar, final);
		if (COMMENT) if (parameters->nvar==2)  cout << "\n final = (" << final[0] << ", " << final[1] << ")" << endl;

		//evaluate tree in final point
		for (int k=0; k<parameters->nvar; k++) 			// assign the right value to all the variables
				(problem->v_list[k])->value = final[k];
		tree_final = tree_value(current_tree, NULL);    // NULL as n_corrections is not a matter here
		if (COMMENT) cout << "\nValue of the tree at the final point : tree_final = " << tree_final << endl;

		// estimate derivative in initial point along conn_v direction
		tree_der[row] = fin_differences(tree_initial, tree_final, dv);
		if (COMMENT) cout << "tree_der[" << row << "] = " << tree_der[row] << endl;

		// check tree expression for errors
		//if (COMMENT) {
			char* expr;
			expr = current_tree->print();
			cout << "\nCurrent tree:" << endl;
			cout << expr;	
			delete [] expr;
		//}
	//// end of the for cycle
	}
	
	// free memory
	delete [] initial;
	delete [] final;
	delete [] conn_v;
}


// function that assess if 1-order inequality constraints are satisfied
// (includes the evaluation of the tree derivative in a given point along a direction declared by a normalised vector)
// INPUT:
// 	ProblemDefinition  pb,
//		address of the node - Binary_Node* c_tree
// OUTPUT:
//		overall score of the tree with regard to 1st order inequality constraints
double Population::get_tree_derivative_given_norm_vector(ProblemDefinition pb, Binary_Node* c_tree)
{
	int COMMENT = 0;

	if (COMMENT) cout << "\nget_tree_derivative_given_norm_vector : enter" << endl;

	cout.precision(10);
	double epsilon = 1.0E-05;    // step for derivative
	double tree_initial, tree_final;
	int n_var = pb.get_n_var();
	double *p_initial = new double [n_var];
	double *p_final = new double [n_var];
	double *vect = new double [n_var];
	double* norm_vect = new double [n_var];
	double* dvect = new double [n_var];
	double magn;
	double tree_der;
	double der_error = 1e+10;
	double const1_error;
	double leak = 0.0;

	// evaluate tree partial derivatives in points and directions given in input file
	for (int row=0; row< pb.n_inequality1; row++) {

		// get coordinates of initial point and of normalised direction
		for (int k=0; k < n_var; k++) {
			p_initial[k] = pb.data_inequality1[row][k];
			vect[k] = pb.data_inequality1[row][k+n_var];
		}
		if (COMMENT) {
			cout << "\n" << row << ")";
			if (n_var==2)  {
				cout << " Computing tree derivative in point (" << p_initial[0] << ", " << p_initial[1] << ")" << endl;
				cout << "\nalong direction (" << vect[0] << ", " << vect[1] << ")" << endl;
			}
		}

		// compute connecting vector magnitude
		magn = vect_magn(vect,n_var);
		if (COMMENT) cout << "\nmagnitude connecting vector = " << magn << endl;

		// normalise connecting vector if its magnitude is not 1
		if (COMMENT) cout << "\nNormalising vector ...";
		vect_div_scal(vect, magn, n_var, norm_vect);
		if (COMMENT)
			if (n_var==2) cout << "\nNormalised vector = (" << norm_vect[0] << ", " << norm_vect[1] << ")" << endl;
		magn = vect_magn(norm_vect,n_var);
		if (COMMENT) cout << "magn = " << magn << endl;

		// multiply normalised connecting vector by an infinitesimal number...
		if (COMMENT) cout << "epsilon = " << scientific << epsilon << endl;
		vect_mult_scal(norm_vect, epsilon, n_var, dvect);
		if (COMMENT)
			if (n_var==2) cout << "dvect = (" << dvect[0] << ", " << dvect[1] << ")" << endl;
		magn = vect_magn(dvect,n_var);
		if (COMMENT) cout << "dvect magnitude = " << magn << endl;

		// get final point
		vect_sum(p_initial, dvect, n_var, p_final);
		if (COMMENT) if (n_var==2)  cout << "\n p_final = (" << p_final[0] << ", " << p_final[1] << ")" << endl;

		//evaluate tree in initial point
		for (int k=0; k<n_var; k++) 			// assign the right value to all the variables for the i-th sample case
				(*pb.v_list[k]).value = p_initial[k];
		tree_initial= tree_value(c_tree, NULL);    // NULL as n_corrections is not a matter here
		if (COMMENT) cout << "\nValue of the tree at the initial point : tree_initial = " << tree_initial << endl;

		//evaluate tree in final point
		for (int k=0; k<n_var; k++) 			// assign the right value to all the variables for the i-th sample case
				(*pb.v_list[k]).value = p_final[k];
		tree_final=  tree_value(c_tree, NULL);       // NULL as n_corrections is not a matter here
		if (COMMENT) cout << "\nValue of the tree at the final point : tree_final = " << tree_final << endl;

		// estimate derivative in initial point along dvect direction
		// ATTENTION! The positive direction in the domain is always the one of dvect!!!
		tree_der = fin_differences(tree_initial, tree_final, vect_magn(dvect,n_var));
		if (COMMENT) cout << "tree_der = " << tree_der << endl;
// OK so far....

		// sum the leaks (compare the tree derivative with the constraint1 (> or <))
		if (pb.constraints1[row] == '>') {
			if (tree_der <= pb.data_inequality1[row][2*n_var])
				// 1st approach: just count the unfeasible points
				//leak++;
				// 2nd approach: measure the distance from the feasible region
				leak = leak + fabs(pb.data_inequality1[row][2*n_var] - tree_der);
			}
		if (pb.constraints1[row] == '<') {
			if (tree_der  >= pb.data_inequality1[row][2*n_var])
				// 1st approach: just count the unfeasible points
				//leak++;
				// 2nd approach: measure the distance from the feasible region
				leak = leak + fabs(pb.data_inequality1[row][2*n_var] - tree_der);
		}

	}

	// free memory
	delete [] p_initial;
	delete [] p_final;
	delete [] vect;
	delete [] norm_vect;
	delete [] dvect;

	cout.precision(5);
	// return error value on first order partial derivative constraint
	if (COMMENT) cout << "\nleak = " << leak << endl;

	if (COMMENT) cout << "\nget_tree_derivative_given_norm_vector : exit" << endl;

 	return leak;
}



// function that recognizes if two trees are identical
// input: pointer to root node 1, pointer to root node 2
// output: int. 1 if trees are identical, 0 if are different
int Population::identical_trees(Binary_Node *p_tree1, Binary_Node* p_tree2)
{
	int COMMENT =0; //1 comments, 0 silent...
	char *expr1, *expr2; 
	int r;
	int size1, size2;
	size1 = p_tree1->count();
	size2 = p_tree2->count();

	if (COMMENT) {
		expr1 = p_tree1->print();
		expr2 = p_tree2->print();
		cout << "\nFirst tree: " << expr1;
		cout << "\nSecond tree: " << expr2;
		delete [] expr1;
		delete [] expr2;
	}

	if (size1 != size2) 
		// trees have different size, so they are different...
		return 0;

	// if you are here, trees have same size...
	r = identical_nodes((Node*)p_tree1, (Node*)p_tree2);
	
	// if r=1 trees are identical, if r=0 trees are different
	return r;
}
 
//function to compare two nodes
//input: pointer to nodes
// output: 1 if nodes are equal (type and content), 0 if they are not
int Population::identical_nodes(Node*p_n1, Node*p_n2)
{
	// check type
	if (p_n1->type != p_n2->type) return 0;
	
	// if you are here, the nodes are of the same type
	switch (p_n1->type) {
		
		case NODE_BINARY: {
			if (((Binary_Node*)p_n1)->get_func() != ((Binary_Node*)p_n2)->get_func())
				// different function 
				return 0;
			// if you are here, binary nodes have the same function: check children
			int rl = identical_nodes(((Binary_Node*)p_n1)->get_left(),((Binary_Node*)p_n2)->get_left());
			int rr = identical_nodes(((Binary_Node*)p_n1)->get_right(),((Binary_Node*)p_n2)->get_right());
			if ((!rl) || (!rr))
				// left children or right children are different
				return 0;
			//if you are here, left and right node are equal	: if (rl && rr)
			return 1;
		}
		break;

		case NODE_UNARY: {
			if (((Unary_Node*)p_n1)->get_func() != ((Unary_Node*)p_n2)->get_func())
				// different function 
				return 0;
			// if you are here, unary nodes have the same function: check children
			int r = identical_nodes(((Unary_Node*)p_n1)->get_child(),((Unary_Node*)p_n2)->get_child());
			
			if (!r)
				// different children
				return 0;	
			// identical children
			return 1;
		}
		break;

		case NODE_VAR: {
			if ( ((Terminal_Var*)p_n1)->val_p() != ((Terminal_Var*)p_n2)->val_p())
				// different variable 
				return 0;	
			//	identical variable
			return 1;
		}
		break;

		case NODE_CONST: {
			if ( ((Terminal_Const*)p_n1)->value(NULL) != ((Terminal_Const*)p_n2)->value(NULL))
				// different value 
				return 0;	
			//	identical value
			return 1;
		}
		
	 //end of switch
	}

}


// function that updates the number of constructive, destructive and neutral
// genetic operations for each genetic operator at the current generation and
// for the given tree.
// counter of destructive/neutral/constractive ops -> *_perf
// variable storing sum of fitness delta for destructive/neutral/constractive ops -> *_av_delta
// (called by Population::evaluate)
// INPUT: address of tree with fitness value already evaluated
void Population::measure_genetic_op_performance(Binary_Node* tree)
{
	int COMMENT = 0;


	// compute variation in error (accuracy)
	double delta = -(tree->fitness-tree->parent_fitness); // positive if fitness decreases

	// son_of => distinguish among:
	// 0 = reproduction, 1 = crossover, 2 = subtree mutation, 3 = point mutation

	switch (tree->son_of) {
		// check if constructive, neutral or destructive on the basis of eps_neutral
		// result of the genetic operations:
		// 0: destructive, 1: neutral, 2: constructive
		case 0:	{ // 0 = reproduction.  This is done because if copies are found new individuals are inserted...
			tot_repr++;
			if (delta>eps_neutral) {
				// constructive
				reproduction_perf[2]++;
				repr_av_delta[2] = repr_av_delta[2] - delta;
			}
			if (fabs(delta)<=eps_neutral) {
				// neutral
				reproduction_perf[1]++;
				repr_av_delta[1] = repr_av_delta[1] - delta;
			}
			if (delta<-eps_neutral) {
				// destructive
				reproduction_perf[0]++;
				repr_av_delta[0] = repr_av_delta[0] - delta;
			}
		}
		break;

		case 1:	{ // 1 = crossover
			tot_cross++;
			if (delta>eps_neutral) {
				// constructive
				crossover_perf[2]++;
				cross_av_delta[2] = cross_av_delta[2] - delta;
			}
			if (fabs(delta)<=eps_neutral) {
				// neutral
				crossover_perf[1]++;
				cross_av_delta[1] = cross_av_delta[1] - delta;
			}
			if (delta<-eps_neutral) {
				// destructive
				crossover_perf[0]++;
				cross_av_delta[0] = cross_av_delta[0] - delta;
			}
		}
		break;

		case 2:	{ // 2 = subtree mutation
			tot_smut++;
			if (delta>eps_neutral) {
				// constructive
				s_mutation_perf[2]++;
				smut_av_delta[2]=smut_av_delta[2]-delta;
			}
			if (fabs(delta)<=eps_neutral) {
				// neutral
				s_mutation_perf[1]++;
				smut_av_delta[1]=smut_av_delta[1]-delta;

			}
			if (delta<-eps_neutral) {
				// destructive
				s_mutation_perf[0]++;
				smut_av_delta[0]=smut_av_delta[0]-delta;
			}
		}
		break;

		case 3:	{ // 3 = point mutation performance
			tot_pmut++;
			if (delta>eps_neutral) {
				// constructive
				p_mutation_perf[2]++;
				pmut_av_delta[2]=pmut_av_delta[2]-delta;
			}
			if (fabs(delta)<=eps_neutral) {
				// neutral
				p_mutation_perf[1]++;
				pmut_av_delta[1]= pmut_av_delta[1]-delta;
			}
			if (delta<-eps_neutral) {
				// destructive
				p_mutation_perf[0]++;
				pmut_av_delta[0]=pmut_av_delta[0]-delta;
			}
		}
		break;
	}

}

// no training: perturbation of genetic operators rates when proportion of constructive crossovers 
// goes below 5%
void Population::adapt_genetic_operators_rates_notused(void)
{
	// version in which mutation is increased only when constructive crossover 
  // is under observed_var_thr for a given number of generations. After 3 generations
  // of alterated genetic operators rates everything goes back to original rates
  
  
  int COMMENT = 0;

	int adaptive_gen_ops = 0;  //switch for adaptive approach: 1 = ON, 0 = off

	if (!adaptive_gen_ops)
		return;

  int ngen_idle=5;
  int ngen_perturb=4;
	double observed_var; //variable used to monitor the evolution and adapt genetic ops rates
	double loc_observed_var_thr; //transformation of observed_var (now average)
	double rate_delta = .05;   // minimum increase of genetic ops rates
  double learning_bound = .9; // threshold that defines end of learning phase
	

	observed_var = (double)crossover_perf[2]/(double)tot_cross; // proportion of constructive crossovers
  observed_var_thr = .05;
  
  if (learning_window<=ngen_idle) {
     if (observed_var < observed_var_thr) {
          learning_window++;
     } else {
          learning_window=0;
     } 
  } else {
		if (COMMENT) {
			cout << "\n\nPopulation::adapt_genetic_operators_rates() : adapting phase";
		}
		// insert here the extreme gen. ops values to perturb the search 
  	repr_rate = 0.0;
		cross_rate = .0;
		mut_rate = 1.0 - (repr_rate+cross_rate);
    window_counter++;
    if (window_counter==ngen_perturb+1) {
       learning_window = 0;
       repr_rate = parameters->repr_rate; // repr_rate;
	     cross_rate = parameters->cross_rate;
     	 mut_rate = parameters->mut_rate;
       window_counter =0;
	  }
  }

	if (COMMENT) {
		cout << "\nlearning_window = " << learning_window;
		cout << "\ntot_cross = " << tot_cross;
		cout << "\nsum_observed_var = " << scientific << sum_observed_var;
		cout << "\nobserved_var = " << scientific << observed_var;
		cout << "\nobserved_var_thr =" << scientific << observed_var_thr;
		cout << "\nrepr_rate = " << repr_rate;
		cout << "\ncross_rate = " << cross_rate;
    cout << "\nmut_rate = " << mut_rate;
	}

}

// training stage followed by adaptive stage
void Population::adapt_genetic_operators_rates_notused2(void)
{
	int COMMENT = 0;

	int adaptive_gen_ops = 1;  //switch for adaptive approach: 1 = ON, 0 = off

	if (!adaptive_gen_ops)
		return;

	double observed_var; //variable used to monitor the evolution and adapt genetic ops rates
	double loc_observed_var_thr; //transformation of observed_var (now average)
	double rate_delta = .05;   // minimum increase of genetic ops rates
	//double learning_var=.0; //variable monitored to define learning phase
  double learning_bound = .95; // threshold that defines end of learning phase
	

	observed_var = (double)crossover_perf[2]/(double)tot_cross; // proportion of constructive crossovers
  //&/&/&/ test
  learning_var=100.0;
  observed_var_thr = .05;
  //&/&/&/
   
	// check if you are still in the initial set of generations used as learning phase
	//if (window_count<learning_window) {
  if (((learning_var<=learning_bound) && (learning_on)) || learning_window==0) {
		// update learning window counter (set to 0 by the Population constructor)
		learning_window++; 
    sum_learning_var = sum_learning_var + (double)crossover_perf[0]/(double)tot_cross; // proportion of destructive crossovers
    learning_var = sum_learning_var/learning_window;
		// rules to adapt genetic operators rates and reproduction rate
		// compute sum of constructive crossovers proportions
		sum_observed_var = sum_observed_var + observed_var; //n. of neutral crossovers. Add cross_av_delta as well? Could be a good idea
		if (learning_window)
			observed_var_thr = sum_observed_var/(double)learning_window;
		else {
			cerr << "\nPopulation::adapt_genetic_operators_rates() : ERROR : learning_window = 0! Exit.";
			exit(-1);
		}		

		if (COMMENT) {
			cout << "\n\nPopulation::adapt_genetic_operators_rates() : learning phase";
		}
	}

	//if (window_count>learning_window) {
	if (learning_var>learning_bound) {
		// here you have passed the learning period - so adaptive phase
		learning_on = 0;
		// insert here the rules to change genetic operators rates

		if (observed_var < observed_var_thr) {
			repr_rate = repr_rate - rate_delta;
			cross_rate = cross_rate - rate_delta;
			//mut_rate = mut_rate + .05;  //not used in new_spawn...
			if (cross_rate<0.05) {
				cross_rate=0.05;
			}
			if (repr_rate < 0.05) {
				repr_rate = 0.05;
			}
      
		}
		if (observed_var >= observed_var_thr) {
			repr_rate = repr_rate + rate_delta;
			cross_rate = cross_rate + rate_delta;
			//mut_rate = mut_rate - .05;  //not used in new_spawn...
			if (cross_rate>0.95) {
				cross_rate=0.95;
				//repr_rate=0.0;
			}
			if (repr_rate>0.95) {
				repr_rate=0.95;
				//cross_rate=0.1;
			}
		
		}

    //cout << "\nPRIMA repr_rate = " << repr_rate;
		//cout << "\nPRIMA cross_rate = " << cross_rate;
    double tot_repr_cross = repr_rate+cross_rate;
    if (tot_repr_cross>1.0) {
				cross_rate = cross_rate/tot_repr_cross; //1.0;
		    repr_rate = repr_rate/tot_repr_cross; //0.0;
    }
    mut_rate = 1.0 - (repr_rate+cross_rate);

		if (COMMENT) {
			cout << "\n\nPopulation::adapt_genetic_operators_rates() : adapting phase";
		}
	}

	if (COMMENT) {
		cout << "\nlearning_window = " << learning_window;
		cout << "\ntot_cross = " << tot_cross;
		cout << "\nsum_observed_var = " << scientific << sum_observed_var;
		cout << "\nobserved_var = " << scientific << observed_var;
		cout << "\nobserved_var_thr =" << scientific << observed_var_thr;
		cout << "\nrepr_rate = " << repr_rate;
		cout << "\ncross_rate = " << cross_rate;
    cout << "\nmut_rate = " << mut_rate;
	}

}




// training stage followed by adaptive stage
void Population::adapt_genetic_operators_rates(void)
{
	int COMMENT = 0;

	int adaptive_gen_ops = 0;  //switch for adaptive approach: 1 = ON, 0 = off

	if (!adaptive_gen_ops)
		return;

	double observed_var; //variable used to monitor the evolution and adapt genetic ops rates
	double loc_observed_var_thr; //transformation of observed_var (now average)
	double rate_delta = .05;   // minimum increase of genetic ops rates
	//double learning_var=.0; //variable monitored to define learning phase
  double learning_bound = .9; // threshold that defines end of learning phase
	

	observed_var = (double)crossover_perf[2]/(double)tot_cross; // proportion of constructive crossovers
  //&/&/&/ test values to remove training stage
  //learning_var=100.0;
  //observed_var_thr = .05;
  //&/&/&/
   
	// check if you are still in the initial set of generations used as learning phase
	//if (window_count<learning_window) {
  if (((learning_var<=learning_bound) && (learning_on)) || learning_window==0) {
		// update learning window counter (set to 0 by the Population constructor)
		learning_window++; 
    sum_learning_var = sum_learning_var + (double)crossover_perf[0]/(double)tot_cross; // proportion of destructive crossovers
    learning_var = sum_learning_var/learning_window;
		// rules to adapt genetic operators rates and reproduction rate
		// compute sum of constructive crossovers proportions
		sum_observed_var = sum_observed_var + observed_var; //n. of constructive crossovers. Add cross_av_delta as well? Could be a good idea
		if (learning_window)
			observed_var_thr = sum_observed_var/(double)learning_window;
		else {
			cerr << "\nPopulation::adapt_genetic_operators_rates() : ERROR : learning_window = 0! Exit.";
			exit(-1);
		}		

		if (COMMENT) {
			cout << "\n\nPopulation::adapt_genetic_operators_rates() : learning phase";
		}
	}

	//if (window_count>learning_window) {
	if (learning_var>learning_bound) {
		// here you have passed the learning period - so adaptive phase
		learning_on = 0;
		// insert here the rules to change genetic operators rates

		if (observed_var < observed_var_thr) {
			repr_rate = repr_rate - rate_delta;
			cross_rate = cross_rate - rate_delta;
			//mut_rate = mut_rate + .05;  //not used in new_spawn...
			if (cross_rate<0.05) {
				cross_rate=0.05;
			}
			if (repr_rate < 0.05) {
				repr_rate = 0.05;
			}
      
		}
		if (observed_var >= observed_var_thr) {
			repr_rate = repr_rate + rate_delta;
			cross_rate = cross_rate + rate_delta;
			//mut_rate = mut_rate - .05;  //not used in new_spawn...
			if (cross_rate>0.95) {
				cross_rate=0.95;
				//repr_rate=0.0;
			}
			if (repr_rate>0.95) {
				repr_rate=0.95;
				//cross_rate=0.1;
			}
		
		}

    //cout << "\nPRIMA repr_rate = " << repr_rate;
		//cout << "\nPRIMA cross_rate = " << cross_rate;
    double tot_repr_cross = repr_rate+cross_rate;
    if (tot_repr_cross>1.0) {
				cross_rate = cross_rate/tot_repr_cross; //1.0;
		    repr_rate = repr_rate/tot_repr_cross; //0.0;
    }
    mut_rate = 1.0 - (repr_rate+cross_rate);

		if (COMMENT) {
			cout << "\n\nPopulation::adapt_genetic_operators_rates() : adapting phase";
		}
	}

	if (COMMENT) {
		cout << "\nlearning_window = " << learning_window;
		cout << "\ntot_cross = " << tot_cross;
		cout << "\nsum_observed_var = " << scientific << sum_observed_var;
		cout << "\nobserved_var = " << scientific << observed_var;
		cout << "\nobserved_var_thr =" << scientific << observed_var_thr;
		cout << "\nrepr_rate = " << repr_rate;
		cout << "\ncross_rate = " << cross_rate;
    cout << "\nmut_rate = " << mut_rate;
	}

}
//&/

void Population::compute_statistics(void)
// function that computes basic ARCHIVE statistics for the CURRENT generation
// Fit (fitness initially meant RMSE): min, max, mean and unb. est. of the pop. variance on the archive
// size: minimum,  average, maximum in the archive
// depth: minimum,  average, maximum in the archive
// F value: minimum,  average, maximum, variance in the archive
{
	int COMMENT = 0;

	int trepr = composition[0]+composition[1];   //archive size (although composition[1] is no. of new individuals generated from scratch... correct to include them?)
	double *af = new double[trepr];

	//compute statistics for Fit: min, max, mean and unb. est. of the pop. variance evaluated on the ARCHIVE!!!
	if (COMMENT) cout << "\n\nError of the archive members";
	for (int i=0; i < trepr; i++) {
			af[i] = complete_trees[i]->fitness;
			if (COMMENT) cout << "\naf[" << i << "] = " << af[i];
	}
	double Fit[4];      //F[0]->Fmin  F[1]->ave   F[2]->max    F[3]->unb. est. of the population variance
	basic_stat_analysis(af,Fit,trepr);
	Fit_max = Fit[2];
	Fit_min = Fit[0];
	Fit_ave = Fit[1];
	Fit_var = Fit[3];
	if (COMMENT) {
		cout << "\nmax Fit_max = "<<Fit_max;
		cout << "\nmin Fit_min = "<<Fit_min;
		cout << "\nMean = " << Fit_ave;
		cout << "\nPop variance = " << Fit_var;
	}

	//size statistics: Smin=minimum    Save= average   Smax=maximum on the archive-------------------
	if (COMMENT) cout << "\n\nSize of the archive members";
	for (int i=0; i<trepr; i++) {
		af[i] = complete_trees[i]->count();
		if (COMMENT) cout << "\naf[" << i << "] = " << af[i];
	}
	double S[4];      //S[0]->Smin S[1]->Save   S[2]->Smax;
	basic_stat_analysis(af,S,trepr);
	S_max = S[2];
	S_min = S[0];
	S_ave = S[1];
	//S_var = S[3];
	if (COMMENT) {
		cout << "\nS_max = "<<S_max;
		cout << "\nS_min = "<<S_min;
		cout << "\nMean = " << S_ave;
		//cout << "\nPop variance = " << S_var;
	}

	//depth statistics: Dmin=minimum    Dave= average   Dmax=maximum   on the archive
	if (COMMENT) cout << "\n\nDepth of the archive members";
	for (int i=0; i<trepr; i++) {
			af[i] = complete_trees[i]->calc_depth();
			if (COMMENT) cout << "\naf[" << i << "] = " << af[i];
	}
	double D[4];      //D[0]->Dmin D[1]->Dave   D[2]->Dmax;
	basic_stat_analysis(af,D,trepr);
	D_max = D[2];
	D_min = D[0];
	D_ave = D[1];
	//D_var = S[3];
	if (COMMENT) {
		cout << "\nD_max = "<<D_max;
		cout << "\nD_min = "<<D_min;
		cout << "\nMean = " << D_ave;
		//cout << "\nPop variance = " << D_var;
	}

	//F value statistics: F_min=minimum    F_ave= average    F_max=maximum-------------------
	if (COMMENT) cout << "\n\nF of the archive members";
	for (int i=0; i<trepr; i++) {
		af[i] = complete_trees[i]->F;
		if (COMMENT) cout << "\naf[" << i << "] = " << af[i];
	}
	double F[4];      //F[0]->Fmin F[1]->ave   F[2]->max    F[3]->unbiased estimate of the population variance;
	basic_stat_analysis(af,F,trepr);
	F_max = F[2];
	F_min = F[0];
	F_ave = F[1];
	F_var = F[3];
	if (COMMENT) {
		cout << "\nF_max = "<<F_max;
		cout << "\nF_min = "<<F_min;
		cout << "\nMean = " << F_ave;
		cout << "\nPop variance = " << F_var;
	}

	delete[] af;

	//pen_ord1 statistics: F_min=minimum    F_ave= average    F_max=maximum-------------------
	double C1_ave = 0.0;

	//compute statistics for C1: min, max, mean and unb. est. of the pop. variance evaluated on the archive
	//mean
	for (int i=0; i<trepr; i++)
		C1_ave += complete_trees[i]->pen_ord1;

	C1_ave = C1_ave/(double)trepr;

	pen_ord1_ave = C1_ave;


	// compute statistics for each single objective of aggregate F evaluated on the ARCHIVE!!!
	trepr=composition[0];
	if (COMMENT) cout << "\n\nError of the archive members";
	double a,b,c;
	T_perc_ave[0]=1.0e+10;
	Fc_perc_ave[0]=1.0e+10;
	Fc_ave[0]=1.0e+10;
	// objective cycle
	for (int j=1; j < 12; j++) {   // objective j
		// trees cycle, same objective
		a=0.0;
		b=0.0;
		c=0.0;
		for (int i=0; i < trepr; i++) {
			a=a+100*(complete_trees[i]->T[j])/(complete_trees[i]->F);  // percentual contribution (includes weights)
			b=b+100*(complete_trees[i]->Fc[j])/(complete_trees[i]->F); // percentual contribution (does NOT include weights)
			c=c+(complete_trees[i]->Fc[j]);  // contribution (does NOT include weights)
		}
		// average percentage contribution relative to objective j throughout the archive trepr
		T_perc_ave[j]=a/(double)trepr;  // includes Fweight[i]
		// average percentage contribution relative to objective j throughout the archive trepr
		Fc_perc_ave[j]=b/(double)trepr;  // does NOT include Fweight[i]
		// average contribution relative to objective j throughout the archive trepr
		Fc_ave[j]=b/(double)trepr;  // does NOT include Fweight[i]
	}
}



/*
void Population::inherit_parameters(Binary_Node *p_tree1, Binary_Node *p_tree2)
{
	int COMMENT=1;   //1 comments, 0 silent...

	// sort new
	qsort(new_trees,size,sizeof(Binary_Node *),tree_comp);
	int r = identical_trees(new_trees[0], trees[0]);
}
*/

// Function that updates the external archive, storing the best tree found in the current generation (and its parameterless version)
void Population::update_ext_archive(void)
{
	int COMMENT=0;   //1 comments, 0 silent...
	static int n_times=0;
	int n_param;
	char *expr;

	if (COMMENT)  cout << "\nPopulation::update_ext_archive" << endl;
	// copy best structure and complete tree.
	// really IMPORTANT: trees and complete_trees must have been already sorted according to F
	// delete old individual
	if (best_tree)
		delete best_tree;
		//extfitness
		// if best_tree exists, update if it is better than the previous one
		//if (trees[0]->fitness < best_tree->fitness) {   //if you want to use fitness, sort according to fitness first!!
		//	if (COMMENT) cout << "\ntrees[0].fitness < best_tree.fitness : SUBSTITUTION... ";

	// copy the parameterless tree
	best_tree = (Binary_Node *)tree_copy(trees[0],NULL);
	best_tree->F = trees[0]->F;
	// copy the full tree (numerical coefficients included)
	best_complete_tree = (Binary_Node *)tree_copy(complete_trees[0],NULL);
	best_complete_tree->F = trees[0]->F;


	// fetch the addresses of the parameters and assign p_par in best_complete_tree
	n_param = best_complete_tree->find_parameters();  //calls int Binary_Node::find_parameters(void)
	
	if (COMMENT) {
		expr = complete_trees[0]->print();	
		cout << "\ncomplete_trees[0] : " << expr;
		delete [] expr;
		expr = trees[0]->print();	
		cout << "\ntrees[0] : " << expr;
		delete [] expr;
		expr = best_complete_tree->print();	
		cout << "\nbest_complete_tree : " << expr;
		delete [] expr;
		expr = best_tree->print();	
		cout << "\nbest_tree : " << expr;
		delete [] expr;
		Val node_value;
		cout << endl;
		for (int i=0; i<n_param; i++) {
			node_value = ((Terminal_Const *)best_complete_tree->p_par[i])->value(NULL);
			cout << "value [" << i << "] = " << node_value << endl;
		}
	}

	n_times++;
}




//// OLD VERSION (before 7/12/20), which finds all the pulsations corresponding to all single sin() and cos() terms
//// function to find the addresses of the terminal_const nodes that are "pulsation" or circular frequencies to be optimised in a specific way
//// This function is called by Population::tuning_individual
//// input: nothing
//// output: no of pulsations (vector_puls is updated with the addresses of the pulsations)
void Population::find_pulsations(Binary_Node *tree)
{
	int COMMENT = 0;

	tree->n_pulsations = 0;  // n_pulsations is a member of Binary_Node

	if (COMMENT) {
		cout << "\n\nPopulation::find_pulsations() : START";
		char *expr = tree->print();
		cout << "\ntree : " << expr;
		delete [] expr;
		cout << "\n index_puls = " << tree->index_puls;
		cout << "\n n_pulsations = " << tree->n_pulsations;
	}

	// exit if sin and cos are not primitives in input file
	if (COMMENT) cout << "\n &Sin = " << &Sin << "  &Cos = " << &Cos;
	if ((&Sin==NULL) && (&Cos==NULL)) {
		tree->n_pulsations = 0;
		if (COMMENT) cout << "\nfind_pulsations : n_pulsations = 0";
	}
	// if you are here, sin and cos are among primitives
	else {
		if (COMMENT) cout << "\n Sin or Cos are used";
		Node *grandpa;
		Node *brother;
		string grandpa_fun;
		vector < int > vector_puls; // vector containing the indexes of the pulsations in p_par
		vector < int > vector_var; // vector containing the indexes of the variables corresponding to pulsations
		int i_var;

		//searching among parameters' addresses in p_par
		for (int i=0; i< tree->p_par.size(); i++) {
				// retrieve grandparent's address: if a is the parameter in sin(a*x), then sin is the grandpa unary function, x is the brother
				grandpa = ((tree->p_par[i])->get_parent())->get_parent();
				// retrieve brother's address (parent's right child...)
				brother = ((Binary_Node*)(tree->p_par[i])->get_parent())->get_right();

				if (grandpa)  {   //make sure that grandpa exists (not the case for additional shift term...)
					if (grandpa->type == NODE_UNARY) {

						// check grandpa function
						grandpa_fun = ((Unary_Node*)grandpa)->get_func()->pre_sign;  //31/10/20 improvement: better to use a internal code for operations
						// check operations of all grandpa's ancestors

						// check that the parameter has all the features to be considered a "pulsation"
						if (((!strcmp(grandpa_fun.c_str(),"sin")) || (!strcmp(grandpa_fun.c_str(),"cos"))) && (brother->type==NODE_VAR)) {
							// check that all the functions upstream the multplication of a1*sinusoidal f. are add,sub,shift,neg
							if ((grandpa->get_parent())->check_nodal_functions_upstream()==1) {
								//&/&/&/&/&/&/&/&/&/&/
								// here you are sure that the pulsation belongs to a term like sin(omega*Zi) ...
								// PULSATION : append the index i of the pulsation to vector_puls
								vector_puls.push_back(i);
								tree->n_pulsations = tree->n_pulsations + 1; // added 7/12/20

								// VARIABLE INDEX in v_list : append the variable no. of the pulsation to vector_var
								// correct the number to get the index in v_list first
								i_var = (((Terminal_Var*)brother)->get_var_number()) - 1;
								//check that i_var< num_vars
								if (i_var>=parameters->nvar) {
									cerr << "\nPopulation::find_pulsations : ERROR : i_var >= num_vars !!!";
									exit (-1);
								}
								// update vector_var with the index
								vector_var.push_back(i_var);

								//&/&/&/&/&/&/&/&/&/&/
							}
						}
					}
				}

		}


		// transform vcout << "Index of pulsations in p_par: " << endl;ector into a standard int array
		tree->n_pulsations = vector_puls.size();
		if (!(tree->n_pulsations)) {
			if (COMMENT) cout << "\nfind_pulsations : n_pulsations = 0";
		}
		else	{  // here you are sure there are pulsations
			// initialise arrays with pulsation information
			tree->index_puls = new int[tree->n_pulsations];   // content deleted at the end of Population::tuning_individual
			tree->index_var = new int[tree->n_pulsations];	// content deleted at the end of Population::tuning_individual
			for (int k=0; k < tree->n_pulsations; k++) {
				tree->index_puls[k] = vector_puls[k];
				tree->index_var[k] = vector_var[k];
			}


			if (COMMENT) {
				cout << "\nn_pulsations = " << tree->n_pulsations << endl;
				cout << "Pulsations found:";
				for (int h=0; h < tree->n_pulsations; h++)
					cout << " " <<  (tree->p_par[tree->index_puls[h]])->value(NULL);
				cout << "\nCorresponding variable's index in v_list:"<< endl;
				for (int h=0; h < tree->n_pulsations; h++)
					cout << " " <<  tree->index_var[h];
			}
		}  // end else

	}  // end else

	if (COMMENT) cout << "\nPopulation::find_pulsations : END" << endl;
	//cin.get();
}


//+-+-+-
// NEW VERSION (after 7/12/20), finds only one pulsation corresponding to a single sin() or cos() term
// function to find the addresses of the terminal_const nodes that are "pulsation" or circular frequencies to be optimised in a specific way
// This function is called by Population::tuning_individual
// input: nothing
// output: no of pulsations (vector_puls is updated with the addresses of the pulsations)
//void Population::find_pulsations(Binary_Node *tree)
//{
//	int COMMENT = 0;
//
//	tree->n_pulsations = 0;  // n_pulsations is a member of Binary_Node
//
//	if (COMMENT) {
//		cout << "\n\nPopulation::find_pulsations";
//		char *expr = tree->print();
//		cout << "\ntree : " << expr;
//		delete [] expr;
//		cout << "\n index_puls = " << tree->index_puls;
//		cout << "\n n_pulsations = " << tree->n_pulsations;
//	}
//
//	// exit if sin and cos are not primitives in input file
//	if (COMMENT) cout << "\n &Sin = " << &Sin << "  &Cos = " << &Cos;
//	if ((&Sin==NULL) && (&Cos==NULL)) {
//		tree->n_pulsations = 0;
//		if (COMMENT) cout << "\nPopulation::find_pulsations : n_pulsations = 0 as sin and cos are not among primitives. Exit";
//		// exit, nothing to find if sin() and cos() are not among primitives
//		return;
//	}
//
//	// search for pulsations in sin() and cos() terms, as sin and cos are among primitives
//	if (COMMENT) cout << "\n Sin or Cos are used";
//	Node *grandpa;
//	Node *brother;
//	string grandpa_fun;
//	vector < int > vector_puls; // vector containing the indexes of the pulsations in p_par
//	vector < int > vector_var; // vector containing the indexes of the variables corresponding to pulsations
//	int i_var;
//
//	//searching for a single pulsation among parameters' addresses in p_par: as soon as it is found leave the loop
//	int i=0;
//	while ( (i< tree->p_par.size()) && (tree->n_pulsations==0) ) {
//
//
//			// retrieve grandparent's address: if a is the parameter in sin(a*x), then sin is the grandpa unary function, x is the brother
//			grandpa = ((tree->p_par[i])->get_parent())->get_parent();
//			// retrieve brother's address (parent's right child...)
//			brother = ((Binary_Node*)(tree->p_par[i])->get_parent())->get_right();
//
//			if (grandpa)  {   //make sure that grandpa exists (not the case for additional shift term...)
//				if (grandpa->type == NODE_UNARY) {
//
//					// check grandpa function
//					grandpa_fun = ((Unary_Node*)grandpa)->get_func()->pre_sign;  //31/10/20 improvement: better to use a internal code for operations
//					// check operations of all grandpa's ancestors
//
//					// check that the parameter has all the features to be considered a "pulsation"
//					if (((!strcmp(grandpa_fun.c_str(),"sin")) || (!strcmp(grandpa_fun.c_str(),"cos"))) && (brother->type==NODE_VAR)) {
//						// check that all the functions upstream the multplication of a1*sinusoidal f. are add,sub,shift,neg
//						if ((grandpa->get_parent())->check_nodal_functions_upstream()==1) {
//							//&/&/&/&/&/&/&/&/&/&/
//							// here you are sure that the pulsation belongs to a term like sin(omega*Zi) ...
//							// PULSATION : append the index i of the pulsation to vector_puls
//							vector_puls.push_back(i);
//							tree->n_pulsations = tree->n_pulsations + 1; // added 7/12/20
//
//							// VARIABLE INDEX in v_list : append the variable no. of the pulsation to vector_var
//							// correct the number to get the index in v_list first
//							i_var = (((Terminal_Var*)brother)->get_var_number()) - 1;
//							// check that i_var< num_vars
//							if (i_var>=parameters->nvar) {
//								cerr << "\nPopulation::find_pulsations : ERROR : i_var >= num_vars !!!";
//								exit (-1);
//							}
//							// update vector_var with the index
//							vector_var.push_back(i_var);
//
//						}
//					}
//				}
//			}
//	i++;
//	} // end while loop
//
//
//	// transform vcout << "Index of pulsations in p_par: " << endl;ector into a standard int array
//	tree->n_pulsations = vector_puls.size();
//	if (!(tree->n_pulsations)) {
//		if (COMMENT) cout << "\nfind_pulsations : n_pulsations = 0";
//	}
//	else	{  // here you are sure there are pulsations
//		// initialise arrays with pulsation information
//		tree->index_puls = new int[tree->n_pulsations];   // content deleted at the end of Population::tuning_individual
//		tree->index_var = new int[tree->n_pulsations];	// content deleted at the end of Population::tuning_individual
//		for (int k=0; k < tree->n_pulsations; k++) {
//			tree->index_puls[k] = vector_puls[k];
//			tree->index_var[k] = vector_var[k];
//		}
//
//
//		if (COMMENT) {
//			cout << "\nn_pulsations = " << tree->n_pulsations << endl;
//			cout << "Pulsations found:";
//			for (int h=0; h < tree->n_pulsations; h++)
//				cout << " " <<  (tree->p_par[tree->index_puls[h]])->value(NULL);
//			cout << "\nCorresponding variable's index in v_list:"<< endl;
//			for (int h=0; h < tree->n_pulsations; h++)
//				cout << " " <<  tree->index_var[h];
//		}
//	}  // end else
//
//
//	if (COMMENT) cout << "\nPopulation::find_pulsations : exit" << endl;
//
//	return;
//}
//+-+-+-


// function that searches the position (number of node) of the node containing a division
// called by Population::evaluate to enable factorisation bonus
void Population::search_first_op(Binary_Node *tree, Node *cur_node, int cur_depth)
{
  int COMMENT = 0;

  // list of operations whose presence has to be checked
  // BINARY
  int found_div =-1; // division
  int found_mult =-1; // multiplication
  // UNARY
  //int found_inv = -1; // reciprocal
  int found_sin = -1; // sin
  int found_cos = -1; // cos 
  int found_exp = -1; // exponential
  //int found_nxp = -1; // negative exponential
  int found_tanh = -1; //tanh
  int found_log = -1; //log

  if (cur_node->type == NODE_BINARY) {
	// check division
	found_div =  strcmp((((Binary_Node*)cur_node)->get_func())->sign,SDiv.sign);
	// check multiplication
	found_mult = strcmp((((Binary_Node*)cur_node)->get_func())->sign,Mult.sign);
   
	if ((found_div==0) || (found_mult==0)) { //or (!found_mult))
		tree->depth_first_op = cur_depth;
     if (COMMENT) {
    	 char *expr;
    	 expr = tree->print();
		 cout << "\n\n Tree: \n" << expr;
		 delete [] expr;
		 cout << "\n Found sought operation at depth : " << cur_depth;
		 if (found_div==0) cout << " - DIVISION";
		 if (found_mult==0) cout << " - MULTIPLICATION";
     }
			return;
		}
		else {
			cur_depth++;
			search_first_op(tree,((Binary_Node*)cur_node)->get_left(),cur_depth);
			search_first_op(tree,((Binary_Node*)cur_node)->get_right(),cur_depth);
		}
	}

	if (cur_node->type == NODE_UNARY) {
		// check reciprocal 
    //found_inv = strcmp((((Unary_Node*)cur_node)->get_func())->sign,Inverse.sign);
		// check sin
    found_sin = strcmp((((Unary_Node*)cur_node)->get_func())->pre_sign,Sin.pre_sign);
    // check cos
    found_cos  = strcmp((((Unary_Node*)cur_node)->get_func())->pre_sign,Cos.pre_sign);
    // check exp
    found_exp = strcmp((((Unary_Node*)cur_node)->get_func())->pre_sign,Exp.pre_sign);
    // check nxp
    //found_nxp = strcmp((((Unary_Node*)cur_node)->get_func())->sign,Negexp.sign);
    // check tanh
    found_tanh = strcmp((((Unary_Node*)cur_node)->get_func())->pre_sign,Tanh.pre_sign);
    // check log   
    found_log = strcmp((((Unary_Node*)cur_node)->get_func())->pre_sign,Logn.pre_sign);
    
    //if ((found_inv==0) || (found_sin==0) || (found_cos==0)) {
		if ((found_sin==0) || (found_cos==0) || (found_exp==0) ||  (found_tanh==0) || (found_log==0)) { //|| (found_nxp==0)
    	tree->depth_first_op = cur_depth;
      if (COMMENT) {
        char *expr;
		expr = tree->print();
		cout << "\n\n Tree: \n" << expr;
		delete [] expr;
		cout << "\n Found sought operation at depth : " << cur_depth;
        if (found_sin==0) cout << " - SIN";
        if (found_cos==0) cout << " - COS";
        if (found_exp==0) cout << " - EXP";
      }
			return;
		}
		else {
			cur_depth++;
			search_first_op(tree,((Unary_Node*)cur_node)->get_child(),cur_depth);
		}
	}

	return;
}

// function that seraches for diverging polynomials (22/10/20 to improve extrapolation properties)
void Population::search_high_level_polynomials(Binary_Node *tree, Node *cur_node)
{
	int COMMENT = 0;
	if (COMMENT) {
		cout << "\n\nPopulation::search_high_level_polynomials() : START";
		char *expr;
		expr = cur_node->print();
		cout << "\n Node: " << expr;
		delete [] expr;
	}

	// Binary node
	int found_add=0;
	int found_sub=0;
	int found_mult=0;
	// Unary node
	int found_shift=0;
	// Terminal node (variable)
	int found_variable=0;

	// check operation of root node
	// Binary node
	if (cur_node->type == NODE_BINARY) {
		// check addition
		if (!strcmp((((Binary_Node*)cur_node)->get_func())->sign,Add.sign)) found_add=1;
		// check subtraction
		if (!strcmp((((Binary_Node*)cur_node)->get_func())->sign,Sub.sign)) found_sub=1;
		// check multiplication
		if (!strcmp((((Binary_Node*)cur_node)->get_func())->sign,Mult.sign)) found_mult=1;
	}
	// Unary node: check only shift: if there are polynomial down below, they are argument of this unary function, and the user should avoid using diverging functions if bounded behaviour is requested...
	if (cur_node->type == NODE_UNARY) {
		// check shift
		if (!strcmp((((Unary_Node*)cur_node)->get_func())->post_sign, Shift.post_sign)) found_shift=1;
	}
	// Terminal Variable
	if (cur_node->type == NODE_VAR) {
		found_variable = 1;
	}
	// Terminal Const: by definition a constant, so no diverging issues


	// check in branches of high level polynomials
	// Binary Node
	if ((found_add) || (found_sub) || (found_mult)) {
		// call recursively offspring function
		search_high_level_polynomials(tree, ((Binary_Node*)cur_node)->get_left());
		search_high_level_polynomials(tree, ((Binary_Node*)cur_node)->get_right());
	}
	// Unary Node
	if (found_shift) {
		// call recursively offspring function
		search_high_level_polynomials(tree, ((Unary_Node*)cur_node)->get_child());
	}
	// Terminal Node - Variable
	if (found_variable) {
		tree->diverging=1;   // here a high level, diverging polynomial has been found
		if (COMMENT) cout << "\n\nPopulation::search_high_level_polynomials() : high level, diverging polynomial has been found";
	}


	if (COMMENT) {
		cout << "\n\nPopulation::search_high_level_polynomials() : END";
	}


	return;

}


// function that purges high level polynomials in trees without parameters (29/11/20)
// each branch of a tree that might contain a diverging term is analysed
// if divergent subtrees are found, two alternative behaviours can be adopted:
// 1. the terminal variable of the subtree is turned into a constant
// 2. the terminal variable of the subtree is replaced by sin(variable) (which is bounded!)
void Population::purge_diverging_terms(Binary_Node *tree, Node *cur_node)
{

	int COMMENT = 0;
	if (COMMENT) {
		cout << "\n\nPopulation::purge_high_level_polynomials() : START";
		char *expr;
		expr = cur_node->print();
		cout << "\n Node: " << expr;
		delete [] expr;
	}

	// Binary node
	int found_add=0;
	int found_sub=0;
	int found_mult=0;
	int found_spow=0;
	// Unary node
	int found_shift=0;
	int found_exp=0;
	int found_square=0;
	int found_cube=0;
	// Terminal node (variable)
	int found_variable=0;

	// check operation of root node
	// Binary node
	if (cur_node->type == NODE_BINARY) {
		// check addition
		if (!strcmp((((Binary_Node*)cur_node)->get_func())->sign,Add.sign)) found_add=1;
		// check subtraction
		if (!strcmp((((Binary_Node*)cur_node)->get_func())->sign,Sub.sign)) found_sub=1;
		// check multiplication
		if (!strcmp((((Binary_Node*)cur_node)->get_func())->sign,Mult.sign)) found_mult=1;
		// check spow
	}
	// Unary node
	if (cur_node->type == NODE_UNARY) {
		// check shift
		if (!strcmp((((Unary_Node*)cur_node)->get_func())->post_sign, Shift.post_sign)) found_shift=1;
		// check exp
		if (!strcmp((((Unary_Node*)cur_node)->get_func())->pre_sign, Exp.pre_sign)) found_exp=1;
		// check square
		if (!strcmp((((Unary_Node*)cur_node)->get_func())->post_sign, Square.post_sign)) found_square=1;
		// check cube
		if (!strcmp((((Unary_Node*)cur_node)->get_func())->post_sign, Cube.post_sign)) found_cube=1;

	}
	// Terminal Variable
	if (cur_node->type == NODE_VAR) {
		found_variable = 1;
	}
	// Terminal Const: by definition a constant, so no diverging issues


	// check in branches of high level polynomials
	// Binary Node
	if ((found_add) || (found_sub) || (found_mult)) {
		// call recursively offspring function
		purge_diverging_terms(tree, ((Binary_Node*)cur_node)->get_left());
		purge_diverging_terms(tree, ((Binary_Node*)cur_node)->get_right());
	}
	// Unary Node
	if (found_shift || found_exp || found_square || found_cube) {
		// call recursively offspring function
		purge_diverging_terms(tree, ((Unary_Node*)cur_node)->get_child());
	}
	// Terminal Node - Variable
	if (found_variable) {
		// 1. here the variable found is the cause of the divergent behaviour: replace it with a constant (1) - TOO CONSERVATIVE!!
//		eliminate_subtree(cur_node, (Val)1.0);
//		if (COMMENT) {
//			cout << "\n\nPopulation::purge_diverging_terms() : diverging polynomial has been found - replaced by 1.0";
//			char *expr;
//			expr = cur_node->print();
//			cout << "\n Node: " << expr;
//			delete [] expr;
//		}
		// 2. here the variable is replaced by a sin function of the same variable
		// create "sin" node and attach child variable randomly selected
		//REPLACE DIVERGENT VARIABLE
		Unary_Func *f=problem->u_func_list[int_rand(problem->num_u_funcs)]; // if this line works, a random function among unary functions is selected
		Unary_Node *new_edsin = new Unary_Node(NULL,&Sin);  //f);  //&Sin);

		Terminal_Var *new_edvar = new Terminal_Var((Node*)new_edsin,problem->v_list[int_rand(parameters->nvar)]);
		new_edsin->set_child((Node*)new_edvar);
		if (COMMENT) {
			cout << "\n\nPopulation::purge_diverging_terms() : new sin() subtree";
			cout << "\n new_edsin = " << new_edsin;
			cout << "\n new_edvar = " << new_edvar;
			char *expr;
			expr = new_edsin->print();
			cout << "\n Subtree: " << expr;
			delete [] expr;
		}
		// replace the divergent cur_node with newly created subtree
		int i=0;
		Node* pa_cur_node = cur_node->get_parent();
		i=insert_subtree((Node**)&tree, (Node*)cur_node, (Node*)new_edsin);
		// ok till here! Mind that here cur_node does not exist any longer, as it is deleted by insert_subtree!
		if (COMMENT) {
			cout << "\n\nPopulation::purge_diverging_terms() : diverging polynomial has been found - replaced by sin()";
			char *expr;
			expr = pa_cur_node->print();
			cout << "\n Node: " << expr;
			delete [] expr;
		}


	}


	if (COMMENT) {
		cout << "\n\npurge_diverging_terms() : END";
	}
//cin.get();
	return;
}


Node* Population::remove_unary_node(Unary_Node *u_node)
{
	// function to bypass a node. It removes a unary node.
    // INPUT		: node to be removed (must be a unary node)
	// OUTPUT	: address of the child of the submitted node if removal is performed
	//			       	: None if removal is not performed
    
    // address of the unary node's child
    Node *r = u_node->get_child();
    
    Node *father = u_node->get_parent();
    
    if (u_node->get_parent()) {
        // update parent's child
        // binary node
        if (father->type == NODE_BINARY) {
        	//(Binary_Node *)father = (Binary_Node*)(u_node->parent);
        	if (((Binary_Node *)father)->get_left() == (Node*)u_node)
            	// left child
            	((Binary_Node *)father)->set_left(r);
            else
                // right child
               ((Binary_Node *)father)->set_right(r);
        }
       	// unary node
        if (father->type == NODE_UNARY) {
            ((Unary_Node*)father)->set_child(r);
        }
        // update child's parent
        r->set_parent(father);
    }        
    else {
        // current node is root node - so far impossible, as unary nodes can't be root nodes 
        // update child's parent
        (u_node->get_child())->set_parent(NULL);
        //new_root = u_node->child;
    }
        
    // cut connections of the current node and delete it
    u_node->set_parent(NULL);
    u_node->set_child(NULL);
    delete u_node;
    
    // return the address of the bypassed node's child     
    return r;	
}



// function to generate a random number (random ephemeral constant, as Koza calls it...)
Val Population::constant_generation(Val vmin, Val vmax, Val step)
{
	// Rational random number is generated dividing two integer random numbers
	// as the sequence of random numbers is not actually random, if you seed the random generator (srand) each time you want
	// a random number you will get the first element of the pseudo-random list, that is you will have every time the same 
	// NUMBER!!! So seed the generator just once at the beginning of the main program!
	
	// to select a point between a (min) and b (max), where the increment is dz follow the linear function:
	// y = min +step * x
	// where step = (max-min)/(n-1), width of the range divided in n-1 parts (n points)
	int n = (int)((vmax-vmin)/step) + 1;
	Val c = vmin + step*int_rand(n);
	//cout << c << " ";
	return c;
}




// function to generate a random number between two limits (overloaded version of constant_generation())
Val Population::constant_generation(Val vmin, Val vmax, Val step, Val a, Val b)
{
	// Rational random number is generated dividing two integer random numbers
	// as the sequence of random numbers is not actually random, if you seed the random generator (srand) each time you want
	// a random number you will get the first element of the pseudo-random list, that is you will have every time the same 
	// NUMBER!!! So seed the generator just once at the beginning of the main program!
	
	// to select a point between a (min) and b (max), where the increment is "step" follow the linear function:
	// y = a +dz * x   
	// where step = (b-a)/(n-1), width of the range divided in n-1 parts (n points)
	
	// force the extremes to be in the range defined in the input file [MINRAND, MAXRAND]
	if (a<vmin)
		a=vmin;
	if (b>vmax)
		b = vmax;
	
	//
	int n = (int)((b-a)/step) + 1;
	Val c = a + step*int_rand(n);
	
	return c;
}



void Population::compute_selected_nodes_statistics(Binary_Node* p_tree, int node_num)
{
	int node_depth = 0;
	int node_type = 0;
	Node* current_parent;
	Node* p_node;
	
	// retrieve the pointer to the node and the node type
	p_node = p_tree->find(node_num);
	node_type = p_node->type;
	
	// compute depth of the selected node
	current_parent = p_node->get_parent();
	while (current_parent !=NULL) {
		node_depth++;
		current_parent = current_parent->get_parent();
	}
	
	// now node depth is known. Update counters:
	total_nodes_selected++;
	selected_nodes_per_depth[node_depth]++;
	selected_nodes_per_type[node_type]++; 
}


//*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/**/*/*/*/*/*/*/*/*/*/*/*/*/

// function that allocate a double array (matrix)
// input: no of rows, no of columns
// output: address of the double array (pointer)
double** Population::zeros(int nrow, int ncol)
{
	double** M;
	M = new double*[nrow]; // dynamic array of pointers: number of rows specified first
	for (int i=0; i<nrow; i++) {
		M[i] = new double[ncol];
		// set each entry equal to 0
		for (int j=0; j<ncol; j++)
			M[i][j] = 0.0;
	}

	return M;
}


// function that allocate a single array (vector)
// input: no of entries
// output: address of the single array (pointer)
double* Population::zeros(int nentries)
{
	double* v;
	// allocate
	v = new double[nentries];
	if (!v) {
		// output error message
	    cerr << "zeros(int nentries) : ERROR ! Can't allocate v pointer to double[nentries]!\n";
	    exit(-1);
	}

	// initialise
	for (int i=0; i<nentries; i++) {
		// set each entry equal to 0
		v[i]=0.0;
	}

	return v;
}



int Population::regenerate_on_boundary(double** p_coord, int spacedim, double** bounds)
{
	/*
		Function that regenerate a point in case it is out of bounds.
		The unfeasible point is regenerated on the boundary.
		Inputs:
		- p_coord    point coordinates - array spacedim x 1
		- spacedim  number of dimensions
		- bounds    array with lower and upper bounds (array spacedim x 2)
		Output:
		- p_coord    new feasible coordinates of the new point on the design space boundary
	*/

	// check current position - array
	for (int i=0; i<spacedim; i++) {
    	// check single coordinate is inside corresponding bounds
        if (*(p_coord[i]) > bounds[i][1])
        	// p_coord[i] exceeds upper bound
        	*(p_coord[i]) = bounds[i][1];
        if (*(p_coord[i]) < bounds[i][0])
			// p_coord is below lower bound
        	*(p_coord[i]) = bounds[i][0];

	}
    return 1;

}




int Population::regenerate_from_swarm_best(double** p_coord, int spacedim, int nparticles, double** bounds, double** zLocal)
{
    /*
        Function that regenerate a particle in case it is out of bounds.
        The unfeasible coordinates of the particle are selected randomly
        from the coordinates of the swarm best particles (see Xu11).
        Inputs:
        - p_coord  	coordinates of the unfeasible particle - array spacedim x 1
        - bounds   	array with lower and upper bounds (array spacedim x 2)
        - zLocal   	matrix containing the coordinates of the personal best of each particle so far (n_var x Npop - 2D even for n_var =1))
        Output:
        - p_coord    new feasible coordinates of the new point on the design space boundary
        - number of particle chosen for the regeneration
    */

	int s_particle=0;

    // check current position - array
	for (int i=0; i<spacedim; i++) {
    	// check single coordinate is inside corresponding bounds
    	if ((*(p_coord[i]) < bounds[i][0]) || (*(p_coord[i]) > bounds[i][1])) {
        	// pick the corresponding coordinate of a randomly selected particle among best so far
            s_particle = rand() % nparticles; // random integer between 0 and nparticles-1
            // assign the new feasible coordinates
            for (int i=0; i<spacedim; i++) *(p_coord[i])=zLocal[i][s_particle];
    	}
    }


	return s_particle;

}



/*
function to perform Particle Swarm Optimisation (minimisation of the given objective function)
Inputs:
- ntree 		pointer to root node of the tree whose numerical coefficients have to be optimised
- p_objfun     	pointer to objective function whose minima have to be found
- spacedim		number of independent variables (space dimension)
- swarmcentre 	centre of the initial swarm or of the design space (array)
- swarmspan   	radius of the sphere in which the particles are generated
- nparticles  	number of particles
- niterations   number of generations
- bounds      	array with lower and upper bounds for each variable (spacedim x 2)
- results		array where results are saved (spacedim+1, that is coordinates of the optimum in the design space and corresponding output)
Outputs assigned to array passed by reference:
- zIncumb  position of the global optimal point (minimum)
- yIncumb  value of the objective function at the global minimum point found zIncumb
- NOT YET: function_output_Incumb value of the metamodel (NOT the objective function) at zIncumb
- NOT YET: constraints_output_Incumb  array containing the values of the output constraints at the found optimum
*/
//int psominimize(Binary_Node *ntree, Population *Pp, int spacedim, double* swarmcentre, double* swarmspan, int nparticles, int niterations, double** bounds, double* results)
int Population::psominimize(Binary_Node *ntree, int spacedim, double* swarmcentre, double* swarmspan, int nparticles, int niterations, double** bounds, double* results)
{
	int COMMENT=0; //1=print comments; 0=silent

	if (COMMENT) cout << "\n\nPopulation::psominimize()\n" << endl;

	// check input
	if (COMMENT) {
		cout << "\npsominimize" << endl;
		cout << "spacedim = " << spacedim << endl;
		cout << "swarmcentre: ";
		for (int i=0; i<spacedim; i++) {
			cout << swarmcentre[i] << ", ";
		}
		cout << endl;
		cout << "swarmspan: ";
		for (int i=0; i<spacedim; i++) {
			cout << swarmspan[i] << ", ";
		}
		cout << endl;
		cout << "nparticles = " << nparticles << endl;
		cout << "niterations = " << niterations << endl;
		for (int i=0; i<spacedim; i++) {
			cout << "Bounds of variable " << i << ": " << "(" << bounds[i][0] << ", " << bounds[i][1] << ")" << endl;
		}
	}

	// random value generator
	int seed = -1;
	if (seed<0) {
		// if seed = -1 use random seed (time)
	    time_t *tp = NULL; // seed the random value generator
		srand((unsigned int) time (tp));    //to be used normally
		seed = time(tp); // attention! two runs that starts with the same seed are identical!
		if (COMMENT) cout << "\n\nRandom value generator : seed =-1 : seed randomly generated = " << seed << endl;
	} else {
		// if seed >0 use the value given in input file
	    srand(seed);
	    if (COMMENT) cout << "\n\n Random value generator : used seed = " << seed << endl;
	}


	// list of variables (all the pointers are dynamically allocated, so remember to free memory!)
	double MAX_VALUE = 1000000.0;
	double** pop;		// matrix of particles' positions - matrix (spacedim x nparticles)
	double** vel;		// matrix of particles' velocities - matrix (spacedim x nparticles)
	double* yValues;	// objective function value associated to a particle - array (1xnparticles)

	double** zLocal;	// archive storing the position of particles personal best  - matrix (spacedim x nparticles)
	double* yLocal;		// best objective values obtained by each particle (particle/personal bests) - array (1 x nparticles)
	double* zIncumb;	// archive storing position of the particle that corresponds to the global best value of objective function
	double yIncumb;		// global best objective value - scalar

	double* function_output_Local; // value of the metamodel at the local (particle) minimum of the objective function in the design space
	double function_output_Incumb; // value of the metamodel at the global minimum (best among particles) of the objective function in the design space

	double* constraints_output_Local; // array of scalars for now - values of constraints at particle best position - matrix (nparticles x n_constraints=1)
	double constraints_output_Incumb; // scalar for now - values of constraints at global best position - array (1 x n_constraints)


	double* thisZ=zeros(spacedim);	// coordinates of current particle, dynamically initialised to a vector of zeros
	double** p_particle=new double*[spacedim]; // array of addresses of particle coordinates (used to modify/correct pop[][])

	// PSO engine tuning parameters
	double wCurr = 0.75;   // inertia weight: regulates the trade off between global and local exploration
	// see dynamic version of wCurr in CSC2011 paper...
	double wLocal = 0.15;  // coefficient of trust in the single particle
	double wGlob = 0.15;   // coefficient of trust in the swarm or global trend

	double r1=0; //random float value
	double r2=0; //random float value
	int v=0;
	double g=0.0;

	int ok=0;

	// --------------------------------------------------------------------



	// ------------------ INITIALISATION --------------------------
	if (COMMENT) {
		cout << "\nStart particles' position INITIALISATION" << endl;
		cout << "\nnparticles=" << nparticles << endl;
	}

	// yValues is an array storing the values of the function whose minima are sought corresponding to the particles
	// create and initialise yValues
	yValues = zeros(nparticles);

	// particles' position INITIALISATION ------------
	// pop is a matrix of dimensions space_dimension x n_particles : each column contains the coordinates of a single particle
	// create pop
	pop = zeros(spacedim, nparticles);
	// randomly initialise pop
	for (int particle=0; particle<nparticles; particle++) {
		for (int i=0; i<spacedim; i++) {
			v = rand() % 200 + 1; // random integer between 1 and 200. ATENTION when you use a normalised space
			g = (double(v)-100.0)/100.0; // random real number between -0.99 and 1.00
			pop[i][particle] = swarmspan[i]*g + swarmcentre[i];
		}

		////// IMPORTANT! Check that the initial swarm is inside the design space and regenerate points in case
		// necessary as swarmCentre and swarmSpan are not determined taking into account the input design space bounds...
		for (int i=0; i<spacedim; i++) {
			p_particle[i] = &(pop[i][particle]);
		}
		ok = regenerate_on_boundary(p_particle, spacedim, bounds);
	}


	// particles' velocity INITIALISATION ------------
	// it's the first set of velocities randomly generated...
    // vel is a matrix of the same dimensions of pop : each column represents a single particle velocity
    // velocity - random in [-swarmSpan,swarmSpan]
	// create vel
	vel = zeros(spacedim, nparticles);
	// initialise vel (FASTER IF IT IS DONE INSIDE THE CYCLE USED TO INITIALISE POP!!!)
	for (int j=0; j<nparticles; j++) {
		for (int i=0; i<spacedim; i++) {
			int v = rand() % 200 + 1; // random integer between 1 and 200
			double g = (double(v)-100.0)/100.0; // random real number between -0.99 and 1.00
			vel[i][j] = swarmspan[i]*g;
		}
	}

	// initialise best score so far in the swarm - global or incumbent
    yIncumb = MAX_VALUE;       			// initial value of the objective function (univariate function)
    zIncumb=zeros(spacedim);   			// initial position of the global minimum in the design space
    function_output_Incumb = MAX_VALUE;  // initial value of the metamodel at the global minimum in the design space
    constraints_output_Incumb = MAX_VALUE;   // scalar for now - it will be an array 1 x n_output_constraints

    // initialise local best score - considering the position of the single particle so far only
    yLocal = zeros(nparticles);
    for (int j=0; j<nparticles; j++) yLocal[j]=MAX_VALUE;
    zLocal = zeros(spacedim, nparticles);
    function_output_Local = zeros(nparticles);
    for (int j=0; j<nparticles; j++) function_output_Local[j]=MAX_VALUE;
    constraints_output_Local = zeros(nparticles); // values of constraints - single scalar for now
    // -------------------------------------------------------------------

	// --------- INITIAL FITNESS EVALUATION -----------------------------------
    if (COMMENT) cout << "\n\nPopulation::psominimize - Iteration = 0" << endl;
    double out = 0.0;
    for (int particle=0; particle<nparticles; particle++) {
		// assign coordinates of current particle
    	for (int i=0; i<spacedim; i++) thisZ[i] = pop[i][particle];

        // evaluate objective function related to particle ii and assign value to yValues[ii]
    	if (COMMENT) cout << "\nPopulation::psominimize - Particle n. " << particle << endl;
    	out = pso_objfunction(thisZ, spacedim, ntree); //p_objfun(thisZ, spacedim, ntree);

        yValues[particle]=out;    // assign objective function value to particle container

        /* print particle position
        cout << "\n Particle=" << particle << " (";
        for (int i=0; i<spacedim; i++) cout << thisZ[i] << " ";
        cout << ")   yValues=" << out;
        */

        function_output_Local[particle] = 0.0;   // assign metamodel value corresponding to particle best to particle container (identical to objective function for now)
        constraints_output_Local[particle] = 0.0; 	// values of the output constraints (scalar and null for now)
        // save minimisation objective and position of each particle in their respective container
        yLocal[particle] = yValues[particle];  // if you want to save the history, turn yLocal into a matrix...
        for (int i=0; i<spacedim; i++) zLocal[i][particle] = thisZ[i];

        // update incumbents (or global best)
		if (yValues[particle]<yIncumb) {
			yIncumb = yValues[particle];
			for (int i=0; i<spacedim; i++) zIncumb[i] = thisZ[i];
			function_output_Incumb = function_output_Local[particle];
			constraints_output_Incumb = constraints_output_Local[particle];
		}

    }

    // average output in the initial swarm
	//lastPopAve = mean(yValues)



	// ------------------ SEARCH THROUGH GENERATIONS -----------------------
	// generations/iterations
    for (int gg=0; gg<niterations; gg++) {

    	if (COMMENT) cout << "\n\nPopulation::psominimize - Iteration : " << gg << endl;

        // particles
        for (int particle=0; particle<nparticles; particle++) {
        	if (COMMENT) cout << "\nPopulation::psominimize - Particle n. " << particle << endl;

        	// GENERATE
        	r1 = double(rand() % 100)/100.0; // random float between 0.00 and 0.99
        	r2 = double(rand() % 100)/100.0; // random float between 0.00 and 0.99
        	for (int i=0; i<spacedim; i++) {
        		// update velocity
        		vel[i][particle] = wCurr*vel[i][particle] + wLocal*r1*(zLocal[i][particle]-pop[i][particle]) + wGlob*r2*(zIncumb[i]-pop[i][particle]);
        		// update position
        		pop[i][particle] = pop[i][particle] + vel[i][particle];
        	}

			//// CHECK that the particle is inside the design space (bounds) and regenerate a feasible point in case
        	// collect position
        	for (int i=0; i<spacedim; i++) {
        		p_particle[i] = &(pop[i][particle]);
        	}
            v = rand() % 100 + 1;     // v2 in the range 1 to 100
            if (v >= 50) {
                // regenerate on boundary
        		ok = regenerate_on_boundary(p_particle, spacedim, bounds);
            } else {
            	ok=regenerate_from_swarm_best(p_particle, spacedim, nparticles, bounds, zLocal);
            }
        	////


            // EVALUATE
        	for (int i=0; i<spacedim; i++) thisZ[i] = pop[i][particle];
        	//out = p_objfun(thisZ, spacedim, ntree);
        	out = pso_objfunction(thisZ, spacedim, ntree);  //out is the result of fitness evaluation. so RMSE error
			yValues[particle]=out;



            // UPDATE
			// if performance improved, update single particle and global performance
			if (yValues[particle]<yLocal[particle]) {
				// update particle best objective value so far
				yLocal[particle] = yValues[particle];
				// update position
				for (int i=0; i<spacedim; i++) zLocal[i][particle]=pop[i][particle];
				// update metamodel value corresponding to particle best
				function_output_Local[particle] = 0.0;
				// update constraints value corresponding to particle best
			    constraints_output_Local[particle] = 0.0; // values of the output constraints;
				// update incumbent (or global best)
                if (yLocal[particle]<yIncumb) {
                	// update global best objective value so far
                	yIncumb = yLocal[particle];
                	// update position
                	for (int i=0; i<spacedim; i++) zIncumb[i]=zLocal[i][particle];
                	// update metamodel value corresponding to global best
                	function_output_Incumb = function_output_Local[particle];
					// update constraints value corresponding to particle best
					constraints_output_Incumb = constraints_output_Local[particle];
                }
			}

		// end particles loop ---------------------------------------
        }


        if (COMMENT) {
        	cout << "\nGlobal best so far: ( ";
        	for (int i=0; i<spacedim; i++) cout << zIncumb[i] << " ";
        	cout << ")   yIncumb=" << yIncumb;
        }

		//wCurr = wCurr*.8   # dynamic reduction of global weight

        // average values and termination criterion
        //lastPopAve = mean(yValues)

        // STOP CRITERIA
        /*
        // if contour lines are targeted
        best_abs_residual = fabs(function_output_Incumb - output_level)
        if best_abs_residual<= min_residual:
            // terminate the search as the found point error is below the desired threshold
            return zIncumb, yIncumb, function_output_Incumb, constraints_output_Incumb
		*/


    // end iterations (generations) loop ---------------------------------------
    }

    // update tree with best parameters found (is this the right position? Maybe it is better to update the tree inside Population)
    update_complete_tree(ntree, zIncumb, spacedim);


    if (COMMENT) {
    	cout << "\n\nRESULT: Global best coordinates: ( ";
    	for (int i=0; i<spacedim; i++) cout << zIncumb[i] << " ";
    	cout << ")   yIncumb=" << yIncumb;
    }

	// store global minimum coordinates in results[0:spacedim-1] and corresponding function output in results[spacedim]
    for (int i=0; i<spacedim; i++) results[i]=zIncumb[i];
    results[spacedim]=yIncumb;
    if (COMMENT) {
    	cout << "\n\nCHECK RESULTS initialisation: Global best coordinates: ( ";
    	for (int i=0; i<spacedim; i++) cout << results[i] << " ";
    	cout << ")   yIncumb=" << results[spacedim];
    }

    // free all dynamically allocated variables (pop, vel, yValues, ..)
    delete[] thisZ; //*
    delete[] yValues; //*
    delete[] zIncumb; //*
    delete[] yLocal; //*

    for (int i=0; i<spacedim; i++) {
    	delete[] zLocal[i];  //*
    	delete[] pop[i];	//*
		delete[] vel[i];	//*
    }
    delete[] zLocal; //*
    delete[] pop; //*
    delete[] vel; //*

    delete[] function_output_Local;
    delete[] constraints_output_Local;
    delete[] p_particle; //*


    // ADD LATER ON: function_output_Incumb, constraints_output_Incumb;

    if (COMMENT) cout << "\n\nPopulation::psominimize() - exit \n" << endl;
    return 1;

}





// interface function called by Population::tuning_individual
//void hypso_launcher(double (*p_objfun)(double*, int, Binary_Node*), Binary_Node *ntree, int spacedim, double* x)
//void hypso_launcher(Population* Pp, Binary_Node *ntree, int spacedim, double* x)
void Population::hypso_launcher(RunParameters* pr, Binary_Node *ntree, int spacedim, double* x)
{
	int COMMENT=0;

	if (COMMENT) cout << "\n\nPopulation::pso_launcher\n" << endl;

	int out = 0;


	// input parameters
	int nparticles = pr->pso_nparticles; //10;  //min 40
	int niterations = pr->pso_niterations; //10;  //min 20
	double* swarmcentre=zeros(spacedim);
	double* swarmspan=zeros(spacedim);
	double** bounds=zeros(spacedim,2); // rows corresponds to variable (1st row= 1st variable) - 1st column= lower bound - 2nd column= upper bound

	// input parameters for swarm bounds: upper and lower bounds for sin/cos terms need to be computed from input data (Nyquist issue)
	for (int i=0;i<spacedim;i++) {
		swarmcentre[i]=(11.0+1.0)/2.0; // {1.0}; //swarmcentre[2] = {2.0, 2.0};  // for nD design spaces
		swarmspan[i]= (11.0-1.0)/2.0; //problem->Ny_omega_max; //0.5*3.14152/(6.40/100.0); //currently defined only for 1D problems  // for nD design spaces
		bounds[i][0] = 1.0; //-(problem->Ny_omega_max); //-0.5*3.14152/(6.40/100.0); //currently defined only for 1D problems // i-th variable lower bound
		bounds[i][1] = 11.0; //problem->Ny_omega_max; //0.5*3.14152/(6.40/100.0); // 200.0;  // i-th variable upper bound
		//(problem->v_list[tree->index_var[i]])->omega_lim
	}

	// optimisation output container (just an array for now, easier to integrate into preexisting code)
	double* res = zeros(spacedim+1);  // minimum coordinates [spacedim], value of the function at minimum [1],

	// launch single PSO search for global minima (function for now)
	out = psominimize(ntree, spacedim, swarmcentre, swarmspan, nparticles, niterations, bounds, res);

	// checks and output
	if (COMMENT) {
		if (out==1) {
			cout << "\n\nOK, psominimize executed correctly.";
			cout << "\nCoordinates of global minimum: ( ";
			for (int i=0; i<spacedim; i++) cout << res[i] << " ";
			cout << ")   Function at global minimum=" << res[spacedim];
		} else {
			cout << "\nERROR in psominimize!";
		}
	}

	// copy the optimised values found to the array provided as input
	for (int i=0; i<spacedim; i++) x[i]=res[i];

	// free dynamically allocated memory
	for (int j=0;j<spacedim;j++) delete[] bounds[j];
	delete[] bounds;
	delete[] swarmspan;
	delete[] swarmcentre;
	delete[] res;

	if (COMMENT) cout << "\n\nPopulation::exit pso_launcher\n" << endl;

}


//*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/**/*/*/*/*/*/*/*/*/*/*/*/*/



