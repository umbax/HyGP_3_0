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


#ifndef CLASS_BINARY_NODE_H_
#define CLASS_BINARY_NODE_H_

#include <vector>

using namespace std;

// dependencies
#include "../modules/Val_type.h"
#include "./class_NODE_base.h"
#include "../modules/func_primitives_prototypes.h"

// derived node class BINARY_NODE

// binary operation node
class Binary_Node : public Node {

  private:
    Node *left;                 // left child
    Node *right;                // right child
    Binary_Func *f;             // the function of this node

  public:
    // constructor takes the parent and the function
    Binary_Node(Node *,Binary_Func *);
    ~Binary_Node(void);         // destructor will delete the children

    // adds a node to the left
    void set_left(Node *);
    // adds a node to the right
    void set_right(Node *);
    // accesses the left child
    Node *get_left(void);
    // accesses the right child
    Node *get_right(void);
    // accesses the function
    Binary_Func *get_func(void);
    //function to modify f (see Population::mutation)
    void change_f (Binary_Func*);

	// variables that define the STATE OF THE INDIVIDUAL (tree) - must be used only for root nodes
    int son_of;				// genetic operation that generated the tree: 0 = reproduction, 1 = crossover, 2 = subtree mutation, 3 = point mutation
    double parent_fitness;	// fitness of best parent (for crossover parents are two...)
    double fitness;         // fitness of the tree, now called "error" (if this is the root). Building data set. IT SHOULD BE OF TYPE Val  !!
    double fitness_test;	// fitness of the tree, now called "error" (if this is the root). Test data set
    Val F;					// value of the aggregative function used for multiobjective approach (see Population::evaluate(void), point 4)
	Val R2; 				// coefficient of determination, useful to compare qualities (see Understanding Statistics pag. 553 in pencil...)
	Val R2_test;
	int hits;				// number of hits of the tree (if this is the root). Building data set
	int hits_test;			// number of hits of the tree (if this is the root). Test data set
	int n_tuning_parameters;   	//number of tuning parameters (holds only for complete trees - with parameters...)
	int n_corrections;						// corrections are due to protected operations, like safe division
	int n_corrections_test;
	int *index_puls; 		//address of the array containing the position of the pulsations in x[]
	int *index_var;	 	//address of the array containing the position of the variable (in v_list) the corresponding pulsation in index_puls refers to
	int n_pulsations;  // see find_pulsations(...)
	int depth_first_op;  //former depth_first_div
	double pen_ord0;		//score related to the no. of times 0-order ineq. constraints are not satisfied
	double pen_ord1;		//score related to the no. of times first order ineq. constraints are not satisfied
	double T1;
	double T2;
	double T3;
	double T4;
	double T5;
	double T6;
	double T7;      // factorisation term

	int twin;										// if twin = 1, the current tree structure is identical to another tree structure (perhaps the best tree in the previous gen)

    // value function just needs to pass the values of the children to f
    Val value(int*);            			// value function: returns the value of the tree (assign the values of the variables first!!!)
    int count(void);            			// count function: returns the number of nodes the tree is composed of
    int calc_depth(void); 		   	//returns the tree (or subtree) depth (the subtree root is the current node)
	Node *find(int);            			// finder function: returns the pointer to the desired node, given its number
    char *print(void);          		// printer function

	// function to allocate parameters
	int op_check(void);
	void check_allocation(Node **, int);

	//parameter tuning (the following variables and functions are meant to be used with the root node trees[i])
	int find_parameters(void); //function that builds p_par
	vector <Node *> p_par;  //array containing the addresses to terminal const nodes in the whole tree

	// function to show tree attributes
	void show_state(void);
};


#endif /* CLASS_BINARY_NODE_H_ */
