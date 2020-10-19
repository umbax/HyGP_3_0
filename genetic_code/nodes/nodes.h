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

/*
 * nodes.h
 *
 *  Created on: May 18, 2017
 *      Author: umba
 */

#ifndef NODES_NODES_H_
#define NODES_NODES_H_


#include <vector>

using namespace std;

// dependencies
#include "../modules/Val_type.h"
//#include "./class_NODE_base.h"
#include "../modules/func_primitives_prototypes.h"
#include "../modules/primitives.h"
#include "../modules/variable.h"




// base NODE class definition (UTF-8 codification)

// base class for the 4 types of nodes
class Node {
  protected:
    Node *parent;               // the parent node
    int n_children;             // number of child nodes
    int n_valid;                // node count validity flag: 1 if the number of nodes is right
	int mult_div;				// operation variable: 1 if the node's operation is * or / (used in param. allocation)
	int multdiv_valid;		// operation control validity flag: 1 if the subtree's operations have been checked
	int depth;						//depth of the subtree rooted in the current node
	int d_valid;						//depth count validity flag: 1 if the subtree's depth is right

  public:
    Node(void);                 	// constructor
    virtual ~Node(void);                // virtual destructor (important to free memory allocate to derived classes)

    int type;                   // indicates the type of node
    int hits;					// number of hits - termination criterion

	void invalidate_count(void); // called when the number of children of this node has changed
    void invalidate_mult_div(void);    // called when the operations in the current node and in the subtrees
											// are not only * and /
	void invalidate_calc_depth(void); //called when the subtree's depth has changed (linked to the number of nodes)

	virtual Val value(int*) = 0; // this is the function that takes care of everything needed to assign a value to this node
    virtual int count(void) = 0; // returns the number of children of the subtree (the subtree root is the current node)
    virtual int calc_depth(void) = 0; //return the subtree depth (the subtree root is the current node)
	virtual Node *find(int) = 0; // used to return a particular node
    virtual char *print(void) = 0; // prints the subtree's expression to a string
    virtual int op_check(void) = 0;      //check the operations in the subtree
	virtual void check_allocation(Node **,int) = 0;   //performs the parameter allocation, if a terminal_var node

	// function that checks the functions of the ancestor nodes linked upstream to the given node.
	// Returns 1 if the node corresponding to the given pointer is a root node: so all checks have been passed up to reach root node
	int check_nodal_functions_upstream(void);
	// resets the parent node
    //void set_parent(Node *n) {parent = n;};
	void set_parent(Node *);

	// accesses the parent node
    //Node *get_parent(void) {return parent;};
	Node *get_parent(void);

};


/*
 * BINARY NODE
 */


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
	double tree_mean;		// mean of the values returned by the complete tree on the building data set
	double tree_variance;	// variance of the values returned by the complete tree on the building data set
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
	double T8;

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


/*
 * UNARY NODE
 */


class Unary_Node : public Node {

  private:
    Node *child;                // the (single) child
    Unary_Func *f;              // the function of this node

  public:
    // constructor takes the parent and the function
    Unary_Node(Node *,Unary_Func *);
    ~Unary_Node(void);          // destructor will delete the child

    // adds a child node
    void set_child(Node *);
    // accesses the child node
    Node *get_child(void);
    // accesses the function
    Unary_Func *get_func(void);
    //function to modify f (see Population::mutation)
    void change_f (Unary_Func*);


    // value function just passes the value of the child to f
    Val value(int*);            // value function
    int count(void);            // counter function
    int calc_depth(void); 		   //return the subtree depth (the subtree root is the current node)
    Node *find(int);            // finder function
    char *print(void);          // printer function


	// function to allocate parameters
	int op_check(void);
	void check_allocation(Node **, int);

};



/*
// TERMINAL_VAR
*/

// variable terminal node.  uses a Val ref for the value
class Terminal_Var : public Node {

  private:
    Variable *var;      // address of the variable value

  public:
    // constructor takes parent and variable pointer
    Terminal_Var(Node *, Variable *);

	int get_var_number(void);
	Variable *get_var_address(void){ return var;};    // NEW
	// value function: returns the current variable value
	//the first var refers to the member of Terminal_Var, the second to the member of Z (see line 247 master.cpp)
    Val value(int*);

    // this returns the value pointer (for copying)
    Variable *val_p(void);

    // counter function is just 1 for this node
    int count(void);

    int calc_depth(void);
   	Node *find(int);            // finder function
    char *print(void);          // printer function

	// function to allocate parameters
	int op_check(void);
	void check_allocation(Node **, int);


};



/*
//derived node class TERMINAL_CONST
*/


// a constant value terminal node
class Terminal_Const : public Node {

  private:
    Val constant;               // the constant value

  public:
    // constructor takes the parent and value
    Terminal_Const(Node *, Val);

    // the value function just needs to return the constant
    Val value(int*);
	// counter function is just 1 for this node
    int count(void);
    int calc_depth(void);
	Node *find(int);            // finder function
    char *print(void);          // printer function

	// function to allocate parameters
	int op_check(void);
	void check_allocation(Node **, int);

	// function to assign the numerical value if the constant
	void assign(Val c);
};


#endif /* NODES_NODES_H_ */
