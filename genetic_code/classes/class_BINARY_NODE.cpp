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


// class BINARY_NODE
// Binary Node function definitions

#include <iostream>  // basic i/o commands: cout, cin, scientific, fixed
#include <stdlib.h>   //exit
#include <cstdlib>   // NULL
#include <cstdio>   // sprintf, sscanf
#include <string.h>    //  string class (C++): strlen, strdup
#include <iomanip>  // manipulators (text format): setw

using namespace std;

// included dependencies
#include "./class_BINARY_NODE.h"
#include "../modules/primitives.h"
#include "../modules/nodes_types.h"   	//just a header, no source file
#include "../tree_functions/tree_operations.h"



// binary node constructor - takes the parent node and the function of this 
// node
Binary_Node::Binary_Node(Node *parent_node, Binary_Func *func)
{
    // just store the args
    parent = parent_node;
    f = func;

    // set the type
    type = NODE_BINARY;

    // set left and right to null initially
    left = right = NULL;

	//initialise state
    son_of = -1;
    parent_fitness = -1.0;
	fitness = .0;  // this is RMSE error, see fitness_func
	fitness_test = .0;  //
	F = .0;
	R2 = .0;
	R2_test = .0;
	hits = 0;
	hits_test = 0;
	n_tuning_parameters = 0;
	n_corrections = 0;
	n_corrections_test = 0;
	index_puls = NULL; 		//address of the array containing the position of the pulsations in x[]
	index_var = NULL;	 	//address of the array containing the position of the variable (in v_list) the corresponding pulsation in index_puls refers to
	n_pulsations = -1;  // see find_pulsations(...)
	depth_first_op = -1;   //depth of the first division found (if found). -1 means not found
	pen_ord0 = 0.0;
	pen_ord1 = 0.0;
	T1 = 1.0e+10;
	T2 = 1.0e+10;
	T3 = 1.0e+10;
	T4 = 1.0e+10;
	T5 = 1.0e+10;
	T6 = 1.0e+10;
	T7 = 1.0e+10;
	twin = 0;
}





// destructor - deletes the right and left subtrees
Binary_Node::~Binary_Node(void)
{
    // delete the left
    delete left;

    // delete the right
    delete right;
}


// adds a node to the left
void Binary_Node::set_left(Node *l) {left = l;}

// adds a node to the right
void Binary_Node::set_right(Node *r) {right = r;}

// accesses the left child
Node* Binary_Node::get_left(void) {return left;}

// accesses the right child
Node* Binary_Node::get_right(void) {return right;}

// accesses the function
Binary_Func* Binary_Node::get_func(void) {return f;}

//function to modify f (see Population::mutation)
void Binary_Node::change_f (Binary_Func *p_fun) {f = p_fun;}


// binary node counter
int Binary_Node::count(void)
{
    // if the current count is valid, return that value
    if (n_valid)
      return n_children;
    
    // otherwise recalculate

    // sum up the left and right (if they exist yet)
    n_children = 0;
    if (left)
      n_children += left->count();
    if (right)
      n_children += right->count();

    // count this node
    n_children++;

    // this is now valid
    n_valid = 1;

    return n_children; //total number of nodes, included the present
}





// binary node finder function
// input: number of the desired node 
// output: pointer to the desired node
Node *Binary_Node::find(int number)
{
    // if the number is 1, return this node
    if (number==1)
      return this;

    // decrement the number, to count this node
    number--;

    // check this value against the number of children in the left subtree.
    // if it is less or equal, the target node is in that subtree, so we 
    // should look there
    if ((left) && (number <= left->count()))
      return (left->find(number));

    // if its not in the left, subtract the number of children of the left
    // from the count
    if (left)
      number -= left->count();

    // check this against the number of children of the right subtree.
    // if it is less or equal, we will find it there
    if ((right) && (number <= right->count()))
      return (right->find(number));

    // if its not in the left or right, we are not going to find it
    return NULL;
}




    
// value function - pretty simple, just evaluate the function using the 
// left and right children's values
Val Binary_Node::value(int* p_n_corrections)
{
    // if the left or right does not exist, generate an error and return
    // zero
    if ((!left) || (!right)) {
        cerr << "\nBinary_Node::value : ERROR! Binary node is missing a child!\n" << endl;
        exit(-1);
    }

    // compute recursively the value of the tree. "This" is useful to update the n of corrections done by the protected operations 
	// in the future implement a check that returns immediately the value if the evaluation has been already done correctly... like in count(void)
    return (f->eval(left->value(p_n_corrections), right->value(p_n_corrections), p_n_corrections));
}





// print function - prints the subtree's expression to a string.
char *Binary_Node::print(void)
{
    // if one of the children is missing, there is an error
    if ((!left) || (!right)) {
        char msg[80];
        sprintf(msg,"BINARY NODE MISSING CHILD");
        return (strdup(msg));
    }

    // get the left and right children's print strings
	char *left_str = left->print();
   	char *right_str = right->print();

	// make a new string to hold the composite string
    // (left length + operator length + right length + 2 for the paren's)
	int newlength = strlen(left_str) + strlen(f->sign) + strlen(right_str) + 2 + 1;
    char *newstr = new (nothrow) char[newlength];

	if (newstr == 0) {
		// write to stderr...
		cerr << "\nBinary_Node::print allocation error : not enough memory to create newstr";
		exit (-1);
	}

    // print into it
    sprintf(newstr,"(%s%s%s)", left_str, f->sign, right_str);

	//delete the 2 child string
	delete [] left_str;
	delete [] right_str;
	
    // return the new string
    return (newstr);
}


int Binary_Node::op_check(void)
{
	if (multdiv_valid)
		//the check on the operation is reliable
		return mult_div;
	
	// otherwise the operations in the subtree must be checked
	int l_check=0, r_check=0, c_check=0;
	mult_div = 0;

	if (left)
		if (left->op_check()==1)
			l_check = 1;
		else
			l_check = 0;   //I know I could not write it...
	if (right)
		if (right->op_check()==1) 
			r_check = 1;
		else
			r_check = 0;

	// check if the operation is either a multiplication or a division
	if ( (!strcmp(f->sign,Mult.sign)) || (!strcmp(f->sign,SDiv.sign)) )
		c_check = 1;
	else
		c_check = 0;
	
	if ((l_check) && (r_check) && (c_check))
		mult_div = 1;
	else
		mult_div = 0;


	// this is now valid
	multdiv_valid = 1;
	return mult_div;
}


void Binary_Node::check_allocation(Node **p_tree, int number_op) 
{
	//Val value;
	Node* p_node = this;

	//----------------- TRY : insert a constant for SHIFT --------------------------
	// if number_op = 1  SHIFT
	if (number_op==1)  {
		// added: keep it, but SHIFT does not give good results!!!
		// 1 - SHIFT insertion
		//choose random value 
		//value = constant_generation();	
		// SHIFT: insert the parameter with the chosen value
		insert_parameter(this, p_tree,  &Add, (Val)1.0);		   //1.0 instead of random value	
		//
		p_node = parent;
	}
	//--------------------------------------------------------------

	// inserts a parameter if *,/ are the only operations in the subtree (included the current binary node)
	if (op_check())      {   //1 if *,/ are the only operations in the subtree  
		// SCALE
		//choose random value
		//value = constant_generation();
		// insert the parameter with the chosen value
		insert_parameter(p_node, p_tree, &Mult, (Val)1.0);        //1.0 instead of random value	
		// try (SHIFT)
		// go to subtrees and perform only SHIFT insertion
		//left->check_allocation(p_tree,1);
		//right->check_allocation(p_tree,1);
		//
	}
	else
	{
		// go to subtrees and perform SCALE insertion
		left->check_allocation(p_tree,3);
		right->check_allocation(p_tree,3);

	}
}


// function to find the addresses of the terminal_const nodes (parameters)
// input: nothing
// output: no of parameters (p_par is updated with the addresses of the parameters)
int Binary_Node::find_parameters(void)
{
	int n_nodes = count();
	int n_const = 0;
	
	// clear vector content to avoid appending new terms to old terms
	p_par.clear();
	
	for (int i=1; i<n_nodes+1; i++) {		//searching through the nodes, following the numeric order
		if ((this->find(i))->type == NODE_CONST) {		//true if the current node is terminal_const
			p_par.push_back (this->find(i));						// appends the pointer to terminal_const to the vector
			//printf ("Binary_Node: find pointers %f", (this->find(i))->value());		
			n_const++;
		}
	}	
	
	
	return n_const;    //returns the number of constant in p_par
}





// depth counter: binary node as subtree's root 
int Binary_Node::calc_depth(void)
{
    // if the current count is valid, return that value
    if (d_valid)
      return depth;
    
    // otherwise recalculate

    // sum up the left and right (if they exist yet)
   depth = 0;
    if ((left) && (right))
      depth = max(left->calc_depth(), right->calc_depth());

    // count this node
    depth++;

    // this is now valid
    d_valid = 1;

    return depth; //overall subtree depth (step to the current node included)
}




// function to show tree attributes
void Binary_Node::show_state(void)
{
	char *expr;
	expr = print();
	cout << "\n\nTree expression:";
	cout << "\n" << expr;
	cout << "\nDepth = " << calc_depth();
	cout << "\nF = " << F;
	cout << "\nR2 = " << R2;
	cout << "\nObjectives : " << endl;
	//
	cout << setw(27) << "O1 : Fitness (error) = " << scientific << setw(12) << fitness;
	cout << setw(6) << "  T1 = " << scientific << setw(12) <<  T1;
	cout << "% = "<< fixed << setw(12)  << 100.0*T1/F << endl;
	//
	cout << setw(27) << "O2 : n_tuning_parameters = " << setw(12) << n_tuning_parameters;
	cout <<  setw(6) << "  T2 = " << scientific << setw(12) << T2;
	cout << "% = " << fixed << setw(12)  << 100.0*T2/F << endl;
	//
	cout << setw(27) << "O3 : n_corrections = " << scientific << setw(12) <<n_corrections;
	cout << setw(6) << "  T3 = " << scientific << setw(12) << T3;
	cout << "% = " << fixed << setw(12) << 100.0*T3/F << endl;
	//
	cout << setw(27) << "O4 : size = " << scientific << setw(12) << count();
	cout << setw(6) << "  T4 = " << scientific << setw(12) << T4;
	cout << "% = " << fixed << setw(12) << 100.0*T4/F << endl;
	//
	cout << setw(27) << "O5 : pen_ord0 = " << scientific << setw(12) << pen_ord0;
	cout << setw(6) << "  T5 = " << scientific << setw(12) << T5;
	cout << "% = " << fixed << setw(12) << 100.0*T5/F << endl;
	//
	cout << setw(27) << "O6 : pen_ord1 = " << scientific << setw(12) << pen_ord1;
	cout << setw(6) << "  T6 = " << scientific << setw(12) << T6;
	cout << "% = " << fixed << setw(12) << 100.0*T6/F << endl;
	//
	cout << setw(27) << "O7 : depth_first_op = " << scientific << setw(12) << depth_first_op;
	cout << setw(6) << "  T7 = " << scientific << setw(12) << T7;
	cout << "% = " << fixed << setw(12) << 100.0*T7/F << endl;

}
