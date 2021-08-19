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
 * nodes.cpp
 *
 *  Created on: May 18, 2017
 *      Author: umba
 */


using namespace std;

// included dependencies
#include <iostream>  // basic i/o commands: cout, cin, scientific, fixed
#include <stdlib.h>   //exit
#include <cstdlib>   // NULL
#include <cstdio>   // sprintf, sscanf
#include <string>    //  string class (C++)
#include <string.h>    //  string class (C++): strlen, strdup
#include <iomanip>  // manipulators (text format): setw

#include "nodes_types.h"   	//just a header, no source file
#include "../modules/primitives.h"
#include "../tree_functions/tree_operations.h"


/*
// base NODE class
// function definitions (virtual functions are defined in the derived classes)
*/


// node class constructor.  right now just initializes some common variables
Node::Node(void)
{
    // set number of children
    n_children = 0;

    // make this number invalid
    n_valid = 0;

    // set the parent ref to null
    parent = NULL;

	// set the operation control to an uncorrect value (will be changed by unary and binary node constructor)
    mult_div = 0;
    multdiv_valid = 0;  //it means that operations have to be checked

	// set subtree's depth
	depth = 0;

	//set the subtree's depth invalid
	d_valid = 0;

	// default type
	type = -1;

	// hits
	hits=0;
}



// node class destructor.  this will reset the parent's pointer to this node
// to NULL, so nothing bad happens in the future.  it also invalidate's
// the parent's child count
Node::~Node(void)
{
    // reset the parent's pointer to this one.

    // if there is no parent, no problem
    if (!parent)
      return;

    // if the parent is unary, this is easy
    if (parent->type == NODE_UNARY)
      // reset it
      ((Unary_Node *)parent)->set_child(NULL);

    else {
        // in this case (binary) we need to figure out if this was the left or
        // right child
        if (((Binary_Node *)parent)->get_left() == this)
          ((Binary_Node *)parent)->set_left(NULL);
        else
          ((Binary_Node *)parent)->set_right(NULL);
    }

   // force a new count
    parent->invalidate_count();
	//invalidate calc_depth: as a result of the change in number of children, depth may have changed
	parent->invalidate_calc_depth();  //try!

}





// this is called when the number of children of a node has changed.  it
// sets the invalid flag of this node and its parent to force a recalculation
// when the time comes
void Node::invalidate_count(void)
{
    // set the flag
    n_valid = 0;

    // if there is a parent, call the function
    if (parent)
      parent->invalidate_count();
}


// called when any node in the subtree changes.
void Node::invalidate_mult_div(void)
{
    // set the flag
    multdiv_valid = 0;

    // if there is a parent, call the function
    if (parent)
      parent->invalidate_mult_div();
}

void Node::invalidate_calc_depth(void)
{
    // set the flag
    d_valid = 0;

    // if there is a parent, call the function
    if (parent)
      parent->invalidate_calc_depth();
}


// function that resets the parent node
void Node::set_parent(Node *n)
{
	parent = n;
}

// function that returns the parent node's address
Node *Node::get_parent(void)
{
	return parent;
}

// function that checks the functions of the ancestor nodes linked upstream to the given node.
// Returns 1 if the node corresponding to the given pointer is a root node: so all checks have been passed up to reach root node
int Node::check_nodal_functions_upstream(void)
{
	int COMMENT=0;
	int pos=-1;

	if (COMMENT) cout << "\n\nBinary_Node::check_nodal_functions_upstream : START" << endl;

	// if root node, return 1 as all check have been passed
	if (parent==NULL) return 1;

	// here p_node is not a root node, so check its parent which is necessarily a functional node and in case call recursively this function
	string parent_fun="";

	if (parent->type == NODE_UNARY) {
		// BEFORE argument
		if (*(((Unary_Node*)parent)->get_func()->pos)==0) parent_fun=((Unary_Node*)parent)->get_func()->pre_sign;
		// AFTER argument
		if (*(((Unary_Node*)parent)->get_func()->pos)==1) parent_fun=((Unary_Node*)parent)->get_func()->post_sign;
		// BEFORE and AFTER argument
		if (*(((Unary_Node*)parent)->get_func()->pos)==2) {
			string pre_s=((Unary_Node*)parent)->get_func()->pre_sign;
			string post_s=((Unary_Node*)parent)->get_func()->post_sign;
			parent_fun=pre_s+post_s;
		}

		if (COMMENT) cout << "\nFunction of parent node: " << parent_fun;

		// check that the parent function belongs to a specific subset=[shift, neg]
		if ( (!strcmp(parent_fun.c_str(),"^1")) || (!strcmp(parent_fun.c_str(),"(-1.0)*")) ) {
			if (COMMENT) cout << "\nOK go upstream!" << endl;
			return parent->check_nodal_functions_upstream();
		} else
			return 0;
	}

	if (parent->type == NODE_BINARY) {
		parent_fun=((Binary_Node*)parent)->get_func()->sign;
		if (COMMENT) cout << "\nFunction of parent node: " << parent_fun;
		// check that the parent function belongs to a specific subset=[add, sub]
		if ( (!strcmp(parent_fun.c_str()," + ")) || (!strcmp(parent_fun.c_str()," - ")) || (!strcmp(parent_fun.c_str()," * "))) {
			if (COMMENT) cout << "\nOK go upstream!" << endl;
			return parent->check_nodal_functions_upstream();
		} else
			return 0;

	}

	return 0; // the function should never reach this point...

	if (COMMENT) cout<<"\n\nBinary_Node::check_nodal_functions_upstream : END";
}


/*
// class BINARY_NODE
// Binary Node function definitions
*/


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
    diverging=0;
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
	tree_mean = 0.0;
	tree_variance = 0.0;
	tree_min=0.0;
	tree_max=0.0;
	r_k = NULL;		// autocorrelation array
	r_k_size =0;
	first_acf_root_tree = 0.0;
	tot_variation_tree = 0.0;
	index_puls = NULL; 		//address of the array containing the position of the pulsations in x[]
	index_var = NULL;	 	//address of the array containing the position of the variable (in v_list) the corresponding pulsation in index_puls refers to
	p_par.clear();		// vector containing the addresses to terminal const nodes in the whole tree
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
	T8 = 1.0e+10;
	T9 = 1.0e+10;
	T10 = 1.0e+10;
	T11 = 1.0e+10;
	twin = 0;
}





// destructor - deletes the right and left subtrees
Binary_Node::~Binary_Node(void)
{
	// delete the left
    delete left;

    // delete the right
    delete right;

    // delete autocorrelation array
    delete[] r_k;
}


// set the left child
void Binary_Node::set_left(Node *l) {left = l;}

// set the right child
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


// function that inserts numerical parameters in selected locations to get complete trees
// (rules defined in Table 5.1 page 132 and at page 142 Armani's thesis)
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
	cout << setw(27) << "O1 : Fitness (RMSE) = " << scientific << setw(12) << fitness;
	cout << setw(6) << "  T1 = " << scientific << setw(12) <<  T1;
	cout << " % = "<< fixed << setw(12)  << 100.0*T1/F << endl;
	//
	cout << setw(27) << "O2 : n_tuning_parameters = " << setw(12) << n_tuning_parameters;
	cout <<  setw(6) << "  T2 = " << scientific << setw(12) << T2;
	cout << " % = " << fixed << setw(12)  << 100.0*T2/F << endl;
	//
	cout << setw(27) << "O3 : n_corrections = " << scientific << setw(12) <<n_corrections;
	cout << setw(6) << "  T3 = " << scientific << setw(12) << T3;
	cout << " % = " << fixed << setw(12) << 100.0*T3/F << endl;
	//
	cout << setw(27) << "O4 : size = " << scientific << setw(12) << count();
	cout << setw(6) << "  T4 = " << scientific << setw(12) << T4;
	cout << " % = " << fixed << setw(12) << 100.0*T4/F << endl;
	//
	cout << setw(27) << "O5 : pen_ord0 = " << scientific << setw(12) << pen_ord0;
	cout << setw(6) << "  T5 = " << scientific << setw(12) << T5;
	cout << " % = " << fixed << setw(12) << 100.0*T5/F << endl;

	cout << setw(27) << "O6 : pen_ord1 = " << scientific << setw(12) << pen_ord1;
	cout << setw(6) << "  T6 = " << scientific << setw(12) << T6;
	cout << " % = " << fixed << setw(12) << 100.0*T6/F << endl;

	cout << setw(27) << "O7 : depth_first_op = " << scientific << setw(12) << depth_first_op;
	cout << setw(6) << "  T7 = " << scientific << setw(12) << T7;
	cout << " % = " << fixed << setw(12) << 100.0*T7/F << endl;

	cout << setw(27) << "O8 : tree mean (train. data) = " << scientific << setw(12) << tree_mean;
	cout << " tree variance (train. data) = " << scientific << setw(12) << tree_variance;
	cout << setw(6) << "  T8 = " << scientific << setw(12) << T8;
	cout << " % = " << fixed << setw(12) << 100.0*T8/F << endl;

	cout << setw(27) << "O9 : high level polynomials : diverging = " << scientific << setw(12) << diverging;
	cout << setw(6) << "  T9 = " << scientific << setw(12) << T9;
	cout << " % = " << fixed << setw(12) << 100.0*T9/F << endl;

	cout << setw(27) << "tree min = " << tree_min << "   tree max = " << tree_max << endl;

	cout << "Autocorrelation array r_k size: " << r_k_size << endl;
//	cout << "\nAutocorrelation r_k:" << endl;
//	for (int i=0; i<r_k_size; i++) cout << r_k[i] << endl;
	cout << setw(27) <<  "10 : Point at which Autocorrelation Function drops by half: " << first_acf_root_tree;
	cout << setw(6) << "  T10 = " << scientific << setw(12) << T10;
	cout << " % = " << fixed << setw(12) << 100.0*T10/F << endl;

	cout << setw(27) <<  "11 : Tree total variation : " << scientific << tot_variation_tree;
	cout << setw(6) << "  T11 = " << scientific << setw(12) << T11;
	cout << " % = " << fixed << setw(12) << 100.0*T11/F << endl;

}




/*
// derived node class UNARY_NODE
*/


// unary node constructor - like the binary one, but a different function type
Unary_Node::Unary_Node(Node *parent_node, Unary_Func *func)
{
    // just store the args
    parent = parent_node;
    f = func;

    // set the child to null initially
    child = NULL;

    // set the type
    type = NODE_UNARY;
}


// destructor deletes the child subtree
Unary_Node::~Unary_Node(void)
{
    // delete the child
    delete child;

}

// adds a child node
void Unary_Node::set_child(Node *c) {child = c;}

// accesses the child node
Node* Unary_Node::get_child(void) {return child;}

// accesses the function
Unary_Func* Unary_Node::get_func(void) {return f;}

//function to modify f (see Population::mutation)
void Unary_Node::change_f (Unary_Func *p_fun) {f = p_fun;}


// unary node counter
int Unary_Node::count(void)
{
    // if the current count is valid, return that value
    if (n_valid)
      return n_children;

    // otherwise recalculate
    n_children = 0;

    // add in the child size (if it exists yet)
    if (child)
      n_children += child->count();

    // count this node
    n_children++;

    // this is valid now
    n_valid = 1;

    return n_children;
}



// unary node finder function
Node *Unary_Node::find(int number)
{
    // if the number is 1, return this node
    if (number==1)
      return this;

    // decrement the number, to count this node
    number--;

    // check this value against the number of nodes in the child subtree.
    // if it is less or equal, the target node is below this node.
    if ((child) && (number <= child->count()))
      return (child->find(number));

    // if its not in the child subtree, we are not going to find it
    return NULL;
}



// unary node value function.  checks for a child, and evaluates f using
// its value
Val Unary_Node::value(int* p_n_corrections)
{
    // check for missing child
    if (!child) {
        cerr << "\n Unary_Node::value : a unary node is missing its child!\n" << endl;
        return (Val)(0.0);
    }

    return (f->eval(child->value(p_n_corrections),p_n_corrections));
}



// print function - prints the subtree's expression to a string.
char *Unary_Node::print(void)
{
    // if one the child is missing, there is an error
    if (!child) {
        char msg[80];
        sprintf(msg,"UNARY NODE MISSING CHILD");
        return (strdup(msg));
    }

    // get the child's print string
	char *child_str = child->print();

    // make a new string to hold the composite string
    // (child length + operator length + 2 for the paren's + 1 for the ending character of the string(\0) )
    int newlength = strlen(child_str) + strlen(f->pre_sign) + strlen(f->post_sign) + 4 + 1;
    char *newstr = new (nothrow) char[newlength];

	//check allocation
	if (newstr == 0) {
		// write to stderr...
		cerr << "\nUnary_Node::print allocation error : not enough memory to create newstr";
		exit (-1);
	}


    // print into it
	// check correct value of *(f->pos) (positioning of unary function sign)
	if ((*(f->pos)<0) || (*(f->pos)>2)) {
		cerr << "\nUnary_Node::print : *(f->pos) error : not 0,1 or 2. Exit";
		exit (-1);
	}
	// sign of the operator BEFORE the argument
	//if (*(f->pos)== 0) sprintf(newstr,"(%s(%s))",f->sign,child_str);

	// sign of the operator AFTER the argument
	//if 	(*(f->pos) == 1) sprintf(newstr,"((%s)%s)", child_str, f->sign);

	// sign of the operator BEFORE AND AFTER the argument
	sprintf(newstr,"(%s(%s)%s)", f->pre_sign, child_str, f->post_sign);

    // delete the child's string
    delete [] child_str;

    // return the new string
    return (newstr);
}


// unary node counter
int Unary_Node::calc_depth(void)
{
    // if the current count is valid, return that value
    if (d_valid)
      return depth;

    // otherwise recalculate
    depth = 0;

    // add in the child size (if it exists yet)
    if (child)
      depth = child->calc_depth();

    // count this node
   depth++;

    // this is valid now
   d_valid = 1;

    return depth;   //depth of the subtree (the step to the current node -subtree root node- is included)
}

// function op_check
int Unary_Node::op_check(void) {return 0;}


// function that inserts numerical parameters in selected locations to get complete trees
// (rules defined in Table 5.1 page 132 and at page 142 Armani's thesis)
void Unary_Node::check_allocation(Node **p_tree, int number_op)
{
	//Val value;
	int child_op;
	Node* p_node = this;

	//----------------- TRY : insert a constant for SHIFT --------------------------
	// if number_op = 1  SHIFT
	if (number_op==1)  {
		// SHIFT: insert the parameter with the chosen value
		insert_parameter(this, p_tree,  &Add, (Val)1.0);      //1.0 instead of random value
		//
		//p_node = parent;
	}
	//--------------------------------------------------------------

	child_op = 0;

	// SCALE
	// if unary function is COS, SIN or EXP insert a parameter for the amplitude
	if ((&Sin) || (&Cos) || (&Exp) || (&RectWave))   //check if these primitives are used... to avoid run-time errors
		//if  ((!strcmp(f->sign,"sin")) || (!strcmp(f->sign,"cos")) || (!strcmp(f->sign,"exp")))  {	// WORKS
		if  ((!strcmp(f->pre_sign,Sin.pre_sign)) || (!strcmp(f->pre_sign,Cos.pre_sign)) || (!strcmp(f->pre_sign,Exp.pre_sign)) || (!strcmp(f->pre_sign,RectWave.pre_sign)))  {
			// insert AMPLITUDE PARAMETER
			insert_parameter(p_node, p_tree, &Mult, (Val)1.0);		   //1.0 instead of random value
		}

	// if unary function is POW1, SQUARE or CUBE insert a parameter for the translation of the argument
	if ((&Shift)) { // || (&Square) || (&Cube))   //check if these primitives are used... to avoid run-time errors
		if  (!strcmp(f->post_sign,"^1")) { // || (!strcmp(f->sign,"^2")) || (!strcmp(f->sign,"^3")))  {			//check on Unary_Func *f;
			// INSERT TRANSLATION PARAMETER
			// set the parameter for check_allocation
			child_op = 1;
		}
	}

	// unary function block HfreqSine: insert amplitude but avoid insertion for child
	if ((&HfreqSine)) {
		if  (!strcmp(f->pre_sign,"hfreqsin")) {  //&& (child->type==NODE_VAR add the condition that the child is a variable node...
			// insert AMPLITUDE PARAMETER
			insert_parameter(p_node, p_tree, &Mult, (Val)1.0);		   //1.0 instead of random value
			// set the parameter child_op to inhibit parameter insertion
			//cout << "\nUnary_Node::check_allocation : HFreqSine detected! No parameter to be inserted for child" << endl;
			child_op = -1;
		}
	}


	// for all the other cases, don't insert a parameter
	child->check_allocation(p_tree, child_op);
}



/*
// derived node class TERMINAL_VAR
*/


// function definitions

// variable terminal node constructor.  like the above, but stores a ref
Terminal_Var::Terminal_Var(Node *parent_node, Variable *v)
{
    // uh, store them
    parent = parent_node;
    var = v;

    // set the type
    type = NODE_VAR;
}


// function that returns the corresponding number of the variable in the node (see v_list)
int  Terminal_Var::get_var_number(void)
{
	string var_name;
	const char *var_no_char;
	int var_no_int;

	int COMMENT = 0;

	// retrieve the complete name
	var_name = var->name;
	if (COMMENT) cout << "\nTerminal_Var::get_var_index : var_name = " << var_name;

	// chop the variable letter off (first character)
	var_name.erase(0,1);
	var_no_char = var_name.c_str();
	if (COMMENT) cout << "\nTerminal_Var::get_var_index : var_no_char = " << var_no_char;

	// cast the remaining string into int
	sscanf(var_no_char, "%d", &var_no_int);
	if (COMMENT) cout << "\nTerminal_Var::get_var_index : var_no_int = " << var_no_int;

	return var_no_int;
}


// value function: returns the current variable value
//the first var refers to the member of Terminal_Var, the second to the member of Z (see line 247 master.cpp)
Val Terminal_Var::value(int* p_n_corrections)
{
	return var->value;
}

// this returns the value pointer (for copying)
Variable* Terminal_Var::val_p(void)
{
	return var;
}

// counter function is just 1 for this node
int Terminal_Var::count(void)
{
	return 1;
}

int Terminal_Var::calc_depth(void)
{
	return 0;
}

// constant terminal node find function
Node *Terminal_Var::find(int i)
{
    if (i==1)
      return this;
    else
      return NULL;
}





// variable terminal node print function
char *Terminal_Var::print(void)
{
    // make a copy of the variable name
    char *newstr = new (nothrow) char[strlen(var->name)+1];

	//check allocation
	if (newstr == 0) {
		// write to stderr...
		cerr << "\nTerminal_Var::print allocation error : not enough memory to create newstr";
		exit (-1);
	}

	// copy
	strcpy(newstr,var->name);

    // return the copy
    return newstr;
}


int Terminal_Var::op_check(void)
{
	return 1;
}


void Terminal_Var::check_allocation(Node **p_tree, int number_op)
{
	//Val value;
	Node* p_node = this;

	/*// ----------------------
	// To avoid inserting a parameter if parent is a sin, cos, exp
	char* op_parent = "ciao";
	 if (parent->type== 1) {
			op_parent = (((Unary_Node*)parent)->get_func())->sign;
		}
	if ((!strcmp(op_parent,"sin")) || (!strcmp(op_parent,"cos")) || (!strcmp(op_parent,"exp")))  {
		// cout << "\n op_parent = " << op_parent << " : parameter NOT introduced" << endl;
	}
	else {
	// -------------------------*/

	// 2 - SCALE insertion


	//choose random value
	//value = constant_generation();
	// SCALE: insert the parameter with the chosen value
	if (number_op!=-1) insert_parameter(p_node, p_tree, &Mult, (Val)1.0);     // number_op=-1 if parent is block HfreqSine

	p_node = parent;
	// if number_op = 1  SHIFT added
	if (number_op==1) {
		// 1 - SHIFT insertion
		//choose random value
		//value = constant_generation();
		// SHIFT: insert the parameter with the chosen value
		insert_parameter(p_node, p_tree,  &Add, (Val)1.0);     //1.0 instead of random value
	}

}




/*
//derived node class TERMINAL_CONST
*/

// function definitions

// constant terminal node constructor.  takes parent and const value
Terminal_Const::Terminal_Const(Node *parent_node, Val c)
{
    // just store them again
    parent = parent_node;
    constant = c;

    // set the type
    type = NODE_CONST;
}


// the value function just needs to return the constant
Val Terminal_Const::value(int* n_corrections)
{
	return constant;
}

// counter function is just 1 for this node
int Terminal_Const::count(void)
{
	return 1;
}

int Terminal_Const::calc_depth(void)
{
	return 0;
}


// constant terminal node find function
Node *Terminal_Const::find(int i)
{
    if (i==1)
      return this;
    else
      return NULL;
}





// constant terminal node print function
char *Terminal_Const::print(void)
{
	int COMMENT = 0;
	char *newstr;
	int dim=0;

	// just print the value to a string
    char prn[30];
 	if (COMMENT)  cout << "\n  prn =  " << prn << endl;
	dim=sprintf(prn,"%1.10e", (double)constant);
	if (COMMENT) {
		cout << "\n  prn =  " << prn << endl;
		cout << "\n  dim = " << dim;
		cout << "\nstrlen = " << strlen(prn);
	}

	// initialise newstr
    newstr = new (nothrow) char[strlen(prn)+1];

	//check allocation
	if (newstr == 0) {
		// write to stderr...
		cerr << "\nTerminal_Const::print allocation error : not enough memory to create newstr";
		exit (-1);
	}
	if (COMMENT) cout << "\n  newstr = " << newstr;


	// copy only the first strlen(prn) characters of prn to newstr
    strncpy(newstr,prn,strlen(prn));
	newstr[strlen(prn)]='\0';

	if (COMMENT) cout << "\n newstr = " << newstr;

	// return the copied string (newstr)
    return newstr;
}


int Terminal_Const::op_check(void)
{
	return 1;
}

void Terminal_Const::check_allocation(Node **, int)
{
	// intentionally left blank - do nothing!
}

// function to assign the numerical value if the constant
void Terminal_Const::assign(Val c)
{
	constant = c;
}




