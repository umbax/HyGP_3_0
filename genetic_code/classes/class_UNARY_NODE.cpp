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


// derived node class UNARY_NODE

#include <iostream>  // basic i/o commands: cout, cin, scientific, fixed
#include <stdlib.h>   //exit
#include <cstdlib>   // NULL

using namespace std;

// included dependencies
#include "./class_UNARY_NODE.h"
#include "../modules/nodes_types.h"   	//just a header, no source file
#include <cstdio>   // sprintf
#include <string.h>    //  string class (C++): strlen, strdup
#include "../tree_functions/tree_operations.h"


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
    int newlength = strlen(child_str) + strlen(f->sign) + 4 + 1;
    char *newstr = new (nothrow) char[newlength];
	
	//check allocation
	if (newstr == 0) {
		// write to stderr...
		cerr << "\nUnary_Node::print allocation error : not enough memory to create newstr";
		exit (-1);
	}


    // print into it
	if (*(f->pos)== 0)
		// sign of the operator BEFORE the argument
    	sprintf(newstr,"(%s(%s))",f->sign,child_str);

	if 	(*(f->pos) == 1)
		// sign of the operator AFTER the argument
		sprintf(newstr,"((%s)%s)", child_str, f->sign);

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


void Unary_Node::check_allocation(Node **p_tree, int number_op) 
{
	//Val value;
	int child_op;
	Node* p_node = this;

	//----------------- TRY : insert a constant for SHIFT --------------------------
	// if number_op = 1  SHIFT
	if (number_op==1)  {
		// added: keep it, but SHIFT does not give good results!!!
		// 1 - SHIFT insertion
		//choose random value 
		//value = constant_generation();
		// SHIFT: insert the parameter with the chosen value
		insert_parameter(this, p_tree,  &Add, (Val)1.0);      //1.0 instead of random value	
		//
		//p_node = parent;
	}
	//--------------------------------------------------------------
	
	child_op = 0;

	// SCALE 
	// if unary function is COS, SIN or EXP insert a parameter for the amplitude
	if ((&Sin) || (&Cos) || (&Exp))   //check if these primitives are used... to avoid run-time errors 
		//if  ((!strcmp(f->sign,"sin")) || (!strcmp(f->sign,"cos")) || (!strcmp(f->sign,"exp")))  {	// WORKS
		if  ((!strcmp(f->sign,Sin.sign)) || (!strcmp(f->sign,Cos.sign)) || (!strcmp(f->sign,Exp.sign)))  {
		// AMPLITUDE PARAMETER
			//choose random value
			//value = constant_generation();
			// insert the parameter with the chosen value
			insert_parameter(p_node, p_tree, &Mult, (Val)1.0);		   //1.0 instead of random value		
		}

		// if unary function is POW1, SQUARE or CUBE insert a parameter for the translation of the argument
	if ((&Shift)) // || (&Square) || (&Cube))   //check if these primitives are used... to avoid run-time errors 
		if  (!strcmp(f->sign,"^1")) { // || (!strcmp(f->sign,"^2")) || (!strcmp(f->sign,"^3")))  {			//check on Unary_Func *f;
			// INSERT TRANSLATION PARAMETER
			// set the parameter for check_allocation 
			child_op = 1;
		}

	
	// for all the other cases, don't insert a parameter
	child->check_allocation(p_tree, child_op);
}

