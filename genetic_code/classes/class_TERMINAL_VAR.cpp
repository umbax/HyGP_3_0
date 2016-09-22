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


// derived node class TERMINAL_VAR

#include <iostream>  // basic i/o commands: cout, cin, scientific, fixed
#include <stdlib.h>   //exit
#include <cstdio>   // sprintf, sscanf
#include <string>    //  string class (C++)
#include <string.h>    //  string class (C++): strlen, strdup

using namespace std;

// included dependencies
#include "./class_TERMINAL_VAR.h"
#include "../modules/primitives.h"
#include "../modules/nodes_types.h"   	//just a header, no source file
#include "../tree_functions/tree_operations.h"



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
		insert_parameter(p_node, p_tree, &Mult, (Val)1.0);     //1.0 instead of random value	
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
