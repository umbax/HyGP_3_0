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

#include <cstdlib>   // NULL

using namespace std;


// included dependencies
#include "./class_NODE_base.h"
#include "../modules/nodes_types.h"   	//just a header, no source file
#include "./class_BINARY_NODE.h"	// risk of recursive loop
#include "./class_UNARY_NODE.h"		// risk of recursive loop


// base NODE class
// function definitions (virtual functions are defined in the derived classes)

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
	multdiv_valid = 0;  //it means that operations have to be checked
	
	// set subtree's depth
	depth = 0;

	//set the subtree's depth invalid
	d_valid = 0;

	// default type
	type = -1;
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
