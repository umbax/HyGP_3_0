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

#ifndef CLASS_NODE_BASE_H_
#define CLASS_NODE_BASE_H_

// dependencies
#include "../modules/Val_type.h"

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
	   
	// resets the parent node
    //void set_parent(Node *n) {parent = n;};
	void set_parent(Node *);

	// accesses the parent node
    //Node *get_parent(void) {return parent;};
	Node *get_parent(void);

};


#endif /* CLASS_NODE_BASE_H_ */


