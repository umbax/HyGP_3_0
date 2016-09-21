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


// UNARY NODE


#ifndef CLASS_UNARY_NODE_H_
#define CLASS_UNARY_NODE_H_

// dependencies
#include "../modules/func_primitives_prototypes.h"

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

#endif /* CLASS_UNARY_NODE_H_ */
