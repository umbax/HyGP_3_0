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


#include "../tree_functions/tree_operations.h"


// functions performing operations on entire trees (or subtrees)


// this recursively copies a tree.  it takes a pointer to the tree to be
// input: pointer to the tree to be copied (src) and a pointer to the parent that the copied tree will have (parent)
// output: pointer to the new, copied tree.
// It copies only the function, variable or value of each node. That means that
// all the data members have the values imposed by the constructor!
Node *tree_copy(Node *src, Node *parent)
{
	int COMMENT =0;

	if (COMMENT)
		cout << "tree_copy : Node* src = " << src << "  Node* parent =  " << parent << endl;

	// make a new node of the correct type, then go to the children if needed
    switch (src->type) {

      Node *dest;

      case NODE_BINARY: {
          // less casting
          Binary_Node *bn = (Binary_Node *)src;
          // make the new node
          dest = (Node *)new Binary_Node(parent,bn->get_func());

          // do the left child
          ((Binary_Node *)dest)->set_left(tree_copy(bn->get_left(),dest));

          // do the right child
          ((Binary_Node *)dest)->set_right(tree_copy(bn->get_right(),dest));

          return dest;
      }

      case NODE_UNARY: {
          // less casting
          Unary_Node *un = (Unary_Node *)src;
          // make the new node
          dest = (Node *)new Unary_Node(parent,un->get_func());

          // do the child
          ((Unary_Node *)dest)->set_child(tree_copy(un->get_child(),dest));

          return dest;
      }

      case NODE_CONST:
        // make the new node
        dest = (Node *)new Terminal_Const(parent,src->value(NULL));  //value(NULL) as here there's no interest in counting the number of corrections

        return dest;


      case NODE_VAR:
        // make the new node
        dest = (Node *)new Terminal_Var(parent,((Terminal_Var *)src)->val_p());

		if (COMMENT)
			cout << "\ntree_copy : Node* src = " << src << "  Node* parent = parent = " << parent << endl;

        return dest;
    }
}





// utility function for sorting the population according to F
// (element 0 with lowest F)
int tree_comp_F(const void *a1, const void *a2)
{
    //REMEMBER! Now the criterion for sorting is F, not pure fitness=RMSE!
	// see void Population::evaluate(void)

	// recast the args
    const Binary_Node *x = *((const Binary_Node **)a1);
    const Binary_Node *y = *((const Binary_Node **)a2);

    // do the comparison
    if (x->F > y->F)
      return 1;  //return -1;
    else
      if (x->F == y->F)
        return 0;
      else
        return -1;  //return 1;
}





// utility function for sorting the population according to RMSE error (wrongly called fitness)
// (element 0 with lowest fitness)
int tree_comp_fitness(const void *a1, const void *a2)
{
	//NOTE! Now the criterion for sorting is pure fitness=RMSE!

	// recast the args
    const Binary_Node *x = *((const Binary_Node **)a1);
    const Binary_Node *y = *((const Binary_Node **)a2);

    // do the comparison
    if (x->fitness > y->fitness)
      return 1;  //return -1;
    else
      if (x->fitness == y->fitness)
        return 0;
      else
        return -1;  //return 1;
}


// Functions for inserting a node_const (parameter)
//-----------------------------------------------------------
// given a pointer to node_var, attaches to it a subtree made of a binary_node (multiplication) and its sons,
// a node_const and the node_var whose pointer is given
// input: pointer p to a node_var, pointer to the array "trees" containing the pointers to the root nodes, pointer to a Binary_Func, parameter's value 
// output:  void
void insert_parameter(Node *cur, Node **p_tree, Binary_Func *opr, Val value) 
{
	int COMMENT =0;

	if (COMMENT)  cout << "\ninsert_parameter(...)" << endl;
	// get the pointer to node_var's parent 
	Node *par = cur->get_parent();
	
	// set node_var's parent to NULL
	cur->set_parent(NULL);
	
	// create new binary_node
	Binary_Node *newbin = new Binary_Node(par, opr);  
	((Node *)newbin)->set_parent(NULL);

	// set binary_node just created as son of cur's parent
	if (par!=NULL)
	{
		if (par->type == NODE_UNARY)
			((Unary_Node *)par)->set_child((Node *)newbin);
		if (par->type == NODE_BINARY)
			{
			if (((Binary_Node *)par)->get_left() == cur)
				// node_var was a left child
				((Binary_Node *)par)->set_left((Node *)newbin);
			else
				// node_var was a right child
				((Binary_Node *)par)->set_right((Node *)newbin);
		}
	
	// set par as parent of the newly created binary_node
	((Node *)newbin)->set_parent(par);
	}
	else
		{
		//the current node is the root node, trees[n] must be updated
		//cout << "\nsubstituting root\n";
		*((Binary_Node **)p_tree)=newbin;
	};

	// create a new terminal_const node
	Terminal_Const *newconst =new Terminal_Const((Node *)newbin, value);
	// link it to the binary node just created, as left child
	newbin->set_left((Node *)newconst);

	//link the old Terminal_Var node
	//as right child of the binary node just created
	newbin->set_right((Node *)cur);
	
	
	// set the newly created binary_node as parent of Terminal_Var node
	cur->set_parent((Node *)newbin);
	
	//INVALIDATE OPERATIONS
	//invalidate count: the number of children has changed
	((Node *)newbin)->invalidate_count();
	//invalidate calc_depth: as a result of the change in number of children, depth may have changed
	((Node *)newbin)->invalidate_calc_depth();
	    
	if (COMMENT) cout << "\nexit insert_parameter(...)" << endl;
}



// function that inserts a subtree, replacing the given p_old_subtree with the selected p_new_subtree
int insert_subtree(Node **p_tree, Node *p_old_subtree, Node *p_new_subtree)             //&trees[i]
{	
	int COMMENT =0;    //1 comments, 0 silent...
	char *expr;
	if (COMMENT) {
		cout << "\nsubtree_mutation : replacing old subtree with new one" << endl;
		cout << "\n	address of tree root node = " << *p_tree;
		cout << "\n	address of address of tree root node = " << p_tree;
		cout << "\n	address of NEW subtree root node = " << p_new_subtree;
	}

	// get the pointer to old subtree's parent 
	Node *parent = p_old_subtree->get_parent();
	if (COMMENT) cout << "\n	address of selected node's parent = " << parent;

	// tell the father (parent) who is the child
	if (parent!=NULL)  { 
		
		//the current node is not the root node...
		if (parent->type == NODE_UNARY) {
			if (COMMENT) cout << "\n	address of old child = " << ((Unary_Node*)parent)->get_child();
			((Unary_Node *)parent)->set_child((Node *)p_new_subtree);
			if (COMMENT) cout << "\n	address of new child = " << ((Unary_Node*)parent)->get_child();
		}
		if (parent->type == NODE_BINARY) {
			if (((Binary_Node *)parent)->get_left() == p_old_subtree) {
				// p_old_subtree was a left child
				if (COMMENT) cout << "\n	address of old left child = " << ((Binary_Node*)parent)->get_left();
				((Binary_Node *)parent)->set_left((Node *)p_new_subtree);
				if (COMMENT) cout << "\n	address of new left child = " << ((Binary_Node*)parent)->get_left();
			}
			else {
				// p_old_subtree was a right child
				if (COMMENT) cout << "	address of old right child = " << ((Binary_Node*)parent)->get_right();
				((Binary_Node *)parent)->set_right((Node *)p_new_subtree);
				if (COMMENT) cout << "	address of new right child = " << ((Binary_Node*)parent)->get_right();
			}
		}
	
		// tell the new child (p_new_subtree) who is his father
		((Node *)p_new_subtree)->set_parent(parent);
		if (COMMENT) cout << "\n	address of new child's parent = " << ((Node *)p_new_subtree)->get_parent();
	
		//INVALIDATE OPERATIONS
		//invalidate count: the number of children has changed
		((Node *)parent)->invalidate_count();
		//invalidate calc_depth: as a result of the change in number of children, depth may have changed
		((Node *)parent)->invalidate_calc_depth();	
	}

	else	{
		//THE CURRENT NODE IS THE ROOT NODE: trees[n] must be updated!!!
		// no invalidate count or calc_depth because the parent is NULL!
		//cout << "\nsubstituting root\n";
		*((Binary_Node **)p_tree) = (Binary_Node*)p_new_subtree;
	}


	
	// set p_old_subtree's parent to NULL: REALLY IMPORTANT OPERATION TO DO BEFORE DELETING A NODE!
	p_old_subtree->set_parent(NULL);
	if (COMMENT) cout << "\n	address of OLD subtree root node parent = " << p_old_subtree->get_parent() << "MUST BE NULL !!!!";
	
	//free memory occupied by p_old_subtree
	delete p_old_subtree;

	expr = (*p_tree)->print();
	if (COMMENT) {
		cout << "\n	address of tree root node after subtree insertion = " << *p_tree << endl;
		cout << "\nNew tree after subtree mutation: " << expr;
	}

	// if everything is ok
	return 1;
}



// function to replace a subtree with a terminal_const node with a given value.
// input: address of the subtree root node, value of the terminal_const node replacing the subtree
// output: void
void eliminate_subtree(Node *subroot, Val value)
{
	int COMMENT =0;  //1 comments, 0 silent...
	Node* parent;
	int type_parent;
	int left = 0;

	//cout << "\n\neliminate_subtree" << endl;

	// retrieve the address of the parent of the subtree root node and store its type
	parent = subroot->get_parent();
	type_parent = parent->type;
	
	if (COMMENT) cout << "Parent of the subtree root node - address " << parent << ", type " << type_parent << endl;
	
	// if parent is a binary node, find out if the subtree root node is a left or a right child
	if (type_parent == NODE_BINARY)
		if (((Binary_Node*)parent)->get_left() == subroot)
			// subroot is a left child
			left = 1;
		else
				// subroot is a right child
			left = 0;

	
	// delete the subtree, deleting the subtree root node (upward and downward connections are safely managed by node deconstructor).
	delete subroot;

	// create a terminal const node with value 0 (upward connection managed by node constructor)
	Terminal_Const *new_const_node = new Terminal_Const(parent, value);

	// tell the parent that the newly created node is its child (take care of the type of the parent node)
	if (type_parent == NODE_BINARY) {
		if (left)
			// left child
			((Binary_Node*)parent)->set_left((Node*)new_const_node);
		else
			// right child
				((Binary_Node*)parent)->set_right((Node*)new_const_node);
	}

	if (type_parent == NODE_UNARY)
		((Unary_Node*)parent)->set_child((Node*)new_const_node);

}


// forse è meglio riscrivere questa funzione: si osserva il nodo corrente e se
// il parent è divisione, il nodo corrente è figlio destro e uguale a 0 (terminal_const)
// o può assumere valore nullo del design space (terminal_var) allora si sostituisce
// Strategy: Ed -> comment out lines followed by //editstruct2
// Strategy: Ed2 -> remove comment from lines followed by //editstruct2
// Strategy: Ed3 -> as Ed2, but the check whether variable domain contains 0 is replaced by random choice
void edit_tree(Node *cur_node, Node** p_tree)
{
	int COMMENT =0;  //1 comments, 0 silent...

	// list of operations whose presence has to be checked (remember: if div is found, found_div is 0 - outcome of strcomp)
	int found_div=-1;
	int c=0;	// for random integer

	// get lower and upper bounds of the design space
	// see v_list[k].lower_b and upper_b

	// check that Sdiv and inv are actually among the primitives

	// check that
	if (cur_node->type == NODE_BINARY) {
		found_div =  strcmp((((Binary_Node*)cur_node)->get_func())->sign,SDiv.sign);
		Node *l_child = ((Binary_Node*)cur_node)->get_left();
		Node *r_child = ((Binary_Node*)cur_node)->get_right();
		//check if if current node is a division and right node is a variable
		if ((found_div==0) && (r_child->type == NODE_VAR)) {
			// check that the variable is zero in the design space
			Variable *var = ((Terminal_Var*)r_child)->get_var_address();
		
			c = int_rand(2); //Ed3
			if (c) {//Ed3: random decision  //Ed, Ed2: CHECK ON VARIABLE DOMAIN: if ((var->lower_b)*(var->upper_b) <= 0.0) {
				// here the variable assumes zero value in the design space:
				// either replace the node with a terminal const = 1.0 or add a terminal_const = 1.0
				c = int_rand(2);       //editstruct2
				if (c) {               //editstruct2
					// replace the variable with a terminal const = 1.0
					if (COMMENT) cout << "\nedit_tree: variable " << var->name << " found: lb*up = " << (var->lower_b)*(var->upper_b);
					if (COMMENT) cout << "\nreplacing it with 1.0";
					   eliminate_subtree(r_child, (Val)1.0);
				} else { ////////////////////////////////////////// NEW PART //editstruct2
					// add a terminal_const = 1.0     //editstruct2
					insert_parameter((Node*)r_child, (Node **)p_tree, &Add, (Val)1.0) ;  //editstruct2
				}     //editstruct2
        
			} else
					// right child well-defined: launch same function in left sibling
					edit_tree(l_child, p_tree);

		} else {
				// launch same function in right child and left child.
				edit_tree(l_child, p_tree);
				edit_tree(r_child, p_tree);
		}

	}
	if (cur_node->type == NODE_UNARY) {
		Node *child = ((Unary_Node*)cur_node)->get_child();
		// for now unary functions give all well-defined results: launch same function in left sibling
		// In the future insert procedure to avoid undefined inv here
		edit_tree(child, p_tree);
	}
	if ((cur_node->type == NODE_VAR) || (cur_node->type == NODE_CONST)) {
		return;
	}

}



