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


// TERMINAL_VAR

#ifndef CLASS_TERMINAL_VAR_H_
#define CLASS_TERMINAL_VAR_H_

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


#endif /* CLASS_TERMINAL_VAR_H_ */
