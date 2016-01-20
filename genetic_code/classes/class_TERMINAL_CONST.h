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


#ifndef CLASS_TERMINAL_CONST_H_
#define CLASS_TERMINAL_CONST_H_

//derived node class TERMINAL_CONST

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


#endif /* CLASS_TERMINAL_CONST_H_ */
