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


//derived node class TERMINAL_CONST

/*
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
	dim=sprintf(prn,"%1.5e", (double)constant);
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


