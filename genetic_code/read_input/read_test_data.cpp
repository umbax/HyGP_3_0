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

// function to read the input file containing the test data set
void read_test_data(string FILE_TEST_DATA,  RunParameters* pr, ProblemDefinition* pb)
{

	int COMMENT=0;

	// declare the name
	string file;
	const char *expr1;
	// string to char conversion
	file =  FILE_TEST_DATA;
	expr1 = file.c_str();

	// open the stream
	ifstream fin;
	fin.open(expr1);

	// check that the file exists and opened in the right mode
	if (!fin.is_open()) {   //open attempt failed
		cerr << "\n\nread_test_data() : ERROR: The test data file " << FILE_TEST_DATA << " doesn't exist!!!" << endl;
		exit(-1);
	}
	else {
		cout << "\n\nread_test_data() : Reading from file " << FILE_TEST_DATA << "  : ..." << endl;
	}


	// in capital letters the parameters to import
	string Input;
	char linea[300];
	char ch;
	int check_header = 0;
	int m = 0;
	int count=0;
	int N_VAR = 0;
	int N_TEST_CASES = 0;

	//----------------------------------------------------------------------------------------------------
	// No of variables and no of test cases acquisition
	//----------------------------------------------------------------------------------------------------
	while (getline(fin, Input) && (!check_header)) {
		if (COMMENT) cout << "result of getline : " << Input << endl;
		// Input contains entire line
		if (Input[0]!='#') {
			stringstream Parse;
			Parse << Input;
			vector<double> Line;
			double Value;
			while (Parse >> Value) {
				// only the values that are double are actually copied... so comment is recognised!
				Line.push_back( Value );
			}
			N_VAR = Line[0];
			N_TEST_CASES = Line[1];
			if (COMMENT) cout << "\n N_VAR = " << N_VAR << "   N_TEST_CASES = " << N_TEST_CASES << endl;
			check_header=1;
		}
	}

	// checks
	if ((!N_VAR) || (!N_TEST_CASES)) {
		cerr << "\nERROR : N_VAR = " << N_VAR << "   N_TEST_CASES = " << N_TEST_CASES << ". Exit.";
		exit(-1);
	}

	if (pr->nvar!=N_VAR) {
		cerr << "\nERROR : the number of variables reported in " << FILE_TEST_DATA << " does not match with the number in the input file";
		exit(-1);
	}

	//----------------------------------------------------------------------------------------------------
	// Dynamic allocation of DATA  (NFITCASES rows, NVAR+1 columns): first allocate the row, then the columns for each row (no other way)
	//-----------------------------------------------------------------------------------------------------
	Val **TEST_DATA = NULL;
	TEST_DATA = new Val*[N_TEST_CASES];   //rows
	if (TEST_DATA==NULL)  {
		cerr << "\nERROR: dynamic allocation of TEST_DATA failed!! Test data can't be imported" << endl;
		exit(-1);
	}
	for (int i=0; i<N_TEST_CASES; i++) {
		TEST_DATA[i] = new Val[N_VAR+1];  // columns
		if (TEST_DATA[i]==NULL)  {
			cerr << "\nERROR: dynamic allocation of TEST_DATA[" << i << "]  failed!! Test data can't be imported" << endl;
			exit(-1);
		}
	}

	// get TEST_DATA line by line
	double Value;
	int row = 0;
	int col = 0;
	while (getline(fin, Input) && (row<N_TEST_CASES)) {
		if (COMMENT) cout << "result of getline : " << Input << endl;
		// Input contains entire line
	    if (Input[0]!='#') {
			stringstream Parse;
	    	Parse << Input;
	    	vector<double> Line;
	    	while (Parse >> Value) {
	    		// only the values that are double are actually copied... so comment is recognised!
	    		Line.push_back( Value );
	    		col++;
	    		if (COMMENT) cout << "filling Line... col = " << col << endl;
	    	}
	    	if (COMMENT) cout << "size of Line = "<< Line.size() << endl;
	    	if (Line.size() != N_VAR+1) {
	    		cerr << "Line " << row << ": expected " << N_VAR+1 << " values, received " << Line.size() << " values, aborting." << endl;
	    		cerr << "CHECK THE NUMBER OF VARIABLES DECLARED IN THE TEST DATA FILE!" << endl;
	    		exit(-1);
	    	}
		    if (COMMENT) cout << "TEST_DATA:" << endl;
		    for (int i=0; i<N_VAR+1; i++) {
		    	TEST_DATA[row][i] = Line[i];
		    	if (COMMENT) cout << scientific << TEST_DATA[row][i] << "  ";
		    }
		    row++;
		    col = 0;
		    if (COMMENT) cout << endl;
		}
		else
		   	if (COMMENT) cout << "comment" << endl;
	}

	// checks (m is the number of record imported)
	if (row < N_TEST_CASES) {
		cerr << "\nERROR: the number of valid records imported (" << m << ") is smaller than N_TEST_CASES (" << N_TEST_CASES << ").";
		cerr << "\nCheck the test data file (see commented lines...)!!" << endl;
		exit(-1);
	}

	// show the results
	if (COMMENT) {
		cout << "\nFile content." << endl;
		cout << "Value of the parameters:" << endl;
		cout << "\nN_VAR = " << N_VAR;
		cout << "\nN_TEST_CASES = " << N_TEST_CASES;

		cout << "\n\nContent of TEST_DATA(" << N_TEST_CASES  << "," << N_VAR+1 << "):" << endl;
		cout << setw(3) << " ";
		for (int j=0; j< N_VAR+1;  j++) cout << left << setw(22) << j ;
			cout << endl;
			for (int i=0; i< N_TEST_CASES; i++)  {//m must be equal to N_TEST_CASES
				cout << left << setw(3) << i;
				for (int j=0; j< N_VAR+1;  j++) {
					cout <<  scientific <<  left << setw(22)  << TEST_DATA[i][j];
				}
				cout << endl;
			}
		}

	//-----------------------------------------------------------------------------------------------------
	// all that follows is operations on read data (operations on file finished...)
	//-----------------------------------------------------------------------------------

	// set ProblemDefinition pointers so that test data can be retrieved
	pb->data_test = TEST_DATA;
	pb->n_test = N_TEST_CASES;

	// compute test data set statistical properties
	// in the future write a class "data set" with the following as a member function
	Val sum_output_test = .0;
	Val y_ave_test = .0;
	for (int k=0; k < pb->n_test; k++) {
		sum_output_test = sum_output_test + pb->data_test[k][pr->nvar]*pb->data_test[k][pr->nvar];
		y_ave_test = y_ave_test + pb->data_test[k][pr->nvar];
	}
	pb->sum_output_test = sum_output_test;
	pb->y_ave_test = y_ave_test/(double)(pb->n_test);

	Val Sy_test = .0;
	for (int k=0; k < pb->n_test; k++)
		Sy_test = Sy_test + (pb->data_test[k][pr->nvar] - pb->y_ave_test)*(pb->data_test[k][pr->nvar] - pb->y_ave_test);
	pb->Sy_test = Sy_test;
	//

	cout << "Ok";
}

