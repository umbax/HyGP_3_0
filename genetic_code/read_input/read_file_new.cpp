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

#include "./read_file_new.h"

// function that
// - reads the main input file, containing the HyGP hyperparameters and the building data set
// (the building data set can be split in tuning data set and validation data set (crossvalidation approach),
// but this is done later with the function kfolds_split (see ProblemDefinition)
// - initialise a few member variables of ProblemDefinition (see Sy, y_ave, etc)
void read_input_file(string FILE_INPUT,  RunParameters* pr, ProblemDefinition* pb)
{

	int COMMENT=0;

	cout << "\n\nread_input_file : enter" << endl;

	// declare the name
	string file;
	const char *expr1;
	// string to char conversion
	file =  FILE_INPUT;
	expr1 = file.c_str();

	// open the stream
	ifstream fin;
	fin.open(expr1);

	// check that the file exists and opened in the right mode
	if (!fin.is_open()) {   //open attempt failed
		cerr << "\nERROR: The input file doesn't exist!!!" << endl;
		exit(-1);
	}
	else {
		cout << "\nReading from file " << FILE_INPUT << "  :" << endl;
	}


	// in capital letters the parameters to import
	char linea[300];
	char ch;
	int NUM_U_FUNCS, NUM_B_FUNCS;
	string func_bin;  //for binary functions' list - Andrey's idea
	string func_uni;  //for unary functions' list - Andrey's idea
	int check_header = 0;
	int m = 0;
	int count=0;
	string word;
	double p_par_values[33];
	int n_par = 33;
	string PAR_LABELS[33]	=		 {	"SEED=",
                                                            "NVAR=",
															"MINRAND=",
															"MAXRAND=",
															"MAX_N_PERIODS=",
															"STEP=",
															"NFITCASES=",
															"METHOD=",
															"DEPTH_MAX=",
															"DEPTH_MIN=",
															"DEPTH_LIM=",
															"p_FULL=",
															"REPR_RATE=",
															"CROSS_RATE=",
															"MUT_RATE=",
															"COMP_RATE=",
															"NEW_RATE=",
															"M=",
															"G=",
															"NORMALISED_ERR=",
															"MINMAX=",
															"W_COMPLEXITY=",
															"W_N_CORRECTIONS=",
															"W_SIZE=",
															"W_FACTORISATION=",
															"N_GUESSES=",
															"CROSSVALIDATION=",
															"FOLDS_N=",
															"THRESHOLD=",
															"N_INEQUALITY0=",
															"W_PEN_ORD0=",
															"N_INEQUALITY1=",
															"W_PEN_ORD1="};
	string BIN_LABEL = "BINARY_FUN=";
	string UNA_LABEL = "UNARY_FUN=";



	//----------------------------------------------------------------------------------------------------
	// Parameters acquisition
	//----------------------------------------------------------------------------------------------------
	// 1st step - READ PARAMETERS, discarding commented lines (#) or not recognised lines
	int max_checks = 5;
	int n_checks = 0;
	while ((!check_header) && (n_checks<max_checks)) {
		cout << "\nReading " << PAR_LABELS[count] << "...";
		fin >> word;
		if (!word.compare(PAR_LABELS[count])) {
			fin >> p_par_values[count];
			count++;
			cout << "Ok.";
			n_checks = 0;
		}
		fin.getline(linea, 300); // go to next line, the rest is comment
		n_checks++;
		if (count==n_par) check_header=1;
	}

	if (!check_header) {
		cout << "\n" << PAR_LABELS[count] << " not found. Exit.";
		cerr << "\nERROR : " << PAR_LABELS[count] << " not found. Exit";
		exit(-1);
	}
	else
		cout << "\nParameters imported." << endl;

	// 2nd step - READ BINARY FUNCTIONS, discarding commented lines (#) or not recognised lines
	check_header = 0;
	while (!check_header) {// && (fin.get(ch))) {     //gets FALSE when end of file is reached or header is read (check_header =1)
		fin >> word;
		if (!word.compare(BIN_LABEL)) {
			getline(fin,func_bin);
			check_header = 1;
		}
		else
			fin.getline(linea, 300); // go to next line, the rest is comment
	}

	if (COMMENT)  cout << "\nBinary functions written in input file:" << func_bin << endl;

	// 3rd step - READ UNARY FUNCTIONS, discarding commented lines (#) or not recognised lines
	check_header = 0;
	while (!check_header) {    //gets FALSE when end of file is reached or header is read (check_header =1)
		fin >> word;
		if (!word.compare(UNA_LABEL)) {
			getline(fin,func_uni);
			check_header = 1;
		}
		else
			fin.getline(linea, 300); // go to next line, the rest is comment
	}
	if (COMMENT)  cout << "\nUnary functions written in input file:" << func_uni << endl;


	// here surely the parameters and the primitives have been imported. Now allocate matrix DATA dynamically...
	pr->seed = (int)(p_par_values[0]);
	pr->nvar = (int)(p_par_values[1]);
	pr->minrand = p_par_values[2];   
	pr->maxrand = p_par_values[3];  
	pr->step = p_par_values[5];         
	pr->max_n_periods = p_par_values[4];
	pr->nfitcases = (int)(p_par_values[6]);
	pr->method = (int)(p_par_values[7]);
	pr->depth_max = (int)(p_par_values[8]);
	pr->depth_min = (int)(p_par_values[9]);
	pr->depth_lim = (int)(p_par_values[10]);
	pr->p_FULL = p_par_values[11];
	pr->repr_rate = p_par_values[12];
	pr->cross_rate = p_par_values[13];
	pr->mut_rate = p_par_values[14];
	pr->comp_rate = p_par_values[15];
	pr->new_rate = p_par_values[16];
	pr->M = (int)(p_par_values[17]);
	pr->G = (int)(p_par_values[18]);
	pr->normalised=(bool)(p_par_values[19]);
	pr->minmax=(bool)(p_par_values[20]);
	pr->w_complexity = p_par_values[21];
	pr->w_n_corrections = p_par_values[22];
	pr->w_size = p_par_values[23];
	pr->w_factorisation = p_par_values[24];  // switch between standard approach (<0) and factorisation bonus (>0)
	pr->n_guesses = (int)(p_par_values[25]);
	pr->crossvalidation = (int)(p_par_values[26]);
	pr->folds_n = (int)(p_par_values[27]);
	pr->threshold = p_par_values[28];
	pr->n_inequality0 = (int)p_par_values[29];
	pr->w_pen_ord0 = p_par_values[30];
	pr->n_inequality1 = (int)p_par_values[31];
	pr->w_pen_ord1 = p_par_values[32];
	

	//----------------------------------------------------------------------------------------------------
	// Dynamic allocation of DATA  (NFITCASES rows, NVAR+1 columns): first allocate the row, then the columns for each row (no other way)
	//-----------------------------------------------------------------------------------------------------
	Val **DATA = NULL;
	DATA = new Val*[pr->nfitcases];
	if (DATA==NULL)  {
		cerr << "\nERROR: dynamic allocation of DATA failed!! Input data can't be imported" << endl;
		exit(-1);
	}
	for (int i=0; i<pr->nfitcases; i++) {
		DATA[i] = new Val[pr->nvar+1];
		if (DATA[i]==NULL)  {
			cerr << "\nERROR: dynamic allocation of DATA[" << i << "]  failed!! Input data can't be imported" << endl;
			exit(-1);
		}
	}


	// get DATA line by line
	string Input;
	double Value;
	int row = 0;
	int col = 0;
	while (getline(fin, Input) && (row<pr->nfitcases)) {
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
	    		if (Line.size() != pr->nvar+1) {
	    			cerr << "Line " << row << ": expected " << (pr->nvar)+1 << " values, received " << Line.size() << " values, aborting." << endl;
	    			cerr << "CHECK THE NUMBER OF VARIABLES DECLARED IN THE INPUT FILE!" << endl;
	    			exit(-1);
	    		}
	    		if (COMMENT) cout << "DATA:" << endl;
	    		for (int i=0; i<(pr->nvar)+1; i++) {
	    			DATA[row][i] = Line[i];
	    			if (COMMENT) cout << scientific << DATA[row][i] << "  ";
	    		}
	    		row++;
	    		col = 0;
	    		if (COMMENT) cout << endl;
	    }
	    else
	    	if (COMMENT) cout << "comment" << endl;
	}



	// checks (m is the number of record imported)
	if (row < pr->nfitcases) {
		cerr << "\nERROR: the number of valid records imported (" << m << ") is smaller than NFITCASES (" << pr->nfitcases << ").";
		cerr << "\nCheck the input file (see commented lines...)!!" << endl;
		exit(-1);
	}

	// show the results
	if (COMMENT) {
		cout << "\nFile content." << endl;
		cout << "Value of the parameters:" << endl;
		for (int k=0; k<n_par; k++)
		 	cout << PAR_LABELS[k] << p_par_values[k] << endl;

		cout << "Content of DATA(" <<  pr->nfitcases << "," << (pr->nvar)+1 << "):" << endl;
		cout << setw(3) << " ";
		for (int j=0; j< (pr->nvar) +1;  j++) cout << left << setw(22) << j ;
		cout << endl;
		for (int i=0; i< pr->nfitcases; i++)  {//m must be equal to NFITCASES
				cout << left << setw(3) << i;
				for (int j=0; j< (pr->nvar) +1;  j++) {
					cout <<  scientific <<  left << setw(22)  << DATA[i][j];
				}
				cout << endl;
		}
	}

	//---------------------------------------------------------------------------------------------------
	// Dynamic allocation of matrix DATA_INEQUALITY0 (FITCASES rows, NVAR+1 columns): first allocate the row, then the columns for each row (no other way)
	//---------------------------------------------------------------------------------------------------
	Val **DATA_INEQUALITY0 = NULL;
	char *CONSTRAINTS0 = NULL;
	if (pr->n_inequality0) {
		DATA_INEQUALITY0 = new Val*[pr->n_inequality0];
		if (DATA_INEQUALITY0 == NULL)  {
			cerr << "\nERROR: dynamic allocation of DATA_INEQUALITY0 failed... Exit" << endl;
			exit(-1);
		}
		for (int i=0; i<pr->n_inequality0; i++) {
			DATA_INEQUALITY0[i] = new Val[pr->nvar+1];
			if (DATA_INEQUALITY0[i]==NULL)  {
				cerr << "\nERROR: dynamic allocation of DATA_INEQUALITY0["<< i << "] failed... Exit" << endl;
				exit(-1);
			}
		}

		CONSTRAINTS0 = new char[pr->n_inequality0];

		// get DATA_INEQUALITY0 and CONSTRAINT0 reading line by line
		char temp;
		row = 0;
		col = 0;
		while (getline(fin, Input) && (row<pr->n_inequality0)) {
			if (COMMENT) cout << "\nResult of getline : " << Input << endl;
			// Input contains entire line
			if (Input[0]!='#') {
					stringstream Parse;
					//cout << Parse.str() << endl;
					Parse << Input;
					vector<double> Line;
					temp = 'n';
					while (Parse >> Value) {
						// only the values that are double are actually copied... so comment is recognised!
						Line.push_back( Value );
						col++;
						if (COMMENT) cout << "filling Line... col = " << col << "  item : " << Value << endl;
						// read constraint after output is read
						if (col == pr->nvar+1)
							Parse >> temp;
					}
					// comment
					if (COMMENT) {
						cout << "size of Line (constraint excluded) = "<< Line.size() << endl;
						cout << "constraint : " << temp << endl;
					}
					// checks
					if (Line.size() != pr->nvar+1) {
						cerr << "\nInequality constraint0 : Line " << row << ": expected " << (pr->nvar)+2 << " values, received " << Line.size()+1 << " values, aborting." << endl;
						cerr << "CHECK THE NUMBER OF VARIABLES DECLARED IN THE INPUT FILE!" << endl;
						exit(-1);
					}
					if (temp=='n') {
						cerr << "ERROR: Relationship not read: temp = " << temp << endl;
						exit(-1);
					}
					if (!((temp=='>') || (temp=='<'))) {
						cerr << "ERROR : constraint not recognised : temp = " << temp << endl;
						cerr << "Check 0-order inequality constraint in input file" << endl;
						exit(-1);
					}
					// show results
					if (COMMENT) cout << "DATA_INEQUALITY0[" << row << "]:" << endl;
					// values
					for (int i=0; i<(pr->nvar)+1; i++) {
						DATA_INEQUALITY0[row][i] = Line[i];
						if (COMMENT) {
							cout << scientific << DATA_INEQUALITY0[row][i] << "  ";
							cout << temp;
						}

					}
					// relationship
					CONSTRAINTS0[row] = temp;

					row++;
					col = 0;
					if (COMMENT) cout << endl;
			}
			else
				if (COMMENT) cout << "comment" << endl;
		}
	}


	//---------------------------------------------------------------------------------------------------
	// Dynamic allocation of matrix DATA_INEQUALITY1 (FITCASES rows, NVAR+1 columns): first allocate the row, then the columns for each row (no other way)
	//---------------------------------------------------------------------------------------------------
	COMMENT = 1;
	Val **DATA_INEQUALITY1 = NULL;
	char *CONSTRAINTS1 = NULL;
	if (pr->n_inequality1) {
		DATA_INEQUALITY1 = new Val*[pr->n_inequality1];
		if (DATA_INEQUALITY1 == NULL)  {
			cerr << "\nERROR: dynamic allocation of DATA_INEQUALITY1 failed... Exit" << endl;
			exit(-1);
		}
		for (int i=0; i < pr->n_inequality1; i++) {
			DATA_INEQUALITY1[i] = new Val[2*(pr->nvar)+1];
			if (DATA_INEQUALITY1[i]==NULL)  {
				cerr << "\nERROR: dynamic allocation of DATA_INEQUALITY1["<< i << "] failed... Exit" << endl;
				exit(-1);
			}
		}
		CONSTRAINTS1 = new char[pr->n_inequality1];
		// get DATA_INEQUALITY1 and CONSTRAINT1 reading line by line
		char temp;
		row = 0;
		col = 0;
		while (getline(fin, Input) && (row<pr->n_inequality1)) {
			if (COMMENT) cout << "\nResult of getline : " << Input << endl;
			// Input contains entire line
			if (Input[0]!='#') {
				stringstream Parse;
				//cout << Parse.str() << endl;
				Parse << Input;
				vector<double> Line;
				temp = 'n';
				while (Parse >> Value) {
					// only the values that are double are actually copied... so comment is recognised!
					Line.push_back( Value );
					col++;
					if (COMMENT) cout << "filling Line... col = " << col << "  item : " << Value << endl;
					// read constraint after output is read
					if (col == 2*(pr->nvar)+1)
						Parse >> temp;
				}
				// comment
				if (COMMENT) {
					cout << "size of Line (constraint excluded) = "<< Line.size() << endl;
					cout << "constraint : " << temp << endl;
				}
				// checks
				if (Line.size() != 2*(pr->nvar)+1) {
					cerr << "\nInequality constraint1 : Line " << row << ": expected " << 2*((pr->nvar))+2 << " values, received " << Line.size()+1 << " values, aborting." << endl;
					cerr << "CHECK THE NUMBER OF COLUMNS DECLARED IN INEQUALITY CONSTRAINT1 MATRIX!" << endl;
					exit(-1);
				}
				if (temp=='n') {
					cerr << "ERROR: Relationship not read: temp = " << temp << endl;
					exit(-1);
				}
				if (!((temp=='>') || (temp=='<'))) {
					cerr << "ERROR : constraint not recognised : temp = " << temp << endl;
					cerr << "Check 1-order inequality constraint in input file" << endl;
					exit(-1);
				}
				// show results
				if (COMMENT) cout << "DATA_INEQUALITY1[" << row << "]:" << endl;
				// values
				for (int i=0; i<2*(pr->nvar)+1; i++) {
					DATA_INEQUALITY1[row][i] = Line[i];
					if (COMMENT) {
						cout << scientific << DATA_INEQUALITY1[row][i] << "  ";
						cout << temp;
					}
				}
				// relationship
				CONSTRAINTS1[row] = temp;

				row++;
				col = 0;
				if (COMMENT) cout << endl;
			}
			else
				if (COMMENT) cout << "comment" << endl;
		}
	}



	//close the stream
	fin.close();




	//-----------------------------------------------------------------------------------------------------
	// all that follows is operations on read data (operations on file finished...)
	//-----------------------------------------------------------------------------------------------------


	// set problem definition attributes
		
//-+-+-+-+-+-+-+-+-+-+-+




//-+-+-+-+-+-+-+-+-+-+-+

	// BINARY FUNCTIONS : the names to be inserted in the input file are between " "
	pb->division = NULL;   //address of the "sdiv" operation, if among primitives
	for (int k=0; k<7; k++)   //the number (7) is the number of dummy operations set in Problem Definition class
		pb->b_func_list[k] = &(pb->dummy_bin[k]);
		
	NUM_B_FUNCS=0;
	if (func_bin.find("add") != string::npos) {
      pb->dummy_bin[NUM_B_FUNCS]=Add;
		NUM_B_FUNCS++;
	}
	if (func_bin.find("sub") != string::npos) {
	  	pb->dummy_bin[NUM_B_FUNCS]=Sub;
	 	NUM_B_FUNCS++;
	}
	if (func_bin.find("mult") != string::npos) {
      pb->dummy_bin[NUM_B_FUNCS]=Mult;
		NUM_B_FUNCS++;
	}
	if (func_bin.find("sdiv") != string::npos) {
      pb->dummy_bin[NUM_B_FUNCS]=SDiv;
      pb->division = pb->b_func_list[NUM_B_FUNCS];
		NUM_B_FUNCS++;
		cout << "\nDIVISION FOUND!" << endl;
	}
	if (func_bin.find("spow") != string::npos) {
      pb->dummy_bin[NUM_B_FUNCS]=Spow;
		NUM_B_FUNCS++;
	}
	pb->num_b_funcs = NUM_B_FUNCS;

	//cout << "\n read_file_new CHECK: pb->num_b_funcs = " << pb->num_b_funcs << endl;


	// UNARY FUNCTIONS : the names to be inserted in the input file are between " "
	for (int k=0; k<15; k++)   //the number (15) is the number of dummy operations set in Problem Definition class
		pb->u_func_list[k] = &(pb->dummy_uni[k]);
		
	NUM_U_FUNCS=0;
	if (func_uni.find("shift") != string::npos) {
    	pb->dummy_uni[NUM_U_FUNCS]=Shift;
		NUM_U_FUNCS++;
 	}
	if (func_uni.find("neg") != string::npos) {
		pb->dummy_uni[NUM_U_FUNCS]=Negation;
		NUM_U_FUNCS++;
 	}
	if (func_uni.find("square") != string::npos) {
     	pb->dummy_uni[NUM_U_FUNCS]=Square;
		NUM_U_FUNCS++;
 	}
	if (func_uni.find("cube") != string::npos) {
     	pb->dummy_uni[NUM_U_FUNCS]=Cube;
		NUM_U_FUNCS++;
	}
	if (func_uni.find("exp") != string::npos) {
     	pb->dummy_uni[NUM_U_FUNCS]=Exp;
		NUM_U_FUNCS++;
	}
	if (func_uni.find("nxp") != string::npos) {
     	pb->dummy_uni[NUM_U_FUNCS]=Negexp;
		NUM_U_FUNCS++;
	}
	if (func_uni.find("sin") != string::npos) {
     	pb->dummy_uni[NUM_U_FUNCS]=Sin;
		NUM_U_FUNCS++;
	}
	if (func_uni.find("cos") != string::npos) {
      	pb->dummy_uni[NUM_U_FUNCS]=Cos;
		NUM_U_FUNCS++;
	}
	if (func_uni.find("inv") != string::npos) {
      	pb->dummy_uni[NUM_U_FUNCS]=Inverse;
		NUM_U_FUNCS++;
	}
	if (func_uni.find("log") != string::npos) {
      	pb->dummy_uni[NUM_U_FUNCS]=Logn;
		NUM_U_FUNCS++;
	}
	if (func_uni.find("abs") != string::npos) {
      	pb->dummy_uni[NUM_U_FUNCS]=Absval;
		NUM_U_FUNCS++;
	}
	if (func_uni.find("cosh") != string::npos) {
      	pb->dummy_uni[NUM_U_FUNCS]=Cosh;
		NUM_U_FUNCS++;
	}
	if (func_uni.find("sinh") != string::npos) {
      	pb->dummy_uni[NUM_U_FUNCS]=Sinh;
		NUM_U_FUNCS++;
	}
	if (func_uni.find("tanh") != string::npos) {
      	pb->dummy_uni[NUM_U_FUNCS]=Tanh;
		NUM_U_FUNCS++;
	}
	if (func_uni.find("rectwave") != string::npos) {
	      	pb->dummy_uni[NUM_U_FUNCS]=RectWave;
			NUM_U_FUNCS++;
	}
	
	pb->num_u_funcs = NUM_U_FUNCS;


	// set ProblemDefinition private attributes (write a constructor?)
	pb->set_data(DATA);
	pb->set_n_data(pr->nfitcases);
	pb->set_n_var(pr->nvar);
	pb->set_n_cols(pr->nvar + 1);
	if (pr->crossvalidation==1) {
		pb->set_n_folds(pr->folds_n);
		pb->set_folds_table();  // dynamically allocate folds_table (see ProblemDefinition.h) for crossvalidation
	}
	if (pr->crossvalidation==0) {
		pb->set_n_folds(0);
		pb->set_folds_table();  // dynamically allocate folds_table (see ProblemDefinition.h) for crossvalidation
	}

	pb->data_tuning = DATA;			// as a default, tuning and evaluation data are set identical
	pb->n_tuning = pr->nfitcases;
	pb->data_validation = DATA;		// as a default, tuning and evaluation data are set identical
	pb->n_validation = pr->nfitcases;

	pb->data_inequality0 = DATA_INEQUALITY0;
	pb->n_inequality0 = pr->n_inequality0;
	pb->constraints0 = CONSTRAINTS0;
	pb->data_inequality1 = DATA_INEQUALITY1;
	pb->n_inequality1 = pr->n_inequality1;
	pb->constraints1 = CONSTRAINTS1;

	cout << "\n\nread_input_file : data correctly imported. Exit \n" << endl;
}



// function to read the secondary input file, containing the test data set
void read_test_data(string FILE_TEST_DATA,  RunParameters* pr, ProblemDefinition* pb)
{

	int COMMENT=0;

	cout << "\n\nread_test_data : enter" << endl;

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
		cerr << "read_test_data() : ERROR: The test data file " << FILE_TEST_DATA << " doesn't exist!!!" << endl;
		exit(-1);
	}
	else {
		cout << "read_test_data() : Reading from file " << FILE_TEST_DATA << "  : ..." << endl;
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

//*/*/*/*/ to be moved to ProblemDefinition
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
//*/*/*/*/

	cout << "read_test_data : points correctly imported = " << row << endl;
	cout << "read_test_file : test data correctly imported. Exit \n" << endl;

}




