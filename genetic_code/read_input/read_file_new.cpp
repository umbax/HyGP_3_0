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

// function to read the input file containing the HyGP hyperparameters and the building data set
// (the building data set can be split in tuning data set and evaluation data set (cross validation data set),
// but this is done later with the function SPLIT

void read_input_file(string FILE_INPUT,  RunParameters* pr, ProblemDefinition* pb)
{

	int COMMENT=0;

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
															"SPLIT=",
															"VALIDATING_LINES=",
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
	pr->split = (int)(p_par_values[26]);
	pr->validating_lines = (int)(p_par_values[27]);
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
	
	pb->num_u_funcs = NUM_U_FUNCS;


	// set ProblemDefinition private attributes (write a constructor?)
	pb->set_data(DATA);
	pb->set_n_data(pr->nfitcases);
	pb->set_n_var(pr->nvar);
	pb->set_n_cols(pr->nvar + 1);

	pb->data_tuning = DATA;			// as a default, tuning and evaluation data are set identical
	pb->n_tuning = pr->nfitcases;
	pb->data_evaluation = DATA;		// as a default, tuning and evaluation data are set identical
	pb->n_evaluation = pr->nfitcases;

	pb->data_inequality0 = DATA_INEQUALITY0;
	pb->n_inequality0 = pr->n_inequality0;
	pb->constraints0 = CONSTRAINTS0;
	pb->data_inequality1 = DATA_INEQUALITY1;
	pb->n_inequality1 = pr->n_inequality1;
	pb->constraints1 = CONSTRAINTS1;

	// compute stats on the whole input data set (tuning+evaluation)
	// in the future write a class "data set" with the following as a member function
	Val sum_output = .0;
	Val y_ave = .0;
	for (int k=0; k < pr->nfitcases; k++) {
		sum_output = sum_output + pb->get_data(k,pr->nvar)*pb->get_data(k,pr->nvar);
		y_ave = y_ave + pb->get_data(k,pr->nvar);
	}
	pb->sum_output = sum_output;
	pb->y_ave = y_ave/(double)(pr->nfitcases);
	
	// Mind that DATA is the whole data set imported as building data set,
	// so if you use cross validation on a subset of DATA (see data_evaluation) you need to recompute Sy
	// on data_evaluation
	Val Sy = 0.;
	for (int k=0; k < pr->nfitcases; k++)
		Sy = Sy + (pb->get_data(k,pr->nvar) - pb->y_ave)*(pb->get_data(k,pr->nvar) - pb->y_ave);
	pb->Sy = Sy;
	
}

