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


// functional primitives definition


// BINARY FUNCTIONS

// function for adding 2 values
Val add(Val a, Val b, int* n_corrections)
{
    return (a+b);                           /// variables have to be defined in the program or class
}
// the character sign for adding
char add_sign[] = " + ";
// the binary function representation
Binary_Func Add = {add, add_sign};

// function for subtraction
Val sub(Val a, Val b, int* n_corrections)
{
    return (a-b);
}
// char sign
char sub_sign[] = " - ";
// func rep
Binary_Func Sub = {sub, sub_sign};

// multiplication
Val mult(Val a, Val b, int* n_corrections)
{
    return (a*b);
}
char mult_sign[] = " * ";
Binary_Func Mult = {mult, mult_sign};

// safe division
Val sdiv(Val a, Val b, int* n_corrections)
{
    // if b is zero or abs(b)<MIN_VAL, add a correction (update counter)
    if ((b==(Val)0.) || (abs(b)<MIN_VAL)) {
		// check that n_corrections is a valid address...
		if (n_corrections)
			*n_corrections = *n_corrections+1;
      	return ((Val)1.0);
	}
	// otherwise
    return (a/b);
}
char sdiv_sign[] = " / ";
Binary_Func SDiv = {sdiv, sdiv_sign};


// power : this is for z^(other node not terminal_const) for example
Val spow(Val a, Val e, int* n_corrections)
{
    // a>0
    if (a>MIN_VAL)
		return (Val)pow(a,e);

	// a=0
	if (abs(a)<MIN_VAL) {
		if (e>MIN_VAL)
			return (Val)0.0;
		if (e<MIN_VAL) {
			// check that n_corrections is a valid address...
			if (n_corrections)
				*n_corrections = *n_corrections+1;
      		return ((Val)1.0);
		}
	}

	// a<0
	if (a<-MIN_VAL) {
		if (abs(e)<MIN_VAL)
			return (Val)1.0;
		if ((abs(e)>MIN_VAL) && (abs(e)<1.0)) {
			// check that n_corrections is a valid address...
			if (n_corrections)
				*n_corrections = *n_corrections+1;
      		return a;
		}
		if (abs(e)>=1)
			return (Val)pow(a,e);
	}
}
char spow_sign[] = "^";
Binary_Func Spow = {spow, spow_sign};




// UNARY FUNCTIONS

Val pow1(Val a, int* n_corrections)
{
	return a;
}
char pow1_sign[] = "^1";
int pow1_pos = 1;   //0 if BEFORE, 1 if AFTER the argument
Unary_Func Shift = {pow1, pow1_sign, &pow1_pos};


// it's strongly suggested not using this function, as can
// interfere with the tuning process
Val negation(Val a, int* n_corrections)
{
	return (-a);
}
char negation_sign[]="(-1.0)*";
int negation_pos = 0; //0 if BEFORE, 1 if AFTER the argument
Unary_Func Negation = {negation, negation_sign, &negation_pos};


Val square(Val a, int* n_corrections)
{
	return (a*a);
}
char square_sign[] = "^2";
int square_pos = 1;   //0 if BEFORE, 1 if AFTER the argument
Unary_Func Square = {square, square_sign, &square_pos};

Val cube(Val a, int* n_corrections)
{
	return (a*a*a);
}
char cube_sign[] = "^3";
int cube_pos = 1;      //0 if BEFORE, 1 if AFTER the argument
Unary_Func Cube = {cube, cube_sign, &cube_pos};


Val exponential(Val a, int* n_corrections)
{
	Val r;
	r = exp(a);

	return r;
}
char exp_sign[] = "exp";
int exp_pos = 0;      //0 if BEFORE, 1 if AFTER the argument
Unary_Func Exp = {exponential, exp_sign, &exp_pos};


Val negexponential(Val a, int* n_corrections)
{
	Val r;
	r = exp(-a);

	return r;
}
char negexp_sign[] = "1/exp";
int negexp_pos = 0;      //0 if BEFORE, 1 if AFTER the argument
Unary_Func Negexp = {negexponential, negexp_sign, &negexp_pos};


Val sine(Val a, int* n_corrections)
{
	// a is expressed in radians!!!;
    return sin(a);
}
char sin_sign[] = "sin";
int sin_pos = 0;      //0 if BEFORE, 1 if AFTER the argument
Unary_Func Sin = {sine, sin_sign, &sin_pos};

Val cosine(Val a, int* n_corrections)
{
	// a is expressed in radians!!!;
    return cos(a);
}
char cos_sign[] = "cos";
int cos_pos = 0;      //0 if BEFORE, 1 if AFTER the argument
Unary_Func Cos = {cosine, cos_sign, &cos_pos};


//Val tan(Val a, int* n_corrections)
//{
//	return tan(a);
//}
//char tan_sign[] = "tan";
//int tan_pos = 0;      //0 if BEFORE, 1 if AFTER the argument
//Unary_Func Tan = {tan, tan_sign, &tan_pos};



Val inverse(Val a, int* n_corrections)
{
	 // if abs(a)< MIN_VAL = 1.8e-19, add a correction (update counter)
    if ((a==(Val)0.) || (abs(a)<MIN_VAL)) {
		// check that n_corrections is a valid address...
		if (n_corrections)
			*n_corrections = *n_corrections+1;
      	return ((Val)1.0);  //quite problematic...
	}
	// otherwise
    return (((Val)1.0)/a);
}
char inverse_sign[] = "^(-1)";
int inverse_pos = 1;      //0 if BEFORE, 1 if AFTER the argument
Unary_Func Inverse = {inverse, inverse_sign, &inverse_pos};



Val logarithm(Val a, int* n_corrections)
{
	 // if a <0 or a=0, add a correction (update counter) and return -MAX_VAL
	 if (a<MIN_VAL) {
	 	// check that n_corrections is a valid address...
		if (n_corrections)
			*n_corrections = *n_corrections+1;
      	return ((Val)(-MAX_VAL));
	 }
	 // otherwise
	 return ((Val)log(a));
}
char logn_sign[] = "log";
int logn_pos = 0;      //0 if BEFORE, 1 if AFTER the argument
Unary_Func Logn = {logarithm, logn_sign, &logn_pos};


Val absolute(Val a, int* n_corrections)
{
	// a<0
	if (a<-MIN_VAL)
		return -a;
	// a>0
	if (a>MIN_VAL)
		return a;
	// a=0
	if (abs(a)<MIN_VAL)
		return 0.0;
}
char abs_sign[] = "abs";
int abs_pos = 0;      //0 if BEFORE, 1 if AFTER the argument
Unary_Func Absval = {absolute, abs_sign, &abs_pos};


Val hypercos(Val a, int* n_corrections)
{
	return cosh(a);
}
char hypercos_sign[] = "cosh";
int hypercos_pos = 0;      //0 if BEFORE, 1 if AFTER the argument
Unary_Func Cosh = {hypercos, hypercos_sign, &hypercos_pos};


Val hypersin(Val a, int* n_corrections)
{
	return sinh(a);
}
char hypersin_sign[] = "sinh";
int hypersin_pos = 0;      //0 if BEFORE, 1 if AFTER the argument
Unary_Func Sinh = {hypersin, hypersin_sign, &hypersin_pos};


Val hypertan(Val a, int* n_corrections)
{
	return tanh(a);
}
char hypertan_sign[] = "tanh";
int hypertan_pos = 0;      //0 if BEFORE, 1 if AFTER the argument
Unary_Func Tanh = {hypertan, hypertan_sign, &hypertan_pos};
