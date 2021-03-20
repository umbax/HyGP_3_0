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

#include <cmath>
//#include <math.h>

using namespace std;

// included dependencies
#include "../modules/Val_type.h"		//just a header, no source file
#include "../modules/func_primitives_prototypes.h"



// BINARY FUNCTIONS (Binary_Func)

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
    // if b is zero or abs(b)<MIN_VAL, add a correction (update counter) - protected approach
    if ((b==(Val)0.) || (fabs(b)<MIN_VAL)) {
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



// UNARY FUNCTIONS (Unary_Func)

Val pow1(Val a, int* n_corrections)
{
	return a;
}
char pow1_sign_pre[] = "";
char pow1_sign_post[] = "^1";
int pow1_pos = 1;   //0 if BEFORE, 1 if AFTER the argument, 2 BEFORE AND AFTER the argument
Unary_Func Shift = {pow1, pow1_sign_pre, pow1_sign_post, &pow1_pos};


// it's strongly suggested not using this function, as can
// interfere with the tuning process
Val negation(Val a, int* n_corrections)
{
	return (-a);
}
char negation_sign_pre[]="(-1.0)*";
char negation_sign_post[]="";
int negation_pos = 0; //0 if BEFORE, 1 if AFTER the argument, 2 BEFORE AND AFTER the argument
Unary_Func Negation = {negation, negation_sign_pre, negation_sign_post, &negation_pos};


Val square(Val a, int* n_corrections)
{
	return (a*a);
}
char square_sign_pre[] = "";
char square_sign_post[] = "^2";
int square_pos = 1;   //0 if BEFORE, 1 if AFTER the argument, 2 BEFORE AND AFTER the argument
Unary_Func Square = {square, square_sign_pre, square_sign_post, &square_pos};

Val cube(Val a, int* n_corrections)
{
	return (a*a*a);
}
char cube_sign_pre[] = "";
char cube_sign_post[] = "^3";
int cube_pos = 1;      //0 if BEFORE, 1 if AFTER the argument, 2 BEFORE AND AFTER the argument
Unary_Func Cube = {cube, cube_sign_pre, cube_sign_post, &cube_pos};


Val exponential(Val a, int* n_corrections)
{
	Val r;
	r = exp(a);

	return r;
}
char exp_sign_pre[] = "exp";
char exp_sign_post[] = "";
int exp_pos = 0;      //0 if BEFORE, 1 if AFTER the argument, 2 BEFORE AND AFTER the argument
Unary_Func Exp = {exponential, exp_sign_pre, exp_sign_post, &exp_pos};


Val negexponential(Val a, int* n_corrections)
{
	Val r;
	r = exp(-a);

	return r;
}
char negexp_sign_pre[] = "1/exp";
char negexp_sign_post[] = "";
int negexp_pos = 0;      //0 if BEFORE, 1 if AFTER the argument, 2 BEFORE AND AFTER the argument
Unary_Func Negexp = {negexponential, negexp_sign_pre, negexp_sign_post, &negexp_pos};


Val sine(Val a, int* n_corrections)
{
	// a is expressed in radians!!!;
    return sin(a);
}
char sin_sign_pre[] = "sin";
char sin_sign_post[] = "";
int sin_pos = 0;      //0 if BEFORE, 1 if AFTER the argument, 2 BEFORE AND AFTER the argument
Unary_Func Sin = {sine, sin_sign_pre, sin_sign_post, &sin_pos};

Val cosine(Val a, int* n_corrections)
{
	// a is expressed in radians!!!;
    return cos(a);
}
char cos_sign_pre[] = "cos";
char cos_sign_post[] = "";
int cos_pos = 0;      //0 if BEFORE, 1 if AFTER the argument, 2 BEFORE AND AFTER the argument
Unary_Func Cos = {cosine, cos_sign_pre, cos_sign_post, &cos_pos};


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
char inverse_sign_pre[] = "";
char inverse_sign_post[] = "^(-1)";
int inverse_pos = 1;      //0 if BEFORE, 1 if AFTER the argument, 2 BEFORE AND AFTER the argument
Unary_Func Inverse = {inverse, inverse_sign_pre, inverse_sign_post, &inverse_pos};



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
char logn_sign_pre[] = "log";
char logn_sign_post[] = "";
int logn_pos = 0;      //0 if BEFORE, 1 if AFTER the argument, 2 BEFORE AND AFTER the argument
Unary_Func Logn = {logarithm, logn_sign_pre, logn_sign_post, &logn_pos};


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
char abs_sign_pre[] = "abs";
char abs_sign_post[] = "";
int abs_pos = 0;      //0 if BEFORE, 1 if AFTER the argument, 2 BEFORE AND AFTER the argument
Unary_Func Absval = {absolute, abs_sign_pre, abs_sign_post, &abs_pos};


Val hypercos(Val a, int* n_corrections)
{
	return cosh(a);
}
char hypercos_sign_pre[] = "cosh";
char hypercos_sign_post[]="";
int hypercos_pos = 0;      //0 if BEFORE, 1 if AFTER the argument, 2 BEFORE AND AFTER the argument
Unary_Func Cosh = {hypercos, hypercos_sign_pre, hypercos_sign_post, &hypercos_pos};


Val hypersin(Val a, int* n_corrections)
{
	return sinh(a);
}
char hypersin_sign_pre[] = "sinh";
char hypersin_sign_post[] = "";
int hypersin_pos = 0;      //0 if BEFORE, 1 if AFTER the argument, 2 BEFORE AND AFTER the argument
Unary_Func Sinh = {hypersin, hypersin_sign_pre, hypersin_sign_post,  &hypersin_pos};


Val hypertan(Val a, int* n_corrections)
{
	return tanh(a);
}
char hypertan_sign_pre[] = "tanh";
char hypertan_sign_post[] = "";
int hypertan_pos = 0;      //0 if BEFORE, 1 if AFTER the argument, 2 BEFORE AND AFTER the argument
Unary_Func Tanh = {hypertan, hypertan_sign_pre, hypertan_sign_post, &hypertan_pos};


// TEST: block defining high frequency sine term
Val hfreqsine(Val a, int* n_corrections)
{
	// a is expressed in radians!!!;
    return sin((2.0*PI/2.0)*a);		// no numerical parameter inserted here (see Unary_Node::check_allocation): frequency is set in an alternative way, it will not be optimised by Omegalim approach
}
char hfreqsine_sign_pre[] = "hfreqsin";
char hfreqsine_sign_post[] = "";
int hfreqsine_pos = 0;      //0 if BEFORE, 1 if AFTER the argument, 2 BEFORE AND AFTER the argument
Unary_Func HfreqSine = {hfreqsine, hfreqsine_sign_pre, hfreqsine_sign_post, &hfreqsine_pos};


Val rect_wave(Val a, int* n_corrections)
{
	// a is expressed in radians!!!;
	return (Val)cbrt(sin(a));		// used cbrt instead of pow() to avoid errors linked to the computation of 1/3
}
char rect_wave_sign_pre[] = "nthroot(sin";    //"rectwave"
char rect_wave_sign_post[] = ",3)";	// nthroot(x,3) command used by Matlab and Octave
int rect_wave_pos = 2;      //0 if BEFORE, 1 if AFTER, 2 BEFORE AND AFTER the argument
Unary_Func RectWave = {rect_wave, rect_wave_sign_pre, rect_wave_sign_post, &rect_wave_pos};


// Heaviside step function 1/(1+exp(-2kx)) ?
Val heavi_step(Val a, int* n_corrections)
{
	return (Val)(1/(1+exp(-2000.0*a)));		// used cbrt instead of pow() to avoid errors linked to the computation of 1/3
}
char heavi_step_sign_pre[] = "1/(1+exp(-2000*";
char heavi_step_sign_post[] = "))";
int heavi_step_pos = 2;      //0 if BEFORE, 1 if AFTER, 2 BEFORE AND AFTER the argument
Unary_Func HeaviStep = {heavi_step, heavi_step_sign_pre, heavi_step_sign_post, &heavi_step_pos};


