/*
 * primitives.h
 *
 *  Created on: Sep 22, 2016
 *      Author: umba
 */

#ifndef GENETIC_CODE_MODULES_PRIMITIVES_H_
#define GENETIC_CODE_MODULES_PRIMITIVES_H_


// included dependencies
#include "../modules/func_primitives_prototypes.h"

// Binary functions
extern Binary_Func Add;
extern Binary_Func Sub;
extern Binary_Func Mult;
extern Binary_Func SDiv;
extern Binary_Func Spow;

// Unary functions
extern Unary_Func Shift;
extern Unary_Func Negation;
extern Unary_Func Square;
extern Unary_Func Cube;
extern Unary_Func Exp;
extern Unary_Func Negexp;
extern Unary_Func Sin;
extern Unary_Func Cos;
extern Unary_Func Inverse;
extern Unary_Func Logn;
extern Unary_Func Absval;
extern Unary_Func Cosh;
extern Unary_Func Sinh;
extern Unary_Func Tanh;
extern Unary_Func RectWave;
extern Unary_Func HfreqSine;

#endif /* GENETIC_CODE_MODULES_PRIMITIVES_H_ */
