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


#ifndef FUNC_PRIMITIVES_PROTOTYPES_H_
#define FUNC_PRIMITIVES_PROTOTYPES_H_

// functional primitives prototypes

// included dependencies
#include "../modules/Val_type.h"


// Binary functions
typedef struct Binary_Func {
    Val (*eval)(Val,Val,int*);               // the function pointer (the int* is used to increment the number of times a correction is done by protected operations)
    char *sign;                         // the infix operation sign
}  Binary_Func;


// Unary functions
typedef struct Unary_Func {
    Val (*eval)(Val, int*);                   // function pointer
    char *pre_sign;                         // the prefix operation sign (before the argument)
    char *post_sign;                        // the suffix operation sign (after the argument)
	int *pos;  							// the position of the operation sign (0 if BEFORE, 1 if after the argument)
}  Unary_Func;


#endif /* FUNC_PRIMITIVES_PROTOTYPES_H_ */
