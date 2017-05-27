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


#ifndef VAL_TYPE_H_
#define VAL_TYPE_H_

typedef double Val;

//// used in problem_definition.h
//typedef struct fold_t {
//  int fold_id;   // fold number
//  Val* p_data_row;  // pointer containing the address of the corresponding row in data
//} fold_t;

// definition of max and minimum absolute values and pi
#define MAX_VAL 1.8e+19
#define MIN_VAL 1.8e-19
#define PI 3.14159265358979323846264338327950288419716939937510582

#endif /* VAL_TYPE_H_ */
