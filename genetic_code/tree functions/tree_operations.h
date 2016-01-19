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


Node *tree_copy(Node *, Node *);
int tree_comp_F(const void *, const void *);
int tree_comp_fitness(const void *, const void *);
void insert_parameter(Node *, Node **, Binary_Func *, Val);
int insert_subtree(Node **, Node *, Node *);    
void eliminate_subtree(Node *, Val);
void edit_tree(Node *,  Node**);
