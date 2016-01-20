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

/*
 * split_file_name.cpp
 *
 *  Created on: Mar 1, 2011
 *      Author: cnua
 */

void split_file_name (const string& str, string& folder, string& file)
{
  size_t found;
  //cout << "Splitting: " << str << endl;
  found = str.find_last_of("/\\");

  folder = str.substr(0,found);
  file = str.substr(found+1);
  //cout << " folder: " << str.substr(0,found) << endl;
  //cout << " file: " << str.substr(found+1) << endl;


}
