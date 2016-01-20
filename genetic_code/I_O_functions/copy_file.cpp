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


int copy_file(string PATH, string PATH_COPY)
{
	cout << "\nCopying input file ....";

	const char *expr, *expr_copy;
	expr = PATH.c_str();
	expr_copy = PATH_COPY.c_str();

	ifstream fin;
    ofstream fout;

    fin.open(expr);
    fout.open(expr_copy, ios_base::out | ios_base::trunc);
      
    if (fin == NULL || fout == NULL) {
    	cerr <<"\nERROR opening " << PATH << " or " << PATH_COPY;
        cerr << "\nExit";
        exit(-1);
    }
      
      // read from the first file then write to the second file
      char c;
      while(!fin.eof())
       {
            fin.get(c);
            fout.put(c);
       }
      
      //close the stream
      fin.close();
      fout.close();

      cout << "Ok";
}
