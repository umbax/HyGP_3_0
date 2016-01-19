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
 * seed_generator.cpp
 *
 *  Created on: Mar 10, 2011
 *      Author: cnua
 */

void seed_generator(int **p_SEED, int N_RUNS, int input_seed)
{
	*p_SEED = new int[N_RUNS];

	if (!p_SEED) {
		cerr << "main : ERROR: Can't allocate SEED array\n";
		cerr << "Exit";
		exit(-1);
	}

	if (input_seed<0) {
		// if seed = -1 use random seed (time)
	    time_t *tp = NULL; // seed the random value generator
		srand((unsigned int) time (tp));    //to be used normally
		//parameters.seed = time(tp);
		cout << "\nSEED= -1 in input file: seeds randomly generated :";
		for (int k=0; k<N_RUNS; k++) {
			(*p_SEED)[k] = (rand() % 100000)*100+ (rand() % 100000);
			cout << "\nrun " << k << " : " << (*p_SEED)[k];
		}
	} else {
		// if seed >0 use the value given in input file
		cout << "\nSEED>0 in input file: used same seed for all runs:";
		for (int k=0; k<N_RUNS; k++) {
			(*p_SEED)[k] =  input_seed; //srand(parameters.seed); use this in the single evolution
			cout << "\nrun " << k << " : " << (*p_SEED)[k];
		}
	}
	cout << endl;
}

/*
int main ()
{
	cout << "\nSetting seeds for random number generator...";
	int* SEED;
	seed_generator(SEED, 5, -1);
	cout << "seed = 12" << endl;
	seed_generator(SEED, 5, 12);
}
*/
