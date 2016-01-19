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

// variable

// variable constructor
Variable::Variable(void)
//old constructor with default values
//(Val initial_val = 0.0, char initial_name[10] = "not set", \
		double initial_range = 0.0, Val initial_omega_lim = MAX_VAL)
{
	value = 0.0; //initial_val;
	strcpy(name, "not_set");
	lower_b = 0.0;
	upper_b = 0.0;
	range = 0.0; //initial_range;
	omega_lim = MAX_VAL; //initial_omega_lim;
}

void Variable::show_status(void)
{
	cout << "\n\nVariable \"" << name << "\" status :";
	cout << "\nvalue = " << value;
	cout << "\nname = " << name;
	cout << "\nlower_b = " << lower_b;
	cout << "\nupper_b = " << upper_b;
	cout << "\nrange = " << range;
	cout << "\nomega_lim = " << omega_lim << endl;
}
