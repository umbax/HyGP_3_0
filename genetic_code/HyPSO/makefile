# Copyright 2021 Dr Umberto Armani
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Basic Makefile 
# all is the default target and is called when just makefile with
# no options is run.

all: hypso

hypso: ../SQP/MI0L2_c/M.o ../SQP/MI0L2_c/T.o  test_main.o hypso_launcher.o ../tree_functions/vector_derivative_functions.o ../tree_functions/tree_operations.o ../classes/problem_definition.o ../classes/run_parameters.o ../classes/class_POPULATION.o ../modules/variable.o  ../nodes/nodes.o ../modules/primitives.o 
	gfortran -std='legacy' -o hypso ../SQP/MI0L2_c/M.o ../SQP/MI0L2_c/T.o test_main.o hypso_launcher.o ../tree_functions/vector_derivative_functions.o ../tree_functions/tree_operations.o ../classes/problem_definition.o ../classes/run_parameters.o ../classes/class_POPULATION.o ../modules/variable.o  ../nodes/nodes.o ../modules/primitives.o -L/usr/lib/gcc/x86_64-redhat-linux/4.0.2/ -lstdc++ #-fopenmp 

test_main.o: test_main.cpp
	g++ -c -g test_main.cpp -o test_main.o

hypso_launcher.o: hypso_launcher.cpp
	g++ -c -g -W hypso_launcher.cpp -o hypso_launcher.o

../SQP/MI0L2_c/M.o: 
	gfortran -std='legacy' -c ../SQP/MI0L2_c/MI0L2.FOR -o ../SQP/MI0L2_c/M.o

../SQP/MI0L2_c/T.o: 
	gfortran -fsecond-underscore -c ../SQP/MI0L2_c/TI0L2.FOR  -o ../SQP/MI0L2_c/T.o

vector_derivative_functions.o: ./genetic_code/tree_functions/vector_derivative_functions.cpp
	g++ -c -g ./genetic_code/tree_functions/vector_derivative_functions.cpp -o ./genetic_code/tree_functions/vector_derivative_functions.o

tree_operations.o: ./genetic_code/tree_functions/tree_operations.cpp
	g++ -c -g ./genetic_code/tree_functions/tree_operations.cpp -o ./genetic_code/tree_functions/tree_operations.o

problem_definition.o: ./genetic_code/classes/problem_definition.cpp
	g++ -c -g ./genetic_code/classes/problem_definition.cpp -o ./genetic_code/classes/problem_definition.o

run_parameters.o: ./genetic_code/classes/run_parameters.cpp
	g++ -c -g ./genetic_code/classes/run_parameters.cpp -o ./genetic_code/classes/run_parameters.o

#read_file_new.o: ./genetic_code/read_input/read_file_new.cpp
#	g++ -c -g ./genetic_code/read_input/read_file_new.cpp -o ./genetic_code/read_input/read_file_new.o

population.o: ./genetic_code/classes/class_POPULATION.cpp
	g++ -c -g ./genetic_code/classes/class_POPULATION.cpp -o ./genetic_code/classes/class_POPULATION.o

#reporter.o:
#	g++ -c -g ./genetic_code/classes/reporter.cpp -o ./genetic_code/classes/reporter.o

#input_check.o: 
#	g++ -c -g ./genetic_code/input_checks/input_check.cpp -o ./genetic_code/input_checks/input_check.o

variable.o: 
	g++ -c -g ./genetic_code/modules/variable.cpp -o ./genetic_code/modules/variable.o

nodes.o:
	g++ -c -g ./genetic_code/nodes/nodes.cpp -o ./genetic_code/nodes/nodes.o

primitives.o:
	g++ -c -g ./genetic_code/modules/primitives.cpp -o ./genetic_code/modules/primitives.o
	# use extern to define Add, Sub, etc also mind the typedef struct renomination... use a class?
	
	
clean:
	rm ./test_main.o
	rm ./hypso_launcher.o
	rm ../SQP/MI0L2_c/M.o
	rm ../SQP/MI0L2_c/T.o
	rm ../tree_functions/vector_derivative_functions.o
	rm ../tree_functions/tree_operations.o
	rm ../classes/problem_definition.o
	rm ../classes/run_parameters.o
	rm ../classes/class_POPULATION.o
	rm ../modules/variable.o
	rm ../nodes/nodes.o
	rm ../modules/primitives.o
	rm ./hypso
	
	