# Copyright 2016 Dr Umberto Armani
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
all: gp

gp: master.o T.o M.o ./genetic_code/input_checks/input_check.o ./genetic_code/modules/variable.o ./genetic_code/nodes/nodes.o ./genetic_code/modules/primitives.o
# for Linux on feng-gps1
	gfortran -std='legacy' -o gp M.o T.o master.o ./genetic_code/input_checks/input_check.o ./genetic_code/modules/variable.o ./genetic_code/nodes/nodes.o ./genetic_code/modules/primitives.o -L/usr/lib/gcc/x86_64-redhat-linux/4.0.2/ -lstdc++ -fopenmp
# for Linux on pre-pc1017 (NOTE the version number! 11/2010)
# g77 -o gp M.o T.o master.o -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2/ -lstdc++ 
# for development under Windows (Eclipse)   D:/MinGW/lib/gcc/mingw32/4.5.0/
#	g77 -o gp M.o T.o master.o -L d:/airbus_epd_py25_v4230201/lib/site-packages/mingw-3.4.5n1-py2.5-win32.egg/EGG-INFO/usr/bin/../lib/gcc/mingw32/3.4.5/ -lstdc++ 

M.o: ./genetic_code/SQP/MI0L2_c/MI0L2.FOR
	gfortran -std='legacy' -c ./genetic_code/SQP/MI0L2_c/MI0L2.FOR -o M.o

T.o: ./genetic_code/SQP/MI0L2_c/TI0L2.FOR
	gfortran -fsecond-underscore -c ./genetic_code/SQP/MI0L2_c/TI0L2.FOR  -o T.o
	
master.o: parallel_master.cpp
	# parallel compilation (OpenMP)
	#g++ -c -g parallel_master.cpp -o master.o -fopenmp
	# normal compilation
	g++ -c -g master.cpp -o master.o

input_check.o: 
	g++ -c -g ./genetic_code/input_checks/input_check.cpp -o ./genetic_code/input_checks/input_check.o

variable.o: 
	g++ -c -g ./genetic_code/modules/variable.cpp -o ./genetic_code/modules/variable.o

nodes.o:
	g++ -c -g ./genetic_code/nodes/nodes.cpp -o ./genetic_code/nodes/nodes.o

primitives.o:
	g++ -c -g ./genetic_code/modules/primitives.cpp -o ./genetic_code/modules/primitives.o
	# use extern to define Add, Sub, etc also mind the typedef struct renomination... use a class?

clean :
	rm T.o
	rm M.o
	rm ./genetic_code/input_checks/input_check.o
	rm ./genetic_code/modules/variable.o
	rm ./genetic_code/nodes/nodes.o
	rm ./genetic_code/modules/primitives.o
	rm ./master.o
	rm ./gp
	rm ./gp_openmp
	
