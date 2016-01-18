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

gp: master.o T.o M.o
# for Linux on feng-gps1
	gfortran -std='legacy' -o gp M.o T.o master.o -L/usr/lib/gcc/x86_64-redhat-linux/4.0.2/ -lstdc++ -fopenmp
# for Linux on pre-pc1017 (NOTE the version number! 11/2010)
# g77 -o gp M.o T.o master.o -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2/ -lstdc++ 
# for development under Windows (Eclipse)   D:/MinGW/lib/gcc/mingw32/4.5.0/
#	g77 -o gp M.o T.o master.o -L d:/airbus_epd_py25_v4230201/lib/site-packages/mingw-3.4.5n1-py2.5-win32.egg/EGG-INFO/usr/bin/../lib/gcc/mingw32/3.4.5/ -lstdc++ 

M.o: ./genetic_code/SQP/MI0L2_c/MI0L2.FOR
	gfortran -std='legacy' -c ./genetic_code/SQP/MI0L2_c/MI0L2.FOR -o M.o

T.o: ./genetic_code/SQP/MI0L2_c/TI0L2.FOR
	gfortran -fsecond-underscore -c ./genetic_code/SQP/MI0L2_c/TI0L2.FOR  -o T.o
	
master.o: master.cpp
	# parallel compilation
	# g++ -c -g parallel_master.cpp -o master.o -fopenmp
	# normal compilation
	g++ -c -g master.cpp -o master.o
	
clean :
	rm M.o
	rm T.o
	rm master.o
	rm gp
