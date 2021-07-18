# Copyright 2017 Dr Umberto Armani
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
all: test

test: test.o ../genetic_code/SQP/MI0L2_c/M.o ../genetic_code/SQP/MI0L2_c/T.o 
# for Linux on feng-gps1
	gfortran -std='legacy' -o test ../genetic_code/SQP/MI0L2_c/M.o ../genetic_code/SQP/MI0L2_c/T.o ./test.o -L/usr/lib/gcc/x86_64-redhat-linux/4.0.2/ -lstdc++ -fopenmp

test.o: ./test_master.cpp
	g++ -c -g ./test_master.cpp -o ./test.o

../genetic_code/SQP/MI0L2_c/M.o: ../genetic_code/SQP/MI0L2_c/MI0L2.FOR
	gfortran -std='legacy' -c ../genetic_code/SQP/MI0L2_c/MI0L2.FOR -o ../genetic_code/SQP/MI0L2_c/M.o

../genetic_code/SQP/MI0L2_c/T.o: ../genetic_code/SQP/MI0L2_c/TI0L2.FOR
	gfortran -fsecond-underscore -c ../genetic_code/SQP/MI0L2_c/TI0L2.FOR  -o ../genetic_code/SQP/MI0L2_c/T.o

clean :
	rm ./test.o
	rm ./test
	
	