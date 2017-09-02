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
 * single_run.cpp
* it is basically the independent thread that can be launched independently
 *  Created on: Mar 2, 2011
 *      Author: cnua
 */
int single_run (int run_id, int cur_run, int* SEED, string OUTPUT_STRING, \
		RunParameters* p_parameters,\
		ProblemDefinition* p_problem, \
		Reporter pop_reporter,\
		time_t start, time_t finish, double delta_t, int argc)
{
	// create directory run_k
	char num_field[10];
	sprintf(num_field, "%d", cur_run+1);
	string num_field_s = num_field;
	string DIR_RUN_K =  OUTPUT_STRING + "/run_" + num_field;
	mkdir(DIR_RUN_K.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );

	// seed the random generator
	srand(SEED[cur_run]);

	// start time
	time(&start);
	int elapsed_time;

	//---------------------------------------------------------------------
	// CREATE THE POPULATION (constructor called)
	//---------------------------------------------------------------------
	Variable* Z;    //to be included in ProblemDefinition...
	p_problem->initialise_variables(&Z, p_parameters->max_n_periods);

	// do you really need new? no, but don't touch it for now
	Population* P = new Population(p_parameters, p_problem);
	if (!P) {
		cerr << "\nmain : Error creating population!!!\n";
		exit(-1);
	}

	//as it's hard to pass Pop as a parameter to fdf_c__ through fortran functions, treat it as a global variable
	Pop = P;
	cout << "\nP = " << P;
	cout << "\nPop = " << Pop;

	///// INITIAL GENERATION (0) ///////////////////////////////////////

	// trees before parameters insertion
			cout << "\nInitialization of the population (generation 0)\n";
			P->print_population_without_parameters(0);


			/// split the data set in tuning set (data_tuning) and validation set (data_validation). See SPLIT and VALIDATING_LINES in input file
			//this function will also allow to increase the number of fitness cases during the run...
			//P->split_data(p_parameters, p_problem, 0,1);
			// here or before parallel section begins?
			// Depends if you want to change fitness cases during the evolution...


			/////////// OTHER GENERATIONS ///////////////////////////////////////
			int check_end = 0;
			int last_gen=0;
			for (int i=0; i<p_parameters->G+1; i++) {   // generations

				if (i) {   //skip generation 0
					// split the data for the current generation (this function will also allow to increase the number of fitness cases during the run...)
					//P->split_data(i,G,split); // not used...
				/*
					if ((i%6)==0) {   //6
						// KILLING and FILLING
						P->kill_and_fill(&problem);
						//cin.get();
					}
					else
				//*/
					// GENETIC OPERATORS: sorting, reproduction, crossover, mutation, pruning
					P->new_spawn(*p_parameters, *p_problem, p_parameters->nfitcases,i);

					// print population WITHOUT parameters after genetic operations
					//if (COMMENT) {
					//	printf("\n\n***********Generation %d after genetic operations (not sorted, new trees marked with f9.999999E+99)************\n", i-1);
					//	P->print_population_without_parameters(i-1);
					//	printf("\n***********************************************************************\n");
					//}
				}

				// evaluate fitness function (in structural GP parameters are added and tuned first, then the evaluation is performed)
				P->evaluate(i,p_parameters->G);

				//----------------------------------------------------------------------------------
				// extfitness - keeping the individual with least fitness value -------
				// sort according to fitness (error)
				//P->sort(i,tree_comp_fitness);

				///printf("\n\n***********Generation %d after genetic operations (not sorted, new trees marked with f9.999999E+99)************\n", i-1);
				//P->print_population_without_parameters(i-1);
				//printf("\n***********************************************************************\n");

				// update the best individual - structure and complete tree (for PARAMETER INHERITANCE)
				//P->update_ext_archive();
				//---------------------------------------------------------------------------------------------

				// sort according to F
				// VITAL! Both populations must be sorted, trees[] and complete_trees[]...
				P->sort(i,tree_comp_F);

				///printf("\n\n***********Generation %d after genetic operations (not sorted, new trees marked with f9.999999E+99)************\n", i-1);
				///P->print_population_without_parameters(i-1);
				///printf("\n***********************************************************************\n");

				// update the best individual - structure and complete tree (for PARAMETER INHERITANCE)
				P->update_ext_archive();

				// compute elapsed time
				elapsed_time = (int)(P->compute_time(start, finish, &delta_t));

				if (VERBOSE) {
					// print elapsed time
					cout << "Elapsed time: " << elapsed_time << " sec";  //total seconds

					// print out the best member - population WITHOUT parameters
					P->print_population_without_parameters(i);

					// print out the best member - population WITH parameters
					P->print_population_with_parameters(i);
				}

				// compute statistical data relating to population (vital if data is shared through populations)
				// IT's REALLY IMPORTANT that this function is executed after evaluate and sort,
				// as evaluate uses statistical data referring to the previous generation
				P->compute_statistics();

				// evaluate termination condition
				check_end=P->terminate(p_parameters->threshold);
				last_gen = i;

				// for the split data set, re-tune and re-evaluate the individuals on the merged data set (2017:REPLACED BY CROSSVALIDATION!!)
//				//if (p_parameters->split) {
//				//	if ((check_end) || (i==p_parameters->G)) {
//				//		cout << "\nBest Individual re-tuning and re-evaluation on the whole dataset" << endl;
//				//		P->split_data(p_parameters, p_problem,last_gen,last_gen);
//				//		P->evaluate(i,p_parameters->G);
//				//	}
//				//}

				// PRINT TO FILE OPERATIONS (in case of crash, data is saved)  -------------------------
				// print to file (if termination criterion met, it closes the stream of data to file )
				pop_reporter.stats2file(p_parameters, P, DIR_RUN_K, i, check_end);
				// write the test-points (training set), and the corresponding values of g_obj and the best individual (only one)
				pop_reporter.points2file(p_parameters, p_problem, P, DIR_RUN_K, i, check_end, start, finish, delta_t,SEED[cur_run]);
				// write best individual's expression ATTENZIONE: cambiato il 7/4/2015... da verificare
				pop_reporter.update_best2file_build(P, DIR_RUN_K, i, check_end); //old function: update_best2file(P, DIR_RUN_K, i, check_end);
				// update the list of the best-so-far individuals (see elite or archive) - truncation! ATTENZIONE: cambiato il 7/4/2015... da verificare
				pop_reporter.archive2file_build(P, DIR_RUN_K, i, check_end);  //old function: objective_table2file(P, DIR_RUN_K, i, check_end);
				// update no of tree evaluations
				pop_reporter.n_tree_eval2file(P, DIR_RUN_K, i, check_end);
				// -------------------------------------------------------------------------------------

				if (check_end)
					break;
			}

			# pragma omp critical
			{
				cout << "Elapsed time to complete run " << cur_run+1 << ": " << elapsed_time << " sec" << endl;  //total seconds
				//termination criterion satisfied
				if (check_end) {
					cout << "Termination criterion satisfied (RMSE < " << p_parameters->threshold << ")." << endl;
					cout << "Possible solution: " << endl;
					P->print_population_with_parameters(last_gen);
					cout << "Check latest_archive.txt for solutions\n" << endl;
				}
				else {
					P->print_population_with_parameters(last_gen);
				}
			}

			// just for test
			//P->get_tree_derivative_given_norm_vector(problem, P->complete_trees[0]);

			//---------------------------------------------------------------------
			// END OPERATIONS
			//---------------------------------------------------------------------
			// write node statistics to file
			pop_reporter.node_stats2file(p_parameters, P, DIR_RUN_K);


			//here operations to evaluate error on TEST DATA SET
			// evaluate fitness (RMSE and R2) on test data set (only if test data has been provided)
			if (argc==4) {
				// show data_test
				cout << "problem->show_data_validation() : show current data_test :" << endl;
				P->problem->show_data_test();
				// evaluate complete individuals on test data set provided by the user
				P->evaluate_complete_trees(); // SET CORRECTLY Mprobl.data_test, n_test, Sy_test after implementing function to read test data set
				// sort according to error (RMSE) - better not to use it to keep order and to recognise performance on building and test data sets...
				//P->sort(last_gen,tree_comp_fitness); // non va: ordina in ordine decrescente e alcune volte pone a 0 RMSE e R2. Perché?
				// find and print best individual on test data set to file
				pop_reporter.best2file_test(P, DIR_RUN_K, last_gen);
				// print archive evaluated on the test data set to file
				pop_reporter.archive2file_test(P, DIR_RUN_K, last_gen);  // insert a function to order in rmse decreasing order, leaving however the name of the run
			}


			// end - free memory allocated to Population
			delete P;


}
