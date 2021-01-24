# HyGP
C++ Hybrid Genetic Programming code for symbolic regression of explicit metamodels from data

Main features

### Memetic approach
The memetic/hybrid approach is implemented using a sequential quadratic programming (SQP) algorithm to tune the numerical coefficients of the individuals. The number of random initial guesses of the numerical coefficients can be set by the user.

### Encoding
Models or individuals are represented by trees with unary and binary operations. The root node is always a binary node.

### Selection
Tournament out of three individuals selected from elite or whole population

### Genetic operators
The implemented genetic operators are:
1. reproduction (copy of the *elite* unchanged to new generation. In case of copies, the reproduction operator replaces copies in the elite with new individuals generated from scratch)
1. crossover
1. mutation (alternatively subtree mutation (even generations) and point mutation (odd generations)). 
The genetic operators are independently applied, that is a model can be subjected to only one of the three operators at each generation.

### Fitness function
Multiobjective definition through weighted approach or MinMax (non Pareto). The objectives used to control the shape and behaviour of the evovled models are:

* individual root mean square error divided by average elite root mean square error
* number of numerical coefficients (squared) 
* number of illegal operations (i.e. division by zero)
* number of nodes (model size)
* variance and average of target function, computed from provided data
* max and min value of provided data
* first root of autocorrelation function (for 1D problems)

The "strategies" that define how the above mentioned objectives are combined to define the fitness function F can be selected from STRAT_STATP parameter in input file:

###### Weighted approach:
Given RMSE(i,t) the Root Mean Square Error of the i-th HyGP model at generation t as evaluated on the building data set (declared in the input file), N_(tuning coeff) the number of numerical coefficient in the model, N_(illegal op) the number of illegal operations found in the model, N_nodes the number of nodes of the model  

	F = a_1 RMSE(i,t)/RMSE(i,t-1) + a_2 N_(tuning coeff)+ a_3 N_(illegal op) + a_4 N_nodes+ a_8 F_8

  * Strategy 6 (STRAT_STATP=6):
  
     F_8=(sqrt(|input_data_variance - tree_variance|))^3 + |input_data_mean – tree_mean|^3
	
  * Strategy 7 (STRAT_STATP=7):
	
	F_8=(sqrt(|input_variance - tree_variance|))^2 + |input_mean – tree_mean|^2
	
  * Strategy 8 (STRAT_STATP=8): 
	
	F_8 = (sqrt(|input_variance - tree_variance|))^3 + |input_mean – tree_mean|^3 + |input_max – tree_max|^3 + |input_min – tree_min|^3
	
  * Strategy 9 (STRAT_STATP=9): 
	
	F_8 =(sqrt(|input_variance - tree_variance|))^3 + |input_mean – tree_mean|^3 + |input_max – tree_max|^3 + |input_min – tree_min|^3 + a9 * diverging (1,0)
	
  * Strategy 10 (STRAT_STATP=10):
   
	as Strategy 8, but with editing of high level polynomials if found before fitness function evaluation. Diverging terms are replaced by numerical coefficients.

  * Strategy 13 (STRAT_STATP=13):
  
     F = a_1 RMSE(i,t)/RMSE(i,t-1) + a_2 N_(tuning coeff)+ a_3 N_(illegal op) + a_4 N_nodes + a_8 F_8 + a_10 F_10

  with:
     
     F8 = (sqrt|input_variance- tree_variance|)^3 + |input_mean - tree_mean|^3
     F_10 = (sqrt|first_acf_root_input-first_acf_root_tree|)^3

###### MinMax approach:

  * Strategy 11 (STRAT_STATP=11):
     
     F = Min{Max[a1 RMSE(i,t)/RMSE(i,t-1), a2 N_(tuning coeff), a3 N_(illegal op), a_4 N_nodes, a5 F8 as in Strategy 8 ]}

  




### Bounded models
The parameter BOUNDED in input file forces HyGP to evolve bounded models (BOUNDED= 1), for extrapolation purposes.


### Termination criteria
Evolution ends when the best individual root mean square error goes below a user defined threshold or the maximum number of generations is reached

### Further measures
Maximal depth restriction implemented to avoid generation of trees of excessive depth.

### Execution
Sequential or parallel (SGE array job)
