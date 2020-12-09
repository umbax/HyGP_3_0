# HyGP
C++ Hybrid Genetic Programming code for symbolic regression of explicit metamodels from data

Main features

Encoding: models or individuals are represented by trees with unary and binary operations

Selection: tournament out of three individuals selected from elite or whole population

Genetic operators: reproduction (elite), subtree crossover and mutation indipendently applied. Mutation is alternatively subtree mutation (even generations) and point mutation (odd generations). The reproduction operator replace copies in the elite with new individuals generated from scratch 

Fitness function: multiobjective definition through weighted approach or MinMax (non Pareto) - Strategy selected from STRAT_STATP parameter in input file. 
Objectives:
- individual root mean square error divided by average elite root mean square error
- number of numerical coefficients (squared) 
- number of illegal operations (i.e. division by zero)
- number of nodes (model size)
- variance and average of target function, computed from provided data
- max and min value of provided data

The parameter BOUNDED in input file allows the used to evolve bounded models (BOUNDED= 1), for extrapolation purposes.

Maximal depth restriction implemented to avoid generation of trees of excessive depth.

Termination criteria:
Evolution ends when the best individual root mean square error goes below a user defined threshold or the maximum number of generations is reached

The hybrid approach is implemented using a sequential quadratic programming (SQP) algorithm to tune the numerical coefficients of the individuals. The number of random initial guesses of the numerical coefficients can be set by the user.

Execution: sequential or parallel (SGE array job)
