%%% Copyright 2021 Dr Umberto Armani
%%% 
%%% Licensed under the Apache License, Version 2.0 (the "License");
%%% you may not use this file except in compliance with the License.
%%% You may obtain a copy of the License at
%%%  
%%%       http://www.apache.org/licenses/LICENSE-2.0
%%%   
%%% Unless required by applicable law or agreed to in writing, software
%%% distributed under the License is distributed on an "AS IS" BASIS,
%%% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%%% See the License for the specific language governing permissions and
%%% limitations under the License.


%%%  HyGPview.m 
%%%  script that:
%%%  - plots the target function and the best tree of a HyGP experiment
%%%  - recomputes error metrics evaluating the HyGP model on the building or test data sets
%%%   conventions: R matrix containing the training data set (i.e. input file data)
%%%                Rtest matrix containing the validation data set
%%% Functions called:
%%% 1) get_experiment
%%% 2) read_INPUT_FILE 
%%% 3) read_TEST_FILE 
%%% 4) assign_f_obj_fun
%%% 5) count_runs
%%% 6) compute_expression
%%% 7) plot_act_vs_est
%%% 8) plot_model
%%% 9) plot autocorrelation function
%%% 10) plot objectives (after reading run_1/inputdata_stats.txt)
%%% Notes: 31/12/20 the script has been tested only for 1D and 2D problems... 

clear all;
close all;
format long e;

% settings and variables ----------------------------------------
plot_char = '-b';
size_title = 11;
pi=3.141592653589793e+00;
sep = filesep;  % set separator
objectives_build=struct('RMSE',0.0,'Corrections',0,'R2',0.0,'Mean',0.0,'Var',0.0,'Min',0.0,'Max',0.0,'ACF',0.0,'Tot_variation',0.0);
objectives_test=struct('RMSE',0.0,'Corrections',0,'R2',0.0,'Mean',0.0,'Var',0.0,'Min',0.0,'Max',0.0,'ACF',0.0,'Tot_variation',0.0);
objectives_target=struct('RMSE',0.0,'Corrections',0,'R2',0.0,'Mean',0.0,'Var',0.0,'Min',0.0,'Max',0.0,'ACF',0.0,'Tot_variation',0.0);
%-------------------------------------------------------

% GET THE NAME OF THE EXPERIMENT AND ITS DIRECTORY (function)
%first_dir = ['E:' sep 'Research_Leeds' sep 'Linux' sep 'code' sep 'C++' sep 'experiments' sep 'structuralGP_OLD' sep 'output' sep]
first_dir=['D:\relocated_users\relocated_umba\Documents\OneDrive - University of Leeds\Research\2020_Research_ML_AI_Vassili\output']
[experiment, directory_name] = get_experiment(first_dir);
%if experiment selection is cancelled, directory_name should be zero
%and nothing should happen
if (directory_name == 0)
  disp('Quitted normally.')
  return
end

% GET THE NUMBER OF RUNS IN THE EXPERIMENT
last_run = count_runs(directory_name);
if (last_run<1)
    disp('No runs in the selected experiment. Exit.')
    exit
end
    
obj_test_array(last_run)=objectives_test;

% READ INPUT FILE 
%[PAR, R, CONSTR0, CONSTR1]=read_INPUT_FILE(directory_name);
[NVAR, NFITCASES, R, CONSTR0, CONSTR1]=read_INPUT_FILE(directory_name); 
%NVAR = int16(PAR(2))
%NFITCASES = int16(PAR(8))

% READ FILE CONTAINING TEST DATA SET
%[NVAR_test, NTESTCASES, Rtest, found_test]=read_TEST_FILE(directory_name, NVAR);
[NVAR_test, NTESTCASES, Rtest, found_test]=read_TEST_FILE(directory_name);


% cycle (end at the last line of the script))
run=1;
while (run<=last_run)

    close all;
    clear R2_now;

    if (run)
        clear GEN F FITNESS R2 HITS EXPR;    
    end

    % check the number of generations
    filename = 'best_gp.txt'; % best model - to access the individuals in the archive use 'objectives.txt';
    file = [directory_name '/run_' num2str(run) '/' filename];
    % extract number of generations from file (for my GP implementations)
    [GEN] = textread(file,'%d %*[^\n]',-1,'bufsize',16382,'commentstyle','shell','endofline','\n');
    tot_gen = size(GEN,1) - 1;       %size(GEN,1) is the number of rows in the text file...
    disp(tot_gen)
    
    % ask for the generation of the best individual
    if (tot_gen==-1)
        disp('No generation in this run');
        break;
    end
  
    % load HyGP model expression
    % reach the right generation on the file
    fid = fopen(file);   % best_gp.txt
    count_gen=0;
    while (count_gen<tot_gen-1)
        % use text scan! text read starts every time from the beginning of the file!!!
        C = textscan(fid,'%u %*[^\n]',1,'commentstyle','#');
        count_gen = C{1};
    end    
    % read tree objectives(TREE_MEAN, TREE_VAR, TREE_MIN, TREE_MAX, etc) at the requested generation
    disp([experiment ': tree in run ' num2str(run) ', generation ' num2str(tot_gen) ])
    %C = textscan(fid,'%u %f %f %f %d %q',1,'bufsize',16382,'commentstyle','#'); %,'endofline','\n')
    C = textscan(fid,'%u %f %f %u %f %f %f %f %f %f %f %f %q',1,'commentstyle','#'); %,'endofline','\n')  % best_gp.txt
    fclose(fid);
    
    %# Gen F Fitness(RMSE) Corrections R2(adim.) Mean Var Min Max First_ACF Tot_variation Expression
    count_gen = C{1};
    F= C{2};
    RMSE_HyGP = C{3};
    CORRECTIONS = C{4};
    R2 = C{5};
    TREE_MEAN = C{6};
    TREE_VAR = C{7};
    TREE_MIN = C{8};
    TREE_MAX = C{9};
    TREE_FIRST_ACF = C{10};
    TREE_TOT_VARIATION = C{11};
    MAXABSERROR = C{12};
    EXPR = char(C{13});

    %
    run
    tree_string = EXPR 
    RMSE_HyGP   
    R2_HyGP = R2  
    disp(' ')

    % evaluate individual: compute the objectives on build and test data sets
    [tree_output_build, objectives_build]=evaluate_individual(NFITCASES, NVAR, R, tree_string);
    [tree_output_test, objectives_test]=evaluate_individual(NTESTCASES, NVAR, Rtest, tree_string);
    
    % store results in an array of struct
    obj_test_array(run)=objectives_test;  %% 27/3/23 ERROR! How to populate an array of struct!?
    
    run=run+1;
end

%clear all 

