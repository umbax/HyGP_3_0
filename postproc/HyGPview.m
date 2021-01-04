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
%%% Notes: 31/12/20 the script has been tested only for 1D and 2D problems... 

clear all;
close all;
format long e;

% settings ---------------------------------------
plot_char = '-b';
size_title = 11;
pi=3.141592653589793e+00;
sep = filesep;  % set separator 
%-------------------------------------------------------

% get the name of the experiment and its directory (function)
%first_dir = ['~' sep 'code' sep 'C++' sep 'experiments' sep 'structuralGP_OLD' sep 'output' sep]
first_dir = ['E:' sep 'Research_Leeds' sep 'Linux' sep 'code' sep 'C++' sep 'experiments' sep 'structuralGP_OLD' sep 'output' sep]
[experiment, directory_name] = get_experiment(first_dir);
%if experiment selection is cancelled, directory_name should be zero
%and nothing should happen
if (directory_name == 0)
  disp('Quitted normally.')
  return
end

% READ INPUT FILE 
[PAR, R, CONSTR0, CONSTR1]=read_INPUT_FILE(directory_name);
NVAR = int16(PAR(2))
NFITCASES = int16(PAR(7))

% READ FILE CONTAINING TEST DATA SET
[NVAR_test, NTESTCASES, Rtest, found_test]=read_TEST_FILE(directory_name, NVAR);

% SELECT OBJECTIVE FUNCTION
f_obj = assign_f_obj_fun();
f_obj


% cycle (end at the last line of the script))
another='0';
while (another~='n')

close all;
clear R2_now;

if (another~='0')
    clear GEN F FITNESS R2 HITS EXPR;    
end

% check the number of runs and ask the run no.
last_run = count_runs(directory_name);
message = [ num2str(last_run) ' runs in the selected experiment. Which? '];
run=-1;
while ((run==-1) || (run>last_run))
    run = input(message);
    if (run>last_run)
        disp('Not available!!')
    end
end

% check the number of generations
filename = 'best_gp.txt'; % to access the individuals in the archive use 'objectives.txt';
file = [directory_name '/run_' num2str(run) '/' filename];
% extract number of generations from file (for my GP implementations)
[GEN] = textread(file,'%d %*[^\n]',-1,'bufsize',16382,'commentstyle','shell','endofline','\n');
tot_gen = size(GEN,1) - 1;       %size(GEN,1) is the number of rows in the text file...

% ask for the generation of the best individual
if (tot_gen==-1)
    disp('No generation in this run');
    break;
end
if (tot_gen==0)
    st = ['Only the initial generation (0) in the current run.'];
else 
   st = [ num2str(tot_gen) ' generations in the current run.'];
end

st = 'Which generation? ';
cur_gen=-1;
while ((cur_gen==-1) || (cur_gen>tot_gen))
    cur_gen = input(st);
    if (cur_gen>tot_gen)
        disp('Not available!!')
    end
end


% load HyGP model expression
% reach the right generation on the file
fid = fopen(file);   % best_gp.txt
count_gen=0;
cur_gen
while (count_gen<cur_gen-1)
    % use text scan! text read starts every time from the beginning of the file!!!
    C = textscan(fid,'%u %*[^\n]',1,'commentstyle','#');
    count_gen = C{1};
end    
% read tree at the requested generation
disp([experiment ': tree in run ' num2str(run) ', generation ' num2str(cur_gen) ])
%C = textscan(fid,'%u %f %f %f %d %q',1,'bufsize',16382,'commentstyle','#'); %,'endofline','\n')
C = textscan(fid,'%u %f %f %f %d %q',1,'commentstyle','#'); %,'endofline','\n')  % best_gp.txt
cur_gen
count_gen = C{1}
F= C{2}
FITNESS = C{3}
R2 = C{4}
HITS = C{5}
EXPR = char(C{6})
tree_string = EXPR 
fitness = FITNESS   
R2_now = R2  
disp(' ')
% close stream fid
fclose(fid);
% % correct tree expression to make it MATLAB readable
% tree = correct_equation(tree_string);  % 29/12/20 not needed if compute_expression is used
% print the latex expression of tree. 
% ATTENTION: Nesting of {, [, and ( cannot exceed a depth of 32. 
% 8/9/20 Turned off to avoid the problem.
%disp(' ')
%disp('Tree expression in LATEX:')
%tree_latex=print_latex_tree(NVAR,'Z',tree_string,'n');
%disp(tree_latex)
%title('FontSize',size_title,'Interpreter','latex');


%---------------------------------------------
% Figure 2 : visualize DoE in the design space
%---------------------------------------------
if (NVAR==1)
    figure(2);  
    plot(R(:,1),zeros(NFITCASES,1), 'or','MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize',4);
    xlabel('Z1')
end
if (NVAR==2) 
    figure(2);  
    plot(R(:,1),R(:,2), 'or','MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize',4);
    xlabel('Z1')
    ylabel('Z2')
end
%legend('validating and building pts',  'Location', 'SouthOutside')
%axis tight equal;
title({experiment; 'Selected points in the domain'},'FontSize',size_title,'Interpreter','none');


%------------------------------------------------------------------------------------------
% Figure 3 : plot objective function and selected HyGP model on building and test data sets
% (RMSE and R2 are read from 'best_gp.txt' and currently only evaluated on building data set
%------------------------------------------------------------------------------------------
% plot model on building data set
plot_model(experiment, run, cur_gen, fitness, R2_now, NVAR, R, tree_string, f_obj, 'Building');
% plot model on test data set
if (found_test) 
    plot_model(experiment, run, cur_gen, -9999, -9999, NVAR_test, Rtest, tree_string, f_obj, 'Test');
end

%------------------------------------------------------------------------
% ACTUAL VS ESTIMATED RESPONSES : bring RMSE, R2, min, max computations
% in this (main) file?
%------------------------------------------------------------------------
% plot actual vs estimated response on building data set
plot_act_vs_est(NFITCASES, NVAR, R, tree_string, experiment, run, cur_gen, 'Building');
% plot actual vs estimated response on test data set (if found)
if (found_test) 
    plot_act_vs_est(NTESTCASES, NVAR, Rtest, tree_string, experiment, run, cur_gen, 'Test');
end
disp(['Shown image: run ' num2str(run) ', generation ' num2str(cur_gen) ]); 


% cycle
another=input('Another operation on the same experiment?(y/n)','s');
end


clear all;  
