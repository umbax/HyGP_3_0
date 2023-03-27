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

% READ INPUT FILE 
%[PAR, R, CONSTR0, CONSTR1]=read_INPUT_FILE(directory_name);
[NVAR, NFITCASES, R, CONSTR0, CONSTR1]=read_INPUT_FILE(directory_name); 
%NVAR = int16(PAR(2))
%NFITCASES = int16(PAR(8))

% READ FILE CONTAINING TEST DATA SET
%[NVAR_test, NTESTCASES, Rtest, found_test]=read_TEST_FILE(directory_name, NVAR);
[NVAR_test, NTESTCASES, Rtest, found_test]=read_TEST_FILE(directory_name);

% SELECT OBJECTIVE FUNCTION
f_obj = assign_f_obj_fun();

% ASK WHICH RUN TO PLOT
message = [ num2str(last_run) ' runs in the selected experiment. Which? '];
run=-1;
while ((run==-1) || (run>last_run))
    run = input(message);
    if (run>last_run)
        disp('Not available!!')
    end
end

% cycle (end at the last line of the script))
while (run)

    close all;
    clear R2_now;

    if (run)
        clear GEN F FITNESS R2 HITS EXPR;    
    end
%     %-----------------------------------------
%     % Figure 1 : plot F evolution of all runs
%     %-----------------------------------------
%     Fave_matrix=[];
%     leg={};         %cell(last_run,1);
%     s_leg=0;
%     ok_run=0;
%     filename = 'data_gp.txt';
%     for i_run=1:last_run; 
%         file = [directory_name '/run_' num2str(i_run) '/' filename]; 
%         fid = fopen(file);
%         if (fid~=-1)
%             C = textscan(fid,'%u %f %f %f %f %f %f %f %f %f %f %f %f %f %f','commentstyle','#'); % 15 columns
%             fclose(fid);
%             s_leg=s_leg+1;
%             ok_run=ok_run+1;
%             Fave_matrix(:,ok_run)=C{12};  % average F->C{13}, min F->C{12}
%             % if you launch HyGPvoew with an experiment with different no
%             % of generations per run it goes into error as dimensions are
%             % different!
%             leg(s_leg)={['Run ' num2str(i_run)]};    % leg(i_run)={['Run ' num2str(i_run)]}
%         end
%         
%     end
%     ok_run;
%     size(Fave_matrix);
%     plot(Fave_matrix,'-o')
%     legend(leg)
%     title({experiment; 'Evolution of minimum aggregate error (F) throughout the runs'},'FontSize',size_title,'Interpreter','none');

    % check the number of generations
    filename = 'best_gp.txt'; % best model - to access the individuals in the archive use 'objectives.txt';
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
        st = [ num2str(tot_gen) ' generations in the current run. Which? '];
    end

   % st = 'Which generation? ';
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
    % read tree objectives(TREE_MEAN, TREE_VAR, TREE_MIN, TREE_MAX, etc) at the requested generation
    disp([experiment ': tree in run ' num2str(run) ', generation ' num2str(cur_gen) ])
    %C = textscan(fid,'%u %f %f %f %d %q',1,'bufsize',16382,'commentstyle','#'); %,'endofline','\n')
    C = textscan(fid,'%u %f %f %u %f %f %f %f %f %f %f %f %q',1,'commentstyle','#'); %,'endofline','\n')  % best_gp.txt
    fclose(fid);
    cur_gen
    %# Gen F Fitness(RMSE) Corrections R2(adim.) Mean Var Min Max First_ACF Tot_variation Expression
    count_gen = C{1}
    F= C{2}
    FITNESS = C{3}
    CORRECTIONS = C{4}
    R2 = C{5}
    TREE_MEAN = C{6}
    TREE_VAR = C{7}
    TREE_MIN = C{8}
    TREE_MAX = C{9}
    TREE_FIRST_ACF = C{10}
    TREE_TOT_VARIATION = C{11}
    MAXABSERROR = C{12}
    EXPR = char(C{13})

    %
    tree_string = EXPR 
    RMSE_HyGP = FITNESS   
    R2_HyGP = R2  
    disp(' ')

    % here you can place the call to the function that computes the
    % objectives on build and test data sets
    
    
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

    %------------------------------------------------------------------------
    % ACTUAL VS ESTIMATED RESPONSES : RMSE, R2, min, max, etc are computed
    % in plot_act_vs_est script. Transfer plotting operations in this (main) file?
    %------------------------------------------------------------------------
    % plot actual vs estimated response on TRAINING data set
    [tree_output_build, objectives_build]=evaluate_individual(NFITCASES, NVAR, R, tree_string);
   
    %------------------------------------------------------------
    % plot estimated against actual response on TRAINING DATA SET
    %------------------------------------------------------------
    figure;
    plot(R(:,NVAR+1),tree_output,'or',R(:,NVAR+1),R(:,NVAR+1),'b')
    line1_title = experiment;
    line2_title =  ['Run ' num2str(run) ', Gen. ' num2str(cur_gen) ' - Training data set'];
    %line3_title = [dataset_type ' data set'];
    line4_title = ['Matlab RMSE = ' num2str(objectives_build.RMSE, '%e')];
    line5_title = ['Matlab R2 = ' num2str(objectives_build.R2, '%e')];
    %line6_title = ['Matlab Max error = ' num2str(error_max, '%e')];  % not yet in objectives_build
    %line7_title = ['Matlab Max rel error (%) = ' num2str(error_rel_max,'%e')]; % not yet in objectives_build
    line8_title = [' '];
    title({line1_title; line2_title; line4_title; line5_title;},'FontSize',size_title,'Interpreter','none');
    ylabel('Estimated','FontSize',size_title)
    xlabel('Actual', 'FontSize',size_title)
    axis equal;
    outbut_bounds=[min(R(:,NVAR+1)) max(R(:,NVAR+1))];
    if (size(outbut_bounds,1)==1 && size(outbut_bounds,2)==2)
        xlim(outbut_bounds);
        ylim(outbut_bounds);
    end
    
    
    % plot actual vs estimated response on TEST data set (if found)
    if (found_test) 
        [tree_output_test, objectives_test]=evaluate_individual(NTESTCASES, NVAR, Rtest, tree_string);
        
        %------------------------------------------------------------
        % plot estimated against actual response on TEST DATA SET
        %------------------------------------------------------------
        figure;
        plot(Rtest(:,NVAR+1),tree_output,'or',Rtest(:,NVAR+1),Rtest(:,NVAR+1),'b')
        line1_title = experiment;
        line2_title =  ['Run ' num2str(run) ', Gen. ' num2str(cur_gen) ' - Test data set'];
        %line3_title = [dataset_type ' data set'];
        line4_title = ['Matlab RMSE = ' num2str(objectives_test.RMSE, '%e')];
        line5_title = ['Matlab R2 = ' num2str(objectives_test.R2, '%e')];
        %line6_title = ['Matlab Max error = ' num2str(error_max, '%e')];  % not yet in objectives_build
        %line7_title = ['Matlab Max rel error (%) = ' num2str(error_rel_max,'%e')]; % not yet in objectives_build
        line8_title = [' '];
        title({line1_title; line2_title; line4_title; line5_title;},'FontSize',size_title,'Interpreter','none');
        ylabel('Estimated','FontSize',size_title)
        xlabel('Actual', 'FontSize',size_title)
        axis equal;
        outbut_bounds=[min(Rtest(:,NVAR+1)) max(Rtest(:,NVAR+1))];
        if (size(outbut_bounds,1)==1 && size(outbut_bounds,2)==2)
            xlim(outbut_bounds);
            ylim(outbut_bounds);
        end
    end
    disp(['Shown image: run ' num2str(run) ', generation ' num2str(cur_gen) ]); 

    
    %------------------------------------------------------------------------------------------
    % Figure 3 : plot objective function and selected HyGP model on building and test data sets
    % (RMSE and R2 are read from 'best_gp.txt' and currently only evaluated on building data set in plot_act_vs_est
    %------------------------------------------------------------------------------------------
    % plot model on building data set
    % RMSE_HyGP has to be equal to objectives_test.RMSE, R2_HyGP has to be equal to objectives_test.R2 
    % as same individual evaluated on the same data sets (training and test)
    plot_model(experiment, run, cur_gen, RMSE_HyGP, R2_HyGP, NVAR, R, tree_string, f_obj, 'Building');
    % plot model on test data set
    if (found_test) 
        % objectives.RMSE, objectives.R2 are the RMSE and the R2 of the current individual calculated by MATLAB
        % - not necessarily identical to RMSE_HyGP_test computed by HyGP.
        % This is because HyGP sorts the Population after evaluating the individuals on the test data set,
        % so the best individual on the train data set may be different from the best one on the test data set.
        plot_model(experiment, run, cur_gen, objectives_test.RMSE, objectives_test.R2, NVAR_test, Rtest, tree_string, f_obj, 'Test');
    end

    

    %------------------------------------------------------------------------------------------
    % Autocorrelation function of input data (building) set and selected HyGP model
    % (relies on the vector tree_output_build returned by plot_act_vs_est()
    %------------------------------------------------------------------------------------------
    if (NVAR==1)
        % Autocorrelation function plot
        % define lag: by default equal to half the number of NFITCASES 
        max_lag = floor(NFITCASES/2)-1;
        % compute autocorrelation function of original data
        [acf_original, ACF_pt_original] = compute_ACF(R(:,2), max_lag);
        % compute autocorrelation function of selected HyGP tree
        [acf_tree_build,ACF_pt_tree_build] = compute_ACF(tree_output_build, max_lag);
        delta=(ACF_pt_tree_build-floor(ACF_pt_tree_build)).*(R(ceil(ACF_pt_tree_build),1)-R(floor(ACF_pt_tree_build),1));
        objectives_build.ACF=R(floor(ACF_pt_tree_build),1)+delta;
        % plot the two autocorrelation functions
        figure;
        hold on
        plot(R(1:max_lag+1,1),acf_original,'-k',R(1:max_lag+1,1),acf_tree_build ,'-r')   
        % the line that follows can be improved.. too much discretised!
        plot(R(floor(ACF_pt_original),1),acf_original(floor(ACF_pt_original)),'dk',R(floor(ACF_pt_tree_build),1),acf_tree_build(floor(ACF_pt_tree_build)) ,'or')   
        hold off
        title('Autocorrelation function: original signal vs HyGP model')
        % captions
        line1_title = experiment;
        line2_title =  ['Run ' num2str(run) ', Gen. ' num2str(cur_gen) ' - Building data set'];
        line3_title = ['Autocorrelation function: original signal vs HyGP model'];
        title({line1_title; line2_title; line3_title},'FontSize',size_title,'Interpreter','none');
        legend('Original signal','HyGP model','Location','best')
    end
    
    %------------------------------------------------------------------------------------------
    % FFT plot of selected HyGP model
    % (relies on the vector tree_output_build returned by plot_act_vs_est()
    %------------------------------------------------------------------------------------------   
    if (NVAR==1)
        % FFT plot on training data set
        disp('Plot GP spectrum on training data set ...')
        HyGPfft(R, tree_output_build, experiment, run, cur_gen, 'Building');
        if (found_test) 
            % FFT plot on test data set
            disp('Plot GP spectrum on test data set ...')
            HyGPfft(Rtest, tree_output_test, experiment, run, cur_gen, 'Test');
        end
    end

    %--------------------------------------------------------------------------
    % Objectives on input data (building) set
    % Fitness(RMSE) Corrections R2(adim.) Mean Var Min Max First_ACF Tot_variation
    %--------------------------------------------------------------------------
    % read target data statistical properties on building data set from run_1/inputdata_stats.txt
    file = [directory_name sep '/run_' num2str(run) sep 'inputdata_stats.txt']  % sometimes run_1 fails with sge experiments.... do not know why
    fileID = fopen(file);
    C = textscan(fileID, '%f %*[^\n]','headerlines', 1); % textscan stops when it finds a different format, so at # First 5 record/fitness cases with target= 
    fclose(fileID);
    C_array=cell2mat(C);
    C_array(3,1) % Sum of target values: sum_output
    y_ave=C_array(4,1) % Average target value: y_ave 
    C_array(5,1) % SStot Total sum of squares of (target- y_ave) : Sy 
    y_var=C_array(6,1) % Target variance: y_var
    y_max=C_array(7,1) % Max value of target: y_max
    y_min=C_array(8,1) % Min value of target: y_min
    y_ACF=C_array(9,1) % First autocorrelation function root
    tot_variation_input=C_array(10,1) % Total variation: tot_variation_input 

    figure;
    marker_size=10;
    % Fitness(RMSE)
    subplot(1,9,1);
    hold on
    plot(0, FITNESS, 'o', 'MarkerSize',marker_size)
    plot(0, objectives_build.RMSE, 'xg', 'MarkerSize',marker_size)
    plot(0, objectives_test.RMSE, 'dr', 'MarkerSize',marker_size)
    hold off
    title({'RMSE'; '   '})
    % Corrections
    subplot(1,9,2);
    plot(0, CORRECTIONS, 'o', 'MarkerSize',marker_size)
    title({'Corrections'; '   '})
    % R2(adim.)
    subplot(1,9,3);
    hold on
    plot(0, R2, 'o', 'MarkerSize',marker_size)
    plot(0, objectives_build.R2, 'xg', 'MarkerSize',marker_size)
    plot(0, objectives_test.R2, 'dr', 'MarkerSize',marker_size)
    hold off
    title({'R2'; '   '})
    ylim([-1 1]);
    % Mean
    subplot(1,9,4);
    hold on
    plot(0, TREE_MEAN, 'o', 'MarkerSize',marker_size)
    plot(0, objectives_build.Mean, 'xg', 'MarkerSize',marker_size)
    plot(0, objectives_test.Mean, 'dr', 'MarkerSize',marker_size)
    plot(0, y_ave, '+k','MarkerSize',marker_size)
    hold off
    title({'Mean'; '   '})
    % Var
    subplot(1,9,5);
    hold on
    plot(0, TREE_VAR, 'o', 'MarkerSize',marker_size)
    plot(0, objectives_build.Var, 'xg', 'MarkerSize',marker_size)
    plot(0, objectives_test.Var, 'dr', 'MarkerSize',marker_size)
    plot(0, y_var, '+k','MarkerSize',marker_size)
    hold off
    title({'Var'; '   '})
    % Min
    subplot(1,9,6);
    hold on
    plot(0, TREE_MIN, 'o', 'MarkerSize',marker_size)
    plot(0, objectives_build.Min, 'xg', 'MarkerSize',marker_size)
    plot(0, objectives_test.Min, 'dr', 'MarkerSize',marker_size)
    plot(0, y_min, '+k','MarkerSize',marker_size)
    hold off
    title({'Min'; '   '})
    % Max
    subplot(1,9,7);
    hold on
    plot(0, TREE_MAX, 'o', 'MarkerSize',marker_size)
    plot(0, objectives_build.Max, 'xg', 'MarkerSize',marker_size)
    plot(0, objectives_test.Max, 'dr', 'MarkerSize',marker_size)
    plot(0, y_max, '+k','MarkerSize',marker_size)
    hold off
    title({'Max'; '   '})
    % First_ACF
    subplot(1,9,8);
    hold on
    plot(0, TREE_FIRST_ACF, 'o', 'MarkerSize',marker_size)
    plot(0,objectives_build.ACF, 'xg', 'MarkerSize',marker_size) 
    plot(0, y_ACF, '+k','MarkerSize',marker_size)
    hold off
    title({'ACF'; '   '})
    % Tot_variation
    subplot(1,9,9);
    hold on
    plot(0, TREE_TOT_VARIATION, 'o', 'MarkerSize',marker_size)
    plot(0, objectives_build.Tot_variation, 'xg', 'MarkerSize',marker_size)
    plot(0, objectives_test.Tot_variation, 'dr', 'MarkerSize',marker_size)
    plot(0, tot_variation_input, '+k','MarkerSize',marker_size)
    hold off
    legend('BUILD','BUILD (Matlab)','TEST (Matlab)','Target BUILD','Location','southoutside')
    title({'Tot variation'; '   '})
    % general title
    line1_title = experiment;
    line2_title =  ['Run ' num2str(run) ', Gen. ' num2str(cur_gen) ' - Objectives'];
    sgtitle({line1_title; line2_title},'FontSize',size_title,'Interpreter','none');

    
    % CYCLE: ASK WHICH RUN TO PLOT
    message = [ num2str(last_run) ' runs in the selected experiment. Which? (0 to exit) '];
    run=-1;
    while ((run==-1) || (run>last_run))
        run = input(message);
        if (run>last_run)
            disp('Not available!!')
        end
    end
end

clear all 

