%%% script to denormalise the output of a HyGP model evovled using normalised data
%%% It currently works only for NVAR=1 !!!!

%%% input: normalised tranining data set (input file) and test data set (test file)
%%% output: peformance of the best tree for the run and generation specified by the user

clear all;
close all;
format long e;

% settings and variables ----------------------------------------
plot_char = '-b';
size_title = 13;
pi=3.141592653589793e+00;
sep = filesep;  % set separator
objectives_build_ORIG=struct('RMSE',0.0,'Corrections',0,'R2',0.0,'Mean',0.0,'Var',0.0,'Min',0.0,'Max',0.0,'ACF',0.0,'Tot_variation',0.0);
objectives_test_ORIG=struct('RMSE',0.0,'Corrections',0,'R2',0.0,'Mean',0.0,'Var',0.0,'Min',0.0,'Max',0.0,'ACF',0.0,'Tot_variation',0.0);
objectives_target_ORIG=struct('RMSE',0.0,'Corrections',0,'R2',0.0,'Mean',0.0,'Var',0.0,'Min',0.0,'Max',0.0,'ACF',0.0,'Tot_variation',0.0);
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
[NVAR, NFITCASES, R, CONSTR0, CONSTR1]=read_INPUT_FILE(directory_name); % ERROR if NVAR>1!
%NVAR = int16(PAR(2))
%NFITCASES = int16(PAR(8))


% READ FILE CONTAINING TEST DATA SET
%[NVAR_test, NTESTCASES, Rtest, found_test]=read_TEST_FILE(directory_name, NVAR);
[NVAR_test, NTESTCASES, Rtest, found_test]=read_TEST_FILE(directory_name);



% % normalising coefficients for independent data
% a_ind = input("Coefficient a for independent data (x): ")
% b_ind = input("Coefficient b for independent data (y): ")
a_ind = 977.826513
b_ind = 0.989990

% % normalising coefficients for dependent data
% a_dep = input("Coefficient a for dependent data (x): ")
% b_dep = input("Coefficient b for dependent data (y): ")
a_dep = 183.833620 
b_dep = 6.003371

% ASK WHICH RUN TO PLOT
message = [ num2str(last_run) ' runs in the selected experiment. Which? '];
run=-1;
while ((run==-1) || (run>last_run))
    run = input(message);
    if (run>last_run)
        disp('Not available!!')
    end
end

% cycle through all the runs
while (run)

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
    
    tree_string = EXPR 
    RMSE_HyGP = FITNESS   
    R2_HyGP = R2  

    % compute error on training data set: data set is denormalised using a_ind and b_ind and then RMSE is computed
    [R_ORIG, output_build_ORIG, objectives_build_ORIG]=evaluate_individual_to_ORIG_1D(a_ind, a_dep, b_ind, b_dep, NFITCASES, NVAR, R, tree_string, 'Building');
    % compute error on test data set (if found): data set is denormalised using a_ind and b_ind and then RMSE is computed
    if (found_test) 
        [Rtest_ORIG, output_test_ORIG, objectives_test_ORIG]=evaluate_individual_to_ORIG_1D(a_ind, a_dep, b_ind, b_dep, NTESTCASES, NVAR, Rtest, tree_string, 'Test');
    end


    %------------------------------------------------------------------------------------------
    % Figure 3 : plot target function and selected HyGP model on test data set
    %------------------------------------------------------------------------------------------
    hold on
    % plot the sampling/target points
    plot(Rtest_ORIG(:,1), Rtest_ORIG(:,2), '-ok', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize',2);
    % visualize selected tree
    plot(Rtest_ORIG(:,1),output_test_ORIG,'r','LineWidth', 1); 
    plot ([(11-b_ind)/a_ind (11-b_ind)/a_ind], [min(output_test_ORIG) max(output_test_ORIG)], 'b', 'Linewidth', 4)  %&%&% line beyond which there is extrapolation
    legend('Original signal','HyGP model','Location','best')    
    % captions
    line1_title = experiment;
    line2_title =  ['Run ' num2str(run) ', Gen. ' num2str(cur_gen) ' - test data set'];
    line3_title = ['RMSE = ' num2str(objectives_test_ORIG.RMSE, '%e')];
    line4_title = ['R2 = ' num2str(objectives_test_ORIG.R2, '%e')];
    title({line1_title; line2_title; line3_title; line4_title},'FontSize',size_title,'Interpreter','none');
    xlim([Rtest_ORIG(1,1) Rtest_ORIG(NTESTCASES,1)]);
    hold off

    %------------------------------------------------------------------------------------------
    % FFT plot of selected HyGP model
    % (relies on the vector tree_output_build returned by plot_act_vs_est()
    %------------------------------------------------------------------------------------------   
    if (NVAR==1)
        % FFT plot on training data set
        disp('Plot GP model spectrum on training data set ...')
        HyGPfft(R_ORIG, output_build_ORIG, experiment, run, cur_gen, 'Build');
        if (found_test) 
            % FFT plot on test data set
            disp('Plot GP model spectrum on test data set ...')
            HyGPfft(Rtest_ORIG, output_test_ORIG, experiment, run, cur_gen, 'Test');
        end
    end

    %------------------------------------------------------------------------------------------
    % write data to file (only referring to test data set for now 6/7/23)
    %------------------------------------------------------------------------------------------
    format long;  
    % write to file test data
    testdata2file=[Rtest_ORIG(:,1), Rtest_ORIG(:,2), output_test_ORIG, output_test_ORIG-Rtest_ORIG(:,2)]
    file_ORIG_name = [directory_name sep 'Denormalised_GPdata_test_r' num2str(run) 'g' num2str(cur_gen) '.txt']
    fileID = fopen(file_ORIG_name,'w');
    fprintf(fileID, '# %s Run=%s Gen=%s\r\n',experiment, num2str(run), num2str(cur_gen));
    fprintf(fileID, '# No. of test cases:\r\n');
    fprintf(fileID, '%d \r\n', size(testdata2file,1));
    fprintf(fileID, '# Objectives on test data set: RMSE  Corrections  R2  Mean  Var  Min  Max  ACF  Tot_variation\r\n'); 
    fprintf(fileID, '%f  %d  %f  ', objectives_test_ORIG.RMSE, objectives_test_ORIG.Corrections, objectives_test_ORIG.R2);
    fprintf(fileID, '%f  %f  %f  %f  ',objectives_test_ORIG.Mean, objectives_test_ORIG.Var, objectives_test_ORIG.Min, objectives_test_ORIG.Max);
    fprintf(fileID, '%f  %f \r\n', objectives_test_ORIG.ACF, objectives_test_ORIG.Tot_variation); 
    fprintf(fileID, '# Output on test data set: x  y_target  y_GPmodel  residual (y_GPmodel-y_target)\r\n');
    fprintf(fileID, '%e  %e  %e  %e \r\n', testdata2file');  % important! Transpose the matrix, otherwise it does not work as intended! 
    fclose(fileID);
    %writematrix(testdata2file,[directory_name sep 'Denormalised_GP_data_test_check.txt'])
    
    % CYCLE: ASK WHICH RUN TO PLOT
    message = [ num2str(last_run) ' runs in the selected experiment. Which? (0 to exit) '];
    run=-1;
    while ((run==-1) || (run>last_run))
        run = input(message);
        if (run>last_run)
            disp('Not available!!')
        end
    end

% end cycle
end





function [R_ORIG, tree_output_ORIG, objectives_ORIG]=evaluate_individual_to_ORIG_1D(a_ind, a_dep, b_ind, b_dep, NFITCASES, NVAR, R, tree_string, dataset_type)
    
    %%% ATTENTION: funciton valid ONLY for NVAR=1!!!
    %%% _ORIG : variable trasnformed into ORIGINAL (physical) space (basically, denormalised) 
    if ((NVAR>1) || (NVAR<1)) 
        disp('This function only works for NVAR = 1. Exit')
        exit
    end

    disp('evaluate_individual_to_ORIG_1D():')
    %-----------------------------------------------
    % data analysis on given data set R (may be building or test data set)
    %-----------------------------------------------
    
    % settings ---------------------------------------
    %plot_char = '-b';
    %size_title = 11;
    pi=3.141592653589793e+00;
    
    tree_output_ORIG = zeros(NFITCASES,1);
    sum_tree_output_ORIG=0.0;
    sum_squared_tree_output_ORIG=0.0;
    tot_variation_ORIG=0.0;
    errorabs_ORIG=0.;
    errorsq_ORIG = 0.;
    error_ORIG = 0.0;
    error_max_ORIG = 0.0;
    error_rel_ORIG = 0.0;
    error_rel_max_ORIG = 0.0;
    norm_errorsq_ORIG=0.;
    norm_RMSE_ORIG = 0.;
    Mean_ORIG = 0.;
    
    % denormalise R  (NVAR must be 1, NVAR>1 not accepted so far)
    % independent variable
    R_ORIG(:,1)=(R(:,1)-b_ind.*ones(NFITCASES,1))./a_ind
    % dependent variable (target)
    R_ORIG(:,2)=(R(:,2)-b_dep.*ones(NFITCASES,1))./a_dep;

    % error computation (MIND that the script assumes NVAR=1)
    for k=1:NFITCASES  % rows of the data matrix (MATLAB comincia a contare gli array da 1!)
        %x=[];
        % assign the value of normalised independent variables for the k-th fitness case
        x=R(k,1);    
        
        % compute tree value in normalised space and denormalise it
        %tree_output(k) = compute_expression(tree_string, x) % normalised version
        tree_output_ORIG(k) = (compute_expression(tree_string, x)-b_dep)/a_dep;
                
        sum_tree_output_ORIG=sum_tree_output_ORIG+tree_output_ORIG(k);
        sum_squared_tree_output_ORIG=sum_squared_tree_output_ORIG+tree_output_ORIG(k)*tree_output_ORIG(k);
        if (k>1) 
            tot_variation_ORIG=tot_variation_ORIG+abs(tree_output_ORIG(k)-tree_output_ORIG(k-1)); 
        end
        error_ORIG = tree_output_ORIG(k) - R_ORIG(k,2);
        error_rel_ORIG = 100.0*(error_ORIG/R_ORIG(k,2));
        
        % update error as sum of the absolute differences
        errorabs_ORIG = errorabs_ORIG + abs(tree_output_ORIG(k) - R_ORIG(k,2));
        % update error as root mean square error
        errorsq_ORIG = errorsq_ORIG + (tree_output_ORIG(k) - R_ORIG(k,2))^2;
        % update sum of the squares of the relative errors
        if (abs(R_ORIG(k,2))< 1e-12)
            norm_errorsq_ORIG = norm_errorsq_ORIG + ((tree_output_ORIG(k) - R_ORIG(k,2)))^2;
        else
            norm_errorsq_ORIG = norm_errorsq_ORIG + ((tree_output_ORIG(k) - R_ORIG(k,2))/R_ORIG(k,2))^2;
        end
        
        % update absolute error
        if (abs(error_ORIG)> abs(error_max_ORIG))
            error_max_ORIG = error_ORIG;    
        end    
        if  (abs(error_rel_ORIG) > abs(error_rel_max_ORIG))
            error_rel_max_ORIG = error_rel_ORIG;
        end
    end   % end NFITCASES cycle
       
    % output values of HyGP tree on the given data set (build or test)
    %--------------------------------------------------------------------
    N = cast(NFITCASES,'double');
    % 1) RMSE
    %rmse_ORIG = sqrt(errorsq_ORIG/N);
    rmse_ORIG = sqrt(mean((R_ORIG(:,2) - tree_output_ORIG(:)).^2));
    objectives_ORIG.RMSE=rmse_ORIG;    
    % 2) Corrections
    objectives_ORIG.Corrections=-1;    
    % 2) R2  
    Serr_ORIG = sum(((R_ORIG(:,2) - tree_output_ORIG(:)).^2));
    y_ave_ORIG = mean(R_ORIG(:,2))*ones(N,1);
    Sy_ORIG =  sum((R_ORIG(:,2)-y_ave_ORIG).^2);
    objectives_ORIG.R2 = 1.- Serr_ORIG/Sy_ORIG;   % R2 is R2_ORIG, obviously....
    % 3) Mean
    Mean_ORIG=sum_tree_output_ORIG/N;
    objectives_ORIG.Mean=Mean_ORIG;  % Mean is Mean_ORIG, obviously....
    % 4) Variance (see reformulation in http://datagenetics.com/blog/november22017/index.html)
    objectives_ORIG.Var=sum_squared_tree_output_ORIG/N - Mean_ORIG*Mean_ORIG;
    % 5) Min
    objectives_ORIG.Min=min(tree_output_ORIG);
    % 6) Max
    objectives_ORIG.Max=max(tree_output_ORIG);
    % 7) ACF
    objectives_ORIG.ACF = NaN;    % 18/11/21 to be completed!
    % 8) Total variation
    objectives_ORIG.Tot_variation=tot_variation_ORIG;

end
