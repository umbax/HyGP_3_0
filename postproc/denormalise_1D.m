%%% script to denormalise the output of a HyGP model evovled using normalised data
%%% It currently works only for NVAR=1 !!!!


clear all;
close all;
format long e;

% settings and variables ----------------------------------------
plot_char = '-b';
size_title = 11;
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
[NVAR, NFITCASES, R, CONSTR0, CONSTR1]=read_INPUT_FILE(directory_name); % ERROR is NVAR>1!
%NVAR = int16(PAR(2))
%NFITCASES = int16(PAR(8))

% READ FILE CONTAINING TEST DATA SET
%[NVAR_test, NTESTCASES, Rtest, found_test]=read_TEST_FILE(directory_name, NVAR);
[NVAR_test, NTESTCASES, Rtest, found_test]=read_TEST_FILE(directory_name);



% normalising coefficients for independent data
a_ind = input("Coefficient a for independent data: ")
b_ind = input("Coefficient b for independent data: ")

% normalising coefficients for dependent data
a_dep = input("Coefficient a for dependent data: ")
b_dep = input("Coefficient b for dependent data: ")

% cycle through all the runs
run=1;
    
    % check the number of generations
    filename = 'best_gp.txt'; % best model - to access the individuals in the archive use 'objectives.txt';
    file = [directory_name '/run_' num2str(run) '/' filename];
    % extract number of generations from file (for my GP implementations)
    [GEN] = textread(file,'%d %*[^\n]',-1,'bufsize',16382,'commentstyle','shell','endofline','\n');
    tot_gen = size(GEN,1) - 1;       %size(GEN,1) is the number of rows in the text file...
    
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
    [R_ORIG, output_build_ORIG, objectives_build_ORIG]=compute_error(a_ind, a_dep, b_ind, b_dep, NFITCASES, NVAR, R, tree_string, 'Building');
    % compute error on test data set (if found): data set is denormalised using a_ind and b_ind and then RMSE is computed
    if (found_test) 
        [Rtest_ORIG, output_test_ORIG, objectives_test_ORIG]=compute_error(a_ind, a_dep, b_ind, b_dep, NTESTCASES, NVAR, Rtest, tree_string, 'Test');
    end

    % write to file points_GP_ORIGINAL.txt
    output_build_ORIG
    objectives_build_ORIG.RMSE
    % R(:,NVAR+1) tree_output(:)

% end cycle

% write to file archives_BEST_ORIGINAL.txt

% write to file archives_BEST_TEST_ORIGINAL.txt




function [R_ORIG, tree_output_ORIG, objectives]=compute_error(a_ind, a_dep, b_ind, b_dep, NFITCASES, NVAR, R, tree_string, dataset_type)
    
    disp('compute_error():')
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
    
    % denormalise R  (NVAR must be 1, NVAR>1 not accepted so far)
    % independent variable
    R_ORIG(:,1)=(R(:,1)-b_ind.*ones(NFITCASES,1))./a_ind;
    % dependent variable
    R_ORIG(:,2)=(R(:,2)-b_dep.*ones(NFITCASES,1))./a_dep;

    % error computation (MIND that the script assumes NVAR=1)
    for k=1:NFITCASES  % rows of the data matrix (MATLAB comincia a contare gli array da 1!)
        %x=[];
        % assign the value of normalised independent variables for the k-th fitness case
        x=R(k,1);    %x=[x R(k,1)];   % Example: Z3 = R(k,3)  for the k-th fitness case (record)
        
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
        if (R_ORIG(k,2)==0)
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
    
    % percentual variation from actual to estimated
    disp(' ')
    disp([dataset_type ' data set quality indicators:'])
    disp('type comp to compare actual and estimated output')
    disp('For the table with estimated and actual output type comp')
    disp('      Actual (Obj. function)         Estimated (Tree)         Variation%      ')
    comp = [R_ORIG(:,2) tree_output_ORIG(:) 100.0*(tree_output_ORIG(:)-R_ORIG(:,2))./R_ORIG(:,2)];
    
    max_perc_var_ORIG = max(abs(comp(:,3)));
    disp(['Maximum absolute % variation = ' num2str(max_perc_var_ORIG)])
    
    % output values of HyGP tree on the given data set (build or test)
    %--------------------------------------------------------------------
    N = cast(NFITCASES,'double');
    % 1) RMSE
    %rmse_ORIG = sqrt(errorsq_ORIG/N);
    rmse_ORIG = sqrt(mean((comp(:,1)-comp(:,2)).^2));
    objectives.RMSE=rmse_ORIG;    
    % 3) R2  
    Serr = sum(((comp(:,1)-comp(:,2)).^2));
    y_ave = mean(R(:,2))*ones(N,1);
    Sy =  sum((R(:,2)-y_ave).^2);
    objectives.R2=1.- Serr/Sy;
    % 4) Mean
    objectives.Mean=sum_tree_output_ORIG/N;
    % 5) Variance (see reformulation in http://datagenetics.com/blog/november22017/index.html)
    objectives.Var=sum_squared_tree_output_ORIG/N - objectives.Mean*objectives.Mean;
    % 6) Min
    objectives.Min=min(tree_output_ORIG);
    % 7) Max
    objectives.Max=max(tree_output_ORIG);
    % 8) ACF
    % 18/11/21 to be completed!
    % 9) Total variation
    objectives.Tot_variation=tot_variation_ORIG;


end
