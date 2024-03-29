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


%%% evaluate_individual.m
%%% MATLAB function to evaluate the individual on the training or test data set
%%% it also plots the actual vs estimated response on the given data set

function [tree_output, objectives]=evaluate_individual(NFITCASES, NVAR, R, tree_string)
    
    disp('start plot_act_vs_est:')
    %-----------------------------------------------
    % data analysis on given data set R (may be building or test data set)
    %-----------------------------------------------
    
    % settings ---------------------------------------
    plot_char = '-b';
    size_title = 11;
    pi=3.141592653589793e+00;
    
    tree_output = zeros(NFITCASES,1);
    sum_tree_output=0.0;
    sum_squared_tree_output=0.0;
    tot_variation=0.0;
    errorabs=0.;
    errorsq = 0.;
    error = 0.0;
    error_max = 0.0;
    error_rel = 0.0;
    error_rel_max = 0.0;
    norm_errorsq=0.0;
    norm_RMSE = 0.0;
    
    % error computation
    for k=1:NFITCASES  % rows of the data matrix (MATLAB comincia a contare gli array da 1!)
        x=[];
        % assign the value of independent variables for the k fitness case
        for j=1:NVAR   % columns of the data matrix
            % collect data is a single vector to pass to "compute_expression"
            x=[x R(k,j)];   % Example: Z3 = R(k,3)  for the k-th fitness case (record)
        end

        % compute tree value
        tree_output(k) = compute_expression(tree_string, x);
        %x
        %tree_output(k)
        
        sum_tree_output=sum_tree_output+tree_output(k);
        sum_squared_tree_output=sum_squared_tree_output+tree_output(k)*tree_output(k);
        if (k>1) 
            tot_variation=tot_variation+abs(tree_output(k)-tree_output(k-1)); 
        end
        error = tree_output(k) - R(k,NVAR+1);
        error_rel = 100.0*(error/R(k,NVAR+1));
        
        % update error as sum of the absolute differences
        errorabs = errorabs + abs(tree_output(k) - R(k,NVAR+1));
        % update error as root mean square error
        errorsq = errorsq + (tree_output(k) - R(k,NVAR+1))^2;
        % update sum of the squares of the relative errors
        if (R(k,NVAR+1)==0)
            norm_errorsq = norm_errorsq + ((tree_output(k) - R(k,NVAR+1)))^2;
        else
            norm_errorsq = norm_errorsq + ((tree_output(k) - R(k,NVAR+1))/R(k,NVAR+1))^2;
        end
        
        % update absolute error
        if (abs(error)> abs(error_max))
            error_max = error;    
        end    
        if  (abs(error_rel) > abs(error_rel_max))
            error_rel_max = error_rel;
        end
    end   % end NFITCASES cycle
    
    % percentual variation from actual to estimated
%     disp(' ')
%     disp('Quality indicators:')
%     disp('type comp to compare actual and estimated output')
%     disp('For the table with estimated and actual output type comp')
%     disp('      Actual (Obj. function)         Estimated (Tree)         Variation%      ')
%     comp = [R(:,NVAR+1) tree_output(:) 100.0*(tree_output(:)-R(:,NVAR+1))./R(:,NVAR+1)];
%     
%     max_perc_var = max(abs(comp(:,3)));
%     disp(['Maximum absolute % variation = ' num2str(max_perc_var)])
    
    % HyGP tree objectives on the given data set (build or test)
    %--------------------------------------------------------------------
    N = cast(NFITCASES,'double');
    % 1) RMSE
    rmse = sqrt(errorsq/N);
    %rmse = sqrt(mean((comp(:,1)-comp(:,2)).^2));
    objectives.RMSE=rmse;    
    % 2) R2  
    %Serr = sum(((comp(:,1)-comp(:,2)).^2));
    Serr = sum(( (R(:,NVAR+1) - tree_output(:)) .^2));
    y_ave = mean(R(:,NVAR+1))*ones(N,1);
    Sy =  sum((R(:,NVAR+1)-y_ave).^2);
    objectives.R2=1.- Serr/Sy;
    % 3) Mean
    objectives.Mean=sum_tree_output/N;
    % 4) Variance (see reformulation in http://datagenetics.com/blog/november22017/index.html)
    objectives.Var=sum_squared_tree_output/N - objectives.Mean*objectives.Mean;
    % 5) Min
    objectives.Min=min(tree_output);
    % 6) Max
    objectives.Max=max(tree_output);
    % 7) ACF
    % 18/11/21 to be completed!
    % 8) Total variation
    objectives.Tot_variation=tot_variation;
    
    % norm_RMSE = sqrt(norm_errorsq/N);
    % error_max; % to be added to objectives instead of hits? Is it L-infinite norm? 
    % error_rel_max;
    % errorabs;  
   
%disp('end plot_act_vs_est:')
return
end

