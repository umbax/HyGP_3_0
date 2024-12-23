% function that plots the evolution of the average of the best 
% (minimum) RMSE (fitness) on training data set


function [r] = plot_average_training_F(experiment_list, directory_list)

disp('Function plot_average_training_F called...')

% settings
plot_char = '-b';
size_title = 12;
size_label = 12;

% variables
sep = filesep; % set separator 
found = zeros(1,15);
MAX_VALUE = 1.e20;
MIN_VALUE = 1.e-20;

dim_experiment_list = size(experiment_list);
n_experiments = dim_experiment_list(1,2);

% set directory and file names
filename = 'F_min.txt';   %RMSE error
n_found = 0; 
S_cell = cell([n_experiments 1]);
Smean_cell = cell([n_experiments 1]);
gen_cell = cell([n_experiments 1]);
%directory_name = directory_list{1};
%experiment_name = experiment_list{1};
%found = search_file(directory_name,filename);

% load the data (rows = generations, columns = runs)
for j=1:n_experiments
    
    % check if fitness_min.txt exists in the experiment directory
    found = 0;
    found = search_file(directory_list{j},filename)
    
    if (found)
        file = [directory_list{j} sep  filename];
        % load the data (rows = generations, columns = runs)
        delimiter = ' ';
        Fit = dlmread(file)
        dim_Fit=size(Fit);
        [Fitmean, nrun_valid] = smart_mean(Fit)
        s = size(Fitmean)
        gen = [0:s(1)-1]
        % update cells
        S_cell{j,1} = Fit;
        Smean_cell{j,1} = Fitmean;
        gen_cell{j,1} = gen;
    else
        m = [filename ' not found in ' directory_name];
        disp(m)
    end
    n_found = n_found + found;
    
end

if (n_found ~= n_experiments)
    disp('Not all files have been found')
end

% plot the minium F found (in the run archive) at each generation for each run
% ATTENTION: REALLY MESSY FOR MORE THAN ONE EXPERIMENT!!
nh = figure();
set(nh,'DefaultAxesLineStyleOrder','-|-.|--|:', 'DefaultAxesColorOrder',[0 0 1])
% this way up to 4 different lien styles can be created: you'll have to find
% a way to incorporate different markers as well...
% {'+','o','*','.','x','s','d','^','v','>','<','p','h'} 
%plot(gen',Fitmean,plot_char)
for j=1:n_experiments
    plot(gen_cell{j},S_cell{j})
    hold all
end
title('Minimum F for each run (archives)','FontSize',size_title);
xlabel('Generation', 'FontSize',size_label);
ylabel('F','FontSize',size_label);
legend(experiment_list, 'interpreter', 'none', 'FontSize', size_label,'Location','best')


% plot the average of the minium F found at each generation throughout all runs
nh = figure();
set(nh,'DefaultAxesLineStyleOrder','-|-.|--|:', 'DefaultAxesColorOrder',[0 0 1])
% this way up to 4 different lien styles can be created: you'll have to find
% a way to incorporate different markers as well...
% {'+','o','*','.','x','s','d','^','v','>','<','p','h'} 
%plot(gen',Fitmean,plot_char)
for j=1:n_experiments
    plot(gen_cell{j},Smean_cell{j})
    hold all
end
title('Average of minimum F through archives','FontSize',size_title);
xlabel('Generation', 'FontSize',size_label);
ylabel('F','FontSize',size_label);
legend(experiment_list, 'interpreter', 'none', 'FontSize', size_label)

r = 1;
return