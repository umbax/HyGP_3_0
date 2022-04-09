% script to plot on the same graph.
% - average size
% - average depth
% as a function of generations for two different experiments

function plot_size_dynamics(experiment_list, directory_list)

%disp('Function plot_size_dynamics called...')

% variables
sep = filesep; % set separator 
found = zeros(1,15);
MAX_VALUE = 1.e20;
MIN_VALUE = 1.e-20;

% settings
%s_fit = 'RMSE';
%plot_char = '-b';
s_F = 'F';
size_title = 12;
size_label = 12;
h = []; % figure handles

dim_experiment_list = size(experiment_list);
n_experiments = dim_experiment_list(1,2);

% verify if all files size_ave.txt are present
filename = 'size_ave.txt'; 
n_found = 0; 
Smean_cell = cell([n_experiments 1]);
gen_cell = cell([n_experiments 1]);
for j=1:n_experiments
    
    % check if size_ave.txt exists in the experiment directory
    found = 0;
    found = search_file(directory_list{j},filename); 
    
    if (found)
        file = [directory_list{j} sep  filename];
        % load the data (rows = generations, columns = runs)
        S = dlmread(file);
        size(S);
        [Smean, k] = smart_mean(S);
        s = size(Smean)
        gen = [0:s(1)-1];
        % update cells
        Smean_cell{j,1} = Smean;
        gen_cell{j,1} = gen;
    else
        m = [filename ' not found in ' directory_list{j}];
        disp(m)
    end
    n_found = n_found + found;
    
end

if (n_found ~= n_experiments)
    disp('Not all files have been found')
end

nh = figure();
% 1 COLOUR, 4 DIFFERENT LINE STYLES 
%set(nh,'DefaultAxesLineStyleOrder','-|-.|--|:', 'DefaultAxesColorOrder',[0 0 1])
% 2 COLOURS, 4 DIFFERENT LINE STYLES 
%set(nh,'DefaultAxesLineStyleOrder','-|-.|--|:', 'DefaultAxesColorOrder',[0 0 1; 0 1 0])
% 1 COLOUR, 8 DIFFERENT LINE STYLES 
set(nh,'DefaultAxesLineStyleOrder','-|-.|--|:|-x|-d|-*|-o', 'DefaultAxesColorOrder',[0 0 1])
% {'+','o','*','.','x','s','d','^','v','>','<','p','h'} 
for j=1:n_experiments
    plot(gen_cell{j},Smean_cell{j})
    hold all
end
title('Average size of archive individuals','FontSize',size_title);
xlabel('Generation', 'FontSize', size_label);
ylabel('Size (nodes)', 'FontSize', size_label);
legend(experiment_list, 'interpreter', 'none', 'FontSize', size_label)




% check cell content
%celldisp(Smean_cell)
%celldisp(gen_cell)  
