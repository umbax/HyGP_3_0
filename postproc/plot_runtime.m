% function to plot the distribution of the runtime required by each run

function plot_runtime(experiment_list, directory_list)

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
filename = 'time_stat.txt'; 
n_found = 0; 
time_collector = cell([1 n_experiments]);
for j=1:n_experiments
    
    % check if size_ave.txt exists in the experiment directory
    found = 0;
    found = search_file(directory_list{j},filename); 
    
    if (found)
        file = [directory_list{j} sep  filename];
        % load the data (rows = generations, columns = runs)
        S = dlmread(file);
        %dim_S = size(S)
        % update time_matrix
        time_collector{1,j} = S;
    else
        m = [ filename ' not found in ' directory_list{j}];
        disp(m)
    end
    n_found = n_found + found;
    
end

if (n_found ~= n_experiments)
    m = ['Not all files ' filename ' files have been found']
    disp(m)
end

nh = figure();
set(nh,'DefaultAxesLineStyleOrder','-|-.|--|:', 'DefaultAxesColorOrder',[0 0 1])
% this way up to 4 different lien styles can be created: you'll have to find
% a way to incorporate different markers as well...
% {'+','o','*','.','x','s','d','^','v','>','<','p','h'} 
set(gca, 'FontSize', size_label) % to change size of the labels' font
boxplot(cell2mat(time_collector),'label',experiment_list,'orientation', 'horizontal','notch','off')
xlabel('time per GP run (s)', 'FontSize', size_label);




% check cell content
%celldisp(Smean_cell)
%celldisp(gen_cell)  
