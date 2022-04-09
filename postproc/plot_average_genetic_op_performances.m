% script to plot on the same graph.
% - average size
% - average depth
% as a function of generations for two different experiments

function plot_average_genetic_op_performances(experiment_list, directory_list, selected_exp_no)

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
filename = 'averaged_adaptation_data.txt'; 
n_found = 0; 

% check if size_ave.txt exists in the experiment directory
n_exp = selected_exp_no;
disp([directory_list{n_exp}  sep  filename])
found = search_file(directory_list{n_exp},filename); 
experiment_name = experiment_list{n_exp};
directory_name = directory_list{n_exp};    
if (found)
    file = [directory_name sep  filename];
    % load the data (rows = generations, columns = runs)
    D = dlmread(file,' ',1,0);   % skip header....
    sD = size(D)
else
    m = [filename ' not found in ' directory_name];
    disp(m)
end

%n_found = n_found + found;
%if (n_found ~= n_experiments)
%disp('Not all files have been found')
%end

% redefinition
gen = D(:,1);
active_runs = D(:,2);
pop_size = D(:,3);
repr_rate = D(:,4); 
cross_rate = D(:,5); 
mut_rate = D(:,6); 
repr_tot = D(:,7);
reproduction_perf = D(:,8:10);
repr_av_delta = D(:,11:13);
cross_tot = D(:,14);
crossover_perf = D(:,15:17);
cross_av_delta = D(:,18:20);
tot_smut = D(:,21);
s_mutation_perf = D(:,22:24);
smut_av_delta = D(:,25:27);
tot_pmut = D(:,28);
p_mutation_perf = D(:,29:31);
pmut_av_delta = D(:,32:34);

%{  
% absolute parameters (not percentages -kept for checks)
nh = figure();
% REPRODUCTION
subplot(4,1,1)
hold on
plot(gen, reproduction_perf(:,1), 'k-.', gen, reproduction_perf(:,2), 'k:', gen, reproduction_perf(:,3), 'k--');
title('Reproduction performance','FontSize',size_title);
xlabel('Generation', 'FontSize',size_label);
ylabel('%', 'FontSize',size_label);

% CROSSOVER
subplot(4,1,2)
hold on
plot(gen,cross_tot,'k-',gen, crossover_perf(:,1), 'k-.', gen, crossover_perf(:,2), 'k:', gen, crossover_perf(:,3), 'k--');
title('Crossover performance','FontSize',size_title);
xlabel('Generation', 'FontSize',size_label);
ylabel('%', 'FontSize',size_label);

% SUBTREE MUTATION
subplot(4,1,3)
hold on
plot(gen,tot_smut,'k-',gen, s_mutation_perf(:,1), 'k-.', gen, s_mutation_perf(:,2), 'k:', gen, s_mutation_perf(:,3), 'k--');
title('Reproduction performance','FontSize',size_title);
xlabel('Generation', 'FontSize',size_label);
ylabel('%', 'FontSize',size_label);

% POINT MUTATION
subplot(4,1,4)
hold on
plot(gen,tot_pmut,'k-',gen, p_mutation_perf(:,1), 'k-.', gen, p_mutation_perf(:,2), 'k:', gen, p_mutation_perf(:,3), 'k--');
s_legend = cell([4,1]);
s_legend{1,1} = 'total';
s_legend{2,1} = 'destructive';
s_legend{3,1} = 'neutral';
s_legend{4,1} = 'constructive';
legend(s_legend, 'interpreter', 'none', 'FontSize', size_label)
title('Point mutation performance','FontSize',size_title);
xlabel('Generation', 'FontSize',size_label);
ylabel('%', 'FontSize',size_label);
%}

% percentages computation
repr_perc = zeros(sD(1),3);
cross_perc = zeros(sD(1),3);
s_mut_perc = zeros(sD(1),3);
p_mut_perc = zeros(sD(1),3);
for k=1:sD(1)
    repr_perc(k,:) = 100*reproduction_perf(k,:)./repr_tot(k);
    cross_perc(k,:) = 100*crossover_perf(k,:)./cross_tot(k);
    s_mut_perc(k,:) = 100*s_mutation_perf(k,:)./tot_smut(k);
    p_mut_perc(k,:) = 100*p_mutation_perf(k,:)./tot_pmut(k);
end
%plot!
nh = figure();
set(nh, 'Position', [10 10 700 800]) % width height])
% REPRODUCTION
subplot(4,1,1)
hold on
plot(gen, repr_perc(:,1), 'bd-.', gen, repr_perc(:,2), 'bo:', gen,repr_perc(:,3), 'b+-');
line1_title = [experiment_name ': average percentages'];
line2_title = 'Reproduction';
title({line1_title; line2_title},'FontSize',size_title,  'interpreter', 'none');
xlabel('Generation', 'FontSize',size_label);
ylabel('%', 'FontSize',size_label);
grid on

% CROSSOVER
subplot(4,1,2)
hold on
plot(gen, cross_perc(:,1), 'bd-.', gen, cross_perc(:,2), 'bo:', gen, cross_perc(:,3), 'b+-');
title('Crossover','FontSize',size_title);
xlabel('Generation', 'FontSize',size_label);
ylabel('%', 'FontSize',size_label);
grid on

% SUBTREE MUTATION
subplot(4,1,3)
hold on
plot(gen, s_mut_perc(:,1), 'bd-.', gen, s_mut_perc(:,2), 'bo:', gen, s_mut_perc(:,3), 'b+-');
title('Subtree mutation','FontSize',size_title);
xlabel('Generation', 'FontSize',size_label);
ylabel('%', 'FontSize',size_label);
grid on

% POINT MUTATION
subplot(4,1,4)
hold on
plot(gen, p_mut_perc(:,1), 'bd-.', gen, p_mut_perc(:,2), 'bo:', gen, p_mut_perc(:,3), 'b+-');
s_legend = cell([3,1]);
s_legend{1,1} = 'destructive';
s_legend{2,1} = 'neutral';
s_legend{3,1} = 'constructive';
legend(s_legend, 'interpreter', 'none', 'FontSize', size_label)
title('Point mutation','FontSize',size_title);
xlabel('Generation', 'FontSize',size_label);
ylabel('%', 'FontSize',size_label);
grid on

% plot average fitness increase/decrease for the genetic operators
nh = figure();
set(nh, 'Position', [10 10 700 800]) % width height])
% REPRODUCTION
subplot(4,1,1)
hold on
plot(gen, log10(abs(repr_av_delta(:,1))), 'bd-.', gen, log10(abs(repr_av_delta(:,2))), 'bo:', gen, log10(abs(repr_av_delta(:,3))), 'b+-');
line1_title = [experiment_name ': average error variations'];
line2_title = 'Reproduction';
title({line1_title; line2_title},'FontSize',size_title,  'interpreter', 'none');
xlabel('Generation', 'FontSize',size_label);
ylabel('log10(abs(dRMSE))', 'FontSize',10);
grid on
ylim([-20, +120]);
% CROSSOVER
subplot(4,1,2)
hold on
plot(gen, log10(abs(cross_av_delta(:,1))), 'bd-.', gen, log10(abs(cross_av_delta(:,2))), 'bo:', gen, log10(abs(cross_av_delta(:,3))), 'b+-');
title('Crossover','FontSize',size_title);
xlabel('Generation', 'FontSize',size_label);
ylabel('log10(abs(dRMSE))', 'FontSize',10);
grid on
ylim([-20, +120]);
% SUBTREE MUTATION
subplot(4,1,3)
hold on
plot(gen, log10(abs(smut_av_delta(:,1))), 'bd-.', gen, log10(abs(smut_av_delta(:,2))), 'bo:', gen, log10(abs(smut_av_delta(:,3))), 'b+-');
title('Subtree mutation','FontSize',size_title);
xlabel('Generation', 'FontSize',size_label);
ylabel('log10(abs(dRMSE))', 'FontSize',10);
grid on
ylim([-20, +120]);
% POINT MUTATION
subplot(4,1,4)
hold on
plot(gen, log10(abs(pmut_av_delta(:,1))), 'bd-.', gen,log10(abs(pmut_av_delta(:,2))), 'bo:', gen, log10(abs(pmut_av_delta(:,3))), 'b+-');
s_legend = cell([3,1]);
s_legend{1,1} = 'destructive';
s_legend{2,1} = 'neutral';
s_legend{3,1} = 'constructive';
legend(s_legend, 'interpreter', 'none', 'FontSize', size_label)
title('Point mutation','FontSize',size_title);
xlabel('Generation', 'FontSize',size_label);
ylabel('log10(abs(dRMSE))', 'FontSize',10);
grid on
ylim([-20, +120]);

% plot reproduction, crossover and mutation rates
figure;
hold on
plot(gen, repr_rate, 'b-', 'LineWidth',2)
plot(gen, cross_rate, 'r-.', 'LineWidth',2)
plot(gen, mut_rate, 'g--', 'LineWidth',2);
%s_legend = cell([3,1]);
s_legend{1,1} = 'reproduction';
s_legend{2,1} = 'crossover';
s_legend{3,1} = 'mutation';
legend(s_legend, 'interpreter', 'none', 'FontSize', size_label)
line1_title = [experiment_name ': genetic operators rates'];
title({line1_title},'FontSize',size_title,  'interpreter', 'none');
xlabel('Generation', 'FontSize',size_label);
ylabel('n/n_{tot}', 'FontSize',size_label);
ylim([0.0 1.0]);
grid on
