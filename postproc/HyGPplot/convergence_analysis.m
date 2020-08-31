% script to plot the points returned by HyGP model and the target points for all HyGP runs of an experiment.
% providing this way a sort of convergence analysis or sensitivity of the evolution on the initial seed.
% The script does not evaluate the tree returned by HyGP, simply plots the points in the two columns
% in file "points_gp.txt" written by HyGP in each run folder 

% get path of selected experiment and its name
sep = filesep;  % set separator 
init_path="/home/umba/Documents/Eclipse/HyGP_3_0/output/200720_Elnaz/Mode1";
exp_path = uigetdir(init_path, "Selet the experiment folder");
split_point = rindex(exp_path, sep)+1;
exp_name = substr(exp_path, split_point)

% define storage
tree_matrix = [];

% cycle through the runs and plot the corresponding output in points_gp.txt 
run_name = "run_1";
current_run = 1;
next_run=1;
while exist ([exp_path sep run_name], "dir")
  current_run=next_run;
  legend_names{current_run}=run_name;
  
  file_path = [exp_path sep run_name sep "points_gp.txt"];
  [n, target, tree, residuals] = textread([file_path], "%n %f %f %f", "headerlines", 12);
  %size(n)
  %size(tree)
  %size(target)
  tree_matrix(:,current_run)=tree;
  
  next_run = current_run+1;
  run_name = ["run_" int2str(next_run)];
  
endwhile

figure(1);
hold on;
plot(n,tree_matrix,'linewidth', 3) %'-x;HyGP model;'
legend_names{current_run+1}="target";
plot(n,target, '-og;target;', 'linewidth', 3)
grid on
legend_names
h = legend(legend_names);
set(h,'Interpreter', 'none')


%msg = ["Found " int2str(current_run) " runs. Execution terminated, no other runs found"];
%msgbox (msg, "End of script")

