% script to plot the points returned by HyGP model and the target points.
% The script does not evaluate the tree returned by HyGP, simply plots the points in the two columns
% in file "points_gp.txt" written by HyGP in each run folder 

sep = filesep;  % set separator 
dir_name="/home/umba/Documents/Eclipse/HyGP_3_0/output/200720_Elnaz/Mode1";
[file,path] = uigetfile('.txt','Select the file points_gp.txt', dir_name)
[n, target, tree, residuals] = textread ([path sep file], "%n %f %f %f", "headerlines", 12);
plot(n,target, '-og;target;', 'linewidth', 3, tree, '-xr;HyGP model;', 'linewidth', 3)
grid on
legend