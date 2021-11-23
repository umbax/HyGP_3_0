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


%%% plot_model.m
%%% function that plots the model

function [outputArg1,outputArg2] = plot_model(experiment, run, cur_gen, fitness, R2_now, NVAR, R, tree_string, f_obj, dataset_type)

size_title = 11;

% define variables
% generate NVAR variables with format Z1, Z2, ... , ZNVAR
NAMES = cell(NVAR,1);
for k=1:NVAR
    NAMES{k,1}=strcat('Z',num2str(k));
end
Z = genvarname(NAMES);
% build the matrix of independent variables
IND_VAR = R(:,1:NVAR);

% set limits to the design space and create mesh for plotting
lim_var=zeros(NVAR, 2);
for k=1:NVAR
    lim_var(k,1) = min(IND_VAR(:,k));   %lower limit variable k   %RatPol2D .05   % Rosenbrock -2.
    lim_var(k,2) = max(IND_VAR(:,k)); %upper limit variable k     %RatPol2D 6.05   % Rosenbrock 2.
end

% upper limit for plots
f_ulim = 1.2*max(R(:,NVAR+1)); %200;%4000;   %1.05  %1.2*max(max(f_obj));   
% lower limit for plots
f_llim = min(R(:,NVAR+1)); %0;  %.95  %&1.2*min(min(f_obj));

% define parameters for objective function visualization
% 1 INPUT VARIABLE 
if (NVAR==1)    
    npoints_1D=size(R,1);     %6000 %1000.0;
    Z1 = [lim_var(1,1) : abs(lim_var(1,2)- lim_var(1,1))/(npoints_1D-1.0) : lim_var(1,2)]
    % compute HyGP model values in Z1
    tree_values=zeros(1,size(Z1,2));
    tree_values_size = size(tree_values,2)
    for i=1:tree_values_size
        tree_values(i)= compute_expression(tree_string, Z1(i));  % not possible to use as input the whole vector x
    end
end
% 2 INPUT VARIABLES
if (NVAR==2)
    npoints_2D=51.0; %number of levels for full factorial (number of points) equal on both axes
    z1_grid1 = [lim_var(1,1) : abs(lim_var(1,2)- lim_var(1,1))/(npoints_2D-1.0) : lim_var(1,2)];
    z2_grid2 = [lim_var(2,1) : abs(lim_var(2,2)-lim_var(2,1))/(npoints_2D-1.0) : lim_var(2,2)];
    [Z1, Z2]  = meshgrid(z1_grid1 , z2_grid2);   % full-factorial mesh in indipendent z1_grid1 and z2_grid2 : Z1 and Z2 are matrices!
    % compute HyGP model values in Z1, Z2  ... how to use meshgrid?
    tree_values=zeros(size(z2_grid2,2), size(z1_grid1,2));
    for col=1:size(z1_grid1,2)
        for row=1:size(z2_grid2,2)
            % compute tree value in position (row, col)
            x(1)=Z1(row,col);
            x(2)=Z2(row,col);
            tree_values(row, col)=compute_expression(tree_string, x); 
        end % end row cycle 
    end % end col cycle
end
% more than 2 INPUT VARIABLES
if (NVAR>2) 
    l_var = [];
    % set fixed variables to define the "slice"
    for k=1:NVAR-2
        st = ['Number of the variable Zi. i= ' ];
        i = input(st);
        st = ['Value of variable Z' num2str(i) '= '];
        val_i= input(st);
        s = ['Z' num2str(i) '=' num2str(val_i)];
        eval(s);
        l_var = [l_var i]; 
    end
    % find out which variables are actually free
    l_var = sort(l_var)
    active_vars = [];
    if (l_var(1)~=1) 
        active_vars = [active_vars 1]; 
    end
    for k=1:size(l_var,2)-1
        if (l_var(k+1)~=l_var(k)+1)
            active_vars = [active_vars l_var(k)+1];
        end
    end
    disp('Active variables:')
    active_vars
    % create the mesh
    grid1 = [lim_var(active_vars(1),1) : abs(lim_var(active_vars(1),2)- lim_var(active_vars(1),1))/50. : lim_var(active_vars(1),2)];
    grid2 = [lim_var(active_vars(2),1) : abs(lim_var(active_vars(2),2)-lim_var(active_vars(2),1))/50. : lim_var(active_vars(2),2)];
    s=['[Z' num2str(active_vars(1)) ' , Z' num2str(active_vars(2)) ' ]  = meshgrid(grid1 , grid2)];'];   % mesh in indipendent variables Z1, Z2
    eval(s) % to create the mesh
    % compute HyGP model values in Zi, Zj
end
%&%&%

%%%
%------------------------------------
% Figure 1 : plot objective function
%------------------------------------
if (f_obj)
    figure(1);
    if (NVAR==1)    
        plot(Z1,eval(f_obj),'k');
        xlabel('Z1')
    end
    if (NVAR==2)    
        mesh(Z1,Z2,eval(f_obj), 'EdgeColor','black');
        xlabel('Z1')
        ylabel('Z2')
        view(-38,30)
    end
end
%%%

% tricks for plotting correctly
if (NVAR<3)
    figure
    set(gca, 'nextplot', 'replacechildren');
end
% 1 INPUT VARIABLE 
if (NVAR==1)    
    if (f_obj)
        % plot the objective function
        plot(Z1,eval(f_obj),'k');
    end
    hold on
    % plot the sample points
    plot(R(:,1), R(:,NVAR+1), '-ok', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize',2);
    % visualize specified tree
      %tree_string_c = correct_equation(tree_string);  %% if there are problems in computing tree_values
      %plot(Z1,eval(tree_string_c),'b','LineWidth', 1);    
    plot(Z1,tree_values,'r','LineWidth', 1); 
    plot ([11 11], [min(tree_values) max(tree_values)], 'b', 'Linewidth', 4)  %&%&% line beyond which there is extrapolation
    hold off
    xlim([lim_var(1,1) lim_var(1,2)]);
    %ylim([f_llim f_ulim]);
    % problem specific
    %line([.3 .92],[0 0], 'LineStyle','--',  'Color', 'k')
    xlabel('Z1')
end
% 2 INPUT VARIABLES 
if (NVAR==2)    
    if (f_obj)
        mesh(Z1,Z2,eval(f_obj), 'EdgeColor','black');
    end
    hold on
    % plot the sample points
    plot3(R(:,1),R(:,2), R(:,NVAR+1), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize',4);
    %zlim([f_llim f_ulim]);
    % visualize specified tree
    %axes(r)
    % plot the objective function
    mesh(Z1,Z2,tree_values, 'EdgeColor','r');
    % Az -38 El 30
    view(-38,30)
    xlabel('Z1')
    ylabel('Z2')
end 
% more than 2 INPUT VARIABLES 
if (NVAR>2)    
    % space left empty on purpose
end

legend('Original signal','HyGP model','Location','best')    
% captions
line1_title = experiment;
line2_title =  ['Run ' num2str(run) ', Gen. ' num2str(cur_gen) ' - ' dataset_type ' data set'];
line3_title = ['HyGP RMSE = ' num2str(fitness, '%e')];
line4_title = ['HyGP R2 = ' num2str(R2_now, '%e')];
title({line1_title; line2_title; line3_title; line4_title},'FontSize',size_title,'Interpreter','none');
%zlim([min(f_llim,f_ulim) max(f_llim,f_ulim)]);  % sometimes when uncommented the plot is blank... don't know why...
hold off


end

