%   script to represent graphically the objective function and
%   the one corresponding to the best tree of the experiment chosen
%   conventions: R matrix containing the training data set (i.e. input file data)
%                Rtest matrix containing the validation data set

clear all;
close all;
format long e;

% settings ---------------------------------------
plot_char = '-b';
size_title = 12;
pi=3.141592653589793e+00;
sep = filesep;  % set separator 
%-------------------------------------------------------

% get the name of the experiment and its directory (function)
%first_dir = ['~' sep 'code' sep 'C++' sep 'experiments' sep 'structuralGP_OLD' sep 'output' sep]
first_dir = ['E:' sep 'Research_Leeds' sep 'Linux' sep 'code' sep 'C++' sep 'experiments' sep 'structuralGP_OLD' sep 'output' sep]
[experiment, directory_name] = get_experiment(first_dir);
%if experiment selection is cancelled, directory_name should be zero
%and nothing should happen
if (directory_name == 0)
  disp('Quitted normally.')
  return
end



% READ INPUT FILE 
[PAR, R, CONSTR0, CONSTR1] = read_INPUT_FILE(directory_name);
NVAR = int16(PAR(2))
NFITCASES = int16(PAR(7));
SPLIT= int16(PAR(27));

% READ FILE CONTAINING TEST DATA SET
%oldFolder = cd(directory_name);
% select a file with TEST data
s = sprintf('\n');
disp(s);
disp('Select the file containing the TEST data set...');
first_path = '~/code/C++/experiments/structuralGP_OLD/output'
%first_path = 'U:\Linux\code\C++\experiments\structuralGP_OLD\output'; % AIRBUS
[file,path] = uigetfile('.txt','Select the file with the TEST data set', directory_name)
% read test data
% open stream 
found_test=0;
if (file~=0)
    fid_test = fopen(strcat(path,file),'r'); 

    if (fid_test~=-1)
        found_test = 1;
        % read header
        C = textscan(fid_test,'%d %d %*[^\n]',1,'commentstyle','#');
        NVAR_test=C{1};
        if (NVAR_test~=NVAR)
            disp('ERROR! NVAR different from NVAR_test:');
            disp('the number of variables in building and test data set are different. Exit');
            return
        end
    
        NTESTCASES=C{2};
        N = cast(NTESTCASES,'double');
        % read TEST data
        Rtest=zeros(NTESTCASES, NVAR+1);
        for k=1:NTESTCASES
            for j=1:NVAR
                C = textscan(fid_test,'%f',1,'commentstyle','#');
                Rtest(k,j)=C{1};
            end
            C = textscan(fid_test,'%f %*[^\n]',1,'commentstyle','#');
            Rtest(k,j+1)=C{1};
        end
        % close stream
        fclose(fid_test);
    end
else
    disp('No test data set selected.')
end


% SELECT OBJECTIVE FUNCTION
f_obj = assign_f_obj_fun();
f_obj


% split validating and building points (split DoE embedded in GP, not the external test data set)
% june 2017: as split has been changed to crossvalidation, check all the
% operations related to SPLIT
%if (SPLIT)
%    R_VAL = R(1:VALIDATING_LINES,:)
%    NVALCASES = size(R_VAL,1)
%    R_BUILD =R(VALIDATING_LINES+1:NFITCASES,:) 
%    NBUILDCASES = size(R_BUILD,1)
%end

%---------------------------------------------------------------------------
% cycle (end at the last line of the script))
another='0';
while (another~='n')

close all;
clear R2_now;

if (another~='0')
    clear GEN F FITNESS R2 HITS EXPR;    
end

% check the number of runs
last_run = count_runs(directory_name);
message = [ num2str(last_run) ' runs in the selected experiment. Which? '];
run=-1;
while ((run==-1) || (run>last_run))
    run = input(message);
    if (run>last_run)
        disp('Not available!!')
    end
end



% check the number of generations
filename = 'best_gp.txt'; % to access the individuals in the archive use 'objectives.txt';
file = [directory_name '/run_' num2str(run) '/' filename];
%
% extract number of generations from file (for my GP implementations)
[GEN] = textread(file,'%d %*[^\n]',-1,'bufsize',16382,'commentstyle','shell','endofline','\n');
dim_GEN = size(GEN);
tot_gen = dim_GEN(1) - 1;       %dim_GEN(1) is the number of rows in the text file...

% ask for the generation of the best individual
if (tot_gen==-1)
    disp('No generation in this run');
    break;
end
if (tot_gen==0)
    st = ['Only the initial generation (0) in the current run. 1 for single image, 0 for movie:'];
else 
   st = [ num2str(tot_gen) ' generations in the current run. 1 for single image, 0 for movie:'];
end

shot = input(st);
if (shot == 1)
    st = 'Which generation? ';
    gen(1)=-1;
    while ((gen(1)==-1) || (gen(1)>tot_gen))
        gen(1) = input(st);
        if (gen(1)>tot_gen)
            disp('Not available!!')
        end
    end
    gen(2) = gen(1);
end
if (shot == 0)
    disp('First and last frame (generation number)');
    gen(1) = input('start ');
    gen(2) = input('end ');
    FRAMES = moviein(gen(2)-gen(1)+1);
end


% generate NVAR variables with format Z1, Z2, ... , ZNVAR
NAMES = cell(NVAR,1);
for k=1:NVAR
    NAMES{k,1}=strcat('Z',num2str(k));
end
Z = genvarname(NAMES);
%

% build the matrix of independent variables
IND_VAR = R(:,1:NVAR);


% set limits to the plot
lim_var=zeros(NVAR, 2);
for k=1:NVAR
    lim_var(k,1) = min(IND_VAR(:,k));   %lower limit variable k   %RatPol2D .05   % Rosenbrock -2.
    lim_var(k,2) = max(IND_VAR(:,k)); %upper limit variable k     %RatPol2D 6.05   % Rosenbrock 2.
end

% plot upper limit
f_ulim = 1.2*max(R(:,NVAR+1)); %200;%4000;   %1.05  %1.2*max(max(f_obj));   
% plot lower limit
f_llim = min(R(:,NVAR+1)); %0;  %.95  %&1.2*min(min(f_obj));


%%%%%%%% 
% parameters for objective function visualization
% 1 INPUT VARIABLE 
if (NVAR==1)    
    Z1 = [lim_var(1,1) : abs(lim_var(1,2)- lim_var(1,1))/350. : lim_var(1,2)];
end
% 2 INPUT VARIABLES
if (NVAR==2)
    grid1 = [lim_var(1,1) : abs(lim_var(1,2)- lim_var(1,1))/50. : lim_var(1,2)];
    grid2 = [lim_var(2,1) : abs(lim_var(2,2)-lim_var(2,1))/50. : lim_var(2,2)];
    [Z1, Z2]  = meshgrid(grid1 , grid2);   % mesh in indipendent variables Z1, Z2
end
% % more than 2 INPUT VARIABLES
% if (NVAR>2) 
%     l_var = [];
%     % set fixed variables to define the "slice"
%     for k=1:NVAR-2
%         st = ['Number of the variable Zi. i= ' ];
%         i = input(st);
%         st = ['Value of variable Z' num2str(i) '= '];
%         val_i= input(st);
%         s = ['Z' num2str(i) '=' num2str(val_i)];
%         eval(s);
%         l_var = [l_var i]; 
%     end
%     % find out which variables are actually free
%     l_var = sort(l_var)
%     active_vars = [];
%     if (l_var(1)~=1) 
%         active_vars = [active_vars 1]; 
%     end
%     for k=1:size(l_var,2)-1
%         if (l_var(k+1)~=l_var(k)+1)
%             active_vars = [active_vars l_var(k)+1];
%         end
%     end
%     disp('Active variables:')
%     active_vars
%     % create the mesh
%     grid1 = [lim_var(active_vars(1),1) : abs(lim_var(active_vars(1),2)- lim_var(active_vars(1),1))/100. : lim_var(active_vars(1),2)];
%     grid2 = [lim_var(active_vars(2),1) : abs(lim_var(active_vars(2),2)-lim_var(active_vars(2),1))/100. : lim_var(active_vars(2),2)];
%     s=['[Z' num2str(active_vars(1)) ' , Z' num2str(active_vars(2)) ' ]  = meshgrid(grid1 , grid2);'];   % mesh in indipendent variables Z1, Z2
%     eval(s)
% end
% %%%%%%%% 


%plot objective function
if (f_obj)
    if (NVAR==1)    
        plot(Z1,eval(f_obj),'k');
    end
    if (NVAR==2)    
        mesh(Z1,Z2,eval(f_obj), 'EdgeColor','black');
    end
    if (NVAR>2)
        s = ['mesh(Z' num2str(active_vars(1)) ',Z' num2str(active_vars(2)) ' ,eval(f_obj));'];
        eval(s);
    end
end

%visualize random points in the domain (R2)
if (NVAR<3)
    figure(3);
    if (SPLIT)
        % validating points and building points 
        if (NVAR==1)
            plot(R_VAL(:,1), zeros(NVALCASES,1), 'or',R_BUILD(:,1), zeros(NBUILDCASES,1), 'xb'); 
        end
        if (NVAR==2) 
            plot(R_VAL(:,1),R_VAL(:,2), 'or',R_BUILD(:,1),R_BUILD(:,2), 'xb'); 
        end
        legend('Validating pts' , 'Building pts', 'Location', 'SouthOutside')
    else
        if (NVAR==1)
            plot(R(:,1),zeros(NFITCASES,1), 'or','MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize',4);
        end
        if (NVAR==2) 
            plot(R(:,1),R(:,2), 'or','MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize',4);
        end
        %legend('validating and building pts',  'Location', 'SouthOutside')
    end
    %axis tight equal;
    title({experiment; 'Selected points in the domain'},'FontSize',size_title,'Interpreter','none');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN MOVIE CYCLE
% reach the right generation on the file
fid = fopen(file);   
count_gen=0;
gen(1)
while (count_gen<gen(1)-1)
    % use text scan! text read starts every time from the beginning of
    % the file!!!
    C = textscan(fid,'%u %*[^\n]',1,'commentstyle','#');
    count_gen = C{1};
end

% set parameters
counter = 1;% for movie
if (NVAR<3)
    figure(4)
    set(gca, 'nextplot', 'replacechildren');
end

% starts cycle (MOVIE OR SINGLE SNAPSHOT)
for cur_gen = gen(1):gen(2)    % cur_gen is the generation counter
    
    if (~shot) 
        disp(' ')
        disp(' ')
        disp(['Generation ' num2str(cur_gen)])  
    end
    
    %---------------------------------
    % EXTRACT TREE EXPRESSION
    %---------------------------------
    disp([experiment ': tree in run ' num2str(run) ', generation ' num2str(cur_gen) ])
    C = textscan(fid,'%u %f %f %f %d %q',1,'bufsize',16382,'commentstyle','#'); %,'endofline','\n')
    cur_gen
    count_gen = C{1}
    F= C{2}
    FITNESS = C{3}
    R2 = C{4}
    HITS = C{5}
    EXPR = char(C{6})
    tree_string = EXPR; 
    fitness = FITNESS   
    R2_now = R2  
    disp(' ')
    disp('Tree expression before correction:')
    disp(tree_string)
    % correct tree expression to make it MATLAB readable
    tree = correct_equation(tree_string);
    disp(' ')
    disp('Tree expression after correction:')
    disp(tree)
    % print the latex expression of tree
    disp(' ')
    disp('Tree expression in LATEX:')
    tree_latex=print_latex_tree(NVAR,'Z',tree_string,'n');
    disp(tree_latex)
    %title('FontSize',size_title,'Interpreter','latex');
    
    
    
    % --------------------------------------------
    % PLOT OBJECTIVE FUNCTION and TREE EXPRESSION
    % --------------------------------------------
    % 1 INPUT VARIABLE 
    if (NVAR==1)    
        if (f_obj)
            % plot the objective function
            plot(Z1,eval(f_obj),'k');
        end
        hold on
        % plot the sample points
        plot(R(:,1), R(:,NVAR+1), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize',4);
        % visualize specified tree
        plot(Z1,eval(tree),'r'); %r
        hold off
        xlim([lim_var(1,1) lim_var(1,2)]);
        ylim([f_llim f_ulim]);
        % problem specific
        line([.3 .92],[1 1], 'LineStyle','--', 'Color', 'k')
        %line([.3 .92],[0 0], 'LineStyle','--',  'Color', 'k')
    end
    % 2 INPUT VARIABLES 
    if (NVAR==2)    
        if (f_obj)
            % plot the objective function
            mesh(Z1,Z2,eval(f_obj), 'EdgeColor','black');
        end
        hold on
        % plot the sample points
        plot3(R(:,1),R(:,2), R(:,NVAR+1), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize',4);
        %zlim([f_llim f_ulim]);
        % visualize specified tree
        %axes(r)
        mesh(Z1,Z2,eval(tree), 'EdgeColor','red');
    end 
    % more than 2 INPUT VARIABLES 
    if (NVAR>2)    
       
           
       
    end
    
    % captions
    line1_title = experiment;
    line2_title =  ['Run ' num2str(run) ', generation ' num2str(cur_gen) ];
    line3_title = ['RMSE = ' num2str(fitness, '%e')];
    line4_title = ['R2 = ' num2str(R2_now, '%e')];
    title({line1_title; line2_title; line3_title; line4_title},'FontSize',size_title,'Interpreter','none');
    f_llim
    f_ulim
    zlim([min(f_llim,f_ulim) max(f_llim,f_ulim)]);
    hold off
    
    if (~shot) 
        FRAMES(:,counter) = getframe(gcf);
        counter = counter+1;
    end
end
if (~shot) 
    h2 = figure;
    movie(h2, FRAMES, 1, 2)   %h2 is really important! Otherwise the axes are shown!
    movie2avi(FRAMES, ['movieHQ_G' num2str(gen(1)) '_' num2str(gen(2))], 'fps', 1, 'Quality', 100);   
end
%%%%%%%%%%%%%%%%% END MOVIE CYCLE
% close stream fid
fclose(fid);

%------------------------------------------------------------------------
% ACTUAL VS ESTIMATED RESPONSES
%------------------------------------------------------------------------
if (shot)
    % call the function to plot actual vs estimated response
    plot_act_vs_est(NFITCASES, NVAR, NAMES, R, tree, experiment, run, cur_gen, 'Building');
    if (found_test) 
        plot_act_vs_est(NTESTCASES, NVAR, NAMES, Rtest, tree, experiment, run, cur_gen, 'Test');
    end
    disp(['Shown image: run ' num2str(run) ', generation ' num2str(cur_gen) ]); 
end

% cycle
another=input('Another operation on the same experiment?(y/n)','s');

end

%cd(oldFolder);
close all;
clear all;  
