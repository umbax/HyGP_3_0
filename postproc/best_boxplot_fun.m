% function to plot the distribution of the best individual of each run
% at the end of an experiment

function best_boxplot_fun(experiment_list, directory_list, dataset)


dim_experiment_list = size(experiment_list);
n_experiments = dim_experiment_list(1,2);

% settings
sep = filesep; % set separator 
fig_handles = [];   %to keep track of figures' handles and close them at the end
s_fit = '';
plot_char = '-b';
s_F = 'F';
size_title = 12;
size_labels = 12;
LIM_VAL = 1000000000; %100.0; for test cases in ASMO-UK paper  %RMSE threshold between good and bad individual

% analysis on train data or test data?
% train data
if strcmp(dataset,'train')
    file_data = 'archives_best.txt';
    mark = 'training data';
end
% test data
if strcmp(dataset,'test')
    file_data = 'archives_best_TEST.txt';
    mark = 'test data';
end

% get the largest number of runs 
N_runs = zeros(n_experiments,1);
for k=1:n_experiments
    % retrieve number of runs per each experiment
    N_runs(k) = count_runs(directory_list{k});
end
N_runs
max_n_runs = max(N_runs)

% initialise the matrix that will contain the best individual's fitness 
% per each run (n columns = n experiments, each column is an experiment)
RMSE = zeros(max_n_runs, n_experiments);

experiment_list
% renaming of experiments: do not use it in normal conditions ------------
%experiment_list= {'  without penalisation  ', '  with penalisation (p=3) '}
%experiment_list{5} = 'omegalim_shift p=2 a5=0.0001'
%experiment_list{6} = 'omegalim_shift p=3 a5=0.01'
%experiment_list{7} = 'omegalim_shift p=3 a5=0.001'
%experiment_list{8} = 'omegalim_shift p=3 a5=0.0001'
%-------------------------------------------------------------------------

% you need to know Sy to compute the limit on R2 corresponding to LIM_VAL--
% set separator 
sep = filesep;
%first_path = ['~' sep 'code' sep 'C++' sep 'experiments'  sep 'structuralGP_OLD' sep 'output'];
%first_path = 'U:\Linux\code\C++\experiments\structuralGP_OLD\output'; % AIRBUS
first_path = 'D:\relocated_users\relocated_umba\Documents\OneDrive - University of Leeds\Research\2020_Research_ML_AI_Vassili\output';
[file,path] = uigetfile('.txt',['Select the file with the ' mark ' set'], first_path);

% read training data set (input file)
[PAR_value, R, CONSTR0, CONSTR1]=read_INPUT_FILE(path); 
NVAR = int16(PAR_value(2));
NCASES = int16(PAR_value(7));

% read test data set (not active yet)

N = cast(NCASES,'double');  % used?
% compute Sy
disp('Data matrix first line: ') 
R(1, :)
ave = ones(NCASES,1)*mean(R(:,NVAR+1));
Sy = sum((R(:,NVAR+1)-ave(:)).^2)
%--------------------------------------------------------------------------



% --------------------------
% RMSE boxplots
% --------------------------
size(RMSE)
for k=1:n_experiments
    % define the file 
    file = [directory_list{k} sep file_data]; 
    % check that file exists
    found = search_file(directory_list{k},file_data); 
    if (found==0)
        msg = ['File ' file ' NOT found! Launch test_eval or build archives_best_TEST first!']
        disp(msg)
        break
    end
    % load the data (rows = generations, columns = runs)
    file
    %fid = fopen(file,'r'); 
    RMSEsingle = textread(file,'%*d %*s %f %*f %*d %*[^\n]',-1,'bufsize', 50000, 'commentstyle','shell','endofline','\n')
    %C = textscan(fid,'%d %d %*[^\n]',1,'commentstyle','#');
    %fclose(fid);
    % initialize the RMSE matrix only the first time
    RMSEsingle
    RMSE(1:N_runs(k),k) = RMSEsingle;
    RMSE(N_runs(k)+1:max_n_runs,k) = nan; 
end

if (found==0)
    disp('Exit')
    return
end


RMSEgood=substitute_larger(RMSE,LIM_VAL)
% plot the distribution of the individuals of the archive of each run
% IMPORTANT! Boxplot does not consider NaN values. It neglects them. 
figure;
set(gca, 'FontSize', size_labels) % to change size of the labels' font
boxplot(RMSEgood,'label',experiment_list,'orientation', 'horizontal','notch','off')

xlabel('RMSE', 'FontSize', 12);
%max_x = min(max(max(RMSEgood)),35);
max_x = max(max(RMSEgood));
xlim([-.05 max_x]);
%xlim([-0.001 .5]);
title({['RMSE on ' mark];'Distribution of best individuals (one per run)'}, 'FontSize',size_title,'Interpreter','none');
h = gcf;
fig_handles = [fig_handles h];


% ---------------------
% R2 boxplots
% ---------------------
R2 = zeros(max_n_runs, n_experiments);
for k=1:n_experiments
    % define the file 
    file = [directory_list{k} '/' file_data]; 
    % read R2
    [r2single] = textread(file,'%*d %*s %*f %f %*d %*[^\n]',-1,'bufsize', 50000,'commentstyle','shell','endofline','\n');
    % load read R2 in DATA(row = run, column = experiment)
    R2(1:N_runs(k),k) = r2single;
    R2(N_runs(k)+1:max_n_runs,k) = nan; 
end
MIN_R2= 1-(N*LIM_VAL^2)/Sy
R2good=substitute_smaller(R2,MIN_R2)
% plot the distribution of the individuals of the archive of each run
% IMPORTANT! Boxplot does not consider NaN values. It neglects them. 
figure;
boxplot(R2good,'label',experiment_list,'orientation', 'horizontal','notch','off')
xlabel('R^2');
min_x = max(min(min(R2good)),-1.0);
xlim([min_x 1.05]);
%xlim([-2 1.05]);
title({['R^2 on ' mark];'Distribution of best individuals (one per run)'}, 'FontSize',size_title);
h = gcf;
fig_handles = [fig_handles h];
% change labels font size
%text_h=findobj(gca,'type','text')
%set(text_h,'FontSize',size_labels)


% find best individual on training or test data set ------------
disp('Best individual found per each experiment on the given dataset:'); 
experiment_list
[Max_value, corresp_run]=max(R2good(:,:))
[best, best_exp] = max(Max_value);
disp('Absolute best individual found in experiment: ');
experiment_list(best_exp)
disp('run: ');
corresp_run(best_exp)
disp('Score:')
best
%corresp_run(best_exp)
%experiment_list(best_exp)
%------------------------------------------------------------




%------------------------------------
% perform statistical analysis on RMSE: in the future put it in a function!
%------------------------------------
% select which parameter to be analysed
disp('Statistical analysis performed on RMSE')
DATA=RMSE;
%disp('Statistical analysis performed on R2')
%DATA=R2;

%--------------------------------------------------------------
N_tot = zeros(n_experiments,1);
N_undef = zeros(n_experiments,1);
P_undef = zeros(n_experiments,1);
N_bad = zeros(n_experiments,1);
P_bad = zeros(n_experiments,1);

% count undefined and bad individuals percentages
for k=1:n_experiments
    for j=1:N_runs(k)
        if (strcmpi(num2str(DATA(j,k)),'NaN'))
            % pathologic individual! (RMSE=NaN)
            N_undef(k)= N_undef(k)+1;
        else    
            if (DATA(j,k)>=LIM_VAL)
                % bad individual
                N_bad(k)= N_bad(k)+1;
            end
        end
        N_tot(k)=N_tot(k)+1;
    end
    P_undef(k)=100.0*N_undef(k)/N_tot(k);
    P_bad(k)=100.0*N_bad(k)/N_tot(k);
end

% print percentages (in the future use uitable)
s = sprintf('Experiment name   Total_no.  Perc undefined Perc RMSE> %f ', LIM_VAL);
disp(s);
for k=1:n_experiments
    %t = sprintf('%s\t%d\t%d\t%.1f\t%d\t%.1f', char(experiment_list(k)), N_tot(k), N_undef(k), P_undef(k), N_bad(k), P_bad(k));
    t= sprintf('%s\t\t%d\t%.1f%%\t%.1f%%', char(experiment_list(k)), N_tot(k), P_undef(k), P_bad(k));
    disp(t);
end


% descriptive statistics of the samples purged of bad and undefined -------
% individuals
med=zeros(1,n_experiments);
mu=zeros(1,n_experiments);
sigma2=zeros(1,n_experiments);
iqrange=zeros(1,n_experiments);
% substitute RMSE values >LIM_VAL with NaN:
% the statistical analysis will be done on the data set purged of outliers
DATA = substitute_larger(DATA, LIM_VAL);  % valid if DATA = RMSE
for k=1:n_experiments
    x=remove_nan(DATA(:,k)); 
    med(1,k) = median(x);
    mu(1,k) = mean(x);  % mean of the best individuals per run per each exp.
    iqrange(1,k)=iqr(x);
    sigma2(1,k) = var(x); % variance of the best individuals per run per each exp.
end

dat = cell(n_experiments, 5);
% set data fot table (see sketch)
for k=1:n_experiments
    dat{k,1} = experiment_list{1,k};            % name
    dat{k,2} = num2str(N_tot(k), '%d');         % total number of individuals
    dat{k,3} = num2str(P_undef(k), '%.0f');    % perc undefined
    dat{k,4} = num2str(P_bad(k), '%.0f');      % perc bad
    dat{k,5} = num2str(med(1,k), '%.3e');         % median
    dat{k,6} = num2str(iqrange(1,k), '%.3e');     % interquartile range
    %dat{k,3} = num2str(mu(1,k), '%e');         % mean
    %dat{k,4} = num2str(sigma2(1,k),'%e');      % variance
end

cnames = {'Experiment','No. tot','P_undef','P_bad','Median (rest)','IQR (rest)'};
cformat = {'char','numeric','numeric', 'numeric','numeric','numeric'}; 
width = 800;
height = 300;
cwidth = {310 45 60 60 110 110};
f = figure('Name', ['RMSE statistics on ' mark],'Position',[100 100 width height]);
% uitable accepted only by recent MATLAB releases?
uitable('Units','Pixels', 'Data',dat,'ColumnName', cnames, 'ColumnFormat', cformat, 'ColumnWidth', cwidth,'Parent',f, 'Position',[1 1 width height]);
fig_handles = [fig_handles f];
%--------------------------------------------------------------



% stop here if experiment_list is made of only a single experiment
if (n_experiments==1)
    return;
end

% NON PARAMETRIC TEST ---------------------------------
% perform Kruskal-Wallis test - (Hp: same distribution, only translated) STILL TO BE TESTED!!!
% H0: all samples have the same median
% H1 : the median is different
alpha = 0.05;
[p_kw,table,stats] = kruskalwallis(DATA, experiment_list); 
h = gcf;
fig_handles = [fig_handles h];

[c,m] = multcompare(stats, 'alpha', alpha, 'estimate', 'kruskalwallis');
title('')
xlabel('')
h = gcf;
fig_handles = [fig_handles h];

disp(['Result of Kruskal-Wallis test on ' mark ': p-value = ' num2str(p_kw)]);
disp(['significance level: ' num2str(alpha)]);

if (p_kw<alpha)
    disp('The null hypothesis (H0) that all the samples come from the same population is rejected');
    disp('Wilcoxon rank-sum test (called also Mann-Whitney U test) will be performed')
    disp('---------------------------')
    p_wil=(-1)*ones(n_experiments);
    for k=1:n_experiments
       for j=k+1:n_experiments
            % perform single Wilcoxon rank-sum test between to samples
            disp(['Comparing ' char(experiment_list(k)) ' with ' char(experiment_list(j))]);
            % IMPORTANT! REMOVE NaN from the vectors...
            x=remove_nan(DATA(:,k));
            y=remove_nan(DATA(:,j));
            [p,h] = ranksum(x,y,'alpha',alpha);
            p_wil(j,k)=p;
            disp(['p_value = ' num2str(p) ' h = ' num2str(h)]);
            if (h==1) 
                disp('NOT EQUAL MEDIAN!'); 
            end
       end
    end
    % print table to ease the reading
    dat = cell(n_experiments, n_experiments+1);
    % set data fot table (see sketch)
    clear cnames;
    clear cformat;
    clear cwidth;
    cnames = cell(1,n_experiments+1);
    cnames{1} = 'Experiments';
    cformat = cell(1,n_experiments+1);
    cformat{1} = 'char';
    cwidth = cell(1,n_experiments+1);
    cwidth{1} = 310;
    for k=1:n_experiments
        dat{k,1} = experiment_list{1,k};            % name
        for j=1:n_experiments
            dat{k,j+1} = num2str(p_wil(k,j), '%.2f');         % p_value
        end
        cnames{1,k+1} = num2str(k);
        cformat{1,k+1} = 'numeric';
        cwidth{1,k+1} = 110;
    end
    width = 800;
    height = 300;
    f = figure('Name', ['Wilcoxon rank-sum test on ' mark ': p_values'],'Position',[100 100 width height]);
    % uitable accepted only by recent MATLAB releases?
    uitable('Units','Pixels', 'Data',dat,'ColumnName', cnames, 'ColumnFormat', cformat, 'ColumnWidth', cwidth,'Parent',f, 'Position',[1 1 width height]);
    %'Units','Pixels',
    fig_handles = [fig_handles f];
else
   disp('Not enough evidence to reject the null hypothesis :')
   disp('all the samples come from the same population. End of the test.')
end
% --------------------------------------------------------

% PARAMETRIC TEST ---------------------------------
% perform ANOVA for mean analysis
% PROBLEM: the observations should come from a normal distribution..
s = sprintf('\n\n');
disp(s);
disp('Parametric test : ONE-WAY ANOVA')
% H0: all samples have the same mean
% H1 : the mean is different
[p_value,table,stats] = anova1(DATA, experiment_list); 
h = gcf;
fig_handles = [fig_handles h];

[c,m] = multcompare(stats, 'alpha', alpha, 'estimate', 'anova1');
h = gcf;
fig_handles = [fig_handles h];

title('')
xlabel('')
disp(['Result of ANOVA test (F statistics): p-value = ' num2str(p_value)]);
disp(['significance level: ' num2str(alpha)]);
if (p_value<alpha)
    disp('The null hypothesis (H0) that all the samples come from the same population is rejected');
    disp('ANOVA tests between single pairs of experiments will be performed')
    disp('---------------------------')
    p_an=(-1)*ones(n_experiments);
    for k=1:n_experiments
       for j=k+1:n_experiments
            % perform single ANOVA test between two samples
            disp(['Comparing ' char(experiment_list(k)) ' with ' char(experiment_list(j))]);
            % IMPORTANT! anova1 seems to neglect NaN...
            x=DATA(:,k);
            y=DATA(:,j);
            w=[x y];
            group{1,1}=char(experiment_list{1,k});
            group{2,1}=char(experiment_list{1,j});
            p = anova1(w,group,'off');
            p_an(j,k)=p;
            disp(['p_value = ' num2str(p)]);
            if (p_value<alpha)
                disp('MEANS DIFFERENT!');
            end
       end
    end
    
    
    % print table to ease the reading
    clear dat;
    dat = cell(n_experiments, n_experiments+1);
    % set data fot table (see sketch)
    clear cnames;
    clear cformat;
    clear cwidth;
    cnames = cell(1,n_experiments+1);
    cnames{1} = 'Experiments';
    cformat = cell(1,n_experiments+1);
    cformat{1} = 'char';
    cwidth = cell(1,n_experiments+1);
    cwidth{1} = 310;
    for k=1:n_experiments
        dat{k,1} = experiment_list{1,k};            % name
        for j=1:n_experiments
            dat{k,j+1} = num2str(p_an(k,j), '%.2f');         % p_value
        end
        cnames{1,k+1} = num2str(k);
        cformat{1,k+1} = 'numeric';
        cwidth{1,k+1} = 110;
    end
    width = 800;
    height = 300;
    f = figure('Name', ['One-way Anova test on ' mark ': p_values'],'Position',[100 100 width height]);
    % uitable accepted only by recent MATLAB releases?
    uitable('Units','Pixels', 'Data',dat,'ColumnName', cnames, 'ColumnFormat', cformat, 'ColumnWidth', cwidth,'Parent',f, 'Position',[1 1 width height]);
    %'Units','Pixels',
    fig_handles = [fig_handles f];
else
   disp('Not enough evidence to reject the null hypothesis :')
   disp('all the samples come from the same population (same mean). End of the test.')
end
% --------------------------------------------------------
s = sprintf('\n\nPress q + return to close all figures and exit\n');
q='';
while (strcmp(q,'q')==0)
    q = input(s,'s');
end
for k=1:size(fig_handles,2)
    disp(fig_handles(k))
    close(fig_handles(k))
end    
clear all;


