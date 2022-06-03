% function to print to file the expressions of the best individuals
% in the list in latex format 

%{
% get the name of the experiment and its directory (function)
%first_dir = '~/code/C++/experiments/structuralGP_OLD/output/'
[experiment, directory_name] = get_experiment('');
%if experiment selection is cancelled, directory_name should be zero
%and nothing should happen
if (directory_name == 0)
  disp('Quitted normally.')
  return
end
%}

%************************ FUNCTION ***********************
function print2file_best_latex(experiment_list, directory_list)
char(directory_list)

disp('print2file')
dim_experiment_list = size(experiment_list);
n_experiments = dim_experiment_list(1,2);

s='y';   % change according to your needs
if (s=='n')
    disp('Expressions will not be simplified')
end
if (s=='y')
    disp('Expressions will be simplified by MATLAB')
end

% open the file 
start_path = char(directory_list{1});
[filename,path] = uiputfile('*.txt','Save Latex expressions to file',start_path);
    
% open stream
fid_w=fopen([path '/' filename],'w');
if (fid_w==-1)
    disp('WRITE: Not able to open the file!')
end
% write header
fprintf(fid_w,'# - Experiment nameLatex expressions \n');
fprintf(fid_w,'# - - Latex expression of best individual of the experiment\n');
for k=1:n_experiments
    disp(['Best tree n. ' num2str(k) ' ...'])
    fprintf(fid_w,['# ' char(experiment_list{k}) '\n']);
    
    % open stream to read expression from file best_tree.txt -----------------
    path_best = char(directory_list{k});
    file_best = [ path_best '/best_tree.txt'];
    fid_r = fopen(file_best,'r');
    if (fid_r==-1)
        disp(['READ: Not able to open ' file_best]);
    end
    % read run no.
    C = textscan(fid_r,'%u %*[^\n]',1,'commentstyle','#');
    % read expression of best tree
    C = textscan(fid_r,'%*u %*f %*f %*f %*d %q %*[^\n]',1,'bufsize',16382,'commentstyle','#');
    % get char expression
    best_tree = char(C{1}); 
    % close stream to fid_r
    fclose(fid_r);
    disp('   read...')
    % --------------------------------------------------------------

    % open stream to read NVAR from file input_file.txt -----------------
    % read test data
    % open stream 
    path_input = char(directory_list{k});
    file_input = [ path_input '/input_file.txt']
    fid_NVAR = fopen(file_input,'r'); 
    if (fid_NVAR==-1)
        disp(['READ: Not able to open ' file_input]);
    end
    % read NVAR
    check = false;
    C = textscan(fid_NVAR,'%s %d',1,'delimiter','=','commentstyle','#');
    while ((check == false) & (~feof(fid_NVAR))) 
        C = textscan(fid_NVAR,'%s %d',1,'delimiter','=','commentstyle','#');
        check = strcmp(char(C{1}),'NVAR');
    end
    if (check == false)
        disp('field NVAR not found. Exit');
        break;
    end    
    
    NVAR=C{2};
    disp([' NVAR = ' num2str(NVAR)])
    % close stream to fid_NVAR
    fclose(fid_NVAR);
    %--------------------------------------------------------------------
    
    % get latex expression (y for simplify, n for as it is) 
    best_tree_latex = print_latex_tree(NVAR,'Z',best_tree,s);
    
    % print latex expression to file
    fprintf(fid_w,'%s\n', best_tree_latex); 
    disp('      ...saved')
end

% close stream to fid_w
fclose(fid_w);

disp('Done.')



