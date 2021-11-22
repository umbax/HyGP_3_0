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


%%% read_INPUT_FILE.m
%%% script to read GP input file
%%% input: name and path of the file
%%% output: cell containing parameters names and their values PAR, record matrix R   

%{
%**************TO BE REMOVED**********************
% get the name of the experiment and its directory (function)
first_dir = '~/code/C++/experiments/structuralGP_OLD/output/'
[experiment, directory_name] = get_experiment(first_dir);
%if experiment selection is cancelled, directory_name should be zero
%and nothing should happen
if (directory_name == 0)
  disp('Quitted normally.')
  return
end
%************************************
%}

function [PAR_value, R, CONSTR0, CONSTR1]=read_INPUT_FILE(directory_name)


% directories definition
exp_dir = [directory_name '/'];
filename = 'input_file.txt'
file = [exp_dir filename]
% open file for reading
fid = fopen(file);   
if (fid==-1) 
    disp('input_file.txt not found. Exit.')
    return
end

% list of labels (cell - to convert to char use char())
PAR_name  ={'SEED',
            'NVAR',
            'MINRAND',
            'MAXRAND',
            'MAX_N_PERIODS',
            'PSO_NPARTICLES',
            'PSO_NITERATIONS',
            'NFITCASES',
            'METHOD',
            'DEPTH_MAX',
            'DEPTH_MIN',
            'DEPTH_LIM',
            'p_FULL',
            'REPR_RATE',
            'CROSS_RATE',
            'MUT_RATE',
            'COMP_RATE',
            'NEW_RATE',
            'M',
            'G',
            'BOUNDED',
            'STRAT_STATP',
            'W_STRAT_STATP',
            'W_ACF',
            'W_TVARIATION',
            'W_COMPLEXITY',
            'W_N_CORRECTIONS',
            'W_SIZE',
            'W_FACTORISATION',
            'N_GUESSES',
            'CROSSVALIDATION',
            'FOLDS_N',
            'THRESHOLD',
            'N_INEQUALITY0',
            'W_PEN_ORD0',
            'N_INEQUALITY1',
            'W_PEN_ORD1'}; 
        
n_PAR = numel(PAR_name);
PAR_value = zeros(n_PAR,1);
k_PAR=1;


   
label = 'ciao';
while ((~strcmp(label,'BINARY_FUN=')) & (k_PAR<n_PAR+1))
        
        disp(['Reading ', char(PAR_name(k_PAR)) ,' ...']);
        
        % = textscan(fid,  %u %u %u %f %f %f %f %u %u %u %f %u',1,'commentstyle','#');
        C = textscan(fid,'%s %f',1,'delimiter','=','commentstyle','#');
        
        % check if the first entry is equal to the name of the 
        %disp(char(C{1}));
        %disp(num2str(C{2}));
        %disp(PAR_name(k_PAR));
        if (strcmp(char(C{1}),PAR_name(k_PAR)))
            PAR_value(k_PAR)=C{2};
            k_PAR=k_PAR+1;
        end
       
end

% define explicitly NVAR and NFTICASES so to avoid using PAR_value(...) in the rest of the script        
NVAR = int16(PAR_value(2));
NFITCASES = int16(PAR_value(8));


%read binary operations
bin_op = '';
disp(['Reading binary operations ...']);
C = textscan(fid,'%s %s',1,'delimiter','=','commentstyle','#');
if (strcmp(char(C{1}),'BINARY_FUN'))
   bin_op = char(C{2});
   k_PAR=k_PAR+1;
end

% read unary operations
un_op = '';
disp(['Reading unary operations ...']);
C = textscan(fid,'%s %s',1,'delimiter','=','commentstyle','#');
%label = fgets(fid);
if (strcmp(char(C{1}),'UNARY_FUN'))
    un_op = char(C{2});
end

% declare record matrix
R = zeros(NFITCASES, NVAR+1);

% read data
format_data = '%f';
for k=2:NVAR+1
    format_data = strcat(format_data, ' %f');
end
format_data = strcat(format_data, ' %*[^\n]');
k=1;

disp(['Reading sample matrix (R)...']);
while ((~feof(fid)) && (k<NFITCASES+1))
    C = textscan(fid, format_data, 1,'commentstyle','#');
    for j=1:NVAR+1
        R(k,j) = C{j};
    end
    k=k+1;
end

% close stream fid
fclose(fid);

CONSTR0 = zeros(1,1);
CONSTR1 = zeros(1,1);

% print data read
for k=1:n_PAR
    disp([char(PAR_name(k)) ' = ' num2str(PAR_value(k))]);
end
disp('');
disp(['Binary operations:' bin_op]);
disp(['Unary operations:' un_op(1:size(un_op,2))]);
disp(' ');

%R;
%CONSTR0;
%CONSTR1;