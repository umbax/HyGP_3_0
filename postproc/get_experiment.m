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


%%% get_experiment.m
%%% FUNCTION TO RETRIEVE THE NAME OF AN EXPERIMENT
%%% AND ITS DIRECTORY THROUGH A SMALL GUI
%%% input: last directory visited (in order to spare time when comparing similar experiments)

function [experiment, directory_name] = get_experiment(last_directory)

% set separator according to operating system
% ispc isunix
sep = filesep;
% go to the directory containing the experiment
dim_dir = size(last_directory);
% extract the name of the experiment
pos = strfind(last_directory,sep);  %find the last occurrence of /
last = max(pos)-1;
last_address = last_directory(1:last);


% get the directory
if (strcmp(last_directory,''))
    address = ['~' sep 'code' sep 'C++' sep 'experiments' sep 'structuralGP_OLD' sep 'output'];
else
    address = last_address;
end

directory_name = uigetdir(address,'Select the experiment');
dim_dir = size(directory_name);
length = dim_dir(2);

% extract the name of the experiment
pos = strfind(directory_name, filesep);  %find the last occurrence of /
start = max(pos)+1;
experiment = directory_name(start:length);