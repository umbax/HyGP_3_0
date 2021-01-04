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


%%% count_runs.m
%%% function to count the number of runs an experiment is composed of

function [last_run] = count_runs(directory_name) 

file_list = dir(directory_name);
dim_file_list = size(file_list);
n_run = 0;
last_run = 0;
for k=3:dim_file_list(1)    %start from 3 because the first entries are . and ..
    name_file = file_list(k).name;
    if (name_file(1:3)=='run')
        length_name_file = size(name_file);
        n_run = name_file(5:length_name_file(2));
        num_run = str2num(n_run);
        if (num_run>last_run)
            last_run = num_run;
        end 
    end
end