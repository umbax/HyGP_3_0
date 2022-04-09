% function that look for a file 
% input: directory where to search, file name (strings)
% output: 1 if the file is found, 0 otherwise

function [found] = search_file(directory, filename)

% files in the directory
file_list = dir(directory);
dim_file_list = size(file_list);
found = 0;

for k=3:dim_file_list(1)    %start from 3 because the first entries are . and ..
    name_file = file_list(k).name;
    if (strcmp(name_file,filename))
        found = 1;
        break
    end
end

%found

