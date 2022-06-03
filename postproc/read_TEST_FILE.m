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


%%% read_TEST_FILE.m
%%% script to read file with test data set for HyGP
%%% input: directory of experiment (may be empty), no of variables of the problem NVAR
%%% output: no. of variables NVAR_test, no of test records NTESTCASES, test record matrix Rtest   


function [NVAR_test, NTESTCASES, Rtest, found_test]=read_TEST_FILE(directory_name, NVAR)
% READ FILE CONTAINING TEST DATA SET
%oldFolder = cd(directory_name);
% select a file with TEST data
s = sprintf('\n');
disp(s);

% SELECT FILE WITH TEST DATA SET
disp('Select the file containing the TEST data set...');
%first_path = ''
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
        disp('Reading test data set...')
        NVAR_test=C{1}
        NTESTCASES=C{2}
        N = cast(NTESTCASES,'double');
        
%        
%         if (NVAR_test~=NVAR)
%             disp('ERROR! NVAR different from NVAR_test:');
%             disp('the number of variables in building and test data set are different. Exit');
%             NVAR_test=-1;
%             NTESTCASES=-1;
%             Rtest=NaN;
%             found_test=0;
%             return
%         end
    
       
        % read TEST data
        Rtest=zeros(NTESTCASES, NVAR_test+1);
        for k=1:NTESTCASES
            for j=1:NVAR_test
                C = textscan(fid_test,'%f',1,'commentstyle','#');
                Rtest(k,j)=C{1};
            end
            C = textscan(fid_test,'%f %*[^\n]',1,'commentstyle','#');
            Rtest(k,NVAR_test+1)=C{1};

        end
        % close stream
        fclose(fid_test);

        disp('First line of test data set')
        disp(Rtest(1,:))

    else
        % no test data file found
        disp('No test data set selected. Exit')
        NVAR_test=-1;
        NTESTCASES=-1;
        Rtest=NaN;
        found_test=0;
    end
else
   % no test data file specified
   disp('No test data set specified. Exit')
   NVAR_test=-1;
   NTESTCASES=-1;
   Rtest=NaN;
   found_test=0; 
end

end
