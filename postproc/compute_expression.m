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


%%% compute_expression.m
%%% Script to evaluate a HyGP expression recursively to overcome the limit on the number of nested brackets.
%%% Limit of 32 in Matlab. For example when evaluating the expression with eval(): 
%%% "(-2.08915e-01 + ((((-1.14654e+02 * (Z1 * Z1)) + (((1.75968e+02 * (cos((-1.88441e+02 * (sin((7.53433e+01 * Z1))))))) * (((1.93809e+02 + (-1.52844e+02 * (sin((-4.01539e+00 * (sin((-5.57157e+01 * Z1)))))))))^1)) + (1.60423e+02 * (Z1 * Z1)))) + ((((-1.39354e-05 * (9.57839e+01 + ((Z1 * Z1) * Z1))))^1) + ((1.08022e+02 * Z1) + (((((((-3.94578e+01 + (9.66158e+01 * Z1)))^1) + ((-1.19353e+02 * Z1) + (-1.88264e+02 * Z1))) + (((-2.03441e+01 * (cos((7.01351e+01 * (sin((-1.61505e+02 * Z1))))))) * (((-1.20601e+02 + (5.60993e+01 * (sin((((4.71973e+01 * Z1) + (-1.24871e+01 * Z1)) + (((-8.61147e+00 + (((((-1.68367e+02 + (8.23163e+01 * (cos((-1.75453e+02 * Z1))))))^1) * (-1.87124e+02 * (cos(((-1.97540e+01 * Z1) + (((-1.24248e+02 + ((((-3.32365e+01 * Z1) + (-9.32993e+01 * Z1)) + (((-7.74957e+01 + (2.67666e+01 * (sin((1.32025e+02 * Z1))))))^1)) + (((1.04254e+02 + (-7.51206e+01 * Z1)))^1))))^1)))))) + (((((-1.10380e+02 * Z1) + (2.10961e+02 * Z1)) * ((-4.43406e+01 * Z1) + (5.72305e+01 * Z1))) + ((1.90894e+02 * (sin((1.36962e+02 * Z1)))) + (-9.26591e+01 * Z1))) + (-1.13912e+02 * (Z1 * Z1))))))^1)))))))^1)) + (-6.61853e+01 * (Z1 * Z1)))) + ((((-1.07870e+02 * (-8.26056e+01 + (Z1 * Z1))))^1) + ((-7.99889e+01 * Z1) + (1.00126e+02 * Z1)))) + (((1.28345e+02 * (1.87209e+02 + (Z1 * Z1))))^1))))) + (1.30254e+02 * (cos((-7.25318e+01 * (Z1 * Z1)))))))"
%%% Error: Nesting of {, [, and ( cannot exceed a depth of 32.
%%% Sin and cos take arguments in radians, as in HyGP

% function that computes the value of a HyGP model (text expression) 
% ATTENTION: to improve precision of the result, increase the number of decimal digits in formatSpec='%.20e'
function [value] = compute_expression(expr, x)
%disp('Start evaluation cycle')
%expr
% define precision
formatSpec='%.20e'; % to match HyGp format

nvar=size(x,2);

% assign input values to variables Z1, Z2, ecc
if (nvar==1)
    Z1=x;    
end
if (nvar==2)
    Z1=x(1);
    Z2=x(2);
end
if (nvar>2)
    % 31/12/2020 to be completed (the id number of the variables to be updated is needed)
    f = msgbox('compute_expression.m : No. of variables >2 : function still not implemented!','ATTENTION!');
end

% find the boundaries of the hierarchical "blocks" that make up the model expression
l_expr = strlength(expr);
map=zeros(2,l_expr);  % row 1 contains the index of the bracket position, row 2 the bracket counter 
bracket_count=1;
map(:,1)=[1;bracket_count];
j=1;
map_size=1;
for i=2:l_expr
    % block start
    if ((expr(i)=='(')) %  && ((expr(i-1)==" ") || (expr(i-1)=="(")))
        bracket_count=bracket_count+1;
        j=j+1;
        map(:,j)=[i; bracket_count];
        map_size=map_size+1;
    end
    % block end
    if ((expr(i)==')')) % && ((expr(i+1)==" ") || (expr(i+1)==")"))) % attenzione! questa condizione "salva" la parentesi di chiusura di sin e cos, quando in realtà non dovrebbe!
        bracket_count=bracket_count-1;
        j=j+1;
        map(:,j)=[i; bracket_count];
        map_size=map_size+1;
    end
    if (map(2,j)==0)
        % reached last closing bracket of the expression
        break    
    end
end

%map_size
map=map(:,1:map_size);
max_bracket=max(map(2,:));

% evaluation completed if the max depth of a bracket is 1
if (max_bracket<=32) %(max_bracket==1) evaluate until the deepest brackets block, slow! -- (max_bracket<=32) % evaluate directly with eval() if nesting of brackets has a depth lower than 32 (Matlab limit)
    value=double(eval(expr));
    return
end

% populate list of brackets at next inner level
replace_map=zeros(1,map_size);  % index in map corresponding the the max depth of the bracket (size might be too big, but it is more efficient like this)
j=0;
for i=2:map_size
    %i
    % count opening bracket
    if (j)
        if ((map(2,i)==2) && (map(2,replace_map(j))<2))
            j=j+1;
            replace_map(j) = i;  % opening bracket
        end
    else
        if (map(2,i)==2)
            j=j+1;
            replace_map(j) = i;  % opening bracket
        end
    end
    % count closing bracket
    if (map(2,i)==1)
        j=j+1;
        replace_map(j) = i;    % corresponding closing bracket
    end
end
%replace_map
replace_map_size=j;
replace_map=replace_map(1:replace_map_size);

% evaluate recursively the blocks one level deeper than current bracket block level
reduced_expr = expr(1:map(1, replace_map(1))-1);  
for k=1:2:replace_map_size
        %k
        % evaluate block at replace_map(k)
        s_to_evaluate = expr(map(1, replace_map(k)): map(1, replace_map(k+1))); 
        s_to_evaluate_value = compute_expression(s_to_evaluate, x);
        % replace value in the original expression you need to evaluate 
        % (not possible to use a vector as independent variable as a single
        % value is needed...)
        reduced_expr=strcat(reduced_expr, ' (', num2str(s_to_evaluate_value, formatSpec), ')');
        
        % add also the operations among blocks (present if replace_map_size>2)
        %disp('Add the joining section of expr....')
        %k
        %replace_map_size
        if ((replace_map_size>=4) && (k<replace_map_size-2))
            reduced_expr=strcat(reduced_expr, expr(map(1, replace_map(k+1))+1: map(1, replace_map(k+2))-1)); % mind that strcat removed the trailing space in char arrays  
        end
        %disp('End joining section of expr....')
end     
reduced_expr=strcat(reduced_expr,expr(map(1, replace_map(replace_map_size))+1:l_expr)); % mind that strcat removed the trailing space in char arrays  


%disp('Computed value:')
value = double(eval(reduced_expr));

%disp('End evaluation cycle')
return
end

