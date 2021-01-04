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


%%% assign_f_obj_fun.m
%%% function to assign the objective function of the GP search (if known)
%%% the function is chosen among a list

function [f_obj] = assign_f_obj_fun()

% cell array containing the strings of the functions
%fun_cell = cell(7,1);

% OBJECTIVE FUNCTION LIST
% 2 INPUT VARIABLES
% 1. Rosenbrock  -2<Z1<2 -2<Z2<2
C{1,1} = '100*((Z2-Z1.^2).^2)+(1-Z1).^2';
% 2. RatPol2D    0<Z1<6  0<Z2<6   (pag. 118 Katia's thesis)
C{2,1} =  '((Z1-3.0).^4.0 + (Z2-3.0).^3.0 - (Z2-3.0))./((Z2-2.0).^4.0 + 10.0)'; 
% 3. Branin-Hoo  -5<Z1<10   0<Z2<15  (Haftka paper)
C{3,1} = '(Z2 - ((5.1*(Z1.^2))./(4*(3.141592653589793).^2)) + ((5.*Z1)./(3.141592653589793)) - 6).^2 + 10.*(1 - 1./(8*3.141592653589793)).*cos(Z1) + 10';
% 4. First Alvarez function (pag. 66 thesis)   0<Z1<5   0<Z2<5
C{4,1} = '(30 + Z1.*sin(Z1)).*(4+exp(-Z2))';
% 5. Kotanchek 0<Z1<4   0<Z2<4  (pag. 118 Katia's thesis)
C{5,1}= 'exp(-(Z1-1).^2)  ./  (1.2 + (Z2 - 2.5).^2 )';
% 6. Kotanchekplus1 0<Z1<4   0<Z2<4
C{6,1}= '1.0 + (exp(-(Z1-1).^2)  ./  (1.2 + (Z2 - 2.5).^2 ))';
% 7. RatPol2D + 50  0<Z1<6  0<Z2<6
C{7,1} =  '50.0 + ((Z1-3).^4 + (Z2-3).^3 - (Z2-3))./((Z2-2).^4 + 10)'; 
% 8. Betts function 4<Z1<10  4<Z2<10
C{8,1} = '((Z2.^3)./(27*sqrt(3.0))).*(((Z1-3.0).^2)-9.0)';
% 9. Peyman function 1<Z1<400 1<Z2<400
C{9,1} = '(9.5.*Z1.^2)./(Z1.^2 + 25000)-(0.01.*Z1)./(Z1./25000 + Z2./50 + 1) - 0.01.*Z1 + 0.05';


% 1 INPUT VARIABLE
% 10. Z1*sinZ1  0<=Z1<=pi
C{10,1} = 'Z1.*sin(Z1)';
% 11. Z1*sin(5Z1)  0<=Z1<=pi
C{11,1} = 'Z1.*sin(5.*Z1)';
% 12. exp(-Z1)*(Z1^3)*cos(Z1)*sin(Z1)*( cos(Z1)((sin(Z1))^2) )-1 
C{12,1} = 'exp(-Z1).*(Z1.^3).*cos(Z1).*sin(Z1).*( cos(Z1).*((sin(Z1)).^2) -1)';

% 4 INPUT VARIABLES
% 13. 0.81+(24.3*(((2.0*Z1)+(3.0*(Z2)^2))/((4.0*(Z3)^3)+(5.0*(Z4)^4))))
C{13,1} = '0.81+(24.3.*(((2.0.*Z1)+(3.0.*((Z2).^2)))./((4.0.*((Z3).^3))+(5.0.*((Z4).^4)))))';

% NO FUNCTION
% 99. ""
C{99,1} = '';
 

% show the functions
disp('Select the objective function:')
disp('2 INPUT VARIABLES: ')
disp(' 1. Rosenbrock')
disp(' 2. RatPol2D')
disp(' 3. Branin-Hoo')
disp(' 4. Alvarez function')
disp(' 5. Kotanchek')
disp(' 6. Kotanchek PLUS 1')
disp(' 7. RatPol2D + 50')
disp(' 8. Betts')
disp(' 9. Peyman function')
disp('1 INPUT VARIABLE:')
disp(' 10. Z1*sin(Z1)')
disp(' 11. Z1*sin(5*Z1)')
disp(' 12. Salustowicz1D')
disp('4 INPUT VARIABLES:')
disp(' 13. 0.81+(24.3*(((2.0*Z1)+(3.0*(Z2)^2))/((4.0*(Z3)^3)+(5.0*(Z4)^4))))')
disp('-------------------')
disp(' 99. NO FUNCTION TO BE EVALUATED')

s = 'Select a number : ';
n_fun = input(s);

% assign the selected function
f_obj = char(C{n_fun,1});
