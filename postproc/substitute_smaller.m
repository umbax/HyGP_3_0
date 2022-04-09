% function to substitute values bigger than max_value with nan in a given matrix 
% Nan will be then removed from single columns using remove_nan.m
% input: matrix or vector, max_value
% output: another (amended) matrix or vector

function [Wclean] = substitute_smaller(W, min_value)

nrow = size(W,1);
ncol = size(W,2);
Wclean = zeros(nrow, ncol);

for k=1:nrow
    for j=1:ncol
        if (W(k,j)<=min_value)
            Wclean(k,j)= nan;
        else
            Wclean(k,j)= W(k,j);
        end
    end
end















end