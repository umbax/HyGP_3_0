% function to remove NaN entries from a vector
% ONLY VECTOR!

function [xclean] = remove_nan(x) 

dimx = size(x);
d1=max(dimx);
d2=min(dimx);


if (d2>1)
    disp('remove_nan: the input is not a vector!!')
    disp('The input is returned with no changes')
    xclean = x;
    return 
end
ok=0;
for k=1:d1
    
    if (~strcmpi(num2str(x(k)),'NaN'))
        if (~ok)
           xclean = x(k); 
           ok=1;
        else
            xclean = [xclean x(k)];
        end
    end 
end

end
