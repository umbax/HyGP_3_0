function [epurated_mean,nrun_valid] = smart_mean(A)
    % this function computes the average along each row of A, neglecting error
    % entries (-1)
    dim_A = size(A);
    n_rows = dim_A(1);
    n_columns = dim_A(2);
    n_runvalid=zeros(n_rows);
    total= 0.;
    n=0;
    for i=1:n_rows
        for j=1:n_columns
            if (A(i,j)~=-1)
                total = total + A(i,j);
                n=n+1;
            end
            nrun_valid(i)= n; 
        end
        if (n)
            epurated_mean(i,1) = total/n;
            total = 0.;
            n=0;
        end

    end
end