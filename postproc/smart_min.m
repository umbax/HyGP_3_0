function [epurated_min, row_min, column_min,last_gen_ord] = smart_min(A)
% this function computes the minimum along each ROW of A, neglecting error
% entries (-1). If all columns contain -1, the computation stops.
% So the vector "epurated_min" contains for each generation
% the minimum in each row of A.
% It also returns the corresponding position:
% (row_min, column_min)
    
    dim_A = size(A);
    n_rows = dim_A(1);
    n_columns = dim_A(2);
    epurated_min = -1;
    r_min = -1;
    c_min = -1;
    minimum=-1;
    n=0;
    last_gen = -1*ones(2, n_columns); % run and the corresponding last_generation
    for i=1:n_rows
        for j=1:n_columns
            if (A(i,j)~= -1)
                if (minimum==-1)
                    minimum=A(i,j);
                    r_min = i;
                    c_min = j;
                end
                if (A(i,j)<minimum);
                    minimum = A(i,j);
                    r_min = i;
                    c_min = j;              
                end
                n=n+1;
            end
            
            % update last_gen 
            if ((A(i,j)== -1) && (i>1))
                if (A(i-1,j)~= -1)
                    last_gen(1,j) = j;
                    last_gen(2,j) = i-2;
                end
            end
            if ((A(i,j)~=-1) && (i==n_rows))
                last_gen(1,j) = j;
                last_gen(2,j) = i-1;
            end
        end
        if (n)
            % update minimum value on the row
            %epurated_min(i,1) = minimum;
            epurated_min = [epurated_min ; minimum];
                
            % update position of global minimum value 
            if (i==1)
                global_min = minimum;
                row_min = r_min;
                column_min = c_min;
            else    
                if (epurated_min(i,1) <= global_min) 
                    global_min = epurated_min(i,1);
                    row_min = r_min;
                    column_min = c_min;
                end
            end    
            n=0;
            minimum=-1;
        else
            % all column contain -1, so the computation stops
            break
        end
    end
    % remove the first element, -1
    %epurated_min
    dim_epurated_min = size(epurated_min);
    %disp(['Size epurated_min BEFORE ' num2str(size(epurated_min))])
    %epurated_min
    epurated_min = epurated_min(2:dim_epurated_min(1));
    %disp(['Size epurated_min AFTER ' num2str(size(epurated_min))])
    row_min = row_min-1;   % this is REALLY IMPORTANT, due to the shift in epurated_min
    last_gen_ord= sort(last_gen(2,:));

